import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import atrail
import bench
import clip
import euler_transform
import foldtp
import origami_slicer
import slice_stl


class PipelineTests(unittest.TestCase):
    def assert_single_trail(self, path):
        report = foldtp.FoldTP(path).validate('test')
        self.assertTrue(report.ok, report.errors)
        self.assertTrue(report.is_eulerian_complete)
        self.assertTrue(report.is_a_trail)
        self.assertTrue(report.is_single_stroke)
        metrics = bench.toolpath_metrics(path)
        self.assertEqual((metrics['strokes'], metrics['travel_moves'], metrics['self_crossings']), (1, 0, 0))

    def test_valid_fixture(self):
        self.assert_single_trail(json.loads((ROOT/'examples/square_tour.fold').read_text()))

    def test_crossing_fixture_is_detected(self):
        raw = json.loads((ROOT/'examples/figure8_crossing.fold').read_text())
        self.assertFalse(foldtp.FoldTP(raw).validate('crossing').is_a_trail)

    def test_clipped_miura_square_and_star(self):
        outlines = [ [[.3,.3],[9.7,.3],[9.7,9.7],[.3,9.7]],
                     [[5+4.6*math.cos(math.radians(90+k*36))*(1 if k%2==0 else .42),
                       5+4.6*math.sin(math.radians(90+k*36))*(1 if k%2==0 else .42)] for k in range(10)]]
        for outline in outlines:
            with self.subTest(vertices=len(outline)):
                pattern = bench.miura(16, 16, a=1, b=1, c=.36)
                for vertex in pattern['vertices_coords']:
                    vertex[0] -= 3; vertex[1] -= 3
                clipped, _ = clip.clip_to_polygon(pattern, outline)
                path, problems = atrail.ATrailGenerator(clipped).generate(patch=True)
                self.assertFalse(problems, problems)
                self.assert_single_trail(path)

    def test_euler_baseline(self):
        path, problems, base_edges = euler_transform.gupta_toolpath(euler_transform.hex_mesh(0,0,10,10,2))
        self.assertFalse(problems, problems)
        self.assertEqual(len(path['edges_vertices']), 4*base_edges)
        self.assert_single_trail(path)

    def test_slicer_cli_outside_repository(self):
        with tempfile.TemporaryDirectory() as tmp:
            mesh = Path(tmp)/'gear.stl'; output = Path(tmp)/'gear.gcode'
            slice_stl.gen_gear_stl(str(mesh), teeth=10, R=12, root=8.4, h=6)
            result = subprocess.run([sys.executable, str(ROOT/'origami_slicer.py'), str(mesh),
                                     '-o', str(output), '--fit-mm', '24', '--spacing', '1.5',
                                     '--max-layers', '2', '--walls', '0'], cwd=tmp,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            stats = json.loads(result.stdout.split('\nWrote ')[0])
            self.assertEqual(stats['layers'], 2)
            self.assertEqual(stats['layer_strokes'], [1,1])
            self.assertEqual(stats['failed_layers'], 0)
            self.assertIn('M83', output.read_text())
            self.assertIn('G1 X', output.read_text())

    def test_cli_rejects_invalid_parameters(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)/'out.gcode'
            result = subprocess.run([sys.executable, str(ROOT/'origami_slicer.py'), 'missing.stl',
                                     '-o', str(output), '--spacing', '0'], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            self.assertIn('positive finite', result.stderr)
            self.assertFalse(output.exists())

    def test_browser_bundle_matches_sources(self):
        result = subprocess.run([sys.executable, str(ROOT/'build_web.py')], cwd=tempfile.gettempdir(),
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        for name in ('atrail.py','clip.py','foldtp.py','slice_stl.py','origami_slicer.py'):
            self.assertEqual((ROOT/name).read_bytes(), (ROOT/'web/py'/name).read_bytes())


if __name__ == '__main__':
    unittest.main()
