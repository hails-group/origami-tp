#!/usr/bin/env python3
"""Local static server for the slicer web app.

Slicing now runs entirely in the browser, so this is just a
static file server for local development. Run `python3 build_web.py` first to populate web/py/.

  python3 server.py        # -> http://localhost:8765
"""
import http.server, functools, sys
from pathlib import Path

WEB = Path(__file__).resolve().parent / "web"
PORT = 8765

if __name__ == "__main__":
    if not (WEB / "py" / "origami_slicer.py").exists():
        sys.exit("web/py is empty — run:  python3 build_web.py")
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(WEB))
    print(f"origami slicer (static, in-browser compute) -> http://localhost:{PORT}   (Ctrl-C to stop)")
    try:
        http.server.ThreadingHTTPServer(("127.0.0.1", PORT), handler).serve_forever()
    except KeyboardInterrupt:
        print("\nbye")
