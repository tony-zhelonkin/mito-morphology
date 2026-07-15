#!/usr/bin/env python3
"""Launch napari with the napari-mcp bridge (inside the viz container).

The bridge binds 127.0.0.1:9999; the container runs with --network host, so that
is the host's localhost:9999 — reachable by Claude Code and by the marimo notebook.
Only ONE napari-mcp call at a time (the bridge is not concurrency-safe).
"""
import napari
from napari_mcp.bridge_server import NapariBridgeServer

PORT = 9999


def main() -> None:
    viewer = napari.Viewer(title=f"mito-viz — MCP bridge :{PORT}")
    server = NapariBridgeServer(viewer, port=PORT)
    if server.start():
        print(f"[napari-mcp] bridge up -> http://localhost:{PORT}/mcp", flush=True)
    else:
        print("[napari-mcp] WARNING: bridge failed to start", flush=True)
    napari.run()


if __name__ == "__main__":
    main()
