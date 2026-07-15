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
    # On a headless Xvfb display the window comes up 1x1; maximize it so VNC/noVNC
    # shows a usable GUI rather than a blank desktop.
    try:
        win = viewer.window._qt_window
        win.resize(1600, 950)
        win.showMaximized()
    except Exception as exc:  # pragma: no cover - display quirks
        print(f"[napari-mcp] could not maximize window: {exc}", flush=True)
    server = NapariBridgeServer(viewer, port=PORT)
    if server.start():
        print(f"[napari-mcp] bridge up -> http://localhost:{PORT}/mcp", flush=True)
    else:
        print("[napari-mcp] WARNING: bridge failed to start", flush=True)
    napari.run()


if __name__ == "__main__":
    main()
