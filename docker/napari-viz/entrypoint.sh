#!/usr/bin/env bash
# Start the full interactive display stack inside the container:
#   Xvfb (virtual display) + software GL -> napari + MCP bridge (:9999)
#   x11vnc (:5901, optional) -> view the real napari GUI
#   marimo (:2720) -> reactive report / notebook
# Container runs with --network host, so all of these are on the host's localhost.
set -euo pipefail

REPO=${MITO_REPO:-/repo}
NOTEBOOK=${MITO_NOTEBOOK:-$REPO/02_analysis/explore.py}
MARIMO_PORT=${MARIMO_PORT:-2720}
VNC_PORT=${VNC_PORT:-5901}
NOVNC_PORT=${NOVNC_PORT:-6080}

cleanup() { pkill -P $$ 2>/dev/null || true; }
trap cleanup EXIT

# 1. virtual display + software GL (mesa llvmpipe)
Xvfb :99 -screen 0 1920x1080x24 -nolisten tcp &
sleep 1
export DISPLAY=:99
fluxbox >/dev/null 2>&1 &

# 2. VNC of the real napari window (localhost only; reach via SSH tunnel)
if [ "${ENABLE_VNC:-1}" = "1" ]; then
  x11vnc -display :99 -localhost -forever -shared -nopw -quiet -rfbport "$VNC_PORT" &
  echo "[viz] VNC on localhost:$VNC_PORT" >&2
  # noVNC: browser-based VNC client (localhost only; reach via SSH tunnel).
  websockify --web=/opt/viz/novnc "127.0.0.1:$NOVNC_PORT" "localhost:$VNC_PORT" >/tmp/novnc.log 2>&1 &
  echo "[viz] noVNC on http://localhost:$NOVNC_PORT/vnc.html" >&2
fi

# 3. napari + napari-mcp bridge (binds 127.0.0.1:9999)
python /opt/viz/napari_bridge.py &
echo "[viz] napari-mcp bridge on localhost:9999" >&2

# 4. marimo edit (write reactive code; localhost only)
echo "[viz] marimo on localhost:$MARIMO_PORT  (notebook: $NOTEBOOK)" >&2
cd "$REPO"
exec marimo edit --headless --host 127.0.0.1 --port "$MARIMO_PORT" --no-token "$NOTEBOOK"
