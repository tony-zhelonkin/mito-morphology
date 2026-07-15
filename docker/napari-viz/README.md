# mito-viz — containerized interactive viewer

A **disposable Docker image** bundling the whole display stack so it never touches the host
or the analysis pixi env:

- **Xvfb + software GL (mesa llvmpipe)** — a virtual display for napari
- **napari + napari-mcp bridge** (`:9999`) — Claude Code and the notebook drive it
- **marimo** (`:2720`) — reactive notebook / report (edit mode, write your own cells)
- **x11vnc** (`:5901`, optional) — view the *real* napari GUI in a VNC client

The container runs with `--network host`, so all three land on the host's **localhost**
(reachable by Claude Code directly, and by you over an SSH tunnel). Nothing binds a public
interface.

## Build (once)
```bash
docker build -t mito-viz:latest .
```

## Run against an experiment repo
```bash
./run.sh /path/to/TBK1i-ISD90-mito-morphology
# ENABLE_VNC=0 ./run.sh ...   # skip the VNC server
```
The repo is mounted at `/repo`; the notebook is `<repo>/02_analysis/explore.py`.

## Connect from your laptop
```bash
ssh -L 2720:localhost:2720 -L 5901:localhost:5901 you@remote
# browser   -> http://localhost:2720      (marimo — reactive report)
# VNC client-> localhost:5901             (live napari GUI)
```

## Notes
- **One napari-mcp call at a time** — the bridge is not concurrency-safe.
- Software GL is used (robust, no GPU setup). The box has an NVIDIA GPU; to use it, add
  `--gpus all` + EGL and drop `LIBGL_ALWAYS_SOFTWARE` — not needed for this workload.
- marimo runs `--no-token` on localhost only; reach it strictly through the SSH tunnel.
- The image pins `mito-morphology` from git, matching the compute package.
