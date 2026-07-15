#!/usr/bin/env bash
# Launch the mito-viz interactive stack against one experiment repo.
#   ./run.sh /abs/path/to/<experiment-repo>
# Then, from your laptop:
#   ssh -L 2720:localhost:2720 -L 5901:localhost:5901 you@remote
#   browser -> http://localhost:2720   (marimo)
#   VNC client -> localhost:5901       (live napari GUI)
set -euo pipefail

REPO="$(realpath "${1:?usage: run.sh <experiment-repo-path>}")"
[ -f "$REPO/experiment.toml" ] || { echo "not an experiment repo (no experiment.toml): $REPO" >&2; exit 1; }

exec docker run --rm -it \
  --name mito-viz \
  --network host \
  -e MITO_REPO=/repo \
  -e ENABLE_VNC="${ENABLE_VNC:-1}" \
  -v "$REPO":/repo \
  mito-viz:latest
