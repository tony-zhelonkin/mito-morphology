# mito-viz — disposable interactive display stack (Xvfb + software GL + napari +
# napari-mcp bridge + marimo + VNC). Nothing here touches the host or the analysis
# pixi env; an experiment repo is mounted at /repo and its 03_results are viewed.
FROM condaforge/miniforge3:latest

# System X / GL / VNC stack (kept inside the image, never on the host).
RUN apt-get update && apt-get install -y --no-install-recommends \
      xvfb x11vnc fluxbox tini procps x11-utils \
      libgl1-mesa-dri libglu1-mesa libegl1 mesa-utils \
      libglib2.0-0 libdbus-1-3 libxkbcommon-x11-0 libxcb-cursor0 \
      libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-randr0 \
      libxcb-render-util0 libxcb-shape0 libxcb-xinerama0 libxcb-util1 \
    && rm -rf /var/lib/apt/lists/*

COPY environment.yml /tmp/environment.yml
RUN mamba env create --yes -f /tmp/environment.yml && mamba clean -afy

# noVNC static files (browser VNC client); served by websockify at runtime.
ADD https://github.com/novnc/noVNC/archive/refs/tags/v1.5.0.tar.gz /tmp/novnc.tgz
RUN mkdir -p /opt/viz/novnc && tar xzf /tmp/novnc.tgz -C /opt/viz/novnc --strip-components=1 && rm /tmp/novnc.tgz

# Make the viz env the default interpreter.
ENV PATH=/opt/conda/envs/viz/bin:$PATH \
    CONDA_DEFAULT_ENV=viz \
    DISPLAY=:99 \
    LIBGL_ALWAYS_SOFTWARE=1 \
    GALLIUM_DRIVER=llvmpipe \
    QT_QPA_PLATFORM=xcb \
    MPLBACKEND=Agg

COPY napari_bridge.py /opt/viz/napari_bridge.py
COPY entrypoint.sh /opt/viz/entrypoint.sh
RUN chmod +x /opt/viz/entrypoint.sh

# --network host is used at runtime, so EXPOSE is documentation only.
EXPOSE 9999 2720 5901 6080
ENTRYPOINT ["tini", "--", "/opt/viz/entrypoint.sh"]
