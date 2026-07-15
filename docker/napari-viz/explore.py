"""mito-viz — reactive exploration notebook (marimo).

Runs inside the mito-viz container against the mounted experiment repo (cwd = /repo).
Left/right dropdowns compare two samples' cell-QC overlays + per-cell tables; a button
pushes the selected sample into the live napari viewer via the napari-mcp bridge.

Edit freely — marimo is reactive, so changing a cell re-runs everything downstream.
"""
import marimo

app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    from pathlib import Path
    import pandas as pd
    from mitomorph.config import load_config

    root = Path.cwd()
    cfg = load_config(root)
    cellqc = cfg.results / "04_cellqc"
    return mo, pd, Path, cfg, cellqc


@app.cell
def _(mo, cfg):
    mo.md(f"# {cfg.experiment} — mito-morphology\nInteractive cell-QC review.")
    return


@app.cell
def _(cellqc):
    # discover samples that have a rendered cell-QC overlay
    samples = sorted(
        p.name for p in cellqc.glob("*/*")
        if (p / "figures").glob("*_cellqc_overlay.png")
        and any((p / "figures").glob("*_cellqc_overlay.png"))
    )
    return (samples,)


@app.cell
def _(mo, samples):
    left = mo.ui.dropdown(samples, value=samples[0] if samples else None, label="Left sample")
    right = mo.ui.dropdown(
        samples, value=samples[1] if len(samples) > 1 else (samples[0] if samples else None),
        label="Right sample",
    )
    mo.hstack([left, right], justify="start", gap=2)
    return left, right


@app.cell
def _(mo, cellqc):
    def overlay_path(sample):
        hits = list((cellqc).glob(f"*/{sample}/figures/{sample}_cellqc_overlay.png"))
        return hits[0] if hits else None

    def panel(sample):
        p = overlay_path(sample)
        if p is None:
            return mo.md(f"**{sample}**: no overlay (run `mito cellqc-viz`)")
        return mo.vstack([mo.md(f"**{sample}**"), mo.image(str(p), width=520)])

    return panel, overlay_path


@app.cell
def _(mo, panel, left, right):
    mo.hstack([panel(left.value), panel(right.value)], widths="equal", gap=2)
    return


@app.cell
def _(mo, pd, cellqc, left, right):
    def cells_table(sample):
        hits = list(cellqc.glob(f"*/{sample}/tables/{sample}_cells.csv"))
        if not hits:
            return mo.md(f"_{sample}: no table_")
        df = pd.read_csv(hits[0])
        cols = ["cell_id", "mito_area_um2", "skeleton_length_um", "branch_count",
                "mito_solidity", "mito_form_factor", "focus_score", "qc_pass"]
        return mo.ui.table(df[cols], selection=None, page_size=8, label=sample)

    mo.hstack([cells_table(left.value), cells_table(right.value)], widths="equal", gap=2)
    return


@app.cell
def _(mo):
    push = mo.ui.button(label="Push LEFT sample to napari")
    mo.md("### Live napari (via napari-mcp bridge on localhost:9999)")
    push
    return (push,)


@app.cell
async def _(mo, cfg, left, push):
    # Push the selected sample's mito MIP + nuclei + territories into the live viewer.
    mo.stop(not push.value, mo.md("_Click the button to load the left sample into napari._"))
    import numpy as np, tifffile
    from fastmcp import Client

    sample = left.value
    group = cfg.parse_sample(sample)[0]
    mip = cfg.results / "02_extract" / group / sample / "figures" / f"{sample}_mito_mip.tif"
    masks = cfg.results / "04_cellqc" / group / sample / "masks"

    code = f'''
import numpy as np, tifffile
viewer.layers.clear()
viewer.add_image(tifffile.imread(r"{mip}"), name="{sample} mito", colormap="gray")
for _n in ("nuclei", "territories"):
    _p = r"{masks}/" + _n + ".npy"
    try:
        viewer.add_labels(np.load(_p), name="{sample} " + _n)
    except Exception as _e:
        print("skip", _n, _e)
viewer.reset_view()
'''
    try:
        async with Client("http://localhost:9999/mcp") as client:
            tools = [t.name for t in await client.list_tools()]
            tool = "execute_code" if "execute_code" in tools else tools[0]
            await client.call_tool(tool, {"code": code})
        status = mo.md(f"✅ pushed **{sample}** to napari (tool: `{tool}`)")
    except Exception as e:
        status = mo.md(f"⚠️ bridge not reachable: {e}")
    status
    return


if __name__ == "__main__":
    app.run()
