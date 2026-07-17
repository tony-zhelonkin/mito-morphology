"""OME-Zarr (OME-NGFF 0.4) substrate: open locked Zeiss .czi in napari/Fiji,
and preserve segmentation masks as versioned, provenanced label groups.

Why NGFF **0.4** (not 0.5): 0.4 is the widely-shipped, stable spec that today's
napari (napari-ome-zarr), Fiji (MoBIE / bio-formats), BigDataViewer and vizarr all
read without fuss. It targets Zarr **v2** on disk, which the current tool ecosystem
expects; 0.5 (Zarr v3) is newer than the software our wet-lab collaborators run, so
we deliberately pin 0.4 for interoperability. The raw Zeiss .czi is locked/opaque,
so converting it once into this open substrate is what makes the data browsable.

The seam
--------
This module is the shared on-disk contract between the two segmentation modules::

    Module A (segmentation) --write_labels--> labels/<version>/ <--read_labels-- Module B (analysis)

`czi_to_zarr` materialises a raw image pyramid once. Module A writes territory /
instance masks into ``labels/<version>/`` with provenance; Module B (and napari,
Fiji) reads whichever version it wants. Label versions are additive and never
clobbered, so retraining or expert curation adds ``1_retrained`` / ``2_expert_final``
alongside ``0_baseline`` instead of destroying the earlier mask.

``ome_zarr`` and ``zarr`` are imported *lazily* inside each function: they live in
the GPU env and may be absent where this module is merely imported (CLI help path,
the CPU compute env), so importing this module must never require them.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .config import ExperimentConfig
from .io import read_manifest, resolve_samples

NGFF_VERSION = "0.4"

# Multiscale pyramid: 4 levels, isotropic 2x XY downsample per level.
_N_SCALES = 4
_DOWNSCALE = 2


def czi_to_zarr(
    czi_path: Path | str,
    out_zarr_path: Path | str,
    px_um: float,
    channel_names: list[str] | None = None,
) -> Path:
    """Read a Zeiss ``.czi`` with bioio and write an OME-NGFF 0.4 multiscale image.

    Channel axes are preserved and the physical XY pixel size (``px_um``, microns)
    is written into the ``coordinateTransformations`` so napari/Fiji show real-world
    scale. A 4-level, 2x XY pyramid is written for smooth multi-resolution viewing.

    Parameters
    ----------
    czi_path : source CZI.
    out_zarr_path : destination ``.zarr`` store (created / overwritten).
    px_um : physical XY pixel size in microns (from the manifest).
    channel_names : optional channel labels; falls back to the CZI's own names.

    Returns the ``Path`` to the written store.
    """
    from bioio import BioImage
    from ome_zarr.io import parse_url
    from ome_zarr.scale import Scaler
    from ome_zarr.writer import write_image

    czi_path = Path(czi_path)
    out_zarr_path = Path(out_zarr_path)
    out_zarr_path.parent.mkdir(parents=True, exist_ok=True)

    img = BioImage(czi_path)
    data = img.get_image_data("CZYX")  # (C, Z, Y, X)
    if channel_names is None:
        try:
            channel_names = list(img.channel_names)
        except Exception:
            channel_names = []

    axes = [
        {"name": "c", "type": "channel"},
        {"name": "z", "type": "space", "unit": "micrometer"},
        {"name": "y", "type": "space", "unit": "micrometer"},
        {"name": "x", "type": "space", "unit": "micrometer"},
    ]
    # One scale transform per pyramid level; XY scale grows with the 2x downsample.
    coordinate_transformations = [
        [{"type": "scale",
          "scale": [1.0, 1.0, px_um * (_DOWNSCALE**lvl), px_um * (_DOWNSCALE**lvl)]}]
        for lvl in range(_N_SCALES)
    ]

    import zarr

    store = parse_url(str(out_zarr_path), mode="w").store
    root = zarr.group(store=store)
    write_image(
        image=data,
        group=root,
        axes=axes,
        coordinate_transformations=coordinate_transformations,
        scaler=Scaler(downscale=_DOWNSCALE, max_layer=_N_SCALES - 1, method="nearest"),
        storage_options={"chunks": (1, 1, data.shape[-2], data.shape[-1])},
    )
    if channel_names:
        root.attrs["omero"] = {
            "version": NGFF_VERSION,
            "channels": [{"label": name, "active": True} for name in channel_names],
        }
    return out_zarr_path


def write_labels(
    zarr_path: Path | str,
    labels: np.ndarray,
    name: str,
    provenance: dict,
) -> str:
    """Write an instance-label array into ``labels/<name>/`` with provenance.

    ``provenance`` (e.g. ``{"algorithm":..., "version":..., "params":..., "source":...}``)
    is stored in the label group's ``.zattrs`` under ``"mitomorph_provenance"`` so the
    mask is self-describing and auditable.

    Overwrite policy: label versions are **immutable**. If ``labels/<name>`` already
    exists this raises ``FileExistsError`` rather than silently clobbering a mask a
    downstream module may depend on — callers pick a fresh version name
    (``0_baseline`` -> ``1_retrained`` -> ``2_expert_final``). Returns ``name``.
    """
    from ome_zarr.io import parse_url
    from ome_zarr.writer import write_labels as _write_labels

    zarr_path = Path(zarr_path)
    if (zarr_path / "labels" / name).exists():
        raise FileExistsError(
            f"label version {name!r} already exists in {zarr_path}; "
            "label versions are immutable — use a new version name"
        )

    import zarr

    labels = np.asarray(labels)
    store = parse_url(str(zarr_path), mode="a").store
    root = zarr.group(store=store)
    _write_labels(labels, group=root, name=name, axes=_label_axes(labels.ndim))
    root["labels"][name].attrs["mitomorph_provenance"] = dict(provenance)
    return name


def read_labels(zarr_path: Path | str, name: str) -> np.ndarray:
    """Read the full-resolution label array from ``labels/<name>/``."""
    import zarr
    from ome_zarr.io import parse_url

    root = zarr.group(store=parse_url(str(Path(zarr_path)), mode="r").store)
    # Level "0" is full resolution within the label multiscale group.
    return np.asarray(root["labels"][name]["0"])


def list_label_versions(zarr_path: Path | str) -> list[str]:
    """List label version names present under ``labels/`` (sorted)."""
    labels_dir = Path(zarr_path) / "labels"
    if not labels_dir.is_dir():
        return []
    return sorted(p.name for p in labels_dir.iterdir() if p.is_dir())


def run(cfg: ExperimentConfig, samples: list[str] | None) -> None:
    """Batch ``.czi`` -> OME-Zarr conversion over requested samples (or all raw CZI).

    Writes one store per sample at ``03_results/00_ome_zarr/<group>/<sample>.zarr``,
    the canonical browsable substrate segmentation label groups are added into.
    """
    manifest_rows = read_manifest(cfg)
    for filestem, czi in resolve_samples(cfg, samples):
        if filestem in cfg.skip_samples:
            print(f"WARNING: skipping {cfg.experiment}/{filestem} (config skip_samples)")
            continue
        if not czi.exists():
            print(f"ERROR: missing raw CZI: {czi}")
            continue
        group, _ = cfg.parse_sample(filestem)
        out = cfg.results / "00_ome_zarr" / group / f"{filestem}.zarr"
        try:
            px_um = float(manifest_rows[(cfg.experiment, czi.name)]["px_x_um"])
            print(f"converting {cfg.experiment}/{czi.name} -> {out} ...", flush=True)
            czi_to_zarr(czi, out, px_um)
        except (KeyError, ValueError) as exc:
            print(f"ERROR: {cfg.experiment}/{filestem}: {exc}")


def _label_axes(ndim: int) -> list[dict]:
    """Spatial axes for a label array (2D yx or 3D zyx)."""
    space = [
        {"name": "z", "type": "space", "unit": "micrometer"},
        {"name": "y", "type": "space", "unit": "micrometer"},
        {"name": "x", "type": "space", "unit": "micrometer"},
    ]
    return space[-ndim:]
