"""Segmentation: locked Li mitochondria mask, nuclei, and cell territories.

`li` is the locked mitochondrial mask (bake-off winner). The MVP threshold
comparison lives in `threshold_masks` for provenance / sensitivity checks.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage
from skimage import feature, filters, measure, morphology, segmentation

from .config import ExperimentConfig
from .io import area_px, radius_px, to_uint8

# remove_small_objects/holes: skimage 0.26 replaced min_size/area_threshold with
# max_size. Verified empirically that max_size=N reproduces the old min_size=N
# exactly, preserving the frozen MVP behaviour (min 32 px). Kept fixed by decision.
_MIN_OBJECT_PX = 32

# --- Cellpose whole-cell footprint (Module A) ---------------------------------
# cellprob < 0 grows masks into dimmer edge mito (fixes under-coverage of big
# cells); flow > default keeps lower-quality masks. MITO_GAMMA < 1 brightens dim
# peripheral mito so Cellpose recalls faint edges / low-signal cells.
CELLPROB_THRESHOLD = -2.0
FLOW_THRESHOLD = 0.7
MITO_GAMMA = 0.5

# LAMBDA (in EDT-pixel units): cost of crossing one non-mito pixel relative to one
# pixel of Euclidean distance. High enough that any mito path beats a dark-gap
# crossing, so watershed ridges fall in the mito-empty gaps *between* cells.
_ELEVATION_LAMBDA = 50.0


def preprocess_li(mip: np.ndarray, px_um: float) -> np.ndarray:
    """Match the MVP Li baseline with a physical top-hat radius."""
    img = filters.gaussian(mip, sigma=1.0, preserve_range=True)
    return morphology.white_tophat(to_uint8(img), morphology.disk(radius_px(1.5, px_um)))


def postprocess_li(mask: np.ndarray) -> np.ndarray:
    """Match the MVP Li baseline postprocessing."""
    mask = morphology.remove_small_objects(mask.astype(bool), max_size=_MIN_OBJECT_PX)
    return morphology.remove_small_holes(mask, max_size=_MIN_OBJECT_PX)


def pipeline_li(mip: np.ndarray, px_um: float) -> tuple[np.ndarray, float]:
    """Locked Li mitochondrial mask; returns (mask, li_threshold)."""
    pre = preprocess_li(mip, px_um)
    threshold = float(filters.threshold_li(pre))
    return postprocess_li(pre > threshold), threshold


def threshold_masks(mip: np.ndarray) -> dict[str, tuple[np.ndarray, float]]:
    """MVP threshold comparison (otsu / li / yen / adaptive) for sensitivity."""
    img = filters.gaussian(mip, sigma=1.0, preserve_range=True)
    img = morphology.white_tophat(img, morphology.disk(15))
    thresholds = {
        "otsu": filters.threshold_otsu(img),
        "li": filters.threshold_li(img),
        "yen": filters.threshold_yen(img),
    }
    masks = {name: (img > value, float(value)) for name, value in thresholds.items()}
    local = filters.threshold_local(img, block_size=81, method="gaussian")
    masks["adaptive"] = (img > local, float(np.mean(local)))
    return {name: (postprocess_li(m), v) for name, (m, v) in masks.items()}


def segment_nuclei(
    nucleus_mip: np.ndarray, cfg: ExperimentConfig, px_um: float
) -> tuple[np.ndarray, np.ndarray]:
    """Segment nuclei with Otsu thresholding and distance watershed splitting."""
    smooth = filters.gaussian(nucleus_mip, sigma=1.0, preserve_range=True)
    if float(smooth.max()) <= float(smooth.min()):
        nucleus_mask = np.zeros(smooth.shape, dtype=bool)
    else:
        nucleus_mask = smooth > filters.threshold_otsu(smooth)
    nucleus_mask = ndimage.binary_fill_holes(nucleus_mask)
    nucleus_mask = morphology.remove_small_objects(
        nucleus_mask,
        max_size=area_px(cfg.min_nucleus_area_um2, px_um),
    )

    dist = ndimage.distance_transform_edt(nucleus_mask)
    seed_dist = filters.gaussian(dist, sigma=2.0, preserve_range=True)
    coords = feature.peak_local_max(
        seed_dist,
        labels=nucleus_mask,
        min_distance=radius_px(cfg.nucleus_seed_min_distance_um, px_um),
        exclude_border=False,
    )
    markers = np.zeros(nucleus_mask.shape, dtype=np.int32)
    if coords.size:
        markers[tuple(coords.T)] = np.arange(1, len(coords) + 1)
    elif nucleus_mask.any():
        markers = measure.label(nucleus_mask).astype(np.int32)

    nuclei = segmentation.watershed(-dist, markers, mask=nucleus_mask)
    return nuclei.astype(np.int32), nucleus_mask.astype(bool)


def build_territories(
    nucleus_labels: np.ndarray, cfg: ExperimentConfig, px_um: float
) -> np.ndarray:
    """LEGACY / DEPRECATED — radius-capped Euclidean Voronoi from nuclei.

    Splits purely by nearest nucleus, so distal mito of elongated/stellate cells
    bleed into a neighbouring nucleus and adjacent cells cannot be separated.
    Superseded by ``build_territories_v2`` (footprint-bounded, marker-controlled
    watershed). Kept only as a fallback when Cellpose is unavailable; do not
    remove (backtracking / provenance)."""
    nucleus_mask = nucleus_labels > 0
    # isotropic_dilation is EDT-based (O(N), exact for a disk); a grey dilation with
    # this ~170 px radius footprint is O(N*K) and blows time/memory on full frames.
    foreground = morphology.isotropic_dilation(
        nucleus_mask, radius_px(cfg.max_cell_radius_um, px_um)
    )
    elevation = ndimage.distance_transform_edt(~nucleus_mask)
    territories = segmentation.watershed(elevation, markers=nucleus_labels, mask=foreground)
    return territories.astype(np.int32)


def _enhance_mito(mito_mip: np.ndarray) -> np.ndarray:
    """Gamma-brighten dim mito so Cellpose recalls faint peripheries / cells."""
    x = to_uint8(mito_mip).astype(np.float32) / 255.0
    return (np.power(x, MITO_GAMMA) * 255.0).astype(np.float32)


def cell_bodies(
    mito_mip: np.ndarray,
    nucleus_mip: np.ndarray,
    cfg: ExperimentConfig,
    px_um: float,
    gpu: bool = False,
) -> np.ndarray:
    """Cellpose instance labels for whole cell bodies from the mito channel.

    Runs ``cyto3`` on a 2-channel image [cytoplasm=enhanced mito, nucleus] and
    returns the INSTANCE label array (one id per detected body) — Cellpose's
    *learned* cell-cell boundaries are the primary segmentation, kept intact by
    ``resolve_instances`` for the 1-nucleus majority. Cellpose is imported lazily
    (the env may lack it at import time)."""
    from cellpose import models

    diameter = 2.0 * radius_px(cfg.max_cell_radius_um, px_um)  # cell ~ 2x radius cap
    rgb = np.zeros((*mito_mip.shape, 3), dtype=np.float32)
    rgb[..., 0] = _enhance_mito(mito_mip)
    rgb[..., 1] = to_uint8(nucleus_mip)
    model = models.Cellpose(gpu=gpu, model_type="cyto3")
    masks, _, _, _ = model.eval(
        rgb, diameter=diameter, channels=[1, 2],
        cellprob_threshold=CELLPROB_THRESHOLD, flow_threshold=FLOW_THRESHOLD,
    )
    return masks.astype(np.int32)


def resolve_instances(
    nucleus_labels: np.ndarray,
    cp_bodies: np.ndarray,
    mito_mask: np.ndarray,
    cfg: ExperimentConfig,
    px_um: float,
) -> np.ndarray:
    """Reconcile Cellpose instance bodies with the nucleus seeds.

    Cellpose's learned boundary is trusted; we intervene only where the nucleus
    count disagrees, so a clean single cell is never carved up by a neighbour's
    basin (the v1 "basin bleed" regression). Per Cellpose body:

    - **1 nucleus**  -> keep the whole body, labelled by that nucleus id.
    - **>=2 nuclei** -> genuine merge: split THIS body only, by a nucleus-seeded
      watershed confined to the body (geodesic elevation: cheap through mito,
      costly across dark gaps).
    - **0 nuclei**   -> drop (unnucleated debris; not a cell).

    A nucleus covered by no body (Cellpose miss) gets a radius-capped disk over
    still-unclaimed pixels, so no nucleus is silently dropped.

    Territories are labelled by nucleus id (cell_id == nucleus_id)."""
    nucleus_mask = nucleus_labels > 0
    out = np.zeros_like(nucleus_labels, dtype=np.int32)
    claimed: set[int] = set()

    # Shared geodesic elevation, reused for any multi-nucleus body split.
    elevation = ndimage.distance_transform_edt(~nucleus_mask).astype(np.float32)
    elevation += _ELEVATION_LAMBDA * (~mito_mask & ~nucleus_mask).astype(np.float32)

    for body_id in range(1, int(cp_bodies.max()) + 1):
        body = cp_bodies == body_id
        if not body.any():
            continue
        nids = np.unique(nucleus_labels[body])
        nids = nids[nids > 0]
        if len(nids) == 0:
            continue  # unnucleated debris -> drop
        if len(nids) == 1:
            out[body] = nids[0]  # trust Cellpose boundary as-is (no bleed)
            claimed.add(int(nids[0]))
            continue
        # >=2 nuclei: genuine merge -> split THIS body only, seeded by its nuclei.
        markers = np.where(body, nucleus_labels, 0)
        sub = segmentation.watershed(elevation, markers=markers, mask=body)
        out[body] = sub[body]
        claimed.update(int(n) for n in nids)

    # Fallback: nuclei Cellpose missed entirely -> radius-capped disk.
    radius = radius_px(cfg.max_cell_radius_um, px_um)
    for nid in range(1, int(nucleus_labels.max()) + 1):
        if nid in claimed:
            continue
        nuc = nucleus_labels == nid
        if not nuc.any():
            continue
        claim = morphology.isotropic_dilation(nuc, radius) & (out == 0)
        out[claim] = nid
    return out


def build_territories_v2(
    mito_mip: np.ndarray,
    nucleus_mip: np.ndarray,
    nucleus_labels: np.ndarray,
    mito_mask: np.ndarray,
    cfg: ExperimentConfig,
    px_um: float,
    gpu: bool = False,
) -> np.ndarray:
    """New Module-A entrypoint: Cellpose bodies -> nucleus-reconciled instances."""
    bodies = cell_bodies(mito_mip, nucleus_mip, cfg, px_um, gpu=gpu)
    return resolve_instances(nucleus_labels, bodies, mito_mask, cfg, px_um)
