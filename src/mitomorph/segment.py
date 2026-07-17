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

# Fallback noise floor for the provenance-only `threshold_masks` tool. The run
# path (`pipeline_li`) instead uses a PHYSICAL floor (cfg.mito_min_object_um2)
# scaled by px_um, so the smallest retained object is the same physical size
# across batches with different pixel sizes.
_MIN_OBJECT_PX_FALLBACK = 32

# Cellpose (cell_bodies), watershed (resolve_instances), and the mito noise floor
# are now config-driven knobs on ExperimentConfig (run-config.yaml), not module
# constants — so a run's parameters are tracked and reproducible.


def preprocess_li(mip: np.ndarray, px_um: float) -> np.ndarray:
    """Match the MVP Li baseline with a physical top-hat radius."""
    img = filters.gaussian(mip, sigma=1.0, preserve_range=True)
    return morphology.white_tophat(to_uint8(img), morphology.disk(radius_px(1.5, px_um)))


def postprocess_li(mask: np.ndarray, min_object_px: int = _MIN_OBJECT_PX_FALLBACK) -> np.ndarray:
    """Match the MVP Li baseline postprocessing (min-object floor in pixels)."""
    mask = morphology.remove_small_objects(mask.astype(bool), max_size=min_object_px)
    return morphology.remove_small_holes(mask, max_size=min_object_px)


def pipeline_li(
    mip: np.ndarray, px_um: float, cfg: ExperimentConfig | None = None
) -> tuple[np.ndarray, float]:
    """Locked Li mitochondrial mask; returns (mask, li_threshold).

    The noise floor is physical (cfg.mito_min_object_um2 -> px via px_um); falls
    back to the frozen 32 px when no cfg is supplied (e.g. provenance tooling)."""
    pre = preprocess_li(mip, px_um)
    threshold = float(filters.threshold_li(pre))
    min_px = (
        area_px(cfg.mito_min_object_um2, px_um) if cfg is not None
        else _MIN_OBJECT_PX_FALLBACK
    )
    return postprocess_li(pre > threshold, min_px), threshold


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


def _enhance_mito(mito_mip: np.ndarray, gamma: float) -> np.ndarray:
    """Gamma-brighten dim mito so Cellpose recalls faint peripheries / cells."""
    x = to_uint8(mito_mip).astype(np.float32) / 255.0
    return (np.power(x, gamma) * 255.0).astype(np.float32)


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
    rgb[..., 0] = _enhance_mito(mito_mip, cfg.mito_gamma)
    rgb[..., 1] = to_uint8(nucleus_mip)
    model = models.Cellpose(gpu=gpu, model_type=cfg.cellpose_model_type)
    masks, _, _, _ = model.eval(
        rgb, diameter=diameter, channels=[1, 2],
        cellprob_threshold=cfg.cellprob_threshold, flow_threshold=cfg.flow_threshold,
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
    elevation += cfg.elevation_lambda * (~mito_mask & ~nucleus_mask).astype(np.float32)

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
