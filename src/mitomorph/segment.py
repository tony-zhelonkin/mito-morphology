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
    """Build nucleus-derived, radius-capped nearest-nucleus territories."""
    nucleus_mask = nucleus_labels > 0
    # isotropic_dilation is EDT-based (O(N), exact for a disk); a grey dilation with
    # this ~170 px radius footprint is O(N*K) and blows time/memory on full frames.
    foreground = morphology.isotropic_dilation(
        nucleus_mask, radius_px(cfg.max_cell_radius_um, px_um)
    )
    elevation = ndimage.distance_transform_edt(~nucleus_mask)
    territories = segmentation.watershed(elevation, markers=nucleus_labels, mask=foreground)
    return territories.astype(np.int32)
