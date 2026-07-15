"""Per-cell mitochondrial metrics: skeleton network, shape, and focus."""
from __future__ import annotations

import math

import numpy as np
from scipy import ndimage
from skimage import measure, morphology

FOCUS_EPS = 1.0e-6


def compute_skeleton_metrics(
    skeleton: np.ndarray, binary: np.ndarray, px_um: float
) -> dict[str, float | int]:
    """Compute skeleton network metrics from a binary mask."""
    total_length_um = float(skeleton.sum()) * px_um

    kernel = np.ones((3, 3), dtype=np.uint8)
    kernel[1, 1] = 0
    neighbor_count = ndimage.convolve(
        skeleton.astype(np.uint8), kernel, mode="constant", cval=0
    )
    neighbor_count = neighbor_count * skeleton

    junctions = neighbor_count >= 3
    endpoints = neighbor_count == 1
    junction_count = int(measure.label(junctions).max())
    endpoint_count = int(endpoints.sum())

    branch_skel = skeleton & ~morphology.dilation(junctions, morphology.disk(1))
    branch_labels = measure.label(branch_skel)
    branch_count = int(branch_labels.max())

    if branch_count > 0:
        branch_lengths = [
            float((branch_labels == label).sum()) * px_um
            for label in range(1, branch_count + 1)
        ]
        mean_branch_length_um = float(np.mean(branch_lengths))
    else:
        mean_branch_length_um = 0.0

    dist = ndimage.distance_transform_edt(binary)
    if skeleton.any():
        mean_diameter_um = float(dist[skeleton].mean()) * 2.0 * px_um
    else:
        mean_diameter_um = 0.0

    return {
        "total_length_um": round(total_length_um, 1),
        "branch_count": branch_count,
        "junction_count": junction_count,
        "endpoint_count": endpoint_count,
        "mean_branch_length_um": round(mean_branch_length_um, 3),
        "mean_diameter_um": round(mean_diameter_um, 3),
    }


def mito_shape_metrics(mito_cell: np.ndarray) -> tuple[float, float]:
    """Return area-weighted solidity and form factor for a cell's mito mask."""
    labels = measure.label(mito_cell)
    props = measure.regionprops(labels)
    if not props:
        return 0.0, 0.0

    total_area = float(sum(prop.area for prop in props))
    solidity = sum(prop.area * prop.solidity for prop in props) / total_area
    form_factor_sum = 0.0
    for prop in props:
        if prop.perimeter > 0:
            form_factor = 4.0 * math.pi * prop.area / (prop.perimeter**2)
        else:
            form_factor = 0.0
        form_factor_sum += prop.area * form_factor
    return float(solidity), float(form_factor_sum / total_area)


def focus_score(nucleus_mip: np.ndarray, laplace: np.ndarray, territory: np.ndarray) -> float:
    """DAPI-based focus: variance of Laplacian normalized by mean intensity^2."""
    lap_vals = laplace[territory]
    nuc_vals = nucleus_mip[territory].astype(float)
    return float(np.var(lap_vals) / (float(nuc_vals.mean()) ** 2 + FOCUS_EPS))
