# Concern 05 — Is Nellie's over-segmentation attributable to voxel anisotropy?

**Question posed:** Is "Nellie prefers isotropic voxels / the ~7.5–14:1 Z-anisotropy of the project's
confocal stacks caused the severe over-segmentation" grounded in anything real, or is it our
speculation? And does Nellie *internally* handle anisotropy (read Z-res, adjust the Frangi
sigma) — or not?

**Method:** read Nellie's local source (`init-draft/01_modules/nellie/`), the *Nature Methods*
2025 paper + PMC full text, and the GitHub repo/docs. Every claim below is bucketed as
**DOCUMENTED** (Nellie's own source/paper says so, with citation), **INFERRED-STRONG**
(well-grounded imaging/sampling principle), or **SPECULATION** (our guess, no source).

---

## Bottom line (read this first)

1. **Nellie DOES internally handle anisotropy.** It reads per-axis physical spacing from
   metadata and builds an *anisotropic* Frangi sigma vector (`sigma_z = sigma / z_ratio`) and
   computes the Hessian with *physical* per-axis spacing. This is DOCUMENTED in both source and
   paper. So the naive form of the claim — "Nellie assumes isotropic voxels / ignores the
   Z-step" — is **FALSE / unsupported.** Nellie explicitly does *not* require isotropic voxels
   and does *not* resample to isotropic.

2. **BUT the over-segmentation is still attributable to anisotropy** — just via a different,
   deeper mechanism than "Nellie forgot about Z." Nellie's anisotropy handling is
   *mathematically correct* (correct physical units) but **correctness cannot manufacture
   information the Z-axis never sampled.** At 7.5–14:1 anisotropy, tubules thinner than the
   Z-step occupy a single Z-plane, so the 3D vesselness + 3D connected-component labelling has
   no Z-neighbourhood to bridge them across slices → fragmentation. This is INFERRED-STRONG
   (a sampling/optics argument), and it is directly corroborated by the observed behaviour.

**Verdict: PARTLY GROUNDED, with the emphasis corrected.** Anisotropy is very likely the
dominant driver of the fragmentation, but *not* because Nellie ignores it — rather because the
Z-axis is undersampled relative to the structures, and no amount of correct sigma-scaling fixes
undersampling. "Nellie prefers isotropic voxels" is true in the weak sense (isotropic, well-Z-
sampled data is where a 3D vesselness pipeline works best) but false in the strong sense (Nellie
does not internally assume or enforce isotropy).

---

## DOCUMENTED — Nellie reads per-axis spacing and adapts the filter

**(a) It reads per-axis physical spacing (X, Y, Z) from image metadata.**
`nellie/im_info/verifier.py` populates a `dim_res = {'X','Y','Z','T'}` dict from ImageJ, OME,
raw TIFF tags, or ND2 metadata:
- OME: `dim_res['Z'] = metadata.images[0].pixels.physical_size_z` (verifier.py:243)
- ImageJ: `dim_res['Z'] = metadata['spacing']` (verifier.py:229)
- The extraction script for this dataset writes `spacing = 0.635 µm` into the ImageJ TIFF, and
  §11 of the init-draft session log confirms Nellie reads it back correctly.

**(b) It builds an ANISOTROPIC Frangi sigma vector — Z is scaled by the anisotropy ratio.**
`nellie/segmentation/filtering.py`:
- `z_ratio = z_res / x_res` (filtering.py:78)
- `min_radius_px = min_radius_um / dim_res['X']` — radii converted µm→px using XY (filtering.py:88-89)
- `_get_sigma_vec`: for 3D, `sigma_vec = (sigma / z_ratio, sigma, sigma)` (filtering.py:285).
  i.e. the Z Gaussian sigma is *divided* by the anisotropy ratio so that the physical smoothing
  scale is the same on all three axes. This is genuine anisotropy compensation, not an
  isotropic-pixel assumption.

**(c) The Hessian is computed in PHYSICAL units (per-axis spacing).**
`_get_spacing(ndim)` returns `(z_res, y_res, x_res)` in µm (filtering.py:265-275), and the
Hessian derivatives use it: `self.xp.gradient(image, *spacing)` (filtering.py:518). So the
second-derivative ridge measure is physically scaled per axis.

**(d) The min-object-size filter is a PHYSICAL volume, anisotropy-aware.**
`labelling.py:_compute_min_area_pixels` (line 209-219): `volume_um3 = 4/3·π·r³`, then
`volume_px = volume_um3 / (x_res·y_res·z_res)` — uses the true anisotropic voxel volume.

**(e) The same `z_ratio` scaling is reused downstream** in mocap marking
(`mocap_marking.py:125,297,326`) and skeleton/networking (`networking.py:73,460`). Anisotropy
awareness is consistent across the whole pipeline, not just the filter.

**(f) The PAPER states this explicitly** (Nature Methods 2025, PMC11978511):
> "Our filter … **automatically adjusts the filter's effective range based on voxel dimensions
> to adapt to various magnifications and anisotropies.**"
> "Our pipeline contrasts with the current SOTA … which are not adaptive to intrinsic image
> metadata … [and] use the same filter parameters across all scales."

So adapting to anisotropy is an advertised design feature. Nellie is presented as *not*
requiring the user to resample to isotropic.

**(g) There is NO resample-to-isotropic code path.** A repo-wide grep for
`resample|isotrop|anisotrop|rescale|zoom(` finds only comments — no upsampling of Z to isotropic
anywhere. Nellie works in native anisotropic voxel space. (The only literal "anisotropic"
mention is a comment in `mocap_marking.py:443` noting KD-tree uses pixel units.)

**What Nellie does NOT document:** a stated minimum Z-sampling, a maximum tolerable anisotropy
ratio, or any warning that sparse-Z confocal data over-segments. The paper/README give no
numeric anisotropy limit (confirmed by web search). The GitHub docs stress only that metadata
must be read correctly (fields turn green), implying "as long as the Z-step is correctly
entered, you're fine" — which is the part our evidence qualifies below.

---

## INFERRED-STRONG — why correct anisotropy handling still fragments this data

Correct physical-unit scaling does not create Z-resolution that was never acquired. Concrete
numbers for the init-draft stack (XY = 0.0852 µm, Z = 0.635 µm, `z_ratio` = 7.45), Nellie
defaults (`min_radius_um=0.25`, `max_radius_um=1.0`):

- Frangi scales (XY, in px): `sigma_min = min_radius_px/2 = 1.47`, `sigma_max = max_radius_px/3 = 3.91`
  (filtering.py:297-300).
- Their **Z components** are `sigma/z_ratio`: **0.20 px to 0.53 px** across *every* scale.
  → The Z Gaussian is **sub-pixel at all scales**; the Hessian's Z second-derivative is taken
  over a ~1-slice neighbourhood. There is effectively no scale-space in Z.
- A tubule of radius 0.25 µm has diameter 0.5 µm = **0.79 Z-pixels** — it lives in a *single*
  z-plane. A 3D vesselness filter needs a coherent ridge across the plane normal; here there is
  no across-Z sample to form one, and 3D connected-component labelling cannot link a structure
  that appears in only one slice to its continuation in the next.

This is the classic, well-grounded failure mode of Frangi/Hessian vesselness + 3D connected
components on **under-sampled (sub-Nyquist in Z) anisotropic stacks**: thin tubules "fall
between slices," 3D connectivity breaks, and the object count explodes into fragments. It is a
sampling-theory consequence (Nyquist), not a bug in Nellie's scaling. Grounding:
- Frangi et al., *Multiscale vessel enhancement filtering*, MICCAI 1998 — vesselness is defined
  on the local Hessian; its discrimination degrades when a scale cannot be resolved on the grid.
- General confocal-anisotropy literature (web search surfaced this as a "recognized failure
  mode": fragmentation of 3D structures on anisotropic/sparse-Z data; e.g. anisotropic neuron-
  tracing and 2D-stack-consensus 3D segmentation papers exist specifically to work *around* it).

The two current experiment repos are **worse** on exactly this axis: TBK1i ~9:1
(XY 0.0706 / Z 0.633) and glutamine ~14:1 (XY 0.0659 / Z 0.923) → the Z sigma components shrink
further and tubule-per-slice occupancy drops, so we should expect *equal or worse* fragmentation
out of the box.

This inference is corroborated by the empirical session-log results (init-draft
`session_log.md`): the only levers that reduced object count were widening the **XY Frangi scale
range** (Preset C, −28%) and **post-hoc morphological closing** (−60% more) — i.e. bridging the
broken 3D connectivity after the fact. `Label.min_radius_um` alone did essentially nothing
(§9). That pattern is exactly what an under-sampled-Z connectivity failure predicts.

---

## SPECULATION (flagged honestly)

- That resampling the stacks to isotropic before Nellie would *fix* it — plausible but unproven;
  Z-upsampling only interpolates, it cannot recover unsampled tubule structure, so it may reduce
  fragmentation cosmetically without recovering the thinnest tubules. Untested here.
- The exact contribution split between anisotropy vs. the other documented difficulties (dim
  8-bit data, dense/"beefy" perinuclear mitochondria with Z-overlap — user note §12) is not
  quantified. Anisotropy is likely dominant for the *peripheral tubular* network; perinuclear
  clumping is a separate problem.
- Whether a 2D-per-slice segmentation + across-Z stitching would outperform Nellie's native 3D
  path on this anisotropy — a reasonable architectural hypothesis, but ours, not sourced.

---

## Answers to the two direct questions

- **Does Nellie internally handle anisotropy (read Z-res, adjust sigma)?** **YES, documented.**
  It reads `dim_res['Z']`, scales the Frangi sigma in Z by `1/z_ratio`, computes the Hessian in
  physical units, and sizes the min-object filter by true voxel volume. It does *not* resample
  to isotropic and does *not* assume isotropic pixels.
- **Is "anisotropy caused the over-segmentation" grounded?** **PARTLY — and the framing must be
  corrected.** Not because Nellie ignores Z (it doesn't), but because at 7.5–14:1 the Z-axis is
  sub-Nyquist for these tubules, so 3D ridge detection + 3D connectivity fragment regardless of
  correct scaling. The fix is not "tell Nellie about the Z-step" (it already knows) but to
  address undersampled-Z connectivity (XY-scale widening, morphological reconnection, or a
  2D-stack/anisotropy-aware strategy).

---

### Citations
- Local source: `init-draft/01_modules/nellie/nellie/segmentation/filtering.py:78,88-89,265-275,285,297-300,518`;
  `.../labelling.py:96,209-219`; `.../im_info/verifier.py:227-244`; `.../mocap_marking.py:125,297,326`;
  `.../networking.py:73,460`.
- Session evidence: `init-draft/02_analysis/nellie_mcp_log/session_log.md` §6-12.
- Paper: Lefebvre et al., *Nellie…*, Nature Methods 2025, doi:10.1038/s41592-025-02612-7 —
  https://www.nature.com/articles/s41592-025-02612-7 ; full text https://pmc.ncbi.nlm.nih.gov/articles/PMC11978511/
- Repo/docs: https://github.com/aelefebv/nellie
- Method grounding: Frangi et al., *Multiscale vessel enhancement filtering*, MICCAI 1998.
