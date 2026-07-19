# The data we work with, its limits, and why 2D-MIP (not 3D)

**Date:** 2026-07-18 · **Status:** design rationale (reference section for the docs).
**Evidence:** `concern-05-nellie-anisotropy-evidence.md` (Nellie source + Nat Methods paper),
the init-draft Nellie MVP (`init-draft/02_analysis/nellie_mcp_log/session_log.md`), and the
per-file pixel-size metadata in each experiment's `03_results/00_inventory/tables/czi_manifest.csv`.

---

## 1. Here's the data we were working with

Confocal z-stacks of BMDCs on a Zeiss system (core facility), 2 channels
(mitochondria + nucleus), 8-bit, 63× objective. The defining property is the
**voxel geometry** — how big one 3D pixel is along each axis:

| Dataset | XY (µm/px) | Z (µm/slice) | **Z : XY anisotropy** | z-slices |
|---|---|---|---|---|
| init-draft (ISD903) | 0.0852 | 0.635 | **~7.5 : 1** | 16 |
| **TBK1i** experiment | 0.0706 | 0.633 | **~9 : 1** | 7–16 |
| **glutamine** experiment | 0.0659 | 0.923 | **~14 : 1** | 14–19 |

The XY sampling is fine (arguably oversampled). The **Z sampling is coarse**, and
it got *coarser* across experiments — the glutamine set steps ~0.92 µm between
slices. That per-session variation in z-step (and one field acquired at coarser XY)
is the fingerprint of **per-session operator settings at a shared scope**, not a
fixed protocol — which is exactly why the pipeline carries pixel size from each
file's metadata rather than a hard-coded constant. *(In init-draft that constant
was a real bug: `PIXEL_SIZE_UM = 0.1` vs the true 0.0852.)*

## 2. Isotropic vs anisotropic — the one concept that governs everything

**Isotropic** = a voxel is the same physical size in X, Y, and Z (a cube).
**Anisotropic** = the voxel is a tall, thin box. Ours are boxes ~7.5–14× deeper
than they are wide.

```
ANISOTROPIC (our confocal stacks)        ISOTROPIC (ideal for 3D analysis)
   XY 0.07 µm                                 all axes equal
  ┌──┐                                       ┌──────┐
  │  │  Z 0.63–0.92 µm  (9–14× taller)       │      │   a cube — 3D neighbourhoods
  │  │                                       │      │   are trustworthy in every axis
  └──┘                                       └──────┘
  a tall box
```

Two causes, different in kind:

- **Optics sets an unavoidable floor (~3:1).** A confocal's resolution is
  intrinsically worse in Z: lateral ≈ λ/(2·NA), but axial ≈ 2·λ·n/(NA²). The
  axial term scales with **NA²**, so even a great objective is ~2.5–3.5× coarser
  in Z. Some anisotropy is baked in by diffraction — not a Zeiss quirk; every
  confocal (Leica/Nikon/Zeiss) has it.
- **The rest is an acquisition choice (the operator's z-step).** To *resolve* a
  ~0.6 µm optical Z-resolution by Nyquist you'd sample Z at **~0.2–0.3 µm**. At
  0.63–0.92 µm the stack is **undersampled in Z by ~3–4×** on top of the optical
  floor — the two multiply to the 9:1 / 14:1 we see. Coarse z-steps are a
  *rational* choice at the scope: faster acquisition, less photobleaching (crucial
  for a dim AF488 mito stain), smaller files, and often the goal was a qualitative
  "is it fragmented" overview, not quantitative 3D reconstruction.

## 3. Here are the limitations of this kind of data

1. **There is no trustworthy 3D structure to recover.** A mitochondrial tubule is
   ~0.3 µm thick — that is **<1 z-slice** tall at our sampling. Thin tubules fall
   *between* slices, so their true 3D connectivity was never captured. This is a
   **sub-Nyquist-Z** limit: it is a property of the acquisition, and **no software
   — Nellie included — can interpolate back information the microscope never
   sampled.** Isotropic resampling does not fix it (interpolation invents, it does
   not recover).
2. **3D object-connectivity fragments.** Any 3D method (vesselness/Frangi + 3D
   connected components) has no across-Z neighbourhood to bridge a tubule that
   lives in one plane, so continuous networks shatter into many small objects.
3. **Between-session variability.** The z-step (and occasionally XY) differs per
   session, so absolute 3D measurements would not even be comparable across images
   without care. Metadata-carried pixel size + relative between-condition contrasts
   sidestep this.

## 4. Why we went 2D-MIP, rather than "3D-Nellie this out"

We tested the 3D route honestly in init-draft. On the 7.5:1 data, **Nellie
over-segmented severely — 2,364 organelle fragments on a single image**, tubules
shattered, mean branch length collapsing to fragment-scale. Tuning the Frangi
scale range helped only partially (−28%), and only *post-hoc morphological closing*
brought the count down ~60%. The current experiments are 9:1 and 14:1 — worse.

**Crucially, this is not a Nellie defect** (a myth worth killing explicitly):

- **Nellie reads per-axis physical spacing and builds an *anisotropic* Frangi
  sigma** (`filtering.py:78,285`): `sigma_z = sigma / z_ratio`, `z_ratio =
  z_res/x_res`; the Hessian is computed in physical units per-axis; the min-object
  filter uses true anisotropic voxel volume. The Nat Methods paper states the
  filter "automatically adjusts based on voxel dimensions to adapt to various
  magnifications and anisotropies." There is **no** hidden "assume isotropic"
  assumption (grep for any resample-to-isotropic path comes back empty).
- So the fragmentation is **not** "Nellie forgot the Z-step." It is that at 7.5:1+,
  the Z sigma component is **sub-pixel at every Frangi scale** and a tubule occupies
  <1 z-plane — correct sigma-scaling cannot manufacture connectivity across planes
  that were never sampled. It is the **data**, not the tool.

Given that, **collapsing to a maximum-intensity projection (MIP) is the honest
response, not a workaround for a tool flaw**: we are not discarding good 3D
information — there was no trustworthy 3D information to discard. In 2D the tubule
is contiguous, connectivity is real, and classical adaptive thresholding (Li here;
Sauvola/Phansalkar explored) is tunable, inspectable, and scalable. The init-draft
classical route *did* initially over-segment too, but that was a **parameter**
problem we corrected (v2 tuning); the 3D route's fragmentation was **structural**.

## 5. Evidence honesty (documented / inferred / speculation)

Per `concern-05`, kept explicit so this doc doesn't overclaim:

- **Documented:** Nellie reads z-spacing and scales sigma anisotropically; no
  isotropy assumption; the paper claims anisotropy adaptation. *(source + paper)*
- **Inferred-strong:** sub-Nyquist-Z undersampling is the dominant driver of the
  fragmentation — a standard sampling/optics principle, corroborated by the MVP log
  (only XY-scale widening + morphological closing helped; the min-radius filter did
  nothing — the signature of broken 3D connectivity, not an object-size problem).
- **Speculation (flagged):** that isotropic resampling would help (it would not
  recover unsampled tubules); and the exact split between anisotropy vs
  perinuclear-clumping contributions is unquantified.

## 6. Consequences carried into the design

- **"3D deferred" is doing real work.** For *this* instrument/protocol, 3D isn't a
  free future upgrade — it needs an **acquisition** change (finer z-step ~0.25 µm,
  or Airyscan / deconvolution toward isotropy), not new software. State this so a
  future reader doesn't assume 3D is one flag away.
- **Nellie is at best an *optional* segmentation back-end**, not the default — for
  hypothetical future isotropic data. On the data we actually have, tuned classical
  2D-MIP segmentation is the correct default. This tempers the earlier
  "anchor on Nellie" language (see `00-SYNTHESIS §2`, `07-positioning-and-moat.md`).
- **All morphometry is on the 2D MIP** — a consistent, systematic bias, so
  *relative* between-condition contrasts remain valid; absolute 3D shape does not.
  This caveat already leads every stats artifact and should stay there.
