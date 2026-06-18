# Directive — M1c-prep: torch→Mathematica self-consistent background pipeline + A≠0 revalidation

You are the **CODER** (Codex codes / Claude reviews). Build the TARGET-BLIND cross-engine machinery for M1c
(decision 06): export the *self-consistent* torch WP1 background, import it into the Mathematica derived chain
(replacing the smoke ρ0), and **re-validate the derived chain with the WP1 gauge fields ACTIVE**. This is
M1c-PREP only — **NO frozen values, NO GATE A, NO physical claim, NO frozen committed packet** (those are
M1c-run). Placeholder params stay; any R_norm computed here is MACHINERY (real self-consistent background but
NOT frozen physical params). Follow decision 06 (M2/M3) + `_scratch/m1c_background_design.log`.

## Part 1 — torch self-consistent WP1 background export (Python, target-blind)
In the torch package `software/stage1_solver/src/stage1_solver/` (NOT the firewalled
`research/pde_audit/simulation/`; do NOT touch `physical_export_permitted` — that guard is unrelated to this
target-blind background export), add an exporter that:
- Runs the VALIDATED `run_branch_continuation` WP1 stationary solve (`coupled_branch.py`) at a TRACTABLE
  resolution (the existing step-3 smoke resolution is fine — the point is self-consistency + the pipeline, not
  high-res; keep each solve under `timeout 600`), with the current PLACEHOLDER (target-blind) params.
- **Confirms the solve is genuinely self-consistent: the coupled stationary residual is SMALL (~1e-8 Newton
  floor), NOT the M1b smoke ~243.** Record the residual norm + convergence evidence. This is the whole point —
  a real stationary solution, not an analytic profile.
- Exports a canonical JSON background bundle: grid (faces/centers/volumes, measure), fields
  `ψ_R0,ψ_I0,A_00,A_r0,A_w0`, derived `ρ0`, `R0(w)`, `Z(w)`, the config/params, source revision, the residual
  norms, and a **content-only hash** (canonical sorted JSON, NO live timestamp in the hashed content — set up
  the M1c byte-reproducibility pattern now). Export to a gitignored `runs/` path.

## Part 2 — Mathematica import + shared-background swap + A≠0 revalidation (.wls)
A new producer (e.g. `mathematica/mt15_06_m1c_prep_crossengine.wls`):
- Imports the torch background JSON; builds interpolation/projection (conservative projection for ρ0 + integral
  weights, point interpolation for smooth fields); records interpolation order + ONE refinement.
- Replaces the hardcoded smoke ρ0/geometry/gauge in the M1b BdG (`mt15_02`) and the Spike-2 transfer
  (`mt15_04`) logic with this ONE shared imported background (reuse/refactor those modules — do not fork a
  third smoke background).
- **A≠0 REVALIDATION (the #1 risk).** The smoke background had ALL gauge fields zero; the real WP1 has the
  scalar `A_00 ≠ 0` (it enters the BdG `q·A0` term — dormant before) and the perturbation gauge response `δA`
  active in the transfer. (Confirm `A_r0/A_w0` on the isotropic WP1 — expected ~0 by isotropic symmetry; state
  what you find.) Re-run the full derived chain on the real background and re-run EVERY validity check with the
  gauge terms now active: BdG eigen-residuals + Stieltjes; the Spike-2 current Fréchet (vs the 8c phasor with
  the now-active gauge terms), pure-gauge zero-transfer, basis-invariance, V2-09 regression, Green residuals,
  N0>0. They must still pass. **Characterize what CHANGED vs the smoke** (which gauge terms activated; any
  spatial-gauge SIGN sensitivity).
- Run the V2-22C/V2-21 direct-derived chain → R_norm, with a LOUD label: **real self-consistent background but
  PLACEHOLDER (unfrozen) params → MACHINERY, NOT physical; proves the cross-engine pipeline only.**

## Acceptance criteria
- Torch WP1 export residual SMALL (self-consistent, ~1e-8) + convergence evidence — NOT smoke-level.
- MMA import faithful: imported fields match torch fields to interpolation order (report order + the one
  refinement); ONE shared background drives both M1b and Spike-2.
- A≠0 revalidation: ALL Spike-1/2 + M1b validity checks PASS on the real background; a written characterization
  of what changed when the gauge terms activated (esp. A_00 in the BdG).
- End-to-end cross-engine pipeline runs → a (machinery) R_norm, loudly labeled non-physical.
- Byte-reproducibility pattern set up: background export + any packet use content-only hashing; two reruns →
  identical hashes (flag any residual nondeterminism).
- NOTHING frozen; NO GATE A; NO physical claim; firewalled `simulation/` untouched; `physical_export_permitted`
  untouched.

## Working rules
Codex codes BOTH the Python exporter + the `.wls`. ≤2 `math`/`wolframscript` kernels; `timeout 600` per script
RUN (incl. the torch solve — pick a resolution that converges within it); never wrap your session in a shell
`timeout`. Files: Python under `src/stage1_solver/`, `.wls` + report under `mathematica/`; run artifacts →
gitignored `runs/`. No `git add`/commit. No network/GPU.

## Report back
The torch WP1 export (residual norm proving self-consistency, resolution, what fields exported, the
content-hash); the MMA import fidelity (interp order + refinement); the A≠0 revalidation results (every check
PASS/with-numbers + the characterization of what activated/changed vs smoke); the machinery R_norm with its
loud label; the byte-reproducibility result. Then Claude runs the transliteration-fidelity audit (torch export
+ MMA import + the gauge-active operators) + adversarial (self-consistency genuine? revalidation checks genuine
not can't-fail? R_norm honestly labeled? target-blind? no guard/firewall touched?) + arbiter.
