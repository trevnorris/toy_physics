---
unit_id: 121
batch: retro
created_at: 2026-06-01T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-01T21:48:42Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 121 (retro-sweep: dual-engine .wl)

> This is a RETRO-SWEEP directive. Stage 121 was audited + verified in batch IV.3 (SymPy-only).
> Under the dual-engine rule (a Mathematica `.wl` is REQUIRED on every stage Mathematica CAN
> independently verify), it is missing its second engine. The ONLY change is to ADD the `.wl`.
> (Any prior IV.3 directive content for this unit is preserved in git history.)

The SymPy audit script for this stage is correct and is the REFERENCE engine. Do NOT modify it. Do NOT touch `paper.tex` or `notes/`. The only required change is the dual-engine gap below.

After creating the script, RUN it (`timeout 600 math -script <path>`) and iterate until it exits 0 with all checks passing. A timeout (exit 124) is a FAILURE — reformulate the math, never raise the cap. The orchestrator independently re-runs afterward.

## F1 — missing_mathematica

**Issue:** Stage 121 ("geometric selection of r") is dual-engine-capable — every claim is elementary symbolic algebra over positive reals (one quadratic-in-`r` solve, substitutions, a nested-radical exact-value identity, and a threshold evaluation), with no integral, BVP, or transcendental root — but it has no Mathematica `.wl`. Under the dual-engine rule an independent second-engine verification is required wherever Mathematica can do the math.

**Required change (you design the route and write the script):**
Create `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage121_geometric_r_selection_mathematica_audit.wl`.
- Independently re-verify EVERY load-bearing assertion in the SymPy script `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage121_geometric_r_selection_sympy_audit.py`. Read that script to enumerate the claims and their target conclusions; the paper card `paper/stages/stage_121.tex` and the stage notes file are the source of truth for the math. (The paper card uses an older section header — anchor to the live numbering and the `.py`; do NOT edit the card or notes.)
- Use Mathematica-NATIVE primitives (`Solve`/`Reduce`, `Simplify`/`FullSimplify`/`RootReduce`/`PowerExpand` under positivity `$Assumptions`, exact `N[...]`, etc.) via a DIFFERENT derivation route than the SymPy script — NOT a line-by-line port with the same variable names and step order. Reference an existing verified `.wl` (e.g. `mathematica/moving_throat_pde_stage187_*_mathematica_audit.wl`) ONLY for house idioms (the `expectZero` helper that `Exit[1]`s on failure, `$Assumptions` positivity declarations, `stripCE` for `ConditionalExpression`, the `math -script` convention).
- Assert cross-engine agreement: each conclusion the `.wl` derives must match the SymPy result.

**Claim manifest** (the `.wl` must independently witness each; match the SymPy PASS-label so the verifier can pair them):
- **M1** — `r_geom closed-form (explicit)`: imposing `L_W = L` on the Stage-99 length law `L_W = (π a/2)·Sqrt((1+r²)/3)` and solving for `r` (positive branch) gives `r_geom² = 12 L²/(π² a²) − 1`.
- **M2** — `tube-length relation`: back-substituting `r_geom` into `L_W = (π a/2)·Sqrt((1+r²)/3)` returns `L_W = L` exactly.
- **M3** — `r_F1 symbolic value at L/a = 37/20`: evaluating `r_geom` at `L/a = 37/20` yields `r_F1 = Sqrt(4107 − 100π²)/(10π)` (numeric ≈ 1.77799353547498).
- **M4** — `r_c(F1) symbolic value`: `r_c = r_F1² = 4107/(100π²) − 1` (numeric ≈ 3.16126101219081).
- **M5** — `Omega_W identification at L_W = L`: with `L_W = L`, `Ω_W = π c_s/(2 L)`. NOTE: in the `.py` this reduces to a definitional substitution of the half-wave pole `Ω_W(L_W) = π c_s/(2 L_W)`; a cross-engine PARITY check (state the half-wave relation, evaluate at `L_W = L`) is acceptable here — it is low-information, not a deep identity. Do not dress it up as an independent derivation.
- **M6** — `r_geom vanishes at existence threshold`: `r_geom = 0` at the exact threshold `L/a = π/(2 Sqrt(3))` (substitute the EXACT symbol, not a decimal, so the residual is an exact 0).

**Route / correctness requirements (acceptance criteria — you choose how to satisfy them):**
- The positive branch of the quadratic solve must be selected explicitly (positivity `$Assumptions` and/or `PowerExpand` of a positive radicand) — do not accidentally grab the negative root.
- The nested-radical equality in M3 (`r_F1` derived vs the boxed surd) must be certified symbolically (`RootReduce`/`FullSimplify`, or compare squares plus a positivity check) — a numeric-only match is NOT sufficient (the audit is a symbolic-zero check); also reproduce the ≈ decimal via exact `N`.
- `L/a = 37/20` is an imported anchor (the carried preferred aspect ratio) and `L_W = (π a/2)·Sqrt((1+r²)/3)` is the Stage-99 length law — treat both as upstream-owned inputs, not as things re-derived here. The stage's content is the CONSEQUENCE of imposing `L_W = L` on that law.

**Anti-transliteration:** a `.wl` that merely re-types the SymPy closed forms (`r_geom = Sqrt(12 L²/(π²a²) − 1)`, `r_F1 = Sqrt(4107 − 100π²)/(10π)`, `r_c = 4107/(100π²) − 1`) and subtracts them against themselves is a transliteration and will be REJECTED at verification. The independent route is to START from the Stage-99 length law, impose `L_W = L`, and let Mathematica DERIVE `r_geom` (positive branch via `Solve`/`Reduce`); then SUBSTITUTE `L/a = 37/20` and let Mathematica PRODUCE `r_F1`/`r_c`, with the boxed surd targets appearing only as the right-hand side of the comparison.

**Comment hygiene:** avoid any `*)` substring inside Mathematica comments (premature comment close — a recurring defect).

**Verification command:** the verifier runs `redteam exec-mathematica 121`, confirms exit 0 with all PASS lines (≥ 6 substantive checks: M1–M6), and reviews that the `.wl` is a genuinely independent route (native primitives, different decomposition) whose conclusions agree with the SymPy engine.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage121_geometric_r_selection_mathematica_audit.wl`
- summary: Added the missing Mathematica audit deriving the positive geometric radius branch from the Stage-99 tube-length law and checking M1-M6 against the SymPy conclusions.
- deviation: none
