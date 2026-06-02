---
unit_id: 122
batch: retro
created_at: 2026-06-01T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-01T21:49:45Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 122 (retro-sweep: dual-engine .wl)

> This is a RETRO-SWEEP directive. Stage 122 was audited + verified in batch IV.3 (SymPy-only;
> its IV.3 traction-ratio de-tautologization is already applied and verified). Under the
> dual-engine rule (a Mathematica `.wl` is REQUIRED on every stage Mathematica CAN independently
> verify), it is missing its second engine. The ONLY change is to ADD the `.wl`.
> (The prior IV.3 directive content for this unit is preserved in git history.)

The SymPy audit script for this stage is correct and is the REFERENCE engine. Do NOT modify it. Do NOT touch `paper.tex` or `notes/`. The only required change is the dual-engine gap below.

After creating the script, RUN it (`timeout 600 math -script <path>`) and iterate until it exits 0 with all checks passing. A timeout (exit 124) is a FAILURE — reformulate the math, never raise the cap. The orchestrator independently re-runs afterward.

## F1 — missing_mathematica

**Issue:** Stage 122 ("mouth source compensation test") is dual-engine-capable — every load-bearing claim is closed-form scalar algebra over surds (`Sqrt(4107−100π²)`, `Sqrt(3)`, `π`): the two roots of a compensation quadratic, a polynomial-defect identity, a non-vanishing check, and a ratio identity — but it has no Mathematica `.wl`. Under the dual-engine rule an independent second-engine verification is required wherever Mathematica can do the math.

**Required change (you design the route and write the script):**
Create `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage122_mouth_source_compensation_test_mathematica_audit.wl`.
- Independently re-verify EVERY load-bearing assertion in the SymPy script `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py`. Read that script to enumerate the claims and their target conclusions; the paper card `paper/stages/stage_122.tex` and the stage notes file are the source of truth for the math. (The paper card / notes use older "Stage 139 / 221 / 223" numbering — anchor to the live `scripts/` owners: stage 119 owns the compensation law + the `g = C/T_m` traction law, stage 121 owns `r_F1`. Do NOT edit the card or notes.)
- Use Mathematica-NATIVE primitives (`Solve`/`Reduce`/`Roots` for the quadratic, `FullSimplify`/`RootReduce`/`Together`/`Cancel` under positivity `$Assumptions`, exact `N[...]`) via a DIFFERENT derivation route than the SymPy script — NOT a line-by-line port. Reference an existing verified `.wl` (e.g. `mathematica/moving_throat_pde_stage187_*_mathematica_audit.wl`) ONLY for house idioms (the `expectZero` helper that `Exit[1]`s on failure, an `expectNonzero` analog, `$Assumptions`, `stripCE`, the `math -script` convention).
- Assert cross-engine agreement: each conclusion the `.wl` derives must match the SymPy result.

**Claim manifest** (match the SymPy PASS-label so the verifier can pair them):
- **M1** — `gminus exact form` / `gplus exact form`: the compensated-branch roots at `r = r_F1` equal the exact surds `g_∓ = (2 Sqrt(4107−100π²) ∓ 37 Sqrt(3))/(20π)` (minus branch carries `−37 Sqrt(3)`; numeric `g_− ≈ 0.758035078944663`, `g_+ ≈ 2.79795199200529`).
- **M2** — `compensation quadratic at gminus` / `compensation quadratic at gplus`: each root satisfies `1 + r_F1² − 4(g − r_F1)² = 0`.
- **M3** — `defect closed form`: the natural branch (`g_nat = 1`) is OFF the compensation surface with exact defect `C_nat := 1 + r_F1² − 4(g_nat − r_F1)² = (−12321 + 80π Sqrt(4107−100π²))/(100π²)` (numeric ≈ 1.74016524722739).
- **M4** — `natural off compensation`: `C_nat` is strictly NONZERO (anchor with the exact closed form AND a precision-`N` check that it is ≈ 1.7402, clearly away from 0).
- **M5** — `traction ratio (-) identity` / `traction ratio (+) identity`: with the stage-119 traction law `g = C/T_m` (`C` the SAME background-fixed constant on every branch), `T_m^(±)/T_m^nat = (C/g_±)/(C/g_nat) = g_nat/g_± = 1/g_±` since `g_nat = 1` (numeric ratios ≈ 1.31920016339112 and ≈ 0.357404273860789).

**Route / correctness requirements (acceptance criteria — you choose how to satisfy them):**
- **M5 is the de-tautologization the IV.3 SymPy fix established — the `.wl` MUST NOT re-introduce the tautology.** Carry `C` as a POSITIVE FREE SYMBOL, build `Tm_nat = C/g_nat`, `Tm_± = C/g_±`, form the ratio, and let Mathematica COMPUTE the cancellation of `C`. Do NOT set `T_ratio := 1/g_±` directly, and do NOT pre-substitute `g_nat = 1` in a way that hides the cancellation. `g_nat = 1` is the equal-normalized natural-branch ansatz (sourced independently), NOT `1/g_±`. The non-tautology property: the residual `T_ratio_minus − 1/g_−` reduces to `(g_nat − 1)/g_−`, which is 0 ONLY because the ansatz sets `g_nat = 1` — a mis-stated natural normalization would FAIL.
- **M1**: derive `g_±` by SOLVING the quadratic `1 + r² = 4(g − r)²` at `r = r_F1` (the M2 relation), then confirm they equal the surd forms — root selection (which is `−` vs `+`) pinned by sign so the minus branch carries `−37 Sqrt(3)`. Do not hard-type the roots and only check membership.
- `r_F1` must be reconstructed from the geometric relation `r² = 12 R²/π² − 1` with `R = 37/20` (anchored to stage 121), not imported as a bare numeric.
- Positivity `$Assumptions`: declare `C > 0` and the surd radicands positive (`4107 − 100π² > 0`) so cancellation and root signs resolve.

**Anti-transliteration:** the `.wl` must NOT simply re-type and self-check these closed forms — derive them: the roots `(2 Sqrt(4107−100π²) ± 37 Sqrt(3))/(20π)` (solve the quadratic), the defect `(−12321 + 80π Sqrt(4107−100π²))/(100π²)` (substitute `g_nat=1` and `r_F1` into `1 + r² − 4(g−r)²` and simplify), and the ratios `1/g_±` (emerge from `C` cancellation). The exact surd targets and numeric anchors may appear only as final comparison RHS, never as the starting definitions.

**Comment hygiene:** avoid any `*)` substring inside Mathematica comments (premature comment close).

**Verification command:** the verifier runs `redteam exec-mathematica 122`, confirms exit 0 with all PASS lines (≥ 8 substantive checks: M1×2, M2×2, M3, M4, M5×2), and reviews that the `.wl` is a genuinely independent route (native primitives, different decomposition; M5 cancellation computed not asserted) whose conclusions agree with the SymPy engine.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage122_mouth_source_compensation_test_mathematica_audit.wl`
- summary: Added an independent Mathematica audit deriving the compensation roots, natural-branch defect, nonzero defect check, and traction-ratio identities from symbolic solve and cancellation.
- deviation: none
