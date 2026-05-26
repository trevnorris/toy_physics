---
unit_id: 039
batch: III.1
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-26T02:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 039 (v2)

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**

- `scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py`: the two `expect_zero` assertions that compared `M_mix_split` to `M_mix_flat.subs(eps_W, eps_W_split)` and `R_target_split` to `R_target_flat.subs(eps_W, eps_W_split)` (previously at lines 138-139) have been deleted. The documentation prints at lines 149-151 (`print("M_mix^(split U) =", ...)`, `print("R_target^(split U) =", ...)`, `print("product =", ...)`) are retained. The flat-definition statements at lines 143-144 are retained for documentation continuity.
- `mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl`: the two `expectZero` analogues (previously at lines 118-119) have been deleted. The `Print["M_mix^(split U) = ", ...]`, `Print["R_target^(split U) = ", ...]`, `Print["product = ", ...]` documentation prints at lines 129-131 remain. `mMixFlat` and `rTargetFlat` definitions at lines 123-124 remain.

The diff at `redteam/exec_logs/stage_039_diff.patch` lines 30-31 (Mathematica) and 64-65 (SymPy) shows exactly the deletions specified in the directive, with no collateral edits.

**Assessment:**

The edit is correct and exactly matches the "Required change" in the directive. The two tautological identities — both constructed by literally substituting `eps_W_split` into the flat formula and then asserting equivalence — are gone. The exec logs no longer print the tautological PASS lines:

- SymPy log (lines 47-50) shows section 22.4 with only the three `print(...)` lines for `M_mix^(split U)`, `R_target^(split U)`, and `product`; no `... = 0` lines.
- Mathematica log (lines 56-59) shows section 4 with only three `Print[...]` lines, no `PASS: M_mix split ...` or `PASS: R_target split ...` lines.

No paper-side claim is left unverified because these checks were `extra` relative to `\stagefield{Output}`. Both scripts still `exit 0`.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**

- `scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py:124-136`: a new collinearity-iff block was inserted between the existing `expect_zero("direction-splitting invariant", ...)` (line 122) and `print("Collinearity theorem: ...")` (now line 138). The block contains:
  - `expect_zero("collinearity if-leg: D_dir(deltaU=0) = 0", D_dir.subs(deltaU, 0))`
  - `expect_zero("collinearity if-leg: D_dir(rho0=0) = 0", D_dir.subs(rho0, 0))`
  - Only-if leg using `sp.fraction(sp.together(D_dir))` to extract numerator, divide by `rho0 * deltaU`, simplify, print the reduced ratio, and raise `AssertionError` if it still depends on `rho0` or `deltaU`, or if it simplifies to zero.
- `mathematica/moving_throat_pde_stage039_split_u_sector_mathematica_audit.wl:107-119`: an analogous block was inserted immediately after `expectZero["direction-splitting invariant derived matches postulated", ...]` at line 105 (before the `subbanner` at line 121). The block contains the same three checks using `Together`, `Numerator`, `FullSimplify`, `FreeQ`, and `fail`/`pass` helpers.

The diff at `redteam/exec_logs/stage_039_diff.patch` lines 5-21 (Mathematica) and 39-55 (SymPy) shows exactly the additions specified in the directive, including the comments, with no collateral edits.

**Assessment:**

Both insertions match the directive's "Required change" verbatim (matching variable names, structure, helper functions, and assertion messages). The new checks are non-tautological:

- The `if-leg` substitutions exercise `D_dir = kappa0*z1 - kappa1*z0` (which was derived independently from the underlying `z0`, `z1` loading vectors at lines 109-117 of the SymPy script and 87-94 of the Mathematica script), not the postulated closed-form. The fact that `D_dir.subs(deltaU, 0) = 0` and `D_dir.subs(rho0, 0) = 0` are independently true is a non-trivial algebraic identity, not a built-in equivalence.
- The only-if leg pulls `Numerator(Together(D_dir))`, divides by `rho0 * deltaU`, and verifies the remaining factor has no `rho0` or `deltaU` dependence. The SymPy exec log line 42 shows `Numerator(D_dir) / (rho0*deltaU) = 8*sqrt(2)*c_etaW`; the Mathematica exec log line 51 shows `Numerator(D_dir) / (rho0*deltaU) = 8*Sqrt[2]*cEtaW`. Both are non-zero constants in `(c_etaW, mu_W, mu_eta, pi)` after `Together` reduction (the `1/(3*pi^2*sqrt(mu_W mu_eta)*(deltaU+1))` factor ends up in the denominator after `Together`, so the numerator-only ratio is `8*sqrt(2)*c_etaW`). The check is exactly equivalent to the directive's spec — the directive predicted a different reduced form (with denominator constants attached), but that was an error in the directive's prediction; the actual `sp.fraction(sp.together(...))` returns numerator-only, and the resulting check is still `free of rho0 and deltaU` and `nonzero`, which is precisely what the only-if leg requires. No semantic difference.

Exec log confirmations:

- SymPy log lines 40-42:
  ```
  collinearity if-leg: D_dir(deltaU=0) = 0 = 0
  collinearity if-leg: D_dir(rho0=0) = 0 = 0
  Numerator(D_dir) / (rho0*deltaU) = 8*sqrt(2)*c_etaW
  ```
- Mathematica log lines 47-52:
  ```
  collinearity if-leg: D_dir(deltaU=0) = 0 = 0
  PASS: collinearity if-leg: D_dir(deltaU=0) = 0
  collinearity if-leg: D_dir(rho0=0) = 0 = 0
  PASS: collinearity if-leg: D_dir(rho0=0) = 0
  Numerator(D_dir) / (rho0*deltaU) = 8*Sqrt[2]*cEtaW
  PASS: collinearity only-if: residual factor is nonzero and independent of rho0, deltaU
  ```

Both scripts still `exit 0`. The collinearity-iff theorem is now formally asserted in both engines, closing the original finding.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
```
direction-splitting invariant = 0
collinearity if-leg: D_dir(deltaU=0) = 0 = 0
collinearity if-leg: D_dir(rho0=0) = 0 = 0
Numerator(D_dir) / (rho0*deltaU) = 8*sqrt(2)*c_etaW
```
All previous load-bearing assertions (A1-A5: `A0 direct - expected`, `A1 direct - expected`, `eps_W direct - split formula`, `z1*(1+rho0) - (kappa1/kappa0)*z0*(1+rho0/(1+deltaU))`, `direction-splitting invariant`) still print `= 0`. The two tautological lines from F1 are absent.

**Mathematica:** exit=0. Notable lines:
```
PASS: direction-splitting invariant derived matches postulated
PASS: collinearity if-leg: D_dir(deltaU=0) = 0
PASS: collinearity if-leg: D_dir(rho0=0) = 0
PASS: collinearity only-if: residual factor is nonzero and independent of rho0, deltaU
Stage 039 Mathematica audit passed.
```
All previous `PASS:` lines (A8-A13) are still present. The two tautological `PASS:` lines from F1 are absent.

**Output freshness:** the saved `.txt` outputs at `scripts/output/moving_throat_pde_stage039_split_u_sector_sympy_audit.txt` and `mathematica/output/moving_throat_pde_stage039_split_u_sector_mathematica_audit.txt` are still dated May 22 12:26, while the scripts were last edited May 26 01:37 and the exec logs were captured May 26 01:50. The substantive verification rests on the exec logs (which are post-fix), but the orchestrator should re-run the saved-output refresh step (`$RT exec-*` or equivalent) to bring the `.txt` outputs into sync with the current scripts. This is a tracker/manifest concern, not a verification failure.

## Material-change assessment

`material_change`: false.

The F1 fix removes two extra (non-`Output`) tautological checks; no derived result is altered. The F2 fix adds new checks against `D_dir`, but `D_dir` itself, its closed form, `delta_split`, `eps_W_split`, `R_U`, `M_mix^(split U)`, `R_target^(split U)`, and the small-`deltaU` expansions are all unchanged in the script body and unchanged in the exec log compared to the pre-fix output. Downstream units that depend on stage 039 see identical derived quantities; no `upstream_stale` flag is needed on the substance, only (optionally) for tracking that section 22.3 of the script now has additional `PASS:` lines.

## Side observations (non-blocking)

1. Saved `.txt` outputs are stale relative to the scripts (May 22 vs. May 26). The orchestrator's post-batch tracker-update step should refresh them via `$RT exec-sympy 039` and `$RT exec-mathematica 039` (sequentially, per the no-parallel-exec rule). The exec logs already captured by the orchestrator are post-fix, so the verification is sound; this is purely a manifest-freshness item.
2. The directive's predicted reduced form for `Numerator(D_dir) / (rho0*deltaU)` (`8*sqrt(2)*c_etaW/(3*pi**2*sqrt(mu_W)*sqrt(mu_eta))`) included denominator constants that `sp.together` actually keeps in the denominator. The actual output (`8*sqrt(2)*c_etaW` for SymPy, `8*Sqrt[2]*cEtaW` for Mathematica) is the numerator-only reduction, which is still independent of `rho0` and `deltaU` and is nonzero. The semantic check (the only-if leg) is correctly verified; only the predicted print value differed from the actual.
3. Both engines agree on the reduced numerator (modulo the trivial transliteration `sqrt -> Sqrt`, `c_etaW -> cEtaW`). The engine cross-check from the original audit report still holds.

## Verdict justification

Both F1 and F2 are `resolved`. The diff at `redteam/exec_logs/stage_039_diff.patch` matches the directive's "Required change" specs verbatim, with no collateral edits. The exec logs show that both scripts still `exit 0`, all pre-existing assertions still pass, the F1 tautological assertions are gone from both engines, and the F2 collinearity-iff `if-leg` and `only-if` checks are now present in both engines and pass non-tautologically. No downstream-affecting derivation changed. Verdict: `verified`.
