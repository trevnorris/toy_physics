---
unit_id: 076
batch: III.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 076

This is the iteration-2 verification for the current 2-finding audit of unit 076
(F1 = tautological_check on `Theta_w` alternative-form assertion; F2 =
hardcoded_result provenance comment on `ref_factor = 1/20`). A prior
iteration-1 verification existed at this path against a 3-finding directive
and has been overwritten.

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- SymPy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.py:61-67`. The block

  ```python
  # Independent route: Theta_w as (2 rho_w mu_star / (hbar c_sw))^2
  Theta_w_alt = sp.simplify((2 * rho_w * mu_star_solved / (hbar * csw))**2)
  print("Theta_w (enthalpy lock) =", Theta_w)
  expect_zero("Theta_w vs alternative-form derivation", Theta_w - Theta_w_alt)
  ```

  is replaced with

  ```python
  # Closed-form target from notes section 2: Theta_w = lambda_mu^2 m^2 rho_w^2 c_sw^2 / (4 hbar^2).
  # This independently states the simplified form; the assertion below exercises the /4 factor
  # in the enthalpy lock (mu_* = lambda_mu * m * c_sw^2 / 4).
  Theta_target = sp.Rational(1, 4) * lambda_mu**2 * m**2 * rho_w**2 * csw**2 / hbar**2
  print("Theta_w (enthalpy lock) =", Theta_w)
  expect_zero("Theta_w under enthalpy lock", Theta_w - Theta_target)
  ```

- Mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.wl:54-58`. The `thetaCanonical = FullSimplify[(2*rhoW*muStarSolved/(hbar*cSw))^2, ...]` line is replaced with `thetaTarget = FullSimplify[(1/4)*lambdaMu^2*mpsi^2*rhoW^2*cSw^2/hbar^2, ...]`, and the `expectZero` label is renamed from "Theta_w vs alternative-form derivation" to "Theta_w under enthalpy lock".

**Assessment:**
The edit matches the directive's "After" blocks byte-for-byte (cross-checked
against `redteam/exec_logs/stage_076_diff.patch` and the current state of both
edited files). The new `Theta_target` is written *independently of
`mu_star_solved`* — it depends only on the primitive symbols
`lambda_mu, m, rho_w, csw, hbar` with the `1/4` factor hand-typed. The LHS
`Theta_w` is still constructed from `4 * rho_w**2 * mu_star_solved**2 /
(hbar**2 * csw**2)` with `mu_star_solved = sp.solve(enthalpy_lock,
mu_star_sym)[0]` and `enthalpy_lock = mu_star_sym - lambda_mu * m * csw**2 /
4`. So the residual `Theta_w - Theta_target` now substantively exercises the
`/4` factor in the enthalpy lock. Per the original report's self-test (notes
section §195(3)), if the enthalpy lock factor were changed to `/5`, the
residual would be `16/25 - 1 != 0`, confirming the check is no longer the
trivial `(2x)^2 = 4 x^2` identity.

Output evidence:
- SymPy `scripts/output/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.txt` line 9: `Theta_w (enthalpy lock) = c_sw**2*lambda_mu**2*m**2*rho_w**2/(4*hbar**2)`; line 10: `Theta_w under enthalpy lock = 0`.
- Mathematica `mathematica/output/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.txt` line 10: `Theta_w (enthalpy lock) = (cSw^2*lambdaMu^2*mpsi^2*rhoW^2)/(4*hbar^2)`; lines 11–12: `Theta_w under enthalpy lock = 0` then `PASS: Theta_w under enthalpy lock`.

No collateral edits beyond the comment text and the renamed assertion label.

### F2 — hardcoded_result

**Classification:** resolved

**What changed:**
- SymPy lines 79–83: the prior block

  ```python
  # Reference-branch convention: ell = a * ref_factor with ref_factor = 1/20.
  # TODO(provenance): cite the upstream stage that fixes ref_factor. This factor is
  # the load-bearing piece of the "25" in the normalized reference identity.
  ref_factor = sp.Rational(1, 20)  # reference-branch convention: ell = a * ref_factor  (see F2 below for provenance)
  ```

  is replaced with

  ```python
  # Reference-branch convention: ell = a * ref_factor with ref_factor = 1/20.
  # Source: Family-1 reference-branch description carried forward as input to this stage
  # (notes/stages/moving_throat_pde_stage076_n5_wall_depth_lock.md section 4).
  # This factor is the load-bearing piece of the "25" in the normalized reference identity.
  ref_factor = sp.Rational(1, 20)
  ```

- Mathematica lines 69–73: the `TODO(provenance)` line inside the `(* ... *)` block is replaced with the `Source: Family-1 reference-branch description ... section 4` text; the trailing inline comment on `refFactor = 1/20;` is removed. Numeric value unchanged.

**Assessment:**
Edit matches the directive's "After" blocks exactly (verified against the diff
and current file state). The `TODO(provenance)` token is removed from both
files. The upstream anchor `notes/stages/moving_throat_pde_stage076_n5_wall_depth_lock.md section 4` is now cited inline. No assertion change, no value change. The transcripts show byte-identical numeric content compared to the pre-fix outputs cited in the audit report:
- SymPy line 14: `Theta_w (reference branch, general a) = 25*lambda_mu**2*rho_w**2/a**2`
- SymPy line 15: `Theta_w (reference branch, normalized wall units) = 25*lambda_mu**2*rho_w**2`
- SymPy line 16: `normalized reference factor = 0`
- Mathematica lines 17, 18, 19–20 mirror these with the matching Mathematica symbols and PASS marker.

Comment-only change, as the directive specified.

## Exec log assessment

**SymPy:** exit=0 (per orchestrator pre-context: "Both engines exit 0"). The
file `/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/stage_076_sympy.log` is not present on disk — only `stage_076_diff.patch` was retained in `exec_logs/` for this iteration. The canonical output transcript at `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.txt` was regenerated post-fix (mtime 2026-05-27 02:15 vs script mtime 2026-05-27 02:13). Notable lines:

```
n=5 enthalpy identity = 0
n=3 residual (should be NONZERO) = 3*K*rho**2/4
Theta_w under enthalpy lock = 0
healing-lock reduction = 0
normalized reference factor = 0
```

The script uses `raise AssertionError` on any nonzero residual, so a clean
output through line 16 implies exit 0.

**Mathematica:** exit=0 (per orchestrator pre-context). Similarly the
`stage_076_mathematica.log` file is not present in `exec_logs/`. Canonical
output transcript at `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.txt` (mtime 2026-05-27 02:19 vs script mtime 2026-05-27 02:13) terminates with:

```
PASS: n=5 enthalpy identity
PASS: Theta_w under enthalpy lock
PASS: healing-lock reduction
PASS: normalized reference factor

Stage 076 Mathematica audit passed.
```

The script's `expectZero` calls `fail[...] -> Exit[1]` on nonzero residual; the
"audit passed" trailer plus four PASS lines plus `Exit[0]` confirms a clean
run.

**Output freshness:** confirmed. SymPy output 02:15 > script 02:13. Mathematica
output 02:19 > script 02:13. Both `.txt` outputs are newer than their
respective sources, so the transcripts reflect the post-fix state.

## Material-change assessment

`material_change`: false.

F1 changes the *form* of one assertion (replacing a tautological check with a
substantive closed-form comparison) but does not change any printed numeric or
symbolic value. The line `Theta_w (enthalpy lock) = c_sw**2*lambda_mu**2*m**2*rho_w**2/(4*hbar**2)` (sympy line 9) and the reference-branch results `25*lambda_mu**2*rho_w**2/a**2` and `25*lambda_mu**2*rho_w**2` (lines 14–15) are byte-identical pre- and post-fix. Only the *label* of one assertion changes ("Theta_w vs alternative-form derivation" → "Theta_w under enthalpy lock"), which no downstream unit can observe. F2 is comment-only. The headline boxed identity `Theta_w = 25 lambda_mu^2 rho_w^2` is unchanged. No downstream unit needs re-auditing on the basis of this stage's edits.

## Side observations (non-blocking)

- Banner strings still read `STAGE 59` (SymPy line 30, docstring line 3) and `STAGE 059` (Mathematica line 26). The auditor flagged this as a cosmetic note (original report line 187, "not blocking"), and the directive did not include it. The Mathematica completion line correctly says "Stage 076 Mathematica audit passed." Mentioned for awareness only.
- Exec-log files at the expected paths `redteam/exec_logs/stage_076_{sympy,mathematica}.log` are not present. The orchestrator's pre-context note guarantees both engines exited 0, and the canonical `.txt` transcripts are fresh and consistent with clean runs, so this does not block verification — but future iterations should preserve those artifacts at the documented paths for auditability.

## Verdict justification

Both findings are resolved. The edits in `redteam/exec_logs/stage_076_diff.patch` match the directive's "After" blocks character-for-character in both the SymPy and Mathematica scripts. F1's new assertion is non-tautological: it compares `Theta_w` (constructed from `mu_star_solved` via the enthalpy lock with `/4`) against an independently-written closed-form target `(1/4) * lambda_mu^2 * m^2 * rho_w^2 * c_sw^2 / hbar^2`; perturbing the enthalpy lock factor to `/5` would make the residual nonzero per the auditor's self-test. F2 is a comment-only edit that removes the `TODO(provenance)` token and cites the upstream notes §4 anchor. Both engines exit 0 (per orchestrator pre-context), the canonical `.txt` outputs are regenerated post-fix and show all assertion residuals at `0` with the renamed labels, and no printed numeric or symbolic value changed. No regressions, no material downstream impact, no collateral edits.
