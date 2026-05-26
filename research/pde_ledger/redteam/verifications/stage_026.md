---
unit_id: 026
batch: II.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26T00:30:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 026 (batch II.1 v2)

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py:203-212` — between the existing `sp.pprint(sp.simplify(K_req))` print (line 201) and the back-substitution check `expect_zero("residual @ K_req", ...)` (now at line 213), Codex inserted the new structural assertion exactly as the directive prescribed:
  ```python
  K_req_paper = (
      B0
      + Q / Delta
      + mhat**2 * kappa**2 * (Omega_U**2 * lambda_W + lambda_R * lambda_U)**2
        / (target * Delta**2)
  )
  expect_zero("K_req - K_req_paper", sp.simplify(K_req - K_req_paper))
  ```
- `mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl:163-173` — between `Print["K_req = ", fmt[kReq]]` (line 161) and the existing `expectZero["residual @ K_req", ...]` (now line 174), Codex inserted:
  ```mathematica
  kReqPaper = FullSimplify[
    b0
    + q/delta
    + mhat^2 * kappa^2 * (omegaU^2*lambdaW + lambdaR*lambdaU)^2
      / (target * delta^2),
    Assumptions -> $Assumptions
  ];
  expectZero["K_req - K_req_paper", kReq - kReqPaper];
  ```

The git diff (`redteam/exec_logs/stage_026_diff.patch`) confirms both edits are local to the prescribed file:line ranges, with no collateral edits. Symbol names match scope (capitalized in SymPy, lowercase in Mathematica), comments match the directive verbatim.

**Assessment:**
The edit is correct and addresses the finding directly. The new assertion is non-tautological: `K_req` comes from `sp.solve(sp.Eq(residual, 0), K)[0]` (resp. `k /. First[Solve[residual == 0, k]]`), while `K_req_paper`/`kReqPaper` is a hand-constructed three-term decomposition built from in-scope symbols `B0, Q, Delta, mhat, kappa, Omega_U, lambda_W, lambda_R, lambda_U, target`. The difference simplifies to zero only if the solver's symbolic output actually matches the paper's three-term structural form; a structurally divergent (but still residual-solving) output would not pass through `sp.simplify`/`FullSimplify` to zero. The assertion is exactly the kind of independent structural check the auditor requested.

The exec logs show both new assertions pass cleanly:
- SymPy log line 127: `K_req - K_req_paper = 0`, followed by line 128 `residual @ K_req = 0`.
- Mathematica log lines 85-88: `K_req - K_req_paper = 0`, `PASS: K_req - K_req_paper`, `residual @ K_req = 0`, `PASS: residual @ K_req`.

One small hygiene note: in the Mathematica file, the new local `kReqPaper` is not added to the `Module[{...}]` declaration at line 148. The directive flagged this conditionally ("if Mathematica complains... add `kReqPaper` to the Module declaration"). The exec log shows no warning and a clean exit 0, so Mathematica is silently treating `kReqPaper` as a Global symbol within the Module. This is a minor scoping hygiene issue but not a regression — the assertion fires correctly with the intended semantics. See side observations.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
```
K_req - K_req_paper = 0
residual @ K_req = 0
On the constant wall branch, bare quadrupole wall stiffness is K = K_eta + 6*T_Omega
```
All prior assertions (`int u0^2 - 1 = 0`, `int f0^2 - 1 = 0`, `kappa_n - expected = 0`, `kappa - 2*sqrt(2)/pi = 0`, `overlap_u0_f0 - kappa = 0`, `overlap_u0_u0 - 1 = 0`) also evaluate to zero. No regressions.

**Mathematica:** exit=0. Notable lines:
```
K_req - K_req_paper = 0
PASS: K_req - K_req_paper
residual @ K_req = 0
PASS: residual @ K_req
```
All prior PASS lines (`int u0^2 - 1`, `int f0^2 - 1`, `kappa_n (analytic) - (fundamental thm)`, `kappa - 2*sqrt(2)/pi`, `overlap_u0_f0 - kappa`, `overlap_u0_u0 - 1`) still fire. The final banner reads `Stage 9 Mathematica audit passed.`, a pre-existing cosmetic template artifact unrelated to F1.

**Output freshness:** Exec logs themselves are fresh (mtime 2026-05-26 00:14, after script mtimes of 2026-05-25 23:35). The orchestrator's exec logs are the authoritative source the verifier consults, and they contain the new `K_req - K_req_paper = 0` line. The committed canonical `.txt` outputs under `scripts/output/` and `mathematica/output/` carry older mtimes (2026-05-21 17:01 / 17:02) and were not refreshed by the red-team exec invocation; those are tracked artifacts the user typically refreshes at commit time and do not affect this verification.

## Material-change assessment

`material_change`: false.

The edits add two new assertions but do not alter any derived expression, output value, or downstream variable. `K_req`, the residual, and all branch-substituted quantities (`C, G_U, G_W, R, Delta, Q, P, B0, Z0, N0, D0, P0`) print identically before and after the fix. No downstream unit consumes the new symbol `K_req_paper`/`kReqPaper`. Downstream stages (notably 032's use of `mhat_-^2 = s_-/kappa_0^2` for the D/N source map, and the Part-II appendix closure for this branch) are unaffected.

## Side observations (non-blocking)

1. In `mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl` line 148, the `Module[{kappa, delta, q, p, b0, z0, n0, d0, target, residual, kReq, kGeom}]` declaration does not include the newly introduced `kReqPaper`. Mathematica accepts this silently (treats it as a Global symbol), and the exec passes cleanly, so the assertion semantics are intact. A future cleanup could add `kReqPaper` to the Module locals list for hygiene, as the directive suggested conditionally. Not a verification blocker.
2. The Mathematica script's final banner reads `Stage 9 Mathematica audit passed.` rather than `Stage 26 ...`. This is a pre-existing template carryover from earlier stages and was not within F1's scope. Flagging only as a note.
3. The committed `scripts/output/.../sympy_audit.txt` and `mathematica/output/.../mathematica_audit.txt` files for stage 026 are older than the script mtimes and do not yet contain the new `K_req - K_req_paper = 0` line. This is expected: the red-team exec pipeline writes to `redteam/exec_logs/`, not to the canonical `output/` directories. The user can refresh canonical outputs at commit time.

## Verdict justification

F1 is fully resolved: both target files received the prescribed structural assertion at the prescribed line ranges with no deviation, the diff is minimal and local, the exec logs show both new checks pass with exit 0, and the new assertion is non-tautological because it compares `Solve`'s symbolic output against an independently-constructed three-term decomposition (`B0 + Q/Delta + mhat^2 * kappa^2 * (...)^2 / (target * Delta^2)`). All previously-passing assertions continue to pass; no closed-form values changed. No regressions, no collateral edits, no downstream material change. Verdict: `verified`.
