---
unit_id: 216
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 216

## Per-finding outcomes

### F1 — missing_verification_script (subtype: script_doesnt_cover_claim)

**Classification:** resolved

**What changed:**
Codex added two helpers to `scripts/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate_sympy_audit.py` — `expect_zero(name, expr)` (lines 5-7, `assert sp.simplify(expr) == 0`) and `require(name, condition)` (lines 10-11, `assert condition`) — both of which raise `AssertionError` on failure. It then inserted failing checks against every previously-computed-but-unasserted quantity, leaving all original `out.append` prints intact:
- M1: line 39 `expect_zero("M1 gradient-optimal unit norm", grad_norm_sq - 1)`.
- M2: line 40 `expect_zero("M2 gradient-optimal max slope", k_grad - k_norm)` (`k_grad = Σ a·k`, `k_norm = sqrt(S)`).
- M3: lines 56-57, looping over the five faces, `expect_zero(... face gap, diff - k_axis**2)` plus `require(... strict positive gap, (k_axis**2).is_positive is True)`.
- M4: lines 87-88, `expect_zero` on `identity_check` and `slack_check`.
- M5: line 96, `expect_zero("M5 barycenter leverage", w_eq - 4)`.
- M6: line 113, `expect_zero("M6 certified bracket quadratic residual", sp.Rational(1,2)*kappa*tau**2 - k*tau + H0)`.

Grep confirms 10 `expect_zero`/`require` call sites; nothing else in the body was altered (no computed expression's form changed, no new symbol/constant introduced — `positive=True` assumptions on the slopes and on `H0,kappa,k` are unchanged).

**Assessment:**
Correct and complete. Every manifest claim M1-M6 now has a falsifiable assertion whose left side is a residual that must `simplify` to exactly `0` (or, for M3 strictness, a boolean that must be `True`). These are genuine: e.g. M5 asserts `w_eq - 4`, so a sign flip in `w_sigma` would flip the exit code. The CRITICAL check — M6 — is **non-tautological**: it does NOT re-print `tau`; it substitutes the formed `tau` into the bracketing quadratic `(1/2)κτ² − kτ + H0` and requires the residual to vanish, which is exactly "the closed form solves the quadratic." If `tau` were the wrong root or had a wrong denominator, the residual would be nonzero and the assert would fire. SymPy exits 0 with all checks active.

### F2 — missing_verification_script (subtype: missing_mathematica)

**Classification:** resolved

**What changed:**
Codex created `mathematica/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate_mathematica_audit.wl` (161 lines). It declares the slopes and `H0,k,kappa` as positive reals via `$Assumptions` (lines 61-66, including the real-branch condition `k^2 > 2 H0 kappa`), uses `expectZero`/`expectTrue` failing-test helpers backed by `fail[...] -> Exit[1]` (lines 21-48), prints a labeled transcript, and `Exit[0]` only when all checks pass. It carries an `expectZero` `ConditionalExpression`-stripping idiom (lines 27-34) consistent with project policy. All six manifest claims have independent failing checks.

**Assessment — independence is LOAD-BEARING, not rubber-stamped:**
This is a genuinely independent derivation, NOT a transliteration of the `.py`. The two engines use materially different decompositions, and the two claims the directive specifically required to be *derived* (M2, M6) are derived:

- **M2 (derived via constrained maximization, not re-stated).** The `.py` simply forms the closed-form ray `a_grad = k/||k||`, computes `k_grad = Σ a·k`, and checks `k_grad − √S = 0`. The `.wl` instead builds a Lagrangian and solves the KKT system:
  - `lagrangian = kVec.aVec - mu (aVec.aVec - 1)` (line 72)
  - `stationarity = Thread[(D[lagrangian, #] & /@ aVec) == 0]` (line 73)
  - `lagrangeRules = Solve[Join[stationarity, {aVec.aVec == 1}], Append[aVec, mu], Reals]` (line 74), then selects the positive-`mu` branch (lines 75-78) and reads off `lagrangeValue = kVec.aVec /. positiveRule`, checking `lagrangeValue - Sqrt[S] == 0` (line 88) plus that this branch dominates the other stationary branch (lines 89-92). So `√S` is *obtained as the constrained maximum*, not asserted by hand. (M1 likewise read from the same Solve, line 87.)

- **M6 (derived via solving the quadratic, not the `.py` substitution route).** The `.py` substitutes the stated `tau` into the residual. The `.wl` instead solves the quadratic for `tau`:
  - `rootRules = Solve[quadraticResidual == 0, tau, Reals]` (line 138), `roots = ... stripConditional[tau /. rootRules]` (line 139)
  - selects the root equal to the stated closed form `statedTau = 2 H0/(k + Sqrt[k^2 - 2 H0 kappa])` (lines 140-145) and identifies the other root (lines 146-150)
  - then `expectZero["M6 solved smaller root minus stated bracket form", solvedTau - statedTau]` (line 153), the residual check (line 154), positivity (line 155), and smaller-root ordering `statedTau < otherTau` (line 156). The log confirms `Solve` returned `{(k - Sqrt[k^2 - 2 H0 kappa])/kappa, (k + Sqrt[...])/kappa}` and that the rationalized stated form matched the smaller root.

- **M4/M5 (different decomposition: matrix quadratic forms vs scalar sums).** The `.py` writes `w_sigma` as an explicit scalar sum of products and `pair_sum` as a Python double-loop. The `.wl` realizes them as quadratic forms: `offDiagonalForm = (all-ones) - IdentityMatrix` and `laplacianForm = n·I - Outer[Times,ones,ones]`, with `wSigma = aVec.offDiagonalForm.aVec` and `pairDifferenceEnergy = Total[(Subtract @@ #)^2 & /@ Subsets[aVec,{2}]]` (lines 111-124), and M5's `w_Sigma(eq)=4` is cross-checked a second way via `Eigenvalues[offDiagonalForm] = {4,-1,-1,-1,-1}` and `Max[...] - 4 == 0` (lines 129-133). This is not a port of the `out.append` choreography.

No new physical constant is introduced; the `.wl` uses only the manifest symbols/values. Mathematica exits 0 with 22 PASS / 0 FAIL.

## Exec log assessment

**SymPy:** exit=0. The log is the printed transcript (assertions raise rather than print on success, which is expected for `assert`-based checks); notable lines confirm the asserted quantities took their paper values:
- `||a_grad||^2 = 1` (M1), `k_grad = sqrt(k_U**2 + k_W**2 + k_c**2 + k_gamma**2 + k_lambda**2)` (M2)
- `(k_5^grad)^2 - (k_Q_hatlambda^grad)^2 = k_lambda**2` etc. (M3), `w_sigma - (...) = 0` and the slack `= 0` (M4), `w_sigma(a_eq) = 4` (M5), `tau(H0,k,kappa) = 2*H0/(k + sqrt(-2*H0*kappa + k**2))` (M6).
`# exit_code: 0`.

**Mathematica:** exit=0, 22 `PASS:` lines, 0 FAIL. Notable lines:
- `M2 positive Lagrange value minus Sqrt[S] = 0` / `PASS` and `M2 positive branch dominates the other stationary branch = True` / `PASS`.
- `quadratic roots from Solve = {(k - Sqrt[k^2 - 2*H0*kappa])/kappa, (k + Sqrt[k^2 - 2*H0*kappa])/kappa}`, then `M6 solved smaller root minus stated bracket form = 0` / `PASS` and `M6 stated bracket quadratic residual = 0` / `PASS`.
- `off-diagonal leverage eigenvalues = {4, -1, -1, -1, -1}`, `M5 top quadratic-form eigenvalue minus 4 = 0` / `PASS`.
- `All Stage 216 Mathematica audit checks passed.` `# exit_code: 0`.

**Output freshness:** confirmed. SymPy `.txt` mtime `2026-06-02 11:32:09` and Mathematica `.txt` mtime `2026-06-02 11:32:09` are both newer than their script mtimes (`11:23:13` for both `.py` and `.wl`) — saved outputs were re-generated post-fix.

## Material-change assessment

`material_change`: false. The edits only (a) added assertions to a script that already computed the correct forms and (b) added a second engine that re-derives the same paper values. No computed expression's symbolic form was changed and no constant was altered, so no downstream-relied result moved.

## Side observations (non-blocking)

- The directive's `## Applied: F2` block names the target as `..._mathematica_audit.wl`, matching the file that exists and that the exec log ran; the original report's F2 "required change" cited a shorter filename (`..._gate.wl` without the `_mathematica_audit` suffix). The actual created file uses the standard `_mathematica_audit.wl` convention and is the one executed — a naming nit only, not a defect.
- The `.py` M3 strict-positivity uses `(k_axis**2).is_positive is True`, which leverages the `positive=True` symbol assumptions; this is sound here and not a weakening (no broader assumption was introduced). Non-blocking.

## Verdict justification

Both findings are resolved. F1: the SymPy script now carries 10 real, falsifiable assertions covering M1-M6 (no longer an unconditional PASS), and the CRITICAL M6 check substitutes `tau` into the bracketing quadratic and requires the residual to vanish — a genuine non-tautological check, not a re-print. F2: a new, genuinely independent `.wl` derives M2 via Lagrange-multiplier constrained maximization and M6 via `Solve` of the quadratic (the two derivations the directive mandated), and uses matrix quadratic forms / eigenvalues for M4/M5 — a different decomposition than the `.py`, so it is not a transliteration. Both engines exit 0 (SymPy 0; Mathematica 22 PASS / 0 FAIL), outputs are fresh, and no new physical constant was introduced. Verdict: verified.
