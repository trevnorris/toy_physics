---
unit_id: 064
batch: III.3
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 064

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py:95-109` inserts a "MATCHED-LAYER INTEGRAL REDUCTION" block that defines a concrete Gaussian `chi_phi_y = exp(-y^2/(2 L^2))` with constant `H = Hw`, integrates over `(-oo, oo)` to obtain `Npp_int = sqrt(pi)*L`, `I1_int = sqrt(pi)*L/H_w`, `I2_int = sqrt(pi)*L/H_w**2`, and asserts `I1_int - Npp_int/Hw == 0` and `I2_int - Npp_int/Hw**2 == 0`.
- `mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl:68-82` adds the analogous block (`chiPhiY = Exp[-y^2/(2 L^2)]`, `nppInt`, `i1Int`, `i2Int` via `Integrate[...]`), then `expectZero` on the same two reductions.

**Assessment:**
The new checks are non-tautological. `Npp_int`, `I1_int`, `I2_int` are each obtained by an independent `sp.integrate`/`Integrate` call; they are not back-substitutions of the answer. The relationship `I1_int = Npp_int/Hw` is now a real consequence of `H` being constant inside the integral, not of a hand-substitution dictionary. Output confirms `sqrt(pi)*L/H_w - sqrt(pi)*L/H_w = 0` reduces algebraically. Both engines pass.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/...sympy_audit.py:156-165` replaces the old tautological block with `Theta_abs = Hw*Nss`, `Lambda_abs = g_phi*Osp`, `soft_abs = Lambda_abs**2/Theta_abs`, prints `soft_abs = I1**2*g_phi**2/(H_w*I2)` (a non-trivial closure-level form), then applies `soft_abs.subs(I2, I1**2/Npp).subs(Npp, I1*Hw)` to land on `g_phi**2 * I1` and asserts `soft_matched - g_phi**2 * I1 == 0`.
- `mathematica/...mathematica_audit.wl:123-129` mirrors the closure-level derivation with the same `softAbs` / `softMatched` flow.

**Assessment:**
The closure-form `soft_abs = I1^2*g_phi^2/(H_w*I2)` is a real algebraic consequence of `Theta = Hw*Nss = Hw*g_phi^2*I2` and `Lambda = g_phi^2*I1`; it is printed before any substitution, demonstrating non-triviality.

Codex deviated by adding a second substitution `Npp -> I1*Hw` beyond the directive's `I2 -> I1^2/Npp`. The deviation is necessary and grounded: with `I2 -> I1^2/Npp` alone, `soft_abs` reduces to `g_phi^2 * Npp/Hw`, not `g_phi^2 * I1`; the directive's target form `g_phi^2 * I1` additionally requires the matched-layer integral identity `Npp = Hw*I1`, which F1's concrete Gaussian computation independently establishes (`Npp_int - Hw*I1_int = sqrt(pi)*L - Hw*(sqrt(pi)*L/Hw) = 0`). So the extra substitution uses a previously verified closure fact, not a hand-wave. The assertion is non-tautological because both substitutions are derived from independently-checked closures rather than asserted as the answer.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
- `scripts/...sympy_audit.py:56-68` adds a "Local static linear-response closure" block: defines `F_loc = (1/2) H(y) sigma^2 - g_phi chi_phi(y) sigma`, calls `sp.solve(sp.diff(F_loc, sigma_loc), sigma_loc)`, asserts a unique solution, then `expect_zero("closure law chi_sigma = g_phi chi_phi/H", chi_sigma_closure - g_phi*chi_phi(y)/H(y))`.
- `scripts/...sympy_audit.py:70-82` adds an "Integral-level overlap invariants" block: concrete Gaussian profile, `Osp_int_check = integrate(chi_sigma_g*chi_phi_g)`, `Nss_int_check = integrate(chi_sigma_g**2)`, then `expect_zero` on `Osp_int_check - g_phi*I1_int_check` and `Nss_int_check - g_phi^2*I2_int_check`.
- `mathematica/...mathematica_audit.wl:33-56` mirrors both blocks (`Solve[D[fLoc,sigmaLoc] == 0, ...]` and the integral overlap checks).

**Assessment:**
The closure law is now derived from minimisation of `F_loc` rather than asserted. The overlap-invariant checks integrate `chi_sigma * chi_phi` and `chi_sigma^2` directly and compare to `g_phi*I1_int_check` and `g_phi^2*I2_int_check` (also defined as integrals), so the equality is the integrand-level identity `chi_sigma = g_phi*chi_phi/H` propagating through `int`. Output shows `closure law chi_sigma = g_phi chi_phi/H = 0`, `overlap O = g_phi * I1 (integral form) = 0`, `overlap N_ss = g_phi^2 * I2 (integral form) = 0`. Both engines pass.

### F4 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- `mathematica/...mathematica_audit.wl:88-89` refactors the matched-layer block so `c2Const = FullSimplify[(gPhi*i1Int)^2/(nppInt*gPhi^2*i2Int), ...]` and `gEqConst = FullSimplify[gPhi^2*i1Int/kX, ...]` — i.e., consume the integral-derived `nppInt`, `i1Int`, `i2Int` from the F1 Mathematica block rather than the abstract `i1, i2` via a `constSubs` dictionary.
- The Mathematica gain assertion now compares to `gPhi^2*nppInt/(kX*hw)` (the integral form), not the abstract `npp/(kX*hw)`.

**Assessment:**
The two engines now reach the matched-layer conclusions via different intermediate quantities: the `.py` script still uses the abstract `const_subs = {I1: Npp/Hw, I2: Npp/Hw**2}` on `C2` and `Geq`, while the `.wl` uses concrete `Integrate[]` outputs. The F1 Gaussian-integral block in each engine independently verifies the substitution, so when sympy substitutes and Mathematica integrates, the agreement on the final values (`1`, `g_phi^2*N_pp/(K_X*H_w)`) is now adversarial. The softening block (F2) still has parallel structure across engines but F2's reformulation removed the original `constSubs`/`const_subs` tautology so the parallelism no longer reaches `PASS` via the same tautological path. F4's specific concern (matched-layer block transliteration) is addressed.

## Exec log assessment

**SymPy:** exit log file `redteam/exec_logs/stage_064_sympy.log` is absent. However, `scripts/output/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.txt` was regenerated at 19:48 (after script mtime 19:46). All `expect_zero` calls printed `... = 0` without raising `AssertionError`, and the final `STAGE 47 AUDIT PASSED` banner appears. Notable lines:

- `closure law chi_sigma = g_phi chi_phi/H = 0`
- `overlap O = g_phi * I1 (integral form) = 0`, `overlap N_ss = g_phi^2 * I2 (integral form) = 0`
- `matched-layer I1 reduction = 0`, `matched-layer I2 reduction = 0`
- `Lambda^2/Theta (closure form) = I1**2*g_phi**2/(H_w*I2)`, `Lambda^2/Theta (matched layer) = I1*g_phi**2`
- `equilibrium softening equals g_phi^2 I1 = 0`

**Mathematica:** exit log file `redteam/exec_logs/stage_064_mathematica.log` is also absent. However, `mathematica/output/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.txt` was regenerated at 19:48 (after script mtime 19:47). All `expectZero` calls printed `PASS: ...` and the final `Stage 064 Mathematica audit passed.` line appears. Notable lines:

- `PASS: closure law chi_sigma = g_phi chi_phi/H`
- `PASS: overlap O = g_phi * I1 (integral form)`, `PASS: overlap N_ss = g_phi^2 * I2 (integral form)`
- `PASS: matched-layer I1 reduction`, `PASS: matched-layer I2 reduction`
- `Lambda^2/Theta (closure form) = (gPhi^2*i1^2)/(hw*i2)`, `Lambda^2/Theta (matched layer) = gPhi^2*i1`
- `PASS: equilibrium softening equals g_phi^2 I1`

**Output freshness:** Both `.txt` outputs (mtime 19:48) are newer than their corresponding scripts (mtimes 19:46 and 19:47). Output regeneration is confirmed even though the orchestrator's `.log` capture files are missing.

## Material-change assessment

`material_change`: false.

The edits add new verifications (closure law, integral reductions) and reformulate existing ones (softening, matched-layer C^2/Geq) but do not change any derived constant, sign, or symbolic identity that downstream stages would import. The final printed forms in the unchanged matched-layer block (`C^2 | H=const = 1`, `G_eq | H=const = N_pp*g_phi**2/(H_w*K_X)`) are identical to the pre-edit values; only the path to those values became non-tautological. No downstream concern.

## Side observations (non-blocking)

- `mathematica/...mathematica_audit.wl:88-89` uses `Assumptions -> lInt > 0 && hw > 0 && ...` even though the actual expression (`(gPhi*i1Int)^2/(nppInt*gPhi^2*i2Int)`) involves only `L`-based integrals (defined at line 73-76 with symbol `L`, not `lInt`). The `lInt > 0` assumption is vestigial — harmless because the expression simplifies without it, but the assumption clause has no effect on these specific integrands. Not blocking.
- The sympy block at lines 70-82 (F3 integral overlap, with `L_int`) and the matched-layer block at lines 95-109 (F1 reduction, with `L`) duplicate the Gaussian integration with different dummy names. This is intentional per the directive's "use distinct names like `Osp_int_check`, `Nss_int_check`" guidance to keep the F1 and F3 blocks independent; just noted that both integrals could have been folded into one. Not blocking.

## Verdict justification

All four findings are resolved. F1 and F2 replaced the original `const_subs` tautologies with derivations grounded in (a) concrete-Gaussian integration that independently produces `Npp/Hw = I1` and `Npp/Hw^2 = I2`, and (b) closure-form `Lambda^2/Theta = g_phi^2*I1^2/(Hw*I2)` whose matched-layer reduction now consumes two independently-verified facts. F3 added a real linear-response derivation of `chi_sigma = g_phi*chi_phi/H` via `solve(diff(F_loc, sigma_loc) == 0, sigma_loc)` and integral-level overlap-invariant checks. F4 differentiated the Mathematica matched-layer block to consume integral-derived quantities while the sympy block continues to use the abstract substitution, restoring adversarial cross-engine value. Codex's single deviation (F2's extra `Npp -> I1*Hw` substitution) is necessary and uses an F1-verified identity rather than hand-encoding the answer. Output files were regenerated post-edit even though the orchestrator's `.log` capture files are missing, and both audits printed every `expect_zero`/`expectZero` as zero plus a final "audit passed" banner. No regressions in the diff.
