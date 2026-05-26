---
unit_id: 025
batch: II.1
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-26T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 025 (batch II.1 v2)

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- SymPy: `scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py:161-163` — the single tautological line `expect_zero("N0 - P^2/Delta^2", N0 - P**2 / Delta**2)` was replaced with a local rebuild from primitive symbols:
  ```python
  P_raw = OmegaU**2 * GW + R * GU
  Delta_raw = OmegaU**2 * OmegaW**2 - R**2
  expect_zero("N0 reconstructed from raw symbols", N0 - P_raw**2 / Delta_raw**2)
  ```
- Mathematica: `mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl:119-123` — the corresponding `expectZero["N0 - P^2/Delta^2", n0 - p^2/delta^2]` was replaced with a scoped `Module[{pRaw, deltaRaw}, ...]` that rebuilds `pRaw` and `deltaRaw` from `omegaU, omegaW, gU, gW, r` and asserts `expectZero["N0 reconstructed from raw symbols", n0 - pRaw^2/deltaRaw^2]`.

**Assessment:**
The edits match the directive's "Required change" code blocks verbatim in both languages. The replacement is substantively non-tautological. `N0` was constructed at sympy line 80 via `N0 = sp.simplify(P**2 / Delta**2)` where `P` and `Delta` are the cached outputs of `zero_frequency_coefficients()`. The new check reconstructs `P_raw` and `Delta_raw` independently from the declared primitive symbols (`OmegaU, OmegaW, GU, GW, R`) and tests `N0 - P_raw**2/Delta_raw**2 == 0`. If a future edit decoupled the cached `P, Delta` from their primitive definitions, this assertion would fail. Exec logs confirm the new assertion is present and passes: sympy log line 73 (`N0 reconstructed from raw symbols = 0`) and mathematica log lines 68-69 (value `0`, then `PASS: N0 reconstructed from raw symbols`). The old assertion name `N0 - P^2/Delta^2` no longer appears in either log. No collateral edits.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
- SymPy: `scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py:140-148` — inserted a new `II.3 — P = 0 corollary (paper Checks item 4)` block at the end of `normalization_formula()` (after the existing II.2 block, as directed). The block sets `P_zero_sub = {GW: -R * GU / OmegaU**2}`, computes `N0_at_Pzero = sp.simplify(N0.subs(P_zero_sub))` and `residual_at_Pzero = sp.simplify((mhat**2 * P0_compact - target).subs(P_zero_sub))`, then runs `expect_zero("N0 vanishes when P=0", N0_at_Pzero)` followed by `expect_zero("(mhat^2*P0 - target) at P=0 equals -target", residual_at_Pzero + target)`.
- Mathematica: `mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl:102-111` — inserted the parallel block before the closing `];` of `normalizationFormula[]`, scoped in `Module[{pZeroSub, n0AtPzero, residualAtPzero}, ...]` with the same substitution `gW -> -r*gU/omegaU^2` and the same two `expectZero` assertions.

**Assessment:**
Both edits match the directive's "Required change" code blocks verbatim. The substitution `GW = -R*GU/OmegaU**2` forces `P = OmegaU**2 * GW + R*GU = 0` symbolically, so `N0 = P^2/Delta^2 = 0` and `mhat^2 * P0_compact = 0`, giving `residual = -target`. Adding `target` then yields 0. Both checks are non-tautological:
- `N0_at_Pzero` is computed by substituting into the *cached* `N0` (which depends on the cached `P`), so it tests that the substitution propagates correctly through `N0`'s construction. If `N0` had been built from something other than `P^2/Delta^2`, this check could detect a mismatch.
- The `residual + target == 0` check is non-trivial: it verifies (a) the LHS of the target equation vanishes when P=0 (because `P0_compact` has `P^2` in its numerator) and (b) `target` is the constant the script declared earlier. Any decoupling of `target` from `54 G c_s^5 / (5 a^5 c^5)` would survive the `+ target` cancellation only if `residual = -target` as an identity, which would not in general be true.

Exec logs confirm both new assertions pass:
- sympy log lines 51-56: `II.3 — P = 0 corollary (paper Checks item 4)`, `N0 at P=0 = 0`, `N0 vanishes when P=0 = 0`, `(mhat^2*P0 - target) at P=0 = -54*G*c_s**5/(5*a**5*c**5)`, `(mhat^2*P0 - target) at P=0 equals -target = 0`.
- mathematica log lines 43-50: `II.3: P = 0 corollary (paper Checks item 4)`, `N0 at P=0 = 0`, `PASS: N0 vanishes when P=0`, `(mhat^2*P0 - target) at P=0 = (-54*cs^5*gConst)/(5*a^5*cSpeed^5)`, `PASS: (mhat^2*P0 - target) at P=0 equals -target`.

The printed residual value `-54*G*c_s**5/(5*a**5*c**5)` (sympy) / `(-54*cs^5*gConst)/(5*a^5*cSpeed^5)` (mathematica) is genuinely `-target`, not zero, so the assertion is meaningful. No collateral edits.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `II.3 — P = 0 corollary (paper Checks item 4)` (line 51) — new subsection appears as directed.
- `N0 vanishes when P=0 = 0` (line 54) — new F2 assertion passes (`expect_zero` prints the simplified value then succeeds because non-zero would raise `AssertionError`).
- `(mhat^2*P0 - target) at P=0 = -54*G*c_s**5/(5*a**5*c**5)` (line 55) — confirms residual = -target, the non-trivial part of F2's check.
- `N0 reconstructed from raw symbols = 0` (line 73) — F1's new substantive check passes.

**Mathematica:** exit=0. Notable lines:
- `II.3: P = 0 corollary (paper Checks item 4)` (line 43) — new subsection appears.
- `PASS: N0 vanishes when P=0` (line 47) and `PASS: (mhat^2*P0 - target) at P=0 equals -target` (line 50) — both new F2 assertions pass.
- `PASS: N0 reconstructed from raw symbols` (line 69) — F1's new substantive check passes.
- The old `N0 - P^2/Delta^2` assertion name appears in neither log, confirming the tautological check was removed in both engines.
- Trailing log line: `Stage 8 Mathematica audit passed.` followed by `# exit_code: 0`.

**Output freshness:** The saved `.txt` outputs under `scripts/output/` and `mathematica/output/` have mtime 2026-05-21 16:58 — older than the edited scripts (mtime 2026-05-25 23:34). They were NOT regenerated post-fix. However, the orchestrator's exec logs at `redteam/exec_logs/stage_025_{sympy,mathematica}.log` are fresh (mtime 2026-05-26 00:14) and contain the post-edit transcripts with the new assertions. Verification is based on those fresh exec logs, as the prompt directs.

## Material-change assessment

`material_change`: false.

Both edits are script-internal verification additions/replacements. F1 replaces one tautological identity (`N0 - P^2/Delta^2`) with a substantively equivalent identity (`N0 - P_raw^2/Delta_raw^2`) that uses primitive-symbol rebuilds. The residual still evaluates to 0; the symbolic expressions consumed by downstream stages (`Delta, Q, P, N0, D0, P0, P0_compact, target, mhat^2`, and the derivatives `dP_0/dK, dP_0/dX`) are unchanged. F2 adds new corollary checks at the `P=0` slice without altering any derived expression. No downstream unit could observe a different result.

## Side observations (non-blocking)

- The saved `.txt` outputs under `scripts/output/` and `mathematica/output/` were not regenerated after Codex's edits and now lag the scripts by one iteration. The fresh `redteam/exec_logs/` files are authoritative for this verification, but a later orchestrator pass that compares saved outputs to current scripts will see staleness.
- The previous v1 verification of stage 025 (different batch / different finding set with 8 findings) is overwritten by this v2 verification. The orchestrator's batch tracking should reflect that the v2 directive (2 findings, both resolved) is the current state.

## Verdict justification

Both findings were applied exactly as the directive prescribed — the new lines match the directive's code blocks verbatim, the old tautological assertion is removed, no collateral edits appear in `stage_025_diff.patch`, and both engines' exec logs exit 0 with the new assertion names reported as `PASS` (Mathematica) or as `... = 0` followed by no raised exception (SymPy). The F1 replacement is non-tautological (rebuilds N0 from declared primitive symbols rather than the cached aliases). The F2 additions are non-tautological (the residual evaluates to `-target`, a non-zero expression, before the `+ target` cancellation). No regressions. Verdict: `verified`.
