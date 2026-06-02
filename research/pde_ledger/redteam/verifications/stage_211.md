---
unit_id: 211
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

# Verification — unit 211

## Per-finding outcomes

### F1 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
Codex created the new Mathematica audit at the registered path
`mathematica/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl`
(untracked, so it does not appear in the captured diff; confirmed present on disk and exec-run to exit 0). The script defines `Phi`, `tau`, the boxed numerators, and verifies the full claim manifest M1–M6 with hard-failing guards (`expectZero`/`expectTrue`, both routed through `fail[...] -> Exit[1]`).

**Assessment:**
The edit is correct, and — critically — it is a GENUINELY INDEPENDENT derivation, not a transliteration of the SymPy `.py`. Three corresponding sections confirm a different decomposition:

1. **M1 stationary numerator (the load-bearing test).** The `.py` only defines `Nr,Ns` from the boxed forms and checks `diff(Phi,r) - Nr/(...) == 0` (`.py:53-58,73-76`); it never extracts the numerator of the differentiated `Phi`. The `.wl` adds exactly the directive-mandated extra route: `directNumerator[var_] := FullSimplify[Numerator[Together[D[Phi, var]]], ...]` (`.wl:82-85`), then checks `directNumerator[r] - Nr == 0` and `directNumerator[s] - Ns == 0` (`.wl:96-97`) BEFORE the derivative-law check (`.wl:98-99`). The log prints the directly-extracted numerators (mathematica log lines 14-15) and matches them to the paper `N_*` (log 16-19). This is "differentiate Phi, common-denominator via Together, read off the numerator, equate to the paper's claimed N_*" — the different decomposition the anti-transliteration guard required.

2. **Total degrees (M2/M3).** The `.py` uses `sp.Poly(Ccross, r, s).total_degree()` (`.py:98-104`). The `.wl` uses a CoefficientRules route: `totalDegreeRS[poly_] := Max[Total /@ CoefficientRules[Expand[poly], {r, s}][[All, 1]]]` (`.wl:49-53`), applied at `.wl:104,114-115`. Not a port of `Poly.total_degree()`.

3. **M4 Bezout.** Both engines compute the product from the MEASURED degrees, not hardcoded `4*6`: `.py:167` `poly_C.total_degree() * poly_Sr.total_degree()`; `.wl:126` `bezoutBound = crossDegree srDegree`. Non-tautological in both.

Manifest coverage in the `.wl`: M1 (`:96-99`), M2 identity + degree==4 (`:106-108`), M3 two identities + two degree==6 (`:117-122`), M4 product==24 (`:126-128`), M5 Delta_iso + tau_iso reductions (`:147-148`), M6 N_r(1,1)=N_s(1,1)=0 (`:166-167`). Every manifest item is an explicit, hard-guarded, non-tautological check. The script uses native Mathematica primitives (`D[]`, `Together`, `Numerator`, `FullSimplify`, `CoefficientRules`, `/.`) throughout. Resolved.

### F2 — stale_output (cosmetic stage-number banner)

**Classification:** resolved

**What changed:**
`scripts/...sympy_audit.py:35` banner string changed from `"STAGE 194 — ..."` to `"STAGE 211 — ..."`. The captured diff (`stage_211_diff.patch`) shows exactly one hunk: line 35 banner swap, nothing else. No assertion, symbol, or math line touched.

**Assessment:**
Correct and minimal. The SymPy exec log line 8 now reads `STAGE 211 — FULL INTERIOR SIMPLEX OPTIMIZER AND FINITE ALGEBRAIC CANDIDATE REDUCTION`, and the closing `All Stage 211 identities verified.` is preserved (log line 234). Script exits 0. Resolved.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 211 — FULL INTERIOR SIMPLEX OPTIMIZER AND FINITE ALGEBRAIC CANDIDATE REDUCTION` (log:8) — banner now corrected.
- `exact dPhi/dr compiler = 0`, `exact dPhi/ds compiler = 0` (log:42-43).
- `deg C_cross = 4`, `deg S_r = 6`, `deg S_s = 6` (log:213-215); `Bezout bound ... = 24` (log:232).
- `All Stage 211 identities verified.` (log:234), `# exit_code: 0` (log:235).

**Mathematica:** exit=0, 15 PASS / 0 FAIL. Notable lines:
- `M1 D[Phi,r] numerator minus paper N_r = 0` → `PASS` (log:16-17); same for `s` (log:18-19); derivative laws PASS (log:20-23).
- `M2 C_cross has total degree 4 = True` → PASS (log:31-32); M3 sextic degrees both PASS (log:43-46).
- `M4 computed Bezout product = 24` → PASS (log:51-53).
- M5 Delta_iso/tau_iso reductions PASS (log:58-61); M6 symmetric N_r(1,1)/N_s(1,1) PASS (log:66-69).
- `All Stage 211 identities verified.` (log:71), `# exit_code: 0` (log:72). PASS count = 15 (4+2+4+1+2+2).

**Output freshness:** confirmed. SymPy `.txt` mtime 2026-06-02 10:54:50 > `.py` mtime 10:44:38. Mathematica `.txt` mtime 2026-06-02 10:54 > `.wl` mtime 10:46:09. Both saved outputs were regenerated post-fix.

## Material-change assessment

`material_change`: false. F1 ADDS an independent second-engine verification (no result changed). F2 changed only a `print` string. No derived quantity that downstream units could depend on was altered.

## Side observations (non-blocking)

- The `.wl` deviates harmlessly from the `.py` symbol names (`Ccoef/Dcoef/Ecoef/Fcoef` vs `C/D/E/F`) to avoid Mathematica's protected `C`, `D`, `E`, `F` — this is good practice and reinforces non-transliteration, not a defect.
- M5's `tau_iso` check correctly runs under tightened assumptions (`DeltaIsoExpected > 0 && krs^2 - 2 H0 u > 0`, `.wl:145,148`) so the nested square root resolves to the intended branch; the residual is `0`, consistent with the `.py`. No action needed.

## Verdict justification

Both findings are `resolved`. F2 is a one-line cosmetic banner fix verified against the diff (no collateral edits). F1 is a substantive new Mathematica engine that independently re-derives M1–M6 via a genuinely different decomposition (numerator-extraction from `D[Phi,var]` for the stationary law, `CoefficientRules`-based total degrees, measured-degree Bezout product), passing all 15 hard guards with exit 0. Both engines exit 0, saved outputs are fresh, and no result downstream units depend on was changed. Verdict: verified.
