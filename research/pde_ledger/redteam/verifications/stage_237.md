---
unit_id: 237
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-03T08:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 237

## Per-finding outcomes

### F1 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
New `mathematica/moving_throat_pde_stage237_..._mathematica_audit.wl` (304 lines) covering M1–M7 with a residual-zero helper (`expectZero` → `normalizeExpr` (FullSimplify under `$Assumptions`, ConditionalExpression-stripped) → `fail`/`Exit[1]`).

**Assessment:**
Genuinely independent route, not a transliteration. Different var names (camelCase `epsEta`, `cEtaU`, `KPhiEff`), native primitives (`FullSimplify`, `Series`/`Normal`, `Coefficient`, `Det`, `D`, `PowerExpand`, `Solve` with positivity `$Assumptions`), and a DIFFERENT decomposition at several checks: M2 uses `Solve[staticExp==1, Rratio, Reals]` then compares the solved root vs. the static curve (vs. SymPy's `.subs`); M5 extracts the first-order term via `Coefficient[Normal[Series], t, 1]` (vs. SymPy series-subtract); M7 codim-two uses `Solve[{-R1-cEta E1==0, E1==0}, {R1,E1}]` (a stronger route than SymPy's weak homogeneous `LUsolve([0,0])` that the auditor flagged as A19). Covers the full M1–M7 manifest; 24 PASS / 0 FAIL, exit 0.

### F2 — insufficient_verification (the load-bearing finding)

**Classification:** resolved

**What changed:**
SymPy §5 (lines 241–294) rebuilt: dressing quantities are now `sp.Function`-wrapped over `support_args=(zeta, M_tr, lambda_phi, K_phi_eff)`, `q_eta_support_frame` is built from them, a guard (line 252) asserts `support_args ⊆ q_eta_support_frame.free_symbols`, and a `impose_support_independence` helper zeroes ONLY `Derivative` atoms whose `.expr.func` is a dressing-sector function before the four `assert_zero(diff(...))`. lambda_phi/K_phi_eff derivatives are taken on `q_eta_support_composite` (zeta/M_tr substituted by their explicit support-variable expressions), so the chain path is live. Mathematica M6 mirrors this and adds an explicit `Not[FreeQ[#, Derivative]]` guard (lines 243–246) that FAILS if the chain rule produced no Derivative terms (i.e. trips on the original vacuous case).

**Assessment:**
GENUINELY NON-VACUOUS. The original trap was `diff(EXPR, VAR)` with `VAR ∉ free_symbols` = 0 for any EXPR. Now the support vars ARE reachable (function args + composite substitution), the chain rule produces nonzero Derivative terms, and only the sector-independence Derivative atoms are zeroed — an explicit support-variable factor leaking into a dressing definition is NOT a Derivative atom and survives, tripping the assert. Confirmed via the negative control: if `K_eta_eff_support` carried a `K_phi_eff` factor, `d/dK_phi_eff` retains `-1/K_phi_eff` (not a Derivative atom) → assert fails. The Mathematica `Not[FreeQ[#, Derivative]]` precondition is an even stronger live-channel guard.

### F3 — tautological_check

**Classification:** resolved

**What changed:**
SymPy line 174–177 now asserts `M_packet*[R1,E1] - [q_nt_linear.subs(t,1), q_eta_linear.subs(t,1)]` against the INDEPENDENTLY series-derived linear forms (lines 156–162), not against a hand-expansion of the same product.

**Assessment:**
Substantive. `q_nt_linear`/`q_eta_linear` come from `sp.series` of the actual log expressions, independent of `M_packet`. A wrong `M_packet` entry (e.g. wrong `c_eta`) now trips the assert. Sign/`t`-factor reconciliation checks out (`q_nt_linear.subs(t,1) = -R1-c_eta E1` = first component; `E1` = second).

### F4 — tautological_check

**Classification:** resolved

**What changed:**
SymPy line 310: `qeta_param.subs(qeta_param,0)` (= 0==0) replaced by `q_eta.subs(eps_eta, eps_eta_ref)` on the real `q_eta = log(eps_eta/eps_eta_ref)`.

**Assessment:**
Substantive. Evaluates `log(1)=0` on the actual `q_eta` expression; a mis-defined `q_eta` would fail. Mathematica M7 carries the matching `qEta /. epsEta -> epsEtaRef` endpoint check.

## Exec log assessment

**SymPy:** exit=0. Notable: `dq_eta/dzeta = 0` … `dq_eta/dK_phi_eff = 0` (the post-`impose_support_independence` results print 0, confirming the asserts run on the rule-reduced expressions), and `All Stage 237 symbolic checks passed.`

**Mathematica:** exit=0. 24 PASS / 0 FAIL — matches the M1–M7 manifest count exactly (M1:2, M2:5, M3:4, M4:1, M5:1, M6:4, M7:7). Notable: `PASS: M6 D[q_eta, K_phi_eff]` and `All Stage 237 Mathematica symbolic checks passed.`

**Output freshness:** confirmed. Both `output/*.txt` mtimes are 2026-06-03 08:16, newer than the `.py` (06-02 22:22) and `.wl` (06-02 22:24). Committed outputs match the exec logs (no silent parser-skip).

**FAIL path:** non-vacuous in both engines — SymPy `assert_zero` raises `AssertionError` on nonzero (with a real `_numeric_zero` two-sample fallback); Mathematica `fail` does `Exit[1]`.

## Material-change assessment

`material_change`: false. All four fixes strengthen verification only — no derived constant, sign, or formula consumed by downstream units (238/239/242) changed. The paper formulas were already correct (`paper_alignment: aligned`); F2/F3/F4 changed how claims are checked, and F1 adds a second engine. No downstream re-audit needed on account of unit 237.

## Side observations (non-blocking)

- The SymPy A19-style homogeneous `M_packet.LUsolve([0,0])` codim-two check (line 313) remains weak (the auditor noted it; it was not a directive finding). The new Mathematica M7 covers the same codim-two point with a stronger explicit `Solve`, so the claim is now well-anchored across the two engines. Not blocking.

## Verdict justification

All four findings are `resolved`. The load-bearing F2 fix is genuinely non-vacuous in both engines: the support variables are now reachable through the chain rule (function-wrapped dressing quantities plus composite substitution), only the sector-independence Derivative atoms are zeroed, and an explicit support-variable leak would survive and trip an assert — confirmed by the negative control and, in Mathematica, by an explicit `Not[FreeQ[#, Derivative]]` live-channel guard. F1's `.wl` is an independent decomposition (different primitives, names, and several genuinely different routes), not a transliteration, and covers the full M1–M7 manifest. F3 and F4 are now substantive checks that a wrong value would fail. Both engines exit 0 with PASS counts matching the manifest, outputs are fresh, and the FAIL paths are live. No regressions in the diff; no material change downstream.
