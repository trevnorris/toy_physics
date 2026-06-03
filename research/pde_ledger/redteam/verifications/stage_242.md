---
unit_id: 242
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-03T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 242 (checkpoint)

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
The whole `.wl` was re-authored (`mathematica/...stage242...mathematica_audit.wl`, diff @@ -56,248 +56,268). All three idiosyncratic SymPy idioms the auditor flagged are gone:
- The artificial "above/below bound" term-split (`sigmaSel - (1/varrho-2) - 1/(3 varrho)`) is **removed entirely**; window membership now via `Resolve[ForAll[...], Reals]` (wl:154-165).
- The abstract-zeta-function device (`RtrSbFn[zeta]` + derivative-substitution rules) is **removed**; q-packet support-blindness now via `D[closedObservablePacket, zeta]` on the actual closed forms (wl:202-212).
- The `Exp[t·d]` first-order parametrization is **removed**; infinitesimal objects now via a `logDrift[expr, vars, drifts]` total-log-differential helper `Sum_v v·D[Log[expr],v]·dv` (wl:89-92, used at 216/227/240).

The placement leg is also independent: the `.wl` builds `sigmaFromRho` and `sigmaFromEpsilon` as two separate closed forms and cross-checks them (wl:109-129), rather than asserting the `.py`'s single residual. Section structure renumbered I–VII → M1–M7, variable names renamed (epsilonSelected/trailObservable/targetObservable). The 3×3 orbit matrix is retained, which the directive explicitly permitted (native `LinearSolve`/`Inverse`), and even there the route differs (`LinearSolve` vs `.inv()*`, packet entries derived from definitions not aliased).

**Assessment:** Genuinely independent now, not a second transliteration. An error in the SymPy choreography (term-split, abstract-fn, `Exp[t·d]`) would no longer be silently echoed because the `.wl` uses different primitives and decompositions. Exits 0, 24 PASS, no FAIL.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
SymPy (py:93-97) adds `ratio = sp.nsimplify(Pi_tr/C_mix)` (= 4/3, `C_mix` cancels), then `assert ratio == 4/3` and `assert ratio > 1 and ratio < 2`, printing `[ok] selected branch strictly inside lowest symmetric twin window`. The two original gap-equalities are kept (documenting the gap form) per directive. The `.wl` (wl:152-165) asserts `demandRatio - 4/3 == 0` (demandRatio derived from `traceLoad/mixedCapacity` on the real closed forms) and certifies `1 < (4/3·X)/X < 2` for all `lambdaWin>0, 0<epsilonWin<1` via `Resolve`, gated by `expectTrue`.

**Assessment (F2_window_nonvacuous = yes):** The strict inclusion `C_mix < Pi_tr < 2 C_mix` is now genuinely exercised, not a relabeled equality. The SymPy `> 1 and < 2` is a concrete bool that trips if the ratio were set to 1 or 2. The `.wl` `Resolve` returns `False` (→ `fail` → `Exit[1]`) if the hardcoded factor were 2; and the separate `demandRatio == 4/3` ties the literal back to the actual `traceLoad = (4/3)·mixedCapacity` construction. A value at/outside the window trips an assert.

### F3 — tautological_check

**Classification:** resolved

**What changed:**
SymPy (py:213-218) replaces `Theta_1 = dln_Rtr` with `Sigma_tr = (1+deltaU)dchi0 + (1+chi0)ddeltaU`, `C_tr_star = chi0·deltaU/((1+chi0)(1+deltaU)(1+chi0+deltaU))`, `Theta_1 = simplify(-C_tr_star·Sigma_tr)`, then `assert_zero("Theta_1 packet form matches dln R_tr", Theta_1 - dln_Rtr)`. The `.wl` mirrors this with `thetaFromPacket = -trailCoefficient·trailSigma` vs `trailLogDrift` (wl:253-258, 268).

**Assessment:** No longer `x - x == 0`. `Theta_1`/`thetaFromPacket` is built from the packet definition, fully independent of `dln_Rtr`/`trailLogDrift` (which come from the differential route at py:189 / wl:227). The residual is a genuine cancellation that would fail if `C_tr_star` were mis-stated (e.g. dropping the `(1+chi0+deltaU)` factor). `Xi_1`/`R_1` checks (real cancellations) were correctly left intact, and `Theta_1` remains consistent in the orbit-packet matrix block.

### F4 — stale_output

**Classification:** resolved

**What changed:** `wl:59` banner now `STAGE 242 — ...` (diff line 10); regenerated `.txt:8` shows `STAGE 242`. In-file footer (wl:342) reads `All Stage 242 ... passed.`

**Assessment:** Grep over the full `.wl` finds no remaining `22x` label. This file only, no batch renumber. Confirmed.

## Exec log assessment

**SymPy:** exit=0. 36 `[ok]` lines, 0 FAIL. Notable: `[ok] selected branch strictly inside lowest symmetric twin window` (F2), `[ok] Theta_1 packet form matches dln R_tr` (F3), `All Stage 242 symbolic checks passed.`

**Mathematica:** exit=0. 24 PASS, 0 FAIL. Banner `STAGE 242` (F4). Notable: `M3 strict twin-window inclusion by Resolve = True` / `PASS` (F1+F2), `M6 Theta_1 packet form equals dln R_tr = 0` / `PASS` (F3), `M5 closed-form q-packet zeta derivatives` (independent route, F1).

Call-site reconciliation: `.wl` has 26 `expectZero`/`expectTrue` tokens − 2 helper definitions (wl:48, 61) = 24 call sites = 24 PASS; no silent parser-skip. Matrix/vector checks (`closedObservablePacket`, `orbit compiler maps`, `LinearSolve`, `formal inverse map`) each emit one MatrixForm + one `allZeroQ`-gated PASS that Flattens over all entries, so a non-zero entry would `fail`→`Exit[1]`. FAIL path non-vacuous on both engines.

**Output freshness:** confirmed. sympy.py mtime 1780497979 < sympy.txt 1780498657; wl mtime 1780498108 < wl.txt 1780498657. Both `.txt` regenerated post-fix.

## Material-change assessment

`material_change`: false. All four edits add or strengthen checks (F2 strict inequality, F3 substantive packet-form derivation) and re-author the second engine (F1) or correct a label (F4). No derived constant changed: the placement coordinate `2(1-eps)/3`, demand ratio `4/3`, threshold rewrites, `R_target` bridge, the compiler drifts, and the product law `R_target·M_mix = C_mix` are all identical to pre-fix. The sample point is unchanged (eps=17/55, varrho=76/165, sigma=17/19, R_tr=19/25, R_target=2749376/54071875). Nothing downstream depends on these script internals beyond the already-stated results, which are preserved. No downstream staleness.

## Side observations (non-blocking)

- The `.wl` `Resolve` window certificate hardcodes the `4/3` factor inside the `ForAll` expression (numerator `(4/3)·X` over denominator `X`), so it literally certifies `1 < 4/3 < 2`. This is acceptable because the separate `demandRatio - 4/3 == 0` check derives the `4/3` from the actual `traceLoad/mixedCapacity` closed forms, anchoring the literal to the construction. The two together exercise the deliverable; not a concern.
- The SymPy script still uses the `Exp[t·d]` route for the infinitesimal compilers — but the auditor scoped F1 to the `.wl` only (the `.py` was deemed faithful), and the directive only required the `.wl` to abandon that idiom for independence. The `.py` `Exp[t·d]` route is the first engine and is fine; the `.wl` now provides the independent total-log-differential cross-check. Not a finding.

## Verdict justification

All four findings are `resolved`. The re-authored `.wl` is genuinely independent — the three flagged transliteration idioms (artificial term-split, abstract-zeta-function support-blind device, `Exp[t·d]` parametrization) are all removed and replaced with native, differently-decomposed routes (`Resolve`, closed-form `D[...,zeta]`, `logDrift` total differential), restoring the checkpoint's dual-engine guarantee. The strict twin-window inclusion is now exercised on both engines and would trip if the branch were at/outside the window. The `Theta_1` check is substantive (packet-form LHS independent of `dln_Rtr`). The banner reads 242 with no residual 22x labels. Both engines exit 0 with non-vacuous FAIL paths, pass counts reconcile exactly (sympy 36, mathematica 24), and outputs are fresh. `material_change: false` — no derived result changed.
