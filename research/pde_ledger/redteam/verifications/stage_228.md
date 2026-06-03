---
unit_id: 228
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 228

## Per-finding outcomes

### F1 — missing_verification_script (subtype: missing_mathematica)

**Classification:** resolved

**What changed:**
New independent Mathematica audit created at `mathematica/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_mathematica_audit.wl` (379 lines), covering M1–M9. The file is untracked (new), which is why `stage_228_diff.patch` is empty (the orchestrator diff capture only sees tracked-file changes); the script itself is present at the Target path and its committed output `.txt` exists.

**Assessment:** Correct and genuinely independent. The `.wl` re-derives every claim through a different primitive family and decomposition than the SymPy script — it is NOT a transliteration:
- Pole census (M7): `NSolve[poly==0, y, WorkingPrecision->60]` on the quartic in `y=omega^2` (wl:246), then `Select`/`Sort` for positive real roots — NOT the forbidden SymPy `np.roots`-on-companion-matrix route.
- Linear responses (M1–M3): `Series[...,{eps,0,1}]` + `Coefficient[...,eps,1]` (wl:97-100) and `Coefficient[Expand[...], var]` for the rows (wl:136) — series-truncation, not SymPy's `sp.diff(...,eps).subs(eps,0)`.
- Rank/nullity (M4): native `MatrixRank`/`NullSpace` (wl:176-183).
- Reduced determinant (M5): `Det[...]` + `Factor`/`FullSimplify` on `piCoeff.transferBasisColumns` (wl:195-199).
- Slopes (M8): analytic implicit differentiation `dy = -(F_eps/F_y)` plus `D[Log[rq], ...]` chain rule (wl:288-301) — method-independent of SymPy's finite-difference slopes (and stronger).
- Generators (M6): `NullSpace` first vector, oriented by `Xi_1>0`, Euclidean-normalized (wl:208-215).

All 9 manifest items map to live, TrueQ-gated assertions. Exec log exits 0 with all checks PASS. M1–M9 are genuine and non-tautological.

### F2 — paper_misalignment (subtype: notes_contradicts_script)

**Classification:** resolved

**What changed:**
Notes lines 151-152 now read `+\frac{196\pi^2}{98\pi^2-25}x_{\Omega_U}` / `x_{\Omega_W}` (previously `247\pi^2`). No script change (per the user-authorized RESOLVED block: script:194/output:14 authoritative).

**Assessment:** Correct. The authorized direction was notes→`196`. No residual `247` anywhere in the notes (grep confirms zero `247` matches; the only `215` matches are the unrelated numeric `0.69293215`). Cross-engine corroboration: the new `.wl` independently computes the `delta_1` row coefficient as `196 Pi^2/(98 Pi^2-25)` (log line 10, and M3 row PASS at log:17-18), matching the SymPy value and confirming the notes' `247` was the typo.

### F3 — paper_misalignment (subtype: notes_contradicts_script)

**Classification:** resolved

**What changed:**
Notes line 196 now reads `\frac{196(200+147\pi^2)(80000+343225\pi^2+43218\pi^4)}{475(...)}` (previously `247(251+215\pi^2)`). No script change.

**Assessment:** Correct. Stale `247`/`251`/`215` factors gone; `196`/`200`/`147` present. Cross-engine corroboration: the new `.wl` independently computes the reduced determinant as `196*(200+147*Pi^2)*(80000+343225*Pi^2+43218*Pi^4)/(475*(8670000+14894275*Pi^2+2117682*Pi^4))` (log line 36, M5 PASS at log:37-40) via native `Det`/`Factor` — matching the SymPy form and confirming the notes' `247(251+215 Pi^2)` was the typo. This is strong independent evidence the F2/F3 typo-direction was right.

## Exec log assessment

**SymPy:** exit=0. Re-run (timestamp 17:25:22) reproduces `delta_1 = ... + 196*pi**2*xOU/(...) + 196*pi**2*xOW/(...)` (log:11) and `det[...] = 196*(200 + 147*pi**2)*(...)` (log:18). Consistent with the authoritative script values.

**Mathematica:** exit=0. 51 PASS, 0 FAIL. Headline lines:
- `M1 split identity = 0` / `PASS` (log:11-12).
- `delta_1 row = {0, 0, 50/(25 - 98*Pi^2), 2 + 50/(-25 + 98*Pi^2), ...}` → `M3 exact delta_1 row = 0 / PASS` (`2 + 50/(-25+98Pi^2)` is algebraically `196 Pi^2/(98 Pi^2-25)`; the M3 expectZero of the squared-difference against the literal `196 Pi^2/...` target confirms equality).
- `det[...] = (196*(200 + 147*Pi^2)*...)/(475*(...))` → `M5 reduced determinant = 0 / PASS` and `M5 ... nonzero = True / PASS` (log:36-40).
- M9 decisive inequalities: `M9 numerator both dynamic exceeds transported static = True / PASS`, `M9 denominator both dynamic exceeds transported static = True / PASS`, `M9 denominator nonempty dynamic exceeds transported static = True / PASS` (log:125-130).

**M9 `dynamic > static` genuineness:** Confirmed a genuine computed comparison. `numDynamic`/`denDynamic` come from `finiteDynamicCeilings[rqBase, ..RqSlopes, threshold]` (wl:335-349, `Log[rqVals/threshold]/(-slope)`), and `numStatic`/`denStatic` from `budget/xiNum`,`budget/xiDen` (wl:357-358). The three `expectTrue[..., a > b]` calls (wl:373-375) compare two independently-computed quantities and are `TrueQ`-gated — not an assumed/echoed relation. The numerator-nonempty case is correctly excluded as `Infinity` (one wall pole improves; `M9 numerator nonempty dynamic is Infinity = True / PASS`), exactly as the SymPy script does.

**Token-convention non-vacuity:** Confirmed genuine. The `fail[name,...]` helper executes `Exit[1]` (wl:8-12). Every `expect*` helper's negative branch routes to `fail`: `If[TrueQ[res===0], pass, fail]` (wl:26), `If[TrueQ[res], pass, fail]` (32), `If[TrueQ[diff<tol], pass, fail]` (38), `If[TrueQ[mag<tol], pass, fail]` (44), `If[TrueQ[maxDiff<tol], pass, fail]` (51). The `TrueQ` wrapper means an unevaluated/symbolic comparison yields `False` → `fail` → `Exit[1]`, so a non-decisive result cannot produce a spurious PASS. Check count fully reconciles: 56 `expect*` occurrences − 5 helper definitions = 51 real invocations = 51 PASS tokens; no assertion was silently skipped.

**Output freshness:** Confirmed. Both outputs re-generated post-fix at 2026-06-02 17:26:42, newer than `.wl` (17:20:02) and `.py` (2026-05-11). SymPy output `.txt` is the only tracked file modified for this stage besides the notes.

## Material-change assessment

`material_change`: **false**. No script logic or derived result changed: F1 adds a new second engine that reproduces (does not alter) the existing SymPy results, and F2/F3 are notes-only typo corrections that bring the notes into agreement with the already-authoritative script. No downstream unit's inputs change. The cross-engine agreement (196 / 196(200+147π²), the pole census, the slopes, and the decisive ceilings) strengthens confidence without modifying any carried-forward quantity.

## Side observations (non-blocking)

- The dead `*0` scaffolding flagged in the original report (SymPy py:268) is not present in the Mathematica route; the Mathematica orientation uses `unitPositiveXi` (wl:208-215) and re-asserts the `Xi_1>0` sign cleanly. No false-pass risk.
- The static budget literals (`0.367930328492646`, `0.737619063660757`, threshold `21.854566296358396`) in the `.wl` (wl:351-353) are genuine upstream carry-forwards: the two budgets match sibling stage 227 (`scripts/...stage227..._sympy_audit.py:271-272`) verbatim. They are declared inputs to the M9 comparison, consistent with the directive.

## Verdict justification

All three findings are `resolved`. F1's new `.wl` is a genuine, independent re-derivation (NSolve quartic, Series+Coefficient, native MatrixRank/NullSpace/Det, analytic implicit-function slopes) of all nine manifest items, not a transliteration; it exits 0 with 51/51 checks passing under a non-vacuous Exit[1]-on-FAIL token convention. The decisive M9 `dynamic > static` inequality is a true computed comparison of independently-derived ceilings. F2/F3 notes typos are corrected to the script-authoritative `196`/`200`/`147` with no stale `247`/`215` residue, and the new Mathematica engine independently corroborates both values — strong cross-engine evidence the typo direction was correct. No regressions, no material change.
