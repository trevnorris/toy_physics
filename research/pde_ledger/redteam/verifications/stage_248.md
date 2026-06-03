---
unit_id: 248
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-03T12:55:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 248

CHECKPOINT stage → higher bar applied. Four applied findings (F0 user-approved notes,
F1 insufficient_verification, F2 mathematica_transliteration via an iter-2 satisfaction
reframe, F3 stale banner). Both engines exit 0; 30 PASS lines in the Mathematica log,
0 FAIL; 27 SymPy `assert`s, log clean. Read the code in both engines; did not execute.

## Per-finding outcomes

### F0 — notes correction (USER-APPROVED)

**Classification:** resolved

**What changed:** `notes/.../stage248_..._sympy_audit.md:506` `\times 168\%` → `\times 100\%`.

**Assessment:** Confirmed ONLY notes:506 changed. The boxed line now reads
`(T_new/T_Coul - 1) × 100%` with the next line `≈ 23.3128%`. No `168` remains anywhere
in the notes file. The published card `paper/stages/stage_248.tex` has no `168` (untouched).
The SymPy script's multiplier `improve_pct = 100.0*(ratio_num-1.0)` (py:259) and its
assert `abs(improve_pct - 23.3128) < 1e-3` (py:307) are intact — script unchanged for F0.
Working-tree diff for this unit contains only the .wl, .py, their two output .txt, and this
one notes .md; no other prose touched. Resolved.

### F1 — insufficient_verification (dynamic diagnostics identity)

**Classification:** resolved

**What changed:** SymPy py:181-200 and Mathematica wl:196-206 add the `lambda_th` identity
check. Both build `lambda_th_def = V(r_+)/V'(r_+)` and `lambda_th_E = E/V'(r_+)` with `V` an
undefined function head (`sp.Function("V")` py:45; `Vfun` an unassigned head in the .wl),
`r_+ = r_plus(E_turn)` / `rpFun[Eturn]`. They compute the raw pre-substitution gap, then
substitute the turning condition `V(r_+) → E`.

**Assessment:** Non-vacuous in BOTH engines and free of the variable-independence trap.
SymPy: `lambda_th_gap_raw` prints as `(-E_turn + V(r_plus(E_turn)))/Derivative(...)`, a
genuinely nonzero symbolic object, so `assert lambda_th_gap_raw != 0` (structural inequality
to S.Zero) holds; after the V=E substitution `lambda_th_gap == 0` (log: "lambda_th gap (V=E) = 0").
Mathematica: `lambdaThGapRaw` prints `(-Eturn + Vfun[rpFun[Eturn]])/Derivative[1][Vfun][...]`,
the non-vacuity guard `If[TrueQ[lambdaThGapRaw == 0], fail, pass]` passes (gap is NOT zero),
and `expectZero["lambda_th identity under V(r_+)=E", lambdaThGap]` = 0. Because `V` is
undefined, `V'(r_+)` is a real symbolic object (not identically zero), so the check is not the
diff-w.r.t.-absent-variable self-test trap. The `Xi_turn` line is a definitional presence
marker (`Xi1(r_+)`/`Xi1fun[rp]-Xi1fun[rp]`), trivially 0 by construction — as the directive
explicitly intended (anchor only; the load-bearing check is the lambda_th pair). Resolved.

### F2 — mathematica_transliteration (CHECKPOINT higher bar; iter-2 reframe)

**Classification:** resolved

**What changed:** Section II of the .wl (wl:93-136) was rewritten. The
`vcritNewSolved`/`vcontactCoulSolved` derivations via `First[Solve[...]] /. ConditionalExpression`
(the SymPy minus-sign-branch mirror) and the withdrawn `Reduce`/`ToRules` attempt are GONE.
In their place: a satisfaction route. The hand-written compiler forms `vcritNew = Sqrt[2(Vpeak-V0)/ms]`
and `vcontactCoul = Sqrt[2(1/rContact-1/r0)/ms]` are substituted into their defining launch-energy
equations and the residuals `EAtVcrit - Vpeak` and `EAtVcontact - 1/rContact` are shown to
`FullSimplify` to 0 (wl:121,123). A non-vacuity guard checks the pre-substitution gaps
`deltaNew = ElaunchNew - Vpeak` and `deltaCoul = ElaunchCoul - 1/rContact` are NOT identically
zero (with an explicit `FreeQ[..., v0]` guard, wl:108-119). A positive-branch guard confirms
`vcritNew > 0` and `vcontactCoul > 0` (wl:125-136) under local physical assumptions
`thresholdAssumptions = $Assumptions && Vpeak > V0 && r0 > rContact` (wl:95).

**Assessment:** Checkpoint higher-bar independence is GENUINELY met. The route does NOT re-run
Solve/Reduce to re-derive the speeds — it independently verifies each hand-written closed form
satisfies its own defining equation via native substitution + FullSimplify, a primitive
distinct from SymPy's `solve(...)`-then-compare. All five iter-2 acceptance criteria hold:
(1) the `Solve/ConditionalExpression` derivations and the `Reduce/ToRules` attempt are gone;
(2) substitution checks present and = 0 in the log; (3) non-vacuity guards present, pre-sub
gaps genuinely depend on `v0` (log: "delta_new pre-substitution = (ms*v0^2)/2 + V0 - Vpeak",
"delta_coul ... = r0^(-1) - rContact^(-1) + (ms*v0^2)/2"), so a wrong closed form would make
the satisfaction check fail; (4) positive-branch guards pass (log: both "= True"), and they
are non-vacuous — without `Vpeak > V0` / `r0 > rContact` the Sqrt sign would be indeterminate,
so the local thresholdAssumptions do real work; (5) substitution+FullSimplify only, no Solve
re-run, Sections I/III/IV/V unchanged, exit 0. The two `Solve` calls remaining at wl:174/181
are the pre-existing Section III transport-law machinery (`Derivative[1][rpFun][Eturn] -> drPlus`),
explicitly out of F2 scope and untouched by the diff. The deviation (using the local domain
inequalities) is justified and disclosed: global `$Assumptions` do not imply the Sqrt branches
are positive. Resolved.

### F3 — stale banner label

**Classification:** resolved

**What changed:** wl:65 banner `"STAGE 231 — ..."` → `"STAGE 248 — DYNAMIC EVENT-CHAIN COMPILER
FROM THE RELAXED STATIONARY BARRIER FRONT END"`.

**Assessment:** Confirmed in source (wl:65) and in the regenerated output header (log line 8
reads "STAGE 248"). Resolved.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `lambda_th gap (raw)   = (-E_turn + V(r_plus(E_turn)))/Derivative(V(r_plus(E_turn)), r_plus(E_turn))` (non-vacuous)
- `lambda_th gap (V=E)   = 0`
- `Xi_turn symbolic      = Xi1(r_plus(E_turn))`
- `All symbolic and numerical checks passed.`

**Mathematica:** exit=0, 30 PASS, 0 FAIL. Notable lines:
- `delta_new pre-substitution = (ms*v0^2)/2 + V0 - Vpeak` → `PASS: delta_new non-vacuous gap`
- `PASS: v_crit,new satisfies defining energy`, `PASS: v_contact,coul satisfies defining energy`
- `v_crit,new positive branch = True` / `v_contact,coul positive branch = True` → both PASS
- `lambda_th gap (raw) = (-Eturn + Vfun[rpFun[Eturn]])/Derivative[1][Vfun][rpFun[Eturn]]` → `PASS: lambda_th identity non-vacuous guard`
- `lambda_th identity under V(r_+)=E = 0` → PASS
- Header reads "STAGE 248". One benign `Limit::alimv` warning in Section III (pre-existing, not introduced by the diff).

**Output freshness:** confirmed. Both .txt outputs (mtime 2026-06-03 12:33:32) are newer than
their scripts (.py 10:12:32, .wl 10:21:22). Outputs regenerated post-fix.

## Material-change assessment

`material_change`: false. No derived result changed. The re-derivation in Section II lands on
the same compiler forms (`vcritNew`, `vcontactCoul`); only the verification surface changed
(solver-echo → satisfaction route). F1 added a new diagnostic identity check exercising an
already-boxed card relation (no new derived quantity). F3 is a cosmetic label. F0 fixes a
notes-only prose typo (the published card and the script multiplier were already correct).
Downstream units (e.g. Stage 249 helicity compiler consuming `lambda_th`) are unaffected;
the boxed values are unchanged.

## Side observations (non-blocking)

- The `Xi_turn` markers in both engines (`Xi_turn_sym == Xi1(r_plus_E)` / `Xi1fun[rp]-Xi1fun[rp]`)
  are tautological by construction. This is by design per the directive (definitional presence
  anchor); the substantive new content is the `lambda_th` identity pair, which is non-vacuous.
  Not a basis to fail verification.
- The `Limit::alimv` warning in Section III is pre-existing (transport/Coulomb endpoint limit),
  not introduced by these edits, and the affected checks still pass at 0.

## Verdict justification

All four findings are resolved with no regressions. F1's lambda_th identity is non-vacuous and
trap-free in both engines (undefined `V`, nonzero raw gap collapsing to 0 only under the turning
condition). F2 meets the checkpoint higher-bar independence requirement: the SymPy `Solve`/
`ConditionalExpression` mirror is gone and replaced by a genuinely independent native
satisfaction route (substitution + FullSimplify), guarded by non-vacuity and positive-branch
checks that would catch a wrong closed form; it does not secretly re-run Solve. F0 is notes:506-only
(×100%) with card and script multiplier intact; F3 banner now reads STAGE 248. Both engines exit
0 with all checks passing and outputs freshly regenerated. material_change is false — the
re-derivation reproduces the same compiler forms. Verdict: verified.
