---
unit_id: 239
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-03T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 239 (checkpoint)

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
The entire derivation body of the `.wl` was re-authored (diff replaces wl:56-373).
Independent-route choices, contrasted with the `.py` it formerly mirrored:

- Compiler matrix: the `.py` hardcodes `S_rm_dep = Matrix([[0,0],[0,-1],[1,-1]])`
  then `C_phys_dep = S_rm_dep * M_phys` (py:174-175). The new `.wl` instead
  *defines* the boxed dependent vector `dependentDelta = {0, -V, U-V}` (wl:135-139)
  and *recovers* the compiler by Jacobian: `compilerJacobian = Table[D[...]]`
  (wl:140-144). The `SrmDep . Mphys` two-step is gone (grep confirms no `SrmDep`,
  `CphysDep`, or "inherited from Stage" string remains).
- Left inverse: the `.py` uses `solve(...)` compared `!=` a literal `expected_solution`
  dict, plus a hardcoded `L_phys_dep = Matrix([[0,-1,1],[0,-1,0]])` (py:197-212).
  The new `.wl` uses `leftNative = PseudoInverse[compilerJacobian]` (wl:150) and
  verifies `leftNative . compilerJacobian == I_2`, `leftNative . dependentDelta == {U,V}`,
  and the coordinate formulas symbolically (wl:152-154). No literal-list compare,
  no `expectedInverse` (grep confirms gone).
- Orbit lock: the `.py` does `solve([Eq(y_dep[1],0),Eq(y_dep[2],0)],[U,V])` compared
  `!= [{U:0,V:0}]` (py:310-316). The new `.wl` uses native
  `Reduce[dependentDelta[[2]]==0 && dependentDelta[[3]]==0, {U,V}, Reals]` then
  `Equivalent[orbitLockReduced, U==0 && V==0]` (wl:229-235).
- Projectors built via `DiagonalMatrix[{1,0}]` / `DiagonalMatrix[{0,1}]` (wl:117-118)
  rather than the `.py`'s explicit `{{1,0},{0,0}}` literals.
- Labels are an entirely new M1–M9 manifest scheme, not the `.py`'s verbatim English.
- First-order form differentiates the full perturbed dependent vector at once
  (wl:242-251) rather than per-coordinate `Ufirst`/`Vfirst` then recompose.

**Assessment:**
Genuinely independent route. The compiler is now obtained *forward* (boxed vector →
Jacobian) where the `.py` went *backward* (decomposition matrix → product); these are
distinct derivations that would not share an algebra slip. `PseudoInverse` and `Reduce`
are native primitives that compute the left inverse and lock set from scratch rather
than comparing to a pre-written answer, so they would catch a wrong target value. The
shared `Derivative[...]→0` support-blindness primitive is the one the directive
explicitly permits as unavoidable; its surrounding choreography (`Outer[D, vec,
supportVars]`, `supportResiduals` helper) is the verifier's own, not the `.py`'s
`applyfunc(lambda x: diff(...).subs(...))`. All nine manifest claims M1–M9 are present
and PASS. No collateral edits beyond the re-author + relabel.

### F2 — stale stage label (script-side, cosmetic)

**Classification:** resolved

**What changed:**
Banner is now `"STAGE 239 — RIGID-MOUTH PHYSICAL NORMAL FORM"` (wl:65). All `Stage221`
suffixes renamed to `Stage238` (`T2Stage238`, `RtargetStage238`, `qNtRigidStage238`,
`UStage238`, `VStage238`, `xPhysStage238`, `yDepStage238`; wl:158-176). The now-unused
`Bcoeff`/`Ccoeff` symbols were dropped from `Clear`/`$Assumptions` (consistent with the
re-author, harmless). `grep -i "STAGE 222\|Stage221"` returns nothing (exit 1).

**Assessment:**
Correct and complete. Single-file fix; no batch renumber. Output header reads
"STAGE 239" (mathematica output line 8). Per directive note, F2 is satisfied by the
F1 re-author. The `.py` "Stage 236/221" comments were correctly left untouched
(explicitly out of scope).

## Exec log assessment

**SymPy:** exit=0. The `.py` was not modified (diff touches only the `.wl`); it is the
independent first engine. 52 `[ok]` lines + "All Stage 239 symbolic checks passed."
Notable: "left inverse of physical compiler", "U = V = 0 cancels the dependent defect",
"first-order dependent correction" all `[ok]`.

**Mathematica:** exit=0. 42 PASS, 0 FAIL. Notable:
`M5 native left inverse matrix → MatrixForm[{{0,0},{0,0}}] PASS`;
`M9 Reduce orbit-lock set = U == 0 && V == 0` followed by `M9 lock set is chart origin
= True PASS`; `M4 Jacobian rebuilds dependent vector → {{0,0,0}} PASS`. Check count
reconciles: 41 `expectZero` call sites + 1 `expectTrue` call site = 42 = PASS lines.
FAIL path is non-vacuous: `expectZero`→`allZeroQ` (`TrueQ[#==0]`) and `expectTrue`
(`TrueQ[res]`) both route a wrong residual to `fail[...]` which prints the residual and
`Exit[1]`; a wrong target would flip exit to 1.

**Output freshness:** confirmed. `.wl` mtime 2026-06-03 08:36:10; `.wl` output mtime
08:42:08 (newer); `.py` output mtime 08:42:08 — both outputs regenerated post-fix.

## Material-change assessment

`material_change`: false. The edit re-authored the second-engine derivation route and
corrected cosmetic stage labels. No derived constant, matrix, or identity changed: the
compiler `(0,-V,U-V)`, `M_phys=I_2`, the target ratio, the left-inverse formulas, and
the orbit-lock equivalence are identical to before and to the `.py`. Nothing downstream
depends on the `.wl`'s internal variable names or route. No downstream staleness.

## Side observations (non-blocking)

- The SymPy emits 52 `[ok]` lines vs the report's approximate "~50/51"; this is the
  unmodified first engine and the count is fine (the report estimate was loose).
- `PseudoInverse[compilerJacobian]` yields the minimum-norm left inverse, which may
  differ from the `.py`'s hand-picked `L_phys_dep = {{0,-1,1},{0,-1,0}}`. That is
  expected and acceptable: the checkpoint claim is that *a* left inverse exists with
  `L.C=I_2` recovering `U,V`, which `PseudoInverse` satisfies — and the difference is
  itself further evidence of route independence, not a defect.

## Verdict justification

Both findings are resolved. F1's re-authored `.wl` is genuinely independent — it builds
the compiler forward by Jacobian (vs the `.py`'s backward `S_rm_dep . M_phys`), obtains
the left inverse via native `PseudoInverse` (vs literal-list compare), and resolves the
orbit lock via `Reduce`+`Equivalent` (vs `Solve`-then-`=!=`), with its own M1–M9 labels;
the only shared primitive (`Derivative→0`) is the one the directive permits. F2's banner
and `Stage221` suffixes are corrected to 239/238 with no stale 22x labels remaining and
no batch renumber. Both engines exit 0; the 42 Mathematica checks reconcile with call
sites; the FAIL/`Exit[1]` path is real; outputs are freshly regenerated. No regressions,
no material change.
