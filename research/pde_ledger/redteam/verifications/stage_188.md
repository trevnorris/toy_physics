---
unit_id: 188
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-01T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 188

## Per-finding outcomes

### F1 — insufficient_verification (Section VII vacuous M*0 checks)

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.py:172-196`. The
diff (`exec_logs/stage_188_diff.patch` lines 30-62) deletes the two `Matrix([0,0,0])`
substitution checks and replaces Section VII with: (a) two bijection round-trips on the
GENERIC `obs` packet — `C_def_to_obs * (C_obs_to_def * obs) - obs` (line 179) and
`C_quot_to_obs * (C_obs_to_quot * obs) - obs` (line 183) — and (b) two
`det*(1/det) - 1` nonzero-determinant checks (lines 191, 195).

**Assessment:**
Genuinely de-tautologized. Both round-trips reference `obs = Matrix([dr, dn, de])`
(line 71), NOT the zero vector. They reduce to `obs - obs = 0` only because
`C_def_to_obs` and `C_quot_to_obs` are true inverses of the forward compilers — a singular
compiler would fail them — so they exercise injectivity/shared-zero-set, which is the
nontrivial content of the iff. The output (lines 256-269) shows all four printing the zero
vector / `0`. No `M*0` remains as a load-bearing check. Not tautological.

### F2 — tautological_check (Cstar := 1/Ctr then assert Cstar*Ctr-1)

**Classification:** resolved

**What changed:**
Same file, lines 46-49 and line 62. The diff (patch lines 9-13, 21-22) replaces
`Cstar = sp.simplify(1/Ctr)` with the independent appendix closed form
`(1+chi0s)(1+deltaUs)(1+chi0s+deltaUs)/(chi0s*deltaUs)`, and replaces the assertion
`Cstar*Ctr - 1` with `Cstar - 1/Ctr`.

**Assessment:**
Genuinely de-tautologized. `Cstar` is now built from the four positive symbols via its own
closed form, and `Ctr` is independently built (lines 42-45) as
`chi0s*deltaUs/((1+chi0s)(1+deltaUs)(1+chi0s+deltaUs))`. The check `Cstar - 1/Ctr`
therefore compares two independently constructed expressions and would fail on any mistype
in either; it is no longer value-into-its-own-defining-equation. Output line 25 shows
`C_* - 1/C_tr,* = 0`. The adjacent load-bearing `Astar - Bstar*Ctr` check is untouched
(line 63). Value of `Cstar` is unchanged, so the downstream `diag(-Cstar,1,1)` compiler
and `det` are unaffected — no collateral edit, confirmed by the tightly scoped diff.

### F3 — missing_mathematica (no .wl; dual-engine rule)

**Classification:** resolved

**What changed:**
New `mathematica/moving_throat_pde_stage188_branch_observables_completion_mathematica_audit.wl`
(277 lines) created. It re-verifies every load-bearing assertion of the CORRECTED SymPy
script: the two coefficient identities (`C_* - 1/C_tr,*` and `A_tr,* - B_* C_tr,*`), both
compiler matrices + inverses, the factorization, the determinants, the complementary-drift
identity, and the corrected Section VII generic-`obs` bijection round-trips.

**Assessment — independent route, NOT a transliteration:**
The `.wl` uses a genuinely different derivation machine than the `.py`:
- The SymPy script hand-TYPES each compiler matrix as a literal
  (`sp.diag(-Cstar,1,1)`, `C_quot_to_def = Matrix([...])`,
  `C_obs_to_def_expected = Matrix([...])`) and then multiplies/inverts.
- The `.wl` does NOT type the compilers to verify them. It writes the underlying *drift
  relations* as equations (`sigTr == firstLogDrift[Exp[-cStar*small*rObs]]`,
  `nObs == nCompositeDrift`, `rDef + xi == etaComplementDrift`, etc.), `Solve`s them for
  the packet variables over `Reals`, and then RECOVERS each compiler matrix as the
  Jacobian `Outer[D, outputs, inputs]` (`linearCompiler`, lines 45-48). The forward,
  quotient, and inverse compilers are each obtained from an independent
  `Solve` + differentiate pipeline (lines 103-129, 133-167, 178-241).
- It uses Mathematica-native primitives absent from the SymPy file: `Series`+`Coefficient`
  on `Log` (`firstLogDrift`, lines 40-43), `Solve[..., Reals]`, `Outer[D,...]`, `Inverse`,
  `Det`, `LinearSolve`-class machinery. The complementary drift is independently derived
  by series-expanding `Log[(1 - eta Exp[small etaObs])/(1 - eta)]` (lines 173-176) rather
  than the SymPy script's hand-typed `-epsetas/(1-epsetas)*de`.
- The SymPy-typed matrices appear ONLY as the right-hand cross-engine comparison targets
  (`sympyObsToQuot`, `sympyQuotToDef`, `sympyObsToDef`, `sympyDefToObs`), each subtracted
  from the *independently derived* Mathematica matrix (lines 126, 153, 201-202, 234). That
  is the intended cross-engine agreement check, not a transliteration of the derivation.

**Substantive checks (no vacuous Mathematica check):**
Every `expectZero`/`expectZeroPacket` operand depends on real derived quantities. The
Section VII bijection round-trips (lines 253-260) use `obsToDefDirect.obsVars` /
`obsToQuot.obsVars` on the generic symbolic `obsVars = {rObs,nObs,etaObs}` — NOT a zero
vector — mirroring the corrected SymPy semantics; they hold only because the solved
inverses are true inverses. No `M*0`, no X−X with X typed twice, no hardcoded literal
checked against itself. The `det*(1/det)-1` lines (264-271) use the
independently-computed `Det[obsToDefDirect]`/`Det[obsToQuot]`.

## Exec log assessment

**SymPy:** exit=0 (`exec_logs/stage_188_sympy.log` final line `# exit_code: 0`).
Notable lines: `C_* - 1/C_tr,* = 0`, `A_tr,* - B_* C_tr,* = 0`, Section VII
`C_obs->def then inverse recovers generic obs (bijection, def side) = [0,0,0]` and the two
`1/det(...) well-defined (nonzero det) = 0` lines — all the corrected checks pass on the
generic packet.

**Mathematica:** exit=0 (`exec_logs/stage_188_mathematica.log` final line `# exit_code: 0`;
ends `Stage 188 Mathematica audit passed.`). Every line is a `PASS:`, including
`PASS: C_* - 1/C_tr,*`, `PASS: observable-to-quotient compiler agrees with SymPy`,
`PASS: factorized compiler - direct compiler`,
`PASS: direct compiler agrees with SymPy expected compiler`,
`PASS: C_obs->def then inverse recovers generic obs (bijection, def side)`.

**Output freshness:** confirmed. `sympy_audit.txt` (mtime 11:23:02) and
`mathematica_audit.txt` (mtime 11:23:02) are both newer than their scripts
(`.py` 11:16:43, `.wl` 11:20:18) — regenerated post-fix.

## Cross-engine agreement

Yes. Both engines produce the identical conclusions for all shared quantities:
`C_tr,*`, `C_*`, `B_*`, `A_tr,*` closed forms match; `det(C_obs->quot) = -C_*`,
`det(C_obs->def) = eta/(eta-1) = -eps/(1-eps)` match; all compiler matrices and both
inverses match; and the `.wl` explicitly subtracts the SymPy-typed compilers from its own
independently-derived ones and gets the zero matrix in all four cases
(`observable-to-quotient`, `quotient-to-defect`, `direct`, `inverse` agreement PASSes).
The corrected Section VII bijection checks agree on both engines.

## Material-change assessment

`material_change`: false. F1 and F2 change only HOW two already-true facts are checked
(`Cstar`'s value is unchanged — same closed form, just written non-circularly; Section VII
swaps vacuous scaffolding for a genuine bijection check). No derived value that downstream
units consume is altered. F3 adds a new verification artifact and changes nothing derived.
No downstream unit needs re-audit on account of this stage.

## Side observations (non-blocking)

- The script's banners/ledger still print "STAGE 171" (e.g. line 35, 198) while the file is
  the stage-188 audit. Pre-existing cosmetic mislabel, outside the scope of these findings;
  flagging only, not blocking.

## Verdict justification

All three findings are `resolved`. F1's Section VII now exercises the bijection on the
generic `obs` packet (load-bearing, fails if either compiler is singular) instead of the
vacuous `M*0`; F2 defines `Cstar` from its independent appendix closed form and cross-checks
it against `1/Ctr` (no longer value-into-its-own-equation); and F3 adds a genuinely
independent Mathematica route that derives every compiler via `Series`/`Solve`/`Outer[D]`
(Jacobian) and cross-checks against the SymPy literals — not a line-by-line port. No new
tautology was introduced on either engine, both exec logs exit 0, outputs are fresh, the
diff is scoped exactly to the named lines plus the new `.wl`, and the two engines agree on
every conclusion. Verdict: verified.
