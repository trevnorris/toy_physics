---
unit_id: 192
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-01T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 192

## Per-finding outcomes

### F1 — missing_mathematica (dual-engine gap)

**Classification:** resolved

**What changed:**
Codex created the new second-engine script
`mathematica/moving_throat_pde_stage192_orbit_quotient_projectors_mathematica_audit.wl`
(8133 bytes, 257 lines). It re-verifies every load-bearing assertion of the SymPy
reference `scripts/moving_throat_pde_stage192_orbit_quotient_projectors_sympy_audit.py`
and prints all results, then `Exit[0]`. No SymPy/paper/notes files were touched
(directive scope honored; diff patch empty, consistent with a single new untracked file).

**Assessment — genuine independent route (NOT a transliteration):**

The two scripts reach the same conclusions by materially different Mathematica-native
machinery. Concrete contrasts:

- **Quotient section S.** SymPy builds it by *embedding × inverse*:
  `Edep` (8×3 unit-column embedding) times `Pdep.inv()`, i.e. `Sdep = Edep * Pdep_inv`
  (.py lines 89–94). The `.wl` instead builds it by a **constrained `LinearSolve`** on
  an augmented 8×8 system: `constrainedSystem = Join[quotientRows, freeSelector]`, then
  for each column `LinearSolve[constrainedSystem, Join[UnitVector[3,col], 0,0,0,0,0]]`
  (.wl lines 76–113). This is the directive's expected "S is the unique vector with zero
  free coordinates whose quotient packet is e_col" construction — a different primitive
  (solve an augmented linear system) than SymPy's explicit 3×3 inverse times an embedding.
  The output (txt line 22) shows `sectionFromSolve - sectionExpected == 0`, so the two
  routes land on the identical S.

- **Quotient-failure and orbit pieces.** SymPy obtains them only by projector
  multiplication (`Qquot*Dx`, `Oorb*Dx`). The `.wl` *cross-validates* each against an
  independent constrained solve: `failureByConstrainedSolve = LinearSolve[constrainedSystem,
  Join[packet, 0…]]` (.wl 149–154) and `orbitByKernelSolve = LinearSolve[constrainedSystem,
  Join[{0,0,0}, drift[[freeSlots]]]]` (.wl 176–181), asserting `Q Δx − constrained solve == 0`
  (txt line 53) and `O Δx − constrained kernel solve == 0` (txt line 64). The corroborating
  quantity is the solve, not the typed literal.

- **Orbit-lock equivalence `Q Δx = 0 ⟺ M_* Δx = 0`.** This is the strongest independence
  evidence. SymPy never solves it directly — it infers `⟺` from the `Q = S M_*` construction
  plus an independent q-packet (`M_*(Sq)=q`, `Q(Sq)=Sq`, `O(Sq)=0`). The `.wl` proves it by
  two *separate* `Solve` calls on different left-hand systems: `lockByMap = Solve[M_* Δx == 0,
  {dt,dk,dm}]` (3-row system) and `lockByProjector = Solve[Q Δx == 0, {dt,dk,dm}]` (8-row
  system), then asserts the two dependent-triple solution laws are identical (.wl 226–247;
  txt line 92 `Solve laws ... agree == 0`) and that each law kills the other operator
  (txt lines 94, 96). A genuinely distinct logical route to the same headline equivalence.
  The `.wl` also keeps the SymPy-style corroboration (Section V, txt 78–83) as belt-and-suspenders.

Variable names also differ (drift `dl,dc,dg,du,dk,dw,dm,dt` and parameters `chi,delU,ee,ff`
vs. SymPy's `Delta_*` / `chi0_star,…`), and the section order differs (Mathematica derives S
column-by-column via solve, SymPy in one `Edep*Pdep_inv` product). This is not a re-typed port.

**Assessment — substantive, no tautology:**

Every load-bearing claim is exercised against an *independently computed* quantity, not
substituted back into its own defining equation:
- Pivot block / inverse: the typed `pivotExpected`/`pivotInverseExpected` are compared to
  `quotientRows[[All,depSlots]]` and `Inverse[pivotBlock]` (bottom-up), plus the
  self-validating `P^{-1}P − I`, `P P^{-1} − I` (txt 17,19).
- `M S − I_3 == 0` (txt 26) confirms S is a true right inverse independent of any literal.
- Six projector identities `Q²−Q, O²−O, QO, OQ, MO, MQ−M` all zero (txt 32–43) — genuine,
  would fail if the solve-built S or projectors were wrong.
- Failure support `Q free-coordinate rows == 0` (txt 45) and the dependent-triple closed forms
  cross-checked two ways (projection vs. constrained solve, txt 53/55).
- The one near-vacuous line `S * 0 == 0` (.wl 221) is the same harmless redundant scaffolding
  the original auditor flagged as SymPy A16 — trivially true, load-bearing on nothing, not a
  finding. It does NOT carry the verification: every substantive claim above is independently
  pinned. No `X − X`, no `M*0==0` masquerading as a real check, no literal-compared-to-itself.

**Cross-engine agreement:** confirmed. Same `M_*` entries, same pivot block, same
`det = 1+χ` (txt 13), same inverse, same section S (txt 22), identical projector-identity
residuals, and the failure/orbit closed forms — though printed in a different algebraic
arrangement — reduce to the SymPy-form `failureExpected`/`orbitExpected` literals with zero
residual (txt 54, 65). The mu-component `ff·q_tr/(1+chi) + q_nt − q_η` matches SymPy's
`F_*/(1+chi0)·q_tr + q_nt − q_η`. All conclusions agree.

## Exec log assessment

**SymPy:** exit=n/a. `redteam/exec_logs/stage_192_sympy.log` is absent — expected, since
Codex made no SymPy edit (the directive forbade it; SymPy is the untouched reference engine).

**Mathematica:** exit=0. `redteam/exec_logs/stage_192_mathematica.log` ends with
`# exit_code: 0` and every line is a PASS. Notable:
- `PASS: M S - I_3` (S from constrained LinearSolve is a genuine right inverse)
- `PASS: Solve laws for Q Delta x == 0 and M Delta x == 0 agree` (independent `⟺` proof)
- `PASS: O Delta x - constrained kernel solve` / `PASS: Q Delta x - constrained failure solve`
  (projector pieces corroborated by independent solves)

**Output freshness:** confirmed. `.wl` mtime 2026-06-01 11:09:55; committed
`output/...stage192...audit.txt` mtime 2026-06-01 11:11:38 (newer than the script).
Output regenerated post-creation. Diff patch empty (single new untracked file, no edits to
tracked files).

## Material-change assessment

`material_change`: false. This finding only ADDS an independent second-engine confirmation of
already-verified SymPy conclusions; it derives no new result and changes no carried quantity.
All numbers (M_*, P, S, Q, O, failure/orbit forms) agree with the existing SymPy reference, so
no downstream unit's inputs change. No `upstream_stale` propagation warranted on numerical grounds.

## Side observations (non-blocking)

- The `.wl` retains the cosmetic pre-renumber lineage references inside comments/ledger text
  (e.g. "Stage 187 single-orbit law" concept), consistent with the SymPy script's "STAGE 175"
  banner that the original auditor already noted as cosmetic-only. The `.wl`'s own banners
  correctly say "STAGE 192", and these strings feed no assertion. Not a finding.
- `expectTrue` is defined (.wl 52–56) but unused; harmless dead helper.

## Verdict justification

`verified`. The single finding (missing second engine) is resolved by a genuinely independent
Mathematica script: it reaches the same conclusions through constrained `LinearSolve` section
construction, solve-based corroboration of the failure/orbit pieces, and a standalone two-`Solve`
proof of the orbit-lock equivalence — none of which mirror the SymPy embedding×inverse / q-packet
route. Every substantive assertion is checked against an independently computed quantity (the one
vacuous `S*0==0` line is non-load-bearing scaffolding only). The committed output exits 0 with all
PASS, is fresh relative to the script, and agrees entry-for-entry with the SymPy reference.
