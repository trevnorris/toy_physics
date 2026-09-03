# P2a slab-row join — decision-leg review record (2 legs, computation-backed, convergent → v2)

Directive `directives/S11c_b_p2a_slab_row_join_directive.md`. Legs Codex + Grok (orchestrator-written, rule 7). Logs
`~/.s11_build/S11c_b_p2a/decision_{codex,grok}.log`. **v1 REJECTED by both (10/10 defects, file:line + literal
probe).** v2 folds both once. What HELD (both): A1 `SOURCE_OPERAND≡EXPANDED`, A7 diagnosis (both engines carry face
⇒ `row_residual:427` WL-only subtraction is the #90 break), A8 family (`MU_THETA_OPERATOR`), U/E partners, A10.

## Convergent defects folded into v2
1. **DOF must be pinned** (both): `row_residual` drops DOF when indexing (`:566`); an unpinned mass/center-face DOF
   lets a one-sided case evade the duplicate guard + collide. ⇒ A1 `DOF=MASS`, A4 `DOF=CENTER_FACE` exactly;
   `SLAB_OBJECTS = {U_MOMENTUM_ROWS, THICKNESS_ROW, MASS_EVOLUTION_ROW, CENTER_FACE_GENERALIZED_ROW}`.
2. **SIGN conventions — TWO, both surfaced not normalized** (both): PY `−kinetic` vs WL `+kinetic` in complete rows
   (`sympy_audit.py:2325/2946` vs `…audit.wl:1345`); PY `+diff` vs WL `−linearVirtualVariation` face sign
   (`:2171`/`:1139`). ⛔ Comparator must not flip a sign to force zero NOR ship `2K` silently ⇒ A5: emit the residual
   with the convention intact; adjudicate representational-vs-real in the STEP RECORD (rule 1/6).
3. **A6 must CONTAIN the closed disposition table** (both enumerated the full PY/WL inventory), with an
   observed-path-set guard. `FACE_FLUX_BOUNDARY_OPERANDS`(19-key)/`FACE_SHAPE_SUBSTRATE`(8-field) = engine-only raw
   T-substrate (joining it = the §3c-forbidden parallel route; the current whole-container `DOF=ALL` join is
   removed). ⛔ 19 keys, not `20` (Codex/Grok both counted).
4. **Origin-family name mismatches** (both): raw same-name join is false — `PY ADVECTIVE={BACKGROUND_ADVECTION}` vs
   `WL ADVECTIVE=div(σ u_t)`; `PY MASS_EVOLUTION` vs `WL ACCUMULATION`. ⇒ A6b adapters
   `PY(MASS_EVOLUTION+ADVECTIVE)↔WL(ACCUMULATION+ADVECTIVE)`, `PY FACE_FLUX↔WL FLUX`, `PY FACE_VIRTUAL_WORK.ROWS↔WL
   FACE.rows`; no un-adapted same-name residual.
5. **μ_θ** (both): adding `MU_THETA_OPERATOR` to `SLAB_OBJECTS` no-ops (`axes.OBJECT` is None → continue) ⇒ A8 a
   SEPARATE `_family_cases("MU_THETA_OPERATOR")` residual-emitting exact-scalar ε¹ comparison, ⛔ not
   `FACE_ATTRIBUTION_ONLY` (emits no residual); not double-joined with `DIVERGENCE_FORM_SOURCE.MU_THETA`.
6. **A7 provenance must be SYMMETRIC** (Codex): current face provenance is WL-only; a builder could delete
   discovery+emission and pass ⇒ role `PROVENANCE_ONLY_ALREADY_IN_COMPLETE_ROW`, both engines, no arithmetic use.
7. **CENTER_FACE currently WL-only attribution** (Grok): `_emit_center_face_rows` emits WL-only, no PY operand/
   residual ⇒ A8 CENTER_FACE flows through `_strong_rows`, `_emit_center_face_rows` is NOT the residual.
8. **A9 freshness + argparse** (both): `row_residual` has NO argparse — always loads stale `C.DEFAULT_PY/WL`
   (`:1086-1096`; committed defaults verified stale — PY lacks `FACE_GENERALIZED_FORCE_ROWS`, WL lacks
   `THETA_VARIATION_RULE`). ⇒ add `--py`/`--wl` SELECTABLE + USED + a CONTENT SHA-256 (a path hash ≠ freshness).
   Guards over EXACT FULL keys, per-engine set+alignment; ⛔ NOT energy families (row_residual doesn't load them).
9. **P2a↔P2b sequencing** (Codex): P2a's row/μθ parse passes through `_extra_basic`'s applied→bare collapse
   (`widthBackground(chi)−widthBackground(x)→0`) which B4 (P2b) fixes ⇒ **P2b lands before P2a's final validation.**

## Status
v2 folded once (rule 7). ⇒ Codex build (after P2b's B4) → 2 build legs. The two sign conventions (kinetic, face) are
the one genuine physics adjudication — resolved by SURFACING in the residual + step-record interpretation, not
comparator normalization.
