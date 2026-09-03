# P2a — slab-row cross-engine join + `row_residual` structural fixes — decision list v2

⚠ **v2 — folded BOTH decision legs (Codex + Grok, computation-backed, convergent), rule 7 one pass.** Record
`directives/_measurements/S11c_b_p2a_slab_row_join_directive.md`. Targets: `scripts/S11c_b_cross_engine_comparator.py`
(`extract_slab`, `extract_term_origins`) AND `scripts/S11c_b_row_residual.py`. Physics-bearing; ⛔ never
blanket-collapse. Inputs = the fresh 4-case primaries `.out` from both engines (P1-WL gate + PY `PRIMARIES_ONLY`).
⚠ **Sequencing: P2b (the §3a bridge, incl. the argument-sensitive applied→bare guard B4) lands BEFORE P2a's final
`row_residual` validation** — P2a's row/μθ parse passes through the `_extra_basic` applied→bare collapse that B4
fixes (`widthBackground(chi) − widthBackground(x) → 0`); until B4 lands, P2a validation can manufacture agreement.

## What HOLDS (both legs agree — do not re-litigate)
PY `THETA_BALANCE.SOURCE_OPERAND ≡ EXPANDED` (same `mass_balance`); A7's diagnosis (both engines carry face in the
complete U/E rows, so `row_residual:427`'s WL-only face subtraction is the #90-stale break); μ_θ is
`MU_THETA_OPERATOR` (not the slab θ-row, not `MU_THETA_FACE_BINDING`); U/E partners
`U_BODY_BALANCE.EXPANDED↔U_MOMENTUM_ROWS`, `E_W_BALANCE.EXPANDED↔THICKNESS_ROW`; `extract_coupling` unchanged.

## extract_slab — what must be TRUE (key = OBJECT **and** DOF; `row_residual` drops DOF when indexing, so an
## unpinned DOF lets a one-sided case evade the duplicate guard and collide — `comparator.py:250`, `row_residual:566`)
A1. Join ONLY `THETA_BALANCE.EXPANDED` → WL `MASS_EVOLUTION_ROW` at **exactly `OBJECT=MASS_EVOLUTION_ROW, DOF=MASS`**.
    `SOURCE_OPERAND` = duplicate-source diagnostic (not a second join); `EVOLUTION_TERM_ORIGINS` → origin family
    only. ⛔ Remove the stale `THETA_BALANCE → ("MU_THETA","THETA")` map (`comparator.py:768`).
A2. Remove/re-key the `ADVECTIVE_MASS_OPERAND → MASS_EVOLUTION_ROW` emit (`comparator.py:777`) — else A1 duplicates
    the axis → `_family_cases` RAISES (`adjudicated_comparison.py:1091`). ADVECTIVE stays an origin/diagnostic.
A3. Re-key or EXCLUDE the WL canonical `MU_THETA` promotion of `DIVERGENCE_FORM_SOURCE.MU_THETA`
    (`comparator.py:797-800`) — it is the un-activated energy variation, NOT the residual's μ_θ. ⛔ After A1 removes
    PY's slab `MU_THETA`, a one-sided WL slab `MU_THETA` makes `row_residual` RAISE. No canonical slab `MU_THETA` on
    EITHER engine.
A4. JOIN WL `CENTER_FACE_GENERALIZED_ROW` ↔ PY `FACE_GENERALIZED_FORCE_ROWS.CENTER_FACE_GENERALIZED_ROW` at exactly
    `OBJECT=CENTER_FACE_GENERALIZED_ROW, DOF=CENTER_FACE` (⛔ NOT optionally excluded — it has a real partner). ⛔ Do
    NOT re-join `FACE_GENERALIZED_FORCE_ROWS.{U,E_W}` (already inside the A5 rows = double-count) or `THETA_FACE_FLUX`
    (already inside the A1 mass row). ⛔ Do NOT sign-normalize (see A5).
A5. **SIGN CONVENTIONS ARE SURFACED, NOT NORMALIZED (rule 1/6, scripts-print-never-assert).** The two engines carry
    OPPOSITE whole-row sign conventions in the complete U/E rows — PY `operator_from_density` has `−kinetic`
    (`sympy_audit.py:2325/2363`, subtracted again post-fold `:2946`), WL has `+kinetic` (`…audit.wl:1345`); likewise
    the face generalized-force sign (PY `+diff(work,δv)` `sympy_audit.py:2171` vs WL `−linearVirtualVariation`
    `…audit.wl:1139`). ⛔ The comparator MUST NOT flip a sign to force zero (that manufactures agreement) NOR ship a
    `2K`/`2·face` residual as physics silently. The joined operands are the complete rows AS EMITTED; the residual is
    emitted with the sign convention intact, and whether a whole-row global negation is representational or a real
    discrepancy is ADJUDICATED IN THE STEP RECORD, not in the comparator.
A6. **The directive CONTAINS the closed disposition table** (⛔ not "classify these"); every path is JOINED(partner,
    key) or ENGINE-ONLY(reason); an exact observed-path-set guard raises on any unlisted path (the extractor
    silently ignores unknown paths today):
    - PY slab: `U_BODY_BALANCE.EXPANDED`→U_MOMENTUM_ROWS(U); `E_W_BALANCE.EXPANDED`→THICKNESS_ROW(E_W);
      `THETA_BALANCE.EXPANDED`→MASS_EVOLUTION_ROW(MASS); `FACE_GENERALIZED_FORCE_ROWS.CENTER_FACE_GENERALIZED_ROW`→
      CENTER_FACE(CENTER_FACE); `{U/E_W_BALANCE}.{LOCAL,DIVERGENCE_FLUX}`, `THETA_BALANCE.{SOURCE_OPERAND,
      EVOLUTION_TERM_ORIGINS}`, `FACE_GENERALIZED_FORCE_ROWS.{U,E_W,THETA_FACE_FLUX,SOURCE_OPERANDS}`,
      `ADVECTIVE_MASS_OPERAND`, `MU_THETA_FACE_BINDING`, `FACE_FLUX_BOUNDARY_OPERANDS` (the 19-key raw T-substrate
      bundle, ⛔ engine-only — joining it is the §3c-forbidden parallel route) = ENGINE-ONLY/diagnostic.
    - WL slab: `ROWS.{U_MOMENTUM_ROWS,THICKNESS_ROW,MASS_EVOLUTION_ROW,CENTER_FACE_GENERALIZED_ROW}` = JOINED (above);
      `DIVERGENCE_FORM_SOURCE.*`, `BACKGROUND_FIRST_JETS`, `VIRTUAL_CONSTRAINT`, `THETA_VARIATION_RULE`,
      `FACE_SHAPE_SUBSTRATE` (8 fields/face, raw substrate) = ENGINE-ONLY. ⛔ The `FACE_SHAPE_SUBSTRATE`/
      `FACE_FLUX_BOUNDARY_OPERANDS` whole-container join (`comparator.py:779-780/804-805`, 19-vs-8 unlike payloads)
      is REMOVED. ⛔ Use the exact 19-name substrate schema, not a `20` count.
A6b. **Origin-family name mismatches get ADAPTERS or are engine-only — NOT a raw same-name join.** `extract_term_origins`
    uses raw names (`comparator.py:809`), but PY `ADVECTIVE`={BACKGROUND_ADVECTION} vs WL `ADVECTIVE`=div(σ u_t)
    (background advection + velocity dilatation); PY `MASS_EVOLUTION`={DENSITY_TIME+VELOCITY_DILATATION} vs WL
    `ACCUMULATION`=density-time. The joined buckets are `PY(MASS_EVOLUTION+ADVECTIVE) ↔ WL(ACCUMULATION+ADVECTIVE)`,
    `PY FACE_FLUX ↔ WL FLUX`, `PY FACE_VIRTUAL_WORK.ROWS.{U,E_W} ↔ WL FACE.{U_MOMENTUM_ROWS,THICKNESS_ROW}`; `KINETIC`
    matches. Any bucket without an adjudicated adapter is diagnostic-only, never a raw same-name residual.

## row_residual — what must be TRUE
A7. Remove the one-sided face subtraction (`row_residual:427`): face stays in BOTH complete rows (compare like
    complete rows). Face-origin provenance is EMITTED SYMMETRICALLY (both PY and WL leaves, role
    `PROVENANCE_ONLY_ALREADY_IN_COMPLETE_ROW`) and ⛔ never used arithmetically; ⛔ a builder may not delete face
    discovery+emission together (the self-ledger would then never miss them).
A8. μ_θ is compared via a SEPARATE `_family_cases("MU_THETA_OPERATOR")` call (⛔ NOT by adding it to `SLAB_OBJECTS` —
    that no-ops: `axes.get("OBJECT")` is None → continue): a residual-EMITTING exact-scalar ε¹ comparison, 4 aligned
    `BRANCH×DENSITY` cases, one PY/WL/residual per case, ⛔ not `FACE_ATTRIBUTION_ONLY` (which emits no residual). Not
    double-joined with `DIVERGENCE_FORM_SOURCE.MU_THETA`. Post-remap `SLAB_OBJECTS = {U_MOMENTUM_ROWS, THICKNESS_ROW,
    MASS_EVOLUTION_ROW, CENTER_FACE_GENERALIZED_ROW}` (μ_θ NOT in it). `_emit_center_face_rows` is NOT the CENTER_FACE
    residual — CENTER_FACE flows through `_strong_rows` like the other slab objects.
A9. Expected-set GUARDS over EXACT FULL KEYS (before any axis is dropped), checking EACH engine's set + alignment
    (the loader only unions, so a case missing from BOTH vanishes): row-bearing objects 4 each (U/E_W/MASS/
    CENTER_FACE), μ_θ 4, origins `KINETIC` 4, coupling (existing `required_coupling_counts`), admissibility. ⛔ Do
    NOT guard energy families here — `row_residual` does not load them (that guard is the comparator's/P2b's).
    `row_residual` gains `--py`/`--wl` (it has NO argparse today, `:1086-1096`) that are SELECTABLE AND USED; the
    emit includes each path AND a SHA-256 CONTENT hash of the parsed bytes (⛔ a path hash does not prove freshness —
    the committed defaults are stale). A11's validation command passes the freshly generated files explicitly.
A10. Output states operands + residual per object; no verdict token, no assert on the residual (rule 5).

## Provenance
A11. S11c-a not regressed. Validate by RUNNING `row_residual --py <fresh PY> --wl <fresh WL>` on the fresh 4-case
     `.out`, AFTER P2b (B4) lands. The residual (incl. any whole-row sign convention) is OUR finding for the step
     record.

## Legs
Codex build ⇒ 2 build legs (fresh Claude + Grok). Mandatory negative ablations that must BITE: wrong `DOF` on the
mass/center-face join → duplicate-guard/projection-collision surfaces; re-introduce the one-sided face subtraction →
residual gains the face term; drop the μ_θ `_family_cases` call → μ_θ residual vanishes; break a disposition entry →
the observed-path-set guard raises; a raw same-name `ADVECTIVE` join → manufactured disagreement; remove the same
case from both `.out` → the expected-set guard raises (not silent).
