# S11c-b adjudication layer — BUILD review, two legs, consolidated (SOUND, NO FALSE MATCH)

**Artifact:** `scripts/S11c_b_adjudicated_comparison.py` (Codex-built; imports the committed comparator + v1
reconcile; Bridge A `bRho↦B_rho_3/W_0`; sort-routing + total-bijection containers + exact-multiset accounting +
jet-conservation). **Legs (Codex-written ⇒ fresh Claude agent + Grok):** prompt
`_legs/S11c_b_adjudication_build_review.md`. Raw: Grok `~/.s11_build/S11c_b_adjudication_build_grok.txt`;
Claude-agent scratch `/tmp/claude-1000/-var-projects-toy-physics/53620ffb-.../scratchpad/adjud_leg/`.
**Both legs: SOUND, NO BLOCKING FINDINGS. The 38 ALGEBRAIC MATCH are all GENUINE — no false MATCH.**

## Route accounting (both legs reproduced the builder's exactly)
`join=147 + py_only=36 + wl_only=48 = 231`; `ALGEBRAIC 38 MATCH / 36 FLAG` + `CONTAINER 0 MATCH / 4 FLAG`
(KINETIC) + `STRUCTURE_INCOMPLETE 69` + `COVERAGE 84` = 231; `CASE_ID_MULTISET_EQUAL=true`; 12
`NAMESPACE_INCOMPLETE`; `JET_CONSERVED=231 / JET_LOST=0`.

## The 38 MATCH — independently located + reduced by BOTH legs (identical)
`ADMISSIBILITY_OPERATOR_OPERAND ×16` (BODY_FORCE 8 + PER_FACE_TRACTION 8) + `ADMISSIBILITY_SUPPORT_OPERAND ×20`
+ `ENERGY_BASIS_COUNT ×2` (Integer(26)==26). SLAB_OPERATOR, COUPLING_KERNEL, COUPLING_KERNEL_TERM_ORIGINS,
MU_THETA_OPERATOR, SLAB_OPERATOR_TERM_ORIGINS carry **ZERO** ALGEBRAIC MATCH.
- Both legs hand-reduced representative MATCHes (BODY_FORCE scalars, a PER_FACE 4-tuple, ENERGY_BASIS_COUNT):
  each reduces to exactly 0 using ONLY {cited renames, Bridge A, integral linearity}; corroborated by random-
  rational numeric zero-tests. The ADMISSIBILITY BODY_FORCE MATCHes are byte-identical raw (0 renames applied).
- `--drop-rename LWidth`: MATCH 38→26 (12 flip to FLAG) — the rename is surgical (a genuine same-variable
  identification: PY `L_W` vs WL `LWidth` are the two different engines, and only the cited rename zeroes it).
- `--collapse-jet w1_profile_d1=w1_profile`: fires `JET_LOST` and changes the MU_THETA residual — jets are
  load-bearing, not pre-collapsed.

## ⭐ KEY honest finding — Bridge A creates ZERO MATCH (report in the step record)
`--drop-bridge-a` leaves MATCH at 38 (both legs): NO MATCH operand contains a `bRho` atom, so the 38 agreements
are RENAME-LEVEL (already latent in v1) and Bridge A cannot manufacture them. Bridge A IS load-bearing
elsewhere — on MU_THETA / SLAB it removes `bRho` and changes the (still-FLAG) residual — and is source-correct
(`B_rho_3=bRho·W_0`: WL:472/1621-1630, PY:1130-1140, spec:102). So Bridge A normalizes FLAG residuals; it
resolves no case to agreement.

## Routing / protection / accounting / hygiene — all cleared (both legs)
- ALGEBRAIC entered only for arithmetic `sp.Expr` (positive test; `Symbol` admitted despite Boolean
  inheritance; `Eq/Ne/Equivalent/Not/And`, `TextAtom` excluded). ADMISSIBILITY_RESIDUAL (Boolean `Not`/`Equivalent`
  vs Expr) and ADJOINTNESS_RELATION (TextAtom vs Equality) → `STRUCTURE_INCOMPLETE`, not folded.
- KINETIC CONTAINER = total bijection over all 4 cited leaves (`U_MOMENTUM_ROWS[0/1/2]`+`THICKNESS_ROW`,
  PY:1573 / WL:851); all 4 are CONTAINER FLAG (no subset MATCH).
- Protected: 0 `gamma_s11cb_*`, 0 `bRho` in any MATCH operand; 07/10 PY-only unjoinable; ENERGY_BASIS reps →
  STRUCTURE_INCOMPLETE. None folded/adapted into a MATCH.
- Accounting: EMITTED=CLASSIFIED=231; duplicate-key probe RAISES (no silent drop); all hook error-paths exit 2.
- Hygiene: no assert/PASS/FAIL/VERDICT/target on measured payloads.

## Non-blocking observations (carry forward)
1. `--residual-leaf-budget` is DECORATIVE (parsed/printed, never passed to the residual computation; routed
   arithmetic uses deterministic `H._residual_value`). ⇒ the layer is REPRODUCIBLE run-to-run (resolves the
   Phase-1 leaf-budget nondeterminism worry). Optional: label the printed line informational-only.
2. The full run is ~797s (>600s foreground ceiling); the agent used the harness's native background+notify (no
   setsid/&/monitor loop) and scoped per-case deltas via the layer's exact functions — no stall.

## Verdict + what remains
Layer SOUND; the 38 MATCH are genuine cross-engine agreements; no false MATCH. The engines AGREE on the
admissibility operator/support/count and DIFFER (genuine post-bridge FLAG) on the coupling kernel + slab
operator rows + mu_theta + kinetic. ⇒ ORCHESTRATOR ADJUDICATION (rule 13): characterize each FLAG residual as a
genuine one-engine finding vs an unbridged representation, then the honest step record. STRUCTURE_INCOMPLETE (69)
+ NAMESPACE (12) are owed (adjointness, ENERGY_BASIS quotient, DIMENSIONS/HOMOGENEITY shapes, controls).
