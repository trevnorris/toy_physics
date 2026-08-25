# S11c-a T7 — feasibility / adjudication MATRIX (step 0 of the reconciliation plan)

This is the systematic per-tag cross-engine audit the revised plan calls for
(`S11c_a_comparator_reemit_plan.md` §0). It **characterises** where the two independently-built engines
diverge and **classifies** each divergence structurally. It does **not** decide which engine is
physically correct — that is step 1 (adjudication vs the spec), after two independent legs. No expected
cross-engine RESULT is asserted anywhere (rule 5).

Every count / axis-set is a COMPUTED object (rule 2). Commands + literal stdout are in
`_measurements/S11c_a_T7_adjudication_matrix.md`; census apparatus + frozen stdout under
`_measurements/s11ca_t7_census/`. **Round-1 fold applied** — both legs (Codex + Grok) reproduced the
primary census and found the classification of the DERIVED/CONTROL families understated; every fold
below was re-verified by the orchestrator against the streams (rule 13). See "Round-1 leg fold" at end.

## Inputs (hash-locked)
- PY tag stream `~/.s11_build/S11c_a_sympy_engine.out`
  sha256 `6386471555b1e99d0aeb0f716eea30f839d59be50c0cedd4677ea7b376b79129` (fresh run of `9b6438fa`).
- WL tag stream `research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out`
  sha256 `82062bd36cfb07b1f18631077f0c63ac1cbce7834967686f680fa9f30019e4ec` (committed `ddecdbc2`).
- Engine sources: PY `scripts/S11c_a_interface_geometry_sympy_audit.py` (2096 ln, `9b6438fa`);
  WL `mathematica/S11c_a_interface_geometry_mathematica_audit.wl` (1814 ln, `ddecdbc2`).
- Join: **39 tag stems exact 1-to-1** (8 PY-local + 1 WL-local excluded).

## The decisive computed object: per-tag LEAF-CASE count + AXIS-SET
Raw top-level key count is NOT comparable (PY nests quantity-grouped tags one level; WL flattens every
axis into one pipe-string). The census flattens both to leaf cases and classifies every key token into a
case AXIS ∈ {QUANTITY, BRANCH, DENSITY, FACE, DOF, VIRTUAL_DOF, DIRECTION, FIELD, ORIGIN}. Full table:
`_measurements/s11ca_t7_census/axis_census.stdout`; leaf keys: `axis_census.keys.txt`.
⚠ **Census classification limit (fold):** the orchestrator census maps PY integer face `±1`→FACE but PY
integer DIRECTION `1/2/3`→`OTHER:2`,`OTHER:3` — so it under-reports a shared DIRECTION axis on the
control tags. Both legs' independent censuses classify PY integer directions as DIRECTION; that is
correct. Counts are unaffected; only the axis-diff labels for Family H were.

## Classification buckets (structural — step 0)
- **AGREE-STRUCTURE** — identical axis-set AND leaf count. (Leaf *semantics/content* still to be checked
  — Family F; an agreeing grid can hide a content divergence.)
- **AXIS-DIVERGENCE** — one engine carries a case axis the other collapses. A genuine case-set
  difference. Step-1 sub-question: is the extra axis REAL (object differs across it ⇒ the other engine is
  MISSING cases) or REDUNDANT (byte-identical across it ⇒ dup-axis)? ⛔ Not decided here.
- **MISSING-COVERAGE** — one engine emits a whole QUANTITY (or face/DOF slice) that the other omits, where
  the spec requires it. A genuine finding, ⛔ not serialization.
- **DECOMPOSITION** — one engine itemises an aggregate into sub-keys (FIELD/ORIGIN), and/or stores extra
  source/provenance leaves. Reconcilable ONLY by a declared algebraic reconstruction identity for the
  itemised part; ⛔ WL-only provenance leaves are NOT reconstructible from a PY value and must be tracked
  separately. Until shown, candidate serialization, unverified.
- **CONTENT-DIVERGENCE** — grid agrees but the stored leaf content differs (representation and/or an
  extra/absent field). Includes the leaf-representation split (Family F).
- **BESPOKE** — a record/premise/global object, not a case grid.
- **CONSEQUENCE** — a control/aggregate whose case set is inherited from its parent primaries; resolves
  once ALL its parents are adjudicated. ⛔ Only valid when no INDEPENDENT divergence is also present.

---

## FAMILY A — AGREE-STRUCTURE grid (11 tags): axis-set + count identical
`BACKGROUND_STATE` (4) · `CLOSURE_SHAPE_DERIV` (16) · `CONORMAL_DERIV` (8) · `EVOLUTION_MASS_BALANCE` (8)
· `FACE_MEASURE_SHAPE_DERIV` (8) · `FACE_NORMAL` (8) · `FACE_VELOCITY` (8) · `TRACTION` (16) ·
`VIRTUAL_CONSTRAINT` (8) · **`FACE_MAP_LAB_HELD` (2)** · **`FACE_MAP_MATERIAL_ADVECTED` (2)**
(the two FACE_MAP tags moved here from BESPOKE — fold: both are 2 leaves keyed by FACE, PY integer `±1`
vs WL `PLUS/MINUS`; the orchestrator flattener returned BESPOKE only because it rejects PY singleton
integer keys. Grid agrees; leaf-schema differs — a Family-F concern, not bespoke).

Bucket: **AGREE-STRUCTURE grid**, pending the Family-F leaf check. ⚠ These certify the CASE GRIDS align,
⛔ NOT that stored objects match. Two members carry a **CONTENT-DIVERGENCE** under the agreeing grid:
- **`BACKGROUND_STATE` (fold, verified; refined in step-1 verdicts):** WL stores fields `BOUNDARY_LOADS,
  AFFINITY_ZERO, FACE_FLUX_ZERO, FACE_VELOCITY_ZERO, THETA_ZERO` (+ W_BG/MU_R_BG/RHOBR_BG). ⚠ PY **already
  emits** the zero-conditions (`θ⁰/V_0/J_0/A_0 = 0`); it is missing only `BOUNDARY_LOADS` (the earlier
  "carries none of the … fields" was corrected in Verdict BG). Likely a cross-tag PLACEMENT difference (PY
  carries the admissibility conditions in `ADMISSIBILITY_PREMISE`), ⛔ not proof PY omits them — step 1
  confirms placement. Classify CONTENT-DIVERGENCE, not 1:1.
- All Family-A tags remain subject to Family F (leaf representation).

## FAMILY B — DENSITY-AXIS divergence (systemic; both directions)
| tag(s) | PY axes | WL axes | axis-diff |
|---|---|---|---|
| `KINEMATIC_BALANCE`, `RELATIVE_FLUX` | B×FACE×DOF (8) | B×DENS×FACE×DOF (16) | **WL carries DENSITY, PY does not** |
| `PROJECTION_SHAPE_DERIV`, `_STATIC_OPERAND`, `_DYNAMIC_OPERAND`, `_RESIDUAL` | B×DENS×DOF (8) | B×DOF (4) | **PY carries DENSITY, WL does not** |

Root cause (PY): `scripts/…sympy_audit.py:916` & `:1466` —
`representatives = DENSITY_REPS if quantity in {"TRACTION","VIRTUAL_WORK_SHAPE_DERIV","CLOSURE_SHAPE_DERIV"} else ("RHO4_CONSTANT",)`
— a hardcoded per-quantity whitelist; projection is handled by a separate path (`:1037` `rho_shape` vs
`rho0` `:1039`) that DOES vary density. PY's density policy is itself non-uniform across quantities.
Root cause (WL): density profiles `mathematica/…audit.wl:846-874`; flux/kinematic emit both reps.

Bucket: **AXIS-DIVERGENCE** (genuine; ⛔ not serialization). Step-1 test (do NOT run here): per law, take
the CAS residual of the object under `RHO4_CONSTANT` vs `RHOBR_CONSTANT`. ⚠ **Narrowed (fold):** a
CANCELLING residual proves the two reps are semantically redundant, but does **not** by itself prove that
OMITTING a spec-required case key is correct — §4/§5 must still say whether the axis is required. A
NON-cancelling residual ⇒ the omitting engine is MISSING cases. ⛔ Outcome not assumed (rule 5).

## FAMILY C — VIRTUAL_DOF divergence (one root cause; controls inherit B **and** C)
| tag | PY | WL | axis-diff |
|---|---|---|---|
| `VIRTUAL_WORK_SHAPE_DERIV` | 8, B×DENS×DOF (virtual DOF tied to physical) | 16, B×DENS×DOF×VIRTUAL_DOF | **WL adds VIRTUAL_DOF (full physical×virtual matrix)** |

Root cause (PY): `:919-924` `cases[(branch, dof, representative)]` — one virtual DOF per physical DOF
(diagonal). WL: `:1051` `{physicalDof, dofNames}, {virtualDof, dofNames}`, key `|VIRTUAL_DOF_` (`:1039`)
— full 2×2 incl. 8 off-diagonal `DOF_x|VIRTUAL_DOF_y` (x≠y) absent in PY.

Bucket: **AXIS-DIVERGENCE / SEMANTIC** (physics). Step-1 vs §4 T-d / §5: full bilinear physical×virtual
form (WL) or diagonal contraction (PY)?

⚠ **Corrected parentage (fold — both legs, verified).** `REP_INVARIANCE_RESIDUAL`,
`REP_INVARIANCE_EULERIAN_OPERAND`, `REP_INVARIANCE_MATERIAL_OPERAND`,
`CONTROL_INDEPENDENCE_{BASE,CORRUPTED}_OPERAND`, `CONTROL_INDEPENDENCE_RESIDUAL` (PY 64 / WL 80) inherit
**BOTH** divergences, not virtual-work alone. Verified per-quantity split of the 64→80 gap
(`REP_INVARIANCE_RESIDUAL`): `RELATIVE_FLUX` 8→16 (**+8 DENSITY, Family B**); `VIRTUAL_WORK_SHAPE_DERIV`
8→16 (**+8 VIRTUAL_DOF, Family C**); the other four quantities match. ⛔ **The virtual-work adjudication
alone does NOT clear these families** — a genuine `RELATIVE_FLUX` density disagreement lives inside them.
Bucket: **CONSEQUENCE(B+C)**; re-census after BOTH are resolved.

## FAMILY D — WL FIELD-explosion (DECOMPOSITION)
| tag | PY | WL | axis-diff |
|---|---|---|---|
| `FACE_SHIFT` | 8, B×FACE×DOF | 80, B×FACE×DOF×FIELD | WL itemises the face-shift into per-FIELD components |

Fold (both legs, verified): PY `VALUE` is 4 groups (scalar pressure + 4-vector velocity + scalar density +
4-vector current) ↔ WL's 10 `FIELD_*` leaves — CONSISTENT with a lossless derivative-level decomposition
(⛔ still requires a declared reconstruction identity, not byte-compare). ⚠ WL FIELD leaves additionally
carry exact-source/provenance operands with no PY counterpart — track those separately (not reconstructible
from a PY value). Bucket: **DECOMPOSITION + WL-only provenance**.

## FAMILY E — ORIGIN repartition (DECOMPOSITION, sets differ)
| tag | PY origins | WL origins | note |
|---|---|---|---|
| `EVOLUTION_TERM_ORIGINS` (PY 8 / WL 24) | **4**: DENSITY_TIME, VELOCITY_DILATATION, BACKGROUND_ADVECTION, TRUE_AREA_FACE_FLUX | **3**: ORIGIN_TIME_DERIVATIVE, ORIGIN_INPLANE_TRANSPORT, ORIGIN_TRUE_AREA_FACE_FLUX | **4↔3 REPARTITION** — PY splits in-plane into dilatation+advection; WL keeps one INPLANE_TRANSPORT |
| `PROJECTION_TERM_ORIGINS` (PY 8 / WL 20) | 5 dynamic/static categories | 5 categories under different names | + Family-B DENSITY divergence (PY carries it, WL not) |

Fold (both legs, verified): ⛔ this is NOT "WL itemises a PY aggregate" — the origin SETS differ (4 vs 3),
and PY is the finer partition. Bucket: **DECOMPOSITION (repartition)** — reconciling needs a declared
mapping of the two origin partitions, ⛔ not relabeling; and WL origin leaves carry extra source operands.

## FAMILY F — leaf-SEMANTICS / representation (orthogonal to grid; applies to A–E)
The census compares CASE GRIDS only. Fold (both legs sampled every Family-A leaf, verified):
- **PY explicit graded `(background, ε·derivative, total)` triple with `epsilon_shape`:** ONLY
  `FACE_NORMAL`, `CONORMAL_DERIV`, `FACE_MEASURE_SHAPE_DERIV`.
- **Other Family-A tags** (`TRACTION`, `FACE_VELOCITY`, `CLOSURE_SHAPE_DERIV`, `EVOLUTION_MASS_BALANCE`,
  `VIRTUAL_CONSTRAINT`, `BACKGROUND_STATE`): PY `VALUE` is a Mul/Matrix with `epsilon_shape` baked in —
  NOT the 3-tuple. (⇒ the earlier blanket "each shape-derivative is a graded tuple" was overstated.)
- **WL** stores the derivative COEFFICIENT (waveOrder stripped, `shapeDerivative` `:42-53`) under a nested
  `EXPRESSION`, order in `MULTIGRADE_EPSILON_ETA_SIGMAW`. Field vocab: PY {VALUE, MULTIGRADE,
  DIMENSION_L_T_M}; WL {EXPRESSION, EXACT_SOURCE, EXACT_TRUE_AREA_SOURCE, MULTIGRADE_EPSILON_ETA_SIGMAW,
  DIMENSION_L_T_M, SHAPE_DERIVATIVE, + tag-specific GRAPH_EVALUATED_SOURCE, ORIENTATION_OBJECT,
  SLAB_VELOCITY, closure operands}.

Bucket: **CONTENT-DIVERGENCE**, split into (i) coefficient↔graded-tuple **normalization** (reconcilable
by a declared algebraic reconstruction of the graded value from coefficient+order — ⛔ not relabeling) and
(ii) **WL-only source/provenance content** (EXACT_SOURCE etc.) that a graded-tuple reconstruction does
NOT produce. ⛔ No Family-A tag is 1:1-serialization on grid agreement alone.

## FAMILY G — BESPOKE records (individual)
| tag | structure | note |
|---|---|---|
| `ADMISSIBILITY_PREMISE` | PY 4-tuple record `(SUPPORT_STABILISED_BACKGROUND, (f_hold_0,t_hold_±_0), STATIONARITY_NOT_TESTED…, S11CB_…_RESERVED)`; WL 2-entry | supplied premise; align as a record. ⚠ may hold the boundary-load conditions WL puts in BACKGROUND_STATE (Family A) — check placement in step 1 |
| `DIMENSIONS` | PY 38 flat quantity→dimension entries; WL 4 field-keyed | different granularity; whole-run dimension object |
| `BACKGROUND_DENSITY_MAP` | PY 4 (BRANCH×DENSITY); WL 2 (DENSITY) | **PY branches the map, WL does not** (small AXIS-DIVERGENCE) |
| `HOMOGENEITY_{BASE,CONTROL,RESIDUAL}_OPERAND` | PY bespoke (non-grid); WL 8 quantity-keyed (CLOSURE/EVOLUTION/PROJECTION/VIRTUAL_WORK) | §6 homogeneity control; describe both structures in step 1 |
(FACE_MAP moved to Family A.)

## FAMILY H — CONTROL families: NOT pure consequence (fold — split three ways)
| tag(s) | PY | WL | verified content of the gap |
|---|---|---|---|
| `CONTROL_FORM_{BASE,ABLATED}_OPERAND`, `CONTROL_FORM_RESIDUAL` | 480 | 744 | see split below |
| `UNIFORM_LIMIT_{S11B,S11CA}_OPERAND`, `UNIFORM_LIMIT_RESIDUAL` | 164 | 232 | see split below |

Fold (Codex, orchestrator-verified). This family is NOT merely CONSEQUENCE; it decomposes into:
1. **MISSING-COVERAGE (genuine finding).** WL form-control and WL uniform-limit each omit **5 whole
   quantities** that PY covers: `EVOLUTION_TERM_ORIGINS`, `PROJECTION_STATIC_OPERAND`,
   `PROJECTION_DYNAMIC_OPERAND`, `PROJECTION_RESIDUAL`, `PROJECTION_TERM_ORIGINS`. Verified: PY form-control
   quantity set = 18, WL = 13, WL-only = ∅ (WL is a strict subset). §5b/§5c require every T-object / S11c-a
   object ⇒ this is a coverage GAP, ⛔ NOT collapsible by density/virtual-DOF adjudication.
2. **FACE-axis divergence.** WL expands FACE where PY uses a singleton `BOTH_FACES`
   (`VIRTUAL_CONSTRAINT`, `EVOLUTION_MASS_BALANCE` 24/48 in form control; virtual work 24/96 combining
   FACE and VIRTUAL_DOF; uniform virtual work 8/32). A step-1 FACE-handling question.
3. **Inherited B + C + D** (density, virtual-dof, field-explosion) via the primaries the controls consume.
   Plus: **DIRECTION is SHARED** — both engines carry it (PY integer `1/2/3`, WL `DIRECTION_1/2/3`); drop
   DIRECTION from any "WL-only" set. `UNIFORM_LIMIT` also carries WL special tokens `SIGMA_E_ZERO`,
   `GRADIENT_SIGMA_E_ZERO` (uniform smoke-test values) and PY `*_TERM_ORIGINS` tokens.

Bucket: **MISSING-COVERAGE + FACE-axis + CONSEQUENCE(B+C+D)**. ⛔ Do not close on B/C alone; re-census
after items 1-3. Producers: PY form-control `:1406`; WL form-control `:1379`; WL uniform `:1537`.

---

## Summary of what step 1 must adjudicate (physics-bearing; ⛔ not decided here)
1. **Density axis (Family B + PROJECTION_TERM_ORIGINS):** RHO4-vs-RHOBR residual per law; §4/§5 fixes the
   normative policy. Both directions occur. (Redundant residual ≠ licence to drop a spec-required key.)
2. **Virtual-work matrix (Family C):** full bilinear (WL) vs diagonal (PY) — §4 T-d / §5.
3. **Control MISSING-COVERAGE (Family H.1):** are the 5 WL-absent quantities required by §5b/§5c? If so WL
   is MISSING cases — a build defect, not serialization.
4. **FACE handling (Family H.2):** BOTH_FACES singleton (PY) vs expanded FACE (WL) in the controls.
5. **Leaf semantics + content (Family F, BACKGROUND_STATE):** graded-tuple↔coefficient reconstruction
   identity; WL-only provenance leaves; boundary-load placement (BACKGROUND_STATE vs ADMISSIBILITY).
6. **Decompositions (D, E):** reconstruction identity (FIELD) / origin-partition mapping (4↔3 ORIGIN).
7. **Bespoke (G):** BACKGROUND_DENSITY_MAP branch axis (PY-only); DIMENSIONS/HOMOGENEITY granularity.
8. **Consequences:** re-census the control families after 1-4; confirm what collapses.

## Round-1 leg fold (provenance — both legs, orchestrator-verified)
Both legs independently re-derived the census (scripts + stdout under `~/.s11_build/`): Grok
`s11ca_t7_independent_census.{py,stdout}`; Codex `s11ca_t7_codex_independent_census.{py,stdout}`. Both
REPRODUCED the primary Family-A/B/C counts and the join (39, none missing). Findings folded above (each
re-verified against the streams — commands in the twin): Family C parentage B+C (both); Family H
MISSING-COVERAGE of 5 quantities + FACE-axis + shared DIRECTION (Codex); Family F three-way leaf split +
WL-only provenance (both); BACKGROUND_STATE content divergence + FACE_MAP-not-bespoke (Codex); EVOLUTION
origins 4↔3 repartition (both); tag-name precision + density-test decidability narrowing (Codex). No
rule-5 leak found by either leg. Rule-2 twin gaps flagged by Codex are repaired in the twin.

## What is NOT yet measured (honest gaps for step 1)
- The density RHO4-vs-RHOBR residual (item 1) is NAMED, ⛔ deliberately NOT run (step 1).
- BACKGROUND_STATE boundary-load PLACEMENT (is it in PY's ADMISSIBILITY_PREMISE?) — described, not yet
  cross-tag-verified.
- HOMOGENEITY PY structure is from raw heads (flattener returned BESPOKE), not fully parsed.
