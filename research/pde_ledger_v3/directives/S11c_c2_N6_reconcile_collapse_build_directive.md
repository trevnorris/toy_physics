# Build directive — S11c-c2 N6 reconcile COLLAPSE instrument (v3)

⚠ **Codex-authored build directive.** It governs the CODE astra writes; astra does NOT decide physics. This
instrument implements the reconcile collapse TEST framed + vetted in
`_measurements/S11c_c2_N6_reconcile_question.md`. ⛔ The instrument PRINTS; it decides NOTHING.

**v3 — Codex-authored (rule 15); folds round-2 MUST findings.**

The frozen energy-basis change and census crosswalk are now written below, the nonexistent density rewrite is
removed, every primitive has a typed scope, and the controls include locality and dropped-primitive sensitivity.

## Deliverable

`research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_collapse.py` (+ a `test_…` alongside).

## Object

Per surfaced cross-engine family, per matched key, per retained grade `(η^i σ_W^j), i,j∈{0,1}`, per fixed
`(anchoring,density)` case: apply the one frozen, typed bridge dictionary of name/structure/normalization
correspondences (§ Dictionary) to the two engines' emitted operands, then PRINT
`bridge(operand_SymPy) − bridge(operand_WL)` and a three-valued collapse witness. ⛔ Never a verdict, ⛔ never an
assertion.

Surfaced families tested:

- `N6RC_CARRIER_EULERIAN`, `N6RC_CARRIER_MATERIAL` (channel a);
- `N6COV_SOURCE_ACTUAL`, `N6COV_SOURCE_PREDICTED`, `N6COV_FROZEN_PHI` (channel b — Φ is a collapse operand;
  see § Φ);
- `N6COV_SOURCE_BASELINE`, computed and labelled the nominal-control duplicate of `SOURCE_ACTUAL` at the
  engine settings `kappa_a=1, kappa_j=0` / `actualAdvectionCoefficient=1, actualJunkCoefficient=0`
  (`S11c_c2_N6_covariance_sympy.py:34-35,137-154`; `S11c_c2_N6_mathematica_audit.wl:18-20,845-849`), not a
  third independent discriminator.

Also emit the control/premise census `N6COV_ACTUAL_CONTROL_PARAMETERS`, `N6COV_PHI_DOMAIN_CENSUS`, and the
separate live-density premise comparison (§ Census).

## Inputs (all committed; content present on disk)

- `DEFAULT_PY` = the three SymPy `.out` files
  (`scripts/out/S11c_c2_N6_{covariance,reconcile,diagnostic}_sympy.out`, `7c0790ab`); `DEFAULT_WL` =
  `mathematica/out/S11c_c2_N6_mathematica_audit.out` (`ae73b884`).
- The two engines + shared modules. All frozen physics maps and scales below are grounded only in these:
  `scripts/S11c_c2_N6_{covariance,reconcile,diagnostic}_sympy.py`,
  `scripts/S11c_c2_selfenergy_fold_sympy_audit.py`, `scripts/S11c_b_brane_operator_sympy_audit.py`, and
  `mathematica/S11c_c2_N6_mathematica_audit.wl`.
- `scripts/S11c_c2_N6_cross_engine_comparator.py` — reuse its extraction primitives.

## Operand-pair extraction (⛔ do NOT re-implement the join/schema bridge; ⛔ do NOT edit the comparator)

`compare_family` materializes operands, computes the residual, emits `CASE`, then releases them
(`S11c_c2_N6_cross_engine_comparator.py:794-851`); `object_work` returns only accounting (`:862-898`). Therefore
do not attempt to recover operand pairs from the `CASE` path. Import the comparator's extraction primitives —
`load_py_jsonl`, `load_wl`, `verified_spelling_maps`/`ACTIVE_NAMES`/`POINT_NAMES`, typed leaf decoders,
`extract_py_numeric`, `extract_wl_numeric`, `extract_meta`, and `materialize` — and, scoped to the families above,
group matched leaves by the comparator's own typed key before any residual. The only new layers are the frozen
dictionary and witness. If a needed primitive cannot be imported without editing the comparator, STOP and
report; a comparator helper is a separate reviewed patch, not an in-band edit.

## THE FROZEN DICTIONARY — exact primitives, direction, and typed scope

⛔ Only the primitive correspondences declared here are allowed. The production bridge first applies the
comparator's verified mechanical WL→SymPy spelling map, then the typed primitives below in their stated
canonical direction. It never solves or matches against a carrier, μ, source, Φ image, or residual. Emit every
declared primitive as `{map_id, primitive_id, wl_domain, canonical_image, scope, stage, engine_sites}` and emit
its pre-rewrite occurrence census.

### Map 1 — grad-θ jet duals

Freeze, for `i=1,2,3`, the canonical atom relation
`thetaJet<i> (WL) → theta_d<i> (comparator spelling) ← grad_theta_<i> (SymPy operator spelling)`.
The SymPy bridge is explicitly emitted by `S11c_c2_N6_reconcile_sympy.py:59-69`; the operator atoms are declared
at `S11c_b_brane_operator_sympy_audit.py:227-233`; `wave_jet` sends `grad_theta_<i>` to `theta_d<i>` at
`S11c_c2_selfenergy_fold_sympy_audit.py:139-169`; WL constructs the same physical jet through
`jet["theta",{i}]` and its bridge at `S11c_c2_N6_mathematica_audit.wl:84-89,945-948`. This is an atom map, never
a field or carrier equality.

### Map 2 — source-wave point-preserving map

For every source wave jet that actually occurs, freeze
`jet[f,I] (WL bare atom) → ∂_I s11cc2Field_f(Y,t) (SymPy applied jet)` at the same point `Y`, preserving `f`, the
sorted spatial multi-index, time-derivative count, and all arguments. WL `sourceMap` retains registered bare jet
atoms (`S11c_c2_N6_mathematica_audit.wl:707-711`); SymPy `source_value` applies `wave_jet(...,Y)`
(`S11c_c2_N6_diagnostic_sympy.py:392-396`), whose field/derivative construction is
`S11c_c2_selfenergy_fold_sympy_audit.py:133-173`. The comparator deliberately preserves applied heads and
arguments (`S11c_c2_N6_cross_engine_comparator.py:287-289,990-993`).

### Map 3 — profile jets and constructor scales

Freeze the multi-index spelling `w1ProfileJet<I> → w1_profile_d<I>` and
`m1ProfileJet<I> → m1_profile_d<I>` (sorted spatial `I`). The exact constructor scale for `r=|I|≥1` is

`D_I WBg = sigmaW·w1ProfileJet<I>/LW^(r-1)` and
`D_I muRBg = sigmaW·(muR/W0)·m1ProfileJet<I>/LW^(r-1)`,

with canonical SymPy spellings `W_bg`, `sigma_W`, `L_W`, `mu_R`, and `W_0`. These formulas are the WL
`profileRules` at `S11c_c2_N6_mathematica_audit.wl:127-132` and the SymPy profile/jet constructors at
`S11c_b_brane_operator_sympy_audit.py:767-778,900-908`. For numeric leaves, which the engines already graded,
apply only the surviving profile-jet atom/scalar spelling; never redo zero-jet profile expansion. For Φ only,
the pre-grade path may apply the engine formulas `WBg→W0(1+etaBg·w1Profile)` and
`muRBg→muR(1+etaBg·m1Profile)` before extracting the four coefficients
(`S11c_c2_N6_mathematica_audit.wl:127-141`; `S11c_b_brane_operator_sympy_audit.py:900-908`); never apply those
zero-jet formulas to an already graded leaf.

### Map 4 — frozen retained energy-basis change (complete table)

This map is physics-bearing normalization, not mechanical spelling. Its canonical production direction is
**WL coefficient atom → the SymPy coefficient expression in the last column**. Preserve all displayed factors;
do not simplify them into a different convention.

Define `R_W := W_0/W_bg` (WL spelling before the comparator fold: `W0/WBg`). The constructors use different
thickness bases: WL contracts `eW` and its jets directly (`S11c_c2_N6_mathematica_audit.wl:182-186,211-226`),
whereas SymPy inserts the local thickness
`E := W_0·e_W/W_bg` and `D_i E = R_W·(e_W_di − e_W·W_bg_di/W_bg)`
(`S11c_b_brane_operator_sympy_audit.py:754-755,782-830,1631-1665`). Collecting the SymPy density on the WL
retained contraction basis, at the engines' retained first-background-jet order, gives the table below. Terms
with two background first jets are outside `σ_W^0,σ_W^1`; that observation is used only to derive this frozen
table — production never runs a truncation/projector rewrite.

Notation: `q=∇θ`, `r=∇e_W`, `G_ai=∂_i u_a`, `gW=∇W_bg`, `gM=∇mu_R_bg`, and `trG=G_ii`.

| ID | retained WL contraction (coefficient factored out) | WL coefficient in density | SymPy contraction/coefficient source | frozen WL→SymPy coefficient image |
|---|---|---|---|---|
| E01 | `WBg·θ²` | `bRho/2` | `θ² : B_rho_3·W_bg/(2·W_0)` | `bRho → B_rho_3/W_0` (equiv. `B_rho_3=bRho·W_0`); both `1/2` factors remain |
| E02 | `WBg·θ·eW` | `cCoupling` | `θ·E : C·W_bg` | `cCoupling → R_W·C` |
| E03 | `θ·(r·gM)` | `energyCoefficient3` | `MU_R_BG/FIRST_JET_CONTRACTION_13 : gamma_s11cb_mu_r_bg_13` | `energyCoefficient3 → R_W·gamma_s11cb_mu_r_bg_13` |
| E04 | `r·q` | `energyCoefficient4` | `q·D(E) : kappa_theta_W` | `energyCoefficient4 → R_W·kappa_theta_W` |
| E05 | `eW·(gM·q)` | `energyCoefficient5` | `MU_R_BG/FIRST_JET_CONTRACTION_14 : gamma_s11cb_mu_r_bg_14` | `energyCoefficient5 → R_W·gamma_s11cb_mu_r_bg_14` |
| E06 | `θ·(gM·q)` | `energyCoefficient6` | `MU_R_BG/FIRST_JET_CONTRACTION_12 : gamma_s11cb_mu_r_bg_12` | `energyCoefficient6 → gamma_s11cb_mu_r_bg_12` |
| E07 | `q·q` | `energyCoefficient7` | `q·q : kappa_theta/2` | `energyCoefficient7 → kappa_theta/2` |
| E08 | `θ·(gM·u)` | `energyCoefficient8` | `MU_R_BG/FIRST_JET_CONTRACTION_04 : gamma_s11cb_mu_r_bg_04` | `energyCoefficient8 → gamma_s11cb_mu_r_bg_04` |
| E09 | `θ·trG` | `energyCoefficient9` | `θ·trG : G_theta_u` | `energyCoefficient9 → G_theta_u` |
| E10 | `gM_a·q_i·G_ai` | `energyCoefficient10` | `MU_R_BG/FIRST_JET_CONTRACTION_06 : gamma_s11cb_mu_r_bg_06` | `energyCoefficient10 → gamma_s11cb_mu_r_bg_06` |
| E11 | `gM_i·q_a·G_ai` | `energyCoefficient11` | `MU_R_BG/FIRST_JET_CONTRACTION_07 : gamma_s11cb_mu_r_bg_07` | `energyCoefficient11 → gamma_s11cb_mu_r_bg_07` |
| E12 | `(gM·q)·trG` | `energyCoefficient12` | `MU_R_BG/FIRST_JET_CONTRACTION_08 : gamma_s11cb_mu_r_bg_08` | `energyCoefficient12 → gamma_s11cb_mu_r_bg_08` |
| E13 | `θ·(r·gW)` | `energyCoefficient13` | `W_BG/FIRST_JET_CONTRACTION_13 : gamma_s11cb_w_bg_13` | `energyCoefficient13 → R_W·gamma_s11cb_w_bg_13` |
| E14 | `eW·(gW·q)` | `energyCoefficient14` | `W_BG/FIRST_JET_CONTRACTION_14` plus the chain-rule part of `q·D(E)` | `energyCoefficient14 → R_W·(gamma_s11cb_w_bg_14 − kappa_theta_W/W_bg)` |
| E15 | `θ·(gW·q)` | `energyCoefficient15` | `W_BG/FIRST_JET_CONTRACTION_12 : gamma_s11cb_w_bg_12` | `energyCoefficient15 → gamma_s11cb_w_bg_12` |
| E16 | `θ·(gW·u)` | `energyCoefficient16` | `W_BG/FIRST_JET_CONTRACTION_04 : gamma_s11cb_w_bg_04` | `energyCoefficient16 → gamma_s11cb_w_bg_04` |
| E17 | `(gW·q)·trG` | `energyCoefficient17` | `W_BG/FIRST_JET_CONTRACTION_08 : gamma_s11cb_w_bg_08` | `energyCoefficient17 → gamma_s11cb_w_bg_08` |
| E18 | `gW_i·q_a·G_ai` | `energyCoefficient18` | `W_BG/FIRST_JET_CONTRACTION_07 : gamma_s11cb_w_bg_07` | `energyCoefficient18 → gamma_s11cb_w_bg_07` |
| E19 | `gW_a·q_i·G_ai` | `energyCoefficient19` | `W_BG/FIRST_JET_CONTRACTION_06 : gamma_s11cb_w_bg_06` | `energyCoefficient19 → gamma_s11cb_w_bg_06` |

The WL order and coefficient allocation are frozen by its contraction generator, quotient, retained filter, and
`known`/`energyCoefficient<i>` assignment (`S11c_c2_N6_mathematica_audit.wl:175-226`). The SymPy uniform rows
E01/E02/E04/E07/E09 come from `uniform_coefficient`
(`S11c_b_brane_operator_sympy_audit.py:1584-1623,1763-1793`); the first-jet contraction numbering and coefficient
allocation come from `NEW_COEFFICIENTS`, `enumerate_new_candidates`, and `construct_energy`
(`S11c_b_brane_operator_sympy_audit.py:355-367,1510-1527,1727-1872`). The `R_W` factors and the E14 triangular
mixing come only from the exact local-thickness substitution and its total derivative (`:754-755,782-830,
1631-1665`). Do not use the constitutive EL at `S11c_c2_N6_diagnostic_sympy.py:318-349` to build or alter this
table: μ is a tested downstream object, not a pairing oracle.

Implement the table as a fixed 19-row triangular basis-change object over `Q(W_0,W_bg)` and verify its symbolic
rank/invertibility there. The expected table is exactly the one printed above. ⛔ No runtime fit, solve, monomial
match, or residual-guided choice; ⛔ never equate the resulting energy densities or μ objects; ⛔ never expand
`R_W` in `η` inside an already graded leaf.

### Map 5 — Φ abstract-jet spelling

Freeze only domain atom naming, sorted multi-index convention, and derivative syntax for Φ. SymPy discovers and
prolongs abstract jet paths without point evaluation (`S11c_c2_N6_covariance_sympy.py:63-129`); WL constructs
the abstract `phiMap` on its registered domain (`S11c_c2_N6_mathematica_audit.wl:84-94,236-243`). Φ images stay
as operands. ⛔ No whole-Φ equality and no `Y` evaluation.

### Map 6 — REMOVED from production rewrites

There is no leftover density atom map. SymPy substitutes the imported density atom with the live expression
`inputs.density[(rho,)][1]` inside `source_terms` (`S11c_c2_N6_diagnostic_sympy.py:378-383`); WL computes
`density3=density4·WBg` and passes it into `sourceBind`, which substitutes `rhoFace→density`
(`S11c_c2_N6_mathematica_audit.wl:380-381,839-840,869-872`). Thus neither local name is an emitted operand atom.
Do not rewrite a constant to a live expression and do not add an inert density rule. Compare these constructions
only in the separately emitted live-density premise object (§ Census).

### Map 7 — literal Jacobian spelling

Freeze `1+tr(grad_u) ↔ 1+Sum_i jet["u"<>i,{i}]` only when that literal factor survives in a scoped operand.
The SymPy frozen relation is emitted at `S11c_c2_N6_reconcile_sympy.py:59-68`; WL uses the factor inside
`materialAmplitude` at `S11c_c2_N6_mathematica_audit.wl:244-247`. Do not run the degree-2 projection: that is
already a construction stage in SymPy (`S11c_b_brane_operator_sympy_audit.py:1970-1981`) and WL (`:244-247`).

### Map 8 — occurrence-conditioned ε/ω primitives

Split ε placement and `omega` realness into separate exact primitives and activate either only if its atom is
present after maps 1-7. The source-side ε normalization is explicit at
`S11c_c2_N6_diagnostic_sympy.py:378-389`; WL binds its source before kernel construction at
`S11c_c2_N6_mathematica_audit.wl:364-385`. `omega` is a retained symbolic source parameter in WL
(`S11c_c2_N6_mathematica_audit.wl:102-113`). Omit an absent primitive. Do not import a default on-shell or
Fourier-of-derivative identity.

### Typed family/stage scope (closed allow-list)

`C_E`, `C_M`, `S_A`, `S_P`, `S_B`, and `PHI` below mean exactly the six surfaced families named in § Object.
“Occurrence-gated” means a primitive may act only on leaves in the listed family/stage whose pre-rewrite operand
contains its exact declared domain. Anything outside the row is forbidden.

| map | allowed family | allowed stage/use | expressly forbidden |
|---|---|---|---|
| 1 grad-θ dual | `C_E,C_M,S_A,S_P,S_B`; `PHI` only as abstract domain/image jet spelling | materialized numeric operand after grade; Φ abstract spelling during `extract_meta` path | source-point evaluation of Φ; census values |
| 2 source-wave at `Y` | `S_A,S_P,S_B` only | materialized source operand, after source extraction | `C_E,C_M,PHI`, every census/premise object |
| 3 profile names/scales | `C_E,C_M,S_A,S_P,S_B`; `PHI` | numeric: surviving atom/scalar spelling after grade; Φ: profile substitution before its four coefficients are extracted | zero-jet profile expansion on a graded leaf |
| 4 energy basis | `C_E,C_M,S_A,S_P,S_B` only | coefficient atoms in materialized numeric operands, occurrence-gated, using the fixed table | `PHI`, census/premise, μ/energy reconstruction, runtime matching |
| 5 Φ spelling | `PHI` only | abstract map-domain keying and derivative syntax; then each graded image coefficient | every carrier/source and any `Y` evaluation |
| 7 Jacobian | `C_M,S_A,S_B` only | literal surviving factor after grade, occurrence-gated | `C_E,S_P,PHI`; degree projection |
| 8 ε placement | `S_A,S_P,S_B` only | exact atom occurrence in materialized source operand | carriers, Φ, census, kernel/on-shell rules |
| 8 ω realness | `C_E,C_M,S_A,S_P,S_B` only | assumption normalization for an occurring `omega`, never an expression identity | Φ, census, on-shell/Fourier rules |

The comparator's own mechanical key/name folds, including pressure-slot spelling, remain available for carriers
and sources (`S11c_c2_N6_cross_engine_comparator.py:287-326`). They do not enlarge this physics-map scope.

⛔ **Excluded:** whole-carrier/energy/μ/source/Φ equality; any construction-stage operation as an operand rewrite
(degree-2 projection, zero-jet profile η re-expansion, whole velocity equality, or whole source-solve-factor
equality); `LAB_HELD↔MATERIAL_ADVECTED` (anchoring is a retained key); `σ_W↔η`, `σ_W→0`, or any cross-grade
binding; default on-shell/Fourier identities. The engines build material velocity and source factors independently
(`S11c_c2_N6_reconcile_sympy.py:82-92`; `S11c_c2_N6_diagnostic_sympy.py:378-389`;
`S11c_c2_N6_mathematica_audit.wl:274,304,366-381,862-872`), so those objects stay in the residual.

## Φ as a graded collapse operand (MUST — Φ is ungraded metadata)

`FROZEN_PHI` is emitted as an ungraded substitution map (SymPy
`S11c_c2_N6_covariance_sympy.py:116-124`; WL `S11c_c2_N6_mathematica_audit.wl:975-976`), and the comparator's
`extract_meta` produces `MAP_VARIABLE`/`FIELD_PATH` keys with no grade axes
(`S11c_c2_N6_cross_engine_comparator.py:688-765`). For Φ: extract map images through `extract_meta`; apply only
the pre-grade profile substitution allowed by map 3; explicitly extract all four `(η^i σ_W^j)` coefficients;
then apply maps 1/3/5 within their Φ scopes and run the witness coefficientwise. Do not compare ungraded images,
reuse the numeric `CASE` grade path, or apply map 2.

## Retained order (⛔ coefficientwise; ⛔ no proxy; ⛔ no η re-introduction)

Operate independently on every `(η^i σ_W^j), i,j∈{0,1}`. Numeric leaves are already grade-keyed; Φ is graded
only by the preceding section. SymPy declares/extracts the rectangle at
`S11c_c2_N6_diagnostic_sympy.py:43,245-248`; WL does so at
`S11c_c2_N6_mathematica_audit.wl:123-141`; comparator keys retain both axes independently
(`S11c_c2_N6_cross_engine_comparator.py:96-97`). ⛔ No grade combining, cross-grade cancellation,
`η·σ_W→η²`, `σ_W→0`, `σ_W↔η`, or map that introduces `η` into a graded leaf. A physical σ/η relation may be
printed as a separate object after all coefficientwise comparisons; never apply it. Enforce this structurally
with the grade-combining tripwire.

## Collapse witness (three-valued; PRINT residual first)

Per leaf, after the typed dictionary: PRINT both operands and the bridged residual, then a bounded witness in
`{collapsed, residual_remains, undecided}`:

- attempt bounded symbolic `cancel`/`together`/`factor` under a per-leaf budget;
- otherwise run a fresh joint finite-field PIT using primes and seeds disjoint from both engines' sets (SymPy
  primes at `S11c_c2_N6_diagnostic_sympy.py:669`; WL primes at
  `S11c_c2_N6_mathematica_audit.wl:733`), reject denominator zeros, derive and emit the degree/exclusion bound,
  and enforce maximum attempts, wall time, and RSS;
- emit `undecided` on parse/budget/resource exhaustion and emit PIT provenance (primes, seeds, bound, draws,
  rejected denominators, and conditional δ).

The witness is not a pass/fail token and carries no interpretation. Never coerce a native Boolean into a CAS
zero witness.

## MANDATORY controls (same pipeline; PRINT, never assert)

Predeclare a defining-relation operand pair for every primitive before reading production residuals. For map 4,
that is the exact table row (E14 uses its displayed two-term triangular relation); for the other maps it is the
exact atom/factor pair stated above. Define a primitive as **active** iff its exact domain occurs in at least one
pre-rewrite surfaced leaf inside its typed scope. Emit inactive primitives honestly in the reachability census;
do not force them active and do not silently discard them.

1. **One-sided corruption + locality.** For every active primitive, replace its image on one engine side with a
   fresh same-type sentinel in a temporary control variant. Run normal extraction, typed rewrite, residual, and
   witness. Emit the baseline/control operands and residuals plus
   `control_delta = baseline_bridged_residual − corrupted_bridged_residual`. Require a nonzero control delta on
   at least one leaf in the primitive's predeclared use-set and require zero control delta on every surfaced leaf
   outside that use-set. This is the locality check against an over-broad rule.
2. **Dropped primitive.** For every active primitive, delete only that primitive in a temporary variant and feed
   its predeclared defining-relation operand through the same extraction/rewrite/witness path. Emit both operands,
   the residual, witness, and deletion delta; require a residual on that defining-relation operand. Do not use an
   already surfaced production residual as this control.
3. **Production reachability.** Emit, for every declared primitive, its typed scope, exact domain occurrence
   count, and exact leaf-key use-set before rewriting. An empty use-set is an honest inert-map record, not a
   control failure and not permission to broaden the rule.
4. **Blanket-collapse ablation.** In a temporary control variant only, add a whole-operand identification and run
   it through the same loader, typed-key grouping, rewrite, residual-print, and witness path. The ablation must
   drive the witness to `collapsed` on every ablation leaf. Emit both extracted operands, the post-ablation
   residual, and witness. This equality is unreachable from production code.
5. **Grade-combining tripwire.** Attempted zero-jet profile expansion on a graded leaf, introduction of `η`,
   binding of `σ_W`, or access to a different grade must be structurally rejected and emitted as a control
   object.

Controls report computed CAS objects and deltas only. They do not emit verdict prose and do not assert a
production residual.

## Census audit (channel d — load-bearing control/premise, not rewrite input)

The semantic crosswalk is closed and contains exactly these six pairs:

| census family | SymPy field | WL field | comparison object | engine sites |
|---|---|---|---|---|
| `ACTUAL_CONTROL_PARAMETERS` | `kappa_a` | `ADVECTION` | scalar operands + residual | `S11c_c2_N6_covariance_sympy.py:137-153`; `S11c_c2_N6_mathematica_audit.wl:977-980` |
| `ACTUAL_CONTROL_PARAMETERS` | `kappa_j` | `JUNK` | scalar operands + residual | `S11c_c2_N6_covariance_sympy.py:137-153`; `S11c_c2_N6_mathematica_audit.wl:977-980` |
| `PHI_DOMAIN_CENSUS` | `imported_theta_e_jets` | `DOMAIN` | sorted abstract-domain atom sets after map-5 spelling | `S11c_c2_N6_covariance_sympy.py:87-115`; `S11c_c2_N6_mathematica_audit.wl:236-243,850-852` |
| `PHI_DOMAIN_CENSUS` | `coverage` | `COVERAGE` | per-domain-atom incidence/path operands + structural residual | `S11c_c2_N6_covariance_sympy.py:99-115`; `S11c_c2_N6_mathematica_audit.wl:236-243,850-852` |
| `PHI_DOMAIN_CENSUS` | `uncovered` | `UNCOVERED` | sorted abstract-domain atom sets | `S11c_c2_N6_covariance_sympy.py:103-115`; `S11c_c2_N6_mathematica_audit.wl:236-243,850-852` |
| `PHI_DOMAIN_CENSUS` | `max_present_jet_rank` | `MAX_RANK` | scalar operands + residual | `S11c_c2_N6_covariance_sympy.py:103-115`; `S11c_c2_N6_mathematica_audit.wl:236-243,850-852` |

No other semantic aliases may be inferred. Existing identical field spelling/case normalization from
`extract_meta` remains mechanical and is emitted separately. Every field not paired by the table or by exact
mechanical spelling stays surfaced one-sided. In particular, WL `THICKNESS` and `MATERIAL_NORMAL` are WL-only;
emit them with `pairing=wl_only` and their raw operands. Do not create a SymPy field or zero. Likewise retain WL
`JUNK_CASE`, `JUNK_SYMBOL`, `JUNK_DIMENSIONS`, `JUNK_ASSUMPTION` and unmatched SymPy fields from the exact
engine records as one-sided unless the comparator's unchanged exact-spelling key already joins them. The source
inventories are `S11c_c2_N6_covariance_sympy.py:105-115,146-153` and
`S11c_c2_N6_mathematica_audit.wl:977-980`.

Emit a separate `LIVE_DENSITY_PREMISE` object containing the two constructor operands and their structural
comparison: SymPy's live replacement target `inputs.density[(rho,)][1]` and WL's
`density3=density4·WBg`, with sites `S11c_c2_N6_diagnostic_sympy.py:378-383` and
`S11c_c2_N6_mathematica_audit.wl:380-381,839-840`. This object is census/premise-only and must never be passed
to the operand rewrite dictionary.

## ⭐⭐⭐ THE THREE CLAUSES (verbatim — non-negotiable)

> **1. The script may PRINT computed objects. It may NOT state conclusions.** An emit payload must be a CAS
> object (an expression, a residual, a symbolic zero-witness), ⛔ never prose describing a result.
> **2. PRINT the residual; do NOT assert it.** ⛔ No `assert residual == 0`, ⛔ no residual-zero exit.
> **3. Interpretation belongs to the reconcile record.** ⛔ The script does not editorialise.

The only hand-combination of physical symbols allowed is the exact frozen dictionary above. Every production
residual and witness must be reached by computation from emitted operands and those tables.

## Value-free / leak discipline

⛔ State no expected collapse outcome, no expected witness count, no production pass condition, and no statement
about whether the shipped run collapses. Prior-run counts are historical facts, not targets; do not reproduce
them. Map-4 factors such as `W_0/W_bg` and `1/2` are constructor facts, not collapse results.

## Builder fence (astra) — build → verify → run → report → STOP

Build the instrument; verify that the deliverable and test exist and are non-empty; run the test; run the
instrument scoped to the surfaced families and four `(anchoring,density)` cases with bounded per-object time/RSS;
write its `.out`; report; STOP. Do not author or run a review of your own work, call another AI, edit the
comparator, or edit outside the deliverable, its test, and its `.out`. Independent fresh-Claude + Grok review is
the orchestrator's next step.

## Definition of Done

- The script imports comparator extraction primitives, groups matched operand pairs by the comparator's typed
  key before residual formation, and does not edit the comparator.
- It emits and applies only maps 1/2/3/4/5/7/8 with the closed scope table; map 4 is the exact 19-row table above;
  map 6 has no production rewrite; no whole-object/construction-stage/grade-combining map is reachable.
- Numeric operands remain coefficientwise; Φ uses `extract_meta` → pre-grade profile expansion → explicit four-
  coefficient extraction → scoped maps → witness.
- Every active primitive has corruption-locality and dropped-primitive objects; every declared primitive has a
  reachability record; blanket ablation and the grade tripwire use the production pipeline.
- The exact six-row census crosswalk, all one-sided fields, and `LIVE_DENSITY_PREMISE` are emitted without a
  synthesized counterpart or operand rewrite.
- Residuals precede three-valued witnesses; PIT disjointness, denominator rejection, derived bounds, finite
  attempt/time/RSS budgets, exhaustion-to-`undecided`, and provenance emission are enforced.
- `test_…` covers the map-4 table/rank and E14 mixing; typed-scope rejection (especially map 2 on Φ/carriers);
  absence of map 6; corruption locality; dropped primitives; reachability; same-path blanket ablation; grade
  tripwire; Φ grading; extraction reuse; native-Boolean rejection; and no assert/verdict/prose result path.
- A scoped `.out` is produced with no expected value encoded anywhere.

## Builder report (return this)

Return deliverable/test/output paths; token usage; the emitted primitive/site/scope table and reachability
census; neutral run accounting by family/key/grade and by witness token; corruption-locality and dropped-
primitive control deltas; census crosswalk/one-sided accounting; peak RSS/runtime; and confirmation that the
builder fence held. ⛔ No physics conclusion.
