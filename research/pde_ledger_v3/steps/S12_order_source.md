# S12 — the source required by the order balance, with separate boundary data

**Status, 2026-10-20:** record and register pass written for review. The constructions are accepted at
`a3da43d9`; the comparator and its filed production output are accepted at `fe485b73`, with the engine-row
note below. This record has not had its own review. No drain law, throat solution or flow is delivered.
The light-sector closeout, including the flow-only oracle-comparison record, is accepted at `e070c7f3`.

The result is conditional: the committed material/order balances require the printed three regional
sources and two face sources on fields satisfying **U1's support condition**, and separately require
the total-material face balances. The second order object is a labelled fast-ordering specialization.
Neither the support identification nor fast ordering is derived. Source/controller forms and boundary
data remain distinct and OPEN. (Retrieve M01–M04, M07, M11.)

## Evidence and retrieval

Paths here are relative to `research/pde_ledger_v3/`. **Retrieve Mnn** means the exact command
`python3 research/pde_ledger_v3/scripts/S12_record_measurements.py --section Mnn`, run from the repository
root. Each command and its literal stdout are in
[`_measurements/S12_record_measurements.md`](_measurements/S12_record_measurements.md).
[`../scripts/S12_record_measurements.py`](../scripts/S12_record_measurements.py) regenerates that file.
It retrieves stored text, JSON fields, named-object shapes and literal-match counts only. It imports
neither engine nor comparator and performs no algebra, probe evaluation, normalization or new check.

| Lookup | Authority retrieved |
|---|---|
| M01 | Step decision list revision 4, cleared at `b08673a2`; premises, domain, inventories and OPEN items. |
| M02 | Build rounds 0–2, comparator rounds 0–2 and their acceptance notes; light-closeout acceptance; record-contract fold. |
| M03 / M04 / M10 | Both accepted transcripts' objects / both source constructions / every dimensional object. |
| M05 | Filed comparator configuration, name table, source correspondences, printed summary and record shape. |
| M06 | Every unjoined occurrence and joined row, both leaf counts and every uncompared leaf path. |
| M07 / M08 / M09 | Every comparison with stored sample counts / every residual other than literal `"0"` and false text comparison / every explicit coverage item. |
| M11 / M12 / M13 | Inventory source rows / register schema, rest-on method and census / identified prior-art sources and historical second-route disagreement. |
| M14 / M15 / M16 | Filed run records and comparator contract / final register pass and status / literal residual and outcome counts. |

## 1. Supplied object, domain and standing

The two balances and signed rate split are **supplied, unfalsifiable in this build**, from the user's
2026-07-24 committed law (decision list, supplied equations; plan S12; v2 stage 006):

```
∂t n + ∇₄·(n u) = 0
∂t(χ_B n) + ∇₄·(χ_B n u + J_χ) = n Γ_B
Γ_B = Γ_return − Γ_drain,       [Γ_B] = T⁻¹.
```

The signed source converts order in one conserved material. Positive means bulk material re-ordering
onto the brane; negative means de-ordering. The ruled-out frozen-wall material sink plus remote return
is not this object. No static limit is posed. (Retrieve M01, M03, M11.)

The domain has real coordinates `(t,x¹,x²,x³,w)` and live faces `w₋(x,t)<w₊(x,t)`. Region 1 is below,
region 2 inside, region 3 above. Outside, order is zero; inside, `f(t,x,w)` is live with `0≤f≤1`.
The faces bound the region in which shear can be held; they are **not** defined by the support of a
chosen order profile. Density, each velocity component, each order-flux component and the interior
order have independent regional restrictions and one-sided traces, may jump, and have no
face-concentrated input part. The source `n Γ_B` can have such a part. Point values at the faces and
jumps are unspecified. (Retrieve M01, M03–M05.)

No regime is chosen and no field is frozen. In particular `J_χ=0` is not adopted. Transfer is only to
where the committed law holds and the shear-holding region has two graph faces over the in-plane
coordinates. The puncture's own interior is excluded, with S22. The S9b/O2 far-zone setting is a
separate restricted inventory item (§5); it is not imposed on these local objects. (Retrieve M01, M11.)

## 2. What the constructions print on each part of the domain

The following is compact notation for the printed expansions, not an additional calculation. Write
`n_r,u_r,J_r` for each region's fields and `D₂f=∂t f+u₂·∇₄f`. At a face `s`, let
`a_{r,s}=u_r^w−∂t w_s−u_r^i∂i w_s` and `b_{r,s}=J_r^w−J_r^i∂i w_s`, always using that
region's one-sided trace at `w=w_s`. These symbols only abbreviate terms in M03. Both deltas below
are normalized in coordinate `w`; they are not surface-area densities. (Retrieve M03, M04, M07.)

| Part | General-order source `n Γ_B` printed there | Total-material balance printed there | Cross-engine evidence |
|---|---|---|---|
| Below `w<w₋` | `∇₄·J₁` | `∂t n₁+∇₄·(n₁u₁)=0` | `general_order/region/1/source` and `material/region/1/balance`: literal exact residual `0`. |
| Inside `w₋<w<w₊` | `n₂D₂f+∇₄·J₂` | `∂t n₂+∇₄·(n₂u₂)=0` | Corresponding `region/2` rows: literal exact residual `0`. |
| Above `w>w₊` | `∇₄·J₃` | `∂t n₃+∇₄·(n₃u₃)=0` | Corresponding `region/3` rows: literal exact residual `0`. |
| Lower face | `C₋ δ(w−w₋)`, `C₋=n₂f₋a_{2,−}+b_{2,−}−b_{1,−}` | `n₂a_{2,−}−n₁a_{1,−}=0` | `face/1/trace_coefficient`, `source_contribution` and material `requirement`: supplied-domain residual `0`; as-computed step difference retained (§4). |
| Upper face | `C₊ δ(w−w₊)`, `C₊=−n₂f₊a_{2,+}+b_{3,+}−b_{2,+}` | `n₃a_{3,+}−n₂a_{2,+}=0` | Corresponding `face/2` rows: supplied-domain residual `0`; as-computed step difference retained (§4). |

**Retrieval for every table cell:** M03 contains both literal expansions; M07 contains their paired
operands, residuals and stored sample counts; M08 retains every as-computed difference. The material
requirements are printed, never solved or substituted at a face (M02, M04). Neither order face
coefficient is required to vanish. No product of a jumping order with the global material-balance
distribution is used: the engines eliminate regional bulk material balances by their different
printed routes (M03, M04; the elimination operands themselves are not a formed comparison).

The complete `general_order/source_distribution`, `regular_distribution_part` and
`source_with_face_traces` compare as three regular parts plus two normalized face parts. Their regular
parts print exact residual `0`; their face parts print supplied-domain residual `0`. This statement
does not extend to the unrestricted coefficient diagnostics or omitted leaves (§4). (Retrieve M07–M09.)

### Item 2 — labelled conditional specialization

The printed premise is **"conditional specialization: ordering keeps pace with the slab; G is general
in both offsets"**. The interior order is `G(w−w₋,w₊−w)`, with both arguments general. The interior
printed source is
`n₂[(u₂^w−∂t w₋−u₂^i∂i w₋)G₁−(u₂^w−∂t w₊−u₂^i∂i w₊)G₂]+∇₄·J₂`,
with `u₂` evaluated at the live interior point `(t,x,w)`; `G₁,G₂` denote the two printed argument derivatives of `G`.
Below and above, the sources remain the respective flux divergences. At the lower face, the printed
`f₋` is replaced by `G(0,w₊−w₋)`; at the upper face, `f₊` by `G(w₊−w₋,0)`. The material
requirements remain separately printed. The `conditional_fast_ordering/region/{1,2,3}/source` rows
have exact residual `0`; the two face traces and delta parts have supplied-domain residual `0`, with
the same step qualification. (Retrieve M03, M07, M08.)

This is **"a labelled conditional specialization"**, **"not a limit earned from item 1"**. The
additional range `0≤G≤1` is supplied. No ordering timescale, relaxation equation or derivative limit
has been earned. (Retrieve M01, M03.)

### Item 3 — dimensional output

The engines use a length/time basis with number count dimensionless. They assign density `(-4,0)`,
order `(0,0)`, velocity `(1,−1)`, faces `(1,0)` and rates `(0,−1)`. They print the inference of
order-flux dimensions as `(-3,−1)`, and dimensional operands/residuals for both balances, both source
objects, their regular parts and face terms, and the material face objects. The printed source and
balance term vectors are `(-4,−1)`; face trace term vectors are `(-3,−1)`. The additive zero has
**"no intrinsic unit"**. (Retrieve M03, M10.)

The flux units are solved from the order balance: homogeneity is **definitional for `J_χ` and
informative for the other terms**, as the acceptance notes state. The pre-inference residuals involving
the unknown flux exponents are not a failed homogeneity result. Cross-engine comparison covers only
the units of term-bound ledger components and other explicitly paired values; the unmatched/nested
ledger limits below remain. (Retrieve M02, M07, M09, M10.)

## 3. Premises and which objects depend on them

The quoted standings are those of the step decision list. It supplies no general "re-tested at knit"
standing for these six items; none is assigned here. (Retrieve M01, M02.)

| Input, source/date | Standing quoted from the list | Printed dependence |
|---|---|---|
| U1, user 2026-10-19, decision list "User premises" | **"an adopted substrate input to a conditional model"** | Its support reading determines the zero exterior order and live interior order in both packages, hence their regional and face source objects and their dimensional reports. |
| U1 support reading, same date/source, "How U1 enters" | **"an identification, flagged"**; **"a condition on which order fields are allowed"**, **"not a formula for the order"** | The shear-holding region is independently bounded by live faces. If order persists outside, or extent is defined by order, this identification fails. It is not a derived shear/material law. |
| U1 re-ordering clause, same user/source/date | **"Not encoded"**; **"an OPEN admissibility requirement"** | Enters no computed object. Support alone does not enforce that settling bulk material re-orders. No local sign or rate follows. |
| U2, user 2026-10-19, same premise section | **"an adopted substrate input to a conditional model"** | No sign restriction is imposed on either source package or its faces. The supplied signed split stays printed; no sign-selection calculation is performed. |
| U3, user 2026-10-19, same section | **"an adopted substrate input to a conditional model"**; **"a statement about the bulk's IR condition"** | Enters no computed local source object. It is carried in the boundary/domain inventory, not a near-face rest condition or an energy budget. |
| U4, user 2026-10-19, same section | **"an adopted substrate input to a conditional model"**; **"a sign premise on a global total, not a rate law"** | Enters no computed object. No global integration, growth rate or spreading-versus-thickening result is printed. |
| "Fields at a face", revision 3 and revision 4, 2026-10-20; user choice and pass-4/5 dispositions named in the list | **"This holds for the fields the source is built from"**; each is **"a live function with its own one-sided value"**; **"none of these fields has a part concentrated on a face itself"** | Function class for regional density, velocity, flux and interior order in both constructions, including material face requirements and both order source packages. **"The rule does not apply to `n Γ_B`"**. It excludes concentrated inputs, not computed face sources. |
| Fast ordering, item 2 of the decision list authored 2026-10-19; the current revision 4 is dated/cleared 2026-10-20 | **"A labelled conditional specialization"**; **"This assumes ordering keeps pace with the slab"** | Only `conditional_fast_ordering/*` and its dimensions. General-order and material objects do not depend on fast ordering. |

The 2026-10-20 user decision **"Leave it out for now"** applies to order carried along a face.
All face-concentrated parts of the four input fields remain an OPEN possibility outside this object.
The source is reported as `n Γ_B`; no division by a possibly jumping or vanishing density is performed.
(Retrieve M01–M04.)

## 4. Cross-engine findings and complete coverage limits

The comparator's configuration prints **"real coordinates; strict lower-face < upper-face; 0 <=
interior order <= 1, supplied by engines"**, forms `affine`, `quadratic`, `rational`, three rational
base points and seven strata: generic, each lower/upper zero-jump stratum, each lower/upper zero-order
stratum, and each lower/upper combined stratum. Its probe note is **"order profiles have a linear zero
at the selected face, remain live and separate argument derivatives there; all printed probe points
carry their supplied-domain values"**. Those are finite probes, not a proof over all live functions.
(Retrieve M05, M07.)

The printed comparison limits are retained literally:

- **"coefficient on support divided by absolute derivative of delta argument in w"**;
- **"unrestricted coefficient representation diagnostic; not a delta part"**;
- **"nested ledgers follow each engine expression tree; no source-defined cross-engine node correspondence"**;
- **"SymPy prints a distributional bulk multiple; Wolfram prints time-derivative substitution rules. These are different kinds of object."**;
- **"the Wolfram payload under this tag contains spatial currents only"**;
- **"no unique identical term binding across these expression trees"**;
- **"repeated tag; unique join rejected"** and **"tag absent in other engine; no differently named pairing"**;
- **"zero for any findings; nonzero only for operational failure"**, and **"no residual target was given"**.

These are quotations from M05–M09. A zero process exit is not agreement. Both original operands remain
in the filed output; M07 is a retrieval projection of comparison fields, not a substitute transcript.

### Every printed row class

| Stored class (number of records) | What the record can use, and its limit |
|---|---|
| `text` (111) | Labels/notes compared as text; 121 component text outcomes, with five false outcomes detailed below. |
| `value` (64) | Scalar/equation/Boolean or distribution parts with explicit operands and exact residuals; supplied-domain values remain beside them. |
| `ledger` (35), `ledger_residual` (26) | Units or unit residuals only for uniquely bound identical terms; term binding keys receive no compared-value credit. |
| `reference` (26) | Explicit reference expression and unit components, within the printed pairing. |
| `nested_ledger` (70) | Includes empty nested ledgers and the 46 rows explicitly not compared; no general node correspondence. |
| `current` (3), `regional_current` (3) | Four spatial current components per object/region; SymPy's time-density components are unpaired. |
| `tree` (26) | Explicit semantic component paths, including fields, traces and requirements; no position-only generalization. |
| `incompatible` (2) | Both material-elimination operands are printed as not formed; no equivalence is claimed. |
| `unit` (1), `inference` (1), `inference_residual` (1) | Flux unit inference and its explicit term-associated comparisons, subject to ledger accounting. |
| `configuration` (1), `test_functions` (21), `probe_domain` (315) | Bindings, finite forms/points/strata and supplied-domain probe values; they are not physical closure choices. |
| `operands` (369), `comparison` (369), `residual_before_domain` (1574) | Original inputs, stored outcomes and emission before domain restriction; repeated residual emissions are not independent evidence. |
| `accounting` (534), `summary` (1), `resources` (1) | Every joined/unjoined occurrence and leaf census; the last is operational metadata only. |

All cardinalities and class coverage above are retrieval of named stored objects (M05–M07). No class
is silently counted as a full physics comparison. In particular syntax heads count in parsed leaves,
while ledger binding keys and matrix wrapper heads receive no compared-value credit (M14).

### Every printed difference

M08 names every differing tag, component path, both operands, as-computed and supplied-domain
residual, and each stored sample whose residual is not literal `"0"`. The following exhausts its
families; Boolean `"False"` is a logical zero, not an arithmetic nonzero. No residual is simplified
or explained away by this record. (Retrieve M07, M08, M16.)

| Difference/family | Printed finding and classification |
|---|---|
| SymPy face-step factor | As-computed differences occur at both faces in the balance/source distributions, traced distributions, face contributions/trace coefficients and material requirements, including `material/face_requirements`. They print `0` on the supplied domain. The build acceptance explicitly settles the representation: SymPy carries a step in `w₊−w₋`, set to one by strict face order; Wolfram does not. **Representational on that supplied domain**, with every raw difference retained. (M02, M07, M08.) |
| Six unrestricted coefficient diagnostics | `material/face/{1,2}/coefficient_operand`, `general_order/face/{1,2}/coefficient_operand`, `conditional_fast_ordering/face/{1,2}/coefficient_operand`: nonzero expressions persist on the supplied domain in `regular/above` for face 1 and `regular/below` for face 2, and their stored rational samples differ. They are explicitly **not delta parts**. Their face-supported comparisons are reported separately, above. No cited disposition settles the unrestricted expressions' difference: **OPEN cross-engine diagnostic differences; no owner named**. (M07–M09.) |
| Domain Boolean representations | `supplied/domain`, `region_domains`, and both packages' `region/{1,2,3}/domain`: extra face-order/real-coordinate content appears in the as-computed XOR operands. Sources explicitly supply strict ordering and real coordinates (SymPy real symbols; Wolfram `Element`); the comparator prints `False` on that domain. **Representational under those supplied assumptions**, not a domain derived here. `conditional_fast_ordering/interior_range` and the supplied interior-range component already print `False` as computed. (M01, M03–M05, M07, M08.) |
| Three scope-label differences | `calculus_scope`: "regional fields and traces" differs at character 11 ("independent live restrictions and one-sided traces" / "independent live functions ... unrestricted one-sided traces"); "input singular parts" at 0 ("no input field has" / "input fields have no"); "source singular parts" at 20 ("can" / "may"). M08 retains the complete strings. **Representational wording**: the revision-4 face rule and accepted build checks explicitly settle their common function class, including source face parts. (M01, M02, M04, M08.) |
| Two `engine` text differences | Runtime text is Python `3.12.3 ...` versus Wolfram `15.0.1 ...`; the second is SymPy `1.14.0` versus platform `Linux-x86-64`. The acceptance says **"Its two parts are not counterparts, so its paired comparison carries no information, and no claim rests on it."** **Accepted pairing defect / identification metadata**, not a physics difference or evidence of equivalence. User: **"Accept, note the engine row"**, 2026-10-20. (M02, M08.) |

The exact literal census is **1514 `"0"` among 1574 as-computed residual fields**, **1556 `"0"`
and 12 `"False"` among 1574 supplied-domain residual fields**, and **98028 `"0"` among 99162
sample residual fields**; every other literal is retained in M08/M16. These counts include diagnostic
and Boolean rows and do not mean that every counted sample is an independent comparison. There are
116 true and five false text outcomes; no `unresolved_coverage` field is printed. (Retrieve M16.)

### Every unresolved coverage class

The comparator's literal summary is
`{"unjoined_py":84,"unjoined_wl":81,"joined":369,"joined_with_not_compared":46,"joined_with_not_formed":2,"joined_below_parsed_count":144}`.
These categories overlap. M06 is the complete per-occurrence/per-row list, with both parsed and
actually-compared counts and **every** uncompared leaf path; M09 is the complete component list with
reasons. Each row in those lists remains bounded by that accounting, whether or not its compared
components have zero residuals. (Retrieve M05, M06, M09.)

| Coverage class | All members and standing |
|---|---|
| Repeated progress tags | `WORKING_ON`: 54 SymPy and 56 Wolfram occurrences, rejected for unique joins. Process text, not a missing physical result. |
| Unjoined SymPy-only tags | Both order packages' `face/{1,2}/delta_operand`, `normalization_operand`, `material_extraction/coefficients`, `material_extraction/local_operands`; `material/bulk_balance_operand`; `order_flux_dimension_inference/solutions`; `primitive_function_dimensions`; each also has its printed `/comparison` note. No differently named pairing is made. |
| Unjoined Wolfram-only tags | For `material` and both order packages: `conserved_density_operand`, `region_density_operands`, with their `/comparison` notes; both order packages' `face/{1,2}/delta`, with notes; `material/time_derivative_rules` and `order_flux_dimension_inference/rules`, with notes; `execution`. |
| Not compared | 46 named nested dimensional operand/residual rows. The full list, rather than an inferred equivalence of expression trees, is M06/M09. **OPEN node correspondence; no owner named.** |
| Not formed | Both packages' `material_elimination_operand`: bulk multiple versus derivative rules. **Different printed object kinds**, as the accepted source constructions settle; no compared elimination claim. |
| Unpaired current components | `material`, `general_order`, `conditional_fast_ordering`: each `conserved_current_operand` lacks a paired SymPy time density; each `region_current_operands` lacks its three regional time densities. Sources settle that Wolfram prints only spatial currents under those tags. **Representation/coverage limit**, no unprinted join. |
| Unpaired dimension terms | M09 lists every operand/residual term with no unique identical term binding, including source, regular-part, balance, face and order-field ledgers. **OPEN term correspondence; no owner named.** The units of paired terms do not certify the missing terms. |
| Below parsed count | All 144 rows in M06 keep their actual counts, including the foregoing ledgers/current components and syntax-wrapper/binding-key omissions. No full-payload comparison is inferred from a joined name. **OPEN coverage wherever a value lacks a counterpart; no owner named.** |

No isolated-distribution coverage is printed in any of the 49 `py_coverage` / `wl_coverage` pairs
(each stored list is empty). That finite printed extraction result closes none of the unjoined,
component, ledger or physical-model OPEN items. (Retrieve M07, M09, M16.)

Two acceptance qualifications travel with all claims: SymPy obtains regional objects by substituting
regional fields, so those objects are **not a test of its support encoding**; and dimensional
homogeneity is definitional for the inferred order-flux units. The acceptance reports reviewed
independent constructions/controls and reproducibility; it is not this record's independent run or
new comparison. The comparator acceptance states that no leg found another pairing of different
objects, not that every leaf was compared. (Retrieve M02, M14.)

## 5. Separate sourced inventories

These are pointers to **each source row**, not replacements for its restrictions. They are not
computed objects, imposed profiles or completed closures. Retrieve M11 for the literal rows.

| Inventory | Source rows, kept with their original scope |
|---|---|
| Source/controller functions | `V3_STEP_PLAN.md:610–611`; `../pde_ledger_v2/notes/stage045_nonvariational_block_prep.md:45` and `:50`, with the return-controller row at `:45` read alongside its re-hosting context. The full local source/controller inventory is also `directives/O2_SHARED_PHYSICS.md`, §5 inventory table. All functional forms remain OPEN. |
| Boundary/domain data | `V3_STEP_PLAN.md:612`; `../../docs/two_throat_simulation_handoff_spec.md:340`, the **entire IR/return row**, including zero modes/moduli and moving/radiative prescriptions; stage 045's boundary row. U3 bears on the bulk IR condition only. No galactic/expansion boundary is imposed. |
| Restricted S9b/O2 domain | `steps/O2_steady_brane_balance.md:74–76`, the complete setting sentence. Its isolated spherical mass at rest with flowing drain, far-field steady profiles, far-field brane-rest lab time, linear optics and leading eikonal are **restrictions**, not S12 results. |
| Momentum/energy partners | `directives/O2_premise_decision_list.md:32–35`; `directives/O2_SHARED_PHYSICS.md:134`, `:243–248`, `:280`, `:584`; `steps/O2_steady_brane_balance.md:396`. Those rows preserve additional partners/reactions, the OPEN bulk-direction reduction and its normal-response/map operands, and conversion/carried-energy partners with reaction/supply systems. Local transported momentum fixes no carried total energy. |
| Native identifications | `steps/O2_steady_brane_balance.md:399`, O6 `𝒥_map`; `directives/O2_SHARED_PHYSICS.md:605`, S12 interface. These are the distinct native exchange/projection/measure and source/boundary handoffs, not the S14a bridge. |
| Separate power obligations | `SUBSTRATE_REQUIREMENTS.md`, `R-S12-01` and `R-O2-01` (retrieve M12). The former's S12 target is a historical register inference, live owner **no owner named**; the latter's target is **no owner named**. Neither whole obligation is assigned to S12. U3's material supply is not either power budget. |

**Return controllers: OPEN.** Stage 045 calls `I_ret/P_ret` projectors enforcing exact global
mass/momentum/energy closure **"despite the local sink"**. Its surrounding text says to re-host
the machinery on dynamical order conversion, while the committed law supplies local conservation
of total material. Those are the two sourced statements; the old mass-sink purpose cannot be
imported as an operative requirement. Whether these controllers have a remaining role, and which
global order/momentum/energy closure would call for it, needs a new assessment beyond retrieval.
It remains with the S12 source/controller inventory (`R-S12-02`); **no separate owner is named for
that assessment**. No controller is chosen or discarded here. (Retrieve M01, M11, M12.)

## 6. What remains OPEN, with its source's owner

| OPEN item | Owner / source limit |
|---|---|
| A13 material branch; phase-carrying content absent from this real-fraction object | S1 gate, propagated through S5/S12 (O2 record's branch row, M11); current decision list retains A13. |
| U1's re-ordering admissibility clause | No owner named; not encoded, with no rate or local sign chosen (M01). |
| Concentrated input fields, including order carried along a face | No owner named; explicitly left out for now, 2026-10-20 (M01). |
| Interior-order law/history/strain; faces' law; form of `J_χ` | No separate owner named for the first two; source/flux forms remain the S12 inventory. Slab material completion is separately S8/S22 (M01, M11). |
| "the projected order-loss source · the definition of `J` · the profile-dependent `J → Q` map · the controlled `ω → 0` return law" | S14a, plan `:636–637` (M01, M11). |
| Drain-strength/mass link | S13 (M01); no `GM` link delivered. |
| Flow profile | Each force states its requirement; S21 compares them (M01). A mismatch reopens S9b. |
| Slab's material law; puncture interior | S8/S22 for the former; S22 for the excluded interior (M01). |
| Controller role; separate boundary laws/data; native O6 identifications; conversion momentum/energy partners | S12 inventory and the normal-response/map/holder owners in the individual O2 rows; no completed native reduction (M11). |
| Non-passive reservoir and budget; net-power supplier and budget | Separate existing `R-S12-01` and `R-O2-01`, live owner no owner named (M12). |
| Unrestricted coefficient differences; unpaired/nested ledger/component coverage | No owner named; every tag and leaf remains in M06/M08/M09. |
| Prior-art applicability and historical objections for each candidate named below | No owner named for a new S12 check; none performed here (M13). |

Questions requiring computation beyond retrieval are therefore left open: controller re-hosting,
any cross-engine extension to unpaired components or nested/term ledgers, adjudication of the six
unrestricted coefficient differences, enforcement of U1's re-ordering clause, constitutive source/
face/flux laws, native reductions, and the already assigned S13/S14a/S21/S22 objects. They are not
prerequisites solved inside this record.

## 7. Register pass, both routes, and source disagreements

The register's schema and rest-on test govern. An adopted value/form is not a derived requirement;
an OPEN future convenience alone generates no entry. The S12 pass merges the printed source and
material-face compatibility objects into `R-S12-02`, whose native conversion/source and separate
boundary-data obligation already owns that content. No means of satisfying those requirements is
selected. The global premises, omitted concentrated fields and unsolved constitutive/mechanical
questions do not acquire a derived status. (Retrieve M01, M12, M15.)

**Census:** 18 entries before (`560c6c51`) and 18 after, all OPEN. **New entries by target: none.**
Four existing entries gain S12 source/provenance: `R-S1-01` (coordinate dimension only), `R-S1-02`
(conditional shear-region identification only), `R-S12-02` (printed compatibility requirements),
and `R-S21-01` (framing decision only, no new optical evidence). `R-S12-01` and `R-O2-01` remain
distinct carried conditions without a computed S12 power consumer; neither is delivered or duplicated.
The full per-entry rest-on census is the S12 register pass (M12, M15).

**Second route:** the register requires it for every sector with identified prior art, including
objects outside a step record's own doubts. Accordingly the pass reads the inherited O2 scoping §5
and the S9b register inventory (M13). Surface/exchange/reference candidates are Gurtin–Murdoch,
Abels–Garcke–Grün, Scriven, Helfrich and Rodriguez–Hoger–McCulloch. Flow candidates include
Hamilton–Lisle, Visser and Unruh, and the inherited Cahill/Sbitnev inflow class. None supplies an
applied S12 oracle or reproduced S12 result: supplied local balances and the support identification
are not an independently established match to a candidate's domain. For **each** candidate the
questions of applicability, transferable historical failures/objections and their answer here remain
OPEN, **no owner named**. The source access limitations remain those reported in scoping §5; no
fresh literature verification is claimed. The other inherited S9b candidates, each individually,
retain their existing OPEN handoffs and no S12 result is attributed to them. No new entry rests on
an unapplied oracle. The four inherited river/acoustic objections remain: driving flow and material
destination with S12/O2; scalar sound versus shear light as S9b's scoped exclusion; kinematics versus
dynamics in the analog-gravity handoff. This source identity is not an answer to them (M13, M15).

**Source disagreements, preserved:**

- Stage 045's controller purpose **"despite the local sink"** and global
  **"`∫(S_drain + S_leakage) = ∫f = ∫w = 0`"** are historical sink language. The current
  decision list instead supplies **"`∂_t n + ∇₄·(n u) = 0` (total material)"** and says
  **"The frozen-wall mass sink with a remote return stays ruled out."** The source also calls for
  re-hosting. The unresolved role question is retained, not answered with the old sink (M01, M11).
- The register's `R-S12-02` asks for **"the native drain/return functions and their separate boundary
  data"**; the 2026-10-19 user framing says S12 **"states what any conversion law must satisfy and
  leaves the forms open"**. The obligation remains OPEN, with the newly printed compatibility
  constraints merged into it; S12 record completion is not delivery of those functions (M01, M12, M15).
- Historical second-route disagreement: the earlier directive says **"Run it only on a prior result
  that a record in the reading list identifies as reproduced by the sector"**; the register says
  **"Run this route for every sector with identified prior art"** and records that its MacCullagh
  objections were absent from the original records. The register wins here; the earlier disagreement
  remains unresolved, quoted in M13. No new operative disagreement with this record contract is found.
- The accepted comparator intended counterpart components; its engine row pairs library-version text
  with a platform identifier. Both sides and the user's acceptance qualification remain in §4 and M02/M08.

The status-only edits record light-closeout acceptance at `e070c7f3` in the register, the S9b oracle
record and repository `STATUS.md`, including the former wait-on-review line. `STATUS.md` states S12's
accepted constructions/comparator and this record's own review-pending standing. No other scientific
record is upgraded. Commits made for this deliverable: **none**. (Retrieve M15.)
