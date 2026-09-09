The directive is not clear to build. No explicit whole-carrier/μ/source equality remains, but the dictionary is still not fully frozen, two admitted rules lack valid production scope, and the mandatory controls are incomplete.

## MUST findings

1. **The energy-basis map still delegates physics-bearing normalization choices to the builder.**

The directive requests an injective coefficient-name table “by contraction identity,” but provides neither the table nor its normalization/direction ([directive:81](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md:81)). Those choices are not mechanical spelling:

- WL assigns `bRho/2` to `WBg·theta²`, `cCoupling` to `WBg·theta·eW`, and otherwise creates `energyCoefficient<i>` in retained-basis order ([WL:211](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:211), [WL:222](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:222)).
- SymPy assigns, for example, `B_rho_3·W_bg/(2W_0)`, `C·W_bg`, and `kappa_theta/2` to its corresponding abstract invariants ([brane:1584](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:1584)); its first-jet coefficients are indexed separately as `gamma_s11cb_*` ([brane:357](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:357), [brane:1820](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:1820)).

Thus rules such as `bRho → B_rho_3/W_0` versus `bRho → B_rho_3`, or `energyCoefficient → kappa_theta/2` versus `kappa_theta`, change the residual. “Construct it from the sites” makes astra decide this during the build, contrary to freezing it before testing. The directive must supply the complete contraction-ID table, canonical orientation, and scale for every pair—or prescribe a deterministic, independently reviewed basis-change algorithm and its expected table—not merely an allowed map type.

2. **Map 6 is not a leftover-name map of the emitted operands and cannot satisfy a load-bearing control.**

SymPy replaces the imported c1 density atom by `inputs.density[(rho,)][1]` before forming source terms ([diagnostic:378](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:378), [diagnostic:382](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:382)). WL likewise substitutes its already-computed `density3` through `rhoFace → density` inside `sourceBind` ([WL:380](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:380)), with `density3=density4·WBg` computed before the source calls ([WL:839](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:839), [WL:869](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:869)).

Consequently, `inputs.density[(rho,)][1]` and `density3` are expressions/local variables, not surviving emitted atom names. Map 6 therefore has two possible implementations:

- inert, in which case its corruption/drop controls cannot bite; or
- a new constant↔live-expression rewrite, reopening exactly the c1 rule-17 hazard the directive says is forbidden.

Remove it from the production rewrite dictionary. The live-density construction can remain a separately emitted premise/census comparison.

3. **The dictionary needs typed scopes; map 2 is presently over-broad.**

The source-wave rule is independently justified only at the source extraction boundary: SymPy inserts applied wave jets at `Y` ([diagnostic:392](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:392), [selfenergy:139](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:139)), whereas WL’s `sourceMap` retains bare registered jets ([WL:707](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:707)).

The directive nevertheless says every rule is applied unchanged to carriers and Φ as well as sources ([directive:63](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md:63)). Φ is an abstract jet map with no `Y` evaluation in either engine ([covariance:63](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:63), [WL:236](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:236)). Applying the source rule there changes the claim from equality of abstract jet maps to equality after point evaluation.

Freeze each rule with a typed family/stage scope. In particular, the bare↔applied-`Y` map belongs only to `SOURCE_*`; Φ receives domain/multi-index spelling and pre-grade profile substitution, not source-point evaluation.

4. **The mandatory controls are incomplete and cannot establish locality or dropped-map sensitivity.**

The directive’s four controls are actually: map corruption, reachability census, blanket collapse, and grade tripwire ([directive:137](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md:137)). It has no dropped-primitive control.

Its corruption control requires a nonzero delta on using leaves, but does not require zero delta on every non-using leaf. A mistakenly broad rule can therefore pass. It also explicitly permits admitted maps with no occurrence ([directive:143](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md:143)), making their corruption controls unable to bite.

Required repair:

- corrupt each active map on one side;
- require a nonzero control delta on at least one declared-use leaf;
- require zero control delta outside that map’s declared use-set;
- delete each active primitive and require a residual on a predeclared defining-relation control operand;
- run blanket collapse through the same extraction/rewrite/witness path;
- retain the grade-combining tripwire.

Synthetic defining-relation operands can make these controls able to fail without leaking the shipped collapse pattern.

5. **The census crosswalk remains illustrative rather than frozen, and one example invents a SymPy counterpart.**

The directive says “e.g.” before several aliases ([directive:154](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md:154)). That leaves the builder to decide which fields are equivalent.

In particular, SymPy emits `kappa_a`, `kappa_j`, baseline parameters, material tag, junk symbol/dimension/wave, and inserted amplitude—no material-normal field ([covariance:146](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:146)). WL separately emits `ADVECTION`, `THICKNESS`, `JUNK`, `JUNK_CASE`, `JUNK_SYMBOL`, `JUNK_DIMENSIONS`, `JUNK_ASSUMPTION`, and `MATERIAL_NORMAL` ([WL:977](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:977)). Therefore `MATERIAL_NORMAL` is one-sided, like `THICKNESS`; it cannot be crosswalked to an implicit SymPy zero without changing the census claim.

The directive must enumerate the exact census crosswalk and derived comparisons, while retaining all other fields as one-sided. The domain-census mappings should likewise explicitly pair `DOMAIN`, `COVERAGE`, `UNCOVERED`, and `MAX_RANK` with their SymPy counterparts.

## Eight-map audit

| Map | Assessment |
|---|---|
| 1. Grad-θ jet duals | Primitive and grounded: the engines use `theta_d{i}` versus `grad_theta_{i}`, with `wave_jet` supplying the dual spelling ([geometry:140](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py:140), [brane:227](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:227), [WL:945](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:945)). |
| 2. Source-wave point map | Primitive for source families only; over-broad if applied to Φ or carriers. |
| 3. Profile-jet names/scales | Primitive. WL and SymPy independently encode the same σ-scaled background jets ([WL:127](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:127), [brane:767](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:767)). Post-grade η re-expansion is correctly forbidden. |
| 4. Energy-basis pairing | Legitimate in principle, but not frozen because the exact normalization table is missing. |
| 5. Φ domain/multi-index spelling | Primitive; the directive correctly retains Φ values as the tested operands. |
| 6. Leftover density name | No surviving emitted-name domain; remove from production rewrites. |
| 7. Jacobian spelling | Primitive only if the literal factor survives. Degree-2 projection is correctly excluded: it is already performed inside both material-amplitude constructions ([brane:1970](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:1970), [WL:244](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:244)). |
| 8. Conditional c1 atoms | Safe only as occurrence-conditioned exact rules. In these scoped source constructions, default Fourier/on-shell maps are correctly excluded because source extraction precedes the kernel ([diagnostic:378](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:378), [WL:364](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:364)). |

## Checks that pass

- **Retained order:** sound. SymPy declares and preserves four independent grades ([diagnostic:43](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:43), [diagnostic:177](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:177), [diagnostic:245](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:245)); WL does likewise ([WL:123](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:123), [WL:439](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:439)). No σ binding or grade combination is permitted.

- **Witness:** sound in concept. Residual precedes the three-valued witness; PIT primes/seeds are disjoint from both engines; denominator rejection, derived bounds, attempt/time/RSS limits, and exhaustion→`undecided` are required ([directive:125](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md:125)). No verdict or residual-zero exit is authorized.

- **BASELINE:** correctly a nominal duplicate, not a third discriminator. SymPy’s actual and baseline differ only through shipped `kappa_a=1`, `kappa_j=0` controls ([covariance:137](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:137)); WL constructs the corresponding baseline separately at the same settings ([WL:845](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:845)).

- **Value-free discipline:** honored. The comparator-run counts are identified as historical facts rather than collapse targets; no shipped witness count or pass threshold appears.

- **Extraction and fence:** importing the comparator’s loaders, spelling maps, typed decoders, numeric extraction, `extract_meta`, and `materialize` is sufficient without editing it ([comparator:631](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:631), [comparator:688](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:688), [comparator:862](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:862)). The build→verify→run→report→STOP fence is explicit ([directive:179](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md:179)).

**DIRECTIVE NOT-SOUND — MUST: freeze the exact normalized energy-coefficient table and typed map scopes; remove the nonexistent leftover-density rewrite; add dropped-map and corruption-locality controls; and enumerate the complete census crosswalk without synthesizing a SymPy material-normal field.**
The directive is not clear to build. No explicit whole-carrier/μ/source equality remains, but the dictionary is still not fully frozen, two admitted rules lack valid production scope, and the mandatory controls are incomplete.

## MUST findings

1. **The energy-basis map still delegates physics-bearing normalization choices to the builder.**

The directive requests an injective coefficient-name table “by contraction identity,” but provides neither the table nor its normalization/direction ([directive:81](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md:81)). Those choices are not mechanical spelling:

- WL assigns `bRho/2` to `WBg·theta²`, `cCoupling` to `WBg·theta·eW`, and otherwise creates `energyCoefficient<i>` in retained-basis order ([WL:211](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:211), [WL:222](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:222)).
- SymPy assigns, for example, `B_rho_3·W_bg/(2W_0)`, `C·W_bg`, and `kappa_theta/2` to its corresponding abstract invariants ([brane:1584](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:1584)); its first-jet coefficients are indexed separately as `gamma_s11cb_*` ([brane:357](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:357), [brane:1820](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:1820)).

Thus rules such as `bRho → B_rho_3/W_0` versus `bRho → B_rho_3`, or `energyCoefficient → kappa_theta/2` versus `kappa_theta`, change the residual. “Construct it from the sites” makes astra decide this during the build, contrary to freezing it before testing. The directive must supply the complete contraction-ID table, canonical orientation, and scale for every pair—or prescribe a deterministic, independently reviewed basis-change algorithm and its expected table—not merely an allowed map type.

2. **Map 6 is not a leftover-name map of the emitted operands and cannot satisfy a load-bearing control.**

SymPy replaces the imported c1 density atom by `inputs.density[(rho,)][1]` before forming source terms ([diagnostic:378](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:378), [diagnostic:382](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:382)). WL likewise substitutes its already-computed `density3` through `rhoFace → density` inside `sourceBind` ([WL:380](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:380)), with `density3=density4·WBg` computed before the source calls ([WL:839](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:839), [WL:869](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:869)).

Consequently, `inputs.density[(rho,)][1]` and `density3` are expressions/local variables, not surviving emitted atom names. Map 6 therefore has two possible implementations:

- inert, in which case its corruption/drop controls cannot bite; or
- a new constant↔live-expression rewrite, reopening exactly the c1 rule-17 hazard the directive says is forbidden.

Remove it from the production rewrite dictionary. The live-density construction can remain a separately emitted premise/census comparison.

3. **The dictionary needs typed scopes; map 2 is presently over-broad.**

The source-wave rule is independently justified only at the source extraction boundary: SymPy inserts applied wave jets at `Y` ([diagnostic:392](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:392), [selfenergy:139](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:139)), whereas WL’s `sourceMap` retains bare registered jets ([WL:707](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:707)).

The directive nevertheless says every rule is applied unchanged to carriers and Φ as well as sources ([directive:63](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md:63)). Φ is an abstract jet map with no `Y` evaluation in either engine ([covariance:63](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:63), [WL:236](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:236)). Applying the source rule there changes the claim from equality of abstract jet maps to equality after point evaluation.

Freeze each rule with a typed family/stage scope. In particular, the bare↔applied-`Y` map belongs only to `SOURCE_*`; Φ receives domain/multi-index spelling and pre-grade profile substitution, not source-point evaluation.

4. **The mandatory controls are incomplete and cannot establish locality or dropped-map sensitivity.**

The directive’s four controls are actually: map corruption, reachability census, blanket collapse, and grade tripwire ([directive:137](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md:137)). It has no dropped-primitive control.

Its corruption control requires a nonzero delta on using leaves, but does not require zero delta on every non-using leaf. A mistakenly broad rule can therefore pass. It also explicitly permits admitted maps with no occurrence ([directive:143](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md:143)), making their corruption controls unable to bite.

Required repair:

- corrupt each active map on one side;
- require a nonzero control delta on at least one declared-use leaf;
- require zero control delta outside that map’s declared use-set;
- delete each active primitive and require a residual on a predeclared defining-relation control operand;
- run blanket collapse through the same extraction/rewrite/witness path;
- retain the grade-combining tripwire.

Synthetic defining-relation operands can make these controls able to fail without leaking the shipped collapse pattern.

5. **The census crosswalk remains illustrative rather than frozen, and one example invents a SymPy counterpart.**

The directive says “e.g.” before several aliases ([directive:154](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md:154)). That leaves the builder to decide which fields are equivalent.

In particular, SymPy emits `kappa_a`, `kappa_j`, baseline parameters, material tag, junk symbol/dimension/wave, and inserted amplitude—no material-normal field ([covariance:146](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:146)). WL separately emits `ADVECTION`, `THICKNESS`, `JUNK`, `JUNK_CASE`, `JUNK_SYMBOL`, `JUNK_DIMENSIONS`, `JUNK_ASSUMPTION`, and `MATERIAL_NORMAL` ([WL:977](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:977)). Therefore `MATERIAL_NORMAL` is one-sided, like `THICKNESS`; it cannot be crosswalked to an implicit SymPy zero without changing the census claim.

The directive must enumerate the exact census crosswalk and derived comparisons, while retaining all other fields as one-sided. The domain-census mappings should likewise explicitly pair `DOMAIN`, `COVERAGE`, `UNCOVERED`, and `MAX_RANK` with their SymPy counterparts.

## Eight-map audit

| Map | Assessment |
|---|---|
| 1. Grad-θ jet duals | Primitive and grounded: the engines use `theta_d{i}` versus `grad_theta_{i}`, with `wave_jet` supplying the dual spelling ([geometry:140](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py:140), [brane:227](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:227), [WL:945](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:945)). |
| 2. Source-wave point map | Primitive for source families only; over-broad if applied to Φ or carriers. |
| 3. Profile-jet names/scales | Primitive. WL and SymPy independently encode the same σ-scaled background jets ([WL:127](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:127), [brane:767](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:767)). Post-grade η re-expansion is correctly forbidden. |
| 4. Energy-basis pairing | Legitimate in principle, but not frozen because the exact normalization table is missing. |
| 5. Φ domain/multi-index spelling | Primitive; the directive correctly retains Φ values as the tested operands. |
| 6. Leftover density name | No surviving emitted-name domain; remove from production rewrites. |
| 7. Jacobian spelling | Primitive only if the literal factor survives. Degree-2 projection is correctly excluded: it is already performed inside both material-amplitude constructions ([brane:1970](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:1970), [WL:244](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:244)). |
| 8. Conditional c1 atoms | Safe only as occurrence-conditioned exact rules. In these scoped source constructions, default Fourier/on-shell maps are correctly excluded because source extraction precedes the kernel ([diagnostic:378](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:378), [WL:364](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:364)). |

## Checks that pass

- **Retained order:** sound. SymPy declares and preserves four independent grades ([diagnostic:43](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:43), [diagnostic:177](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:177), [diagnostic:245](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:245)); WL does likewise ([WL:123](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:123), [WL:439](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:439)). No σ binding or grade combination is permitted.

- **Witness:** sound in concept. Residual precedes the three-valued witness; PIT primes/seeds are disjoint from both engines; denominator rejection, derived bounds, attempt/time/RSS limits, and exhaustion→`undecided` are required ([directive:125](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md:125)). No verdict or residual-zero exit is authorized.

- **BASELINE:** correctly a nominal duplicate, not a third discriminator. SymPy’s actual and baseline differ only through shipped `kappa_a=1`, `kappa_j=0` controls ([covariance:137](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:137)); WL constructs the corresponding baseline separately at the same settings ([WL:845](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:845)).

- **Value-free discipline:** honored. The comparator-run counts are identified as historical facts rather than collapse targets; no shipped witness count or pass threshold appears.

- **Extraction and fence:** importing the comparator’s loaders, spelling maps, typed decoders, numeric extraction, `extract_meta`, and `materialize` is sufficient without editing it ([comparator:631](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:631), [comparator:688](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:688), [comparator:862](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:862)). The build→verify→run→report→STOP fence is explicit ([directive:179](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md:179)).

**DIRECTIVE NOT-SOUND — MUST: freeze the exact normalized energy-coefficient table and typed map scopes; remove the nonexistent leftover-density rewrite; add dropped-map and corruption-locality controls; and enumerate the complete census crosswalk without synthesizing a SymPy material-normal field.**
