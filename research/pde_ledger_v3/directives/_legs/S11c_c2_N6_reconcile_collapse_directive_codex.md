The directive is not clear to build. It contains several computation-changing ambiguities, including one rule-5/blanket-collapse hazard.

## MUST findings

1. **The “frozen dictionary” is not an executable frozen dictionary, and map 6 can assume part of the source question.**

   Entries 2, 3, 5, 6, 7, and 8 are descriptions of construction stages or families of possible rules, not explicit domain→image substitutions with direction, scope, and normalization. Most seriously, “face-velocity + source-solve factors” permits equating independently computed subobjects:

   - SymPy independently builds material velocity and extracts source factors from the imported response at [reconcile:82–92](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_sympy.py:82) and [diagnostic:378–389](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:378).
   - WL independently derives its algebraic source solve and geometric velocities at [WL:366–381](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:366) and [WL:855–872](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:855).
   - The adopted resolution explicitly keeps physical face-velocity correctness open at [RESOLVED:41–43](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_RESOLVED.md:41).

   Therefore map 6 is not primitive as written. Mapping whole velocity or source-solve factors would suppress a component of `SOURCE_*` by assumption. It must be replaced by explicit symbol, jet, parameter, and normalization rules; computed velocity/source factors must remain in the residual.

   Map 5 likewise needs a predeclared injective coefficient table derived solely from contraction IDs—not a run-time match or solve against μ/source residuals. WL assigns generated coefficients to its quotient basis at [WL:211–228](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:211), while SymPy selects and coefficients a separately generated basis at [brane engine:1727–1872](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:1727). The directive supplies no pair table, scale factors, or injectivity/completeness check.

2. **A needed primitive source-wave representation map is missing.**

   SymPy’s source value inserts an applied wave jet at `Y` ([diagnostic:392–396](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:392)); `wave_jet` constructs applied fields and derivatives ([selfenergy engine:139–169](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:139)). WL’s `sourceMap` instead retains registered bare jet atoms ([WL:707–711](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:707)).

   The reused comparator deliberately does **not** identify bare fields with applied fields and preserves applied heads/arguments ([comparator:287–289](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:287), [comparator:990–993](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:990)). Map 1 only bridges `a.grad_theta[i]` and `b.grad_theta[i]`; it does not bridge WL’s bare source jets to SymPy’s point-labelled applied jets. An exact, point-preserving map such as `jet[f,I] ↔ ∂I s11cc2Field_f(Y,t)` is required.

3. **The coefficientwise contract is not implemented for `FROZEN_PHI`.**

   Carrier and source numeric leaves really are keyed by independent `ETA,SIGMA` axes ([comparator:107–113](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:107), [comparator:143–169](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:143)). But `FROZEN_PHI` is metadata, and `extract_meta` produces `MAP_VARIABLE`/`FIELD_PATH` keys without grade axes ([comparator:59–60](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:59), [comparator:688–765](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:688)).

   The engines emit full ungraded Φ rules at [covariance:116–124](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:116) and [WL:975–976](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:975). Thus the directive’s claim that “operands are already graded by the comparator key” is false for Φ. It must prescribe profile expansion followed by explicit extraction of all four coefficients before the witness.

4. **The four controls do not all bite as specified.**

   - The blanket control is tautological if implemented by setting one operand equal to the other. It cannot prove the normal path lacks such a replacement.
   - The one-sided corruption control only requires the resulting residual to be nonzero. A leaf already nonzero before corruption satisfies that condition. It must require a nonzero **control delta** between baseline and corrupted bridged residuals.
   - Maps 2, 3, and 8 may have no surviving production-leaf occurrence because those operations were already performed inside the engines. They therefore cannot each satisfy “corruption moves a leaf.”
   - The dropped-map control is conditional on a production leaf collapsing, so it can vacuously pass.
   - “Every leaf collapsed” is incompatible with legitimate parse/budget `undecided` leaves unless the blanket control bypasses extraction, which makes it still more tautological.

   Each map needs an independent defining-relation control plus a separate production reachability census. The grade tripwire is otherwise sound once Φ is actually graded.

5. **The census audit lacks the semantic schema needed to compute its claimed equivalences.**

   SymPy emits `kappa_a`, `kappa_j`, `baseline_parameters`, etc. ([covariance:146–153](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:146)); WL emits `ADVECTION`, `THICKNESS`, `JUNK`, `MATERIAL_NORMAL`, etc. ([WL:977–980](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:977)). The comparator lowercases field names but explicitly adds no aliases such as `kappa_a=ADVECTION` ([comparator:683–685](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:683)).

   The same issue affects `max_present_jet_rank` versus `MAX_RANK`. Raw “matched keys” therefore cannot perform the promised production-control/max-rank equivalence. The directive must freeze a census-only semantic crosswalk while retaining unmatched and one-sided fields such as WL `THICKNESS`.

6. **The prescribed extraction path does not expose operand pairs before residual construction.**

   `compare_family` materializes operands, computes `difference`, then emits `CASE` ([comparator:794–825](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:794)); it subsequently releases them. `object_work` returns accounting, not pairs ([comparator:862–898](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:862)). Consequently, the builder cannot simultaneously:

   - reuse that path verbatim,
   - obtain CAS operands before any residual,
   - avoid reimplementing pairing, and
   - obey the no-outside-edit fence.

   A comparator API returning matched `Leaf` pairs is needed, or the directive must explicitly authorize and specify the adapter.

## Eight-map audit

| Map | Assessment |
|---|---|
| 1. Jet vocabulary | Primitive and grounded at [reconcile:46–69](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_sympy.py:46) and [WL:939–948](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:939), but incomplete for source applied-field jets. |
| 2. Jacobian/degree two | Genuine shared construction identity: SymPy applies `1+tr(∇u)` and extracts wave degree two at [brane engine:1970–1981](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:1970); WL does so at [WL:244–247](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:244). It is not normally a surviving operand rewrite and belongs in construction/census validation. |
| 3. Live density | Primitive only with explicit case and point formulas. The cited SymPy line is wrong: the substitution is at [diagnostic:382–389](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:382), and WL defines `density4,density3` at [WL:839–840](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:839). The c1 precedent made c2 re-adjudication mandatory; it did not authorize an unconditional constant↔live-field equality ([c1 precedent:137–163](/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_reconcile.md:137)). |
| 4. Background/profile jets | Primitive and needed, but must enumerate all retained jet ranks and X/Y locations. Exact formulas are at [WL:127–132](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:127) and [brane engine:767–821](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:767). The directive’s SymPy `:204,281` citation is incomplete. |
| 5. Energy basis | Potentially primitive only as a frozen contraction-ID coefficient table. No such table or normalization is supplied; current wording is unsafe. |
| 6. Velocity/source factors | Not primitive as written; it can equate independently computed subobjects and partially assume the source answer. MUST replace. |
| 7. Φ spelling | Primitive only for domain names, sorted multi-indices, and derivative syntax. Φ values are independently constructed at [covariance:63–129](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:63) and [WL:236–247](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:236); value equality must remain the residual. |
| 8. c1 maps | Not one frozen map and has no two-engine citations. Name folding is already mechanical; Fourier/on-shell rules are absent from the scoped source construction and conflict with the default-Fourier exclusion. Split into exact conditional rules or omit inactive rules. |

## Parts that are sound

The independent retained rectangle is real in both engines ([diagnostic:43,245–248](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:43), [WL:123–141](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:123)); the σ-binding, cross-anchoring, and default on-shell exclusions are correct. `SOURCE_BASELINE` is correctly classified as the shipped nominal duplicate ([covariance:137–154](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:137), [WL:845–849](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:845)). The three witness values, residual-before-witness order, no-verdict clauses, value-free discipline, and builder fence are conceptually correct.

The PIT clause still needs enforceable bounds: disjoint fresh primes/seeds, denominator rejection, a derived degree/exclusion bound, maximum attempts/time/RSS, and `undecided` on exhaustion. Emitting only primes, draws, and a claimed δ is insufficient to validate the probabilistic witness.

**DIRECTIVE NOT-SOUND — MUST: freeze explicit primitive rewrite tables; remove/replace the velocity/source-factor blanket-collapse path; add the point-preserving source-wave map and coefficientwise Φ extraction; define census schema mappings; repair the controls, PIT bounds, and operand-pair extraction API.**
tokens used
The directive is not clear to build. It contains several computation-changing ambiguities, including one rule-5/blanket-collapse hazard.

## MUST findings

1. **The “frozen dictionary” is not an executable frozen dictionary, and map 6 can assume part of the source question.**

   Entries 2, 3, 5, 6, 7, and 8 are descriptions of construction stages or families of possible rules, not explicit domain→image substitutions with direction, scope, and normalization. Most seriously, “face-velocity + source-solve factors” permits equating independently computed subobjects:

   - SymPy independently builds material velocity and extracts source factors from the imported response at [reconcile:82–92](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_sympy.py:82) and [diagnostic:378–389](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:378).
   - WL independently derives its algebraic source solve and geometric velocities at [WL:366–381](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:366) and [WL:855–872](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:855).
   - The adopted resolution explicitly keeps physical face-velocity correctness open at [RESOLVED:41–43](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_RESOLVED.md:41).

   Therefore map 6 is not primitive as written. Mapping whole velocity or source-solve factors would suppress a component of `SOURCE_*` by assumption. It must be replaced by explicit symbol, jet, parameter, and normalization rules; computed velocity/source factors must remain in the residual.

   Map 5 likewise needs a predeclared injective coefficient table derived solely from contraction IDs—not a run-time match or solve against μ/source residuals. WL assigns generated coefficients to its quotient basis at [WL:211–228](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:211), while SymPy selects and coefficients a separately generated basis at [brane engine:1727–1872](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:1727). The directive supplies no pair table, scale factors, or injectivity/completeness check.

2. **A needed primitive source-wave representation map is missing.**

   SymPy’s source value inserts an applied wave jet at `Y` ([diagnostic:392–396](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:392)); `wave_jet` constructs applied fields and derivatives ([selfenergy engine:139–169](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:139)). WL’s `sourceMap` instead retains registered bare jet atoms ([WL:707–711](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:707)).

   The reused comparator deliberately does **not** identify bare fields with applied fields and preserves applied heads/arguments ([comparator:287–289](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:287), [comparator:990–993](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:990)). Map 1 only bridges `a.grad_theta[i]` and `b.grad_theta[i]`; it does not bridge WL’s bare source jets to SymPy’s point-labelled applied jets. An exact, point-preserving map such as `jet[f,I] ↔ ∂I s11cc2Field_f(Y,t)` is required.

3. **The coefficientwise contract is not implemented for `FROZEN_PHI`.**

   Carrier and source numeric leaves really are keyed by independent `ETA,SIGMA` axes ([comparator:107–113](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:107), [comparator:143–169](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:143)). But `FROZEN_PHI` is metadata, and `extract_meta` produces `MAP_VARIABLE`/`FIELD_PATH` keys without grade axes ([comparator:59–60](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:59), [comparator:688–765](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:688)).

   The engines emit full ungraded Φ rules at [covariance:116–124](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:116) and [WL:975–976](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:975). Thus the directive’s claim that “operands are already graded by the comparator key” is false for Φ. It must prescribe profile expansion followed by explicit extraction of all four coefficients before the witness.

4. **The four controls do not all bite as specified.**

   - The blanket control is tautological if implemented by setting one operand equal to the other. It cannot prove the normal path lacks such a replacement.
   - The one-sided corruption control only requires the resulting residual to be nonzero. A leaf already nonzero before corruption satisfies that condition. It must require a nonzero **control delta** between baseline and corrupted bridged residuals.
   - Maps 2, 3, and 8 may have no surviving production-leaf occurrence because those operations were already performed inside the engines. They therefore cannot each satisfy “corruption moves a leaf.”
   - The dropped-map control is conditional on a production leaf collapsing, so it can vacuously pass.
   - “Every leaf collapsed” is incompatible with legitimate parse/budget `undecided` leaves unless the blanket control bypasses extraction, which makes it still more tautological.

   Each map needs an independent defining-relation control plus a separate production reachability census. The grade tripwire is otherwise sound once Φ is actually graded.

5. **The census audit lacks the semantic schema needed to compute its claimed equivalences.**

   SymPy emits `kappa_a`, `kappa_j`, `baseline_parameters`, etc. ([covariance:146–153](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:146)); WL emits `ADVECTION`, `THICKNESS`, `JUNK`, `MATERIAL_NORMAL`, etc. ([WL:977–980](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:977)). The comparator lowercases field names but explicitly adds no aliases such as `kappa_a=ADVECTION` ([comparator:683–685](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:683)).

   The same issue affects `max_present_jet_rank` versus `MAX_RANK`. Raw “matched keys” therefore cannot perform the promised production-control/max-rank equivalence. The directive must freeze a census-only semantic crosswalk while retaining unmatched and one-sided fields such as WL `THICKNESS`.

6. **The prescribed extraction path does not expose operand pairs before residual construction.**

   `compare_family` materializes operands, computes `difference`, then emits `CASE` ([comparator:794–825](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:794)); it subsequently releases them. `object_work` returns accounting, not pairs ([comparator:862–898](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:862)). Consequently, the builder cannot simultaneously:

   - reuse that path verbatim,
   - obtain CAS operands before any residual,
   - avoid reimplementing pairing, and
   - obey the no-outside-edit fence.

   A comparator API returning matched `Leaf` pairs is needed, or the directive must explicitly authorize and specify the adapter.

## Eight-map audit

| Map | Assessment |
|---|---|
| 1. Jet vocabulary | Primitive and grounded at [reconcile:46–69](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_sympy.py:46) and [WL:939–948](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:939), but incomplete for source applied-field jets. |
| 2. Jacobian/degree two | Genuine shared construction identity: SymPy applies `1+tr(∇u)` and extracts wave degree two at [brane engine:1970–1981](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:1970); WL does so at [WL:244–247](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:244). It is not normally a surviving operand rewrite and belongs in construction/census validation. |
| 3. Live density | Primitive only with explicit case and point formulas. The cited SymPy line is wrong: the substitution is at [diagnostic:382–389](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:382), and WL defines `density4,density3` at [WL:839–840](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:839). The c1 precedent made c2 re-adjudication mandatory; it did not authorize an unconditional constant↔live-field equality ([c1 precedent:137–163](/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_reconcile.md:137)). |
| 4. Background/profile jets | Primitive and needed, but must enumerate all retained jet ranks and X/Y locations. Exact formulas are at [WL:127–132](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:127) and [brane engine:767–821](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py:767). The directive’s SymPy `:204,281` citation is incomplete. |
| 5. Energy basis | Potentially primitive only as a frozen contraction-ID coefficient table. No such table or normalization is supplied; current wording is unsafe. |
| 6. Velocity/source factors | Not primitive as written; it can equate independently computed subobjects and partially assume the source answer. MUST replace. |
| 7. Φ spelling | Primitive only for domain names, sorted multi-indices, and derivative syntax. Φ values are independently constructed at [covariance:63–129](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:63) and [WL:236–247](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:236); value equality must remain the residual. |
| 8. c1 maps | Not one frozen map and has no two-engine citations. Name folding is already mechanical; Fourier/on-shell rules are absent from the scoped source construction and conflict with the default-Fourier exclusion. Split into exact conditional rules or omit inactive rules. |

## Parts that are sound

The independent retained rectangle is real in both engines ([diagnostic:43,245–248](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:43), [WL:123–141](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:123)); the σ-binding, cross-anchoring, and default on-shell exclusions are correct. `SOURCE_BASELINE` is correctly classified as the shipped nominal duplicate ([covariance:137–154](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:137), [WL:845–849](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:845)). The three witness values, residual-before-witness order, no-verdict clauses, value-free discipline, and builder fence are conceptually correct.

The PIT clause still needs enforceable bounds: disjoint fresh primes/seeds, denominator rejection, a derived degree/exclusion bound, maximum attempts/time/RSS, and `undecided` on exhaustion. Emitting only primes, draws, and a claimed δ is insufficient to validate the probabilistic witness.

**DIRECTIVE NOT-SOUND — MUST: freeze explicit primitive rewrite tables; remove/replace the velocity/source-factor blanket-collapse path; add the point-preserving source-wave map and coefficientwise Φ extraction; define census schema mappings; repair the controls, PIT bounds, and operand-pair extraction API.**
