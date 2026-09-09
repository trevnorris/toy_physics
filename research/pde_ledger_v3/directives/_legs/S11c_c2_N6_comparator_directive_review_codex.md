he DAG is reconstructable, but the literal schema is incompletely specified.

Evidence: The four operations are `number`, `variable`, `add`, `mul`, and `pow` as constructed in the diagnostic source ([node constructors](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:115)); references are emitted topologically ([emitter](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_sympy.py:151)). But applied functions and derivatives are encoded as `(sp.srepr(e), "formal_jet")` variable literals ([compiler](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:237)); those strings occur in the committed DAG ([actual DAG](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_c2_N6_reconcile_sympy.out:45)).

Why wrong: A builder believing all literals are simple names may turn formal jets into opaque symbols and manufacture symbolic disagreements.

Minimal fix: Specify all literal variants—`global`, local point, `formal_jet`, `algebraic_i`, and integral-guard—and require a synthetic round-trip containing a real `formal_jet` srepr literal.

### 8. Rule 5 is explicitly breached, and the DoD can accept an all-deferred/all-mismatched instrument

Decision: The directive is not value-free and its coverage gate is non-failing.

Evidence: It states that `REP_INVARIANCE_RESIDUAL` is a “wrong-object — **nonzero is EXPECTED**” ([directive](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_comparator_build_directive.md:70)). The engine instead withholds every expected residual/control value ([engine header](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:12)).

The DoD allows every family to terminate as `axis_set_mismatch` or `deferred_oversize` and expressly refuses any positive-join gate ([DoD](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_comparator_build_directive.md:191)). Combined with the unspecified budget and missing bridge, zero actual comparisons can satisfy the prose.

Why wrong: The expected-nonzero statement contaminates later adjudication; the coverage rule cannot distinguish a cautious comparator from a nonfunctional extractor.

Minimal fix: Delete the quoted expectation. Require a source-derived synthetic fixture for every joined schema and treat zero extracted leaves from two nonempty, declared-shared containers as an operational extractor failure—not a physics disagreement and not a residual target.

## Decisions that are sound

The directive correctly preserves the inherited three-valued residual machinery and bans per-case verdict tokens. The precedent really imports `residual` and its four typed outcomes ([c1 imports](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c1_cross_engine_comparator.py:30)), provides the WL→SymPy canonicalizer ([translator](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c1_cross_engine_comparator.py:462)), and has lazy parse/release objects ([lazy cases](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c1_cross_engine_comparator.py:519)).

The named SymPy-only residuals are also correctly excluded from the cross-engine join, and the raw PIT and digest seals should remain. The directive becomes sound only after the eight corrections above, especially the completed join set and pre-registered schema bridge.
hook: Stop
hook: Stop Completed
tokens used
223,353
# Overall decision: NOT-SOUND

The raw PIT seal is correct, but the directive would still permit a comparator that misses four shared covariance families, fails to join the real schemas, mishandles valid repeated SymPy tags, or defers every symbolic comparison. It also contains one explicit expected-value leak.

## Independent N6 capability assessment

The committed PIT tables cannot be cross-subtracted:

- WL uses `{1000000009, 998244353, 1004535809}` with eight draws and its own seeded Mersenne-Twister streams ([WL sampler](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:723), [actual `LOCAL_PROBE`](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S11c_c2_N6_mathematica_audit.out:39)).
- SymPy uses `{1000000009, 998244353, 1000000033}` and hash-oracle sampling ([SymPy sampler](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:668)). The committed covariance, reconcile, and diagnostic runs have respectively 10, 11, and 12 draws and seeds `110603`/`110602` ([covariance](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_c2_N6_covariance_sympy.out:38), [reconcile](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_c2_N6_reconcile_sympy.out:61), [diagnostic](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_c2_N6_diagnostic_sympy.out:60)).
- Thus the sets are different, though not disjoint: the first two primes coincide. The sample assignments do not.

Symbolic operands genuinely exist:

- WL places them in each leaf’s `ARITHMETIC` field ([emitter](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:802), [actual carrier payload](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S11c_c2_N6_mathematica_audit.out:3)).
- SymPy emits a shared DAG plus per-family `columns`/`root_ids` ([DAG emitter](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_sympy.py:151), [actual `R_N6_NODES`](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_c2_N6_reconcile_sympy.out:57)).

The digest seal is also correct: WL’s `Hash[..., "SHA256"]` is emitted as a large integer ([WL emitter](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:957), [actual output](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S11c_c2_N6_mathematica_audit.out:18)); SymPy emits hexadecimal SHA-256 text ([SymPy helper](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:91), [actual output](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_c2_N6_reconcile_sympy.out:16)).

Therefore: sealing the committed raw PIT residues is right. A new comparator-generated common-point evaluation could be separately specified, but it would be a new measurement contract, not a legitimate subtraction of the existing residue tables.

## Ranked findings

### 1. The JOIN SET drops four genuinely shared, computed covariance families

Decision: The covariance join list is incomplete.

Evidence: The directive includes only seven covariance families ([directive](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_comparator_build_directive.md:84)), but both emitters explicitly include:

- `SOURCE_CONTROL_DELTA`
- `R_COV_CONTROL_DELTA`
- `PHI_DOMAIN_CENSUS`
- `ACTUAL_CONTROL_PARAMETERS`

See the complete WL `covNames` list ([WL source](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:30)) and the actual WL tags ([control deltas](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S11c_c2_N6_mathematica_audit.out:25), [census/configuration](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S11c_c2_N6_mathematica_audit.out:29)). SymPy constructs the same control deltas ([SymPy source](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:194)) and emits all four ([census](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_c2_N6_covariance_sympy.out:2), [parameters](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_c2_N6_covariance_sympy.out:16), [deltas](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_c2_N6_covariance_sympy.out:94)).

Why wrong: It removes the covariance independence control, its residual, the substitution-domain coverage guard, and the actual control settings. A comparator could report covariance agreement while never checking that both engines ran the same control or covered the same field-map domain.

Minimal fix: Add all four to the covariance join. Keep `PROVENANCE` separately surfaced with digest subfields sealed.

### 2. The “column↔slot bridge is post-run adjudication” decision confuses schema normalization with physics reconciliation

Decision: The necessary coordinate bridge must be fixed from the emitters before values are compared.

Evidence:

- SymPy carrier columns are `(row, slot, face, grade)` ([SymPy source](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_sympy.py:95)); WL uses `(row, face, slot, {1,η,σ})` and adds face `"SUM"` ([WL source](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:697), [carrier construction](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:712)). SymPy uses face `0` for the computed sum ([SymPy source](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:508)).
- SymPy reconcile sources are already summed over waves and keyed only by `(face,grade)` ([SymPy source](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_sympy.py:105)); WL retains `(face,wave,grade)` columns ([WL source](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:707)).
- SymPy weak columns use `(direction,target)`, numeric signature `0/6/9/12`, and face `-1/0/1` ([SymPy construction](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:491)). WL uses `TARGET_FROM_TRIAL`, named Fourier family, and component vectors ([WL construction](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:552)).

Why wrong: These are deterministic emitter-schema relations, not findings to infer from the measured values. Leaving them to “post-run adjudication” permits either zero joins or value-fitted/manufactured joins. Directly subtracting WL component tuples from SymPy’s formal-basis scalar DAG also produces structural mismatches rather than a measurement.

Minimal fix: Pre-register a source-derived, injective bridge covering axis order, `{1,η,σ}↔(η,σ)`, `"SUM"↔0`, signature names↔`0/6/9/12`, block-name mapping, WL wave-column aggregation, and component/formal-basis reconstruction. Keep only physics-bearing representation differences sealed.

### 3. The primary dimensions channel is extracted from the wrong SymPy container

Decision: `N6RC_DIMENSIONS` alone cannot supply dimensions for every joined family.

Evidence: The directive describes only `data.eulerian/material` in the standalone dimensions tag ([directive](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_comparator_build_directive.md:156)). In reality SymPy attaches a `dimension` array to every numeric emission, aligned with `data.columns` ([SymPy PIT emitter](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:769)). There is no standalone `N6COV_DIMENSIONS` family. WL, conversely, places `DIMENSIONS` beside each leaf’s `ARITHMETIC` ([WL emitter](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:802)); its standalone `N6RC_DIMENSIONS` is a much broader nested census ([WL source](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:986)).

Why wrong: Following the stated trap would omit covariance dimensions and compare unlike global containers for reconcile.

Minimal fix: For every SymPy numeric family, zip the top-level `dimension[i]` with `data.columns[i]`; pair that with the corresponding WL leaf/component `DIMENSIONS`. Treat standalone `N6RC_DIMENSIONS` as an additional guard, not the sole dimensions source.

### 4. The reused loader grammar and duplicate rule reject the real streams

Decision: The c1 `load_wl` cannot be reused as stated, and duplicate detection lacks the required `probe` axis.

Evidence: c1 recognizes colon syntax and splits on `:` ([c1 regex](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c1_cross_engine_comparator.py:55), [loader](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c1_cross_engine_comparator.py:189)). N6 WL uses `WL_S11CC2_... = <|...|>` ([actual first tag](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S11c_c2_N6_mathematica_audit.out:1)).

Moreover, valid SymPy guard rows repeat the same object/anchoring/density with `probe="EULERIAN"` and `"MATERIAL"` ([actual rows](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_c2_N6_diagnostic_sympy.out:124), [second route](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_c2_N6_diagnostic_sympy.out:127)). The route is deliberately supplied as the emitter’s `probe` label ([SymPy source](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:837)).

Why wrong: A literal reuse fails before extraction; duplicate rejection keyed only by object/anchoring/density rejects valid guard routes.

Minimal fix: Specify an N6 WL grammar accepting `=` and define SymPy duplicate identity over all typed top-level labels, including `probe`, then map guard `probe` to the corresponding WL embedded route axis.

### 5. PIT-derived support is not an ordinary two-valued structural fact

Decision: The support channel lacks the required evidence lattice.

Evidence: The directive calls nonzero support a “STRUCTURAL fact” ([directive](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_comparator_build_directive.md:50)). The SymPy engine explicitly says a nonzero modular numerator is one-sided evidence, while all-zero samples mean only “no nonzero found” ([engine contract](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:8)); its emitted support is merely `any(row[i][0] != 0)` ([emitter](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:764)). WL likewise emits an explicitly conditional all-zero semantics ([WL source](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:817)).

Why wrong: Treating `false` as “structurally zero” can create both false agreement and false disagreement across independent probes.

Minimal fix: Type each observation as `NONZERO_WITNESSED` or `NO_NONZERO_FOUND`; the latter must propagate as undecided. Never subtract the native support booleans.

### 6. §N8/T7 is being applied beyond its current scope, and its ≥64 GB deferral is overgeneralized

Decision: The generic behavioral clauses are reusable, but the substantive N6 scope needs its own contract clause.

Evidence: §N8 makes the c2 self-energy increment load-bearing and requires surfacing the §3d representation questions ([governing contract](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:477)); §5c separately defines the fixed-axis N6 two-route construction ([§5c](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:303)). The ≥64 GB entries concern the S11c-b full residual ([deferred record](/var/projects/toy_physics/research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md:16)) and four named c1 families ([c1 entry](/var/projects/toy_physics/research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md:114)), not arbitrary N6 objects.

The N6 WL emitter’s complete family list contains no direct `t_s`, DtN whole-form, or flat-symbol family ([WL family list](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:25)), while the directive excludes provenance/bookkeeping ([directive](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_comparator_build_directive.md:149)).

Why wrong: The builder cannot fulfill the stated §3d surfacing duty from these N6 streams, and the unspecified “fixed per-object size/time budget” can reclassify any N6 symbolic channel as deferred without authority.

Minimal fix: Add an N6-specific contract clause defining its exact families, schema bridges, support semantics, and concrete resource limits. Explicitly place the main-c2 §3d representations outside this N6-only instrument unless new operands are emitted.

### 7. The `_NODES` reader is necessary, but “there is no srepr here” is materially false

Decision: The DAG is reconstructable, but the literal schema is incompletely specified.

Evidence: The four operations are `number`, `variable`, `add`, `mul`, and `pow` as constructed in the diagnostic source ([node constructors](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:115)); references are emitted topologically ([emitter](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_sympy.py:151)). But applied functions and derivatives are encoded as `(sp.srepr(e), "formal_jet")` variable literals ([compiler](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:237)); those strings occur in the committed DAG ([actual DAG](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_c2_N6_reconcile_sympy.out:45)).

Why wrong: A builder believing all literals are simple names may turn formal jets into opaque symbols and manufacture symbolic disagreements.

Minimal fix: Specify all literal variants—`global`, local point, `formal_jet`, `algebraic_i`, and integral-guard—and require a synthetic round-trip containing a real `formal_jet` srepr literal.

### 8. Rule 5 is explicitly breached, and the DoD can accept an all-deferred/all-mismatched instrument

Decision: The directive is not value-free and its coverage gate is non-failing.

Evidence: It states that `REP_INVARIANCE_RESIDUAL` is a “wrong-object — **nonzero is EXPECTED**” ([directive](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_comparator_build_directive.md:70)). The engine instead withholds every expected residual/control value ([engine header](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:12)).

The DoD allows every family to terminate as `axis_set_mismatch` or `deferred_oversize` and expressly refuses any positive-join gate ([DoD](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_comparator_build_directive.md:191)). Combined with the unspecified budget and missing bridge, zero actual comparisons can satisfy the prose.

Why wrong: The expected-nonzero statement contaminates later adjudication; the coverage rule cannot distinguish a cautious comparator from a nonfunctional extractor.

Minimal fix: Delete the quoted expectation. Require a source-derived synthetic fixture for every joined schema and treat zero extracted leaves from two nonempty, declared-shared containers as an operational extractor failure—not a physics disagreement and not a residual target.

## Decisions that are sound

The directive correctly preserves the inherited three-valued residual machinery and bans per-case verdict tokens. The precedent really imports `residual` and its four typed outcomes ([c1 imports](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c1_cross_engine_comparator.py:30)), provides the WL→SymPy canonicalizer ([translator](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c1_cross_engine_comparator.py:462)), and has lazy parse/release objects ([lazy cases](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c1_cross_engine_comparator.py:519)).

The named SymPy-only residuals are also correctly excluded from the cross-engine join, and the raw PIT and digest seals should remain. The directive becomes sound only after the eight corrections above, especially the completed join set and pre-registered schema bridge.
