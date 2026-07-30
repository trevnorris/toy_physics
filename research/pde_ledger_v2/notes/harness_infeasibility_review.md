# DESIGN REVIEW — CROSS-ENGINE DIMENSION HARNESS

## 1. VERDICT

**NOT_FEASIBLE_AS_SPECIFIED.**

The stop condition in §7 applies. The current Wolfram programs do not expose a complete per-quantity value surface, and the combination of R2, R7, and “do not modify `.wl`” leaves no
completeness-verifiable way to create one. A smaller preservation gate is feasible and better matched
to the rewrite.

## 2. FEASIBILITY

**Answer: no—not for every quantity under the stated rules.**

- Corpus census: 44 `.wl`/44 `.out` files = 43 numbered stages plus one midway codimension audit.
  Of the 30 dimension-bearing stages, only these 9 outputs render values produced by dimension
  computation: 004, 011, 012, 013, 016, 018, 021, 023, 027
  (`mathematica/out/ledger_stage004_gnls_action_dimensional_foundation_mathematica_audit.out:9`;
  `mathematica/out/ledger_stage021_dimensional_closure_mathematica_audit.out:17`). The remaining 21 reduce the actual tooth to
  `PASS`; some add hardcoded human prose, which is not a runtime value payload
  (`mathematica/out/ledger_stage031_puncture_deflection_field_identity_source_mathematica_audit.out:48`;
  `mathematica/out/ledger_stage042_charge_coupled_scalar_mathematica_audit.out:68`).
- Thus the “9 of 30” measurement is right as a computed-output census, but even those nine do not
  define a common schema or generally emit every dimension-bearing quantity. The 13 other numbered
  stages have no Python dimensional machinery (5 PARTIAL, 8 UNCOVERED); the 44th file computes
  algebraic codimension, not physical exponents (`mathematica/midway_knob_audit_codimension_mathematica.wl:80`).
- Reading existing `.out` harder cannot recover omitted data. Stage038 computes eight four-axis
  values, including `muDim={2,1,-4,-1}` (`mathematica/ledger_stage038_sealed_landing_electric_bc_r1_mathematica_audit.wl:675`,
  `:697`), but its output prints a fixed summary that omits `mu_R`
  (`mathematica/out/ledger_stage038_sealed_landing_electric_bc_r1_mathematica_audit.out:31`).
- A normal driver cannot `Get[file]; export[...]`: all 44 files call `Exit`. Representative terminal
  exits are `mathematica/ledger_stage004_gnls_action_dimensional_foundation_mathematica_audit.wl:556`
  and `mathematica/ledger_stage042_charge_coupled_scalar_mathematica_audit.wl:2017`.
- I tested the strongest plausible wrapper. Temporarily localising `System\`Exit` lets `Get` return,
  and stage004's exposed `deriveDictionary[]` then yields exact vectors
  (probe output, ⛔ **not retained** — it lived in gitignored `_scratch/dim_harness/review_scratch/PROBE_RESULTS.txt`
  and no copy survives, so this and the next probe result are recorded here, not auditable). This proves a wrapper
  works **where a public constructor already returns the values**.
- It does not generalise. Stage031's 21-row `unitRows` is local to the terminal `Module`
  (`mathematica/ledger_stage031_puncture_deflection_field_identity_source_mathematica_audit.wl:121`,
  `:444`); stage042's 14-entry base map is local to `dimensionGuard`
  (`mathematica/ledger_stage042_charge_coupled_scalar_mathematica_audit.wl:819`,
  `:833`) and the function returns only aggregates, not that map (`:898`). The probe accordingly
  recovered aggregate dimensions but not `B,C,K,qL,...` (same unretained probe output).
- Recovering such locals requires transforming held source/DownValues, evaluator instrumentation, or
  30 stage-specific replicas. That is in-memory modification or reimplementation, not a passive
  `Get`; it perturbs execution, has no independent completeness oracle, and can silently read literals
  instead of values “the engine actually computes.” Treating it as “the file was not edited” is a
  loophole, not satisfaction of R2/R7.

**Honest options:** (1) permit each of the 30 covered `.wl` files to emit a versioned, exact,
axis-labelled value map, accepting a new Wolfram-output baseline; (2) retain the no-edit rule and
scope cross-engine comparison to the values already exposed by the nine stages; or (3) use a
pre-rewrite Python named-value golden baseline as the rewrite-preservation gate, keep the existing
Wolfram/PASS gates, and reserve dual-engine work for targeted high-risk stages. Option 3 buys most of
the rewrite protection at much lower cost.

## 3. FACT CHECK

- **CONFIRMED — 9/30 computed outputs.** The exact nine are listed above. “The other 21 only print
  PASS” needs the qualification that some print hardcoded dimension prose; they still expose no
  complete computed payload.
- **PARTLY — 6 conventions / 4 axis sets / fractions.** The six storage/basis conventions and extremes
  are real: stage008 is two-axis `(L,T)` (`scripts/ledger_stage008_monopole_dipole_return_spec_sympy_audit.py:461`,
  `:484`), stage038 is `(M,L,T,E-charge)`
  (`scripts/ledger_stage038_sealed_landing_electric_bc_r1_sympy_audit.py:692`), and stage042 uses
  `(stiffness,length,time)` (`scripts/ledger_stage042_charge_coupled_scalar_sympy_audit.py:819`,
  `:825`). Eight scripts use fractional dimension
  arithmetic. The phrase “four within `{L,T,M}`” is misleading: seven of the eight use the ordinary
  L/M/T axis set; four are specifically the `(L,M,T)`-ordered cases with non-integral values/probes
  (012, 021, 023, 027; `scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:473`;
  `scripts/ledger_stage027_port_checks_closure_sympy_audit.py:208`).
- **PARTLY — 22/30 transpose.** `22` is 14 `(L,T,M)` plus 8 `(L,M,T)` only if an unstated target order
  `(M,L,T)` is chosen. The committed interchange decision instead makes named maps canonical and says
  no positional order is canonical (`manifests/DIM_ORDER_DECISION.md:10`, `:14`). Stages008/038/042
  need arity/basis policy, not mere transposition. This is a contingent design count, not a corpus fact.
- **CONFIRMED — COVERED 30 · PARTIAL 5 · UNCOVERED 8**, under the report's stated criterion. An AST
  recount gives 3,752 dimensional lines in exactly 30 Python files. The PARTIAL Wolfram files really
  restate dimension-dependent monomials without Python dimension maps (e.g. stage020's `aPower`,
  `mathematica/ledger_stage020_provenance_partition_mathematica_audit.wl:93`; stage024's `gBase`,
  `mathematica/ledger_stage024_density_port_derivation_mathematica_audit.wl:141`). Two sampled
  UNCOVERED stages confirm the distinction: stage017 merely binds `cited016DimensionalOk=True`
  (`mathematica/ledger_stage017_grouped_p2_lane_isotropy_mathematica_audit.wl:34`;
  `scripts/ledger_stage017_grouped_p2_lane_isotropy_sympy_audit.py:62`), while stage043's
  “dimension” is CAD/Krull dimension and explicitly says physical homogeneity is N/A
  (`mathematica/ledger_stage043_irreducible_count_range_mathematica_audit.wl:572`, `:1224`).
- **PARTLY — load-bearing claim.** The hole is real, but the explanation is too universal. Stage038
  compares `unit_state()` with a separately written `EXPECTED_UNIT_STATE`
  (`scripts/ledger_stage038_sealed_landing_electric_bc_r1_sympy_audit.py:708`, `:741`); moving both
  together would pass, and no runner compares the Wolfram values. Stage031 is worse already: 20 of 21
  rows use identical actual/target expressions
  (`scripts/ledger_stage031_puncture_deflection_field_identity_source_sympy_audit.py:554`, `:577`;
  `mathematica/ledger_stage031_puncture_deflection_field_identity_source_mathematica_audit.wl:447`,
  `:470`). But stage004 derives values and checks several independent compositions/absolute bases
  (`scripts/ledger_stage004_gnls_action_dimensional_foundation_sympy_audit.py:185`, `:288`), so not
  every module error necessarily
  moves every target. Strictly, the `.wl` reveals **no** Python error automatically because there is no
  cross-run comparison. A pre-rewrite Python value baseline closes the rewrite-preservation hole
  without requiring full Wolfram extraction.

## 4. REQUIREMENT CRITIQUE

**R3 — semantic comparison.** Named axes are necessary but not sufficient: quantity identity also needs
a stable `(stage, quantity_id)` mapping. Refusal is correct for stage038 versus `{L,T,M}` unless an
explicit physical conversion for `E-charge` is supplied; dropping the extra axis can turn inequivalent
units into a false match. Likewise, `(L,T)` may embed into `{L,T,M}` with `M=0` only if omission is
defined to mean zero, which stage008 does not state. A shared-axis projection may be reported as a
non-verdict diagnostic, but never as agreement. The same rule applies to stage042's stiffness axis.

**R5 — able to fail.** A one-exponent mutation proves the comparator can notice one exposed value; it
does not prove extraction completeness. A harness that extracts one quantity and omits the other 225
can pass that demonstration. It also misses quantity deletion, axis-label swaps applied consistently,
semantic alias errors, and compensating changes. Require at least a real one-sided value mutation, an
axis-label transposition, and a quantity-deletion test against a fixed inventory.

**R6 — categories.** The four categories are neither exhaustive nor safely exclusive. A fifth terminal
case is **indeterminate/extraction error**: the engine computes a value but the adapter cannot recover,
parse, identify, or axis-label it. That must not be called absent coverage. Engine failure/malformed
output and ambiguous identity should be reason codes beneath this case. Define precedence so one item
cannot be both absent and refused.

**R7 — non-perturbation.** It is achievable only with fresh-process extraction and stronger evidence,
not naive imports/repeated calls. Fifteen scripts capture mutation state at import
(`scripts/ledger_stage040_cone_lock_readjudication_sympy_audit.py:29`); stage037 accumulates the
order-sensitive `BUILD_LOG`
(`scripts/ledger_stage037_route_b_boost_structural_relation_sympy_audit.py:242`, `:1094`); stage044
writes a verdict file (`scripts/ledger_stage044_parent_action_reconciliation_sympy_audit.py:1345`).
The directive/measurement is wrong about stage043 having `lru_cache`: current stage043 imports no
cache and `dimension_record` has no decorator
(`scripts/ledger_stage043_irreducible_count_range_sympy_audit.py:15`, `:589`); stage040 does cache
mutation-sensitive work (`scripts/ledger_stage040_cone_lock_readjudication_sympy_audit.py:502`,
`:548`, `:593`).
PASS multisets plus unchanged Wolfram output do not prove C7 teeth or stage044 side effects survived.

**Acceptance §4.4 — cross-stage conflicts.** The claim that this harness would catch them all is false.
The surviving inventory in `notes/rewrite_reference_table.md` §3 catalogues the conflicts; the corpus
context is in `manifests/DIMENSION_REWRITE.md` §7. Same-stage Python-vs-Wolfram values cannot diagnose:

- stage031's self-comparison structure or stage017's hardcoded `True`;
- `r_BA`, which is dimensionless in both incompatible systems (`manifests/DIMENSION_REWRITE.md` §7);
- `epsilon*` versus `Z*_ret`, which needs an external alias assertion (`:97`);
- register locus mis-attribution, which is a register-provenance audit (`:99`);
- `K_eta`, `mu_eta/T_w`, and `M0`, which are cross-**stage**, not cross-engine, disagreements;
- `A_E` across stage037/038, which R3 must refuse absent an `E-charge` conversion.

The raw baseline could support separate structural/register/cross-stage analyses, but none is required.
The “or explain out of scope” escape therefore makes §4.4 weak evidence, not the strongest evidence.

## 5. STANDING QUESTIONS + SCOPE

- **Can acceptance be satisfied without the work? Yes.** Report all unrecovered values as absent, show
  a synthetic one-record mutation, explain every known conflict out of scope, touch no scripts (making
  R7 vacuous), and finish quickly. There is no normative quantity inventory or minimum recovered count.
- **False/synthesised claims?** R2 plus no Wolfram export pressures a builder to mislabel parsed literals
  or reimplemented values as engine-computed. R6 correctly forbids filling gaps. R4 also requests
  symbolic exponents although the corpus has none (`notes/measure_register_sufficiency.md:71`), so its
  proof must be explicitly a synthetic test fixture, not baseline data. The “all nine” benefit claim is
  false as written.
- **Spellings instead of properties?** R8 specifies “YAML or markdown, not JSON”; machine readability is
  a versioned-schema/validator property, and free-form Markdown does not provide it. §4.4 also needs
  stable conflict IDs and expected detection modes, not an ambiguous number.
- **Pre-design?** R8 chooses serialisation, §6 chooses placement, and R5 chooses one demonstration route.
  State lossless exact rationals, stable identity, validation, isolation, and failure properties; leave
  format and adapter architecture to the builder.
- **Worth/correct size?** A full 30-stage Wolfram adapter layer is over-built for preserving a 3,752-line
  Python rewrite, yet under-built for correctness because it cannot see assertion vacuity, aliases,
  incompatible unit systems, or cross-stage conflicts. Use two gates: (A) snapshot every current Python
  quantity into reviewed named-axis records before rewrite and compare after rewrite, plus existing
  PASS-multiset/normalised-Wolfram-output/C7 checks; (B) adjudicate known semantic conflicts separately,
  with targeted dual-engine exports only where they add independent physics evidence. This directly
  catches symmetric module and transposition errors without pretending it proves the old values correct.

## 6. REQUIRED AMENDMENTS

1. Choose explicitly between relaxing the `.wl` no-edit rule for a stable export protocol, or reducing
   R2 to the nine stages/exposed quantities. Reject in-memory source rewriting as a no-edit loophole.
2. Define a versioned, independently reviewed `(stage, quantity_id, axes, exact exponents, source locus)`
   inventory with expected per-stage counts; extraction failure must fail completeness, not become absence.
3. Separate rewrite preservation, same-stage cross-engine agreement, cross-stage identity checking,
   structural assertion audit, and register provenance into distinct scopes and counts.
4. Replace `22/30` with a target-independent named-axis rule. Permit cross-basis conversion only through
   an explicit, injective conversion contract; otherwise refuse and optionally show a non-verdict projection.
5. Add `indeterminate/extraction-error` and reason codes, plus mutually exclusive classification precedence.
6. Strengthen R5 with real-corpus one-sided value, axis-swap, and deletion mutations; require restoration
   and unchanged inventory cardinality after each test.
7. Strengthen R7 to use fresh subprocesses, rerun affected C7 mutations, verify stage037 build-log
   cardinality/order, and hash stage044's side-effect content. Correct the false stage043 cache warning.
8. Replace “nine known conflicts” with stable conflict IDs and, for each, the responsible detector or an
   explicit non-goal. Remove the claim that a same-stage value comparator catches all of them.
9. Replace the YAML/Markdown enumeration with a validated schema for exact rationals, axes, identities,
   provenance, status, and error reasons; generate a separate human summary.

## 7. STRONGEST OBJECTION

The harness spends almost the cost of a second dimension rewrite to manufacture observability that the
independent engine deliberately does not expose, while its advertised acceptance oracle consists mostly
of defects it cannot detect by construction. A reviewed pre-rewrite Python value snapshot catches the
actual rewrite risk—symmetric movement and transposition—more directly; targeted Wolfram changes should
be justified by specific physics disputes, not by an unattainable “every quantity/all nine” promise.
