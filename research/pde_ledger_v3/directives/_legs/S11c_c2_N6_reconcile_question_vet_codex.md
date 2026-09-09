The core question is correct: the open gap is cross-engine agreement of the operands, not agreement of already-cancelled residuals. But the framing needs revision before it can safely govern an instrument.

## Findings

1. **MUST — §4 includes the desired conclusions as candidate “identities,” making collapse potentially vacuous.**

   Framing: “blind-graph-Eulerian-slab-face ↔ imported-slab identity” and “constitutive-source blind-derivation ↔ imported-μ identity” ([question:101](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_question.md:101)).

   Source: those are precisely the independently built operands being tested: SymPy imports and differentiates the slab ([reconcile:190](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_sympy.py:190)), whereas WL constructs and differentiates blind rows ([WL:859](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:859)); likewise SymPy builds actual/predicted sources at [covariance:184](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:184), while WL builds them independently at [WL:843](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:843).

   Directly substituting either whole-object “identity” would assume the answer. Replace them with primitive, independently derivable mappings:

   - Background/profile jets: WL explicitly uses `WBg→W0(1+ηw1)`, `muRBg→muR(1+ηm1)` and their σ-scaled derivatives ([WL:127](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:127)); SymPy implements the corresponding profile and derivative rules at [fold:204](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:204) and [fold:281](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:281).
   - Live density: SymPy rebinds the imported c1 density before source extraction ([diagnostic:378](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:378)); WL defines `density4`, `density3=density4 WBg` ([WL:839](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:839)). This is especially necessary because the c1 precedent made the live-density issue mandatory for c2 ([c1 reconcile:133](/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_reconcile.md:133)).
   - Energy-basis coefficients/contractions: WL generates a quotient basis and coefficient for each retained contraction ([WL:211](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:211)); SymPy varies its independently constructed energy termwise ([diagnostic:318](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:318)). Mapping coefficients by contraction identity is legitimate; equating the resulting μ objects is not.
   - Face velocity and source-solve factors: both enter the emitted source operands directly ([covariance:189](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:189), [WL:859](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:859), [WL:364](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:364)) but are missing from §4.

   The anchoring map may derive formulas within each fixed case, but must not identify `LAB_HELD` with `MATERIAL_ADVECTED`: anchoring is a retained comparator key ([comparator:86](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:86)). On-shell and Fourier-kernel identities appear irrelevant to these particular strong carrier/source/Φ families: SymPy removes the DtN/resolvent before source extraction ([diagnostic:379](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:379)), while WL’s kernel construction begins only after `sourceBind` ([WL:380](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:380)). They should be admitted only if an inspected residual actually contains such atoms.

2. **MUST — §3(b)’s common-bridge idea is right, but its failure inference is too strong and BASELINE is not independent.**

   Framing: failure to collapse means the engines “compute physically different sources” and “would be a genuine cross-engine finding” ([question:74](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_question.md:74)). Yet §6 correctly says non-collapse remains `UNDECIDED`, not disagreement ([question:131](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_question.md:131)).

   Source: at shipped settings, SymPy uses `kappa_a=1`, `kappa_j=0` ([covariance:34](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:34)), making ACTUAL equal to its tag-1 BASELINE ([covariance:137](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:137)). WL likewise gives actual and baseline the same `(1,1,0)` coefficients ([WL:18](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:18), [WL:845](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:845)). The committed tally confirms `SOURCE_CONTROL_DELTA` is 160/0 cross-engine ([tally:14](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_comparator_run_tally.txt:14)).

   Therefore:

   - Require one fixed, predeclared **bridge dictionary**, not literally one identity.
   - Apply it unchanged to ACTUAL and PREDICTED separately.
   - Check BASELINE too, but label it a nominal-control duplicate rather than a third independent discriminator.
   - Collapse earns agreement for the matched emitted source operands.
   - Non-collapse only establishes that the declared bridge did not reconcile them. It remains `UNDECIDED` unless bridge completeness is independently proven.

3. **MUST — retained-order §5 is right, but the proposed “σ_W binding” could itself become a forbidden projection.**

   Framing correctly forbids `σ_W→0` ([question:116](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_question.md:116)), but §4 also proposes “the σ_W binding” ([question:105](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_question.md:105)).

   Source: both engines retain four independent coefficients. SymPy defines `GRADES={(0,0),(0,1),(1,0),(1,1)}` ([diagnostic:43](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:43)) and extracts them independently ([diagnostic:245](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:245)); WL independently truncates η and σ to first order ([WL:123](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:123)). SymPy also deliberately excludes a σ equality from its ordinary profile substitutions ([fold:204](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:204)).

   The collapse must therefore compare all four grades coefficientwise. A physical σ/η relation may be reported after that comparison, but must not combine grades, permit cross-grade cancellation, or turn `ησ` into a discarded `η²` term.

4. **The §2 structural argument is sound, but “zero” needs the PIT qualification.**

   Framing: the within-engine residuals “are 0,” hence their cross-engine differences are trivially zero ([question:54](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_question.md:54)).

   Source confirms the algebraic definitions: SymPy forms `C_E−C_M` at [reconcile:210](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_sympy.py:210) and `ACTUAL−PREDICTED` at [covariance:194](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:194); WL forms the corresponding differences at [WL:874](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:874) and [WL:878](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:878). Thus their equality says nothing about cross-engine operand equality.

   However, the covariance instrument explicitly says all-zero samples mean only “no nonzero found” ([covariance:8](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:8)), and the adopted resolution preserves that qualification ([RESOLVED:23](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_RESOLVED.md:23)). Replace exact “is 0” with “no nonzero found at the retained rectangle, under the adopted PIT-qualified disposition.”

   Also, §1 reverses the comparator sign: it says `WL−SymPy` ([question:16](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_question.md:16)), while the comparator passes SymPy as operand A and WL as operand B ([comparator:794](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:794)), and `residual` computes `py_value−wl_value` ([base comparator:823](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_cross_engine_comparator.py:823)). This does not affect zero/nonzero classifications but should be corrected.

5. **§3(d) should call these load-bearing control/premise records, not merely metadata whose “physics” is to be dismissed.**

   Framing calls them “Provisionally bookkeeping” and asks whether they “carry no physics” ([question:89](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_question.md:89)).

   Source: `ACTUAL_CONTROL_PARAMETERS` records coefficients that directly select the material amplitude and inserted junk ([covariance:137](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:137)); WL’s corresponding parameters directly enter `materialAmplitude`, and `MATERIAL_NORMAL` controls the carrier geometry ([WL:845](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:845), [WL:977](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:977)). `PHI_DOMAIN_CENSUS` carries coverage, uncovered atoms and maximum rank ([covariance:105](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:105)); failures abort construction ([covariance:125](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py:125), [WL:850](/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl:850)).

   Audit these eight nonzero leaves semantically as production-control equivalence and domain-coverage equivalence. They are not source/carrier operands, but neither are they non-load-bearing bookkeeping.

## Exact reframing

Use this as the operative question:

> For each fixed `(anchoring,density)` case and each retained coefficient `(η^i σ_W^j)`, `i,j∈{0,1}`, does one predeclared, source-derived primitive bridge dictionary—containing field/profile-jet, live-density, energy-basis, face-velocity, source-solve and ε-placement conventions, but no whole-carrier, whole-μ or whole-source equality—make the matched cross-engine residuals vanish separately for `C_E`, `C_M`, SOURCE_ACTUAL, SOURCE_PREDICTED and Φ? Apply the same dictionary unchanged to ACTUAL and PREDICTED; report BASELINE as the nominal-control duplicate. Collapse earns agreement only for those matched emitted operands. Failure to collapse remains UNDECIDED unless completeness of the bridge has independently been established.

The authorship/review architecture is otherwise sound: the comparator deliberately omits physics identifications ([comparator:4](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py:4)); the primitive dictionary and instrument must be frozen and reviewed before testing.

For scope, choose **A**. The complete comparator finished all 34 families with zero deferrals, ~330 MB peak RSS and an ~80 MB largest object ([run data:14](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_comparator_run_data.md:14)). Run the collapse family-by-family and grade-by-grade under bounded workers; only an actually oversized individual channel should fall back to B.

**QUESTION NEEDS-WORK — A:** a faithful coefficientwise retained-rectangle test is available on a normal box, but the bridge must first be rewritten as primitive identities and non-collapse must remain UNDECIDED.
