# S11c-c2 SymPy builder report

Status: the full detached production audit and export publication completed. Mechanical checks and scoped reduction-tool evidence are recorded below; physics dispositions remain for the step record and independent reviews.

## Authority and boundary

`directives/S11c_c2_SHARED_PHYSICS.md` matches the bytes at commit `16849fc6`. All of spec §1, its AGREE/UNDECIDED dispositions, the inherited weak-restriction convention, the outgoing-bulk ansatz, the profile/density relations, and the retained grading are SUPPLIED and unfalsifiable within this build. Running this script does not independently validate those inputs. The face-force convention and the #90 closure-fold convention remain in the inherited slot carriers. No cross-engine normalization or disposition is made here.

The requested implementation is `scripts/S11c_c2_selfenergy_fold_sympy_audit.py`; the generated delta is `scripts/S11c_c2_exports.py`. No Wolfram engine, comparator, reviewer, additional builder, or commit was launched. The parent owns the fresh Claude and Grok reviews.

## Frozen symbol map

The source is the `DELTA_P` quantity under `s11c_c1_face_response / CASES / (anchoring, face, density) / VALUE`. The build checks that the two prefixed response roots have producer `S11c-c1`. Neither bare S11b response key is bound. The coefficient export's `PRESSURE` quantity is not a pressure construction operand.

The face representation is `DELTA_W` for both faces of both anchorings. `build_face` (line 477) applies these stages:

1. Bind `rho_br_bg_rho4_constant` to the matching `background_density_map` VALUE before composing the pressure source. The live representative is the imported field relation; the other representative is independently read from its own density-map case.
2. Read `DELTA_P`, divide its source by `epsilon_shape`, and evaluate its ordered resolvent–DtN product using the kernel bridge. The slab slot already supplies its external epsilon.
3. Map each face's c1 `mu_theta` identifier to the matching anchoring/density `mu_theta_operator` amplitude. Both faces use the anchoring's single constitutive operand. That imported operand itself contains epsilon, so its amplitude is obtained by division by epsilon.
4. Map each c1 velocity identifier to the matching `face_velocity / (anchoring, face, DELTA_W)` amplitude, also divided by its imported epsilon. `kinematic_balance` is not used as the velocity operand.
5. Map the c1 whole-form `Z` identifier and resolvent identifier to the computed kernel-algebra matrices. The construction uses `dtn_kernel`, including its two legs. The raw `dtn_operator` is not silently substituted as the construction operator.
6. Compute the normal jet by differentiating the outgoing continuation ansatz with respect to the actual normal coordinate and evaluating at the reference face. Substitute both the pressure and normal-jet slots per face. No extra closed `J` is added.

Every case emits its full pressure, normal jet, density substitution, four identifier-class bindings, inverse operand, inverse matrix, pressure matrix, and second-scattering kernel under `PY_S11CC2_FOLD_SYMBOL_MAP_*`.

Exact application lines in the final script: live-density binding **482–484**; c1 epsilon normalization **487**; mu amplitude **493**, velocity amplitude **494**, and their simultaneous identification **499–500**; the kernel bridge **366–417**, including inverse evaluation **396–398** and ordered intermediate-leg composition **402–416**; pressure integral application **501**; normal differentiation **506–509**; slot map **510–512**; substitution into the full rows **553–555**. `REPRESENTATION = 'DELTA_W'` is fixed at line 47. All eight source cases are bound explicitly by `build_face`; the table below spells out the identifier names. Both faces of each anchoring map to that anchoring's single slab mu operand, resolved through the matching density case of `mu_theta_operator`.

| Anchoring / density / face | c1 mu identifier → slab mu | c1 velocity identifier → face-velocity case | Whole-form identifier / resolvent identifier → kernel bridge |
|---|---|---|---|
| `LAB_HELD / RHO4_CONSTANT / 1` | `s11cc1_mu_theta_lab_held_plus` → `mu_theta_L` | `s11cc1_V_lab_held_plus` → `(LAB_HELD, 1, DELTA_W)` | `s11cc1_dtn_operator_lab_held_plus` / `s11cc1_response_resolvent_lab_held_plus_rho4_constant` → computed per-case matrices |
| `LAB_HELD / RHO4_CONSTANT / -1` | `s11cc1_mu_theta_lab_held_minus` → `mu_theta_L` | `s11cc1_V_lab_held_minus` → `(LAB_HELD, -1, DELTA_W)` | `s11cc1_dtn_operator_lab_held_minus` / `s11cc1_response_resolvent_lab_held_minus_rho4_constant` → computed per-case matrices |
| `LAB_HELD / RHOBR_CONSTANT / 1` | `s11cc1_mu_theta_lab_held_plus` → `mu_theta_L` | `s11cc1_V_lab_held_plus` → `(LAB_HELD, 1, DELTA_W)` | `s11cc1_dtn_operator_lab_held_plus` / `s11cc1_response_resolvent_lab_held_plus_rhobr_constant` → computed per-case matrices |
| `LAB_HELD / RHOBR_CONSTANT / -1` | `s11cc1_mu_theta_lab_held_minus` → `mu_theta_L` | `s11cc1_V_lab_held_minus` → `(LAB_HELD, -1, DELTA_W)` | `s11cc1_dtn_operator_lab_held_minus` / `s11cc1_response_resolvent_lab_held_minus_rhobr_constant` → computed per-case matrices |
| `MATERIAL_ADVECTED / RHO4_CONSTANT / 1` | `s11cc1_mu_theta_material_advected_plus` → `mu_theta_M` | `s11cc1_V_material_advected_plus` → `(MATERIAL_ADVECTED, 1, DELTA_W)` | `s11cc1_dtn_operator_material_advected_plus` / `s11cc1_response_resolvent_material_advected_plus_rho4_constant` → computed per-case matrices |
| `MATERIAL_ADVECTED / RHO4_CONSTANT / -1` | `s11cc1_mu_theta_material_advected_minus` → `mu_theta_M` | `s11cc1_V_material_advected_minus` → `(MATERIAL_ADVECTED, -1, DELTA_W)` | `s11cc1_dtn_operator_material_advected_minus` / `s11cc1_response_resolvent_material_advected_minus_rho4_constant` → computed per-case matrices |
| `MATERIAL_ADVECTED / RHOBR_CONSTANT / 1` | `s11cc1_mu_theta_material_advected_plus` → `mu_theta_M` | `s11cc1_V_material_advected_plus` → `(MATERIAL_ADVECTED, 1, DELTA_W)` | `s11cc1_dtn_operator_material_advected_plus` / `s11cc1_response_resolvent_material_advected_plus_rhobr_constant` → computed per-case matrices |
| `MATERIAL_ADVECTED / RHOBR_CONSTANT / -1` | `s11cc1_mu_theta_material_advected_minus` → `mu_theta_M` | `s11cc1_V_material_advected_minus` → `(MATERIAL_ADVECTED, -1, DELTA_W)` | `s11cc1_dtn_operator_material_advected_minus` / `s11cc1_response_resolvent_material_advected_minus_rhobr_constant` → computed per-case matrices |

The exact code operations are in `build_face`, `kernel_bridge`, and `kernel_apply`; the source index records their spans.

## Computation and representation

| Object / operation | Computation location |
|---|---|
| Input binding / producer verification | `bind_inputs`, line 177 |
| Generated-coefficient dimensions | `infer_dimensions`, line 212 |
| Whole-field trial/jet ansatz | `wave_jet`, line 140 |
| Continuous source-point field map | `at_source`, line 246 |
| Two-leg inverse and intermediate-momentum composition | `kernel_bridge`, line 366 |
| Solved outgoing branch bindings | `outgoing_spectral`, line 438 |
| Fourier profile bindings | `profile_bindings`, line 732 |
| Closed pressure and computed normal jets | `build_face`, line 477 |
| Close before weak extraction | `build_case`, line 532 |
| Both weak off-diagonal blocks | `extract`, line 326 |
| Same-extract subtraction | `difference`, line 100 |
| Rectangular retained-grade computation | `shape_coefficients`, line 579 |
| Native traction / independent face-work pairing | `traction_pairing`, line 772 |
| All import-level controls | `control`, line 765 |
| Modulus-profile form control | `modulus_form`, line 959 |
| Anchoring field pullback and jets | `representation_pullback`, line 980 |
| Imported coupling regression coordinate bridge | `regression_coordinates`, line 1018 |
| Own-row bind-closure publication | `publish`, line 807 |

The row representation contains the assembled strong `U`, `THETA`, and `E_W` balances, retaining all derivatives in the imported expanded rows. Per-face pressure-slot contributions are computed from those same rows. The reported parity quantities are the sum and difference of these per-face contributions. The full imported open coupling root is used only in the ordering regression.

The continuous kernels act through explicit Fourier integrals on arbitrary source fields. Their integrands and ordered kernel compositions are computed. The arbitrary profile/test fields are not assigned a profile or evaluated numerically. The outgoing roots are solved with SymPy and retained as named branch functions with their computed definitions attached to every exported primary. The profile Fourier bindings are also attached. These bindings are part of the serialized object and are required when interpreting or comparing the named functions.

The inverse is evaluated through the inherited rectangular retained grade: the two-leg triangular solve supplies the diagonal and first-scattering pieces; a three-leg triangular solve supplies the mixed second-scattering term, with the intermediate momentum integrated explicitly. Higher shape grades are outside this representation. The normal extension used for the w-jet is explicit in `build_face`. This is a retained-order inverse kernel, not an all-orders inversion on a prescribed numerical background. The exported named branch functions and profile transforms have attached computed definitions; their defining integrals remain on arbitrary fields.

There is no independent adjointness construction in this implementation, so no `SELF_ENERGY_ADJOINTNESS_RESIDUAL` is emitted. Both off-diagonal blocks are emitted. The traction pairing compares the slab's inherited generalized-force contribution to independently composed native covector work; its kinetic/stored operand is the full slab mechanical power with the inherited face generalized-force rows removed. The test velocities use the imported velocity amplitudes, with the external wave epsilon divided out consistently with the pressure map; physical quadratic mechanical power restores that test-amplitude epsilon. This is a comparison of face work, not an independent derivation of all kinetic and stored energies. It does not adjudicate c1's separate far-field ENERGY audit.

The N6 comparison constructs both imported anchorings and applies the density/thickness field map to the full material rows before re-extraction. Its two one-sided probes each read the anchoring's native face-normal covector, flip the upper face's first slope component, and propagate the computed component ratio into the corresponding first Fourier tilt of that anchoring's imported DtN kernel. The original and corrupted native normals are emitted with each probe. The companion anchoring is emitted from its separately constructed baseline; no difference of an object with itself is presented as a companion test. These precise probe scopes are exposed for review; they are not a claim that the unresolved c1 density, traction, whole-form, leg-labeling, or ENERGY dispositions have closed.

## Import guard evidence

The positional two-parent fold measured `2485` rows, with source counts `[['/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_exports.py', 2441], ['/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c1_exports.py', 44]]` and overwrites `[]`. The actual access-recording call is `assert_lookups_equal_manifest(bind_inputs, fold, IMPORT_KEYS)`, and its observed lookup set equals the manifest on these files. `check_consumer` resolves `543` keys. The literal evidence is `_measurements/S11c_c2_sympy_guard_evidence.json`.

Final `IMPORT_KEYS` (34 keys):

- `L_W`
- `Lambda_A_0`
- `Lambda_V_0`
- `Lambda_X_0`
- `W_0`
- `W_bg`
- `background_density_map`
- `c_s0`
- `closure_shape_deriv`
- `conormal_deriv`
- `coupling_kernel`
- `dtn_flat_symbol`
- `dtn_kernel`
- `dtn_operator`
- `epsilon_shape`
- `eta_bg`
- `face_measure_shape_deriv`
- `face_normal`
- `face_shift`
- `face_velocity`
- `kinematic_balance`
- `mu_R`
- `mu_theta_operator`
- `omega`
- `relative_flux`
- `rho_br`
- `rho_br_bg_rho4_constant`
- `rho_m`
- `s11c_c1_face_response`
- `s11c_c1_face_response_coeffs`
- `sigma_W`
- `slab_operator`
- `slab_operator_term_origins`
- `traction`

## Export and measurement evidence

The declared outgoing roots are `s11cc2ClosedSlabOperator`, `s11cc2ClosedCouplingKernel`, and `s11cc2SelfEnergyIncrement`. Publication computes their own bind-closure over the real parent fold, adds referenced new coordinate/function declarations and their dimension rows, checks exact-key collisions, runs `assert_delta_is_minimal`, and checks serialization against the live values. The generated module pins the five required input digests. The runtime evidence records the exact exported key set and guard operands.

Production stdout is `/tmp/S11c_c2_selfenergy_fold_sympy_audit.out`. Each §5 control has its own literal operand and residual tags there. No residual is asserted against a physics target. The completion marker is `_measurements/S11c_c2_sympy_completion.txt`; `/usr/bin/time -v` measures the Python process in `_measurements/S11c_c2_sympy_resources.txt`; `_measurements/S11c_c2_sympy_monitor.tsv` samples that process's PID, elapsed seconds, RSS, VSZ, and CPU percentage. `_measurements/S11c_c2_sympy_progress.json` records the current case/control.

`_measurements/S11c_c2_sympy_source_index.json` records the script hash, line count, Python lexical-token count, and exact function spans. The runtime tag index provides the exact computation/emit location and fresh write-key for each emitted object. Prior development attempts and the single-case smoke run are separate `/tmp/S11c_c2_attempt*` / `/tmp/s11cc2_smoke*` files; they are not the production completion evidence.

## Final measured evidence

The production completion marker contains `0`. The audit contains 1053 lines, 82,702 bytes, and 22,953 Python lexical tokens (a source-size plausibility measure, not a model-token billing count). Its SHA-256 is `6360721815b0b8c4a52393c388c37b20c858bca04615d8e671964af49a87c7dd`. The plain, non-symlink export contains 71 own rows and 60,516,900 bytes. Publication reloaded the serialized module and emitted its literal equality tuple before the round-trip guard.

The actual lookup manifest has 34 keys, exactly matching the recorded accesses; its consumer closure has 543 keys. The three-parent fold has 2556 rows. Own delta membership equals the computed own bind-closure, imported-key collisions are `[]`, and the independent file reread records all five digest comparisons in `_measurements/S11c_c2_postbuild_verify.json`. Closed primary pressure/normal-slot references are `[]`; raw c1 whole-form/resolvent references are `[]`.

The stdout index records 154 unique tags across 498,811,405 bytes, duplicate tags `[]`, and 0 untagged records. Every record carries a computed multigrade and dimension object. The complete literal metadata, branch definitions, and residual values remain in the indexed stdout; they are not converted into physics dispositions here.

The streaming dimension-metadata scan found missing dimension fields on `[]` and literal `nan` entries on `[]`. This records the dimension walker's output; the independent corruption probe below tests its sensitivity.

Python PID(s) sampled: `['2439638']`. Runtime self-measurement: 1159.807 seconds and 2,689,720 KiB peak RSS. All requested in-script controls ran serially; none was deferred. The inherited c1 giant-family and ≥64 GB cross-engine comparisons remain outside this build. The final `/usr/bin/time -v` record is:

```text
	Command being timed: "python -u /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py"
	User time (seconds): 1159.48
	System time (seconds): 2.33
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 19:22.18
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 2689720
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 5
	Minor (reclaiming a frame) page faults: 1784272
	Voluntary context switches: 24
	Involuntary context switches: 14767
	Swaps: 0
	File system inputs: 184
	File system outputs: 1093544
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

The read-only verification script emits the base and deliberately corrupted dimensional operands and their equality boolean. It separately compares the retained-grade implementation to SymPy series on the same inherited live-density pressure operand, then changes only that operand’s density binding. Literal probes are in `_measurements/S11c_c2_postbuild_verify.out`; these method checks do not independently verify the supplied physics.

## Reduction tool scope and literal outcomes

The reduction tools are absent from this active checkout. The repository copies used are `/var/projects/toy_physics_review_w3/baseline/HEAD_tree/research/pde_ledger_v3/reduction/derived_or_declared.py` and `engine_output_checks.py`. Their literal diagnostics, resource records, and completion markers are under `_measurements/S11c_c2_*`. No output was declared a supplied premise to force a triage result.

`derived_or_declared.py` was run on the deliverable with `S11CC2_PACKAGE=TRIAGE`. That package computes the real LAB_HELD/RHO4 closed slab, closed kernel, and increment, plus dimensional bindings; it selects four emitted quantities and does not publish an export. This is a scoped dependency triage, not coverage of the complete control suite. The legacy tool has a 590-second child timeout; the full production run exceeds that window. Literal diagnostics (large constant payloads remain in the tool stdout):

```text
BASELINE: RAN tags=4 symbols=2046
PERTURBATION collapse(s11cc2X1=s11cc2X2): SKIPPED reason=exit=1; last_output="TypeError: unsupported operand type(s) for -: 'StrictGreaterThan' and 'StrictLessThan'"
PERTURBATION collapse(s11cc2X3=s11cc2Y1): RAN tags=4 comparable=4
PERTURBATION collapse(s11cc2Y2=s11cc2Y3): SKIPPED reason=exit=1; last_output="TypeError: unsupported operand type(s) for -: 'StrictGreaterThan' and 'StrictLessThan'"
PERTURBATION collapse(s11cc2Time=s11cc2NormalCoordinate): RAN tags=4 comparable=4
PERTURBATION collapse(s11cc2MiddleMomentum1=s11cc2MiddleMomentum2): RAN tags=4 comparable=4
PERTURBATION collapse(s11cc2MiddleNormalMomentum=A_T_s11cb_1): RAN tags=4 comparable=4
PERTURBATIONS_RAN: 4/6
CONSTANT_COUNT: 1
CONSTANT_ADJUDICATION: this list is for adjudication; a tag may legitimately not depend on any collapsed pair, and CONSTANT is not proof of literal source text
DERIVED_TAGS: count=3 tags=PY_S11CC2_CLOSED_SLAB_OPERATOR_LAB_HELD_RHO4_CONSTANT,PY_S11CC2_CLOSED_COUPLING_KERNEL_LAB_HELD_RHO4_CONSTANT,PY_S11CC2_SELF_ENERGY_INCREMENT_LAB_HELD_RHO4_CONSTANT
SUPPLIED_PREMISES: path=/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.premises status=MISSING (treated as an empty declaration) declared=0
UNDECLARED_CONSTANTS: count=1 tags=PY_S11CC2_DIMENSION_COEFFICIENT_BINDINGS
DERIVED_OR_DECLARED: FAIL
```

The legacy output checker was run against the full production stdout with the c2 declaration file. Its parser is restricted to `S` followed by digits and an underscore, so it does not recognize the required `S11CC2` prefix. The full output-checker result is an operational limitation, not successful validation. The separate naming ablation uses an unchanged literal payload and changes only the prefix to measure that parser boundary. Summary:

```text
stdout: /tmp/S11c_c2_legacy_output_checks.out
OPERATIONAL_FAILURE: cannot parse tagged output /tmp/S11c_c2_selfenergy_fold_sympy_audit.out: engine emitted no valid tagged records
```

Literal naming-ablation evidence:

```json
{
  "source": "/tmp/S11c_c2_selfenergy_fold_sympy_audit.out",
  "literal_line": "PY_S11CC2_DTN_WHOLEFORM_DEPENDENCE_LAB_HELD_RHO4_CONSTANT: Tuple(Tuple(Str('VALUE'), Tuple()), Tuple(Str('MULTIGRADE'), Tuple()), Tuple(Str('DIMENSION_L_T_M'), Tuple()))",
  "literal_keys": [],
  "literal_parser_error": "engine emitted no valid tagged records",
  "renamed_keys": [
    "PY_S11_DTN_WHOLEFORM_DEPENDENCE_LAB_HELD_RHO4_CONSTANT"
  ],
  "renamed_ignored_count": 0,
  "normalization_types": {
    "PY_S11_DTN_WHOLEFORM_DEPENDENCE_LAB_HELD_RHO4_CONSTANT": "ParsedValue"
  }
}
```

## Per-object computation and literal control index

All computation/emit numbers below refer to `scripts/S11c_c2_selfenergy_fold_sympy_audit.py` at the SHA above. The computation column identifies the algebra or imported-operand rerun that produces the object, rather than merely its print call. The byte location is zero-based offset / byte count in `/tmp/S11c_c2_selfenergy_fold_sympy_audit.out`; every §5 residual is referenced literally here. Full JSON with hashes and actual multigrades is `_measurements/S11c_c2_sympy_object_provenance.json`. Each per-case emit write-key is fresh and injective; the delta’s three aggregate root keys are separately listed above.

| Literal stdout tag | Computation lines | Emit line | Fresh write-key | Stdout offset / bytes |
|---|---|---|---|---|
| `PY_S11CC2_DIMENSION_COEFFICIENT_BINDINGS` | 234 | 866 | `s11cc2DimensionCoefficientBindings` | 0 / 11735 |
| `PY_S11CC2_CLOSED_SLAB_OPERATOR_LAB_HELD_RHO4_CONSTANT` | 553, 559 | 893 | `s11cc2ClosedSlabOperatorLabHeldRho4Constant` | 11735 / 4013128 |
| `PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_LAB_HELD_RHO4_CONSTANT` | 550, 560 | 893 | `s11cc2ClosedSlabOperatorTermOriginsLabHeldRho4Constant` | 4024863 / 4396607 |
| `PY_S11CC2_CLOSED_SLAB_OPERATOR_PARITY_BLOCKS_LAB_HELD_RHO4_CONSTANT` | 884, 885 | 893 | `s11cc2ClosedSlabOperatorParityBlocksLabHeldRho4Constant` | 8421470 / 3533252 |
| `PY_S11CC2_CLOSED_COUPLING_KERNEL_LAB_HELD_RHO4_CONSTANT` | 561 | 893 | `s11cc2ClosedCouplingKernelLabHeldRho4Constant` | 11954722 / 4366525 |
| `PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_LAB_HELD_RHO4_CONSTANT` | 887 | 893 | `s11cc2ClosedCouplingKernelTermOriginsLabHeldRho4Constant` | 16321247 / 4189623 |
| `PY_S11CC2_SELF_ENERGY_CLOSED_EXTRACT_OPERAND_LAB_HELD_RHO4_CONSTANT` | 561 | 893 | `s11cc2SelfEnergyClosedExtractOperandLabHeldRho4Constant` | 20510870 / 4356956 |
| `PY_S11CC2_SELF_ENERGY_OPEN_EXTRACT_OPERAND_LAB_HELD_RHO4_CONSTANT` | 562 | 893 | `s11cc2SelfEnergyOpenExtractOperandLabHeldRho4Constant` | 24867826 / 1812958 |
| `PY_S11CC2_SELF_ENERGY_INCREMENT_LAB_HELD_RHO4_CONSTANT` | 563 | 893 | `s11cc2SelfEnergyIncrementLabHeldRho4Constant` | 26680784 / 2575252 |
| `PY_S11CC2_FOLD_SYMBOL_MAP_LAB_HELD_RHO4_CONSTANT` | 397, 484, 499, 501, 509 | 893 | `s11cc2FoldSymbolMapLabHeldRho4Constant` | 29256036 / 3535040 |
| `PY_S11CC2_ORDERING_EXTRACT_FIRST_OPERAND_LAB_HELD_RHO4_CONSTANT` | 901, 904 | 905 | `s11cc2OrderingExtractFirstOperandLabHeldRho4Constant` | 32791076 / 4227179 |
| `PY_S11CC2_ORDERING_COMMUTATOR_LAB_HELD_RHO4_CONSTANT` | 906 | 906 | `s11cc2OrderingCommutatorLabHeldRho4Constant` | 37018255 / 431555 |
| `PY_S11CC2_TRACTION_MECHANICAL_CONTRIB_LAB_HELD_RHO4_CONSTANT` | 780, 784 | 908 | `s11cc2TractionMechanicalContribLabHeldRho4Constant` | 37449810 / 3820491 |
| `PY_S11CC2_TRACTION_SLAB_POWER_PAIRING_LAB_HELD_RHO4_CONSTANT` | 792, 798, 799, 800, 801 | 909 | `s11cc2TractionSlabPowerPairingLabHeldRho4Constant` | 41270301 / 13730136 |
| `PY_S11CC2_TRACTION_SLAB_POWER_PAIRING_RESIDUAL_LAB_HELD_RHO4_CONSTANT` | 803 | 910 | `s11cc2TractionSlabPowerPairingResidualLabHeldRho4Constant` | 55000437 / 11622816 |
| `PY_S11CC2_TRACTION_SIGN_OPERAND_LAB_HELD_RHO4_CONSTANT` | 779, 792, 801 | 912 | `s11cc2TractionSignOperandLabHeldRho4Constant` | 66623253 / 13730299 |
| `PY_S11CC2_TRACTION_SIGN_RESIDUAL_LAB_HELD_RHO4_CONSTANT` | 913 | 913 | `s11cc2TractionSignResidualLabHeldRho4Constant` | 80353552 / 7918643 |
| `PY_S11CC2_DTN_WHOLEFORM_DEPENDENCE_LAB_HELD_RHO4_CONSTANT` | 915, 916 | 916 | `s11cc2DtnWholeformDependenceLabHeldRho4Constant` | 88272195 / 170 |
| `PY_S11CC2_FLAT_SYMBOL_USAGE_LAB_HELD_RHO4_CONSTANT` | 918, 919 | 919 | `s11cc2FlatSymbolUsageLabHeldRho4Constant` | 88272365 / 804 |
| `PY_S11CC2_ROUTING_MECHANICAL_ONLY_OPERAND_LAB_HELD_RHO4_CONSTANT` | 563, 766, 920 | 767 | `s11cc2RoutingMechanicalOnlyOperandLabHeldRho4Constant` | 88273169 / 1709943 |
| `PY_S11CC2_ROUTING_MECHANICAL_ONLY_RESIDUAL_LAB_HELD_RHO4_CONSTANT` | 768, 920 | 768 | `s11cc2RoutingMechanicalOnlyResidualLabHeldRho4Constant` | 89983112 / 856943 |
| `PY_S11CC2_ROUTING_TRACTION_CHANNEL_OPERAND_LAB_HELD_RHO4_CONSTANT` | 563, 766, 921 | 767 | `s11cc2RoutingTractionChannelOperandLabHeldRho4Constant` | 90840055 / 1709413 |
| `PY_S11CC2_ROUTING_TRACTION_CHANNEL_RESIDUAL_LAB_HELD_RHO4_CONSTANT` | 768, 921 | 768 | `s11cc2RoutingTractionChannelResidualLabHeldRho4Constant` | 92549468 / 857626 |
| `PY_S11CC2_ZERO_DTN_OPERAND_LAB_HELD_RHO4_CONSTANT` | 563, 766, 922 | 767 | `s11cc2ZeroDtnOperandLabHeldRho4Constant` | 93407094 / 6282 |
| `PY_S11CC2_ZERO_DTN_RESIDUAL_LAB_HELD_RHO4_CONSTANT` | 768, 922 | 768 | `s11cc2ZeroDtnResidualLabHeldRho4Constant` | 93413376 / 2560452 |
| `PY_S11CC2_LAMBDA_A_LIMIT_OPERAND_LAB_HELD_RHO4_CONSTANT` | 563, 766, 923 | 767 | `s11cc2LambdaALimitOperandLabHeldRho4Constant` | 95973828 / 4406 |
| `PY_S11CC2_LAMBDA_A_LIMIT_RESIDUAL_LAB_HELD_RHO4_CONSTANT` | 768, 923 | 768 | `s11cc2LambdaALimitResidualLabHeldRho4Constant` | 95978234 / 2562371 |
| `PY_S11CC2_IMPERMEABLE_LIMIT_OPERAND_LAB_HELD_RHO4_CONSTANT` | 563, 766, 924 | 767 | `s11cc2ImpermeableLimitOperandLabHeldRho4Constant` | 98540605 / 4409 |
| `PY_S11CC2_IMPERMEABLE_LIMIT_RESIDUAL_LAB_HELD_RHO4_CONSTANT` | 768, 924 | 768 | `s11cc2ImpermeableLimitResidualLabHeldRho4Constant` | 98545014 / 2562374 |
| `PY_S11CC2_UNIFORM_LIMIT_OPERAND_LAB_HELD_RHO4_CONSTANT` | 563, 766, 925 | 767 | `s11cc2UniformLimitOperandLabHeldRho4Constant` | 101107388 / 16655 |
| `PY_S11CC2_UNIFORM_LIMIT_RESIDUAL_LAB_HELD_RHO4_CONSTANT` | 768, 925 | 768 | `s11cc2UniformLimitResidualLabHeldRho4Constant` | 101124043 / 2550016 |
| `PY_S11CC2_MU_R_FORM_OPERAND_LAB_HELD_RHO4_CONSTANT` | 563, 766, 926 | 767 | `s11cc2MuRFormOperandLabHeldRho4Constant` | 103674059 / 3288349 |
| `PY_S11CC2_MU_R_FORM_RESIDUAL_LAB_HELD_RHO4_CONSTANT` | 768, 926 | 768 | `s11cc2MuRFormResidualLabHeldRho4Constant` | 106962408 / 5319734 |
| `PY_S11CC2_CLOSED_SLAB_OPERATOR_LAB_HELD_RHOBR_CONSTANT` | 553, 559 | 893 | `s11cc2ClosedSlabOperatorLabHeldRhobrConstant` | 112282142 / 3131875 |
| `PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_LAB_HELD_RHOBR_CONSTANT` | 550, 560 | 893 | `s11cc2ClosedSlabOperatorTermOriginsLabHeldRhobrConstant` | 115414017 / 3056958 |
| `PY_S11CC2_CLOSED_SLAB_OPERATOR_PARITY_BLOCKS_LAB_HELD_RHOBR_CONSTANT` | 884, 885 | 893 | `s11cc2ClosedSlabOperatorParityBlocksLabHeldRhobrConstant` | 118470975 / 2863428 |
| `PY_S11CC2_CLOSED_COUPLING_KERNEL_LAB_HELD_RHOBR_CONSTANT` | 561 | 893 | `s11cc2ClosedCouplingKernelLabHeldRhobrConstant` | 121334403 / 3417565 |
| `PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_LAB_HELD_RHOBR_CONSTANT` | 887 | 893 | `s11cc2ClosedCouplingKernelTermOriginsLabHeldRhobrConstant` | 124751968 / 2904370 |
| `PY_S11CC2_SELF_ENERGY_CLOSED_EXTRACT_OPERAND_LAB_HELD_RHOBR_CONSTANT` | 561 | 893 | `s11cc2SelfEnergyClosedExtractOperandLabHeldRhobrConstant` | 127656338 / 3407996 |
| `PY_S11CC2_SELF_ENERGY_OPEN_EXTRACT_OPERAND_LAB_HELD_RHOBR_CONSTANT` | 562 | 893 | `s11cc2SelfEnergyOpenExtractOperandLabHeldRhobrConstant` | 131064334 / 1506401 |
| `PY_S11CC2_SELF_ENERGY_INCREMENT_LAB_HELD_RHOBR_CONSTANT` | 563 | 893 | `s11cc2SelfEnergyIncrementLabHeldRhobrConstant` | 132570735 / 1930219 |
| `PY_S11CC2_FOLD_SYMBOL_MAP_LAB_HELD_RHOBR_CONSTANT` | 397, 484, 499, 501, 509 | 893 | `s11cc2FoldSymbolMapLabHeldRhobrConstant` | 134500954 / 2221101 |
| `PY_S11CC2_ORDERING_EXTRACT_FIRST_OPERAND_LAB_HELD_RHOBR_CONSTANT` | 901, 904 | 905 | `s11cc2OrderingExtractFirstOperandLabHeldRhobrConstant` | 136722055 / 3387667 |
| `PY_S11CC2_ORDERING_COMMUTATOR_LAB_HELD_RHOBR_CONSTANT` | 906 | 906 | `s11cc2OrderingCommutatorLabHeldRhobrConstant` | 140109722 / 402185 |
| `PY_S11CC2_TRACTION_MECHANICAL_CONTRIB_LAB_HELD_RHOBR_CONSTANT` | 780, 784 | 908 | `s11cc2TractionMechanicalContribLabHeldRhobrConstant` | 140511907 / 2315371 |
| `PY_S11CC2_TRACTION_SLAB_POWER_PAIRING_LAB_HELD_RHOBR_CONSTANT` | 792, 798, 799, 800, 801 | 909 | `s11cc2TractionSlabPowerPairingLabHeldRhobrConstant` | 142827278 / 9957583 |
| `PY_S11CC2_TRACTION_SLAB_POWER_PAIRING_RESIDUAL_LAB_HELD_RHOBR_CONSTANT` | 803 | 910 | `s11cc2TractionSlabPowerPairingResidualLabHeldRhobrConstant` | 152784861 / 8392937 |
| `PY_S11CC2_TRACTION_SIGN_OPERAND_LAB_HELD_RHOBR_CONSTANT` | 779, 792, 801 | 912 | `s11cc2TractionSignOperandLabHeldRhobrConstant` | 161177798 / 9957746 |
| `PY_S11CC2_TRACTION_SIGN_RESIDUAL_LAB_HELD_RHOBR_CONSTANT` | 913 | 913 | `s11cc2TractionSignResidualLabHeldRhobrConstant` | 171135544 / 4774734 |
| `PY_S11CC2_DTN_WHOLEFORM_DEPENDENCE_LAB_HELD_RHOBR_CONSTANT` | 915, 916 | 916 | `s11cc2DtnWholeformDependenceLabHeldRhobrConstant` | 175910278 / 171 |
| `PY_S11CC2_FLAT_SYMBOL_USAGE_LAB_HELD_RHOBR_CONSTANT` | 918, 919 | 919 | `s11cc2FlatSymbolUsageLabHeldRhobrConstant` | 175910449 / 805 |
| `PY_S11CC2_ROUTING_MECHANICAL_ONLY_OPERAND_LAB_HELD_RHOBR_CONSTANT` | 563, 766, 920 | 767 | `s11cc2RoutingMechanicalOnlyOperandLabHeldRhobrConstant` | 175911254 / 1280020 |
| `PY_S11CC2_ROUTING_MECHANICAL_ONLY_RESIDUAL_LAB_HELD_RHOBR_CONSTANT` | 768, 920 | 768 | `s11cc2RoutingMechanicalOnlyResidualLabHeldRhobrConstant` | 177191274 / 641860 |
| `PY_S11CC2_ROUTING_TRACTION_CHANNEL_OPERAND_LAB_HELD_RHOBR_CONSTANT` | 563, 766, 921 | 767 | `s11cc2RoutingTractionChannelOperandLabHeldRhobrConstant` | 177833134 / 1279564 |
| `PY_S11CC2_ROUTING_TRACTION_CHANNEL_RESIDUAL_LAB_HELD_RHOBR_CONSTANT` | 768, 921 | 768 | `s11cc2RoutingTractionChannelResidualLabHeldRhobrConstant` | 179112698 / 642417 |
| `PY_S11CC2_ZERO_DTN_OPERAND_LAB_HELD_RHOBR_CONSTANT` | 563, 766, 922 | 767 | `s11cc2ZeroDtnOperandLabHeldRhobrConstant` | 179755115 / 6283 |
| `PY_S11CC2_ZERO_DTN_RESIDUAL_LAB_HELD_RHOBR_CONSTANT` | 768, 922 | 768 | `s11cc2ZeroDtnResidualLabHeldRhobrConstant` | 179761398 / 1915445 |
| `PY_S11CC2_LAMBDA_A_LIMIT_OPERAND_LAB_HELD_RHOBR_CONSTANT` | 563, 766, 923 | 767 | `s11cc2LambdaALimitOperandLabHeldRhobrConstant` | 181676843 / 4407 |
| `PY_S11CC2_LAMBDA_A_LIMIT_RESIDUAL_LAB_HELD_RHOBR_CONSTANT` | 768, 923 | 768 | `s11cc2LambdaALimitResidualLabHeldRhobrConstant` | 181681250 / 1917364 |
| `PY_S11CC2_IMPERMEABLE_LIMIT_OPERAND_LAB_HELD_RHOBR_CONSTANT` | 563, 766, 924 | 767 | `s11cc2ImpermeableLimitOperandLabHeldRhobrConstant` | 183598614 / 4410 |
| `PY_S11CC2_IMPERMEABLE_LIMIT_RESIDUAL_LAB_HELD_RHOBR_CONSTANT` | 768, 924 | 768 | `s11cc2ImpermeableLimitResidualLabHeldRhobrConstant` | 183603024 / 1917367 |
| `PY_S11CC2_UNIFORM_LIMIT_OPERAND_LAB_HELD_RHOBR_CONSTANT` | 563, 766, 925 | 767 | `s11cc2UniformLimitOperandLabHeldRhobrConstant` | 185520391 / 16656 |
| `PY_S11CC2_UNIFORM_LIMIT_RESIDUAL_LAB_HELD_RHOBR_CONSTANT` | 768, 925 | 768 | `s11cc2UniformLimitResidualLabHeldRhobrConstant` | 185537047 / 1905009 |
| `PY_S11CC2_MU_R_FORM_OPERAND_LAB_HELD_RHOBR_CONSTANT` | 563, 766, 926 | 767 | `s11cc2MuRFormOperandLabHeldRhobrConstant` | 187442056 / 2460094 |
| `PY_S11CC2_MU_R_FORM_RESIDUAL_LAB_HELD_RHOBR_CONSTANT` | 768, 926 | 768 | `s11cc2MuRFormResidualLabHeldRhobrConstant` | 189902150 / 3882204 |
| `PY_S11CC2_CLOSED_SLAB_OPERATOR_MATERIAL_ADVECTED_RHO4_CONSTANT` | 553, 559 | 893 | `s11cc2ClosedSlabOperatorMaterialAdvectedRho4Constant` | 193784354 / 4437985 |
| `PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_MATERIAL_ADVECTED_RHO4_CONSTANT` | 550, 560 | 893 | `s11cc2ClosedSlabOperatorTermOriginsMaterialAdvectedRho4Constant` | 198222339 / 5050165 |
| `PY_S11CC2_CLOSED_SLAB_OPERATOR_PARITY_BLOCKS_MATERIAL_ADVECTED_RHO4_CONSTANT` | 884, 885 | 893 | `s11cc2ClosedSlabOperatorParityBlocksMaterialAdvectedRho4Constant` | 203272504 / 4040716 |
| `PY_S11CC2_CLOSED_COUPLING_KERNEL_MATERIAL_ADVECTED_RHO4_CONSTANT` | 561 | 893 | `s11cc2ClosedCouplingKernelMaterialAdvectedRho4Constant` | 207313220 / 6000541 |
| `PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_MATERIAL_ADVECTED_RHO4_CONSTANT` | 887 | 893 | `s11cc2ClosedCouplingKernelTermOriginsMaterialAdvectedRho4Constant` | 213313761 / 9173840 |
| `PY_S11CC2_SELF_ENERGY_CLOSED_EXTRACT_OPERAND_MATERIAL_ADVECTED_RHO4_CONSTANT` | 561 | 893 | `s11cc2SelfEnergyClosedExtractOperandMaterialAdvectedRho4Constant` | 222487601 / 5990972 |
| `PY_S11CC2_SELF_ENERGY_OPEN_EXTRACT_OPERAND_MATERIAL_ADVECTED_RHO4_CONSTANT` | 562 | 893 | `s11cc2SelfEnergyOpenExtractOperandMaterialAdvectedRho4Constant` | 228478573 / 1895584 |
| `PY_S11CC2_SELF_ENERGY_INCREMENT_MATERIAL_ADVECTED_RHO4_CONSTANT` | 563 | 893 | `s11cc2SelfEnergyIncrementMaterialAdvectedRho4Constant` | 230374157 / 4204540 |
| `PY_S11CC2_FOLD_SYMBOL_MAP_MATERIAL_ADVECTED_RHO4_CONSTANT` | 397, 484, 499, 501, 509 | 893 | `s11cc2FoldSymbolMapMaterialAdvectedRho4Constant` | 234578697 / 3579897 |
| `PY_S11CC2_ORDERING_EXTRACT_FIRST_OPERAND_MATERIAL_ADVECTED_RHO4_CONSTANT` | 901, 904 | 905 | `s11cc2OrderingExtractFirstOperandMaterialAdvectedRho4Constant` | 238158594 / 4551148 |
| `PY_S11CC2_ORDERING_COMMUTATOR_MATERIAL_ADVECTED_RHO4_CONSTANT` | 906 | 906 | `s11cc2OrderingCommutatorMaterialAdvectedRho4Constant` | 242709742 / 1788273 |
| `PY_S11CC2_TRACTION_MECHANICAL_CONTRIB_MATERIAL_ADVECTED_RHO4_CONSTANT` | 780, 784 | 908 | `s11cc2TractionMechanicalContribMaterialAdvectedRho4Constant` | 244498015 / 3862584 |
| `PY_S11CC2_TRACTION_SLAB_POWER_PAIRING_MATERIAL_ADVECTED_RHO4_CONSTANT` | 792, 798, 799, 800, 801 | 909 | `s11cc2TractionSlabPowerPairingMaterialAdvectedRho4Constant` | 248360599 / 15232502 |
| `PY_S11CC2_TRACTION_SLAB_POWER_PAIRING_RESIDUAL_MATERIAL_ADVECTED_RHO4_CONSTANT` | 803 | 910 | `s11cc2TractionSlabPowerPairingResidualMaterialAdvectedRho4Constant` | 263593101 / 12698484 |
| `PY_S11CC2_TRACTION_SIGN_OPERAND_MATERIAL_ADVECTED_RHO4_CONSTANT` | 779, 792, 801 | 912 | `s11cc2TractionSignOperandMaterialAdvectedRho4Constant` | 276291585 / 15232613 |
| `PY_S11CC2_TRACTION_SIGN_RESIDUAL_MATERIAL_ADVECTED_RHO4_CONSTANT` | 913 | 913 | `s11cc2TractionSignResidualMaterialAdvectedRho4Constant` | 291524198 / 7573948 |
| `PY_S11CC2_DTN_WHOLEFORM_DEPENDENCE_MATERIAL_ADVECTED_RHO4_CONSTANT` | 915, 916 | 916 | `s11cc2DtnWholeformDependenceMaterialAdvectedRho4Constant` | 299098146 / 179 |
| `PY_S11CC2_FLAT_SYMBOL_USAGE_MATERIAL_ADVECTED_RHO4_CONSTANT` | 918, 919 | 919 | `s11cc2FlatSymbolUsageMaterialAdvectedRho4Constant` | 299098325 / 813 |
| `PY_S11CC2_ROUTING_MECHANICAL_ONLY_OPERAND_MATERIAL_ADVECTED_RHO4_CONSTANT` | 563, 766, 920 | 767 | `s11cc2RoutingMechanicalOnlyOperandMaterialAdvectedRho4Constant` | 299099138 / 3317086 |
| `PY_S11CC2_ROUTING_MECHANICAL_ONLY_RESIDUAL_MATERIAL_ADVECTED_RHO4_CONSTANT` | 768, 920 | 768 | `s11cc2RoutingMechanicalOnlyResidualMaterialAdvectedRho4Constant` | 302416224 / 879140 |
| `PY_S11CC2_ROUTING_TRACTION_CHANNEL_OPERAND_MATERIAL_ADVECTED_RHO4_CONSTANT` | 563, 766, 921 | 767 | `s11cc2RoutingTractionChannelOperandMaterialAdvectedRho4Constant` | 303295364 / 2520370 |
| `PY_S11CC2_ROUTING_TRACTION_CHANNEL_RESIDUAL_MATERIAL_ADVECTED_RHO4_CONSTANT` | 768, 921 | 768 | `s11cc2RoutingTractionChannelResidualMaterialAdvectedRho4Constant` | 305815734 / 1676102 |
| `PY_S11CC2_ZERO_DTN_OPERAND_MATERIAL_ADVECTED_RHO4_CONSTANT` | 563, 766, 922 | 767 | `s11cc2ZeroDtnOperandMaterialAdvectedRho4Constant` | 307491836 / 18456 |
| `PY_S11CC2_ZERO_DTN_RESIDUAL_MATERIAL_ADVECTED_RHO4_CONSTANT` | 768, 922 | 768 | `s11cc2ZeroDtnResidualMaterialAdvectedRho4Constant` | 307510292 / 4177656 |
| `PY_S11CC2_LAMBDA_A_LIMIT_OPERAND_MATERIAL_ADVECTED_RHO4_CONSTANT` | 563, 766, 923 | 767 | `s11cc2LambdaALimitOperandMaterialAdvectedRho4Constant` | 311687948 / 231335 |
| `PY_S11CC2_LAMBDA_A_LIMIT_RESIDUAL_MATERIAL_ADVECTED_RHO4_CONSTANT` | 768, 923 | 768 | `s11cc2LambdaALimitResidualMaterialAdvectedRho4Constant` | 311919283 / 4394335 |
| `PY_S11CC2_IMPERMEABLE_LIMIT_OPERAND_MATERIAL_ADVECTED_RHO4_CONSTANT` | 563, 766, 924 | 767 | `s11cc2ImpermeableLimitOperandMaterialAdvectedRho4Constant` | 316313618 / 182114 |
| `PY_S11CC2_IMPERMEABLE_LIMIT_RESIDUAL_MATERIAL_ADVECTED_RHO4_CONSTANT` | 768, 924 | 768 | `s11cc2ImpermeableLimitResidualMaterialAdvectedRho4Constant` | 316495732 / 4345114 |
| `PY_S11CC2_UNIFORM_LIMIT_OPERAND_MATERIAL_ADVECTED_RHO4_CONSTANT` | 563, 766, 925 | 767 | `s11cc2UniformLimitOperandMaterialAdvectedRho4Constant` | 320840846 / 16664 |
| `PY_S11CC2_UNIFORM_LIMIT_RESIDUAL_MATERIAL_ADVECTED_RHO4_CONSTANT` | 768, 925 | 768 | `s11cc2UniformLimitResidualMaterialAdvectedRho4Constant` | 320857510 / 4179304 |
| `PY_S11CC2_MU_R_FORM_OPERAND_MATERIAL_ADVECTED_RHO4_CONSTANT` | 563, 766, 926 | 767 | `s11cc2MuRFormOperandMaterialAdvectedRho4Constant` | 325036814 / 4917637 |
| `PY_S11CC2_MU_R_FORM_RESIDUAL_MATERIAL_ADVECTED_RHO4_CONSTANT` | 768, 926 | 768 | `s11cc2MuRFormResidualMaterialAdvectedRho4Constant` | 329954451 / 5452871 |
| `PY_S11CC2_CLOSED_SLAB_OPERATOR_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 553, 559 | 893 | `s11cc2ClosedSlabOperatorMaterialAdvectedRhobrConstant` | 335407322 / 3489829 |
| `PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 550, 560 | 893 | `s11cc2ClosedSlabOperatorTermOriginsMaterialAdvectedRhobrConstant` | 338897151 / 3576380 |
| `PY_S11CC2_CLOSED_SLAB_OPERATOR_PARITY_BLOCKS_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 884, 885 | 893 | `s11cc2ClosedSlabOperatorParityBlocksMaterialAdvectedRhobrConstant` | 342473531 / 3303824 |
| `PY_S11CC2_CLOSED_COUPLING_KERNEL_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 561 | 893 | `s11cc2ClosedCouplingKernelMaterialAdvectedRhobrConstant` | 345777355 / 4856918 |
| `PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 887 | 893 | `s11cc2ClosedCouplingKernelTermOriginsMaterialAdvectedRhobrConstant` | 350634273 / 7145727 |
| `PY_S11CC2_SELF_ENERGY_CLOSED_EXTRACT_OPERAND_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 561 | 893 | `s11cc2SelfEnergyClosedExtractOperandMaterialAdvectedRhobrConstant` | 357780000 / 4847349 |
| `PY_S11CC2_SELF_ENERGY_OPEN_EXTRACT_OPERAND_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 562 | 893 | `s11cc2SelfEnergyOpenExtractOperandMaterialAdvectedRhobrConstant` | 362627349 / 1613350 |
| `PY_S11CC2_SELF_ENERGY_INCREMENT_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 563 | 893 | `s11cc2SelfEnergyIncrementMaterialAdvectedRhobrConstant` | 364240699 / 3330109 |
| `PY_S11CC2_FOLD_SYMBOL_MAP_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 397, 484, 499, 501, 509 | 893 | `s11cc2FoldSymbolMapMaterialAdvectedRhobrConstant` | 367570808 / 2265958 |
| `PY_S11CC2_ORDERING_EXTRACT_FIRST_OPERAND_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 901, 904 | 905 | `s11cc2OrderingExtractFirstOperandMaterialAdvectedRhobrConstant` | 369836766 / 3703268 |
| `PY_S11CC2_ORDERING_COMMUTATOR_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 906 | 906 | `s11cc2OrderingCommutatorMaterialAdvectedRhobrConstant` | 373540034 / 1573417 |
| `PY_S11CC2_TRACTION_MECHANICAL_CONTRIB_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 780, 784 | 908 | `s11cc2TractionMechanicalContribMaterialAdvectedRhobrConstant` | 375113451 / 2350450 |
| `PY_S11CC2_TRACTION_SLAB_POWER_PAIRING_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 792, 798, 799, 800, 801 | 909 | `s11cc2TractionSlabPowerPairingMaterialAdvectedRhobrConstant` | 377463901 / 11218386 |
| `PY_S11CC2_TRACTION_SLAB_POWER_PAIRING_RESIDUAL_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 803 | 910 | `s11cc2TractionSlabPowerPairingResidualMaterialAdvectedRhobrConstant` | 388682287 / 9309719 |
| `PY_S11CC2_TRACTION_SIGN_OPERAND_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 779, 792, 801 | 912 | `s11cc2TractionSignOperandMaterialAdvectedRhobrConstant` | 397992006 / 11218497 |
| `PY_S11CC2_TRACTION_SIGN_RESIDUAL_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 913 | 913 | `s11cc2TractionSignResidualMaterialAdvectedRhobrConstant` | 409210503 / 4549679 |
| `PY_S11CC2_DTN_WHOLEFORM_DEPENDENCE_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 915, 916 | 916 | `s11cc2DtnWholeformDependenceMaterialAdvectedRhobrConstant` | 413760182 / 180 |
| `PY_S11CC2_FLAT_SYMBOL_USAGE_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 918, 919 | 919 | `s11cc2FlatSymbolUsageMaterialAdvectedRhobrConstant` | 413760362 / 814 |
| `PY_S11CC2_ROUTING_MECHANICAL_ONLY_OPERAND_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 563, 766, 920 | 767 | `s11cc2RoutingMechanicalOnlyOperandMaterialAdvectedRhobrConstant` | 413761176 / 2657765 |
| `PY_S11CC2_ROUTING_MECHANICAL_ONLY_RESIDUAL_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 768, 920 | 768 | `s11cc2RoutingMechanicalOnlyResidualMaterialAdvectedRhobrConstant` | 416418941 / 664057 |
| `PY_S11CC2_ROUTING_TRACTION_CHANNEL_OPERAND_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 563, 766, 921 | 767 | `s11cc2RoutingTractionChannelOperandMaterialAdvectedRhobrConstant` | 417082998 / 1977211 |
| `PY_S11CC2_ROUTING_TRACTION_CHANNEL_RESIDUAL_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 768, 921 | 768 | `s11cc2RoutingTractionChannelResidualMaterialAdvectedRhobrConstant` | 419060209 / 1344727 |
| `PY_S11CC2_ZERO_DTN_OPERAND_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 563, 766, 922 | 767 | `s11cc2ZeroDtnOperandMaterialAdvectedRhobrConstant` | 420404936 / 18457 |
| `PY_S11CC2_ZERO_DTN_RESIDUAL_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 768, 922 | 768 | `s11cc2ZeroDtnResidualMaterialAdvectedRhobrConstant` | 420423393 / 3303251 |
| `PY_S11CC2_LAMBDA_A_LIMIT_OPERAND_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 563, 766, 923 | 767 | `s11cc2LambdaALimitOperandMaterialAdvectedRhobrConstant` | 423726644 / 231336 |
| `PY_S11CC2_LAMBDA_A_LIMIT_RESIDUAL_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 768, 923 | 768 | `s11cc2LambdaALimitResidualMaterialAdvectedRhobrConstant` | 423957980 / 3519930 |
| `PY_S11CC2_IMPERMEABLE_LIMIT_OPERAND_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 563, 766, 924 | 767 | `s11cc2ImpermeableLimitOperandMaterialAdvectedRhobrConstant` | 427477910 / 182115 |
| `PY_S11CC2_IMPERMEABLE_LIMIT_RESIDUAL_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 768, 924 | 768 | `s11cc2ImpermeableLimitResidualMaterialAdvectedRhobrConstant` | 427660025 / 3470709 |
| `PY_S11CC2_UNIFORM_LIMIT_OPERAND_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 563, 766, 925 | 767 | `s11cc2UniformLimitOperandMaterialAdvectedRhobrConstant` | 431130734 / 16665 |
| `PY_S11CC2_UNIFORM_LIMIT_RESIDUAL_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 768, 925 | 768 | `s11cc2UniformLimitResidualMaterialAdvectedRhobrConstant` | 431147399 / 3304899 |
| `PY_S11CC2_MU_R_FORM_OPERAND_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 563, 766, 926 | 767 | `s11cc2MuRFormOperandMaterialAdvectedRhobrConstant` | 434452298 / 3859984 |
| `PY_S11CC2_MU_R_FORM_RESIDUAL_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 768, 926 | 768 | `s11cc2MuRFormResidualMaterialAdvectedRhobrConstant` | 438312282 / 4015341 |
| `PY_S11CC2_DENSITY_LIVE_MINUS_FROZEN_LAB_HELD` | 930 | 930 | `s11cc2DensityLiveMinusFrozenLabHeld` | 442327623 / 1524219 |
| `PY_S11CC2_DENSITY_LIVE_MINUS_FROZEN_MATERIAL_ADVECTED` | 930 | 930 | `s11cc2DensityLiveMinusFrozenMaterialAdvected` | 443851842 / 1940349 |
| `PY_S11CC2_REP_INVARIANCE_LAB_OPERAND_RHO4_CONSTANT` | 563 | 939 | `s11cc2RepInvarianceLabOperandRho4Constant` | 445792191 / 2565667 |
| `PY_S11CC2_REP_INVARIANCE_MATERIAL_OPERAND_RHO4_CONSTANT` | 934, 935, 936, 937, 938 | 940 | `s11cc2RepInvarianceMaterialOperandRho4Constant` | 448357858 / 4739926 |
| `PY_S11CC2_REP_INVARIANCE_RESIDUAL_RHO4_CONSTANT` | 941 | 941 | `s11cc2RepInvarianceResidualRho4Constant` | 453097784 / 5582453 |
| `PY_S11CC2_REP_INDEPENDENCE_MATERIAL_NORMAL_OPERANDS_RHO4_CONSTANT` | 420, 423, 424 | 943 | `s11cc2RepIndependenceMaterialNormalOperandsRho4Constant` | 458680237 / 1554 |
| `PY_S11CC2_REP_INDEPENDENCE_MATERIAL_OPERAND_RHO4_CONSTANT` | 563, 942, 947 | 944 | `s11cc2RepIndependenceMaterialOperandRho4Constant` | 458681791 / 4780166 |
| `PY_S11CC2_REP_INDEPENDENCE_MATERIAL_RESIDUAL_RHO4_CONSTANT` | 945 | 945 | `s11cc2RepIndependenceMaterialResidualRho4Constant` | 463461957 / 1237382 |
| `PY_S11CC2_REP_INDEPENDENCE_LAB_OPERAND_RHO4_CONSTANT` | 563 | 946 | `s11cc2RepIndependenceLabOperandRho4Constant` | 464699339 / 2565669 |
| `PY_S11CC2_REP_INDEPENDENCE_LAB_NORMAL_OPERANDS_RHO4_CONSTANT` | 420, 423, 424 | 948 | `s11cc2RepIndependenceLabNormalOperandsRho4Constant` | 467265008 / 1549 |
| `PY_S11CC2_REP_INDEPENDENCE_LAB_CORRUPTED_OPERAND_RHO4_CONSTANT` | 563, 947 | 949 | `s11cc2RepIndependenceLabCorruptedOperandRho4Constant` | 467266557 / 3134251 |
| `PY_S11CC2_REP_INDEPENDENCE_LAB_RESIDUAL_RHO4_CONSTANT` | 950 | 950 | `s11cc2RepIndependenceLabResidualRho4Constant` | 470400808 / 1204095 |
| `PY_S11CC2_REP_INDEPENDENCE_MATERIAL_COMPANION_OPERAND_RHO4_CONSTANT` | 563 | 951 | `s11cc2RepIndependenceMaterialCompanionOperandRho4Constant` | 471604903 / 4194963 |
| `PY_S11CC2_REP_INVARIANCE_LAB_OPERAND_RHOBR_CONSTANT` | 563 | 939 | `s11cc2RepInvarianceLabOperandRhobrConstant` | 475799866 / 1920634 |
| `PY_S11CC2_REP_INVARIANCE_MATERIAL_OPERAND_RHOBR_CONSTANT` | 934, 935, 936, 937, 938 | 940 | `s11cc2RepInvarianceMaterialOperandRhobrConstant` | 477720500 / 4231716 |
| `PY_S11CC2_REP_INVARIANCE_RESIDUAL_RHOBR_CONSTANT` | 941 | 941 | `s11cc2RepInvarianceResidualRhobrConstant` | 481952216 / 5653817 |
| `PY_S11CC2_REP_INDEPENDENCE_MATERIAL_NORMAL_OPERANDS_RHOBR_CONSTANT` | 420, 423, 424 | 943 | `s11cc2RepIndependenceMaterialNormalOperandsRhobrConstant` | 487606033 / 1555 |
| `PY_S11CC2_REP_INDEPENDENCE_MATERIAL_OPERAND_RHOBR_CONSTANT` | 563, 942, 947 | 944 | `s11cc2RepIndependenceMaterialOperandRhobrConstant` | 487607588 / 3417286 |
| `PY_S11CC2_REP_INDEPENDENCE_MATERIAL_RESIDUAL_RHOBR_CONSTANT` | 945 | 945 | `s11cc2RepIndependenceMaterialResidualRhobrConstant` | 491024874 / 260560 |
| `PY_S11CC2_REP_INDEPENDENCE_LAB_OPERAND_RHOBR_CONSTANT` | 563 | 946 | `s11cc2RepIndependenceLabOperandRhobrConstant` | 491285434 / 1920636 |
| `PY_S11CC2_REP_INDEPENDENCE_LAB_NORMAL_OPERANDS_RHOBR_CONSTANT` | 420, 423, 424 | 948 | `s11cc2RepIndependenceLabNormalOperandsRhobrConstant` | 493206070 / 1550 |
| `PY_S11CC2_REP_INDEPENDENCE_LAB_CORRUPTED_OPERAND_RHOBR_CONSTANT` | 563, 947 | 949 | `s11cc2RepIndependenceLabCorruptedOperandRhobrConstant` | 493207620 / 2017410 |
| `PY_S11CC2_REP_INDEPENDENCE_LAB_RESIDUAL_RHOBR_CONSTANT` | 950 | 950 | `s11cc2RepIndependenceLabResidualRhobrConstant` | 495225030 / 260555 |
| `PY_S11CC2_REP_INDEPENDENCE_MATERIAL_COMPANION_OPERAND_RHOBR_CONSTANT` | 563 | 951 | `s11cc2RepIndependenceMaterialCompanionOperandRhobrConstant` | 495485585 / 3320532 |
| `PY_S11CC2_EXPORT_ROUNDTRIP` | 848 | 849 | `s11cc2ExportRoundtrip` | 498806117 / 5288 |
