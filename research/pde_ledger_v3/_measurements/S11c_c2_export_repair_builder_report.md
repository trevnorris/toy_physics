# S11c-c2 export repair — builder report

Build → verify → report completed. Publication-only; physics construction and emissions are unchanged.

The generated `scripts/S11c_c2_exports.py` is **22,441,522 bytes**, down from **60,516,900 bytes**: **38,075,378 bytes removed (62.92%; 2.70× smaller)**.

## Exact edit sites

Only `scripts/S11c_c2_selfenergy_fold_sympy_audit.py` was edited as source:

- `EXPORT_ROOTS`, current line 48: removed `s11cc2SelfEnergyIncrement`; retained both declared deliverables.
- New adjacent publication helper `publication_compact`, current line 806.
- `publish`, current line 852: compact only case `VALUE` trees; bounded root-name displays; retain declaration schema, F9, closure, minimality, digest structure, and structural roundtrip; add strict emitted-versus-restored semantics and reciprocal guards, literal evidence and size measurements; install only after passing.
- The sole `export_key` dictionary inside `run`, current line 1051: removed increment export routing. The increment's existing emit remains.

AST comparison against the pre-edit source passed after excluding only these authorized sites. All other construction, extraction, control, grading, emission, and physics-loop code is identical. No construction change was needed.

## Representation and size

Transform: exact `collect` followed by `factor_terms(fraction=False)` on per-case VALUE expressions, with reciprocal/calculus atoms protected locally and restored immediately; Integral boundaries and container types retained; any subexpression whose reciprocal-power set would change keeps its original form. No expansion is used to produce export values, no CSE is used, and no temporary/hold symbol survives in the delivered values. Expansion occurs only in the separate semantic guard.

`s11cc2ClosedCouplingKernel` VALUE totals: 18,599,401 → 13,917,507 bytes (25.17% smaller). `s11cc2ClosedSlabOperator` VALUE totals: 15,031,549 → 8,422,000 bytes (43.97% smaller). `s11cc2SelfEnergyIncrement` VALUE totals: 11,997,906 → 0 bytes (100.00% smaller; absent from delta, still emitted).

Sizes below are UTF-8 byte lengths of each case's `srepr(VALUE)`, excluding its unchanged metadata. Baseline sizes were measured from the original 60,516,900-byte module; the final build's emitted sizes match those baseline sizes exactly. Zero means absent from the delta.

| Object | Case | Before bytes | After bytes | Reduction |
| --- | --- | ---: | ---: | ---: |
| `s11cc2ClosedCouplingKernel` | `LAB_HELD/RHO4_CONSTANT` | 4,355,993 | 3,375,460 | 22.51% |
| `s11cc2ClosedCouplingKernel` | `LAB_HELD/RHOBR_CONSTANT` | 3,407,032 | 2,596,031 | 23.80% |
| `s11cc2ClosedCouplingKernel` | `MATERIAL_ADVECTED/RHO4_CONSTANT` | 5,990,000 | 4,422,717 | 26.16% |
| `s11cc2ClosedCouplingKernel` | `MATERIAL_ADVECTED/RHOBR_CONSTANT` | 4,846,376 | 3,523,299 | 27.30% |
| `s11cc2ClosedSlabOperator` | `LAB_HELD/RHO4_CONSTANT` | 4,002,816 | 2,271,711 | 43.25% |
| `s11cc2ClosedSlabOperator` | `LAB_HELD/RHOBR_CONSTANT` | 3,121,562 | 1,718,913 | 44.93% |
| `s11cc2ClosedSlabOperator` | `MATERIAL_ADVECTED/RHO4_CONSTANT` | 4,427,664 | 2,502,613 | 43.48% |
| `s11cc2ClosedSlabOperator` | `MATERIAL_ADVECTED/RHOBR_CONSTANT` | 3,479,507 | 1,928,763 | 44.57% |
| `s11cc2SelfEnergyIncrement` | `LAB_HELD/RHO4_CONSTANT` | 2,564,686 | 0 | 100.00% |
| `s11cc2SelfEnergyIncrement` | `LAB_HELD/RHOBR_CONSTANT` | 1,919,652 | 0 | 100.00% |
| `s11cc2SelfEnergyIncrement` | `MATERIAL_ADVECTED/RHO4_CONSTANT` | 4,194,000 | 0 | 100.00% |
| `s11cc2SelfEnergyIncrement` | `MATERIAL_ADVECTED/RHOBR_CONSTANT` | 3,319,568 | 0 | 100.00% |

The two operator displays are exactly their root names (24 and 26 characters). The row field sets and all retained declaration fields/values/displays match the original module by AST comparison. `value_kind`, `class`, `step`, `route`, and declaration `dimension_key` are preserved.

## Build and deliverable verification

Final full command:

```sh
python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py > /tmp/S11c_c2_export_repair.out 2>&1
```

Exit status 0; elapsed 1447.073 seconds; peak RSS 1,917,760 KiB (1.829 GiB). Only one heavy CAS process ran at a time; no timeout or kill was used. The unchanged harness wrote `_measurements/S11c_c2_sympy_guard_evidence.json` at its checkpoints and `_measurements/S11c_c2_sympy_progress.json` from the existing `build_case` writes (current lines 533–564). Both are run side effects, not builder-authored review artifacts.

The first full attempt stopped at the new pole guard before replacing the target: collecting a kernel coefficient to zero erased reciprocal expressions. The helper was corrected within publication to retain any affected original subexpression. The same `publish` function then passed against all eight emitted payloads, followed by the final full run above. All literal results below are from that final full run.

Normal Python import of the delivered module passed. Its `_restore(srepr(value))` roundtrip passed for all 70 rows. The restored delta folded over the two parents passes `check_consumer` and `assert_delta_is_minimal`, with no F9 collisions. The five `BUILD_INPUT_DIGESTS` entries retain their schema and all match current files, including the edited audit script.

Both roots bind through `fold[key]['value']` → `cases(...)` to exactly all four cases, with the five named payload fields preserved. A direct consumer API probe differentiated all 44 scalar VALUE leaves with respect to their present `epsilon_shape` symbol (literal zero channels handled directly); all returned valid SymPy expressions with no Dummy or hold wrapper. This was a binding/differentiation check, not a downstream physics step.

Tag/emit-key identity: **154/154 names match** the prior stdout/provenance index, with no duplicates, additions, or omissions. This includes the retained publication roundtrip tag and all **four SELF_ENERGY_INCREMENT emits**. No byte comparison of the large physics payloads was used. New publication evidence uses separate `EXPORT_*` diagnostic lines.

## R1–R4 confirmation

- **R1 PASS:** Exactly 70 own rows: the two operator roots plus their new-symbol/dimension declaration closure over the two-parent fold. Increment, origins, parity, §3d, §5 controls, and operand rows are absent. No concrete S11c-d manifest is assumed; both declared prospective deliverables are retained.
- **R2 PASS:** Cased nested payloads and metadata survive. Structural roundtrip is separate and retained. The added generated-module semantic guard checks case sets, mapping keys/uniqueness, pair and tuple arities, labels and matrix shapes, exact metadata, and Integral-aware scalar expansion. All **44/44 expanded differences are literal 0**; all **44/44 reciprocal-power sets match exactly**, including bases and exponents. No singular factor is cancelled in the delivered compaction.
- **R3 PASS:** Bounded operator displays replace the oversized duplicate; other row fields and declarations are unchanged.
- **R4 PASS:** Generated module imports; all 70 rows roundtrip; binding and differentiation pass; restored closure/minimality/F9 pass; all five build digests match.

## Literal guard outputs from the final full build

The retained structural roundtrip output (70 true entries):

```text
PY_S11CC2_EXPORT_ROUNDTRIP: Tuple(Tuple(Str('VALUE'), Tuple(true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true)), Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)))), Tuple(Str('DIMENSION_L_T_M'), Tuple(ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]), ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)]]))))
```

The separate semantic, pole, size, closure, minimality, and collision outputs are reproduced verbatim below. Every container/metadata/equality predicate is True and every expanded scalar difference is 0. These checks ran against `namespace['LEDGER']` produced by executing the actual generated module source, before its atomic installation.

```text
EXPORT_F9_COLLISIONS = []
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel case_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT payload_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE mapping_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS mapping_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/THETA algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/THETA reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/THETA expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/THETA expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/E_W algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/E_W reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/E_W expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/E_W expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/DIV_U algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/DIV_U reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/DIV_U expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/DIV_U expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE mapping_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/THETA algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/THETA reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/THETA expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/THETA expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/E_W algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/E_W reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/E_W expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/E_W expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/DIV_U algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/DIV_U reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/DIV_U expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/DIV_U expanded_difference_is_zero = True
EXPORT_VALUE_BYTES s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT = {'emitted_srepr_bytes': 4355993, 'compact_srepr_bytes': 3375460}
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/MULTIGRADE metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/DIMENSION_L_T_M metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/COMPUTED_BRANCH_BINDINGS metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/FOURIER_PROFILE_BINDINGS metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT payload_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE mapping_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS mapping_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/THETA algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/THETA reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/THETA expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/THETA expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/E_W algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/E_W reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/E_W expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/E_W expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/DIV_U algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/DIV_U reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/DIV_U expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/DIV_U expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE mapping_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/THETA algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/THETA reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/THETA expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/THETA expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/E_W algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/E_W reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/E_W expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/E_W expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/DIV_U algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/DIV_U reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/DIV_U expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/DIV_U expanded_difference_is_zero = True
EXPORT_VALUE_BYTES s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT = {'emitted_srepr_bytes': 3407032, 'compact_srepr_bytes': 2596031}
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/MULTIGRADE metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/DIMENSION_L_T_M metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/COMPUTED_BRANCH_BINDINGS metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/LAB_HELD/RHOBR_CONSTANT/FOURIER_PROFILE_BINDINGS metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT payload_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE mapping_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS mapping_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/THETA algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/THETA reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/THETA expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/THETA expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/E_W algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/E_W reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/E_W expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/E_W expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/DIV_U algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/DIV_U reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/DIV_U expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/DIV_U expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE mapping_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/THETA algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/THETA reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/THETA expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/THETA expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/E_W algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/E_W reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/E_W expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/E_W expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/DIV_U algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/DIV_U reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/DIV_U expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/DIV_U expanded_difference_is_zero = True
EXPORT_VALUE_BYTES s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT = {'emitted_srepr_bytes': 5990000, 'compact_srepr_bytes': 4422717}
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/MULTIGRADE metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/DIMENSION_L_T_M metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/COMPUTED_BRANCH_BINDINGS metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHO4_CONSTANT/FOURIER_PROFILE_BINDINGS metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT payload_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE mapping_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS mapping_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/THETA algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/THETA reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/THETA expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/THETA expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/E_W algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/E_W reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/E_W expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/E_W expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/DIV_U algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/DIV_U reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/DIV_U expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/DIV_U expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE mapping_keys = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/THETA algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/THETA reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/THETA expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/THETA expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/E_W algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/E_W reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/E_W expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/E_W expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/DIV_U algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/DIV_U reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/DIV_U expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THICKNESS_TO_TRANSVERSE/DIV_U expanded_difference_is_zero = True
EXPORT_VALUE_BYTES s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT = {'emitted_srepr_bytes': 4846376, 'compact_srepr_bytes': 3523299}
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/MULTIGRADE metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/DIMENSION_L_T_M metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/COMPUTED_BRANCH_BINDINGS metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedCouplingKernel/MATERIAL_ADVECTED/RHOBR_CONSTANT/FOURIER_PROFILE_BINDINGS metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator case_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT payload_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE mapping_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/U tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/U tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/U[0] algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/U[0] reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/U[0] expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/U[0] expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/U[1] algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/U[1] reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/U[1] expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/U[1] expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/U[2] algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/U[2] reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/U[2] expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/U[2] expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/THETA algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/THETA reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/THETA expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/THETA expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/E_W algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/E_W reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/E_W expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/E_W expanded_difference_is_zero = True
EXPORT_VALUE_BYTES s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT = {'emitted_srepr_bytes': 4002816, 'compact_srepr_bytes': 2271711}
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/MULTIGRADE metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/DIMENSION_L_T_M metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/COMPUTED_BRANCH_BINDINGS metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/FOURIER_PROFILE_BINDINGS metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT payload_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE mapping_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/U tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/U tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/U[0] algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/U[0] reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/U[0] expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/U[0] expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/U[1] algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/U[1] reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/U[1] expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/U[1] expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/U[2] algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/U[2] reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/U[2] expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/U[2] expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/THETA algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/THETA reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/THETA expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/THETA expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/E_W algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/E_W reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/E_W expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/VALUE/E_W expanded_difference_is_zero = True
EXPORT_VALUE_BYTES s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT = {'emitted_srepr_bytes': 3121562, 'compact_srepr_bytes': 1718913}
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/MULTIGRADE metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/DIMENSION_L_T_M metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/COMPUTED_BRANCH_BINDINGS metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/LAB_HELD/RHOBR_CONSTANT/FOURIER_PROFILE_BINDINGS metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT payload_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE mapping_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/U tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/U tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/U[0] algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/U[0] reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/U[0] expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/U[0] expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/U[1] algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/U[1] reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/U[1] expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/U[1] expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/U[2] algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/U[2] reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/U[2] expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/U[2] expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THETA algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THETA reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THETA expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/THETA expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/E_W algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/E_W reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/E_W expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/VALUE/E_W expanded_difference_is_zero = True
EXPORT_VALUE_BYTES s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT = {'emitted_srepr_bytes': 4427664, 'compact_srepr_bytes': 2502613}
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/MULTIGRADE metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/DIMENSION_L_T_M metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/COMPUTED_BRANCH_BINDINGS metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHO4_CONSTANT/FOURIER_PROFILE_BINDINGS metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT payload_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/emitted mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/emitted pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/emitted unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/restored mapping_tuple = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/restored pair_arities = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/restored unique_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE mapping_keys = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/U tuple_type = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/U tuple_arity = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/U[0] algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/U[0] reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/U[0] expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/U[0] expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/U[1] algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/U[1] reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/U[1] expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/U[1] expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/U[2] algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/U[2] reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/U[2] expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/U[2] expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THETA algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THETA reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THETA expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/THETA expanded_difference_is_zero = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/E_W algebraic_leaf = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/E_W reciprocal_powers_unchanged = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/E_W expanded_difference = 0
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/VALUE/E_W expanded_difference_is_zero = True
EXPORT_VALUE_BYTES s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT = {'emitted_srepr_bytes': 3479507, 'compact_srepr_bytes': 1928763}
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/MULTIGRADE metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/DIMENSION_L_T_M metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/COMPUTED_BRANCH_BINDINGS metadata_exact = True
EXPORT_SEMANTIC s11cc2ClosedSlabOperator/MATERIAL_ADVECTED/RHOBR_CONSTANT/FOURIER_PROFILE_BINDINGS metadata_exact = True
EXPORT_RESTORED_CLOSURE = ['s11cc2ClosedCouplingKernel', 's11cc2ClosedSlabOperator', 's11cc2Coefficientm1Profile', 's11cc2Coefficientm1ProfileDimension', 's11cc2Coefficientw1Profile', 's11cc2Coefficientw1ProfileDimension', 's11cc2FieldeW', 's11cc2FieldeWDimension', 's11cc2Fieldtheta', 's11cc2FieldthetaDimension', 's11cc2Fieldu1', 's11cc2Fieldu1Dimension', 's11cc2Fieldu2', 's11cc2Fieldu2Dimension', 's11cc2Fieldu3', 's11cc2Fieldu3Dimension', 's11cc2FourierW1ProfileHatTransfer', 's11cc2FourierW1ProfileHatTransferDimension', 's11cc2FourierW1ProfileJetHat1', 's11cc2FourierW1ProfileJetHat1Dimension', 's11cc2FourierW1ProfileJetHat2', 's11cc2FourierW1ProfileJetHat2Dimension', 's11cc2FourierW1ProfileJetHat3', 's11cc2FourierW1ProfileJetHat3Dimension', 's11cc2MiddleMomentum1', 's11cc2MiddleMomentum1Dimension', 's11cc2MiddleMomentum2', 's11cc2MiddleMomentum2Dimension', 's11cc2MiddleMomentum3', 's11cc2MiddleMomentum3Dimension', 's11cc2OutgoingNormalMomentum', 's11cc2OutgoingNormalMomentumDimension', 's11cc2TestA0', 's11cc2TestA0Dimension', 's11cc2TestA1', 's11cc2TestA1Dimension', 's11cc2TestA2', 's11cc2TestA2Dimension', 's11cc2TestE', 's11cc2TestEDimension', 's11cc2TestPhi', 's11cc2TestPhiDimension', 's11cc2TestTheta', 's11cc2TestThetaDimension', 's11cc2Time', 's11cc2TimeDimension', 's11cc2TrialA0', 's11cc2TrialA0Dimension', 's11cc2TrialA1', 's11cc2TrialA1Dimension', 's11cc2TrialA2', 's11cc2TrialA2Dimension', 's11cc2TrialE', 's11cc2TrialEDimension', 's11cc2TrialPhi', 's11cc2TrialPhiDimension', 's11cc2TrialTheta', 's11cc2TrialThetaDimension', 's11cc2X1', 's11cc2X1Dimension', 's11cc2X2', 's11cc2X2Dimension', 's11cc2X3', 's11cc2X3Dimension', 's11cc2Y1', 's11cc2Y1Dimension', 's11cc2Y2', 's11cc2Y2Dimension', 's11cc2Y3', 's11cc2Y3Dimension']
EXPORT_RESTORED_MINIMALITY = {'exported_keys': ['s11cc2ClosedCouplingKernel', 's11cc2ClosedSlabOperator', 's11cc2Coefficientm1Profile', 's11cc2Coefficientm1ProfileDimension', 's11cc2Coefficientw1Profile', 's11cc2Coefficientw1ProfileDimension', 's11cc2FieldeW', 's11cc2FieldeWDimension', 's11cc2Fieldtheta', 's11cc2FieldthetaDimension', 's11cc2Fieldu1', 's11cc2Fieldu1Dimension', 's11cc2Fieldu2', 's11cc2Fieldu2Dimension', 's11cc2Fieldu3', 's11cc2Fieldu3Dimension', 's11cc2FourierW1ProfileHatTransfer', 's11cc2FourierW1ProfileHatTransferDimension', 's11cc2FourierW1ProfileJetHat1', 's11cc2FourierW1ProfileJetHat1Dimension', 's11cc2FourierW1ProfileJetHat2', 's11cc2FourierW1ProfileJetHat2Dimension', 's11cc2FourierW1ProfileJetHat3', 's11cc2FourierW1ProfileJetHat3Dimension', 's11cc2MiddleMomentum1', 's11cc2MiddleMomentum1Dimension', 's11cc2MiddleMomentum2', 's11cc2MiddleMomentum2Dimension', 's11cc2MiddleMomentum3', 's11cc2MiddleMomentum3Dimension', 's11cc2OutgoingNormalMomentum', 's11cc2OutgoingNormalMomentumDimension', 's11cc2TestA0', 's11cc2TestA0Dimension', 's11cc2TestA1', 's11cc2TestA1Dimension', 's11cc2TestA2', 's11cc2TestA2Dimension', 's11cc2TestE', 's11cc2TestEDimension', 's11cc2TestPhi', 's11cc2TestPhiDimension', 's11cc2TestTheta', 's11cc2TestThetaDimension', 's11cc2Time', 's11cc2TimeDimension', 's11cc2TrialA0', 's11cc2TrialA0Dimension', 's11cc2TrialA1', 's11cc2TrialA1Dimension', 's11cc2TrialA2', 's11cc2TrialA2Dimension', 's11cc2TrialE', 's11cc2TrialEDimension', 's11cc2TrialPhi', 's11cc2TrialPhiDimension', 's11cc2TrialTheta', 's11cc2TrialThetaDimension', 's11cc2X1', 's11cc2X1Dimension', 's11cc2X2', 's11cc2X2Dimension', 's11cc2X3', 's11cc2X3Dimension', 's11cc2Y1', 's11cc2Y1Dimension', 's11cc2Y2', 's11cc2Y2Dimension', 's11cc2Y3', 's11cc2Y3Dimension'], 'required_keys': ['s11cc2ClosedCouplingKernel', 's11cc2ClosedSlabOperator', 's11cc2Coefficientm1Profile', 's11cc2Coefficientm1ProfileDimension', 's11cc2Coefficientw1Profile', 's11cc2Coefficientw1ProfileDimension', 's11cc2FieldeW', 's11cc2FieldeWDimension', 's11cc2Fieldtheta', 's11cc2FieldthetaDimension', 's11cc2Fieldu1', 's11cc2Fieldu1Dimension', 's11cc2Fieldu2', 's11cc2Fieldu2Dimension', 's11cc2Fieldu3', 's11cc2Fieldu3Dimension', 's11cc2FourierW1ProfileHatTransfer', 's11cc2FourierW1ProfileHatTransferDimension', 's11cc2FourierW1ProfileJetHat1', 's11cc2FourierW1ProfileJetHat1Dimension', 's11cc2FourierW1ProfileJetHat2', 's11cc2FourierW1ProfileJetHat2Dimension', 's11cc2FourierW1ProfileJetHat3', 's11cc2FourierW1ProfileJetHat3Dimension', 's11cc2MiddleMomentum1', 's11cc2MiddleMomentum1Dimension', 's11cc2MiddleMomentum2', 's11cc2MiddleMomentum2Dimension', 's11cc2MiddleMomentum3', 's11cc2MiddleMomentum3Dimension', 's11cc2OutgoingNormalMomentum', 's11cc2OutgoingNormalMomentumDimension', 's11cc2TestA0', 's11cc2TestA0Dimension', 's11cc2TestA1', 's11cc2TestA1Dimension', 's11cc2TestA2', 's11cc2TestA2Dimension', 's11cc2TestE', 's11cc2TestEDimension', 's11cc2TestPhi', 's11cc2TestPhiDimension', 's11cc2TestTheta', 's11cc2TestThetaDimension', 's11cc2Time', 's11cc2TimeDimension', 's11cc2TrialA0', 's11cc2TrialA0Dimension', 's11cc2TrialA1', 's11cc2TrialA1Dimension', 's11cc2TrialA2', 's11cc2TrialA2Dimension', 's11cc2TrialE', 's11cc2TrialEDimension', 's11cc2TrialPhi', 's11cc2TrialPhiDimension', 's11cc2TrialTheta', 's11cc2TrialThetaDimension', 's11cc2X1', 's11cc2X1Dimension', 's11cc2X2', 's11cc2X2Dimension', 's11cc2X3', 's11cc2X3Dimension', 's11cc2Y1', 's11cc2Y1Dimension', 's11cc2Y2', 's11cc2Y2Dimension', 's11cc2Y3', 's11cc2Y3Dimension'], 'allowed_infra_keys': []}
EXPORT_SEMANTIC delta closure_unchanged = True
```

Builder task ends here. No review, independent comparator, or downstream step was launched; no self-review, derived-or-declared, output-checker, or finalize artifact was created, and `.claude/skills/*` was not read.
