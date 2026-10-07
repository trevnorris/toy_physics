# Measurements — O2 comparator build r0 review dispositions (generated 2026-10-07 15:03)

Generator: `_scratch/s9b_build/gen/o2_comparator_build_r0_lookups.sh` (sed/grep/sha256sum only; lines cut at 300 characters).

```
$ (cd /var/projects/toy_physics && sha256sum -c _scratch/s9b_build/o2_comparator_build_review_baseline_r0.sha256)
research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py: OK
research/pde_ledger_v3/scripts/test_O2_cross_engine_comparator.py: OK
_scratch/s11c/o2-comparator-build/comparison.jsonl: OK
_scratch/s11c/o2-comparator-build/accounting.jsonl: OK
_scratch/s11c/o2-comparator-build/catalog.json: OK
_scratch/s11c/o2-comparator-build/summary.json: OK
research/pde_ledger_v3/_measurements/O2_comparator_builder_report.md: OK
research/pde_ledger_v3/_measurements/O2_comparator_builder_measurements.txt: OK
```

## R1 — the mass_loss_orientation join and the two engine lines it cites
```
$ sed -n '424p' scripts/O2_cross_engine_comparator.py
    J('mass_loss_orientation','MASS_INPUT/outward_loss','MASS_INPUT/RHS',378,244,'§3.1,5'),
```

```
$ sed -n '266p;377,378p' scripts/O2_live_balance_sympy_audit.py
    mass_equation = sp.Eq(mass_divergence, -jn, evaluate=False)
        'MASS_INPUT': record(current=mass_current, divergence=mass_divergence,
                             outward_loss=jn, supplied_equation=mass_equation,
```

```
$ sed -n '244,247p' mathematica/O2_live_balance_mathematica_audit.wl
massRHS = -profiles["j_n"];
emit["MASS_INPUT", <|"DensityOperand" -> massDensity, "VelocityOperand" -> vPlane,
  "Flux" -> massFlux, "Divergence" -> massDivergence, "RHS" -> massRHS,
  "Equation" -> (massDivergence == massRHS), "Residual" -> (massDivergence - massRHS),
```

```
$ grep -n 'mass_equation\|profile_j_n' scripts/O2_cross_engine_comparator.py | head -8
384:    J('mass_equation','MASS_INPUT/supplied_equation','MASS_INPUT/Equation',266,247,'§3.1'),
```

## R2 — the signed-action census and the face-support entry in both engines
```
$ sed -n '672,691p' scripts/O2_cross_engine_comparator.py
    def signed_actions(v,sign=1):
        if v.head == 'Mul':
            for child in v.args:
                if child.head == 'Number' and child.value.startswith('-'):
                    sign *= -1
            for child in v.args:
                if child.head != 'Number':
                    signed_actions(child,sign)
        elif v.head == 'Apply':
            head = v.args[0]
            if head.value.startswith('py::OPEN_'):
                orientations[(head.value,sign)] += 1
            elif head.value == 'wl::OpenAction' and len(v.args) > 1:
                orientations[(json.dumps(data(v.args[1]),separators=(',',':')),sign)] += 1
            # An action's arguments are dependencies, not signed balance terms.
            for child in v.args[1:]:
                signed_actions(child,1)
        else:
            for child in v.args:
                signed_actions(child,sign if v.head in ('Add','Derivative') else 1)
```

```
$ sed -n '259,260p' scripts/O2_live_balance_sympy_audit.py
    entries = ((1, storage), (1, transport), (1, carry), (1, partners),
               (-1, internal), (-1, reduced_load)) + body_entries
```

```
$ grep -n 'MechanicalFaceSupport' mathematica/O2_live_balance_mathematica_audit.wl
231:  entry[MechanicalFaceSupport,-1,mechanicalLoad,origin["4,5","7"]],
```

## R3 — what the controls compare
```
$ sed -n '30,40p' scripts/test_O2_cross_engine_comparator.py
            c.run(p,w,output,accounting,joins=(row,),names=names)
        parsed = [json.loads(line) for line in output.getvalue().splitlines()]
        residual = [line['comparison'] for line in parsed if line['kind']=='residual']
        return residual,parsed,[json.loads(line) for line in accounting.getvalue().splitlines()]

    def changes(self,py,wl,mutation):
        before = self.fixture(py,wl)[0]
        after = self.fixture(py,mutation)[0]
        self.assertNotEqual(before,after,'mutation must move the printed residual, not just operand text')

    def test_one_sided_operand(self):
```

```
$ grep -c 'structure' scripts/test_O2_cross_engine_comparator.py
1
```

```
$ grep -n "'py_value'\|'wl_value'" scripts/O2_cross_engine_comparator.py | head -6
756:    return {'outcome':'exact','py_value':str(left),'wl_value':str(right),
```

```
$ grep -n 'def test' scripts/test_O2_cross_engine_comparator.py
40:    def test_one_sided_operand(self):
43:    def test_form_change(self):
46:    def test_each_name_binding_repoint(self):
67:    def test_applied_argument_stripped(self):
70:    def test_live_profile_frozen(self):
73:    def test_open_head(self):
77:    def test_open_named_operand(self):
81:    def test_open_argument(self):
85:    def test_open_orientation(self):
89:    def test_derivative_order(self):
93:    def test_binder_structure(self):
97:    def test_nested_sibling_removed(self):
101:    def test_moved_tag_and_key_stays_joined(self):
109:    def test_duplicate_join_rejected(self):
114:    def test_parent_child_join_rejected(self):
119:    def test_duplicate_name_rejected_both_sides(self):
125:    def test_boolean_does_not_hide_algebraic_sibling(self):
134:    def test_lossless_function_metadata_survives(self):
141:    def test_empty_containers_are_accounted(self):
146:    def test_matching_text_and_names_not_zero_evidence(self):
152:    def test_grammatical_unknown_head_is_coverage(self):
156:    def test_malformed_input_rejected(self):
162:    def test_declared_tables_are_injective(self):
165:    def test_independent_census_accounts_for_nested_sibling(self):
174:    def test_profile_derivative_reaches_scalar_subtraction(self):
180:    def test_unpaired_open_subtree_mutation_moves_residual(self):
```

## R4 — the four OPEN operands: join rows, unjoined reasons, name table
```
$ sed -n '430,431p;804p;808p' scripts/O2_cross_engine_comparator.py
    J('material_identifications','COUPLED_INPUTS/operands/29','MATERIAL_MOMENTUM/DifferentiatedSection/OtherDependence/@6',162,133,'§3.2'),
    J('energy_overlap','COUPLED_INPUTS/operands/30','B_E_STEADY/Entries/2/Object/@2/1',162,289,'§6'),
        'COUPLED_INPUTS/operands/24':'PY separately registers the S12 reaction system; WL carries it within additional momentum-partner actions, with no second register occurrence (§5).',
        'COUPLED_INPUTS/operands/31':'PY separately registers face-support partition; WL carries unresolved support partition within complete-hold actions, no separate register occurrence (§3.3).',
```

```
$ grep -n 'material_action_compatibility\|UnresolvedStressInertiaNormalIdentifications\|energy_accounting_overlap\|UnresolvedEnergyOccurrenceIdentifications\|S12_reaction_system\|OPENReactionSystem\|face_support_partition\|UnresolvedSupportPartition' scripts/O2_cross_engine_comparator.py
```

```
$ grep -c 'material_action_compatibility' scripts/O2_live_balance_sympy_audit.py; grep -c 'energy_accounting_overlap' scripts/O2_live_balance_sympy_audit.py; grep -c 'S12_reaction_system' scripts/O2_live_balance_sympy_audit.py; grep -c 'face_support_partition' scripts/O2_live_balance_sympy_audit.py
2
2
2
2
```

```
$ grep -c 'UnresolvedStressInertiaNormalIdentifications' mathematica/O2_live_balance_mathematica_audit.wl; grep -c 'UnresolvedEnergyOccurrenceIdentifications' mathematica/O2_live_balance_mathematica_audit.wl; grep -c 'OPENReactionSystem' mathematica/O2_live_balance_mathematica_audit.wl; grep -c 'UnresolvedSupportPartition' mathematica/O2_live_balance_mathematica_audit.wl
4
1
1
1
```

