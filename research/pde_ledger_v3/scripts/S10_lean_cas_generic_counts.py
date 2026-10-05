"""Bind the 42 generic D3 integer emissions to full matrices on their chart.

The N1/N3 matrices are assembled from the existing imported scalar trees.
All seven counts have explicit semantic obligations; their generic chart is
never extended to the parallel or perpendicular exceptional strata.
"""
import S10_lean_cas_counts as counts


def fields(engine, root):
    name = f'root{root}'
    M, S, B = [f'({engine}GenericCountReference.{name}{label} rho mu sigma z k)'
               for label in ('Matrix', 'Stack', 'Basis')]
    transverse = int(root == 1)
    return [('N2_RANK', 'rank', 2, f'{M}.rank', 'ℕ'),
            ('N2_NULLITY', 'nullity', 1, f'matrixNullity {M}', 'ℕ'),
            ('N3_STACKED_RANK', 'stacked_rank', 3-transverse, f'{S}.rank', 'ℕ'),
            ('N3_TRANSVERSE_NULLITY', 'transverse_nullity', transverse, f'matrixNullity {S}', 'ℕ'),
            ('N4_NULLITY_DIFFERENCE', 'difference', 1-transverse, f'nullityDifference {M} {S}', 'ℤ'),
            ('N7_BASIS_COUNT', 'basis_count', 1, f'basisCount {B}', 'ℕ'),
            ('N7_BASIS_COUNT_RESIDUAL' if engine == 'PY' else 'N7_COUNT_RESIDUAL',
             'basis_residual', 0, f'basisCountResidual {B} {M}', 'ℤ')]


def build(bridge, engine, path):
    b = bridge.Builder(engine, path, namespace=engine+'GenericCounts',
                       support='GenericCountReference', tactic='count_eval')
    for root in range(3):
        for tag, label, expected, semantic, number_type in fields(engine, root):
            b.record(f'ROOT{root+1}_'+tag, lambda x: counts.count_shape(bridge, x),
                     [(str(expected), (0,0,0), [])])
            b.records[-1]['count_semantics'] = {
                'quantity': label, 'case': f'root{root}', 'number_type': number_type,
                'reference': semantic,
                'lean_semantic_claim': engine+f'GenericCounts.root{root}_'+label,
                'coefficient_domain': ['rho != 0', 'mu != 0', '0 < sigma', 'sigma != 1'],
                'coordinate_domain': ['k0 != 0', 'k1 != 0', 'k2 != 0'],
                'lean_chart': 'GenericChart sigma k',
                'subtraction': 'signed integer' if number_type == 'ℤ' else None}
    return b


def reference_module(builders):
    lines = ['import S10Audit.CAS.GenericCountSupport', '',
             'set_option backward.isDefEq.respectTransparency false', '',
             'namespace S10Audit.CAS', 'open S10Pilot', 'noncomputable section', '']
    audits = ['genericChart_mode_dimensions']
    params = '(rho mu sigma z : ℝ) (k : Vec 3)'
    args = params+' (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k)'
    val, call = 'rho mu sigma z k', 'rho mu sigma z k hr hm h'
    for engine, b in builders.items():
        ns = engine+'GenericCountReference'
        lines += ['namespace '+ns, '']
        for root in range(3):
            name, transverse = f'root{root}', int(root == 1)
            M, S, B = [f'({name}{label} {val})' for label in ('Matrix', 'Stack', 'Basis')]
            target = f'((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma k {root}) k)'
            for label, suffix, rows in [('Matrix', 'N1_MATRIX', 3), ('Stack', 'N3_STACKED_MATRIX', 4)]:
                record = next(rec for rec in b.records if rec['tag'].endswith(f'_ROOT{root+1}_'+suffix))
                names = [engine+'.'+c['lean_tree'] for c in record['cells']]
                tree = '!['+', '.join('!['+', '.join(names[3*i:3*i+3])+']' for i in range(rows))+']'
                def cell_proof(cell):
                    factors = cell['raw_denominator_factors'] + [['atom', x] for x in cell['reference_domain_symbols']]
                    assert all(d[0] == 'atom' and d[1] in ('rho', 'sigma') for d in factors)
                    symbols = dict.fromkeys(d[1] for d in factors)
                    extra = ''.join(' _hr' if x == 'rho' else ' _hs' for x in symbols)
                    return 'exact '+engine+'.'+cell['lean_claim']+' '+val+extra
                ref = target if label == 'Matrix' else 'appendConstraint '+target+' k'
                lines += [f'def {name}{label} {params} : Matrix (Fin {rows}) (Fin 3) ℝ :=',
                          f'  fun i j => Expr.eval (values {val})',
                          f'    (({tree} : Fin {rows} → Fin 3 → Expr Symbol) i j)', '',
                          f'theorem {name}{label}_reference {params} (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :',
                          f'    {name}{label} {val} = {ref} := by',
                          '  ext i j', '  fin_cases i', '  all_goals fin_cases j',
                          *['  · '+cell_proof(c) for c in record['cells']], '']
            lines += [f'def {name}Basis {params} : Fin 1 → Vec 3 :=',
                      f'  fun _ => {engine}.basis {val} {root}', '',
                      f'theorem {name}_basis_independent {params} (h : GenericChart sigma k) :',
                      f'    LinearIndependent ℝ {B} := by',
                      f'  exact {engine}.basis_independent {val} {root} h', '',
                      f'theorem {name}_basis_complete {args} (a : Vec 3) :',
                      f'    {M}.mulVec a = 0 ↔ a ∈ Submodule.span ℝ (Set.range {B}) := by',
                      f'  rw [{name}Matrix_reference {val} hr (ne_of_gt h.1),',
                      '    Matrix.smul_mulVec, smul_eq_zero, or_iff_right (by norm_num : (1/2 : ℝ) ≠ 0),',
                      f'    referenceMatrix_basis_complete rho mu sigma k a {root} hr hm h]',
                      f'  have hb : Set.range {B} = {{referenceBasis sigma k {root}}} := by',
                      f'    change Set.range (fun _ : Fin 1 => {engine}.basis {val} {root}) = _',
                      f'    rw [Set.range_const, {engine}.basis_reference rho mu sigma z k {root} h]',
                      '  rw [hb]', '',
                      f'theorem {name}_kernel {args} :',
                      f'    LinearMap.ker {M}.mulVecLin = rootMode sigma k {root} := by',
                      f'  rw [{name}Matrix_reference {val} hr (ne_of_gt h.1)]',
                      '  ext a', f'  exact coordinate_kernel_mode rho mu sigma k a {root} hr hm', '',
                      f'theorem {name}_stack_kernel {args} :',
                      f'    LinearMap.ker {S}.mulVecLin = rootMode sigma k {root} ⊓ transverseSpace k := by',
                      f'  rw [{name}Stack_reference {val} hr (ne_of_gt h.1),',
                      f'    ← {name}Matrix_reference {val} hr (ne_of_gt h.1),',
                      f'    appendConstraint_kernel, {name}_kernel {call}]', '',
                      f'theorem {name}_nullity {args} : matrixNullity {M} = 1 := by',
                      f'  rw [matrixNullity, {name}_kernel {call}]',
                      f'  exact (genericChart_mode_dimensions sigma k {root} h).1', '',
                      f'theorem {name}_transverse_nullity {args} : matrixNullity {S} = {transverse} := by',
                      f'  rw [matrixNullity, {name}_stack_kernel {call}]',
                      f'  exact (genericChart_mode_dimensions sigma k {root} h).2', '',
                      f'theorem {name}_rank {args} : {M}.rank = 2 := by',
                      f'  exact matrix_rank_of_nullity {M} ({name}_nullity {call})', '',
                      f'theorem {name}_stacked_rank {args} : {S}.rank = {3-transverse} := by',
                      f'  exact matrix_rank_of_nullity {S} ({name}_transverse_nullity {call})', '']
            for _, label, expected, semantic, number_type in fields(engine, root)[4:]:
                lines += [f'theorem {name}_{label} '+(params if label == 'basis_count' else args)+
                          f' : {semantic} = ({expected} : {number_type}) := by']
                if label == 'basis_count':
                    lines += ['  rfl', '']
                elif label == 'difference':
                    lines += [f'  rw [nullityDifference, {name}_nullity {call}, {name}_transverse_nullity {call}]',
                              '  norm_num', '']
                else:
                    lines += [f'  rw [basisCountResidual, {name}_basis_count, {name}_nullity {call}]',
                              '  norm_num', '']
            audits += [ns+'.'+name+suffix for suffix in ('Matrix_reference', 'Stack_reference',
                       '_basis_independent', '_basis_complete', '_kernel', '_stack_kernel',
                       '_nullity', '_transverse_nullity', '_rank', '_stacked_rank',
                       '_difference', '_basis_count', '_basis_residual')]
        lines += ['end '+ns, '']
    return '\n'.join(lines+['end', 'end S10Audit.CAS', '']), audits


def binding_module(builders):
    lines = ['import S10Audit.CAS.PYGenericCounts', 'import S10Audit.CAS.WLGenericCounts', '',
             'namespace S10Audit.CAS', 'open S10Pilot', 'noncomputable section', '']
    audits = []
    for engine, b in builders.items():
        lines += ['namespace '+engine+'GenericCounts', '']
        for record in b.records:
            cell, semantics = record['cells'][0], record['count_semantics']
            label = semantics['quantity']
            claim = semantics['case']+'_'+label
            lines += [f'theorem {claim} (rho mu sigma z : ℝ) (k : Vec 3)',
                      '    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_h : GenericChart sigma k) :',
                      f'    Expr.eval (values rho mu sigma z k) {cell["lean_tree"]} =',
                      f'      (({semantics["reference"]} : {semantics["number_type"]}) : ℝ) := by',
                      f'  rw [{cell["lean_claim"]} rho mu sigma z k,',
                      f'    {engine}GenericCountReference.{claim} rho mu sigma z k'+
                      ('' if label == 'basis_count' else ' _hr _hm _h')+']',
                      '  norm_num', '']
            audits.append(engine+'GenericCounts.'+claim)
        lines += ['end '+engine+'GenericCounts', '']
    return '\n'.join(lines+['end', 'end S10Audit.CAS', '']), audits


def generate(bridge, generic_builders):
    outputs, engines, builders = {}, {}, {}
    reference, audits = reference_module(generic_builders)
    outputs[bridge.GENERATED/'GenericCountReference.lean'] = reference
    for engine, path in bridge.INPUTS.items():
        b = build(bridge, engine, path)
        builders[engine] = b
        outputs[bridge.GENERATED/(engine+'GenericCounts.lean')] = b.output()
        engines[engine] = {'path': str(path.relative_to(bridge.BASE)), 'sha256': bridge.sha(path.read_bytes()),
                           'records': b.records, 'unique_trees': len(b.nodes)}
        audits += [engine+'GenericCounts.'+name for name in b.audits]
    bindings, selected = binding_module(builders)
    outputs[bridge.GENERATED/'GenericCountBindings.lean'] = bindings
    return outputs, engines, audits+selected
