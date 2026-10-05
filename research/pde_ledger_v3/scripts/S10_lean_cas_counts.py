"""Connect exceptional N2/N3/N4/N7 integer emissions to actual matrix counts.

Counts are parsed as bounded integral scalar expressions. Arithmetic equality
and units are checked first; separate semantic theorems connect each value to
matrix rank, kernel dimension, the printed basis index set, or a signed residual.
"""
import S10_lean_cas_reruns as reruns


def count_shape(bridge, value):
    bridge.require(isinstance(value, bridge.Node), 'expected scalar count')
    number = bridge.numeric(value)
    bridge.require(number is not None and number.denominator == 1 and abs(number) <= 32,
                   'count must be a bounded integral constant')
    return [value]


def fields(engine, name, n, transverse):
    M, S, B = [f'({engine}Rerun.{name}{label} rho mu sigma)' for label in ('Matrix','Stack','Basis')]
    return [('N2_RANK', 'rank', 3-n, f'{M}.rank', 'ℕ'),
            ('N2_NULLITY', 'nullity', n, f'matrixNullity {M}', 'ℕ'),
            ('N3_STACKED_RANK', 'stacked_rank', 3-transverse, f'{S}.rank', 'ℕ'),
            ('N3_TRANSVERSE_NULLITY', 'transverse_nullity', transverse, f'matrixNullity {S}', 'ℕ'),
            ('N4_NULLITY_DIFFERENCE', 'difference', n-transverse, f'nullityDifference {M} {S}', 'ℤ'),
            ('N7_BASIS_COUNT', 'basis_count', n, f'basisCount {B}', 'ℕ'),
            ('N7_BASIS_COUNT_RESIDUAL' if engine == 'PY' else 'N7_COUNT_RESIDUAL',
             'basis_residual', 0, f'basisCountResidual {B} {M}', 'ℤ')]


def build(bridge, engine, path):
    b = bridge.Builder(engine, path, namespace=engine+'Counts', support='CountReference', tactic='count_eval')
    for geometry, root, vectors, _ in reruns.CASES[engine]:
        name, n = reruns.case_name(geometry, root), len(vectors)
        prefix = reruns.prefix(engine, geometry)+f'ROOT{root+1}_'
        transverse = 0 if root == 0 else n
        for tag, label, expected, semantic, number_type in fields(engine, name, n, transverse):
            b.record(prefix+tag, lambda x: count_shape(bridge,x), [(str(expected),(0,0,0),[])])
            b.records[-1]['count_semantics'] = {
                'quantity': label, 'case': name, 'number_type': number_type, 'reference': semantic,
                'lean_semantic_claim': engine+'Counts.'+name+'_'+label,
                'coefficient_domain': ['rho != 0','mu != 0','0 < sigma','sigma != 1'],
                'subtraction': 'signed integer' if number_type == 'ℤ' else None}
    return b


def reference_module():
    lines = ['import S10Audit.CAS.CountSupport', '',
             'set_option backward.isDefEq.respectTransparency false', '',
             'namespace S10Audit.CAS', 'open S10Pilot', 'noncomputable section', '']
    audits = []
    args = '(rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1)'
    call = 'rho mu sigma hr hm hs hs1'
    for engine, cases in reruns.CASES.items():
        lines += ['namespace '+engine+'CountReference', '']
        for geometry, root, vectors, _ in cases:
            name, n = reruns.case_name(geometry, root), len(vectors)
            transverse = 0 if root == 0 else n
            p = engine+'Loci.'+geometry+'Point'
            V = f'rootMode sigma {p} {root}'
            M, S = [f'({engine}Rerun.{name}{label} rho mu sigma)' for label in ('Matrix','Stack')]
            basis_ref = engine+'Reference.'+name
            lines += [f'theorem {name}_kernel {args} :', f'    LinearMap.ker {M}.mulVecLin = {V} := by',
                      '  ext a', f'  change {M}.mulVec a = 0 ↔ a ∈ {V}',
                      f'  rw [{engine}Rerun.{name}_basis_complete rho mu sigma a hr hm hs hs1,',
                      f'    {engine}Rerun.{name}Basis_reference rho mu sigma hr (ne_of_gt hs),',
                      f'    {basis_ref}_complete sigma hs hs1]', '',
                      f'theorem {name}_stack_kernel {args} :',
                      f'    LinearMap.ker {S}.mulVecLin = {V} ⊓ transverseSpace {p} := by',
                      f'  have he : {S} = appendConstraint {M} {p} := by',
                      f'    rw [{engine}Rerun.{name}Stack_reference rho mu sigma hr (ne_of_gt hs),',
                      f'      {engine}Rerun.{name}Matrix_reference rho mu sigma hr (ne_of_gt hs)]',
                      '    rfl', f'  rw [he, appendConstraint_kernel, {name}_kernel {call}]', '',
                      f'theorem {name}_nullity {args} : matrixNullity {M} = {n} := by',
                      f'  rw [matrixNullity, {name}_kernel {call}]',
                      f'  exact ({basis_ref}_dimensions sigma hs hs1).1', '',
                      f'theorem {name}_transverse_nullity {args} : matrixNullity {S} = {transverse} := by',
                      f'  rw [matrixNullity, {name}_stack_kernel {call}]',
                      f'  exact ({basis_ref}_dimensions sigma hs hs1).2', '',
                      f'theorem {name}_rank {args} : {M}.rank = {3-n} := by',
                      f'  exact matrix_rank_of_nullity {M} ({name}_nullity {call})', '',
                      f'theorem {name}_stacked_rank {args} : {S}.rank = {3-transverse} := by',
                      f'  exact matrix_rank_of_nullity {S} ({name}_transverse_nullity {call})', '']
            for _, label, expected, semantic, number_type in fields(engine,name,n,transverse)[4:]:
                # Basis cardinality has no coefficient hypotheses. Residuals
                # additionally use the proved nullities of the actual matrices.
                lines += [f'theorem {name}_{label} '+('(rho mu sigma : ℝ)' if label == 'basis_count' else args)+
                          f' : {semantic} = ({expected} : {number_type}) := by']
                if label == 'basis_count':
                    lines += ['  rfl', '']
                elif label == 'difference':
                    lines += [f'  rw [nullityDifference, {name}_nullity {call}, {name}_transverse_nullity {call}]',
                              '  norm_num', '']
                else:
                    lines += [f'  rw [basisCountResidual, {name}_basis_count, {name}_nullity {call}]',
                              '  norm_num', '']
            # Both physical matrix blocks restore their distinct row scales.
            P = f'((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma (c • {p}) {root}) (c • {p}))'
            base = f'((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma {p} {root}) {p})'
            lines += [f'theorem {name}_physical_counts {args} (c : ℝ) (hc : c ≠ 0) :',
                      f'    matrixNullity {P} = {n} ∧ {P}.rank = {3-n} ∧',
                      f'      matrixNullity (appendConstraint {P} (c • {p})) = {transverse} ∧',
                      f'      (appendConstraint {P} (c • {p})).rank = {3-transverse} := by',
                      f'  have hM := matrix_counts_of_kernel _ _',
                      f'    (physical_modal_kernel rho mu sigma (referenceRoot rho mu sigma {p} {root}) c {p} hc)',
                      f'  have hS := matrix_counts_of_kernel _ _',
                      f'    (physical_constraint_kernel rho mu sigma (referenceRoot rho mu sigma {p} {root}) c {p} hc)',
                      f'  have hB : {base} = {M} := ({engine}Rerun.{name}Matrix_reference rho mu sigma hr (ne_of_gt hs)).symm',
                      f'  have hA : appendConstraint {base} {p} = {S} :=',
                      f'    ({engine}Rerun.{name}Stack_reference rho mu sigma hr (ne_of_gt hs)).symm',
                      '  rw [hB] at hM', '  rw [hA] at hS', '  rw [referenceRoot_scale]',
                      f'  exact ⟨hM.1.trans ({name}_nullity {call}), hM.2.trans ({name}_rank {call}),',
                      f'    hS.1.trans ({name}_transverse_nullity {call}), hS.2.trans ({name}_stacked_rank {call})⟩', '']
            audits += [engine+'CountReference.'+name+'_'+label for label in
                       ('kernel','stack_kernel','nullity','transverse_nullity','rank','stacked_rank',
                        'difference','basis_count','basis_residual','physical_counts')]
        lines += ['end '+engine+'CountReference', '']
    lines += ['end','end S10Audit.CAS','']
    return '\n'.join(lines), audits


def binding_module(builders):
    lines = ['import S10Audit.CAS.PYCounts', 'import S10Audit.CAS.WLCounts', '',
             'namespace S10Audit.CAS', 'open S10Pilot', 'noncomputable section', '']
    audits = []
    for engine, b in builders.items():
        lines += ['namespace '+engine+'Counts','']
        for record in b.records:
            cell, semantics = record['cells'][0], record['count_semantics']
            name, label = semantics['case'], semantics['quantity']
            claim = name+'_'+label
            lines += [f'theorem {claim} (rho mu sigma z : ℝ) (k : Vec 3)',
                      '    (_hr : rho ≠ 0) (_hm : mu ≠ 0) (_hs : 0 < sigma) (_hs1 : sigma ≠ 1) :',
                      f'    Expr.eval (values rho mu sigma z k) {cell["lean_tree"]} =',
                      f'      (({semantics["reference"]} : {semantics["number_type"]}) : ℝ) := by',
                      f'  rw [{cell["lean_claim"]} rho mu sigma z k,',
                      f'    {engine}CountReference.{claim} rho mu sigma'+
                      ('' if label == 'basis_count' else ' _hr _hm _hs _hs1')+']',
                      '  norm_num','']
            audits.append(engine+'Counts.'+claim)
        lines += ['end '+engine+'Counts','']
    lines += ['end','end S10Audit.CAS','']
    return '\n'.join(lines), audits


def generate(bridge):
    outputs, engines, builders, audits = {}, {}, {}, []
    ref, selected = reference_module()
    outputs[bridge.GENERATED/'CountReference.lean'] = ref
    audits += selected
    for engine, path in bridge.INPUTS.items():
        b = build(bridge,engine,path)
        builders[engine] = b
        outputs[bridge.GENERATED/(engine+'Counts.lean')] = b.output()
        engines[engine] = {'path': str(path.relative_to(bridge.BASE)), 'sha256': bridge.sha(path.read_bytes()),
                           'records': b.records, 'unique_trees': len(b.nodes)}
        audits += [engine+'Counts.'+name for name in b.audits]
    bindings, selected = binding_module(builders)
    outputs[bridge.GENERATED/'CountBindings.lean'] = bindings
    audits += selected
    audits += ['appendConstraint_kernel','matrix_rank_nullity','matrix_rank_of_nullity',
               'matrix_counts_of_kernel','physical_modal_kernel','physical_constraint_kernel']
    return outputs, engines, audits
