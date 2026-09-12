"""Exceptional rerun arithmetic, complete finite bases, and physical rescaling.

Reference bases are specified independently below. The generated proofs check
membership, independence, and equality with the full theoretical kernel.
The transcript parser preserves the printed basis ordering, including WL's
reversed pair at the parallel double root.
"""

# geometry, root index, ordered reference vectors, pivot coordinates
CASES = {
    'PY': [('parallel', 0, ['![1, 0, 0]'], [0]),
           ('parallel', 1, ['![0, 1, 0]', '![0, 0, 1]'], [1, 2]),
           ('perpendicular', 0, ['![0, 1, 0]'], [1]),
           ('perpendicular', 1, ['![0, 0, 1]'], [2]),
           ('perpendicular', 2, ['![1, 0, 0]'], [0])],
    'WL': [('parallel', 0, ['![1, 0, 0]'], [0]),
           ('parallel', 1, ['![0, 0, 1]', '![0, 1, 0]'], [2, 1]),
           ('perpendicular', 0, ['![0, -54, 1]'], [2]),
           ('perpendicular', 1, ['![0, 1/54, 1]'], [2]),
           ('perpendicular', 2, ['![1, 0, 0]'], [0])],
}


def case_name(geometry, root):
    return geometry+str(root)


def prefix(engine, geometry):
    return ('Q8_' if engine == 'PY' else '')+f'STRATUM{1 if geometry == "parallel" else 2}_'


def build(bridge, engine, path):
    dims = {**bridge.DIMS, 'z': (2, -2, 0), **{f'k{i}': (0, 0, 0) for i in range(3)}}
    def zero(dim):
        monomial = {(-1, -2, 1): 'mu_R', (-3, -6, 3): 'mu_R**3',
                    (2, -2, 0): 'mu_R/rho_br'}[dim]
        return bridge.Node('mul', (bridge.Node('num', (0,)), bridge.Parser(monomial, 'PY').parse()))
    b = bridge.Builder(engine, path, namespace=engine+'Rerun', support='RerunReference',
                       tactic='rerun_eval', units='coordinateUnits', dims=dims, zero_adapter=zero)
    md, zd, bd = (-1, -2, 1), (2, -2, 0), (0, 0, 0)
    for geometry in ('parallel', 'perpendicular'):
        p = engine+'Loci.'+geometry+'Point'
        pre = prefix(engine, geometry)
        nroots = 2 if geometry == 'parallel' else 3
        b.record(pre+'Q3_DETERMINANT', bridge.scalar_shape,
                 [(f'Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma z {p})', (-3, -6, 3), [])])
        b.record(pre+'ROOT_ORDERING', lambda x: x,
                 [(f'referenceRoot rho mu sigma {p} {r}', zd,
                   [] if r == 0 else ['rho']+(['sigma'] if r == 2 else [])) for r in range(nroots)])
    for geometry, r, vectors, _ in CASES[engine]:
        p = engine+'Loci.'+geometry+'Point'
        name = engine+'Reference.'+case_name(geometry, r)
        pre = prefix(engine, geometry)+f'ROOT{r+1}_'
        matrix = f'((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma {p} {r}) {p})'
        domain = [] if r == 0 else ['rho']+(['sigma'] if r == 2 else [])
        def basis_shape(x):
            bridge.require(isinstance(x, list) and len(x) == len(vectors), 'wrong exceptional basis cardinality')
            return [cell for vector in x for cell in bridge.vector_shape(vector, engine)]
        b.record(pre+'N1_MATRIX', lambda x: bridge.matrix_shape(x, 3, 3),
                 [(f'{matrix} {i} {j}', md, domain) for i in range(3) for j in range(3)])
        b.record(pre+'N3_STACKED_MATRIX', lambda x: bridge.matrix_shape(x, 4, 3),
                 [(f'{matrix} {i} {j}' if i < 3 else f'{p} {j}', md if i < 3 else bd,
                   domain if i < 3 else []) for i in range(4) for j in range(3)])
        b.record(pre+('N5_MATRIX_TIMES_K' if engine == 'PY' else 'N5_WAVEVECTOR_PRODUCT'),
                 lambda x: bridge.vector_shape(x, engine),
                 [(f'{matrix}.mulVec {p} {i}', md, domain) for i in range(3)])
        b.record(pre+'N6_NULLSPACE_BASIS', basis_shape,
                 [(f'{name} {j} {i}', bd, []) for j in range(len(vectors)) for i in range(3)])
        b.record(pre+('N6_BASIS_DOT_K' if engine == 'PY' else 'N6_BASIS_DOTS'), lambda x: x,
                 [(f'dot {p} ({name} {j})', bd, []) for j in range(len(vectors))])
        b.record(pre+('N6_BASIS_VECTOR_RESIDUALS' if engine == 'PY' else 'N6_BASIS_RESIDUALS'), basis_shape,
                 [(f'normSq {p} * {name} {j} {i} - dot {p} ({name} {j}) * {p} {i}', bd, [])
                  for j in range(len(vectors)) for i in range(3)])
    return b


def reference_module():
    lines = ['import S10Audit.CAS.RerunSupport', '',
             'set_option backward.isDefEq.respectTransparency false', '',
             'namespace S10Audit.CAS', 'open S10Pilot S10Anisotropic', 'noncomputable section', '']
    audits, definitions = [], []
    for engine, cases in CASES.items():
        lines += ['namespace '+engine+'Reference', '']
        for geometry, r, vectors, pivots in cases:
            name, n = case_name(geometry, r), len(vectors)
            p = engine+'Loci.'+geometry+'Point'
            V = f'rootMode sigma {p} {r}'
            lines += [f'def {name} : Fin {n} → Vec 3 := !['+', '.join(vectors)+']', '',
                      f'theorem {name}_independent : LinearIndependent ℝ {name} := by',
                      '  rw [Fintype.linearIndependent_iff]', '  intro c hc j', '  fin_cases j']
            for pivot in pivots:
                lines += [f'  · have h := congrFun hc {pivot}',
                          f'    simpa [{name}, Fin.sum_univ_succ, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] using h']
            lines += ['', f'theorem {name}_members (sigma : ℝ) (_hs : sigma ≠ 0) (j : Fin {n}) :',
                      f'    {name} j ∈ {V} := by',
                      f'  change normalizedOperator 0 sigma ({"0" if r == 0 else "normSq "+p if r == 1 else "extraValue 0 sigma "+p}) {p} ({name} j) = 0',
                      '  ext i', '  fin_cases j', '  all_goals fin_cases i', '  all_goals',
                      f'    norm_num [{name}, {p}, normalizedOperator, extraValue, extraNumerator, perpSq,',
                      '      normSq, dot, unit, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, _hs]',
                      *(['  all_goals cas_equal'] if r == 2 else []), '',
                      f'theorem {name}_dimensions (sigma : ℝ) ({"_hs" if r == 0 else "hs"} : 0 < sigma) ({"_hs1" if r == 0 or geometry == "parallel" else "hs1"} : sigma ≠ 1) :',
                      f'    Module.finrank ℝ ({V}) = {n} ∧',
                      f'      Module.finrank ℝ ({V} ⊓ transverseSpace {p} : Submodule ℝ (Vec 3)) = {0 if r == 0 else n} := by']
            if r == 0:
                lines += [f'  exact zero_counts (sigma := sigma) {p}_nonzero']
            elif geometry == 'parallel':
                lines += [f'  apply parallel_counts (ne_of_gt hs) {p}_nonzero',
                          f'  norm_num [perpSq, normSq, dot, Fin.sum_univ_three, {p}, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]']
            else:
                lines += [f'  have h := perpendicular_counts hs hs1 {p}_nonzero (by norm_num [{p}] : {p} 0 = 0)',
                          '  exact '+('⟨h.1, h.2.1⟩' if r == 1 else '⟨h.2.2.1, h.2.2.2⟩')]
            lines += ['', f'theorem {name}_complete (sigma : ℝ) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :',
                      f'    Submodule.span ℝ (Set.range {name}) = {V} := by',
                      f'  exact complete_of_dimension {name} _ {name}_independent',
                      f'    ({name}_members sigma (ne_of_gt hs)) ({name}_dimensions sigma hs hs1).1', '']
            audits += [engine+'Reference.'+name+'_'+suffix for suffix in ('independent','members','dimensions','complete')]
            definitions.append(engine+'Reference.'+name)
        lines += ['end '+engine+'Reference', '']
    definitions += [e+'Loci.'+g+'Point' for e in CASES for g in ('parallel','perpendicular')]
    lines += ['open Lean Elab Tactic in', 'elab "rerun_equal" : tactic => do',
              '  evalTactic (← `(tactic| norm_num ['+', '.join(definitions)+',',
              '    referenceMatrix, referenceRoot, Matrix.det_fin_three, Matrix.mulVec, dotProduct,',
              '    extraConeValue, extraValue, extraNumerator, perpSq, coneValue, normSq, dot, Fin.sum_univ_three,',
              '    Fin.ext_iff, Fin.coe_ofNat_eq_mod, -Fin.val_eq_zero_iff,',
              '    Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]))',
              '  if !(← getGoals).isEmpty then', '    evalTactic (← `(tactic| field_simp))',
              '  if !(← getGoals).isEmpty then', '    evalTactic (← `(tactic| ring_nf))',
              '  if !(← getGoals).isEmpty then', '    evalTactic (← `(tactic| norm_num))', '',
              'open Lean Elab Tactic in', 'elab "rerun_eval" h:Lean.Parser.Tactic.rwRule : tactic => do',
              '  evalTactic (← `(tactic| rw [$h]))', '  if !(← getGoals).isEmpty then',
              '    evalTactic (← `(tactic| rerun_equal))', '', 'end', 'end S10Audit.CAS', '']
    return '\n'.join(lines), audits


def binding_module(builders):
    lines = ['import S10Audit.CAS.PYRerun', 'import S10Audit.CAS.WLRerun', '',
             'set_option backward.isDefEq.respectTransparency false', '',
             'namespace S10Audit.CAS', 'open S10Pilot', 'noncomputable section', '']
    audits = []
    for engine, b in builders.items():
        lines += ['namespace '+engine+'Rerun', '']
        for geometry, r, vectors, _ in CASES[engine]:
            name, n = case_name(geometry, r), len(vectors)
            ref = engine+'Reference.'+name
            p = engine+'Loci.'+geometry+'Point'
            pre = prefix(engine, geometry)+f'ROOT{r+1}_'
            def get(suffix):
                return next(rec for rec in b.records if rec['tag'].endswith('_'+pre+suffix))
            for label, suffix, rows in [('Matrix','N1_MATRIX',3),('Stack','N3_STACKED_MATRIX',4),
                                        ('Basis','N6_NULLSPACE_BASIS',n)]:
                record = get(suffix)
                names = [c['lean_tree'] for c in record['cells']]
                tree = '!['+', '.join('!['+', '.join(names[3*i:3*i+3])+']' for i in range(rows))+']'
                target = f'((1/2 : ℝ) • referenceMatrix rho mu sigma (referenceRoot rho mu sigma {p} {r}) {p})'
                if label == 'Stack':
                    target = '!['+', '.join(f'{target} {i}' for i in range(3))+f', {p}]'
                if label == 'Basis':
                    target = ref
                def cell_proof(cell):
                    factors = cell['raw_denominator_factors'] + [['atom', x] for x in cell['reference_domain_symbols']]
                    assert all(d[0] == 'atom' and d[1] in ('rho', 'sigma') for d in factors)
                    symbols = dict.fromkeys(d[1] for d in factors)
                    arguments = ''.join(' _hr' if x == 'rho' else ' _hs' for x in symbols)
                    return 'exact '+cell['lean_claim']+f' rho mu sigma 0 {p}'+arguments
                lines += [f'def {name}{label} (rho mu sigma : ℝ) : '+(f'Fin {rows} → Vec 3' if label == 'Basis' else f'Matrix (Fin {rows}) (Fin 3) ℝ')+' :=',
                          f'  fun i j => Expr.eval (values rho mu sigma 0 {p})',
                          f'    (({tree} : Fin {rows} → Fin 3 → Expr Symbol) i j)', '',
                          f'theorem {name}{label}_reference (rho mu sigma : ℝ) (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :',
                          f'    {name}{label} rho mu sigma = {target} := by', '  ext i j',
                          '  fin_cases i', '  all_goals fin_cases j',
                          *['  · '+cell_proof(c) for c in record['cells']], '']
                audits.append(engine+'Rerun.'+name+label+'_reference')
            lines += [f'theorem {name}_basis_independent (rho mu sigma : ℝ) (hr : rho ≠ 0) (hs : sigma ≠ 0) :',
                      f'    LinearIndependent ℝ ({name}Basis rho mu sigma) := by',
                      f'  rw [{name}Basis_reference rho mu sigma hr hs]', f'  exact {ref}_independent', '',
                      f'theorem {name}_basis_complete (rho mu sigma : ℝ) (a : Vec 3)',
                      '    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :',
                      f'    ({name}Matrix rho mu sigma).mulVec a = 0 ↔',
                      f'      a ∈ Submodule.span ℝ (Set.range ({name}Basis rho mu sigma)) := by',
                      f'  rw [{name}Matrix_reference rho mu sigma hr (ne_of_gt hs),',
                      f'    {name}Basis_reference rho mu sigma hr (ne_of_gt hs),',
                      f'    {ref}_complete sigma hs hs1, coordinate_kernel_mode rho mu sigma {p} a {r} hr hm]', '',
                      f'theorem {name}_physical_complete (rho mu sigma c : ℝ) (a : Vec 3)',
                      '    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1) (hc : c ≠ 0) :',
                      f'    ((1/2 : ℝ) • referenceMatrix rho mu sigma',
                      f'      (referenceRoot rho mu sigma (c • {p}) {r}) (c • {p})).mulVec a = 0 ↔',
                      f'      a ∈ Submodule.span ℝ (Set.range ({name}Basis rho mu sigma)) := by',
                      f'  rw [referenceRoot_scale, coordinate_kernel_scale rho mu sigma _ c {p} a hc,',
                      f'    ← {name}Matrix_reference rho mu sigma hr (ne_of_gt hs)]',
                      f'  exact {name}_basis_complete rho mu sigma a hr hm hs hs1', '']
            audits += [engine+'Rerun.'+name+'_'+s for s in ('basis_independent','basis_complete','physical_complete')]
        lines += ['end '+engine+'Rerun', '']
    lines += ['end', 'end S10Audit.CAS', '']
    return '\n'.join(lines), audits


def generate(bridge):
    outputs, engines, builders, audits = {}, {}, {}, []
    ref, selected = reference_module()
    outputs[bridge.GENERATED/'RerunReference.lean'] = ref
    audits += selected
    for engine, path in bridge.INPUTS.items():
        b = build(bridge, engine, path)
        builders[engine] = b
        outputs[bridge.GENERATED/(engine+'Rerun.lean')] = b.output()
        engines[engine] = {'path': str(path.relative_to(bridge.BASE)), 'sha256': bridge.sha(path.read_bytes()),
                           'unit_convention': 'k = κ p; omegaSquared = κ² z; [κ] = L^-1; p dimensionless; [z] = L² T^-2',
                           'records': b.records, 'unique_trees': len(b.nodes)}
        audits += [engine+'Rerun.'+name for name in b.audits]
    bindings, selected = binding_module(builders)
    outputs[bridge.GENERATED/'RerunBindings.lean'] = bindings
    audits += selected
    audits += ['coordinate_'+s+'_units' for s in ('frequency','matrix','determinant','product','wavevector','residual')]
    audits += ['normSq_scale','referenceMatrix_scale','referenceRoot_scale']
    audits += ['coordinate_'+s+'_scale' for s in ('determinant','kernel','product','dot','residual')]
    audits += ['coordinate_kernel_mode','complete_of_dimension']
    return outputs, engines, audits
