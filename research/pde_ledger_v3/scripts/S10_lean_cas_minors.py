"""D3 minor selections; no CAS simplification is used to establish coverage.

Maps record the printed order (SymPy deduplicates, Wolfram does not). Every
selection-map entry is independently checked against its determinant in Lean.
"""
from itertools import combinations

FAMILIES = [
    ('staticRank', 0, False, 2, 6, [0, 1, 2, 1, 3, 4, 2, 4, 5]),
    ('staticTransverse', 0, True, 3, 4, [0, 1, 2, 3]),
    ('ordinaryRank', 1, False, 2, 4, [0, 1, 2, 1, 3, 2, 2, 2, 2]),
    ('ordinaryTransverse', 1, True, 2, 6,
     [0, 1, 2, 1, 3, 2, 4, 5, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]),
    ('extraRank', 2, False, 2, 6, [0, 1, 2, 1, 3, 4, 2, 4, 5]),
    ('extraTransverse', 2, True, 3, 4, [0, 1, 2, 3]),
]


def family_data(family):
    label, root, stacked, order, count, mapping = family
    rows = list(combinations(range(4 if stacked else 3), order))
    cols = list(combinations(range(3), order))
    assert len(mapping) == len(rows)*len(cols) and set(mapping) == set(range(count))
    suffix = f'ROOT{root+1}_Q8_'+('TRANSVERSE_' if stacked else '')+'RANK_DROP_MINORS'
    reference = label+'Minors mu '+('sigma ' if root else '')+'k'
    matrix = f'{"rootStack" if stacked else "rootMatrix"} rho mu sigma k {root}'
    row_table = f'choose{order}{4 if stacked else 3}'
    col_table = f'choose{order}3'
    return label, root, order, count, mapping, rows, cols, suffix, reference, matrix, row_table, col_table


def build(bridge, engine, path):
    b = bridge.Builder(engine, path, namespace=engine+'Minors', support='MinorSupport', tactic='minor_eval')
    for family in FAMILIES:
        label, root, order, count, mapping, rows, cols, suffix, ref, *_ = family_data(family)
        selected = list(range(len(mapping))) if engine == 'WL' else [mapping.index(i) for i in range(count)]
        cells = []
        provenance = []
        for position in selected:
            rr, cc = rows[position//len(cols)], cols[position % len(cols)]
            dim = tuple(sum(((-1, 0, 0) if r == 3 else (-3, -2, 1))[i] for r in rr)
                        for i in range(3))
            cells.append((f'{ref} {mapping[position]}', dim, []))
            provenance.append({'rows': rr, 'columns': cc, 'order': order,
                               'canonical_index': mapping[position]})
        b.record(suffix, lambda x: x, cells)
        b.records[-1]['minor_selections'] = provenance
        b.records[-1]['full_selection_map'] = mapping
        b.records[-1]['canonical_count'] = count
    return b


def reference_module():
    lines = ['import S10Audit.CAS.MinorSupport', '',
             'set_option backward.isDefEq.respectTransparency false', '',
             'namespace S10Audit.CAS', 'open S10Pilot', 'noncomputable section', '']
    audits = []
    for family in FAMILIES:
        label, root, order, count, mapping, rows, cols, suffix, ref, matrix, rt, ct = family_data(family)
        nr, nc = len(rows), len(cols)
        table = '!['+', '.join('!['+', '.join(str(n) for n in mapping[i*nc:(i+1)*nc])+']'
                              for i in range(nr))+']'
        lines += [f'def {label}Map : Fin {nr} → Fin {nc} → Fin {count} := {table}',
                  f'theorem {label}Map_onto : ∀ c, ∃ i j, {label}Map i j = c := by decide', '',
                  f'theorem {label}_det (rho mu sigma : ℝ) (k : Vec 3)',
                  '    (_hr : rho ≠ 0) (_hs : sigma ≠ 0) :',
                  f'    ∀ i j, {ref} ({label}Map i j) =',
                  f'      Matrix.det (({matrix}).submatrix ({rt} i) ({ct} j)) := by',
                  '  intro i j', '  fin_cases i']
        for i in range(nr):
            lines.append('  · fin_cases j')
            for j in range(nc):
                def entry(a, b):
                    row, col = rows[i][a], cols[j][b]
                    return f'(k {col})' if row == 3 else f'(rootMatrix rho mu sigma k {root} {row} {col})'
                if order == 2:
                    determinant = f'{entry(0,0)} * {entry(1,1)} - {entry(0,1)} * {entry(1,0)}'
                else:
                    terms = [[(0,0),(1,1),(2,2)],[(0,0),(1,2),(2,1)],
                             [(0,1),(1,0),(2,2)],[(0,1),(1,2),(2,0)],
                             [(0,2),(1,0),(2,1)],[(0,2),(1,1),(2,0)]]
                    signs = ['', ' - ', ' - ', ' + ', ' + ', ' - ']
                    determinant = ''.join(sign+' * '.join(entry(a,b) for a,b in term)
                                          for sign,term in zip(signs,terms))
                lines.append(f'    · rw [Matrix.det_fin_{"two" if order == 2 else "three"}]')
                lines.append(f'      change {ref} {mapping[i*nc+j]} = {determinant}')
                lines.append('      minor_equal')
        lines += ['', f'theorem {label}_complete (rho mu sigma : ℝ) (k : Vec 3)',
                  '    (hr : rho ≠ 0) (hs : sigma ≠ 0) :',
                  f'    (∀ c, {ref} c = 0) ↔ AllOrderedMinorsZero ({matrix}) {order} := by',
                  f'  rw [← selection_table_complete ({matrix}) {rt} {ct}',
                  f'    {rt}_complete {ct}_complete {rt}_strict {ct}_strict]',
                  '  constructor', '  · intro h i j',
                  f'    rw [← {label}_det rho mu sigma k hr hs i j]',
                  '    exact h _', '  · intro h c',
                  f'    obtain ⟨i, j, rfl⟩ := {label}Map_onto c',
                  f'    rw [{label}_det rho mu sigma k hr hs i j]',
                  '    exact h i j', '']
        audits += [label+'_det', label+'_complete']
    lines += ['end', 'end S10Audit.CAS', '']
    return '\n'.join(lines), audits


def binding_module(builders):
    lines = ['import S10Audit.CAS.PYMinors', 'import S10Audit.CAS.WLMinors',
             'import S10Audit.CAS.MinorReference', '',
             'set_option backward.isDefEq.respectTransparency false', '',
             'namespace S10Audit.CAS', 'open S10Pilot', 'noncomputable section', '']
    audits = []
    for engine, builder in builders.items():
        lines.append(f'namespace {engine}Minors\n')
        for family, record in zip(FAMILIES, builder.records):
            label, root, order, count, mapping, rows, cols, suffix, ref, matrix, *_ = family_data(family)
            cells = record['cells']
            indices = [s['canonical_index'] for s in record['minor_selections']]
            n = len(cells)
            trees = '!['+', '.join(c['lean_tree'] for c in cells)+']'
            lines += [f'def {label} (rho mu sigma z : ℝ) (k : Vec 3) : Fin {n} → ℝ :=',
                      f'  fun i => Expr.eval (values rho mu sigma z k) (({trees} : Fin {n} → Expr Symbol) i)',
                      f'def {label}Index : Fin {n} → Fin {count} := !['+', '.join(map(str, indices))+']',
                      f'theorem {label}Index_onto : Function.Surjective {label}Index := by decide', '',
                      f'theorem {label}_reference (rho mu sigma z : ℝ) (k : Vec 3) (_hs : sigma ≠ 0) :',
                      f'    ∀ i, {label} rho mu sigma z k i = {ref} ({label}Index i) := by',
                      '  intro i', '  fin_cases i']
            for cell in cells:
                args = []
                for d in cell['raw_denominator_factors']:
                    assert d == ['atom', 'sigma']
                    args.append('_hs')
                lines.append('  · exact '+cell['lean_claim']+' rho mu sigma z k'+''.join(' '+a for a in args))
            lines += ['', f'theorem {label}_zero_iff (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :',
                      f'    (∀ i, {label} rho mu sigma z k i = 0) ↔ (∀ c, {ref} c = 0) := by',
                      '  constructor', '  · intro h c',
                      f'    obtain ⟨i, rfl⟩ := {label}Index_onto c',
                      f'    rw [← {label}_reference rho mu sigma z k hs i]', '    exact h i',
                      f'  · intro h i; rw [{label}_reference rho mu sigma z k hs i]; exact h _', '',
                      f'theorem {label}_complete (rho mu sigma z : ℝ) (k : Vec 3)',
                      '    (hr : rho ≠ 0) (hs : sigma ≠ 0) :',
                      f'    (∀ i, {label} rho mu sigma z k i = 0) ↔ AllOrderedMinorsZero ({matrix}) {order} := by',
                      f'  rw [{label}_zero_iff rho mu sigma z k hs, CAS.{label}_complete rho mu sigma k hr hs]', '']
            audits += [f'{engine}Minors.{label}_{s}' for s in ['reference', 'zero_iff', 'complete']]
        lines.append(f'end {engine}Minors\n')
    for family in FAMILIES:
        label, *_ = family
        lines += [f'theorem {label}_zero_cross_engine (rho mu sigma z : ℝ) (k : Vec 3) (hs : sigma ≠ 0) :',
                  f'    (∀ i, PYMinors.{label} rho mu sigma z k i = 0) ↔',
                  f'      (∀ i, WLMinors.{label} rho mu sigma z k i = 0) := by',
                  f'  rw [PYMinors.{label}_zero_iff rho mu sigma z k hs, WLMinors.{label}_zero_iff rho mu sigma z k hs]', '']
        audits.append(label+'_zero_cross_engine')
    lines += ['end', 'end S10Audit.CAS', '']
    return '\n'.join(lines), audits


def generate(bridge):
    outputs, engines, builders, audits = {}, {}, {}, []
    for engine, path in bridge.INPUTS.items():
        b = build(bridge, engine, path)
        builders[engine] = b
        outputs[bridge.GENERATED/(engine+'Minors.lean')] = b.output()
        engines[engine] = {'path': str(path.relative_to(bridge.BASE)),
                           'sha256': bridge.sha(path.read_bytes()), 'records': b.records,
                           'unique_trees': len(b.nodes)}
        audits += [engine+'Minors.'+a for a in b.audits]
    reference, ra = reference_module()
    bindings, ba = binding_module(builders)
    outputs[bridge.GENERATED/'MinorReference.lean'] = reference
    outputs[bridge.GENERATED/'MinorBindings.lean'] = bindings
    audits += ra+ba+[f'choose{q}{n}_{s}' for q, n in [(2,3),(2,4),(3,3),(3,4)]
                     for s in ['complete','strict']]+['selection_table_complete']
    return outputs, engines, audits
