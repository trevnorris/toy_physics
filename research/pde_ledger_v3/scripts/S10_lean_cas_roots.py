"""D3 solution lists, distinct roots, multiplicities and syntactic filter counts.

Only single-variable arithmetic solution rules and the two named count fields
are accepted. Raw container text is retained alongside the extracted trees.
No solver implementation or conditional-expression grammar is trusted here.
"""
import re
import S10_lean_cas_counts as counts
import S10_lean_cas_reruns as reruns

CASES = ('generic', 'parallel', 'perpendicular')


def container_parser(bridge, raw, engine):
    class Parser(bridge.Parser):
        token = re.compile(r'->|:|<\||\|>|"FilteredCandidateCount"|"DistinctRootCount"|'+bridge.Parser.token.pattern)
    return Parser(raw, engine)


def solutions(bridge, raw, engine):
    p = container_parser(bridge, raw, engine)
    opening, closing = ('[', ']') if engine == 'PY' else ('{', '}')
    p.take(opening)
    roots = []
    while p.peek() != closing:
        if roots:
            p.take(',')
        p.take('{')
        p.take('omegaSquared')
        p.take(':' if engine == 'PY' else '->')
        root = p.expression()
        bridge.require(isinstance(root, bridge.Node), 'expected one scalar root per solution rule')
        roots.append(root)
        bridge.require(len(roots) <= 3, 'too many D3 root candidates')
        p.take('}')
    p.take(closing)
    bridge.require(p.peek() is None and roots, 'empty or trailing solution payload')
    return roots


def list_counts(bridge, raw, engine):
    bridge.require(engine == 'WL', 'root-list association is Wolfram-only')
    p = container_parser(bridge, raw, engine)
    p.take('<|')
    values = []
    for index, key in enumerate(('FilteredCandidateCount', 'DistinctRootCount')):
        if index:
            p.take(',')
        p.take('"'+key+'"')
        p.take('->')
        values += counts.count_shape(bridge, p.expression())
    p.take('|>')
    bridge.require(p.peek() is None, 'trailing root-list association')
    return values


def build(bridge, engine, path, case):
    generic = case == 'generic'
    ns = engine+'Roots'+case.title()
    dims = bridge.DIMS if generic else {**bridge.DIMS, 'z': (2,-2,0), **{f'k{i}':(0,0,0) for i in range(3)}}
    dim = (0,-2,0) if generic else (2,-2,0)
    def zero(d):
        return bridge.zero_lift(d) if generic else bridge.Node('mul',
            (bridge.Node('num',(0,)),bridge.Parser('mu_R/rho_br','PY').parse()))
    b = bridge.Builder(engine, path, namespace=ns, support='RootSupport' if generic else 'RerunReference',
                       tactic='cas_eval' if generic else 'rerun_eval',
                       units='units' if generic else 'coordinateUnits', dims=dims, zero_adapter=zero)
    pre = '' if generic else reruns.prefix(engine, case)
    p = 'k' if generic else engine+'Loci.'+case+'Point'
    distinct = [0,1] if case == 'parallel' else [0,1,2]
    candidates = [0,1,1] if case == 'parallel' and engine == 'WL' else distinct
    for label, suffix, indices, parser in [
        ('candidates', 'Q3_ROOT_SOLUTIONS_RAW' if engine == 'PY' else 'Q3_SOLUTIONS', candidates,
         lambda raw,e: solutions(bridge,raw,e)),
        ('distinct', 'Q3_ROOTS_DISTINCT' if engine == 'PY' else 'Q3_DISTINCT_ROOTS', distinct, None)]:
        b.record(pre+suffix, lambda x:x,
                 [(f'referenceRoot rho mu sigma {p} {r}', dim,
                   [] if r == 0 else ['rho']+(['sigma'] if r == 2 else [])) for r in indices],
                 parse_payload=parser)
        b.records[-1]['root_semantics'] = {'case':case,'kind':label,'root_indices':indices,
            'container':'single-variable solution rules' if label == 'candidates' else 'arithmetic list',
            'lean_semantic_claims':[ns+'.'+label+'_reference', ns+'.'+label+'_complete']}
        if label == 'candidates':
            raw = b.records[-1]['payload']
            for cell,node in zip(b.records[-1]['cells'],solutions(bridge,raw,engine)):
                cell['lean_raw_tree'] = b.intern(node)
        if label == 'distinct':
            b.records[-1]['root_semantics']['lean_semantic_claims'] += [ns+'.distinct_nodup',ns+'.multiplicities']
        elif not (case == 'parallel' and engine == 'PY'):
            b.records[-1]['root_semantics']['lean_semantic_claims'] += [ns+'.candidate_multiset']
    numeric_fields = [('Q3_ROOT_COUNT','distinct_count',len(distinct),'distinct')]
    if engine == 'WL':
        numeric_fields += [('Q3_ROOT_CANDIDATE_COUNT_BEFORE_FILTER','candidate_count',len(candidates),'candidates'),
                           ('Q3_ROOT_CANDIDATE_COUNT_AFTER_FILTER','filtered_count',len(candidates),'filtered')]
    for suffix,label,expected,operand in numeric_fields:
        b.record(pre+suffix,lambda x:counts.count_shape(bridge,x),[(str(expected),(0,0,0),[])])
        b.records[-1]['root_semantics'] = {'case':case,'kind':label,'operand':operand,
            'lean_semantic_claims':[ns+'.'+label]}
    if engine == 'WL':
        b.record(pre+'Q3_ROOT_LIST_COUNTS',lambda x:x,
                 [(str(len(candidates)),(0,0,0),[]),(str(len(distinct)),(0,0,0),[])],
                 parse_payload=lambda raw,e:list_counts(bridge,raw,e))
        b.records[-1]['root_semantics'] = {'case':case,'kind':'list_counts',
            'fields':['FilteredCandidateCount','DistinctRootCount'],
            'lean_semantic_claims':[ns+'.list_counts']}
    for r in b.records:
        r['root_semantics']['coefficient_domain'] = ['rho != 0','mu != 0','0 < sigma','sigma != 1']
        r['root_semantics']['coordinate_domain'] = 'GenericChart sigma k' if generic else p
    return b


def binding_module(bridge, builders):
    lines = ['import S10Audit.CAS.RootSupport',
             *['import S10Audit.CAS.'+b.namespace for b in builders.values()], '',
             'set_option backward.isDefEq.respectTransparency false', '',
             'namespace S10Audit.CAS', 'open S10Pilot S10Anisotropic Polynomial', 'noncomputable section', '']
    audits, controls = [], []
    for (engine, case), b in builders.items():
        ns = b.namespace
        generic = case == 'generic'
        p = 'k' if generic else engine+'Loci.'+case+'Point'
        pre = '' if generic else reruns.prefix(engine,case)
        params, val = '(rho mu sigma z : ℝ) (k : Vec 3)', 'rho mu sigma z k'
        domain = '(hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1)'
        if generic:
            domain += ' (hg : GenericChart sigma k)'
        args, call = params+' '+domain, val+' hr hm hs hs1'+(' hg' if generic else '')
        refargs, refcall = params+' (_hr : rho ≠ 0) (_hs : sigma ≠ 0)', val+' hr (ne_of_gt hs)'
        records = {r['root_semantics']['kind']:r for r in b.records}
        lines += ['namespace '+ns, '']
        for label in ('candidates','distinct'):
            r = records[label]
            cells = r['cells']
            trees = '['+', '.join(c.get('lean_raw_tree',c['lean_tree']) for c in cells)+']'
            refs = '['+', '.join(f'referenceRoot rho mu sigma {p} {i}' for i in r['root_semantics']['root_indices'])+']'
            lines += [f'def {label}Trees : List (Expr Symbol) := {trees}',
                      f'def {label} {params} : List ℝ := {label}Trees.map (Expr.eval (values {val}))', '',
                      f'theorem {label}_reference {refargs} : {label} {val} = {refs} := by',
                      f'  unfold {label} {label}Trees', '  simp only [List.map_cons, List.map_nil]']
            rules = {}
            for c in cells:
                factors = c['raw_denominator_factors']+[['atom',x] for x in c['reference_domain_symbols']]
                assert all(x[0]=='atom' and x[1] in ('rho','sigma') for x in factors)
                extra = ''.join(' _hr' if x=='rho' else ' _hs' for x in dict.fromkeys(x[1] for x in factors))
                tree = c.get('lean_raw_tree',c['lean_tree'])
                rules[tree] = (tree+'_eval '+val if label=='candidates' and c['literal_zero_lift'] else
                               c['lean_claim']+' '+val+extra)
            lines += ['  rw ['+', '.join(rules.values())+']']
            if label == 'candidates':
                lines += ['  rfl']
            lines += ['']
            audits += [ns+'.'+label+'_reference']
        # Full determinant root multiset. The repeated nonzero factor survives.
        setup = [f'  have hk : {p} ≠ 0 := '+('(genericChart_split sigma k hg).1' if generic else engine+'Loci.'+case+'Point_nonzero'),
                 f'  have hq : perpSq 0 {p} '+('= 0' if case=='parallel' else '≠ 0')+' := '+
                 ('(genericChart_split sigma k hg).2' if generic else 'by\n'+
                  f'    norm_num [{p}, perpSq, normSq, dot, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]')]
        if case == 'parallel':
            setup += [f'  have he := (root_coincidence rho mu sigma {p} hr hm (ne_of_gt hs) hs1).mpr hq']
        for label in ('candidates','distinct'):
            lines += [f'theorem {label}_complete {args} (w : ℝ) :',
                      f'    Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma w {p}) = 0 ↔ w ∈ {label} {val} := by',
                      *(['  have hq : perpSq 0 '+p+' = 0 := by',
                         f'    norm_num [{p}, perpSq, normSq, dot, Fin.sum_univ_three, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]',
                         f'  have he := (root_coincidence rho mu sigma {p} hr hm (ne_of_gt hs) hs1).mpr hq'] if case=='parallel' else []),
                      f'  rw [rootPolynomial_complete rho mu sigma w {p} hr (ne_of_gt hs),',
                      f'    rootPolynomial_roots rho mu sigma {p} hr (ne_of_gt hs), {label}_reference {refcall}]',
                      *(['  rw [he]'] if case=='parallel' else []),
                      '  simp [referenceRoot]', '']
            audits += [ns+'.'+label+'_complete']
        operands = bridge.build(engine,bridge.INPUTS[engine]) if generic else reruns.build(bridge,engine,bridge.INPUTS[engine])
        determinant = next(r for r in operands.records if r['tag'].endswith('_'+pre+'Q3_DETERMINANT'))['cells'][0]
        operand_ns = engine if generic else engine+'Rerun'
        lines += [f'theorem emitted_determinant_complete {args} (w : ℝ) :',
                  f'    Expr.eval (values rho mu sigma w k) {operand_ns}.{determinant["lean_tree"]} = 0 ↔ w ∈ distinct {val} := by',
                  f'  rw [{operand_ns}.{determinant["lean_claim"]} rho mu sigma w k]',
                  f'  exact distinct_complete {call} w', '']
        records['distinct']['root_semantics']['lean_semantic_claims'] += [ns+'.emitted_determinant_complete']
        audits += [ns+'.emitted_determinant_complete']
        lines += [f'theorem distinct_nodup {args} : (distinct {val}).Nodup := by',
                  *setup[:2], f'  rw [distinct_reference {refcall}]',
                  f'  exact {"parallel" if case=="parallel" else "split"}_root_nodup rho mu sigma {p} hr hm hs'+
                  (' hk' if case=='parallel' else ' hs1 hk hq'), '']
        audits += [ns+'.distinct_nodup']
        lines += [f'theorem multiplicities {args} :',
                  f'    rootMultiplicity 0 (rootPolynomial rho mu sigma {p}) = 1 ∧',
                  f'      rootMultiplicity (referenceRoot rho mu sigma {p} 1) (rootPolynomial rho mu sigma {p}) = '+('2' if case=='parallel' else '1 ∧'),
                  *([] if case=='parallel' else [f'      rootMultiplicity (referenceRoot rho mu sigma {p} 2) (rootPolynomial rho mu sigma {p}) = 1']),
                  '    := by', *setup[:2],
                  f'  exact {"parallel" if case=="parallel" else "split"}_root_multiplicities rho mu sigma {p} hr hm hs hs1 hk hq', '']
        audits += [ns+'.multiplicities']
        if not (case=='parallel' and engine=='PY'):
            lines += [f'theorem candidate_multiset {args} :',
                      f'    (rootPolynomial rho mu sigma {p}).roots = (candidates {val} : Multiset ℝ) := by',
                      *(setup[1:] if case=='parallel' else []),
                      f'  rw [rootPolynomial_roots rho mu sigma {p} hr (ne_of_gt hs), candidates_reference {refcall}]',
                      *(['  rw [he]'] if case=='parallel' else []), '  rfl', '']
            audits += [ns+'.candidate_multiset']
        # An actual finite-set cardinality, rather than length alone.
        cell = records['distinct_count']['cells'][0]
        lines += [f'theorem distinct_count {args} :',
                  f'    Expr.eval (values {val}) {cell["lean_tree"]} = ((distinct {val}).toFinset.card : ℝ) := by',
                  '  classical', f'  rw [List.toFinset_card_of_nodup (distinct_nodup {call}),',
                  f'    {cell["lean_claim"]} {val}, distinct_reference {refcall}]', '  norm_num', '']
        audits += [ns+'.distinct_count']
        if engine=='WL':
            lines += ['def filteredTrees : List (Expr Symbol) := candidatesTrees.filter frequencyFree',
                      f'def filtered {params} : List ℝ := filteredTrees.map (Expr.eval (values {val}))', '',
                      'theorem filter_keeps_all : filteredTrees = candidatesTrees := by rfl', '',
                      'theorem discarded_empty : candidatesTrees.filter (fun e => !frequencyFree e) = [] := by rfl', '',
                      f'theorem filtered_reference {params} : filtered {val} = candidates {val} := by',
                      '  rw [filtered, filter_keeps_all]; rfl', '']
            audits += [ns+'.filter_keeps_all', ns+'.discarded_empty', ns+'.filtered_reference']
            tag = bridge.PREFIX+pre+'Q3_ROOTS_DISCARDED_BY_FILTER'
            line, raw = b.rows[tag]
            bridge.require(bridge.Parser(raw,engine).parse() == [], 'expected empty discarded-root list')
            controls.append({'tag':engine+'_'+tag,'line':line,'payload':raw,'payload_sha256':bridge.sha(raw.encode()),
                'case':case,'meaning':'candidate trees containing the squared-frequency symbol',
                'lean_semantic_claims':[ns+'.discarded_empty',ns+'.filter_keeps_all']})
            for label, operand in [('candidate_count','candidates'),('filtered_count','filtered')]:
                c = records[label]['cells'][0]
                lines += [f'theorem {label} {params} :',
                          f'    Expr.eval (values {val}) {c["lean_tree"]} = (({operand} {val}).length : ℝ) := by',
                          f'  rw [{c["lean_claim"]} {val}]',
                          *(['  rw [filtered_reference]'] if operand=='filtered' else []), '  rfl', '']
                audits += [ns+'.'+label]
            a,c = records['list_counts']['cells']
            lines += [f'theorem list_counts {args} :',
                      f'    Expr.eval (values {val}) {a["lean_tree"]} = ((filtered {val}).length : ℝ) ∧',
                      f'      Expr.eval (values {val}) {c["lean_tree"]} = ((distinct {val}).toFinset.card : ℝ) := by',
                      '  classical', '  constructor',
                      f'  · rw [{a["lean_claim"]} {val}, filtered_reference]; rfl',
                      f'  · rw [List.toFinset_card_of_nodup (distinct_nodup {call}),',
                      f'      {c["lean_claim"]} {val}, distinct_reference {refcall}]', '    norm_num', '']
            audits += [ns+'.list_counts']
        lines += ['end '+ns, '']
    source = '\n'.join(lines+['end','end S10Audit.CAS',''])
    source = re.sub(r'\b(hr|hm|hs|hs1|hg|z|k)\b', r'_\1', source)
    return source, audits, controls


def generate(bridge):
    outputs, builders, engines, audits = {}, {}, {}, []
    for engine,path in bridge.INPUTS.items():
        for case in CASES:
            b = build(bridge,engine,path,case)
            builders[engine,case] = b
            outputs[bridge.GENERATED/(b.namespace+'.lean')] = b.output()
            engines[b.namespace] = {'engine':engine,'case':case,'path':str(path.relative_to(bridge.BASE)),
                'sha256':bridge.sha(path.read_bytes()),'records':b.records,'unique_trees':len(b.nodes)}
            audits += [b.namespace+'.'+name for name in b.audits]
    bindings, selected, controls = binding_module(bridge,builders)
    outputs[bridge.GENERATED/'RootBindings.lean'] = bindings
    audits += selected
    audits += ['rootPolynomial_eval','rootPolynomial_ne_zero','rootPolynomial_roots','rootPolynomial_complete',
               'rootPolynomial_card','root_nonzero','root_coincidence','split_root_nodup','parallel_root_nodup',
               'split_root_multiplicities','parallel_root_multiplicities','genericChart_split']
    return outputs, engines, controls, audits
