#!/usr/bin/env python3
"""Independent parser comparison and isolated, proof-rejecting bridge mutations.

Only the frozen, selected arithmetic is sent to the older CAS parsers. Malformed
payloads are handled solely by the non-evaluating bridge parser.
"""
from pathlib import Path
import json
import os
import re
import subprocess
import sys
import tempfile
from itertools import product

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import S10_lean_cas_bridge as bridge
import S10_cross_engine_comparator as comparator
import S10_lean_cas_minors as minors
import S10_lean_cas_loci as loci
import S10_lean_cas_reruns as reruns
import S10_lean_cas_counts as counts
import S10_lean_cas_generic_counts as generic_counts
import S10_lean_cas_roots as roots
import S10_lean_cas_coincidence as coincidence
import S10_lean_cas_records as records
import sympy as sp

SCRATCH = bridge.LEAN / 's10/_scratch/cas'


def flatten(value):
    if isinstance(value, (tuple, list, sp.Tuple)):
        return [x for child in value for x in flatten(child)]
    return [value]


def to_sympy(ast):
    op, *args = ast
    if op == 'num':
        return sp.Integer(args[0])
    if op == 'atom':
        return sp.Symbol(args[0])
    if op == 'pow':
        return to_sympy(args[0]) ** args[1]
    a, b = map(to_sympy, args)
    return {'add': lambda: a+b, 'sub': lambda: a-b, 'mul': lambda: a*b,
            'div': lambda: a/b}[op]()


def independent_comparison(builders):
    checked = 0
    for builder in builders.values():
        engine = builder.engine
        parse = comparator.parse_sympy_payload if engine == 'PY' else comparator.parse_wolfram_payload
        for record in builder.records:
            if record.get('record_semantics',{}).get('kind')=='aggregate':
                parsed=independent_coincidence_parse(record['payload'],engine)
                fields=[{str(rule.args[0].args[0]):rule.args[1] for rule in obj.args} for obj in parsed]
                parsed=[x for obj in fields for x in (*obj['Pair'],obj['Equation'].lhs)]
            elif builder.engine=='WL' and record.get('coincidence_arithmetic',{}).get('kind')=='equation':
                parsed = independent_coincidence_parse(record['payload'],engine)[0].lhs
            else:
                parsed = comparator.normalize(parse(record['payload']))
            kind = record.get('root_semantics',{}).get('kind')
            if kind == 'candidates':
                old = []
                for rule in parsed:
                    if engine == 'PY':
                        assert isinstance(rule,dict) and len(rule) == 1
                        key,value = next(iter(rule.items()))
                        assert str(key) == 'omegaSquared'
                    else:
                        assert len(rule) == 1 and rule[0].func.__name__ == 'Rule'
                        key,value = rule[0].args
                        assert str(key) == 'omegaSquared'
                    old.append(comparator.normalize(value))
            elif kind == 'list_counts':
                assert parsed.func.__name__ == 'Association' and len(parsed.args) == 2
                old = []
                for rule,key in zip(parsed.args,('FilteredCandidateCount','DistinctRootCount')):
                    assert rule.func.__name__ == 'Rule' and str(rule.args[0]) == 'Str('+key+')'
                    old.append(rule.args[1])
            else:
                old = flatten(parsed)
            assert len(old) == len(record['cells']), record['tag']
            for expression, cell in zip(old, record['cells']):
                expression = sp.sympify(expression)
                renamed = expression.xreplace({s: sp.Symbol(bridge.NAMES['WL'][s.name])
                                                for s in expression.free_symbols})
                assert sp.cancel(renamed - to_sympy(cell['ast'])) == 0, (record['tag'], cell['index'])
                checked += 1
    return checked


def must_reject(name, operation):
    try:
        operation()
    except bridge.BridgeError as error:
        return {'name': name, 'outcome': 'REJECTED', 'reason': str(error)}
    raise AssertionError('accepted malformed input: ' + name)


def parser_checks():
    checks = []
    for engine, raw in [('PY', '__import__("os")'), ('PY', 'sqrt(k1)'),
                        ('WL', 'Sqrt[k1]'), ('PY', 'unknown'), ('PY', '1.0'),
                        ('PY', 'k1**(-1)'), ('PY', 'k1**33'), ('PY', 'k1^2'),
                        ('WL', 'k1**2'), ('PY', 'k1 k2'), ('PY', '[k1]+[k2]'),
                        ('PY', 'Matrix([[k1],[k2,k3]])'), ('PY', '1;2'),
                        ('PY', '9'*31)]:
        checks.append(must_reject(engine+':'+raw, lambda e=engine, r=raw: bridge.Parser(r, e).parse()))
    for raw in ['1/0', 'mu_R/(2-2)', 'mu_R*(1/(3-3))']:
        checks.append(must_reject(raw, lambda r=raw: bridge.numeric(bridge.Parser(r, 'PY').parse())))
    for raw in ['rho_br*omegaSquared + mu_R*k1', 'k1-k1+1']:
        checks.append(must_reject(raw, lambda r=raw: bridge.dimension(bridge.Parser(r, 'PY').parse())))
    assert bridge.node_json(bridge.Parser('-k1**2', 'PY').parse()) == [
        'mul', ['num', -1], ['pow', ['atom', 'k0'], 2]]
    assert bridge.node_json(bridge.Parser('k1-k2-k3', 'PY').parse()) == [
        'sub', ['sub', ['atom', 'k0'], ['atom', 'k1']], ['atom', 'k2']]
    assert bridge.node_json(bridge.Parser('k1/k2/k3', 'WL').parse()) == [
        'div', ['div', ['atom', 'k0'], ['atom', 'k1']], ['atom', 'k2']]
    source = bridge.INPUTS['PY'].read_text()
    first = source.splitlines()[0]
    with tempfile.TemporaryDirectory(prefix='s10_bridge_') as directory:
        path = Path(directory)/'input.out'
        path.write_text(source+'\n'+first+'\n')
        checks.append(must_reject('duplicate_tag', lambda: bridge.read_rows(path, 'PY')))
        path.write_text(source.replace('PY_S10_XFORM_ANISO_D3_Q2_MATRIX_A:',
                                       'PY_S10_XFORM_ANISO_D3_REMOVED_MATRIX_A:', 1))
        checks.append(must_reject('missing_selected_record', lambda: bridge.build('PY', path)))
        path.write_text(source.replace('PY_', 'WL_', 1))
        checks.append(must_reject('wrong_engine', lambda: bridge.read_rows(path, 'PY')))
        tag = bridge.PREFIX+'Q8_STRATUM1_ROOT2_N6_NULLSPACE_BASIS'
        raw = bridge.read_rows(bridge.INPUTS['PY'],'PY')[tag][1]
        path.write_text(source.replace('PY_'+tag+': '+raw,
                                       'PY_'+tag+': (Matrix([[0], [1], [0]]),)',1))
        checks.append(must_reject('missing_exceptional_basis_vector', lambda: reruns.build(bridge,'PY',path)))
    for raw in ['1/2', 'rho_br', '33', '[1]']:
        checks.append(must_reject('invalid_count:'+raw,
                      lambda r=raw: counts.count_shape(bridge,bridge.Parser(r,'PY').parse())))
    return checks


def root_parser_checks():
    checks = []
    for engine, raw in [
        ('PY','[{wrong: 0}]'), ('PY','[{omegaSquared: 0, k1: 1}]'),
        ('PY','[{omegaSquared: 0}] trailing'), ('PY','[{omegaSquared: [0]}]'),
        ('PY','[]'), ('WL','{{omegaSquared -> ConditionalExpression[0, rhoBr > 0]}}'),
        ('WL','{{omegaSquared -> 0, omegaSquared -> 1}}'),
        ('WL','{{k1 -> 0}}'), ('WL','{omegaSquared -> 0}'),
        ('WL','{{omegaSquared -> 0.0}}')]:
        checks.append(must_reject('root_rules:'+engine+':'+raw,
                                 lambda e=engine,r=raw:roots.solutions(bridge,r,e)))
    for raw in ['<|"FilteredCandidateCount" -> 3, "DistinctRootCount" -> 3, "Extra" -> 1|>',
                '<|"DistinctRootCount" -> 3, "DistinctRootCount" -> 3|>',
                '<|"FilteredCandidateCount" -> 3, "DistinctRootCount" -> 1/2|>',
                '<|"FilteredCandidateCount" -> 3|>']:
        checks.append(must_reject('root_counts:'+raw,lambda r=raw:roots.list_counts(bridge,r,'WL')))
    return checks


def lean_run(name, source, expected_exit):
    target = SCRATCH / (name+'.lean')
    target.write_text(source)
    command = ['lake', 'env', 'lean', '-DwarningAsError=true', str(target.relative_to(bridge.LEAN))]
    result = subprocess.run(command, cwd=bridge.LEAN, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=240,
                            env={**os.environ, 'LAKE_CACHE_DIR': '.lake/cache'})
    (SCRATCH/(name+'.log')).write_text(result.stdout)
    assert result.returncode == expected_exit, (name, result.returncode, result.stdout)
    if expected_exit:
        assert re.search('unsolved goals|tactic|Type mismatch: After simplification',
                         result.stdout, re.I), (name, result.stdout)
        assert not re.search(r'unknown|unexpected|failed to synthesize|not found|unused', result.stdout, re.I), result.stdout
    return {'name': name, 'outcome': 'REJECTED' if expected_exit else 'PASS',
            'exit_status': result.returncode, 'source_sha256': bridge.sha(source.encode()),
            'log_sha256': bridge.sha(result.stdout.encode()), 'command': command,
            'diagnostic': result.stdout.strip()[:3000]}


def lean_mutations(builders):
    checks = []
    for engine, builder in builders.items():
        # Feed an altered transcript cell through the real parser and generator.
        suffix = 'Q2_MATRIX_B'
        tag = bridge.PREFIX+suffix
        line, raw = builder.rows[tag]
        # A dimensionless coefficient error must survive typing and fail equality.
        b = bridge.Builder(engine, bridge.INPUTS[engine])
        assert '/2' in raw
        b.rows[tag] = (line, raw.replace('/2', '/3', 1))
        b.record(suffix, lambda x: bridge.matrix_shape(x, 3, 3),
                 [(f'(1/2) * referenceMatrix rho mu sigma z k {i} {j}', (-3, -2, 1), [])
                  for i in range(3) for j in range(3)])
        checks.append(lean_run(engine+'_wrong_coefficient', b.output(), 1))
    # Value and unit checks must reject distinct errors; use positive controls
    # for each proof before removing an assumption or changing a unit.
    header = 'import S10Audit.CAS.Support\nopen S10Audit S10Audit.CAS S10Pilot S10Anisotropic\n'
    denominator = header + '''example (sigma : ℝ) (k : Vec 3) (h : BasisDomain sigma k 0) :
    referenceBasis sigma k 0 2 = 1 := by
  exact referenceBasis_last sigma k 0 h
'''
    checks.append(lean_run('denominator_control', denominator, 0))
    checks.append(lean_run('missing_chart_denominator', header + '''example (sigma : ℝ) (k : Vec 3) :
    referenceBasis sigma k 0 2 = 1 := by
  cas_equal
''', 1))
    checks.append(lean_run('wrong_basis_normalization', header + '''example (sigma : ℝ) (k : Vec 3)
    (h : BasisDomain sigma k 0) : referenceBasis sigma k 0 2 = 2 := by
  rw [referenceBasis_last sigma k 0 h]
  norm_num
''', 1))
    typed = header + '''example : Expr.HasDim units (.atom .mu) (dimensions (-1) (-2) 1) := by
  apply castDim (Expr.HasDim.atom (u := units) .mu)
  norm_num [units, dimensions, Dimension.Exponent.ringEquivRat, Dimension.Exponent.equivRat]
'''
    checks.append(lean_run('unit_control', typed, 0))
    checks.append(lean_run('wrong_slot_unit', typed.replace('dimensions (-1) (-2) 1',
                                                          'dimensions (-2) (-2) 1'), 1))
    return checks


def coordinate_checks():
    checks, evaluated, points = [], 0, 0
    def logical_ast(tree):
        op, *args = tree
        if op in ('eq0', 'lt0', 'gt0'):
            x = sp.Symbol('k'+str(args[0]+1))
            return {'eq0': lambda: sp.Eq(x,0), 'lt0': lambda: x<0, 'gt0': lambda: x>0}[op]()
        return (sp.And if op == 'and' else sp.Or)(*map(logical_ast,args))
    def old_rule(rule):
        if isinstance(rule, sp.Equality):
            return rule
        assert rule.func.__name__ == 'Rule'
        left,right = rule.args
        if right.func.__name__ == 'ConditionalExpression':
            return sp.And(sp.Eq(left,right.args[0]),right.args[1])
        return sp.Eq(left,right)
    for engine,path in bridge.INPUTS.items():
        rows = bridge.read_rows(path,engine)
        parse = comparator.parse_sympy_payload if engine == 'PY' else comparator.parse_wolfram_payload
        for family in minors.FAMILIES:
            tag = bridge.PREFIX+minors.family_data(family)[7].removesuffix('MINORS')+'LOCUS'
            raw = rows[tag][1]
            old = comparator.normalize(parse(raw))
            old_prop = sp.Or(*(sp.And(*map(old_rule,branch)) for branch in old))
            new_prop = logical_ast(loci.CoordinateParser(bridge,raw,engine).locus())
            # The restricted locus grammar compares coordinates only with zero;
            # these 27 sign assignments exhaust its distinct truth patterns.
            for values in product((-1,0,1),repeat=3):
                sub = {sp.Symbol('k'+str(i+1)):v for i,v in enumerate(values)}
                assert bool(old_prop.subs(sub)) == bool(new_prop.subs(sub)), (engine,tag,values)
                evaluated += 1
        for index in (1,2):
            suffix = f'Q8_STRATUM{index}_POINT' if engine == 'PY' else f'STRATUM{index}_Q8_POINT'
            raw = rows[bridge.PREFIX+suffix][1]
            old = comparator.normalize(parse(raw))
            expected = {str(rule.args[0]):rule.args[1] for rule in old}
            point = loci.CoordinateParser(bridge,raw,engine).point()
            assert all(sp.Rational(x.numerator,x.denominator) == expected['k'+str(i+1)]
                       for i,x in enumerate(point))
            points += 1
    for engine,raw,method in [
        ('PY','((Eq(k4, 0),),)','locus'),
        ('PY','((Eq(k1, 1),),)','locus'),
        ('PY','((Eq(k1, 0),),) trailing','locus'),
        ('PY','((),)','locus'),
        ('WL','{{k1 -> ConditionalExpression[0, k2 >= 0]}}','locus'),
        ('WL','{{k1 -> ConditionalExpression[0, unknown > 0]}}','locus'),
        ('WL','{{k1 -> ConditionalExpression[1, k2 > 0]}}','locus'),
        ('WL','{{k1 -> Unsafe[0]}}','locus'),
        ('WL','{k1 -> 1, k1 -> 0, k3 -> 0}','point'),
        ('PY','(Eq(k1, 1), Eq(k2, 0))','point'),
        ('WL','{k1 -> 1/0, k2 -> 0, k3 -> 0}','point'),
        ('WL','{k1 -> ConditionalExpression[0, k2 > 0], k2 -> 0, k3 -> 1}','point')]:
        checks.append(must_reject(engine+':'+raw,
                      lambda e=engine,r=raw,m=method:getattr(loci.CoordinateParser(bridge,r,e),m)()))
    return evaluated, points, checks


def minor_locus_mutations():
    checks = []
    ref = (bridge.GENERATED/'MinorReference.lean').read_text()
    # Check a row/column selection error, even though the set of polynomial
    # values and its dimensions are unchanged by this permutation.
    ref = ref.split('def staticTransverseMap')[0]+'\nend\nend S10Audit.CAS\n'
    assert ref.count('![![0, 1, 2],') == 1
    checks.append(lean_run('swapped_minor_selection',ref.replace('![![0, 1, 2],','![![0, 2, 1],',1),1))
    header = 'import S10Audit.CAS.MinorLoci\nopen S10Audit.CAS S10Pilot\n'
    for name,tree in [
        ('missing_perpendicular_branch',['or',['and',['eq0',1],['eq0',2]]]),
        ('missing_parallel_branch',['or',['and',['eq0',0]]])]:
        text = header + f'''def mutantLocus (k : Vec 3) : Prop := {loci.proposition(tree)}
example (k : Vec 3) : mutantLocus k ↔ (k 0 = 0 ∨ (k 1 = 0 ∧ k 2 = 0)) := by
  classical
  locus_eval mutantLocus
'''
        checks.append(lean_run(name,text,1))
    rows = bridge.read_rows(bridge.INPUTS['WL'],'WL')
    raw = rows[bridge.PREFIX+'ROOT3_Q8_TRANSVERSE_RANK_DROP_LOCUS'][1]
    assert 'k2 > 0 || k2 < 0' in raw
    mutated = raw.replace('k2 > 0 || k2 < 0','k2 > 0 && k2 < 0',1)
    tree = loci.CoordinateParser(bridge,mutated,'WL').locus()
    checks.append(lean_run('corrupt_conditional_guard',header+f'''def mutantLocus (k : Vec 3) : Prop := {loci.proposition(tree)}
example (k : Vec 3) : mutantLocus k ↔ (k 0 = 0 ∨ (k 1 = 0 ∧ k 2 = 0)) := by
  classical
  locus_eval mutantLocus
''',1))
    return checks



def rerun_mutations():
    checks = []
    header = ('import S10Audit.CAS.RerunReference\n'
              'set_option backward.isDefEq.respectTransparency false\n'
              'open S10Audit S10Audit.CAS S10Pilot S10Anisotropic\n')
    ref = (bridge.GENERATED/'RerunReference.lean').read_text()
    # Reuse the exact generated independence proof, with both basis entries
    # still present but made equal. Shape and membership alone cannot catch it.
    start = ref.index('def parallel1 : Fin 2')
    stop = ref.index('theorem parallel1_members',start)
    fragment = ref[start:stop]
    checks.append(lean_run('parallel_basis_control',header+fragment,0))
    assert fragment.count('!['+'![0, 1, 0], ![0, 0, 1]]') == 1
    duplicate = fragment.replace('!['+'![0, 1, 0], ![0, 0, 1]]','!['+'![0, 1, 0], ![0, 1, 0]]',1)
    checks.append(lean_run('dependent_parallel_basis',header+duplicate,1))
    # The full-kernel dimension is two: a one-dimensional replacement must
    # fail even when its retained vector is a valid null vector.
    incomplete = header + """example (sigma : ℝ) (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    Module.finrank ℝ (rootMode sigma PYLoci.parallelPoint 1) = 1 := by
  rw [(PYReference.parallel1_dimensions sigma hs hs1).1]
  norm_num
"""
    checks.append(lean_run('incomplete_parallel_kernel',incomplete,1))
    support = (bridge.GENERATED/'RerunSupport.lean').read_text()
    start = support.index('theorem referenceMatrix_scale')
    stop = support.index('theorem referenceRoot_scale')
    scaling = support[start:stop]
    checks.append(lean_run('coordinate_scale_control',header+scaling,0))
    checks.append(lean_run('wrong_frequency_scale',header+scaling.replace('(c ^ 2 * z)', '(c * z)',1),1))
    start = support.index('theorem coordinate_frequency_units')
    stop = support.index('theorem coordinate_matrix_units')
    units = support[start:stop].replace('^ (2 : ℕ)','^ (1 : ℕ)',1)
    checks.append(lean_run('wrong_coordinate_unit_scale',header+units,1))
    return checks


def count_mutations():
    checks = []
    for engine, suffix, replacement, expected in [
        ('PY', 'Q8_STRATUM1_ROOT2_N2_NULLITY', '1', 2),
        ('WL', 'STRATUM1_ROOT1_N3_TRANSVERSE_NULLITY', '1', 0)]:
        b = bridge.Builder(engine,bridge.INPUTS[engine],namespace=engine+'Counts',
                           support='CountReference',tactic='count_eval')
        tag = bridge.PREFIX+suffix
        line,raw = b.rows[tag]
        assert raw == str(expected)
        b.rows[tag] = (line,replacement)
        b.record(suffix,lambda x:counts.count_shape(bridge,x),[(str(expected),(0,0,0),[])])
        checks.append(lean_run(engine+'_wrong_count',b.output(),1))
    header = ('import S10Audit.CAS.CountReference\n'
              'set_option backward.isDefEq.respectTransparency false\n'
              'open S10Audit S10Audit.CAS S10Pilot\n')
    args = '(rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1)'
    call = 'rho mu sigma hr hm hs hs1'
    rank = header+f"""example {args} :
    (PYRerun.parallel0Stack rho mu sigma).rank + matrixNullity (PYRerun.parallel0Stack rho mu sigma) = 4 := by
  rw [PYCountReference.parallel0_stacked_rank {call},
    PYCountReference.parallel0_transverse_nullity {call}]
  norm_num
"""
    checks.append(lean_run('wrong_stack_domain_dimension',rank,1))
    signed = header+f"""example {args} :
    basisCountResidual (fun _ : Fin 1 => PYReference.parallel1 0)
      (PYRerun.parallel1Matrix rho mu sigma) = (-1 : ℤ) := by
  norm_num [basisCountResidual, basisCount, PYCountReference.parallel1_nullity {call}]
"""
    checks.append(lean_run('signed_count_residual_control',signed,0))
    checks.append(lean_run('truncated_count_residual',signed.replace('= (-1 : ℤ)','= (0 : ℤ)',1),1))
    return checks

def generic_count_mutations():
    checks = []
    for engine, suffix, replacement, expected in [
        ('PY', 'ROOT1_N2_NULLITY', '2', 1),
        ('WL', 'ROOT3_N3_TRANSVERSE_NULLITY', '1', 0)]:
        b = bridge.Builder(engine, bridge.INPUTS[engine], namespace=engine+'GenericCounts',
                           support='GenericCountReference', tactic='count_eval')
        tag = bridge.PREFIX+suffix
        line, raw = b.rows[tag]
        assert raw == str(expected)
        b.rows[tag] = (line, replacement)
        b.record(suffix, lambda x: counts.count_shape(bridge, x), [(str(expected),(0,0,0),[])])
        checks.append(lean_run(engine+'_wrong_generic_count', b.output(), 1))
    # The generic matrix formulas still make sense at this exceptional point,
    # but the generic extra-root transverse count does not extend there.
    header = ('import S10Audit.CAS.GenericCountReference\n'
              'import S10Audit.CAS.CountReference\n'
              'set_option backward.isDefEq.respectTransparency false\n'
              'open S10Audit S10Audit.CAS S10Pilot\n')
    boundary = header+"""example (rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0)
    (hs : 0 < sigma) (hs1 : sigma ≠ 1) :
    matrixNullity (PYGenericCountReference.root2Stack rho mu sigma 0 PYLoci.perpendicularPoint) = 1 := by
  have outside : ¬ GenericChart sigma PYLoci.perpendicularPoint := by
    intro h
    exact h.2.2.1 (by norm_num [PYLoci.perpendicularPoint])
  have same : PYGenericCountReference.root2Stack rho mu sigma 0 PYLoci.perpendicularPoint =
      PYRerun.perpendicular2Stack rho mu sigma := by
    rw [PYGenericCountReference.root2Stack_reference rho mu sigma 0 PYLoci.perpendicularPoint hr (ne_of_gt hs),
      PYRerun.perpendicular2Stack_reference rho mu sigma hr (ne_of_gt hs)]
    rfl
  norm_num [same, PYCountReference.perpendicular2_transverse_nullity rho mu sigma hr hm hs hs1]
"""
    checks.append(lean_run('generic_perpendicular_boundary_control', boundary, 0))
    checks.append(lean_run('generic_count_beyond_chart', boundary.replace('perpendicularPoint) = 1',
                                                                           'perpendicularPoint) = 0', 1), 1))
    return checks


def root_mutations():
    checks = []
    header = ('import S10Audit.CAS.RootBindings\n'
              'set_option backward.isDefEq.respectTransparency false\n'
              'open S10Audit S10Audit.CAS S10Pilot Polynomial\n')
    args = '(rho mu sigma : ℝ) (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : 0 < sigma) (hs1 : sigma ≠ 1)'
    multiplicity = header+f"""example {args} :
    rootMultiplicity (referenceRoot rho mu sigma WLLoci.parallelPoint 1)
      (rootPolynomial rho mu sigma WLLoci.parallelPoint) = 2 := by
  rw [(WLRootsParallel.multiplicities rho mu sigma 0 0 hr hm hs hs1).2]
"""
    checks.append(lean_run('root_multiplicity_control', multiplicity, 0))
    wrong = multiplicity.replace('parallelPoint) = 2','parallelPoint) = 1',1)+'  norm_num\n'
    checks.append(lean_run('collapsed_parallel_multiplicity', wrong, 1))
    incomplete = header+"""example (rho mu sigma : ℝ) (k : Vec 3)
    (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    (rootPolynomial rho mu sigma k).roots.card = 2 := by
  rw [rootPolynomial_card rho mu sigma k hr hs]
  norm_num
"""
    checks.append(lean_run('incomplete_cubic_root_multiset', incomplete, 1))
    b = bridge.Builder('WL',bridge.INPUTS['WL'],namespace='WrongRootCount',support='RootSupport',tactic='count_eval')
    suffix = 'STRATUM1_Q3_ROOT_COUNT'
    line,raw = b.rows[bridge.PREFIX+suffix]
    assert raw == '2'
    b.rows[bridge.PREFIX+suffix] = (line,'3')
    b.record(suffix,lambda x:counts.count_shape(bridge,x),[('2',(0,0,0),[])])
    checks.append(lean_run('confused_candidate_and_distinct_count',b.output(),1))
    b = roots.build(bridge,'WL',bridge.INPUTS['WL'],'generic')
    zero = b.records[0]['cells'][0]
    assert zero['literal_zero_lift'] and zero['lean_raw_tree'] != zero['lean_tree']
    control = header+f"""example : frequencyFree WLRootsGeneric.{zero['lean_raw_tree']} = true ∧
    frequencyFree (.sub (.atom .z) (.atom .z)) = false := by decide
"""
    checks.append(lean_run('raw_root_filter_control',control,0))
    checks.append(lean_run('unit_adapted_zero_filter',control.replace(zero['lean_raw_tree'],zero['lean_tree'],1),1))
    return checks


def coincidence_symbol(name):
    name={'rho_br':'rhoBr','mu_R':'muR','s_rho':'sRho'}.get(name,name)
    return {'True':sp.true,'False':sp.false,'EmptySet':sp.EmptySet}.get(name,sp.Symbol(comparator.mechanical_lower_camel(name)))


def coincidence_call(head,args):
    operations={'Plus':lambda:sp.Add(*args),'Times':lambda:sp.Mul(*args),'Power':lambda:args[0]**args[1],
                'Equal':lambda:sp.Eq(*args),'Unequal':lambda:sp.Ne(*args),'Less':lambda:sp.Lt(*args),
                'Greater':lambda:sp.Gt(*args),'LessEqual':lambda:sp.Le(*args),'GreaterEqual':lambda:sp.Ge(*args),
                'And':lambda:sp.And(*args),'Or':lambda:sp.Or(*args),'Sqrt':lambda:sp.sqrt(*args),
                'Q.integer':lambda:sp.Q.integer(*args),'Q.real':lambda:sp.Q.real(*args),'Q.positive':lambda:sp.Q.positive(*args)}
    if head in operations: return operations[head]()
    if head=='List': return tuple(args)
    if head=='Element': return sp.Symbol('Element('+','.join(map(str,args))+')')
    if head=='Inequality':
        assert len(args)==5 and args[1]==args[3]==sp.Symbol('Less')
        return sp.And(sp.Lt(args[0],args[2]),sp.Lt(args[2],args[4]))
    assert head in ('ConditionalExpression','Rule','Association','_Str'),head
    return sp.Function(head)(*args)


def coincidence_to_sympy(t):
    op,*args=t
    if op=='num': return sp.Integer(args[0])
    if op=='symbol': return coincidence_symbol(args[0])
    if op=='seq': return tuple(map(coincidence_to_sympy,args))
    if op=='assoc': return coincidence_call('Association',list(map(coincidence_to_sympy,args)))
    if op=='key': return coincidence_call('_Str',[sp.Symbol(args[0])])
    if op=='call':
        head={'Eq':'Equal','Ne':'Unequal'}.get(args[0],args[0])
        return coincidence_call(head,list(map(coincidence_to_sympy,args[1:])))
    a=coincidence_to_sympy(args[0])
    if op=='neg': return -a
    b=coincidence_to_sympy(args[1])
    if op=='sub': return a-b
    if op=='div': return a/b
    return coincidence_call({'add':'Plus','mul':'Times','pow':'Power','eq':'Equal','ne':'Unequal',
        'lt':'Less','gt':'Greater','le':'LessEqual','ge':'GreaterEqual','and':'And','or':'Or','rule':'Rule'}[op],[a,b])


def independent_coincidence_parse(raw,engine):
    if engine=='WL' and raw=='{}': return comparator.parse_wolfram_payload(raw)
    if engine=='PY': return comparator.normalize(comparator.parse_sympy_payload(raw))
    # Use SymPy's independent Mathematica tokenizer and full-form parser. Its final
    # Boolean conversion cannot accept Element inside And; retain membership as
    # an opaque Boolean atom here instead of changing the production comparator.
    from sympy.parsing.mathematica import MathematicaParser
    parser=MathematicaParser()
    def walk(t):
        if isinstance(t,str):
            return sp.Integer(t) if re.fullmatch(r'-?[0-9]+',t) else coincidence_symbol(t)
        return coincidence_call(t[0],[walk(x) for x in t[1:]])
    prepared=raw.replace('<|','Association[').replace('|>',']')
    return walk(parser._from_tokens_to_fullformlist(parser._from_mathematica_to_tokens(prepared)))


def coincidence_parser_checks():
    _,_,records,_=coincidence.generate(bridge)
    for r in records:
        def evaluate(t):
            if isinstance(t,(tuple,list,sp.Tuple)): return tuple(map(evaluate,t))
            return t.doit() if isinstance(t,sp.Basic) else t
        old=evaluate(independent_coincidence_parse(r['payload'],r['tag'][:2]))
        fresh=coincidence_to_sympy(r['mixed_ast'])
        assert old==fresh,(r['tag'],old,fresh)
    malformed=[]
    for engine,raw in [('PY','__import__("os")'),('WL','Run[1]'),('WL','k1 -> ConditionalExpression[0]'),
                       ('WL','k1 -> ConditionalExpression[0, Element[k1,Integers]]'),
                       ('PY','Q.real(unknown)'),('WL','k1 -> Sqrt[k2, k3]'),('PY','Eq(k1,0) junk'),
                       ('WL','k1 -> ConditionalExpression[0, sRho < 0]; 1'),
                       ('PY','Q.positive(k1**33)'),('WL','k1 -> ConditionalExpression[0,Maybe]'),
                       ('PY','Q.positive(1/0)'),('WL','k1 -> Sqrt[1/k2]'),
                       ('PY','True && False'),('WL','True & False')]:
        def check(e=engine,raw=raw):
            t=coincidence.parse(bridge,raw,e)
            if t[0]=='rule': coincidence.assignment(bridge,t)
            else: coincidence.prop(bridge,t)
        malformed.append(must_reject('coincidence:'+raw,check))
    return len(records),malformed


def coincidence_mutations():
    source=(bridge.GENERATED/'CoincidenceBindings.lean').read_text()
    header='import S10Audit.CAS.CoincidenceBindings\nnamespace S10Audit.CAS.Coincidence\nopen S10Pilot S10Anisotropic\nnoncomputable section\nset_option backward.isDefEq.respectTransparency false\n'
    checks=[]
    def block(name):
        start=source.index('def '+name+' ')
        end=source.index('\ndef ',start+4) if '\ndef ' in source[start+4:] else source.index('\nend\n',start)
        text=source[start:end]
        return re.sub(r'\b'+re.escape(name)+r'\b','probe',text).replace(name+'_correct','probe_correct').replace(name+'_sound','probe_sound').replace(name+'_exists','probe_exists')
    controls=[('guarded_locus','WL_generic_locus1'),('allowed_banks','WL_generic_allowed2'),('printed_witness','PY_generic_witness2')]
    for label,name in controls:
        text=block(name)
        checks.append(lean_run(label+'_control',header+text,0))
    guards=block('WL_generic_locus1')
    assert '(_sigma < 0)' in guards
    checks.append(lean_run('admit_negative_sigma_branches',header+guards.replace('(_sigma < 0)','(_sigma > 0)'),1))
    banks=block('WL_generic_allowed2')
    assert '((_k 0 > 0)' in banks
    checks.append(lean_run('omit_positive_parallel_bank',header+banks.replace('((_k 0 > 0)','((_k 0 < 0)',1),1))
    status=block('PY_generic_status2')
    checks.append(lean_run('coincidence_wrong_nonempty_status',header+status.replace(': Prop := True',': Prop := False',1),1))
    witness=block('PY_generic_witness2')
    checks.append(lean_run('coincidence_zero_witness',header+witness.replace('(_k 0 = 1)','(_k 0 = 0)',1),1))
    geom=block('PY_generic_locus2')
    checks.append(lean_run('coincidence_omit_axis_constraint',header+geom.replace('((_k 1 = 0) ∧ (_k 2 = 0))','(_k 1 = 0)',1),1))
    # A scalar mutation must still pass the grammar and fail the root-difference proof.
    b=coincidence.build(bridge,'PY',bridge.INPUTS['PY'],'generic')
    r=b.records[0];c=r['cells'][0]
    claim=f'''import S10Audit.CAS.CoincidenceBindings
namespace S10Audit.CAS
open S10Pilot S10Anisotropic
noncomputable section
example (rho mu sigma z : ℝ) (k : Vec 3) (hr : rho ≠ 0) (hs : sigma ≠ 0) :
    Expr.eval (values rho mu sigma z k) {b.namespace}.{c['lean_tree']} =
      referenceRoot rho mu sigma k 1 - referenceRoot rho mu sigma k 0 := by
  cas_eval {b.namespace}.{c['lean_tree']}_eval
'''
    checks.append(lean_run('coincidence_reversed_difference',claim,1))
    return checks


def record_parser_checks():
    _,_,data,_=records.generate(bridge)
    def evaluate(t):
        if isinstance(t,(tuple,list,sp.Tuple)): return tuple(map(evaluate,t))
        return t.doit() if isinstance(t,sp.Basic) else t
    for rec in data:
        old=evaluate(independent_coincidence_parse(rec['payload'],rec['tag'][:2]))
        if 'mixed_ast' in rec:
            fresh=coincidence_to_sympy(rec['mixed_ast'])
        elif rec['kind']=='root_sign':
            fresh={'zero':sp.Integer(0),'positive':sp.Integer(1),'negative':sp.Integer(-1),
                   'undecided':sp.Symbol('undecidedUnderJointAssumptions')}[rec['reported_sign']]
        elif rec['kind']=='skipped_decision':
            fresh=sp.false
        else:
            names={'spectrum_status':'rootsReturned','stratum_status':'notSkippedAllowedBranch',
                   'skipped_reason':'locusConflictsWithPositiveWavevectorNorm'}
            fresh=sp.Symbol(names[rec['kind']])
        assert old==fresh,(rec['tag'],old,fresh)
    malformed=[]
    def agg(raw): return records.associations(bridge,coincidence.parse(bridge,raw,'WL'),'generic')
    raw=bridge.read_rows(bridge.INPUTS['WL'],'WL')[bridge.PREFIX+'Q3_ROOT_COINCIDENCE_LOCI'][1]
    for name,bad in [
        ('duplicate_key',raw.replace('"Equation" ->','"Pair" ->',1)),
        ('missing_key',raw.replace('"Intersects" -> False','',1).replace(', |>','|>')),
        ('noninteger_pair',raw.replace('"Pair" -> {1, 2}','"Pair" -> {1/2, 2}',1)),
        ('out_of_range_pair',raw.replace('"Pair" -> {1, 2}','"Pair" -> {1, 4}',1)),
        ('nonzero_equation_rhs',raw.replace(' == 0',' == 1',1)),
        ('unknown_outcome',raw.replace('decidedEmpty','unresolved',1)),
        ('missing_pair', '{'+raw[raw.index('<|',raw.index('|>')+2):]),
    ]:
        malformed.append(must_reject('aggregate:'+name,lambda bad=bad:agg(bad)))
    for raw in ('2','1.0','positive','undecided','0;1'):
        malformed.append(must_reject('sign:'+raw,lambda raw=raw:records.reported_sign(bridge,raw)))
    # Exact solve-variable extraction must reject a same-unit composite expression.
    with tempfile.TemporaryDirectory(prefix='s10_record_') as directory:
        path=Path(directory)/'input.out'
        source=bridge.INPUTS['PY'].read_text()
        tag=bridge.PREFIX+'Q3_SPECTRUM_SOLVE_CONDITION_OPERANDS'
        payload=bridge.read_rows(bridge.INPUTS['PY'],'PY')[tag][1]
        replacement=payload.removesuffix(', omegaSquared)')+', 2*omegaSquared)'
        assert replacement!=payload
        path.write_text(source.replace(payload,replacement,1))
        malformed.append(must_reject('composite_solve_variable',lambda:records.build(bridge,'PY',path,'generic')))
    return len(data),malformed


def record_mutations():
    record_source=(bridge.GENERATED/'RecordBindings.lean').read_text()
    status_source=(bridge.GENERATED/'StatusBindings.lean').read_text()
    header='import S10Audit.CAS.RecordBindings\nimport S10Audit.CAS.StatusBindings\nnamespace S10Audit.CAS.Records\nopen S10Pilot S10Anisotropic Polynomial\nnoncomputable section\nset_option backward.isDefEq.respectTransparency false\n'
    checks=[]
    def block(source,name):
        start=source.index('def '+name+' ')
        end=source.index('\ndef ',start+4) if '\ndef ' in source[start+4:] else source.index('\nend\n',start)
        text=source[start:end]
        return text.replace(name,'probe')
    alias=block(record_source,'PY_q8root3_status2')
    checks.append(lean_run('q8_duplicate_control',header+alias,0))
    wrong=alias.replace(': Prop := True',': Prop := False',1)
    assert wrong!=alias
    # Change both the copied field and its proposed reference proof: the wrong
    # status must fail the universal semantic theorem, not just an identifier lookup.
    wrong=wrong.replace('  rfl','  simp [probe, Coincidence.PY_generic_status2]',1)
    checks.append(lean_run('q8_duplicate_wrong_status',header+wrong,1))
    pairs=block(record_source,'WL_generic_aggregate0_pairs')
    wrong=pairs.replace('[(1, 2), (1, 3), (2, 3)]','[(1, 2), (1, 3), (1, 3)]',1)
    checks.append(lean_run('aggregate_wrong_pair_index',header+wrong,1))
    # One changed coefficient passes the mixed parser and unit checker but fails equality.
    with tempfile.TemporaryDirectory(prefix='s10_record_') as directory:
        path=Path(directory)/'input.out'
        source=bridge.INPUTS['WL'].read_text()
        tag=bridge.PREFIX+'Q3_ROOT_COINCIDENCE_LOCI'
        payload=bridge.read_rows(bridge.INPUTS['WL'],'WL')[tag][1]
        changed=payload.replace('"Equation" -> -(', '"Equation" -> -2*(',1)
        assert changed!=payload
        path.write_text(source.replace(payload,changed,1))
        b=records.build(bridge,'WL',path,'generic')
        checks.append(lean_run('aggregate_wrong_root_difference',b.output(),1))
    sign=block(status_source,'PY_generic_root2_sign')
    checks.append(lean_run('undecided_sign_with_positive_proof_control',header+sign,0))
    wrong=sign.replace('(classifiedSign 2).Holds','ReportedSign.negative.Holds',1)
    # Ask the same algebraic proof to certify the changed sign, allowing simplification
    # to expose a mathematical goal rather than treating elaboration as the rejection.
    line='  exact reference_root_sign _rho _mu _sigma _k 2 hr hm hs _hk'
    assert line in wrong
    wrong=wrong.replace(line,'  have hp := reference_root_sign _rho _mu _sigma _k 2 hr hm hs _hk\n  simp only [classifiedSign, ReportedSign.Holds] at hp ⊢\n  norm_num at hp ⊢\n  try linarith',1)
    checks.append(lean_run('wrong_extra_root_sign',header+wrong,1))
    # Removing the nonzero-wavevector premise admits k=0 and must fail positivity.
    root_positive=header+'''example (rho mu sigma : ℝ) (k : Vec 3)
    (hr : 0 < rho) (hm : 0 < mu) (hs : 0 < sigma) :
    0 < referenceRoot rho mu sigma k 2 := by
  change 0 < (mu / rho) * extraValue 0 sigma k
  apply mul_pos (div_pos hm hr)
  apply extraValue_pos hs
  simp
'''
    checks.append(lean_run('positive_root_without_nonzero_wavevector',root_positive,1))
    skipped=header+'''example : PY_skipped_test = true := by
  simp [PY_skipped_test]
'''
    checks.append(lean_run('admitting_skipped_origin',skipped,1))
    solver=header+'''example (rho mu sigma lam z : ℝ) (k : Vec 3)
    (hd : coincidenceDomain rho mu sigma lam) :
    Expr.eval (values rho mu sigma z k) PYRecordsParallel.n16 = 0 ↔
      z ∈ PYRootsParallel.candidates rho mu sigma z k := by
  exact PY_parallel_spectrum_returned_complete rho mu sigma lam z k hd
'''
    # Resolve the actual imported determinant tree; do not hard-code its interned id.
    b=records.build(bridge,'PY',bridge.INPUTS['PY'],'parallel')
    actual=b.records[0]['cells'][0]['lean_tree']
    solver=solver.replace('PYRecordsParallel.n16','PYRecordsParallel.'+actual)
    checks.append(lean_run('solver_operand_complete_control',solver,0))
    return checks


def main():
    SCRATCH.mkdir(parents=True, exist_ok=True)
    tracked = [*bridge.INPUTS.values(), *bridge.GENERATED.glob('*.lean'),
               Path(bridge.__file__),Path(minors.__file__),Path(loci.__file__),Path(reruns.__file__),Path(counts.__file__),
               Path(generic_counts.__file__),Path(roots.__file__),Path(coincidence.__file__),Path(records.__file__),ROOT/'scripts/S10_exports.py']
    before = {str(p.relative_to(ROOT)): bridge.sha(p.read_bytes()) for p in tracked}
    bridge.generate(check=True)
    builders = {e: bridge.build(e, p) for e, p in bridge.INPUTS.items()}
    minor_builders = {e+'Minors':minors.build(bridge,e,p) for e,p in bridge.INPUTS.items()}
    rerun_builders = {e+'Rerun':reruns.build(bridge,e,p) for e,p in bridge.INPUTS.items()}
    count_builders = {e+'Counts':counts.build(bridge,e,p) for e,p in bridge.INPUTS.items()}
    generic_count_builders = {e+'GenericCounts':generic_counts.build(bridge,e,p) for e,p in bridge.INPUTS.items()}
    root_builders = {e+case:roots.build(bridge,e,p,case) for e,p in bridge.INPUTS.items() for case in roots.CASES}
    coincidence_builders = {e+'Coincidence'+case:coincidence.build(bridge,e,p,case) for e,p in bridge.INPUTS.items() for case in coincidence.CASES}
    record_builders = {e+'Records'+case:records.build(bridge,e,p,case) for e,p in bridge.INPUTS.items() for case in records.CASES}
    metadata_records,metadata_rejections = record_parser_checks()
    coincidence_records,coincidence_rejections = coincidence_parser_checks()
    cells = independent_comparison({**builders,**minor_builders,**rerun_builders,**count_builders,**generic_count_builders,**root_builders,**coincidence_builders,**record_builders})
    malformed = parser_checks()
    malformed += root_parser_checks()
    malformed += coincidence_rejections
    malformed += metadata_rejections
    _,_,filters,_ = roots.generate(bridge)
    for record in filters:
        assert comparator.normalize(comparator.parse_wolfram_payload(record['payload'])) == ()
    coordinate_patterns,points,coordinate_rejections = coordinate_checks()
    mutations = lean_mutations(builders)
    mutations += minor_locus_mutations()
    mutations += rerun_mutations()
    mutations += count_mutations()
    mutations += generic_count_mutations()
    mutations += root_mutations()
    mutations += coincidence_mutations()
    mutations += record_mutations()
    after = {str(p.relative_to(ROOT)): bridge.sha(p.read_bytes()) for p in tracked}
    assert before == after
    bridge.generate(check=True)
    report = {'status': 'PASS', 'independent_parser_cells': cells,
              'strict_parser_rejections': malformed, 'lean_checks': mutations,
              'coordinate_sign_pattern_comparisons': coordinate_patterns,
              'independently_parsed_points': points,'coordinate_parser_rejections':coordinate_rejections,
              'independently_parsed_root_filters':len(filters),
              'independently_parsed_coincidence_records':coincidence_records,
              'independently_parsed_metadata_records':metadata_records,
              'canonical_sha256_before': before, 'canonical_sha256_after': after,
              'instrument_sha256': bridge.sha(Path(__file__).read_bytes())}
    (ROOT/'_measurements/S10_lean_cas_bridge_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(f'PASS: {cells} independently parsed cells, {len(malformed)} parser rejections, '
          f'{coordinate_patterns} locus sign patterns, {points} points, '
          f'{len(coordinate_rejections)} coordinate-parser rejections, '
          f'{len(mutations)} Lean controls/mutations; canonical hashes unchanged.')


if __name__ == '__main__':
    main()
