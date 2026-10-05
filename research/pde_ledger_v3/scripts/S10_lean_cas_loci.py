"""Strict parsing of emitted coordinate loci and targeted stratum points.

The grammar preserves unions, conjunctions and Wolfram conditional-rule guards.
It accepts only the coordinate/rational forms used by this D3 certificate.
"""
import re


class CoordinateParser:
    token = re.compile(r'\s+|->|\|\||&&|[A-Za-z][A-Za-z0-9]*|[0-9]+|[(),{}\[\]<>/+-]')

    def __init__(self, bridge, raw, engine):
        self.b, self.engine, self.tokens, self.pos = bridge, engine, [], 0
        bridge.require(engine in ('PY', 'WL') and len(raw) < 100000, 'coordinate parser domain')
        end = 0
        for match in self.token.finditer(raw):
            bridge.require(match.start() == end, f'unsupported coordinate syntax at {end}')
            end = match.end()
            if not match.group().isspace():
                self.tokens.append(match.group())
        bridge.require(end == len(raw), 'unsupported coordinate trailing syntax')

    def peek(self):
        return self.tokens[self.pos] if self.pos < len(self.tokens) else None

    def take(self, expected=None):
        value = self.peek()
        self.b.require(value is not None and (expected is None or value == expected),
                       f'coordinate token: expected {expected}, got {value}')
        self.pos += 1
        return value

    def coordinate(self):
        name = self.take()
        self.b.require(name in ('k1', 'k2', 'k3'), 'expected wavevector coordinate')
        return int(name[1])-1

    def rational(self):
        parts = []
        if self.peek() in ('-', '+'):
            parts.append(self.take())
        number = self.take()
        self.b.require(number.isdigit() and len(number) < 20, 'expected bounded integer')
        parts.append(number)
        if self.peek() == '/':
            parts.append(self.take())
            denominator = self.take()
            self.b.require(denominator.isdigit() and len(denominator) < 20, 'expected denominator')
            parts.append(denominator)
        node = self.b.Parser(''.join(parts), self.engine).parse()
        return self.b.numeric(node)

    def condition_atom(self):
        if self.peek() == '(':
            self.take('(')
            term = self.condition()
            self.take(')')
            return term
        coord = self.coordinate()
        op = self.take()
        self.b.require(op in ('<', '>'), 'only strict coordinate comparisons supported')
        self.b.require(self.rational() == 0, 'comparison requires literal zero')
        return ['lt0' if op == '<' else 'gt0', coord]

    def conjunction(self):
        result = self.condition_atom()
        while self.peek() == '&&':
            self.take('&&')
            result = ['and', result, self.condition_atom()]
        return result

    def condition(self):
        result = self.conjunction()
        while self.peek() == '||':
            self.take('||')
            result = ['or', result, self.conjunction()]
        return result

    def assignment(self, allow_guard):
        if self.engine == 'PY':
            self.take('Eq'); self.take('(')
            coord = self.coordinate(); self.take(',')
            value = self.rational(); self.take(')')
            return coord, value, None
        coord = self.coordinate(); self.take('->')
        if self.peek() != 'ConditionalExpression':
            return coord, self.rational(), None
        self.b.require(allow_guard, 'a point cannot contain conditional coordinates')
        self.take('ConditionalExpression'); self.take('[')
        value = self.rational(); self.take(',')
        condition = self.condition(); self.take(']')
        return coord, value, condition

    def sequence(self, item):
        opening, closing = ('(', ')') if self.engine == 'PY' else ('{', '}')
        self.take(opening)
        items = []
        while self.peek() != closing:
            if items:
                self.take(',')
                if self.peek() == closing:
                    break
            items.append(item())
        self.take(closing)
        return items

    def locus(self):
        branches = self.sequence(lambda: self.sequence(lambda: self.assignment(True)))
        self.b.require(self.peek() is None and branches, 'empty or trailing locus')
        result = []
        for branch in branches:
            self.b.require(branch, 'empty locus branch')
            terms = []
            for coord, value, guard in branch:
                self.b.require(value == 0, 'locus requires zero coordinate assignments')
                terms.append(['eq0', coord])
                if guard is not None:
                    terms.append(guard)
            result.append(['and', *terms])
        return ['or', *result]

    def point(self):
        assignments = self.sequence(lambda: self.assignment(False))
        self.b.require(self.peek() is None and len(assignments) == 3 and
                       {c for c, _, _ in assignments} == {0, 1, 2}, 'point needs each coordinate exactly once')
        return [next(value for c, value, _ in assignments if c == i) for i in range(3)]


def proposition(tree):
    op, *args = tree
    if op == 'eq0':
        return f'k {args[0]} = 0'
    if op == 'lt0':
        return f'k {args[0]} < 0'
    if op == 'gt0':
        return f'0 < k {args[0]}'
    return '('+(' ∧ ' if op == 'and' else ' ∨ ').join(proposition(a) for a in args)+')'


def rational_lean(value):
    return f'({value.numerator} / {value.denominator} : ℝ)'


def generate(bridge, families):
    outputs, engines, audits = {}, {}, []
    for engine, path in bridge.INPUTS.items():
        rows = bridge.read_rows(path, engine)
        records, points = [], []
        ns = engine+'Loci'
        lines = ['import S10Audit.CAS.MinorLoci', 'import S10Audit.CAS.MinorBindings', '', 'namespace S10Audit.CAS.'+ns,
                 'open S10Pilot', 'noncomputable section', '']
        for label, root, stacked, *_ in families:
            suffix = f'ROOT{root+1}_Q8_'+('TRANSVERSE_' if stacked else '')+'RANK_DROP_LOCUS'
            tag = bridge.PREFIX+suffix
            bridge.require(tag in rows, 'missing selected locus: '+tag)
            line, raw = rows[tag]
            tree = CoordinateParser(bridge, raw, engine).locus()
            expected = 'k = 0' if root == 0 else ('k 0 = 0 ∨ (k 1 = 0 ∧ k 2 = 0)'
                                                if root == 2 and stacked else 'k 1 = 0 ∧ k 2 = 0')
            lines += [f'def {label} (k : Vec 3) : Prop := {proposition(tree)}', '',
                      f'theorem {label}_geometry (k : Vec 3) : {label} k ↔ ({expected}) := by',
                      '  classical', f'  locus_eval {label}', '',
                      f'theorem {label}_minors (rho mu sigma z : ℝ) (k : Vec 3)',
                      '    (hm : mu ≠ 0) (hs : sigma ≠ 0) (_hs1 : sigma ≠ 1) :',
                      f'    {label} k ↔ (∀ i, {engine}Minors.{label} rho mu sigma z k i = 0) := by',
                      f'  rw [{label}_geometry, {engine}Minors.{label}_zero_iff rho mu sigma z k hs,',
                      f'    CAS.{label}_locus mu '+('sigma ' if root else '')+'k hm'+
                      (' hs _hs1' if root == 2 else ' _hs1' if root else '')+']', '',
                      f'theorem {label}_all_minors (rho mu sigma : ℝ) (k : Vec 3)',
                      '    (hr : rho ≠ 0) (hm : mu ≠ 0) (hs : sigma ≠ 0) (hs1 : sigma ≠ 1) :',
                      f'    {label} k ↔ AllOrderedMinorsZero',
                      f'      ({"rootStack" if stacked else "rootMatrix"} rho mu sigma k {root}) '+
                      str(3 if stacked and root != 1 else 2)+' := by',
                      f'  rw [{label}_minors rho mu sigma 0 k hm hs hs1,',
                      f'    {engine}Minors.{label}_complete rho mu sigma 0 k hr hs]', '']
            records.append({'tag': engine+'_'+tag, 'line': line, 'payload': raw,
                            'payload_sha256': bridge.sha(raw.encode()), 'predicate_ast': tree,
                            'lean_locus': ns+'.'+label, 'geometry': expected})
            audits += [ns+'.'+label+'_'+s for s in ['geometry','minors','all_minors']]
        for index, geometry in [(1, 'parallel'), (2, 'perpendicular')]:
            suffix = f'Q8_STRATUM{index}_POINT' if engine == 'PY' else f'STRATUM{index}_Q8_POINT'
            tag = bridge.PREFIX+suffix
            bridge.require(tag in rows, 'missing targeted point: '+tag)
            line, raw = rows[tag]
            point = CoordinateParser(bridge, raw, engine).point()
            name = geometry+'Point'
            vector = '!['+', '.join(map(rational_lean, point))+']'
            lines += [f'def {name} : Vec 3 := {vector}',
                      f'theorem {name}_nonzero : {name} ≠ 0 := by',
                      '  intro h',
                      f'  have hv := congrFun h {0 if index == 1 else 1}',
                      f'  norm_num [{name}] at hv',
                      f'theorem {name}_target : extraTransverse {name} := by',
                      f'  rw [extraTransverse_geometry]; norm_num [{name}, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]',
                      f'theorem {name}_separate : '+
                      (f'ordinaryRank {name} ∧ {name} 0 ≠ 0' if index == 1 else f'¬ ordinaryRank {name} ∧ {name} 0 = 0')+' := by',
                      f'  rw [ordinaryRank_geometry]; norm_num [{name}, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]', '']
            points.append({'tag': engine+'_'+tag, 'line': line, 'payload': raw,
                           'payload_sha256': bridge.sha(raw.encode()),
                           'point': [[x.numerator,x.denominator] for x in point], 'geometry': geometry})
            audits += [ns+'.'+name+'_'+s for s in ['nonzero','target','separate']]
        lines += ['end', 'end S10Audit.CAS.'+ns, '']
        outputs[bridge.GENERATED/(ns+'.lean')] = '\n'.join(lines)
        engines[engine] = {'path': str(path.relative_to(bridge.BASE)),
                           'sha256': bridge.sha(path.read_bytes()), 'loci': records, 'points': points}
    lines = ['import S10Audit.CAS.PYLoci', 'import S10Audit.CAS.WLLoci', '',
             'namespace S10Audit.CAS', 'open S10Pilot', '']
    for label, *_ in families:
        lines += [f'theorem {label}_locus_cross_engine (k : Vec 3) : PYLoci.{label} k ↔ WLLoci.{label} k := by',
                  f'  rw [PYLoci.{label}_geometry, WLLoci.{label}_geometry]', '']
        audits.append(label+'_locus_cross_engine')
    lines += ['end S10Audit.CAS', '']
    outputs[bridge.GENERATED/'LocusBindings.lean'] = '\n'.join(lines)
    audits += [label+'_locus' for label,*_ in families]
    audits += ['wavevector_zero_iff','normSq_zero_iff','perp_pair_zero_iff','positive_or_negative']
    audits += ['exceptional_strata_disjoint','extraTransverse_nonzero_partition',
               'perpendicular_drops_only_transverse','parallel_drops_both','static_no_nonzero_rank_drop']
    return outputs, engines, audits
