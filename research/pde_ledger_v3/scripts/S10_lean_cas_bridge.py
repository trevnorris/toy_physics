#!/usr/bin/env python3
"""Preserve selected CAS arithmetic as Expr trees and emit Lean obligations.

The parser accepts only exact integer/rational arithmetic, natural powers,
explicitly mapped symbols, containers, and SymPy's Matrix container spelling.
It never evaluates a transcript as Python or asks a CAS to simplify its trees.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
import sys

BASE = Path(__file__).resolve().parents[1]
LEAN = BASE / 'lean'
INPUTS = {'PY': BASE / 'scripts/out/S10_anisotropic_strata_sympy_audit.out',
          'WL': BASE / 'mathematica/out/S10_anisotropic_strata_mathematica_audit.out'}
GENERATED = LEAN / 's10/S10Audit/CAS'
PREFIX = 'S10_XFORM_ANISO_D3_'
NAMES = {
    'PY': {'rho_br': 'rho', 'mu_R': 'mu', 's_rho': 'sigma', 'omegaSquared': 'z'},
    'WL': {'rhoBr': 'rho', 'muR': 'mu', 'sRho': 'sigma', 'omegaSquared': 'z'},
}
for mapping in NAMES.values():
    mapping.update({f'k{i+1}': f'k{i}' for i in range(3)})
DIMS = {'rho': (-3, 0, 1), 'mu': (-1, -2, 1), 'sigma': (0, 0, 0),
        'z': (0, -2, 0), **{f'k{i}': (-1, 0, 0) for i in range(3)}}


class BridgeError(ValueError):
    pass


@dataclass(frozen=True)
class Node:
    op: str
    args: tuple


def require(test, message):
    if not test:
        raise BridgeError(message)


class Parser:
    token = re.compile(r'\s+|\*\*|[A-Za-z_][A-Za-z_0-9]*|[0-9]+|[+*/^(),{}\[\]-]')

    def __init__(self, raw: str, engine: str):
        require(engine in NAMES, 'unsupported engine')
        require(len(raw) <= 100000, 'payload exceeds parser bound')
        self.engine, self.tokens, self.pos = engine, [], 0
        end = 0
        for match in self.token.finditer(raw):
            require(match.start() == end, f'unsupported syntax at byte {end}')
            end = match.end()
            if not match.group().isspace():
                self.tokens.append(match.group())
        require(end == len(raw), f'unsupported trailing syntax at byte {end}')

    def peek(self):
        return self.tokens[self.pos] if self.pos < len(self.tokens) else None

    def take(self, expected=None):
        value = self.peek()
        require(value is not None and (expected is None or value == expected),
                f'expected {expected or "token"}, got {value}')
        self.pos += 1
        return value

    def parse(self):
        result = self.expression()
        require(self.peek() is None, f'unconsumed token {self.peek()}')
        return result

    def expression(self):
        result = self.product()
        while self.peek() in ('+', '-'):
            op = self.take()
            other = self.product()
            require(isinstance(result, Node) and isinstance(other, Node), 'container arithmetic')
            result = Node('add' if op == '+' else 'sub', (result, other))
        return result

    def product(self):
        result = self.unary()
        while self.peek() in ('*', '/'):
            op = self.take()
            other = self.unary()
            require(isinstance(result, Node) and isinstance(other, Node), 'container arithmetic')
            result = Node('mul' if op == '*' else 'div', (result, other))
        return result

    def unary(self):
        if self.peek() in ('+', '-'):
            op = self.take()
            term = self.unary()
            require(isinstance(term, Node), 'container unary operation')
            return term if op == '+' else Node('mul', (Node('num', (-1,)), term))
        return self.power()

    def power(self):
        result = self.atom()
        if self.peek() in ('**', '^'):
            op = self.take()
            require(op == ('**' if self.engine == 'PY' else '^'), 'wrong engine power syntax')
            exponent = self.unary()
            require(isinstance(exponent, Node) and exponent.op == 'num' and
                    0 <= exponent.args[0] <= 32, 'only bounded natural powers are supported')
            require(isinstance(result, Node), 'container power')
            result = Node('pow', (result, exponent.args[0]))
        return result

    def container(self, opening, closing, parens=False):
        self.take(opening)
        if self.peek() == closing:
            self.take(closing)
            return []
        first = self.expression()
        if parens and self.peek() == closing:
            self.take(closing)
            return first
        result = [first]
        while self.peek() != closing:
            self.take(',')
            if self.peek() == closing:
                break
            result.append(self.expression())
        self.take(closing)
        return result

    def atom(self):
        token = self.peek()
        if token == '(':
            return self.container('(', ')', parens=True)
        if token == '[':
            require(self.engine == 'PY', 'Wolfram function syntax is unsupported')
            return self.container('[', ']')
        if token == '{':
            require(self.engine == 'WL', 'Python sets/dicts are unsupported')
            return self.container('{', '}')
        if token == 'Matrix':
            require(self.engine == 'PY', 'unexpected Matrix function')
            self.take('Matrix')
            self.take('(')
            result = self.container('[', ']')
            self.take(')')
            require(result and all(isinstance(row, list) for row in result), 'matrix row shape')
            require(len({len(row) for row in result}) == 1, 'ragged matrix')
            return result
        token = self.take()
        if token.isdigit():
            require(len(token) <= 30, 'integer literal exceeds parser bound')
            return Node('num', (int(token),))
        require(token in NAMES[self.engine], f'unknown symbol {token}')
        return Node('atom', (NAMES[self.engine][token],))


def dimension(node, dims=DIMS):
    if node.op == 'num':
        return (0, 0, 0)
    if node.op == 'atom':
        return dims[node.args[0]]
    if node.op == 'pow':
        return tuple(x * node.args[1] for x in dimension(node.args[0], dims))
    a, b = (dimension(arg, dims) for arg in node.args)
    if node.op in ('add', 'sub'):
        require(a == b, f'inhomogeneous additive branches: {a} versus {b}')
        return a
    return tuple(x + y if node.op == 'mul' else x - y for x, y in zip(a, b))


def numeric(node):
    if node.op == 'num':
        return Fraction(node.args[0])
    if node.op == 'atom':
        return None
    if node.op == 'pow':
        a = numeric(node.args[0])
        return None if a is None else a ** node.args[1]
    a, b = map(numeric, node.args)
    if node.op == 'div':
        require(b != 0, 'literal zero denominator')
    if a is None or b is None:
        return None
    return {'add': lambda: a+b, 'sub': lambda: a-b, 'mul': lambda: a*b,
            'div': lambda: a/b}[node.op]()


def denominators(node):
    """Return factors whose nonvanishing suffices for every raw division."""
    def factors(term):
        value = numeric(term)
        if value is not None:
            require(value != 0, 'zero denominator')
            return []
        if term.op in ('mul', 'div'):
            return factors(term.args[0]) + factors(term.args[1])
        if term.op == 'pow' and term.args[1] > 0:
            return factors(term.args[0])
        return [term]
    found = []
    if node.op == 'div':
        found.extend(factors(node.args[1]))
    for arg in node.args:
        if isinstance(arg, Node):
            found.extend(denominators(arg))
    return list(dict.fromkeys(found))


def node_json(node):
    return [node.op, *[node_json(a) if isinstance(a, Node) else a for a in node.args]]


def lean_symbol(name):
    return f'(.k {name[1:]})' if name.startswith('k') else f'.{name}'


def lean_number(n):
    return f'({n})' if n < 0 else str(n)


def arithmetic(node):
    if node.op == 'num':
        return lean_number(node.args[0])
    if node.op == 'atom':
        a = node.args[0]
        return f'(k {a[1:]})' if a.startswith('k') else a
    if node.op == 'pow':
        return f'({arithmetic(node.args[0])} ^ {node.args[1]})'
    op = {'add': '+', 'sub': '-', 'mul': '*', 'div': '/'}[node.op]
    return f'({arithmetic(node.args[0])} {op} {arithmetic(node.args[1])})'


def zero_lift(dim):
    """Typed-zero adapter: only a literal zero may receive its declared slot unit."""
    l, t, m = dim
    require(t % 2 == 0, 'zero slot requires unsupported half power')
    terms = [('rho', m), ('z', -t//2), ('k0', -l-3*m)]
    unit = Node('num', (1,))
    for name, power in terms:
        if power:
            factor = Node('pow', (Node('atom', (name,)), abs(power)))
            unit = Node('mul' if power > 0 else 'div', (unit, factor))
    return Node('mul', (Node('num', (0,)), unit))


def read_rows(path, engine):
    rows = {}
    for number, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip():
            continue
        match = re.fullmatch(r'(PY|WL)_([A-Z][A-Z0-9_]*): (.*)', line)
        if not match:
            require(engine == 'WL' and line ==
                    'Solve::svars: Equations may not give solutions for all "solve" variables.',
                    f'{path}:{number}: unrecognized transcript line')
            continue
        require(match[1] == engine, f'{path}:{number}: wrong engine tag')
        require(match[2] not in rows, f'{path}:{number}: duplicate tag {match[2]}')
        rows[match[2]] = (number, match[3])
    return rows


def sha(data):
    return hashlib.sha256(data).hexdigest()


def dim_lean(dim):
    return 'dimensions ' + ' '.join(lean_number(x) for x in dim)


class Builder:
    def __init__(self, engine, path, *, namespace=None, support='Support', tactic='cas_eval',
                 units='units', dims=DIMS, zero_adapter=zero_lift):
        self.engine, self.path = engine, path
        self.namespace = namespace or engine
        self.support, self.tactic = support, tactic
        self.units, self.dims, self.zero_adapter = units, dims, zero_adapter
        self.rows = read_rows(path, engine)
        self.nodes, self.definitions, self.claims, self.records, self.audits = {}, [], [], [], []

    def intern(self, node):
        if node in self.nodes:
            return self.nodes[node]
        children = [self.intern(a) for a in node.args if isinstance(a, Node)]
        name = f'n{len(self.nodes)}'
        self.nodes[node] = name
        dim = dimension(node, self.dims)
        if node.op == 'num':
            expr = f'.scalar {lean_number(node.args[0])}'
            proof = f'  exact castDim (Expr.HasDim.scalar {lean_number(node.args[0])}) dimensions_zero.symm'
        elif node.op == 'atom':
            expr = f'.atom {lean_symbol(node.args[0])}'
            proof = f'  exact Expr.HasDim.atom (u := {self.units}) {lean_symbol(node.args[0])}'
        elif node.op == 'pow':
            expr = f'.pow {children[0]} {node.args[1]}'
            proof = (f'  apply castDim (Expr.HasDim.pow {node.args[1]} {children[0]}_dim)\n'
                     '  norm_num [dimensions_pow, dimensions_zero]')
        else:
            expr = f'.{node.op} {children[0]} {children[1]}'
            term = f'Expr.HasDim.{node.op} {children[0]}_dim {children[1]}_dim'
            proof = ('  exact ' + term if node.op in ('add', 'sub') else
                     f'  apply castDim ({term})\n  norm_num [dimensions_{node.op}, dimensions_zero]')
        evaluation = '  rfl' if not children else (
            '  simp only [' + ', '.join([name, 'Expr.eval', *dict.fromkeys(x+'_eval' for x in children)]) + ']')
        self.definitions.append(f'''def {name} : Expr Symbol := {expr}
theorem {name}_dim : Expr.HasDim {self.units} {name} ({dim_lean(dim)}) := by
{proof}
theorem {name}_eval (rho mu sigma z : ℝ) (k : Vec 3) :
    Expr.eval (values rho mu sigma z k) {name} = {arithmetic(node)} := by
{evaluation}
''')
        return name

    def record(self, suffix, shape, cells, *, parse_payload=None):
        """cells supplies independent reference operands and declared slot units."""
        tag = PREFIX + suffix
        require(tag in self.rows, f'{self.engine}: missing {tag}')
        line, raw = self.rows[tag]
        parsed = Parser(raw, self.engine).parse() if parse_payload is None else parse_payload(raw, self.engine)
        payload = shape(parsed)
        require(len(payload) == len(cells), f'{tag}: wrong cell count')
        names, dims, equations, domains, entries = [], [], [], [], []
        for index, (node, (expected, dim, extra_domain)) in enumerate(zip(payload, cells)):
            require(isinstance(node, Node), f'{tag}: non-scalar cell {index}')
            numeric(node)  # reject zero numeric denominators even in nested arithmetic
            raw_domain = denominators(node)
            domain = list(dict.fromkeys(raw_domain + [Node('atom', (x,)) for x in extra_domain]))
            lifted = node.op == 'num' and node.args[0] == 0 and dim != (0, 0, 0)
            tree = self.zero_adapter(dim) if lifted else node
            require(dimension(tree, self.dims) == dim, f'{tag}[{index}]: {dimension(tree, self.dims)} != {dim}')
            name = self.intern(tree)
            claim = suffix.lower() + f'_cell{index}'
            hypotheses = ' '.join(f'(_h{i} : {arithmetic(d)} ≠ 0)' for i, d in enumerate(domain))
            equation = f'Expr.eval (values rho mu sigma z k) {name} = {expected}'
            self.claims.append(f'''theorem {claim} (rho mu sigma z : ℝ) (k : Vec 3) {hypotheses} :
    {equation} := by
  {self.tactic} {name}_eval
''')
            names.append(claim)
            dims.append((name, dim))
            equations.append(equation)
            domains.append(domain)
            entries.append({'index': index, 'ast': node_json(node), 'dimension': list(dim),
                            'literal_zero_lift': lifted,
                            'raw_denominator_factors': [node_json(d) for d in raw_domain],
                            'reference_domain_symbols': extra_domain,
                            'lean_tree': name, 'lean_claim': claim, 'reference': expected})
        union_domain = list(dict.fromkeys(d for domain in domains for d in domain))
        hd = {d: f'h{i}' for i, d in enumerate(union_domain)}
        hypotheses = ' '.join(f'({hd[d]} : {arithmetic(d)} ≠ 0)' for d in union_domain)
        aggregate = suffix.lower()
        proofs = [f'{name} rho mu sigma z k' + ''.join(' '+hd[d] for d in domain)
                  for name, domain in zip(names, domains)]
        # A singleton has no conjunction constructor.
        value_proof = proofs[0] if len(proofs) == 1 else '⟨'+', '.join(proofs)+'⟩'
        dim_proof = dims[0][0]+'_dim' if len(dims) == 1 else '⟨'+', '.join(n+'_dim' for n, _ in dims)+'⟩'
        value_prop = ' ∧\n    '.join(equations)
        dim_prop = ' ∧\n    '.join(f'Expr.HasDim {self.units} {name} ({dim_lean(dim)})' for name, dim in dims)
        self.claims.append(f'''theorem {aggregate}_values (rho mu sigma z : ℝ) (k : Vec 3) {hypotheses} :
    {value_prop} :=
  {value_proof}

theorem {aggregate}_dimensions :
    {dim_prop} :=
  {dim_proof}
''')
        self.audits += [aggregate+'_values', aggregate+'_dimensions']
        self.records.append({'tag': self.engine+'_'+tag, 'line': line, 'payload': raw,
                             'payload_sha256': sha(raw.encode()), 'cells': entries,
                             'lean_aggregate': aggregate})

    def output(self):
        return (f'''-- Generated by scripts/S10_lean_cas_bridge.py; use --check to verify provenance.
-- Input SHA-256: {sha(self.path.read_bytes())}
import S10Audit.CAS.{self.support}

set_option backward.isDefEq.respectTransparency false

namespace S10Audit.CAS.{self.namespace}
open S10Pilot S10Anisotropic
noncomputable section

''' + '\n'.join(self.definitions+self.claims) + f'\nend\nend S10Audit.CAS.{self.namespace}\n')


def scalar_shape(value):
    require(isinstance(value, Node), 'expected scalar')
    return [value]


def matrix_shape(value, rows, cols):
    require(isinstance(value, list) and len(value) == rows, f'expected {rows} rows')
    require(all(isinstance(row, list) and len(row) == cols for row in value), f'expected {cols} columns')
    return [x for row in value for x in row]


def vector_shape(value, engine):
    if engine == 'PY':
        return matrix_shape(value, 3, 1)
    require(isinstance(value, list) and len(value) == 3, 'expected three-vector')
    return value


def build(engine, path):
    b = Builder(engine, path)
    md, kd, zd, bd = (-3, -2, 1), (-1, 0, 0), (0, -2, 0), (0, 0, 0)
    for suffix, factor in [('Q2_MATRIX_A', '(-1)'), ('Q2_MATRIX_B', '(1/2)'),
                           ('Q2_MATRIX_RESIDUAL', '(-3/2)')]:
        b.record(suffix, lambda x: matrix_shape(x, 3, 3),
                 [(f'{factor} * referenceMatrix rho mu sigma z k {i} {j}', md, [])
                  for i in range(3) for j in range(3)])
    b.record('Q2_MATRIX_ENTRY_RATIO', scalar_shape, [('(-2)', bd, [])])
    b.record('Q3_DETERMINANT', scalar_shape,
             [('Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma z k)', (-9, -6, 3), [])])
    b.record('ROOT_ORDERING', lambda x: x,
             [(f'referenceRoot rho mu sigma k {r}', zd, [] if r == 0 else ['rho'] + (['sigma'] if r == 2 else []))
              for r in range(3)])
    for r in range(3):
        root = f'(referenceRoot rho mu sigma k {r})'
        matrix = f'((1/2 : ℝ) • referenceMatrix rho mu sigma {root} k)'
        basis = f'(referenceBasis sigma k {r})'
        root_domain = [] if r == 0 else ['rho'] + (['sigma'] if r == 2 else [])
        basis_domain = [['k2'], ['k1'], ['sigma', 'k0', 'k2']][r]
        rp = f'ROOT{r+1}_'
        b.record(rp+'N1_MATRIX', lambda x: matrix_shape(x, 3, 3),
                 [(f'{matrix} {i} {j}', md, root_domain) for i in range(3) for j in range(3)])
        b.record(rp+'N3_STACKED_MATRIX', lambda x: matrix_shape(x, 4, 3),
                 [(f'{matrix} {i} {j}' if i < 3 else f'k {j}', md if i < 3 else kd,
                   root_domain if i < 3 else []) for i in range(4) for j in range(3)])
        b.record(rp+('N5_MATRIX_TIMES_K' if engine == 'PY' else 'N5_WAVEVECTOR_PRODUCT'), lambda x: vector_shape(x, engine),
                 [(f'{matrix}.mulVec k {i}', (-4, -2, 1), root_domain) for i in range(3)])
        def basis_shape(x):
            require(isinstance(x, list) and len(x) == 1, 'expected one basis vector in this generic D3 chart')
            return vector_shape(x[0], engine)
        b.record(rp+'N6_NULLSPACE_BASIS', basis_shape,
                 [(f'{basis} {i}', bd, basis_domain) for i in range(3)])
        dot_suffix = 'N6_BASIS_DOT_K' if engine == 'PY' else 'N6_BASIS_DOTS'
        b.record(rp+dot_suffix, lambda x: x, [(f'dot k {basis}', kd, basis_domain)])
        res_suffix = 'N6_BASIS_VECTOR_RESIDUALS' if engine == 'PY' else 'N6_BASIS_RESIDUALS'
        b.record(rp+res_suffix, basis_shape,
                 [(f'normSq k * {basis} {i} - dot k {basis} * k {i}', (-2, 0, 0), basis_domain)
                  for i in range(3)])
    return b


def matrix_bindings(builders):
    text = ['import S10Audit.CAS.PY', 'import S10Audit.CAS.WL',
            'import S10Audit.CAS.BasisCompletion', '',
            'set_option backward.isDefEq.respectTransparency false', '',
            'namespace S10Audit.CAS', 'open S10Pilot', 'noncomputable section', '']
    audits = []
    for engine, builder in builders.items():
        text.append(f'namespace {engine}\n')
        for label, suffix, factor in [('matrixA', 'Q2_MATRIX_A', '(-1)'),
                                      ('matrixB', 'Q2_MATRIX_B', '(1/2)'),
                                      ('matrixResidual', 'Q2_MATRIX_RESIDUAL', '(-3/2)')]:
            record = next(r for r in builder.records if r['tag'].endswith('_'+suffix))
            names = [c['lean_tree'] for c in record['cells']]
            matrix = '!['+', '.join('!['+', '.join(names[3*i:3*i+3])+']' for i in range(3))+']'
            text.append(f'''def {label} (rho mu sigma z : ℝ) (k : Vec 3) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => Expr.eval (values rho mu sigma z k)
    (({matrix} : Fin 3 → Fin 3 → Expr Symbol) i j)

theorem {label}_reference (rho mu sigma z : ℝ) (k : Vec 3) :
    {label} rho mu sigma z k = ({factor} : ℝ) • referenceMatrix rho mu sigma z k := by
  ext i j
  fin_cases i <;> fin_cases j
'''+''.join(f'  · exact {suffix.lower()}_cell{i} rho mu sigma z k\n' for i in range(9)))
            audits.append(f'{engine}.{label}_reference')
        text.append(f'''theorem routes (rho mu sigma z : ℝ) (k : Vec 3) :
    matrixA rho mu sigma z k = (-2 : ℝ) • matrixB rho mu sigma z k := by
  rw [matrixA_reference, matrixB_reference, smul_smul]
  norm_num

theorem matrixB_action (rho mu sigma omega : ℝ) (k : Vec 3) :
    matrixB rho mu sigma (omega ^ 2) k =
      (1/2 : ℝ) • actionMatrix .anisotropic 0 rho mu sigma 1 omega k := by
  rw [matrixB_reference, referenceMatrix_action]

theorem matrixB_kernel (rho mu sigma z : ℝ) (k a : Vec 3) :
    (matrixB rho mu sigma z k).mulVec a = 0 ↔
      (referenceMatrix rho mu sigma z k).mulVec a = 0 := by
  rw [matrixB_reference, Matrix.smul_mulVec, smul_eq_zero]
  norm_num

end {engine}
''')
        audits += [f'{engine}.{name}' for name in ['routes', 'matrixB_action', 'matrixB_kernel']]
        text.append(f'namespace {engine}\n')
        basis_records = [next(record for record in builder.records
                              if record['tag'].endswith(f'_ROOT{r+1}_N6_NULLSPACE_BASIS'))
                         for r in range(3)]
        basis_trees = '!['+', '.join('!['+', '.join(c['lean_tree'] for c in record['cells'])+']'
                                    for record in basis_records)+']'
        text.append(f'''def basis (rho mu sigma z : ℝ) (k : Vec 3) : Fin 3 → Vec 3 :=
  fun r j => Expr.eval (values rho mu sigma z k)
    (({basis_trees} : Fin 3 → Fin 3 → Expr Symbol) r j)

theorem basis_reference (rho mu sigma z : ℝ) (k : Vec 3) (r : Fin 3)
    (h : GenericChart sigma k) : basis rho mu sigma z k r = referenceBasis sigma k r := by
  rcases h with ⟨hs, _hs1, h0, h1, h2⟩
  have hsn : sigma ≠ 0 := ne_of_gt hs
  ext j
  fin_cases r <;> fin_cases j
''')
        hypothesis_names = {'sigma': 'hsn', 'k0': 'h0', 'k1': 'h1', 'k2': 'h2'}
        for record in basis_records:
            for cell in record['cells']:
                assert all(d[0] == 'atom' for d in cell['raw_denominator_factors'])
                domain = list(dict.fromkeys([d[1] for d in cell['raw_denominator_factors']]
                                             + cell['reference_domain_symbols']))
                args = ''.join(' '+hypothesis_names[d] for d in domain)
                text.append(f'  · exact {cell["lean_claim"]} rho mu sigma z k{args}')
        text.append('''
theorem basis_independent (rho mu sigma z : ℝ) (k : Vec 3) (r : Fin 3)
    (h : GenericChart sigma k) :
    LinearIndependent ℝ (fun _ : Fin 1 => basis rho mu sigma z k r) := by
  rw [basis_reference rho mu sigma z k r h]
  exact referenceBasis_independent sigma k r (genericChart_basisDomain sigma k r h)

theorem basis_complete (rho mu sigma : ℝ) (k a : Vec 3) (r : Fin 3)
    (hr : rho ≠ 0) (hm : mu ≠ 0) (h : GenericChart sigma k) :
    (matrixB rho mu sigma (referenceRoot rho mu sigma k r) k).mulVec a = 0 ↔
      a ∈ Submodule.span ℝ {basis rho mu sigma (referenceRoot rho mu sigma k r) k r} := by
  rw [matrixB_kernel, basis_reference rho mu sigma _ k r h]
  exact referenceMatrix_basis_complete rho mu sigma k a r hr hm h
''')
        text.append(f'end {engine}\n')
        audits += [f'{engine}.{name}' for name in ['basis_reference', 'basis_independent', 'basis_complete']]
    for label in ['matrixA', 'matrixB', 'matrixResidual']:
        text.append(f'''theorem {label}_cross_engine (rho mu sigma z : ℝ) (k : Vec 3) :
    PY.{label} rho mu sigma z k = WL.{label} rho mu sigma z k := by
  rw [PY.{label}_reference, WL.{label}_reference]
''')
        audits.append(label+'_cross_engine')
    text.append('''theorem basis_cross_engine (rho mu sigma z : ℝ) (k : Vec 3) (r : Fin 3)
    (h : GenericChart sigma k) : PY.basis rho mu sigma z k r = WL.basis rho mu sigma z k r := by
  rw [PY.basis_reference rho mu sigma z k r h, WL.basis_reference rho mu sigma z k r h]
''')
    audits.append('basis_cross_engine')
    text += ['end', 'end S10Audit.CAS', '']
    return '\n'.join(text), audits


def generate(check=False):
    import S10_lean_cas_minors as minor_bridge
    import S10_lean_cas_loci as locus_bridge
    import S10_lean_cas_reruns as rerun_bridge
    import S10_lean_cas_counts as count_bridge
    import S10_lean_cas_generic_counts as generic_count_bridge
    import S10_lean_cas_roots as root_bridge
    import S10_lean_cas_coincidence as coincidence_bridge
    import S10_lean_cas_records as record_bridge
    outputs, engines, builders = {}, {}, {}
    audits = []
    for engine, path in INPUTS.items():
        builder = build(engine, path)
        builders[engine] = builder
        target = GENERATED / f'{engine}.lean'
        outputs[target] = builder.output()
        engines[engine] = {'path': str(path.relative_to(BASE)), 'sha256': sha(path.read_bytes()),
                           'records': builder.records, 'unique_trees': len(builder.nodes)}
        audits += [f'#print axioms S10Audit.CAS.{engine}.{name}' for name in builder.audits]
    bindings, binding_audits = matrix_bindings(builders)
    outputs[GENERATED/'Bindings.lean'] = bindings
    minor_outputs, minor_engines, minor_audits = minor_bridge.generate(sys.modules[__name__])
    outputs.update(minor_outputs)
    locus_outputs, locus_engines, locus_audits = locus_bridge.generate(sys.modules[__name__], minor_bridge.FAMILIES)
    outputs.update(locus_outputs)
    rerun_outputs, rerun_engines, rerun_audits = rerun_bridge.generate(sys.modules[__name__])
    outputs.update(rerun_outputs)
    count_outputs, count_engines, count_audits = count_bridge.generate(sys.modules[__name__])
    outputs.update(count_outputs)
    generic_count_outputs, generic_count_engines, generic_count_audits = generic_count_bridge.generate(sys.modules[__name__], builders)
    outputs.update(generic_count_outputs)
    root_outputs, root_engines, root_controls, root_audits = root_bridge.generate(sys.modules[__name__])
    outputs.update(root_outputs)
    coincidence_outputs, coincidence_engines, coincidence_records, coincidence_audits = coincidence_bridge.generate(sys.modules[__name__])
    outputs.update(coincidence_outputs)
    record_outputs, record_engines, record_metadata, record_audits = record_bridge.generate(sys.modules[__name__])
    outputs.update(record_outputs)
    audit_imports = ['import S10Audit.CAS.RecordBindings','import S10Audit.CAS.StatusBindings','import S10Audit.CAS.CoincidenceBindings','import S10Audit.CAS.Bindings', 'import S10Audit.CAS.LocusBindings',
                     'import S10Audit.CAS.RerunBindings', 'import S10Audit.CAS.CountBindings',
                     'import S10Audit.CAS.GenericCountBindings', 'import S10Audit.CAS.RootBindings']
    audits += [f'#print axioms S10Audit.CAS.{name}' for name in record_audits]
    audits += [f'#print axioms S10Audit.CAS.{name}' for name in coincidence_audits]
    audits += [f'#print axioms S10Audit.CAS.{name}' for name in root_audits]
    audits += [f'#print axioms S10Audit.CAS.{name}' for name in generic_count_audits]
    audits += [f'#print axioms S10Audit.CAS.{name}' for name in count_audits]
    audits += [f'#print axioms S10Audit.CAS.{name}' for name in rerun_audits]
    audits += [f'#print axioms S10Audit.CAS.{name}' for name in minor_audits]
    audits += [f'#print axioms S10Audit.CAS.{name}' for name in locus_audits]
    audits += [f'#print axioms S10Audit.CAS.{name}' for name in binding_audits]
    audits += [f'#print axioms S10Audit.CAS.{name}' for name in [
        'rho_units', 'mu_units', 'referenceMatrix_action', 'referenceMatrix_normalized',
        'referenceBasis_last', 'referenceBasis_independent', 'staticBasis_span',
        'ordinaryBasis_span', 'extraBasis_span', 'genericChart_basisDomain',
        'referenceRoot_normalized', 'referenceMatrix_basis_complete']]
    outputs[GENERATED/'Audit.lean'] = '\n'.join(audit_imports)+'\n\n'+'\n'.join(audits)+'\n'
    manifest = {'scope': 'XFORM_ANISO D3 generic expressions, complete rank-drop minors, coordinate loci, targeted points, exceptional reruns, generic/exceptional rank/count records, root lists, multiplicities, syntactic filter counts coincidence predicates/decisions/witnesses, aggregate/Q8 records, root signs, spectrum operands/statuses and skipped-stratum decisions from both focused transcripts',
                'parser': 'strict arithmetic and selected logical grammars; no CAS evaluation or simplification',
                'squared_frequency': 'generic/minor expressions: independent real z = omega^2; fixed-point reruns: reduced z with omega^2 = κ^2 z, as declared in each rerun engine record',
                'zero_policy': 'literal zero only may be lifted to its explicitly declared slot unit',
                'generator_sha256': sha(Path(__file__).read_bytes()),
                'minor_generator_sha256': sha(Path(minor_bridge.__file__).read_bytes()),
                'locus_generator_sha256': sha(Path(locus_bridge.__file__).read_bytes()),
                'rerun_generator_sha256': sha(Path(rerun_bridge.__file__).read_bytes()),
                'count_generator_sha256': sha(Path(count_bridge.__file__).read_bytes()),
                'generic_count_generator_sha256': sha(Path(generic_count_bridge.__file__).read_bytes()),
                'record_generator_sha256': sha(Path(record_bridge.__file__).read_bytes()),
                'coincidence_generator_sha256': sha(Path(coincidence_bridge.__file__).read_bytes()),
                'root_generator_sha256': sha(Path(root_bridge.__file__).read_bytes()),
                'generated_sha256': {str(p.relative_to(BASE)): sha(s.encode()) for p, s in outputs.items()},
                'engines': engines, 'minor_engines': minor_engines, 'locus_engines': locus_engines,
                'rerun_engines': rerun_engines, 'count_engines': count_engines,
                'generic_count_engines': generic_count_engines, 'root_engines': root_engines,
                'root_filter_records': root_controls, 'coincidence_engines': coincidence_engines,
                'coincidence_records': coincidence_records, 'record_engines': record_engines, 'record_metadata': record_metadata}
    outputs[BASE/'_measurements/S10_lean_cas_bridge_manifest.json'] = json.dumps(manifest, indent=2)+'\n'
    for path, text in outputs.items():
        if check:
            require(path.exists() and path.read_text() == text, f'stale generated artifact: {path}')
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)
    arithmetic_records = [r for group in (engines, minor_engines, rerun_engines, count_engines, generic_count_engines, root_engines, coincidence_engines, record_engines) for data in group.values() for r in data['records']]
    print(('CHECKED' if check else 'GENERATED') +
          f': {len(arithmetic_records)} arithmetic records/{sum(len(r["cells"]) for r in arithmetic_records)} cells; '+
          f'{sum(len(d["loci"]) for d in locus_engines.values())} loci/'+
          f'{sum(len(d["points"]) for d in locus_engines.values())} targeted points')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Reject stale generated files without writing.')
    args = parser.parse_args()
    generate(args.check)


if __name__ == '__main__':
    main()
