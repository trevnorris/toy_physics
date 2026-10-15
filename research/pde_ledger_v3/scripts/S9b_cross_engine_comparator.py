#!/usr/bin/env python3
"""S9b cross-engine instrument: lossless input, explicit bindings, no target.

Readers adapted from O2_cross_engine_comparator.py (constructor AST/InputForm
Pratt parser). No engine code is imported or executed. Unsupported mathematics
is coverage, malformed input is operational failure. Outputs are JSON lines.
"""
from __future__ import annotations
import argparse
import faulthandler
import signal
import ast
from collections import Counter, defaultdict, OrderedDict
from dataclasses import dataclass, asdict
from functools import lru_cache
from itertools import product
import hashlib
import io
import json
from pathlib import Path
import re
import resource
import sys
import time
import mpmath as mp
import sympy as sp
from sympy.core.function import AppliedUndef
from sympy.core.relational import Relational
from sympy.logic.boolalg import Boolean
from sympy.logic.boolalg import BooleanFunction

HERE = Path(__file__).resolve().parent
PY_SOURCE = 'research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py'
WL_SOURCE = 'research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl'
SPEC = 'research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md'
DEFAULT_PY = HERE / 'out/S9b_light_bending_sympy_audit.out'
DEFAULT_WL = HERE.parent / 'mathematica/out/S9b_light_bending_mathematica_audit.out'

class InputError(ValueError):
    """Missing, ambiguous, or ungrammatical input/configuration."""


@dataclass(frozen=True, slots=True)
class Node:
    head: str
    args: tuple = ()
    value: str = ''
    options: tuple = ()


def atom(head, value):
    return Node(head, value=str(value))


def seq(*items):
    return Node('Sequence', tuple(items))


def record(items):
    keys = [key for key, _ in items]
    if len(keys) != len(set(keys)):
        raise InputError('duplicate record key')
    return Node('Record', tuple(Node('Entry', (value,), key) for key, value in items))


def py_parse(source):
    """Parse lossless constructor grammar without eval or engine execution."""
    def walk(n):
        if isinstance(n, ast.Constant):
            return atom('Boolean' if isinstance(n.value, bool) else
                        'Text' if isinstance(n.value, str) else 'Number', n.value)
        if isinstance(n, ast.Name):
            return atom('Boolean' if n.id in ('true', 'false', 'True', 'False')
                        else 'Name', n.id)
        if isinstance(n, (ast.Tuple, ast.List)):
            return seq(*(walk(a) for a in n.elts))
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub):
            v = walk(n.operand)
            if v.head == 'Number':
                return atom('Number', '-' + v.value)
            return Node('Mul', (atom('Number', -1), v))
        if not isinstance(n, ast.Call):
            raise InputError('outside lossless constructor grammar: ' + type(n).__name__)
        options = []
        for k in n.keywords:
            if k.arg is None:
                d = ast.literal_eval(k.value)
                if not isinstance(d, dict):
                    raise InputError('constructor options are not a dictionary')
                options.extend((str(key), repr(value)) for key, value in d.items())
            else:
                options.append((k.arg, repr(ast.literal_eval(k.value))))
        options = tuple(sorted(options))
        a = tuple(walk(x) for x in n.args)
        if isinstance(n.func, ast.Call):
            return Node('Apply', (walk(n.func), *a), options=options)
        if not isinstance(n.func, ast.Name):
            raise InputError('non-constructor Python call')
        h = n.func.id
        if h in ('Symbol', 'Dummy', 'Str', 'Function'):
            if len(a) != 1 or a[0].head != 'Text':
                raise InputError('invalid named constructor ' + h)
            return Node({'Symbol': 'Symbol', 'Dummy': 'Dummy', 'Str': 'Text',
                         'Function': 'FunctionName'}[h], value=a[0].value, options=options)
        if h == 'Integer':
            if len(a) != 1 or a[0].head != 'Number':
                raise InputError('invalid Integer')
            return a[0]
        if h == 'Tuple':
            if a and all(v.head == 'Sequence' and len(v.args) == 2
                         and v.args[0].head == 'Text' for v in a) and len({v.args[0].value for v in a}) == len(a):
                return record([(v.args[0].value, v.args[1]) for v in a])
            return Node('Sequence', a)
        if h in ('ImmutableDenseMatrix', 'MutableDenseMatrix'):
            if len(a) != 1 or a[0].head != 'Sequence':
                raise InputError('invalid matrix constructor')
            return Node('Matrix', a[0].args)
        return Node(h, a, options=options)
    try:
        return walk(ast.parse(source, mode='eval').body)
    except (SyntaxError, TypeError, ValueError) as exc:
        raise InputError('SymPy grammar: ' + str(exc)) from exc


# A small InputForm parser: arithmetic, associations, lists and arbitrary held
# application heads. No Mathematica evaluator, automatic head rewriting or alarm.
TOKEN = re.compile(r'\s+|"(?:\\.|[^"\\])*"|<\||\|>|->|==|!=|>=|<=|&&|\|\||'
                   r'[A-Za-z$][A-Za-z0-9$`_]*|(?:\d+(?:\.\d*)?|\.\d+)(?:\*\^[+-]?\d+)?|'
                   r'[\[\]{}(),+*/^<>!|\-\']')


class WLParser:
    def __init__(self, source):
        self.tokens = []
        end = 0
        for match in TOKEN.finditer(source):
            if match.start() != end:
                raise InputError('Wolfram grammar near ' + repr(source[end:end+60]))
            end = match.end()
            if not match.group().isspace():
                self.tokens.append(match.group())
        if source[end:].strip():
            raise InputError('Wolfram trailing grammar: ' + repr(source[end:end+60]))
        self.tokens.append('<END>')
        self.i = 0

    def peek(self):
        return self.tokens[self.i]

    def take(self, expected=None):
        v = self.peek()
        if expected is not None and v != expected:
            raise InputError(f'Wolfram expected {expected!r}, found {v!r}')
        self.i += 1
        return v

    def separated(self, closing):
        values = []
        if self.peek() != closing:
            while True:
                values.append(self.expr())
                if self.peek() != ',':
                    break
                self.take(',')
        self.take(closing)
        return tuple(values)

    def expr(self, minimum=0):
        tok = self.take()
        if tok == '(':
            left = self.expr()
            self.take(')')
        elif tok == '{':
            left = Node('Sequence', self.separated('}'))
        elif tok == '<|':
            entries = self.separated('|>')
            if any(e.head != 'Rule' or e.args[0].head != 'Text' for e in entries):
                raise InputError('association requires text-key rules')
            left = record([(e.args[0].value, e.args[1]) for e in entries])
        elif tok.startswith('"'):
            left = atom('Text', ast.literal_eval(tok))
        elif tok == '!':
            left = Node('Not', (self.expr(6),))
        elif tok in ('+', '-'):
            v = self.expr(35)
            left = v if tok == '+' else Node('Mul', (atom('Number', -1), v))
        elif re.fullmatch(r'\d+(?:\.\d*)?|\.\d+', tok):
            left = atom('Number', tok)
        elif tok in ('True', 'False'):
            left = atom('Boolean', tok)
        elif re.fullmatch(r'[A-Za-z$][A-Za-z0-9$`_]*', tok):
            left = atom('Symbol', tok)
        else:
            raise InputError('unexpected Wolfram token ' + repr(tok))
        precedence = {'|': 3, '->': 2, '||': 4, '&&': 5, '==': 10, '!=': 10,
                      '>': 10, '<': 10, '>=': 10, '<=': 10,
                      '+': 20, '-': 20, '*': 30, '/': 30, '^': 40}
        heads = {'|':'Alternatives','->': 'Rule', '==': 'Equality', '!=': 'Unequality', '>': 'StrictGreaterThan',
                 '<': 'StrictLessThan', '>=': 'GreaterThan', '<=': 'LessThan',
                 '&&': 'And', '||': 'Or', '+': 'Add', '*': 'Mul', '^': 'Pow'}
        while True:
            op = self.peek()
            if op == '[' and 50 >= minimum:
                self.take('[')
                left = Node('Apply', (left, *self.separated(']')))
                continue
            if op == "'" and 50 >= minimum:
                self.take()
                left = Node('Prime', (left,))
                continue
            implicit = (op in ('(', '{', '<|') or op.startswith('"') or
                        re.match(r'^[A-Za-z$\d]', op) is not None)
            prec = precedence.get(op, 30 if implicit else -1)
            if prec < minimum:
                break
            if not implicit:
                self.take()
            else:
                op = '*'
            right = self.expr(prec if op in ('^', '->') else prec + 1)
            if op == '-':
                right = Node('Mul', (atom('Number', -1), right))
                op = '+'
            elif op == '/':
                right = Node('Pow', (right, atom('Number', -1)))
                op = '*'
            head = heads[op]
            # Flatten associative syntax, without evaluating operands.
            args = (*left.args, right) if head in ('Add', 'Mul') and left.head == head else (left, right)
            left = Node(head, args)
        return left


def wl_parse(source):
    p = WLParser(source)
    result = p.expr()
    p.take('<END>')
    return result


@dataclass(frozen=True)
class NameRow:
    py: str
    wl: str
    spec: str
    py_line: int
    wl_line: int
    kind: str


# Explicit injective physical-name bindings. Bound integration coordinates are
# scoped by their own binders, not identified with a second physical coordinate.
# No Part-D OPEN actions occur in this Parts-A--C vocabulary.
NAME_TABLE = tuple(NameRow(*row) for row in [
    ('s9b_r', 'r', 'radial coordinate; 141,175', 668, 45, 'coordinate'),
    ('b', 'b', 'flyby impact parameter; 177', 668, 45, 'coordinate'),
    ('Z_E', 'ZE', 'emitter endpoint; 178-180', 668, 45, 'coordinate'),
    ('Z_R', 'ZR', 'reflector endpoint; 178-180', 668, 45, 'coordinate'),
    ('r_E', 'rE', 'emitter radius; 180', 669, 69, 'coordinate'),
    ('r_R', 'rR', 'reflector radius; 180', 669, 69, 'coordinate'),
    ('s9b_b_min', 'bFar', 'far-zone lower bound; 172-176', 669, 45, 'coordinate'),
    ('GM', 'GM', 'independent orbital mass parameter; 202-204', 670, 46, 'parameter'),
    ('c_0', 'c0', 'asymptotic optical speed; 116', 668, 45, 'parameter'),
    ('K', 'K', 'bulk pressure coefficient; 219-220', 670, 46, 'parameter'),
    ('n', 'n', 'bulk pressure exponent; 219-220', 670, 46, 'parameter'),
    ('m', 'm', 'bulk particle mass; 219-220', 671, 46, 'parameter'),
    ('s', 's', 'optical power exponent; 246', 670, 46, 'parameter'),
    ('rho_0', 'rho0', 'asymptotic bulk number density; 223', 671, 46, 'parameter'),
    ('Phi', 'Phi', 'outward sphere flux; repair directive item 10', 1041, 46, 'parameter'),
    ('s9b_gamma', 'gamma', 'reference-family gamma; 208-217', 849, 68, 'parameter'),
    ('rho_bulk', 'rho', 'local bulk number density; 219,248', 672, 92, 'parameter'),
    ('s9b_fractional_bulk_change', 'fLocal', 'local fractional bulk change; 222', 1019, 94, 'parameter'),
    ('s9b_delta', 'deltaProfile', 'live optical speed shift; 124,141', 397, 50, 'profile'),
    ('s9b_V', 'velocityProfile', 'live radial material velocity; 60,141,175', 398, 50, 'profile'),
    ('s9b_xi', 'xiProfile', 'live graph displacement; 83-87,141', 399, 51, 'profile'),
    ('rho_br_profile', 'rhoBr', 'live brane mass density; 57,69,80', 400, 63, 'profile'),
    ('s9b_f', 'fProfile', 'live fractional bulk density change; 222', 401, 53, 'profile'),
    ('mu_perp', 'muPerp', 'live transverse stiffness; 52-57', 402, 63, 'profile'),
    ('j_n', 'jn', 'live normal mass exchange; 69,80', 673, 67, 'profile'),
    ('s9b_mass_V', 'vForward', 'forward no-loss velocity; repair item 10', 1042, 651, 'profile'),
    ('s9b_velocity_squared', 'yProfile', 'velocity-square reduction; 236-240', 447, 504, 'profile'),
])


def vocabulary():
    grades = ['G%d%d%d' % g for g in product(range(2), range(3), range(2))]
    logs = ('G100', 'G010', 'G020', 'G001')
    obs = ('DEFLECTION', 'RADAR')
    responses = ('CONSTANT', 'FIXED_RATIO', 'POWER')
    branch = ['BRANCH_EXISTENCE', 'PATH_TRAVERSAL', 'BRANCH_TYPE']
    a = ['A_'+o+'_'+g for o in ('DEFLECTION','ROUND_TRIP','ONE_WAY_ER','ONE_WAY_RE','NONRECIPROCAL') for g in grades]
    a += ['A_RADAR_LOG_'+g for g in logs]
    b = ['B_'+o+'_'+v for o in ('GAMMA','RESIDUAL','CONDITION','IMPLIED_JN') for v in obs]+['B_GAMMA_DIFFERENCE']
    c = ['C_'+r+'_'+o+'_'+v for r in responses for o in ('CONDITION','IMPLIED_JN') for v in obs]+['C_'+r+'_N_DEPENDENCE' for r in responses]
    names = branch+a+['A_NONRECIPROCAL_PATH_DEPENDENCE']+b+c
    names += ['R_'+r+'_'+v for r in ('FLOW_ONLY','SPEED_ONLY','TILT_ONLY') for v in ['A_DEFLECTION_'+g for g in grades]+['A_RADAR_LOG_'+g for g in logs]+b]
    names += ['R_BULK_ONLY_C_'+r+'_'+o+'_'+v for r in responses[1:] for o in ('CONDITION','IMPLIED_JN') for v in obs]
    names += ['F_'+r+'_'+v for r in ('FLOW','LIVE') for v in branch+a+[x for x in b if 'IMPLIED_JN' not in x]]
    names += ['F_LIVE_C_'+r+'_CONDITION_'+v for r in responses for v in obs]+['F_MASS_SOLUTION']
    return tuple(names)


SHARED = vocabulary()


def validate_tables(names, joins):
    for field in ('py', 'wl'):
        values = [getattr(row, field) for row in names]
        if len(values) != len(set(values)):
            raise InputError('duplicate name-table '+field+' row')
    if len(joins) != len(set(joins)):
        raise InputError('duplicate join row')


def data(n):
    return [n.head, n.value, list(n.options), [data(a) for a in n.args]]


def census(n):
    """Independent parse census, before role extraction; empty objects count.

    Each terminal and each constructor assumption/option is a separate leaf.
    Arithmetic coverage does not claim to compare assumption options.
    """
    return (sum(census(a) for a in n.args) if n.args else 1) + len(n.options) + (1 if n.head=='Entry' else 0)


@dataclass(frozen=True)
class Ref:
    node: Node
    path: tuple = ()

    def child(self, i):
        return Ref(self.node.args[i], self.path+(i,))

    def fields(self):
        return {a.value: self.child(i).child(0) for i,a in enumerate(self.node.args)} if self.node.head == 'Record' else {}


def leaf_paths(ref, options=True):
    n = ref.node
    found = set().union(*(leaf_paths(ref.child(i), options) for i in range(len(n.args)))) if n.args else {ref.path}
    if n.head=='Entry': found.add(ref.path+('key',))
    if options:
        found.update(ref.path+('option', k) for k,_ in n.options)
    return found


TAG = re.compile(rb'^(PY|WL)_(LOCAL_)?S9B_([A-Z0-9_]+): ?(.*)\r?\n?$')


class Stream:
    """Index one-line framing; retain duplicate occurrences for accounting.

    Payload parsing is lazy per tag; giant repeated ray domains do not require
    keeping the entire production stream's AST in memory at once.
    """
    def __init__(self, path, engine):
        self.path, self.engine = Path(path), engine
        self.index = defaultdict(list)
        self.digest = hashlib.sha256()
        try:
            with self.path.open('rb') as f:
                line = 0
                while True:
                    offset = f.tell()
                    raw = f.readline()
                    if not raw:
                        break
                    line += 1
                    self.digest.update(raw)
                    m = TAG.fullmatch(raw)
                    if not m or m[1].decode().lower() != engine:
                        raise InputError(f'{path}:{line}: unframed/wrong-engine input')
                    key = ('LOCAL/' if m[2] else '')+m[3].decode()
                    self.index[key].append((offset, line))
        except OSError as exc:
            raise InputError(str(exc)) from exc
        if not self.index:
            raise InputError('empty stream '+str(path))

    def read(self, key, occurrence=0):
        with self.path.open('rb') as f:
            f.seek(self.index[key][occurrence][0])
            raw = f.readline()
        text = TAG.fullmatch(raw)[4].decode().strip()
        try:
            return py_parse(text) if self.engine == 'py' else wl_parse(text)
        except (InputError, RecursionError) as exc:
            raise InputError(f'{self.path}:{self.index[key][occurrence][1]}: {exc}') from exc


class Unsupported(Exception):
    """Parsed content for which no mathematical comparison is implemented."""


class UnresolvedEvaluation(Exception):
    """The evaluator has not established a value; this is not undefinedness."""


class UndefinedIntegral(Exception):
    """A proved nonintegrable pole, not a quadrature failure."""


class Quantifier(BooleanFunction):
    """Held scoped ForAll/Exists; never inferred from finitely many probes."""
    @property
    def free_symbols(self):
        return (self.args[2].free_symbols | self.args[3].free_symbols)-set(self.args[1])


class Bindings:
    def __init__(self, rows):
        self.maps = {'py': {r.py:r.py for r in rows}, 'wl': {r.wl:r.py for r in rows}}
        self.kinds = {r.py:r.kind for r in rows}

    def name(self, value, engine):
        return self.maps[engine].get(value, engine+'::'+value)


def call(n, name):
    return n.head == 'Apply' and n.args and n.args[0].head == 'Symbol' and n.args[0].value == name


def wl_head(n):
    if n.head != 'Apply':
        return None, ()
    if n.args[0].head == 'Symbol':
        return n.args[0].value, n.args[1:]
    if call(n.args[0], 'Inactive') and len(n.args[0].args) == 2:
        return 'Inactive/'+n.args[0].args[1].value, n.args[1:]
    return None, ()


RELATIONS = {'Equality':sp.Eq, 'Unequality':sp.Ne, 'StrictGreaterThan':sp.Gt,
             'StrictLessThan':sp.Lt, 'GreaterThan':sp.Ge, 'LessThan':sp.Le}
@dataclass(frozen=True)
class FunctionRow:
    py: str
    wl: tuple[str, ...]
    meaning: str


# Native syntax aliases share ONE row; distinct rows have distinct targets.
# These are mathematical constructors, not additional physical-name bindings.
# Inactive aliases retain their binder and are never executed by a WL engine.
FUNCTION_TABLE = tuple(FunctionRow(*r) for r in [
    ('Add',('Plus',),'sum'), ('Mul',('Times',),'product'),
    ('Pow',('Power','densityPower'),'power, including the source-defined zero exponent convention'),
    ('Equality',('Equal',),'equality'), ('Unequality',('Unequal',),'inequality'),
    ('StrictGreaterThan',('Greater',),'strict greater-than'),
    ('StrictLessThan',('Less',),'strict less-than'),
    ('GreaterThan',('GreaterEqual',),'greater-than or equal'),
    ('LessThan',('LessEqual',),'less-than or equal'),
    ('Sequence',('List',),'ordered tuple'),
    ('Integral',('Integrate','Inactive/Integrate'),'held integral with scoped limits'),
    ('Limit',('Inactive/Limit',),'held limit'),
    ('Contains',('Element',),'set membership'),
    ('ForAll',('ForAll','Inactive/ForAll'),'universal quantifier over a scoped real binder'),
    ('Exists',('Exists','Inactive/Exists'),'existential quantifier over a scoped real binder'),
    ('And',('And',),'logical conjunction'), ('Or',('Or',),'logical disjunction'),
    ('Not',('Not',),'logical negation'), ('Piecewise',('Piecewise',),'ordered conditional branches'),
    ('Rational',('Rational',),'exact rational number'),
    ('Min',('Min',),'minimum'), ('Max',('Max',),'maximum'),
    ('Abs',('Abs',),'absolute value'), ('sin',('Sin',),'sine'),
    ('cos',('Cos',),'cosine'), ('tan',('Tan',),'tangent'),
    ('sec',('Sec',),'secant'), ('acos',('ArcCos',),'principal inverse cosine'),
    ('log',('Log',),'principal logarithm'), ('exp',('Exp',),'exponential'),
    ('sqrt',('Sqrt',),'principal square root'),
    ('re',('Re',),'real part'), ('im',('Im',),'imaginary part'),
])


def function_bindings(rows=None):
    rows=FUNCTION_TABLE if rows is None else rows
    targets=[r.py for r in rows]; aliases=[a for r in rows for a in r.wl]
    if len(set(targets))!=len(targets) or len(set(aliases))!=len(aliases):
        raise InputError('noninjective function-name table')
    return {alias:row.py for row in rows for alias in row.wl}


def mapped_function(head):
    return function_bindings().get(head)


def math_node(n, bindings, engine, scope=None):
    """Translate supported CAS syntax, retaining functions and binder scopes.

    No assumptions are added to real symbols. Source assumptions remain in the
    printed parse and census. Unmapped free names stay engine-qualified.
    Set/quantifier representations are deliberately not reconciled.
    """
    scope = {} if scope is None else scope
    h,a = n.head,n.args
    def cv(x, local=None):
        return math_node(x, bindings, engine, scope if local is None else local)
    if h == 'Number': return sp.Rational(n.value)
    if h == 'Boolean': return sp.true if n.value.lower() == 'true' else sp.false
    if h in ('Name','Symbol','Dummy'):
        constants = {'oo':sp.oo,'Infinity':sp.oo,'pi':sp.pi,'Pi':sp.pi,'E':sp.E,
                     'I':sp.I,'true':sp.true,'false':sp.false, 'Reals':sp.S.Reals,
                     'EmptySet':sp.EmptySet,'ComplexInfinity':sp.zoo,'zoo':sp.zoo,
                     'Undefined':sp.nan,'nan':sp.nan}
        if n.value in scope: return scope[n.value]
        if n.value in constants: return constants[n.value]
        key = bindings.name(n.value,engine)
        if bindings.kinds.get(key) == 'profile':
            raise Unsupported('unapplied live function '+key)
        return sp.Symbol(key)
    if h == 'Apply':
        wh,wa = wl_head(n)
        if mapped_function(wh) is not None:
            h,a = mapped_function(wh),wa
        else:
            head = a[0]
            # Derivative[k][f][x] and f'[x]; neither loses its argument.
            order = None
            if head.head == 'Apply' and call(head.args[0], 'Derivative'):
                if len(head.args[0].args) == 2 and len(head.args) == 2:
                    order,head = int(head.args[0].args[1].value),head.args[1]
            elif head.head == 'Prime':
                order = 0
                while head.head == 'Prime': order,head = order+1,head.args[0]
            if order is not None:
                if head.head != 'Symbol' or len(a) != 2: raise Unsupported('derivative arity')
                key = bindings.name(head.value,engine)
                if bindings.kinds.get(key) != 'profile': raise Unsupported('unmapped derivative function '+key)
                u = sp.Dummy('profile_argument')
                return sp.Subs(sp.Derivative(sp.Function(key)(u),(u,order),evaluate=False),u,cv(a[1]))
            if head.head in ('FunctionName','Symbol'):
                key = bindings.name(head.value,engine)
                if bindings.kinds.get(key) == 'profile':
                    if len(a) != 2: raise Unsupported('profile arity differs')
                    return sp.Function(key)(cv(a[1]))
            raise Unsupported('unsupported applied head '+str(wh or head.head))
    if h in ('Add','Mul'):
        values = [cv(x) for x in a]
        if not all(isinstance(x,sp.Expr) for x in values): raise Unsupported('non-scalar arithmetic')
        return (sp.Add if h == 'Add' else sp.Mul)(*values)
    if h == 'Pow':
        if len(a)!=2: raise Unsupported('power arity')
        base,exponent = cv(a[0]),cv(a[1])
        # Both native definitions make zero powers one, including at zero.
        return sp.Pow(base,exponent)
    if h == 'Rational': return sp.Rational(*(cv(x) for x in a))
    if h in RELATIONS:
        if len(a)!=2: raise Unsupported('nonbinary relation')
        if a[0].head in RELATIONS:
            return sp.And(cv(a[0]),RELATIONS[h](cv(a[0].args[-1]),cv(a[1]),evaluate=False))
        return RELATIONS[h](*(cv(x) for x in a),evaluate=False)
    if h in ('And','Or','Not'):
        values = [cv(x) for x in a]
        if not all(isinstance(x,Boolean) for x in values): raise Unsupported('nonboolean connective')
        return getattr(sp,h)(*values,evaluate=False)
    if h == 'Contains':
        left,right = cv(a[0]),cv(a[1])
        if isinstance(left,sp.Tuple): return sp.And(*(sp.Contains(x,right,evaluate=False) for x in left))
        return sp.Contains(left,right,evaluate=False)
    if h in ('ForAll','Exists'):
        variables=a[0].args if a[0].head=='Sequence' else (a[0],)
        local=dict(scope)
        for i,v in enumerate(variables): local[v.value]=sp.Symbol('_quant_'+str(len(scope)+i),real=True)
        restriction=cv(a[1],local) if len(a)==3 else sp.true
        return Quantifier(sp.Integer(h=='ForAll'),sp.Tuple(*(local[v.value] for v in variables)),restriction,cv(a[-1],local))
    if h == 'ConditionSet':
        variables=a[0].args if a[0].head=='Sequence' else (a[0],)
        local=dict(scope)
        for i,v in enumerate(variables): local[v.value]=sp.Symbol('_set_'+str(len(scope)+i),real=True)
        bound=tuple(local[v.value] for v in variables)
        return sp.ConditionSet(bound[0] if len(bound)==1 else sp.Tuple(*bound),cv(a[1],local),cv(a[2]))
    if h in ('ProductSet','Intersection','Union','Complement'):
        return getattr(sp,h)(*(cv(x) for x in a),**({'evaluate':False} if h in ('Intersection','Union') else {}))
    if h in ('Sequence','Alternatives','FiniteSet'): return (sp.Tuple if h in ('Sequence','Alternatives') else sp.FiniteSet)(*(cv(x) for x in a))
    if h == 'Interval': return sp.Interval(*(cv(x) for x in a),**{k:ast.literal_eval(v) for k,v in n.options})
    if h in ('sin','cos','tan','sec','log','exp','sqrt','Abs','Min','Max','re','im','acos','cbrt'):
        return getattr(sp,h)(*(cv(x) for x in a))
    if h == 'Limit':
        if engine=='wl':
            if len(a)!=2 or a[1].head!='Rule': raise Unsupported('limit binder grammar')
            variable,destination=a[1].args
        else:
            if len(a) not in (3,4): raise Unsupported('limit arity')
            variable,destination=a[1:3]
        local=dict(scope); bound=sp.Symbol('_limit_'+str(len(scope)))
        local[variable.value]=bound
        direction=a[3].value if engine=='py' and len(a)==4 else '+'
        return sp.Limit(cv(a[0],local),bound,cv(destination),dir=direction)
    if h == 'Integral':
        if len(a)<2: raise Unsupported('missing integral limits')
        local = dict(scope); limits = []
        for spec in a[1:]:
            if spec.head != 'Sequence' or len(spec.args)!=3 or spec.args[0].head not in ('Symbol','Dummy'):
                raise Unsupported('unsupported integral binder')
            rawvar,lo,hi = spec.args
            var = sp.Symbol('_bound_'+str(len(local)))
            lower,upper = cv(lo,local),cv(hi,local)
            local[rawvar.value] = var
            limits.append((var,lower,upper))
        return sp.Integral(cv(a[0],local),*limits)
    if h == 'Derivative':
        specs=[]
        for x in a[1:]:
            parts=x.args if x.head=='Sequence' else (x,)
            specs.append((cv(parts[0]),int(cv(parts[1]))) if len(parts)==2 else cv(parts[0]))
        return sp.Derivative(cv(a[0]),*specs,evaluate=False)
    if h == 'Subs':
        def members(x): return x.args if x.head=='Sequence' else (x,)
        variables,points = members(a[1]),members(a[2]); local=dict(scope)
        if len(variables)!=len(points): raise Unsupported('substitution arity differs')
        for i,v in enumerate(variables): local[v.value]=sp.Symbol('_sub_'+str(len(scope)+i))
        return sp.Subs(cv(a[0],local),tuple(local[v.value] for v in variables),tuple(cv(x) for x in points))
    if h == 'Piecewise':
        pairs=[]
        if engine=='wl':
            if not a or a[0].head!='Sequence': raise Unsupported('piecewise grammar')
            branches=a[0].args
            fallback = cv(a[1]) if len(a)>1 else sp.S.Zero
        else: branches=a; fallback=None
        for pair in branches:
            if pair.head not in ('Sequence','ExprCondPair') or len(pair.args)!=2: raise Unsupported('piecewise branch')
            pairs.append((cv(pair.args[0]),cv(pair.args[1])))
        if fallback is not None: pairs.append((fallback,sp.true))
        return sp.Piecewise(*pairs,evaluate=False)
    raise Unsupported('unsupported constructor '+h)
@dataclass
class View:
    ref: Ref
    role: str = 'value'
    children: dict | None = None


def is_text(n, value): return n.head=='Text' and n.value==value

def attached(ref):
    return ((ref.node.head=='Sequence' and len(ref.node.args)==2 and is_text(ref.node.args[0],'ATTACHED_DOMAIN'))
            or (ref.node.head=='Record' and set(ref.fields())=={'ATTACHED_DOMAIN'}))


def domain_ref(ref):
    return ref.fields()['ATTACHED_DOMAIN'] if ref.node.head=='Record' else ref.child(1)


def gate_pair(ref):
    if ref.node.head=='Sequence' and len(ref.node.args)==2 and attached(ref.child(0)):
        return domain_ref(ref.child(0)),ref.child(1)
    # Tuple((ATTACHED_DOMAIN, condition), (text, value), ...) is a record.
    fields=ref.fields()
    if 'ATTACHED_DOMAIN' in fields:
        rest={k:v for k,v in fields.items() if k!='ATTACHED_DOMAIN'}
        if len(rest)==1:
            # A gate's second element can itself be a text/value pair. Keep
            # that pair's original record entry and key in the value view.
            return fields['ATTACHED_DOMAIN'],View(ref,'value',{k:extract(v,'py') for k,v in rest.items()})
    return None


def extract(ref, engine, role='value'):
    """Source-defined wrappers only; no algebraic representation reconciliation.

    PY 76-93,605-611; WL 270-272,618-623,661 define these roles.
    All wrapper metadata remains in raw operands and parsed leaf counts. Only
    compared descendants contribute to comparison coverage.
    """
    n=ref.node; fields=ref.fields()
    if attached(ref): return extract(domain_ref(ref),engine,'domain')
    if engine=='py' and n.head=='Sequence':
        a=n.args
        if len(a)==3 and is_text(a[0],'APPLICABLE_OBJECT') and attached(ref.child(1)):
            return View(ref,role,{'value':extract(ref.child(2),engine),
                                  'applicability':extract(ref.child(1),engine)})
        if len(a)==2 and all(gate_pair(ref.child(i)) is not None for i in range(2)):
            result={}
            for i,name in enumerate(('on','outside')):
                d,v=gate_pair(ref.child(i))
                result[name]=View(ref.child(i),'case',{'domain':extract(d,engine,'domain'),
                    'value':v if isinstance(v,View) else extract(v,engine)})
            return View(ref,role,result)
        if len(a)==4 and is_text(a[0],'LOCAL_BULK_DENSITY_PARTS'):
            parts={}
            for i,name in ((1,'positive'),(2,'nonpositive')):
                labelref=ref.child(i)
                if labelref.node.head!='Sequence' or len(labelref.node.args)!=2:
                    raise Unsupported('malformed density label structure')
                objectref=ref.child(3)
                if objectref.node.head!='Sequence' or len(objectref.node.args)!=2:
                    raise Unsupported('malformed density object structure')
                parts[name]=View(labelref,'density_part',{
                    'label':View(labelref.child(0),'label'),
                    'label_condition':extract(labelref.child(1),engine,'domain'),
                    'object':extract(objectref.child(1),engine)})
            return View(ref,role,parts)
    if engine=='wl' and fields:
        if 'Object' in fields and 'Premise' in fields:
            inner=extract(fields['Object'],engine,role)
            # The complete forward premise is separately visible and unpaired.
            return View(ref,role,{'value':inner,'premise':extract(fields['Premise'],engine,'domain')})
        if 'OnDomain' in fields and 'OutsideDomain' in fields:
            result={}
            for source,dest in (('OnDomain','on'),('OutsideDomain','outside')):
                item=fields[source]
                if not call(item.node,'ConditionalExpression') or len(item.node.args)!=3:
                    return View(ref,role,{'representation':View(ref,role)})
                result[dest]=View(item,'case',{'value':extract(item.child(1),engine),
                                              'domain':extract(item.child(2),engine,'domain')})
            for key,item in fields.items():
                if key not in ('OnDomain','OutsideDomain'): result['metadata/'+key]=extract(item,engine)
            return View(ref,role,result)
        if 'DensityParts' in fields:
            parts={}; rows=fields['DensityParts']
            if rows.node.head!='Sequence' or len(rows.node.args)!=2:
                raise Unsupported('density part count differs')
            for i,name in enumerate(('positive','nonpositive')):
                item=rows.child(i); f=item.fields()
                if not {'Label','LabelCondition','Object'}<=f.keys(): raise Unsupported('density fields absent')
                parts[name]=View(item,'density_part',{'label':View(f['Label'],'label'),
                    'label_condition':extract(f['LabelCondition'],engine,'domain'),
                    'object':extract(f['Object'],engine)})
            for key,item in fields.items():
                if key!='DensityParts': parts['metadata/'+key]=extract(item,engine,'domain' if key=='SuppliedBulkDensityPremise' else 'value')
            return View(ref,role,parts)
    if engine=='wl' and call(n,'ConditionalExpression') and len(n.args)==3:
        return View(ref,role,{'value':extract(ref.child(1),engine),'domain':extract(ref.child(2),engine,'domain')})
    if n.head=='Record':
        return View(ref,role,{k:extract(v,engine,'domain' if k in DOMAIN_KEYS else role) for k,v in fields.items()})
    if n.head in ('Sequence','Matrix'):
        return View(ref,role,{str(i):extract(ref.child(i),engine,role) for i in range(len(n.args))})
    return View(ref,role)


DOMAIN_KEYS={'Domain','BranchDomain','AlongRayDomain','ConstantDomain','GeneralConstantDomain',
             'NoRealVelocityDomain','LabelCondition','SuppliedBulkDensityPremise'}
# Labels are role bindings, compared along with the conditions attached to them.
LABEL_TABLE={'RHO_POSITIVE':'rho > 0','RHO_NONPOSITIVE':'rho <= 0'}


# Rational points are probes, not premises. Distinct endpoint distances and
# independent GM/Phi values are deliberate; no orbital/profile identification.
POINTS=tuple({k:sp.Rational(v) for k,v in row.items()} for row in [
    dict(s9b_r='11/3',b='7/3',Z_E='101/3',Z_R='137/3',c_0='19/8',GM='1/101',
         Phi='2/103',K='5/7',n='2',m='3/2',s='1/2',rho_0='11/7',s9b_gamma='2/3',s9b_b_min='1/7'),
    dict(s9b_r='17/5',b='11/5',Z_E='151/5',Z_R='173/5',c_0='29/12',GM='2/107',
         Phi='3/109',K='7/11',n='3',m='5/3',s='2',rho_0='13/11',s9b_gamma='3/5',s9b_b_min='1/11'),
    dict(s9b_r='23/7',b='13/7',Z_E='211/7',Z_R='263/7',c_0='31/16',GM='3/113',
         Phi='5/127',K='11/13',n='1/2',m='7/5',s='-1',rho_0='17/13',s9b_gamma='5/7',s9b_b_min='1/13')])
PROFILE_CITATIONS = SPEC+':52-57,80,124-149,175-176,219-224; repair directive item 9 (density signs; no numeric f bound)'
PROFILE_NAMES=('inverse','shifted','exponential','optical_zero','delta_zero','velocity_zero','tilt_zero',
               'bulk_zero_at_point','bulk_negative_at_point','brane_zero_at_point','stiffness_zero_at_point')


def profile_formulas(family,point):
    return _profile_formulas(family,tuple(sorted(point.items())))


@lru_cache(maxsize=256)
def _profile_formulas(family,point_items):
    point=dict(point_items)
    x=sp.Symbol('_profile_x'); c=point['c_0']; anchor=point['s9b_r']
    if family=='exponential':
        d,v,xi,f=sp.exp(-x)/107,c*sp.exp(-x)/23,sp.exp(-x)/29,sp.exp(-x)/31
    elif family=='shifted':
        d,v,xi,f=1/(103*(x+1)),c/(19*sp.sqrt(x+1)),sp.log(1+x)/23,-1/(29*(x+1))
    else:
        d,v,xi,f=1/(101*x),c/(17*sp.sqrt(x)),2*sp.sqrt(x)/19,1/(23*x)
    rho=1+1/(13*x)
    if family=='optical_zero': d,v,xi,f=map(sp.sympify,(0,0,0,0))
    if family=='delta_zero': d=sp.S.Zero
    if family=='velocity_zero': v=sp.S.Zero
    if family=='tilt_zero': xi=sp.S.Zero
    if family=='bulk_zero_at_point': f=-anchor/x
    if family=='bulk_negative_at_point': f=-2*anchor/x
    if family=='brane_zero_at_point': rho=(x-anchor)**2/(1+(x-anchor)**2)
    if family=='stiffness_zero_at_point': d=(1+d)*(x-anchor)**2/(1+(x-anchor)**2)-1
    mu=rho*c**2*(1+d)**2
    exchange=-sp.diff(x*x*rho*v,x)/(x*x)
    return x,dict(s9b_delta=d,s9b_V=v,s9b_xi=xi,rho_br_profile=rho,s9b_f=f,mu_perp=mu,
                   j_n=exchange,s9b_mass_V=point['Phi']/(4*sp.pi*x*x*rho),s9b_velocity_squared=v*v)


def configurations(expressions=()):
    # Generic points plus individual zero/sign strata for every emitted free
    # parameter used by this object. Uncovered joint strata are reported below.
    for i,p in enumerate(POINTS):
        for profile in PROFILE_NAMES:
            yield f'p{i}/{profile}',p,profile
    free=set().union(*(e.free_symbols for e in expressions if isinstance(e,sp.Basic)))
    p=POINTS[0]
    for name in ('GM','Phi','K','n','s','py::s9b_velocity_integration_constant','wl::velocityIntegrationConstant','rho_bulk','s9b_fractional_bulk_change'):
        if sp.Symbol(name) in free:
            for val in (sp.S.Zero,sp.S.One,-sp.S.One,sp.Rational(1,2)):
                yield name+'='+str(val),distinct_probe(p,name,val),'inverse'


def distinct_probe(point,name,value):
    result=dict(point,**{name:value})
    occupied={value}
    for key in sorted(result):
        if key==name: continue
        candidate=result[key]
        while candidate in occupied: candidate+=sp.Rational(1,211)
        result[key]=candidate; occupied.add(candidate)
    return result


def substitute_profiles(expr,point,profile):
    x,forms=profile_formulas(profile,point)
    def replace_fn(a):
        key=a.func.__name__
        if key not in forms: raise Unsupported('test profile unavailable for '+key)
        if len(a.args)!=1: raise Unsupported('profile argument arity')
        return forms[key].subs(x,a.args[0])
    @lru_cache(maxsize=None)
    def prepare(node):
        # Finish each derivative and its Subs before rebuilding the parent.
        # Separate global passes temporarily constructed products of Subs(0,
        # ..., point) at different points, whose SymPy ordering raises while
        # comparing tuple-valued internal content. No such intermediate is
        # needed to evaluate the original sampled expression.
        args=tuple(prepare(a) for a in node.args)
        if isinstance(node,AppliedUndef):
            return replace_fn(node.func(*args))
        if isinstance(node,sp.Derivative):
            value=sp.Derivative(*args,evaluate=False)
            return value if value.has(sp.Integral) else evaluate_profile_derivative(value)
        if isinstance(node,sp.Subs):
            value=sp.Subs(*args)
            return evaluate_substitution(value)
        return node.func(*args) if args!=node.args else node
    value=prepare(expr)
    prepare.cache_clear()
    substitutions={sp.Symbol(k):v for k,v in point.items()}
    substitutions[sp.Symbol('r_E')]=sp.sqrt(point['b']**2+point['Z_E']**2)
    substitutions[sp.Symbol('r_R')]=sp.sqrt(point['b']**2+point['Z_R']**2)
    # SymPy ConditionSet._eval_subs can leave external parameters unchanged.
    # Set binders were alpha-renamed by the reader, so this free-parameter
    # replacement cannot capture their variables or rewrite derivatives in
    # scalar operands outside the sets.
    value=value.replace(lambda n:isinstance(n,sp.ConditionSet),lambda n:n.xreplace(substitutions))
    return value.subs(substitutions, simultaneous=True)


def evaluate_profile_derivative(value):
    return value.doit()


def evaluate_substitution(value):
    return value if value.has(sp.Derivative,sp.Integral) else value.doit(deep=False)


def math_text(value):
    """Print stored mathematical content without evaluating it to order terms.

    SymPy's default Add printer can call complex() on numerical integrals and
    run quadrature while formatting an undefined sample. Stored argument order
    retains the full expression without that accidental second computation.
    """
    return sp.sstr(value,order='none')


def numeric_value(expr,point,profile):
    return _numeric_value(expr,tuple(sorted(point.items())),profile)


def sampled_truth(value):
    """Decide only established truth, including scoped sets and quantifiers."""
    if value in (sp.true,sp.false): return {'truth':bool(value)}
    if isinstance(value,(sp.And,sp.Or)):
        results=[sampled_truth(v) for v in value.args]
        decisive=False if isinstance(value,sp.And) else True
        if any(r.get('truth') is decisive for r in results): return {'truth':decisive}
        if all('truth' in r for r in results): return {'truth':not decisive}
        return {'unresolved':'; '.join(sorted({r['unresolved'] for r in results if 'unresolved' in r}))}
    if isinstance(value,sp.Not):
        result=sampled_truth(value.args[0])
        return {'truth':not result['truth']} if 'truth' in result else result
    if isinstance(value,Quantifier):
        universal,variables,restriction,predicate=value.args
        counter=sp.And(restriction,sp.Not(predicate)) if universal else sp.And(restriction,predicate)
        try:
            empty=set_empty(sp.ConditionSet(variables[0] if len(variables)==1 else variables,counter,sp.S.Reals if len(variables)==1 else sp.ProductSet(*([sp.S.Reals]*len(variables)))))
        except UnresolvedEvaluation as exc: return {'unresolved':str(exc)}
        return {'truth':empty if universal else not empty} if empty is not None else {'unresolved':'quantified condition could not be decided over its full bound domain'}
    if isinstance(value,(sp.Equality,sp.Unequality)) and (isinstance(value.lhs,sp.Set) or isinstance(value.rhs,sp.Set)):
        if value.lhs==value.rhs: return {'truth':isinstance(value,sp.Equality)}
        other=value.rhs if value.lhs==sp.EmptySet else value.lhs if value.rhs==sp.EmptySet else None
        try:
            empty=set_empty(other) if other is not None else None
        except UnresolvedEvaluation as exc: return {'unresolved':str(exc)}
        if empty is not None: return {'truth':empty if isinstance(value,sp.Equality) else not empty}
        return {'unresolved':'set equality/emptiness could not be decided'}
    if value.free_symbols: return {'unresolved':'free names '+','.join(sorted(map(str,value.free_symbols)))}
    if value.has(sp.Integral,sp.Derivative,sp.Subs,sp.Limit):
        return {'unresolved':'predicate contains an unevaluated functional'}
    if value.has(sp.nan,sp.zoo): return {'unresolved':'predicate contains undefined scalar content'}
    reduced=sp.simplify(value)
    if reduced in (sp.true,sp.false): return {'truth':bool(reduced)}
    return {'unresolved':'predicate truth not decided: '+math_text(reduced)}


def set_empty(value):
    if value==sp.EmptySet: return True
    if isinstance(value,(sp.Interval,sp.FiniteSet)) or value==sp.S.Reals:
        return value.is_empty
    if isinstance(value,sp.ProductSet):
        states=[set_empty(v) for v in value.args]
        if True in states: return True
        return False if all(v is False for v in states) else None
    if isinstance(value,sp.ConditionSet):
        if set_empty(value.base_set) is True: return True
        variables=value.sym if isinstance(value.sym,sp.Tuple) else (value.sym,)
        condition=value.condition
        if condition is sp.false: return True
        if condition is sp.true: return set_empty(value.base_set)
        if isinstance(value.base_set,sp.FiniteSet):
            truths=[sampled_truth(condition.subs(value.sym,x)) for x in value.base_set]
            if any(t.get('truth') is True for t in truths): return False
            if all(t.get('truth') is False for t in truths): return True
        # Exact univariate polynomial/rational decision, not finite sampling.
        # More general quantified ray constructions remain explicit coverage.
        if len(variables)==1 and not (condition.free_symbols-set(variables)) and not condition.has(AppliedUndef,sp.Integral,sp.Derivative,sp.Subs,Quantifier,sp.Contains):
            x=variables[0]
            try:
                if all((r.lhs-r.rhs).is_rational_function(x) for r in condition.atoms(Relational)):
                    # as_set's general-coefficient solver can expand huge
                    # radicals and launch numerical assumption checks. This
                    # decision procedure supports rational coefficients;
                    # other coefficient fields get a sampled unknown, never
                    # a truth inferred from a finite set of point probes.
                    for relation in condition.atoms(Relational):
                        numerator,denominator=sp.together(relation.lhs-relation.rhs).as_numer_denom()
                        for polynomial in (numerator,denominator):
                            if not all(coefficient.is_Rational for coefficient in sp.Poly(polynomial,x).all_coeffs()):
                                raise UnresolvedEvaluation('quantified/set truth requires root isolation with non-rational coefficients; not implemented')
                    solution=condition.as_set()
                    intersection=sp.Intersection(solution,value.base_set)
                    return intersection.is_empty
            except (NotImplementedError,ValueError,TypeError): pass
    return None


@lru_cache(maxsize=4096)
def _numeric_value(expr,point_items,profile):
    try:
        value=substitute_profiles(expr,dict(point_items),profile)
    except TypeError as exc:
        # SymPy refuses to construct an order comparison with NaN or a
        # non-real number. That is an undecidable sampled predicate, not an
        # operational input error and not a truth value. Other TypeErrors
        # remain bugs and propagate.
        if str(exc)=='Invalid NaN comparison' or str(exc).startswith('Invalid comparison of non-real'):
            return {'unresolved':'sampled predicate: '+str(exc)}
        raise
    if isinstance(value,Boolean) and not isinstance(value,sp.Expr): return sampled_truth(value)
    if value.has(sp.Limit):
        return {'unresolved':'held limit: limiting value not established'}
    if value.free_symbols:
        return {'unresolved':'free names '+','.join(sorted(map(str,value.free_symbols)))}
    if value.has(sp.nan,sp.zoo): return {'undefined':math_text(value),'reason':'exact substitution is undefined'}
    if isinstance(value,sp.Set): return {'unresolved':'set is not a scalar truth value'}
    if value.has(sp.Integral,sp.Derivative,sp.Subs):
        try:
            low=numerical_expression(value,18); high=numerical_expression(value,30)
        except (UnresolvedEvaluation,ValueError,NotImplementedError,ZeroDivisionError) as exc:
            return {'unresolved':'quadrature '+type(exc).__name__+': '+str(exc)}
        except UndefinedIntegral as exc:
            return {'undefined':str(exc)}
        with mp.workdps(50):
            discrepancy=abs(high-low)
            if not mp.isfinite(high) or discrepancy>mp.mpf('1e-14')*abs(high):
                return {'unresolved':'independent quadrature precision check did not establish 12 printed digits',
                        'evaluations':evaluation_pair(low,high)}
            return {'evaluations':evaluation_pair(low,high),'approx_12':mp.nstr(high,12),'precision_difference':mp.nstr(discrepancy,5),
                    'method':'convergence screening; quadrature precision/panel checks; complete-functional derivative stencil/step checks where needed',
                    'reliable_digits':12,'reliability':'numerical evidence, not a certified error bound'}
    if not isinstance(value,sp.Expr): return {'unresolved':'nonscalar sample '+math_text(value)}
    if value.is_finite is False: return {'undefined':math_text(value),'reason':'nonfinite exact scalar'}
    return {'exact':math_text(value),'approx_30':math_text(value.evalf(30)),
            'evaluations':evaluation_pair(*[scalar_number(value,d) for d in (48,60)])}


def scalar_number(value,digits):
    with mp.workdps(digits):
        return mp.mpc(str(sp.re(value).evalf(digits)),str(sp.im(value).evalf(digits)))


def evaluation_pair(low,high):
    return [[mp.nstr(mp.re(x),65),mp.nstr(mp.im(x),65)] for x in (low,high)]


def residual_scale(a,b):
    return max(abs(x) for x in (*a,*b))


def numerical_residual(left,right):
    """Subtract computed operand evaluations; never gate on residual size.

    Values are estimates at two working precisions. The reliability statement
    is ABSOLUTE at operand scale, not significant digits of a near-zero result.
    Even an unreliable computed difference retains both estimates.
    """
    if 'evaluations' not in left or 'evaluations' not in right:
        return {'unresolved':'operand evaluations unavailable for subtraction'}
    with mp.workdps(70):
        a=[mp.mpc(*v) for v in left['evaluations']]
        b=[mp.mpc(*v) for v in right['evaluations']]
        low,high=(x-y for x,y in zip(a,b))
        scale=residual_scale(a,b); discrepancy=abs(high-low)
        reliable=all(mp.isfinite(x) for x in (*a,*b,low,high)) and discrepancy<=mp.mpf('1e-14')*scale
        return {'computed_18':mp.nstr(low,50),'computed_30':mp.nstr(high,60),
                'precision_difference':mp.nstr(discrepancy,12),'operand_scale':mp.nstr(scale,30),
                'absolute_tolerance':mp.nstr(mp.mpf('1e-14')*scale,30),
                'reliability_established':reliable,
                'reliability':'absolute precision consistency at operand scale; numerical evidence, not a certified bound'}


def trig_poles(function,x,lo,hi):
    """Locate the real poles of affine-argument trig functions exactly.

    The floor/ceiling bounds enumerate the entire closed interval, not a set
    of numerical probes. Non-affine arguments remain an explicit limitation.
    Also inspect reciprocal sin/cos written as powers instead of tan/sec.
    """
    roots=set()
    for node in function.atoms(sp.tan,sp.sec,sp.cot,sp.csc,sp.sin,sp.cos):
        if node.func in (sp.sin,sp.cos) and not any(
                p.exp.is_negative and p.base==node for p in function.atoms(sp.Pow)):
            continue
        argument=node.args[0]
        if x not in argument.free_symbols: continue
        slope=sp.diff(argument,x); intercept=argument-slope*x
        if x in slope.free_symbols or x in intercept.free_symbols or slope.is_real is not True or slope.is_zero is not False:
            raise UnresolvedEvaluation('trigonometric pole locations: affine real argument not established for '+math_text(node))
        a,b=(argument.subs(x,v) for v in (lo,hi))
        if slope.is_negative: a,b=b,a
        offset=sp.Rational(1,2) if node.func in (sp.tan,sp.sec,sp.cos) else sp.S.Zero
        first=sp.ceiling(a/sp.pi-offset); last=sp.floor(b/sp.pi-offset)
        if not (first.is_Integer and last.is_Integer):
            raise UnresolvedEvaluation('trigonometric pole index bounds not established for '+math_text(node))
        roots.update((sp.pi*(k+offset)-intercept)/slope for k in range(int(first),int(last)+1))
    return roots


def endpoint_power_bound(expr,t):
    """Return q when |expr| = O(t**q), t -> 0+, or leave it unresolved.

    A nonzero leading coefficient may contain exp/log of constants. An
    exponential with a finite limit is bounded even when its own series
    implementation cannot supply a leading term. Products/sums of these
    bounds cover algebraic singularities multiplied by analytic profiles.
    """
    if not expr.has(t):
        if expr.is_finite is True: return sp.S.Zero
        raise UnresolvedEvaluation('endpoint coefficient finiteness not established')
    if expr.is_Mul: return sum(endpoint_power_bound(a,t) for a in expr.args)
    if expr.is_Add: return min(endpoint_power_bound(a,t) for a in expr.args)
    try:
        coefficient,power=expr.leadterm(t,cdir=1)
        if not coefficient.has(t) and coefficient.is_finite is True and coefficient.is_zero is False:
            return power
    except (NotImplementedError,sp.PoleError,ValueError):
        pass
    try:
        limit=sp.limit(expr,t,0,dir='+')
    except (NotImplementedError,sp.PoleError,ValueError) as exc:
        raise UnresolvedEvaluation('endpoint asymptotic bound unavailable for '+math_text(expr)) from exc
    if limit.is_finite is True: return sp.S.Zero
    raise UnresolvedEvaluation('endpoint power or boundedness not established for '+math_text(expr))


def endpoint_convergence(function,x,endpoint,direction):
    t=sp.Dummy('endpoint',positive=True)
    transformed=(function.subs(x,1/t)/t**2 if endpoint==sp.oo else
                 function.subs(x,-1/t)/t**2 if endpoint==-sp.oo else function.subs(x,endpoint+direction*t))
    # First try the complete expression: cancellation and exponential decay
    # can establish boundedness even when individual algebraic factors grow.
    try:
        bound=endpoint_power_bound(transformed,t)
        if (bound>-1) is sp.true: return
        detail='absolute power bound '+math_text(bound)+' is not greater than -1'
    except UnresolvedEvaluation as exc:
        detail=str(exc)
    try:
        if sp.limit(transformed,t,0,dir='+').is_finite is True: return
    except (NotImplementedError,sp.PoleError,ValueError):
        pass
    raise UnresolvedEvaluation('endpoint convergence at '+math_text(endpoint)+': '+detail)


def real_range_bound(expr,x,lo,hi):
    """Elementary interval enclosure, with exact symbolic endpoints.

    Used only as a sufficient nonzero-denominator certificate. In particular,
    shifted profiles composed with sec on a pole-free first-quadrant range
    do not become unsupported just because their denominators are not polynomials.
    Dependency overestimation can refuse a certificate; it cannot admit a pole.
    """
    if expr==x: return lo,hi
    if not expr.has(x):
        if expr.is_real is True and expr.is_finite is True: return expr,expr
        raise UnresolvedEvaluation('real interval coefficient not established')
    ranges=[real_range_bound(a,x,lo,hi) for a in expr.args]
    if expr.is_Add: return sum(a for a,b in ranges),sum(b for a,b in ranges)
    if expr.is_Mul:
        a=b=sp.S.One
        for c,d in ranges:
            values=(a*c,a*d,b*c,b*d); a,b=sp.Min(*values),sp.Max(*values)
        return a,b
    if expr.is_Pow:
        a,b=ranges[0]; power=expr.exp
        if not power.has(x) and power.is_real and ((a>0) is sp.true or (a==0 and power.is_positive)):
            values=(a**power,b**power)
            return sp.Min(*values),sp.Max(*values)
        if power.is_Integer and power.is_nonnegative:
            values=[a**power,b**power]
            if (a<=0) is sp.true and (b>=0) is sp.true: values.append(sp.S.Zero)
            return sp.Min(*values),sp.Max(*values)
    if len(ranges)==1:
        a,b=ranges[0]
        if expr.func==sp.exp: return sp.exp(a),sp.exp(b)
        if expr.func==sp.log and (a>0) is sp.true: return sp.log(a),sp.log(b)
        if (a>=0) is sp.true and (b<sp.pi/2) is sp.true:
            if expr.func in (sp.sin,sp.tan,sp.sec): return expr.func(a),expr.func(b)
            if expr.func==sp.cos: return sp.cos(b),sp.cos(a)
    raise UnresolvedEvaluation('real interval enclosure not established for '+math_text(expr))


def integral_screen(integral):
    """Conservative absolute-convergence screen before any quadrature.

    All algebraic denominator roots are inspected, including interior roots.
    Unknown singular sets or endpoint orders are unresolved, never finite values.
    No Cauchy principal value replaces an improper integral.
    """
    if len(integral.limits)!=1: raise UnresolvedEvaluation('multidimensional integral')
    (x,lo,hi),=integral.limits
    f=integral.function
    if f.has(sp.Integral): raise UnresolvedEvaluation('nested integral convergence not established')
    if lo==hi: return
    if (lo>hi) is sp.true: lo,hi=hi,lo
    roots=trig_poles(f,x,lo,hi)
    for power in f.atoms(sp.Pow):
        if power.exp.is_negative:
            base=power.base
            if x not in base.free_symbols: continue
            if not base.is_polynomial(x):
                # exp is entire and has no zeros; other denominators require
                # a separate proof, not successful quadrature nodes.
                if base.func==sp.exp: continue
                if base.func in (sp.sin,sp.cos): continue # zeros screened above
                if base.is_rational_function(x):
                    # Isolate numerator zeros of a rational denominator.
                    # Its own poles occur among the nested negative powers
                    # and are screened independently; no value is changed.
                    base=sp.cancel(base).as_numer_denom()[0]
                else:
                    try:
                        a,b=real_range_bound(base,x,lo,hi)
                        if (a>0) is sp.true or (b<0) is sp.true: continue
                    except UnresolvedEvaluation:
                        pass
                    raise UnresolvedEvaluation('nonpolynomial denominator zero set not established for '+math_text(base))
            try:
                polynomial=sp.Poly(base,x)
                roots.update(polynomial.real_roots())
            except (sp.PolynomialError,NotImplementedError):
                raise UnresolvedEvaluation('denominator real roots not established')
    for root in roots:
        inside=sp.And(root>lo,root<hi)
        if inside is sp.true:
            # Cancellation across summands and fractional singularities can
            # be removable/integrable. We do not call them divergent without
            # proving it; neither do we print an unvalidated finite value.
            raise UnresolvedEvaluation('interior denominator zero at '+math_text(root)+'; local integrability not established')
        if inside not in (sp.true,sp.false): raise UnresolvedEvaluation('singularity location undecided')
    for endpoint,sign in ((lo,1),(hi,-1)):
        if endpoint not in roots and endpoint not in (sp.oo,-sp.oo): continue
        endpoint_convergence(f,x,endpoint,sign)


@lru_cache(maxsize=2048)
def numerical_integral(integral,digits):
    if integral.function.has(sp.Integral):
        evaluated=integral.doit()
        if evaluated.has(sp.Integral): raise UnresolvedEvaluation("nested integral convergence not established")
        with mp.workdps(digits+30): return mp.mpc(str(sp.re(evaluated.evalf(digits+30))),str(sp.im(evaluated.evalf(digits+30))))
    integral_screen(integral)
    (variable,lower,upper),=integral.limits
    function=sp.lambdify(variable,integral.function,modules='mpmath',cse=True)
    bounds=sp.lambdify((),(lower,upper),modules='mpmath',cse=True)
    with mp.workdps(2*digits+30):
        lo,hi=bounds()
        first,error=mp.quad(function,[lo,hi],method='tanh-sinh',error=True)
    with mp.workdps(2*digits+50):
        lo,hi=bounds()
        middle=(lo+hi)/2 if mp.isfinite(lo) and mp.isfinite(hi) else lo+1 if mp.isfinite(lo) else hi-1 if mp.isfinite(hi) else mp.mpf(0)
        second,error2=mp.quad(function,[lo,middle,hi],method='tanh-sinh',error=True)
        tolerance=mp.power(10,-digits+4)*abs(second)
        if not quadrature_consistent(first,second,error,error2,tolerance):
            raise UnresolvedEvaluation('quadrature error/precision/panel consistency check did not converge')
        return +second


def quadrature_consistent(first,second,error,error2,tolerance):
    return all(mp.isfinite(v) for v in (first,second,error,error2)) and max(error,error2,abs(second-first))<=tolerance


@lru_cache(maxsize=32)
def derivative_stencil(order):
    radius=order+2
    nodes=tuple(range(-radius,radius+1))
    weights=tuple(sp.finite_diff_weights(order,nodes,0)[order][-1])
    rate=2*radius+2-order if order%2==0 else 2*radius+1-order
    return nodes,weights,rate


def derivative_consistent(coarse,fine,scale,digits):
    return all(mp.isfinite(v) for v in (coarse,fine,scale)) and abs(fine-coarse)<=mp.power(10,-digits+4)*scale


def translation_independent(value,variable):
    """Recognize an exactly parameter-independent complete functional.

    Translate each finite integral's dummy coordinate by its lower bound.
    This tests independence within one operand; it never substitutes another
    engine's representation or differentiates singular boundary terms.
    """
    replacements={}
    for integral in value.atoms(sp.Integral):
        if len(integral.limits)!=1: return False
        x,lo,hi=integral.limits[0]
        if lo in (sp.oo,-sp.oo) or hi in (sp.oo,-sp.oo): return False
        t=sp.Dummy('translated',real=True)
        replacements[integral]=sp.Integral(integral.function.subs(x,t+lo),(t,0,hi-lo))
    return not value.xreplace(replacements).has(variable)


@lru_cache(maxsize=1024)
def numerical_derivative(value,point_items,digits):
    """Differentiate the complete functional, including its moving limits.

    Never apply Leibniz's rule to singular boundary terms separately. Exact
    rational, symmetric finite-difference stencils evaluate the original
    functional at every node through the same convergence/precision checks.
    Two step sizes and the caller's independent working precisions check the
    estimate. This is numerical evidence, not a certified derivative bound.
    No information from the other operand enters this evaluation.
    """
    points=dict(point_items)
    variable,order=value.variable_count[0]
    order=int(order)
    if variable not in points:
        raise UnresolvedEvaluation('held derivative: evaluation point unavailable for '+math_text(variable))
    center=points[variable]
    if center.is_real is not True or center.is_finite is not True:
        raise UnresolvedEvaluation('held derivative: finite real evaluation point not established for '+math_text(variable))
    if not 1<=order<=4:
        raise UnresolvedEvaluation('held derivative: stencil order not implemented: '+str(order))
    remainder=(sp.Derivative(value.expr,*value.variable_count[1:],evaluate=False)
               if len(value.variable_count)>1 else value.expr)
    if translation_independent(remainder,variable):
        # Still establish that the undifferentiated integral exists at this
        # point before returning the derivative of a constant functional.
        numerical_expression(remainder,digits,points)
        return mp.mpf(0)
    nodes,weights,rate=derivative_stencil(order)
    exponent=max(3,(digits+6+rate-1)//rate)
    precision=digits+order*exponent+12
    estimates=[]; scales=[]
    for divisor in (1,2):
        step=sp.Rational(1,divisor*10**exponent)
        values=[]
        for node in nodes:
            shifted=dict(points); shifted[variable]=center+node*step
            progress_event('derivative node',variable=math_text(variable),order=int(order),
                           digits=digits,step=str(step),node=node)
            try:
                values.append(numerical_expression(remainder,precision,shifted))
            except (UnresolvedEvaluation,UndefinedIntegral,ZeroDivisionError) as exc:
                raise UnresolvedEvaluation('held derivative: complete functional unavailable at '+
                    math_text(variable)+'='+math_text(shifted[variable])+': '+str(exc)) from exc
        with mp.workdps(precision+30):
            h=scalar_number(step,precision+30)
            baseline=values[len(nodes)//2]
            estimates.append(mp.fsum(scalar_number(w,precision+30)*(v-baseline)
                                      for w,v in zip(weights,values))/h**order)
            # An absolute step-consistency check at the functional's scale
            # also permits derivatives that vanish; no residual-size test.
            scales.append(max(abs(v) for v in values))
    with mp.workdps(precision+30):
        scale=max(*scales,*(abs(v) for v in estimates))
        if not derivative_consistent(*estimates,scale,digits):
            raise UnresolvedEvaluation('held derivative: independent step-size consistency not established for '+math_text(variable))
        return +estimates[-1]


def numerical_functional(value,digits,points):
    if isinstance(value,sp.Subs):
        locations=tuple(p.subs(points,simultaneous=True) for p in value.point)
        local=dict(points); local.update(zip(value.variables,locations))
        return numerical_expression(value.expr,digits,local)
    if isinstance(value,sp.Derivative):
        return numerical_derivative(value,tuple(sorted(points.items(),key=lambda p:sp.default_sort_key(p[0]))),digits)
    if isinstance(value,sp.Integral):
        return numerical_integral(value.subs(points,simultaneous=True),digits)
    raise UnresolvedEvaluation('functional evaluator unavailable for '+type(value).__name__)


def numerical_expression(value,digits,points=None):
    points={} if points is None else points
    # Only outer integrals are independent quadrature terms. An inner integral
    # can depend on its outer binder and must remain inside that evaluator.
    integrals=set()
    traversal=sp.preorder_traversal(value)
    for node in traversal:
        if isinstance(node,(sp.Integral,sp.Derivative,sp.Subs)):
            integrals.add(node)
            traversal.skip()
    integrals=sorted(integrals,key=sp.default_sort_key)
    slots=tuple(sp.Dummy('quadrature_'+str(i)) for i in range(len(integrals)))
    expression=value.xreplace(dict(zip(integrals,slots))).subs(points,simultaneous=True)
    if expression.free_symbols-set(slots):
        raise UnresolvedEvaluation('functional evaluation has unassigned names '+','.join(sorted(map(str,expression.free_symbols-set(slots)))))
    function=sp.lambdify(slots,expression,modules='mpmath',cse=True)
    numbers=[numerical_functional(integral,digits,points) for integral in integrals]
    with mp.workdps(digits+15):
        return +function(*numbers)


CURRENT_OBJECT=None
SHOW_PROGRESS=False


def progress_event(stage,**fields):
    if SHOW_PROGRESS:
        print('S9B_COMPARATOR_PROGRESS '+json.dumps(dict(object=CURRENT_OBJECT,stage=stage,monotonic_seconds=time.monotonic(),**fields)),
              file=sys.stderr,flush=True)

def sample_difference(a,b,residual=None):
    if 'undefined' in a or 'undefined' in b:
        return {'coverage':'undefined operand; comparison not formed'}
    if 'unresolved' in a or 'unresolved' in b:
        return {'coverage':'unresolved coverage'}
    if residual is not None:
        if 'undefined' in residual: return {'coverage':'undefined residual; comparison not formed'}
        if 'unresolved' in residual: return {'coverage':'unresolved coverage'}
        return {'coverage':'residual formed','residual':residual}
    if 'truth' in a and 'truth' in b: return {'boolean_difference':a['truth']!=b['truth']}
    if 'exact' in a and 'exact' in b:
        # Re-evaluation uses already-created exact SymPy objects in compare_math;
        # this field carries no untrusted eval of printed text.
        return {'same_exact_representation':a['exact']==b['exact']}
    return {'coverage':'unresolved coverage','reason':'no arithmetic comparison formed'}


def condition_pieces(expr):
    """Disjoint truth strata of each union and its complement, syntactically.

    No exponential DNF expansion is needed: each Or's ith piece is its operand
    with earlier operands false. Nested unions get their own coverage entries.
    """
    seen=set(); pieces=[]
    def add(label,value):
        key=sp.srepr(value,order='none')
        if key not in seen: seen.add(key); pieces.append((label,value))
    add('condition',expr); add('complement',sp.Not(expr))
    for sub in sp.preorder_traversal(expr):
        if isinstance(sub,sp.Or):
            prior=[]
            for i,arg in enumerate(sub.args):
                add('union piece '+str(i),sp.And(arg,*(sp.Not(x) for x in prior)))
                prior.append(arg)
            add('union complement',sp.And(*(sp.Not(x) for x in prior)))
        if isinstance(sub,Relational):
            add('relation true',sub); add('relation false',sp.Not(sub))
            if isinstance(sub.lhs,sp.Expr) and isinstance(sub.rhs,sp.Expr):
                add('relation boundary',sp.Eq(sub.lhs,sub.rhs))
    return pieces


def compare_math(left,right,role):
    boolean = isinstance(left,Boolean) and not isinstance(left,sp.Expr)
    boolean_r = isinstance(right,Boolean) and not isinstance(right,sp.Expr)
    if boolean != boolean_r:
        return {'outcome':'not_formed','reason':'boolean/scalar representation differs'},False
    if isinstance(left,sp.Set) or isinstance(right,sp.Set):
        return {'outcome':'not_formed','reason':'set representation retained; no set equivalence implemented'},False
    if boolean:
        residual=sp.Xor(left,right)
        result={'outcome':'predicate','py_condition':math_text(left),'wl_condition':math_text(right),
                'boolean_residual':math_text(residual)}
    elif isinstance(left,sp.Expr) and isinstance(right,sp.Expr):
        residual=symbolic_difference(left,right)
        result={'outcome':'symbolic','py_value':math_text(left),'wl_value':math_text(right),'residual':math_text(residual),
                'zero_class': 'zero' if residual.is_zero is True else 'nonzero' if residual.is_zero is False else 'unresolved'}
    else: return {'outcome':'not_formed','reason':'unsupported value type'},False
    # Functionals and predicates get numerical evidence even when their formal
    # subtraction simplifies. Live profiles are always identical on both sides.
    if not boolean and not (left.has(AppliedUndef,sp.Integral,sp.Derivative,sp.Subs,sp.Piecewise) or right.has(AppliedUndef,sp.Integral,sp.Derivative,sp.Subs,sp.Piecewise)):
        result['samples']=[]
        result['sample_policy']='exact scalar subtraction'
        return result,True
    samples=[]; covered=set(); pieces=[]
    if boolean:
        pieces=[(eng,label,p) for eng,e in (('py',left),('wl',right)) for label,p in condition_pieces(e)]
    for label,point,profile in configurations((left,right)):
        progress_event('sample',configuration=label,role=role)
        lv=numeric_value(left,point,profile); rv=numeric_value(right,point,profile)
        sampled_residual=None if boolean else numerical_residual(lv,rv)
        record={'configuration':label,'point':{k:str(v) for k,v in point.items()},'profile':profile,
                'profile_formulas':{k:str(v) for k,v in profile_formulas(profile,point)[1].items()},
                'py':lv,'wl':rv,'comparison':sample_difference(lv,rv,sampled_residual)}
        if not boolean:
            record['residual']=sampled_residual
        samples.append(record)
        for i,(_,_,piece) in enumerate(pieces):
            if i not in covered and numeric_value(piece,point,profile).get('truth') is True: covered.add(i)
    result['samples']=samples
    if pieces:
        result['strata']=[{'engine':eng,'piece':label,'condition':math_text(p),
                           'coverage':stratum_coverage(p,i in covered)} for i,(eng,label,p) in enumerate(pieces)]
    formed=any(sample_formed(row['comparison']) for row in samples)
    result['comparison_formed']=formed
    return result,formed


def sample_formed(finding):
    return (finding.get('coverage')=='residual formed' or 'boolean_difference' in finding
            or 'same_exact_representation' in finding)


def stratum_coverage(piece,sampled):
    return 'empty' if piece is sp.false else 'sampled' if sampled else 'unresolved coverage'


def symbolic_difference(left,right):
    residual=left-right
    if not residual.has(sp.Integral,sp.Derivative,sp.Subs,sp.Piecewise): residual=sp.cancel(residual)
    return residual


def text_comparison(left,right,role):
    lv=LABEL_TABLE.get(left,left) if role=='label' else left
    return lv==right


COMPARISON_KINDS=('value','domain','marker','label')


def comparison_kind(result,role):
    if result.get('outcome') in ('marker','label'): return result['outcome']
    return 'label' if role=='label' else 'domain' if role=='domain' else 'value'


def compare_domains(left,right,bindings):
    translated=[]; errors=[]
    for node,engine in ((left,'py'),(right,'wl')):
        if node is None:
            translated.append(None); errors.append('no bound counterpart'); continue
        try:
            translated.append(math_node(node,bindings,engine)); errors.append(None)
        except Unsupported as exc:
            translated.append(None); errors.append(str(exc))
    if not any(errors) and all(isinstance(e,Boolean) and not isinstance(e,sp.Expr) for e in translated):
        return compare_math(*translated,'domain')
    samples=[]
    for label,point,profile in configurations(tuple(e for e in translated if e is not None)):
        values=[{'unresolved':'domain translation: '+error} if error else numeric_value(expr,point,profile)
                for expr,error in zip(translated,errors)]
        samples.append({'configuration':label,'point':{k:str(v) for k,v in point.items()},'profile':profile,
            'profile_formulas':{k:str(v) for k,v in profile_formulas(profile,point)[1].items()},
            'py':values[0],'wl':values[1],'comparison':sample_difference(*values)})
    return {'outcome':'predicate','samples':samples,'strata':[
        {'engine':eng,'piece':'condition and complement','coverage':'unresolved coverage','reason':error}
        for eng,error in zip(('py','wl'),errors) if error]},False
def structural_delta(a,b,bindings):
    """Serialization diagnostic only; never a mathematical zero."""
    def signature(n,engine):
        if n.head in ('Symbol','Dummy','FunctionName'):
            value=bindings.name(n.value,engine)
        else: value=n.value
        return (n.head,value,n.options,tuple(signature(x,engine) for x in n.args))
    sa,sb=signature(a,'py'),signature(b,'wl')
    return {'same_bound_serialization':sa==sb,
            'py_tree_digest':hashlib.sha256(repr(sa).encode()).hexdigest(),
            'wl_tree_digest':hashlib.sha256(repr(sb).encode()).hexdigest(),
            'scope':'representation only; not a mathematical residual'}


class Comparator:
    def __init__(self,names):
        self.bindings=Bindings(names)
        self.cache=OrderedDict()
        self.attempted=[set(),set()]
        self.formed_leaves={kind:[set(),set()] for kind in COMPARISON_KINDS}

    def leaf(self,a,b,role):
        left,right=a.node,b.node
        if role=='domain':
            key=(left,right,role)
            if key not in self.cache:
                self.cache[key]=compare_domains(left,right,self.bindings)
                if len(self.cache)>256: self.cache.popitem(last=False)
            return self.cache[key]
        if left.head=='Text' or right.head=='Text':
            if left.head!=right.head:
                return {'outcome':'not_formed','reason':'marker/text versus other value'},False
            comparison=text_comparison(left.value,right.value,role)
            return {'outcome':'label' if role=='label' else 'marker',
                    'comparison':comparison,'py_text':left.value,'wl_text':right.value},comparison is not None
        try:
            lv=math_node(left,self.bindings,'py'); rv=math_node(right,self.bindings,'wl')
        except Unsupported as exc:
            return {'outcome':'not_formed','reason':str(exc),
                    'representation':structural_delta(left,right,self.bindings)},False
        key=(lv,rv,role)
        if key not in self.cache:
            self.cache[key]=compare_math(lv,rv,role)
            if len(self.cache)>256: self.cache.popitem(last=False)
        self.cache.move_to_end(key)
        return self.cache[key]

    def compare(self,a,b):
        role='domain' if a.role=='domain' or b.role=='domain' else a.role
        # Source wrappers carry extra metadata, not changes to object identity.
        # Extract the declared Object/value, and account for every unmatched
        # premise/applicability field separately; never conjoin unlike gates.
        for side,view,other in (('py',a,b),('wl',b,a)):
            if view.children and 'value' in view.children and ('applicability' in view.children or 'premise' in view.children):
                result,counts=self.compare(view.children['value'],other) if side=='py' else self.compare(other,view.children['value'])
                extras={k:self.absent(v,side) for k,v in view.children.items() if k!='value'}
                return {'outcome':'wrapper','comparison':result,'unpaired_metadata':extras},counts
        if a.children is not None or b.children is not None:
            ca,cb=a.children or {},b.children or {}
            if a.children is None or b.children is None:
                return {'outcome':'not_formed','reason':'container/scalar representation differs',
                        'representation':structural_delta(a.ref.node,b.ref.node,self.bindings)},(set(),set())
            result={}; used=[set(),set()]
            for key in sorted(ca.keys()|cb.keys()):
                if key not in ca or key not in cb:
                    result[key]=self.absent(ca.get(key,cb.get(key)),'py' if key in ca else 'wl')
                else:
                    result[key],counts=self.compare(ca[key],cb[key])
                    used[0].update(counts[0]); used[1].update(counts[1])
            if not ca and not cb:
                # Empty containers are values whose constructor is observable.
                used=(leaf_paths(a.ref),leaf_paths(b.ref))
                kind=comparison_kind({},role)
                for i in range(2): self.formed_leaves[kind][i].update(used[i])
                return {'outcome':'empty_container','comparison':a.ref.node.head==b.ref.node.head,
                        'comparison_kind':kind,'comparison_formed':True},used
            return {'outcome':'container','role':role,'children':result},tuple(used)
        for i,ref in enumerate((a.ref,b.ref)):
            self.attempted[i].update(leaf_paths(ref,False))
        result,consumed=self.leaf(a.ref,b.ref,role)
        kind=comparison_kind(result,role)
        result=dict(result,role=role,py_path=a.ref.path,wl_path=b.ref.path,
                    py_operand=data(a.ref.node),wl_operand=data(b.ref.node),
                    representation=structural_delta(a.ref.node,b.ref.node,self.bindings),
                    comparison_kind=kind,comparison_formed=bool(consumed))
        used=(leaf_paths(a.ref,False),leaf_paths(b.ref,False)) if consumed else (set(),set())
        for i in range(2): self.formed_leaves[kind][i].update(used[i])
        return result,used

    def absent(self,view,engine):
        result={'outcome':'not_formed','reason':'nested sibling has no bound counterpart','engine':engine,
                'path':view.ref.path,'role':view.role,'parsed_leaves':census(view.ref.node)}
        if view.children:
            result['children']={k:self.absent(v,engine) for k,v in view.children.items()}
        elif view.role=='domain':
            self.attempted[0 if engine=='py' else 1].update(leaf_paths(view.ref,False))
            # A missing counterpart must not suppress truth sampling of the
            # attached condition that is present on the other engine.
            result['domain_comparison']=compare_domains(view.ref.node if engine=='py' else None,
                view.ref.node if engine=='wl' else None,self.bindings)[0]
        return result


def coverage_inventory(root,engine):
    """Enumerate emitted branches independently of paired extraction.

    Identical repetitions are occurrence-listed, not silently dropped. Every
    condition is printed with a complement; sample coverage is only asserted
    by the paired comparison's strata entries, never by a raw parse match.
    """
    result={}
    def add(ref,role):
        payload=data(ref.node); encoded=json.dumps(payload,separators=(',',':'))
        key=hashlib.sha256(encoded.encode()).hexdigest()
        if key not in result:
            result[key]={'role':role,'condition':payload,'complement':['Not',payload],
                         'coverage':'unresolved coverage unless sampled in paired comparison', 'paths':[]}
        result[key]['paths'].append(ref.path)
    def visit(ref):
        n=ref.node; h,a=n.head,n.args; wh,wa=wl_head(n)
        if attached(ref): add(ref.child(1),'attached domain')
        if h=='Record':
            for key,field in ref.fields().items():
                if key in DOMAIN_KEYS: add(field,'attached domain '+key)
        if wh=='ConditionalExpression' and len(a)==3: add(ref.child(2),'conditional domain')
        if h=='Or' or wh=='Or':
            start=1 if wh else 0
            for i in range(start,len(a)): add(ref.child(i),'union piece (overlap retained)')
            add(ref,'union and complement')
        if h=='Piecewise':
            for i,x in enumerate(a):
                if len(x.args)==2: add(ref.child(i).child(1),'piecewise stratum')
        if wh=='Piecewise' and len(a)>1 and a[1].head=='Sequence':
            for i,x in enumerate(a[1].args):
                if len(x.args)==2: add(ref.child(1).child(i).child(1),'piecewise stratum')
        # Relational boundaries expose zero-parameter/profile/density strata,
        # including inside unsupported quantified/set constructions.
        if h in RELATIONS or wh in ('Equal','Unequal','Greater','Less','GreaterEqual','LessEqual'):
            add(ref,'predicate and boundary audit')
        for i in range(len(a)): visit(ref.child(i))
    visit(Ref(root))
    return list(result.values())


def emit(out,value):
    out.write(json.dumps(value,separators=(',',':'),ensure_ascii=False)+'\n'); out.flush()


def run(py_path,wl_path,out,accounting,*,names=NAME_TABLE,joins=SHARED,progress=True):
    global SHOW_PROGRESS,CURRENT_OBJECT
    SHOW_PROGRESS=progress
    started=time.monotonic(); validate_tables(names,joins); function_bindings()
    streams={'py':Stream(py_path,'py'),'wl':Stream(wl_path,'wl')}
    comparator=Comparator(names)
    emit(out,{'kind':'configuration','name_table':[asdict(x) for x in names],
              'function_table':[asdict(x) for x in FUNCTION_TABLE],
              'sources':{'py':PY_SOURCE,'wl':WL_SOURCE,'spec':SPEC},
              'input_sha256':{e:s.digest.hexdigest() for e,s in streams.items()},
              'profiles':[{'name':p,'admissibility':PROFILE_CITATIONS,
                 'qualification':'boundary probe, outside propagating subset where its condition is false' if p in ('brane_zero_at_point','stiffness_zero_at_point') else 'radial test profile; optical tails within supplied counting',
                 'formulas':{k:str(v) for k,v in profile_formulas(p,POINTS[0])[1].items()}} for p in PROFILE_NAMES],
              'points':[{k:str(v) for k,v in p.items()} for p in POINTS],
              'numeric_scope':'exact substitutions; absolute convergence screening and tanh-sinh error, precision and panel checks; complete-functional derivatives checked at two rational step sizes and working precisions; 12 printed digits only after reliability checks; numerical evidence not a certified bound; unresolved strata are coverage findings',
              'target':'no residual target was supplied'})
    summary=Counter(rows_without_formed_comparison=0,rows_without_formed_value_comparison=0)
    sample_outcomes=Counter(); sample_outcomes_by_kind={kind:Counter() for kind in COMPARISON_KINDS}
    rows_without_value=[]
    allkeys=sorted(set(joins)|set(streams['py'].index)|set(streams['wl'].index))
    for key in allkeys:
        CURRENT_OBJECT=key
        row_started=time.monotonic()
        if progress: print('S9B_COMPARATOR_WORK '+key,file=sys.stderr,flush=True)
        occurrences={e:len(s.index.get(key,())) for e,s in streams.items()}
        if key.startswith('LOCAL/'):
            for e,s in streams.items():
                for i in range(occurrences[e]):
                    n=s.read(key,i) # Local input grammar is validated too.
                    row={'name':key,'engine':e,'accounting':'local','reason':'engine-local tag; not joined',
                         'line':s.index[key][i][1],'occurrence':i,'parsed_leaves':census(n)}
                    emit(accounting,row); emit(out,dict(kind='accounting',**row)); summary['local_'+e]+=1
            continue
        reason = ('duplicate shared-name occurrence; injective join rejected' if any(x>1 for x in occurrences.values()) else
                  'outside shared vocabulary' if key not in joins else
                  'missing from '+','.join(e for e,c in occurrences.items() if c==0) if not all(occurrences.values()) else None)
        if reason:
            row={'name':key,'accounting':'unjoined','reason':reason,'occurrences':occurrences}
            for e,s in streams.items():
                row[e+'_parsed_leaves']=[census(s.read(key,i)) for i in range(occurrences[e])]
            emit(accounting,row); emit(out,dict(kind='accounting',**row)); summary['unjoined']+=1
            continue
        progress_event('parse')
        nodes={e:s.read(key) for e,s in streams.items()}
        parsed={e:census(n) for e,n in nodes.items()}
        # Raw operands are emitted before extraction/comparison and any outcome.
        emit(out,{'kind':'operands','name':key,**{e:data(n) for e,n in nodes.items()}})
        comparator.attempted=[set(),set()]
        comparator.formed_leaves={kind:[set(),set()] for kind in COMPARISON_KINDS}
        try:
            views={e:extract(Ref(n),e) for e,n in nodes.items()}
            progress_event('compare')
            result,used=comparator.compare(views['py'],views['wl'])
        except Unsupported as exc:
            result={'outcome':'not_formed','reason':'extraction '+str(exc)}; used=(set(),set())
        emit(out,{'kind':'comparison','name':key,'comparison':result})
        def count_samples(value,kind='value'):
            if isinstance(value,dict):
                kind=value.get('comparison_kind', 'domain' if value.get('role')=='domain' else kind)
                for sample in value.get('samples',[]):
                    finding=sample['comparison']
                    outcome=finding.get('coverage','boolean comparison formed' if 'boolean_difference' in finding else 'exact representation comparison formed')
                    sample_outcomes[outcome]+=1
                    sample_outcomes_by_kind[kind][outcome]+=1
                for k,v in value.items():
                    if k!='samples': count_samples(v,'domain' if k=='domain_comparison' else kind)
            elif isinstance(value,list):
                for v in value: count_samples(v,kind)
        count_samples(result)
        # Census is computed independently from parse, not from extracted views.
        row={'name':key,'accounting':'joined',
             **{e:{'parsed_leaves':parsed[e],'compared_leaves':len(used[i]),
                    'compared_leaves_by_kind':{kind:len(comparator.formed_leaves[kind][i]) for kind in COMPARISON_KINDS},
                    'uncompared_leaves':parsed[e]-len(used[i]),
                    'attempted_leaves':len(comparator.attempted[i] | used[i]),
                    'attempted_unformed_leaves':len(comparator.attempted[i]-used[i]),
                    'unattempted_leaves':parsed[e]-len(comparator.attempted[i] | used[i])} for i,e in enumerate(('py','wl'))}}
        row['comparison_kinds_formed']={kind:any(comparator.formed_leaves[kind]) for kind in COMPARISON_KINDS}
        row['comparison_formed']=any(row['comparison_kinds_formed'].values())
        row['value_comparison_formed']=row['comparison_kinds_formed']['value']
        if not row['value_comparison_formed']:
            summary['rows_without_formed_value_comparison']+=1
            rows_without_value.append(key)
        if not row['comparison_formed']: summary['rows_without_formed_comparison']+=1
        for kind,formed in row['comparison_kinds_formed'].items():
            if formed: summary['rows_with_formed_'+kind+'_comparison']+=1
        if any(row[e]['attempted_leaves']>parsed[e] for e in nodes):
            raise RuntimeError('comparison coverage exceeds independent parse census')
        emit(accounting,row); emit(out,dict(kind='accounting',**row)); summary['joined']+=1
        if any(row[e]['uncompared_leaves'] for e in nodes): summary['rows_with_uncompared_leaves']+=1
        for e,n in nodes.items():
            emit(out,{'kind':'stratum_inventory','name':key,'engine':e,'strata':coverage_inventory(n,e)})
        # Bounded memoization reuses only identical mathematical inputs. Every
        # row, configuration and coverage record is still emitted in full.
        progress_event('completed',seconds=time.monotonic()-row_started,
                       samples_cache=_numeric_value.cache_info()._asdict(),
                       integrals_cache=numerical_integral.cache_info()._asdict())
    summary.update({'wall_seconds':time.monotonic()-started,
                    'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    summary['sample_outcomes']=dict(sample_outcomes)
    summary['sample_outcomes_by_kind']={kind:dict(counts) for kind,counts in sample_outcomes_by_kind.items()}
    summary['rows_without_formed_value_comparison_names']=rows_without_value
    emit(out,{'kind':'completion',**summary})
    return dict(summary)


def main():
    faulthandler.register(signal.SIGUSR1,all_threads=True)
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--py',type=Path,default=DEFAULT_PY)
    parser.add_argument('--wl',type=Path,default=DEFAULT_WL)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--accounting',type=Path,required=True)
    args=parser.parse_args()
    for path in (args.output,args.accounting):
        if 'out' in path.parts: parser.error('write outputs in a private scratch directory')
    try:
        with args.output.open('x') as output,args.accounting.open('x') as accounting:
            result=run(args.py,args.wl,output,accounting)
        print(json.dumps(result),flush=True)
        return 0
    except InputError as exc:
        print('operational input error: '+str(exc),file=sys.stderr,flush=True)
        return 2


if __name__=='__main__':
    sys.exit(main())
