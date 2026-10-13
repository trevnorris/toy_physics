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
    ('s9b_velocity_integration_constant', 'velocityIntegrationConstant',
     'velocity reduction constant; repair item 6', 470, 524, 'parameter'),
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
WL_BUILTINS = {'Plus':'Add','Times':'Mul','Power':'Pow','Equal':'Equality',
               'Unequal':'Unequality','Greater':'StrictGreaterThan','Less':'StrictLessThan',
               'GreaterEqual':'GreaterThan','LessEqual':'LessThan', 'List':'Sequence',
               'Integrate':'Integral', 'Inactive/Integrate':'Integral',
               'Inactive/Limit':'Limit','Element':'Contains'}


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
        if wh in WL_BUILTINS or wh in ('And','Or','Not','Piecewise','ConditionalExpression','Rational','Min','Max','Abs','Sin','Cos','Tan','Sec','ArcCos','Log','Exp','Sqrt','Re','Im','densityPower'):
            h,a = WL_BUILTINS.get(wh,wh),wa
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
    if h in ('Pow','densityPower'):
        if len(a)!=2: raise Unsupported('power arity')
        base,exponent = cv(a[0]),cv(a[1])
        # Both native definitions make zero powers one, including at zero.
        return sp.Pow(base,exponent)
    if h == 'Rational': return sp.Rational(*(cv(x) for x in a))
    if h in RELATIONS:
        if len(a)!=2 or any(x.head in ('Equality','StrictGreaterThan','StrictLessThan','LessThan','GreaterThan') for x in a):
            raise Unsupported('chained or nonbinary relation')
        return RELATIONS[h](*(cv(x) for x in a),evaluate=False)
    if h in ('And','Or','Not'):
        values = [cv(x) for x in a]
        if not all(isinstance(x,Boolean) for x in values): raise Unsupported('nonboolean connective')
        return getattr(sp,h)(*values,evaluate=False)
    if h == 'Contains':
        left,right = cv(a[0]),cv(a[1])
        if isinstance(left,sp.Tuple): return sp.And(*(sp.Contains(x,right,evaluate=False) for x in left))
        return sp.Contains(left,right,evaluate=False)
    if h in ('Sequence','Alternatives','FiniteSet'): return (sp.Tuple if h in ('Sequence','Alternatives') else sp.FiniteSet)(*(cv(x) for x in a))
    if h == 'Interval': return sp.Interval(*(cv(x) for x in a))
    if h in ('sin','cos','tan','log','exp','sqrt','Abs','Min','Max','re','im','acos','Sin','Cos','Tan','Log','Exp','Sqrt','Re','Im','ArcCos','Sec'):
        function = {'Sin':sp.sin,'Cos':sp.cos,'Tan':sp.tan,'Log':sp.log,'Exp':sp.exp,'Sqrt':sp.sqrt,'Re':sp.re,'Im':sp.im,'ArcCos':sp.acos,'Sec':sp.sec}.get(h,getattr(sp,h,None))
        return function(*(cv(x) for x in a))
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
    return ref.node.head=='Sequence' and len(ref.node.args)==2 and is_text(ref.node.args[0],'ATTACHED_DOMAIN')


def extract(ref, engine, role='value'):
    """Source-defined wrappers only; no algebraic representation reconciliation.

    PY 76-93,605-611; WL 270-272,618-623,661 define these roles.
    All wrapper metadata remains in raw operands and parsed leaf counts. Only
    compared descendants contribute to comparison coverage.
    """
    n=ref.node; fields=ref.fields()
    if attached(ref): return extract(ref.child(1),engine,'domain')
    if engine=='py' and n.head=='Sequence':
        a=n.args
        if len(a)==3 and is_text(a[0],'APPLICABLE_OBJECT') and attached(ref.child(1)):
            return View(ref,role,{'value':extract(ref.child(2),engine),
                                  'applicability':extract(ref.child(1),engine)})
        if len(a)==2 and all(x.head=='Sequence' and len(x.args)==2 and attached(ref.child(i).child(0)) for i,x in enumerate(a)):
            return View(ref,role,{name: View(ref.child(i), 'case', {
                'domain':extract(ref.child(i).child(0),engine),
                'value':extract(ref.child(i).child(1),engine)}) for i,name in enumerate(('on','outside'))})
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
    dict(s9b_r='11/3',b='7/3',Z_E='101/3',Z_R='137/3',c_0='7/3',GM='1/101',
         Phi='2/103',K='5/7',n='2',m='3/2',s='1/2',rho_0='11/7',s9b_gamma='2/3',s9b_b_min='1/7'),
    dict(s9b_r='17/5',b='11/5',Z_E='151/5',Z_R='173/5',c_0='11/5',GM='2/107',
         Phi='3/109',K='7/11',n='3',m='5/3',s='2',rho_0='13/11',s9b_gamma='3/5',s9b_b_min='1/11'),
    dict(s9b_r='23/7',b='13/7',Z_E='211/7',Z_R='263/7',c_0='13/7',GM='3/113',
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
    for name in ('GM','Phi','K','n','s','s9b_velocity_integration_constant','rho_bulk','s9b_fractional_bulk_change'):
        if sp.Symbol(name) in free:
            for val in (sp.S.Zero,sp.S.One,-sp.S.One,sp.Rational(1,2)):
                yield name+'='+str(val),dict(p,**{name:val}),'inverse'


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
            return value if value.has(sp.Integral) else value.doit()
        if isinstance(node,sp.Subs):
            return sp.Subs(*args).doit(deep=False)
        return node.func(*args) if args!=node.args else node
    value=prepare(expr)
    prepare.cache_clear()
    substitutions={sp.Symbol(k):v for k,v in point.items()}
    substitutions[sp.Symbol('r_E')]=sp.sqrt(point['b']**2+point['Z_E']**2)
    substitutions[sp.Symbol('r_R')]=sp.sqrt(point['b']**2+point['Z_R']**2)
    return value.subs(substitutions, simultaneous=True)


def math_text(value):
    """Print stored mathematical content without evaluating it to order terms.

    SymPy's default Add printer can call complex() on numerical integrals and
    run quadrature while formatting an undefined sample. Stored argument order
    retains the full expression without that accidental second computation.
    """
    return sp.sstr(value,order='none')


def numeric_value(expr,point,profile):
    return _numeric_value(expr,tuple(sorted(point.items())),profile)


@lru_cache(maxsize=4096)
def _numeric_value(expr,point_items,profile):
    value=substitute_profiles(expr,dict(point_items),profile)
    if value.free_symbols:
        return {'unresolved':'free names '+','.join(sorted(map(str,value.free_symbols)))}
    if value.has(sp.Derivative,sp.Limit): return {'unresolved':'unevaluated derivative or limit'}
    if value.has(sp.nan,sp.zoo): return {'undefined':math_text(value)}
    if value in (sp.true,sp.false): return {'truth':bool(value)}
    if isinstance(value,Boolean) and not isinstance(value,sp.Expr):
        reduced=sp.simplify(value)
        if reduced in (sp.true,sp.false): return {'truth':bool(reduced)}
        return {'unresolved':'predicate '+math_text(reduced)}
    if isinstance(value,sp.Set): return {'exact':math_text(value)}
    # Independent precision evaluations retain their discrepancy; no tolerance
    # becomes an equality claim. evalf's precision is not a certified bound.
    if value.has(sp.Integral):
        for integral in value.atoms(sp.Integral):
            if len(integral.limits)!=1: return {'unresolved':'multidimensional integral'}
        try:
            low=numerical_expression(value,18); high=numerical_expression(value,30)
        except (ValueError,NotImplementedError) as exc:
            return {'unresolved':'quadrature '+type(exc).__name__+': '+str(exc)}
        if not mp.isfinite(high): return {'undefined':str(high)}
        with mp.workdps(40):
            difference=mp.nstr(high-low,30)
        return {'approx_18':mp.nstr(low,18),'approx_30':mp.nstr(high,30),
                'precision_difference':difference,
                'method':'compiled mpmath tanh-sinh; independent integral evaluation; fixed working precision',
                'error_bound':'not certified; numerical evidence only'}
    if not isinstance(value,sp.Expr): return {'unresolved':'nonscalar sample '+math_text(value)}
    return {'exact':math_text(value),'approx_30':math_text(value.evalf(30))}



@lru_cache(maxsize=2048)
def numerical_integral(integral,digits):
    """Evaluate the original integral, with no cross-engine rewriting.

    Compile once per distinct substituted integral/precision. SymPy evalf walks
    the symbolic integrand at every node and repeatedly reintegrates a difference
    to chase relative accuracy lost in cancellation. Instead evaluate each
    integral independently, retain both precisions, and do arithmetic on those
    numbers. This is numerical evidence, never a certified zero or error bound.
    """
    (variable,lower,upper),=integral.limits
    if integral.function.has(sp.Integral):
        # Retain the original evaluator's nested-integral support. The fast
        # path applies to single integrals; it never removes a supported case.
        evaluated=integral.evalf(digits)
        if evaluated.has(sp.Integral):
            raise NotImplementedError('quadrature not formed')
        real,imag=evaluated.as_real_imag()
        with mp.workdps(digits+15):
            return mp.mpc(str(real),str(imag))
    function=sp.lambdify(variable,integral.function,modules='mpmath',cse=True)
    bounds=sp.lambdify((),(lower,upper),modules='mpmath',cse=True)
    with mp.workdps(digits+15):
        lo,hi=bounds()
        return +mp.quad(function,[lo,hi],method='tanh-sinh')


def numerical_expression(value,digits):
    # Only outer integrals are independent quadrature terms. An inner integral
    # can depend on its outer binder and must remain inside that evaluator.
    integrals=set()
    traversal=sp.preorder_traversal(value)
    for node in traversal:
        if isinstance(node,sp.Integral):
            integrals.add(node)
            traversal.skip()
    integrals=sorted(integrals,key=sp.default_sort_key)
    slots=tuple(sp.Dummy('quadrature_'+str(i)) for i in range(len(integrals)))
    expression=value.xreplace(dict(zip(integrals,slots)))
    function=sp.lambdify(slots,expression,modules='mpmath',cse=True)
    numbers=[numerical_integral(integral,digits) for integral in integrals]
    with mp.workdps(digits+15):
        return +function(*numbers)


CURRENT_OBJECT=None
SHOW_PROGRESS=False


def progress_event(stage,**fields):
    if SHOW_PROGRESS:
        print('S9B_COMPARATOR_PROGRESS '+json.dumps(dict(object=CURRENT_OBJECT,stage=stage,monotonic_seconds=time.monotonic(),**fields)),
              file=sys.stderr,flush=True)

def sample_difference(a,b):
    if 'truth' in a and 'truth' in b: return {'boolean_difference':a['truth']!=b['truth']}
    if 'exact' in a and 'exact' in b:
        # Re-evaluation uses already-created exact SymPy objects in compare_math;
        # this field carries no untrusted eval of printed text.
        return {'same_exact_representation':a['exact']==b['exact']}
    return {'coverage':'unresolved' if 'unresolved' in a or 'unresolved' in b else 'numeric operands printed'}


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
        residual=left-right
        if not residual.has(sp.Integral,sp.Derivative,sp.Subs,sp.Piecewise): residual=sp.cancel(residual)
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
        record={'configuration':label,'point':{k:str(v) for k,v in point.items()},'profile':profile,
                'profile_formulas':{k:str(v) for k,v in profile_formulas(profile,point)[1].items()},
                'py':lv,'wl':rv,'comparison':sample_difference(lv,rv)}
        if not boolean:
            record['residual']=numeric_value(residual,point,profile)
        samples.append(record)
        for i,(_,_,piece) in enumerate(pieces):
            if i not in covered and numeric_value(piece,point,profile).get('truth') is True: covered.add(i)
    result['samples']=samples
    if pieces:
        result['strata']=[{'engine':eng,'piece':label,'condition':math_text(p),
                           'coverage':'sampled' if i in covered else 'unresolved coverage'} for i,(eng,label,p) in enumerate(pieces)]
    return result,True
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

    def leaf(self,a,b,role):
        left,right=a.node,b.node
        if left.head=='Text' or right.head=='Text':
            if left.head!=right.head:
                return {'outcome':'not_formed','reason':'marker/text versus other value'},False
            lv=LABEL_TABLE.get(left.value,left.value) if role=='label' else left.value
            return {'outcome':'label' if role=='label' else 'marker',
                    'comparison':lv==right.value,'py_text':left.value,'wl_text':right.value},True
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
                return {'outcome':'empty_container','comparison':a.ref.node.head==b.ref.node.head},(leaf_paths(a.ref),leaf_paths(b.ref))
            return {'outcome':'container','role':role,'children':result},tuple(used)
        result,consumed=self.leaf(a.ref,b.ref,role)
        result=dict(result,role=role,py_path=a.ref.path,wl_path=b.ref.path,
                    py_operand=data(a.ref.node),wl_operand=data(b.ref.node),
                    representation=structural_delta(a.ref.node,b.ref.node,self.bindings))
        used=(leaf_paths(a.ref,False),leaf_paths(b.ref,False)) if consumed else (set(),set())
        return result,used

    @staticmethod
    def absent(view,engine):
        return {'outcome':'not_formed','reason':'nested sibling has no bound counterpart','engine':engine,
                'path':view.ref.path,'role':view.role,'parsed_leaves':census(view.ref.node)}


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
    started=time.monotonic(); validate_tables(names,joins)
    streams={'py':Stream(py_path,'py'),'wl':Stream(wl_path,'wl')}
    comparator=Comparator(names)
    emit(out,{'kind':'configuration','name_table':[asdict(x) for x in names],
              'sources':{'py':PY_SOURCE,'wl':WL_SOURCE,'spec':SPEC},
              'input_sha256':{e:s.digest.hexdigest() for e,s in streams.items()},
              'profiles':[{'name':p,'admissibility':PROFILE_CITATIONS,
                 'qualification':'boundary probe, outside propagating subset where its condition is false' if p in ('brane_zero_at_point','stiffness_zero_at_point') else 'radial test profile; optical tails within supplied counting',
                 'formulas':{k:str(v) for k,v in profile_formulas(p,POINTS[0])[1].items()}} for p in PROFILE_NAMES],
              'points':[{k:str(v) for k,v in p.items()} for p in POINTS],
              'numeric_scope':'exact substitutions plus independent 18/30 digit compiled tanh-sinh quadrature with 15 guard digits; no certified integration error bound; residual arithmetic does not chase cancellation; unresolved strata are coverage findings',
              'target':'no residual target was supplied'})
    summary=Counter()
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
        try:
            views={e:extract(Ref(n),e) for e,n in nodes.items()}
            progress_event('compare')
            result,used=comparator.compare(views['py'],views['wl'])
        except Unsupported as exc:
            result={'outcome':'not_formed','reason':'extraction '+str(exc)}; used=(set(),set())
        emit(out,{'kind':'comparison','name':key,'comparison':result})
        # Census is computed independently from parse, not from extracted views.
        row={'name':key,'accounting':'joined',
             **{e:{'parsed_leaves':parsed[e],'compared_leaves':len(used[i]),
                    'uncompared_leaves':parsed[e]-len(used[i])} for i,e in enumerate(('py','wl'))}}
        if any(row[e]['compared_leaves']>parsed[e] for e in nodes):
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
