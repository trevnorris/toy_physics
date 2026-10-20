#!/usr/bin/env python3
"""S12 cross-engine measurement instrument; no residual target.

Native sstr and InputForm readers, explicit semantic bindings, distribution
parts, exact residuals before supplied-domain restriction, and rational probes.
Only the InputForm tokenizer/Pratt reader and parse-census design are adapted
from the accepted S9b comparator; no source engine is imported or executed.
"""
from __future__ import annotations
import argparse
import ast
from collections import Counter, defaultdict
from dataclasses import dataclass, asdict
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import re
import resource
import sys
import time
import sympy as s
from sympy.core.function import AppliedUndef
from sympy.core.relational import Relational
from sympy.logic.boolalg import Boolean

class InputError(ValueError):
    """Operational input/configuration error."""

class Unsupported(Exception):
    """Parsed content without an implemented value comparison."""

@dataclass(frozen=True)
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


def py_parse(source):
    """Lossless sstr grammar, not eval/sympify; native syntax stays in the AST."""
    def walk(n):
        if isinstance(n, ast.Constant):
            return atom('Text' if isinstance(n.value, str) else
                        'Boolean' if isinstance(n.value, bool) else 'Number', n.value)
        if isinstance(n, ast.Name):
            return atom('Symbol', n.id)
        if isinstance(n, (ast.Tuple, ast.List, ast.Set)):
            return Node('FiniteSet' if isinstance(n, ast.Set) else 'Sequence',
                        tuple(walk(x) for x in n.elts))
        if isinstance(n, ast.Dict):
            return Node('Dictionary', tuple(Node('Pair', (walk(k), walk(v)))
                                           for k, v in zip(n.keys, n.values)))
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, (ast.USub, ast.UAdd)):
            v = walk(n.operand)
            return Node('Mul', (atom('Number', -1), v)) if isinstance(n.op, ast.USub) else v
        if isinstance(n, ast.BinOp):
            a, b = walk(n.left), walk(n.right)
            if isinstance(n.op, ast.Sub):
                return Node('Add', (a, Node('Mul', (atom('Number', -1), b))))
            if isinstance(n.op, ast.Div):
                return Node('Mul', (a, Node('Pow', (b, atom('Number', -1)))))
            heads = {ast.Add:'Add', ast.Mult:'Mul', ast.Pow:'Pow', ast.BitAnd:'And', ast.BitOr:'Or'}
            if type(n.op) in heads:
                return Node(heads[type(n.op)], (a, b))
        if isinstance(n, ast.Compare):
            heads = {ast.Lt:'StrictLessThan', ast.LtE:'LessThan', ast.Gt:'StrictGreaterThan',
                     ast.GtE:'GreaterThan', ast.Eq:'Equality', ast.NotEq:'Unequality'}
            parts = [n.left, *n.comparators]
            rels = tuple(Node(heads[type(op)], (walk(parts[i]), walk(parts[i+1])))
                         for i, op in enumerate(n.ops))
            return rels[0] if len(rels) == 1 else Node('And', rels)
        if isinstance(n, ast.Call) and not n.keywords:
            return Node('Apply', (walk(n.func), *(walk(x) for x in n.args)))
        raise InputError('outside sstr grammar: '+type(n).__name__)
    try:
        return walk(ast.parse(source, mode='eval').body)
    except (SyntaxError, ValueError, KeyError) as exc:
        raise InputError('SymPy grammar: '+str(exc)) from exc


@dataclass(frozen=True)
class NameRow:
    py: str
    wl: str
    object: str
    py_line: str
    wl_line: str
    arity: int = 0


# Each row is a substitution naming an object, never a residual condition.
# Citations refer to the two accepted S12_order_source engines at a3da43d9.
NAME_TABLE = (
    NameRow('t','t','time coordinate','131','92'),
    NameRow('x1','x1','first in-plane coordinate','132','92'),
    NameRow('x2','x2','second in-plane coordinate','132','92'),
    NameRow('x3','x3','third in-plane coordinate','132','92'),
    NameRow('w','w','bulk coordinate','132-134','92'),
    NameRow('n','n','generic material density','136,141','93',5),
    NameRow('chi','chi','generic ordered fraction','136,142','96',5),
    NameRow('f','f','live interior order','136,163','111',5),
    NameRow('G','g','conditional order in two offsets','136,164','112',2),
    NameRow('u1','u1','generic first velocity component','137,143','94',5),
    NameRow('u2','u2','generic second velocity component','137,143','94',5),
    NameRow('u3','u3','generic third velocity component','137,143','94',5),
    NameRow('uw','um','generic bulk velocity component','137,143','94',5),
    NameRow('j1','j1','generic first order flux component','138,144','95',5),
    NameRow('j2','j2','generic second order flux component','138,144','95',5),
    NameRow('j3','j3','generic third order flux component','138,144','95',5),
    NameRow('jw','jm','generic bulk order flux component','138,144','95',5),
    NameRow('wminus','wm','lower live face','139,146','97',4),
    NameRow('wplus','wp','upper live face','139,146','97',4),
    NameRow('GammaB','gb','signed conversion rate','140,145','103-104',5),
    NameRow('GammaReturn','gr','return rate','140,145','104',5),
    NameRow('GammaDrain','gd','drain rate','140,145','104',5),
    NameRow('n_region1','nRegion1','density below slab','179,188','116,124',5),
    NameRow('n_region2','nRegion2','density inside slab','179,188','116,124',5),
    NameRow('n_region3','nRegion3','density above slab','179,188','116,124',5),
    NameRow('u1_region1','u1Region1','first velocity below slab','180-183,189','117-118,125',5),
    NameRow('u2_region1','u2Region1','second velocity below slab','180-183,189','117-118,125',5),
    NameRow('u3_region1','u3Region1','third velocity below slab','180-183,189','117-118,125',5),
    NameRow('u4_region1','u4Region1','bulk velocity below slab','180-183,189','117-118,125',5),
    NameRow('u1_region2','u1Region2','first velocity inside slab','180-183,189','117-118,125',5),
    NameRow('u2_region2','u2Region2','second velocity inside slab','180-183,189','117-118,125',5),
    NameRow('u3_region2','u3Region2','third velocity inside slab','180-183,189','117-118,125',5),
    NameRow('u4_region2','u4Region2','bulk velocity inside slab','180-183,189','117-118,125',5),
    NameRow('u1_region3','u1Region3','first velocity above slab','180-183,189','117-118,125',5),
    NameRow('u2_region3','u2Region3','second velocity above slab','180-183,189','117-118,125',5),
    NameRow('u3_region3','u3Region3','third velocity above slab','180-183,189','117-118,125',5),
    NameRow('u4_region3','u4Region3','bulk velocity above slab','180-183,189','117-118,125',5),
    NameRow('j1_region1','j1Region1','first flux below slab','184-187,190','119-120,126',5),
    NameRow('j2_region1','j2Region1','second flux below slab','184-187,190','119-120,126',5),
    NameRow('j3_region1','j3Region1','third flux below slab','184-187,190','119-120,126',5),
    NameRow('j4_region1','j4Region1','bulk flux below slab','184-187,190','119-120,126',5),
    NameRow('j1_region2','j1Region2','first flux inside slab','184-187,190','119-120,126',5),
    NameRow('j2_region2','j2Region2','second flux inside slab','184-187,190','119-120,126',5),
    NameRow('j3_region2','j3Region2','third flux inside slab','184-187,190','119-120,126',5),
    NameRow('j4_region2','j4Region2','bulk flux inside slab','184-187,190','119-120,126',5),
    NameRow('j1_region3','j1Region3','first flux above slab','184-187,190','119-120,126',5),
    NameRow('j2_region3','j2Region3','second flux above slab','184-187,190','119-120,126',5),
    NameRow('j3_region3','j3Region3','third flux above slab','184-187,190','119-120,126',5),
    NameRow('j4_region3','j4Region3','bulk flux above slab','184-187,190','119-120,126',5),
    NameRow('j_length','jLength','unknown flux length exponent','228-230','150'),
    NameRow('j_time','jTime','unknown flux time exponent','228-230','150'),
)


def validate_tables(names, correspondences):
    for field in ('py','wl'):
        keys = [getattr(r,field) for r in names]
        if len(keys) != len(set(keys)):
            raise InputError('duplicate name row: '+field)
    keys = [c.tag for c in correspondences]
    if len(keys) != len(set(keys)):
        raise InputError('duplicate join/component row')


def census(n):
    """Independent parse census; every terminal, including empty lists, counts."""
    return sum(census(x) for x in n.args) if n.args else 1


def leaf_paths(n, prefix=()):
    if not n.args:
        return {prefix}
    return set().union(*(leaf_paths(x,prefix+(i,)) for i,x in enumerate(n.args)))


def node_at(n, path):
    for i in path:
        n = n.args[i]
    return n


def data(n):
    return [n.head,n.value,[data(a) for a in n.args]]


LABELS = {'engine','governing_input','calculus_scope','dimension_basis','completion','execution'}


def is_label(tag):
    return tag in LABELS or tag == 'WORKING_ON' or tag.endswith(
        ('/comparison','/dimension_comparison','/premise','_rhs/dimension'))


@dataclass
class Entry:
    tag: str
    raw: str
    line: int
    node: Node


class Stream:
    def __init__(self,path,engine):
        self.path, self.engine = Path(path), engine
        self.index = defaultdict(list)
        try:
            raw = self.path.read_bytes()
        except OSError as exc:
            raise InputError(str(exc)) from exc
        self.sha256 = hashlib.sha256(raw).hexdigest()
        entries = []
        for i,line in enumerate(raw.decode().splitlines(),1):
            m = re.fullmatch(r'([A-Za-z0-9_/]+): ?(.*)',line)
            if m:
                entries.append([m[1],m[2],i])
            elif entries and (line.startswith((' ','[',']','(' ,')','{','}'))):
                entries[-1][1] += '\n'+line
            else:
                raise InputError(f'{path}:{i}: unframed input')
        if not entries:
            raise InputError('empty input '+str(path))
        for tag,text,line in entries:
            try:
                if is_label(tag) and engine == 'py':
                    node = atom('Text',text)
                else:
                    node = (py_parse if engine == 'py' else wl_parse)(text)
            except InputError as exc:
                raise InputError(f'{path}:{line}:{tag}: {exc}') from exc
            self.index[tag].append(Entry(tag,text,line,node))


class Bindings:
    def __init__(self,names):
        self.maps = {'py':{r.py:r.py for r in names}, 'wl':{r.wl:r.py for r in names}}
        self.arities = {r.py:r.arity for r in names}

    def name(self,name,engine):
        return self.maps[engine].get(name,engine+'::'+name)


RELATIONS = {'Equality':s.Eq,'Unequality':s.Ne,'StrictLessThan':s.Lt,
             'StrictGreaterThan':s.Gt,'LessThan':s.Le,'GreaterThan':s.Ge}
BUILTINS = {'Eq':'Equality','Ne':'Unequality','Lt':'StrictLessThan','Le':'LessThan',
            'Gt':'StrictGreaterThan','Ge':'GreaterThan','HeavisideTheta':'Heaviside',
            'Sin':'sin','Cos':'cos','Exp':'exp','Sqrt':'sqrt','Rational':'Rational'}


def call(n,name):
    return n.head=='Apply' and n.args[0].head=='Symbol' and n.args[0].value==name


def math_node(n,b,engine,scope=None):
    scope = {} if scope is None else scope
    def cv(x,local=None):
        return math_node(x,b,engine,scope if local is None else local)
    h,a = n.head,n.args
    if h=='Number':
        return s.Rational(n.value)
    if h=='Boolean':
        return s.true if n.value.lower()=='true' else s.false
    if h=='Symbol':
        if n.value in scope:
            return scope[n.value]
        constants = {'True':s.true,'False':s.false,'true':s.true,'false':s.false,
                     'Reals':s.S.Reals,'pi':s.pi,'Pi':s.pi,'E':s.E}
        if n.value in constants:
            return constants[n.value]
        key = b.name(n.value,engine)
        if b.arities.get(key,0):
            raise Unsupported('unapplied live function '+key)
        return s.Symbol(key,real=True)
    if h=='Text':
        raise Unsupported('text in mathematical payload')
    if h=='Apply':
        head,*args = a
        # Native multiargument partial derivative, with formal slots before
        # substituting composite arguments. No function becomes a bare symbol.
        if head.head=='Apply' and call(head.args[0],'Derivative'):
            orders = head.args[0].args[1:]
            fun = head.args[1]
            if fun.head!='Symbol' or len(orders)!=len(args):
                raise Unsupported('native derivative arity')
            key = b.name(fun.value,engine)
            if b.arities.get(key)!=len(args):
                raise Unsupported('unbound derivative function/arity '+key)
            slots = tuple(s.Dummy('arg'+str(i),real=True) for i in range(len(args)))
            deriv = s.Derivative(s.Function(key,real=True)(*slots),
                                 *((v,int(cv(k))) for v,k in zip(slots,orders) if int(cv(k))),
                                 evaluate=False)
            points = tuple(cv(x) for x in args)
            # Coordinate-only partials use their original coordinate variables.
            if all(isinstance(x,s.Symbol) for x in points) and len(set(points))==len(points):
                return s.Derivative(s.Function(key,real=True)(*points),
                       *((v,int(cv(k))) for v,k in zip(points,orders) if int(cv(k))),evaluate=False)
            return s.Subs(deriv,slots,points)
        if head.head!='Symbol':
            raise Unsupported('unsupported composite application head')
        name = head.value
        if name in b.maps[engine] and b.arities.get(b.name(name,engine),0):
            key = b.name(name,engine)
            if len(args)!=b.arities[key]:
                raise Unsupported('live function arity '+key)
            return s.Function(key,real=True)(*(cv(x) for x in args))
        h,a = BUILTINS.get(name,name),tuple(args)
    if h in ('Add','Mul','Pow'):
        return getattr(s,h)(*(cv(x) for x in a))
    if h in RELATIONS:
        if len(a)!=2:
            raise Unsupported('relation arity')
        if a[0].head in RELATIONS:
            return s.And(cv(a[0]),RELATIONS[h](cv(a[0].args[-1]),cv(a[1]),evaluate=False),evaluate=False)
        return RELATIONS[h](*(cv(x) for x in a),evaluate=False)
    if h in ('And','Or','Not'):
        return getattr(s,h)(*(cv(x) for x in a),evaluate=False)
    if h=='Element':
        def alternatives(x):
            return sum((alternatives(y) for y in x.args),()) if x.head=='Alternatives' else (x,)
        variables = alternatives(a[0])
        return s.And(*(s.Contains(cv(v),cv(a[1]),evaluate=False) for v in variables),evaluate=False)
    if h in ('Sequence','Alternatives'):
        return s.Tuple(*(cv(x) for x in a))
    if h in ('Heaviside','DiracDelta','Abs','sin','cos','exp','sqrt','Rational'):
        return getattr(s,h)(*(cv(x) for x in a))
    if h=='Derivative':
        specs = []
        for x in a[1:]:
            specs.append((cv(x.args[0]),int(cv(x.args[1]))) if x.head=='Sequence' else cv(x))
        return s.Derivative(cv(a[0]),*specs,evaluate=False)
    if h=='Subs':
        def members(x):
            return x.args if x.head=='Sequence' else (x,)
        variables,points = members(a[1]),members(a[2])
        if len(variables)!=len(points):
            raise Unsupported('substitution arity')
        local=dict(scope)
        for i,v in enumerate(variables):
            local[v.value]=s.Symbol('_slot'+str(len(scope)+i),real=True)
        return s.Subs(cv(a[0],local),tuple(local[v.value] for v in variables),tuple(cv(x) for x in points))
    if h=='Piecewise':
        branches=a[0].args if engine=='wl' else a
        pairs=[]
        for pair in branches:
            if pair.head!='Sequence' or len(pair.args)!=2:
                raise Unsupported('piecewise branch grammar')
            pairs.append((cv(pair.args[0]),cv(pair.args[1])))
        if engine=='wl':
            pairs.append((cv(a[1]) if len(a)>1 else s.S.Zero,s.true))
        return s.Piecewise(*pairs,evaluate=False)
    raise Unsupported('unsupported mathematical head '+h)


# Explicit semantic schemas for container components. A schema is checked in
# full before any child comparison, so missing siblings cannot shift a pairing.
SCALAR = None
SPACE = ('x1','x2','x3','w')
REGIONS = ('below','inside','above')
UNIT = ('length exponent','time exponent')
FIELDS = (('density',SCALAR),('velocity',SPACE),('order flux',SPACE),('order',SCALAR))


@dataclass(frozen=True)
class Correspondence:
    tag: str
    kind: str
    object: str
    py_line: str
    wl_line: str
    schema: object = None
    # An explicit override is also used for synthetic component-repoint controls.
    paths: tuple = ()


def contracts():
    rows=[]
    def add(tag,kind,obj,py,wl,schema=None):
        rows.append(Correspondence(tag,kind,obj,str(py),str(wl),schema))
    for tag,py,wl in [('engine',122,84),('governing_input',123,85),('calculus_scope','124-130','86-91'),
                      ('dimension_basis',232,149),('completion',371,285)]:
        add(tag,'text',tag,py,wl)
    for tag,py,wl in [('supplied/material_balance',203,133),('supplied/order_balance',204,134),
                      ('supplied/conversion_identity',205,135)]:
        add(tag,'value',tag,py,wl)
    add('supplied/faces','tree','lower and upper live faces',206,136,('lower face','upper face'))
    add('supplied/domain','tree','face order/coordinate assumptions; interior range',207,137,
        ('face order and real coordinates','interior range'))
    add('live_fields','tree','generic density, spatial velocity/flux, fraction, faces',208,138,
        (('density',None),('velocity',SPACE),('order flux',SPACE),('fraction',None),('faces',('lower','upper'))))
    for tag,obj,py,wl,sch in [
        ('region_domains','three open region domains',209,139,REGIONS),
        ('fields/region_indicators','indicators of the three open regions',210,140,REGIONS),
        ('density/region_fields','regional densities',211,141,REGIONS),
        ('velocity/region_fields','regional spatial velocities',212,142,tuple((r,SPACE) for r in REGIONS)),
        ('flux/region_fields','regional spatial fluxes',213,143,tuple((r,SPACE) for r in REGIONS)),
        ('density/distribution_field','density distribution',214,144,None),
        ('velocity/distribution_field','spatial velocity distributions',215,145,SPACE),
        ('flux/distribution_field','spatial flux distributions',216,146,SPACE)]:
        add(tag,'value' if sch is None else 'tree',obj,py,wl,sch)
    add('order_flux_dimension_inference/operands','inference','expanded order terms and reference with units',237,155)
    add('order_flux_dimension_inference/residual','inference_residual','term-associated unit residuals',239,157)
    add('order_flux_dimension_inference/object','unit','inferred flux units',244,161)
    add('supplied/material_balance_rhs','value','material balance right hand side',247,163)
    add('supplied/material_balance_rhs/dimension','text','dimension note',248,164)
    add('supplied/order_balance_rhs/dimension_comparison','text','dimension comparison note',
        '111-113,250','63-64,166')
    add('material/region_current_operands','regional_current','regional conserved mass currents: spatial components only',258,176)
    add('material/region_balance_operands','tree','regional material equations',259,177,REGIONS)
    add('material/conserved_current_operand','current','conserved mass current: spatial components only',260,182)
    for tag,py,wl in [('balance_distribution_operand',262,185),('balance_distribution',266,187),
                     ('face_requirements',291,212),('balance_with_face_traces',292,213)]:
        add('material/'+tag,'tree' if tag=='face_requirements' else 'value','material '+tag,py,wl,
            ('lower','upper') if tag=='face_requirements' else None)
    for r in (1,2,3):
        add(f'material/region/{r}/balance','value',f'material equation in region {r}',270,191)
    for k in (1,2):
        for tag,kind,py,wl,sch in [('coefficient_operand','value',280,201,None),('trace_coefficient','value',281,202,None),
            ('balance_operands','tree',282,203,('trace','material rhs')),('balance_residual','value',283,204,None),
            ('requirement','value',285,206,None),('balance_contribution','value',286,207,None)]:
            add(f'material/face/{k}/'+tag,kind,f'material face {k} '+tag,py,wl,sch)
    for package in ('general_order','conditional_fast_ordering'):
        for tag,kind,py,wl,sch in [
            ('order_field','value',297,221,None),('region_order_operands','tree',298,222,REGIONS),
            ('region_current_operands','regional_current',307,233,None),
            ('conserved_current_operand','current',312,237,None),
            ('material_elimination_operand','incompatible',323,242,None),
            ('balance_source_operand','value',326,241,None),('source_distribution','value',328,244,None),
            ('regular_distribution_part','value',330,246,None),('source_with_face_traces','value',370,283,None)]:
            add(package+'/'+tag,kind,package+' '+tag,py,wl,sch)
        for r in (1,2,3):
            for tag,kind,py,wl,sch in [('domain','value',338,255,None),('field_operands','tree',339,256,FIELDS),
                                     ('source','value',341,259,None)]:
                add(f'{package}/region/{r}/'+tag,kind,f'region {r} '+tag,py,wl,sch)
        for k in (1,2):
            for tag,kind,py,wl,sch in [
                ('locus','value',359,272,None),('field_trace_operands','tree',360,273,
                        tuple((r,FIELDS) for r in REGIONS[k-1:k+1])),
                ('coefficient_operand','value',363,275,None),('trace_coefficient','value',364,276,None),
                ('source_contribution','value',365,277,None),('material_requirement','value',366,278,None)]:
                add(f'{package}/face/{k}/'+tag,kind,f'face {k} '+tag,py,wl,sch)
    add('conditional_fast_ordering/premise','text','conditional specialization label',300,224)
    add('conditional_fast_ordering/interior_range','value','conditional interior range',301,225)
    return tuple(rows)


CORRESPONDENCES = contracts()


def correspondence(tag,rows):
    bytag={r.tag:r for r in rows}
    if tag in bytag:
        return bytag[tag]
    suffix=tag.rsplit('/',1)[-1]
    base=tag.rsplit('/',1)[0]
    if suffix in ('comparison','dimension_comparison'):
        parent=correspondence(base,rows)
        if parent:
            return Correspondence(tag,'text','comparison note for '+parent.object,
                                  '35-37; '+parent.py_line,'21-23; '+parent.wl_line)
    if suffix in ('dimension_operands','dimension_residual','dimension_reference_operand',
                  'nested_dimension_operands','nested_dimension_residual'):
        # Only dimension reports actually emitted by the sources are bound.
        roots={'supplied/material_balance_lhs':('246','162'), 'supplied/order_balance_lhs':('249','165'),
               'supplied/order_balance_rhs':('250','166'), 'supplied/conversion_identity_rhs':('251','167')}
        parent=correspondence(base,rows)
        if base.endswith('/balance_lhs'):
            parent=correspondence(base[:-4],rows)
        if parent or base in roots:
            py,wl=roots.get(base,(parent.py_line,parent.wl_line) if parent else ('',''))
            kind={'dimension_operands':'ledger','dimension_residual':'ledger_residual',
                  'dimension_reference_operand':'reference','nested_dimension_operands':'nested_ledger',
                  'nested_dimension_residual':'nested_ledger'}[suffix]
            return Correspondence(tag,kind,'dimension report '+base,'97-117; '+py,'56-68; '+wl)
    return None


COORDS = s.symbols('t x1 x2 x3 w',real=True)
T,X1,X2,X3,W = COORDS
WM = s.Function('wminus',real=True)(*COORDS[:-1])
WP = s.Function('wplus',real=True)(*COORDS[:-1])
FORMS = ('affine','quadratic','rational')
STRATA = ('generic','jumps_zero_lower','jumps_zero_upper','order_zero_lower','order_zero_upper',
          'jumps_and_order_zero_lower','jumps_and_order_zero_upper')
BASE_POINTS = ((s.Rational(1,7),s.Rational(2,9),s.Rational(-1,5),s.Rational(3,11)),
               (s.Rational(-2,7),s.Rational(1,6),s.Rational(3,10),s.Rational(-1,8)),
               (s.Rational(3,8),s.Rational(-2,11),s.Rational(1,9),s.Rational(2,7)))
FORMAL = s.symbols('a0:5',real=True)


def basic_form(form,args,seed):
    weighted = sum((i+1)*(seed+1)*a for i,a in enumerate(args))
    if form=='affine':
        return s.Rational(seed+1,7)+weighted/31
    if form=='quadratic':
        return s.Rational(seed+1,7)+sum((i+1)**2*(seed+1)*a*a+(i+1)*a
                                        for i,a in enumerate(args))/37
    if form=='rational':
        return s.Rational(seed+1,7)+sum((i+1)*(seed+1)*a/(1+a*a)
                                        for i,a in enumerate(args))/41
    raise InputError('unknown test form '+form)


@lru_cache(None)
def face_forms(form):
    lower = -1+basic_form(form,COORDS[:-1],1)/100
    upper = lower+3+(1+sum((i+1)*x*x for i,x in enumerate(COORDS[:-1])))/100
    return lower,upper


@lru_cache(None)
def profile(name,arity,form,stratum):
    """Live, argument-separating probes. Zero strata retain nonconstant fields."""
    args=FORMAL[:arity]
    lower,upper=face_forms(form)
    toargs=dict(zip(COORDS,args))
    lo,hi=lower.xreplace(toargs),upper.xreplace(toargs)
    if name in ('wminus','wplus'):
        return (lo,hi)[name=='wplus']
    seed=1+int(hashlib.sha256(name.encode()).hexdigest()[:4],16)%13
    base=basic_form(form,args,seed)
    if name in ('f','G','chi'):
        # These finite probes are checked against the supplied ranges at every
        # printed point. Bounds are not asserted globally for polynomial forms.
        base=s.Rational(1,3)+(base-s.Rational(seed+1,7))/10000
        if name in ('f','G') and 'order_zero' in stratum:
            side='lower' if stratum.endswith('lower') else 'upper'
            distance=(args[4]-lo if side=='lower' else hi-args[4]) if name=='f' else args[side=='upper']
            # Linear zero at the physical face, positive on its interior side.
            # The additional quadratic term makes the finite exterior probes
            # stay in the printed range too. Argument derivatives at the face
            # remain distinct; no live field is made constant on a stratum.
            base=distance*(1+8*distance)*(1+base)/1000
        return base
    if '_region' in name and ('jumps_zero' in stratum or 'jumps_and_order_zero' in stratum):
        stem,region=name.split('_region')
        common_seed=1+int(hashlib.sha256(stem.encode()).hexdigest()[:4],16)%13
        base=basic_form(form,args,common_seed)
        distance=args[4]-(lo if stratum.endswith('lower') else hi)
        base+=int(region)*distance**2*(1+sum((i+1)*a*a for i,a in enumerate(args)))/47
    return base


@lru_cache(None)
def test_expression(expr,form,stratum):
    """Substitute entire functions before derivatives and Subs are evaluated."""
    def replace_function(fun):
        name=fun.func.__name__
        # math_node constructs undefined functions only from explicit table
        # entries. Synthetic table rows use the same live probe generator.
        if not fun.args:
            raise Unsupported('unapplied function in evaluation '+name)
        return profile(name,len(fun.args),form,stratum).xreplace(dict(zip(FORMAL,fun.args)))
    value=expr.replace(lambda x:isinstance(x,AppliedUndef),replace_function)
    return value.doit()


@lru_cache(None)
def point_for(form,index,location):
    lower,upper=face_forms(form)
    env=dict(zip(COORDS[:-1],BASE_POINTS[index]))
    lo,hi=lower.subs(env),upper.subs(env)
    positions={'below':lo-s.Rational(index+1,4),'inside':lo+(hi-lo)*s.Rational(index+1,4),
               'above':hi+s.Rational(index+1,4),'lower':lo,'upper':hi}
    env[W]=positions[location]
    # Dimension unknowns are not solved/fixed by the comparator.
    env[s.Symbol('j_length',real=True)]=s.Rational(2+index,7)
    env[s.Symbol('j_time',real=True)]=s.Rational(-3-index,11)
    return env


@lru_cache(None)
def evaluate(expr,form,stratum,index,location):
    value=test_expression(expr,form,stratum).subs(point_for(form,index,location)).doit()
    if isinstance(value,s.Expr):
        value=s.cancel(value)
    if value.has(AppliedUndef,s.Derivative,s.Subs) or value.free_symbols:
        raise Unsupported('evaluation retains a live function, derivative or free symbol: '+str(value))
    if value.has(s.zoo,s.nan,s.oo,-s.oo):
        raise Unsupported('non-finite evaluation')
    return value


def supplied_domain_value(expr):
    """Only strict live-face order and real coordinates; no field equations.

    This operates on a residual already emitted as computed. It never changes
    either operand or identifies two representations before subtraction.
    """
    gap=WP-WM
    replacements={s.Lt(WM,WP,evaluate=False):s.true,s.Gt(WP,WM,evaluate=False):s.true}
    replacements.update((s.Contains(z,s.S.Reals,evaluate=False),s.true) for z in COORDS)
    for h in expr.atoms(s.Heaviside):
        arg=s.expand(h.args[0])
        if s.expand(arg-gap)==0:
            replacements[h]=s.S.One
        elif s.expand(arg+gap)==0:
            replacements[h]=s.S.Zero
    for a in expr.atoms(s.Abs):
        if s.expand(a.args[0]-gap)==0 or s.expand(a.args[0]+gap)==0:
            replacements[a]=gap
    value=expr.xreplace(replacements)
    return s.expand(value) if isinstance(value,s.Expr) else value


def exact_difference(a,b):
    if isinstance(a,Boolean) and not isinstance(a,s.Expr):
        return s.Xor(a,b,evaluate=False)
    return s.expand(a-b)


def text_difference(a,b):
    return {'equal_text':a==b,'first_different_character':next((i for i,(x,y) in enumerate(zip(a,b)) if x!=y),
                                                             min(len(a),len(b)) if a!=b else None)}


def region_part(expr,region):
    """Strict open-region restrictions of printed steps, not face point values."""
    replacements={}
    unknown=[]
    for step in expr.atoms(s.Heaviside):
        arg=step.args[0]
        found=False
        for face_index,face in enumerate((WM,WP)):
            distance=W-face
            for sign in (1,-1):
                ratio=s.cancel(arg/(sign*distance))
                if ratio.is_number and ratio.is_positive:
                    positive=region>face_index+1
                    replacements[step]=s.Integer(positive if sign==1 else not positive)
                    found=True
                    break
            if found:
                break
        if not found:
            unknown.append(str(step))
    value=s.expand(expr.xreplace(replacements))
    return value,unknown


def distribution_parts(expr):
    """Extract from this operand alone; include 1/|d argument/dw| on its face.

    Higher deltas, nonlinear roots, products of singular atoms and unknown
    step supports are visible coverage outcomes, never silently dropped.
    """
    expr=s.expand(expr)
    deltas=expr.atoms(s.DiracDelta)
    regular=s.expand(expr.xreplace(dict.fromkeys(deltas,s.S.Zero)))
    parts={}
    coverage=[]
    for r in (1,2,3):
        val,unknown=region_part(regular,r)
        parts['regular/'+REGIONS[r-1]]=val
        if unknown:
            coverage.append({'part':'regular/'+REGIONS[r-1],'unisolated_steps':unknown})
    singular=expr-regular
    removed=s.S.Zero
    for d in sorted(deltas,key=str):
        if len(d.args)>1 and d.args[1]!=0:
            coverage.append({'singular_atom':str(d),'reason':'higher delta derivative is not isolated'})
            continue
        arg=d.args[0]
        slope=s.diff(arg,W)
        if slope.has(W) or slope==0:
            coverage.append({'singular_atom':str(d),'reason':'not a nonzero affine argument in w'})
            continue
        root=s.cancel(W-arg/slope)
        face=next((i for i,f in enumerate((WM,WP)) if s.expand(root-f)==0),None)
        if face is None:
            coverage.append({'singular_atom':str(d),'reason':'support is not either supplied live face'})
            continue
        coefficient=expr.coeff(d)
        if coefficient.has(s.DiracDelta) or expr.has(d**2):
            coverage.append({'singular_atom':str(d),'reason':'nonlinear singular product'})
            continue
        trace=s.Subs(coefficient,W,(WM,WP)[face]).doit()/s.Abs(slope)
        key='delta/'+('lower','upper')[face]
        parts[key]=parts.get(key,s.S.Zero)+trace
        removed+=coefficient*d
    for key in ('delta/lower','delta/upper'):
        parts.setdefault(key,s.S.Zero)
    leftover=s.expand(singular-removed)
    if leftover!=0:
        coverage.append({'unisolated_content':str(leftover)})
    return parts,coverage


def locations(tag,part):
    if part.startswith('regular/'):
        return (part.split('/')[-1],)
    if part.startswith('delta/'):
        return (part.split('/')[-1],)
    face=re.search(r'/face/([12])/',tag)
    if face:
        return (('lower','upper')[int(face[1])-1],)
    region=re.search(r'/region/([123])/',tag)
    if region:
        return (REGIONS[int(region[1])-1],)
    return ('inside',)


def sample_values(a,b,tag,part,*,forms=FORMS,points=range(3),strata=STRATA):
    rows=[]
    for form in forms:
        for stratum in strata:
            for index in points:
                for location in locations(tag,part):
                    row={'form':form,'stratum':stratum,'point_index':index,'location':location,
                         'point':{str(k):str(v) for k,v in point_for(form,index,location).items()}}
                    try:
                        av,bv=evaluate(a,form,stratum,index,location),evaluate(b,form,stratum,index,location)
                        residual=exact_difference(av,bv)
                        row.update(py_value=str(av),wl_value=str(bv),residual=str(residual))
                    except Unsupported as exc:
                        row.update(unresolved_coverage=str(exc))
                    rows.append(row)
    return rows


def compare_values(a,b,tag,part='value',sample_options=None):
    """No target check and no verdict; both computed residuals remain present."""
    row={'part':part,'py_operand':str(a),'wl_operand':str(b)}
    if isinstance(a,s.Equality) or isinstance(b,s.Equality):
        if not isinstance(a,s.Equality) or not isinstance(b,s.Equality):
            return dict(row,not_formed='different kinds of relation')
        return dict(row,relation_head=type(a).__name__,components=[
            compare_values(a.lhs,b.lhs,tag,part+'/lhs',sample_options),
            compare_values(a.rhs,b.rhs,tag,part+'/rhs',sample_options)])
    if (isinstance(a,s.Expr)) != (isinstance(b,s.Expr)):
        return dict(row,not_formed='different kinds of mathematical object')
    try:
        residual=exact_difference(a,b)
    except (TypeError,ValueError,NotImplementedError) as exc:
        # A genuinely unsupported symbolic subtraction still gets probes.
        row['exact_unresolved']=str(exc)
    else:
        row['as_computed_residual']=str(residual)
        if ACTIVE_OUTPUT is not None:
            emit(ACTIVE_OUTPUT,{'kind':'residual_before_domain','tag':tag,'part':part,
                                'py_operand':str(a),'wl_operand':str(b),'as_computed_residual':str(residual)})
        row['on_supplied_domain_residual']=str(supplied_domain_value(residual))
    row['samples']=sample_values(a,b,tag,part,**(sample_options or {}))
    return row


def compare_scalar(a,b,tag,binding,sample_options=None):
    # Equations are structural components; distributions beneath each side are
    # still separated. This is also the path for standalone dimension values.
    if isinstance(a,s.Equality) or isinstance(b,s.Equality):
        if not isinstance(a,s.Equality) or not isinstance(b,s.Equality):
            return {'not_formed':'different relation kinds','py_operand':str(a),'wl_operand':str(b)}
        return {'relation_head':type(a).__name__,'components':[
            compare_scalar(a.lhs,b.lhs,tag,binding+'/lhs',sample_options),
            compare_scalar(a.rhs,b.rhs,tag,binding+'/rhs',sample_options)]}
    if a.has(s.DiracDelta) or b.has(s.DiracDelta) or any(h.args[0].has(W) for h in a.atoms(s.Heaviside)|b.atoms(s.Heaviside)):
        pa,ca=distribution_parts(a)
        pb,cb=distribution_parts(b)
        return {'distribution_parts':[
            compare_values(pa[key],pb[key],tag,key,sample_options) for key in pa],
            'py_coverage':ca,'wl_coverage':cb,
            'normalization':'coefficient on support divided by absolute derivative of delta argument in w'}
    return compare_values(a,b,tag,binding,sample_options)


def schema_children(schema):
    return tuple((x,None) if isinstance(x,str) else x for x in schema)


def shaped(n,schema):
    if schema is None:
        return n.head not in ('Sequence','Dictionary') and not call(n,'Matrix')
    return n.head=='Sequence' and len(n.args)==len(schema) and all(
        shaped(child,sub) for child,(_,sub) in zip(n.args,schema_children(schema)))


def schema_paths(schema,prefix=(),name=''):
    if schema is None:
        return [(name,prefix,prefix)]
    if not schema:
        return [(name,prefix,prefix)]
    result=[]
    for i,(label,sub) in enumerate(schema_children(schema)):
        result+=schema_paths(sub,prefix+(i,),name+'/'+label)
    return result


def matrix_unit(n):
    if call(n,'Matrix'):
        seqnode=n.args[1]
        if seqnode.head!='Sequence' or len(seqnode.args)!=2 or any(
                row.head!='Sequence' or len(row.args)!=1 for row in seqnode.args):
            raise Unsupported('dimension matrix is not a two-entry column')
        return [(node_at(n,(1,i,0)),(1,i,0)) for i in range(2)]
    if n.head=='Sequence' and len(n.args)==2 and all(x.head!='Sequence' for x in n.args):
        return [(node_at(n,(i,)),(i,)) for i in range(2)]
    raise Unsupported('dimension value is not a length/time pair')


class Comparator:
    def __init__(self,names,sample_options=None):
        self.bindings=Bindings(names)
        self.sample_options=sample_options
        self.compared=[set(),set()]

    def scalar(self,p,w,tag,name,pp=(),wp=()):
        try:
            a=math_node(p,self.bindings,'py')
            b=math_node(w,self.bindings,'wl')
        except Unsupported as exc:
            return {'component':name,'not_compared':str(exc),'py_parse':data(p),'wl_parse':data(w)}
        if isinstance(a,s.Tuple) or isinstance(b,s.Tuple):
            return {'component':name,'not_formed':'container found where scalar component is bound'}
        ans=compare_scalar(a,b,tag,name,self.sample_options)
        if not has_key(ans,'not_formed') and not ans.get('py_coverage') and not ans.get('wl_coverage'):
            self.compared[0].update(leaf_paths(p,pp))
            self.compared[1].update(leaf_paths(w,wp))
        return dict(component=name,**ans)

    def unit(self,p,w,tag,name,pp=(),wp=()):
        try:
            pa,wa=matrix_unit(p),matrix_unit(w)
        except Unsupported as exc:
            return {'component':name,'not_formed':str(exc)}
        return {'component':name,'components':[self.scalar(x,y,tag,name+'/'+UNIT[i],pp+px,wp+wx)
                   for i,((x,px),(y,wx)) in enumerate(zip(pa,wa))]}

    def tree(self,p,w,tag,c):
        if not shaped(p,c.schema) or not shaped(w,c.schema):
            return {'not_formed':'component counts or kinds differ from the cited semantic schema',
                    'schema':c.schema,'py_parse':data(p),'wl_parse':data(w)}
        if c.schema==():
            self.compared[0].update(leaf_paths(p))
            self.compared[1].update(leaf_paths(w))
            return {'empty_container_comparison':True}
        return {'components':[self.scalar(node_at(p,pp),node_at(w,wp),tag,name,pp,wp)
                              for name,pp,wp in (c.paths or schema_paths(c.schema))]}

    def currents(self,p,w,tag,regional):
        ps=tuple((r,('density',)+SPACE) for r in REGIONS) if regional else ('density',)+SPACE
        ws=tuple((r,SPACE) for r in REGIONS) if regional else SPACE
        if not shaped(p,ps) or not shaped(w,ws):
            return {'not_formed':'conserved-current component counts/kinds differ from the cited schemas'}
        paths=[]
        extras=[]
        for r in range(3) if regional else (None,):
            prefix=(r,) if r is not None else ()
            stem=REGIONS[r]+'/' if r is not None else ''
            for i,label in enumerate(SPACE):
                paths.append((stem+'spatial/'+label,prefix+(i+1,),prefix+(i,)))
            extras.append({'engine':'py','path':prefix+(0,), 'object':stem+'time-density',
                           'reason':'the Wolfram payload under this tag contains spatial currents only',
                           'parse':data(node_at(p,prefix+(0,)))})
        return {'components':[self.scalar(node_at(p,pp),node_at(w,wp),tag,name,pp,wp) for name,pp,wp in paths],
                'unpaired_components':extras}

    def ledger(self,p,w,tag,residual=False,operands=None,pp=(),wp=()):
        """Term identity is the binding; CAS ordering is never a binding.

        Only exactly identical translated term objects identify ledger rows.
        A changed/missing term stays unpaired. No fuzzy/algebraic pairing.
        """
        if p.head!='Sequence' or w.head!='Sequence':
            return {'not_formed':'dimension ledger is not a list'}
        refs=operands or (p,w)
        if any(r.head!='Sequence' for r in refs):
            return {'not_formed':'dimension residual lacks a term ledger'}
        maps=[]
        for side,(n,ref,engine) in enumerate(zip((p,w),refs,('py','wl'))):
            if len(n.args)!=len(ref.args):
                return {'not_formed':'unit count differs from its term count'}
            mapping=defaultdict(list)
            for i,pair in enumerate(ref.args):
                if pair.head!='Sequence' or len(pair.args)!=2:
                    return {'not_formed':'term/unit record has different kinds or counts'}
                try:
                    key=math_node(pair.args[0],self.bindings,engine)
                except Unsupported as exc:
                    key=(engine,i,str(exc))
                mapping[key].append(i)
            maps.append(mapping)
        components=[]
        extras=[]
        for key in sorted(set(maps[0])|set(maps[1]),key=str):
            inds=[m.get(key,[]) for m in maps]
            if any(len(x)!=1 for x in inds):
                extras.append({'term':str(key),'py_indices':inds[0],'wl_indices':inds[1],
                               'reason':'no unique identical term binding across these expression trees'})
                continue
            i,j=inds[0][0],inds[1][0]
            if residual:
                components.append(self.unit(p.args[i],w.args[j],tag,'unit residual of '+str(key),pp+(i,),wp+(j,)))
            else:
                components.append(self.scalar(p.args[i].args[0],w.args[j].args[0],tag,'term '+str(key),pp+(i,0),wp+(j,0)))
                components.append(self.unit(p.args[i].args[1],w.args[j].args[1],tag,'units of '+str(key),pp+(i,1),wp+(j,1)))
        if not p.args and not w.args:
            self.compared[0].add(pp)
            self.compared[1].add(wp)
        return {'components':components,'unpaired_components':extras}

    def compare(self,p,w,tag,c,streams):
        kind=c.kind
        if kind=='text':
            a=label_text(p,streams['py'].index[tag][0].raw)
            b=label_text(w,streams['wl'].index[tag][0].raw)
            result=text_difference(a,b)
            self.compared=[leaf_paths(p),leaf_paths(w)]
            return {'py_text':a,'wl_text':b,'text_difference':result}
        if kind=='incompatible':
            return {'not_formed':'SymPy prints a distributional bulk multiple; Wolfram prints time-derivative substitution rules. These are different kinds of object.'}
        if kind=='value':
            if not shaped(p,None) or not shaped(w,None):
                return {'not_formed':'scalar correspondence has a container operand'}
            return self.scalar(p,w,tag,c.object)
        if kind=='tree':
            return self.tree(p,w,tag,c)
        if kind=='current' or kind=='regional_current':
            return self.currents(p,w,tag,kind=='regional_current')
        if kind=='unit':
            return self.unit(p,w,tag,c.object)
        if kind=='reference':
            schema=('expression','unit')
            if p.head!='Sequence' or w.head!='Sequence' or len(p.args)!=2 or len(w.args)!=2:
                return {'not_formed':'dimension reference is not an expression/unit pair'}
            return {'components':[self.scalar(p.args[0],w.args[0],tag,'reference expression',(0,),(0,)),
                                  self.unit(p.args[1],w.args[1],tag,'reference unit',(1,),(1,))]}
        if kind=='nested_ledger':
            if p.head=='Sequence' and w.head=='Sequence' and not p.args and not w.args:
                self.compared=[leaf_paths(p),leaf_paths(w)]
                return {'empty_container_comparison':True}
            return {'not_compared':'nested ledgers follow each engine expression tree; no source-defined cross-engine node correspondence'}
        if kind in ('ledger','ledger_residual'):
            ops=None
            if kind=='ledger_residual':
                base=tag.rsplit('/',1)[0]+'/dimension_operands'
                if any(len(st.index.get(base,[]))!=1 for st in streams.values()):
                    return {'not_formed':'no unique term ledger for dimension residual'}
                ops=tuple(streams[e].index[base][0].node for e in ('py','wl'))
            return self.ledger(p,w,tag,kind=='ledger_residual',ops)
        if kind in ('inference','inference_residual'):
            if kind=='inference_residual':
                refs=tuple(st.index['order_flux_dimension_inference/operands'][0].node.args[0]
                           for st in (streams['py'],streams['wl']))
                return self.ledger(p,w,tag,True,refs)
            if p.head!='Sequence' or w.head!='Sequence' or len(p.args)!=2 or len(w.args)!=2:
                return {'not_formed':'inference payload is not a term ledger/reference pair'}
            results=[self.ledger(p.args[0],w.args[0],tag,pp=(0,),wp=(0,))]
            for side,n in enumerate((p.args[1],w.args[1])):
                if n.head!='Sequence' or len(n.args)!=2:
                    return {'not_formed':'inference reference is not an expression/unit pair'}
            results.extend([self.scalar(p.args[1].args[0],w.args[1].args[0],tag,'reference expression',(1,0),(1,0)),
                            self.unit(p.args[1].args[1],w.args[1].args[1],tag,'reference unit',(1,1),(1,1))])
            return {'components':results}
        raise InputError('unknown correspondence kind '+kind)


def label_text(n,raw):
    return n.value if n.head=='Text' else raw


def has_key(obj,key):
    if isinstance(obj,dict):
        return key in obj or any(has_key(v,key) for v in obj.values())
    if isinstance(obj,list):
        return any(has_key(v,key) for v in obj)
    return False


def emit(out,obj):
    out.write(json.dumps(obj,separators=(',',':'),ensure_ascii=False)+'\n')
    out.flush()


ACTIVE_OUTPUT = None


def run(py_path,wl_path,out,*,names=NAME_TABLE,correspondences=CORRESPONDENCES,progress=True,sample_options=None,probe_metadata=True):
    global ACTIVE_OUTPUT
    started=time.monotonic()
    validate_tables(names,correspondences)
    streams={e:Stream(path,e) for e,path in (('py',py_path),('wl',wl_path))}
    cmp=Comparator(names,sample_options)
    emit(out,{'kind':'configuration','name_table':[asdict(r) for r in names],
              'sources':{'py':'scripts/S12_order_source.py at a3da43d9',
                         'wl':'mathematica/S12_order_source.wl at a3da43d9'},
              'input_sha256':{e:st.sha256 for e,st in streams.items()},
              'domain':'real coordinates; strict lower-face < upper-face; 0 <= interior order <= 1, supplied by engines',
              'forms':FORMS,'strata':STRATA,'base_points':[[str(x) for x in p] for p in BASE_POINTS],
              'zero_stratum_probe_note':'order profiles have a linear zero at the selected face, remain live and separate argument derivatives there; all printed probe points carry their supplied-domain values',
              'target':'no residual target was given'})
    # Print every live-function probe and the complete domain check at every
    # location, form, stratum and rational point, once, beside the comparisons.
    for form in FORMS if probe_metadata else ():
        for stratum in STRATA:
            emit(out,{'kind':'test_functions','form':form,'stratum':stratum,
                      'formal_arguments':[str(x) for x in FORMAL],
                      'functions':{r.py:str(profile(r.py,r.arity,form,stratum)) for r in names if r.arity}})
            for index in range(3):
                for loc in ('below','inside','above','lower','upper'):
                    env=point_for(form,index,loc)
                    checks={'face_order':evaluate(s.Lt(WM,WP,evaluate=False),form,stratum,index,loc),
                            'f_range':evaluate(s.And(s.Ge(s.Function('f',real=True)(*COORDS),0),
                                                    s.Le(s.Function('f',real=True)(*COORDS),1)),form,stratum,index,loc),
                            'G_range':evaluate(s.And(s.Ge(s.Function('G',real=True)(W-WM,WP-W),0),
                                                    s.Le(s.Function('G',real=True)(W-WM,WP-W),1)),form,stratum,index,loc)}
                    emit(out,{'kind':'probe_domain','form':form,'stratum':stratum,'point_index':index,'location':loc,
                              'point':{str(k):str(v) for k,v in env.items()},
                              'supplied_domain_values':{k:str(v) for k,v in checks.items()}})
    summary=Counter()
    for tag in sorted(set(streams['py'].index)|set(streams['wl'].index)):
        if progress:
            print('WORKING_ON: '+tag,file=sys.stderr,flush=True)
        entries={e:st.index.get(tag,[]) for e,st in streams.items()}
        counts={e:len(v) for e,v in entries.items()}
        if any(v!=1 for v in counts.values()):
            reason='repeated tag; unique join rejected' if any(v>1 for v in counts.values()) else 'tag absent in other engine; no differently named pairing'
            for e,values in entries.items():
                for i,entry in enumerate(values):
                    emit(out,{'kind':'accounting','tag':tag,'accounting':'unjoined','engine':e,'occurrence':i,
                              'reason':reason,'line':entry.line,'parsed_leaves':census(entry.node),
                              'original_operand':entry.raw,'parse':data(entry.node)})
                    summary['unjoined_'+e]+=1
            continue
        p,w=entries['py'][0],entries['wl'][0]
        # Original tagged operands precede extraction, any guard or residual.
        emit(out,{'kind':'operands','tag':tag,'py':{'tag':p.tag,'line':p.line,'raw':p.raw,'parse':data(p.node)},
                  'wl':{'tag':w.tag,'line':w.line,'raw':w.raw,'parse':data(w.node)}})
        c=correspondence(tag,correspondences)
        cmp.compared=[set(),set()]
        if c is None:
            result={'not_compared':'no source-cited component correspondence for this joined name'}
        else:
            emit(out,{'kind':'component_correspondence',**asdict(c)})
            ACTIVE_OUTPUT=out
            result=cmp.compare(p.node,w.node,tag,c,streams)
            ACTIVE_OUTPUT=None
        if tag.endswith('/coefficient_operand'):
            result['role']='unrestricted coefficient representation diagnostic; not a delta part'
        emit(out,{'kind':'comparison','tag':tag,'comparison':result})
        parsed={'py':census(p.node),'wl':census(w.node)}
        actual={e:len(cmp.compared[i]) for i,e in enumerate(('py','wl'))}
        missing={e:sorted([list(x) for x in leaf_paths(entry.node)-cmp.compared[i]],key=str)
                 for i,(e,entry) in enumerate((('py',p),('wl',w)))}
        emit(out,{'kind':'accounting','tag':tag,'accounting':'joined','parsed_leaves':parsed,
                  'actually_compared_leaves':actual,'uncompared_leaf_paths':missing})
        summary['joined']+=1
        summary['joined_with_not_compared']+=has_key(result,'not_compared')
        summary['joined_with_not_formed']+=has_key(result,'not_formed')
        summary['joined_below_parsed_count']+=any(actual[e]<parsed[e] for e in parsed)
    emit(out,{'kind':'resources','runtime_seconds':time.monotonic()-started,
              'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    emit(out,{'kind':'summary','counts':dict(summary),'exit_semantics':'zero for any findings; nonzero only for operational failure'})
    return dict(summary)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--py',type=Path,required=True)
    ap.add_argument('--wl',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    try:
        with args.output.open('x') as out:
            run(args.py,args.wl,out)
    except (InputError,OSError) as exc:
        print('operational_error: '+str(exc),file=sys.stderr,flush=True)
        return 2
    return 0


if __name__=='__main__':
    sys.exit(main())


