#!/usr/bin/env python3
"""S16 native-stream comparator. No residual target was given.

The S12/S9b lossless constructor and InputForm readers and S12 partial-slot
translation are adapted here. Neither accepted engine is imported or run.
Bindings name objects, never outcomes. JSONL contains the original operands,
exact residuals before domain restrictions, rational probes and parse census.
Parser and partial-slot machinery adapted from de81b9fc; S16 contracts, integral
evaluation and probe profiles are new. No accepted engine is imported or run.
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
            if any(e.head != 'Rule' for e in entries):
                raise InputError('association requires rules')
            left = Node('Dictionary', tuple(Node('Pair', e.args) for e in entries))
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
            if head in ('StrictLessThan','LessThan','StrictGreaterThan','GreaterThan') and left.head in ('StrictLessThan','LessThan','StrictGreaterThan','GreaterThan'):
                left = Node('And', (left, Node(head, (left.args[-1], right))))
                continue
            # Flatten associative syntax, without evaluating operands.
            args = (*left.args, right) if head in ('Add', 'Mul', 'Alternatives') and left.head == head else (left, right)
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


# Source citations below refer to accepted 1d990ee6, without importing engines.
from itertools import product
from io import StringIO

@dataclass(frozen=True)
class NameRow:
    py: str
    wl: str
    object: str
    py_line: str
    wl_line: str
    arity: int = 0
    alias: str = ''

PARTS = (('bulk.faceMinus','bulkSource1'), ('bulk.facePlus','bulkSource2'),
         ('bulk.reference.centerMotion','bulkReference.centerMotion'),
         ('bulk.reference.radiusMotion','bulkReference.radiusMotion'),
         ('bulk.weight','bulkWeight'),
         *((f'surface{k}.{o}',f'surface{k+1}.{o}') for k in range(2)
           for o in ('rest','centerMotion','radiusMotion')))
SURFACE_TABLE=({'object':'inner sphere','py_index':0,'wl_index':1,'py_line':'39,51-54,108-110','wl_line':'21,38,47-53'},
               {'object':'outer sphere','py_index':1,'wl_index':2,'py_line':'39,51-54,108-110','wl_line':'21,38,47-53'})
FAMILIES = (('Nrate','nrate'), *((f'Drate{i}',f'drate{i+1}') for i in range(3)))

def name_contract():
    rows=[]
    def add(p,w,obj,pl,wl,arity=0,alias=''):
        rows.append(NameRow(p,w,obj,pl,wl,arity,alias))
    for n in ('x1','x2','x3','t','rho','theta','phi'):
        add(n,n,'region coordinate '+n,'33-36','18-19,37-43')
    for i in range(1,4):
        add(f'Xc{i}',f'c{i}',f'region centre component {i}','38','20',1)
        add(f'D{i}',f'dip[{i}]',f'item 1 dipole component {i}','102','121',1)
        add(f'P{i}',f'mom[{i}]',f'item 1 momentum component {i}','103','122',1)
    for p,w in (('r_in','rin'),('R','rout')):
        add(p,w,'region inner radius' if p=='r_in' else 'region outer radius','39','21',1)
    add('N','count','item 1 material count','101','120',1)
    add('Sigma','Sigma','supplied material density','40','22',4)
    for stem,pl,wl,obj in (('F','41','23','material flux'),('G','42','24','momentum density'),('f','44','25','force density')):
        for i in range(1,4): add(f'{stem}{i}',f'{stem}{i}',f'supplied {obj} component {i}',pl,wl,4)
    for i,j in product(range(1,4),repeat=2):
        add(f'Pi{i}{j}',f'Pi{i}{j}',f'supplied momentum flux ordered indices ({i},{j})','43','26-27',4)
    for p,w in (('n1','n1'),('n3','n3'),('a_minus','aminus'),('a_plus','aplus')):
        add(p,w,'supplied material face quantity '+p,'45','28',4)
    # q has two explicitly printed arities, for the volume and surface templates.
    add('q','q','integral definition test integrand','105,109','52-53',-1)
    for p,w in (('unitLength','"length"'),('unitTime','"time"'),('unitCount','"count"'),('unitMomentum','"momentum"')):
        add(p,w,'item 5 dimension basis '+p,'342,356','205,217')
    for pf,wf in FAMILIES:
        for pp,wp in PARTS:
            add(pf+'_'+pp.replace('.','_'),wf+wp.replace('.',''),
                'item 2 named rate term '+pf+'.'+pp,'309-313,320-321','178-184',1,pf+'.term.'+pp)
    return tuple(rows)

NAME_TABLE=name_contract()
ALIAS_NAMES={r.py for r in NAME_TABLE if r.alias}

@dataclass(frozen=True)
class End:
    tag: str
    path: tuple=()

@dataclass(frozen=True)
class Binding:
    key: str
    rule: str
    py: End
    wl: End
    object: str
    py_line: str
    wl_line: str
    kind: str='value'
    normalized: bool=False

# Shapes are semantic layouts fixed from cited source constructors, never
# inferred by comparing two payloads. Keys name dictionary entries, not values.
ROOT_SHAPES={}
RULES={}

def binding_contract():
    rows=[]
    def add(rule,p,w,obj,pl,wl,kind='value',pp=(),wp=(),normalized=False):
        citation=(obj.split(':')[0],pl,wl)
        if rule not in RULES: RULES[rule]=citation
        rows.append(Binding(f'{rule}:{len(rows)}',rule,End(p,pp),End(w,wp),obj,pl,wl,kind,normalized))
    def group(rule,p,w,shape,obj,pl,wl,kind='value',normalized=False):
        ROOT_SHAPES['py',p]=shape; ROOT_SHAPES['wl',w]=shape
        for path in shape_paths(shape):
            add(rule,p,w,obj+': '+str(path),pl,wl,kind,path,path,normalized)
    for mode in ('unnormalized','normalized'):
        add('domains','domain.'+mode,'domain.'+mode,'region and item 2 domain: '+mode,'100,104','40,50-51','condition')
    for p,w in (('Volume','V'),('Surface.0','B.1'),('Surface.1','B.2')):
        add('integral_definitions','integral.'+p,'integral.'+w,'region integral definition: '+p,'105-110','47-53')
    add('supplied_material','supplied.material','supplied.material','supplied material equation','49,111','30,54')
    group('supplied_momentum','supplied.momentum','supplied.momentum',('v',)*3,'supplied momentum equation','50,112','31,55')
    for p,w,shape in (('surfaces','surfaces.chart',(('v',)*3,)*2),('normal','surfaces.normal',(('v',)*3,)*2),
                      ('areaMeasure','surfaces.measure',('v',)*2),('volumeMeasure','volume.measure','v')):
        group('geometry_'+p,'geometry.'+p,w,shape,'region '+p,'51-54,64-68,74,114-117','37-46,56-58')
    indices=[(),*((i,) for i in range(3)),*product(range(3),repeat=2)]
    for k,idx in enumerate(indices):
        add('moments','moment.'+str(idx),'moment.'+str(k),'item 4 moment: '+str(idx),'56-57,197-200','116-119')
    defs=[('count.definition','count.definition')]+[(f'{obj}.definition.{i}',f'{obj}.definition.{i+1}') for obj in ('dipole','momentum') for i in range(3)]
    for p,w in defs:
        group('definitions',p,w,{'lhs':'v','rhs':'v'},'item 1 definition: '+p,'101-103,197,201-204','120-122')
    for pf,wf,is_momentum in [('count','count',False)]+[(f'{obj}.{i}',f'{obj}.{i+1}',obj=='momentum') for obj in ('dipole','momentum') for i in range(3)]:
        parts=((('bulk.force','bulkSource1'),)+PARTS[2:]) if is_momentum else PARTS
        for pp,wp in (*parts,('rate','rate')):
            add('rates_'+pp,pf+'.'+pp,wf+'.'+wp,'item 1 rate piece: '+pf+'.'+pp,'159-193,205-209','99-125')
        group('force_companion',pf+'.integratedForceCompanion',wf+'.integratedForceCompanion',('v',)*3,'item 1 integrated force companion','193,205-209','113')
        # The local evolution equation is split-vs-grouped, not a positional join.
        add('local_rate',pf+'.localRate',wf+'.derivation','item 1 solved local rate: '+pf,'161,170,205-209','86,92,101',wp=('localRate',0,'rhs'))
        ROOT_SHAPES['wl',wf+'.derivation']={'localRate':({'lhs':'v','rhs':'v'},),'windowDerivative':'v','afterParts':'v'}
    for p in ('momentumReceived','momentumRate'):
        group(p,p,p,('v',)*3,'item 3 '+p,'211-214','127-129')
    for i in range(3):
        for k in range(2):
            for origin in ('centerMotion','radiusMotion'):
                add('momentum_surface_motion','momentumSurfaceMotion','momentumSurfaceMotion','item 1 momentum surface motion','215-216','130',pp=(i,f'surface{k}.{origin}'),wp=(i,f'surface{k+1}.{origin}'))
    ROOT_SHAPES['py','momentumSurfaceMotion']=tuple({p:'v' for p,_ in PARTS if p.endswith('Motion')} for i in range(3))
    ROOT_SHAPES['wl','momentumSurfaceMotion']=tuple({w:'v' for _,w in PARTS if w.startswith('surface') and w.endswith('Motion')} for i in range(3))
    for pf,wf in FAMILIES:
        for pp,wp in PARTS:
            for sub in ('term','termRate'):
                group('aliases_'+sub,pf+'.'+sub+'.'+pp,wf+'.'+sub+'.'+wp,{'lhs':'v','rhs':'v'},'item 2 '+sub+': '+pf+'.'+pp,'309-315,320-321','178-184')
            origins=[('localEvolution','localEvolution')]
            if pp.startswith('bulk.'):
                origins += [('bulk.centerMotion','bulk.centerMotion'),('bulk.radiusMotion','bulk.radiusMotion')]
                origins += [(f'surface{k}.{o}',f'surface{k+1}.{o}') for k in range(2) for o in ('centerMotion','radiusMotion')]
            else:
                origins += [(pp.split('.')[0]+'.'+o,wp.split('.')[0]+'.'+o) for o in ('centerMotion','radiusMotion')]
            for po,wo in origins:
                add('derivative_piece_'+po.split('.')[-1],pf+'.termRatePart.'+pp+'.'+po,wf+'.termRatePart.'+wp+'.'+wo,
                    'item 2 derivative origin: '+pf+'.'+pp+'.'+po,'243-259,315-317,320-321','147-157,181-183')
            for suffix in ('direct','partsSum'):
                add('partition_'+suffix,pf+'.termRatePartition.'+pp+'.'+suffix,wf+'.termRatePartition.'+wp+'.'+suffix,
                    'item 2 same direct derivative versus same origin sum: '+suffix,'221-238,281-287,315-318','168-173,183')
            for pk,wk in [('outsideIntegrals','outsideIntegrals'),('volumeIntegrand','volumeIntegrand'),('surface0Integrand','surface1Integrand'),('surface1Integrand','surface2Integrand')]:
                add('partition_residual',pf+'.termRatePartition.'+pp+'.residual',wf+'.termRatePartition.'+wp+'.residual',
                    'item 2 partition residual component: '+pk,'266-279,288,318','161-167,172','value',(pk,),(wk,))
            add('partition_routes',pf+'.termRatePartition.'+pp+'.routes',wf+'.termRatePartition.'+wp+'.routes',
                'item 2 partition route labels','289,318','173','label')
    for family in ('relative','absolute'):
        for obj in ('definition','rate','acceleration'):
            p='center.'+family+'.'+obj
            group('center_'+obj,p,p,('v',)*3,'item 2 '+p,'322-333','189-196',normalized=True)
    for obj in ('rateForNormalization','accelerationForNormalization'):
        add('count_normalization','count.'+obj,'count.'+obj,'item 2 unnormalized count '+obj,'322-323,334-335','197-198')
    # Dimensional leaves follow each source's independent dimensional walker.
    # Scalar set-of-units versus set-of-exponent-vectors is one dimension value.
    originals=tuple(rows)
    for row in originals:
        if row.kind!='value' or row.rule in ('domains','integral_definitions','local_rate','partition_residual'):
            continue
        ppath=row.py.path; wpath=row.wl.path
        ptag='.'.join((row.py.tag,*(str(v) for v in ppath),'dimensions'))
        # Py dictionary dimension keys are verbatim, including embedded dots.
        wtag=row.wl.tag+'.dimensions'
        add('dimensions',ptag,wtag,'item 5 dimension of '+row.object,'359-407,425-428','230-243','dimension',(),wpath)
    # Labels stay text; no paraphrase changes their content.
    for p,w in [('Volume','V'),('Surface','B'),('X','X'),('premises','physicalPremises'),('shapeChoice','shapeChoice'),('momentOrderChoice','momentOrderChoice'),('f','forceScope'),('Pi','PiScope')]:
        add('schema','schema','schema','schema entry: '+p,'92-99','49','label',(p,),(w,))
    add('dimension_basis','dimensions.basis','dimensions.basis','item 5 dimension basis','342,356','217','basis')
    for n in NAME_TABLE:
        if n.arity and n.py!='q':
            if n.alias: continue  # WL prints fields before aliases are inserted (217 vs 234).
            add('dimension_fields','dimensions.fields','dimensions.fields','item 5 field dimension: '+n.object,'344-355,357','207-217','dimension_field',(n.py,),(n.wl,))
    for p in ('objects','ungovernedRates','fluxDecomposition'):
        add('coverage','coverage.'+p,'coverage.'+p,'coverage label: '+p,'435-437','244-246','label')
    # Comparison notes only for two fully corresponding emitted objects.
    roots={}
    for r in originals:
        if r.rule in ('local_rate','partition_residual','partition_routes','partition_direct','partition_partsSum','integral_definitions','domains'): continue
        roots[r.py.tag,r.wl.tag]=r
    for (p,w),r in roots.items():
        add('comparison_notes',p+'.comparison',w+'.comparison','engine comparison note for '+r.object,'23-26','15','label')
    return tuple(rows)

def shape_paths(shape,path=()):
    if shape=='v': yield path
    elif isinstance(shape,dict):
        for k,v in shape.items(): yield from shape_paths(v,path+(k,))
    else:
        for i,v in enumerate(shape): yield from shape_paths(v,path+(i,))

BINDING_TABLE=binding_contract()

def validate_tables(names,bindings,rules=None):
    for field in ('py','wl'):
        vals=[getattr(n,field) for n in names]
        if len(set(vals))!=len(vals): raise InputError('duplicate name row '+field)
        seen=defaultdict(list)
        for row in bindings:
            end=getattr(row,field)
            if field=='wl' and row.kind=='dimension':
                end=End(end.tag,tuple({'lhs':0,'rhs':1}.get(k,k) for k in end.path))
            for old in seen[end.tag]:
                if old==end.path[:len(old)] or end.path==old[:len(end.path)]:
                    raise InputError('duplicate binding end '+field+' '+str(end))
            seen[end.tag].append(end.path)
    keys=[r.key for r in bindings]
    if len(set(keys))!=len(keys): raise InputError('duplicate binding rule instance')
    if rules is not None and len(set(rules))!=len(rules): raise InputError('duplicate binding rule')


def call(n,name):
    return n.head=='Apply' and n.args[0].head=='Symbol' and n.args[0].value==name

def application(n):
    if n.head!='Apply': return None,()
    head=n.args[0]
    if head.head=='Symbol': return head.value,n.args[1:]
    if call(head,'Inactive'): return head.args[1].value,n.args[1:]
    return None,n.args[1:]

def census(n):
    return sum(census(a) for a in n.args) if n.args else 1

def leaf_paths(n,p=()):
    return set().union(*(leaf_paths(a,p+(i,)) for i,a in enumerate(n.args))) if n.args else {p}

def data(n): return [n.head,n.value,[data(a) for a in n.args]]

def key_text(n):
    if n.head in ('Symbol','Text','Number'): return n.value
    if n.head=='Apply': return key_text(n.args[0])+'['+','.join(key_text(q) for q in n.args[1:])+']'
    raise InputError('nonliteral dictionary key')

# Native sstr removes string quotes, including dictionary values and dot keys.
# These two readers are deliberately limited to engine-emitted label containers.
def py_dictionary(text,labels=False):
    body=text.strip()[1:-1]
    marks=list(re.finditer(r'(?:^|, )([A-Za-z][A-Za-z0-9_.]*): ',body))
    if not marks and body: raise InputError('sstr dictionary key grammar')
    pairs=[]
    for i,m in enumerate(marks):
        val=body[m.end():marks[i+1].start() if i+1<len(marks) else len(body)]
        pairs.append(Node('Pair',(atom('Text',m[1]),atom('Text',val) if labels else parse_payload(val,'py',''))))
    return Node('Dictionary',tuple(pairs))

def parse_payload(text,engine,tag):
    if engine=='wl': return wl_parse(text)
    if tag.endswith('.comparison') or tag=='coverage.fluxDecomposition' or tag=='WORKING_ON': return atom('Text',text)
    if tag=='schema': return py_dictionary(text,True)
    if tag in ('coverage.objects','coverage.ungovernedRates'):
        if not (text.startswith('[') and text.endswith(']')): raise InputError('coverage list grammar')
        # Object names may contain tuple commas. Their delimiter is unambiguous
        # because the engine's object names start with an alphabetic character.
        return seq(*(atom('Text',v) for v in re.split(r', (?=[A-Za-z])',text[1:-1])))
    if tag.endswith('.routes'):
        if not (text.startswith('(') and text.endswith(')')): raise InputError('route tuple grammar')
        return seq(*(atom('Text',v) for v in text[1:-1].split(', ')))
    if text.startswith('{') and ':' in text: return py_dictionary(text)
    if text.startswith('[{'):
        return seq(*(py_dictionary(v) for v in re.findall(r'\{[^{}]*\}',text)))
    return py_parse(text)

@dataclass
class Entry:
    tag: str
    raw: str
    line: int
    node: Node

class Stream:
    def __init__(self,path,engine):
        self.path=Path(path); self.engine=engine; self.index=defaultdict(list); self.progress=[]
        try: raw=self.path.read_bytes(); lines=raw.decode().splitlines()
        except (OSError,UnicodeError) as exc: raise InputError(str(exc)) from exc
        self.sha256=hashlib.sha256(raw).hexdigest(); self.line_count=len(lines)
        framed=[]
        for i,line in enumerate(lines,1):
            match=re.fullmatch(r'([A-Za-z][^=\n]*?) = (.*)',line)
            if match: framed.append([match[1],match[2],i])
            elif framed and framed[-1][0]!='WORKING_ON' and line.startswith((' ','[',']',')','}')):
                framed[-1][1]+='\n'+line
            else: raise InputError(f'{path}:{i}: unassigned line {line[:80]!r}')
        if not framed: raise InputError('empty stream')
        for tag,txt,line in framed:
            try: node=parse_payload(txt,engine,tag)
            except InputError as exc: raise InputError(f'{path}:{line}:{tag}: {exc}') from exc
            e=Entry(tag,txt,line,node)
            (self.progress if tag=='WORKING_ON' else self.index[tag]).append(e)

@dataclass(frozen=True)
class Ref:
    node: Node
    path: tuple=()
    def child(self,k):
        try: return Ref(self.node.args[k],self.path+(k,))
        except IndexError as exc: raise Unsupported('component missing '+str(self.path+(k,))) from exc

def semantic_children(ref):
    n=ref.node
    if n.head=='Sequence': return {i:ref.child(i) for i in range(len(n.args))}
    if n.head=='Dictionary':
        result={}
        for i,p in enumerate(n.args):
            k=key_text(p.args[0])
            if k in result: raise InputError('duplicate payload key '+k)
            result[k]=ref.child(i).child(1)
        return result
    if call(n,'Matrix'):
        rows=ref.child(1)
        if any(x.head!='Sequence' or len(x.args)!=1 for x in rows.node.args): raise Unsupported('column matrix layout')
        return {i:rows.child(i).child(0) for i in range(len(rows.node.args))}
    if n.head=='Rule': return {'lhs':ref.child(0),'rhs':ref.child(1)}
    if call(n,'Eq') or call(n,'Equal'): return {'lhs':ref.child(1),'rhs':ref.child(2)}
    return None

def get_ref(ref,path):
    for k in path:
        children=semantic_children(ref)
        if children is None or k not in children: raise Unsupported('no semantic component '+str(k)+' at '+str(ref.path))
        ref=children[k]
    return ref

def shaped(ref,shape):
    children=semantic_children(ref)
    if shape=='v': return children is None
    expected=shape if isinstance(shape,dict) else dict(enumerate(shape))
    return children is not None and set(children)==set(expected) and all(shaped(children[k],v) for k,v in expected.items())

def shaped_dimension(ref,shape):
    if shape=='v':
        return ref.node.head=='Sequence' and all(q.head=='Sequence' and len(q.args)==4 for q in ref.node.args)
    if isinstance(shape,dict) and set(shape)=={'lhs','rhs'}:
        shape=(shape['lhs'],shape['rhs'])
    expected=shape if isinstance(shape,dict) else dict(enumerate(shape))
    children=semantic_children(ref)
    return children is not None and set(children)==set(expected) and all(shaped_dimension(children[k],v) for k,v in expected.items())


RELATIONS={'Equality':s.Eq,'Unequality':s.Ne,'StrictLessThan':s.Lt,'StrictGreaterThan':s.Gt,'LessThan':s.Le,'GreaterThan':s.Ge}
BUILTINS={'Eq':'Equality','Equal':'Equality','Ne':'Unequality','Unequal':'Unequality','Lt':'StrictLessThan','Less':'StrictLessThan',
          'Le':'LessThan','LessEqual':'LessThan','Gt':'StrictGreaterThan','Greater':'StrictGreaterThan','Ge':'GreaterThan',
          'GreaterEqual':'GreaterThan','HeavisideTheta':'Heaviside','Sin':'sin','Cos':'cos','Exp':'exp','Sqrt':'sqrt',
          'Integrate':'Integral','D':'Derivative'}


def partial_signature(expr):
    if isinstance(expr,s.Subs):
        sig=partial_signature(expr.expr)
        if sig is None:
            return None
        head,orders,points=sig
        return head,orders,tuple(s.Subs(p,expr.variables,expr.point).doit() for p in points)
    if not isinstance(expr,s.Derivative) or not isinstance(expr.expr,AppliedUndef):
        return None
    fun=expr.expr
    counts=dict(expr.variable_count)
    if any(not isinstance(v,s.Symbol) or fun.args.count(v)!=1 or any(p!=v and p.has(v) for p in fun.args) for v in counts):
        return None
    return fun.func,tuple(int(counts.get(p,0)) for p in fun.args),fun.args


def partial_derivative(head,orders,points):
    if not any(orders):
        return head(*points)
    direct=all(isinstance(p,s.Symbol) for p in points) and len(set(points))==len(points)
    slots=points if direct else s.symbols('_s14partial0:'+str(len(points)),**{})
    derivative=s.Derivative(head(*slots),*((v,k) for v,k in zip(slots,orders) if k),evaluate=False)
    return derivative if direct else s.Subs(derivative,slots,points)


def normalize_partial(expr):
    sig=partial_signature(expr)
    if sig is not None:
        return partial_derivative(*sig)
    if isinstance(expr,s.Derivative):
        expanded=expr.doit()
        if expanded!=expr:
            return expanded.replace(lambda q:isinstance(q,(s.Derivative,s.Subs)),normalize_partial)
    return expr


def native_relation(head,*args):
    if len(args)>2:
        return s.And(*(native_relation(head,a,b) for a,b in zip(args,args[1:])),evaluate=False)
    left,right=args
    if head in ('StrictGreaterThan','GreaterThan'):
        head={'StrictGreaterThan':'StrictLessThan','GreaterThan':'LessThan'}[head]
        left,right=right,left
    if head in ('Equality','Unequality') and s.default_sort_key(right)<s.default_sort_key(left):
        left,right=right,left
    return RELATIONS[head](left,right,evaluate=False)


def alpha_integrals(expr):
    """Name integral binders by active lexical depth, without capture.

    Conversion first uses fresh Dummies. Only the completed expression is
    renamed: an inner integral cannot give an outer variable its own name.
    Flattened multi-limit integrals and nested integral syntax share the
    same outside-in binder order. Sibling scopes may reuse names safely.
    """
    def walk(value,env,depth):
        if isinstance(value,s.Integral):
            local=dict(env)
            for i,var in enumerate(reversed(value.variables)):
                local[var]=s.Symbol('_int'+str(depth+i),**{})
            inner_depth=depth+len(value.limits)
            limits=[(local[var],walk(lo,local,inner_depth),walk(hi,local,inner_depth))
                    for var,lo,hi in value.limits]
            return s.Integral(walk(value.function,local,inner_depth),*limits)
        if not value.args:
            return env.get(value,value)
        args=tuple(walk(arg,env,depth) for arg in value.args)
        if args==value.args:
            return value
        if isinstance(value,(s.Derivative,Relational)):
            return value.func(*args,evaluate=False)
        return value.func(*args)
    return walk(expr,{},0)


def math_node(n,b,engine,scope=None):
    return _math_node(n,b,engine,scope)


def _math_node(n,b,engine,scope=None):
    scope={} if scope is None else scope
    def cv(x,local=None):
        return _math_node(x,b,engine,scope if local is None else local)
    h,a=n.head,n.args
    if h=='Number':
        if n.value=='None':
            return s.Symbol(b.name('None',engine),**{})
        return s.Rational(n.value)
    if h=='Text' and ('"'+n.value+'"') in b.maps[engine]:
        return s.Symbol(b.name('"'+n.value+'"',engine))
    if h=='Boolean':
        return s.true if n.value.lower()=='true' else s.false
    if h=='Symbol':
        if n.value in scope:
            return scope[n.value]
        if n.value in ('pi','Pi'):
            return s.pi
        if n.value in ('True','true','False','false'):
            return s.true if n.value.lower()=='true' else s.false
        key=b.name(n.value,engine)
        if b.arities.get(key,0):
            raise Unsupported('unapplied live function '+key)
        return s.Symbol(key,**{})
    if h=='Apply':
        head,*args=a
        if head.head=='Apply' and call(head.args[0],'Derivative'):
            orders=head.args[0].args[1:]
            if len(head.args)!=2:
                raise Unsupported('derivative head structure')
            key=b.name(key_text(head.args[1]),engine)
            if b.arities.get(key)!=len(args) or len(orders)!=len(args):
                raise Unsupported('derivative function arity '+key)
            return partial_derivative(s.Function(key,**{}),tuple(int(cv(q)) for q in orders),tuple(cv(q) for q in args))
        name,args=application(n)
        if name is None and head.head=='Apply':
            name=key_text(head)
        if name is None:
            raise Unsupported('composite application head')
        if name in b.maps[engine] and b.arities.get(b.name(name,engine),0):
            key=b.name(name,engine)
            if b.arities[key]!=-1 and len(args)!=b.arities[key]:
                raise Unsupported('live function arity '+key)
            return s.Function(key,**{})(*(cv(q) for q in args))
        if name in ('Volume','Surface','V','B'):
            return b.integral(engine,name,tuple(cv(q) for q in args))
        h,a=BUILTINS.get(name,name),args
    if h in ('Add','Mul','Pow'):
        return getattr(s,h)(*(cv(q) for q in a))
    if h=='Rule':
        return s.Tuple(*(cv(q) for q in a))
    if h=='Element':
        members=a[0].args if a[0].head in ('Alternatives','Sequence') else (a[0],)
        if a[1].value!='Reals': raise Unsupported('domain set '+a[1].value)
        return s.And(*(s.Contains(cv(q),s.S.Reals,evaluate=False) for q in members),evaluate=False)
    if h in RELATIONS:
        return native_relation(h,*(cv(q) for q in a))
    if h in ('And','Or','Not'):
        return getattr(s,h)(*(cv(q) for q in a),evaluate=False)
    if h=='Sequence':
        return s.Tuple(*(cv(q) for q in a))
    if h in ('Heaviside','DiracDelta','Abs','sin','cos','exp','sqrt','Rational'):
        return getattr(s,h)(*(cv(q) for q in a))
    if h=='Integral':
        if len(a)<2:
            raise Unsupported('integral missing binder')
        # Multi-limit integrals bind inside out. Bounds of an inner integral
        # may use its outer variable; neither name is a free coordinate.
        local=dict(scope)
        variables=[]
        for spec in a[1:]:
            if spec.head!='Sequence' or len(spec.args)!=3 or spec.args[0].head!='Symbol':
                raise Unsupported('definite integral binder grammar')
            var=s.Symbol(b.name(spec.args[0].value,engine))
            local[spec.args[0].value]=var
            variables.append(var)
        limits=[]
        for spec,var in zip(a[1:],variables):
            bounds_scope=dict(local)
            bounds_scope.pop(spec.args[0].value,None)
            bounds_scope.update({k:v for k,v in scope.items() if k==spec.args[0].value})
            limits.append((var,cv(spec.args[1],bounds_scope),cv(spec.args[2],bounds_scope)))
        return s.Integral(cv(a[0],local),*limits)
    if h=='Derivative':
        specs=[]
        for q in a[1:]:
            specs.append((cv(q.args[0]),int(cv(q.args[1]))) if q.head=='Sequence' else cv(q))
        return normalize_partial(s.Derivative(cv(a[0]),*specs,evaluate=False))
    if h=='Subs':
        def members(q):
            return q.args if q.head=='Sequence' else (q,)
        variables,points=members(a[1]),members(a[2])
        if len(variables)!=len(points) or any(v.head!='Symbol' for v in variables):
            raise Unsupported('substitution arity or variable')
        local=dict(scope)
        for i,var in enumerate(variables):
            local[var.value]=s.Symbol('_slot'+str(len(scope)+i),**{})
        return normalize_partial(s.Subs(cv(a[0],local),tuple(local[v.value] for v in variables),tuple(cv(q) for q in points)))
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
    raise Unsupported('mathematical head '+h)


class Names:
    def __init__(self,rows):
        self.rows=rows
        self.maps={'py':{r.py:r.py for r in rows},'wl':{r.wl:r.py for r in rows}}
        self.arities={r.py:r.arity for r in rows}
        self.templates={}; self.aliases={'py':{},'wl':{}}
    def name(self,name,engine): return self.maps[engine].get(name,engine+'::'+name)
    def integral(self,engine,head,args):
        key='volume' if head in ('Volume','V') else ('inner' if int(args[0])==(0 if engine=='py' else 1) else 'outer')
        if head in ('Surface','B') and int(args[0]) not in ((0,1) if engine=='py' else (1,2)):
            raise Unsupported('unknown surface '+str(args[0]))
        if (engine,key) not in self.templates: raise Unsupported('missing printed integral definition '+engine+' '+key)
        template=self.templates[engine,key]
        qcalls=[f for f in template.atoms(AppliedUndef) if f.func.__name__==self.name('q',engine)]
        if len(qcalls)!=1: raise Unsupported('integral definition requires one live q call')
        q=qcalls[0]; integrand=args[-1]
        if key=='volume':
            coords=tuple(s.Symbol(self.name(n,engine)) for n in ('x1','x2','x3','t'))
            # x is bound Cartesian position (Py 82; WL 49 integral schema).
            # Sequential substitution retains Subs around partial derivatives;
            # simultaneous=True instead differentiates the time-dependent chart.
            # The chart contains no Cartesian x, so replacements cannot cascade.
            integrand=integrand.subs(list(zip(coords,q.args)))
        return template.xreplace({q:integrand})
    def load(self,streams):
        for engine,stream in streams.items():
            tags=(('integral.Volume','volume'),('integral.Surface.0','inner'),('integral.Surface.1','outer')) if engine=='py' else (('integral.V','volume'),('integral.B.1','inner'),('integral.B.2','outer'))
            for tag,key in tags:
                entries=stream.index.get(tag,[])
                if len(entries)==1: self.templates[engine,key]=math_node(entries[0].node,self,engine)
            for row in self.rows:
                if not row.alias: continue
                tag=row.alias
                if engine=='wl':
                    pf,rest=tag.split('.term.',1)
                    tag=dict(FAMILIES)[pf]+'.term.'+dict(PARTS)[rest]
                for kind in ('term','termRate'):
                    entries=stream.index.get(tag.replace('.term.','.'+kind+'.'),[])
                    if len(entries)!=1: continue
                    refs=semantic_children(Ref(entries[0].node))
                    if refs is None or set(refs)!={'lhs','rhs'}: raise InputError('alias definition structure '+tag)
                    left=math_node(refs['lhs'].node,Names(NAME_TABLE),engine)
                    right=math_node(refs['rhs'].node,self,engine)
                    if left in self.aliases[engine]: raise InputError('duplicate alias lhs '+str(left))
                    self.aliases[engine][left]=right
    def resolve(self,expr,engine):
        definitions=self.aliases[engine]
        expr=expr.xreplace(definitions)
        if not any(f.func.__name__ in ALIAS_NAMES for f in expr.atoms(AppliedUndef)): return expr
        t=s.Symbol('t')
        def walk(v,active=frozenset()):
            signature=partial_signature(v)
            if signature is not None and signature[0].__name__ in ALIAS_NAMES:
                fun,orders,points=signature
                key=partial_derivative(fun,orders,(t,))
            elif isinstance(v,AppliedUndef) and v.func.__name__ in ALIAS_NAMES:
                fun=v.func;points=v.args;key=fun(t)
            else:
                if not v.args:return v
                args=tuple(walk(q,active) for q in v.args)
                if args==v.args:return v
                return v.func(*args,evaluate=False) if isinstance(v,(s.Derivative,Relational)) else v.func(*args)
            if key not in definitions: raise Unsupported('alias has no printed definition for '+str(key))
            if key in active: raise Unsupported('cyclic printed alias definition '+str(key))
            return walk(definitions[key].subs(t,points[0],simultaneous=True),active|{key})
        return walk(expr)

FORMS=('mixed_quadratic','mixed_cubic','mixed_quartic')
STRATA=('generic','center_velocity_zero','radius_velocity_zero','all_velocities_zero','count_zero')
POINTS=((s.Rational(1,5),s.Rational(2,7),s.Rational(-1,4),s.Rational(1,6)),
        (s.Rational(-2,9),s.Rational(3,8),s.Rational(4,11),s.Rational(-1,7)),
        (s.Rational(3,10),s.Rational(-1,3),s.Rational(2,5),s.Rational(2,9)))
FORMAL=s.symbols('u0:4')
CODES={r.py:i+1 for i,r in enumerate(NAME_TABLE)}

def profile(name,arity,form,stratum,point,options=None):
    options=options or {}; a=FORMAL[:arity]; code=CODES.get(name,101)
    if options.get('duplicate_functions') and name in ('F1','F2'): code=CODES['F1']
    degree=FORMS.index(form)+2
    def coefficient(slot):
        # Independent deterministic coefficients avoid imposing affine relations
        # among successive field/component names. The generated polynomials are
        # printed in full; these hashes contain no measured-stream information.
        raw=hashlib.sha256(f'{code}:{form}:{slot}'.encode()).digest()
        return s.Integer(1+int.from_bytes(raw[:4],'big')%97)
    linear=sum((1000*(j+1)+coefficient(j))*v for j,v in enumerate(a))
    coupling=sum(coefficient(10+i*arity+j)*a[i]*a[j] for i in range(arity) for j in range(i,arity))
    if options.get('symmetric_arguments'):
        linear=code*sum(a);coupling=code*sum(a)**2
    if options.get('uncoupled_arguments'):
        coupling=sum(coefficient(10+j)*v**2 for j,v in enumerate(a))
    result=s.Integer(code)+linear+coupling
    if degree>2:
        result+=sum(coefficient(40 if options.get('symmetric_arguments') else 40+j)*v**degree for j,v in enumerate(a))/7
    if degree==4 and not options.get('uncoupled_arguments'):
        result+=coefficient(45)*s.prod(1+(1 if options.get('symmetric_arguments') else j+1)*v for j,v in enumerate(a))/11
    if arity==1 and (name.startswith('Xc') or name in ('r_in','R','N')):
        z=a[0]-POINTS[point][3]
        slope=(code+coefficient(51)/101)/97
        if name.startswith('Xc') and stratum in ('center_velocity_zero','all_velocities_zero'): slope=s.S.Zero
        if name in ('r_in','R') and stratum in ('radius_velocity_zero','all_velocities_zero'): slope=s.S.Zero
        base=(code+coefficient(50)/101)/103
        if name=='r_in': base=s.Rational(3+point,4)
        elif name=='R': base=s.Rational(11+point,4)
        elif name=='N': base=s.S.Zero if stratum=='count_zero' else s.Rational(code,13)
        result=base+slope*z+(code+1+coefficient(52)/103)*z**2/103+(code+3+coefficient(53)/107)*z**3/107
        if degree==4: result+=(code+5+coefficient(54)/109)*z**4/109
    if name=='Sigma' and arity==4 and stratum=='count_zero':
        # Subtract a constant spatial mean at this stratum's anchor time.
        # All spatial and time derivatives remain those of the live profile.
        result-=zero_count_mean(form,point,tuple(sorted(options.items())))
    return s.expand(result)

@lru_cache(maxsize=128)
def zero_count_mean(form,point,opts=()):
    """Exact mean of the generic density over the concentric spherical shell.

    Odd centred Cartesian moments vanish. For even powers p_i, the angular
    mean is prod((p_i-1)!!)/(sum(p_i)+1)!!; the radial mean supplies the rest.
    This chooses test fields only, and applies no supplied balance equation.
    """
    options=dict(opts); anchor=POINTS[point][3]
    def at(name):
        return profile(name,1,form,'generic',point,options).subs(FORMAL[0],anchor)
    inner,outer=at('r_in'),at('R')
    shifted=profile('Sigma',4,form,'generic',point,options).subs(
        {**{FORMAL[i]:FORMAL[i]+at('Xc'+str(i+1)) for i in range(3)},
         FORMAL[3]:anchor},simultaneous=True)
    mean=s.S.Zero
    for powers,coefficient in s.Poly(shifted,*FORMAL[:3]).terms():
        if any(p%2 for p in powers): continue
        degree=sum(powers)
        angular=s.prod(s.factorial2(p-1) for p in powers)/s.factorial2(degree+1)
        radial=3*(outer**(degree+3)-inner**(degree+3))/((degree+3)*(outer**3-inner**3))
        mean+=coefficient*angular*radial
    return s.cancel(mean)

@lru_cache(maxsize=8192)
def substituted(expr,form,stratum,point,opts=()):
    options=dict(opts)
    def walk(v):
        if isinstance(v,AppliedUndef):
            name=v.func.__name__
            if name in ALIAS_NAMES: raise Unsupported('alias has no resolved printed definition: '+name)
            if name not in CODES: raise Unsupported('no test function for '+name)
            if len(v.args)>4: raise Unsupported('test arity')
            f=profile(name,len(v.args),form,stratum,point,options)
            return f.xreplace(dict(zip(FORMAL,tuple(walk(q) for q in v.args))))
        if isinstance(v,s.Derivative):
            return s.diff(walk(v.expr),*v.variable_count)
        if isinstance(v,s.Subs):
            return walk(v.expr).subs(list(zip(v.variables,tuple(walk(q) for q in v.point))),simultaneous=True)
        if not v.args: return v
        if isinstance(v,s.Integral):
            return s.Integral(walk(v.function),*[(lim[0],*(walk(q) for q in lim[1:])) for lim in v.limits])
        args=tuple(walk(q) for q in v.args)
        if isinstance(v,Relational) or v.func==s.Contains: return v.func(*args,evaluate=False)
        return v.func(*args)
    return walk(expr)

@lru_cache(maxsize=None)
def angular_moment(p,q,full):
    if full:
        if p%2 or q%2: return s.S.Zero
        return 2*s.gamma(s.Rational(p+1,2))*s.gamma(s.Rational(q+1,2))/s.gamma(s.Rational(p+q+2,2))
    if q%2: return s.S.Zero
    return s.gamma(s.Rational(p+1,2))*s.gamma(s.Rational(q+1,2))/s.gamma(s.Rational(p+q+2,2))

NUMERICAL_INTEGRALS={}

def integrate_numeric(integral):
    """Fixed two-precision quadrature; the difference is an error estimate."""
    import mpmath as mp
    variables=tuple(lim[0] for lim in integral.limits)
    evaluator=s.lambdify(variables,integral.function,'mpmath')
    values=[]
    try:
        for precision in (40,60):
            with mp.workdps(precision):
                def walk(level,env):
                    if level<0: return evaluator(*(env[v] for v in variables))
                    var,lo,hi=integral.limits[level]
                    def number(v):
                        replacements={k:s.Float(str(z),precision) for k,z in env.items()}
                        return mp.mpf(str(s.N(v.subs(replacements),precision)))
                    return mp.quad(lambda z:walk(level-1,{**env,var:z}),[number(lo),number(hi)])
                v=walk(len(variables)-1,{})
                if not mp.isfinite(v): raise Unsupported('nonfinite quadrature')
                values.append(mp.nstr(v,precision))
        with mp.workdps(60): error=mp.nstr(abs(mp.mpf(values[1])-mp.mpf(values[0])),60)
    except (ValueError,TypeError,ZeroDivisionError) as exc:
        raise Unsupported('quadrature could not evaluate integral: '+str(exc)) from exc
    key=hashlib.sha256(str(integral).encode()).hexdigest()
    NUMERICAL_INTEGRALS[key]={'integral':str(integral),'working_precisions':[40,60],
        'values':values,'error_estimate':error,'estimate_kind':'absolute difference of two precision evaluations, not a rigorous bound'}
    return s.Float(values[1],60)

@lru_cache(maxsize=16384)
def integrate_exact(expr):
    if not isinstance(expr,s.Integral): return expr
    f=expr.function
    for var,lo,hi in expr.limits:
        # Polynomial Cartesian test functions produce finite trigonometric
        # polynomials. Integrate their monomials with beta integrals, exactly.
        if (lo==0 and hi in (s.pi,2*s.pi)):
            sn,cs=s.Dummy('sin'),s.Dummy('cos')
            polyexpr=s.cancel(f).subs({s.sin(var):sn,s.cos(var):cs})
            if not polyexpr.has(var):
                try: pol=s.Poly(polyexpr,sn,cs)
                except s.PolynomialError: pol=None
                if pol is not None:
                    f=s.Add(*(coef*angular_moment(p,q,hi==2*s.pi) for (p,q),coef in pol.terms()))
                    continue
        try:
            poly=s.Poly(f,var)
        except s.PolynomialError: poly=None
        if poly is not None:
            f=s.Add(*(coef*(hi**(power[0]+1)-lo**(power[0]+1))/(power[0]+1) for power,coef in poly.terms()))
        else:
            f=s.integrate(f,(var,lo,hi))
            if f.has(s.Integral): return integrate_numeric(expr)
    return s.cancel(f)

@lru_cache(maxsize=32768)
def evaluate_cached(expr,form,stratum,index,opts=()):
    # Evaluate integral values compositionally. Derivative and Subs nodes are
    # deliberately excluded: their operations precede evaluation at a point.
    if expr.has(s.Integral) and expr.func in (s.Add,s.Mul,s.Pow):
        values=tuple(evaluate_cached(a,form,stratum,index,opts) for a in expr.args)
        return s.simplify(expr.func(*values))
    value=substituted(expr,form,stratum,index,opts)
    env=dict(zip(s.symbols('x1 x2 x3 t'),POINTS[index]))
    env.update({s.Symbol('rho'):s.Rational(3,2),s.Symbol('theta'):s.Rational(2+index,5),s.Symbol('phi'):s.Rational(1+index,3)})
    env.update(dict(zip(s.symbols('unitLength unitTime unitCount unitMomentum'),map(s.Integer,(2,3,5,7)))))
    # Substitute free coordinates only. Integral variables remain bound.
    value=value.subs(env,simultaneous=True)
    if value.has(s.Integral):
        value=value.replace(lambda q:isinstance(q,s.Integral),integrate_exact)
    value=s.simplify(value)
    if value.has(AppliedUndef,s.Derivative,s.Subs,s.Integral): raise Unsupported('unevaluated live object '+str(value))
    if value.free_symbols: raise Unsupported('unevaluated symbols '+str(value.free_symbols))
    return value

def exact_difference(a,b):
    if a in (s.true,s.false) or isinstance(a,Relational) or a.func in (s.And,s.Or,s.Not,s.Xor,s.Contains):
        return s.Xor(a,b,evaluate=False)
    # No positivity assumption, supplied balance, face law or domain guard.
    residual=s.expand(a-b)
    groups=defaultdict(list)
    for integral in residual.atoms(s.Integral):
        coefficient=residual.coeff(integral)
        if coefficient and not coefficient.has(s.Integral) and not (coefficient.free_symbols & set(integral.variables)):
            residual-=coefficient*integral
            groups[integral.limits].append(coefficient*integral.function)
    for limits,terms in groups.items():
        integrand=s.expand(s.Add(*terms))
        if not integrand.has(AppliedUndef,s.Derivative,s.Subs): integrand=s.trigsimp(integrand)
        if integrand!=0: residual+=s.Integral(integrand,*limits)
    return s.expand(residual)

def resolved_exact_difference(residual):
    """Evaluate canonical substitution wrappers without changing derivative jets.

    A partial derivative is rebuilt from its slot orders and evaluation point,
    never by inserting the point into its differentiated function. Other Subs
    nodes use SymPy's substitution semantics, after resolving their children.
    Integrals remain integrals; only zero integrands and linearity are reduced.
    """
    def walk(value):
        signature=partial_signature(value)
        if signature is not None:
            head,orders,points=signature
            return partial_derivative(head,orders,tuple(walk(p) for p in points))
        if not value.args: return value
        args=tuple(walk(arg) for arg in value.args)
        rebuilt=value.func(*args,evaluate=False) if isinstance(value,(s.Derivative,Relational)) else value.func(*args)
        if isinstance(rebuilt,s.Subs):
            evaluated=rebuilt.doit(deep=False)
            return walk(evaluated) if evaluated!=rebuilt else rebuilt
        if isinstance(rebuilt,s.Integral) and rebuilt.function==0: return s.S.Zero
        return rebuilt
    resolved=walk(residual)
    if resolved in (s.true,s.false) or isinstance(resolved,Relational) or resolved.func in (s.And,s.Or,s.Not,s.Xor,s.Contains):
        return resolved
    return exact_difference(resolved,s.S.Zero)

def label(n):
    if n.head in ('Text','Symbol','Number'): return n.value
    if call(n,'Missing'):
        # Missing's reason code is metadata, its second argument is the label.
        if len(n.args)>2: return label(n.args[2])
        return 'Missing['+label(n.args[1])+']'
    if n.head=='Sequence': return [label(q) for q in n.args]
    raise Unsupported('not a native label')

def dim_value(n,b,engine,single=False):
    basis=tuple(s.Symbol(b.name('"'+name+'"','wl')) for name in ('length','time','count','momentum'))
    def unit(q):
        if engine=='py': return math_node(q,b,engine)
        if q.head!='Sequence' or len(q.args)!=4: raise Unsupported('dimension vector length')
        return s.prod(v**math_node(e,b,engine) for v,e in zip(basis,q.args))
    if single: return unit(n)
    if n.head!='Sequence': raise Unsupported('dimension set layout')
    return s.FiniteSet(*(unit(q) for q in n.args))

def emit(out,obj):
    out.write(json.dumps(obj,sort_keys=True,default=str)+'\n'); out.flush()


def compare_value(a,b,row,names,conditions,options,output=None):
    if options.get('disable_comparison'): return {'not_compared':'comparison disabled'}
    if row.kind in ('label','basis'):
        aa,bb=label(a.node),label(b.node)
        if row.kind=='basis':
            # Literal basis names translate by the explicit basis table.
            bb=[names.name('"'+x+'"','wl') for x in bb]
        return {'py_value':aa,'wl_value':bb,'text_difference':{'py_only':aa,'wl_only':bb} if aa!=bb else []}
    if row.kind in ('value','condition') and (semantic_children(a) is not None or semantic_children(b) is not None):
        raise Unsupported('structural mismatch: scalar binding has a container without a source-cited component correspondence')
    if row.kind in ('dimension','dimension_field'):
        aa,bb=[dim_value(ref.node,names,e,row.kind=='dimension_field') for ref,e in ((a,'py'),(b,'wl'))]
        if row.kind=='dimension':
            return {'py_value':str(aa),'wl_value':str(bb),'dimension_difference':str(s.SymmetricDifference(aa,bb)),
                    'dimension_evaluation':{'basis_values':{'unitLength':2,'unitTime':3,'unitCount':5,'unitMomentum':7},
                    'py':sorted(str(v.subs(dict(zip(s.symbols('unitLength unitTime unitCount unitMomentum'),(2,3,5,7))))) for v in aa),
                    'wl':sorted(str(v.subs(dict(zip(s.symbols('unitLength unitTime unitCount unitMomentum'),(2,3,5,7))))) for v in bb)}}
    else:
        aa=math_node(a.node,names,'py'); bb=math_node(b.node,names,'wl')
    original=(aa,bb)
    aa=names.resolve(aa,'py'); bb=names.resolve(bb,'wl')
    predicates=lambda z:z in (s.true,s.false) or isinstance(z,Relational) or z.func in (s.And,s.Or,s.Not,s.Xor,s.Contains)
    if predicates(aa)!=predicates(bb): raise Unsupported('structural mismatch: predicate against a non-predicate')
    residual=exact_difference(aa,bb)
    result={'py_value':str(original[0]),'wl_value':str(original[1]),'py_alias_resolved':str(aa),'wl_alias_resolved':str(bb),
            'exact_residual':str(residual),'exact_method':'subtraction, expanded; predicates use symmetric difference',
            'resolved_exact_residual':str(resolved_exact_difference(residual)),
            'resolution_method':'evaluate substitution wrappers preserving partial derivative slots, orders and points; expand and group integrals; no supplied equation or domain assumption',
            'domain_context':[str(v) for v in conditions],
            'condition_residual':str(s.Xor(*conditions,evaluate=False)),
            'one_sided_domains':{'py_only':str(s.And(conditions[0],s.Not(conditions[1]))),'wl_only':str(s.And(conditions[1],s.Not(conditions[0])))},
            'residual_on_supplied_domain':{'value':str(residual),'restricted_to':str(s.And(*conditions)),
                'meaning':'symbolic restriction; evaluated values are printed per sample without replacing the as-computed residual'}}
    if output is not None:
        emit(output,{'kind':'exact_comparison','binding':row.key,**result})
    if row.kind=='condition':
        result['exclusive_domains']={'py_only':str(s.And(aa,s.Not(bb))), 'wl_only':str(s.And(bb,s.Not(aa)))}
        probes=[]
        for angle in (s.Rational(-1,3),s.Rational(1,3),s.Rational(7,2)):
            env={s.Symbol('theta'):angle,s.Symbol('phi'):s.Rational(1,4)}
            values=[evaluate_cached(v.subs(env),FORMS[0],'generic',0) for v in (aa,bb)]
            probes.append({'theta':str(angle),'phi':'1/4','point':0,'conditions':list(map(str,values)),'one_sided':values[0]!=values[1]})
        result['domain_boundary_probes']=probes
    samples=[]
    forms=options.get('forms',FORMS); strata=options.get('strata',STRATA); points=options.get('points',range(len(POINTS)))
    opts=tuple(sorted((k,v) for k,v in options.items() if k in ('duplicate_functions','symmetric_arguments','uncoupled_arguments')))
    for form in forms:
        for stratum in strata:
            for index in points:
                rec={'form':form,'stratum':stratum,'point':index}
                normalized=row.normalized or options.get('normalize_all',False)
                try:
                    domains=[evaluate_cached(v,form,stratum,index,opts) for v in conditions]
                    rec['conditions']=list(map(str,domains))
                    rec['condition_residual']=str(s.Xor(*domains,evaluate=False))
                    if stratum=='count_zero' and normalized:
                        rec.update(coverage='outside normalized domain',count='0',one_sided=domains[0]!=domains[1]); samples.append(rec); continue
                    if any(v==s.false for v in domains):
                        rec.update(coverage='outside intersection',one_sided=domains[0]!=domains[1]); samples.append(rec); continue
                    values=[evaluate_cached(v,form,stratum,index,opts) for v in (aa,bb)]
                    rec.update(py=str(values[0]),wl=str(values[1]),residual=str(exact_difference(*values)),coverage='evaluated')
                except Unsupported as exc: rec.update(coverage='unresolved',reason=str(exc))
                samples.append(rec)
    result['samples']=samples
    return result

def unpaired(ref,used,reason):
    all_paths=leaf_paths(ref.node,ref.path)
    if all_paths<=used: return []
    if not all_paths&used:
        return [{'path':ref.path,'parse':data(ref.node),'leaves':len(all_paths),
                 'reason':reason}]
    return [v for i in range(len(ref.node.args)) for v in unpaired(ref.child(i),used,reason)]

def unpaired_reason(tag,engine):
    if engine=='py' and tag=='momentumSurfaceMotion':
        return 'Py 215-216 includes bulk reference-motion entries; WL 130 prints only the four surface-motion entries. Those extra values stay unpaired. Dictionary keys and constructor leaves are not separately compared values.'
    if engine=='wl' and tag.endswith('.derivation'):
        return 'WL 86,92,101 prints the solved-law lhs and window-weighted windowDerivative/afterParts intermediates; Py 170 prints only the solved rhs and Py 171-173 prints unwindowed intermediates. These are different components. Dictionary keys are not compared values.'
    if tag=='dimensions.fields':
        return 'Py 354-357 prints alias dimensions with the field table; WL prints its table at 217 before adding aliases at 234, so the alias entries have no printed counterpart. Field-name keys establish identities through the name table and are not compared values.'
    if tag.endswith('.comparison') or tag=='coverage.fluxDecomposition':
        return 'Missing constructor and reason-code metadata (WL 15,246) has no separate Py text component (Py 23-26,437). The message text itself is compared; metadata remains unpaired.'
    return 'No source-cited value selection covers this component. Constructor and key leaves are not separately compared values; any unselected content remains printed. The binding table supplies the selected component identities and native layouts.'


def unjoined_reason(tag,engine):
    if '.dimensionSum.' in tag:
        return 'Engine-local dimension-sum audit: sums and reference addends are enumerated from different native expression trees (py 409-434; wl 238-243); no term/reference correspondence by tree ordinal.'
    if tag.startswith('calculus.'):
        return 'SymPy-only Leibniz and generic pullback diagnostic (py 60-63,122-128); Wolfram window differentiation is a different intermediate (wl 74-95).'
    if tag in ('geometry.chart','geometry.centerVelocity','geometry.radiusVelocity','geometry.centerNormalSpeed','geometry.radiusNormalSpeed','surfaces.velocity','surfaces.normalSpeed') or any(tag.startswith(k+'.') for k in ('geometry.chart','geometry.centerVelocity','geometry.radiusVelocity','geometry.centerNormalSpeed','geometry.radiusNormalSpeed','surfaces.velocity','surfaces.normalSpeed')):
        return 'Different printed geometry objects: free-rho chart or split centre/radius velocities/speeds (py 51-54,69-73,113,118-121), versus surface charts and total velocities/speeds (wl 56-58). Surface charts themselves are bound; totals cannot bind individual parts.'
    if tag=='moment.indices': return 'Wolfram-only explicit moment index labels (wl 119); SymPy encodes indices in tags (py 198-200).'
    if any(k in tag for k in ('.weightedLocalRate','.afterParts.bulk','.afterParts.current')):
        return 'SymPy unwindowed IBP intermediates (py 171-173); Wolfram prints window-weighted intermediates (wl 87-95,101), different objects.'
    if tag.endswith('.comparison'):
        return 'Engine-local note on an unjoined construction object (py 23-26; wl 15); no same object/counterpart binding.'
    if tag.endswith('.dimensions'):
        return 'Dimension of an engine-only or differently packaged construction object; native dimension walkers py 359-407,425-434 and wl 230-243. No unused source-cited counterpart.'
    return 'No declared source-cited counterpart; retained for accounting, not a completeness claim.'


def run(py_path,wl_path,out,*,names=NAME_TABLE,bindings=BINDING_TABLE,options=None,progress=True,metadata=True):
    started=time.monotonic(); options=options or {}
    validate_tables(names,bindings)
    streams={e:Stream(p,e) for e,p in (('py',py_path),('wl',wl_path))}
    b=Names(names); b.load(streams)
    emit(out,{'kind':'configuration','sources':'S16 engines accepted at 1d990ee6; decision list revision 4 at 0464b80a',
              'input_sha256':{e:st.sha256 for e,st in streams.items()},'binding_rules':RULES,'surface_table':SURFACE_TABLE,
              'binding_table':[asdict(r) for r in bindings],'name_table':[asdict(n) for n in names],
              'table_sizes':{'rules':len(RULES),'bindings':len(bindings),'names':len(names)},
              'root_shapes':[{'engine':e,'tag':tag,'shape':shape} for (e,tag),shape in ROOT_SHAPES.items()],
              'forms':FORMS,'strata':STRATA,'points':[[str(v) for v in p] for p in POINTS],
              'count_zero_construction':'Sigma is its generic live polynomial minus its exact spatial mean on the shell at the probe time; N vanishes at that time. The subtracted mean is constant in x and t, preserving all density derivatives. No supplied balance is imposed.',
              'residual_target':'No residual target was given.'})
    if metadata:
        for form in FORMS:
            for stratum in STRATA:
                for index in range(len(POINTS)):
                    emit(out,{'kind':'test_functions','form':form,'stratum':stratum,'point':index,
                              'formal_arguments':list(map(str,FORMAL)),
                              'functions':{r.py:str(profile(r.py,r.arity,form,stratum,index)) for r in names if r.arity>0 and not r.alias},
                              'q_forms':{arity:str(profile('q',arity,form,stratum,index)) for arity in (3,4)},
                              'geometry_at_point':{r:str(evaluate_cached(s.Function(r)(s.Symbol('t')),form,stratum,index)) for r in ('r_in','R','N')},
                              'free_chart_point':{'rho':'3/2','theta':str(s.Rational(2+index,5)),'phi':str(s.Rational(index+1,3))}})
    for e,st in streams.items():
        for entry in st.progress: emit(out,{'kind':'progress_marker','engine':e,'line':entry.line,'tag':entry.tag,'raw':entry.raw})
    domains={}
    for mode in ('normalized','unnormalized'):
        vals=[]
        for e,st in streams.items():
            entries=st.index.get('domain.'+mode,[])
            vals.append(math_node(entries[0].node,b,e) if len(entries)==1 else s.true)
        domains[mode]=tuple(vals)
    used={e:defaultdict(set) for e in streams}; bound={e:defaultdict(list) for e in streams}; total=Counter(); below=[]
    unresolved=Counter(); unjoined=Counter(); notformed=[]
    for row in bindings:
        if progress: print('WORKING_ON = '+row.key+' '+row.py.tag+' '+str(row.py.path),file=sys.stderr,flush=True)
        es={e:streams[e].index.get(getattr(row,e).tag,[]) for e in streams}
        emit(out,{'kind':'operands','binding':row.key,**{e:[{'tag':v.tag,'line':v.line,'raw':v.raw} for v in vs] for e,vs in es.items()}})
        refs={}; computed={}; parsed={e:0 for e in streams}; actual={e:0 for e in streams}
        try:
            if any(len(v)!=1 for v in es.values()):
                raise Unsupported('binding end missing or repeated: '+str({e:len(v) for e,v in es.items()}))
            for e,vs in es.items():
                end=getattr(row,e); bound[e][end.tag].append(row.key)
                shape=ROOT_SHAPES.get((e,end.tag))
                if shape is not None and not shaped(Ref(vs[0].node),shape):
                    raise Unsupported('structural mismatch with source-defined full layout: '+e+' '+end.tag)
                if e=='wl' and row.kind=='dimension':
                    dimension_shape=ROOT_SHAPES.get((e,end.tag.removesuffix('.dimensions')),'v')
                    if not shaped_dimension(Ref(vs[0].node),dimension_shape):
                        raise Unsupported('structural mismatch with native dimension-container layout: '+end.tag)
                path=end.path
                if row.kind=='dimension' and e=='wl': path=tuple({'lhs':0,'rhs':1}.get(k,k) for k in path)
                refs[e]=get_ref(Ref(vs[0].node),path)
                parsed[e]=census(refs[e].node)
            computed=compare_value(refs['py'],refs['wl'],row,b,
                                   (s.true,s.true) if row.kind=='condition' else domains['normalized' if row.normalized else 'unnormalized'],options,out)
            if 'not_compared' not in computed:
                for e,ref in refs.items():
                    leaves=leaf_paths(ref.node,ref.path)
                    if row.kind=='label' and call(ref.node,'Missing') and len(ref.node.args)>2:
                        leaves=leaf_paths(ref.node.args[2],ref.path+(2,))
                    if not options.get('omit_leaf_credit'): used[e][getattr(row,e).tag].update(leaves); actual[e]=len(leaves)
            for rec in computed.get('samples',[]):
                if rec['coverage']=='unresolved': unresolved[rec['reason']]+=1
            total['value_comparisons']+=1
        except Unsupported as exc:
            computed={'not_formed':str(exc)}; notformed.append({'binding':row.key,'reason':str(exc)})
            for e,vs in es.items():
                if e not in refs: parsed[e]=sum(census(v.node) for v in vs)
        emit(out,{'kind':'comparison','binding':row.key,**computed})
        emit(out,{'kind':'binding_accounting','binding':row.key,'parsed_leaves':parsed,'actually_compared_leaves':actual})
        total['bindings']+=1
        for e in streams:
            total['bound_parsed_leaves_'+e]+=parsed[e]; total['compared_leaves_'+e]+=actual[e]
        if any(actual[e]<parsed[e] for e in streams): below.append(row.key)
    for e,st in streams.items():
        total['lines_'+e]=st.line_count
        for tag,entries in st.index.items():
            for occurrence,v in enumerate(entries):
                keys=bound[e][tag] if len(entries)==1 else []
                count=census(v.node); actual=len(used[e][tag]) if len(entries)==1 else 0
                rec={'kind':'tag_accounting','engine':e,'tag':tag,'line':v.line,'occurrence':occurrence,
                     'parsed_leaves':count,'actually_compared_leaves':actual,'accounting':'bound' if keys else 'unjoined'}
                total['tags_'+e]+=1; total['parsed_leaves_'+e]+=count; total['tag_compared_leaves_'+e]+=actual
                if keys:
                    rec.update(bindings=keys,unpaired_components=unpaired(Ref(v.node),used[e][tag],unpaired_reason(tag,e) if actual else 'No leaves credited by the active comparisons; see their not-formed/not-compared reasons and source citations.'))
                else:
                    reason='Repeated tag; no comparison formed.' if len(entries)>1 else unjoined_reason(tag,e)
                    rec.update(reason=reason,raw=v.raw,parse=data(v.node)); unjoined[e+' '+reason]+=1; total['unjoined_tags_'+e]+=1
                emit(out,rec)
    emit(out,{'kind':'resources','runtime_seconds':time.monotonic()-started,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    summary={'kind':'summary','counts':dict(total),'bindings_below_parsed_count':below,'not_formed':notformed,
             'unjoined_reasons':dict(unjoined),'unresolved_coverage':dict(unresolved),
             'numerical_integrals':list(NUMERICAL_INTEGRALS.values()),'integral_method':'exact polynomial and trigonometric moments, exact CAS integration, then 40/60-digit quadrature with a difference error estimate; inability remains unresolved',
             'exit_semantics':'zero for findings; nonzero only for operational failure'}
    emit(out,summary); return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--py',required=True,type=Path); parser.add_argument('--wl',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    try:
        with args.output.open('x') as out: run(args.py,args.wl,out)
    except InputError as exc:
        print('Operational input error: '+str(exc),file=sys.stderr); return 2
    return 0

if __name__=='__main__': sys.exit(main())
