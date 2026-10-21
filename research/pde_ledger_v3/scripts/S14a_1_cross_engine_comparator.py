#!/usr/bin/env python3
"""S14a part 1 native-stream comparator. No residual target was given.

The S12/S9b lossless constructor and InputForm readers and S12 partial-slot
translation are adapted here. Neither accepted engine is imported or run.
Bindings name objects, never outcomes. JSONL contains the original operands,
exact residuals before domain restrictions, rational probes and parse census.
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


def name_contract():
    rows = [NameRow(p, w, obj, pl, wl) for p, w, obj, pl, wl in (
        ('x','x1','first Cartesian coordinate','45','34'),
        ('y','x2','second Cartesian coordinate','45','34'),
        ('z','x3','third Cartesian coordinate','45','34'),
        ('w','w','bulk coordinate','45','29,34'),
        ('t','t','time','45','29,34'),
        ('radius','r','spherical radius','46','29'),
        ('inner_radius','ri','inner shell radius','46,263','196,222'),
        ('outer_radius','ro','outer shell radius','46,263','196,222'),
        ('theta','theta','polar angle','47,226','55'),
        ('phi','phi','azimuth','47,226','55'),
        ('None','dimensionOfZero','additive zero has no assigned dimension','282-283,308','229'),
    )]
    for p, spherical, ncomp, arity in (('General',False,4,5),('SphericalSpecialization',True,2,3)):
        stem = 's' if spherical else ''
        wl_lines = '29-33,39-43' if spherical else '34-43'
        for face, wface in (('lower','wLower'),('upper','wUpper')):
            rows.append(NameRow(p+'_face_'+face, stem+wface, p+' '+face+' live face',
                                '130; 222,237',wl_lines,arity-1))
        for k, region in enumerate(('below','inside','above'),1):
            fields = [('n',stem+'n'+str(k),'material number density','138')]
            fields += [('u'+str(i),stem+'u'+str(k)+'c'+str(i+1),
                        ('radial' if spherical and i==0 else 'bulk' if i==ncomp-1 else 'Cartesian '+str(i))+' velocity','139')
                       for i in range(ncomp)]
            fields += [('Jchi'+str(i),stem+'j'+str(k)+'c'+str(i+1),
                        ('radial' if spherical and i==0 else 'bulk' if i==ncomp-1 else 'Cartesian '+str(i))+' extra order flux','140')
                       for i in range(ncomp)]
            fields += [('GammaB',stem+'gamma'+str(k),'net conversion rate','142'),
                       ('GammaReturn',stem+'return'+str(k),'return rate','143'),
                       ('GammaDrain',stem+'drain'+str(k),'drain rate','144')]
            if k==2:
                fields += [('chi',stem+'chi','interior order fraction','141')]
            for field, wl, obj, line in fields:
                rows.append(NameRow(p+'_'+region+'_'+field,wl,p+' '+region+' '+obj,line,wl_lines,arity))
    return tuple(rows)


NAME_TABLE = name_contract()
REGIONS = ('below','inside','above')
V = 'value'
DIM = (V,V,V)  # {material-number, length, time}, source py39-43 / wl17-20


@dataclass(frozen=True)
class End:
    tag: str
    path: tuple = ()


@dataclass(frozen=True)
class Component:
    object: str
    py_path: tuple = ()
    wl_path: tuple = ()
    kind: str = 'value'


@dataclass(frozen=True)
class Binding:
    key: str
    py: End
    wl: End
    object: str
    py_line: str
    wl_line: str
    py_shape: object = V
    wl_shape: object = V
    components: tuple = (Component('value'),)
    chart: str = 'General'
    # External validity contexts constrain selected value components. None
    # records that this engine printed an unconditional standalone value.
    # Contexts are domain diagnostics, not additional bound value leaves.
    condition_context: tuple = ()


def vector_components(labels, kind='value'):
    # Labels enumerate the axes explicitly specified by the sources, not
    # components discovered in measured payloads. Each row is independently tested.
    return tuple(Component(label,(i,),(i,),kind) for i,label in enumerate(labels))


def binding_contract():
    rows=[]
    def add(p,w,obj,pl,wl,*,ps=V,ws=V,parts=None,chart='General',context=()):
        pe=p if isinstance(p,End) else End(p)
        we=w if isinstance(w,End) else End(w)
        rows.append(Binding('b'+str(len(rows)+1),pe,we,obj,str(pl),str(wl),ps,ws,
                            tuple(parts) if parts is not None else (Component(obj),),chart,context))
    for p,w,ncomp in (('General','general projection',4),('SphericalSpecialization','spherical specialization',2)):
        chart=p
        axes=('x','y','z','w') if ncomp==4 else ('radius','w')
        planar=axes[:-1]
        add(p+'/faces/domain',w+'/face ordering','supplied U1 face ordering','149','116',chart=chart)
        add(p+'/U1/masks',w+'/region masks','supplied U1 regional masks','132-134,150','112,119',
            ps=(V,)*3,ws=(V,)*3,parts=vector_components(REGIONS),chart=chart)
        for k,region in enumerate(REGIONS):
            for kind in ('material','order'):
                # Split SymPy fields against the explicitly constructed spacetime
                # current components in WL129,135, with time component first.
                pc=p+'/supplied/'+region+'/'+kind
                wc=w+'/'+kind+'/regional currents'
                add(pc+'_density',End(wc,(k,0)),f'supplied {kind} density, {region}',
                    '109-118,138-148','44-45,129,135',chart=chart)
                for i,axis in enumerate(axes):
                    add(End(pc+'_flux',(i,)),End(wc,(k,i+1)),f'supplied {kind} {axis} current, {region}',
                        '112-118,139-148','45,129,135',chart=chart)
                for side,role in enumerate(('differential operand','right hand side')):
                    add(End(pc+'_equation',(side,)),End(w+'/'+kind+'/regional balance operands',(side,k)),
                        f'supplied {kind} balance {role}, {region}',
                        '120-121,145-148','46,130,136',chart=chart)
            add(p+'/supplied/'+region+'/gamma_equation',End(w+'/supplied gamma identity',(k,)),
                'supplied Gamma_B = Gamma_return - Gamma_drain, '+region,'122,142-148','47,117',
                ps=(V,V),ws='equation',parts=(Component('net rate',(0,),('lhs',)),
                                                           Component('return minus drain',(1,),('rhs',))),chart=chart)
        for kind in ('material','order'):
            pc=p+'/'+kind
            wc=w+'/'+kind
            add(pc+'/differentiated distribution',wc+'/distribution divergence',
                'item 1: differentiated distribution of '+kind,'157-166','131-132,138',chart=chart)
            add(pc+'/concentrated face parts',wc+'/face concentrated coefficients',
                'item 1: '+kind+' concentrated coefficients on their faces','160,168','133,139',
                ps=(V,V),ws=(V,V),parts=vector_components(('lower face trace','upper face trace')),chart=chart)
            add(pc+'/closed integral of differential operand',wc+'/closed integral operand',
                'item 1: closed integral of '+kind+' differential operand','177-178,189','160,169',chart=chart)
            add(pc+'/projected density',wc+'/projected density','item 1: projected '+kind+' density',
                '179,190','148,166',chart=chart)
            add(pc+'/projected in-plane flux',wc+'/projected in-plane flux','item 1: projected '+kind+' in-plane flux',
                '180,191','149,167',ps=(V,)*len(planar),ws=(V,)*len(planar),
                parts=vector_components(planar),chart=chart)
            for i,face in enumerate(('lower','upper')):
                add(pc+f'/face {i}/interior leaving term',End(wc+'/relative interior transport at each face',(i,)),
                    'item 1: '+kind+' interior transport at '+face+' face','182-186,194','152,155',chart=chart)
                add(pc+f'/face {i}/closed-projection leaving term',End(wc+'/each closed-projection face term',(i,)),
                    'item 1: '+kind+' closed-projection transport at '+face+' face','187,195','157,168',chart=chart)
                if kind=='material':
                    add(End(pc+f'/face {i}/supplied material face equation/operands',(0,)),
                        End(wc+'/face balance operands',(i,)),
                        'supplied material face equation differential coefficient, '+face,'197','172',chart=chart)
            add(pc+'/closed projected balance/operands',wc+'/projected balance',
                'item 1: projected '+kind+' balance','199-200; 33-35','159,163-179',
                ps=(V,V),ws='equation',parts=(Component('left including faces',(0,),('lhs',)),
                                              Component('source',(1,),('rhs',))),chart=chart)
            # Both construction checks compare the same two operands; their
            # subtraction orders differ, so their residual tags are NOT bound.
            add(End(pc+'/projection integration/operands',(0,)),wc+'/closed integration/operand_b',
                'integration construction check: projected derivative plus faces','199,201','159,162',chart=chart)
            add(End(pc+'/projection integration/operands',(1,)),wc+'/closed integration/operand_a',
                'integration construction check: integrated differential operand','177-178,201','161-162',chart=chart)
            if kind=='order':
                add(pc+'/projected source',wc+'/projected order source','item 1: projected order source with face parts',
                    '188,198','163,165,177',chart=chart)
            # Dimensional exponent values only. Operand expression trees and
            # dimension-ledger term lists do not acquire bindings by identity.
            pd=p+'/dimensions/'+kind
            wd=w+'/dimensions/'+kind
            add(pd+'/density/dimension',wd+'/density','item 4: projected '+kind+' density dimensions',
                '305-308,321','249-250',ps=DIM,ws=DIM,parts=vector_components(('number exponent','length exponent','time exponent')),chart=chart)
            for i,axis in enumerate(planar):
                add(pd+f'/flux {i}/dimension',End(wd+'/plane_flux',(i,)),
                    'item 4: projected '+kind+' '+axis+' flux dimensions','322-323','249-250',
                    ps=DIM,ws=DIM,parts=vector_components(('number exponent','length exponent','time exponent')),chart=chart)
            for i,face in enumerate(('lower','upper')):
                add(pd+f'/face {i}/dimension',End(wd+'/faces',(i,)),
                    'item 4: '+kind+' '+face+' transport dimensions','324-325','249-250',
                    ps=DIM,ws=DIM,parts=vector_components(('number exponent','length exponent','time exponent')),chart=chart)
            add(pd+'/balance left/dimension',wd+'/left','item 4: '+kind+' balance left dimensions',
                '326-327','249-250',ps=DIM,ws=DIM,parts=vector_components(('number exponent','length exponent','time exponent')),chart=chart)
            add(pd+'/balance right/dimension',wd+'/supplied_source','item 4: '+kind+' balance source dimensions',
                '326-327','249-250',ps=DIM if kind=='order' else V,ws=DIM if kind=='order' else V,
                parts=vector_components(('number exponent','length exponent','time exponent')) if kind=='order' else None,chart=chart)
    shell='spherical specialization / shell/'
    for name,other,pl,wl in (
        ('shell domain','domain','263','222'),('face channel','face channel','251-252','201'),
        ('core channel','core channel','253','204-206'),('J','J','254','208'),
        ('order lost','order lost','255','209'),('outward outer flux','outer material flux','256','207'),
        ('storage','storage','257','211-213'),('outer flux radius derivative','outer flux radius derivative','261','217')):
        add('SphericalSpecialization/'+name,shell+other,'items 2-3: '+name,pl,wl+'; 222',chart='SphericalSpecialization')
    add('SphericalSpecialization/shell balance/operands',shell+'balance','item 2: shell balance','258-269','214,222',
        ps=(V,V),ws='equation',parts=(Component('shell left',(0,),('lhs',)),Component('shell source',(1,),('rhs',))),chart='SphericalSpecialization')
    add('SphericalSpecialization/shell balance/residual',shell+'balance residual','item 2: shell balance residual',
        '258-269;35','214,222',chart='SphericalSpecialization')
    add('SphericalSpecialization/face channel minus order lost/residual',shell+'face channel minus order lost',
        'item 2: face channel minus order lost','260,270;35','215,222',chart='SphericalSpecialization')
    for name,other,pl,wl in (('outer projected velocity','projected velocity','275-276','218-219'),
                             ('outer volume flux','projected velocity volume flux','273-277','218-220')):
        add('SphericalSpecialization/'+name,shell+other,'item 3: '+name,pl,wl+';222',
            ps='conditioned',ws='conditioned',parts=(Component('nonzero projected density',('condition',),('condition',),'condition'),
                                                     Component('conditional value',('value',),('value',),'conditional_value')),chart='SphericalSpecialization')
    add('Spherical specialization/sphere area','spherical specialization / sphere measure',
        'items 2-3: spherical area measure','226-235','55-62',chart='SphericalSpecialization')
    for name,other,pl,wl in (('face channel','face channel','329-335','201'),('core channel','core channel','329-335','206'),
        ('J','J','329-335','208'),('order lost','order lost','329-335','209'),
        ('outward outer flux','outer material flux','329-335','207'),('storage','storage','329-335','213'),
        ('shell balance','balance residual','329-335','214'),('face channel minus order lost','face channel minus order lost','329-335','215'),
        ('outer flux radius derivative','outer flux radius derivative','329-335','217')):
        add('SphericalSpecialization/dimensions/'+name+'/dimension','spherical specialization / dimensions/'+other,
            'item 4: '+name+' dimensions',pl,wl+';255-256',ps=DIM,ws=DIM,
            parts=vector_components(('number exponent','length exponent','time exponent')),chart='SphericalSpecialization')
    for name,other in (('projected velocity','projected velocity'),('volume flux','projected velocity volume flux')):
        add('SphericalSpecialization/dimensions/'+name+'/dimension',
            End('spherical specialization / dimensions/'+other,('value',)),
            'item 4: '+name+' dimension value; py308 emits it unconditionally, wl238 attaches a condition','336-337;308','238,255-256',
            ps=DIM,ws=DIM,parts=vector_components(('number exponent','length exponent','time exponent')),chart='SphericalSpecialization',
            context=(None,End('spherical specialization / dimensions/'+other,('condition',))))
    add(End('SphericalSpecialization/dimensions/volume flux/dimension operand',('condition',)),
        End('spherical specialization / dimensions/projected velocity volume flux',('condition',)),
        'item 4: volume-flux dimension domain, nonzero outer projected density','274,305-306,337','218,238,256',
        parts=(Component('printed nonzero-density condition',kind='condition'),),chart='SphericalSpecialization')
    # No-engine-own-check labels are compared only for whole-tag object bindings.
    # Split/grouped label occurrences have no injective counterpart, and remain
    # accounted as unjoined. No labels or values are matched by string identity.
    for row in tuple(rows):
        if not row.py.path and not row.wl.path and '/operands' not in row.py.tag and not row.py.tag.endswith(('/residual','/dimension')):
            add(row.py.tag+'/comparison',row.wl.tag+'/comparison','comparison note about '+row.object,
                '27-30; '+row.py_line,'14; '+row.wl_line,parts=(Component('comparison note',kind='text'),),chart=row.chart)
    return tuple(rows)


BINDING_TABLE = binding_contract()


def validate_tables(names, rows):
    for field in ('py','wl'):
        values=[getattr(r,field) for r in names]
        if len(values)!=len(set(values)):
            raise InputError('duplicate name row: '+field)
    if len({r.key for r in rows})!=len(rows):
        raise InputError('duplicate binding identifier')
    for field in ('py','wl'):
        ends=[]
        for row in rows:
            end=getattr(row,field)
            for old in ends:
                if old.tag==end.tag and (old.path==end.path[:len(old.path)] or end.path==old.path[:len(end.path)]):
                    raise InputError('duplicate or overlapping binding end: '+field+' '+str(end))
            ends.append(end)
        for row in rows:
            paths=[getattr(c,field+'_path') for c in row.components]
            for i,path in enumerate(paths):
                if any(path==other[:len(path)] or other==path[:len(other)] for other in paths[:i]):
                    raise InputError('duplicate component correspondence: '+row.key+' '+field)


def census(n):
    """Independent count from the parser, before extraction or joining."""
    return sum(census(x) for x in n.args) if n.args else 1


def leaf_paths(n,prefix=()):
    if not n.args:
        return {prefix}
    return set().union(*(leaf_paths(x,prefix+(i,)) for i,x in enumerate(n.args)))


def data(n):
    return [n.head,n.value,[data(x) for x in n.args]]


def call(n,name):
    return n.head=='Apply' and n.args[0].head=='Symbol' and n.args[0].value==name


def application(n):
    if n.head!='Apply':
        return None,()
    head=n.args[0]
    if head.head=='Symbol':
        return head.value,n.args[1:]
    if call(head,'Inactive') and len(head.args)==2 and head.args[1].head=='Symbol':
        return head.args[1].value,n.args[1:]
    return None,n.args[1:]


def is_comparison_note(tag):
    # The extra spelling is emitted literally at WL71, outside alone().
    return tag.endswith('/comparison') or tag=='transverse integration / comparison'


def is_label(tag):
    return tag in ('contract','engine','WORKING_ON') or is_comparison_note(tag) or tag.endswith('/representation')


def read_label(text,engine,tag):
    if engine=='wl':
        result=wl_parse(text)
        if result.head not in ('Text','Sequence') or (result.head=='Sequence' and any(n.head!='Text' for n in result.args)):
            raise InputError('label must be a string or string list')
        return result
    if tag=='engine':
        if text[:1] not in ('(','[') or text[-1:] not in (')',']'):
            raise InputError('engine label list framing')
        parts=text[1:-1].rsplit(', ',1) # py221: runtime string, SymPy version
        return seq(*(atom('Text',part) for part in parts))
    # sstr omits quotes even for plain labels. Tests and measured reads share it.
    return atom('Text',text)


@dataclass
class Entry:
    tag: str
    raw: str
    line: int
    node: Node


class Stream:
    def __init__(self,path,engine):
        self.path,self.engine=Path(path),engine
        self.index=defaultdict(list)
        self.progress=[]
        try:
            raw=self.path.read_bytes()
            lines=raw.decode('utf-8').splitlines()
        except (OSError,UnicodeError) as exc:
            raise InputError(str(exc)) from exc
        self.sha256=hashlib.sha256(raw).hexdigest()
        framed=[]
        for i,line in enumerate(lines,1):
            m=re.fullmatch(r'([A-Za-z][^:]*): (.*)',line)
            if m:
                framed.append([m[1],m[2],i])
            elif framed and framed[-1][0]!='WORKING_ON' and line.startswith((' ','[',']',')','}')):
                framed[-1][1]+='\n'+line
            else:
                raise InputError(f'{path}:{i}: unassigned line')
        if not framed:
            raise InputError('empty stream '+str(path))
        for tag,text,line in framed:
            try:
                node=read_label(text,engine,tag) if is_label(tag) else (py_parse if engine=='py' else wl_parse)(text)
            except InputError as exc:
                raise InputError(f'{path}:{line}:{tag}: {exc}') from exc
            entry=Entry(tag,text,line,node)
            (self.progress if tag=='WORKING_ON' else self.index[tag]).append(entry)
        self.line_count=len(lines)


@dataclass(frozen=True)
class Ref:
    node: Node
    path: tuple=()

    def child(self,index):
        if index<0 or index>=len(self.node.args):
            raise Unsupported('component missing at '+str(self.path+(index,)))
        return Ref(self.node.args[index],self.path+(index,))


def conditioned_refs(ref):
    n=ref.node
    h,a=application(n)
    if h=='ConditionalExpression' and len(a)==2:
        return {'value':ref.child(1),'condition':ref.child(2)}
    if n.head=='Dictionary':
        pairs={}
        for i,pair in enumerate(n.args):
            key=pair.args[0]
            if key.head not in ('Symbol','Text') or key.value not in ('object','domain') or key.value in pairs:
                return None
            pairs[key.value]=ref.child(i).child(1)
        if set(pairs)=={'object','domain'}:
            return {'value':pairs['object'],'condition':pairs['domain']}
    return None


def get_ref(ref,path):
    for part in path:
        if isinstance(part,int):
            if ref.node.head!='Sequence':
                raise Unsupported('expected sequence at '+str(ref.path))
            ref=ref.child(part)
        elif part in ('condition','value'):
            values=conditioned_refs(ref)
            if values is None:
                raise Unsupported('expected conditioned object at '+str(ref.path))
            ref=values[part]
        elif part in ('lhs','rhs'):
            h,a=application(ref.node)
            if h not in ('Eq','Equal') or len(a)!=2:
                raise Unsupported('expected equation at '+str(ref.path))
            ref=ref.child(1 if part=='lhs' else 2)
        else:
            raise InputError('unknown component selector '+str(part))
    return ref


def shaped(ref,schema):
    n=ref.node
    if schema==V:
        return n.head not in ('Sequence','Dictionary') and conditioned_refs(ref) is None
    if schema=='equation':
        h,a=application(n)
        return h in ('Eq','Equal') and len(a)==2
    if schema=='conditioned':
        return conditioned_refs(ref) is not None
    return n.head=='Sequence' and len(n.args)==len(schema) and all(shaped(ref.child(i),sub) for i,sub in enumerate(schema))



def root_contract():
    """Source-defined full layouts checked before selecting any binding end.

    A place denotes a semantic component only within its complete native
    container. These layouts are not learned from either measured stream.
    """
    schemas={}
    for p,w,ncomp in (('General','general projection',4),('SphericalSpecialization','spherical specialization',2)):
        schemas[('wl',w+'/supplied gamma identity')]=('equation',)*3
        for region in REGIONS:
            for kind in ('material','order'):
                base=p+'/supplied/'+region+'/'+kind
                schemas[('py',base+'_flux')]=(V,)*ncomp
                schemas[('py',base+'_equation')]=(V,V)
            schemas[('py',p+'/supplied/'+region+'/gamma_equation')]=(V,V)
        for kind in ('material','order'):
            wc=w+'/'+kind
            pc=p+'/'+kind
            schemas[('wl',wc+'/regional currents')]=((V,)*(ncomp+1),)*3
            schemas[('wl',wc+'/regional balance operands')]=((V,)*3,(V,)*3)
            for name in ('relative interior transport at each face','each closed-projection face term','face balance operands'):
                schemas[('wl',wc+'/'+name)]=(V,V)
            schemas[('py',pc+'/projection integration/operands')]=(V,V)
            for i in (0,1):
                schemas[('py',pc+f'/face {i}/supplied material face equation/operands')]=(V,V)
            schemas[('wl',w+'/dimensions/'+kind+'/plane_flux')]=(DIM,)*(ncomp-1)
            schemas[('wl',w+'/dimensions/'+kind+'/faces')]=(DIM,DIM)
    for name in ('projected velocity','projected velocity volume flux'):
        schemas[('wl','spherical specialization / dimensions/'+name)]='conditioned'
    return schemas


ROOT_SCHEMAS=root_contract()


class Names:
    def __init__(self,rows):
        self.maps={'py':{r.py:r.py for r in rows},'wl':{r.wl:r.py for r in rows}}
        self.arities={r.py:r.arity for r in rows}

    def name(self,name,engine):
        return self.maps[engine].get(name,engine+'::'+name)


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
    slots=points if direct else s.symbols('_s14partial0:'+str(len(points)),real=True)
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
                local[var]=s.Symbol('_int'+str(depth+i),real=True)
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
    return alpha_integrals(_math_node(n,b,engine,scope))


def _math_node(n,b,engine,scope=None):
    scope={} if scope is None else scope
    def cv(x,local=None):
        return _math_node(x,b,engine,scope if local is None else local)
    h,a=n.head,n.args
    if h=='Number':
        if n.value=='None':
            return s.Symbol(b.name('None',engine),real=True)
        return s.Rational(n.value)
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
        return s.Symbol(key,real=True)
    if h=='Apply':
        head,*args=a
        if head.head=='Apply' and call(head.args[0],'Derivative'):
            orders=head.args[0].args[1:]
            if len(head.args)!=2 or head.args[1].head!='Symbol':
                raise Unsupported('derivative head structure')
            key=b.name(head.args[1].value,engine)
            if b.arities.get(key)!=len(args) or len(orders)!=len(args):
                raise Unsupported('derivative function arity '+key)
            return partial_derivative(s.Function(key,real=True),tuple(int(cv(q)) for q in orders),tuple(cv(q) for q in args))
        name,args=application(n)
        if name is None:
            raise Unsupported('composite application head')
        if name in b.maps[engine] and b.arities.get(b.name(name,engine),0):
            key=b.name(name,engine)
            if len(args)!=b.arities[key]:
                raise Unsupported('live function arity '+key)
            return s.Function(key,real=True)(*(cv(q) for q in args))
        h,a=BUILTINS.get(name,name),args
    if h in ('Add','Mul','Pow'):
        return getattr(s,h)(*(cv(q) for q in a))
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
            var=s.Dummy('integral_bound',real=True)
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
            local[var.value]=s.Symbol('_slot'+str(len(scope)+i),real=True)
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

# Probe formulas are mathematical inputs, not solutions of supplied balances.
FORMS=('affine','quadratic','cubic')
STRATA=('generic','jump_lower','jump_upper','order_zero_lower','order_zero_upper','density_zero')
LOCATIONS=('below','inside','above','lower','upper')
BASE_POINTS=((s.Rational(2,9),s.Rational(-1,5),s.Rational(3,11),s.Rational(1,7)),
             (s.Rational(1,6),s.Rational(3,10),s.Rational(-1,8),s.Rational(-2,7)),
             (s.Rational(-2,11),s.Rational(1,9),s.Rational(2,7),s.Rational(3,8)))
FORMAL=s.symbols('a0:5',real=True)
X,Y,Z,W,T=s.symbols('x y z w t',real=True)
R,A,B=s.symbols('radius inner_radius outer_radius',real=True)
NAME_CODES={row.py:i+1 for i,row in enumerate(NAME_TABLE)}
SCALAR_NAMES={row.py for row in NAME_TABLE if not row.arity}


def seed(name):
    if name in NAME_CODES:
        return s.Integer(NAME_CODES[name])
    code=int.from_bytes(b'\x01'+name.encode(),'big')
    return len(NAME_CODES)+1+s.Rational(code,code+1)


def argument_coupling(args,code):
    """Every first partial depends on every other argument in this form."""
    return (code+1)*sum((i+1)*(j+1)*args[i]*args[j]
                        for i in range(len(args)) for j in range(i+1,len(args)))/211


def basic_form(form,args,code):
    if form=='affine':
        return (code+1)/7+sum((i+1)*(code+1)*q/101 for i,q in enumerate(args))
    if form=='quadratic':
        return (code+1)/7+sum(((i+1)**2*(code+1)*q*q+(i+1)*q)/103 for i,q in enumerate(args))+argument_coupling(args,code)
    if form=='cubic':
        return (code+1)/7+sum(((i+1)**3*(code+1)*q**3+(i+1)**2*q)/107 for i,q in enumerate(args))
    raise InputError('unknown probe form '+form)


def face_arguments(chart):
    return (X,Y,Z,T) if chart=='General' else (R,T)


@lru_cache(None)
def face_forms(chart,form):
    args=face_arguments(chart)
    lower=-1+basic_form(form,args,s.Integer(1))/100
    upper=lower+3+(1+sum((i+1)*q*q for i,q in enumerate(args)))/100
    return lower,upper


def face_objects(chart):
    return tuple(s.Function(chart+'_face_'+side,real=True)(*face_arguments(chart)) for side in ('lower','upper'))


@lru_cache(None)
def profile(name,arity,form,stratum):
    args=FORMAL[:arity]
    chart='SphericalSpecialization' if name.startswith('SphericalSpecialization_') else 'General'
    face_args=args if '_face_' in name else args[:-2]+(args[-1],)
    lo,hi=face_forms(chart,form)
    replace=dict(zip(face_arguments(chart),face_args))
    lo,hi=lo.xreplace(replace),hi.xreplace(replace)
    if name==chart+'_face_lower':
        return lo
    if name==chart+'_face_upper':
        return hi
    base=basic_form(form,args,seed(name))
    if name.endswith('_chi'):
        base=s.Rational(1,3)+base/10000
        if stratum in ('order_zero_lower','order_zero_upper'):
            distance=args[-2]-(lo if stratum.endswith('lower') else hi)
            base=distance*(1+base)/100
    if stratum.startswith('jump_') and any('_'+region+'_' in name for region in REGIONS):
        region=next(region for region in REGIONS if '_'+region+'_' in name)
        field=name.split('_'+region+'_')[1]
        # U1 has no exterior order field. Other one-sided native fields share
        # only their trace at the selected face; their full functions differ.
        if field!='chi':
            common=basic_form(form,args,s.Integer(200+sum((i+1)*ord(q) for i,q in enumerate(field))))
            distance=args[-2]-(lo if stratum.endswith('lower') else hi)
            base=common+seed(name)*distance**2*(1+sum((i+1)*q*q for i,q in enumerate(args)))/109
    if stratum=='density_zero' and name.endswith('_inside_n'):
        # A nonconstant profile whose integral between the moving faces
        # vanishes, with no replacement of a live function by a constant.
        offset=args[-2]-(lo+hi)/2
        planar=args[:-2]+(args[-1],)
        live=basic_form(form,planar,seed(name))
        base=offset*(1+live)+(offset**2-(hi-lo)**2/12)*(seed(name)+1)/113
        if form=='cubic':
            base+=offset**3*(1+live)/127
    return base


@lru_cache(maxsize=4096)
def test_expression(expr,form,stratum):
    def substitute(fun):
        return profile(fun.func.__name__,len(fun.args),form,stratum).xreplace(dict(zip(FORMAL,fun.args)))
    value=expr.replace(lambda q:isinstance(q,AppliedUndef),substitute)
    # Evaluate derivatives after complete function substitution; held integrals
    # remain held until the point's free coordinates have been substituted.
    value=value.replace(lambda q:isinstance(q,s.Derivative),lambda q:q.doit(deep=False))
    value=value.replace(lambda q:isinstance(q,s.Subs),lambda q:q.doit())
    return value


@lru_cache(None)
def point_for(chart,form,index,location):
    x,y,z,t=BASE_POINTS[index]
    env={X:x,Y:y,Z:z,T:t,R:s.Rational(5+index,4),A:s.Rational(1,2),B:s.Rational(3+index,1),
         s.Symbol('theta',real=True):s.Rational(index+1,5),s.Symbol('phi',real=True):s.Rational(index+1,7)}
    # Additional dimension markers and synthetic mapped coordinates get
    # distinct rational probes; radii and native coordinates keep the above.
    for name,code in NAME_CODES.items():
        key=s.Symbol(name,real=True)
        if name in SCALAR_NAMES and key not in env and key!=W:
            env[key]=s.Rational(code+index+1,19)
    lo,hi=(s.cancel(v.xreplace(env)) for v in face_forms(chart,form))
    env[W]={'below':lo-s.Rational(index+1,4),'inside':lo+(hi-lo)*s.Rational(index+1,4),
            'above':hi+s.Rational(index+1,4),'lower':lo,'upper':hi}[location]
    return env


def free_atoms(expr,kind):
    if isinstance(expr,s.Integral):
        return set()
    result={expr} if isinstance(expr,kind) else set()
    for arg in expr.args:
        result.update(free_atoms(arg,kind))
    return result


def delta_integral(integral,*,symbolic=False):
    """Contract's closed measure: every affine endpoint delta has full weight.

    Return the regular integral, all singular contributions, and a ledger.
    Never delegate an endpoint delta to a CAS integration convention.
    """
    if len(integral.limits)!=1:
        raise Unsupported('delta integral with multiple limits')
    var,lo,hi=integral.limits[0]
    body=s.expand(integral.function)
    atoms=free_atoms(body,s.DiracDelta)
    if not atoms:
        return integral,[]
    regular=body.xreplace(dict.fromkeys(atoms,s.S.Zero))
    contributions=[]
    counted=[]
    removed=s.S.Zero
    for delta in sorted(atoms,key=str):
        if len(delta.args)>1 and delta.args[1]!=0:
            raise Unsupported('higher derivative of delta in integral')
        arg=delta.args[0]
        slope=s.diff(arg,var)
        if slope==0 or slope.has(var):
            raise Unsupported('nonaffine delta integral support')
        root=s.cancel(var-arg/slope)
        coefficient=body.coeff(delta)
        if coefficient.has(s.DiracDelta) or body.has(delta**2):
            raise Unsupported('nonlinear delta product in integral')
        at_endpoint=s.expand(root-lo)==0 or s.expand(root-hi)==0
        if symbolic:
            if not at_endpoint:
                raise Unsupported('symbolic delta integral membership not isolated')
            weight=s.S.One
        else:
            if lo.free_symbols or hi.free_symbols or root.free_symbols:
                raise Unsupported('delta integral retains free membership variables')
            weight=s.Integer(1 if bool(s.Le(lo,root)) and bool(s.Le(root,hi)) else 0)
        contribution=coefficient.subs(var,root)/s.Abs(slope)*weight
        contributions.append(contribution)
        removed+=coefficient*delta
        counted.append({'integral':str(integral),'support':str(root),'endpoint':at_endpoint,
                        'closed_weight':str(weight),'contribution':str(contribution)})
    if s.expand(body-regular-removed)!=0:
        raise Unsupported('delta integral has unisolated singular content')
    rest=s.Integral(regular,(var,lo,hi)) if regular!=0 else s.S.Zero
    return rest+s.Add(*contributions),counted


def exact_endpoint_integrals(expr):
    ledger=[]
    def translate(integral):
        try:
            value,rows=delta_integral(integral,symbolic=True)
        except Unsupported:
            return integral
        ledger.extend(rows)
        return value
    return expr.replace(lambda q:isinstance(q,s.Integral) and q.has(s.DiracDelta),translate),ledger


def evaluated_integrals(expr):
    ledger=[]
    def one(integral):
        if integral.function.has(s.DiracDelta):
            value,counted=delta_integral(integral)
            ledger.extend(dict(method='closed delta measure',**r) for r in counted)
            if value==integral:
                raise Unsupported('delta integral did not isolate')
            evaluated,extra=evaluated_integrals(value)
            ledger.extend(extra)
            return evaluated
        result=integral.function
        for var,lo,hi in integral.limits:
            if result.has(s.Integral):
                raise Unsupported('nested integral remains unresolved')
            try:
                poly=s.Poly(result,var)
            except s.PolynomialError:
                exact=s.Integral(result,(var,lo,hi)).doit()
                if exact.has(s.Integral):
                    if exact.free_symbols:
                        raise Unsupported('numerical integral retains free parameters')
                    import mpmath as mp
                    f=s.lambdify(var,result,'mpmath')
                    with mp.workdps(50):
                        low,errlow=mp.quad(f,[mp.mpf(str(s.N(lo,55))),mp.mpf(str(s.N(hi,55)))],error=True)
                    with mp.workdps(80):
                        high,errhigh=mp.quad(f,[mp.mpf(str(s.N(lo,85))),mp.mpf(str(s.N(hi,85)))],error=True)
                        error=max(abs(high-low),errlow,errhigh)
                        ledger.append({'integral':str(integral),'method':'numerical quadrature','working_precision':[50,80],
                                       'error_estimate':mp.nstr(error,15)})
                        exact=s.Float(mp.nstr(high,80),80)
                result=exact
            else:
                primitive=poly.integrate().as_expr()
                result=primitive.subs(var,hi)-primitive.subs(var,lo)
        return result
    value=expr.replace(lambda q:isinstance(q,s.Integral),one)
    return value,ledger


@lru_cache(maxsize=8192)
def evaluate(expr,chart,form,stratum,index,location):
    value=test_expression(expr,form,stratum).xreplace(point_for(chart,form,index,location))
    value,ledger=evaluated_integrals(value)
    value=value.doit()
    if isinstance(value,s.Expr):
        value=s.cancel(value)
    if value.has(AppliedUndef,s.Derivative,s.Subs,s.Integral) or value.free_symbols:
        raise Unsupported('evaluation retains functions, derivatives, integrals or free symbols: '+str(value))
    if value.has(s.zoo,s.nan,s.oo,-s.oo):
        raise Unsupported('nonfinite evaluation')
    return value,ledger


def supplied_domain_value(expr,chart):
    lo,hi=face_objects(chart)
    gap=hi-lo
    replacements={}
    for step in expr.atoms(s.Heaviside):
        arg=s.expand(step.args[0])
        if s.expand(arg-gap)==0:
            replacements[step]=s.S.One
        elif s.expand(arg+gap)==0:
            replacements[step]=s.S.Zero
    for absolute in expr.atoms(s.Abs):
        if s.expand(absolute.args[0]-gap)==0 or s.expand(absolute.args[0]+gap)==0:
            replacements[absolute]=gap
    return s.expand(expr.xreplace(replacements)) if isinstance(expr,s.Expr) else expr.xreplace(replacements)


def exact_difference(a,b):
    if isinstance(a,Boolean) and not isinstance(a,s.Expr):
        if not isinstance(b,Boolean) or isinstance(b,s.Expr):
            raise Unsupported('Boolean and value have different kinds')
        return s.Xor(a,b)
    if isinstance(a,Relational) or isinstance(b,Relational):
        return s.Xor(a,b)
    return s.expand(a-b)


def region_part(expr,region,chart):
    replacements={}
    unknown=[]
    for step in free_atoms(expr,s.Heaviside):
        found=False
        for i,face in enumerate(face_objects(chart)):
            for sign in (1,-1):
                ratio=s.cancel(step.args[0]/(sign*(W-face)))
                if ratio.is_number and ratio.is_positive:
                    replacements[step]=s.Integer((region>i) if sign==1 else (region<=i))
                    found=True
                    break
            if found:
                break
        if not found and step.args[0].has(W):
            unknown.append(str(step))
    return s.expand(expr.xreplace(replacements)),unknown


def distribution_parts(expr,chart):
    expr=s.expand(expr)
    deltas=free_atoms(expr,s.DiracDelta)
    regular=expr.xreplace(dict.fromkeys(deltas,s.S.Zero))
    parts={}
    coverage=[]
    for i,region in enumerate(REGIONS):
        part,unknown=region_part(regular,i,chart)
        parts['regular/'+region]=part
        if unknown:
            coverage.append({'part':'regular/'+region,'unisolated_steps':unknown})
    removed=s.S.Zero
    for delta in sorted(deltas,key=str):
        if len(delta.args)>1 and delta.args[1]!=0:
            coverage.append({'atom':str(delta),'reason':'higher delta derivative'})
            continue
        slope=s.diff(delta.args[0],W)
        if slope==0 or slope.has(W):
            coverage.append({'atom':str(delta),'reason':'not affine in w'})
            continue
        root=s.cancel(W-delta.args[0]/slope)
        face=next((i for i,f in enumerate(face_objects(chart)) if s.expand(root-f)==0),None)
        coefficient=expr.coeff(delta)
        if face is None or coefficient.has(s.DiracDelta) or expr.has(delta**2):
            coverage.append({'atom':str(delta),'reason':'unknown support or nonlinear singular product'})
            continue
        key='delta/'+('lower','upper')[face]
        trace=s.Subs(coefficient,W,face_objects(chart)[face]).doit()/s.Abs(slope)
        parts[key]=parts.get(key,s.S.Zero)+trace
        removed+=coefficient*delta
    for key in ('delta/lower','delta/upper'):
        parts.setdefault(key,s.S.Zero)
    leftover=s.expand(expr-regular-removed)
    if leftover!=0:
        coverage.append({'unisolated_content':str(leftover)})
    return parts,coverage


def text_difference(a,b):
    return {'equal_text':a==b,'first_different_character':next((i for i,(x,y) in enumerate(zip(a,b)) if x!=y),
                                                              min(len(a),len(b)) if a!=b else None)}


def sample_values(a,b,chart,part,*,conditions=None,forms=FORMS,points=range(3),strata=STRATA,locations=None):
    result=[]
    if locations is None:
        locations=(part.split('/')[-1],) if part.startswith(('regular/','delta/')) else ('inside',)
    for form in forms:
        for stratum in strata:
            for index in points:
                for loc in locations:
                    row={'form':form,'stratum':stratum,'point_index':index,'location':loc,
                         'point':{str(k):str(v) for k,v in point_for(chart,form,index,loc).items() if str(k) in SCALAR_NAMES}}
                    try:
                        env=point_for(chart,form,index,loc)
                        lo,hi=[f.xreplace(env) for f in face_forms(chart,form)]
                        row['supplied_domain']={'face_order':str(s.Lt(lo,hi)),
                            'shell_order':str(s.And(s.Gt(env[A],0),s.Lt(env[A],env[B]))),'positive_radius':str(s.Gt(env[R],0))}
                        if conditions is not None:
                            ac,al=evaluate(conditions[0],chart,form,stratum,index,loc)
                            bc,bl=evaluate(conditions[1],chart,form,stratum,index,loc)
                            row['conditions']={'py':str(ac),'wl':str(bc),'exclusive':str(s.Xor(ac,bc))}
                            if ac!=s.true or bc!=s.true:
                                row['value_coverage']='outside conjunction of printed conditions'
                                result.append(row)
                                continue
                        av,al=evaluate(a,chart,form,stratum,index,loc)
                        bv,bl=evaluate(b,chart,form,stratum,index,loc)
                        row.update(py_value=str(av),wl_value=str(bv),residual=str(exact_difference(av,bv)))
                        if al or bl:
                            row['integral_evaluations']={'py':al,'wl':bl}
                    except Unsupported as exc:
                        row['unresolved_coverage']=str(exc)
                    result.append(row)
    return result


ACTIVE_OUTPUT=None


def emit(out,obj):
    out.write(json.dumps(obj,separators=(',',':'),ensure_ascii=False)+'\n')
    out.flush()


def compare_values(a,b,row,part,options=None,conditions=None):
    result={'part':part,'py_operand':str(a),'wl_operand':str(b)}
    if (isinstance(a,Boolean) and not isinstance(a,s.Expr)) != (isinstance(b,Boolean) and not isinstance(b,s.Expr)):
        return dict(result,not_formed='different kinds of mathematical object')
    if isinstance(a,s.Equality) or isinstance(b,s.Equality):
        if not isinstance(a,s.Equality) or not isinstance(b,s.Equality):
            return dict(result,not_formed='different relation kinds')
        return dict(result,components=[compare_values(a.lhs,b.lhs,row,part+'/lhs',options),
                                       compare_values(a.rhs,b.rhs,row,part+'/rhs',options)])
    a,la=exact_endpoint_integrals(a)
    b,lb=exact_endpoint_integrals(b)
    if la or lb:
        result['closed_endpoint_integrals']={'py':la,'wl':lb}
        result['closed_measure_operands']={'py':str(a),'wl':str(b)}
    try:
        residual=exact_difference(a,b)
    except (Unsupported,TypeError,ValueError,NotImplementedError) as exc:
        result['exact_unresolved']=str(exc)
    else:
        result['as_computed_residual']=str(residual)
        if ACTIVE_OUTPUT is not None:
            emit(ACTIVE_OUTPUT,dict(kind='residual_before_domain',binding=row.key,**result))
        result['on_supplied_domain_residual']=str(supplied_domain_value(residual,row.chart))
    if conditions is not None:
        result['value_restriction']={'py_condition':str(conditions[0]),'wl_condition':str(conditions[1]),
                                     'common_condition':str(s.And(*conditions))}
    sample_options=dict(options or {})
    if 'locations' not in sample_options and not part.startswith(('regular/','delta/')):
        match=re.search(r'/supplied/(below|inside|above)/',row.py.tag)
        face=re.search(r'/face ([01])/',row.py.tag)
        if match:
            sample_options['locations']=(match[1],)
        elif face:
            sample_options['locations']=(('lower','upper')[int(face[1])],)
        elif part.startswith(('lower face','upper face')):
            sample_options['locations']=(part.split()[0],)
    result['samples']=sample_values(a,b,row.chart,part,conditions=conditions,**sample_options)
    return result


def compare_scalar(a,b,row,part,options=None,conditions=None):
    deltas=free_atoms(a,s.DiracDelta)|free_atoms(b,s.DiracDelta)
    steps=free_atoms(a,s.Heaviside)|free_atoms(b,s.Heaviside)
    if deltas or any(q.args[0].has(W) for q in steps):
        pa,ca=distribution_parts(a,row.chart)
        pb,cb=distribution_parts(b,row.chart)
        return {'distribution_parts':[compare_values(pa[key],pb[key],row,key,options,conditions) for key in pa],
                'py_coverage':ca,'wl_coverage':cb,'normalization':'face trace divided by absolute delta-argument slope in w'}
    return compare_values(a,b,row,part,options,conditions)


def has_key(obj,key):
    if isinstance(obj,dict):
        return key in obj or any(has_key(v,key) for v in obj.values())
    if isinstance(obj,(list,tuple)):
        return any(has_key(v,key) for v in obj)
    return False


class Comparator:
    def __init__(self,names,options=None):
        self.names=Names(names)
        self.options=options
        self.compared=[set(),set()]

    def component(self,p,w,row,component,conditions=None):
        if component.kind=='text':
            if p.node.head!='Text' or w.node.head!='Text':
                return {'not_formed':'text component has a different kind'}
            result={'part':component.object,'py_text':p.node.value,'wl_text':w.node.value,
                    'text_difference':text_difference(p.node.value,w.node.value)}
        else:
            try:
                a=math_node(p.node,self.names,'py')
                b=math_node(w.node,self.names,'wl')
            except Unsupported as exc:
                return {'part':component.object,'not_compared':str(exc)}
            if isinstance(a,s.Tuple) or isinstance(b,s.Tuple):
                return {'not_formed':'scalar component contains a sequence'}
            result=compare_scalar(a,b,row,component.object,self.options,conditions)
        if not has_key(result,'not_compared') and not has_key(result,'not_formed') and not result.get('py_coverage') and not result.get('wl_coverage'):
            self.compared[0].update(leaf_paths(p.node,p.path))
            self.compared[1].update(leaf_paths(w.node,w.path))
        return result

    def compare(self,p,w,row,context=None):
        if not shaped(p,row.py_shape) or not shaped(w,row.wl_shape):
            return {'not_formed':'component counts or kinds differ from the cited source schema',
                    'py_shape':row.py_shape,'wl_shape':row.wl_shape}
        parts=[]
        domain_diagnostic=None
        if context is not None:
            domain_diagnostic=compare_values(*context,row,'external printed domain contexts',self.options)
        for component in row.components:
            try:
                pc,wc=get_ref(p,component.py_path),get_ref(w,component.wl_path)
                conditions=context
                if component.kind=='conditional_value':
                    conditions=tuple(math_node(get_ref(ref,('condition',)).node,self.names,engine)
                                     for ref,engine in ((p,'py'),(w,'wl')))
            except Unsupported as exc:
                parts.append({'part':component.object,'not_formed':str(exc)})
                continue
            parts.append(self.component(pc,wc,row,component,conditions))
        # Constructor/key leaves of equations or conditioned objects are
        # credited only when all their children have actually been compared.
        # Arbitrarily omitted children cannot earn credit through this rule.
        for side,ref in enumerate((p,w)):
            if (row.py_shape if side==0 else row.wl_shape) in ('equation','conditioned'):
                childrefs=([get_ref(ref,('lhs',)),get_ref(ref,('rhs',))] if shaped(ref,'equation') else
                           [get_ref(ref,('condition',)),get_ref(ref,('value',))])
                if all(leaf_paths(r.node,r.path)<=self.compared[side] for r in childrefs):
                    self.compared[side].update(leaf_paths(ref.node,ref.path))
        result={'components':parts}
        if domain_diagnostic is not None:
            result['domain_context_diagnostic']=domain_diagnostic
            result['domain_context_census']='External validity contexts are not bound value endpoints and receive no additional value-leaf credit.'
        return result


def curated_unjoined_rationale(tag,engine):
    """Source/table curation only; not evidence about active bindings."""
    if '/projection integration/' in tag or '/closed integration/' in tag:
        return 'Construction residuals subtract the two integration routes in opposite orders (py201-205; wl161-162,15); no residual or note binding.'
    if '/coordinate ' in tag or '/pulled current' in tag or '/transverse ' in tag or tag.startswith('transverse integration'):
        return 'Wolfram-only fixed-slab coordinate construction; SymPy prints native-w integration objects (py171-185; wl64-99,108-125,143-155). Different printed objects; no coordinate reconciliation.'
    if '/bulk derivative integral' in tag:
        return 'Native-w derivative integral or its endpoints (py171-176); Wolfram prints the pulled z-current derivative and relative endpoints (wl143,151-155). Different printed integrands/current objects.'
    if '/supplied material factor' in tag or '/supplied material face factor' in tag:
        return 'SymPy prints the one-sided trace w=face (py213-216); Wolfram prints unrestricted w-dependent factors (wl48,118). Different objects; faces are not substituted to reconcile them.'
    if 'dimension' in tag or '/sum ' in tag:
        return 'Engine-specific dimension audit occurrence: SymPy operand/sum/reference expression trees (py305-313,316-337), or Wolfram assignments/atoms/derivative/source/term ledgers and extra shell objects (wl244-256). No unused source-cited counterpart; term identity and position do not create bindings.'
    if tag.endswith('/projected velocity') and engine=='py':
        return 'SymPy prints velocity at free radius (py271,275); Wolfram prints it only at outer radius ro (wl218-219). The outer-radius conditional counterpart is bound separately.'
    if tag.endswith('/regular parts'):
        return 'Regional differential objects are already compared through supplied equation lhs components (py120-121,163,167; wl130,136); this redundant standalone occurrence has no unused counterpart.'
    if '/regional balance residuals' in tag or '/source law residual' in tag:
        return 'Wolfram construction checks subtract native supplied rhs (wl137,178); no SymPy check against that same counterpart is printed.'
    if '/encoded ' in tag:
        return 'Only SymPy prints the encoded global current before differentiation (py157-158,164-165); Wolfram prints regional currents and their differentiated distribution (wl129-138).'
    if tag.endswith('/projected differential expansion'):
        return 'Only SymPy prints the standalone projected differential expansion (py181,192); Wolfram prints it combined with face terms (wl158-170).'
    if tag.endswith('/face channel by face'):
        return 'Only SymPy prints shell-integrated transport by individual face (py251,263); Wolfram shell objects group interior/concentrated contributions across both faces (wl201-203,222).'
    if 'face channel minus order lost/operands' in tag:
        return 'Only SymPy prints the two operands separately (py270,33-34); Wolfram prints their difference (wl215,222), which is bound to the SymPy residual.'
    if tag.endswith('/projected source'):
        return 'Material source literal zero is already compared in the projected balance rhs (py188,198-200; wl163,173); this redundant occurrence has no unused counterpart.'
    if '/supplied/' in tag or '/supplied material face equation/' in tag or tag.endswith('/velocity'):
        return 'SymPy redundant supplied-field/equation occurrence (py109-148,197); corresponding current and regional equation components are bound where both are printed (wl129-136,172). No unused counterpart for this occurrence.'
    return 'No source-cited end for this tag is declared in the reference table; no semantic counterpart absence or completeness claim is made.'


def unjoined_accounting(entry,engine,streams,bindings,bound_tags,used):
    """Separate curated meaning from computed current-table accounting.

    Reference candidates are declared ends, never inferred from payloads or
    similar names. Their occurrence counts and active binding membership are
    obtained from this run, including when component bindings are removed.
    """
    tag=entry.tag
    other='wl' if engine=='py' else 'py'
    candidates=[]
    for curated in BINDING_TABLE:
        if getattr(curated,engine).tag!=tag:
            continue
        counterpart=getattr(curated,other)
        candidates.append({'curated_binding':curated.key,
            'local_end':asdict(getattr(curated,engine)),
            'counterpart_engine':other,'counterpart_end':asdict(counterpart),
            'py_line':curated.py_line,'wl_line':curated.wl_line,
            'counterpart_occurrences':len(streams[other].index.get(counterpart.tag,[])),
            'active_binding_rows':[r.key for r in bindings if r.py==curated.py and r.wl==curated.wl],
            'counterpart_tag_accounted_bindings':list(bound_tags[other][counterpart.tag]),
            'counterpart_tag_compared_leaves':len(used[other][counterpart.tag])})
    if is_comparison_note(tag):
        result={'reason':'Engine-local comparison note; this occurrence has no unique join under the active table.',
                'classification':'engine-local comparison note','read_as_label':is_label(tag),
                'reason_basis':'reader label classification and current-run accounting; py27-36,201-205; wl14-15,71,146,162'}
    elif candidates:
        result={'reason':'No unique join formed under the active table; curated candidate endpoints and computed current-run facts are listed separately.',
                'reason_basis':'candidate meanings rest on reference-table curation; occurrence counts and active rows are computed'}
    else:
        result={'reason':'Reference-table/source curation (not computed active-table or counterpart-completeness evidence): '+curated_unjoined_rationale(tag,engine),
                'reason_basis':'source/table curation only; completeness of named-object comparison is not computed'}
    if candidates:
        result['curated_binding_candidates']=candidates
    return result


def unpaired_subtrees(ref,covered,reason='No unused source-cited component correspondence; the subtree is retained, not counted as compared.'):
    paths=leaf_paths(ref.node,ref.path)
    if paths<=covered:
        return []
    if not paths & covered:
        return [{'path':list(ref.path),'parse':data(ref.node),'parsed_leaves':census(ref.node),
                 'reason':reason}]
    result=[]
    for i in range(len(ref.node.args)):
        result.extend(unpaired_subtrees(ref.child(i),covered,reason))
    return result


def run(py_path,wl_path,out,*,names=NAME_TABLE,bindings=BINDING_TABLE,progress=True,options=None,probe_metadata=True):
    global ACTIVE_OUTPUT
    started=time.monotonic()
    validate_tables(names,bindings)
    streams={engine:Stream(path,engine) for engine,path in (('py',py_path),('wl',wl_path))}
    cmp=Comparator(names,options)
    emit(out,{'kind':'configuration','sources':{'py':'scripts/S14a_1_projection.py at 94fddf2f','wl':'mathematica/S14a_1_projection.wl at 94fddf2f'},
              'input_sha256':{e:st.sha256 for e,st in streams.items()},'name_table':[asdict(r) for r in names],
              'binding_table':[asdict(r) for r in bindings],
              'root_schemas':[{'engine':e,'tag':tag,'shape':shape} for (e,tag),shape in ROOT_SCHEMAS.items()], 'forms':FORMS,'strata':STRATA,
              'base_points':[[str(x) for x in p] for p in BASE_POINTS],
              'domain':'real coordinates; strict face ordering; positive spherical radius; 0<inner<outer; 0<theta<pi',
              'residual_target':'No residual target was given.'})
    # Table rows and original source references appear before any comparison.
    for engine,stream in streams.items():
        for entry in stream.progress:
            emit(out,{'kind':'progress_marker','engine':engine,'line':entry.line,'tag':entry.tag,'raw':entry.raw})
    if probe_metadata:
        for form in FORMS:
            for stratum in STRATA:
                emit(out,{'kind':'test_functions','form':form,'stratum':stratum,'formal_arguments':list(map(str,FORMAL)),
                          'functions':{r.py:str(profile(r.py,r.arity,form,stratum)) for r in names if r.arity}})
                for chart in ('General','SphericalSpecialization'):
                    for index in range(3):
                        for loc in LOCATIONS:
                            env=point_for(chart,form,index,loc)
                            lo,hi=(q.xreplace(env) for q in face_forms(chart,form))
                            emit(out,{'kind':'probe_domain','chart':chart,'form':form,'stratum':stratum,'point_index':index,'location':loc,
                                      'point':{str(k):str(v) for k,v in env.items() if str(k) in SCALAR_NAMES},
                                      'face_order':str(s.Lt(lo,hi)),'shell_order':str(s.And(s.Gt(env[A],0),s.Lt(env[A],env[B]))),
                                      'radius_positive':str(s.Gt(env[R],0))})
    total=Counter()
    used={e:defaultdict(set) for e in streams}
    bound_tags={e:defaultdict(list) for e in streams}
    below=[]
    for row in bindings:
        if progress:
            print('WORKING_ON: '+row.key+' '+row.object,file=sys.stderr,flush=True)
        entries={e:streams[e].index.get(getattr(row,e).tag,[]) for e in streams}
        if any(len(v)!=1 for v in entries.values()):
            original={}
            for e,values in entries.items():
                payloads=[{'tag':v.tag,'line':v.line,'raw':v.raw,'parse':data(v.node)} for v in values]
                original[e]=payloads[0] if len(payloads)==1 else {'occurrences':payloads}
            emit(out,{'kind':'operands','binding':row.key,**original})
            emit(out,{'kind':'comparison','binding':row.key,'not_compared':'binding end missing or repeated',
                      'occurrences':{e:len(v) for e,v in entries.items()}})
            parsed={}
            for e,values in entries.items():
                count=0
                for entry in values:
                    try:
                        ref=get_ref(Ref(entry.node),getattr(row,e).path)
                    except Unsupported:
                        ref=Ref(entry.node)
                    count+=census(ref.node)
                parsed[e]=count
                total['bound_parsed_leaves_'+e]+=count
            emit(out,{'kind':'binding_accounting','binding':row.key,'parsed_leaves':parsed,
                      'actually_compared_leaves':{e:0 for e in streams},
                      'parse_scope':'all parsed occurrences at each rejected end'})
            below.append(row.key)
            total['bindings']+=1
            total['bindings_not_compared']+=1
            continue
        emit(out,{'kind':'operands','binding':row.key,**{e:{'tag':v[0].tag,'line':v[0].line,'raw':v[0].raw,'parse':data(v[0].node)} for e,v in entries.items()}})
        refs={}
        for e in streams:
            bound_tags[e][entries[e][0].tag].append(row.key)
        try:
            for e,v in entries.items():
                schema=ROOT_SCHEMAS.get((e,v[0].tag))
                if schema is not None and not shaped(Ref(v[0].node),schema):
                    raise Unsupported('full tagged payload differs from cited container schema: '+e+' '+v[0].tag)
            refs={e:get_ref(Ref(v[0].node),getattr(row,e).path) for e,v in entries.items()}
        except Unsupported as exc:
            emit(out,{'kind':'comparison','binding':row.key,'not_formed':str(exc)})
            parse_refs={}
            for e,v in entries.items():
                try:
                    parse_refs[e]=get_ref(Ref(v[0].node),getattr(row,e).path)
                except Unsupported:
                    parse_refs[e]=Ref(v[0].node)
            parsed={e:census(ref.node) for e,ref in parse_refs.items()}
            emit(out,{'kind':'binding_accounting','binding':row.key,'parsed_leaves':parsed,
                      'actually_compared_leaves':{e:0 for e in streams},
                      'parse_scope':'selected end if present, otherwise full tagged parse; no comparison formed'})
            below.append(row.key)
            for e,count in parsed.items():
                total['bound_parsed_leaves_'+e]+=count
            total['bindings']+=1
            total['bindings_not_formed']+=1
            continue
        cmp.compared=[set(),set()]
        ACTIVE_OUTPUT=out
        try:
            context=None
            if row.condition_context:
                domain_values=[]
                for engine,end in zip(('py','wl'),row.condition_context):
                    if end is None:
                        domain_values.append(s.true)
                    else:
                        candidates=streams[engine].index.get(end.tag,[])
                        if len(candidates)!=1:
                            raise InputError('missing or repeated explicit domain context '+end.tag)
                        domain_values.append(math_node(get_ref(Ref(candidates[0].node),end.path).node,cmp.names,engine))
                context=tuple(domain_values)
            comparison=cmp.compare(refs['py'],refs['wl'],row,context)
        finally:
            ACTIVE_OUTPUT=None
        emit(out,{'kind':'comparison','binding':row.key,'comparison':comparison})
        parsed={e:census(ref.node) for e,ref in refs.items()}
        actual={e:len(cmp.compared[i]) for i,e in enumerate(streams)}
        missing={e:sorted([list(p) for p in leaf_paths(ref.node,ref.path)-cmp.compared[i]],key=str) for i,(e,ref) in enumerate(refs.items())}
        emit(out,{'kind':'binding_accounting','binding':row.key,'parsed_leaves':parsed,'actually_compared_leaves':actual,
                  'uncompared_leaf_paths':missing})
        for i,e in enumerate(streams):
            used[e][entries[e][0].tag].update(cmp.compared[i])
            total['bound_parsed_leaves_'+e]+=parsed[e]
            total['compared_leaves_'+e]+=actual[e]
        if any(actual[e]<parsed[e] for e in parsed):
            below.append(row.key)
        total['bindings']+=1
        total['bindings_not_formed']+=has_key(comparison,'not_formed')
        total['bindings_not_compared']+=has_key(comparison,'not_compared')
    for e,stream in streams.items():
        total['tags_'+e]=sum(len(v) for v in stream.index.values())
        total['lines_'+e]=stream.line_count
        for tag,entries in stream.index.items():
            for occurrence,entry in enumerate(entries):
                keys=bound_tags[e][tag] if len(entries)==1 else []
                leaves=census(entry.node)
                actual=len(used[e][tag]) if len(entries)==1 else 0
                record={'kind':'tag_accounting','engine':e,'tag':tag,'line':entry.line,'occurrence':occurrence,
                        'parsed_leaves':leaves,'actually_compared_leaves':actual,'accounting':'bound' if keys else 'unjoined'}
                if keys:
                    reason='No unused source-cited component correspondence; retained without value-leaf credit.'
                    if e=='wl' and tag.startswith('spherical specialization / dimensions/projected velocity'):
                        reason='Conditional wrapper/domain context (wl238,256): domain diagnostics constrain dimension values; repeated outer-density condition has no unused same-point SymPy predicate occurrence (py274-277,305-308,336-337). Constructor heads are syntax, not standalone values.'
                    elif e=='py' and tag.endswith('/supplied material face equation/operands'):
                        reason='Supplied equation rhs zero (py197,33-34); wl172 prints only the differential coefficient, not this rhs component.'
                    elif e=='py' and tag.endswith('/dimensions/volume flux/dimension operand'):
                        reason='Py305-306,337 prints the native volume-flux expression; wl238,256 prints its dimension vector and condition. Only the predicate is a corresponding component here.'
                    record.update(bindings=keys,unpaired_components=unpaired_subtrees(Ref(entry.node),used[e][tag],reason))
                else:
                    if len(entries)>1:
                        record.update(reason='Repeated tag; binding rejected.',reason_basis='computed occurrence count')
                    else:
                        record.update(unjoined_accounting(entry,e,streams,bindings,bound_tags,used))
                    record.update(raw=entry.raw,parse=data(entry.node))
                    total['unjoined_tags_'+e]+=1
                total['parsed_leaves_'+e]+=leaves
                emit(out,record)
    emit(out,{'kind':'resources','runtime_seconds':time.monotonic()-started,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    emit(out,{'kind':'summary','counts':dict(total),'bindings_below_parsed_count':below,
              'exit_semantics':'zero for any findings; nonzero only for operational input failure'})
    return dict(total)


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
    raise SystemExit(main())
