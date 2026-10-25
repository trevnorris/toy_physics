#!/usr/bin/env python3
"""S15 cross-engine comparator. No residual target was given.

Native readers, derivative-slot translation, leaf census, and integral evaluation
adapted from S16 at f08db34a. S15 bindings cite engines at 3192123b and the
S15 decision list revision 3 at 91248b67. Neither engine is imported or run.
"""
from __future__ import annotations
import argparse
import ast
from collections import Counter, defaultdict
from dataclasses import dataclass, asdict
from functools import lru_cache
import hashlib
import json
import multiprocessing
import gc
import os
import threading
from contextlib import contextmanager
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
            if head in ('StrictLessThan','LessThan','StrictGreaterThan','GreaterThan') and (left.head in ('StrictLessThan','LessThan','StrictGreaterThan','GreaterThan') or (left.head=='And' and left.args[-1].head in ('StrictLessThan','LessThan','StrictGreaterThan','GreaterThan'))):
                left = Node('And', (left, Node(head, ((left.args[-1].args[-1] if left.head=='And' else left.args[-1]), right))))
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
    memo={}
    def walk(n): return memo[id(n)]
    def arithmetic(head,args):
        return Node(head,tuple(q for a in args for q in (a.args if a.head==head and head in ("Add","Mul") else (a,))))
    def convert(n):
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
                return arithmetic('Add', (a, arithmetic('Mul', (atom('Number', -1), b))))
            if isinstance(n.op, ast.Div):
                return Node('Mul', (a, Node('Pow', (b, atom('Number', -1)))))
            heads = {ast.Add:'Add', ast.Mult:'Mul', ast.Pow:'Pow', ast.BitAnd:'And', ast.BitOr:'Or'}
            if type(n.op) in heads:
                return arithmetic(heads[type(n.op)], (a, b))
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
        root=ast.parse(source, mode='eval').body
        # Python's AST stores long sums as a left-deep binary tree. Convert
        # children before parents without consuming the Python call stack;
        # associative arithmetic is retained as a flat, unevaluated Node.
        for node in reversed(list(ast.walk(root))):
            if isinstance(node,ast.expr): memo[id(node)]=convert(node)
        return memo[id(root)]
    except (SyntaxError, ValueError, KeyError) as exc:
        raise InputError('SymPy grammar: '+str(exc)) from exc


from itertools import product, combinations

@dataclass(frozen=True)
class NameRow:
    py: str
    wl: str
    object: str
    py_line: str
    wl_line: str
    arity: int = 0
    alias: str = ''

@dataclass(frozen=True)
class End:
    tag: str
    path: tuple = ()

@dataclass(frozen=True)
class Binding:
    key: str
    rule: str
    py: End
    wl: End
    object: str
    py_line: str
    wl_line: str
    kind: str = 'value'


def name_contract():
    rows=[]
    def add(p,w,obj,pl,wl,n=0): rows.append(NameRow(p,w,obj,pl,wl,n))
    for p,w in zip(('x_1','x_2','x_3','t','s','theta','phi'),('x','y','z','t','s','theta','phi')):
        add(p,w,'supplied coordinates / item 1-2 chart variable '+p,'37-41,129-141','18,75-86')
    for k,i in product((1,2),(1,2,3)):
        add(f'X_{k}_{i}',f'xx[{k},{i}]',f'drains: centre {k}, Cartesian component {i}','48','20',1)
    for p,w,obj in [('r_in','rin','inner'),('R','rout','outer')]:
        add(p,w,'region: '+obj+' radius','42-43,231','38,156',1)
    for k in (1,2):
        for p,w,pl,wl in [('rho','rr','50,52','23'),('u','uu','51,53','22')]:
            add(f'{p}_{k}',f'{w}[{k}]',f'drains: {p} radial profile of drain {k}',pl,wl,2)
    add('varrho','vr','drains: remainder density','55','25',4)
    for i in (1,2,3):
        add(f'W_{i}',f'ww[{i}]',f'drains: remainder velocity component {i}','56','24',4)
        add(f'f_{i}',f'ff[{i}]',f'supplied force density component {i}','62','30',4)
    for i,j in product((1,2,3),repeat=2):
        add(f'Pi_o_{i}_{j}',f'oo[{i},{j}]',f'supplied other current ordered indices {i},{j}','61','29',4)
    for p in ('eta','zeta'): add(p,p,'supplied P6 viscosity '+p,'60','28',4)
    add('p_br','pb','supplied P6 pressure as function of density','59','33',1)
    add('c_comp','ccomp','free parameter outside supplied balance','76,170','45,105',4)
    for p,w in zip(('label_1','label_2','label_r'),('a1','a2','ar')):
        add(p,w,'dimensionless label marker '+p,'44-47','18-19')
    for p,w in zip(('Density','Length','Time'),('Dens','Length','Time')):
        add(p,w,'item 3 dimension basis '+p,'145-148,224','94,118')
    # Bound dummies are alpha-renamed lexically by Subs/Derivative/Integral;
    # they are not free-name substitutions. No printed term aliases occur in
    # these sources. The generic printed-definition resolver is retained for
    # synthetic alias controls, without inventing aliases in these streams.
    return tuple(rows)

NAME_TABLE=name_contract()
@dataclass(frozen=True)
class SlotRow:
    py_index: int
    wl_index: int
    py: str
    wl: str
    object: str
    py_line: str
    wl_line: str

# These names are lexical, not additional free fields. SymPy prints _xi_j
# inside a Subs scope; Wolfram prints the derivative slot instead of naming a
# dummy. Alpha-renaming a binder does not change the bound point or derivative.
SLOT_TABLE=tuple(SlotRow(j,j,f'_xi_{j+1} in its lexical Subs scope',
                        f'Derivative order slot {j+1}',f'formal function argument {j+1}',
                        '63,79,116,137,284','32,36,63,84,88-90') for j in range(4))

SURFACE_TABLE=({'object':'inner sphere about X_2, radius r_in, inward normal',
                'py':'inner','wl':'inner','py_line':'231-239','wl_line':'125-131,156'},
               {'object':'outer sphere about X_2, radius R, outward normal',
                'py':'outer','wl':'outer','py_line':'231-239','wl_line':'125-131,156'})
GROUPS=tuple(c for n in (1,2,3) for c in combinations(('drain1','drain2','remainder'),n))
PIECES=('carried','pressure','viscous','other','transport_centre','transport_radius')
ROOT_SHAPES={}
ROOT_KINDS={}
# One declaration per binding rule; repeated parameter instances are not new rules.
RULE_NAMES=('balance_dimensions', 'balance_reduction', 'comparison_notes', 'current_partition', 'dimension_basis', 'field_signatures', 'free_parameter', 'global_domain', 'global_domain_label', 'homogeneity_dimension_operands', 'homogeneity_dimension_pair_residuals', 'local_balance', 'local_balance_equation', 'pressure_condition', 'pressure_partition', 'pressure_pressure_expansion_at_order_one', 'pressure_pressure_whole', 'pressure_pressure_whole_minus_expansion', 'prose', 'region_balance', 'region_balance_equation', 'region_balance_parts', 'region_chart', 'region_group_partition', 'region_integrated_force', 'region_momentum', 'region_momentum_groups', 'region_momentum_rate', 'region_rate_groups', 'region_volume_measure', 'supplied_P5_G', 'supplied_P6_T', 'supplied_Pi', 'supplied_density_profiles_and_remainder', 'supplied_distances', 'supplied_total_density', 'supplied_total_velocity', 'supplied_velocity_profiles_and_remainder', 'surface_area_measure', 'surface_centre_velocity', 'surface_chart', 'surface_group_partition', 'surface_groups', 'surface_outward_normal', 'surface_piece', 'surface_pressure_condition', 'surface_pressure_first_order_choice', 'surface_pressure_groups', 'surface_pressure_whole_minus_expansion', 'surface_radius_velocity', 'surface_surface_velocity', 'surface_velocity_partition')
RULES={}
VEC=('v',)*3
MAT=(VEC,)*3

def shape_paths(shape,path=()):
    if isinstance(shape,str) and shape in ('v','predicate'): yield path
    elif isinstance(shape,dict):
        for k,v in shape.items(): yield from shape_paths(v,path+(k,))
    else:
        for i,v in enumerate(shape): yield from shape_paths(v,path+(i,))


def binding_contract():
    rows=[]
    def add(rule,p,w,obj,pl,wl,kind='value',pp=(),wp=()):
        if rule not in RULE_NAMES: raise InputError('undeclared binding rule '+rule)
        RULES.setdefault(rule,{'object':obj,'py_line':pl,'wl_line':wl})
        rows.append(Binding(f'{rule}:{len(rows)}',rule,End(p,pp),End(w,wp),obj,pl,wl,kind))
    def obj(rule,p,w,description,pl,wl,shape='v',kind='value',note=False,dim=False):
        if isinstance(shape,dict) and set(shape)=={'lhs','rhs'}:
            ROOT_KINDS['py',p]='equation'; ROOT_KINDS['wl',w]='equation'
        if shape!='v':
            ROOT_SHAPES['py',p]=shape; ROOT_SHAPES['wl',w]=shape
        for path in shape_paths(shape):
            add(rule,p,w,description+' component '+str(path),pl,wl,kind,path,path)
        if note: add('comparison_notes',p+'/comparison',w+'/comparison',description+' engine-local absence-of-counterpart note','26-28,'+pl,'13,'+wl,'label')
        if dim:
            for suffix in ('dimension_operands','dimension_pair_residuals'):
                add('homogeneity_'+suffix,p+'/'+suffix,w+'/'+suffix,'item 3 '+description+' '+suffix,'220-223,'+pl,'120,'+wl,'dimension')
    for tag in ('PRESSURE_DOMAIN','REGULARITY','PREMISES','CHOICES','BOOKKEEPING','balance_construction','scope','END_OF_TRANSCRIPT','dimension_convention'):
        add('prose',tag,tag,'domain / supplied premises / scope: '+tag,'71-75,226,311-313','40-44,119,177-179','label')
    ROOT_SHAPES['py','DOMAIN_FOR_EVERY_OBJECT']=('predicate','v'); ROOT_SHAPES['wl','DOMAIN_FOR_EVERY_OBJECT']=('predicate','v')
    add('global_domain','DOMAIN_FOR_EVERY_OBJECT','DOMAIN_FOR_EVERY_OBJECT','region domain relation','68-70','37-39','condition',(0,),(0,))
    add('global_domain_label','DOMAIN_FOR_EVERY_OBJECT','DOMAIN_FOR_EVERY_OBJECT','region domain prose','70','39','label',(1,),(1,))
    add('dimension_basis','dimension_basis','dimension_basis','item 3 printed dimension basis','224','118','basis')
    obj('free_parameter','c_comp_not_present_in_supplied_equations','LIVE_PARAMETER_NOT_IN_SUPPLIED_BALANCE','free c_comp parameter','76','45')
    for p,w,shape,desc in [
        ('distances','radial_distances',('v','v'),'drains radial distances'),
        ('density_profiles_and_remainder','drain_and_remainder_densities',VEC,'drains density profiles and remainder'),
        ('velocity_profiles_and_remainder','drain_and_remainder_velocities',(VEC,)*3,'drains velocity profiles and remainder'),
        ('total_density','total_density','v','drains total density'),
        ('total_velocity','total_velocity',VEC,'drains total velocity'),
        ('P5_G','P5_momentum_density',VEC,'supplied P5 momentum density'),
        ('P6_T','P6_stress',MAT,'supplied P6 stress'),
        ('Pi','momentum_current',MAT,'supplied total current')]:
        obj('supplied_'+p,p,w,desc,'48-66,77-78','20-35,46-48',shape,note=True)
    def check(rule,p,w,desc,pl,wl,shape='v'):
        for suffix in ('operand_a','operand_b','residual'):
            obj(rule,p+'/'+suffix,w+'/'+suffix,desc+' '+suffix,pl,wl,shape)
    check('local_balance','supplied_local_balance','supplied_local_balance','supplied local balance: time derivative of G plus divergence of Pi versus f','79-82','36,49-50',VEC)
    obj('local_balance_equation','supplied_local_balance/equation','supplied_local_balance/equation','supplied local balance equation','82','50',{'lhs':VEC,'rhs':VEC})
    check('current_partition','Pi_channel_partition','current_channel_partition','supplied Pi versus sum of extracted current channels with markers','95-99','53-57',MAT)
    for p,w in [('pressure_whole','pressure_whole'),('pressure_expansion_at_order_one','pressure_expansion'),('pressure_whole_minus_expansion','pressure_whole_minus_expansion')]:
        obj('pressure_'+p,p,w,'item 1 '+p,'111-125','63-71',note=True)
    add('pressure_condition','pressure_expansion/order_and_condition','pressure_expansion/order_and_condition','item 1 expansion order and condition','124','70','label')
    check('pressure_partition','pressure_expansion_partition','pressure_group_partition','item 1 pressure expansion versus sum of support groups','126','72')
    for surface in ('inner','outer'):
        for key in ('chart','outward_normal','surface_velocity','centre_velocity','radius_velocity','area_measure'):
            obj('surface_'+key,surface+'/'+key,surface+'/'+key,'item 1 '+surface+' '+key,'231-239','125-131,156', 'v' if key=='area_measure' else VEC,note=True)
        check('surface_velocity_partition',surface+'/surface_velocity_partition',surface+'/velocity_partition','item 1 '+surface+' total surface velocity versus centre plus radius velocities','240','132',VEC)
        for piece,i in product(PIECES,(1,2,3)):
            tag=f'{surface}/{piece}/component_{i}'
            obj('surface_piece',tag,tag,f'item 1 {surface} {piece}, Cartesian component {i}','241-250','133-139',note=True,dim=True)
            if piece not in ('pressure','other'):
                for labels in GROUPS:
                    p=tag+'/labels_'+'+'.join(labels); w=tag+'/labels_{'+', '.join(labels)+'}'
                    obj('surface_groups',p,w,f'item 1 {surface} {piece}, component {i}, exact label set {labels}','252-258','140-143',note=True,dim=True)
                check('surface_group_partition',tag+'/group_partition_integrand',tag+'/group_partition_integrand',f'item 1 {surface} {piece} component {i}, whole integrand versus sum of support groups','258','143')
            if piece=='pressure':
                for suffix in ('first_order_choice','whole_minus_expansion'):
                    obj('surface_pressure_'+suffix,tag+'/'+suffix,tag+'/'+suffix,f'item 1 {surface} pressure {suffix} component {i}','259-268','145-150',note=True,dim=True)
                add('surface_pressure_condition',tag+'/first_order_domain',tag+'/first_order_domain',f'item 1 {surface} pressure condition component {i}','265','149','label')
                for labels in GROUPS:
                    obj('surface_pressure_groups',tag+'/expansion_labels_'+'+'.join(labels),tag+'/expansion_labels_{'+', '.join(labels)+'}',f'item 1 {surface} pressure expansion component {i}, exact label set {labels}','269-275','151-153',note=True,dim=True)
    for key in ('chart','volume_measure'):
        obj('region_'+key,'region/'+key,'region/'+key,'item 2 region '+key,'135,279','82-83,159',VEC if key=='chart' else 'v',note=True)
    for i in (1,2,3):
        tag=f'region/component_{i}'
        for suffix in ('momentum','momentum_rate','integrated_force'):
            obj('region_'+suffix,tag+'/'+suffix,tag+'/'+suffix,f'item 2 {suffix} component {i}','282-287','162-164,168',note=True,dim=True)
        for labels,suffix in product(GROUPS,('momentum','rate')):
            obj('region_'+suffix+'_groups',tag+'/'+suffix+'_labels_'+'+'.join(labels),tag+'/'+suffix+'_labels_{'+', '.join(labels)+'}',f'item 2 {suffix} component {i}, exact label set {labels}','288-292','165-166',note=True)
        check('region_group_partition',tag+'/momentum_group_partition_integrand',tag+'/momentum_group_partition_integrand',f'item 2 momentum integrand component {i}, whole versus sum of support groups','293','167')
        obj('region_balance_parts',tag+'/balance_parts',tag+'/balance_parts',f'item 2 inner and outer surface pieces and undivided force component {i}','294-296','170',({p:'v' for p in PIECES},{p:'v' for p in PIECES},'v'))
        check('region_balance',tag+'/balance',tag+'/balance',f'item 2 momentum rate versus sum of both surfaces and integrated force component {i}','295-297','169-171')
        obj('region_balance_equation',tag+'/balance_equation',tag+'/balance_equation',f'item 2 region balance equation component {i}','299','171',{'lhs':'v','rhs':'v'})
        add('balance_reduction',tag+'/balance/reduction',tag+'/balance/reduction',f'item 2 component {i} engine reduction label','298','172','label')
        for suffix in ('balance_dimension_operand_a','balance_dimension_operand_b','balance_dimension_residuals','balance_rhs/dimension_operands','balance_rhs/dimension_pair_residuals'):
            add('balance_dimensions',tag+'/'+suffix,tag+'/'+suffix,f'item 3 region balance component {i} '+suffix,'300-304','173-175','dimension')
    # Split-vs-grouped packaging: both field output and argument dimensions are
    # compared as individual triples, joined through the explicit name table.
    for name in NAME_TABLE:
        if not name.arity: continue
        p='dimension_field_signatures'
        for component,wt in ((0,'field_dimension_signatures'),(1,'field_argument_dimension_signatures')):
            paths=[()] if component==0 else [(j,) for j in range(name.arity)]
            for suffix in paths:
                add('field_signatures',p,wt,'item 3 field signature '+name.object,'154-170,225','96-105,118','dimension_field',(name.py,component)+suffix,(name.wl,)+suffix)
    return tuple(rows)

BINDING_TABLE=binding_contract()

def validate_tables(names,bindings,rules=None):
    rules=RULE_NAMES if rules is None else rules
    for field in ('py_index','wl_index'):
        slots=[getattr(r,field) for r in SLOT_TABLE]
        if len(set(slots))!=len(slots): raise InputError('duplicate bound-name row '+field)
    for field in ('py','wl'):
        vals=[getattr(n,field) for n in names]
        if len(set(vals))!=len(vals): raise InputError('duplicate name row '+field)
        seen=defaultdict(list)
        for row in bindings:
            end=getattr(row,field)
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
    if n.head in ('Symbol','Text','Number'): return re.sub(r'\s+', '', n.value) if re.fullmatch(r'[a-z]+\[.*\]',n.value) else n.value
    if n.head=='Apply': return key_text(n.args[0])+'['+','.join(key_text(q) for q in n.args[1:])+']'
    raise InputError('nonliteral dictionary key')

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
        if any(x.head!='Sequence' for x in rows.node.args): raise Unsupported('matrix row layout')
        column=all(len(x.args)==1 for x in rows.node.args)
        return {i:rows.child(i).child(0) if column else rows.child(i) for i in range(len(rows.node.args))}
    if n.head in ('Rule','Equality'): return {'lhs':ref.child(0),'rhs':ref.child(1)}
    if call(n,'Eq') or call(n,'Equal'): return {'lhs':ref.child(1),'rhs':ref.child(2)}
    return None

def get_ref(ref,path):
    for k in path:
        children=semantic_children(ref)
        if children is None or k not in children: raise Unsupported('no semantic component '+str(k)+' at '+str(ref.path))
        ref=children[k]
    return ref

def shaped(ref,shape):
    if shape=='predicate':
        return ref.node.head in ('Boolean','And','Or','Not','Equality','Unequality','StrictLessThan','LessThan','StrictGreaterThan','GreaterThan') or any(call(ref.node,name) for name in ('Eq','Equal','Ne','Unequal','And','Or','Not','Lt','Le','Gt','Ge'))
    children=semantic_children(ref)
    if shape=='v': return children is None
    expected=shape if isinstance(shape,dict) else dict(enumerate(shape))
    return children is not None and set(children)==set(expected) and all(shaped(children[k],v) for k,v in expected.items())



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
    slots=points if direct else s.symbols('_s15partial0:'+str(len(points)),**{})
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
    if h=='Text' and n.value in b.maps[engine]:
        return s.Symbol(b.name(n.value,engine))
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
            mapped_orders=[0]*len(args)
            for row in SLOT_TABLE[:len(args)]:
                if row.py_index>=len(args) or row.wl_index>=len(args): raise Unsupported('bound-name slot outside derivative arity')
                mapped_orders[row.py_index]=int(cv(orders[row.wl_index]))
            return partial_derivative(s.Function(key,**{}),tuple(mapped_orders),tuple(cv(q) for q in args))
        if head.head=='Prime':
            order=0
            while head.head=='Prime': order+=1; head=head.args[0]
            key=b.name(key_text(head),engine)
            if b.arities.get(key)!=1 or len(a)!=2: raise Unsupported('prime derivative arity '+key)
            return partial_derivative(s.Function(key),(order,),(cv(a[1]),))
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



TEXT_TAGS={'PRESSURE_DOMAIN','REGULARITY','PREMISES','CHOICES','BOOKKEEPING','chart_domain',
           'dimension_convention','balance_construction','scope','END_OF_TRANSCRIPT',
           'pressure_expansion/order_and_condition'}

def parse_payload(text,engine,tag):
    if tag=='WORKING_ON':
        return wl_parse(text) if engine=='wl' and text.startswith('"') else atom('Text',text)
    if engine=='py' and (tag in TEXT_TAGS or tag.endswith(('/comparison','/first_order_domain','/balance/reduction'))):
        return atom('Text',text)
    return wl_parse(text) if engine=='wl' else py_parse(text)

@dataclass
class Entry:
    tag: str
    raw: str
    line: int
    node: Node
    original_tag: str = ''

class Stream:
    def __init__(self,path,engine,progress=False):
        self.path=Path(path); self.engine=engine; self.index=defaultdict(list); self.progress=[]
        try: raw=self.path.read_bytes(); lines=raw.decode().splitlines()
        except (OSError,UnicodeError) as exc: raise InputError(str(exc)) from exc
        self.sha256=hashlib.sha256(raw).hexdigest(); self.line_count=len(lines)
        framed=[]
        tag_pattern=r'^([A-Za-z_]([^ :{]|\{[^}]*\})*): (.*)$'
        for i,line in enumerate(lines,1):
            match=re.fullmatch(tag_pattern,line)
            if match: framed.append([match[1],match[3],i])
            elif framed and framed[-1][0]!='WORKING_ON' and line.startswith((' ','[',']',')','}')):
                framed[-1][1]+='\n'+line
            else: raise InputError(f'{path}:{i}: unassigned line {line[:80]!r}')
        if not framed: raise InputError('empty stream')
        record_key=None; partition_context=None
        for tag,txt,line in framed:
            if progress: print('WORKING_ON: reading '+engine+' '+tag+' line '+str(line),file=sys.stderr,flush=True)
            try: node=parse_payload(txt,engine,tag)
            except InputError as exc: raise InputError(f'{path}:{line}:{tag}: {exc}') from exc
            original=tag
            if engine=='py' and tag=='dimension_argument/name':
                record_key=json.dumps(data(node),separators=(',',':'))
            if engine=='py' and tag.startswith('dimension_argument/'):
                if record_key is None: raise InputError('dimension record without identity')
                tag=tag+'@'+record_key
            elif engine=='py' and tag=='partition/unlabelled_coefficient':
                # Source groups(raw) is called immediately following the whole
                # piece's dimension audit (250,253), or region force (287,288).
                # The preceding named object identifies this occurrence; stream
                # ordinal is not an identity.
                if partition_context is None: raise InputError('partition record without source object')
                tag=tag+'@'+partition_context
            elif engine=='py' and tag.endswith('/dimension_pair_residuals'):
                partition_context=tag.removesuffix('/dimension_pair_residuals')
                if partition_context.endswith('/integrated_force'):
                    partition_context=partition_context.removesuffix('/integrated_force')+'/momentum'
            entry=Entry(tag,txt,line,node,original)
            (self.progress if tag=='WORKING_ON' else self.index[tag]).append(entry)

class Names:
    def __init__(self,rows):
        self.rows=rows
        self.maps={'py':{r.py:r.py for r in rows},'wl':{r.wl:r.py for r in rows}}
        self.arities={r.py:r.arity for r in rows}; self.aliases={'py':{},'wl':{}}
    def name(self,name,engine): return self.maps[engine].get(name,engine+'::'+name)
    def load(self,streams):
        for engine,stream in streams.items():
            for row in self.rows:
                if not row.alias: continue
                entries=stream.index.get(row.alias,[])
                if len(entries)!=1: continue
                refs=semantic_children(Ref(entries[0].node))
                if refs is None or set(refs)!={'lhs','rhs'}: raise InputError('alias definition structure '+row.alias)
                left=math_node(refs['lhs'].node,self,engine)
                right=math_node(refs['rhs'].node,self,engine)
                if left in self.aliases[engine]: raise InputError('duplicate alias definition '+str(left))
                self.aliases[engine][left]=right
    def resolve(self,expr,engine):
        definitions=self.aliases[engine]; aliases={r.py for r in self.rows if r.alias}
        def walk(v,active=frozenset()):
            if v in definitions:
                if v in active: raise Unsupported('cyclic alias '+str(v))
                return walk(definitions[v],active|{v})
            if isinstance(v,AppliedUndef) and v.func.__name__ in aliases:
                raise Unsupported('unresolved printed alias '+str(v))
            if not v.args:return v
            args=tuple(walk(q,active) for q in v.args)
            if args==v.args:return v
            return v.func(*args,evaluate=False) if isinstance(v,(s.Derivative,Relational)) else v.func(*args)
        return walk(expr)

FORMS=('mixed_quadratic','mixed_cubic','mixed_quartic')
STRATA=('generic','X2_velocity_zero','radii_velocity_zero','X2_and_radii_velocity_zero','X1_velocity_zero')
POINTS=((s.Rational(3,2),s.Rational(2,7),s.Rational(-1,4),s.Rational(1,6)),
        (s.Rational(7,4),s.Rational(3,8),s.Rational(4,11),s.Rational(-1,7)),
        (s.Integer(2),s.Rational(-1,3),s.Rational(2,5),s.Rational(2,9)))
FORMAL=s.symbols('u0:4')
CODES={r.py:i+1 for i,r in enumerate(NAME_TABLE)}
ALIAS_NAMES=set()

@lru_cache(maxsize=4096)
def profile_cached(name,arity,form,stratum,point,opts=()):
    options=dict(opts); a=FORMAL[:arity]; code=CODES.get(name,101)
    if options.get('duplicate_functions') and name in ('W_1','W_2'): code=CODES['W_1']
    degree=FORMS.index(form)+2
    def coef(slot):
        raw=hashlib.sha256(f'{code}:{form}:{slot}'.encode()).digest()
        return s.Integer(1+int.from_bytes(raw[:4],'big')%97)
    linear=sum((1000*(j+1)+coef(j))*v for j,v in enumerate(a))
    coupling=sum(coef(10+i*arity+j)*a[i]*a[j] for i in range(arity) for j in range(i,arity))
    if options.get('symmetric_arguments'): linear=code*sum(a); coupling=code*sum(a)**2
    if options.get('uncoupled_arguments'): coupling=sum(coef(10+j)*v**2 for j,v in enumerate(a))
    result=code+linear+coupling
    if degree>2: result+=sum(coef(40 if options.get('symmetric_arguments') else 40+j)*v**degree for j,v in enumerate(a))/7
    if degree==4 and not options.get('uncoupled_arguments'):
        result+=coef(45)*s.prod(1+(1 if options.get('symmetric_arguments') else j+1)*v for j,v in enumerate(a))/11
        result+=sum(coef(60 if options.get('symmetric_arguments') else 60+i*arity+j)*a[i]**2*a[j]**2 for i in range(arity) for j in range(i+1,arity))/13
    # Smooth radial forms make an exact polynomial integral feasible without
    # changing any operand. Each has coupled radial/time dependence. Velocity
    # profiles include a radial factor; no profile vanishes on either sphere.
    if name.startswith(('u_','rho_')) and arity==2:
        result=result.subs(a[0],a[0]**2)
        if name.startswith('u_'): result*=a[0]
    if arity==1 and (name.startswith('X_') or name in ('r_in','R')):
        z=a[0]-POINTS[point][3]
        slope=(code+coef(51)/101)/97
        if name.startswith('X_2_') and stratum in ('X2_velocity_zero','X2_and_radii_velocity_zero'): slope=s.S.Zero
        if name.startswith('X_1_') and stratum=='X1_velocity_zero': slope=s.S.Zero
        if name in ('r_in','R') and stratum in ('radii_velocity_zero','X2_and_radii_velocity_zero'): slope=s.S.Zero
        base=(code+coef(50)/101)/103
        if name=='r_in':base=s.Rational(3+point,4)
        elif name=='R':base=s.Rational(11+point,4)
        elif name.startswith('X_1_'):base+=10+int(name[-1])
        result=base+slope*z+(code+1+coef(52)/103)*z**2/103+(code+3+coef(53)/107)*z**3/107
        if degree==4:result+=(code+5+coef(54)/109)*z**4/109
    if options.get('vanishing_form') and name=='W_1':
        # Deliberate bad form for a coverage control, never a default.
        z=a[3]-POINTS[point][3]; result*=z
    return s.expand(result)

def profile(name,arity,form,stratum,point,options=None):
    return profile_cached(name,arity,form,stratum,point,tuple(sorted((options or {}).items())))

@lru_cache(maxsize=8192)
def substituted(expr,form,stratum,point,opts=()):
    options=dict(opts)
    @lru_cache(maxsize=None)
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

def sparse_moments(expr):
    """Exact sparse rational algebra over all integration variables at once.

    Avoid expanding a trigonometric coefficient as a giant SymPy expression
    during one-variable cancellation. Coefficients remain rational numbers;
    monomials are integrated independently with the same exact moments.
    This path accepts only a polynomial after exact division. Other integrands
    continue through the general exact/quadrature path without loss of coverage.
    """
    from sympy.polys.rings import ring
    generators=[];moments=[]
    for var,lo,hi in expr.limits:
        if lo==0 and hi in (s.pi,2*s.pi):
            generators.extend((s.sin(var),s.cos(var)))
            moments.append(('angle',len(generators)-2,hi==2*s.pi))
        else:
            if not (lo.is_Rational and hi.is_Rational):return None
            generators.append(var);moments.append(('power',len(generators)-1,lo,hi))
    # A radial integral can occur inside an angular integral. Those outer
    # variables are parameters here, not numerical arguments to quadrature.
    # Retain their polynomial/trigonometric dependence until the outer binder
    # is integrated. This changes neither the chart nor its evaluation point.
    integrated_count=len(generators)
    parameters=expr.function.free_symbols-set(expr.variables)
    parameter_atoms={v for v in expr.function.atoms(s.sin,s.cos)
                     if v.free_symbols and v.free_symbols<=parameters}
    parameter_atoms.update(parameters)
    if expr.function.has(s.pi):parameter_atoms.add(s.pi)
    generators.extend(sorted(parameter_atoms-set(generators),key=s.default_sort_key))
    parameter_end=len(generators)
    # Derivatives of radial profiles can leave square roots inside sums even
    # when the complete integrand is polynomial. Keep each root as an algebraic
    # generator and reduce z**2 to its base. Never choose a root sign, collapse
    # a function, or replace a derivative point. Different bases stay distinct.
    root_bases=sorted({s.expand(v.base) for v in expr.function.atoms(s.Pow)
                       if v.exp.is_Rational and v.exp.q==2},key=s.default_sort_key)
    generators.extend(s.sqrt(base) for base in root_bases)
    algebra,*variables=ring('a:'+str(len(generators)),s.QQ)
    atoms=dict(zip(generators,variables))
    root_indices={base:parameter_end+j for j,base in enumerate(root_bases)}
    @lru_cache(maxsize=None)
    def root_polynomial(base):
        n,d=convert(base)
        if d!=algebra.one or any(any(p[parameter_end:]) for p in n):
            raise s.PolynomialError('nested or non-polynomial radical base')
        return n
    def reduce_radicals(polynomial):
        if not root_bases or not any(any(v>=2 for v in p[parameter_end:]) for p in polynomial):return polynomial
        reduced=algebra.zero
        for powers,coefficient in polynomial.items():
            powers=list(powers);factor=algebra.one
            for base,index in root_indices.items():
                exponent,powers[index]=divmod(powers[index],2)
                if exponent:factor*=root_polynomial(base)**exponent
            reduced+=algebra.from_dict({tuple(powers):coefficient})*factor
        return reduced
    def quotient(numerator,denominator):
        if denominator==algebra.one:return numerator
        return numerator.exquo(denominator)
    @lru_cache(maxsize=4096)
    def convert(v):
        if v in atoms:return atoms[v],algebra.one
        if v.is_Rational:return algebra.ground_new(s.QQ.convert(v)),algebra.one
        if v.is_Add:
            terms=[convert(q) for q in v.args]
            numerator,denominator=algebra.zero,algebra.one
            for n,d in terms:
                if d==denominator:numerator+=n
                elif d==algebra.one:numerator+=n*denominator
                elif denominator==algebra.one:numerator=numerator*d+n;denominator=d
                else:
                    common=denominator.gcd(d)
                    left=quotient(d,common);right=quotient(denominator,common)
                    numerator=numerator*left+n*right;denominator*=left
            return reduce_fraction(numerator,denominator)
        if v.is_Mul:
            numerator,denominator=algebra.one,algebra.one
            for q in v.args:
                n,d=convert(q);numerator*=n;denominator*=d
            return reduce_fraction(numerator,denominator)
        if v.is_Pow and v.exp.is_Rational and v.exp.q==2:
            base=s.expand(v.base);polynomial=root_polynomial(base)
            root=variables[root_indices[base]];power=int(v.exp//1)
            return (root*polynomial**power,algebra.one) if power>=0 else (root,polynomial**(-power))
        if v.is_Pow and v.exp.is_Integer:
            n,d=convert(v.base);power=int(v.exp)
            if power<0:n,d=d,n;power=-power
            return reduce_radicals(n**power),reduce_radicals(d**power)
        raise s.PolynomialError('non-rational-polynomial subexpression')
    def reduce_fraction(n,d):
        n,d=reduce_radicals(n),reduce_radicals(d)
        if d==algebra.one:return n,d
        try:return quotient(n,d),algebra.one
        except s.ExactQuotientFailed:return n,d
    try:
        n,d=convert(expr.function)
        polynomial=quotient(n,d)
        if any(any(p[parameter_end:]) for p in polynomial):return None
    except (s.PolynomialError,s.ExactQuotientFailed):return None
    finally:convert.cache_clear();root_polynomial.cache_clear()
    evaluation_tick('sparse integral monomials '+str(len(polynomial)))
    coefficients=defaultdict(lambda:s.S.Zero)
    for powers,coefficient in polynomial.items():
        factor=s.QQ.to_sympy(coefficient)
        for item in moments:
            if item[0]=='angle':
                _,index,full=item;factor*=angular_moment(powers[index],powers[index+1],full)
            else:
                _,index,lo,hi=item;power=powers[index]+1
                factor*=(hi**power-lo**power)/power
            if factor==0:break
        if factor:coefficients[powers[integrated_count:parameter_end]]+=factor
    return s.Add(*(coefficient*s.prod(generator**power for generator,power in
                    zip(generators[integrated_count:parameter_end],powers))
                   for powers,coefficient in coefficients.items()))


@lru_cache(maxsize=16384)
def integrate_exact(expr):
    if not isinstance(expr,s.Integral): return expr
    sparse=sparse_moments(expr)
    if sparse is not None:return sparse
    evaluation_tick("general integral algebra")
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
            return integrate_numeric(expr)
    return s.cancel(f)

# Evaluation follows the expression DAG. In particular, differentiate a test
# profile in formal argument slots, then evaluate the argument points; never
# differentiate after replacing the moving chart by its value.
EVALUATION_PROGRESS=None
EVALUATION_NODES=0
EVALUATION_LAST=0.0

def evaluation_tick(stage):
    global EVALUATION_NODES,EVALUATION_LAST
    EVALUATION_NODES+=1
    now=time.monotonic()
    if EVALUATION_PROGRESS and now-EVALUATION_LAST>=5:
        print(f'WORKING_ON: {EVALUATION_PROGRESS} stage={stage} evaluated_nodes={EVALUATION_NODES}',file=sys.stderr,flush=True)
        EVALUATION_LAST=now

@lru_cache(maxsize=16384)
def profile_jet(name,arity,orders,form,stratum,index,opts):
    if name not in CODES: raise Unsupported('no test function for '+name)
    value=profile(name,arity,form,stratum,index,dict(opts))
    jet=s.diff(value,*[(v,n) for v,n in zip(FORMAL,orders) if n]) if any(orders) else value
    # Extract common radial factors before inserting square-root distances.
    # This is polynomial factorization, not a field equation or profile change.
    return s.factor_terms(jet) if name.startswith(('u_','rho_')) else jet

@lru_cache(maxsize=4096)
def integral_linearity(expr):
    """Evaluate a sum of integrals by linearity over identical variable ranges.

    Original parsed operands and the as-computed exact residual stay printed.
    This evaluation identity applies no supplied equation, pairs no objects,
    and does not identify surfaces: each integrand retains its chart, normal,
    measure and field arguments. Products of integrals are not combined.
    """
    if not expr.is_Add:return expr
    groups=defaultdict(list);remaining=[]
    for term in expr.args:
        factors=s.Mul.make_args(term)
        integrals=[q for q in factors if isinstance(q,s.Integral)]
        if len(integrals)!=1:remaining.append(term);continue
        integral=integrals[0]
        coefficient=s.Mul(*(q for q in factors if q is not integral))
        if coefficient.free_symbols & set(integral.variables):remaining.append(term);continue
        groups[integral.limits].append((term,coefficient*integral.function))
    if not any(len(terms)>1 for terms in groups.values()):return expr
    for limits,terms in groups.items():
        if len(terms)==1:remaining.append(terms[0][0])
        else:remaining.append(s.Integral(s.Add(*(integrand for _,integrand in terms)),*limits))
    return s.Add(*remaining)


@lru_cache(maxsize=131072)
def evaluated_node(expr,form,stratum,index,opts,environment):
    evaluation_tick('expression DAG')
    env=dict(environment)
    def ev(q,bindings=environment):
        return evaluated_node(q,form,stratum,index,opts,bindings)
    if not expr.args:return env.get(expr,expr)
    if expr.is_Add:
        linear=integral_linearity(expr)
        if linear!=expr:
            evaluation_tick('integral linearity before sampling')
            return ev(linear)
    sig=partial_signature(expr) if isinstance(expr,(s.Derivative,s.Subs)) else None
    if isinstance(expr,AppliedUndef):sig=(expr.func,(0,)*len(expr.args),expr.args)
    if sig is not None:
        head,orders,points=sig
        jet=profile_jet(head.__name__,len(points),orders,form,stratum,index,opts)
        values=tuple(ev(q) for q in points)
        return jet.xreplace(dict(zip(FORMAL,values)))
    if isinstance(expr,s.Subs):
        replacements={v:ev(q) for v,q in zip(expr.variables,expr.point)}
        return ev(expr.expr,tuple(sorted({**env,**replacements}.items(),key=lambda q:s.default_sort_key(q[0]))))
    if isinstance(expr,s.Derivative):
        # General derivatives (including chain/product rules) are performed
        # before the environment is inserted. This is not a field equation.
        replaced=substituted(expr,form,stratum,index,opts)
        return ev(replaced)
    if isinstance(expr,s.Integral):
        variables=set(expr.variables)
        local=tuple((k,v) for k,v in environment if k not in variables)
        integrand=ev(expr.function,local)
        limits=[]
        for j,lim in enumerate(expr.limits):
            outer={q[0] for q in expr.limits[j+1:]}
            bound_env=tuple((k,v) for k,v in environment if k not in outer)
            limits.append((lim[0],*(ev(q,bound_env) for q in lim[1:])))
        evaluation_tick('integral moments')
        return integrate_exact(s.Integral(integrand,*limits))
    args=tuple(ev(q) for q in expr.args)
    return expr.func(*args)

@lru_cache(maxsize=32768)
def evaluate_cached(expr,form,stratum,index,opts=()):
    env=dict(zip(s.symbols('x_1 x_2 x_3 t'),POINTS[index]))
    env.update({s.Symbol('s'):s.Rational(3,2),s.Symbol('theta'):s.Rational(2+index,5),s.Symbol('phi'):s.Rational(1+index,3)})
    env.update(dict(zip(s.symbols('Density Length Time label_1 label_2 label_r'),map(s.Integer,(2,3,5,7,11,13)))))
    environment=tuple(sorted(env.items(),key=lambda q:s.default_sort_key(q[0])))
    value=evaluated_node(expr,form,stratum,index,opts,environment)
    if value.has(AppliedUndef,s.Derivative,s.Subs,s.Integral): raise Unsupported('unevaluated live object '+str(value))
    if value.free_symbols: raise Unsupported('unevaluated symbols '+str(value.free_symbols))
    # All free-point arithmetic has already been evaluated. Simplifying this
    # scalar does not build the original multi-megabyte symbolic polynomial.
    return s.simplify(value)

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

def label(n):
    if n.head in ('Text','Symbol','Number'): return n.value
    if call(n,'Missing'):
        # Missing's reason code is metadata, its second argument is the label.
        if len(n.args)>2: return label(n.args[2])
        return 'Missing['+label(n.args[1])+']'
    if n.head=='Sequence': return [label(q) for q in n.args]
    raise Unsupported('not a native label')

def dim_value(n,b,engine,single=False):
    basis=tuple(s.Symbol(b.name(name,engine)) for name in (('Density','Length','Time') if engine=='py' else ('Dens','Length','Time')))
    def unit(q):
        if q.head!='Sequence' or len(q.args)!=3: raise Unsupported('dimension triple structure')
        return s.prod(v**math_node(e,b,engine) for v,e in zip(basis,q.args))
    if single: return unit(n)
    if n.head!='Sequence': raise Unsupported('dimension alternatives structure')
    return s.Tuple(*(unit(q) for q in n.args))

def emit(out,obj):
    out.write(json.dumps(obj,sort_keys=True,default=str)+'\n'); out.flush()



def compare_value(a,b,row,names,conditions,options,output=None,progress=False):
    global EVALUATION_PROGRESS
    if options.get('disable_comparison'): return {'not_compared':'comparison disabled'}
    if row.kind in ('label','basis'):
        if row.kind=='label' and (a.node.head!='Text' or b.node.head!='Text'):
            raise Unsupported('structural mismatch: source-defined prose label is not text')
        aa,bb=label(a.node),label(b.node)
        if row.kind=='basis':
            aa=[names.name(x,'py') for x in aa]; bb=[names.name(x,'wl') for x in bb]
        return {'py_value':aa,'wl_value':bb,'text_difference':{'py_only':aa,'wl_only':bb} if aa!=bb else []}
    if row.kind=='value' and (semantic_children(a) is not None or semantic_children(b) is not None):
        raise Unsupported('structural mismatch: scalar binding has a container without source-cited component correspondence')
    if row.kind in ('dimension','dimension_field'):
        aa,bb=[dim_value(ref.node,names,e,row.kind=='dimension_field') for ref,e in ((a,'py'),(b,'wl'))]
        if row.kind=='dimension':
            if len(aa)!=len(bb): raise Unsupported('structural mismatch: different numbers of dimension alternatives; no ordinal pairing')
            # Dimensions are alternatives of additive terms, not term ordinals.
            # Retain multiplicity while comparing by dimension value.
            ca,cb=Counter(aa),Counter(bb)
            env=dict(zip(s.symbols('Density Length Time'),map(s.Integer,(2,3,5))))
            return {'py_value':str(aa),'wl_value':str(bb),
                    'dimension_difference':{'py_only':{str(k):v for k,v in (ca-cb).items()},'wl_only':{str(k):v for k,v in (cb-ca).items()}},
                    'dimension_evaluation':{'basis_values':{str(k):str(v) for k,v in env.items()},'py':sorted(str(v.subs(env)) for v in aa),'wl':sorted(str(v.subs(env)) for v in bb)}}
    else:
        aa=math_node(a.node,names,'py'); bb=math_node(b.node,names,'wl')
    original=(aa,bb); aa=names.resolve(aa,'py'); bb=names.resolve(bb,'wl')
    pred=lambda v:v in (s.true,s.false) or isinstance(v,Relational) or v.func in (s.And,s.Or,s.Not,s.Xor,s.Contains)
    if row.kind=='condition' and (not pred(aa) or not pred(bb)): raise Unsupported('structural mismatch: a condition requires a relation or Boolean')
    if pred(aa)!=pred(bb): raise Unsupported('structural mismatch: predicate against non-predicate')
    residual=exact_difference(aa,bb)
    result={'py_value':str(original[0]),'wl_value':str(original[1]),'py_alias_resolved':str(aa),'wl_alias_resolved':str(bb),
            'exact_residual':str(residual),'exact_method':'symbolic subtraction; integral linearity; predicates use symmetric difference; no supplied equation',
            'domain_context':[str(v) for v in conditions]}
    if output is not None:emit(output,{'kind':'exact_comparison','binding':row.key,**result})
    samples=[]
    forms=options.get('forms',FORMS); strata=options.get('strata',STRATA); points=options.get('points',range(len(POINTS)))
    opts=tuple(sorted((k,v) for k,v in options.items() if k in ('duplicate_functions','symmetric_arguments','uncoupled_arguments','vanishing_form')))
    for form,stratum,index in product(forms,strata,points):
        if progress:
            print(f'WORKING_ON: {row.key} form={form} stratum={stratum} point={index}',file=sys.stderr,flush=True)
        EVALUATION_PROGRESS=f'{row.key} form={form} stratum={stratum} point={index}' if progress else None
        sample_started=time.monotonic()
        rec={'form':form,'stratum':stratum,'point':index}
        try:
            # Only the supplied geometry relation guards a value, after the
            # exact residual has been emitted. Conditions are never read as true.
            domains=[evaluate_cached(v,form,stratum,index,opts) for v in conditions]
            rec['conditions']=list(map(str,domains))
            if any(v!=s.true for v in domains):
                rec.update(coverage='unresolved',reason='test geometry outside or unresolved printed domain')
            else:
                vals=[evaluate_cached(v,form,stratum,index,opts) for v in (aa,bb)]
                rec.update(py=str(vals[0]),wl=str(vals[1]),residual=str(exact_difference(*vals)),coverage='evaluated')
        except Unsupported as exc:rec.update(coverage='unresolved',reason=str(exc))
        if output is not None:emit(output,{'kind':'sample','binding':row.key,'elapsed_seconds':time.monotonic()-sample_started,**rec})
        samples.append(rec)
    EVALUATION_PROGRESS=None
    result['samples']=samples
    return result


def unpaired(ref,used,reason):
    all_paths=leaf_paths(ref.node,ref.path)
    if all_paths<=used:return []
    if not all_paths&used:return [{'path':ref.path,'parse':data(ref.node),'leaves':len(all_paths),'reason':reason}]
    return [v for i in range(len(ref.node.args)) for v in unpaired(ref.child(i),used,reason)]


def unjoined_reason(tag,engine):
    if tag.startswith('dimension_argument/'):
        return 'SymPy-only argument-vs-signature/bound dimension check (py 197,209,214,306-310); Wolfram prints field signatures (wl 118) and object dimensions (wl 120), different checks, no printed argument-check counterpart.'
    if tag.startswith('partition/unlabelled_coefficient'):
        return 'SymPy-only empty-support coefficient from polynomial partition (py 108); Wolfram prints nonempty support groups and their sum (wl 61-62,143,167), not an empty-support coefficient.'
    if tag=='chart_domain':return 'SymPy-only printed angular chart-domain prose (py 227); Wolfram uses assumptions internally (wl 81) without a printed chart-domain line.'
    if tag=='c_comp_not_present_in_supplied_equations/comparison':return 'SymPy prints a local comparison note (py 26-28,76); Wolfram prints the free parameter with emit, without a local comparison note (wl 45).'
    return 'No declared source-cited counterpart; retained and listed for accounting; this is an unresolved binding-coverage question.'


def compare_binding(row,streams,b,conditions,options,out,progress):
    used={e:defaultdict(set) for e in streams};bound={e:defaultdict(list) for e in streams}
    total=Counter();unresolved=Counter();notformed=[]
    if progress:print('WORKING_ON: '+row.key+' '+row.py.tag+' '+str(row.py.path),file=sys.stderr,flush=True)
    es={e:streams[e].index.get(getattr(row,e).tag,[]) for e in streams}
    emit(out,{'kind':'operands','binding':row.key,**{e:[{'tag':v.original_tag,'record_identity':v.tag,'line':v.line,'raw':v.raw} for v in vs] for e,vs in es.items()}})
    refs={}; parsed={e:0 for e in streams};actual={e:0 for e in streams}
    try:
        if any(len(v)!=1 for v in es.values()):raise Unsupported('binding end missing or repeated: '+str({e:len(v) for e,v in es.items()}))
        for e,vs in es.items():
            end=getattr(row,e);bound[e][end.tag].append(row.key)
            if ROOT_KINDS.get((e,end.tag))=='equation' and not (call(vs[0].node,'Eq') or call(vs[0].node,'Equal') or vs[0].node.head=='Equality'):
                raise Unsupported('structural mismatch: source-defined equation relation is absent')
            shape=ROOT_SHAPES.get((e,end.tag))
            if shape is not None and not shaped(Ref(vs[0].node),shape):
                raise Unsupported('structural mismatch with source-defined full layout: '+e+' '+end.tag)
            refs[e]=get_ref(Ref(vs[0].node),end.path);parsed[e]=census(refs[e].node)
        computed=compare_value(refs['py'],refs['wl'],row,b,(s.true,s.true) if row.kind=='condition' else conditions,options,out,progress=progress)
        if 'not_compared' not in computed:
            for e,ref in refs.items():
                leaves=leaf_paths(ref.node,ref.path)
                if not options.get('omit_leaf_credit'):used[e][getattr(row,e).tag].update(leaves);actual[e]=len(leaves)
        for rec in computed.get('samples',[]):
            if rec['coverage']=='unresolved':unresolved[rec['reason']]+=1
        total['value_comparisons']+=1
    except Unsupported as exc:
        computed={'not_formed':str(exc)};notformed.append({'binding':row.key,'reason':str(exc)})
        for e,vs in es.items():
            if e not in refs:parsed[e]=sum(census(v.node) for v in vs)
    emit(out,{'kind':'comparison','binding':row.key,**computed})
    emit(out,{'kind':'binding_accounting','binding':row.key,'parsed_leaves':parsed,'actually_compared_leaves':actual})
    return used,bound,total,unresolved,notformed,parsed,actual

# Workers inherit only read-only input/translation state through fork. Each
# binding has its own durable spool; the parent publishes in table order.
# Completion order never selects, joins, truncates or changes an object.
WORKER_CONTEXT=None

def worker_binding(index):
    streams,b,conditions,options,progress,bindings,spool=WORKER_CONTEXT
    row=bindings[index]
    global EVALUATION_PROGRESS
    EVALUATION_PROGRESS=row.key+' extracting and comparing'
    # Drop expression DAGs from prior objects; coverage is unaffected by cache
    # lifetime. The fork-inherited native streams remain available read-only.
    for cache in (evaluated_node,evaluate_cached,substituted,integrate_exact,integral_linearity):cache.cache_clear()
    s.core.cache.clear_cache();gc.collect()
    before=set(NUMERICAL_INTEGRALS)
    path=spool/f'{index:06d}.jsonl'
    with path.open('x') as out:
        result=compare_binding(row,streams,b,conditions,options,out,progress)
    EVALUATION_PROGRESS=None
    if progress: print('WORKING_ON: completed '+row.key,file=sys.stderr,flush=True)
    return index,result,{k:v for k,v in NUMERICAL_INTEGRALS.items() if k not in before}

def worker_heartbeat():
    # Observation only: no deadline, cancellation, or change of coverage.
    while True:
        threading.Event().wait(15)
        if EVALUATION_PROGRESS:
            print(f'WORKING_ON: worker={os.getpid()} {EVALUATION_PROGRESS} nodes={EVALUATION_NODES}',file=sys.stderr,flush=True)

def worker_init():
    if WORKER_CONTEXT[4]:threading.Thread(target=worker_heartbeat,daemon=True).start()

@contextmanager
def binding_results(streams,b,conditions,options,progress,bindings,out,workers,spool,order):
    global WORKER_CONTEXT
    if workers==1:
        def serial():
            for index,row in enumerate(bindings):
                yield index,compare_binding(row,streams,b,conditions,options,out,progress),{}
        yield serial()
        return
    if workers<1 or workers>len(os.sched_getaffinity(0)):
        raise InputError('worker count exceeds inherited guarded CPU allocation')
    spool=Path(spool)
    spool.mkdir(parents=True,exist_ok=False)
    WORKER_CONTEXT=(streams,b,conditions,options,progress,bindings,spool)
    indices=list(range(len(bindings)))
    if order=='reverse':indices.reverse()
    elif order!='forward':raise InputError('unknown work submission order')
    with multiprocessing.get_context('fork').Pool(workers,initializer=worker_init) as pool:
        pending={}
        def ordered():
            next_index=0
            for result in pool.imap_unordered(worker_binding,indices,chunksize=1):
                pending[result[0]]=result
                while next_index in pending:
                    value=pending.pop(next_index)
                    with (spool/f'{next_index:06d}.jsonl').open() as source:
                        for line in source:out.write(line)
                    out.flush()
                    yield value
                    next_index+=1
        yield ordered()
    WORKER_CONTEXT=None


def run(py_path,wl_path,out,*,names=NAME_TABLE,bindings=BINDING_TABLE,options=None,progress=True,metadata=True,workers=1,spool=None,submission_order="forward"):
    started=time.monotonic(); options=options or {}; validate_tables(names,bindings)
    streams={e:Stream(p,e,progress=progress) for e,p in (('py',py_path),('wl',wl_path))}
    b=Names(names);b.load(streams)
    emit(out,{'kind':'configuration','sources':'S15 engines 3192123b; decision list revision 3 91248b67; comparator contract e305f6d4',
              'input_sha256':{e:st.sha256 for e,st in streams.items()},'binding_rules':RULES,'rule_declarations':RULE_NAMES,
              'surface_table':SURFACE_TABLE,'binding_table':[asdict(r) for r in bindings],'name_table':[asdict(r) for r in names],'bound_name_table':[asdict(r) for r in SLOT_TABLE],
              'table_sizes':{'rules':len(RULES),'bindings':len(bindings),'names':len(names),'bound_names':len(SLOT_TABLE)},
              'root_shapes':[{'engine':e,'tag':tag,'shape':shape,'relation_kind':ROOT_KINDS.get((e,tag))} for (e,tag),shape in ROOT_SHAPES.items()],
              'dummy_rule':'Subs and integral dummies are lexical binders, alpha-renamed without capture; py 137,269 and SymPy differentiation at 63,79,116,284; wl derivative slots and binders at 32,36,63,85-90.',
              'alias_rule':'No printed term-alias symbols/definitions in accepted S15 sources; explicit synthetic aliases resolve only through their printed defining equations.',
              'forms':FORMS,'strata':STRATA,'point_variables':['x_1','x_2','x_3','t'],'points':[[str(v) for v in p] for p in POINTS],
              'constant_probe_values':{'Density':'2','Length':'3','Time':'5','label_1':'7','label_2':'11','label_r':'13'},
              'residual_target':'No residual target was given.'})
    if metadata:
        for form,stratum,index in product(FORMS,STRATA,range(len(POINTS))):
            geometry={r.py:str(evaluate_cached(s.Function(r.py)(s.Symbol('t')),form,stratum,index)) for r in names if r.py.startswith('X_') or r.py in ('r_in','R')}
            emit(out,{'kind':'test_functions','form':form,'stratum':stratum,'point':index,
                      'formal_arguments':list(map(str,FORMAL)),
                      'functions':{r.py:str(profile(r.py,r.arity,form,stratum,index)) for r in names if r.arity>0 and not r.alias},
                      'geometry_at_point':geometry,
                      'pressure_condition':'Polynomial pressure is defined and differentiable at all tested densities; supplied momentum balance is not imposed.',
                      'free_chart_point':{'s':'3/2','theta':str(s.Rational(2+index,5)),'phi':str(s.Rational(index+1,3))}})
    for e,st in streams.items():
        for entry in st.progress: emit(out,{'kind':'progress_marker','engine':e,'line':entry.line,'tag':entry.original_tag,'raw':entry.raw})
    conditions=[]
    for e,st in streams.items():
        entries=st.index.get('DOMAIN_FOR_EVERY_OBJECT',[])
        if not entries:
            conditions.append(s.true)
            continue
        if len(entries)!=1:
            conditions.append(s.Symbol('unresolved_domain_'+e))
            emit(out,{'kind':'unresolved_domain','engine':e,'tag':'DOMAIN_FOR_EVERY_OBJECT','reason':'repeated domain tag; no unique relation','raw':[v.raw for v in entries]})
            continue
        try:
            if not shaped(Ref(entries[0].node),ROOT_SHAPES[e,'DOMAIN_FOR_EVERY_OBJECT']):
                raise Unsupported('source-defined domain payload structure')
            condition=math_node(get_ref(Ref(entries[0].node),(0,)).node,b,e)
            if not (condition in (s.true,s.false) or isinstance(condition,Relational) or condition.func in (s.And,s.Or,s.Not,s.Xor,s.Contains)):
                raise Unsupported('printed domain is not a relation or Boolean')
            conditions.append(condition)
        except Unsupported as exc:
            conditions.append(s.Symbol('unresolved_domain_'+e))
            emit(out,{'kind':'unresolved_domain','engine':e,'tag':'DOMAIN_FOR_EVERY_OBJECT','reason':str(exc),'raw':entries[0].raw})
    conditions=tuple(conditions)
    used={e:defaultdict(set) for e in streams};bound={e:defaultdict(list) for e in streams};total=Counter();below=[]
    unresolved=Counter();unjoined=Counter();notformed=[]
    with binding_results(streams,b,conditions,options,progress,bindings,out,workers,spool,submission_order) as results:
      for index,result,numerical in results:
        row=bindings[index]
        row_used,row_bound,row_total,row_unresolved,row_notformed,parsed,actual=result
        for e in streams:
            for tag,paths in row_used[e].items():used[e][tag].update(paths)
            for tag,keys in row_bound[e].items():bound[e][tag].extend(keys)
        total.update(row_total);unresolved.update(row_unresolved);notformed.extend(row_notformed)
        NUMERICAL_INTEGRALS.update(numerical)
        total['bindings']+=1
        for e in streams:total['bound_parsed_leaves_'+e]+=parsed[e];total['compared_leaves_'+e]+=actual[e]
        if any(actual[e]<parsed[e] for e in streams):below.append(row.key)
    for e,st in streams.items():
        total['lines_'+e]=st.line_count
        for tag,entries in st.index.items():
            for occurrence,v in enumerate(entries):
                keys=bound[e][tag] if len(entries)==1 else [];count=census(v.node);actual=len(used[e][tag]) if len(entries)==1 else 0
                rec={'kind':'tag_accounting','engine':e,'tag':v.original_tag,'record_identity':tag,'line':v.line,'occurrence':occurrence,
                     'parsed_leaves':count,'actually_compared_leaves':actual,'accounting':'bound' if keys else 'unjoined'}
                total['tags_'+e]+=1;total['parsed_leaves_'+e]+=count;total['tag_compared_leaves_'+e]+=actual
                if keys:
                    rec.update(bindings=keys,unpaired_components=unpaired(Ref(v.node),used[e][tag],
                        'Constructor or key leaf is not a compared value, or no active source-cited component correspondence covers this content; see binding accounting.'))
                else:
                    reason='Repeated tag or record identity; no comparison formed.' if len(entries)>1 else unjoined_reason(tag,e)
                    rec.update(reason=reason,raw=v.raw,parse=data(v.node));unjoined[e+' '+reason]+=1;total['unjoined_tags_'+e]+=1
                emit(out,rec)
    emit(out,{'kind':'resources','runtime_seconds':time.monotonic()-started,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    summary={'kind':'summary','counts':dict(total),'bindings_below_parsed_count':below,'not_formed':notformed,
             'unjoined_reasons':dict(unjoined),'unresolved_coverage':dict(unresolved),
             'numerical_integrals':[NUMERICAL_INTEGRALS[k] for k in sorted(NUMERICAL_INTEGRALS)],
             'integral_method':'Exact integral linearity and sparse polynomial/trigonometric moments with outer parameters retained, otherwise 40/60-digit quadrature with difference error estimate; inability remains unresolved',
             'exit_semantics':'zero for findings; nonzero only for operational failure'}
    emit(out,summary);return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--py',required=True,type=Path);parser.add_argument('--wl',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--workers',type=int,default=1)
    parser.add_argument('--spool',type=Path)
    args=parser.parse_args()
    try:
        with args.output.open('x') as out:run(args.py,args.wl,out,workers=args.workers,spool=args.spool or args.output.with_suffix('.bindings'))
    except InputError as exc:
        print('Operational input error: '+str(exc),file=sys.stderr);return 2
    return 0

if __name__=='__main__':sys.exit(main())
