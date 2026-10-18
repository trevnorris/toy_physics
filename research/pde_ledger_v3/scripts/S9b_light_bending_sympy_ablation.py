#!/usr/bin/env python3
"""S9b amendment 12: literal marks and finite searches of emitted site parts.

Run through scripts/s11c_guarded_run.py --pool s9b --memory-gib 22 (user authorization for this round).
No deadlines, retries, full K11 worker, or mathematical equality classifications.
Raw engine payloads are preserved. Only the per-knife witness search evaluates
selected site parts; failure to decide a condition produces no witness.
"""
import argparse
import ast
from collections import Counter
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import os
import subprocess
import sys
import time
sys.dont_write_bytecode = True
import sympy as sp
from sympy.core.function import AppliedUndef
from sympy.core.symbol import Str
from sympy.functions.elementary.piecewise import ExprCondPair
from sympy.printing.repr import ReprPrinter

HERE = Path(__file__).resolve()
ENGINE = HERE.with_name('S9b_light_bending_sympy_audit.py')
ROOT = HERE.parent.parent
REPOSITORY = ROOT.parent.parent
INVENTORY = ROOT/'directives/_measurements/S9b_repair_sympy_tag_inventory.json'
KNIVES = {
    'K1': ('kinetic = (omega-advected_velocity.dot(momenta))**2', 'kinetic = omega**2'),
    'K2': ('inverse_metric = sp.diag(1/A,1,1)', 'inverse_metric = sp.eye(3)'),
    'K3': ('local_speed_squared = chi', 'local_speed_squared = c0**2'),
    'K4a': ('RADIAL_FREEZE = None', 'RADIAL_FREEZE = 0'),
    'K4b': ('RADIAL_FREEZE = None', 'RADIAL_FREEZE = 1'),
    'K4c': ('RADIAL_FREEZE = None', 'RADIAL_FREEZE = 2'),
    'K5': ('return_orientation = -1', 'return_orientation = 1'),
    'K6': ('log_source = roundtrip_time', 'log_source = forward_time'),
    'K7': ('fixed_impact = None', "fixed_impact = sp.Symbol('s9b_fixed_b',positive=True)"),
    'K8': ('mass_density = rho', "mass_density = sp.Symbol('s9b_constant_density',positive=True)"),
    'K9a': ('local_sound_speed = sp.sqrt(local_square)', 'local_sound_speed = sp.sqrt(anchor_square)'),
    'K9b': ('power_density = density', 'power_density = rho0'),
    'K10': ('observable_gate = real_domain', 'observable_gate = sp.Not(real_domain)'),
    'K11': ('azimuthal_component = sp.S.Zero', "azimuthal_component = sp.Symbol('s9b_azimuthal',real=True)"),
}
COMPACT_PARAMETERS = {'d':sp.Rational(1,10), 'v':sp.Rational(1,8),
                      'w':sp.Rational(1,6), 'p':sp.Integer(2),
                      'q':sp.Integer(1), 'h':sp.Integer(1),
                      'L':sp.Integer(1), 'c0':sp.Integer(1)}
COMPACT_AZIMUTHAL = sp.Rational(1,12)
COMPACT_POINTS = ((3,4,0),(4,0,3),(0,6,8))
COMPACT_TARGETS = {'PY_LOCAL_S9B_NONRECIPROCITY_ONEFORM':'oneform',
                   'PY_LOCAL_S9B_NONRECIPROCITY_EXTERIOR_DERIVATIVE':'curl',
                   'PY_LOCAL_S9B_NONRECIPROCITY_PATH_PREDICATE':'predicate'}
RELATIONALS = {name:(lambda a,b,f=f:f(a,b,evaluate=False)) for name,f in
    (('Equality',sp.Eq),('Unequality',sp.Ne),('StrictGreaterThan',sp.Gt),
     ('StrictLessThan',sp.Lt),('GreaterThan',sp.Ge),('LessThan',sp.Le))}


def progress(name):
    print('S9B_HARNESS_WORK '+str(name),file=sys.stderr,flush=True)


@contextmanager
def work_part(name):
    progress('begin '+name)
    start=time.monotonic()
    try:yield
    finally:progress('end '+name+' seconds='+format(time.monotonic()-start,'.6f'))


class PayloadPrinter(ReprPrinter):
    def _print_FunctionClass(self,expr):
        if issubclass(expr,AppliedUndef):
            return 'Function(%r, **%r)' % (expr.__name__,dict(sorted(expr._kwargs.items())))
        return super()._print_FunctionClass(expr)

    def _print_Intersection(self,expr):
        return 'Intersection('+', '.join(self._print(a) for a in expr.args)+', evaluate=False)'

    def _print_Union(self,expr):
        return 'Union('+', '.join(self._print(a) for a in expr.args)+', evaluate=False)'

    def _print_SymmetricDifference(self,expr):
        return 'SymmetricDifference('+', '.join(self._print(a) for a in expr.args)+', evaluate=False)'


def render(value):
    return PayloadPrinter().doprint(value)


def cas(value):
    if isinstance(value,str):return Str(value)
    if isinstance(value,dict):return sp.Tuple(*(sp.Tuple(cas(k),cas(v)) for k,v in value.items()))
    if isinstance(value,(list,tuple)):return sp.Tuple(*(cas(v) for v in value))
    return sp.sympify(value)


def emit(name,value):
    print('PY_LOCAL_S9B_ABLATION_'+name+': '+render(cas(value)),flush=True)


STRUCTURAL_METHOD = ('Literal Python expression AST, preserving node types, field order, '
    'argument order, keyword order, and the original token of every Constant; '
    'no eval, simplification, alpha-renaming, arithmetic or canonicalization. '
    'Whitespace outside constants is syntax, not an expression node. '
    'Different float tokens or precisions cannot be identical.')


def literal_key(payload):
    """Read syntax only. In particular Python never rounds a Float argument."""
    tree=ast.parse(payload,mode='eval')
    # Positions affect neither expression structure nor literal spelling.
    # Keep Constant source segments: ast.Constant.value alone rounds floats.
    lines=[line.encode() for line in payload.splitlines(keepends=True)]
    def source_segment(node):
        if node.lineno!=node.end_lineno:
            return ast.get_source_segment(payload,node)
        return lines[node.lineno-1][node.col_offset:node.end_col_offset].decode()
    def visit(node):
        if isinstance(node,ast.Constant):
            return ('Constant',type(node.value).__name__,source_segment(node))
        if isinstance(node,ast.AST):
            return (type(node).__name__,tuple((k,visit(v)) for k,v in ast.iter_fields(node)))
        if isinstance(node,list):return tuple(visit(x) for x in node)
        return node
    return visit(tree)


def structural_mark(left,right):
    return 'identical' if literal_key(left)==literal_key(right) else 'not compared'


def read_stream(path):
    result={}
    with path.open() as stream:
        for number,line in enumerate(stream,1):
            tag,sep,payload=line.rstrip('\n').partition(': ')
            if not sep or not tag.startswith(('PY_S9B_','PY_LOCAL_S9B_')):
                raise ValueError(('untagged line',str(path),number))
            if tag in result:raise ValueError(('duplicate tag',tag))
            result[tag]=payload
    if not result:raise ValueError(('empty stream',str(path)))
    return result


def decode(payload):
    # ONLY witness candidates and K11 baseline reduction use this. Structural
    # marks always precede and bypass reconstruction by the CAS.
    return eval(payload,{'__builtins__':{},**vars(sp),'Str':Str,
                         'ExprCondPair':ExprCondPair,**RELATIONALS})


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def source_location(source,text):
    if source.count(text)!=1:raise ValueError(('source-location ambiguity',text))
    return {'file':str(ENGINE.relative_to(REPOSITORY)),
            'line':source[:source.index(text)].count('\n')+1,'text':text}


# The witness plans select downstream emitted construction objects, not the
# input digest or a pre-mutation object. They introduce no engine mutation.
PLANS = {
    'K1':('PY_LOCAL_S9B_NONRECIPROCITY_ONEFORM','matrix00',"components = [evaluate_jet(component,mapping).tidy() for component in source['centre']]"),
    'K2':('PY_LOCAL_S9B_NONRECIPROCITY_ONEFORM','matrix00',"components = [evaluate_jet(component,mapping).tidy() for component in source['centre']]"),
    'K3':('PY_LOCAL_S9B_NONRECIPROCITY_ONEFORM','matrix00',"components = [evaluate_jet(component,mapping).tidy() for component in source['centre']]"),
    'K4a':('PY_LOCAL_S9B_OPTICAL_Q',(1,2,0),'optical_q=((arrj**sp.Rational(1,2)*c0)/optical_r.diff(r)).tidy()'),
    'K4b':('PY_LOCAL_S9B_OPTICAL_Q',(0,2,0),'optical_q=((arrj**sp.Rational(1,2)*c0)/optical_r.diff(r)).tidy()'),
    'K4c':('PY_LOCAL_S9B_NONRECIPROCITY_ONEFORM','matrix00','slope = radial_diff(xi,r)'),
    'K5':('PY_LOCAL_S9B_NONRECIPROCITY_ONEFORM','matrix00','path=nonreciprocity(r,c0,delta,velocity,xi,return_orientation)'),
    'K6':('PY_LOCAL_S9B_RADAR_SLOPE_SOURCE_G010','whole','source=log_source[g].subs(endpoint_map)'),
    'K7':('PY_S9B_B_CONDITION_DEFLECTION','reduced_condition','local_mass=substitute(mass,replacement)'),
    'K8':('PY_LOCAL_S9B_MASS_BALANCE','rhs','exchange=sp.simplify(-sp.diff(area*mass_density*velocity,r)/area)'),
    'K9a':('PY_LOCAL_S9B_C_FIXED_RATIO_RESPONSE','response','responses=bulk_responses(rho0,K,n,m,s,f_symbol)'),
    'K9b':('PY_LOCAL_S9B_C_POWER_RESPONSE','response','responses=bulk_responses(rho0,K,n,m,s,f_symbol)'),
    'K10':('PY_S9B_A_DEFLECTION_G000','gate','branch_domain=sp.And(observable_gate,travel_domain)'),
    'K11':('PY_LOCAL_S9B_NONRECIPROCITY_ONEFORM','compact_sample0_matrix00',
            "components = [evaluate_jet(component,mapping).tidy() for component in source['centre']]"),
}


def at_path(value,path):
    for index in path:value=value.args[index]
    return value


def select_part(whole,selector):
    if isinstance(selector,tuple):
        grade=sp.Tuple(*selector)
        for index,row in enumerate(whole):
            if row[0]==grade:return row[1],(index,1)
        # Sparse Jet dictionaries omit zero coefficients. The object is the
        # emitted mapping; an absent entry is not an emitted value witness.
        return None,None
    if selector=='matrix00':return whole[0,0],(2,0)
    if selector=='whole':return whole,()
    if selector=='rhs':return whole.rhs,(1,)
    if selector=='response':return whole[3][1][0][1],(3,1,0,1)
    if selector=='reduced_condition':return whole[4],(4,)
    if selector=='gate':return whole[0][0][1],(0,0,1)
    raise ValueError(selector)


def attached_conditions(whole,path):
    """Collect domain siblings at every containing level, with case logic.

    Searches below deliberately use only local ungated intermediates, except
    the explicit reduced-condition and gate searches. Those two searches
    retain all quantifiers and cannot promote unresolved ray sets to truth.
    Numeric grade keys are dictionary entries, never ordered cases. Density
    sign labels partition the same object and are not attached predicates.
    """
    found=[]
    for depth in range(len(path)+1):
        parent=at_path(whole,path[:depth])
        if not isinstance(parent,sp.Tuple):continue
        for j,child in enumerate(parent):
            if (isinstance(child,sp.Tuple) and len(child)==2 and child[0]==Str('ATTACHED_DOMAIN')):
                location=path[:depth]+(j,1)
                if location not in [p for p,c in found]:found.append((location,child[1]))
    return found


RADIUS=sp.Symbol('s9b_witness_radius',positive=True)
# Each entry is a specific fixed function. These never enter an engine run.
PROFILES = {'s9b_delta':1/(10*RADIUS**2), 's9b_V':1/(8*RADIUS**2),
            's9b_xi':1/(6*RADIUS), 'rho_br_profile':1+1/RADIUS,
            's9b_f':1/(10*RADIUS**2), 'mu_perp':(1+1/RADIUS)*(1+1/(10*RADIUS**2))**2}
POINT = {'s9b_X':2,'s9b_r':2,'c_0':1,'b':12,'Z_E':5,'Z_R':9,
         's9b_b_min':1,'s9b_fixed_b':12,'s9b_constant_density':1,
         'rho_0':1,'K':1,'m':1,'n':3,'s':1,'GM':1,
         'x_1':3,'x_2':4,'x_3':0}


def functions_in(value):
    return {a.func for a in value.atoms(AppliedUndef)}


def specialize_profiles(value,profiles):
    value=value.replace(lambda x:isinstance(x,AppliedUndef) and x.func.__name__ in profiles,
        lambda x:profiles[x.func.__name__].subs(RADIUS,x.args[0]))
    # Differentiate before any free point coordinate is substituted. Subs
    # variables, Integral variables and ConditionSet variables stay bound.
    return value.doit(deep=True)


def evaluate(value,profiles,point):
    value=specialize_profiles(value,profiles)
    value=value.subs(point,simultaneous=True)
    return sp.simplify(value.doit(deep=True))


def finite_exact(value):
    if isinstance(value,(sp.Tuple,sp.MatrixBase)):
        return all(finite_exact(x) for x in value)
    return value in (sp.true,sp.false) or (isinstance(value,sp.Expr) and
        not value.free_symbols and not value.has(AppliedUndef,sp.Integral,sp.Limit,sp.Derivative,sp.Subs)
        and value.is_real is True and value.is_finite is True)


def different_exact(a,b):
    if isinstance(a,sp.Tuple) and isinstance(b,sp.Tuple):
        return len(a)!=len(b) or any(different_exact(x,y) for x,y in zip(a,b))
    if a in (sp.true,sp.false) and b in (sp.true,sp.false):return a!=b
    if not all(isinstance(x,sp.Expr) for x in (a,b)):return False
    return sp.simplify(a-b).is_zero is False


def exact_separation(a,b,path=()):
    """A finite-value inequality certificate, not a structural change claim."""
    if isinstance(a,sp.Tuple) and isinstance(b,sp.Tuple):
        if len(a)!=len(b):return ('finite tuple cardinalities',path,len(a),len(b),sp.Ne(len(a),len(b)))
        for index,(x,y) in enumerate(zip(a,b)):
            if different_exact(x,y):return exact_separation(x,y,path+(index,))
        return ('no exact separation',)
    if all(isinstance(x,sp.Expr) for x in (a,b)):
        delta=sp.simplify(a-b)
        return ('exact component inequality',path,a,b,delta,sp.Ne(delta,0))
    return ('exact truth values',path,a,b,sp.sympify(a!=b))


def binders(value):
    rows=[]
    for node in sp.preorder_traversal(value):
        if isinstance(node,sp.Integral):rows.append(('Integral',node.variables,node.limits))
        elif isinstance(node,sp.ConditionSet):rows.append(('ConditionSet',(node.sym,),node.base_set))
        elif isinstance(node,sp.Limit):rows.append(('Limit',(node.args[1],),node.args[2:]))
        elif isinstance(node,sp.Subs):rows.append(('Subs',node.variables,node.point))
    return rows


def input_reading(parts,conditions,profiles,point,fixings=()):
    whole=sp.Tuple(*parts,*(condition for group in conditions for path,condition in group),*fixings)
    free=whole.free_symbols
    roles=[]
    fixed={eq.lhs for eq in fixings if isinstance(eq,sp.Equality)}
    for symbol in sorted(whole.atoms(sp.Symbol),key=sp.default_sort_key):
        occurrences=[]
        if symbol in free:
            if symbol in fixed:occurrences.append(('fixed by construction',tuple(eq for eq in fixings if eq.lhs==symbol)))
            elif symbol in point:occurrences.append(('assigned at point',point[symbol]))
            else:occurrences.append(('unresolved free input',symbol))
        for kind,variables,ranges in binders(whole):
            if symbol in variables:occurrences.append(('bound',kind,ranges))
        roles.append((symbol,occurrences))
    function_rows=[]
    for function in sorted(functions_in(whole),key=str):
        name=function.__name__
        if name not in profiles:raise ValueError(('unassigned witness profile',name))
        expression=profiles[name]
        required=max([0]+[int(order) for d in whole.atoms(sp.Derivative)
            if any(a.func==function for a in d.atoms(AppliedUndef)) for variable,order in d.variable_count])
        # Rational profiles here have no poles on r>0; log(r) in K11 is
        # differentiated on that same open interval. Print the expressions.
        derivatives=tuple((order,sp.diff(expression,RADIUS,order)) for order in range(required+1))
        function_rows.append((Str(name),sp.Lambda(RADIUS,expression),
            ('range of profile argument',sp.Interval.open(0,sp.oo)),
            ('required derivatives',derivatives),
            ('singularities on profile range',sp.Intersection(
                sp.singularities(expression,RADIUS),sp.Interval.open(0,sp.oo)))))
    return {'symbol_roles':roles,'profile_inputs':function_rows,
            'point':tuple(sorted(point.items(),key=lambda x:sp.default_sort_key(x[0]))),
            'bound_occurrences':binders(whole),
            'reading_scope':'Selected local emitted site part and all containing domains. No on-shell dispersion equation, ray-root choice, matching relation, or order-counting entry is attached to these local intermediate coefficients.',
            'worded_attached_conditions':(),
            'object_labels':(), 'ordered_case_exclusions':()}


# These are dependencies of constructions, not tests on emitted expressions.
# A zero coefficient remains a site part. Source entries below are checked
# against the pinned engine. Export rows inherit their originating object.
FORM={'K1','K2','K3','K11'}
RAY=FORM|{'K4b','K4c'}
OPTICAL=FORM|{'K4a','K4b','K4c'}


def construction_dependencies(name):
    """Source-level slice of build(), including its nested emit functions."""
    name=name.removeprefix('PY_LOCAL_S9B_').removeprefix('PY_S9B_')
    if name.startswith('s9b_'):
        return construction_dependencies(name[4:].upper())
    local={
        'VECTOR_DISPERSION':(FORM,'dispersion = kinetic-local_speed_squared*(momenta.T*inverse_metric*momenta)[0]'),
        'VECTOR_ADVECTION':({'K11'},'advected_velocity = sp.Matrix([U,azimuthal_component,0])'),
        'RAY_RADIAL_DOMAINS':(RAY,'rays=ray_domains(geometry,r,c0,b,ZE,ZR,branch_map)'),
        'BRANCH_EXISTENCE':(RAY,'real_domain = everywhere(r,existence.subs(branch_map),traversed)'),
        'PATH_TRAVERSAL':(RAY,'travel_domain = sp.And('),
        'BRANCH_TYPE':(FORM|{'K4c'},"emit('BRANCH_TYPE',cas(branch_type).subs(branch_map))"),
        'OPTICAL_RADIUS':(FORM|{'K4c'},'optical_r=(aphij**sp.Rational(1,2)*c0).tidy()'),
        'OPTICAL_Q':(OPTICAL,'optical_q=compose_radial(optical_q,r,x,inverse).tidy()'),
        'NONRECIPROCITY_ONEFORM':(FORM|{'K4c','K5'},"oneform = (oneform*(1-return_orientation)/2).applyfunc(sp.factor)"),
        'NONRECIPROCITY_EXTERIOR_DERIVATIVE':(OPTICAL|{'K5'},'curl = sp.ImmutableMatrix('),
        'NONRECIPROCITY_PATH_PREDICATE':(OPTICAL|{'K5'},'closed=sp.And(*(sp.Eq(component,0) for component in curl))'),
        'RADAR_ANGULAR_OPERATOR_RATIO':(OPTICAL|{'K5','K6'},'slope_ratio=sp.simplify(slope_density/(qsum*angular_kernel.subs(B,b)))'),
        'ABEL_REFERENCE_INVERSES':(OPTICAL|{'K5','K6'},'inverse_targets[name]=sp.simplify(-r*r*integral/sp.pi)'),
        'MASS_BALANCE':({'K8'},'exchange=sp.simplify(-sp.diff(area*mass_density*velocity,r)/area)'),
        'F_MASS_SOLUTION':({'K8'},'solution=sp.dsolve(eq(sp.diff(area*mass_density*liveV,r),0),liveV)'),
    }
    if name in local:return local[name]
    if name.startswith('RADAR_'):
        return OPTICAL|{'K5','K6'},'source=log_source[g].subs(endpoint_map)'
    forward=name.startswith(('F_FLOW_','F_LIVE_'))
    if forward:
        name=name.split('_',2)[2]
        if name in ('BRANCH_EXISTENCE','PATH_TRAVERSAL','BRANCH_TYPE'):
            deps,entry=construction_dependencies(name)
            return deps|{'K8'},entry
    # Restriction substitutions are downstream computations; do not discard
    # their zero results on the ground that a profile was set to zero.
    if name.startswith('R_BULK_ONLY_'):name=name[len('R_BULK_ONLY_'):]
    for prefix in ('R_FLOW_ONLY_','R_SPEED_ONLY_','R_TILT_ONLY_'):
        if name.startswith(prefix):name=name[len(prefix):]
    deps=({'K8'} if forward else set())
    response=({'K9a'} if name.startswith('C_FIXED_RATIO_') else
              {'K9b'} if name.startswith('C_POWER_') else set())
    if name.startswith('C_') and name.endswith('_RESPONSE'):
        return response,"responses=bulk_responses(rho0,K,n,m,s,f_symbol)"
    if '_CONDITION_' in name or '_IMPLIED_JN_' in name or name.endswith('_N_DEPENDENCE'):
        deps|=OPTICAL|{'K5','K6','K7','K10'}|response
        if name.endswith('_DEFLECTION'):deps-={'K5','K6'}
        if '_IMPLIED_JN_' in name:deps.add('K8')
        return deps,'local_mass=substitute(mass,replacement)'
    if name.startswith('B_'):
        deps|=OPTICAL|{'K5','K6','K10'}
        if name.endswith('_DEFLECTION'):deps-={'K5','K6'}
        return deps,"solved[name]=gamma_cases(references[name],substitute(comparisons[name],replacement),reference_gamma,GM)"
    if name.startswith('A_'):
        deps|=RAY|{'K10'}
        if name=='A_NONRECIPROCAL_PATH_DEPENDENCE':deps|=OPTICAL|{'K5'}
        elif 'RADAR_LOG' in name:deps|=OPTICAL|{'K5','K6'}
        elif any(word in name for word in ('ROUND_TRIP','ONE_WAY_RE','NONRECIPROCAL')):deps|=OPTICAL|{'K5'}
        else:deps|=OPTICAL
        return deps,'physical_branch=substitute_domain(branch_domain,replacement)'
    return set(),None


def engine_line(source,text):
    # A prefix is used only for a single multi-line construction statement.
    return source_location(source,text)


def mutation_entry(knife,name,source):
    entries={
        'K1':'dispersion = kinetic-local_speed_squared*(momenta.T*inverse_metric*momenta)[0]',
        'K2':'dispersion = kinetic-local_speed_squared*(momenta.T*inverse_metric*momenta)[0]',
        'K3':'dispersion = kinetic-local_speed_squared*(momenta.T*inverse_metric*momenta)[0]',
        'K4a':"name=('s9b_delta','s9b_V','s9b_xi')[RADIAL_FREEZE]",
        'K4b':"name=('s9b_delta','s9b_V','s9b_xi')[RADIAL_FREEZE]",
        'K4c':"name=('s9b_delta','s9b_V','s9b_xi')[RADIAL_FREEZE]",
        'K5':('oneform = (oneform*(1-return_orientation)/2).applyfunc(sp.factor)' if 'NONRECIPROCITY' in name else 'reverse_time=even_time+return_orientation*odd_time'),
        'K6':'source=log_source[g].subs(endpoint_map)',
        'K7':'condition_mass={name:sp.solve(eq(comparisons[name].subs(b,fixed_impact),references[name].subs({b:fixed_impact,reference_gamma:1})),GM)[0] for name in comparisons}',
        'K8':('solution=sp.dsolve(eq(sp.diff(area*mass_density*liveV,r),0),liveV)' if 'F_' in name else 'exchange=sp.simplify(-sp.diff(area*mass_density*velocity,r)/area)'),
        'K9a':'fixed_ratio_response=local_sound_speed/anchor_speed',
        'K9b':'power_response=(power_density/rho0)**s',
        'K10':'branch_domain=sp.And(observable_gate,travel_domain)',
        'K11':'advected_velocity = sp.Matrix([U,azimuthal_component,0])',
    }
    return engine_line(source,entries[knife])


def site_parts(knife,streams,source,prefix=''):
    """Read engine schemas, using semantic identities, never storage offsets.

    AST argument paths are printed as addresses only. Mapping keys, grade
    keys, named fields, matrix coordinates and source-defined branch roles
    establish identities across copies. No CAS or structural mark is used
    to choose the union. Arithmetic trees and set-valued solutions are one
    mathematical part; they are not recursively split into summands/roots.
    """
    copies=[{},{}]
    exclusions=[]
    readings=[]
    def call(node,name):return isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id==name
    def label(node):return node.args[0].value if call(node,'Str') else None
    def tup(node):return node.args if call(node,'Tuple') else None
    def ad(node):return tup(node) is not None and len(node.args)==2 and label(node.args[0])=='ATTACHED_DOMAIN'
    for copy,stream in enumerate(streams):
        external={}
        label_tag='PY_LOCAL_S9B_OBJECT_LABELS'
        if label_tag in stream:
            label_payload=stream[label_tag]
            label_tree=ast.parse(label_payload,mode='eval').body
            for row in label_tree.args:
                external[label(row.args[0])]=[(('object_label',label(row.args[0])),label_payload[n.args[1].col_offset:n.args[1].end_col_offset])
                    for n in ast.walk(row.args[1]) if ad(n)]
        for tag,payload in stream.items():
            progress(knife+' reading '+tag+' copy '+str(copy))
            name=tag.removeprefix('PY_LOCAL_S9B_').removeprefix('PY_S9B_')
            if name in ('EXPORT_CLOSURE','LOCAL_NAMES','SHARED_VOCABULARY','IMPORT_MANIFEST'):
                exclusions.append((copy,tag,(),'expression-extracted symbol, function or row names'));continue
            if name=='BUILD_INPUT_DIGESTS':
                exclusions.append((copy,tag,(),'build-input digest'));continue
            special=name in ('EXPORT_ROUNDTRIP','EXPORT_ROUTES','VELOCITY_REDUCTIONS')
            deps,entry=construction_dependencies(name)
            if not special and knife not in deps:
                exclusions.append((copy,tag,(),'computed before or independently of the mutated result'));continue
            tree=ast.parse(payload,mode='eval').body
            if special and not tree.args:
                readings.append((copy,tag,(),'empty engine mapping: no emitted value/domain entries',engine_line(source,'def cas(value):')))
                continue
            def segment(node):return payload[node.col_offset:node.end_col_offset]
            def add(node,path,identity,role,conditions,partition,origin,ranges=()):
                key=(tag,tuple(identity),partition)
                if key in copies[copy]:raise ValueError(('ambiguous semantic part',key))
                origin_name=next((part[4:].upper() for part in identity if part.startswith('key=s9b_')),name)
                copies[copy][key]={'payload_location':path,'expression':segment(node),'role':role,
                    'conditions':conditions,'density_part':partition,'source_entry':mutation_entry(knife,origin_name,source),
                    'emission_computation_source':engine_line(source,origin),'ranged_inputs':ranges}
            def walk(node,path,identity,conditions,partition,origin,context='',reading_name=None,ranges=()):
                object_name=name if reading_name is None else reading_name
                args=tup(node)
                radial_node=node if args and label(args[0])=='RADIAL_DIFFERENTIAL_IDENTITY' else None
                if args and len(args)==3 and 'IMPLIED_JN_' in object_name:
                    first=tup(args[0])
                    if first and label(first[0])=='RADIAL_DIFFERENTIAL_IDENTITY':radial_node=args[0]
                if radial_node is not None:
                    variable=segment(radial_node.args[1])
                    domain=segment(radial_node.args[2].args[1])
                    ranges=tuple(dict.fromkeys((*ranges,(variable,domain))))
                    readings.append((copy,tag,path,'full radial range from RADIAL_DIFFERENTIAL_IDENTITY; never a point substitution',variable,domain,
                        engine_line(source,"condition=cas(('RADIAL_DIFFERENTIAL_IDENTITY',r,attached_domain(sp.Interval.open(bmin,sp.oo)),attached_domain(domain),eq(GM,local_mass)))")))
                if args is None:
                    if label(node) is not None and context in ('reason','tail_domain'):
                        exclusions.append((copy,tag,path,'ineligible text: worded construction description, not a literal status or label without a mathematical reading'))
                        readings.append((copy,tag,path,'worded construction description read through its engine computation, never as a literal-difference witness',segment(node),
                            engine_line(source,'if not mass.has(Velocity):'),engine_line(source,'if mass.has(sp.Integral):'),
                            engine_line(source,'tail=sp.Integral(forcing.subs(r,t),(t,r,sp.oo))'),
                            'Supplied counting is not an attached condition; ordinary existence of the printed tail is required for any finite-value witness.'))
                        return
                    if label(node) in KNIVES[knife]:
                        exclusions.append((copy,tag,path,"text whose literal difference is the knife's replacement source"));return
                    if object_name=='MASS_BALANCE' and call(node,'Equality'):
                        add(node.args[1],path+(1,),identity+('exchange_fixed_by_mass_balance',),'value',conditions,partition,origin,ranges);return
                    if call(node,'ImmutableDenseMatrix') or call(node,'MutableDenseMatrix'):
                        entries=node.args[0].elts
                        columns=len(entries[0].elts) if entries else 0
                        readings.append((copy,tag,path,'matrix coordinates from srepr nested rows; payload address uses Matrix.args[2] row-major elements',origin))
                        for row,values in enumerate(entries):
                            for column,value in enumerate(values.elts):
                                if object_name=='VECTOR_ADVECTION' and row!=1:
                                    exclusions.append((copy,tag,path+(2,row*columns+column),'component computed independently of the azimuthal result'));continue
                                add(value,path+(2,row*columns+column),identity+(f'matrix[{row},{column}]',),'value',conditions,partition,origin,ranges)
                    else:add(node,path,identity,'text' if label(node) is not None else 'value',conditions,partition,origin,ranges)
                    return
                if args and label(args[0])=='LOCAL_BULK_DENSITY_PARTS':
                    readings.append((copy,tag,path,'density partition: one shared object, two mathematical sign trials',engine_line(source,'def density_parts(value, density):')))
                    for index,sign in ((1,'RHO_POSITIVE'),(2,'RHO_NONPOSITIVE')):
                        sign_expr=segment(args[index].args[1])
                        walk(args[3].args[1],path+(3,1),identity+('shared_density_object',),conditions,
                             (sign,sign_expr),origin,context,object_name,ranges)
                    return
                if ad(node):
                    add(args[1],path+(1,),identity+('domain',),'domain',conditions,partition,origin,ranges);return
                # Engine gated() is an exclusive branch set, not Piecewise.
                if len(args)==2 and all(tup(a) is not None and len(a.args)==2 and ad(a.args[0]) for a in args):
                    readings.append((copy,tag,path,'exclusive branch set: physical branch / complement; no first-match reading',engine_line(source,'def gated(value, domain, *, applicability=None):')))
                    for index,branch in enumerate(args):
                        role=('physical_branch','branch_complement')[index]
                        cond=(path+(index,0,1),segment(branch.args[0].args[1]))
                        add(branch.args[0].args[1],cond[0],identity+(role,'domain'),'domain',conditions,partition,origin,ranges)
                        walk(branch.args[1],path+(index,1),identity+(role,'value'),conditions+[cond],partition,origin,context,object_name,ranges)
                    return
                # Engine dictionaries: cas(Mapping), including sparse Jets.
                mapping=bool(args) and all(tup(a) is not None and len(a.args)==2 and
                    (label(a.args[0]) is not None or (tup(a.args[0]) is not None and
                     all(call(n,'Integer') for n in a.args[0].args))) for a in args)
                if mapping:
                    readings.append((copy,tag,path,'mapping: keys identify entries; missing sparse grades have no emitted-value reading',engine_line(source,'def cas(value):')))
                    for index,row in enumerate(args):
                        key=label(row.args[0]) or segment(row.args[0])
                        next_origin=origin
                        child_name=object_name
                        if object_name=='RAY_RADIAL_DOMAINS' and ((context=='FLYBY' and key in ('radial_action','angular_condition')) or
                            (knife=='K4b' and key in ('radial_action','angular_condition'))):
                            exclusions.append((copy,tag,path+(index,),'ray entry computed before the mutated derivative/result enters'));continue
                        if object_name=='ABEL_REFERENCE_INVERSES' and key=='DEFLECTION':
                            exclusions.append((copy,tag,path+(index,),'reference-only inverse; mutated result does not enter'));continue
                        if special and not identity:
                            child_name=key.removeprefix('PY_LOCAL_S9B_').removeprefix('PY_S9B_')
                            if child_name.startswith('s9b_'):child_name=child_name[4:].upper()
                            row_deps,next_origin=construction_dependencies(key)
                            if name=='VELOCITY_REDUCTIONS':
                                row_deps,next_origin=construction_dependencies(key)
                                if not any(restriction in key for restriction in ('R_SPEED_ONLY_','R_TILT_ONLY_','R_BULK_ONLY_')):
                                    row_deps=row_deps-{'K8','K10'}
                                next_origin='velocity_reductions[tag]=value[-1]'
                            if knife not in row_deps:
                                exclusions.append((copy,tag,path+(index,),'row independent of mutated result, or declaration/name row'));continue
                        walk(row.args[1],path+(index,1),identity+('key='+key,),conditions,partition,next_origin,context=key,reading_name=child_name,ranges=ranges)
                    return
                # Each direct ATTACHED_DOMAIN qualifies the container. The
                # domain part itself is tested with its other conditions.
                local=[(path+(i,1),segment(a.args[1])) for i,a in enumerate(args) if ad(a)]
                header=label(args[0]) if args else None
                schemas={
                    'APPLICABLE_OBJECT':('schema','applicability','object'),
                    'RESPONSE_ON_DEFINITION_DOMAIN':('schema','response','definition_domain'),
                    'RADIAL_DIFFERENTIAL_IDENTITY':('schema','ranged_radius','radial_interval','applicability','identity'),
                    'GAMMA_RELATION':('schema','unknown_gamma','fixing_equation','solution_branches'),
                    'GAMMA_DIFFERENCE_SOLUTION_SETS':('schema','solution_branches'),
                    'VELOCITY_FIXED_BY_RESTRICTION':('status','velocity_equation','applicability','exchange_equation'),
                    'V_UNDETERMINED':('status','reason','free_velocity','material_density','applicability','exchange_equation'),
                    'ON_EACH_CONNECTED_COMPONENT_OF_POSITIVE_V_SQUARED':('branch_label','applicability','velocity_equation','exchange_equation'),
                    'ZERO_LOCUS_DIFFERENTIABLE_JOIN_INCLUDING_ZERO_INTERVALS':('branch_label','applicability','joining_signs','zero_velocity','derivative_equation','exchange_equation'),
                }
                fields=schemas.get(header)
                if fields:readings.append((copy,tag,path,'fixed field meanings from the engine constructor',fields,engine_line(source,origin)))
                if name=='EXPORT_ROUNDTRIP' and len(identity)==1:fields=('live_value','restored_value','serialization_status')
                if object_name.startswith('RADAR_ENDPOINT_COUNTING_CERTIFICATE_') and not identity:fields=('endpoint_bound','endpoint_limit')
                if object_name=='F_MASS_SOLUTION' and not identity:fields=('zero_exchange_premise','flux_equation','velocity_equation')
                if 'IMPLIED_JN_' in object_name and header is None and len(args)==3 and context not in ('solution_branches','velocity_branches'):
                    if tup(args[0]) is not None and label(args[0].args[0])=='RADIAL_DIFFERENTIAL_IDENTITY':fields=('matching_condition','velocity_branches','velocity_reduction')
                if object_name in ('OPTICAL_Q','OPTICAL_RADIUS'):fields=None
                if object_name.startswith('C_') and object_name.endswith('_RESPONSE') and header is None and len(args)==3:
                    fields=('exact_response','linearized_delta','propagation_domain')
                if object_name.endswith('_N_DEPENDENCE') and header is None and len(args)==2 and ad(args[1]):fields=('observable_derivatives','applicability')
                if context=='observable_derivatives':fields=('deflection_derivative','radar_derivative')
                if context=='constants' and len(args)==2:fields=('anchored_constant_equation','excluded_constant_branch')
                if context=='excluded_constant_branch' and len(args)==3:fields=('excluded_domain','exclusion_status','excluded_limit_equation')
                # Source-defined branch sets use their generating equations,
                # not their list offsets, for cross-copy location identities.
                branches=bool(args) and all(tup(a) is not None and any(ad(b) for b in a.args) for a in args)
                gamma_roles=('regular mass stratum','singular compatible stratum','singular incompatible stratum')
                if len(args)==5:gamma_roles=('regular mass stratum',*(f'singular left {left} right {right}' for left in ('compatible','incompatible') for right in ('compatible','incompatible')))
                if context=='solution_branches':
                    readings.append((copy,tag,path,'branch-set identities computed from gamma_cases regular + singular_rows truth generator, or gamma_difference_cases product of left/right strata',
                        engine_line(source,'rows=[regular]+[row for group in singular_rows.values() for row in group]'),gamma_roles))
                if not fields and not branches and context!='solution_branches' and not any(ad(a) for a in args):
                    # An engine-emitted mathematical tuple is a single tuple
                    # value. Do not invent numerical field/position readings.
                    if all(label(a) is not None for a in args):
                        for index,child in enumerate(args):add(child,path+(index,),identity+('literal label='+label(child),),'text',conditions,partition,origin,ranges)
                    else:add(node,path,identity+('mathematical_tuple',),'value',conditions,partition,origin,ranges)
                    readings.append((copy,tag,path,'mathematical tuple value, not a mapping or ordered case list',origin));return
                for index,child in enumerate(args):
                    field=(fields[index] if fields and index<len(fields) else None)
                    if field in ('schema','ranged_radius','joining_signs','unknown_gamma','free_velocity','material_density'):
                        readings.append((copy,tag,path+(index,),'container identity / ranged variable, not an emitted value part',origin));continue
                    if branches:
                        first=label(child.args[0])
                        if object_name.endswith('BRANCH_TYPE'):
                            status=child.args[1]
                            field='branch type '+(label(status) or '+'.join(label(n) for n in status.args))
                        elif first=='ON_EACH_CONNECTED_COMPONENT_OF_POSITIVE_V_SQUARED':
                            rhs=child.args[2].args[1]
                            field='positive-component '+('negative-root' if segment(rhs).startswith('Mul(Integer(-1),') else 'positive-root')
                        elif first=='ZERO_LOCUS_DIFFERENTIABLE_JOIN_INCLUDING_ZERO_INTERVALS':field='zero-join signs='+segment(child.args[2])
                        elif context=='solution_branches':field=gamma_roles[index]
                        else:field='branch predicate='+segment(next(b.args[1] for b in child.args if ad(b)))
                    if field is None:
                        if context=='solution_branches':field=gamma_roles[index]
                        elif len(args)==2 and ad(args[0]):field=('branch_domain','branch_value')[index]
                        elif ad(child):field='container_domain'
                        else:field='branch_object'
                    child_conditions=conditions+[c for c in local if c[0]!=path+(index,1)]
                    if name=='EXPORT_ROUNDTRIP' and len(identity)==1 and field in ('live_value','restored_value'):
                        label_conditions=external.get('PY_S9B_'+object_name,[])
                        if label_conditions:
                            readings.append((copy,tag,path+(index,),'same object-label conditions as the shared object',object_name,label_conditions,engine_line(source,'labels={tag:')))
                            child_conditions+=label_conditions
                    walk(child,path+(index,),identity+(field,),child_conditions,partition,origin,field,object_name,ranges)
            # Large export containers inherit each source row separately.
            origin=entry or "emit('EXPORT_ROUNDTRIP', checks, local=True)"
            if name=='VELOCITY_REDUCTIONS':origin='velocity_reductions[tag]=value[-1]'
            readings.append((copy,tag,(),'parts are mathematical values/domains or literal statuses; tuples/mappings are containers',engine_line(source,origin)))
            if tag in external:
                readings.append((copy,tag,(),'object-label conditions attached to this whole object',external[tag],engine_line(source,'labels={tag:')))
            walk(tree,(),(),external.get(tag,[]),None,origin)
    union=[]
    for key in sorted(set(copies[0])|set(copies[1]),key=repr):
        pair=[copy.get(key) for copy in copies]
        union.append((key,pair))
    condition_catalog={}
    def condition_ref(raw):
        token=hashlib.sha256(raw.encode()).hexdigest()
        condition_catalog[token]=raw
        return token
    descriptions=[]
    for key,pair in union:
        descriptions.append((key,[({k:(tuple((location,condition_ref(raw)) for location,raw in v) if k=='conditions' else v)
            for k,v in p.items() if k!='expression'} if p else ('not evaluated' if knife=='K11' else 'not emitted')) for p in pair]))
    # References de-duplicate printed syntax only, not mathematical objects.
    # Every referenced condition is printed in full in this same tag.
    emit(prefix+knife+'_SITE_PARTS',{'parts':descriptions,'condition_catalog':condition_catalog,
        'exclusions':exclusions,'readings':readings,
        'identity_rule':'engine mapping keys, fixed field meanings, matrix coordinates and generating branch roles; AST paths are addresses only',
        'absent_entry_reading':('not evaluated: compact method does not reach the tag; full-copy emission is not inferred' if knife=='K11' else 'not emitted; no zero or value is invented'),
        'absent_entry_sources':(engine_line(source,'def cas(value):'),engine_line(source,'self.c = {g: v for g, v in self.c.items() if v != 0}'))})
    return union


def direct_search(knife,left_stream,right_stream,source):
    tag,selector,entry=PLANS[knife]
    a_whole,b_whole=decode(left_stream[tag]),decode(right_stream[tag])
    a,pa=select_part(a_whole,selector)
    b,pb=select_part(b_whole,selector)
    if a is None or b is None:
        return {'status':'undecided','tag':tag,'tried':('selected grade is absent from a sparse emitted mapping',selector),
                'source_entry':source_location(source,entry)}
    conditions=[attached_conditions(a_whole,pa),attached_conditions(b_whole,pb)]
    parts=(a,b)
    profiles=dict(PROFILES)
    density_conditions=tuple(whole[index][1] for whole in (a_whole,b_whole) for index in (1,2)) if selector=='response' else ()
    expressions=sp.Tuple(a,b,*(x for group in conditions for p,x in group),*density_conditions)
    fixings=[]
    # f_symbol is the response's formal local argument. The engine explicitly
    # evaluates it on f(r) in `local=response.subs(f_symbol,f)` (both copies).
    if selector=='response':
        fs=next(x for x in expressions.free_symbols if x.name=='s9b_fractional_bulk_change')
        r=next(x for x in expressions.free_symbols if x.name=='s9b_r')
        f=next(f for f in functions_in(expressions) if f.__name__=='s9b_f')
        fixings=[sp.Eq(fs,f(r),evaluate=False)]
    fixed_symbols={eq.lhs for eq in fixings}
    point={symbol:sp.sympify(POINT[symbol.name]) for symbol in expressions.free_symbols-fixed_symbols}
    # A reduced equation is not an attached condition on its own result.
    # MASS_BALANCE fixes j_n; it is handled explicitly as a solution variable.
    reading=input_reading(parts+density_conditions,conditions,profiles,point,fixings)
    reading['profile_definition_checks']={name:sp.limit(expression,RADIUS,sp.oo) for name,expression in profiles.items() if name in ('s9b_f','s9b_delta','s9b_V')}
    reading['asymptotic_bulk_density']=('rho_0 is fixed positive and f tends to zero; density is rho_0*(1+f)',sp.limit(profiles['s9b_f'],RADIUS,sp.oo))
    if knife=='K8':
        jn=a_whole.lhs
        unknown=sp.Symbol('s9b_witness_exchange_value',real=True)
        equations=[sp.Eq(unknown,evaluate(whole.rhs,profiles,point),evaluate=False) for whole in (a_whole,b_whole)]
        solutions=[sp.solveset(equation,unknown,domain=sp.S.Reals) for equation in equations]
        intersection=sp.Intersection(*solutions)
        selected=[next(iter(solution)) for solution in solutions]
        # Complete j_n profiles are fixed by the emitted equation, not by a
        # separately selected mass law or independently chosen jet values.
        fixed_profiles=[specialize_profiles(whole.rhs,profiles).subs(
            {symbol:value for symbol,value in point.items() if symbol.name!='s9b_r'},simultaneous=True)
            for whole in (a_whole,b_whole)]
        reading['symbol_roles'].append((unknown,(('fixed by construction',equations),)))
        reading['derived_profile_inputs']=tuple((jn,expression) for expression in fixed_profiles)
        reading['fixing_equations_source']=(source_location(source,'emit(\'MASS_BALANCE\',eq(jn,exchange))'),)
        emitted_evaluations=[sp.Eq(whole.lhs.subs(point),value,evaluate=False)
            for whole,value in zip((a_whole,b_whole),selected)]
        reading['solution_coordinate_identification']=sp.Eq(unknown,jn.subs(point),evaluate=False)
        reading['derived_profile_roles']='j_n is fixed by each copy of the emitted mass balance; the two derived profiles are not independent inputs.'
        values=selected
        construction_equations=(a_whole,b_whole)
    else:
        unknowns=tuple(eq.lhs for eq in fixings)
        solutions=[]
        equations=[]
        for copy in range(2):
            eqs=tuple(sp.Eq(eq.lhs,evaluate(eq.rhs,profiles,point),evaluate=False) for eq in fixings)
            equations.append(eqs)
            # Exact singleton solution tuple, or {()} for a direct construction.
            solved=tuple(sp.solveset(eq,x,domain=sp.S.Reals) for x,eq in zip(unknowns,eqs))
            solutions.append(sp.ProductSet(*solved) if solved else sp.FiniteSet(sp.Tuple()))
        intersection=sp.Intersection(*solutions)
        if intersection is sp.EmptySet:
            selected=[next(iter(x)) for x in solutions]
        else:
            common=next(iter(intersection))
            selected=[common,common]
        values=[]
        emitted_evaluations=[]
        condition_evaluations=[]
        for value,group,solution in zip(parts,conditions,selected):
            assignments=point|dict(zip(unknowns,solution))
            values.append(evaluate(value,profiles,assignments))
            emitted_evaluations.append(evaluate((a_whole,b_whole)[len(values)-1],profiles,assignments))
            condition_evaluations.append(tuple((path,condition,evaluate(condition,profiles,assignments)) for path,condition in group))
        construction_equations=tuple(sp.Eq(sp.Symbol('s9b_selected_part_value'),value,evaluate=False) for value in parts)
        reading['fixing_equations_source']=((source_location(source,'local=response.subs(f_symbol,f)'),) if fixings else ())
    if knife=='K8':
        condition_evaluations=[(),()]
    certified_conditions=all(result is sp.true for group in condition_evaluations for path,condition,result in group)
    complete_roles=all(not any(role[0]=='unresolved free input' for role in roles) for symbol,roles in reading['symbol_roles'])
    exists=all(finite_exact(value) for value in values)
    changed=different_exact(*values) if exists else False
    density_evaluations=[]
    if selector=='response':
        for whole,solution in zip((a_whole,b_whole),selected):
            assignments=point|dict(zip(unknowns,solution))
            density_evaluations.append(tuple((str(whole[index][0]),whole[index][1],
                evaluate(whole[index][1],profiles,assignments)) for index in (1,2)))
    witness=certified_conditions and complete_roles and exists and changed
    # No mathematical claim is made by the marks. Only this independently
    # computed certificate may furnish a witness, and only one is printed.
    return {'status':'witness' if witness else 'undecided','site_part':'value','tag':tag,
        'payload_locations':(pa,pb),'source_entry':source_location(source,entry),
        'construction_site':source_location(source,KNIVES[knife][0]),
        'reading':reading,'copy_fixing_equations':equations,
        'implicit_solution_variables':((unknown,) if knife=='K8' else unknowns),
        'emitted_construction_equations':construction_equations,
        'copy_solution_sets':solutions,'solution_intersection':intersection,
        'selected_solutions':selected,'solution_selection':('common solution' if intersection is not sp.EmptySet else 'intersection decided empty; each copy uses its own solution'),
        'condition_evaluations':condition_evaluations,
        'density_partition_evaluations':density_evaluations,
        'selected_parts':parts,'emitted_objects_at_inputs':emitted_evaluations,
        'exact_separation_computation':exact_separation(*values),
        'evaluations':values,'finite_exact_evaluations':exists,
        'arithmetic':'EXACT; finite ordinary integrals and derivatives evaluated by SymPy; no numerical approximation, regularization, principal value, or analytic continuation',
        'tried':{'profiles':profiles,'point':point,'exact_difference':sp.simplify(values[0]-values[1]) if all(isinstance(x,sp.Expr) for x in values) else (),
                 'all_attached_conditions_true':certified_conditions}}


def universal_coverage(pair):
    """Coverage of a symbolic identity, with no test-profile restriction."""
    expressions=sp.Tuple(*(decode(p['expression']) for p in pair if p),
        *(decode(c) for p in pair if p for path,c in p['conditions']),
        *(decode(value) for p in pair if p for row in p.get('ranged_inputs',()) for value in row))
    return {'inputs':'All independent real inputs admitted by S9b_SHARED_PHYSICS.md Parts A-C, its symbol definitions, supplied premises and this part\'s attached conditions; no additional restriction.',
        'symbol_inputs':[(s,dict(sorted(s.assumptions0.items()))) for s in sorted(expressions.free_symbols,key=sp.default_sort_key)],
        'live_profile_inputs':sorted((str(f),tuple(sorted({a.args for a in expressions.atoms(AppliedUndef) if a.func==f},key=str))) for f in functions_in(expressions)),
        'bound_occurrences':binders(expressions),'stated_intervals':sorted(expressions.atoms(sp.Interval),key=sp.default_sort_key),
        'declared_ranged_inputs_by_copy':[[(decode(variable),decode(domain)) for variable,domain in p.get('ranged_inputs',())] if p else 'not emitted' for p in pair],
        'profiles':'Every live general radial delta, V, xi, rho_br, f, mu_perp, j_n admitted by those definitions, with required derivatives and ordinary integrals/limits; rho_0 is the positive asymptotic bulk density and f tends to zero.',
        'intervals':'The intervals in the printed part and attached conditions, including the entire r>b_min range of a radial differential identity. Bound variables and quantifiers remain active.',
        'parts_and_conditions':[(p['expression'],p['conditions']) if p else 'not emitted' for p in pair],
        'K11':'No sample-domain or truncation restriction on a no-change decision.'}


def finite_attempt(key,pair,source,knife):
    """One finite symbolic decision and one concrete-profile trial per part.

    A symbolic reflexivity decision covers general inputs. A failed or
    unresolved concrete trial cannot decide no change or impossibility.
    No joint ray solver or new physics construction is introduced here.
    """
    tag,identity,partition=key
    missing='not evaluated' if knife=='K11' else 'not emitted'
    base={'tag':tag,'location_identity':identity,'density_partition':partition,
        'locations':[p['payload_location'] if p else missing for p in pair],
        'roles':[p['role'] if p else missing for p in pair],
        'source_entries':[p['source_entry'] if p else missing for p in pair]}
    if knife=='K11' and tag not in COMPACT_TARGETS:
        return base|{'status':'undecided','CAS_operations':('compact method does not reach this tag; no mathematical comparison',),
            'inputs':'not evaluated','evaluations':('not evaluated','not evaluated'),
            'unmet_requirements':('Both compact operands are not evaluated. No general no-change decision and no witness can be supplied.')}
    values=[decode(p['expression']) if p else None for p in pair]
    conditions=[[(path,decode(raw)) for path,raw in p['conditions']] if p else [] for p in pair]
    range_rows=[[(decode(variable),decode(domain)) for variable,domain in p.get('ranged_inputs',())] if p else [] for p in pair]
    base['declared_ranged_inputs_by_copy']=range_rows
    # Literal engine statuses have no mathematical reading. In particular,
    # do not turn a status spelling into a numerical or physical predicate.
    if knife!='K11' and all(p and p['role']=='text' for p in pair):
        result=sp.sympify(values[0]==values[1])
        return base|{'status':'no change' if result is sp.true else 'witness',
            'CAS_operations':('SymPy Str literal equality',),'texts':values,'CAS_result':result,
            'statement':('literal status equality',values),
            'full_input_coverage':'These emitted literal statuses have no mathematical reading or mathematical inputs. Their containing mathematical domains are separate enumerated parts.',
            'containing_condition_locations':[[path for path,condition in group] for group in conditions],
            'counts_as_density_sign_trial':False,'inputs':'No mathematical inputs for these literal statuses.',
            'evaluations':values,'unmet_requirements':()}
    statement=None
    decision=Str('UNRESOLVED')
    # This CAS equality compares mathematical objects AND every containing
    # condition. It is separate from, and does not inspect, literal marks.
    if knife!='K11' and all(p and p['role']!='text' for p in pair):
        operands=[sp.Tuple(value,*(condition for path,condition in group),sp.Tuple(*(sp.Tuple(*row) for row in ranges)))
            for value,group,ranges in zip(values,conditions,range_rows)]
        statement=sp.Eq(*operands,evaluate=False)
        decision=sp.Eq(*operands)
        if decision is sp.true and not partition:
            return base|{'status':'no change','statement':statement,'CAS_result':decision,
                'CAS_operations':('SymPy Eq on part and attached-condition tuples; general symbolic inputs',),
                'full_input_coverage':universal_coverage(pair),'unmet_requirements':()}
    # A concrete, fixed, asymptotically admissible test, not a family. All
    # quantifiers stay active; partial substitutions are explicitly partial.
    profiles=dict(PROFILES)
    if knife=='K11':
        profiles.update({'s9b_delta':COMPACT_PARAMETERS['d']/RADIUS**COMPACT_PARAMETERS['p'],
            's9b_V':COMPACT_PARAMETERS['v']/RADIUS**COMPACT_PARAMETERS['q'],
            's9b_xi':COMPACT_PARAMETERS['w']*sp.log(RADIUS)})
    if partition and partition[0]=='RHO_NONPOSITIVE':profiles['s9b_f']=-8/RADIUS**2
    wholes=sp.Tuple(*(v for v in values if v is not None),*(c for group in conditions for path,c in group),
        *(value for rows in range_rows for row in rows for value in row),
        *((decode(partition[1]),) if partition else ()))
    free=wholes.free_symbols
    # Reduced identities range over the whole engine interval; their radius
    # cannot be converted into a single test coordinate.
    ranged={variable for rows in range_rows for variable,domain in rows}
    input_point=POINT|({'x_1':COMPACT_POINTS[0][0],'x_2':COMPACT_POINTS[0][1],'x_3':COMPACT_POINTS[0][2],
        's9b_azimuthal':COMPACT_AZIMUTHAL} if knife=='K11' else {})
    point={s:sp.sympify(input_point[s.name]) for s in free-ranged if s.name in input_point}
    profile_rules=lambda x:profiles[x.func.__name__].subs(RADIUS,x.args[0])
    partial=[]
    for value,group in zip(values,conditions):
        def prepare(expr):
            # Substitution evaluates algebra, but deliberately does not
            # claim that unevaluated integrals/ConditionSets were decided.
            return expr.replace(lambda x:isinstance(x,AppliedUndef) and x.func.__name__ in profiles,profile_rules).subs(point,simultaneous=True)
        partial.append({'part':prepare(value) if value is not None else missing,
            'conditions':[(path,c,prepare(c)) for path,c in group]})
    density=[]
    if partition:
        for copy,p in enumerate(pair):
            if p is None:density.append(missing);continue
            condition=decode(partition[1])
            tested=condition.replace(lambda x:isinstance(x,AppliedUndef) and x.func.__name__ in profiles,profile_rules)
            symbols=tested.free_symbols
            density_point={s:sp.sympify(POINT[s.name]) for s in symbols if s.name in POINT and s.name!='s9b_r'}
            tested=tested.subs(density_point,simultaneous=True)
            r=next((s for s in tested.free_symbols if s.name=='s9b_r'),None)
            density_value=condition.lhs.replace(lambda x:isinstance(x,AppliedUndef) and x.func.__name__ in profiles,profile_rules).subs(density_point,simultaneous=True)
            if r is None or r not in ranged:
                assignments=({r:point[r]} if r is not None else {})
                density.append({'condition':condition,'point':density_point|assignments,
                    'density_evaluation':density_value.subs(assignments,simultaneous=True),
                    'CAS_operations':'subs at the same printed concrete input',
                    'sign_met':sp.simplify(tested.subs(assignments,simultaneous=True))})
            else:
                domain=sp.Interval.open(1,2) if partition[0]=='RHO_NONPOSITIVE' else sp.Interval.open(3,sp.oo)
                for variable,stated in range_rows[copy]:
                    if variable==r:domain=sp.Intersection(domain,stated.subs(point,simultaneous=True))
                solution=sp.solveset(tested,r,domain=domain)
                density.append({'condition':condition,'at_profile':tested,'subrange':domain,'CAS_operations':'solveset(sign, r, domain=subrange)',
                    'sign_solution_set':solution,'covers_subrange':sp.Eq(solution,domain),
                    'nonempty_subrange':sp.sympify(domain.is_empty is False),
                    'density':density_value,'rho_0':sp.S.One})
    unassigned=sorted(free-set(point)-ranged,key=sp.default_sort_key)
    unfixed=sorted({f.__name__ for f in functions_in(wholes)}-set(profiles))
    no_change=decision is sp.true
    return base|{'status':'no change' if no_change else 'undecided','statement':statement,'CAS_result':decision,
        'full_input_coverage':universal_coverage(pair) if no_change else 'No full-input decision.',
        'CAS_operations':('Eq on general part and condition tuples (when both mathematical parts emitted)',
            'replace fixed profile functions; simultaneous free-input subs; no bound-variable substitution; active quantifiers retained'),
        'inputs':{'point':point,'ranged':[[ (variable,domain.subs(point,simultaneous=True)) for variable,domain in rows] for rows in range_rows],
            'fixed_profiles':[(name,sp.Lambda(RADIUS,expr)) for name,expr in profiles.items()],
            'f_asymptotic_limit':sp.limit(profiles['s9b_f'],RADIUS,sp.oo),'rho_0':1,
            'bound_occurrences':binders(wholes),'unassigned_inputs':unassigned,'unfixed_profiles':unfixed},
        'evaluations':partial,'density_sign_evaluations':density,
        'unmet_requirements':() if no_change else ('No agreement decision at every admissible input. No impossibility decision.',
            'This partial trial does not solve the joint fixing equations and their cross-copy intersection, decide every active quantified condition, or certify ordinary finite unequal values. No witness is inferred.')}


def finite_search(knife,left,right,source,folder):
    parts=site_parts(knife,(left,right),source)
    if knife=='K7':
        # A reduction kind is an engine-emitted literal status, not a
        # mathematical statement about the solution or the ray gates.
        parts.sort(key=lambda item:(not (item[0][0]=='PY_LOCAL_S9B_VELOCITY_REDUCTIONS' and 'key=kind' in item[0][1]),repr(item[0])))
    # The preferred construction is only a search ordering, never a filter.
    # It is attempted as an individual part, not as a tuple/container value.
    if knife not in ('K7','K10'):
        trial=direct_search(knife,left,right,source)
        signs=set.intersection(*(set(row[0] for row in rows if row[2] is sp.true) for rows in trial.get('density_partition_evaluations',()))) if trial.get('density_partition_evaluations') else set()
        identities=[key for key,pair in parts if key[0]==trial.get('tag') and
            len(trial.get('payload_locations',()))==2 and all(p and p['payload_location']==trial['payload_locations'][i] for i,p in enumerate(pair))
            and (not key[2] or key[2][0] in signs)]
        trial['enumerated_site_parts']=identities
        emit(knife+'_ATTEMPT_PREFERRED',trial)
        if trial['status']=='witness':
            if not identities:raise ValueError(('witness part not enumerated',knife))
            emit(knife+'_WITNESS',trial)
            return {'status':'witness','tag':trial['tag'],'site_part':'value','site_part_count':len(parts),'attempt_count':1}
    unresolved=[]
    counts=Counter()
    for number,(key,pair) in enumerate(parts):
        progress(knife+' finite search '+key[0]+' '+repr(key[1]))
        trial=finite_attempt(key,pair,source,knife)
        emit(knife+'_ATTEMPT_'+str(number),trial)
        counts[trial['status']]+=1
        if trial['status']=='witness':
            emit(knife+'_WITNESS',trial)
            return {'status':'witness','tag':key[0],'site_part':'text','site_part_count':len(parts),'attempt_count':sum(counts.values())}
        if trial['status']!='no change':unresolved.append((key,trial['status'],trial['unmet_requirements']))
    outcome=search_outcome(['undecided'] if unresolved else ['no change'])
    emit(knife+'_SEARCH_OUTCOME',{'status':outcome,'parts_without_no_change_decision':unresolved,
        'attempt_counts':dict(counts),'site_part_count':len(parts),
        'statement':'none found requires a no-change decision at every site part; samples and structural marks decide neither agreement nor impossibility'})
    return {'status':outcome,'site_part_count':len(parts),'attempt_count':sum(counts.values()),
            'parts_without_no_change_decision':len(unresolved)}


def search_outcome(decisions):
    return 'none found' if all(decision=='no change' for decision in decisions) else 'undecided'


def compact_copy(path):
    namespace={'__file__':str(path),'__name__':'s9b_compact_copy'}
    source=path.read_text()
    exec(compile(source,str(path),'exec'),namespace)
    r=sp.Symbol('s9b_r',positive=True)
    config=COMPACT_PARAMETERS
    templates={namespace['Delta']:config['d']*(config['L']/r)**config['p'],
               namespace['Velocity']:config['c0']*config['v']*(config['L']/r)**config['q'],
               namespace['Embedding']:config['w']*sp.log(r/config['L'])}
    replacements={function:(lambda x,value=value:value.subs(r,x)) for function,value in templates.items()}
    result=namespace['nonreciprocity'](r,config['c0'],templates[namespace['Delta']],
        templates[namespace['Velocity']],templates[namespace['Embedding']])
    parameters={sp.Symbol('c_0',positive=True):config['c0'],
                sp.Symbol('s9b_azimuthal',real=True):COMPACT_AZIMUTHAL}
    reached={}
    for tag,name in COMPACT_TARGETS.items():
        rows=[]
        for point in COMPACT_POINTS:
            progress('K11 '+path.name+' '+tag+' sample '+str(point))
            rows.append(sp.simplify(namespace['cas'](result[name]).subs(parameters).subs(
                dict(zip(result['xyz'],map(sp.Integer,point))))))
        reached[tag]=sp.Tuple(*rows)
    return reached,(namespace['profile_substitute'],replacements,parameters),result,templates


def compact_reduce(value,reducer,coordinates):
    apply,replacements,parameters=reducer
    reduced=apply(value,replacements).subs(parameters)
    return sp.Tuple(*(sp.simplify(reduced.subs(dict(zip(coordinates,map(sp.Integer,point)))))
                      for point in COMPACT_POINTS))


def compact_k11(baseline,source,mutant,folder):
    unmutated=mutant.with_name('S9b_compact_unmutated.py')
    unmutated.write_text(source)
    original,reducer,result,templates=compact_copy(unmutated)
    corrupted,other_reducer,other_result,other_templates=compact_copy(mutant)
    emit('K11_CONFIGURATION_INPUTS',{'unmutated':(result['xyz'],result['grades'],templates),
        'corrupted':(other_result['xyz'],other_result['grades'],other_templates)})
    if result['xyz']!=other_result['xyz'] or result['grades']!=other_result['grades'] or templates!=other_templates:
        raise ValueError('K11 configuration mismatch')
    sites=site_parts('K11',(baseline,{tag:render(cas(other_result[name])) for tag,name in COMPACT_TARGETS.items()}),source)
    coverage={'method':'engine vector_dispersion/nonreciprocity; differentiate before coordinate sampling',
        'profile_parameters':COMPACT_PARAMETERS,'profile_functions':templates,
        'azimuthal_parameter':COMPACT_AZIMUTHAL,'sample_points':COMPACT_POINTS,
        'retained_grades':result['grades'],'additional_truncation':False,
        'arithmetic':'EXACT; no rounding, seed, or numerical error',
        'sample_domain_is_coverage_not_condition':True}
    emit('K11_CONFIGURATION',coverage)
    reached=[]
    not_evaluated=[]
    flags=[]
    for tag,full in sorted(baseline.items()):
        progress('K11 '+tag)
        if tag in original:
            reduced=compact_reduce(decode(full),reducer,result['xyz'])
            before,after,reduced_text=render(original[tag]),render(corrupted[tag]),render(reduced)
            # Ordered operands, never subtraction or a zero/change category.
            delta='Tuple('+render(Str('COMPACT_CORRUPTED'))+', '+after+', '+render(Str('COMPACT_UNMUTATED'))+', '+before+', '+render(Str(structural_mark(after,before)))+', '+render(cas(coverage))+')'
            residual='Tuple('+render(Str('COMPACT_UNMUTATED'))+', '+before+', '+render(Str('REDUCED_FULL_BASELINE'))+', '+reduced_text+', '+render(Str(structural_mark(before,reduced_text)))+', '+render(cas(coverage))+')'
            reached.append(tag)
            flags.append({'tag':tag,'evaluated':True,'difference_pair_mark':structural_mark(after,before),
                'method_residual_pair_mark':structural_mark(before,reduced_text),
                'difference_roles':domain_literal_flags(after,before),
                'method_residual_roles':domain_literal_flags(before,reduced_text)})
        else:
            missing=render(Str('NOT_EVALUATED'))
            delta='Tuple('+missing+', '+missing+', '+render(Str('not compared'))+', '+render(cas(coverage))+')'
            residual=delta
            not_evaluated.append(tag)
            flags.append({'tag':tag,'evaluated':False,'difference_pair_mark':'not compared',
                'method_residual_pair_mark':'not compared'})
        print('PY_LOCAL_S9B_ABLATION_K11_'+tag+': Tuple('+full+', Tuple('+render(Str('DIFFERENCE_PAIR'))+', '+delta+'), Tuple('+render(Str('METHOD_RESIDUAL_PAIR'))+', '+residual+'))',flush=True)
    tag=PLANS['K11'][0]
    # Witness only the corrupted/unmutated pair, never the method residual.
    values=(original[tag][0][0,0],corrupted[tag][0][0,0])
    profiles={f.__name__:value.subs(sp.Symbol('s9b_r',positive=True),RADIUS) for f,value in templates.items()}
    symbolic=(decode(baseline[tag])[0,0],other_result['oneform'][0,0])
    expressions=sp.Tuple(*symbolic)
    point={symbol:sp.sympify({'x_1':3,'x_2':4,'x_3':0,'c_0':1,'s9b_azimuthal':COMPACT_AZIMUTHAL}[symbol.name])
        for symbol in expressions.free_symbols}
    reading=input_reading(symbolic,[[],[]],profiles,point)
    reading['configuration_fixed_quantities']=coverage
    reading['copy_azimuthal_constructions']=(sp.S.Zero,COMPACT_AZIMUTHAL)
    # Re-evaluate the full baseline component at the same inputs; this is a
    # reading computation, not the method residual used as a witness.
    evaluations=(evaluate(symbolic[0],profiles,point),evaluate(symbolic[1],profiles,point))
    witness=(all(finite_exact(x) for x in values+evaluations) and different_exact(*values)
        and all(sp.simplify(a-b)==0 for a,b in zip(values,evaluations)))
    search={'status':'witness' if witness else 'undecided',
        'site_part':'value','tag':tag,'payload_location':'DIFFERENCE_PAIR compact operand: sample[0], matrix[0,0]',
        'source_entry':source_location(source,PLANS['K11'][2]),
        'construction_site':source_location(source,KNIVES['K11'][0]),
        'reading':reading,'copy_fixing_equations':((),()),
        'copy_solution_sets':(sp.FiniteSet(sp.Tuple()),sp.FiniteSet(sp.Tuple())),
        'solution_intersection':sp.FiniteSet(sp.Tuple()),'selected_solutions':((),()),
        'solution_selection':'common empty tuple of implicit unknowns; both direct constructions',
        'condition_evaluations':((),()),'evaluations':values,
        'reading_computations':evaluations,'coverage':coverage,
        'arithmetic':'EXACT; ordinary smooth functions on r>0 and transverse radius>0',
        'enumerated_site_parts':[key for key,pair in sites if key[0]==tag and 'matrix[0,0]' in key[1]],
        'tried':('same compact inputs in both copies; first sample and matrix component',)}
    emit('K11_ATTEMPT_PREFERRED',search)
    if witness:emit('K11_WITNESS',search)
    else:
        undecided=[]
        for number,(key,pair) in enumerate(sites):
            trial=finite_attempt(key,pair,source,'K11')
            # Fixed compact profiles cannot certify every admissible input.
            trial['status']='undecided'
            emit('K11_ATTEMPT_'+str(number),trial)
            undecided.append(key)
        emit('K11_SEARCH_OUTCOME',{'status':'undecided','parts_without_no_change_decision':undecided})
    (folder/'marks.json').write_text(json.dumps({'tag_count':len(flags),'tags':flags},indent=2)+'\n')
    receipt={'evaluated_tags':reached,'not_evaluated_tags':not_evaluated,
             'source_sha256':hashlib.sha256(source.encode()).hexdigest(),
             'mutant_sha256':sha(mutant),'harness_sha256':sha(HERE),
             'witness_summary':{'status':'witness' if witness else 'undecided','tag':tag,'site_part':'value','site_part_count':len(sites)}}
    (folder/'compact-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


def run(script,folder):
    folder.mkdir(parents=True,exist_ok=False)
    output=folder/'stdout'
    with output.open('w') as out,(folder/'stderr').open('w') as err:
        child=subprocess.Popen([sys.executable,'-B',str(script)],stdout=out,stderr=subprocess.PIPE,
            text=True,stdin=subprocess.DEVNULL,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
        for line in child.stderr:
            err.write(line);err.flush()
            progress(str(folder)+': '+line.rstrip('\n'))
        code=child.wait()
        child.stderr.close()
    (folder/'exit-code').write_text(str(code)+'\n')
    if code:raise RuntimeError(('worker stopped',str(folder),code))
    return output


def self_checks():
    # Every fixture's actual computed mark/value precedes its assertion.
    checks={}
    def check(name,computed,predicate):
        emit('SOFTWARE_FIXTURE_'+name,computed)
        passed=bool(predicate(computed))
        emit('SOFTWARE_RESULT_'+name,passed)
        checks[name]=passed
        if not passed:raise ValueError(('software fixture',name))
    fixtures=(
        ('syntax_whitespace','Tuple(Integer(1), Integer(2))','Tuple( Integer(1),Integer(2) )','identical'),
        ('argument_order','Tuple(Integer(1), Integer(2))','Tuple(Integer(2), Integer(1))','not compared'),
        ('no_algebra','Add(Integer(1),Integer(1))','Integer(2)','not compared'),
        ('float_decimal_tail',"Float('1.00000000000000001', precision=100)","Float('1.00000000000000002', precision=100)",'not compared'),
        ('float_precision',"Float('1.0', precision=53)","Float('1.0', precision=54)",'not compared'),
        ('python_float_token','f(1.00000000000000001)','f(1.00000000000000002)','not compared'),
        ('binder_names_literal',"Integral(Symbol('x'), (Symbol('x'),0,1))","Integral(Symbol('y'), (Symbol('y'),0,1))",'not compared'))
    for name,left,right,expected in fixtures:
        check(name,{'left':left,'right':right,'mark':structural_mark(left,right)},lambda value: value['mark']==expected)
    x=sp.Symbol('self_check_x',real=True)
    function=sp.Function('s9b_f',real=True)
    outer=sp.Tuple(Str('LOCAL_BULK_DENSITY_PARTS'),sp.Tuple(Str('RHO_POSITIVE'),sp.Gt(x,0)),
        sp.Tuple(Str('RHO_NONPOSITIVE'),sp.Le(x,0)),sp.Tuple(Str('OBJECT_ON_BOTH_PARTS'),
            sp.Tuple(sp.Tuple(Str('ATTACHED_DOMAIN'),sp.Ge(x,0)),x)))
    check('density_labels_are_not_conditions',attached_conditions(outer,(3,1,1)),lambda value:len(value)==1)
    check('derivative_before_point',evaluate(sp.diff(function(x),x),PROFILES,{x:2}),lambda value:value==-sp.Rational(1,40))
    check('integral_bound_not_substituted',evaluate(sp.Integral(function(x),(x,1,2)),PROFILES,{x:99}),lambda value:value==sp.Rational(1,20))
    check('empty_unknown_tuple_is_a_solution',sp.Intersection(sp.FiniteSet(sp.Tuple()),sp.FiniteSet(sp.Tuple())),lambda value:value==sp.FiniteSet(sp.Tuple()))
    check('finite_tuple_witness',{'values':(sp.Tuple(1),sp.Tuple(2)),'finite':finite_exact(sp.Tuple(1,sp.Tuple(2,3))),
        'different':different_exact(sp.Tuple(1),sp.Tuple(2))},lambda value:value['finite'] and value['different'])
    check('unevaluated_not_finite',{'value':sp.Integral(function(x),(x,1,sp.oo)),
        'finite':finite_exact(sp.Integral(function(x),(x,1,sp.oo)))},lambda value:not value['finite'])
    source=ENGINE.read_text()
    r=sp.Symbol('s9b_r',positive=True)
    fs=sp.Symbol('s9b_fractional_bulk_change',real=True)
    f=function(r)
    def response_fixture(response):
        return sp.Tuple(Str('LOCAL_BULK_DENSITY_PARTS'),
            sp.Tuple(Str('RHO_POSITIVE'),sp.Gt(1+f,0)),
            sp.Tuple(Str('RHO_NONPOSITIVE'),sp.Le(1+f,0)),
            sp.Tuple(Str('OBJECT_ON_BOTH_PARTS'),sp.Tuple(
                sp.Tuple(Str('RESPONSE_ON_DEFINITION_DOMAIN'),response,
                    sp.Tuple(Str('ATTACHED_DOMAIN'),sp.Ge(1+fs,0))),
                f,sp.Tuple(Str('ATTACHED_DOMAIN'),sp.Gt(1+f,0)))))
    tag=PLANS['K9a'][0]
    check('response_solution_intersection',direct_search('K9a',{tag:render(response_fixture(1+fs))},
        {tag:render(response_fixture(sp.S.One))},source),lambda value:value['status']=='witness'
        and value['solution_intersection'].is_empty is False and value['selected_solutions'][0]==value['selected_solutions'][1]
        and all(row[2] is sp.true for group in value['condition_evaluations'] for row in group))
    tag=PLANS['K4a'][0]
    grade=sp.Tuple(*PLANS['K4a'][1])
    left={tag:render(sp.Tuple(sp.Tuple(grade,f)))}
    right={tag:render(sp.Tuple(sp.Tuple(grade,sp.S.Zero)))}
    check('numeric_mapping_part_reading',direct_search('K4a',left,right,source),lambda value:value['status']=='witness' and value['selected_parts'][0]==f)
    tag=PLANS['K8'][0]
    jn=sp.Function('j_n',real=True)(r)
    check('decided_disjoint_solutions',direct_search('K8',{tag:render(sp.Eq(jn,f,evaluate=False))},
        {tag:render(sp.Eq(jn,0,evaluate=False))},source),lambda value:value['status']=='witness'
        and value['solution_intersection'] is sp.EmptySet and all(solution.contains(v) is sp.true for solution,v in zip(value['copy_solution_sets'],value['selected_solutions'])))
    parts=site_parts('K1',(left,right),source,prefix='SOFTWARE_FIXTURE_')
    check('mapping_semantic_identity',[(key,[p['payload_location'] if p else 'not emitted' for p in pair]) for key,pair in parts],
        lambda value:len(value)==1 and value[0][0][1][0].startswith('key='))
    equal=site_parts('K1',(left,left),source,prefix='SOFTWARE_EQUAL_')
    check('symbolic_no_change_coverage',finite_attempt(*equal[0],source,'K1'),lambda value:value['status']=='no change' and 'full_input_coverage' in value)
    absent=site_parts('K1',(left,{}),source,prefix='SOFTWARE_ABSENT_')
    check('one_copy_is_not_a_witness',finite_attempt(*absent[0],source,'K1'),lambda value:value['status']=='undecided' and value['evaluations'][1]['part']=='not emitted')
    sign=site_parts('K9a',({PLANS['K9a'][0]:render(response_fixture(1+fs))},
        {PLANS['K9a'][0]:render(response_fixture(1+fs))}),source,prefix='SOFTWARE_SIGNS_')
    check('both_density_signs_enumerated',[key[2][0] for key,pair in sign],lambda value:set(value)=={'RHO_POSITIVE','RHO_NONPOSITIVE'})
    negative=next((key,pair) for key,pair in sign if key[2][0]=='RHO_NONPOSITIVE' and all(p['role']=='value' for p in pair))
    check('negative_sign_trial_is_computed',finite_attempt(*negative,source,'K9a'),lambda value:value['status']=='no change'
        and all(row.get('sign_met') is sp.true or row.get('covers_subrange') is sp.true for row in value['density_sign_evaluations']))
    compact_key=('PY_LOCAL_S9B_NONRECIPROCITY_ONEFORM',('matrix[0,0]',),None)
    check('compact_samples_do_not_decide_no_change',finite_attempt(compact_key,equal[0][1],source,'K11'),
        lambda value:value['status']=='undecided' and value['CAS_result']==Str('UNRESOLVED'))
    check('unreached_compact_part_is_not_evaluated',finite_attempt(*equal[0],source,'K11'),
        lambda value:value['status']=='undecided' and value['evaluations']==('not evaluated','not evaluated'))
    check('outcome_requires_every_no_change',[(statuses,search_outcome(statuses)) for statuses in
        (['no change'],['no change','undecided'],['no change','cannot be tried'])],lambda value:[row[1] for row in value]==['none found','undecided','undecided'])
    check('negative_density_profile_asymptote',{'f':-8/RADIUS**2,'limit':sp.limit(-8/RADIUS**2,RADIUS,sp.oo),
        'local_density':(1-8/RADIUS**2).subs(RADIUS,2)},lambda value:value['limit']==0 and value['local_density']<0)
    # The same mathematical object keeps its field readings and external
    # conditions when publish() embeds it in a serialization-check row.
    gm=sp.Symbol('GM',real=True)
    material=sp.Function('rho_br_profile',positive=True)(r)
    velocity=sp.Function('s9b_V',real=True)(r)
    matching=cas(('RADIAL_DIFFERENTIAL_IDENTITY',r,('ATTACHED_DOMAIN',sp.Interval.open(1,sp.oo)),
        ('ATTACHED_DOMAIN',sp.true),sp.Eq(gm,f,evaluate=False)))
    implied=cas((matching,('V_UNDETERMINED','software fixture description',velocity,material,
        ('ATTACHED_DOMAIN',sp.true),sp.Eq(jn,0,evaluate=False)),{'kind':'UNDETERMINED'}))
    shared_tag='PY_S9B_B_IMPLIED_JN_DEFLECTION'
    forward_tag='PY_S9B_F_LIVE_B_CONDITION_DEFLECTION'
    export_tag='PY_LOCAL_S9B_EXPORT_ROUNDTRIP'
    fixture_stream={shared_tag:render(implied),forward_tag:render(matching),
        export_tag:render(cas({'s9b_b_implied_jn_deflection':(implied,implied,sp.true),
            's9b_f_live_b_condition_deflection':(matching,matching,sp.true)})),
        'PY_LOCAL_S9B_OBJECT_LABELS':render(cas({forward_tag:('FORWARD_NO_FAR_ZONE_LOSS',('ATTACHED_DOMAIN',sp.Eq(jn,0,evaluate=False)))}))}
    exported=site_parts('K7',(fixture_stream,fixture_stream),source,prefix='SOFTWARE_EXPORTED_')
    def shared_view(shared_name,export_key,operand):
        shared=[]
        recovered=[]
        for key,pair in exported:
            tag,identity,partition=key
            value=pair[0]
            summary=(partition,value['role'],value['expression'],tuple(raw for path,raw in value['conditions']),value['ranged_inputs'])
            if tag==shared_name:shared.append((identity,summary))
            if tag==export_tag and identity[:2]==('key='+export_key,operand):recovered.append((identity[2:],summary))
        return (sorted(shared,key=repr),sorted(recovered,key=repr))
    check('exported_implied_object_keeps_field_readings',{
        operand:shared_view(shared_tag,'s9b_b_implied_jn_deflection',operand) for operand in ('live_value','restored_value')},
        lambda value:all(original==recovered and original for original,recovered in value.values()))
    check('exported_object_keeps_label_conditions',{
        operand:shared_view(forward_tag,'s9b_f_live_b_condition_deflection',operand) for operand in ('live_value','restored_value')},
        lambda value:all(original==recovered and original for original,recovered in value.values()))
    check('serialization_status_has_no_physics_condition',[
        pair[0]['conditions'] for key,pair in exported if key[0]==export_tag and key[1][-1]=='serialization_status'],
        lambda value:value and all(not conditions for conditions in value))
    ranged_key,ranged_pair=next((key,pair) for key,pair in exported if key[0]==forward_tag and key[1][-1]=='identity')
    other=dict(ranged_pair[1],expression=render(sp.Eq(gm,0,evaluate=False)))
    check('radial_identity_is_never_point_substituted',finite_attempt(ranged_key,(ranged_pair[0],other),source,'K7'),
        lambda value:value['status']=='undecided' and r not in value['inputs']['point']
        and all(rows and rows[0]==(r,sp.Interval.open(1,sp.oo)) for rows in value['declared_ranged_inputs_by_copy']))
    literal={'payload_location':(),'expression':render(Str('fixture status A')),'role':'text',
        'conditions':[((0,),render(sp.Gt(x,0)))],'source_entry':source_location(source,KNIVES['K7'][0]),'ranged_inputs':()}
    text_key=('PY_LOCAL_S9B_SOFTWARE_TEXT',('literal_status',),None)
    check('text_equality_is_literal_even_with_containers',finite_attempt(text_key,(literal,literal),source,'K7'),
        lambda value:value['status']=='no change' and value['CAS_result'] is sp.true and not value['counts_as_density_sign_trial'])
    changed_literal=dict(literal,expression=render(Str('fixture status B')))
    check('text_difference_is_literal_even_with_containers',finite_attempt(text_key,(literal,changed_literal),source,'K7'),
        lambda value:value['status']=='witness' and value['CAS_result'] is sp.false and not value['counts_as_density_sign_trial'])
    emit('SOFTWARE_CHECKS',checks)


def raw_record(path,source_path,run_id):
    rows=[]
    marks=Counter()
    seen=set()
    with path.open() as stream:
        for number,line in enumerate(stream,1):
            tag,sep,payload=line.rstrip('\n').partition(': ')
            progress('inventory '+run_id+' '+tag)
            if not sep or not tag.startswith(('PY_S9B_','PY_LOCAL_S9B_')) or tag in seen:raise ValueError(('invalid tag',path,number))
            seen.add(tag)
            tree=ast.parse(payload,mode='eval')
            strings=[node.value for node in ast.walk(tree) if isinstance(node,ast.Constant) and isinstance(node.value,str)]
            markers=Counter(word for word in strings if word in ('NOT_ESTABLISHED','NOT_EVALUATED','ATTACHED_DOMAIN',
                'identical','not compared','witness','none found','undecided','no change','cannot be tried','not emitted','not evaluated','UNRESOLVED'))
            marks.update(markers)
            rows.append({'name':tag,'line':number,'markers':dict(sorted(markers.items()))})
    return {'id':run_id,'output':str(path.relative_to(REPOSITORY)),
        'output_sha256':sha(path),'source':str(source_path.relative_to(REPOSITORY)),
        'source_sha256':sha(source_path),'tag_count':len(rows),'marker_counts':dict(sorted(marks.items())),
        'tags':rows}


def make_inventory(scratch,runs,pins,k11,command):
    records=[raw_record(path,source,run_id) for run_id,path,source in runs]
    baseline=read_stream(runs[0][1])
    comparison_counts={}
    # Sidecars contain only structural marks and names. No payloads enter the
    # inventory. Computed witness details remain in the run's stdout stream.
    for knife in KNIVES:
        comparison_counts[knife]=json.loads((scratch/knife/'marks.json').read_text()) if knife!='K11' else {
            'tag_count':len(baseline),'evaluated_count':len(k11['evaluated_tags']),
            'not_evaluated_count':len(k11['not_evaluated_tags']),
            'tags':json.loads((scratch/knife/'marks.json').read_text())['tags']}
    inventory={'status':'BUILD_RUN_COMPLETE_STOPPED_FOR_HANDOFF','commits':[],
        'authority':'S9b repair build directive amendment 12, repair round 9; Parts A-C shared physics',
        'builder':'Codex gpt-6.1-sol at xhigh effort; fresh harness author',
        'inventory_format':'Names, provenance, markers, structural flags and counts only; no tag payload copied.',
        'structural_comparison':{'method_id':'literal_ast_with_original_constant_tokens','source':str(HERE.relative_to(REPOSITORY)),'function':'literal_key','source_sha256':sha(HERE)},
        'flag_qualification':'identical is literal emitted-expression structure only; not compared makes no equality or change claim. Value/domain flags are literal marks, not mathematical differences. Finite site-part search decisions are separate.',
        'source_pins':pins,'runs':records,'command':command,
        'comparison_counts_by_knife':comparison_counts,
        'shared_tags':sorted(tag for tag in baseline if tag.startswith('PY_S9B_')),
        'local_tags':sorted(tag for tag in baseline if tag.startswith('PY_LOCAL_S9B_')),
        'not_established_tags':sorted(tag for tag,payload in baseline.items() if "Str('NOT_ESTABLISHED')" in payload),
        'missing_shared':[str(name) for name in decode(baseline['PY_LOCAL_S9B_SHARED_VOCABULARY'])[2]],
        'knives':[{'knife':knife,'construction_site':{k:v for k,v in source_location(ENGINE.read_text(),pair[0]).items() if k!='text'},
                   'source':str((scratch/knife/'replica/scripts'/ENGINE.name).relative_to(REPOSITORY)),
                   'source_sha256':sha(scratch/knife/'replica/scripts'/ENGINE.name),
                   'witness_search':json.loads((scratch/knife/'witness-summary.json').read_text())} for knife,pair in KNIVES.items()],
        'item13':{'every_b':'Unchanged engine: Abel inversion to radial identities and velocity-square branches.',
            'component_calculus':'Unchanged engine: coordinate spherical Jacobian and divergence on d3x.',
            'profiles':'Unchanged engine: general radial functions, derivatives, integrals and quantified ray sets.'},
        'bindings':[{'fold_key':'c_s0','spec_line':223,'source':'S9b_SHARED_PHYSICS.md','construction':'computed asymptotic bulk sound-speed square'}],
        'shared_profile_families':[],
        'K11':k11,'guard_refusals_or_kills_this_round':[],'item12_events_this_round':[],
        'f9_collisions_baseline_count':len(decode(baseline['PY_LOCAL_S9B_EXPORT_ROUTES'])),
        'open_inputs':['Phi versus GM','unsupplied material relations'],
        'engine_unchanged':sha(ENGINE)==pins[str(ENGINE.relative_to(REPOSITORY))],
        'export_delta_unchanged':sha(ENGINE.with_name('S9b_exports.py'))==pins[str(ENGINE.with_name('S9b_exports.py').relative_to(REPOSITORY))]}
    inventory['shared_count']=len(inventory['shared_tags'])
    inventory['local_count']=len(inventory['local_tags'])
    INVENTORY.write_text(json.dumps(inventory,indent=2)+'\n')


def domain_literal_flags(left,right):
    """Flags on explicit AD carriers without rebuilding any expression."""
    def project(text):
        root=ast.parse(text,mode='eval').body
        lines=[line.encode() for line in text.splitlines(keepends=True)]
        def segment(node):
            if node.lineno!=node.end_lineno:return ast.get_source_segment(text,node)
            # All engine payloads are one line; slicing avoids reparsing the
            # entire large export payload for every literal and AD carrier.
            return lines[node.lineno-1][node.col_offset:node.end_col_offset].decode()
        domains=[]
        def visit(node,path):
            if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id=='Tuple' and len(node.args)==2:
                first=node.args[0]
                if (isinstance(first,ast.Call) and isinstance(first.func,ast.Name) and first.func.id=='Str'
                    and len(first.args)==1 and isinstance(first.args[0],ast.Constant) and first.args[0].value=='ATTACHED_DOMAIN'):
                    domains.append((path,literal_key(segment(node.args[1]))))
                    return ('DOMAIN_SLOT',)
            if isinstance(node,ast.Constant):return literal_key(segment(node))
            if isinstance(node,ast.AST):
                return (type(node).__name__,tuple((field,visit(value,path+(field,))) for field,value in ast.iter_fields(node)))
            if isinstance(node,list):return tuple(visit(value,path+(index,)) for index,value in enumerate(node))
            return node
        value=visit(root,())
        return value,domains
    a,ad=project(left);b,bd=project(right)
    return {'value_structure':'identical' if a==b else 'not compared',
            'explicit_domain_structure':'identical' if ad==bd else 'not compared',
            'explicit_domain_carriers':(len(ad),len(bd))}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scratch',type=Path,required=True)
    parser.add_argument('--checks-only',action='store_true')
    parser.add_argument('--guard-log-directory',type=Path)
    args=parser.parse_args()
    scratch=args.scratch.resolve()
    scratch.relative_to(REPOSITORY/'_scratch')
    scratch.mkdir(parents=True,exist_ok=False)
    source=ENGINE.read_text()
    export=ENGINE.with_name('S9b_exports.py')
    inputs=(HERE,ENGINE,export,ROOT/'directives/S9b_SHARED_PHYSICS.md',ROOT/'directives/S9b_repair_build_directive.md',
            *(ROOT/'scripts'/name for name in ('ledger_fold.py','S11c_b_exports.py','S11c_c1_exports.py','S11c_c2_exports.py')))
    pins={str(path.relative_to(REPOSITORY)):sha(path) for path in inputs}
    (scratch/'source-pins.json').write_text(json.dumps(pins,indent=2)+'\n')
    snapshots=scratch/'sources'
    snapshots.mkdir()
    for path in inputs:(snapshots/path.name).write_bytes(path.read_bytes())
    self_checks()
    if args.checks_only:return
    emit('STRUCTURAL_COMPARISON',STRUCTURAL_METHOD)
    # Baseline source is canonical and unmodified. Its export is deterministic;
    # verify bytes after completion rather than replacing engine publication.
    with work_part('baseline engine'):output=run(ENGINE,scratch/'baseline')
    baseline=read_stream(output)
    runs=[('baseline',output,ENGINE)]
    k11=None
    for knife,(before,after) in KNIVES.items():
        count=source.count(before)
        emit(knife+'_CONSTRUCTION',(before,after,count,source_location(source,before)))
        if count!=1:raise ValueError(('construction-site count',knife,count))
        replica=scratch/knife/'replica'
        scripts=replica/'scripts'
        scripts.mkdir(parents=True)
        (replica/'directives').symlink_to(ROOT/'directives',target_is_directory=True)
        for name in ('ledger_fold.py','S11c_b_exports.py','S11c_c1_exports.py','S11c_c2_exports.py'):
            (scripts/name).symlink_to(ROOT/'scripts'/name)
        mutant=scripts/ENGINE.name
        mutant.write_text(source.replace(before,after,1))
        if knife=='K11':
            with work_part('K11 compact'):k11=compact_k11(baseline,source,mutant,scratch/knife)
            summary=k11['witness_summary']
        else:
            with work_part(knife+' engine'):worker_output=run(mutant,scratch/knife/'run')
            runs.append((knife,worker_output,mutant))
            corrupted=read_stream(worker_output)
            counts=Counter()
            flags=[]
            for tag in sorted(set(baseline)|set(corrupted)):
                progress(knife+' '+tag+' literal structure')
                a=baseline.get(tag,render(Str('MISSING_TAG')))
                b=corrupted.get(tag,render(Str('MISSING_TAG')))
                mark=structural_mark(a,b)
                counts[mark]+=1
                role_flags=domain_literal_flags(a,b)
                flags.append({'tag':tag,'mark':mark,**role_flags})
                print('PY_LOCAL_S9B_ABLATION_'+knife+'_'+tag+': Tuple('+a+', '+b+', '+render(Str(mark))+')',flush=True)
            (scratch/knife/'marks.json').write_text(json.dumps({'tag_count':len(flags),'counts':dict(counts),'tags':flags},indent=2)+'\n')
            if set(baseline)!=set(corrupted):raise ValueError(('tag inventory changed',knife))
            with work_part(knife+' finite witness search'):
                summary=finite_search(knife,baseline,corrupted,source,scratch/knife)
        (scratch/knife/'witness-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    # Nothing can change a pinned source while this run is being reported.
    changed=[path for path,digest in pins.items() if sha(REPOSITORY/path)!=digest]
    if changed:raise ValueError(('pinned input changed during run',changed))
    (scratch/'completed.json').write_text(json.dumps({'pins':pins,'runs':[(i,str(p),str(s)) for i,p,s in runs],
        'K11':k11,'command':sys.argv},indent=2)+'\n')
    if args.guard_log_directory:
        # All tags have flushed. Guard outcome metadata is added mechanically
        # after the guard exits; this step reads no scientific payload values.
        runs.append(('full_harness',args.guard_log_directory.resolve()/'stdout',HERE))
        make_inventory(scratch,runs,pins,k11,sys.argv)


if __name__=='__main__':
    main()
