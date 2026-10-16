#!/usr/bin/env python3
"""S9b amendment 11: literal marks for all tags; one witness search per knife.

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
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_location(source,text):
    if source.count(text)!=1:raise ValueError(('source-location ambiguity',text))
    return {'file':str(ENGINE.relative_to(REPOSITORY)),
            'line':source[:source.index(text)].count('\n')+1,'text':text}


# The witness plans select downstream emitted construction objects, not the
# input digest or a pre-mutation object. They introduce no engine mutation.
PLANS = {
    'K1':('PY_LOCAL_S9B_OPTICAL_Q','whole','arrj=evaluate_jet(arr,mapping).tidy()'),
    'K2':('PY_LOCAL_S9B_OPTICAL_Q','whole','arrj=evaluate_jet(arr,mapping).tidy()'),
    'K3':('PY_LOCAL_S9B_OPTICAL_Q','whole','arrj=evaluate_jet(arr,mapping).tidy()'),
    'K4a':('PY_LOCAL_S9B_OPTICAL_Q','whole','optical_q=((arrj**sp.Rational(1,2)*c0)/optical_r.diff(r)).tidy()'),
    'K4b':('PY_LOCAL_S9B_OPTICAL_Q','whole','optical_q=((arrj**sp.Rational(1,2)*c0)/optical_r.diff(r)).tidy()'),
    'K4c':('PY_LOCAL_S9B_OPTICAL_Q','whole','slope = radial_diff(xi,r)'),
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


def direct_search(knife,left_stream,right_stream,source):
    tag,selector,entry=PLANS[knife]
    a_whole,b_whole=decode(left_stream[tag]),decode(right_stream[tag])
    a,pa=select_part(a_whole,selector)
    b,pb=select_part(b_whole,selector)
    if a is None or b is None:
        return {'status':'none found','tag':tag,'tried':('selected grade is absent from a sparse emitted mapping',selector),
                'source_entry':source_location(source,entry)}
    conditions=[attached_conditions(a_whole,pa),attached_conditions(b_whole,pb)]
    parts=(a,b)
    profiles=dict(PROFILES)
    expressions=sp.Tuple(a,b,*(x for group in conditions for p,x in group))
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
    reading=input_reading(parts,conditions,profiles,point,fixings)
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
    return {'status':'witness' if witness else 'none found','site_part':'value','tag':tag,
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


def cautious_condition_search(knife,left_stream,right_stream,source):
    """Try a real profile/point, but retain unresolved quantified ray gates.

    The candidate cannot be a witness until every containing condition is
    decided on its whole range. No arbitrary ray turning radius, momentum,
    bound radius, or every-b radius is ever assigned as a point input.
    """
    tag,selector,entry=PLANS[knife]
    wholes=[decode(stream[tag]) for stream in (left_stream,right_stream)]
    parts_paths=[select_part(whole,selector) for whole in wholes]
    attempts=[]
    for number,whole in enumerate(wholes):
        part,path=parts_paths[number]
        conditions=attached_conditions(whole,path)
        if knife=='K10':
            # Outer gated payloads are exclusive explicit cases. This is the
            # first case's gate after the K10 result has entered branch_domain.
            # Its travel_domain conjunct is required even when existence is
            # negated. The value applicability carries speed identification.
            conditions.append(((0,1,1,1),whole[0][1][1][1]))
        # Do not call doit on the ray angular closure integrals with unknown
        # roots: they are a solved-construction obligation, not free inputs.
        candidate=part.replace(lambda x:isinstance(x,AppliedUndef) and x.func.__name__ in PROFILES,
            lambda x:PROFILES[x.func.__name__].subs(RADIUS,x.args[0]))
        free=whole.free_symbols
        ranged={symbol for symbol in free if knife=='K7' and symbol.name=='s9b_r'}
        assignments={symbol:sp.sympify(POINT[symbol.name]) for symbol in free-ranged if symbol.name in POINT}
        prepared=candidate.subs(assignments,simultaneous=True)
        # A literal equality decision is sufficient only for this candidate's
        # local equation. It says nothing about its still quantified gates.
        decided=sp.true if prepared is sp.true else sp.false if prepared is sp.false else Str('UNRESOLVED')
        attempts.append({'copy':number,'payload_location':path,'part_at_inputs':prepared,
            'point':assignments,'profiles':PROFILES,
            'ranged':[(symbol,sp.Interval.open(1,sp.oo)) for symbol in ranged],
            'bound_variables':binders(whole),
            'attached_conditions':conditions,
            'partial_evaluation':decided,
            'unresolved':'Exact joint turning-point/angular-closure solution sets and their intersection were not reduced by this witness search. Whole ray gates remain active; no sample radius substitutes a bound or ranged variable.'})
    return {'status':'none found','site_part':'domain' if knife=='K10' else 'value',
        'tag':tag,'source_entry':source_location(source,entry),
        'construction_site':source_location(source,KNIVES[knife][0]),
        'tried':attempts,'qualification':'No witness: unresolved is not false and structural difference is not a change claim.'}


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
    if result['xyz']!=other_result['xyz'] or result['grades']!=other_result['grades'] or templates!=other_templates:
        raise ValueError('K11 configuration mismatch')
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
    emit('K11_WITNESS',{'status':'witness' if witness else 'none found',
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
        'tried':('same compact inputs in both copies; first sample and matrix component',)})
    (folder/'marks.json').write_text(json.dumps({'tag_count':len(flags),'tags':flags},indent=2)+'\n')
    receipt={'evaluated_tags':reached,'not_evaluated_tags':not_evaluated,
             'source_sha256':hashlib.sha256(source.encode()).hexdigest(),
             'mutant_sha256':sha(mutant),'harness_sha256':sha(HERE),
             'witness_summary':{'status':'witness' if witness else 'none found','tag':tag,'site_part':'value'}}
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
    # Software checks, with no physics expected-value acceptance test.
    checks={
        'syntax_whitespace':structural_mark('Tuple(Integer(1), Integer(2))','Tuple( Integer(1),Integer(2) )')=='identical',
        'argument_order':structural_mark('Tuple(Integer(1), Integer(2))','Tuple(Integer(2), Integer(1))')=='not compared',
        'no_algebra':structural_mark('Add(Integer(1),Integer(1))','Integer(2)')=='not compared',
        'float_decimal_tail':structural_mark("Float('1.00000000000000001', precision=100)","Float('1.00000000000000002', precision=100)")=='not compared',
        'float_precision':structural_mark("Float('1.0', precision=53)","Float('1.0', precision=54)")=='not compared',
        'python_float_token':structural_mark('f(1.00000000000000001)','f(1.00000000000000002)')=='not compared',
        'binder_names_literal':structural_mark("Integral(Symbol('x'), (Symbol('x'),0,1))","Integral(Symbol('y'), (Symbol('y'),0,1))")=='not compared',
    }
    x=sp.Symbol('self_check_x',real=True)
    dummy=sp.Symbol('self_check_dummy',positive=True)
    function=sp.Function('s9b_f',real=True)
    outer=sp.Tuple(Str('LOCAL_BULK_DENSITY_PARTS'),sp.Tuple(Str('RHO_POSITIVE'),sp.Gt(x,0)),
        sp.Tuple(Str('RHO_NONPOSITIVE'),sp.Le(x,0)),sp.Tuple(Str('OBJECT_ON_BOTH_PARTS'),
            sp.Tuple(sp.Tuple(Str('ATTACHED_DOMAIN'),sp.Ge(x,0)),x)))
    checks['density_labels_are_not_conditions']=len(attached_conditions(outer,(3,1,1)))==1
    checks['derivative_before_point']=evaluate(sp.diff(function(x),x),PROFILES,{x:2})==-sp.Rational(1,40)
    checks['integral_bound_not_substituted']=evaluate(sp.Integral(function(x),(x,1,2)),PROFILES,{x:99})==sp.Rational(1,20)
    checks['empty_unknown_tuple_is_a_solution']=sp.Intersection(sp.FiniteSet(sp.Tuple()),sp.FiniteSet(sp.Tuple()))==sp.FiniteSet(sp.Tuple())
    checks['finite_tuple_witness']=finite_exact(sp.Tuple(1,sp.Tuple(2,3))) and different_exact(sp.Tuple(1),sp.Tuple(2))
    checks['unevaluated_not_finite']=not finite_exact(sp.Integral(function(x),(x,1,sp.oo)))
    source=ENGINE.read_text()
    # Fixtures exercise the complete witness output path, including the
    # singleton solution intersection and density partition carriers. Their
    # arbitrary functions are software inputs, not expected engine results.
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
    fixture=direct_search('K9a',{tag:render(response_fixture(1+fs))},
        {tag:render(response_fixture(sp.S.One))},source)
    checks['response_solution_intersection']=(fixture['status']=='witness'
        and fixture['solution_intersection'].is_empty is False
        and fixture['selected_solutions'][0]==fixture['selected_solutions'][1]
        and all(row[2] is sp.true for group in fixture['condition_evaluations'] for row in group))
    tag=PLANS['K1'][0]
    fixture=direct_search('K1',{tag:render(sp.Tuple(sp.Tuple(sp.Tuple(0,0,0),f)))},
        {tag:render(sp.Tuple(sp.Tuple(sp.Tuple(0,0,0),sp.S.Zero)))},source)
    checks['whole_numeric_mapping_reading']=fixture['status']=='witness'
    tag=PLANS['K8'][0]
    jn=sp.Function('j_n',real=True)(r)
    fixture=direct_search('K8',{tag:render(sp.Eq(jn,f,evaluate=False))},
        {tag:render(sp.Eq(jn,0,evaluate=False))},source)
    checks['decided_disjoint_solutions']=(fixture['status']=='witness'
        and fixture['solution_intersection'] is sp.EmptySet
        and all(solution.contains(value) is sp.true for solution,value in
                zip(fixture['copy_solution_sets'],fixture['selected_solutions'])))
    emit('SOFTWARE_CHECKS',checks)
    if not all(checks.values()):raise ValueError(('software checks',checks))


def raw_record(path,source_path,run_id):
    stream=read_stream(path)
    rows=[]
    marks=Counter()
    for number,(tag,payload) in enumerate(stream.items(),1):
        tree=ast.parse(payload,mode='eval')
        strings=[node.value for node in ast.walk(tree) if isinstance(node,ast.Constant) and isinstance(node.value,str)]
        markers=Counter(word for word in strings if word in ('NOT_ESTABLISHED','NOT_EVALUATED','ATTACHED_DOMAIN',
            'identical','not compared','witness','none found','UNRESOLVED'))
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
        'authority':'S9b repair build directive amendment 11, repair round 8; Parts A-C shared physics',
        'builder':'Codex gpt-6.1-sol at xhigh effort; fresh harness author',
        'inventory_format':'Names, provenance, markers, structural flags and counts only; no tag payload copied.',
        'structural_comparison':STRUCTURAL_METHOD,
        'flag_qualification':'identical is literal emitted-expression structure only; not compared makes no equality or change claim. Value/domain flags are literal marks, not mathematical differences. One witness search per knife.',
        'source_pins':pins,'runs':records,'command':command,
        'comparison_counts_by_knife':comparison_counts,
        'shared_tags':sorted(tag for tag in baseline if tag.startswith('PY_S9B_')),
        'local_tags':sorted(tag for tag in baseline if tag.startswith('PY_LOCAL_S9B_')),
        'not_established_tags':sorted(tag for tag,payload in baseline.items() if "Str('NOT_ESTABLISHED')" in payload),
        'missing_shared':[str(name) for name in decode(baseline['PY_LOCAL_S9B_SHARED_VOCABULARY'])[2]],
        'knives':[{'knife':knife,'construction_site':source_location(ENGINE.read_text(),pair[0]),
                   'source':str((scratch/knife/'replica/scripts'/ENGINE.name).relative_to(REPOSITORY)),
                   'source_sha256':sha(scratch/knife/'replica/scripts'/ENGINE.name),
                   'witness_search':json.loads((scratch/knife/'witness-summary.json').read_text())} for knife,pair in KNIVES.items()],
        'item13':{'every_b':'Unchanged engine: Abel inversion to radial identities and velocity-square branches.',
            'component_calculus':'Unchanged engine: coordinate spherical Jacobian and divergence on d3x.',
            'profiles':'Unchanged engine: general radial functions, derivatives, integrals and quantified ray sets.'},
        'bindings':[{'fold_key':'c_s0','spec_line':223,'source':'S9b_SHARED_PHYSICS.md','construction':'computed asymptotic bulk sound-speed square'}],
        'shared_profile_families':[],
        'K11':k11,'guard_refusals_or_kills_this_round':([json.loads((scratch.parent/'admission-refusal.json').read_text())]
            if (scratch.parent/'admission-refusal.json').exists() else []),'item12_events_this_round':[],
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
            with work_part(knife+' witness search'):
                search=(cautious_condition_search if knife in ('K7','K10') else direct_search)(knife,baseline,corrupted,source)
            emit(knife+'_WITNESS',search)
            summary={key:search[key] for key in ('status','tag','site_part') if key in search}
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
