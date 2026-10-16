#!/usr/bin/env python3
"""S9b single-site copy mutations from repair build directive item 14.

Run this harness through s11c_guarded_run.py --pool s9b --memory-gib 16.
Workers are serial children of that same cgroup. No worker deadline or retry.
Baseline is the unchanged canonical engine. Mutants write only in their own
scratch replica. Every baseline/mutant tag, including metadata, is differenced.
--baseline-output permits reuse of a completed unchanged baseline, provided its
source digest sidecar matches. It never silently replays a completed worker.
"""
import argparse
import ast
import hashlib
import os
from functools import lru_cache
from fractions import Fraction
from math import isqrt
from itertools import product
from pathlib import Path
import subprocess
import sys
import threading
import time
from contextlib import contextmanager
from collections import Counter
sys.dont_write_bytecode = True
import sympy as sp
from sympy.core.symbol import Str
from sympy.core.function import AppliedUndef
from sympy.printing.repr import ReprPrinter
from sympy.functions.elementary.piecewise import ExprCondPair

HERE = Path(__file__).resolve()
ENGINE = HERE.with_name('S9b_light_bending_sympy_audit.py')
ROOT = HERE.parent.parent
REPOSITORY = ROOT.parent.parent

# Exactly one replacement in exactly one construction per copy. The mutation
# manifest contains no expected payload or bite criterion.
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

RELATIONALS = {name: (lambda a,b, f=f: f(a,b,evaluate=False)) for name,f in
    (('Equality',sp.Eq),('Unequality',sp.Ne),('StrictGreaterThan',sp.Gt),
     ('StrictLessThan',sp.Lt),('GreaterThan',sp.Ge),('LessThan',sp.Le))}

FORWARD_CONDITION_TAGS = frozenset(
    ['PY_S9B_F_'+stage+'_B_CONDITION_'+quantity
     for stage in ('FLOW','LIVE') for quantity in ('DEFLECTION','RADAR')]+
    ['PY_S9B_F_LIVE_C_'+response+'_CONDITION_'+quantity
     for response in ('CONSTANT','FIXED_RATIO','POWER') for quantity in ('DEFLECTION','RADAR')])


_progress_lock=threading.Lock()


def progress(object_name):
    with _progress_lock:
        print('S9B_HARNESS_WORK '+str(object_name), file=sys.stderr, flush=True)


_work_stack=[]
_work_counts=Counter()
_work_detail='idle'
_work_candidate=''
_work_last=0.0


def work_detail(name):
    global _work_detail,_work_last
    _work_detail=name
    _work_counts[name.split(':',1)[0]]+=1
    now=time.monotonic()
    if now-_work_last>=1:
        progress('ADVANCE '+ '/'.join(_work_stack)+' '+_work_candidate+' '+name)
        _work_last=now


@contextmanager
def work_part(name):
    """Operational telemetry only; never enters a mathematical payload."""
    global _work_candidate,_work_detail
    parent_candidate=_work_candidate
    parent_detail=_work_detail
    _work_candidate=''
    _work_detail='enter '+name
    start=time.monotonic()
    _work_stack.append(name)
    progress('PART begin '+ '/'.join(_work_stack))
    try:yield
    finally:
        _work_counts[name]+=1
        _work_counts['parts_finished']+=1
        progress('PART end '+ '/'.join(_work_stack)+' seconds='+format(time.monotonic()-start,'.6f'))
        _work_stack.pop()
        _work_candidate=parent_candidate
        _work_detail=parent_detail


def start_progress_monitor():
    """Show current exact suboperation and completed work, without a deadline."""
    main_thread=threading.get_ident()
    stopped=threading.Event()
    def monitor():
        while not stopped.wait(30):
            frame=sys._current_frames().get(main_thread)
            stack=[]
            while frame is not None and len(stack)<8:
                stack.append(frame.f_code.co_name+':'+str(frame.f_lineno))
                frame=frame.f_back
            counts={name:_work_counts[name] for name in
                    ('parts_finished','parse','equality','witness','evaluate','membership','output')}
            progress('ACTIVE '+ '/'.join(_work_stack)+' '+_work_candidate+' '+_work_detail+' counters='+str(counts)+
                     ' stack='+str(stack))
    threading.Thread(target=monitor,daemon=True).start()
    return stopped

# Amendment 2: finite exact samples, not a universal path-independence claim.
# These dimensionless coordinate/unit choices are evaluation inputs only.
# In particular both copies receive identical nonzero profile amplitudes.
COMPACT_PARAMETERS = {'d':sp.Rational(1,10), 'v':sp.Rational(1,8),
                      'w':sp.Rational(1,6), 'p':sp.Integer(2),
                      'q':sp.Integer(1), 'h':sp.Integer(1),
                      'L':sp.Integer(1), 'c0':sp.Integer(1)}
COMPACT_AZIMUTHAL = sp.Rational(1,12)
COMPACT_POINTS = ((3,4,0),(4,0,3),(0,6,8))
COMPACT_TARGETS = {'PY_LOCAL_S9B_NONRECIPROCITY_ONEFORM':'oneform',
                   'PY_LOCAL_S9B_NONRECIPROCITY_EXTERIOR_DERIVATIVE':'curl',
                   'PY_LOCAL_S9B_NONRECIPROCITY_PATH_PREDICATE':'predicate'}


def compact_copy(path):
    """Call the evaluated copy's own vector-dispersion/one-form construction."""
    source=path.read_text()
    namespace={'__file__':str(path),'__name__':'s9b_compact_copy'}
    exec(compile(source,str(path),'exec'),namespace)
    r=sp.Symbol('s9b_r',positive=True)
    c0=sp.Symbol('c_0',positive=True)
    config=COMPACT_PARAMETERS
    templates={namespace['Delta']:config['d']*(config['L']/r)**config['p'],
               namespace['Velocity']:config['c0']*config['v']*(config['L']/r)**config['q'],
               namespace['Embedding']:config['w']*sp.log(r/config['L'])}
    # The xi template differentiates to the configured h=1 slope. It is a
    # compact-evaluation input only, never the engine's shared-profile ansatz.
    replacements={function:(lambda x,value=value:value.subs(r,x)) for function,value in templates.items()}
    with work_part('construct '+path.name):
        result=namespace['nonreciprocity'](r,config['c0'],templates[namespace['Delta']],
            templates[namespace['Velocity']],templates[namespace['Embedding']])
    parameters={c0:config['c0'],sp.Symbol('s9b_azimuthal',real=True):COMPACT_AZIMUTHAL}
    reached={}
    for tag,name in COMPACT_TARGETS.items():
        rows=[]
        for point in COMPACT_POINTS:
            with work_part('sample '+tag+' '+str(point)):
                rows.append(sp.simplify(namespace['cas'](result[name]).subs(parameters).subs(
                    dict(zip(result['xyz'],map(sp.Integer,point))))))
        reached[tag]=sp.Tuple(*rows)
    functions=[n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef)
               and n.name in ('vector_dispersion','nonreciprocity')]
    lines=tuple((node.lineno,node.end_lineno) for node in functions)
    reducer=(namespace['profile_substitute'],replacements,parameters)
    return reached,reducer,result['xyz'],result['grades'],lines


def compact_reduce(value, reducer, coordinates):
    apply,replacements,parameters=reducer
    reduced=apply(value,replacements).subs(parameters)
    return sp.Tuple(*(sp.simplify(reduced.subs(dict(zip(coordinates,map(sp.Integer,point)))))
                      for point in COMPACT_POINTS))


def _compact_difference(a,b,context=sp.true):
    # The same equality-first, exact-witness rule applies to compact values
    # and method residuals. Sampling coverage is supplied by compact_k11.
    return _value_difference(a,b,context)


def compact_k11(baseline, source, mutant, folder):
    unmutated = mutant.with_name('S9b_compact_unmutated.py')
    unmutated.write_text(source)
    progress('K11 compact unmutated dispersion/nonreciprocity')
    original,parameters,coordinates,grades,lines = compact_copy(unmutated)
    progress('K11 compact mutated dispersion/nonreciprocity')
    corrupted,other_parameters,other_coordinates,other_grades,other_lines = compact_copy(mutant)
    coverage = sp.Tuple(
        sp.Tuple(Str('method'),Str('engine vector_dispersion and nonreciprocity calls; exact profile specialization; derivatives before coordinate sampling')),
        sp.Tuple(Str('profile_parameters'),sp.Tuple(*(sp.Tuple(Str(k),v) for k,v in COMPACT_PARAMETERS.items()))),
        sp.Tuple(Str('azimuthal_parameter'),COMPACT_AZIMUTHAL),
        sp.Tuple(Str('sample_points'),sp.Tuple(*(sp.Tuple(*point) for point in COMPACT_POINTS))),
        sp.Tuple(Str('radial_retained_grades'),sp.Tuple(*(sp.Tuple(*grade) for grade in grades))),
        sp.Tuple(Str('additional_truncation'),Str('NONE beyond the engine retained grades; all vector velocity components enter the same graded dispersion')),
        sp.Tuple(Str('arithmetic'),Str('EXACT; no seed, rounding, or numerical error')),
        sp.Tuple(Str('boolean_difference'),Str('corrupted minus unmutated truth indicators; pointwise predicate only')),
        sp.Tuple(Str('extracted_line_ranges'),sp.Tuple(*(sp.Tuple(*pair) for pair in lines))))
    print('PY_LOCAL_S9B_ABLATION_K11_CONFIGURATION: '+render(coverage),flush=True)
    if coordinates != other_coordinates or grades != other_grades or lines != other_lines:
        raise ValueError('compact extraction/configuration mismatch')
    reached,not_evaluated = [],[]
    for tag,full in sorted(baseline.items()):
        progress('K11 '+tag+': compact evaluation, difference and method residual')
        if tag in original and tag in corrupted:
            with work_part(tag+' reduce baseline'):
                reduced = compact_reduce(full,parameters,coordinates)
            with work_part(tag+' difference'):
                delta = compact_difference(original[tag],corrupted[tag])
            with work_part(tag+' method residual'):
                residual = compact_difference(reduced,original[tag])
            payload = sp.Tuple(full,original[tag],corrupted[tag],
                sp.Tuple(Str('DIFFERENCE'),delta,coverage),
                sp.Tuple(Str('METHOD_RESIDUAL'),residual,coverage),
                sp.Tuple(Str('REDUCED_FULL_BASELINE'),reduced))
            reached.append(tag)
        else:
            missing = Str('NOT_EVALUATED')
            payload = sp.Tuple(full,missing,missing,
                sp.Tuple(Str('DIFFERENCE'),missing,coverage),
                sp.Tuple(Str('METHOD_RESIDUAL'),missing,coverage))
            not_evaluated.append(tag)
        print('PY_LOCAL_S9B_ABLATION_K11_'+tag+': '+render(payload),flush=True)
    # Receipt is run state, not an engine output or an acceptance verdict.
    import json
    (folder/'compact-receipt.json').write_text(json.dumps({
        'source_sha256':hashlib.sha256(source.encode()).hexdigest(),
        'mutant_sha256':hashlib.sha256(mutant.read_bytes()).hexdigest(),
        'harness_sha256':hashlib.sha256(HERE.read_bytes()).hexdigest(),
        'evaluated_tags':reached,'not_evaluated_tags':not_evaluated,
        'configuration':render(coverage)},indent=2)+'\n')


class PayloadPrinter(ReprPrinter):
    def _print_FunctionClass(self, expr):
        if issubclass(expr, AppliedUndef):
            return 'Function(%r, **%r)' % (expr.__name__, dict(sorted(expr._kwargs.items())))
        return super()._print_FunctionClass(expr)

    def _print_Intersection(self, expr):
        return 'Intersection('+', '.join(self._print(a) for a in expr.args)+', evaluate=False)'

    def _print_Union(self, expr):
        return 'Union('+', '.join(self._print(a) for a in expr.args)+', evaluate=False)'

    def _print_SymmetricDifference(self, expr):
        return 'SymmetricDifference('+', '.join(self._print(a) for a in expr.args)+', evaluate=False)'


def render(value):
    return PayloadPrinter().doprint(value)


def parse(path):
    namespace = {'__builtins__': {}, **vars(sp), 'Str':Str,
                 'ExprCondPair':ExprCondPair, **RELATIONALS}
    result = {}
    for line in path.read_text().splitlines():
        tag, sep, payload = line.partition(': ')
        if not sep or not (tag.startswith('PY_S9B_') or tag.startswith('PY_LOCAL_S9B_')):
            raise ValueError(('untagged line',line[:100]))
        if tag in result:
            raise ValueError(('duplicate tag',tag))
        progress("parse "+str(path)+" "+tag)
        work_detail('parse:'+tag)
        result[tag] = eval(payload,namespace)
    if not result:
        raise ValueError(('empty stream',str(path)))
    missing = FORWARD_CONDITION_TAGS-set(result)
    if missing:
        raise ValueError(('missing amended forward conditions', sorted(missing)))
    return result


def _primitive_integral(integral):
    """Extract only factors independent of every integration variable."""
    kernel=sp.factor_terms(integral.function)
    if kernel.is_zero is True:
        return sp.S.Zero
    coefficient,kernel=kernel.as_independent(*integral.variables,as_Add=False)
    return coefficient*sp.Integral(kernel,*integral.limits)


def _exact_relation(head,left,right):
    """Do not let an undecided exact functional trigger numerical quadrature.

    SymPy's eager Eq constructor can call evalf on a definite Integral.
    Retain that exact relation after algebraic cancellation, rather than use
    its numerical heuristic as an exact Boolean witness.
    """
    if left==right and head in (sp.Eq,sp.Ne):
        return sp.true if head is sp.Eq else sp.false
    if isinstance(left,sp.Tuple) and isinstance(right,sp.Tuple) and head in (sp.Eq,sp.Ne):
        equal=(sp.false if len(left)!=len(right) else
               sp.And(*(_exact_relation(sp.Eq,a,b) for a,b in zip(left,right))))
        return equal if head is sp.Eq else sp.Not(equal)
    inert=(sp.Integral,sp.Derivative,sp.Limit,sp.Subs)
    if any(isinstance(x,sp.Basic) and x.has(*inert) for x in (left,right)):
        if isinstance(left,sp.Expr) and isinstance(right,sp.Expr):
            delta=_arithmetic_reduce(left-right)
            if delta==0:return head(sp.S.Zero,sp.S.Zero)
        return head(left,right,evaluate=False)
    if all(isinstance(x,sp.Rational) or x in (sp.oo,-sp.oo) for x in (left,right)):
        return head(left,right)
    return head(left,right,evaluate=False)


def undecided(a,b,residual):
    return sp.Tuple(Str('UNDECIDED'),sp.Tuple(Str('BASELINE'),a),
                    sp.Tuple(Str('CORRUPTED'),b),sp.Tuple(Str('RESIDUAL'),residual))


def _arithmetic_reduce(value):
    # Treat live functionals as independent atoms during algebra. This
    # changes no integration method and never samples a live profile.
    value=value.replace(lambda node:isinstance(node,sp.Integral),_primitive_integral)
    functionals=value.atoms(sp.Integral,sp.Derivative,sp.Limit,sp.Subs)
    opaque={node:sp.Dummy('s9b_comparison_functional') for node in sorted(functionals,key=sp.default_sort_key)}
    reduced=sp.factor_terms(sp.expand(value.xreplace(opaque), power_base=False, power_exp=False, log=False))
    return reduced.xreplace({v:k for k,v in opaque.items()})


def _known_zero(value):
    functionals=value.atoms(sp.Integral,sp.Derivative,sp.Limit,sp.Subs)
    opaque={x:sp.Dummy('s9b_exact_functional') for x in functionals}
    return value.xreplace(opaque).is_zero is True


@lru_cache(maxsize=256)
def _normalized_predicate(value):
    """Canonicalize conditions even inside sets and quantified predicates.

    Only identities established by CAS algebra are used. In particular, an
    unknown multiplicative factor is never divided out of a zero equation.
    """
    def visit(node):
        if not isinstance(node,sp.Basic) or not node.args:
            return node
        if isinstance(node,(sp.ConditionSet,sp.Integral,sp.Limit,sp.Lambda,sp.Subs)):
            node=node.as_dummy()
        if isinstance(node,sp.Integral) and node.function.is_zero is True:
            return sp.S.Zero
        args=tuple(visit(arg) for arg in node.args)
        if isinstance(node,sp.core.relational.Relational):
            left,right=args
            # Held set equalities can become reflexive after their bound
            # conditions reduce. Let the CAS discharge that identity first.
            if isinstance(node,(sp.Equality,sp.Unequality)) and left==right:
                return node.func(left,right)
            if all(isinstance(x,sp.Expr) for x in args):
                delta=_arithmetic_reduce(left-right)
                factors=sp.Mul.make_args(delta)
                divisor=sp.S.One
                for factor in factors:
                    if factor.has(sp.Integral,sp.Derivative,sp.Limit,sp.Subs):continue
                    if isinstance(node,(sp.Equality,sp.Unequality)):
                        if factor.is_zero is False:divisor*=factor
                    elif factor.is_positive is True:
                        divisor*=factor
                if divisor!=1:delta=sp.cancel(delta/divisor)
                return _exact_relation(node.func,delta,sp.S.Zero)
            return node.func(left,right,evaluate=False)
        if isinstance(node,(sp.Intersection,sp.Union,sp.SymmetricDifference)):
            return node.func(*args,evaluate=False)
        return node.func(*args) if args!=node.args else node
    return visit(value)


# Witness profiles are test inputs only. Each entry is a *specific* function
# with fixed rational coefficients, not independent values of radial jets.
# The same function is used at every argument and differentiated before the
# witness point is substituted. No order-counting or matching condition is
# imposed on these tests (amendment 10).
WITNESS_RADIAL = sp.Symbol('s9b_witness_radius', positive=True)
WITNESS_LITERAL_INFINITY = sp.Dummy('witness_literal_infinity',positive=True)
WITNESS_VALUES = (sp.S.Zero, sp.S.One, -sp.S.One, sp.Integer(2),
                  sp.Rational(1,2), sp.Integer(-2))


def _profile_bank(index, functions):
    r = WITNESS_RADIAL
    ordinary = (sp.S.Zero, sp.S.One, r, r*r+1, 1/r, -1/r)
    positive = (sp.S.One, sp.Integer(2), r, r*r+1, 1/r, 1+r)
    return {f:sp.Lambda(r, (positive if f(r).is_positive else ordinary)[(index+offset)%6])
            for offset,f in enumerate(functions)}


def _point(symbols, index):
    result = {}
    for offset, symbol in enumerate(symbols):
        rotation=(index+offset)%len(WITNESS_VALUES)
        candidates = WITNESS_VALUES[rotation:]+WITNESS_VALUES[:rotation]
        for candidate in candidates:
            # Membership in the symbol's declared assumptions is checked;
            # no value is assigned to a bound variable.
            if all(getattr(candidate, 'is_'+key) == value
                   for key,value in symbol.assumptions0.items()
                   if key in ('positive','negative','nonnegative','nonpositive',
                              'nonzero','integer','real') and value is True):
                result[symbol] = candidate
                break
        else:
            return None
    return result


def _exact_sets(node):
    """Evaluate statements about sets, retaining their binders.

    A sample element can disprove emptiness, but cannot establish emptiness.
    Here only exact set reduction can decide an emptiness statement. A set
    *value*, in contrast, is tested by contains(element) by the caller.
    """
    if isinstance(node, sp.ConditionSet):
        variable, condition, base = node.args
        if (isinstance(variable,sp.Symbol) and base.is_subset(sp.S.Reals) is True
                and condition.free_symbols <= {variable}):
            if condition is sp.true:return base
            if condition is sp.false:return sp.EmptySet
            if isinstance(condition,sp.core.relational.Relational):
                try:
                    arithmetic=_IntervalArithmetic()
                    numerator,denominator=arithmetic.rational_function(condition.lhs-condition.rhs,variable)
                    if (max(len(numerator),len(denominator))>3 or
                        any(max(c.numerator.bit_length(),c.denominator.bit_length())>128
                            for c in numerator+denominator)):
                        return node
                except _EvaluationUnresolved:
                    return node
                # Only rational linear/quadratic polynomials reach the exact
                # root/inequality reducer. Higher-degree conditions remain
                # explicit and their element memberships are still tested.
                numerator_expr=sum(_fraction_sympy(c)*variable**i for i,c in enumerate(numerator))
                denominator_expr=sum(_fraction_sympy(c)*variable**i for i,c in enumerate(denominator))
                if denominator_expr==0:return node
                delta=numerator_expr/denominator_expr
                poles=sp.solveset(denominator_expr,variable,domain=base)
                if isinstance(condition,sp.Equality):
                    return sp.Complement(sp.solveset(numerator_expr,variable,domain=base),poles)
                reduced=sp.reduce_inequalities([condition.func(delta,0)],variable)
                return sp.Complement(sp.Intersection(base,reduced.as_set()),poles)
        return node
    if isinstance(node,(sp.Intersection,sp.Union,sp.Complement,sp.SymmetricDifference)):
        return node.func(*node.args)
    if isinstance(node,sp.core.relational.Relational):
        left,right=node.args
        if isinstance(left,sp.Set) and isinstance(right,sp.Set):
            if left==right:return sp.true if isinstance(node,sp.Equality) else sp.false
            # Never infer unequal merely because a held set has a new tree.
            if not left.has(sp.ConditionSet) and not right.has(sp.ConditionSet):
                if left is sp.EmptySet and right.is_empty is False:
                    return sp.false if isinstance(node,sp.Equality) else sp.true
                if right is sp.EmptySet and left.is_empty is False:
                    return sp.false if isinstance(node,sp.Equality) else sp.true
                return _exact_relation(node.func,left,right)
            return node
        return _exact_relation(node.func,left,right)
    return node


def _has_function_class(value):
    if isinstance(value,sp.FunctionClass):return True
    return isinstance(value,sp.Basic) and any(_has_function_class(x) for x in value.args)


WITNESS_NATIVE_TREE_NODES=512
WITNESS_NATIVE_SET_NODES=32


@lru_cache(maxsize=4096)
def _native_evaluation_supported(value):
    pending=[value];nodes=sets=0
    while pending:
        node=pending.pop();nodes+=1
        if nodes>WITNESS_NATIVE_TREE_NODES:return False
        if isinstance(node,sp.Set):
            sets+=1
            if sets>WITNESS_NATIVE_SET_NODES:return False
        if isinstance(node,sp.Basic):pending.extend(node.args)
    return True


def _evaluate_at(value, profiles_tuple, point_tuple):
    # Parameters absent from this operand cannot affect its evaluation. Drop
    # them from the memoization key, not from the witness point or coverage.
    # Every emitted witness still carries the full comparison's assignments.
    if not _native_evaluation_supported(value):
        return Str('EVALUATION_UNRESOLVED_NATIVE_TREE_BUDGET')
    if isinstance(value,sp.Basic):
        functions={node.func for node in value.atoms(AppliedUndef)}
        profiles_tuple=tuple((f,p) for f,p in profiles_tuple if f in functions)
        point_tuple=tuple((x,v) for x,v in point_tuple if x in value.free_symbols)
    return _evaluate_at_cached(value,profiles_tuple,point_tuple)


@lru_cache(maxsize=4096)
def _evaluate_at_cached(value, profiles_tuple, point_tuple):
    # Export inventories contain bare FunctionClass entries, not radial
    # function values. The evaluator must not call Basic.replace on classes.
    # This candidate is unresolved; no structural inequality is called change.
    if _has_function_class(value):return Str('UNEVALUATED_FUNCTION_CLASS')
    try:
        return _evaluate_at_raw(value,profiles_tuple,point_tuple)
    except TypeError as error:
        if not any(message in str(error) for message in ('Invalid NaN comparison','Invalid comparison of non-real')):raise
        return sp.nan  # Undefined candidate, not a witness.
    except ZeroDivisionError:
        return sp.nan


class _WitnessPredicate(sp.logic.boolalg.BooleanFunction):
    """Hold predicate evaluation until substituted operands are defined."""


@lru_cache(maxsize=128)
def _hold_predicates(value):
    return value.replace(lambda x:isinstance(x,sp.logic.boolalg.Boolean)
        and not isinstance(x,sp.Expr) and x not in (sp.true,sp.false)
        and not isinstance(x,_WitnessPredicate),
        lambda x:_WitnessPredicate(Str(x.func.__name__),*x.args))


def _restore_predicate(node):
    if node.has(sp.nan,sp.zoo):return sp.nan
    head=getattr(sp,str(node.args[0]))
    if head is sp.Contains:
        return _membership(node.args[2],node.args[1])
    if issubclass(head,sp.core.relational.Relational):
        return _exact_relation(head,*node.args[1:])
    return head(*node.args[1:])


@lru_cache(maxsize=256)
def _membership(value,element):
    # Hold the condition before ConditionSet.contains substitutes its bound
    # variables. Otherwise that substitution can eagerly evaluate an Eq and
    # launch numerical integration before the exact evaluator sees it.
    def formula(domain,item):
        if isinstance(domain,sp.FiniteSet):
            return sp.Or(*(_exact_relation(sp.Eq,item,x) for x in domain))
        if isinstance(domain,(sp.Union,sp.Intersection,sp.SymmetricDifference)):
            operation={sp.Union:sp.Or,sp.Intersection:sp.And,sp.SymmetricDifference:sp.Xor}[type(domain)]
            return operation(*(formula(x,item) for x in domain.args))
        if isinstance(domain,sp.Complement):
            return sp.And(formula(domain.args[0],item),sp.Not(formula(domain.args[1],item)))
        if isinstance(domain,sp.ProductSet):
            if not isinstance(item,sp.Tuple) or len(item)!=len(domain.args):return sp.false
            return sp.And(*(formula(x,y) for x,y in zip(domain.args,item)))
        if isinstance(domain,sp.ConditionSet):
            def same_shape(a,b):
                if isinstance(a,sp.Tuple)!=isinstance(b,sp.Tuple):return False
                return (not isinstance(a,sp.Tuple) or
                        (len(a)==len(b) and all(same_shape(x,y) for x,y in zip(a,b))))
            if not same_shape(domain.sym,item):return sp.false
            base=formula(domain.base_set,item)
            if base is sp.false:return sp.false
            # This is ConditionSet's own Lambda substitution, with its
            # predicates held and its base membership evaluated exactly.
            condition=sp.Lambda((domain.sym,),_hold_predicates(domain.condition))(item)
            return sp.And(base,condition)
        if isinstance(domain,sp.Interval):
            if isinstance(item,sp.Tuple):return sp.false
            lower=sp.Gt if domain.left_open else sp.Ge
            upper=sp.Lt if domain.right_open else sp.Le
            return sp.And(_exact_relation(lower,item,domain.start),
                          _exact_relation(upper,item,domain.end))
        if domain.has(sp.Integral,sp.Derivative,sp.Limit,sp.Subs):
            return sp.Contains(item,domain,evaluate=False)
        return domain.contains(item)
    condition=formula(value,element)
    if condition.has(sp.nan,sp.zoo):return sp.nan
    condition=condition.replace(lambda x:isinstance(x,_WitnessPredicate),_restore_predicate)
    return _exact_sets(condition)


@lru_cache(maxsize=1024)
def _prepare_profiles(value,profiles_tuple):
    """Reuse exact profile substitution/derivatives across coordinate points."""
    profiles=dict(profiles_tuple)
    value=_hold_predicates(value)
    value=value.xreplace({sp.oo:WITNESS_LITERAL_INFINITY,-sp.oo:-WITNESS_LITERAL_INFINITY})
    # Eq(zoo, 0) would collapse to False and hide an undefined operand.
    # Hold all predicate/Boolean heads through substitution and integration;
    # reject a pole before restoring their actual mathematical operations.
    work_detail('evaluate:profiles')
    def specialize(node):
        if isinstance(node,AppliedUndef) and node.func in profiles:
            return profiles[node.func](*node.args)
        if isinstance(node,sp.Derivative) and not node.has(AppliedUndef,sp.Integral,sp.Limit,sp.Subs):
            return node.doit()
        if isinstance(node,sp.Subs):return _resolve_subs(node)
        return node
    # Resolve each specialized derivative and Subs before rebuilding its
    # enclosing arithmetic. This preserves native simultaneous substitution
    # and prevents unevaluated Subs nodes from breaking Basic.compare.
    value=value.replace(lambda x:isinstance(x,(AppliedUndef,sp.Derivative,sp.Subs)),specialize)
    return value


def _resolve_subs(node):
    if node.expr.has(sp.Derivative):return node
    return node.expr.subs(dict(zip(node.variables,node.point)),simultaneous=True)


def _rational_integral(value):
    # Never launch general symbolic integration in a witness. The original
    # exact integral remains the operand of the finite certified evaluator.
    return sp.S.Zero if value.function==0 else value


def _evaluate_at_raw(value, profiles_tuple, point_tuple):
    if not isinstance(value,sp.Basic):return value
    value=_prepare_profiles(value,profiles_tuple)
    work_detail('evaluate:point')
    value=value.subs(dict(point_tuple), simultaneous=True)
    value=value.replace(lambda x:isinstance(x,sp.Subs),lambda x:_resolve_subs(x))
    value=value.replace(lambda x:isinstance(x,sp.Integral) and x.function==0,
                        lambda x:sp.S.Zero)
    if value.has(sp.nan,sp.zoo,sp.oo,-sp.oo):return sp.nan
    value=value.xreplace({WITNESS_LITERAL_INFINITY:sp.oo})
    # Retain exact integral operands for certified finite-work evaluation.
    # Unsupported or exhausted evaluations remain explicitly unresolved.
    work_detail('evaluate:retain_integrals')
    value=value.replace(lambda x:isinstance(x,sp.Integral),_rational_integral)
    if value.has(sp.nan,sp.zoo):return sp.nan
    work_detail('evaluate:restore_predicates')
    value=value.replace(lambda x:isinstance(x,_WitnessPredicate),_restore_predicate)
    work_detail('evaluate:exact_sets')
    value=value.replace(lambda x:isinstance(x,(sp.ConditionSet,sp.Intersection,
        sp.Union,sp.Complement,sp.SymmetricDifference,sp.core.relational.Relational)),_exact_sets)
    return value


def _decided_value(value):
    if isinstance(value,sp.Tuple):return all(_decided_value(x) for x in value)
    if isinstance(value,Str):return value!=Str('UNEVALUATED_FUNCTION_CLASS')
    if value in (sp.true,sp.false):return True
    if not isinstance(value,sp.Basic):return False
    return (not value.free_symbols and not value.has(AppliedUndef,sp.Integral,
            sp.Derivative,sp.Limit,sp.Subs,sp.Set,sp.Float,sp.nan,sp.zoo,sp.oo,-sp.oo)
            and isinstance(value,sp.Expr) and value.is_finite is True)


def _element_dimensions(value):
    """The element of a set with tuple binders is an ordered tuple."""
    if isinstance(value,sp.ConditionSet):
        return {len(value.sym)} if isinstance(value.sym,sp.Tuple) else _element_dimensions(value.base_set)
    if isinstance(value,sp.ProductSet):return {len(value.args)}
    if isinstance(value,sp.FiniteSet):
        return {len(x) if isinstance(x,sp.Tuple) else 1 for x in value} or {1}
    if isinstance(value,(sp.Union,sp.Intersection,sp.Complement,sp.SymmetricDifference)):
        return set().union(*(_element_dimensions(x) for x in value.args))
    return {1}


@lru_cache(maxsize=256)
def _set_elements(left,right):
    candidates=set()
    for dimension in _element_dimensions(left)|_element_dimensions(right):
        if dimension==1:candidates.update(WITNESS_VALUES)
        else:candidates.update(sp.Tuple(*row) for row in product(WITNESS_VALUES,repeat=dimension))
    for side in (left,right):
        if isinstance(side,sp.FiniteSet):candidates.update(side)
        if isinstance(side,sp.Interval):
            candidates.update((side.start,side.end,(side.start+side.end)/2))
    # An unbounded interval's midpoint can be NaN. Retain that candidate
    # for explicit definedness reporting, without ordering NaN numerically.
    return sorted(candidates,key=lambda x:(x.has(sp.nan),
        sp.srepr(x) if x.has(sp.nan) else sp.default_sort_key(x)))


# Finite-work witness arithmetic. These are evaluation limits, not job clocks.
# Exhausting a limit preserves UNDECIDED; it never establishes equality.
WITNESS_PRECISION_BITS=80
WITNESS_QUADRATURE_PANELS=128
WITNESS_ARITHMETIC_STEPS=50000
WITNESS_POLYNOMIAL_DEGREE=32


class _EvaluationUnresolved(Exception):
    pass


class _IntervalArithmetic:
    def __init__(self):self.steps=0
    def tick(self,amount=1):
        self.steps+=amount
        if self.steps>WITNESS_ARITHMETIC_STEPS:
            raise _EvaluationUnresolved('ARITHMETIC_STEP_BUDGET')
    def pair(self,lo,hi=None):
        self.tick()
        if hi is None:hi=lo
        lo,hi=Fraction(lo),Fraction(hi)
        if lo>hi:raise _EvaluationUnresolved('INVALID_INTERVAL')
        if lo==hi and max(lo.numerator.bit_length(),lo.denominator.bit_length())<=4096:
            return lo,hi
        scale=1<<WITNESS_PRECISION_BITS
        return (Fraction((lo*scale).__floor__(),scale),
                Fraction((hi*scale).__ceil__(),scale))
    def add(self,a,b):return self.pair(a[0]+b[0],a[1]+b[1])
    def neg(self,a):return -a[1],-a[0]
    def mul(self,a,b):
        values=[x*y for x in a for y in b]
        return self.pair(min(values),max(values))
    def div(self,a,b):
        if b[0]<=0<=b[1]:raise _EvaluationUnresolved('DENOMINATOR_NOT_SEPARATED_FROM_ZERO')
        return self.mul(a,self.pair(1/b[1],1/b[0]))
    def power(self,a,n):
        if abs(n)>64:raise _EvaluationUnresolved('POWER_BUDGET')
        if n<0:return self.div(self.pair(1),self.power(a,-n))
        result=self.pair(1)
        while n:
            if n&1:result=self.mul(result,a)
            n//=2
            if n:a=self.mul(a,a)
        return result
    def sqrt(self,a):
        if a[0]<0:raise _EvaluationUnresolved('SQRT_DOMAIN_NOT_ESTABLISHED')
        scale=1<<WITNESS_PRECISION_BITS
        def lower(q):return Fraction(isqrt((q.numerator*scale*scale)//q.denominator),scale)
        lo=lower(a[0]);hi=lower(a[1])
        if hi*hi<a[1]:hi+=Fraction(1,scale)
        return self.pair(lo,hi)
    def constant_pi(self):
        # Machin's formula; alternating-series remainder for each atan.
        def atan_inverse(q):
            total=Fraction(0)
            for k in range(64):
                self.tick();total+=Fraction((-1)**k,(2*k+1)*q**(2*k+1))
            error=Fraction(1,129*q**129)
            return self.pair(total-error,total+error)
        return self.add(self.mul(self.pair(16),atan_inverse(5)),
                        self.mul(self.pair(-4),atan_inverse(239)))
    def logarithm(self,a):
        if a[0]<=0:raise _EvaluationUnresolved('LOG_DOMAIN_NOT_ESTABLISHED')
        def point(q):
            k=0
            while q>2:
                self.tick();q/=2;k+=1
            while q<1:
                self.tick();q*=2;k-=1
            def series(t):
                total=Fraction(0)
                for j in range(64):
                    self.tick();total+=2*t**(2*j+1)/(2*j+1)
                tail=2*abs(t)**129/(129*(1-t*t))
                return self.pair(total-tail,total+tail)
            return self.add(series((q-1)/(q+1)),self.mul(self.pair(k),series(Fraction(1,3))))
        lo,hi=point(a[0]),point(a[1])
        return self.pair(lo[0],hi[1])
    def numeric(self,node):
        self.tick()
        if isinstance(node,sp.Rational):return self.pair(Fraction(int(node.p),int(node.q)))
        if node is sp.pi:return self.constant_pi()
        if isinstance(node,sp.Add):
            result=self.pair(0)
            for arg in node.args:result=self.add(result,self.numeric(arg))
            return result
        if isinstance(node,sp.Mul):
            result=self.pair(1)
            for arg in node.args:result=self.mul(result,self.numeric(arg))
            return result
        if isinstance(node,sp.Pow):
            base,exponent=node.args
            if isinstance(exponent,sp.Integer):return self.power(self.numeric(base),int(exponent))
            if exponent==sp.Rational(1,2):return self.sqrt(self.numeric(base))
            if exponent==sp.Rational(-1,2):return self.div(self.pair(1),self.sqrt(self.numeric(base)))
        if node.func is sp.Abs:
            a=self.numeric(node.args[0])
            return self.pair(0 if a[0]<=0<=a[1] else min(abs(x) for x in a),max(abs(x) for x in a))
        if node.func is sp.log:return self.logarithm(self.numeric(node.args[0]))
        if isinstance(node,sp.Integral):return self.integral(node)
        raise _EvaluationUnresolved('UNSUPPORTED_OR_UNASSIGNED_'+node.func.__name__)
    def polynomial(self,coefficients):
        while len(coefficients)>1 and not coefficients[-1]:coefficients.pop()
        if len(coefficients)>WITNESS_POLYNOMIAL_DEGREE+1:
            raise _EvaluationUnresolved('POLYNOMIAL_DEGREE_BUDGET')
        if any(max(x.numerator.bit_length(),x.denominator.bit_length())>4096 for x in coefficients):
            raise _EvaluationUnresolved('COEFFICIENT_BIT_BUDGET')
        return coefficients
    def padd(self,a,b):
        self.tick(max(len(a),len(b)))
        return self.polynomial([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0)
                                for i in range(max(len(a),len(b)))])
    def pmul(self,a,b):
        if len(a)+len(b)-2>WITNESS_POLYNOMIAL_DEGREE:
            raise _EvaluationUnresolved('POLYNOMIAL_DEGREE_BUDGET')
        self.tick(len(a)*len(b));out=[Fraction(0)]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):out[i+j]+=x*y
        return self.polynomial(out)
    def rational_function(self,node,variable):
        self.tick()
        one=[Fraction(1)]
        if node==variable:return [Fraction(0),Fraction(1)],one
        if isinstance(node,sp.Rational):return [Fraction(int(node.p),int(node.q))],one
        if isinstance(node,sp.Add):
            n,d=[Fraction(0)],one
            for arg in node.args:
                an,ad=self.rational_function(arg,variable)
                n,d=self.padd(self.pmul(n,ad),self.pmul(an,d)),self.pmul(d,ad)
            return n,d
        if isinstance(node,sp.Mul):
            n,d=one,one
            for arg in node.args:
                an,ad=self.rational_function(arg,variable);n,d=self.pmul(n,an),self.pmul(d,ad)
            return n,d
        if isinstance(node,sp.Pow) and isinstance(node.exp,sp.Integer) and abs(node.exp)<=16:
            an,ad=self.rational_function(node.base,variable)
            if node.exp<0:an,ad=ad,an
            n,d=one,one
            for _ in range(abs(int(node.exp))):n,d=self.pmul(n,an),self.pmul(d,ad)
            return n,d
        raise _EvaluationUnresolved('RATIONAL_COEFFICIENT_INTEGRAND_REQUIRED')
    def horner(self,coefficients,x):
        result=self.pair(0)
        for coefficient in reversed(coefficients):result=self.add(self.mul(result,x),self.pair(coefficient))
        return result
    def quadrature(self,n,d,a,b):
        if a==b:return self.pair(0)
        if a>b:return self.neg(self.quadrature(n,d,b,a))
        width=(b-a)/WITNESS_QUADRATURE_PANELS;total=self.pair(0)
        for i in range(WITNESS_QUADRATURE_PANELS):
            x=self.pair(a+i*width,a+(i+1)*width)
            bound=self.div(self.horner(n,x),self.horner(d,x))
            total=self.add(total,self.mul(self.pair(width),bound))
        return total
    def positive_tail(self,n,d,a):
        # x=1/t, including its Jacobian. No root finding or factoring.
        exponent=len(d)-len(n)-2
        if exponent<0:raise _EvaluationUnresolved('IMPROPER_INTEGRAL_CONVERGENCE_NOT_ESTABLISHED')
        numerator=[Fraction(0)]*exponent+list(reversed(n));denominator=list(reversed(d))
        return self.quadrature(numerator,denominator,Fraction(0),1/a)
    def integral(self,node):
        if len(node.limits)!=1 or len(node.limits[0])!=3:
            raise _EvaluationUnresolved('DEFINITE_UNIVARIATE_INTEGRAL_REQUIRED')
        x,a,b=node.limits[0];n,d=self.rational_function(node.function,x)
        if not any(d):raise _EvaluationUnresolved('ZERO_DENOMINATOR')
        if not any(n):return self.pair(0)
        if a is -sp.oo:
            reflected=[v*(-1)**i for i,v in enumerate(n)],[v*(-1)**i for i,v in enumerate(d)]
            if b is sp.oo:
                return self.add(self.add(self.positive_tail(*reflected,Fraction(1)),
                                        self.quadrature(n,d,Fraction(-1),Fraction(1))),
                                self.positive_tail(n,d,Fraction(1)))
            if not isinstance(b,sp.Rational):raise _EvaluationUnresolved('RATIONAL_ENDPOINT_REQUIRED')
            endpoint=Fraction(int(b.p),int(b.q));cut=min(endpoint,Fraction(-1))
            return self.add(self.positive_tail(*reflected,-cut),self.quadrature(n,d,cut,endpoint))
        if not isinstance(a,sp.Rational):raise _EvaluationUnresolved('RATIONAL_ENDPOINT_REQUIRED')
        left=Fraction(int(a.p),int(a.q))
        if b is sp.oo:
            cut=max(left,Fraction(1))
            return self.add(self.quadrature(n,d,left,cut),self.positive_tail(n,d,cut))
        if not isinstance(b,sp.Rational):raise _EvaluationUnresolved('RATIONAL_ENDPOINT_REQUIRED')
        return self.quadrature(n,d,left,Fraction(int(b.p),int(b.q)))


def _fraction_sympy(q):return sp.Rational(q.numerator,q.denominator)


@lru_cache(maxsize=4096)
def _certified_value(value):
    arithmetic=_IntervalArithmetic();certificates=[]
    def numeric(node):
        pair=arithmetic.numeric(node)
        if pair[0]!=pair[1]:certificates.append(sp.Tuple(_fraction_sympy(pair[0]),_fraction_sympy(pair[1])))
        return pair
    def visit(node):
        if node in (sp.true,sp.false):return node,sp.S.Zero
        if isinstance(node,sp.core.relational.Relational):
            a,b=numeric(node.lhs),numeric(node.rhs)
            if isinstance(node,(sp.Equality,sp.Unequality)):
                if a[0]==a[1]==b[0]==b[1]:truth=isinstance(node,sp.Equality)
                elif a[1]<b[0] or b[1]<a[0]:truth=isinstance(node,sp.Unequality)
                else:raise _EvaluationUnresolved('OVERLAPPING_EQUALITY_ENCLOSURES')
            else:
                if isinstance(node,(sp.GreaterThan,sp.StrictGreaterThan)):a,b=b,a
                strict=isinstance(node,(sp.StrictLessThan,sp.StrictGreaterThan))
                if a[1]<b[0] or (not strict and a[1]<=b[0]):truth=True
                elif a[0]>b[1] or (strict and a[0]>=b[1]):truth=False
                else:raise _EvaluationUnresolved('OVERLAPPING_ORDER_ENCLOSURES')
            return sp.true if truth else sp.false,sp.S.Zero
        if isinstance(node,(sp.And,sp.Or,sp.Not,sp.Xor)):
            values=[visit(arg)[0] for arg in node.args]
            if not all(v in (sp.true,sp.false) for v in values):raise _EvaluationUnresolved('UNDECIDED_BOOLEAN_OPERAND')
            return node.func(*values),sp.S.Zero
        if isinstance(node,Str):
            if str(node).startswith(('UNEVALUATED_','EVALUATION_UNRESOLVED_')):
                raise _EvaluationUnresolved(str(node))
            return node,sp.S.Zero
        lo,hi=numeric(node)
        return _fraction_sympy((lo+hi)/2),_fraction_sympy((hi-lo)/2)
    try:
        result,error=visit(value)
        return result,error,sp.Tuple(*certificates),None
    except _EvaluationUnresolved as failure:
        return None,None,None,str(failure)


def _certified_element(element):
    items=list(element) if isinstance(element,sp.Tuple) else [element]
    certificates=[]
    for item in items:
        if isinstance(item,sp.Tuple):
            valid,certificate,reason=_certified_element(item)
        else:
            value,error,certificate,reason=_certified_value(item)
            valid=reason is None
        if not valid:return False,sp.Tuple(),reason
        certificates.extend(certificate)
    return True,sp.Tuple(*certificates),None


# Repeated appearances of the same comparison share the same exact search.
# This cache changes neither the candidates nor the emitted objects.
@lru_cache(maxsize=4096)
def _witness(a,b,context):
    """Finite exact search. Agreement is never evidence of equality.

    Sets get a separately recorded element. Quantified statements only get
    assignments to free variables; their bound variables never enter point.
    The search is intentionally incomplete, so unresolved evaluation is kept.
    """
    global _work_candidate
    whole=sp.Tuple(a,b,context)
    symbols=sorted((x for x in whole.free_symbols
                    if not (x.is_Dummy and x.name=='unresolved_definition_domain')),
                   key=sp.default_sort_key)
    functions=sorted({f.func for f in whole.atoms(AppliedUndef)},key=str)
    if any(len(f.args)!=1 for f in whole.atoms(AppliedUndef)):
        return None,sp.Tuple(Str('SEARCH_UNRESOLVED_NONRADIAL_FUNCTION'))
    attempted=0
    unresolved=Counter()
    for profile_index in range(6 if functions else 1):
        profiles=tuple(_profile_bank(profile_index,functions).items())
        for point_index in range(6 if symbols else 1):
            candidate='profile='+str(profile_index)+' point='+str(point_index)
            _work_candidate=candidate
            work_detail('witness:profile='+str(profile_index)+' point='+str(point_index))
            point=_point(symbols,point_index)
            if point is None:continue
            points=tuple(point.items())
            attempted+=1
            allowed=_evaluate_at(context,profiles,points)
            allowed,allowed_error,allowed_certificate,reason=_certified_value(allowed)
            if reason:unresolved['COMMON_DOMAIN:'+reason]+=1
            if allowed is not sp.true:continue
            left=_evaluate_at(a,profiles,points)
            right=_evaluate_at(b,profiles,points)
            elements=(None,)
            if isinstance(left,sp.Set) and isinstance(right,sp.Set):
                elements=_set_elements(left,right)
            elif isinstance(a,sp.Set) and isinstance(b,sp.Set):
                # Failed set evaluation does not erase its element candidates.
                # Each original candidate is visited and recorded unresolved.
                elements=_set_elements(a,b)
            for element_index,element in enumerate(elements):
                _work_candidate=candidate+' element='+str(element_index)+'/'+str(len(elements))
                work_detail('membership:profile='+str(profile_index)+' point='+str(point_index)+
                            ' element='+str(element_index)+'/'+str(len(elements)))
                x,y=left,right
                element_certificate=sp.Tuple()
                if element is not None:
                    valid,element_certificate,reason=_certified_element(element)
                    if not valid:
                        unresolved['SET_ELEMENT:'+reason]+=1
                        continue
                    try:
                        x=_membership(left,element) if isinstance(left,sp.Set) else left
                        y=_membership(right,element) if isinstance(right,sp.Set) else right
                    except TypeError as error:
                        if not any(message in str(error) for message in ('Invalid NaN comparison','Invalid comparison of non-real')):raise
                        continue
                    except ZeroDivisionError:
                        continue
                x,x_error,x_certificate,x_reason=_certified_value(x)
                y,y_error,y_certificate,y_reason=_certified_value(y)
                if x_reason:unresolved['BASELINE:'+x_reason]+=1
                if y_reason:unresolved['CORRUPTED:'+y_reason]+=1
                if x_reason or y_reason:continue
                if x in (sp.true,sp.false) and y in (sp.true,sp.false):
                    different=x!=y
                elif isinstance(x,Str) and isinstance(y,Str):different=x!=y
                elif isinstance(x,sp.Rational) and isinstance(y,sp.Rational):
                    different=abs(y-x)>x_error+y_error
                else:continue
                interval_used=bool(x_certificate or y_certificate or allowed_certificate or element_certificate)
                if different:
                    return sp.Tuple(Str('CHANGE'),
                        sp.Tuple(Str('POINT'),*(sp.Tuple(k,v) for k,v in points)),
                        sp.Tuple(Str('PROFILES'),*(sp.Tuple(Str(str(k)),v) for k,v in profiles)),
                        sp.Tuple(Str('SET_ELEMENT'),element if element is not None else Str('NOT_APPLICABLE')),
                        sp.Tuple(Str('BASELINE_EVALUATION'),x),
                        sp.Tuple(Str('CORRUPTED_EVALUATION'),y),
                        sp.Tuple(Str('COMMON_DOMAIN_EVALUATION'),allowed),
                        sp.Tuple(Str('ARITHMETIC'),Str('CERTIFIED_RATIONAL_INTERVAL' if interval_used else 'EXACT')),
                        sp.Tuple(Str('PRECISION_BITS'),sp.Integer(WITNESS_PRECISION_BITS) if interval_used else Str('EXACT')),
                        sp.Tuple(Str('ERROR_BOUND_BOTH_EVALUATIONS'),x_error+y_error),
                        sp.Tuple(Str('BASELINE_ERROR_BOUND'),x_error),
                        sp.Tuple(Str('CORRUPTED_ERROR_BOUND'),y_error),
                        sp.Tuple(Str('EVALUATION_ENCLOSURES'),x_certificate,y_certificate,allowed_certificate,element_certificate)),None
    return None,sp.Tuple(Str('WITNESS_SEARCH'),Str('fixed rational points and six explicit radial functions; exact or certified outward rational interval arithmetic'),
        sp.Tuple(Str('ATTEMPTS'),attempted),
        sp.Tuple(Str('EVALUATION_UNRESOLVED'),*(sp.Tuple(Str(reason),sp.Integer(count)) for reason,count in sorted(unresolved.items()))),
        Str('NO_WITNESS_FOUND; agreement or unresolved evaluation proves nothing'))


def _value_difference(a,b,context=sp.true):
    # Equality is always tried before evaluations, even for compact K11.
    if a==b or context is sp.false:return sp.S.Zero
    if isinstance(a,sp.Tuple) and isinstance(b,sp.Tuple) and len(a)==len(b):
        parts=sp.Tuple(*(_value_difference(x,y,context) for x,y in zip(a,b)))
        return sp.S.Zero if all(part==0 for part in parts) else parts
    if isinstance(a,sp.MatrixBase) and isinstance(b,sp.MatrixBase) and a.shape==b.shape:
        parts=sp.Tuple(*(_value_difference(x,y,context) for x,y in zip(a,b)))
        return sp.S.Zero if all(part==0 for part in parts) else sp.Tuple(Str('MATRIX_ENTRY_DIFFERENCES'),parts)
    if isinstance(a,sp.Expr) and isinstance(b,sp.Expr):
        work_detail('equality:algebra')
        residual=_arithmetic_reduce(b-a)
        if _known_zero(residual):return sp.S.Zero
    elif isinstance(a,sp.Basic) and isinstance(b,sp.Basic):
        work_detail('equality:normalize')
        left,right=_normalized_predicate(a),_normalized_predicate(b)
        if left==right:return sp.S.Zero
        if isinstance(left,sp.logic.boolalg.Boolean) and isinstance(right,sp.logic.boolalg.Boolean):
            residual=sp.simplify_logic(sp.Xor(left,right),deep=False)
            if residual is sp.false:return sp.S.Zero
        elif isinstance(left,sp.Set) and isinstance(right,sp.Set):
            residual=sp.SymmetricDifference(left,right,evaluate=False)
        else:residual=sp.Tuple(Str('EQUALITY_QUERY'),a,b)
    else:residual=sp.Tuple(Str('EQUALITY_QUERY'),a,b)
    change,search=_witness(a,b,context)
    return change if change is not None else undecided(a,b,sp.Tuple(residual,search))



def payload_roles(value):
    """Separate explicitly attached domains, never Boolean-valued objects.

    Integer argument paths identify every domain slot, including nested
    condition/velocity cases and copies inside export metadata. The sample
    configuration is supplied outside this function and is not a domain.
    """
    domains={}
    def visit(node,path):
        if isinstance(node,sp.Tuple) and len(node)==2 and node[0]==Str('ATTACHED_DOMAIN'):
            domains[path]=node[1]
            return Str('ATTACHED_DOMAIN_SLOT')
        if isinstance(node,sp.Basic):
            args=tuple(visit(child,path+(index,)) for index,child in enumerate(node.args))
            if args!=node.args:
                return node.func(*args)
        # Export metadata also contains FunctionClass objects, which are
        # classes, not Basic instances. Do not invoke their unbound xreplace
        # method while projecting the enclosing tuple.
        return node
    return visit(value,()),domains


def _attached_predicate(domain,parent):
    """Adapt the engine's domain carriers without treating a set as a bool."""
    if isinstance(domain,sp.logic.boolalg.Boolean) and not isinstance(domain,sp.Expr):
        return domain
    if isinstance(domain,sp.Tuple):
        return sp.And(*(_attached_predicate(part,parent) for part in domain))
    if (isinstance(domain,sp.Set) and isinstance(parent,sp.Tuple) and len(parent)>1
            and parent[0]==Str('RADIAL_DIFFERENTIAL_IDENTITY')):
        # This carrier explicitly names its radial coordinate before the set.
        return sp.Contains(parent[1],domain)
    # A carrier without an explicit point association remains opaque for
    # value witnesses. Its actual set/condition is still compared separately.
    return sp.Dummy('unresolved_definition_domain')


def role_difference(a,b,value_method):
    left,left_domains=payload_roles(a)
    right,right_domains=payload_roles(b)
    def compare(x,y,path=(),context=sp.true):
        if x==y:return sp.S.Zero
        local=[_attached_predicate(domain,parent) for table,parent in ((left_domains,x),(right_domains,y))
               for slot,domain in table.items() if slot[:-1]==path]
        context=sp.And(context,*local)
        if isinstance(x,sp.Tuple) and isinstance(y,sp.Tuple) and len(x)==len(y):
            parts=sp.Tuple(*(compare(a,b,path+(j,),context) for j,(a,b) in enumerate(zip(x,y))))
            return sp.S.Zero if all(part==0 for part in parts) else parts
        with work_part('value'+str(path)):
            return value_method(x,y,context)
    value_delta=compare(left,right)
    if left_domains==right_domains:
        domain_delta=sp.S.Zero
    else:
        rows=[]
        for path in sorted(left_domains.keys()|right_domains.keys()):
            before=left_domains.get(path,Str('NO_ATTACHED_DOMAIN'))
            after=right_domains.get(path,Str('NO_ATTACHED_DOMAIN'))
            with work_part('domain'+str(path)):
                rows.append(sp.Tuple(sp.Tuple(*path),_value_difference(before,after)))
        domain_delta=sp.S.Zero if all(row[1]==0 for row in rows) else sp.Tuple(*rows)
    return sp.Tuple(sp.Tuple(Str('VALUE_DIFFERENCE'),value_delta),
                    sp.Tuple(Str('DOMAIN_DIFFERENCE'),domain_delta))


def difference(a,b):
    return role_difference(a,b,_value_difference)


def compact_difference(a,b):
    return role_difference(a,b,_compact_difference)


def run(script, folder):
    folder.mkdir(parents=True,exist_ok=False)
    output = folder/'stdout'
    with output.open('w') as out, (folder/'stderr').open('w') as err:
        progress(str(folder)+': engine start')
        child = subprocess.Popen([sys.executable,'-B',str(script)],stdout=out,stderr=subprocess.PIPE,
                                 text=True,stdin=subprocess.DEVNULL,
                                 env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
        # Relay operational progress immediately, keeping stdout purely the
        # payload grammar and preserving a separate worker stderr transcript.
        for line in child.stderr:
            err.write(line)
            err.flush()
            print(str(folder)+': '+line.rstrip('\n'),file=sys.stderr,flush=True)
        code = child.wait()
        child.stderr.close()
    (folder/'exit-code').write_text(str(code)+'\n')
    progress(str(folder)+': engine exit '+str(code))
    if code:
        raise RuntimeError(('worker stopped',str(folder),code))
    return output


def main():
    start_progress_monitor()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scratch',type=Path,required=True)
    parser.add_argument('--baseline-output',type=Path)
    parser.add_argument('--knife',choices=tuple(KNIVES),action='append')
    parser.add_argument('--reuse-worker-root',type=Path,help='Reuse only successful, byte-identical completed mutant engine jobs; comparisons are always recomputed.')
    args = parser.parse_args()
    scratch = args.scratch.resolve()
    scratch.relative_to(REPOSITORY/'_scratch')
    scratch.mkdir(parents=True,exist_ok=False)
    source = ENGINE.read_text()
    digest = hashlib.sha256(ENGINE.read_bytes()).hexdigest()
    if args.baseline_output:
        output = args.baseline_output.resolve()
        if output.with_name(output.name+'.source-sha256').read_text().strip()!=digest:
            raise ValueError('baseline source pin mismatch')
    else:
        with work_part('baseline engine'):
            output = run(ENGINE,scratch/'baseline')
        output.with_name(output.name+'.source-sha256').write_text(digest+'\n')
    progress('baseline: parse '+str(output))
    with work_part('parse baseline'):
        baseline = parse(output)
    for knife in args.knife or KNIVES:
        before,after = KNIVES[knife]
        count = source.count(before)
        print('PY_LOCAL_S9B_ABLATION_'+knife+'_CONSTRUCTION: '+render(sp.Tuple(Str(before),Str(after),sp.Integer(count))),flush=True)
        if count != 1:
            raise ValueError(('construction-site count',knife,count))
        replica = scratch/knife/'replica'
        scripts = replica/'scripts'
        scripts.mkdir(parents=True)
        (replica/'directives').symlink_to(ROOT/'directives',target_is_directory=True)
        for name in ('ledger_fold.py','S11c_b_exports.py','S11c_c1_exports.py','S11c_c2_exports.py'):
            (scripts/name).symlink_to(ROOT/'scripts'/name)
        mutant = scripts/ENGINE.name
        mutant.write_text(source.replace(before,after,1))
        if knife == 'K11':
            with work_part('K11 compact'):
                compact_k11(baseline,source,mutant,scratch/knife)
            continue
        previous = args.reuse_worker_root/knife if args.reuse_worker_root else None
        if previous is not None and (previous/'run/exit-code').exists():
            if (previous/'run/exit-code').read_text().strip()!='0':
                raise ValueError(('cannot reuse unsuccessful worker',knife))
            if (previous/'replica/scripts'/ENGINE.name).read_bytes()!=mutant.read_bytes():
                raise ValueError(('cannot reuse different worker source',knife))
            worker_output=previous/'run/stdout'
            progress(knife+': reuse completed byte-identical engine '+str(worker_output))
            (scratch/knife/'reused-worker-output').write_text(str(worker_output.resolve())+'\n')
        else:
            with work_part(knife+' engine'):
                worker_output = run(mutant,scratch/knife/'run')
        progress(knife+': parse '+str(worker_output))
        with work_part('parse '+knife):
            corrupted = parse(worker_output)
        for tag in sorted(set(baseline)|set(corrupted)):
            progress(knife+' '+tag+': baseline, corrupted and difference')
            a = baseline.get(tag,Str('MISSING_TAG'))
            b = corrupted.get(tag,Str('MISSING_TAG'))
            with work_part(knife+' '+tag):
                delta=difference(a,b)
                work_detail('output:render')
                print('PY_LOCAL_S9B_ABLATION_'+knife+'_'+tag+': '+render(sp.Tuple(a,b,delta)),flush=True)
        if set(baseline)!=set(corrupted):
            raise ValueError(('tag inventory changed',knife))


if __name__=='__main__':
    main()
