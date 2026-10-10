#!/usr/bin/env python3
"""S9b Parts A-C: general radial profile functionals, finite optical grades,
fixed-endpoint radar logarithmic slope, and Abel-reduced radial identities.
Only physics input: S9b_SHARED_PHYSICS.md. Emit computed objects, not verdicts.
Only writes: stdout and S9b_exports.py. Run through the 16 GiB s9b guard.
"""
from __future__ import annotations

import argparse
from collections.abc import Mapping
from itertools import product
from pathlib import Path
import hashlib
import sys

sys.dont_write_bytecode = True
import sympy as sp
from sympy.core.function import AppliedUndef
from sympy.core.relational import Relational
from sympy.core.symbol import Str
from sympy.logic.boolalg import Boolean
from sympy.printing.repr import ReprPrinter

HERE = Path(__file__).resolve()
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE.parent))
from ledger_fold import (load_model, check_consumer,
                         assert_lookups_equal_manifest, assert_delta_is_minimal)

IMPORT_KEYS = ('c_s0',)
INPUTS = (HERE, ROOT/'directives/S9b_SHARED_PHYSICS.md',
          *(ROOT/('scripts/'+s+'_exports.py') for s in
            ('S11c_b', 'S11c_c1', 'S11c_c2')), ROOT/'scripts/ledger_fold.py')
GRADES = tuple(product(range(2), range(3), range(2)))
ZERO = (0, 0, 0)
RADIAL_FREEZE = None  # K4a-c: derivative construction
RADIAL_DATA = None
TAGS = {}
LOCAL = []


def cas(value):
    if isinstance(value, str):
        return Str(value)
    if isinstance(value, Mapping):
        return sp.Tuple(*(sp.Tuple(cas(k), cas(v)) for k, v in value.items()))
    if isinstance(value, (tuple, list, set, frozenset)):
        return sp.Tuple(*(cas(v) for v in value))
    return sp.sympify(value)


def emit(name, value, *, local=False):
    local = local or name not in SHARED
    tag = ('PY_LOCAL_S9B_' if local else 'PY_S9B_') + name
    if tag in TAGS:
        raise ValueError(('duplicate tag', tag))
    obj = cas(value)
    TAGS[tag] = obj
    if local:
        LOCAL.append(tag)
    print(tag+': '+export_srepr(obj), flush=True)
    return obj


def eq(a, b):
    return sp.Eq(a, b, evaluate=False)


def everywhere(variable, condition, domain):
    # An empty counterexample set is universal quantification, not sampling.
    return eq(sp.ConditionSet(variable, sp.Not(condition), domain), sp.EmptySet)


class Jet:
    """R[d,V,H]/(d**2,V**3,H**2); discarded grades never get evaluated."""
    def __init__(self, values=0):
        self.c = (dict(values) if isinstance(values, dict) else
                  {ZERO: sp.sympify(values)})
        self.c = {g: v for g, v in self.c.items() if v != 0}

    def __getitem__(self, g):
        return self.c.get(g, sp.S.Zero)

    def __add__(self, other):
        other = asjet(other)
        return Jet({g: self[g]+other[g] for g in GRADES})
    __radd__ = __add__

    def __neg__(self):
        return Jet({g: -v for g, v in self.c.items()})

    def __sub__(self, other):
        return self+-asjet(other)

    def __rsub__(self, other):
        return asjet(other)+-self

    def __mul__(self, other):
        other = asjet(other)
        out = {}
        for g, a in self.c.items():
            for h, b in other.c.items():
                k = tuple(x+y for x, y in zip(g, h))
                if k in GRADES:
                    out[k] = out.get(k, 0)+a*b
        return Jet(out)
    __rmul__ = __mul__

    def __pow__(self, exponent):
        exponent = sp.sympify(exponent)
        if exponent.is_Integer and exponent >= 0:
            out = Jet(1)
            for _ in range(int(exponent)):
                out = out*self
            return out
        base = self[ZERO]
        if base == 0:
            raise ValueError(('nonanalytic jet power', exponent))
        tail = (self-base)*(1/base)
        out, power = Jet(1), Jet(1)
        for k in range(1, 5):
            power = power*tail
            out += sp.binomial(exponent, k)*power
        return out*base**exponent

    def __truediv__(self, other):
        return self*asjet(other)**-1

    def map(self, fn):
        return Jet({g: fn(v) for g, v in self.c.items()})

    def diff(self, variable, order=1):
        return self.map(lambda v: radial_diff(v, variable, order))

    def tidy(self):
        return self.map(lambda v: sp.factor(sp.expand(v)))


def radial_diff(value, variable, order=1):
    if RADIAL_FREEZE is None:
        return sp.diff(value,variable,order)
    name=('s9b_delta','s9b_V','s9b_xi')[RADIAL_FREEZE]
    frozen={node:sp.Dummy('s9b_frozen_profile') for node in value.atoms(AppliedUndef)
            if node.func.__name__==name}
    held=value.xreplace(frozen)
    held=held.replace(lambda a:isinstance(a,sp.Derivative),lambda a:a.doit())
    return sp.diff(held,variable,order).xreplace({v:k for k,v in frozen.items()})


def asjet(value):
    return value if isinstance(value, Jet) else Jet(value)


def evaluate_jet(expression, mapping):
    if expression in mapping:
        return mapping[expression]
    if not any(expression.has(k) for k in mapping):
        return Jet(expression)
    if expression.is_Add:
        return sum((evaluate_jet(a, mapping) for a in expression.args), Jet())
    if expression.is_Mul:
        out = Jet(1)
        for a in expression.args:
            out *= evaluate_jet(a, mapping)
        return out
    if expression.is_Pow:
        return evaluate_jet(expression.base, mapping)**expression.exp
    raise ValueError(('unsupported analytic jet expression', expression))


def compose_radial(jet, old, new, radial):
    shift = radial-new
    out = Jet()
    for g, value in jet.c.items():
        grade = Jet({g: 1})
        power = Jet(1)
        for k in range(5):
            term = grade*power
            if term.c:
                out += term*radial_diff(value, old, k).subs(old, new)/sp.factorial(k)
            power *= shift
    return out.tidy()


def compare(left, right):
    """Total three-valued equality; predicates never enter subtraction."""
    if left == right:
        return sp.true
    if isinstance(left, Mapping) or isinstance(right, Mapping):
        if not isinstance(left, Mapping) or not isinstance(right, Mapping) or left.keys() != right.keys():
            return sp.false
        values = [compare(left[k], right[k]) for k in left]
    elif isinstance(left, (tuple, list, sp.Tuple, sp.MatrixBase)) or isinstance(right, (tuple, list, sp.Tuple, sp.MatrixBase)):
        if type(left) is not type(right) or len(left) != len(right):
            return sp.false
        if isinstance(left, sp.MatrixBase) and left.shape != right.shape:
            return sp.false
        values = [compare(a, b) for a, b in zip(left, right)]
    elif isinstance(left, (str, Str)) or isinstance(right, (str, Str)):
        return sp.false
    elif isinstance(left, (Relational, Boolean)) or isinstance(right, (Relational, Boolean)):
        if isinstance(left, Relational) and isinstance(right, Relational) and left.func == right.func:
            values = [compare(a, b) for a, b in zip(left.args, right.args)]
        else:
            result = sp.simplify(sp.Equivalent(left, right)) if isinstance(left, Boolean) and isinstance(right, Boolean) else None
            return result if result in (sp.true, sp.false) else Str('UNDECIDED')
    elif isinstance(left, sp.Expr) and isinstance(right, sp.Expr):
        result = sp.simplify(left-right).is_zero
        return sp.sympify(result) if result is not None else Str('UNDECIDED')
    else:
        return Str('UNDECIDED')
    if sp.false in values:
        return sp.false
    return sp.true if all(v == sp.true for v in values) else Str('UNDECIDED')


def row(value, class_tag='DERIVED'):
    return {'value': value, 'display': str(value), 'value_kind': 'COMPUTED_OBJECT',
            'class': class_tag, 'step': 'S9b',
            'description': 'S9b_SHARED_PHYSICS.md Parts A-C; general-profile radial differential conditions',
            'evidence': 'S9b_light_bending_sympy_audit.py; emitted construction and condition cases'}


class ExportReprPrinter(ReprPrinter):
    """Preserve undefined-function identity, including its assumptions.

    SymPy's default srepr prints only Function(name), dropping _kwargs.
    These kwargs participate in UndefinedFunction equality (and its pickle
    reconstruction); dropping them changes the object revived by D3.
    """
    def _print_FunctionClass(self, expr):
        if issubclass(expr, AppliedUndef):
            return 'Function(%r, **%r)' % (expr.__name__, dict(sorted(expr._kwargs.items())))
        return super()._print_FunctionClass(expr)


def export_srepr(value):
    return ExportReprPrinter().doprint(value)


def roundtrip_equal(live, decoded):
    # D3 checks faithful serialization, not mathematical equivalence. Exact
    # reconstruction also checks ConditionSet binders and function assumptions.
    # A mismatch must be emitted and rejected, never sent to Boolean simplify.
    return sp.sympify(live == decoded)


def publish(fold, objects, declarations):
    candidates = {k: row(v) for k, v in objects.items()}
    candidates.update({k: row(v, c) for k, (v, c) in declarations.items()})
    routed, routes, roots = {}, {}, []
    for key, record in candidates.items():
        actual = key
        if key in fold:
            prior = fold[key]['value']
            same_object = key in IMPORT_KEYS
            outcome = compare(prior, record['value']) if same_object else Str('UNDECIDED')
            route = 'F9B_EQUAL' if outcome == sp.true else 'F9C_PRESERVE'
            if outcome != sp.true:
                actual = 's9b_'+key
            record['f9_operands'] = sp.Tuple(prior, record['value'])
            record['f9_comparison'] = outcome
            record['corroborated_steps'] = (fold[key].get('step'), 'S9b') if outcome == sp.true else ()
            routes[key] = (actual, outcome, record['f9_operands'])
        else:
            route = 'F9A_ABSENT'
        record['route'] = route
        if actual in routed or (actual != key and actual in fold):
            raise ValueError(('routed key collision', actual))
        routed[actual] = record
        if key in objects:
            roots.append(actual)
    emit('EXPORT_ROUTES', routes, local=True)
    combined = dict(fold) | routed
    closure = check_consumer(combined, roots)['closure']
    own = set(routed).intersection(closure)
    delta = {k: routed[k] for k in sorted(own)}
    emit('EXPORT_CLOSURE', sorted(own), local=True)
    assert_delta_is_minimal(delta, own)
    digests = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in INPUTS}
    relational_code = {
        'Equality': 'Eq', 'Unequality': 'Ne', 'StrictGreaterThan': 'Gt',
        'StrictLessThan': 'Lt', 'GreaterThan': 'Ge', 'LessThan': 'Le'}
    lines = ['# Generated S9b own-rows delta. General-profile conditional objects.',
             'from types import MappingProxyType', 'import sympy as sp',
             'from sympy.core.symbol import Str',
             'from sympy.functions.elementary.piecewise import ExprCondPair',
             '_RELATIONALS = {']
    lines += [repr(k)+': lambda a, b: sp.'+v+'(a, b, evaluate=False),' for k, v in relational_code.items()]
    lines += ['}', 'def _restore(s):',
              "    return eval(s, {'__builtins__': {}, **vars(sp), 'Str': Str, 'ExprCondPair': ExprCondPair, **_RELATIONALS})",
              'IMPORT_KEYS = '+repr(IMPORT_KEYS), 'EXPORT_ROOTS = '+repr(tuple(roots)),
              'BUILD_INPUT_DIGESTS = MappingProxyType('+repr(digests)+')', '_LEDGER = {']
    for k, record in delta.items():
        fields = []
        for field, value in record.items():
            encoded = '_restore('+repr(export_srepr(value))+')' if isinstance(value, sp.Basic) or isinstance(value, sp.FunctionClass) else repr(value)
            fields.append(repr(field)+': '+encoded)
        lines.append(repr(k)+': {'+', '.join(fields)+'},')
    lines += ['}', 'LEDGER = MappingProxyType({k: MappingProxyType(v) for k,v in _LEDGER.items()})', 'del _LEDGER', '']
    source = '\n'.join(lines)
    namespace = {}
    exec(compile(source, str(ROOT/'scripts/S9b_exports.py'), 'exec'), namespace)
    restored = namespace['LEDGER']
    checks = {}
    for k in delta:
        live, decoded = delta[k]['value'], restored[k]['value']
        # Serialization check, not a physics residual or an F9b residual.
        checks[k] = (live, decoded, roundtrip_equal(live, decoded))
    emit('EXPORT_ROUNDTRIP', checks, local=True)
    if any(v[2] != sp.true for v in checks.values()):
        raise ValueError('D3 serialization round-trip')
    restored_closure = check_consumer(dict(fold) | dict(restored), roots)['closure']
    assert_delta_is_minimal(restored, set(restored).intersection(restored_closure))
    (ROOT/'scripts/S9b_exports.py').write_text(source)
    emit('BUILD_INPUT_DIGESTS', digests, local=True)


def shared_vocabulary():
    grades = ['G%d%d%d' % g for g in GRADES]
    loggrades = ('G100','G010','G020','G001')
    obs = ('DEFLECTION','RADAR')
    responses = ('CONSTANT','FIXED_RATIO','POWER')
    branch = ['BRANCH_EXISTENCE','PATH_TRAVERSAL','BRANCH_TYPE']
    a = ['A_'+name+'_'+g for name in ('DEFLECTION','ROUND_TRIP','ONE_WAY_ER','ONE_WAY_RE','NONRECIPROCAL') for g in grades]
    a += ['A_RADAR_LOG_'+g for g in loggrades]
    b = ['B_'+name+'_'+o for name in ('GAMMA','RESIDUAL','CONDITION','IMPLIED_JN') for o in obs]+['B_GAMMA_DIFFERENCE']
    c = ['C_'+response+'_'+name+'_'+o for response in responses for name in ('CONDITION','IMPLIED_JN') for o in obs]+['C_'+response+'_N_DEPENDENCE' for response in responses]
    names = branch+a+['A_NONRECIPROCAL_PATH_DEPENDENCE']+b+c
    ar = ['A_DEFLECTION_'+g for g in grades]+['A_RADAR_LOG_'+g for g in loggrades]
    names += ['R_'+kind+'_'+name for kind in ('FLOW_ONLY','SPEED_ONLY','TILT_ONLY') for name in ar+b]
    names += ['R_BULK_ONLY_C_'+response+'_'+name+'_'+o for response in ('FIXED_RATIO','POWER') for name in ('CONDITION','IMPLIED_JN') for o in obs]
    names += ['F_'+stage+'_'+name for stage in ('FLOW','LIVE') for name in branch+a+[n for n in b if 'IMPLIED_JN' not in n]]
    names += ['F_LIVE_C_'+response+'_CONDITION_'+o for response in responses for o in obs]+['F_MASS_SOLUTION']
    return frozenset(names)


SHARED = shared_vocabulary()
Delta = sp.Function('s9b_delta',real=True)
Velocity = sp.Function('s9b_V',real=True)
Embedding = sp.Function('s9b_xi',real=True)
Density = sp.Function('rho_br_profile',positive=True)
BulkChange = sp.Function('s9b_f',real=True)


def profile_substitute(value, replacements):
    """Substitute inside integrals and quantified gates, preserving binders."""
    obj = cas(value)
    if not replacements:
        return obj
    obj = obj.replace(lambda a:isinstance(a,AppliedUndef) and a.func in replacements,
                      lambda a:replacements[a.func](a.args[0]))
    obj = obj.replace(lambda a:isinstance(a,sp.Integral) and a.function==0,lambda a:sp.S.Zero)
    obj = obj.replace(lambda a:isinstance(a,sp.Derivative) and not a.has(sp.Integral),lambda a:a.doit())
    obj = obj.replace(lambda a:isinstance(a,sp.Subs) and not a.has(sp.Integral),lambda a:a.doit())
    return obj


def vector_dispersion(c0,r):
    U,A,chi,omega = sp.symbols('s9b_U s9b_A s9b_chi s9b_omega',real=True)
    momenta = sp.Matrix(sp.symbols('s9b_k_r s9b_k_t s9b_k_z',real=True))
    azimuthal_component = sp.S.Zero  # K11: vector advection construction
    advected_velocity = sp.Matrix([U,azimuthal_component,0])
    kinetic = (omega-advected_velocity.dot(momenta))**2  # K1
    inverse_metric = sp.diag(1/A,1,1)  # K2
    local_speed_squared = chi  # K3
    dispersion = kinetic-local_speed_squared*(momenta.T*inverse_metric*momenta)[0]
    quadratic = sp.hessian(dispersion,momenta)/2
    linear = sp.Matrix([sp.diff(dispersion,k).subs(dict.fromkeys(momenta,0))/(2*omega) for k in momenta])
    centre = (quadratic.inv()*(-linear)).applyfunc(sp.factor)
    J = sp.Symbol('s9b_J',real=True)
    radial = dispersion.subs({momenta[1]:omega*J/r,momenta[2]:0})
    aa,bb,cc = sp.Poly(radial,momenta[0]).all_coeffs()
    beta = sp.factor(-bb/(2*aa*omega))
    even_squared = sp.factor((bb**2-4*aa*cc)/(4*aa**2*omega**2))
    arr = sp.factor(even_squared.subs(J,0))
    aphi = sp.factor(-arr/(sp.diff(even_squared,J,2)/2))
    return dict(U=U,A=A,chi=chi,omega=omega,k=momenta,J=J,
                advection=advected_velocity,dispersion=dispersion,radial=radial,
                centre=centre,beta=beta,even_squared=even_squared,arr=arr,aphi=aphi)


def nonreciprocity(r,c0,delta,velocity,embedding):
    """Engine construction also called by amendment-2 compact evaluation."""
    source = vector_dispersion(c0,r)
    xyz = sp.symbols('x_1 x_2 x_3',real=True)
    radius = sp.sqrt(sum(x*x for x in xyz))
    transverse = sp.sqrt(xyz[0]**2+xyz[1]**2)
    radial_unit = sp.Matrix([sp.diff(radius,x) for x in xyz])
    azimuthal_unit = sp.Matrix([-xyz[1],xyz[0],0])/transverse
    dj = Jet({(1,0,0):delta.subs(r,radius)})
    vj = Jet({(0,1,0):velocity.subs(r,radius)})
    hj = Jet({(0,0,1):radial_diff(embedding,r).subs(r,radius)**2})
    mapping = {source['U']:vj,source['A']:1+hj,source['chi']:(c0*(1+dj))**2}
    for symbol in source['advection'].free_symbols:
        if symbol.name == 's9b_azimuthal':
            mapping[symbol] = Jet({(0,1,0):symbol*transverse/radius**2})
    components = [evaluate_jet(component,mapping).tidy() for component in source['centre']]
    oneform = sp.ImmutableMatrix(radial_unit*sum(components[0].c.values(),sp.S.Zero)+
                                 azimuthal_unit*sum(components[1].c.values(),sp.S.Zero))
    oneform = oneform.applyfunc(sp.factor)
    curl = sp.ImmutableMatrix([sp.simplify(radial_diff(oneform[(i+2)%3],xyz[(i+1)%3])-
                                         radial_diff(oneform[(i+1)%3],xyz[(i+2)%3])) for i in range(3)])
    return dict(oneform=oneform,curl=curl,predicate=sp.Ne(curl,sp.zeros(3,1)),
                xyz=xyz,source=source,grades=GRADES)


def gated(value, domain):
    return cas(((domain,value),(sp.Not(domain),'NOT_ESTABLISHED')))


def build(branch_demo=None):
    r,x,B,b,c0,ZE,ZR = sp.symbols('s9b_r s9b_X s9b_B b c_0 Z_E Z_R',positive=True)
    RE,RR,rt,bmin = sp.symbols('r_E r_R s9b_r_turn s9b_b_min',positive=True)
    GM,K,n,s = sp.symbols('GM K n s',real=True)
    rho0,m = sp.symbols('rho_0 m',positive=True)
    rho_bulk = sp.Symbol('rho_bulk',positive=True)
    jn = sp.Function('j_n',real=True)(r)
    delta,velocity,xi,rho,f = Delta(r),Velocity(r),Embedding(r),Density(r),BulkChange(r)
    geometry = vector_dispersion(c0,r)
    U,A,chi,omega,J = (geometry[k] for k in ('U','A','chi','omega','J'))
    kr = geometry['k'][0]
    beta,arr,aphi = (geometry[k] for k in ('beta','arr','aphi'))
    slope = radial_diff(xi,r)
    radial = geometry['radial']
    radial_omega = sp.solve(radial.subs(J,0),omega)
    speeds = [sp.diff(root,kr) for root in radial_omega]
    discriminant = sp.factor(sp.discriminant(radial.subs(J,0),omega)/(4*kr**2))
    existence = sp.And(sp.Gt(discriminant,0),sp.Ne(rho,0))
    traversal = sp.Lt(sp.factor(sp.prod(speeds)),0)
    profile_map = {U:velocity,A:1+slope**2,chi:(c0*(1+delta))**2}
    endpoint_map = {RE:sp.sqrt(b*b+ZE*ZE),RR:sp.sqrt(b*b+ZR*ZR)}
    real_domain = everywhere(r,existence.subs(profile_map),sp.Interval(rt,sp.oo))
    travel_domain = everywhere(r,traversal.subs(profile_map),sp.Interval(rt,sp.oo))
    branch_type = ((sp.Lt(discriminant,0),('GROWING','DECAYING')),
                   (sp.Eq(discriminant,0),'ABSENT'),
                   (sp.And(existence,sp.Not(traversal)),'UNABLE_TO_TRAVERSE_REQUIRED_DIRECTION'))
    emit('SUPPLIED_DEPENDENCIES',('DISPERSION','SPEED_IDENTIFICATION','ADVECTION','METRIC','COORDINATE_MASS_BALANCE','LAB_HELD','ORDER_COUNTING','REFERENCES','PART_C_RESPONSES'))
    emit('GENERAL_PROFILES',(delta,velocity,xi,rho,f))
    emit('COUNTING',(eq(sp.Function('epsilon')(r),GM/(c0*c0*r)),
                     tuple((g,sp.Rational(g[0])+sp.Rational(g[1],2)+g[2]) for g in GRADES)))
    emit('VECTOR_DISPERSION',geometry['dispersion'])
    emit('VECTOR_ADVECTION',geometry['advection'])
    emit('SPEED_IDENTIFICATION',eq((c0*(1+delta))**2,sp.Function('mu_perp')(r)/rho))
    emit('BRANCH_EXISTENCE',real_domain)
    emit('PATH_TRAVERSAL',travel_domain)
    emit('BRANCH_TYPE',cas(branch_type).subs(profile_map))
    if branch_demo is not None:
        values={'elliptic':{chi:-c0*c0,U:0,A:1},'degenerate':{chi:0,U:0,A:1},'blocked':{chi:c0*c0,U:-2*c0,A:1}}[branch_demo]
        domain=sp.And(existence,traversal).subs(values)
        emit('DEMONSTRATION_BRANCH',(existence.subs(values),traversal.subs(values),domain))
        emit('DEMONSTRATION_BRANCH_TYPE',cas(branch_type).subs(values))
        if domain == sp.false:
            emit('DEMONSTRATION_OBSERVABLES',('NOT_ESTABLISHED',domain))
            emit('LOCAL_NAMES',LOCAL+['PY_LOCAL_S9B_LOCAL_NAMES'])
            return
    fold,_ = load_model(*(ROOT/('scripts/'+s+'_exports.py') for s in ('S11c_b','S11c_c1','S11c_c2')))
    manifest=check_consumer(fold,IMPORT_KEYS)
    bound=assert_lookups_equal_manifest(lambda ledger:{k:ledger[k]['value'] for k in IMPORT_KEYS},fold,IMPORT_KEYS)
    emit('IMPORT_MANIFEST',(IMPORT_KEYS,sorted(bound['lookups']),sorted(manifest['closure'])))
    pressure=K*rho_bulk**n
    cs_squared=sp.diff(pressure,rho_bulk)/m
    cs0_squared=cs_squared.subs(rho_bulk,rho0)
    emit('BULK_SOUND_ANCHOR',eq(bound['result']['c_s0']**2,cs0_squared))
    emit('IMPORT_BINDING_SOURCE','S9b_SHARED_PHYSICS.md:217; asymptotic bulk sound speed')

    dj,vj,hj=Jet({(1,0,0):delta}),Jet({(0,1,0):velocity}),Jet({(0,0,1):slope**2})
    mapping={chi:(c0*(1+dj))**2,U:vj,A:1+hj}
    arrj=evaluate_jet(arr,mapping).tidy()
    aphij=evaluate_jet(aphi,mapping).tidy()
    betaj=evaluate_jet(beta,mapping).tidy()
    optical_r=(aphij**sp.Rational(1,2)*c0).tidy()
    optical_q=((arrj**sp.Rational(1,2)*c0)/optical_r.diff(r)).tidy()
    inverse=Jet(x)
    inverse_factor=(Jet(r)/optical_r).tidy()
    for _ in range(3):
        inverse=compose_radial(inverse_factor,r,x,inverse)*x
    optical_q=compose_radial(optical_q,r,x,inverse).tidy()
    emit('OPTICAL_RADIUS',optical_r.c)
    emit('OPTICAL_Q',optical_q.c)
    turning=sp.factor(geometry['even_squared']/arr).subs(profile_map).subs(J,B/c0)
    turning_domain=sp.And(eq(turning.subs(r,rt),0),sp.Gt(radial_diff(turning,r).subs(r,rt),0),
                          everywhere(r,sp.Gt(turning,0),sp.Interval.open(rt,sp.oo)))
    observable_gate = real_domain  # K10: existence only
    gate=sp.And(observable_gate,travel_domain,turning_domain,sp.Lt(rt,sp.Min(RE,RR))).subs(endpoint_map).subs(B,b)
    path=nonreciprocity(r,c0,delta,velocity,xi)
    emit('NONRECIPROCITY_ONEFORM',path['oneform'])
    emit('NONRECIPROCITY_EXTERIOR_DERIVATIVE',path['curl'])
    emit('NONRECIPROCITY_PATH_PREDICATE',path['predicate'])
    emit('A_NONRECIPROCAL_PATH_DEPENDENCE',gated(path['predicate'],gate))

    kernel=sp.sqrt(1-B*B/x**2)
    angular_kernel=-2*sp.diff(kernel,B)
    theta={}
    action=Jet()
    for g in GRADES:
        perturbation=optical_q[g]-(1 if g==ZERO else 0)
        theta[g]=sp.Integral(perturbation*angular_kernel,(x,B,sp.oo)) if perturbation!=0 else sp.S.Zero
        if perturbation!=0:
            action+=Jet({g:sum(sp.Integral(perturbation*kernel/c0,(x,B,R)) for R in (RE,RR))})
    full_integrand=optical_q*kernel/c0
    for R in (RE,RR):
        shift=optical_r.map(lambda e:e.subs(r,R))-R
        for k in range(1,4):
            action+=full_integrand.diff(x,k-1).map(lambda e:e.subs(x,R))*shift**k/sp.factorial(k)
    # The flat action comes from the same radial kernel after regularization.
    t=sp.Symbol('s9b_tangent',positive=True)
    regular=B*sp.sqrt(1+t*t)
    flat_primitive=sp.integrate(sp.simplify(kernel.subs(x,regular)*sp.diff(regular,t)),t)
    flat_action=sp.simplify(sum((flat_primitive.subs(t,sp.sqrt(R*R/B**2-1))-flat_primitive.subs(t,0))/c0 for R in (RE,RR)))
    H0=sp.simplify(sp.diff(flat_action,B,2))
    H1=sp.simplify(sp.diff(H0,B))
    # Differentiate complete action functionals, without splitting singular
    # moving-lower-limit boundary terms of their higher derivatives.
    def action_diff(jet,order):
        return jet.map(lambda e:sp.Derivative(e,B,order,evaluate=False))
    first_shift=-action_diff(action,1)/H0
    shift=-(action_diff(action,1)+action_diff(action,2)*first_shift+H1*first_shift**2/2)/H0
    even_time=action+action_diff(action,1)*shift+action_diff(action,2)*shift**2/2+H0*shift**2/2+H1*shift**3/6
    even_time=even_time.map(lambda e:e.subs(B,b))
    odd_time=Jet({g:sp.Integral(value,(r,RE,RR)) for g,value in betaj.c.items()})
    return_orientation = -1  # K5
    forward_time=even_time+odd_time
    reverse_time=even_time+return_orientation*odd_time
    roundtrip_time=forward_time+reverse_time
    nonreciprocal_time=(forward_time-reverse_time)/2
    first_grades=((1,0,0),(0,1,0),(0,2,0),(0,0,1))
    log_source = roundtrip_time  # K6
    radar={}
    for g in first_grades:
        source=log_source[g].subs(endpoint_map)
        derivative=sp.expand(-b*radial_diff(source,b)/2)
        terms=sp.Add.make_args(derivative)
        bulk=sum((term for term in terms if term.has(sp.Integral)),sp.S.Zero)
        remainder=sp.simplify(derivative-bulk)
        # A counting certificate, not a profile ansatz: endpoint values are
        # bounded by independent constants times the supplied radial powers.
        # The triangle inequality makes this independent of profile signs.
        envelope_d,envelope_v,envelope_h=sp.symbols('s9b_bound_delta s9b_bound_V s9b_bound_xi',positive=True)
        envelope=profile_substitute(remainder,{
            Delta:lambda z:envelope_d/z,
            Velocity:lambda z:envelope_v/sp.sqrt(z),
            Embedding:lambda z:2*envelope_h*sp.sqrt(z)})
        bound=sum((sp.Abs(term) for term in sp.Add.make_args(sp.expand(envelope))),sp.S.Zero)
        endpoint_limit=sp.limit(sp.limit(bound,ZE,sp.oo),ZR,sp.oo)
        radar[g]=bulk.replace(lambda a:isinstance(a,sp.Integral),
            lambda a:sp.Integral(a.function,*(tuple(limit[:-1])+(sp.oo,) for limit in a.limits)))
        emit('RADAR_SLOPE_SOURCE_G%d%d%d'%g,source)
        emit('RADAR_FINITE_SLOPE_G%d%d%d'%g,derivative)
        emit('RADAR_ENDPOINT_REMAINDER_G%d%d%d'%g,remainder)
        emit('RADAR_ENDPOINT_COUNTING_CERTIFICATE_G%d%d%d'%g,(bound,endpoint_limit))
    reference_gamma=sp.Symbol('s9b_gamma',real=True)
    reference_theta=(1+reference_gamma)*2*GM/(b*c0**2)
    reference_time=2*(1+reference_gamma)*GM/c0**3*sp.log(4*RE*RR/b**2)
    reference_slope=sp.limit(sp.limit(-b*sp.diff(reference_time.subs(endpoint_map),b)/2,ZE,sp.oo),ZR,sp.oo)
    deflection=sum((theta[g].subs(B,b) for g in first_grades),sp.S.Zero)
    radar_slope=sum(radar.values(),sp.S.Zero)
    references={'DEFLECTION':reference_theta,'RADAR':reference_slope}
    comparisons={'DEFLECTION':deflection,'RADAR':radar_slope}
    gammas={name:sp.solve(eq(references[name],obj),reference_gamma)[0] for name,obj in comparisons.items()}
    residuals={name:obj-references[name].subs(reference_gamma,1) for name,obj in comparisons.items()}
    emit('REFERENCE_OBJECTS',(reference_theta,reference_time,reference_slope))
    qsum=sum((optical_q[g] for g in first_grades),sp.S.Zero)
    def integrand_of(value):
        total=sp.S.Zero
        for term in sp.Add.make_args(sp.expand(value)):
            integrals=term.atoms(sp.Integral)
            if not integrals:
                if term!=0: raise ValueError(('nonintegral slope remainder',term))
                continue
            if len(integrals)!=1: raise ValueError('nonlinear comparison functional')
            integral=next(iter(integrals))
            total+=term/integral*integral.function
        return sp.factor(total)
    slope_density=integrand_of(radar_slope)
    slope_ratio=sp.simplify(slope_density/(qsum*angular_kernel.subs(B,b)))
    emit('RADAR_ANGULAR_OPERATOR_RATIO',slope_ratio)
    # Abel inverse: composing the two square-root kernels gives this beta
    # integral. Decay anchoring removes the integration constant.
    eta=sp.Symbol('s9b_abel_eta',positive=True)
    emit('ABEL_KERNEL_COMPOSITION',sp.integrate(1/sp.sqrt(eta*(1-eta)),(eta,0,1)))
    inverse_targets={}
    for name,ratio in (('DEFLECTION',sp.S.One),('RADAR',slope_ratio)):
        target=references[name].subs(reference_gamma,1)/ratio
        transformed=sp.diff(target/b,b)/sp.sqrt(b*b-r*r)
        substitution=r*sp.sqrt(1+t*t)
        integral=sp.integrate(sp.simplify(transformed.subs(b,substitution)*sp.diff(substitution,t)),(t,0,sp.oo))
        inverse_targets[name]=sp.simplify(-r*r*integral/sp.pi)
    emit('ABEL_REFERENCE_INVERSES',inverse_targets)
    qradial=qsum.subs(x,r)
    condition_mass={name:sp.solve(eq(qradial,target),GM)[0] for name,target in inverse_targets.items()}
    fixed_impact = None  # K7
    if fixed_impact is not None:
        condition_mass={name:sp.solve(eq(comparisons[name].subs(b,fixed_impact),references[name].subs({b:fixed_impact,reference_gamma:1})),GM)[0] for name in comparisons}
    # Supplied coordinate mass measure, computed through the spherical map.
    polar,azimuth=sp.symbols('s9b_polar s9b_azimuth',real=True)
    coords=sp.Matrix([r*sp.sin(polar)*sp.cos(azimuth),r*sp.sin(polar)*sp.sin(azimuth),r*sp.cos(polar)])
    jacobian=sp.trigsimp(coords.jacobian((r,polar,azimuth)).det())
    area=sp.integrate(jacobian,(polar,0,sp.pi),(azimuth,0,2*sp.pi))
    mass_density = rho  # K8
    exchange=sp.simplify(-sp.diff(area*mass_density*velocity,r)/area)
    emit('MASS_BALANCE',eq(jn,exchange))
    emit('COORDINATE_SPHERE_MEASURE',area)
    objects={}
    def substitute(value,replacement):
        return profile_substitute(value,replacement)
    def emit_condition(prefix,replacement,own_domain=sp.true,implied=True):
        domain=sp.And(substitute(gate,replacement),own_domain)
        for name,mass in condition_mass.items():
            local_mass=substitute(mass,replacement)
            condition=cas(('RADIAL_DIFFERENTIAL_IDENTITY',r,sp.Interval.open(bmin,sp.oo),domain,eq(GM,local_mass)))
            tag=prefix+'CONDITION_'+name
            emit(tag,condition)
            objects['s9b_'+tag.lower()]=condition
            if implied:
                tag=prefix+'IMPLIED_JN_'+name
                value=cas((condition,eq(jn,substitute(exchange,replacement))))
                emit(tag,value)
                objects['s9b_'+tag.lower()]=value
    def emit_comparisons(prefix,replacement,domain):
        for name in comparisons:
            emit(prefix+'B_GAMMA_'+name,gated(substitute(gammas[name],replacement),sp.And(domain,sp.Ne(GM,0))))
            emit(prefix+'B_RESIDUAL_'+name,gated(substitute(residuals[name],replacement),domain))
        emit(prefix+'B_GAMMA_DIFFERENCE',gated(substitute(gammas['DEFLECTION']-gammas['RADAR'],replacement),sp.And(domain,sp.Ne(GM,0))))
    def emit_a(prefix,replacement,domain,all_times=True):
        for g in GRADES:
            suffix='_G%d%d%d'%g
            emit(prefix+'A_DEFLECTION'+suffix,gated(substitute(theta[g].subs(B,b),replacement),domain))
            if all_times:
                for name,jet in (('ROUND_TRIP',roundtrip_time),('ONE_WAY_ER',forward_time),('ONE_WAY_RE',reverse_time),('NONRECIPROCAL',nonreciprocal_time)):
                    emit(prefix+'A_'+name+suffix,gated(substitute(jet[g].subs(endpoint_map),replacement),domain))
        for g,value in radar.items():
            emit(prefix+'A_RADAR_LOG_G%d%d%d'%g,gated(substitute(value,replacement),domain))
    emit_a('',{},gate)
    emit_comparisons('',{},gate)
    emit_condition('B_',{})
    restrictions={'FLOW_ONLY':{Delta:lambda x:sp.S.Zero,Embedding:lambda x:sp.S.Zero},
                  'SPEED_ONLY':{Velocity:lambda x:sp.S.Zero,Embedding:lambda x:sp.S.Zero},
                  'TILT_ONLY':{Delta:lambda x:sp.S.Zero,Velocity:lambda x:sp.S.Zero}}
    for label,replacement in restrictions.items():
        prefix='R_'+label+'_'
        domain=substitute(gate,replacement)
        emit_a(prefix,replacement,domain,False)
        emit_comparisons(prefix,replacement,domain)
        emit_condition(prefix+'B_',replacement)
    f_symbol=sp.Symbol('s9b_fractional_bulk_change',real=True)
    cs_ratio=sp.powsimp(sp.powdenest(cs_squared.subs(rho_bulk,rho0*(1+f_symbol))/cs0_squared,force=True),force=True)
    fixed_ratio_response = sp.sqrt(cs_ratio)  # K9a
    power_response = (1+f_symbol)**s  # K9b
    responses=(sp.S.One,fixed_ratio_response,power_response)
    response_rows=[]
    for label,response in zip(('CONSTANT','FIXED_RATIO','POWER'),responses):
        coefficient=sp.simplify(sp.diff(response,f_symbol).subs(f_symbol,0))
        replacement={Delta:lambda x,a=coefficient:a*BulkChange(x)}
        local=response.subs(f_symbol,f)
        domain=sp.And(sp.Contains(local,sp.S.Reals),sp.Gt(local,0))
        if label!='CONSTANT': domain=sp.And(domain,sp.Gt(rho0*(1+f),0))
        if label=='FIXED_RATIO': domain=sp.And(domain,sp.Gt(cs0_squared,0),sp.Gt(cs_squared.subs(rho_bulk,rho0*(1+f)),0))
        emit('C_'+label+'_RESPONSE',(response,coefficient*f,domain))
        emit_condition('C_'+label+'_',replacement,domain)
        emit('C_'+label+'_N_DEPENDENCE',(tuple(sp.diff(substitute(value,replacement),n) for value in condition_mass.values()),domain))
        response_rows.append((label,replacement,domain))
        if label!='CONSTANT':
            emit_condition('R_BULK_ONLY_C_'+label+'_',replacement|{Velocity:lambda x:sp.S.Zero,Embedding:lambda x:sp.S.Zero},domain)
    Phi=sp.Symbol('Phi',real=True)
    liveV=sp.Function('s9b_mass_V')(r)
    solution=sp.dsolve(eq(sp.diff(area*mass_density*liveV,r),0),liveV)
    constant=next(a for a in solution.free_symbols if a.name=='C1')
    flux=sp.simplify(area*mass_density*solution.rhs)
    flux_constant=sp.solve(eq(Phi,flux),constant)[0]
    forwardV=solution.rhs.subs(constant,flux_constant)
    emit('F_MASS_SOLUTION',(eq(jn,0),eq(Phi,flux),eq(liveV,forwardV)))
    for label,extra in (('FLOW',restrictions['FLOW_ONLY']),('LIVE',{})):
        replacement={Velocity:lambda x:forwardV.subs(r,x)}|extra
        prefix='F_'+label+'_'
        domain=substitute(gate,replacement)
        emit(prefix+'BRANCH_EXISTENCE',substitute(real_domain,replacement))
        emit(prefix+'PATH_TRAVERSAL',substitute(travel_domain,replacement))
        emit(prefix+'BRANCH_TYPE',substitute(cas(branch_type).subs(profile_map),replacement))
        emit_a(prefix,replacement,domain)
        emit_comparisons(prefix,replacement,domain)
        emit_condition(prefix+'B_',replacement,implied=False)
        if label=='LIVE':
            for response_label,response_map,response_domain in response_rows:
                emit_condition(prefix+'C_'+response_label+'_',replacement|response_map,response_domain,False)
    labels={tag:('FORWARD_NO_FAR_ZONE_LOSS',eq(jn,0),'Phi is the live outward coordinate-sphere mass flux')
            for tag in TAGS if tag.startswith('PY_S9B_F_')}
    for label,replacement in restrictions.items():
        for tag in TAGS:
            if tag.startswith('PY_S9B_R_'+label+'_'):
                labels[tag]=(label,tuple(eq(profile(r),value(r)) for profile,value in replacement.items()))
    for tag in TAGS:
        if tag.startswith('PY_S9B_R_BULK_ONLY_'):
            labels[tag]=('BULK_ONLY',eq(velocity,0),eq(xi,0))
    emit('OBJECT_LABELS',labels)
    emitted={tag.removeprefix('PY_S9B_') for tag in TAGS if tag.startswith('PY_S9B_')}
    emit('SHARED_VOCABULARY',(sorted(SHARED),sorted(emitted),sorted(SHARED-emitted)))
    if emitted!=SHARED: raise ValueError(('shared vocabulary mismatch',SHARED-emitted,emitted-SHARED))
    declarations={}
    for value in objects.values():
        for symbol in value.atoms(sp.Symbol):
            declarations[symbol.name]=(symbol,'COORDINATE' if symbol in (r,x,b,B,rt,bmin) else 'KNOB')
        for function in value.atoms(AppliedUndef):
            declarations[function.func.__name__]=(function.func,'PREMISE')
    publish(fold,objects,declarations)
    emit('LOCAL_NAMES',LOCAL+['PY_LOCAL_S9B_LOCAL_NAMES'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--branch-demo',choices=('elliptic','degenerate','blocked'))
    build(parser.parse_args().branch_demo)
