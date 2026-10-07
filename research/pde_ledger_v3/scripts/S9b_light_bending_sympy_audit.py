#!/usr/bin/env python3
"""S9b SymPy builder. Sole physics input: S9b_SHARED_PHYSICS.md.

1. The script may PRINT computed objects. It may NOT state conclusions.
2. PRINT the residual; do NOT assert it. Compute -> emit -> then guard.
3. Interpretation belongs to the STEP RECORD.

The ONLY place the physical symbols may be combined by hand is in
CONSTRUCTING THE ACTION and the ANSATZ. Every other expression involving them
must be REACHED BY COMPUTATION. Every control re-enters the chain at the ACTION,
never at a result.

Profile family (not a general-profile claim): delta=d*(L/r)**p,
V^r=c0*v*(L/r)**q, xi'=w*(L/r)**h, f=F*(L/r)**p.
The density rho_br(r) is arbitrary. V denotes the coordinate velocity dr/dt,
as in the supplied contravariant advection term. All amplitudes and exponents
remain symbolic. The asymptotically anchored, simple outer turning branch is
used, conditional on the printed existence/traversal domain.

Method: solve the radial Hamilton-Jacobi quadratic; change to optical
circumference radius X=c0*sqrt(a_phi_phi). Formal inversion in the finite
multigraded ring keeps the turning point at X=B. Integrate its regular action
kernel. For finite endpoints, vary B so that the angular separation remains
fixed; Taylor stationary elimination retains the bent-path mixed terms.
The odd radial one-form is integrated between endpoints, before ray variation.
No singular moving-lower-limit Taylor integrals are used.

Only writes: stdout and S9b_exports.py. Run with python -B through the guard.
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

IMPORT_KEYS = ('rho_br', 'c_s0', 'transverse_dispersion')
INPUTS = (HERE, ROOT/'directives/S9b_SHARED_PHYSICS.md',
          *(ROOT/('scripts/'+s+'_exports.py') for s in
            ('S11c_b', 'S11c_c1', 'S11c_c2')), ROOT/'scripts/ledger_fold.py')
GRADES = tuple(product(range(2), range(3), range(2)))
ZERO = (0, 0, 0)
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
    tag = ('PY_LOCAL_S9B_' if local else 'PY_S9B_') + name
    if tag in TAGS:
        raise ValueError(('duplicate tag', tag))
    obj = cas(value)
    TAGS[tag] = obj
    if local:
        LOCAL.append(tag)
    print(tag+': '+sp.srepr(obj), flush=True)
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
        return self.map(lambda v: sp.diff(v, variable, order))

    def tidy(self):
        return self.map(lambda v: sp.factor(sp.expand(v)))


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
                out += term*sp.diff(value, old, k).subs(old, new)/sp.factorial(k)
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
            'class': class_tag, 'step': 'S9b'}


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
            outcome = compare(prior, record['value'])
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
    lines = ['# Generated S9b own-rows delta. Conditional symbolic profile family.',
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


def build(branch_demo=None):
    # SUPPLIED dispersion, profiles, balance, bulk responses and references.
    r, x, B, b, L, c0, ZE, ZR = sp.symbols('r X B b L c_0 Z_E Z_R', positive=True)
    RE, RR, rt, bmin = sp.symbols('r_E r_R r_turn b_min', positive=True)
    d, v, w, F = sp.symbols('d v w F', real=True)
    p, q, h = sp.symbols('p q h', positive=True)
    GM, rho0, K, m = sp.symbols('GM rho_0 K m', positive=True)
    n, s = sp.symbols('n s', real=True)
    kr, J, omega, U, chi = sp.symbols('k_r J omega U chi', real=True)
    A = sp.Symbol('A', positive=True)
    mu = sp.Function('mu_perp', real=True)(r)
    rho = sp.Function('rho_br_profile', positive=True)(r)
    jn = sp.Function('j_n', real=True)(r)
    delta, velocity, slope, f = d*(L/r)**p, c0*v*(L/r)**q, w*(L/r)**h, F*(L/r)**p
    metric = sp.diag(1+slope**2, r**2, r**2*sp.sin(sp.Symbol('theta', real=True))**2)
    dispersion = (omega-U*kr)**2-chi*(kr**2/A+(omega*J)**2/r**2)
    speed_relation = eq(chi, mu/rho)
    reference_theta = (1+sp.Symbol('gamma', real=True))*2*GM/(b*c0**2)
    reference_log = 2*(1+sp.Symbol('gamma', real=True))*GM/c0**3
    responses = (sp.S.One, None, (1+sp.Symbol('f', real=True))**s)
    bulk_rho = sp.Symbol('rho_bulk', positive=True)
    bulk_pressure = K*bulk_rho**n
    cs_squared = sp.diff(bulk_pressure, bulk_rho)/m

    emit('SUPPLIED_DEPENDENCIES', ('DISPERSION', 'SPEED_IDENTIFICATION', 'ADVECTION',
         'INDUCED_METRIC', 'MASS_BALANCE', 'LAB_HELD', 'ORDER_COUNTING', 'GR_REFERENCES', 'PART_C_RESPONSES'))
    emit('PROFILE_FAMILY', (eq(sp.Function('delta')(r), delta), eq(sp.Function('V')(r), velocity),
         eq(sp.Derivative(sp.Function('xi_w')(r), r), slope), eq(sp.Function('f')(r), f)))
    emit('PROFILE_DOMAIN', sp.And(sp.Ge(p, 1), sp.Ge(q, sp.Rational(1, 2)), sp.Ge(h, sp.Rational(1, 2))))
    emit('ENDPOINTS', (eq(RE, sp.sqrt(b*b+ZE*ZE)), eq(RR, sp.sqrt(b*b+ZR*ZR))))
    emit('DISPERSION', eq(dispersion, 0))
    emit('SPEED_RELATION', speed_relation)

    polynomial = sp.Poly(dispersion, kr)
    aa, bb, cc = polynomial.all_coeffs()
    beta = sp.factor(-bb/(2*aa*omega))
    even_squared = sp.factor((bb**2-4*aa*cc)/(4*aa**2*omega**2))
    arr = sp.factor(even_squared.subs(J, 0))
    # Coefficient extraction is independent of rational-expression rendering.
    aphi = sp.factor(-arr/(sp.diff(even_squared, J, 2)/2))
    radical = sp.factor(even_squared/arr)
    radial_omega = sp.solve(dispersion.subs(J, 0), omega)
    radial_speeds = tuple(sp.diff(root, kr) for root in radial_omega)
    existence = sp.And(sp.Gt(chi, 0), sp.Gt(A, 0), sp.Ne(rho, 0))
    traversal = sp.Gt(chi/A-U**2, 0)
    emit('RADIAL_FREQUENCIES', radial_omega)
    emit('RADIAL_GROUP_VELOCITIES', radial_speeds)
    emit('BRANCH_EXISTENCE', existence)
    emit('PATH_TRAVERSAL', traversal)
    emit('BRANCH_TYPES', ((sp.Lt(chi, 0), ('GROWING', 'DECAYING')),
         (eq(chi, 0), 'ABSENT'), (sp.And(existence, sp.Not(traversal)), 'UNABLE_TO_TRAVERSE_REQUIRED_DIRECTION')))
    emit('OBSERVABLE_DOMAIN', ((sp.And(existence, traversal), 'CONDITIONAL_RAY_FUNCTIONALS'),
         (sp.Not(sp.And(existence, traversal)), ('NOT_ESTABLISHED', 'REAL_PROPAGATING_BIDIRECTIONAL_BRANCH_REQUIRED'))))
    if branch_demo is not None:
        supplied = {'elliptic': {chi: -c0**2, U: 0, A: 1},
                    'degenerate': {chi: 0, U: 0, A: 1},
                    'blocked': {chi: c0**2, U: -2*c0, A: 1}}[branch_demo]
        emit('DEMONSTRATION_INPUT', supplied, local=True)
        emit('DEMONSTRATION_BRANCH', (existence.subs(supplied), traversal.subs(supplied),
             tuple(z.subs(supplied) for z in radial_omega)), local=True)
        emit('OBSERVABLES', ('NOT_ESTABLISHED', 'REAL_PROPAGATING_BIDIRECTIONAL_BRANCH_REQUIRED'))
        emit('LOCAL_NAMES', LOCAL)
        return

    fold, fold_audit = load_model(*(ROOT/('scripts/'+s+'_exports.py') for s in ('S11c_b', 'S11c_c1', 'S11c_c2')))
    manifest = check_consumer(fold, IMPORT_KEYS)
    bound = assert_lookups_equal_manifest(lambda ledger: {k: ledger[k]['value'] for k in IMPORT_KEYS}, fold, IMPORT_KEYS)
    inputs = bound['result']
    anchor = inputs['transverse_dispersion'][0]
    anchor_poly = anchor.lhs-anchor.rhs
    anchor_k = next(a for a in anchor_poly.free_symbols if a.name == 'k')
    stiffness = sp.diff(anchor_poly, anchor_k, 2)/2
    emit('ANCHOR', (anchor, eq(c0**2, stiffness/inputs['rho_br'])))
    emit('IMPORT_MANIFEST', (IMPORT_KEYS, sorted(bound['lookups']), sorted(manifest['closure'])), local=True)
    emit('OPTICAL_RADIAL_ACTION', (beta, even_squared, arr, aphi))

    dj = Jet({(1, 0, 0): delta})
    vj = Jet({(0, 1, 0): velocity})
    hj = Jet({(0, 0, 1): slope**2})
    mapping = {chi: (c0*(1+dj))**2, U: vj, A: 1+hj}
    arrj = evaluate_jet(arr, mapping).tidy()
    aphij = evaluate_jet(aphi, mapping).tidy()
    betaj = evaluate_jet(beta, mapping).tidy()
    optical_r = (aphij**sp.Rational(1, 2)*c0).tidy()
    optical_q = ((arrj**sp.Rational(1, 2)*c0)/optical_r.diff(r)).tidy()
    inverse = Jet(x)
    inverse_factor = (Jet(r)/optical_r).tidy()
    # Three nonzero even generators => nilpotence index four.
    for _ in range(3):
        inverse = compose_radial(inverse_factor, r, x, inverse)*x
    optical_q = compose_radial(optical_q, r, x, inverse).tidy()
    emit('OPTICAL_RADIUS', optical_r.c)
    emit('INVERSE_OPTICAL_RADIUS', inverse.c)
    emit('OPTICAL_RADIAL_COEFFICIENT', optical_q.c)

    profile_map = {chi: (c0*(1+delta))**2, U: velocity, A: metric[0, 0]}
    real_domain = everywhere(r, existence.subs(profile_map), sp.Interval(rt, sp.oo))
    travel_domain = everywhere(r, traversal.subs(profile_map), sp.Interval(rt, sp.oo))
    turning = radical.subs(profile_map).subs(J, B/c0)
    turning_domain = sp.And(eq(turning.subs(r, rt), 0), sp.Gt(sp.diff(turning, r).subs(r, rt), 0),
                           everywhere(r, sp.Gt(turning, 0), sp.Interval.open(rt, sp.oo)))
    emit('RAY_DOMAIN', (real_domain, travel_domain, turning_domain, sp.Lt(rt, sp.Min(RE, RR))))

    # Integrals after y=1-B**2/X**2. betainc uses a finite upper endpoint,
    # including second parameter zero; no Gamma poles are introduced there.
    def action_primitive(power, upper):
        return B**(1-power)*sp.betainc(sp.Rational(3, 2), (power-1)/2, 0, 1-B**2/upper**2)/2

    def monomial(value):
        # The Euler derivative extracts the degree even when SymPy groups
        # positive radial ratios, e.g. (L**2/X**2)**h, into one power base.
        power = sp.simplify(-x*sp.diff(value, x)/value)
        coefficient = sp.simplify(sp.expand_power_base(value, force=True)*x**power)
        if coefficient.has(x):
            raise ValueError(('radial coefficient is not a power', value))
        return coefficient, power

    qdata = {g: monomial(value) for g, value in optical_q.c.items()}
    theta = {}
    action = Jet()
    kernel = sp.sqrt(1-B**2/x**2)
    for g in GRADES:
        coefficient, power = qdata.get(g, (sp.S.Zero, sp.S.Zero))
        # The deflection is -dW/dJ-pi; this convergent angular integral also
        # covers powers <=1, where the infinite radial action itself diverges.
        angular = coefficient*B**(-power)*sp.expand_func(sp.beta((power+1)/2, sp.Rational(1, 2)))
        theta[g] = sp.simplify(angular-(sp.pi if g == ZERO else 0))
        if g != ZERO and coefficient != 0:
            action += Jet({g: coefficient*sum(action_primitive(power, R) for R in (RE, RR))/c0})
    # Flat action is obtained by integrating the same radial kernel.
    tangent = sp.Symbol('tangent', positive=True)
    radial_substitution = B*sp.sqrt(1+tangent**2)
    transformed = sp.simplify(kernel.subs(x, radial_substitution)*sp.diff(radial_substitution, tangent))
    flat_primitive = sp.integrate(transformed, tangent)
    flat_action = sp.simplify(sum((flat_primitive.subs(tangent, sp.sqrt(R**2/B**2-1))-
                                  flat_primitive.subs(tangent, 0))/c0 for R in (RE, RR)))
    # Moving optical endpoint X(R), expanded at fixed physical R.
    full_integrand = optical_q*kernel/c0
    for R in (RE, RR):
        shift = optical_r.map(lambda e: e.subs(r, R))-R
        for k in range(1, 4):
            boundary = full_integrand.diff(x, k-1).map(lambda e: e.subs(x, R))
            action += boundary*shift**k/sp.factorial(k)
    action = action.tidy()
    emit('RADIAL_ACTION_COEFFICIENTS', action.c)
    emit('FLAT_RADIAL_ACTION', flat_action)
    angle = sp.acos(B/RE)+sp.acos(B/RR)
    emit('ENDPOINT_ANGLE', angle.subs(B, b))
    H0 = sp.simplify(sp.diff(flat_action, B, 2))
    H1 = sp.simplify(sp.diff(H0, B))
    first_shift = -action.diff(B)/H0
    shift = -(action.diff(B)+action.diff(B, 2)*first_shift+H1*first_shift**2/2)/H0
    even_time = action+action.diff(B)*shift+action.diff(B, 2)*shift**2/2+H0*shift**2/2+H1*shift**3/6
    # Differentiate before substituting B=b; physical endpoints stay fixed.
    even_time = even_time.map(lambda e: e.subs(B, b))
    emit('ENDPOINT_IMPACT_SHIFT', shift.map(lambda e: e.subs(B, b)).c)

    odd_time = {}
    for g in GRADES:
        integrand = betaj[g]
        if integrand == 0:
            odd_time[g] = sp.S.Zero
        else:
            coef, power = monomial(integrand.subs(r, x))
            odd_time[g] = coef*sp.Piecewise((sp.log(RR/RE), eq(power, 1)),
                             ((RR**(1-power)-RE**(1-power))/(1-power), True))
    # Exactness of the radial one-form is computed in Cartesian coordinates.
    xyz = sp.symbols('x_1 x_2 x_3', real=True)
    radius = sp.sqrt(sum(a*a for a in xyz))
    beta_general = sp.Function('beta_r')(radius)
    oneform = [beta_general*sp.diff(radius, a) for a in xyz]
    curl = sp.ImmutableMatrix([sp.simplify(sp.diff(oneform[(i+2)%3], xyz[(i+1)%3])-
                   sp.diff(oneform[(i+1)%3], xyz[(i+2)%3])) for i in range(3)])
    emit('NONRECIPROCITY_EXTERIOR_DERIVATIVE', curl)
    emit('NONRECIPROCITY_PATH_DEPENDENCE', sp.Eq(curl, sp.zeros(3, 1)))
    emit('NONRECIPROCITY_PRIMITIVE', sp.Integral(beta.subs(profile_map), r))
    endpoint_map = {RE: sp.sqrt(b*b+ZE*ZE), RR: sp.sqrt(b*b+ZR*ZR)}
    first_grades = ((1, 0, 0), (0, 1, 0), (0, 2, 0), (0, 0, 1))
    for g in GRADES:
        suffix = 'D%d_V%d_H%d' % g
        emit('ORDER_'+suffix, sp.Rational(g[0])+sp.Rational(g[1], 2)+g[2])
        emit('DEFLECTION_'+suffix, theta[g].subs(B, b))
        emit('TIME_FORWARD_'+suffix, (even_time[g]+odd_time[g]).subs(endpoint_map))
        emit('TIME_REVERSE_'+suffix, (even_time[g]-odd_time[g]).subs(endpoint_map))
        emit('TIME_RADAR_'+suffix, (2*even_time[g]).subs(endpoint_map))
        emit('TIME_NONRECIPROCAL_'+suffix, odd_time[g].subs(endpoint_map))

    # Coefficient of ln(1/b**2) in the large-endpoint first-order radar time.
    # The beta primitive's large-X integrand is coefficient*X**(-power)/c0.
    # Endpoint motion has only powers (no log); stationary path corrections
    # start at higher order. Two endpoints and two legs give the conversion.
    log_terms = {}
    for g in first_grades:
        coef, power = qdata.get(g, (sp.S.Zero, sp.S.Zero))
        log_terms[g] = 2*coef/c0*sp.Piecewise((1, eq(power, 1)), (0, True))
        emit('RADAR_LOG_D%d_V%d_H%d' % g, log_terms[g])
    deflection = sum(theta[g].subs(B, b) for g in first_grades)
    radar_log = sum(log_terms.values())
    gamma = next(a for a in reference_theta.free_symbols if a.name == 'gamma')
    gamma_theta = sp.solve(eq(reference_theta, deflection), gamma)[0]
    gamma_radar = sp.solve(eq(reference_log, radar_log), gamma)[0]
    residual_theta = deflection-reference_theta.subs(gamma, 1)
    residual_radar = radar_log-reference_log.subs(gamma, 1)
    emit('REFERENCES', (reference_theta, reference_log, eq(gamma, 1)))
    emit('GAMMA_DEFLECTION', gamma_theta)
    emit('GAMMA_RADAR', gamma_radar)
    emit('GAMMA_DIFFERENCE', gamma_theta-gamma_radar)
    emit('RESIDUAL_DEFLECTION', (deflection, reference_theta.subs(gamma, 1), residual_theta))
    emit('RESIDUAL_RADAR', (radar_log, reference_log.subs(gamma, 1), residual_radar))

    # Divergence from the induced volume element, with rho_br and V^r live.
    volume = sp.sqrt(metric.det())
    exchange = sp.simplify(-sp.diff(volume*rho*velocity, r)/volume)
    flat_volume = volume.subs(w, 0)
    flat_exchange = sp.simplify(-sp.diff(flat_volume*rho*velocity, r)/flat_volume)
    emit('MASS_BALANCE', eq(jn, exchange))
    emit('MASS_BALANCE_FLAT', eq(jn, flat_exchange))
    emit('MASS_BALANCE_METRIC_CORRECTION', sp.factor(exchange-flat_exchange))
    emit('MASS_BALANCE_VOLUME', volume)
    mu_profile = sp.solve(speed_relation.subs(chi, (c0*(1+delta))**2), mu)[0]
    emit('STIFFNESS_PROFILE', eq(mu, mu_profile))
    emit('DENSITY_ANCHOR', eq(sp.Limit(rho, r, sp.oo), inputs['rho_br']))

    cs0_squared = cs_squared.subs(bulk_rho, rho0)
    f_symbol = sp.Symbol('f', real=True)
    cs_ratio = sp.powsimp(sp.powdenest(cs_squared.subs(bulk_rho, rho0*(1+f_symbol))/cs0_squared, force=True), force=True)
    responses = (responses[0], sp.sqrt(cs_ratio), responses[2])
    response_delta = [sp.diff(response, f_symbol).subs(f_symbol, 0)*f for response in responses]
    emit('BULK_PRESSURE', eq(sp.Function('P')(bulk_rho), bulk_pressure))
    emit('BULK_SOUND_SPEED', cs_squared)
    emit('BULK_SOUND_ANCHOR', eq(inputs['c_s0']**2, cs0_squared))
    emit('BULK_RESPONSE_DELTAS', response_delta)
    objects = {}
    for name, res in (('deflection', residual_theta), ('radar', residual_radar)):
        condition = everywhere(b, eq(res, 0), sp.Interval.open(bmin, sp.oo))
        objects['light_'+name+'_condition'] = condition
        objects['light_'+name+'_exchange'] = sp.And(condition, eq(jn, exchange))
        emit('B_'+name.upper()+'_CONDITION', condition)
        emit('B_'+name.upper()+'_EXCHANGE', objects['light_'+name+'_exchange'])
        for label, response in zip(('constant', 'fixed_ratio', 'power'), response_delta):
            amplitude = sp.simplify(response/(L/r)**p)
            changed = res.subs(d, amplitude)
            ccondition = everywhere(b, eq(changed, 0), sp.Interval.open(bmin, sp.oo))
            key = 'light_'+label+'_'+name
            objects[key+'_condition'] = ccondition
            objects[key+'_exchange'] = sp.And(ccondition, eq(jn, exchange))
            emit('C_'+label.upper()+'_'+name.upper()+'_CONDITION', ccondition)
            emit('C_'+label.upper()+'_'+name.upper()+'_EXCHANGE', objects[key+'_exchange'])
            emit('C_'+label.upper()+'_'+name.upper()+'_N_DEPENDENCE', sp.diff(changed, n))
    # Explicit conditions carry their quantifiers; applicability is exported as
    # their conjunction, so a later consumer cannot lose the family/ray domain.
    applicability = sp.And(real_domain, travel_domain, turning_domain,
                           sp.Ge(p, 1), sp.Ge(q, sp.Rational(1, 2)), sp.Ge(h, sp.Rational(1, 2)))
    objects = {k: sp.And(applicability, value) for k, value in objects.items()}
    declarations = {}
    for value in objects.values():
        for symbol in value.atoms(sp.Symbol):
            declarations[symbol.name] = (symbol, 'COORDINATE' if symbol in (r, B, b, rt, bmin) else 'KNOB')
        for function in value.atoms(AppliedUndef):
            declarations[function.func.__name__] = (function.func, 'PREMISE')
    publish(fold, objects, declarations)
    emit('LOCAL_NAMES', LOCAL)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--branch-demo', choices=('elliptic', 'degenerate', 'blocked'))
    build(parser.parse_args().branch_demo)
