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
The density rho_br(r) is arbitrary in Parts B/C mass balance; the forward
observable evaluation additionally declares rho_br=R*(L/r)**a, after solving
the no-loss equation for general rho_br. This excludes non-power density
profiles from the forward observable evaluation. All parameters stay symbolic.
The density rho_br(r) is arbitrary in the general mass equation. V denotes the coordinate velocity dr/dt,
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

IMPORT_KEYS = ('c_s0',)
INPUTS = (HERE, ROOT/'directives/S9b_SHARED_PHYSICS.md',
          *(ROOT/('scripts/'+s+'_exports.py') for s in
            ('S11c_b', 'S11c_c1', 'S11c_c2')), ROOT/'scripts/ledger_fold.py')
GRADES = tuple(product(range(2), range(3), range(2)))
ZERO = (0, 0, 0)
RADIAL_FREEZE = None  # K4a-b: derivative construction
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
    if RADIAL_FREEZE is None or RADIAL_DATA is None or variable != RADIAL_DATA[0]:
        return sp.diff(value,variable,order)
    _, amplitudes, exponents = RADIAL_DATA
    index = RADIAL_FREEZE
    result = value
    for _ in range(order):
        result = sp.diff(result,variable)+exponents[index]*amplitudes[index]*sp.diff(result,amplitudes[index])/variable
    return result


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
            'description': 'S9b_SHARED_PHYSICS.md Parts A-C; reduced power-family conditions',
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


def partitions(items):
    """All set partitions, including all exponent coincidences."""
    if not items:
        yield []
        return
    first, *rest = items
    for blocks in partitions(rest):
        yield [[first]]+blocks
        for i in range(len(blocks)):
            yield blocks[:i]+[[first]+blocks[i]]+blocks[i+1:]


def reduced_cases(terms, exponents, variables, mass, flow_amplitude,
                  exchange, jn, domain, quantifier_variable, substitutions=None):
    """Finite power independence on an open half-line, then pivot-complete
    linear elimination in GM and squared signed flow amplitude. Every zero
    pivot is a separate case. No ConditionSet encodes the optical matching.
    terms[0] is the supplied reference contribution, exponent[0] its degree.
    """
    substitutions = substitutions or {}
    fixed_impact = None  # K7: quantifier construction
    output = []
    squared = sp.Symbol('s9b_flow_square', real=True)
    for blocks in partitions(list(range(len(exponents)))):
        equalities = [z for block in blocks for j in block[1:] if (z := exponents[j]-exponents[block[0]]) != 0]
        solutions = sp.solve(equalities, variables, dict=True) if equalities else [{}]
        for mapping in solutions:
            distinct = [sp.Ne(exponents[a[0]], exponents[b[0]])
                        for i,a in enumerate(blocks) for b in blocks[i+1:]]
            case_domain = sp.And(domain, *[sp.Eq(z, 0) for z in equalities], *distinct)
            # Preserve the original exponent stratum in the domain; use the
            # solved map only to normalize its coefficient equations.
            normalized = [sp.factor(t.subs(substitutions).subs(mapping)) for t in terms]
            if fixed_impact is None:
                equations = [sp.factor(sum(normalized[i] for i in block)) for block in blocks]
            else:
                equations = [sum(t*fixed_impact**(-e.subs(mapping)) for t,e in zip(normalized,exponents))]
            mass_row = next(e for e in equations if e.has(mass))
            mass_value = sp.solve(mass_row, mass)[0]
            remaining = [sp.factor(e.subs(mass, mass_value)) for e in equations if e != mass_row]
            if flow_amplitude is None:
                # The no-loss case has already eliminated V through the live
                # flux. Reduce the remaining coefficient equations on its
                # inputs, without introducing an independent V or j_n row.
                output.append((sp.And(case_domain, *[sp.Eq(e, 0) for e in remaining]),
                               eq(mass, mass_value)))
                continue
            remaining = [e.subs(flow_amplitude**2, squared) for e in remaining]
            # GM may itself depend on the flow even when the other groups do
            # not. Then the signed flow remains free; that is a solved graph.
            pivots = [(sp.diff(e, squared), e.subs(squared, 0)) for e in remaining]
            earlier_zero = []
            for a, c in pivots:
                if a == 0:
                    earlier_zero.append(sp.Eq(c, 0))
                    continue
                root = sp.factor(-c/a)
                branch_domain = sp.And(case_domain, *earlier_zero, sp.Ne(a, 0), sp.Ge(root, 0),
                    *[sp.Eq(e.subs(squared, root), 0) for e in remaining])
                for sign in (-1, 1):
                    signed = sign*sp.sqrt(root)
                    solved_mass = sp.factor(mass_value.subs(flow_amplitude, signed))
                    implication = sp.factor(exchange.subs(substitutions).subs(mapping).subs(flow_amplitude, signed))
                    output.append((branch_domain, eq(mass, solved_mass),
                                   eq(flow_amplitude, signed), eq(jn, implication)))
                earlier_zero.extend((sp.Eq(a, 0), sp.Eq(c, 0)))
            free_domain = sp.And(case_domain, *[sp.Eq(a, 0) for a,c in pivots],
                                 *[sp.Eq(c, 0) for a,c in pivots])
            output.append((free_domain, eq(mass, mass_value),
                           eq(jn, exchange.subs(substitutions).subs(mapping))))
    return cas(output)


def gated(value, domain):
    # The complement carries no evaluated observable. Values here are formal
    # conditional functionals, evaluated only inside their stated branch.
    return cas(((domain, value), (sp.Not(domain), 'NOT_ESTABLISHED')))


def build(branch_demo=None):
    # SUPPLIED dispersion, profiles, balance, bulk responses and references.
    r, x, B, b, L, c0, ZE, ZR = sp.symbols('s9b_r s9b_X s9b_B b s9b_L c_0 Z_E Z_R', positive=True)
    RE, RR, rt, bmin = sp.symbols('r_E r_R s9b_r_turn s9b_b_min', positive=True)
    d, v, w, F = sp.symbols('s9b_d s9b_v s9b_w s9b_F', real=True)
    p, q, h = sp.symbols('s9b_p s9b_q s9b_h', positive=True)
    rho0, m = sp.symbols('rho_0 m', positive=True)
    GM, K = sp.symbols('GM K', real=True)
    n, s = sp.symbols('n s', real=True)
    kr, J, omega, U, chi = sp.symbols('s9b_k_r s9b_J s9b_omega s9b_U s9b_chi', real=True)
    A = sp.Symbol('s9b_A', positive=True)
    mu = sp.Function('mu_perp', real=True)(r)
    rho = sp.Function('rho_br_profile', positive=True)(r)
    jn = sp.Function('j_n', real=True)(r)
    delta, velocity, slope, f = d*(L/r)**p, c0*v*(L/r)**q, w*(L/r)**h, F*(L/r)**p
    global RADIAL_DATA
    RADIAL_DATA = (r,(d,v),(p,q))
    slope = slope  # K4c: embedding derivative construction
    metric = sp.diag(1+slope**2, r**2, r**2*sp.sin(sp.Symbol('theta', real=True))**2)
    kinetic = (omega-U*kr)**2  # K1: kinetic construction
    inverse_metric = sp.diag(1/A, 1/r**2)  # K2: metric construction
    local_speed_squared = chi  # K3: local speed construction
    covector = sp.Matrix([kr, omega*J])
    dispersion = kinetic-local_speed_squared*(covector.T*inverse_metric*covector)[0]
    speed_relation = eq(chi, mu/rho)
    reference_theta = (1+sp.Symbol('gamma', real=True))*2*GM/(b*c0**2)
    reference_radar = 2*(1+sp.Symbol('gamma', real=True))*GM/c0**3*sp.log(4*RE*RR/b**2)
    reference_log = sp.simplify(-b*sp.diff(reference_radar,b)/2)
    responses = (sp.S.One, None, (1+sp.Symbol('f', real=True))**s)
    bulk_rho = sp.Symbol('rho_bulk', positive=True)
    bulk_pressure = K*bulk_rho**n
    cs_squared = sp.diff(bulk_pressure, bulk_rho)/m

    emit('SUPPLIED_DEPENDENCIES', ('DISPERSION', 'SPEED_IDENTIFICATION', 'ADVECTION',
         'INDUCED_METRIC', 'MASS_BALANCE', 'LAB_HELD', 'ORDER_COUNTING', 'GR_REFERENCES', 'PART_C_RESPONSES'))
    emit('PROFILE_FAMILY', (eq(sp.Function('delta')(r), delta), eq(sp.Function('V')(r), velocity),
         eq(sp.Derivative(sp.Function('xi_w')(r), r), slope), eq(sp.Function('f')(r), f)))
    emit('PROFILE_DOMAIN', sp.And(sp.Or(sp.Eq(d,0),sp.Ge(p,1)),sp.Or(sp.Eq(v,0),sp.Ge(q,sp.Rational(1,2))),sp.Or(sp.Eq(w,0),sp.Ge(h,sp.Rational(1,2)))))
    emit('ENDPOINTS', (eq(RE, sp.sqrt(b*b+ZE*ZE)), eq(RR, sp.sqrt(b*b+ZR*ZR))))
    emit('DISPERSION', eq(dispersion, 0))
    emit('SPATIAL_METRIC', metric)
    emit('INVERSE_SPATIAL_METRIC', metric.inv())
    light = sp.Function('c_gamma')(r)
    emit('REFERENCE_SPEED',eq(c0,sp.Limit(light,r,sp.oo)))
    emit('POLARIZATION_IDENTIFICATION',(eq(sp.Function('c_gamma_1')(r),light),eq(sp.Function('c_gamma_2')(r),light)))
    emit('SPEED_PROFILE',eq(light,c0*(1+delta)))
    ell = sp.Symbol('s9b_reduction_scale',positive=True)
    emit('EMBEDDING_IDENTIFICATION',eq(sp.Function('xi_w')(r),ell*sp.Function('h')(r)))
    time_coordinate = sp.Symbol('s9b_time',real=True)
    emit('LAB_HELD_ANCHOR',eq(sp.Function('c_gamma_L')(r,time_coordinate),light))
    emit('SPEED_RELATION', speed_relation)
    uniform_density, uniform_stiffness = sp.symbols('s9b_uniform_brane_density s9b_uniform_transverse_stiffness',positive=True)
    emit('UNIFORM_ANCHOR', (eq(uniform_density*omega**2,uniform_stiffness*kr**2),eq(U,0),eq(slope,0)))

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
    radial_discriminant = sp.factor(sp.discriminant(dispersion.subs(J, 0), omega)/4/kr**2)
    existence = sp.And(sp.Gt(radial_discriminant, 0), sp.Ne(rho, 0))
    traversal = sp.Lt(sp.factor(radial_speeds[0]*radial_speeds[1]),0)
    emit('RADIAL_FREQUENCIES', radial_omega)
    emit('RADIAL_GROUP_VELOCITIES', radial_speeds)
    emit('BRANCH_EXISTENCE', existence)
    emit('PATH_TRAVERSAL', traversal)
    emit('BRANCH_TYPES', ((sp.Lt(radial_discriminant, 0), ('GROWING', 'DECAYING')),
         (eq(radial_discriminant, 0), 'ABSENT'), (sp.And(existence, sp.Not(traversal)), 'UNABLE_TO_TRAVERSE_REQUIRED_DIRECTION')))
    emit('OBSERVABLE_DOMAIN', ((sp.And(existence, traversal), 'CONDITIONAL_RAY_FUNCTIONALS'),
         (sp.Not(sp.And(existence, traversal)), ('NOT_ESTABLISHED', 'REAL_PROPAGATING_BIDIRECTIONAL_BRANCH_REQUIRED'))))
    if branch_demo is not None:
        supplied = {'elliptic': {chi: -c0**2, U: 0, A: 1},
                    'degenerate': {chi: 0, U: 0, A: 1},
                    'blocked': {chi: c0**2, U: -2*c0, A: 1}}[branch_demo]
        emit('DEMONSTRATION_INPUT', supplied, local=True)
        emit('DEMONSTRATION_BRANCH', (existence.subs(supplied), traversal.subs(supplied),
             tuple(z.subs(supplied) for z in radial_omega)), local=True)
        demonstration_domain = sp.And(existence,traversal).subs(supplied)
        if demonstration_domain == sp.false:
            emit('OBSERVABLES', ('NOT_ESTABLISHED',demonstration_domain))
            emit('LOCAL_NAMES', LOCAL)
            return

    fold, fold_audit = load_model(*(ROOT/('scripts/'+s+'_exports.py') for s in ('S11c_b', 'S11c_c1', 'S11c_c2')))
    manifest = check_consumer(fold, IMPORT_KEYS)
    bound = assert_lookups_equal_manifest(lambda ledger: {k: ledger[k]['value'] for k in IMPORT_KEYS}, fold, IMPORT_KEYS)
    inputs = bound['result']
    emit('IMPORT_BINDINGS', (('c_s0', inputs['c_s0'],
         'S9b_SHARED_PHYSICS.md:205-207; asymptotic bulk sound speed'),), local=True)
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
    ray_domain = sp.And(real_domain, travel_domain, turning_domain, sp.Lt(rt, sp.Min(RE, RR)))
    observable_gate = ray_domain  # K10: branch gate construction
    emit('RAY_DOMAIN', ray_domain)
    def observable(name, value):
        use_domain = sp.And(observable_gate,sp.Ne(GM,0)) if name.startswith('GAMMA') else observable_gate
        emit(name, gated(value, use_domain))

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
    beta_general = sum(betaj.c.values(),sp.S.Zero).subs(r, radius)
    azimuthal_velocity = sp.zeros(3, 1)  # K11: advection construction
    advected_velocity = sp.Matrix([velocity.subs(r, radius)*sp.diff(radius, a) for a in xyz])+azimuthal_velocity
    spatial_metric = sp.eye(3)+slope.subs(r, radius)**2*sp.Matrix([sp.diff(radius, a) for a in xyz])*sp.Matrix([sp.diff(radius, a) for a in xyz]).T
    # Randers odd form is obtained by completing the same advected quadratic.
    speed_cartesian = (c0*(1+delta.subs(r, radius)))**2
    odd_cartesian = -spatial_metric*advected_velocity/(speed_cartesian-(advected_velocity.T*spatial_metric*advected_velocity)[0])
    oneform = [beta_general*sp.diff(radius, a)+sp.factor(odd_cartesian[i]-odd_cartesian[i].subs({z: 0 for z in azimuthal_velocity.free_symbols if z.name == 's9b_azimuthal'})) for i,a in enumerate(xyz)]
    curl = sp.ImmutableMatrix([sp.simplify(sp.diff(oneform[(i+2)%3], xyz[(i+1)%3])-
                   sp.diff(oneform[(i+1)%3], xyz[(i+2)%3])) for i in range(3)])
    parameter = sp.Symbol('s9b_path_parameter',real=True)
    curve = [sp.Function('s9b_path_'+str(i))(parameter) for i in range(3)]
    pullback = sum(component.subs(dict(zip(xyz,curve)))*sp.diff(coordinate,parameter) for component,coordinate in zip(oneform,curve))
    path_functional = sp.Integral(pullback,(parameter,sp.Symbol('s9b_path_start',real=True),sp.Symbol('s9b_path_end',real=True)))
    observable('TIME_NONRECIPROCAL_PATH_FUNCTIONAL',path_functional)
    emit('NONRECIPROCITY_ONEFORM', oneform)
    emit('NONRECIPROCITY_EXTERIOR_DERIVATIVE', curl)
    emit('NONRECIPROCITY_PATH_DEPENDENCE', sp.Eq(curl, sp.zeros(3, 1)))
    emit('NONRECIPROCITY_PRIMITIVE', sp.Integral(beta.subs(profile_map), r))
    endpoint_map = {RE: sp.sqrt(b*b+ZE*ZE), RR: sp.sqrt(b*b+ZR*ZR)}
    first_grades = ((1, 0, 0), (0, 1, 0), (0, 2, 0), (0, 0, 1))
    return_orientation = -1  # K5: return direction construction
    forward_time = even_time+Jet(odd_time)
    reverse_time = even_time+return_orientation*Jet(odd_time)
    roundtrip_time = forward_time+reverse_time
    nonreciprocal_time = (forward_time-reverse_time)/2
    for g in GRADES:
        suffix = 'D%d_V%d_H%d' % g
        emit('ORDER_'+suffix, sp.Rational(g[0])+sp.Rational(g[1], 2)+g[2])
        observable('DEFLECTION_'+suffix, theta[g].subs(B, b))
        observable('TIME_FORWARD_'+suffix, forward_time[g].subs(endpoint_map))
        observable('TIME_REVERSE_'+suffix, reverse_time[g].subs(endpoint_map))
        observable('TIME_RADAR_'+suffix, roundtrip_time[g].subs(endpoint_map))
        observable('TIME_NONRECIPROCAL_'+suffix, nonreciprocal_time[g].subs(endpoint_map))

    # For a power tail, log terms are the exponent-one stratum. On that
    # stratum differentiate the ACTUAL finite-endpoint time with respect to a
    # common endpoint dilation. Half its limiting dilation derivative is the
    # coefficient of log(1/b**2). Pure endpoint powers and endpoint shifts are
    # included in the source, not dropped by a manually supplied coefficient.
    log_source = roundtrip_time  # K6: radar coefficient source construction
    dilation = sp.Symbol('s9b_endpoint_dilation', positive=True)
    log_terms = {}
    log_coefficients = {}
    for g in first_grades:
        source = log_source[g]
        coef, power = qdata.get(g, (sp.S.Zero, sp.S.Zero))
        if source == 0 or power == 0:
            coefficient = sp.S.Zero
            resonance = sp.false
        else:
            resonance = sp.Eq(power, 1)
            exponent_symbols = [a for a in (p, q, h) if power.has(a)]
            resonant_map = sp.solve(power-1, exponent_symbols[0], dict=True)[0]
            dilated = source.subs(resonant_map).subs({RE: dilation*RE, RR: dilation*RR})
            derivative = sp.diff(dilated, dilation)*dilation/2
            coefficient = sp.simplify(sp.limit(derivative, dilation, sp.oo))
        log_coefficients[g] = coefficient
        log_terms[g] = sp.Piecewise((coefficient, resonance), (0, True))
        emit('RADAR_LOG_SOURCE_D%d_V%d_H%d' % g, source)
        observable('RADAR_LOG_D%d_V%d_H%d' % g, log_terms[g])
    deflection = sum(theta[g].subs(B, b) for g in first_grades)
    radar_log = sum(log_terms.values())
    gamma = next(a for a in reference_theta.free_symbols if a.name == 'gamma')
    gamma_theta = sp.solve(eq(reference_theta, deflection), gamma)[0]
    gamma_radar = sp.solve(eq(reference_log, radar_log), gamma)[0]
    residual_theta = deflection-reference_theta.subs(gamma, 1)
    residual_radar = radar_log-reference_log.subs(gamma, 1)
    emit('REFERENCES', (reference_theta, reference_radar, reference_log, eq(gamma, 1)))
    observable('GAMMA_DEFLECTION', gamma_theta)
    observable('GAMMA_RADAR', gamma_radar)
    observable('GAMMA_DIFFERENCE', gamma_theta-gamma_radar)
    observable('RESIDUAL_DEFLECTION', (deflection, reference_theta.subs(gamma, 1), residual_theta))
    observable('RESIDUAL_RADAR', (radar_log, reference_log.subs(gamma, 1), residual_radar))

    # Supplied coordinate d^3x mass law; no induced-measure object.
    volume = r**2
    mass_density = rho  # K8: mass-balance density construction
    exchange = sp.simplify(-sp.diff(volume*mass_density*velocity, r)/volume)
    emit('MASS_BALANCE', eq(jn, exchange))
    emit('MASS_BALANCE_COORDINATE_JACOBIAN', volume)
    mu_profile = sp.solve(speed_relation.subs(chi, (c0*(1+delta))**2), mu)[0]
    emit('STIFFNESS_PROFILE', eq(mu, mu_profile))

    cs0_squared = cs_squared.subs(bulk_rho, rho0)
    f_symbol = sp.Symbol('f', real=True)
    cs_ratio = sp.powsimp(sp.powdenest(cs_squared.subs(bulk_rho, rho0*(1+f_symbol))/cs0_squared, force=True), force=True)
    fixed_ratio_response = sp.sqrt(cs_ratio)  # K9a: local sound-speed response
    power_response = responses[2]  # K9b: local bulk-density response
    responses = (responses[0], fixed_ratio_response, power_response)
    response_delta = [sp.diff(response, f_symbol).subs(f_symbol, 0)*f for response in responses]
    emit('BULK_PRESSURE', eq(sp.Function('P')(bulk_rho), bulk_pressure))
    emit('BULK_SOUND_SPEED', cs_squared)
    emit('BULK_SOUND_ANCHOR', eq(inputs['c_s0']**2, cs0_squared))
    emit('BULK_RESPONSE_DELTAS', response_delta)
    family_domain = sp.And(sp.Or(sp.Eq(d,0),sp.Ge(p,1)),
                           sp.Or(sp.Eq(v,0),sp.Ge(q,sp.Rational(1,2))),
                           sp.Or(sp.Eq(w,0),sp.Ge(h,sp.Rational(1,2))))
    applicability = sp.And(real_domain, travel_domain, turning_domain, family_domain)
    objects = {}
    # Normalize deflection to a sum of independent powers of b. All
    # coefficients are extracted from the computed observable and reference.
    exponents = (sp.S.One, p, 2*q, 2*h)
    deflection_terms = [-reference_theta.subs(gamma, 1)*b]
    for grade, exponent in zip(((1,0,0),(0,2,0),(0,0,1)), exponents[1:]):
        deflection_terms.append(sp.simplify(theta[grade].subs(B,b)*b**exponent))
    # Radar uses the same exponent stratification; its single coefficient
    # equation retains precisely the resonant contributions. A common dummy
    # degree combines them, while the outer strata resolve all Piecewises.
    def condition_rows(label, replacement, own_domain):
        dterms = [sp.factor(t.subs(replacement)) for t in deflection_terms]
        dc = reduced_cases(dterms, exponents, (p,q,h), GM, v, exchange, jn,
                           sp.And(applicability.subs(replacement), own_domain), b,
                           replacement)
        radar_cases = []
        for blocks in partitions(list(range(4))):
            equalities = [z for block in blocks for j in block[1:] if (z := exponents[j]-exponents[block[0]]) != 0]
            maps = sp.solve(equalities, (p,q,h), dict=True) if equalities else [{}]
            for mapping in maps:
                distinct = [sp.Ne(exponents[a[0]],exponents[z[0]]) for i,a in enumerate(blocks) for z in blocks[i+1:]]
                edomain = sp.And(*[sp.Eq(e,0) for e in equalities],*distinct)
                refblock = next(block for block in blocks if 0 in block)
                # Membership in the reference exponent block evaluates the
                # resonance predicates of the already extracted coefficients.
                changed = sum(log_coefficients[grade] for i,grade in
                    enumerate(((1,0,0),(0,2,0),(0,0,1)),1) if i in refblock)
                changed = sp.sympify(changed).subs(replacement).subs(mapping)
                rterms = [changed-reference_log.subs(gamma,1)]
                radar_cases.extend(reduced_cases(rterms, (sp.S.One,), (), GM, v,
                    exchange.subs(mapping), jn,
                    sp.And(applicability.subs(replacement), own_domain, edomain), b,
                    replacement))
        rc = sp.Tuple(*radar_cases)
        for quantity, cases in (('DEFLECTION',dc),('RADAR',rc)):
            name = label+'_'+quantity+'_CONDITION'
            emit(name,cases)
            objects['s9b_'+name.lower()] = cases
        return dc, rc

    condition_rows('B', {}, sp.true)
    restriction_maps = {'FLOW_ONLY_DELTA_XI_ZERO': {d:0,w:0},
                        'SPEED_ONLY_V_XI_ZERO': {v:0,w:0},
                        'TILT_ONLY_DELTA_V_ZERO': {d:0,v:0}}
    computed = {'DEFLECTION':deflection, 'RADAR_LOG':radar_log,
                'GAMMA_DEFLECTION':gamma_theta,'GAMMA_RADAR':gamma_radar,
                'GAMMA_DIFFERENCE':gamma_theta-gamma_radar,
                'RESIDUAL_DEFLECTION':residual_theta,'RESIDUAL_RADAR':residual_radar}
    for label, replacement in restriction_maps.items():
        for name, obj in computed.items():
            emit(label+'_'+name, gated(obj.subs(replacement), sp.And(observable_gate.subs(replacement),sp.Ne(GM,0)) if name.startswith('GAMMA') else observable_gate.subs(replacement)))
        condition_rows(label,replacement,sp.true)

    response_rows = []
    for label, response, exact_response in zip(('CONSTANT','FIXED_RATIO','POWER'), response_delta, responses):
        amplitude = sp.simplify(response/(L/r)**p)
        replacement = {d: amplitude}
        # Each predicate is obtained from THIS row's response, including its
        # own coefficient-zero strata. It is not inherited from the d family.
        response_local = exact_response.subs(f_symbol,f)
        row_domain = sp.And(sp.Contains(response_local,sp.S.Reals),sp.Gt(response_local,0),
                            sp.Or(sp.Eq(amplitude,0),sp.Ge(p,1)))
        if label == 'FIXED_RATIO':
            row_domain = sp.And(row_domain,sp.Gt(cs0_squared,0),
                                 sp.Gt(cs_squared.subs(bulk_rho,rho0*(1+f)),0))
        response_rows.append((label, replacement, row_domain))
        emit('C_'+label+'_RESPONSE', (exact_response,response,row_domain))
        emit('C_'+label+'_N_DEPENDENCE',sp.diff(response,n))
        condition_rows('C_'+label,replacement,row_domain)
        if label != 'CONSTANT':
            condition_rows('C_'+label+'_BULK_DENSITY_ONLY_V_XI_ZERO',
                           replacement|{v:0,w:0},row_domain)

    # Item 10 forward premise: first solve with a general live density, then
    # evaluate the explicitly declared density power family. No equation here
    # matches GM or chooses Phi's magnitude/sign.
    Phi, constant = sp.symbols('Phi s9b_mass_constant', real=True)
    live_V = sp.Function('s9b_V_live')(r)
    live_balance = -sp.diff(volume*mass_density*live_V,r)/volume
    mass_solution = sp.dsolve(eq(live_balance,0),live_V)
    integration_constant = next(a for a in mass_solution.free_symbols if a.name == 'C1')
    conserved = mass_solution.rhs.subs(integration_constant,constant)
    sphere_flux = sp.integrate(conserved*mass_density*r**2*sp.sin(sp.Symbol('s9b_polar',real=True)),
                             (sp.Symbol('s9b_polar',real=True),0,sp.pi))*2*sp.pi
    forward_V = conserved.subs(constant,sp.solve(eq(Phi,sphere_flux),constant)[0])
    emit('FORWARD_NO_FAR_ZONE_LOSS_MASS_SOLUTION',
         (eq(jn,0),eq(live_balance,0),eq(Phi,sphere_flux),eq(live_V,forward_V)))
    density_scale = sp.Symbol('s9b_density_scale',positive=True)
    density_power = sp.Symbol('s9b_density_power',real=True)
    density_family = density_scale*(L/r)**density_power
    family_V = sp.factor(forward_V.subs(rho,density_family))
    forward_q = sp.simplify(-r*sp.diff(family_V,r)/family_V)
    forward_v = sp.simplify(family_V/(c0*(L/r)**forward_q))
    emit('FORWARD_NO_FAR_ZONE_LOSS_DENSITY_FAMILY',eq(rho,density_family))

    def forward_substitute(obj, replacement):
        # This SymPy version's ConditionSet.subs leaves a condition unchanged
        # for replacements that are not differentiation variables. Rebuild
        # each quantified radial predicate, keeping its bound radius intact.
        obj = sp.sympify(obj)
        zeros = {key:value for key,value in replacement.items() if value == 0}
        def scalar(value):
            return value.subs(zeros).subs(replacement)
        sets = {subset:sp.ConditionSet(subset.sym,scalar(subset.condition),
                                      scalar(subset.base_set))
                for subset in obj.atoms(sp.ConditionSet)}
        return scalar(obj.xreplace(sets))

    def forward_condition_rows(label, extra, own_domain=sp.true):
        # Substitute the computed no-loss velocity BEFORE exponent
        # stratification. In particular q is no longer an independent input.
        # Evaluate zero flux first, before a density exponent can encounter a
        # removable pole in a coefficient whose amplitude has vanished.
        collected = {'DEFLECTION': [], 'RADAR': []}
        for flux_domain, flow in ((sp.Eq(Phi,0),sp.S.Zero),
                                  (sp.Ne(Phi,0),forward_v)):
            replacement = {q:forward_q, v:flow, rho:density_family}|extra
            def substitute(obj):
                return forward_substitute(obj,replacement)
            domain = sp.And(substitute(applicability), substitute(own_domain), flux_domain)
            degrees = tuple(substitute(e) for e in exponents)
            terms = [sp.factor(substitute(t)) for t in deflection_terms]
            # An identically absent term has no exponent constraint. Removing
            # it here is algebraic zero elimination, not a family restriction.
            active = [i for i,t in enumerate(terms) if i == 0 or t != 0]
            variables = (p,density_power,h)
            collected['DEFLECTION'].extend(reduced_cases(
                [terms[i] for i in active], [degrees[i] for i in active],
                variables, GM, None, None, None, domain, b))
            coefficients = [None]+[sp.factor(substitute(log_coefficients[g]))
                                   for g in ((1,0,0),(0,2,0),(0,0,1))]
            active = [i for i,c in enumerate(coefficients) if i == 0 or c != 0]
            for blocks in partitions(active):
                equalities = [z for block in blocks for i in block[1:]
                              if (z := degrees[i]-degrees[block[0]]) != 0]
                maps = sp.solve(equalities,variables,dict=True) if equalities else [{}]
                distinct = [sp.Ne(degrees[a[0]],degrees[z[0]])
                            for i,a in enumerate(blocks) for z in blocks[i+1:]]
                stratum = sp.And(*[sp.Eq(e,0) for e in equalities],*distinct)
                reference_block = next(block for block in blocks if 0 in block)
                changed = sum((coefficients[i] for i in reference_block if i != 0),sp.S.Zero)
                for mapping in maps:
                    residual = (changed-reference_log.subs(gamma,1)).subs(mapping)
                    collected['RADAR'].extend(reduced_cases(
                        [residual], (sp.S.One,), (), GM, None, None, None,
                        sp.And(domain,stratum), b))
        for quantity, cases in collected.items():
            name = label+'_'+quantity+'_CONDITION'
            value = sp.Tuple(*cases)
            emit(name,value)
            objects['s9b_'+name.lower()] = value

    for label, extra in (('DELTA_XI_ZERO',{d:0,w:0}),('DELTA_XI_LIVE',{})):
        substitution = {q:forward_q,v:forward_v,rho:density_family}|extra
        label = 'FORWARD_NO_FAR_ZONE_LOSS_'+label
        forward_domain = forward_substitute(observable_gate,substitution)
        forward_domain = sp.And(forward_domain,family_domain.subs(substitution))
        emit(label+'_BRANCH_EXISTENCE',forward_substitute(real_domain,substitution))
        emit(label+'_PATH_TRAVERSAL',forward_substitute(travel_domain,substitution))
        emit(label+'_DOMAIN',forward_domain)
        for name,obj in computed.items():
            emit(label+'_'+name,gated(obj.subs(substitution),sp.And(forward_domain,sp.Ne(GM,0)) if name.startswith('GAMMA') else forward_domain))
        for name,jet in (('TIME_RADAR',roundtrip_time),('TIME_FORWARD',forward_time),
                         ('TIME_REVERSE',reverse_time),('TIME_NONRECIPROCAL',nonreciprocal_time)):
            emit(label+'_'+name,gated(sum(jet.c.values(),sp.S.Zero).subs(substitution).subs(endpoint_map),forward_domain))
        forward_condition_rows(label,extra)
        if label.endswith('DELTA_XI_LIVE'):
            for response_label, replacement, row_domain in response_rows:
                forward_condition_rows(label+'_C_'+response_label,replacement,row_domain)
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
