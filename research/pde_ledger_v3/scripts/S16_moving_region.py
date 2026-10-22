#!/usr/bin/env python3
"""S16 independent SymPy construction: balance IBP and normal-speed transport.

No data or expressions are read from the Wolfram engine. The geometry is a
spherical integration chart, with completely general Cartesian live fields.
Product-rule integration by parts reduces the supplied PDEs to boundary
operators. Named integral/rate functions form an explicitly printed expression
graph, including their derivatives; they are never numerical substitutions.
"""
from pathlib import Path
from itertools import product
import sympy as s

OUT = Path(__file__).resolve().parent / 'out' / 'S16_moving_region.out'
OUT.parent.mkdir(parents=True, exist_ok=True)
log = OUT.open('w')
objects = {}

def write(name, value):
    log.write(f'{name} = {s.sstr(value)}\n')
    log.flush()

def emit(name, value):
    write(name, value)
    objects[name] = value
    write(name + '.comparison', 'nothing to compare it with; no second route')

def working(name):
    print('WORKING_ON = ' + name, flush=True)
    write('WORKING_ON', name)

# Supplied balances, live fields and allowed control-surface geometry.
x = s.symbols('x1:4', real=True)
t = s.Symbol('t', real=True)
rho = s.Symbol('rho', positive=True)
theta, phi = s.symbols('theta phi', real=True)
args = (*x, t)
c = s.Matrix([s.Function(f'Xc{i}')(t) for i in range(1, 4)])
a, b = (s.Function(n)(t) for n in ('r_in', 'R'))
Sigma = s.Function('Sigma')(*args)
F = s.Matrix([s.Function(f'F{i}')(*args) for i in range(1, 4)])
G = s.Matrix([s.Function(f'G{i}')(*args) for i in range(1, 4)])
Pi = s.Matrix(3, 3, lambda i, j: s.Function(f'Pi{i+1}{j+1}')(*args))
f = s.Matrix([s.Function(f'f{i}')(*args) for i in range(1, 4)])
n1, n3, am, ap = [s.Function(n)(*args) for n in ('n1', 'n3', 'a_minus', 'a_plus')]
face_balance = (-n1*am, n3*ap)
def divergence(vector):
    return sum(s.diff(vector[j], x[j]) for j in range(3))
material = s.diff(Sigma, t) + divergence(F) + sum(face_balance)
momentum = [s.diff(G[i], t) + divergence(Pi.row(i)) - f[i] for i in range(3)]
e = s.Matrix([s.sin(theta)*s.cos(phi), s.sin(theta)*s.sin(phi), s.cos(theta)])
chart = c + rho*e
radii = (a, b)
surfaces = [chart.subs(rho, radius) for radius in radii]
# Generic requested definitions are applied to the input coordinate objects.
y = s.Matrix(x) - c
weights = [s.S.One, *y, *(y[i]*y[j] for i, j in product(range(3), repeat=2))]

# Derive the moving-endpoint orientation from Leibniz, not a supplied sign.
u = s.Function('radialIntegrand')(rho, t)
leibniz_operand = s.Integral(u, (rho, a, b))
leibniz_rate = s.diff(leibniz_operand, t)
endpoint_coefficients = [s.expand(leibniz_rate).coeff(s.diff(r, t)).coeff(u.subs(rho, r)) for r in radii]
jac = s.simplify(chart.jacobian((rho, theta, phi)).det())
area_vector = s.simplify(chart.diff(theta).cross(chart.diff(phi)))
# Area magnitude on 0<theta<pi; orientation is fixed by the radial chart.
area = s.trigsimp(area_vector.dot(chart.diff(rho)))
normal = [s.simplify(sign*area_vector/area) for sign in endpoint_coefficients]
velocity_c = chart.diff(t)
velocity_surfaces = [surface.diff(t) for surface in surfaces]
velocity_r = [s.simplify(velocity_surfaces[k] - velocity_c) for k in range(2)]
speed_c = [s.simplify(normal[k].dot(velocity_c)) for k in range(2)]
speed_r = [s.simplify(normal[k].dot(velocity_r[k])) for k in range(2)]
areas = [area.subs(rho, r) for r in radii]

class Volume(s.Function):
    nargs = 1

class Surface(s.Function):
    nargs = 2

# Surface arguments are chart-pulled; Volume arguments retain bound x.
def at(expr, k):
    return expr.subs(dict(zip(x, surfaces[k])), simultaneous=True)

def vol(expr):
    return Volume(s.expand(expr))

def surf(k, expr):
    return Surface(s.Integer(k), s.expand(expr))

write('schema', {'Volume': 'shell integral in x',
                 'Surface': 'angular integral of chart-pulled integrand times area',
                 'X': 'dipole/N, displacement relative to Xc',
                 'shapeChoice': 'concentric spheres; no field symmetry',
                 'momentOrderChoice': 2,
                 'premises': 'supplied and unfalsifiable',
                 'Pi': 'total supplied flux, without a constitutive split',
                 'f': 'undivided supplied force density, including slab-face exchange'})
write('domain.unnormalized', s.And(a > 0, b > a))
N = s.Function('N')(t)
D = [s.Function(f'D{i}')(t) for i in range(1, 4)]
P = [s.Function(f'P{i}')(t) for i in range(1, 4)]
write('domain.normalized', s.And(a > 0, b > a, s.Ne(N, 0)))
q = s.Function('q')(*args)
write('integral.Volume', s.Integral(q.subs(dict(zip(x, chart)), simultaneous=True)*jac,
                                  (rho, a, b), (theta, 0, s.pi), (phi, 0, 2*s.pi)))
for k in range(2):
    write(f'integral.Surface.{k}', s.Integral(s.Function('q')(theta, phi, t)*areas[k],
                                           (theta, 0, s.pi), (phi, 0, 2*s.pi)))
emit('supplied.material', material)
emit('supplied.momentum', momentum)
emit('geometry.chart', chart)
emit('geometry.surfaces', surfaces)
emit('geometry.volumeMeasure', jac)
emit('geometry.areaMeasure', areas)
emit('geometry.normal', normal)
emit('geometry.centerVelocity', velocity_c)
emit('geometry.radiusVelocity', velocity_r)
emit('geometry.centerNormalSpeed', speed_c)
emit('geometry.radiusNormalSpeed', speed_r)
write('calculus.Leibniz.operand', leibniz_operand)
write('calculus.Leibniz.derivative', leibniz_rate)
write('calculus.endpointOrientation', endpoint_coefficients)
# The following generic pullback derivative is printed as a diagnostic.
# Surface motion below is constructed from the full integrand and normal speed.
write('calculus.pullback.operand', q.subs(dict(zip(x, chart)), simultaneous=True))
write('calculus.pullback.derivative', s.diff(q.subs(dict(zip(x, chart)), simultaneous=True), t))

# Geometric differentiation coordinates are independent only while taking
# partial derivatives. Restoring them through Subs preserves every live jet,
# including a velocity coefficient which has not itself been differentiated.
def geometry_partial_rates(expr):
    center_heads = {z.func for z in c}
    radius_heads = {z.func for z in radii}
    jets = sorted((jet for jet in expr.atoms(s.Derivative)
                   if jet.expr.func in center_heads | radius_heads
                   and all(v == t for v in jet.variables)), key=s.default_sort_key)
    center = (*c, *(jet for jet in jets if jet.expr.func in center_heads))
    radial = (*radii, *(jet for jet in jets if jet.expr.func in radius_heads))
    geometry = center + radial
    slots = tuple(s.Dummy(f'geometry{j}') for j in range(len(geometry)))
    # A separate time coordinate prevents Subs from interpreting a field's
    # partial time derivative as a total derivative of its restored chart.
    field_time = s.Dummy('field_time')
    frozen = expr.xreplace(dict(zip(geometry, slots))).xreplace({t: field_time})
    def restore(value):
        return s.Subs(value, (*slots, field_time), (*geometry, t))
    result = {'localEvolution': restore(s.diff(frozen, field_time))}
    offset = 0
    for label, group in (('centerMotion', center), ('radiusMotion', radial)):
        result[label] = s.Add(*(restore(s.diff(frozen, slots[offset+j]))*s.diff(jet, t)
                               for j, jet in enumerate(group)))
        offset += len(group)
    return result

# Generic product-rule integration by parts, component by component. The
# coefficient of each flux derivative is extracted from the solved PDE.
def build_rate(name, weight, density, flux, balance, sources):
    working(name)
    local = s.solve(balance, s.diff(density, t))[0]
    raw = s.expand(s.diff(weight*density, t).subs(s.diff(density, t), local))
    bulk = raw
    currents = []
    for j in range(3):
        jet = s.diff(flux[j], x[j])
        coefficient = bulk.coeff(jet)
        currents.append(coefficient*flux[j])
        bulk = s.expand(bulk - coefficient*jet - s.diff(coefficient, x[j])*flux[j])
    write(name+'.localRate', local)
    write(name+'.weightedLocalRate', raw)
    write(name+'.afterParts.bulk', bulk)
    write(name+'.afterParts.current', currents)
    parts = {}
    remainder = bulk
    for source_name, source in sources:
        selected = sum(term for term in s.Add.make_args(bulk) if term.has(source))
        parts['bulk.'+source_name] = vol(selected)
        remainder = s.expand(remainder-selected)
    reference_rates = geometry_partial_rates(weight)
    for origin in ('centerMotion', 'radiusMotion'):
        reference = s.expand(density*reference_rates[origin].doit())
        parts['bulk.reference.'+origin] = vol(reference)
        remainder = s.expand(remainder-reference)
    parts['bulk.weight'] = vol(remainder)
    for k in range(2):
        parts[f'surface{k}.rest'] = surf(k, at(s.Matrix(currents).dot(normal[k]), k))
        parts[f'surface{k}.centerMotion'] = surf(k, at(weight*density, k)*speed_c[k])
        parts[f'surface{k}.radiusMotion'] = surf(k, at(weight*density, k)*speed_r[k])
    for part, value in parts.items():
        emit(name+'.'+part, value)
    emit(name+'.rate', s.Add(*parts.values()))
    emit(name+'.integratedForceCompanion', [vol(v) for v in f])
    return parts

working('moments')
moments = [vol(w*Sigma) for w in weights]
indices = [(), *((i,) for i in range(3)), *product(range(3), repeat=2)]
for index, moment in zip(indices, moments):
    emit('moment.'+str(index), moment)
emit('count.definition', s.Eq(N, moments[0], evaluate=False))
for i in range(3):
    emit(f'dipole.definition.{i}', s.Eq(D[i], moments[i+1], evaluate=False))
    emit(f'momentum.definition.{i}', s.Eq(P[i], vol(G[i]), evaluate=False))
count_parts = build_rate('count', s.S.One, Sigma, F, material, [('faceMinus', n1), ('facePlus', n3)])
dipole_parts = [build_rate(f'dipole.{i}', y[i], Sigma, F, material,
                           [('faceMinus', n1), ('facePlus', n3)]) for i in range(3)]
momentum_parts = [build_rate(f'momentum.{i}', s.S.One, G[i], list(Pi.row(i)), momentum[i],
                             [('force', f[i])]) for i in range(3)]
working('momentumReceived')
received = [s.Add(*(parts[key] for key in ('bulk.force', 'surface0.rest', 'surface1.rest')))
            for parts in momentum_parts]
emit('momentumReceived', received)
emit('momentumRate', [s.Add(*parts.values()) for parts in momentum_parts])
emit('momentumSurfaceMotion', [{key: value for key, value in parts.items() if key.endswith('Motion')}
                               for parts in momentum_parts])

# Time differentiation uses normal-speed volume transport and differentiation
# of each surface pullback and area measure. All surface pullbacks are live,
# so differentiation also retains spatial jets and the changing area.
def integral_time(expr):
    if isinstance(expr, Volume):
        density = expr.args[0]
        return vol(s.diff(density, t)) + s.Add(*(surf(k, at(density, k)*speed)
                    for k in range(2) for speed in (speed_c[k], speed_r[k])))
    if isinstance(expr, Surface):
        k, density = expr.args
        k = int(k)
        return surf(k, s.diff(density*areas[k], t)/areas[k])
    if expr.is_Add:
        return s.Add(*(integral_time(term) for term in expr.args))
    if expr.is_Mul:
        factors = expr.args
        return s.Add(*(integral_time(factor)*s.Mul(*(factors[:i]+factors[i+1:]))
                       for i, factor in enumerate(factors)))
    if expr.is_Pow:
        return s.diff(expr, expr.base)*integral_time(expr.base)
    return s.diff(expr, t)

# The local derivative holds all geometric positions and jets fixed. Volume
# transport is then applied to the full integrand as a coefficient; products
# of existing velocities with boundary speed stay with the moving boundary.
def integral_time_parts(expr):
    if isinstance(expr, Volume):
        density = expr.args[0]
        rates = geometry_partial_rates(density)
        pieces = {'localEvolution': vol(rates['localEvolution']),
                  'bulk.centerMotion': vol(rates['centerMotion']),
                  'bulk.radiusMotion': vol(rates['radiusMotion'])}
        for k in range(2):
            pieces[f'surface{k}.centerMotion'] = surf(k, at(density, k)*speed_c[k])
            pieces[f'surface{k}.radiusMotion'] = surf(k, at(density, k)*speed_r[k])
        return pieces
    if isinstance(expr, Surface):
        k, density = expr.args
        k = int(k)
        rates = geometry_partial_rates(density*areas[k])
        return {('localEvolution' if origin == 'localEvolution' else f'surface{k}.'+origin):
                surf(k, value/areas[k]) for origin, value in rates.items()}
    raise TypeError(f'unsupported integral operator {expr.func}')

# A comparison within this engine: direct differentiation versus the sum of
# origin parts. Linearity groups residuals under their original measures.
# Canonicalizing evaluated Subs after differentiation reconciles dummy
# differentiation coordinates without prescribing a value for the residual.
def integral_residual(left, right):
    delta = s.expand(left-right)
    integrals = sorted(delta.atoms(Volume, Surface), key=s.default_sort_key)
    result = {'outsideIntegrals': delta,
              'volumeIntegrand': s.S.Zero,
              **{f'surface{k}Integrand': s.S.Zero for k in range(2)}}
    for integral in integrals:
        coefficient = delta.coeff(integral)
        result['outsideIntegrals'] -= coefficient*integral
        if isinstance(integral, Volume):
            result['volumeIntegrand'] += coefficient*integral.args[0]
        else:
            result[f'surface{int(integral.args[0])}Integrand'] += coefficient*integral.args[1]
    return {key: s.simplify(value.doit()) for key, value in result.items()}

def partition_report(tag, direct, pieces):
    working(tag)
    total = s.Add(*pieces.values())
    write(tag+'.direct', direct)
    objects[tag+'.direct'] = direct
    write(tag+'.partsSum', total)
    objects[tag+'.partsSum'] = total
    write(tag+'.residual', integral_residual(direct, total))
    write(tag+'.routes', ('direct product derivative', 'partial geometric jets and normal-speed transport'))

material_rhs = s.solve(material, s.diff(Sigma, t))[0]
def governed(expr):
    # SymPy represents a field jet evaluated on a moving chart with Subs.
    # Substitute the governed partial time derivative before evaluating Subs;
    # do not replace total time derivatives along a surface by partial ones.
    def replace_derivative(node):
        if isinstance(node, s.Derivative) and node.expr.func == Sigma.func and node.variables == (node.expr.args[-1],):
            return material_rhs.subs(dict(zip(args, node.expr.args)), simultaneous=True)
        return node
    expr = expr.replace(lambda node: isinstance(node, s.Derivative), replace_derivative)
    return expr

working('centerKinematics')
alias_definitions = {}
alias_families = {}
def name_parts(parts, prefix, family):
    named = {}
    for key, value in parts.items():
        alias = s.Function(prefix+'_'+key.replace('.', '_'))(t)
        alias_definitions[alias] = value
        alias_families[alias.func] = family
        named[key] = alias
        emit(prefix+'.term.'+key, s.Eq(alias, value, evaluate=False))
        emit(prefix+'.termRate.'+key, s.Eq(s.diff(alias, t), governed(integral_time(value)), evaluate=False))
        pieces = {origin: governed(contribution) for origin, contribution in integral_time_parts(value).items()}
        for origin, contribution in pieces.items():
            emit(prefix+'.termRatePart.'+key+'.'+origin, contribution)
        partition_report(prefix+'.termRatePartition.'+key, governed(integral_time(value)), pieces)
    return named
named_n = name_parts(count_parts, 'Nrate', 'count')
named_d = [name_parts(parts, f'Drate{i}', 'dipole') for i, parts in enumerate(dipole_parts)]
replacements = {s.diff(N, t): s.Add(*named_n.values()),
                s.diff(N, t, 2): s.Add(*(s.diff(v, t) for v in named_n.values()))}
for i in range(3):
    replacements[s.diff(D[i], t)] = s.Add(*named_d[i].values())
    replacements[s.diff(D[i], t, 2)] = s.Add(*(s.diff(v, t) for v in named_d[i].values()))
relative = [s.cancel(moment/N) for moment in D]
absolute = [s.Add(c[i], relative[i]) for i in range(3)]
emit('center.relative.definition', relative)
emit('center.absolute.definition', absolute)
for label, center in (('relative', relative), ('absolute', absolute)):
    emit('center.'+label+'.rate', [s.expand(s.diff(z, t).xreplace(replacements)) for z in center])
    emit('center.'+label+'.acceleration', [s.expand(s.diff(z, t, 2).xreplace(replacements)) for z in center])
emit('count.rateForNormalization', replacements[s.diff(N, t)])
emit('count.accelerationForNormalization', replacements[s.diff(N, t, 2)])

# Dimensional algebra on all emitted objects. Four independent units: length,
# time, count, momentum. Units of fluxes and source terms follow the supplied
# balances. Face factors have velocity units; their number-density dimensions
# follow the projected material source. No mass unit is introduced.
working('dimensions')
L, T, C, K = s.symbols('unitLength unitTime unitCount unitMomentum', positive=True)
coordinate_units = {**{v: L for v in x}, t: T, rho: L, theta: s.S.One, phi: s.S.One}
field_units = {Sigma.func: C/L**3, N.func: C}
field_units.update({z.func: L for z in (*c, a, b)})
field_units.update({z.func: K/L**3 for z in G})
field_units.update({z.func: field_units[Sigma.func]/T*L for z in F})
field_units.update({z.func: field_units[G[0].func]/T*L for z in Pi})
field_units.update({z.func: field_units[G[0].func]/T for z in f})
field_units.update({z.func: L/T for z in (am, ap)})
field_units.update({z.func: field_units[Sigma.func]/T/field_units[am.func] for z in (n1, n3)})
field_units.update({z.func: field_units[N.func]*L for z in D})
field_units.update({z.func: field_units[G[0].func]*L**3 for z in P})
field_units.update({head: (field_units[N.func] if family == 'count' else field_units[D[0].func])/T
                    for head, family in alias_families.items()})
write('dimensions.basis', (L, T, C, K))
write('dimensions.fields', field_units)

def units(expr, env=None):
    env = coordinate_units if env is None else env
    if expr == 0:
        return frozenset()  # an additive zero has no unique physical dimension
    if expr.is_Number or expr == s.pi:
        return frozenset([s.S.One])
    if expr in env:
        return frozenset([env[expr]])
    if isinstance(expr, Volume):
        return frozenset(s.cancel(u*L**3) for u in units(expr.args[0], env))
    if isinstance(expr, Surface):
        return frozenset(s.cancel(u*L**2) for u in units(expr.args[1], env))
    if isinstance(expr, s.Subs):
        mapping = dict(env)
        for var, point in zip(expr.variables, expr.point):
            mapping[var] = next(iter(units(point, env)))
        return units(expr.expr, mapping)
    if isinstance(expr, s.Derivative):
        dims = units(expr.expr, env)
        for var, order in expr.variable_count:
            dims = frozenset(s.cancel(dim/env[var]**order) for dim in dims)
        return dims
    if expr.is_Add:
        return frozenset().union(*(units(term, env) for term in expr.args))
    if expr.is_Mul:
        choices = [units(term, env) for term in expr.args]
        return frozenset(s.cancel(s.prod(parts)) for parts in product(*choices))
    if expr.is_Pow:
        return frozenset(s.cancel(unit**expr.exp) for unit in units(expr.base, env))
    if expr.func in (s.sin, s.cos):
        return frozenset([s.S.One])
    if expr.func in field_units:
        return frozenset([field_units[expr.func]])
    raise TypeError(f'unregistered dimensional node: {s.sstr(expr)}')

def leaves(value):
    if isinstance(value, dict):
        for key, part in value.items():
            for path, item in leaves(part):
                yield (str(key), *path), item
    elif isinstance(value, (list, tuple, s.MatrixBase)):
        for i, part in enumerate(value):
            for path, item in leaves(part):
                yield (str(i), *path), item
    elif isinstance(value, s.Equality):
        yield ('lhs',), value.lhs
        yield ('rhs',), value.rhs
    else:
        yield (), value

def dimensional_sums(expr, env=None):
    env = coordinate_units if env is None else env
    if isinstance(expr, s.Subs):
        context = dict(env)
        for var, point in zip(expr.variables, expr.point):
            context[var] = next(iter(units(point, env)))
        yield from dimensional_sums(expr.expr, context)
        for point in expr.point:
            yield from dimensional_sums(point, env)
        return
    if expr.is_Add:
        yield expr, env
    for child in expr.args:
        if isinstance(child, s.Basic):
            yield from dimensional_sums(child, env)

for name, value in objects.items():
    for path, expr in leaves(value):
        prefix = '.'.join((name, *path))
        write(prefix+'.dimensions', sorted(units(expr), key=str))
        for i, (summation, context) in enumerate(dimensional_sums(expr)):
            operands = [sorted(units(term, context), key=str) for term in summation.args]
            write(prefix+f'.dimensionSum.{i}.operands', operands)
            residuals = [[s.expand(left-right) for left, right in product(operand, operands[0])]
                         for operand in operands]
            write(prefix+f'.dimensionSum.{i}.residuals', residuals)
write('coverage.objects', list(objects))
write('coverage.ungovernedRates', [z.func for z in (*F, n1, n3, am, ap)])
write('coverage.fluxDecomposition', 'no supplied separate stress/advection constituents of Pi')
log.close()
print('TRANSCRIPT = '+str(OUT), flush=True)
