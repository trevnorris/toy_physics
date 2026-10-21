#!/usr/bin/env python3
"""Direct distributional SymPy construction for S14a part 1.

Only the supplied balances, U1 and face factor form native physical products.
All projections are calculus on those inputs. There are no Wolfram reads,
exports, expected-value tests, constitutive forms, radius limits or freezes.
"""
from pathlib import Path
import sys
import sympy as sp

OUT = Path(__file__).resolve().parent / "out" / "S14a_1_projection.out"
OUT.parent.mkdir(parents=True, exist_ok=True)
transcript = OUT.open("w", buffering=1)


def emit(tag, expression):
    line = f"{tag}: {sp.sstr(expression)}"
    print(line, flush=True)
    transcript.write(line + "\n")


def working(tag):
    emit("WORKING_ON", tag)


def object_(tag, expression, domain=None):
    payload = expression if domain is None else {"domain": domain, "object": expression}
    emit(tag, payload)
    emit(tag + "/comparison", "No second route; nothing to compare it with.")


def balance(tag, left, right):
    emit(tag + "/operands", (left, right))
    emit(tag + "/residual", sp.expand(left - right))
    emit(tag + "/comparison", "Balance residual, not an independent-route comparison.")


# Seed dimensions only, in {material number, length, time} order.
N = sp.ImmutableMatrix([1, 0, 0])
L = sp.ImmutableMatrix([0, 1, 0])
T = sp.ImmutableMatrix([0, 0, 1])
ONE = sp.ImmutableMatrix([0, 0, 0])
unit = {}
x, y, z, w, t = sp.symbols("x y z w t", real=True)
r, a, R = sp.symbols("radius inner_radius outer_radius", positive=True)
theta, phi = sp.symbols("theta phi", real=True)
for coordinate in (x, y, z, w, r, a, R):
    unit[coordinate] = L
unit[t] = T
unit[theta] = unit[phi] = ONE


def native(name, arguments, dimension):
    f = sp.Function(name, real=True)
    unit[f] = dimension
    return f(*arguments)


def integral(expression, variable, lower, upper):
    """Linear definite integration; pull only factors independent of its variable."""
    terms = []
    for term in sp.Add.make_args(sp.expand(expression)):
        coefficient, dependence = term.as_independent(variable, as_Add=False)
        if dependence == 1:
            terms.append(coefficient * (upper - lower))
        else:
            terms.append(coefficient * sp.Integral(dependence, (variable, lower, upper)))
    return sp.Add(*terms)


def divergence(density, flux, coordinates, clock, measure):
    return sp.diff(density, clock) + sp.Add(*[
        sp.diff(measure * component, coordinate) / measure
        for component, coordinate in zip(flux, coordinates)
    ])


def boundary_differentiation(integrand, variable, faces, coordinate):
    """Obtain endpoint contributions by differentiating a generic definite integral."""
    lower, upper = sp.symbols("bound_lower bound_upper", real=True)
    abstract = sp.Integral(integrand, (variable, lower, upper))
    return tuple(
        sp.diff(abstract, bound).subs({lower: faces[0], upper: faces[1]}, simultaneous=True)
        * sp.diff(face, coordinate)
        for bound, face in zip((lower, upper), faces)
    )


def bulk_direction_integration(integrand, variable, faces):
    """Integrate a total derivative with an opaque primitive, then apply it."""
    lower, upper = sp.symbols("integration_lower integration_upper", real=True)
    primitive = sp.Function("integration_primitive")(variable)
    evaluated = sp.Integral(sp.Derivative(primitive, variable), (variable, lower, upper)).doit()
    return tuple(
        evaluated.coeff(primitive.subs(variable, bound)) * integrand.subs(variable, face)
        for bound, face in zip((lower, upper), faces)
    )


def linearize_integrals(expression):
    """Apply definite-integral linearity for the projection-route residual."""
    return expression.replace(
        lambda q: isinstance(q, sp.Integral) and len(q.limits) == 1,
        lambda q: integral(q.function, *q.limits[0]),
    )


def supplied(n, velocity, chi, extra, gamma, gamma_return, gamma_drain, coordinates, clock, measure):
    # SUPPLIED BALANCES. These are native physical products, not derived answers.
    material_density = n
    material_flux = tuple(n * v for v in velocity)
    order_density = chi * n
    order_flux = tuple(chi * n * v + j for v, j in zip(velocity, extra))
    order_source = n * gamma
    return {
        "material_density": material_density, "material_flux": material_flux,
        "order_density": order_density, "order_flux": order_flux,
        "order_source": order_source,
        "material_equation": (divergence(material_density, material_flux, coordinates, clock, measure), sp.S.Zero),
        "order_equation": (divergence(order_density, order_flux, coordinates, clock, measure), order_source),
        "gamma_equation": (gamma, gamma_return - gamma_drain),
        "velocity": velocity,
    }


def construct(label, coordinates, clock, arguments, face_arguments, measure):
    working(label + "/native fields and U1")
    bulk = coordinates[-1]
    faces = tuple(native(label + "_face_" + name, face_arguments, L) for name in ("lower", "upper"))
    lower, upper = faces
    steps = tuple(sp.Heaviside(bulk - f) for f in faces)
    # SUPPLIED U1 ENCODING. Exterior order vanishes by premise, not by a field freeze.
    masks = (1 - steps[0], steps[0] - steps[1], steps[1])
    data = []
    for index, region in enumerate(("below", "inside", "above")):
        prefix = label + "_" + region + "_"
        n = native(prefix + "n", arguments, N - 4 * L)
        velocity = tuple(native(prefix + "u" + str(i), arguments, L - T) for i in range(len(coordinates)))
        extra = tuple(native(prefix + "Jchi" + str(i), arguments, N - 3 * L - T) for i in range(len(coordinates)))
        chi = native(prefix + "chi", arguments, ONE) if index == 1 else sp.S.Zero
        gamma = native(prefix + "GammaB", arguments, -T)
        gamma_return = native(prefix + "GammaReturn", arguments, -T)
        gamma_drain = native(prefix + "GammaDrain", arguments, -T)
        current = supplied(n, velocity, chi, extra, gamma, gamma_return, gamma_drain, coordinates, clock, measure)
        data.append(current)
        for key, value in current.items():
            object_(label + "/supplied/" + region + "/" + key, value)
    object_(label + "/faces/domain", sp.Lt(lower, upper))
    object_(label + "/U1/masks", masks)
    delta = tuple(sp.DiracDelta(bulk - face) for face in faces)
    result = {"faces": faces, "data": data}
    for kind in ("material", "order"):
        working(label + "/" + kind + "/direct distributional derivative")
        densities = tuple(d[kind + "_density"] for d in data)
        fluxes = tuple(d[kind + "_flux"] for d in data)
        global_density = sp.Add(*[m * d for m, d in zip(masks, densities)])
        global_flux = tuple(sp.Add(*[masks[i] * fluxes[i][j] for i in range(3)]) for j in range(len(coordinates)))
        differentiated = sp.expand(divergence(global_density, global_flux, coordinates, clock, measure))
        atoms = tuple(sp.expand(differentiated.coeff(d)).subs(bulk, face) for d, face in zip(delta, faces))
        regular_distribution = differentiated.xreplace(dict.fromkeys(delta, sp.S.Zero))
        memberships = ((0, 0), (1, 0), (1, 1))
        regular = tuple(regular_distribution.xreplace(dict(zip(steps, membership))) for membership in memberships)
        object_(label + "/" + kind + "/encoded density", global_density)
        object_(label + "/" + kind + "/encoded flux", global_flux)
        object_(label + "/" + kind + "/differentiated distribution", differentiated)
        object_(label + "/" + kind + "/regular parts", regular)
        object_(label + "/" + kind + "/concentrated face parts", atoms)
        # All concentrated parts at both closed endpoints receive their full weight.
        working(label + "/" + kind + "/bulk-direction integration")
        bulk_integrand = measure * fluxes[1][-1]
        bulk_differential = sp.diff(bulk_integrand, bulk) / measure
        bulk_integral = sp.Integral(sp.Derivative(bulk_integrand, bulk, evaluate=False), (bulk, lower, upper)) / measure
        bulk_parts = tuple(part / measure for part in bulk_direction_integration(bulk_integrand, bulk, faces))
        object_(label + "/" + kind + "/bulk derivative integral", bulk_integral)
        object_(label + "/" + kind + "/bulk derivative integral by face", bulk_parts)
        integrated_regular = integral(sp.expand(regular[1] - bulk_differential), bulk, lower, upper) + sp.Add(*bulk_parts)
        integrated_distribution = sp.Add(integrated_regular, *atoms)
        projected_density = integral(densities[1], bulk, lower, upper)
        projected_flux = tuple(integral(f, bulk, lower, upper) for f in fluxes[1][:-1])
        projected_derivative = divergence(projected_density, projected_flux, coordinates[:-1], clock, measure)
        endpoint_parts = list(boundary_differentiation(densities[1], bulk, faces, clock))
        for component, coordinate in zip(fluxes[1][:-1], coordinates[:-1]):
            pieces = boundary_differentiation(measure * component, bulk, faces, coordinate)
            endpoint_parts = [old + new / measure for old, new in zip(endpoint_parts, pieces)]
        interior_leaving = tuple(sp.expand(part - endpoint) for part, endpoint in zip(bulk_parts, endpoint_parts))
        leaving = tuple(sp.expand(atom + inside) for atom, inside in zip(atoms, interior_leaving))
        source = integral(data[1]["order_source"], bulk, lower, upper) + sp.Add(*atoms) if kind == "order" else sp.S.Zero
        object_(label + "/" + kind + "/closed integral of differential operand", integrated_distribution)
        object_(label + "/" + kind + "/projected density", projected_density)
        object_(label + "/" + kind + "/projected in-plane flux", projected_flux)
        object_(label + "/" + kind + "/projected differential expansion", projected_derivative)
        for index, (edge, inside, atom) in enumerate(zip(leaving, interior_leaving, atoms)):
            object_(label + f"/{kind}/face {index}/interior leaving term", inside)
            object_(label + f"/{kind}/face {index}/closed-projection leaving term", edge)
            if kind == "material":
                balance(label + f"/{kind}/face {index}/supplied material face equation", atom, sp.S.Zero)
        object_(label + "/" + kind + "/projected source", source)
        left = sp.Add(projected_derivative, *leaving)
        balance(label + "/" + kind + "/closed projected balance", left, source)
        emit(label + "/" + kind + "/projection integration/operands", (left, integrated_distribution))
        emit(label + "/" + kind + "/projection integration/residual",
             sp.expand(linearize_integrals(left) - linearize_integrals(integrated_distribution)))
        emit(label + "/" + kind + "/projection integration/comparison",
             "Differentiation of projected integrals versus integration of the regional differential operand with full face parts.")
        result[kind] = {"density": projected_density, "flux": projected_flux,
                        "derivative": projected_derivative, "faces": leaving,
                        "source": source, "balance": (left, source), "atoms": atoms}
    # SUPPLIED FACE FACTOR. It is printed, never used as a derived comparison target.
    for face_index, face in enumerate(faces):
        for region, d in zip(("below", "inside", "above"), data):
            velocity = d["velocity"]
            factor = velocity[-1] - sp.diff(face, clock) - sp.Add(*[
                velocity[i] * sp.diff(face, coordinates[i]) for i in range(len(coordinates) - 1)
            ])
            object_(label + f"/face {face_index}/{region}/supplied material factor", factor.subs(bulk, face))
    return result


emit("contract", "S14a_1_decision_list.md revision 3; supplied balances, U1 and face factor are premises")
emit("engine", (sys.version, sp.__version__))
general = construct("General", (x, y, z, w), t, (x, y, z, w, t), (x, y, z, t), sp.S.One)

working("Spherical specialization/geometry")
# Geometry is mathematical input, independent of any physical field.
embedding = sp.Matrix([r * sp.sin(theta) * sp.cos(phi), r * sp.sin(theta) * sp.sin(phi), r * sp.cos(theta)])
jacobian = embedding.jacobian((r, theta, phi))
metric = sp.simplify(jacobian.T * jacobian)
volume_measure = sp.refine(sp.sqrt(sp.trigsimp(metric.det())), sp.Q.positive(sp.sin(theta)))
area = sp.integrate(volume_measure, (theta, 0, sp.pi), (phi, 0, 2 * sp.pi))
object_("Spherical specialization/geometry domain", sp.And(r > 0, theta > 0, theta < sp.pi))
object_("Spherical specialization/embedding", embedding)
object_("Spherical specialization/metric", metric)
object_("Spherical specialization/volume Jacobian", volume_measure)
object_("Spherical specialization/sphere area", area)
# Native radial and bulk components, with full (radius,w,time) dependence.
spherical = construct("SphericalSpecialization", (r, w), t, (r, w, t), (r, t), area)
emit("SphericalSpecialization/representation", "u0,Jchi0 are radial; u1,Jchi1 are bulk; all native fields have radius,w,time arguments")


def shell(expression):
    return integral(area * expression, r, a, R)


def surface(expression, radius):
    return (area * expression).subs(r, radius)


working("SphericalSpecialization/shell objects and radius dependence")
material, order = spherical["material"], spherical["order"]
face_parts = tuple(shell(f) for f in material["faces"])
face_channel = sp.Add(*face_parts)
core_channel = -surface(material["flux"][0], a)
J = sp.Add(face_channel, core_channel)
order_lost = -shell(order["source"])
outer_flux = surface(material["flux"][0], R)
storage = sp.diff(shell(material["density"]), t)
shell_left = sp.Add(storage, outer_flux, J)
shell_right = shell(material["balance"][1])
face_minus_order = sp.expand(face_channel - order_lost)
radius_derivative = sp.diff(outer_flux, R)
for tag, value in (
    ("shell domain", sp.And(a > 0, a < R)), ("face channel by face", face_parts),
    ("face channel", face_channel), ("core channel", core_channel), ("J", J),
    ("order lost", order_lost), ("outward outer flux", outer_flux), ("storage", storage),
    ("outer flux radius derivative", radius_derivative),
):
    object_("SphericalSpecialization/" + tag, value)
balance("SphericalSpecialization/shell balance", shell_left, shell_right)
balance("SphericalSpecialization/face channel minus order lost", face_channel, order_lost)
velocity = material["flux"][0] / material["density"]
velocity_domain = sp.Ne(material["density"], 0)
volume_flux = surface(velocity, R)
outer_domain = velocity_domain.subs(r, R)
object_("SphericalSpecialization/projected velocity", velocity, velocity_domain)
object_("SphericalSpecialization/outer projected velocity", velocity.subs(r, R), outer_domain)
object_("SphericalSpecialization/outer volume flux", volume_flux, outer_domain)


def dimension(expression):
    e = sp.sympify(expression)
    if e == 0:
        return None  # An additive zero alone fixes no dimension.
    if e.is_number:
        return ONE
    if e in unit:
        return unit[e]
    if isinstance(e, sp.Add):
        return dimension(e.args[0])
    if isinstance(e, sp.Mul):
        return sum((dimension(q) for q in e.args), ONE)
    if isinstance(e, sp.Pow):
        return e.exp * dimension(e.base)
    if isinstance(e, sp.Derivative):
        return dimension(e.expr) - sum((multiplicity * dimension(var) for var, multiplicity in e.variable_count), ONE)
    if isinstance(e, sp.Subs):
        return dimension(e.expr)
    if isinstance(e, sp.Integral):
        return dimension(e.function) + sum((dimension(lim[0]) for lim in e.limits), ONE)
    if e.func in unit:
        return unit[e.func]
    raise TypeError(f"Unregistered dimensional node: {sp.sstr(e)}")


def audit_dimensions(tag, expression, domain=None):
    object_(tag + "/dimension operand", expression, domain)
    d = dimension(expression)
    emit(tag + "/dimension", None if d is None else tuple(d))
    for index, addition in enumerate(q for q in sp.preorder_traversal(expression) if isinstance(q, sp.Add)):
        operands = tuple((q, tuple(dimension(q))) for q in addition.args)
        emit(tag + f"/sum {index}/operands", operands)
        reference = dimension(addition.args[0])
        emit(tag + f"/sum {index}/dimension residuals", tuple(tuple(dimension(q) - reference) for q in addition.args))


working("dimensional objects")
for label, data in (("General", general), ("SphericalSpecialization", spherical)):
    for kind in ("material", "order"):
        current = data[kind]
        prefix = label + "/dimensions/" + kind
        audit_dimensions(prefix + "/density", current["density"])
        for index, flux in enumerate(current["flux"]):
            audit_dimensions(prefix + f"/flux {index}", flux)
        for index, term in enumerate(current["faces"]):
            audit_dimensions(prefix + f"/face {index}", term)
        for side, expression in zip(("left", "right"), current["balance"]):
            audit_dimensions(prefix + "/balance " + side, expression)
        audit_dimensions(prefix + "/balance residual", current["balance"][0] - current["balance"][1])
for name, value in (
    ("face channel", face_channel), ("core channel", core_channel), ("J", J),
    ("order lost", order_lost), ("outward outer flux", outer_flux), ("storage", storage),
    ("shell balance", shell_left - shell_right), ("face channel minus order lost", face_minus_order),
    ("outer flux radius derivative", radius_derivative),
):
    audit_dimensions("SphericalSpecialization/dimensions/" + name, value)
audit_dimensions("SphericalSpecialization/dimensions/projected velocity", velocity, velocity_domain)
audit_dimensions("SphericalSpecialization/dimensions/volume flux", volume_flux, outer_domain)
transcript.close()
