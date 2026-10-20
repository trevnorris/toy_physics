#!/usr/bin/env python3
"""S12: SymPy construction from the supplied balances and U1.

No imports of the Wolfram engine, its transcript, or any other derived result.
The scalar material-balance coefficient is extracted before substituting U1;
the two-factor support is differentiated by SymPy. Every requested quantity
is emitted, without an expected-value guard or cross-engine comparison.

Distribution calculus retains the original differentiable n and u. J has
independent restrictions on all three open regions, independent one-sided
traces, and arbitrary distributions supported at the faces. Derivatives of
the supported distributions remain formal. No form, continuity, exterior
vanishing, boundary condition, or pointwise face value is selected for J.
The order can jump. No live field is frozen.
"""

from pathlib import Path
import sys

import sympy as s
from sympy.core.function import AppliedUndef


TRANSCRIPT = Path(
    "/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S12_order_source.out"
)


class Printer:
    def __init__(self, path):
        self.stream = path.open("w", buffering=1)

    def emit(self, name, value):
        line = f"{name}: {s.sstr(value)}"
        print(line, flush=True)
        self.stream.write(line + "\n")

    def obj(self, name, value):
        self.emit(name, value)
        self.emit(name + "/comparison", "there is nothing to compare it with")

    def working(self, name):
        self.emit("WORKING_ON", name)


class Dimensions:
    """Structural units, {length,time}, including chain-rule Subs nodes.

    Unit assignments concern primitive meanings only. Additive subexpressions
    retain their operand dimensions and residuals in a printed ledger.
    """

    def __init__(self, coordinates, heads, units):
        self.coordinate_units = dict(zip(coordinates, units))
        self.head_units = heads
        self.ledger = []

    def of(self, expr, local=None):
        local = {} if local is None else local
        if expr in local:
            return local[expr]
        if expr in self.coordinate_units:
            return self.coordinate_units[expr]
        if expr.is_Number:
            return s.zeros(2, 1)
        if isinstance(expr, s.Subs):
            env = dict(local)
            env.update(
                (var, self.of(point, local))
                for var, point in zip(expr.variables, expr.point)
            )
            return self.of(expr.expr, env)
        if isinstance(expr, s.Derivative):
            result = self.of(expr.expr, local)
            for var, count in expr.variable_count:
                result = result - count * self.of(var, local)
            return result
        if isinstance(expr, AppliedUndef):
            return self.head_units[expr.func]
        if expr.func == s.Heaviside:
            self.of(expr.args[0], local)
            return s.zeros(2, 1)
        if expr.func == s.DiracDelta:
            degree = expr.args[1] if len(expr.args) > 1 else s.S.Zero
            return -(degree + 1) * self.of(expr.args[0], local)
        if expr.is_Add:
            operands = list(expr.args)
            vectors = [self.of(x, local) for x in operands]
            reference = vectors[0]
            self.ledger.append(
                (operands, vectors, [v - reference for v in vectors])
            )
            return reference
        if expr.is_Mul:
            return sum((self.of(x, local) for x in expr.args), s.zeros(2, 1))
        if expr.is_Pow and expr.exp.is_Number:
            return expr.exp * self.of(expr.base, local)
        raise TypeError(f"unassigned dimension: {s.srepr(expr)}")

    def report(self, printer, label, expr, reference=None):
        printer.working(label + "/dimensions")
        self.ledger = []
        terms = s.Add.make_args(s.expand(expr))
        vectors = [self.of(term) for term in terms]
        printer.emit(label + "/dimension_operands", list(zip(terms, vectors)))
        printer.emit(
            label + "/nested_dimension_operands",
            [(row[0], row[1]) for row in self.ledger],
        )
        printer.emit(
            label + "/nested_dimension_residual", [row[2] for row in self.ledger]
        )
        if reference is None:
            printer.emit(
                label + "/dimension_comparison", "there is nothing to compare it with"
            )
        else:
            ref_vector = self.of(reference)
            printer.emit(label + "/dimension_reference_operand", (reference, ref_vector))
            printer.emit(label + "/dimension_residual", [v - ref_vector for v in vectors])


def main(p):
    p.working("supplied_inputs")
    p.emit("engine", (sys.version, s.__version__))
    p.emit(
        "governing_input",
        "S12_decision_list.md; supplied equations and U1 encoding are unfalsifiable inputs",
    )
    p.emit(
        "calculus_scope",
        (
            "distributional order; n and u retain differentiable bulk fields with traces",
            "J has independent restrictions and one-sided traces in all three open regions",
            "J also has arbitrary distributions supported at each face; their derivatives remain formal and distributional",
            "no continuity, exterior vanishing, or boundary condition on J",
            "pointwise face values unspecified",
            "interior range and strict face order are supplied domain conditions",
        ),
    )
    time = s.Symbol("t", real=True)
    spatial = s.symbols("x1 x2 x3 w", real=True)
    coordinates = (time,) + spatial
    w = spatial[-1]
    slab_coordinates = coordinates[:-1]
    n_head, chi_head, f_head, g_head = map(s.Function, ("n", "chi", "f", "G"))
    u_heads = list(map(s.Function, ("u1", "u2", "u3", "uw")))
    j_heads = list(map(s.Function, ("j1", "j2", "j3", "jw")))
    face_heads = list(map(s.Function, ("wminus", "wplus")))
    rate_heads = list(map(s.Function, ("GammaB", "GammaReturn", "GammaDrain")))
    density = n_head(*coordinates)
    chi = chi_head(*coordinates)
    velocity = [fun(*coordinates) for fun in u_heads]
    flux = [fun(*coordinates) for fun in j_heads]
    rates = [fun(*coordinates) for fun in rate_heads]
    faces = [fun(*slab_coordinates) for fun in face_heads]

    # BEGIN SUPPLIED BALANCES: the only authored physical combinations besides U1.
    material = s.Eq(
        s.diff(density, time)
        + sum(s.diff(density * v, x) for v, x in zip(velocity, spatial)),
        0,
    )
    ordered = s.Eq(
        s.diff(chi * density, time)
        + sum(
            s.diff(chi * density * v + j, x)
            for v, j, x in zip(velocity, flux, spatial)
        ),
        density * rates[0],
    )
    conversion = s.Eq(rates[0], rates[1] - rates[2])
    # END SUPPLIED BALANCES.
    # BEGIN U1 ENCODING: two offsets, each nonnegative in the region.
    offsets = (w - faces[0], faces[1] - w)
    steps = tuple(s.Heaviside(offset) for offset in offsets)
    support = s.Mul(*steps)
    interior = f_head(*coordinates)
    conditional_interior = g_head(*offsets)
    packages = (
        ("general_order", interior * support, (s.S.Zero, interior, s.S.Zero)),
        ("conditional_fast_ordering", conditional_interior * support,
         (s.S.Zero, conditional_interior, s.S.Zero)),
    )
    domains = (
        s.And(faces[0] < faces[1], w < faces[0]),
        s.And(faces[0] < faces[1], faces[0] < w, w < faces[1]),
        s.And(faces[0] < faces[1], w > faces[1]),
    )
    # END U1 ENCODING.
    # Restrict J independently to each of U1's open strata. The masks only
    # locate the restrictions, and impose no constitutive profile. Keep an
    # independent arbitrary supported distribution at each face as well;
    # this is a distribution field with the same units as J, not a trace or
    # a chosen delta profile. No live field is set to zero in this encoding.
    regional_masks = (
        s.Heaviside(s.Mul(-1, offsets[0])), support,
        s.Heaviside(s.Mul(-1, offsets[1])),
    )
    region_j_heads = [
        [s.Function(f"j{c + 1}_region{r + 1}") for c in range(len(spatial))]
        for r in range(len(domains))
    ]
    face_j_heads = [
        [s.Function(f"j{c + 1}_face{k + 1}") for c in range(len(spatial))]
        for k in range(len(faces))
    ]
    regional_flux = [[fun(*coordinates) for fun in row] for row in region_j_heads]
    supported_flux = [[fun(*coordinates) for fun in row] for row in face_j_heads]
    distributed_flux = [
        s.Add(
            *(s.Mul(mask, row[c]) for mask, row in zip(regional_masks, regional_flux)),
            *(row[c] for row in supported_flux),
        )
        for c in range(len(spatial))
    ]
    flux_substitution = dict(zip(flux, distributed_flux))
    p.obj("supplied/material_balance", material)
    p.obj("supplied/order_balance", ordered)
    p.obj("supplied/conversion_identity", conversion)
    p.obj("supplied/faces", faces)
    p.obj("supplied/domain", (faces[0] < faces[1], s.And(interior >= 0, interior <= 1)))
    p.obj("live_fields", (density, velocity, flux, chi, faces))
    p.obj("region_domains", domains)
    p.obj("flux/region_indicators", regional_masks)
    p.obj("flux/region_fields", regional_flux)
    p.obj("flux/face_supported_fields", supported_flux)
    p.obj("flux/face_support_domains", list(zip(supported_flux, [s.Eq(w, face) for face in faces])))
    p.obj("flux/distribution_field", distributed_flux)

    p.working("material_balance_extraction")
    expanded_order = s.expand(ordered.lhs)
    material_weight = expanded_order.coeff(s.diff(density, time))
    material_multiple = s.expand(s.Mul(material_weight, material.lhs))
    transport_expression = s.expand(s.Add(expanded_order, -material_multiple))
    p.obj("material_extraction/order_operand", expanded_order)
    p.obj("material_extraction/coefficient", material_weight)
    p.obj("material_extraction/balance_multiple", material_multiple)
    p.obj("material_extraction/source_expression", transport_expression)

    # Extract the linear flux differential operator from the supplied law.
    # Apply its computed jet coefficients to each supported distribution.
    # No face-source formula is inserted and no J is set to zero to isolate it.
    flux_operator = [
        (c, x, s.diff(field, x), transport_expression.coeff(s.diff(field, x)))
        for c, field in enumerate(flux) for x in spatial
    ]
    p.obj("flux/differential_operator_operands", flux_operator)
    supported_sources = [
        s.Add(*(weight * s.diff(row[c], x) for c, x, jet, weight in flux_operator))
        for row in supported_flux
    ]
    p.obj("flux/face_supported_sources", supported_sources)

    p.working("order_flux_dimensions")
    length_unit, time_unit = s.Matrix([1, 0]), s.Matrix([0, 1])
    dimensionless = s.zeros(2, 1)
    primitive_units = {n_head: s.Matrix([-4, 0]), chi_head: dimensionless,
                       f_head: dimensionless, g_head: dimensionless}
    primitive_units.update((head, s.Matrix([1, -1])) for head in u_heads)
    primitive_units.update((head, length_unit) for head in face_heads)
    primitive_units.update((head, s.Matrix([0, -1])) for head in rate_heads)
    unknown_units = s.symbols("j_length j_time")
    all_j_heads = j_heads + [head for row in region_j_heads + face_j_heads for head in row]
    primitive_units.update((head, s.Matrix(unknown_units)) for head in all_j_heads)
    dimensions = Dimensions(
        coordinates, primitive_units, [time_unit] + [length_unit] * len(spatial)
    )
    p.emit("dimension_basis", ("length", "time", "number count dimensionless"))
    p.obj("primitive_function_dimensions", primitive_units)
    terms = s.Add.make_args(expanded_order)
    vectors = [dimensions.of(term) for term in terms]
    reference_unit = dimensions.of(ordered.rhs)
    p.emit("order_flux_dimension_inference/operands", (list(zip(terms, vectors)), (ordered.rhs, reference_unit)))
    residuals = [v - reference_unit for v in vectors]
    p.emit("order_flux_dimension_inference/residual", residuals)
    solution = s.linsolve([r for v in residuals for r in v], unknown_units)
    p.obj("order_flux_dimension_inference/solutions", solution)
    unit_rules = dict(zip(unknown_units, next(iter(solution), unknown_units)))
    flux_unit = s.Matrix(unknown_units).subs(unit_rules)
    p.obj("order_flux_dimension_inference/object", flux_unit)
    primitive_units.update((head, flux_unit) for head in all_j_heads)
    dimensions.report(p, "supplied/material_balance_lhs", material.lhs, s.diff(density, time))
    p.obj("supplied/material_balance_rhs", material.rhs)
    p.emit("supplied/material_balance_rhs/dimension", "additive identity has no intrinsic unit")
    dimensions.report(p, "supplied/order_balance_lhs", ordered.lhs, ordered.rhs)
    dimensions.report(p, "supplied/order_balance_rhs", ordered.rhs)
    dimensions.report(p, "supplied/conversion_identity_rhs", conversion.rhs, conversion.lhs)

    # Dirac factors and their normalization are generated by differentiating
    # the actual U1 steps, with one slot retained for each supplied face.
    delta_factors = [
        s.expand(s.diff(step, w) / s.diff(offset, w), diracdelta=True, wrt=w)
        for step, offset in zip(steps, offsets)
    ]
    for label, order_field, regional_orders in packages:
        p.working(label + "/source")
        p.obj(label + "/order_field", order_field)
        if label == "conditional_fast_ordering":
            p.emit(label + "/premise", "conditional specialization: ordering keeps pace with the slab; G is general in both offsets")
            p.obj(label + "/interior_range", s.And(conditional_interior >= 0, conditional_interior <= 1))
        source_operand = transport_expression.subs(chi, order_field).subs(flux_substitution, simultaneous=True)
        p.obj(label + "/flux_substitution_operand", source_operand)
        source = s.expand(source_operand.doit(), diracdelta=True, wrt=w)
        p.obj(label + "/source_distribution", source)
        regular = s.expand(source - s.Add(*supported_sources)).xreplace(dict.fromkeys(delta_factors, s.S.Zero))
        p.obj(label + "/regular_distribution_part", regular)
        dimensions.report(p, label + "/order_field", order_field)
        dimensions.report(p, label + "/source_distribution", source, ordered.rhs)
        dimensions.report(p, label + "/regular_distribution_part", regular, ordered.rhs)
        for k, (domain, regional_order, regional_j) in enumerate(zip(domains, regional_orders, regional_flux), 1):
            name = label + f"/region/{k}"
            p.working(name)
            p.obj(name + "/domain", domain)
            p.obj(name + "/order_operand", regional_order)
            p.obj(name + "/flux_operand", regional_j)
            restricted = s.expand(transport_expression.subs(chi, regional_order).subs(dict(zip(flux, regional_j)), simultaneous=True).doit())
            p.obj(name + "/source", restricted)
            dimensions.report(p, name + "/source", restricted, ordered.rhs)
        face_terms = []
        for k, (face, delta) in enumerate(zip(faces, delta_factors), 1):
            name = label + f"/face/{k}"
            p.working(name)
            coefficient = source.coeff(delta)
            jacobian = s.Abs(s.diff(delta.args[0], w))
            trace = s.expand(coefficient.subs(w, face) / jacobian)
            normalized_delta = s.DiracDelta(
                s.expand(delta.args[0] / s.diff(delta.args[0], w))
            )
            term = s.Add(s.Mul(trace, normalized_delta), supported_sources[k - 1])
            p.obj(name + "/locus", s.Eq(w, face))
            p.obj(name + "/delta_operand", delta)
            p.obj(name + "/normalization_operand", jacobian)
            p.obj(name + "/coefficient_operand", coefficient)
            p.obj(name + "/trace_coefficient", trace)
            p.obj(name + "/supported_flux_operand", supported_flux[k - 1])
            p.obj(name + "/supported_flux_source", supported_sources[k - 1])
            p.obj(name + "/source_contribution", term)
            dimensions.report(p, name + "/trace_coefficient", trace)
            dimensions.report(p, name + "/supported_flux_source", supported_sources[k - 1], ordered.rhs)
            dimensions.report(p, name + "/source_contribution", term, ordered.rhs)
            face_terms.append(term)
        p.obj(label + "/source_with_face_traces", s.Add(regular, *face_terms))
    p.emit("completion", ("general_order", "conditional_fast_ordering", "dimensional_objects"))


if __name__ == "__main__":
    printer = Printer(TRANSCRIPT)
    try:
        main(printer)
    finally:
        printer.stream.close()
