#!/usr/bin/env python3
"""S12: SymPy construction from the supplied balances and U1.

No imports of the Wolfram engine, its transcript, or any other derived result.
Revision 4 gives n, u and J independent live restrictions and one-sided traces
in three open regions. Order is live inside and zero outside by U1. Input
fields have no face-concentrated parts. The source and total-material balance
are differentiated distributionally, with conserved products formed region
by region. Material-balance coefficients are extracted locally before their
bulk contributions are removed. No jump or point value is selected; no live
field is frozen. There is no expected-value guard or cross-engine comparison.
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
    p.emit("governing_input", "S12_decision_list.md revision 4, b08673a2; supplied equations and face rules are unfalsifiable inputs")
    p.emit("calculus_scope", (
        "n, u and J have independent live restrictions and one-sided traces in all three open regions",
        "order is a live interior field, zero outside by U1",
        "no input field has a face-concentrated part; those possibilities remain OPEN outside this object",
        "the computed source can have face-concentrated parts",
        "pointwise face values and jumps are unspecified",
    ))
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

    # BEGIN SUPPLIED BALANCES. The spacetime current components are retained
    # as undifferentiated operands, so discontinuous products can be formed
    # separately in their open regions before distributional differentiation.
    material_current = (density,) + tuple(density * v for v in velocity)
    order_current = (chi * density,) + tuple(
        chi * density * v + j for v, j in zip(velocity, flux)
    )
    material = s.Eq(sum(s.diff(c, z) for c, z in zip(material_current, coordinates)), 0)
    ordered = s.Eq(sum(s.diff(c, z) for c, z in zip(order_current, coordinates)), density * rates[0])
    conversion = s.Eq(rates[0], rates[1] - rates[2])
    # END SUPPLIED BALANCES.
    # BEGIN U1 ENCODING.
    offsets = (w - faces[0], faces[1] - w)
    steps = tuple(s.Heaviside(offset) for offset in offsets)
    support = s.Mul(*steps)
    interior = f_head(*coordinates)
    conditional_interior = g_head(*offsets)
    packages = (
        ("general_order", (s.S.Zero, interior, s.S.Zero)),
        ("conditional_fast_ordering", (s.S.Zero, conditional_interior, s.S.Zero)),
    )
    domains = (
        s.And(faces[0] < faces[1], w < faces[0]),
        s.And(faces[0] < faces[1], faces[0] < w, w < faces[1]),
        s.And(faces[0] < faces[1], w > faces[1]),
    )
    # END U1 ENCODING.
    regional_masks = (
        s.Heaviside(s.Mul(-1, offsets[0])), support,
        s.Heaviside(s.Mul(-1, offsets[1])),
    )
    region_n_heads = [s.Function(f"n_region{r + 1}") for r in range(len(domains))]
    region_u_heads = [
        [s.Function(f"u{c + 1}_region{r + 1}") for c in range(len(spatial))]
        for r in range(len(domains))
    ]
    region_j_heads = [
        [s.Function(f"j{c + 1}_region{r + 1}") for c in range(len(spatial))]
        for r in range(len(domains))
    ]
    regional_density = [fun(*coordinates) for fun in region_n_heads]
    regional_velocity = [[fun(*coordinates) for fun in row] for row in region_u_heads]
    regional_flux = [[fun(*coordinates) for fun in row] for row in region_j_heads]
    generic_fields = [density] + velocity + flux
    region_maps = [
        dict(zip(generic_fields, [n] + u + j))
        for n, u, j in zip(regional_density, regional_velocity, regional_flux)
    ]
    # Only these smooth regional products are multiplied by a partition mask.
    mass_currents = [tuple(c.xreplace(mapping) for c in material_current) for mapping in region_maps]
    mass_local = [s.expand(sum(s.diff(c, z) for c, z in zip(row, coordinates))) for row in mass_currents]
    mass_current_field = tuple(
        s.Add(*(mask * row[c] for mask, row in zip(regional_masks, mass_currents)))
        for c in range(len(coordinates))
    )
    p.obj("supplied/material_balance", material)
    p.obj("supplied/order_balance", ordered)
    p.obj("supplied/conversion_identity", conversion)
    p.obj("supplied/faces", faces)
    p.obj("supplied/domain", (faces[0] < faces[1], s.And(interior >= 0, interior <= 1)))
    p.obj("live_fields", (density, velocity, flux, chi, faces))
    p.obj("region_domains", domains)
    p.obj("fields/region_indicators", regional_masks)
    p.obj("density/region_fields", regional_density)
    p.obj("velocity/region_fields", regional_velocity)
    p.obj("flux/region_fields", regional_flux)
    p.obj("density/distribution_field", s.Add(*(mask * n for mask, n in zip(regional_masks, regional_density))))
    p.obj("velocity/distribution_field", [s.Add(*(mask * row[c] for mask, row in zip(regional_masks, regional_velocity))) for c in range(len(spatial))])
    p.obj("flux/distribution_field", [s.Add(*(mask * row[c] for mask, row in zip(regional_masks, regional_flux))) for c in range(len(spatial))])

    p.working("order_flux_dimensions")
    length_unit, time_unit = s.Matrix([1, 0]), s.Matrix([0, 1])
    dimensionless = s.zeros(2, 1)
    primitive_units = {n_head: s.Matrix([-4, 0]), chi_head: dimensionless,
                       f_head: dimensionless, g_head: dimensionless}
    primitive_units.update((head, s.Matrix([-4, 0])) for head in region_n_heads)
    all_u_heads = u_heads + [head for row in region_u_heads for head in row]
    primitive_units.update((head, s.Matrix([1, -1])) for head in all_u_heads)
    primitive_units.update((head, length_unit) for head in face_heads)
    primitive_units.update((head, s.Matrix([0, -1])) for head in rate_heads)
    unknown_units = s.symbols("j_length j_time")
    all_j_heads = j_heads + [head for row in region_j_heads for head in row]
    primitive_units.update((head, s.Matrix(unknown_units)) for head in all_j_heads)
    dimensions = Dimensions(coordinates, primitive_units, [time_unit] + [length_unit] * len(spatial))
    p.emit("dimension_basis", ("length", "time", "number count dimensionless"))
    p.obj("primitive_function_dimensions", primitive_units)
    terms = s.Add.make_args(s.expand(ordered.lhs))
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

    delta_factors = [
        s.expand(s.diff(step, w) / s.diff(offset, w), diracdelta=True, wrt=w)
        for step, offset in zip(steps, offsets)
    ]
    p.working("material/balance_distribution")
    p.obj("material/region_current_operands", mass_currents)
    p.obj("material/region_balance_operands", [s.Eq(lhs, material.rhs) for lhs in mass_local])
    p.obj("material/conserved_current_operand", mass_current_field)
    mass_raw = s.expand(sum(s.diff(c, z) for c, z in zip(mass_current_field, coordinates)), diracdelta=True, wrt=w)
    p.obj("material/balance_distribution_operand", mass_raw)
    mass_bulk_operand = s.expand(s.Add(*(mask * lhs for mask, lhs in zip(regional_masks, mass_local))))
    p.obj("material/bulk_balance_operand", mass_bulk_operand)
    mass_required = s.expand(mass_raw - mass_bulk_operand)
    p.obj("material/balance_distribution", mass_required)
    dimensions.report(p, "material/balance_distribution_operand", mass_raw, s.diff(density, time))
    dimensions.report(p, "material/balance_distribution", mass_required, s.diff(density, time))
    for r, lhs in enumerate(mass_local, 1):
        p.obj(f"material/region/{r}/balance", s.Eq(lhs, material.rhs))
        dimensions.report(p, f"material/region/{r}/balance_lhs", lhs, s.diff(density, time))
    mass_face_equations, mass_face_terms = [], []
    for k, (face, delta) in enumerate(zip(faces, delta_factors), 1):
        name = f"material/face/{k}"
        p.working(name)
        coefficient = mass_raw.coeff(delta)
        jacobian = s.Abs(s.diff(delta.args[0], w))
        trace = s.expand(coefficient.subs(w, face) / jacobian)
        term = trace * s.DiracDelta(s.expand(delta.args[0] / s.diff(delta.args[0], w)))
        p.obj(name + "/coefficient_operand", coefficient)
        p.obj(name + "/trace_coefficient", trace)
        p.emit(name + "/balance_operands", (trace, material.rhs))
        p.emit(name + "/balance_residual", trace - material.rhs)
        equation = s.Eq(trace, material.rhs)
        p.obj(name + "/requirement", equation)
        p.obj(name + "/balance_contribution", term)
        dimensions.report(p, name + "/trace_coefficient", trace)
        dimensions.report(p, name + "/balance_contribution", term, s.diff(density, time))
        mass_face_equations.append(equation)
        mass_face_terms.append(term)
    p.obj("material/face_requirements", mass_face_equations)
    p.obj("material/balance_with_face_traces", s.Add(*mass_face_terms))

    for label, regional_orders in packages:
        p.working(label + "/conserved_operands")
        order_field = s.Add(*(mask * q for mask, q in zip(regional_masks, regional_orders)))
        p.obj(label + "/order_field", order_field)
        p.obj(label + "/region_order_operands", regional_orders)
        if label == "conditional_fast_ordering":
            p.emit(label + "/premise", "conditional specialization: ordering keeps pace with the slab; G is general in both offsets")
            p.obj(label + "/interior_range", s.And(conditional_interior >= 0, conditional_interior <= 1))
        order_currents = []
        for mapping, q in zip(region_maps, regional_orders):
            field_map = dict(mapping)
            field_map[chi] = q
            order_currents.append(tuple(c.xreplace(field_map) for c in order_current))
        p.obj(label + "/region_current_operands", order_currents)
        current_field = tuple(
            s.Add(*(mask * row[c] for mask, row in zip(regional_masks, order_currents)))
            for c in range(len(coordinates))
        )
        p.obj(label + "/conserved_current_operand", current_field)
        # Extract the coefficient of the regional material density's time
        # derivative from each local order equation. Remove only this local
        # bulk balance. Multiplying a global material delta by a jumping chi
        # would choose an undefined distribution product and is never used.
        local_order = [s.expand(sum(s.diff(c, z) for c, z in zip(row, coordinates))) for row in order_currents]
        weights = [lhs.coeff(s.diff(n, time)) for lhs, n in zip(local_order, regional_density)]
        local_multiples = [s.expand(weight * mass) for weight, mass in zip(weights, mass_local)]
        p.obj(label + "/material_extraction/coefficients", weights)
        p.obj(label + "/material_extraction/local_operands", list(zip(local_order, local_multiples)))
        bulk_multiple = s.expand(s.Add(*(mask * part for mask, part in zip(regional_masks, local_multiples))))
        p.obj(label + "/material_elimination_operand", bulk_multiple)
        p.working(label + "/source")
        raw = s.expand(sum(s.diff(c, z) for c, z in zip(current_field, coordinates)), diracdelta=True, wrt=w)
        p.obj(label + "/balance_source_operand", raw)
        source = s.expand(raw - bulk_multiple)
        p.obj(label + "/source_distribution", source)
        regular = source.xreplace(dict.fromkeys(delta_factors, s.S.Zero))
        p.obj(label + "/regular_distribution_part", regular)
        dimensions.report(p, label + "/order_field", order_field)
        dimensions.report(p, label + "/balance_source_operand", raw, ordered.rhs)
        dimensions.report(p, label + "/source_distribution", source, ordered.rhs)
        dimensions.report(p, label + "/regular_distribution_part", regular, ordered.rhs)
        for r, (domain, q, lhs, part) in enumerate(zip(domains, regional_orders, local_order, local_multiples), 1):
            name = label + f"/region/{r}"
            p.working(name)
            p.obj(name + "/domain", domain)
            p.obj(name + "/field_operands", (regional_density[r - 1], regional_velocity[r - 1], regional_flux[r - 1], q))
            restricted = s.expand(lhs - part)
            p.obj(name + "/source", restricted)
            dimensions.report(p, name + "/source", restricted, ordered.rhs)
        face_terms = []
        for k, (face, delta) in enumerate(zip(faces, delta_factors), 1):
            name = label + f"/face/{k}"
            p.working(name)
            coefficient = source.coeff(delta)
            jacobian = s.Abs(s.diff(delta.args[0], w))
            trace = s.expand(coefficient.subs(w, face) / jacobian)
            normalized_delta = s.DiracDelta(s.expand(delta.args[0] / s.diff(delta.args[0], w)))
            term = s.Mul(trace, normalized_delta)
            trace_fields = [
                (regional_density[r].subs(w, face),
                 [v.subs(w, face) for v in regional_velocity[r]],
                 [j.subs(w, face) for j in regional_flux[r]],
                 regional_orders[r].subs(w, face))
                for r in (k - 1, k)
            ]
            p.obj(name + "/locus", s.Eq(w, face))
            p.obj(name + "/field_trace_operands", trace_fields)
            p.obj(name + "/delta_operand", delta)
            p.obj(name + "/normalization_operand", jacobian)
            p.obj(name + "/coefficient_operand", coefficient)
            p.obj(name + "/trace_coefficient", trace)
            p.obj(name + "/source_contribution", term)
            p.obj(name + "/material_requirement", mass_face_equations[k - 1])
            dimensions.report(p, name + "/trace_coefficient", trace)
            dimensions.report(p, name + "/source_contribution", term, ordered.rhs)
            face_terms.append(term)
        p.obj(label + "/source_with_face_traces", s.Add(regular, *face_terms))
    p.emit("completion", ("material_face_requirements", "general_order", "conditional_fast_ordering", "dimensional_objects"))


if __name__ == "__main__":
    printer = Printer(TRANSCRIPT)
    try:
        main(printer)
    finally:
        printer.stream.close()
