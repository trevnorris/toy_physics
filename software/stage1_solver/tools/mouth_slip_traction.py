"""Linear response of the U1 tangential slip option on one boundary component."""

from pathlib import Path

import sympy as sp
import yaml


def main():
    source = Path(__file__).resolve().parents[2] / "em_charge_attribute"
    body = yaml.safe_load((source / "u1_body_dynamics_inputs.yaml").read_text())
    mechanics = yaml.safe_load((source / "u1_body_mechanics_inputs.yaml").read_text())
    endpoint = mechanics["endpoint_functionals"]["E5"]
    wall = body["endpoints"]["E5"]
    normal, tangent = wall["condition"].split("; ")
    lhs, rhs = tangent.replace("v.tangent", "v").replace("V.tangent", "V").split("=")
    v, V, U, T, T_U, T_V = sp.symbols("v V U T T_U T_V", real=True)
    K, gamma = sp.symbols("K gamma", positive=True)
    names = {"tangentDtN": K, "gammaSigma": gamma, "v": v, "V": V,
             "v_tangent": v, "V_tangent": V}
    operator = sp.sympify(lhs, locals=names) - sp.sympify(rhs, locals=names)
    loss = sp.sympify(endpoint["rayleigh_integrand"], locals=names)
    # T is the omitted inhomogeneous exterior traction in the scalar equation.
    # It stays independent of the physical incident speed until the last step.
    trace, = sp.solve(operator - T, v)
    surface_force = -sp.diff(loss, V)
    force = sp.factor(surface_force.subs(v, trace))
    linear = force.subs(T, T_U * U + T_V * V)
    A, B = (sp.factor(sp.diff(linear, speed)) for speed in (U, V))
    force_on_fluid = -sp.diff(loss, v)
    powers = sp.factor(surface_force * V + force_on_fluid * v)
    zero_trace = sp.factor(trace.subs(T, 0) / V)
    zero_force = sp.factor(force.subs(T, 0) / V)
    slope_sum = sp.factor(A + B)

    # Independently assemble the traction balance from the mechanics input.
    exterior = sp.sympify(mechanics["boundary_operator"]["tangent_DtN"], locals=names)
    stationarity = exterior * v + sp.diff(loss, v)
    matrix_data = wall["trace_system"]
    matrix = sp.Matrix([[sp.sympify(cell, locals=names) for cell in row]
                        for row in matrix_data["matrix"]])
    rhs_vector = sp.Matrix([sp.sympify(cell, locals=names) for cell in matrix_data["rhs"]])
    reference_trace = matrix.LUsolve(rhs_vector)[wall["trace_unknowns"].index("tangent")]
    exterior_force = T - exterior * trace
    # Differentiate the stationarity equation implicitly along each preparation.
    def exterior_slope(speed, applied_slope):
        trace_slope = (applied_slope - sp.diff(stationarity, speed)) / sp.diff(stationarity, v)
        return sp.factor(applied_slope - exterior * trace_slope)

    incident = exterior_slope(U, T_U)
    translation = exterior_slope(V, T_V)
    domain_checks = {
        "Rayleigh surface": (endpoint["rayleigh_surface"], "Sigma"),
        "normal motion on wall": (normal, "v.normal=V.normal"),
        "normal trace imposed": (endpoint["essential"].get("normal"), 1),
    }
    checks = {
        "tangential operator": (operator, stationarity),
        "forced trace in traction equation": (operator.subs(v, trace), T),
        "wall force before solving": (surface_force, -force_on_fluid),
        "solved wall force": (force, exterior_force),
        "wall and material power": (powers, -2 * loss),
        "zero forcing trace per V": (zero_trace, reference_trace),
        "zero forcing wall force per V": (zero_force, -exterior * reference_trace),
        "wall incident slope": (A, incident),
        "wall translation slope": (B, translation),
        "wall slope sum": (slope_sum, incident + translation),
    }
    print("Rayleigh density =", loss, " applied traction =", T, " solved trace =", trace)
    failures = []
    for label, (actual, expected) in domain_checks.items():
        print(f"{label}: supplied={actual} required={expected}")
        if actual != expected:
            failures.append(label)
    for label, (actual, expected) in checks.items():
        residual = sp.factor(actual - expected)
        print(f"{label}: result={sp.factor(actual)} reference={sp.factor(expected)} residual={residual}")
        if residual != 0:
            failures.append(label)
    if failures:
        raise ArithmeticError("Inconsistent E5 wall calculation: " + ", ".join(failures))


if __name__ == "__main__":
    main()
