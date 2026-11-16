"""Linear response of the U1 tangential slip option on one boundary component."""

from pathlib import Path

import sympy as sp
import yaml


def main():
    source = Path(__file__).resolve().parents[2] / "em_charge_attribute"
    body = yaml.safe_load((source / "u1_body_dynamics_inputs.yaml").read_text())
    mechanics = yaml.safe_load((source / "u1_body_mechanics_inputs.yaml").read_text())
    endpoint = mechanics["endpoint_functionals"]["E5"]
    normal, tangent = body["endpoints"]["E5"]["condition"].split("; ")
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
    balance = sp.simplify(operator.subs(v, trace) - T)
    power_balance = sp.simplify(powers + 2 * loss)
    print("surface =", endpoint["rayleigh_surface"], " normal condition =", normal)
    print("tangential operator =", operator, " forcing =", T)
    print("Rayleigh density =", loss, " summed force power =", powers)
    print("solved trace =", trace, " equation residual =", balance)
    print("force on surface =", surface_force, " solved force =", force)
    print("zero forcing: trace per V =", sp.factor(trace.subs(T, 0) / V),
          " force per V =", sp.factor(force.subs(T, 0) / V))
    print("surface force slopes: A =", A, " B =", B, " A+B =", sp.factor(A+B))
    print("opposite-force residual =", sp.simplify(surface_force + force_on_fluid),
          " power residual =", power_balance)
    if balance != 0 or power_balance != 0 or sp.simplify(surface_force + force_on_fluid) != 0:
        raise ArithmeticError("Surface trace, force or power balance failed")


if __name__ == "__main__":
    main()
