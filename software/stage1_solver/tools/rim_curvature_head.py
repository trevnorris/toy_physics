"""Virtual work at a spherical rim and steady power for the bare 0PN order layer."""

import sympy as sp


def verify(label, value):
    residual = sp.factor(value)
    print(f"{label} = {residual}", flush=True)
    if residual != 0:
        raise AssertionError(label)


def spherical_rim():
    r, sigma, tau = sp.symbols("a sigma tau", positive=True)
    dr, dw = sp.symbols("r_s w_s", real=True)
    speed = sp.sqrt(dr**2 + dw**2)
    area = 4 * sp.pi * r**2
    sheet_energy_per_parameter = sigma * area * speed
    sheet = sp.Matrix([sp.diff(sheet_energy_per_parameter, d) / area for d in (dr, dw)])
    collar = sp.diff(tau * area, r) / area
    print("sheet endpoint energy gradient per rim area =", list(sheet), flush=True)
    print("collar radial energy gradient per rim area =", collar, flush=True)
    verify("radial endpoint normalization", sheet[0] - sigma * dr / speed)
    verify("axial endpoint normalization", sheet[1] - sigma * dw / speed)
    reverse = sheet.xreplace({dr: -dr, dw: -dw})
    verify("smooth equal-tension join, summed endpoint gradients squared", (sheet + reverse).dot(sheet + reverse))
    verify("collar curvature coefficient", collar - 2 * tau / r)
    vr, vw = sp.symbols("v_rim_r v_rim_w", real=True)
    rim_power = -area * ((sheet[0] + collar) * vr + sheet[1] * vw)
    verify("mechanical power on a stationary rim", rim_power.subs({vr: 0, vw: 0}))


def steady_layer():
    z = sp.Symbol("z", real=True)
    q = sp.Function("q")(z)
    n, area, mob = [sp.Function(s, positive=True)(z) for s in ("n", "A", "M")]
    k, rate = sp.symbols("k I", real=True)
    f = sp.Function("f")(q)
    qp = q.diff(z)
    u = rate / (n * area)
    # Vary the order energy with its actual normal-tube volume measure.
    integrand = area * (f + k * qp**2 / 2)
    mu = (integrand.diff(q) - integrand.diff(qp).diff(z)) / area
    normal_stress = k * qp**2 / 2 - f
    H = area.diff(z) / area
    drive = normal_stress.diff(z) + k * H * qp**2
    verify("curvature plus endpoint stress gradient plus material order force", drive + mu * qp)

    # For a steady layer that only relaxes, D_t q = -M mu. Carrier force is mu q'.
    force = (-u * qp / mob) * qp
    print("relaxation contribution to material power per tube length =", sp.factor(area * u * force), flush=True)
    Rq = qp**2 / (area * n**2 * mob)
    verify("power plus local dissipation", area * u * force + rate**2 * Rq)
    verify("carrier head derivative plus I times resistance density", force / n + rate * Rq)

    # The endpoint term by itself depends on where the tube is cut. In equilibrium
    # d(normal stress)/dz = -k H q'^2 for every smooth stationary extremum.
    equilibrium_qpp = f.diff(q) / k - H * qp
    verify("static curvature and endpoint cancellation",
           drive.subs(q.diff(z, 2), equilibrium_qpp))


if __name__ == "__main__":
    spherical_rim()
    steady_layer()
