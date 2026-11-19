"""Axial force of the existing conservative cylindrical cavity at free length."""

import sympy as sp


def check(label, result, expected):
    print(f"{label} = {sp.factor(result)}", flush=True)
    if sp.simplify(result-expected) != 0:
        raise AssertionError(label)


R, L, action = sp.symbols("R L J", positive=True)
radial, axial, gap = sp.symbols("v_r v_z omega_g", positive=True)
sigma, pressure, curvature = sp.symbols("sigma p k_b", nonnegative=True)
omega = sp.sqrt(gap**2+(radial/R)**2+(axial/L)**2)
wave = action*omega
wall = sigma*(4*sp.pi*R**2*L+8*sp.pi*R**3/3)
volume = pressure*4*sp.pi*R**3*L/3
bending = 16*sp.pi*curvature*L/9
total = wave + wall + volume + bending
push = -wave.diff(L)
beta = (axial/L)**2/omega**2
check("wave axial force divided by wave energy", push/wave, beta/L)

# Only imposed condition is axial stationarity. The two-sided stress and its
# transfer to the exterior are independent of this internal load balance.
j_star, = sp.solve(total.diff(L), action)
force = push.subs(action, j_star)
check("axially balanced wave push", force,
      16*sp.pi*curvature/9+4*sp.pi*sigma*R**2+4*sp.pi*pressure*R**3/3)
check("balanced force slope at common material parameters", force.diff(R),
      8*sp.pi*sigma*R+4*sp.pi*pressure*R**2)
check("doubling radius changes the balanced push by", force.subs(R,2*R)-force,
      12*sp.pi*sigma*R**2+28*sp.pi*pressure*R**3/3)
check("free-radius outward drive with curvature alone",
      -total.diff(R).subs({pressure:0,sigma:0}), action*radial**2/(R**3*omega))
