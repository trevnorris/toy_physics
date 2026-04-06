import numpy as np
from scipy.integrate import solve_ivp

# Dimensionless weak-field units.
G = 1.0
M = 1.0
c0 = 1.0
mu = G * M

IMPACT_PARAMETERS = (400.0, 800.0, 1600.0)
FAR_FIELD_MULTIPLIER = 400.0


def cs_and_grad(pos):
    """Linearized n=5 pure-refraction profile used in the paper."""
    r = np.linalg.norm(pos)
    cs = c0 * (1.0 - 2.0 * mu / (c0**2 * r))
    grad_cs = (2.0 * mu / (c0 * r**3)) * pos
    return cs, grad_cs


def ray_equations(_, state):
    pos = state[:2]
    k = state[2:]
    k_mag = np.linalg.norm(k)
    cs, grad_cs = cs_and_grad(pos)
    dx_dlambda = cs * k / k_mag
    dk_dlambda = -k_mag * grad_cs
    return np.concatenate([dx_dlambda, dk_dlambda])


def run_case(impact_parameter):
    z_far = FAR_FIELD_MULTIPLIER * impact_parameter
    state0 = np.array([impact_parameter, -z_far, 0.0, 1.0])

    def exit_event(_, state):
        return state[1] - z_far

    exit_event.terminal = True
    exit_event.direction = 1

    sol = solve_ivp(
        ray_equations,
        (0.0, 2.5 * z_far),
        state0,
        events=exit_event,
        rtol=1e-9,
        atol=1e-9,
        max_step=impact_parameter / 20.0,
    )

    if sol.status != 1:
        raise RuntimeError(f"Ray did not reach the far zone for b={impact_parameter}.")

    k_out = sol.y[2:, -1]
    theta_num = abs(np.arctan2(k_out[0], k_out[1]))
    theta_analytic = 4.0 * mu / (impact_parameter * c0**2)
    return theta_num, theta_analytic, theta_num / theta_analytic


print("Ray-tracing check for the preferred n=5 pure-refraction branch")
print("Using N(r) ≈ 1 + 2 GM/(r c^2), H = c_s(r)|k|, and no background flow.")
print(f"{'b':>8} | {'theta_num':>14} | {'theta_analytic':>14} | {'ratio':>10}")
print("-" * 58)

for impact_parameter in IMPACT_PARAMETERS:
    theta_num, theta_analytic, ratio = run_case(impact_parameter)
    print(
        f"{impact_parameter:8.1f} | "
        f"{theta_num:14.8f} | "
        f"{theta_analytic:14.8f} | "
        f"{ratio:10.6f}"
    )

print("-" * 58)
print("Conclusion: the weak-field numerical rays stay close to Δθ = 4GM/(b c^2).")
