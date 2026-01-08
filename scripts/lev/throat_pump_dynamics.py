import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def throat_pump_simulation():
    """
    Simulates the nonlinear 1D radial flow of an n=5 superfluid
    driven by an oscillating boundary at the throat (r=a).
    Calculates the net 'Dark Thrust' (Reaction Force) over time.
    """

    # ==========================================
    # 1. Physics & Grid Configuration
    # ==========================================
    # Parameters for the Stiff Vacuum
    n_poly = 5.0          # Stiff polytrope index (Paper VI)
    gamma_eos = n_poly    # Adiabatic index

    # Normalized Units:
    # Throat Radius a = 1.0
    # Background Density rho_inf = 1.0
    # Sound Speed c_inf = 1.0
    # Mass M ~ 1.0
    a = 1.0
    rho_inf = 1.0
    c_inf = 1.0

    # Stiffness constant K: c^2 = n * K * rho^(n-1)
    # 1^2 = 5 * K * 1^4  -> K = 1/5
    K_eos = c_inf**2 / n_poly

    # Grid (Radial coordinate r)
    nr = 100
    r_max = 20.0 * a
    r = np.linspace(a, r_max, nr)
    dr = r[1] - r[0]

    # ==========================================
    # 2. Initial State (Steady Inflow)
    # ==========================================
    # We start with a mild steady inflow (gravity) to establish the baseline.
    # Conservation of mass: u(r) ~ 1/r^2
    flux_base = -0.05  # Negative = Inflow

    rho_init = np.ones(nr) * rho_inf
    u_init = flux_base / (r**2)

    # State Vector: [Density_0, ..., Density_N, Velocity_0, ..., Velocity_N]
    y0 = np.concatenate([rho_init, u_init])

    # ==========================================
    # 3. The Pumping Drive (Actuation)
    # ==========================================
    # We modulate the velocity at the throat: u(a, t) = u_base + u_pump * sin(wt)
    drive_freq = 2.0 * np.pi * 1.0  # f = 1.0 (Normalized)
    drive_amp = 0.5 * np.abs(flux_base) / (a**2) # 50% modulation amplitude

    def boundary_velocity(t):
        return (flux_base / a**2) + drive_amp * np.sin(drive_freq * t)

    # ==========================================
    # 4. The Euler Solver (Method of Lines)
    # ==========================================
    def time_derivative(t, y):
        # Unpack state
        rho = y[:nr]
        u   = y[nr:]

        # Enforce Boundary Conditions (Ghost point method approx)
        u[0] = boundary_velocity(t)
        # Outer boundary (Non-reflecting / Neumann approx)
        u[-1] = u[-2]
        rho[-1] = rho_inf # Fix far-field density

        # Gradients (Central Difference)
        # d/dr
        drho = np.zeros_like(rho)
        du   = np.zeros_like(u)

        drho[1:-1] = (rho[2:] - rho[:-2]) / (2*dr)
        du[1:-1]   = (u[2:]   - u[:-2])   / (2*dr)

        # One-sided at boundaries for stability
        drho[0] = (rho[1] - rho[0]) / dr
        du[0]   = (u[1]   - u[0])   / dr

        # Pressure Gradient: P = K * rho^n
        # grad(P)/rho = K * n * rho^(n-1) * grad(rho) / rho
        #             = c_s^2(rho) * grad(rho) / rho
        cs_sq = n_poly * K_eos * rho**(n_poly - 1)
        term_pressure = (cs_sq / rho) * drho

        # --- Governing Equations ---
        # 1. Continuity: d_t(rho) = - [ u * d_r(rho) + rho * d_r(u) + 2*rho*u/r ]
        #    (The 2/r term comes from spherical symmetry in 3D/4D)
        dt_rho = - (u * drho + rho * du + 2 * rho * u / r)

        # 2. Momentum: d_t(u) = - [ u * d_r(u) + (1/rho)*d_r(P) ]
        dt_u   = - (u * du + term_pressure)

        # Artificial Viscosity (Smoothing for stability)
        # Adds a small diffusion term: nu * d^2(u)/dr^2
        nu = 0.01 * dr
        dt_u[1:-1] += nu * (u[2:] - 2*u[1:-1] + u[:-2]) / dr**2
        dt_rho[1:-1] += nu * (rho[2:] - 2*rho[1:-1] + rho[:-2]) / dr**2

        return np.concatenate([dt_rho, dt_u])

    # ==========================================
    # 5. Run Simulation
    # ==========================================
    t_span = (0, 20.0) # Run for 20 cycles
    t_eval = np.linspace(t_span[0], t_span[1], 500)

    print("Running solver (this may take a moment)...")
    sol = solve_ivp(time_derivative, t_span, y0, t_eval=t_eval, method='RK45')

    # ==========================================
    # 6. Compute Thrust (Reaction Force)
    # ==========================================
    # Force on Throat Wall (Reaction to Fluid)
    # F = Integral [ P(a) + rho(a)*u(a)^2 ] dA
    # We care about the CHANGE in force relative to the undriven baseline.

    forces = []
    baselines = []

    for i, t in enumerate(sol.t):
        rho_now = sol.y[:nr, i]
        u_now   = sol.y[nr:, i]

        # Pressure at throat wall
        P_wall = K_eos * rho_now[0]**n_poly
        # Momentum flux at throat
        Mom_flux = rho_now[0] * u_now[0]**2

        # Total "Push" into the bulk
        # Area = 4*pi*a^2 (Spherical throat approximation)
        Area = 4 * np.pi * a**2
        F_total = (P_wall + Mom_flux) * Area

        # Calculate what the force WOULD be without the AC drive (Baseline)
        # Baseline u is just the steady flux
        u_base = flux_base / a**2
        rho_base = rho_inf # Approx
        P_base = K_eos * rho_base**n_poly
        F_base = (P_base + rho_base * u_base**2) * Area

        forces.append(F_total)
        baselines.append(F_base)

    forces = np.array(forces)
    baselines = np.array(baselines)

    # Net Thrust = Time Average of (Force - Baseline)
    # Note: If Force > Baseline, we are pushing HARDER into the bulk.
    # Reaction force is UP (Levitation).
    net_thrust_profile = forces - baselines

    # Average over the last 50% of the run (steady state)
    cutoff = len(forces) // 2
    mean_lift = np.mean(net_thrust_profile[cutoff:])

    print(f"\nRESULTS:")
    print(f"Mean 'Dark Thrust' (Lift): {mean_lift:.6f} [Normalized Force]")
    print(f"  (Positive = Levitation, Negative = Suction)")

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(sol.t, net_thrust_profile, label="Instantaneous Lift")
    plt.axhline(mean_lift, color='r', linestyle='--', label=f"Mean Lift = {mean_lift:.4f}")
    plt.title(f"Dark Thrust Validation (n={n_poly})")
    plt.xlabel("Time (Normalized)")
    plt.ylabel("Net Lift Force (Force - Baseline)")
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    throat_pump_simulation()
