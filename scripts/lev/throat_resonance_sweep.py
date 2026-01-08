import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def run_frequency_sweep():
    """
    Sweeps drive frequencies to determine the scaling law of Dark Thrust.
    Also looks for resonant enhancement in the n=5 superfluid.
    """

    # --- Physics Parameters ---
    n_poly = 5.0
    K_eos = 1.0 / n_poly  # Normalized stiffness (c=1, rho=1)
    a = 1.0               # Throat radius

    # We will test normalized frequencies relative to f_natural = c/a
    # f_norm = 1.0 corresponds to the Gamma Ray limit (~10^23 Hz)
    frequencies = [0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]

    mean_lifts = []

    print(f"Starting Resonance Sweep for n={n_poly} Superfluid...")
    print(f"{'Freq (norm)':<12} | {'Mean Lift (norm)':<20} | {'Status'}")
    print("-" * 50)

    for freq in frequencies:
        # Run the simulation for this frequency
        lift = simulate_throat_lift(freq, n_poly, K_eos, a)
        mean_lifts.append(lift)
        print(f"{freq:<12.2f} | {lift:<20.6f} | {'Done'}")

    # --- Analysis & Extrapolation ---
    freqs = np.array(frequencies)
    lifts = np.array(mean_lifts)

    # Fit Power Law: Lift = k * f^alpha
    # log(Lift) = log(k) + alpha * log(f)
    log_f = np.log10(freqs)
    log_L = np.log10(lifts)

    # Linear regression
    coeffs = np.polyfit(log_f, log_L, 1)
    alpha = coeffs[0] # The exponent
    k_const = 10**coeffs[1]

    print("-" * 50)
    print(f"SCALING LAW DETECTED: Lift ~ freq^{alpha:.4f}")

    # Extrapolation to "Non-Gamma" frequencies
    # Target: X-ray (10^-4 norm) or Optical (10^-8 norm)?
    f_target = 1e-4
    lift_projected = k_const * (f_target ** alpha)

    print(f"Projected Lift at f={f_target} (X-ray range): {lift_projected:.6e}")

    # --- Plotting ---
    plt.figure(figsize=(10, 6))
    plt.loglog(freqs, lifts, 'bo-', label=f'Simulation Data (n={int(n_poly)})')

    # Plot Trendline
    f_smooth = np.logspace(np.log10(min(freqs)), np.log10(max(freqs)), 100)
    plt.loglog(f_smooth, k_const * f_smooth**alpha, 'r--', label=f'Scaling: $\propto f^{{{alpha:.2f}}}$')

    plt.title(f"Dark Thrust Frequency Response\nScaling Law: Lift $\propto f^{{{alpha:.2f}}}$")
    plt.xlabel("Frequency (Normalized to $c/a$)")
    plt.ylabel("Mean Lift Force (Normalized)")
    plt.grid(True, which="both", ls="-")
    plt.legend()
    plt.show()

def simulate_throat_lift(drive_freq_Hz, n_poly, K_eos, a):
    """
    Solves the 1D Radial Euler equations for a specific drive frequency.
    Returns the time-averaged Lift Force.
    """
    # Grid Setup
    nr = 100
    r_max = 10.0 * a
    r = np.linspace(a, r_max, nr)
    dr = r[1] - r[0]

    # Initial Conditions (Steady Inflow approx)
    rho_inf = 1.0
    flux_base = -0.1
    u_init = flux_base / (r**2)
    rho_init = np.ones(nr) * rho_inf
    y0 = np.concatenate([rho_init, u_init])

    # Drive Function
    w = 2.0 * np.pi * drive_freq_Hz
    drive_amp = 0.2 * abs(flux_base) / a**2 # 20% modulation

    # Solver Definitions
    def boundary_velocity(t):
        return (flux_base / a**2) + drive_amp * np.sin(w * t)

    def time_derivative(t, y):
        rho = y[:nr]
        u   = y[nr:]

        # BCs
        u[0] = boundary_velocity(t)
        u[-1] = u[-2] # Neumann
        rho[-1] = 1.0 # Fixed far-field

        # Spatial Derivatives
        drho = np.gradient(rho, dr, edge_order=2)
        du   = np.gradient(u, dr, edge_order=2)

        # One-sided correction at throat (index 0)
        drho[0] = (rho[1] - rho[0]) / dr
        du[0]   = (u[1] - u[0]) / dr

        # Pressure Term: grad(P)/rho = n*K*rho^(n-2) * grad(rho)
        # Note: c_s^2 = n*K*rho^(n-1)
        # term = c_s^2 / rho * drho
        cs_sq = n_poly * K_eos * rho**(n_poly - 1)
        press_term = (cs_sq / rho) * drho

        # PDEs
        # dt_rho = - (u*drho + rho*du + 2*rho*u/r)
        dt_rho = - (u * drho + rho * du + (2 * rho * u / r))

        # dt_u = - (u*du + press_term)
        dt_u = - (u * du + press_term)

        # Numerical Dissipation (Viscosity) to capture shock rectification
        visc = 0.005 * dr
        dt_rho[1:-1] += visc * (rho[2:] - 2*rho[1:-1] + rho[:-2]) / dr**2
        dt_u[1:-1]   += visc * (u[2:]   - 2*u[1:-1]   + u[:-2])   / dr**2

        return np.concatenate([dt_rho, dt_u])

    # Run Simulation
    # Ensure we simulate enough cycles to settle (at least 5 cycles)
    period = 1.0 / drive_freq_Hz
    t_end = 10.0 * period
    # To save time in sweep, we analyze the last 3 cycles

    t_span = (0, t_end)
    t_eval = np.linspace(0, t_end, 200 * int(t_end)) # High res for averaging

    sol = solve_ivp(time_derivative, t_span, y0, t_eval=t_eval, method='RK23', rtol=1e-3)

    # Calculate Lift Profile
    # Force ~ P_wall + Momentum_flux
    forces = []

    # Use last 40% of data for steady state average
    idx_start = int(len(sol.t) * 0.6)

    for i in range(idx_start, len(sol.t)):
        rho_val = sol.y[0, i] # rho at throat
        u_val   = sol.y[nr, i] # u at throat

        P_wall = K_eos * rho_val**n_poly
        flux = rho_val * u_val**2
        forces.append(P_wall + flux)

    mean_force = np.mean(forces)

    # Subtract baseline (approximate vacuum pressure 1/n + flux)
    baseline = (K_eos * 1.0**n_poly) + (1.0 * (flux_base/a**2)**2)

    return mean_force - baseline

if __name__ == "__main__":
    run_frequency_sweep()
