import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def mass_engineering_sweep():
    """
    Simulates the time-averaged density at the throat (Effective Mass)
    under 'Suction Mode' drive conditions.
    Checks if we are creating 'Mass Reduction' (Inertial Damping) or
    'Mass Amplification' (Artificial Gravity).
    """

    # Physics Parameters
    n_poly = 5.0
    K_eos = 1.0 / n_poly
    a = 1.0
    rho_inf = 1.0

    # Drive Frequency: f=0.1 (The "Suction" regime found earlier)
    freq = 0.1
    w = 2.0 * np.pi * freq

    # Amplitudes to sweep (Normalized to flux_base)
    # flux_base was -0.1. We sweep drive from weak to strong.
    # drive_amp_ratio = u_osc / u_steady
    amp_ratios = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0]

    results_rho = []
    results_force = []

    print(f"--- Mass Engineering Simulation (n={int(n_poly)}, f={freq}) ---")
    print(f"{'Drive Amp':<12} | {'Mean Density (M_eff)':<22} | {'Mean Force'}")
    print("-" * 55)

    for ratio in amp_ratios:
        rho_mean, force_mean = simulate_throat_density(freq, ratio, n_poly, K_eos, a)

        results_rho.append(rho_mean)
        results_force.append(force_mean)

        # M_eff relative to M_0 (rho_inf)
        m_eff_ratio = rho_mean / rho_inf
        print(f"{ratio:<12.1f} | {rho_mean:<22.4f} ({m_eff_ratio*100:.1f}%) | {force_mean:.4f}")

    # --- Plotting ---
    fig, ax1 = plt.subplots(figsize=(10, 6))

    color = 'tab:red'
    ax1.set_xlabel('Drive Amplitude (Ratio to Steady Flow)')
    ax1.set_ylabel('Effective Mass ($M/M_0$)', color=color)
    ax1.plot(amp_ratios, np.array(results_rho)/rho_inf, 'o-', color=color, linewidth=2, label='Effective Mass')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.axhline(1.0, color='gray', linestyle='--') # Baseline Mass
    ax1.grid(True)

    # Twin axis for Force
    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel('Net Force (Suction/Lift)', color=color)
    ax2.plot(amp_ratios, results_force, 's--', color=color, label='Force')
    ax2.tick_params(axis='y', labelcolor=color)

    plt.title(f"Mass Engineering: Density & Force vs Drive Amplitude\n(Frequency f={freq} c/a)")
    fig.tight_layout()
    plt.show()

def simulate_throat_density(drive_freq, amp_ratio, n_poly, K_eos, a):
    # Grid
    nr = 100
    r_max = 10.0 * a
    r = np.linspace(a, r_max, nr)
    dr = r[1] - r[0]

    # Initial Conditions
    flux_base = -0.1
    u_init = flux_base / (r**2)
    rho_init = np.ones(nr)
    y0 = np.concatenate([rho_init, u_init])

    # Drive
    w = 2.0 * np.pi * drive_freq
    # Amplitude is relative to the base inflow velocity at the throat
    u_steady_throat = flux_base / a**2
    drive_amp = amp_ratio * abs(u_steady_throat)

    def boundary_vel(t):
        return u_steady_throat + drive_amp * np.sin(w * t)

    def pde_system(t, y):
        rho = y[:nr]
        u   = y[nr:]

        # BCs
        u[0] = boundary_vel(t)
        u[-1] = u[-2]
        rho[-1] = 1.0

        # Gradients
        drho = np.gradient(rho, dr, edge_order=2)
        du   = np.gradient(u, dr, edge_order=2)
        drho[0] = (rho[1] - rho[0])/dr
        du[0]   = (u[1] - u[0])/dr

        # Physics
        cs_sq = n_poly * K_eos * rho**(n_poly - 1)

        dt_rho = - (u * drho + rho * du + 2*rho*u/r)
        dt_u   = - (u * du + (cs_sq / rho) * drho)

        # Viscosity
        visc = 0.02 * dr
        dt_rho[1:-1] += visc * np.diff(rho, 2)/dr**2
        dt_u[1:-1]   += visc * np.diff(u, 2)/dr**2

        return np.concatenate([dt_rho, dt_u])

    # Run Simulation
    T = 1.0 / drive_freq
    # Run long enough to settle
    t_end = 4.0 * T
    t_span = (0, t_end)
    t_eval = np.linspace(3.0*T, 4.0*T, 200) # Average last cycle

    # Use max_step to capture wave
    sol = solve_ivp(pde_system, t_span, y0, t_eval=t_eval, method='RK23', max_step=T/50)

    # Analyze Density & Force
    rhos = []
    forces = []

    # Baseline for Force (Steady State)
    rho_base = 1.0 # approx
    u_base = u_steady_throat
    P_base = K_eos * rho_base**n_poly
    F_base = (P_base + rho_base * u_base**2) * (4*np.pi*a**2)

    for i in range(len(sol.t)):
        rho_val = sol.y[0, i]
        u_val   = sol.y[nr, i]

        rhos.append(rho_val)

        P_wall = K_eos * rho_val**n_poly
        F_inst = (P_wall + rho_val * u_val**2) * (4*np.pi*a**2)
        forces.append(F_inst - F_base)

    mean_rho = np.mean(rhos)
    mean_force = np.mean(forces)

    return mean_rho, mean_force

if __name__ == "__main__":
    mass_engineering_sweep()
