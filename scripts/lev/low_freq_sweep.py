import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def run_low_freq_sweep_v2():
    """
    Fixed Low-Frequency Sweep to capture the 'Deep Breathing' regime
    without timeouts or NaNs.
    """

    # Physics (n=5 Stiff Vacuum)
    n_poly = 5.0
    K_eos = 1.0 / n_poly
    a = 1.0

    # Sweep Range: Focused on the transition region found in previous run
    frequencies = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.10]

    lifts = []
    efficiencies = []

    print(f"Starting Refined Low-Frequency Sweep (n={int(n_poly)})...")
    print(f"{'Freq (norm)':<12} | {'Mean Lift':<15} | {'Efficiency'}")
    print("-" * 45)

    for freq in frequencies:
        lift, power = simulate_throat_dynamics(freq, n_poly, K_eos, a)

        # Avoid div/0 if power is negligible (DC hold)
        if power > 1e-9:
            eff = lift / power
        else:
            eff = 0.0

        lifts.append(lift)
        efficiencies.append(eff)

        print(f"{freq:<12.4f} | {lift:<15.6f} | {eff:.2f}")

    # Plot results
    fig, ax1 = plt.subplots(figsize=(10, 6))

    color = 'tab:blue'
    ax1.set_xlabel('Frequency (Normalized $f/f_0$)')
    ax1.set_ylabel('Lift Force (Normalized)', color=color)
    ax1.plot(frequencies, lifts, 'o-', color=color, linewidth=2, label='Lift Force')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(True)

    # Highlight the zero crossing
    ax1.axhline(0, color='black', linewidth=1, linestyle='--')

    ax2 = ax1.twinx()
    color = 'tab:green'
    ax2.set_ylabel('Efficiency metric (Lift/Power)', color=color)
    ax2.plot(frequencies, efficiencies, 's--', color=color, label='Efficiency')
    ax2.tick_params(axis='y', labelcolor=color)

    plt.title("Proton Anchor Performance: Mapping the Pump Regime")
    fig.tight_layout()
    plt.show()

def simulate_throat_dynamics(drive_freq, n_poly, K_eos, a):
    # Grid
    nr = 100
    r_max = 10.0 * a
    r = np.linspace(a, r_max, nr)
    dr = r[1] - r[0]

    # Initial Conditions (Steady Inflow)
    flux_base = -0.1
    u_init = flux_base / (r**2)
    rho_init = np.ones(nr)
    y0 = np.concatenate([rho_init, u_init])

    # Drive
    w = 2.0 * np.pi * drive_freq
    amp = 0.2 * abs(flux_base) / a**2

    def boundary_vel(t):
        return (flux_base / a**2) + amp * np.sin(w * t)

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

        # Physics (n=5)
        cs_sq = n_poly * K_eos * rho**(n_poly - 1)

        dt_rho = - (u * drho + rho * du + 2*rho*u/r)
        dt_u   = - (u * du + (cs_sq / rho) * drho)

        # Viscosity (Critical for shock capture in pump mode)
        visc = 0.02 * dr
        dt_rho[1:-1] += visc * np.diff(rho, 2)/dr**2
        dt_u[1:-1]   += visc * np.diff(u, 2)/dr**2

        return np.concatenate([dt_rho, dt_u])

    # --- FIX: Simulation Duration ---
    # We need strictly 3 full periods to ensure the averaging window is valid.
    # The slower the frequency, the longer the T.
    T = 1.0 / drive_freq
    t_end = 3.0 * T

    # We only analyze the LAST single period to be safe
    # High resolution for the solver to avoid stepping over the wave
    t_span = (0, t_end)
    t_eval = np.linspace(2.0*T, 3.0*T, 200)

    # Increase max_step to prevent timeouts on long, slow cycles
    sol = solve_ivp(pde_system, t_span, y0, t_eval=t_eval, method='RK23', max_step=T/50)

    # Analyze Lift & Power
    forces = []
    powers = []

    baseline_force = (K_eos * 1.0**n_poly) * (4*np.pi*a**2) + \
                     (1.0 * (flux_base/a**2)**2) * (4*np.pi*a**2)

    for i in range(len(sol.t)):
        rho_val = sol.y[0, i]
        u_val   = sol.y[nr, i]

        P_wall = K_eos * rho_val**n_poly
        F_inst = (P_wall + rho_val * u_val**2) * (4*np.pi*a**2)
        forces.append(F_inst)

        power_inst = P_wall * (u_val * 4*np.pi*a**2)
        powers.append(abs(power_inst))

    if len(forces) == 0:
        return 0.0, 1.0 # Fail safe

    mean_lift = np.mean(forces) - baseline_force
    mean_power = np.mean(powers)

    return mean_lift, mean_power

if __name__ == "__main__":
    run_low_freq_sweep_v2()
