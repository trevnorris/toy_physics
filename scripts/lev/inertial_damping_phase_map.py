import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def inertial_damping_phase_map():
    """
    Sweeps Frequency vs. Amplitude to map the 'Agility Mode' (Mass Reduction).
    Includes a 'Structural Integrity' check to ensure the particle doesn't collapse.
    """
    print("--- Phase Space Mapping: Inertial Damping vs. Stability ---")

    # Grid Resolution (Keep moderate for speed)
    freqs = np.linspace(0.05, 0.5, 6)   # Normalized Frequency f/f0
    amps  = np.linspace(0.5, 4.0, 6)    # Drive Amplitude Ratio (u_osc / u_steady)

    # Storage
    mass_map = np.zeros((len(amps), len(freqs)))
    stability_map = np.zeros((len(amps), len(freqs)))

    # Physics (n=5 Stiff Vacuum)
    n_poly = 5.0
    K_eos = 1.0 / n_poly
    a = 1.0

    print(f"Sweeping {len(freqs)}x{len(amps)} grid ({len(freqs)*len(amps)} simulations)...")

    for i, amp_ratio in enumerate(amps):
        for j, freq in enumerate(freqs):

            # Run simulation
            m_eff, min_rho = simulate_point(freq, amp_ratio, n_poly, K_eos, a)

            # Store data
            mass_map[i, j] = m_eff
            stability_map[i, j] = min_rho

            # Status update
            status = "STABLE" if min_rho > 0.5 else "CRITICAL" if min_rho > 0.1 else "COLLAPSE"
            print(f"  f={freq:.2f}, Amp={amp_ratio:.1f} -> Mass: {m_eff:.3f}, MinRho: {min_rho:.2f} [{status}]")

    # --- Plotting the Landscape ---
    X, Y = np.meshgrid(freqs, amps)

    fig, ax = plt.subplots(figsize=(10, 7))

    # Contour Plot for Mass (The Terrain)
    # Levels: 0.90 to 1.10 (Mass reduction to Mass increase)
    levels = np.linspace(0.85, 1.15, 15)
    cp = ax.contourf(X, Y, mass_map, levels=levels, cmap='coolwarm', extend='both')
    cbar = fig.colorbar(cp)
    cbar.set_label('Effective Mass Ratio ($M/M_0$)')

    # Overlay: The "neutral" line (Mass = 1.0)
    ax.contour(X, Y, mass_map, levels=[1.0], colors='white', linewidths=2, linestyles='solid')

    # Overlay: Stability Warning (The Cliff)
    # Hatching or red line where min_rho < 0.5 (Risk of destabilizing particle)
    ax.contour(X, Y, stability_map, levels=[0.5], colors='red', linewidths=2, linestyles='dashed')
    ax.text(freqs[0], amps[-1], 'Red Dashed Line = Stability Limit (Throat Collapse)', color='red', fontsize=9, verticalalignment='top')

    ax.set_title("Inertial Damping Phase Space\n(Blue = Mass Reduction, Red = Mass Increase)")
    ax.set_xlabel("Drive Frequency ($f/f_0$)")
    ax.set_ylabel("Drive Amplitude ($u_{osc}/u_{steady}$)")

    plt.show()

def simulate_point(freq, amp_ratio, n_poly, K_eos, a):
    # Reduced grid for speed
    nr = 50
    r_max = 8.0 * a
    r = np.linspace(a, r_max, nr)
    dr = r[1] - r[0]

    # Initial State
    flux_base = -0.1
    u_init = flux_base / (r**2)
    rho_init = np.ones(nr)
    y0 = np.concatenate([rho_init, u_init])

    # Drive
    w = 2.0 * np.pi * freq
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

        # Derivatives
        drho = np.gradient(rho, dr, edge_order=2)
        du   = np.gradient(u, dr, edge_order=2)
        drho[0] = (rho[1] - rho[0])/dr
        du[0]   = (u[1] - u[0])/dr

        cs_sq = n_poly * K_eos * (np.abs(rho)**(n_poly - 1)) # abs protection

        dt_rho = - (u * drho + rho * du + 2*rho*u/r)
        dt_u   = - (u * du + (cs_sq / (rho+1e-6)) * drho) # div/0 protection

        # Viscosity
        visc = 0.05 * dr # Higher viscosity for stability in coarse grid
        dt_rho[1:-1] += visc * np.diff(rho, 2)/dr**2
        dt_u[1:-1]   += visc * np.diff(u, 2)/dr**2

        return np.concatenate([dt_rho, dt_u])

    # Time integration
    T = 1.0 / freq
    t_span = (0, 3.0*T) # 3 cycles
    # Eval last cycle
    t_eval = np.linspace(2.0*T, 3.0*T, 50)

    try:
        sol = solve_ivp(pde_system, t_span, y0, t_eval=t_eval, method='RK23', max_step=T/20)

        # Stats
        rhos_at_throat = sol.y[0, :]
        mean_rho = np.mean(rhos_at_throat)
        min_rho = np.min(rhos_at_throat)

        return mean_rho, min_rho

    except:
        return 1.0, 0.0 # Fail case

if __name__ == "__main__":
    inertial_damping_phase_map()
