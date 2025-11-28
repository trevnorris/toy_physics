import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time

# Path Setup
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

from superfluid_lib.core import UnitScaler, Grid3D
from superfluid_lib.solvers import PoissonSolverFFT, WaveEquationSolver
from superfluid_lib.dynamics import DyonEnsemble, ForceAccumulator, SymplecticStepper

try:
    import cupy as cp
    to_cpu = cp.asnumpy
except ImportError:
    import numpy as cp
    to_cpu = lambda x: x


class InterpolatedForceAccumulator(ForceAccumulator):
    def _sample_field(self, field_grid, positions):
        # 1. Normalize coordinates to Grid Index Space
        u = (positions + self.grid.L / 2.0) / self.grid.dx

        # 2. Get lower-left corner indices
        i = cp.floor(u).astype(int)
        i = cp.clip(i, 0, self.grid.N - 2)  # Safety clip

        # 3. Fractional part for weighting (0.0 to 1.0)
        f = u - i

        # 4. Trilinear Interpolation
        fx, fy, fz = f[:, 0], f[:, 1], f[:, 2]
        ix, iy, iz = i[:, 0], i[:, 1], i[:, 2]

        # Gather 8 Corners
        c000 = field_grid[ix, iy, iz]
        c100 = field_grid[ix + 1, iy, iz]
        c010 = field_grid[ix, iy + 1, iz]
        c001 = field_grid[ix, iy, iz + 1]
        c110 = field_grid[ix + 1, iy + 1, iz]
        c101 = field_grid[ix + 1, iy, iz + 1]
        c011 = field_grid[ix, iy + 1, iz + 1]
        c111 = field_grid[ix + 1, iy + 1, iz + 1]

        # Blend X
        c00 = c000 * (1 - fx) + c100 * fx
        c01 = c001 * (1 - fx) + c101 * fx
        c10 = c010 * (1 - fx) + c110 * fx
        c11 = c011 * (1 - fx) + c111 * fx

        # Blend Y
        c0 = c00 * (1 - fy) + c10 * fy
        c1 = c01 * (1 - fy) + c11 * fy

        # Blend Z
        return c0 * (1 - fz) + c1 * fz


def compute_sound_speed(mu, a_sim, e, epsilon):
    return np.sqrt(mu / (epsilon * a_sim * (1 - e**2)))


def find_perihelia(radii):
    perihelion_indices = []
    for k in range(1, len(radii) - 1):
        if radii[k] < radii[k - 1] and radii[k] < radii[k + 1]:
            perihelion_indices.append(k)
    return perihelion_indices


def analyze_orbit(radii, times, xs, ys, scaler, target_period_days):
    radii = np.array(radii)
    times = np.array(times)

    perihelion_indices = find_perihelia(radii)
    if len(perihelion_indices) < 2:
        return {
            "period_days": None,
            "period_error": None,
            "precession_deg": None,
            "perihelia": perihelion_indices,
        }

    p_times = times[perihelion_indices]
    p_diffs = np.diff(p_times)
    avg_period_sim = np.mean(p_diffs)
    avg_period_days = scaler.to_phys_time(avg_period_sim) / 86400.0
    err_period = abs(avg_period_days - target_period_days)

    angles = []
    for idx in perihelion_indices:
        px = xs[idx]
        py = ys[idx]
        angles.append(np.arctan2(py, px))

    angles = np.unwrap(np.array(angles))
    shifts = np.diff(angles)
    precession_deg = np.degrees(np.mean(shifts))

    return {
        "period_days": avg_period_days,
        "period_error": err_period,
        "precession_deg": precession_deg,
        "perihelia": perihelion_indices,
    }


def run_mode(
    mode_name,
    use_wave,
    enable_beta,
    beta,
    c_s_mode,
    grid,
    scaler,
    r_peri_sim,
    v_peri_sim,
    T_total_sim,
    gamma,
    sample_stride,
    c_s_scalar_ref,
    target_period_days,
):
    poisson = PoissonSolverFFT(grid)
    wave = WaveEquationSolver(grid, c_s=c_s_mode, gamma=gamma)
    engine = InterpolatedForceAccumulator(
        grid,
        poisson,
        wave,
        beta=beta,
        enable_beta_inertia=enable_beta,
    )
    stepper = SymplecticStepper()

    system = DyonEnsemble(2)
    SUN, MERC = 0, 1
    system.add_particle(SUN, pos=[0, 0, 0], vel=[0, 0, 0], mass=1.0)
    system.add_particle(MERC, pos=[r_peri_sim, 0, 0], vel=[0, v_peri_sim, 0], mass=1e-6)

    dt_cfl = grid.dx / ((c_s_mode if use_wave else c_s_scalar_ref) * 1.732)
    dt = dt_cfl * 0.2  # tighter for accuracy
    steps = int(T_total_sim / dt)

    print(f"\n[Run:{mode_name}] steps={steps} dt={dt:.6e} c_s={c_s_mode} beta={beta} use_wave={use_wave}")

    rho = engine.deposit_mass(system)
    if use_wave:
        phi_curr = poisson.solve(rho, scaler.G_sim)
        phi_prev = phi_curr.copy()
    else:
        phi_curr = None
        phi_prev = None

    orbit_trace_x = []
    orbit_trace_y = []
    radii = []
    times = []

    start_time = time.time()
    progress_interval = max(steps // 20, 1)

    for i in range(steps):
        rho = engine.deposit_mass(system)
        phi_wave = None
        if use_wave:
            phi_next = wave.step(phi_curr, phi_prev, rho, dt, scaler.G_sim)
            phi_prev, phi_curr = phi_curr, phi_next
            phi_wave = phi_curr

        stepper.clear_acc(system)
        engine.compute_gravity_field(
            system,
            rho,
            scaler.G_sim,
            phi_wave=phi_wave,
            use_wave=use_wave,
        )
        stepper.half_kick(system, dt)
        stepper.drift(system, dt)

        rho = engine.deposit_mass(system)
        stepper.clear_acc(system)
        engine.compute_gravity_field(
            system,
            rho,
            scaler.G_sim,
            phi_wave=phi_wave,
            use_wave=use_wave,
        )
        stepper.half_kick(system, dt)

        if i % sample_stride == 0:
            sun_pos = to_cpu(system.pos[SUN])
            merc_pos = to_cpu(system.pos[MERC])
            rel = merc_pos - sun_pos
            orbit_trace_x.append(rel[0])
            orbit_trace_y.append(rel[1])
            radii.append(np.linalg.norm(rel))
            times.append(i * dt)

        if i % progress_interval == 0:
            pct = 100.0 * i / steps
            elapsed = time.time() - start_time
            print(f"[Run:{mode_name}] {pct:5.1f}% complete ({i}/{steps}), elapsed {elapsed:.1f}s")

    wall = time.time() - start_time
    print(f"[Run:{mode_name}] wall {wall:.2f}s sampled {len(radii)} points")

    results = analyze_orbit(
        radii,
        times,
        orbit_trace_x,
        orbit_trace_y,
        scaler,
        target_period_days=target_period_days,
    )

    plot_name = f"mercury_calibration_result_{mode_name}.png"
    plt.figure(figsize=(8, 8))
    plt.plot(orbit_trace_x, orbit_trace_y, lw=0.5, label="Mercury (relative)")
    px = [orbit_trace_x[i] for i in results["perihelia"]]
    py = [orbit_trace_y[i] for i in results["perihelia"]]
    plt.scatter(px, py, c="r", s=15, zorder=5, label="Perihelia")
    plt.scatter([0], [0], c="y", s=80, label="Sun (relative frame)")
    limit_sim = scaler.to_sim_len(0.6 * (1.496e11))
    plt.xlim(-limit_sim, limit_sim)
    plt.ylim(-limit_sim, limit_sim)

    period_str = "n/a"
    prec_str = "n/a"
    if results["period_days"] is not None:
        period_str = f"{results['period_days']:.2f}d"
    if results["precession_deg"] is not None:
        prec_str = f"{results['precession_deg']:.4f} deg/orbit"

    plt.title(
        f"{mode_name} | period {period_str} vs {target_period_days:.2f}d | precession {prec_str}"
    )
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(plot_name)
    print(f"[Run:{mode_name}] saved {plot_name}")

    return {
        "mode": mode_name,
        "period_days": results["period_days"],
        "period_error": results["period_error"],
        "precession_deg": results["precession_deg"],
        "plot": plot_name,
        "steps": steps,
        "dt": dt,
    }


def run_experiment():
    print("==========================================================")
    print("   EXPERIMENT 2.1: MERCURY CALIBRATION & PRECESSION")
    print("==========================================================")

    AU = 1.496e11           # meters
    M_SUN_SI = 1.989e30     # kg
    YEAR_SI = 31557600.0    # seconds

    A_MERCURY = 0.387098 * AU
    E_MERCURY = 0.205630
    T_MERCURY_DAYS = 87.969

    L_REAL = 6.0 * AU
    N_GRID = 192

    scaler = UnitScaler(L0=L_REAL, T0=YEAR_SI, M0=M_SUN_SI)

    r_peri_si = A_MERCURY * (1 - E_MERCURY)
    v_peri_si = np.sqrt(scaler.G_PHYS_SI * M_SUN_SI * (2 / r_peri_si - 1 / A_MERCURY))

    r_peri_sim = scaler.to_sim_len(r_peri_si)
    v_peri_sim = scaler.to_sim_vel(v_peri_si)

    mu = scaler.G_sim * 1.0
    a_sim = scaler.to_sim_len(A_MERCURY)
    # EPSILON_TARGET = 1e-3  # small parameter mu / (c_s^2 a (1-e^2))
    EPSILON_TARGET = 0.02
    C_S_SCALAR = compute_sound_speed(mu, a_sim, E_MERCURY, EPSILON_TARGET)
    C_S_FULL = np.sqrt(6.0) * C_S_SCALAR
    BETA_FULL = 2.5
    GAMMA = 0.0

    grid = Grid3D(N_GRID, 1.0)

    print(f"\n[Setup] Perihelion {r_peri_si/AU:.3f} AU -> {r_peri_sim:.4f} sim")
    print(f"[Setup] v_peri {v_peri_si/1000:.2f} km/s -> {v_peri_sim:.4f} sim")
    print(f"[Setup] a_sim={a_sim:.5f} epsilon={EPSILON_TARGET} c_s_scalar={C_S_SCALAR:.2f} c_s_full={C_S_FULL:.2f}")

    T_total_si = 1.0 * YEAR_SI
    T_total_sim = scaler.to_sim_time(T_total_si)
    sample_stride = 50

    modes = [
        {
            "name": "newtonian",
            "use_wave": False,
            "enable_beta": False,
            "beta": 0.0,
            "c_s": C_S_SCALAR,
        },
        {
            "name": "scalar_only",
            "use_wave": True,
            "enable_beta": False,
            "beta": 0.0,
            "c_s": C_S_SCALAR,
        },
        {
            "name": "full_1pn_toy",
            "use_wave": True,
            "enable_beta": True,
            "beta": BETA_FULL,
            "c_s": C_S_FULL,
        },
    ]

    summaries = []
    for mode in modes:
        result = run_mode(
            mode_name=mode["name"],
            use_wave=mode["use_wave"],
            enable_beta=mode["enable_beta"],
            beta=mode["beta"],
            c_s_mode=mode["c_s"],
            grid=grid,
            scaler=scaler,
            r_peri_sim=r_peri_sim,
            v_peri_sim=v_peri_sim,
            T_total_sim=T_total_sim,
            gamma=GAMMA,
            sample_stride=sample_stride,
            c_s_scalar_ref=C_S_SCALAR,
            target_period_days=T_MERCURY_DAYS,
        )
        summaries.append(result)

    print("\n[Summary]")
    for s in summaries:
        period = "n/a" if s["period_days"] is None else f"{s['period_days']:.3f}d (err {s['period_error']:.3f}d)"
        prec = "n/a" if s["precession_deg"] is None else f"{s['precession_deg']:.5f} deg/orbit"
        print(f"  {s['mode']:14s} | period {period} | precession {prec} | plot {s['plot']}")


if __name__ == "__main__":
    run_experiment()
