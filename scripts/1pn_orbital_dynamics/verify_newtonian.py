import time
import argparse
import numpy as np
import matplotlib.pyplot as plt

# --- 1. Backend Selection (CPU vs GPU) ---
try:
    import cupy as cp
    HAS_GPU = True
    print(" [System] CuPy detected. Running in GPU mode.")
except ImportError:
    import numpy as cp
    HAS_GPU = False
    print(" [System] CuPy not found. Running in CPU mode.")

# --- 2. Core Grid & Physics Classes ---

class SimulationGrid:
    """Manages the 3D spatial domain and FFT k-space."""
    def __init__(self, N, L):
        self.N = N
        self.L = L
        self.dx = L / N

        # Spatial coordinate arrays (centered)
        lin = cp.linspace(-L/2, L/2, N, endpoint=False)
        self.X, self.Y, self.Z = cp.meshgrid(lin, lin, lin, indexing='ij')

        # K-space for FFT
        k_freq = cp.fft.fftfreq(N, d=self.dx) * 2 * cp.pi
        KX, KY, KZ = cp.meshgrid(k_freq, k_freq, k_freq, indexing='ij')
        self.K_sq = KX**2 + KY**2 + KZ**2

        # Inverse Laplacian operator (avoid divide-by-zero at k=0)
        with np.errstate(divide='ignore', invalid='ignore'):
            self.inv_k_sq = -1.0 / self.K_sq
        self.inv_k_sq[0, 0, 0] = 0.0  # Jeans swindle / Zero-mean gauge

class FieldSolver:
    """Handles Poisson (Gravity) and Wave (Scalar Lag) equations."""
    def __init__(self, grid, cs, G=1.0):
        self.grid = grid
        self.cs = cs
        self.G = G

    def solve_poisson(self, rho):
        """Returns Phi_P obeying del^2 Phi = 4*pi*G*rho"""
        # FFT Poisson inherently removes the mean density (Jeans Swindle)
        rho_k = cp.fft.fftn(rho)
        phi_k = (4 * cp.pi * self.G) * rho_k * self.grid.inv_k_sq
        return cp.fft.ifftn(phi_k).real

    def wave_step(self, phi_curr, phi_prev, rho, dt, damping=0.0):
        """Verlet step for wave eq: d^2Phi/dt^2 = cs^2(Laplacian(Phi) - 4piG*(rho - rho_bar))"""
        # Laplacian in k-space is -k^2
        phi_k = cp.fft.fftn(phi_curr)
        lap_phi = cp.fft.ifftn(phi_k * (-self.grid.K_sq)).real

        # --- THE FIX: Explicit Jeans Swindle ---
        # We subtract the mean density so the source matches the
        # periodic Poisson solution.
        rho_background = cp.mean(rho)
        source = 4 * cp.pi * self.G * (rho - rho_background)

        # Acceleration of the field
        # Damping term: -gamma * dPhi/dt approx -gamma * (phi - phi_old)/dt
        vel_term = (phi_curr - phi_prev) / dt
        accel = (self.cs**2) * (lap_phi - source) - damping * vel_term

        phi_next = 2 * phi_curr - phi_prev + accel * (dt**2)
        return phi_next

class ParticleEnsemble:
    """Manages particles. 0 is Central Mass, 1 is Test Particle."""
    def __init__(self, N_particles):
        self.pos = cp.zeros((N_particles, 3), dtype=np.float64)
        self.vel = cp.zeros((N_particles, 3), dtype=np.float64)
        self.mass = cp.ones(N_particles, dtype=np.float64) # Simulation mass

    def deposit_density(self, grid, sigma):
        """Gaussian deposition of mass onto the grid."""
        rho = cp.zeros((grid.N, grid.N, grid.N), dtype=np.float64)

        # Optimization: Only compute Gaussian near particles?
        # For this script, we do full grid (slower but accurate for validation)
        for i in range(len(self.mass)):
            dist_sq = (grid.X - self.pos[i,0])**2 + \
                      (grid.Y - self.pos[i,1])**2 + \
                      (grid.Z - self.pos[i,2])**2

            # Normalize so integral is mass
            norm = self.mass[i] / ((2*np.pi)**1.5 * sigma**3)
            rho += norm * cp.exp(-0.5 * dist_sq / (sigma**2))

        return rho

    def get_accel_from_grid(self, phi, grid):
        """Trilinear interpolation of gradients (-nabla Phi)."""
        # Compute gradients on grid
        grad_x, grad_y, grad_z = cp.gradient(phi, grid.dx)

        # Flatten for mapping
        gx = grad_x.ravel()
        gy = grad_y.ravel()
        gz = grad_z.ravel()

        # Coordinate to index conversion
        grid_coords = (self.pos + (grid.L/2.0)) / grid.dx

        # Map coordinates (using CuPy/SciPy ndimage for simplicity)
        if HAS_GPU:
            from cupyx.scipy.ndimage import map_coordinates
        else:
            from scipy.ndimage import map_coordinates

        # map_coordinates expects (3, N_particles)
        coords = grid_coords.T

        ax = -map_coordinates(grad_x, coords, order=1, mode='wrap')
        ay = -map_coordinates(grad_y, coords, order=1, mode='wrap')
        az = -map_coordinates(grad_z, coords, order=1, mode='wrap')

        return cp.stack((ax, ay, az), axis=1)

# --- 3. Main Simulation Logic ---

def run_simulation(args):
    # Setup
    N = args.res
    L = args.box_size
    grid = SimulationGrid(N, L)
    solver = FieldSolver(grid, cs=args.cs, G=1.0)

    # Particles
    # Index 0: Heavy Central Mass (Static)
    # Index 1: Light Test Particle (Orbiting)
    particles = ParticleEnsemble(2)

    # 1. Central Mass Setup
    particles.pos[0] = cp.array([0, 0, 0])
    particles.vel[0] = cp.array([0, 0, 0])
    particles.mass[0] = 1.0 # mu = G*M = 1.0

    # 2. Test Particle Setup (Circular Orbit at R=1)
    # v_circ = sqrt(GM/r) = sqrt(1/1) = 1.0
    orbit_radius = 2.0
    v_circ = np.sqrt(1.0 / orbit_radius)

    particles.pos[1] = cp.array([orbit_radius, 0, 0])
    particles.vel[1] = cp.array([0, v_circ, 0])
    particles.mass[1] = 0.0 # Test particle (no backreaction on field)

    # Density Width (sigma) - needs to be resolved by grid
    sigma = 1.5 * grid.dx

    print(f"\n--- Initializing Simulation ---")
    print(f"Grid: {N}^3, Box: {L}, dx: {grid.dx:.3f}")
    print(f"Orbit Radius: {orbit_radius} (approx {orbit_radius/grid.dx:.1f} grid cells)")
    print(f"Sound Speed (cs): {args.cs}")

    # --- PHASE 1: Field Relaxation (Static Limit Theorem) ---
    print("\n[Phase 1] Relaxing Scalar Field to Poisson Solution...")

    # Compute Density (Only central mass contributes, test mass is 0)
    rho = particles.deposit_density(grid, sigma)

    # Initial Guess: Pure Poisson (Instantaneous)
    phi_P = solver.solve_poisson(rho)
    phi_curr = phi_P.copy()
    phi_prev = phi_P.copy()

    # Run damped steps
    dt_wave = 0.2 * grid.dx / args.cs # CFL condition
    for _ in range(1000):
        phi_next = solver.wave_step(phi_curr, phi_prev, rho, dt_wave, damping=2.0)
        phi_prev = phi_curr
        phi_curr = phi_next

    # Check Stability
    resid = cp.linalg.norm(phi_curr - phi_P) / cp.linalg.norm(phi_P)
    print(f"Field Relaxation Complete. Residual vs Poisson: {resid:.2e}")

    # --- PHASE 2: Dynamic Orbit ---
    print("\n[Phase 2] Running Orbit with Live Wave Equation...")

    T_orbit = 2 * np.pi * np.sqrt(orbit_radius**3 / 1.0)
    total_time = args.orbits * T_orbit
    dt = T_orbit / args.steps_per_orbit

    # Sub-cycling calculation
    n_wave_substeps = int(np.ceil(dt / dt_wave))
    dt_wave_actual = dt / n_wave_substeps

    print(f"Orbital Period: {T_orbit:.2f}")
    print(f"Particle dt: {dt:.4f}")
    print(f"Wave sub-steps per particle step: {n_wave_substeps}")

    history_r = []
    history_E = []
    history_t = []

    t = 0
    acc = particles.get_accel_from_grid(phi_curr, grid)

    start_time = time.time()

    import tqdm
    pbar = tqdm.tqdm(total=int(args.orbits * args.steps_per_orbit))

    step = 0
    while t < total_time:
        # A. Field Update (Sub-cycling)
        for _ in range(n_wave_substeps):
            phi_next = solver.wave_step(phi_curr, phi_prev, rho, dt_wave_actual, damping=0.0)
            phi_prev = phi_curr
            phi_curr = phi_next

        # B. Particle Update (Kick-Drift-Kick)
        particles.vel[1] += 0.5 * acc[1] * dt
        particles.pos[1] += particles.vel[1] * dt
        acc = particles.get_accel_from_grid(phi_curr, grid)
        particles.vel[1] += 0.5 * acc[1] * dt

        # C. Diagnostics
        if step % 10 == 0:
            pos_cpu = particles.pos[1].get() if HAS_GPU else particles.pos[1]
            vel_cpu = particles.vel[1].get() if HAS_GPU else particles.vel[1]

            r = np.linalg.norm(pos_cpu)
            v_sq = np.linalg.norm(vel_cpu)**2

            # Simple Total Energy = KE + PE_exact (to check orbit mechanics alone)
            # PE_exact = -1/r
            tot_energy = 0.5 * v_sq - (1.0/r)

            history_r.append(r)
            history_E.append(tot_energy)
            history_t.append(t)

        t += dt
        step += 1
        pbar.update(1)

    pbar.close()

    # --- Analysis & Stats Table ---
    history_r = np.array(history_r)
    history_E = np.array(history_E)
    history_t = np.array(history_t)

    # --- SAVE RAW DATA (Crucial for Paper) ---
    np.savez("newtonian_hero_run.npz", t=history_t, r=history_r, E=history_E)
    print("Raw data saved to newtonian_hero_run.npz")

    # Stats
    mean_r = np.mean(history_r)
    max_dev_r = np.max(np.abs(history_r - orbit_radius))

    E0 = history_E[0]
    frac_err_E = np.abs((history_E - E0)/E0)
    max_err_E = np.max(frac_err_E)
    drift_E = (history_E[-1] - E0)/E0

    print("\n" + "="*45)
    print("      SIMULATION DIAGNOSTICS SUMMARY")
    print("="*45)
    print(f"{'Metric':<25} | {'Value':<15}")
    print("-" * 45)
    print(f"{'Target Radius':<25} | {orbit_radius:.6f}")
    print(f"{'Mean Radius':<25} | {mean_r:.6f}")
    print(f"{'Radius Error (%)':<25} | {abs(mean_r-orbit_radius)/orbit_radius*100:.4f} %")
    print(f"{'Max Radius Dev':<25} | {max_dev_r:.6f}")
    print("-" * 45)
    print(f"{'Energy Max Err (%)':<25} | {max_err_E*100:.4f} %")
    print(f"{'Energy Net Drift (%)':<25} | {drift_E*100:.4f} %")
    print("-" * 45)
    print(f"{'Field Relaxation Resid':<25} | {resid:.2e}")
    print(f"{'Grid Resolution (N)':<25} | {args.res}")
    print(f"{'Orbit/Cell Ratio':<25} | {orbit_radius/grid.dx:.2f} px")
    print("="*45 + "\n")

    # Plotting
    plt.figure(figsize=(10, 8))

    plt.subplot(2, 1, 1)
    plt.plot(history_t, history_r, label='Grid Orbit')
    plt.axhline(orbit_radius, color='r', linestyle='--', label='Newtonian Exact')
    plt.ylabel('Orbital Radius')
    plt.title(f'Orbit Stability (Box={L}, N={N}, cs={args.cs})')
    plt.legend()
    plt.grid(True)

    plt.subplot(2, 1, 2)
    plt.plot(history_t, frac_err_E)
    plt.ylabel('Fractional Energy Error')
    plt.xlabel('Time')
    plt.grid(True)

    plt.tight_layout()
    plt.savefig('verify_newtonian_results.png')
    print("Results saved to verify_newtonian_results.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Verify Newtonian Limit of Superfluid Gravity")
    parser.add_argument("--res", type=int, default=64, help="Grid Resolution (per side)")
    parser.add_argument("--box_size", type=float, default=20.0, help="Simulation Box Size")
    parser.add_argument("--cs", type=float, default=100.0, help="Sound Speed")
    parser.add_argument("--orbits", type=float, default=2.0, help="Number of orbits to simulate")
    parser.add_argument("--steps_per_orbit", type=int, default=200, help="Particle timesteps per orbit")

    args = parser.parse_args()

    run_simulation(args)
