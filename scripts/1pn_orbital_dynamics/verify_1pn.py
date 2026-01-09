import time
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# --- 1. Backend Selection ---
try:
    import cupy as cp
    HAS_GPU = True
    print(" [System] CuPy detected. Running in GPU mode.")
except ImportError:
    import numpy as cp
    HAS_GPU = False
    print(" [System] CuPy not found. Running in CPU mode.")

# --- 2. Core Grid & Physics ---

class SimulationGrid:
    def __init__(self, N, L):
        self.N = N
        self.L = L
        self.dx = L / N
        lin = cp.linspace(-L/2, L/2, N, endpoint=False)
        self.X, self.Y, self.Z = cp.meshgrid(lin, lin, lin, indexing='ij')
        k_freq = cp.fft.fftfreq(N, d=self.dx) * 2 * cp.pi
        KX, KY, KZ = cp.meshgrid(k_freq, k_freq, k_freq, indexing='ij')
        self.K_sq = KX**2 + KY**2 + KZ**2
        with np.errstate(divide='ignore', invalid='ignore'):
            self.inv_k_sq = -1.0 / self.K_sq
        self.inv_k_sq[0, 0, 0] = 0.0

class FieldSolver:
    def __init__(self, grid, cs, G=1.0):
        self.grid = grid
        self.cs = cs
        self.G = G

    def solve_poisson(self, rho):
        rho_k = cp.fft.fftn(rho)
        phi_k = (4 * cp.pi * self.G) * rho_k * self.grid.inv_k_sq
        return cp.fft.ifftn(phi_k).real

    def wave_step(self, phi_curr, phi_prev, rho, dt, damping=0.0):
        phi_k = cp.fft.fftn(phi_curr)
        lap_phi = cp.fft.ifftn(phi_k * (-self.grid.K_sq)).real

        # Fixed Source (Jeans Swindle)
        rho_bg = cp.mean(rho)
        source = 4 * cp.pi * self.G * (rho - rho_bg)

        vel_term = (phi_curr - phi_prev) / dt
        accel = (self.cs**2) * (lap_phi - source) - damping * vel_term
        return 2 * phi_curr - phi_prev + accel * (dt**2)

class ParticleEnsemble:
    def __init__(self, N_particles):
        self.pos = cp.zeros((N_particles, 3), dtype=np.float64)
        self.vel = cp.zeros((N_particles, 3), dtype=np.float64)
        self.mass = cp.ones(N_particles, dtype=np.float64)

    def deposit_density(self, grid, sigma):
        rho = cp.zeros((grid.N, grid.N, grid.N), dtype=np.float64)
        for i in range(len(self.mass)):
            dist_sq = (grid.X - self.pos[i,0])**2 + \
                      (grid.Y - self.pos[i,1])**2 + \
                      (grid.Z - self.pos[i,2])**2
            norm = self.mass[i] / ((2*np.pi)**1.5 * sigma**3)
            rho += norm * cp.exp(-0.5 * dist_sq / (sigma**2))
        return rho

    def sample_field_and_grad(self, phi, grid):
        """Returns Phi and Gradient at particle positions."""
        grad_x, grad_y, grad_z = cp.gradient(phi, grid.dx)

        # Manual Trilinear Interpolation Helper
        def interpolate(field):
            coords = (self.pos + (grid.L/2.0)) / grid.dx
            base = cp.floor(coords).astype(cp.int64)
            frac = coords - base

            # Wrap indices
            idx_x = base[:, 0] % grid.N
            idx_y = base[:, 1] % grid.N
            idx_z = base[:, 2] % grid.N

            idx_x1 = (idx_x + 1) % grid.N
            idx_y1 = (idx_y + 1) % grid.N
            idx_z1 = (idx_z + 1) % grid.N

            # 8 corners
            c000 = field[idx_x, idx_y, idx_z]
            c100 = field[idx_x1, idx_y, idx_z]
            c010 = field[idx_x, idx_y1, idx_z]
            c001 = field[idx_x, idx_y, idx_z1]
            c110 = field[idx_x1, idx_y1, idx_z]
            c101 = field[idx_x1, idx_y, idx_z1]
            c011 = field[idx_x, idx_y1, idx_z1]
            c111 = field[idx_x1, idx_y1, idx_z1]

            # Lerp X
            c00 = c000*(1-frac[:,0]) + c100*frac[:,0]
            c01 = c001*(1-frac[:,0]) + c101*frac[:,0]
            c10 = c010*(1-frac[:,0]) + c110*frac[:,0]
            c11 = c011*(1-frac[:,0]) + c111*frac[:,0]

            # Lerp Y
            c0 = c00*(1-frac[:,1]) + c10*frac[:,1]
            c1 = c01*(1-frac[:,1]) + c11*frac[:,1]

            # Lerp Z
            return c0*(1-frac[:,2]) + c1*frac[:,2]

        # Interpolate Phi and Gradient components
        phi_val = interpolate(phi)
        gx = interpolate(grad_x)
        gy = interpolate(grad_y)
        gz = interpolate(grad_z)

        # Gradient is -Force, so Accel = -Gradient
        # But for beta-calc we need the actual gradient of Phi
        grad_vec = cp.stack((gx, gy, gz), axis=1)

        # Base Acceleration (Newtonian part) is -Gradient
        acc_base = -grad_vec

        return phi_val, grad_vec, acc_base

def compute_1pn_accel(acc_base, phi_val, grad_phi, vel, cs, beta):
    """
    Applies the beta-inertia correction:
    sigma = -beta * Phi / cs^2
    a_full = 1/(1+sigma) * [ g - (v.grad_sigma)v + 0.5 v^2 grad_sigma ]
    """
    if beta == 0.0:
        return acc_base

    sigma = -beta * phi_val / (cs**2)
    grad_sigma = -beta * grad_phi / (cs**2)

    # Dot products
    v_sq = cp.sum(vel**2, axis=1)
    v_dot_grad_sigma = cp.sum(vel * grad_sigma, axis=1)

    # Expand dimensions for vector math
    sigma = sigma[:, None]
    v_sq = v_sq[:, None]
    v_dot_grad_sigma = v_dot_grad_sigma[:, None]

    term1 = acc_base
    term2 = -v_dot_grad_sigma * vel
    term3 = 0.5 * v_sq * grad_sigma

    acc_full = (1.0 / (1.0 + sigma)) * (term1 + term2 + term3)
    return acc_full

# --- 3. Simulation Logic ---

def run_1pn_test(args):
    # Setup
    N = args.res
    L = args.box_size
    grid = SimulationGrid(N, L)
    solver = FieldSolver(grid, cs=args.cs, G=1.0)
    particles = ParticleEnsemble(2)

    # 1. Central Mass
    particles.pos[0] = cp.array([0, 0, 0])
    particles.vel[0] = cp.array([0, 0, 0])
    particles.mass[0] = 1.0 # mu = 1.0

    # 2. Test Particle (Eccentric Orbit)
    # Start at Perihelion
    a = 2.0
    e = 0.2
    r_p = a * (1 - e)
    v_p = np.sqrt((1.0/a) * ((1+e)/(1-e)))

    particles.pos[1] = cp.array([r_p, 0, 0])
    particles.vel[1] = cp.array([0, v_p, 0])
    particles.mass[1] = 0.0

    sigma = 1.5 * grid.dx

    print(f"\n--- 1PN Precession Test ---")
    print(f"Beta: {args.beta} (Target: 3.0)")
    print(f"Sound Speed: {args.cs}")
    print(f"Eccentricity: {e}")

    # Theoretical Prediction
    # D_phi = 6 * pi * mu / (cs^2 * a * (1-e^2)) * (beta/3)
    # Note: The formula usually assumes beta=3 gives the GR value (6pi)
    # If beta=0, precession should be 0.
    precess_theory = (2 * args.beta * np.pi * 1.0) / (args.cs**2 * a * (1 - e**2))
    print(f"Theoretical Precession: {precess_theory:.6f} rad/orbit")

    # Relax Field
    print("Relaxing field...")
    rho = particles.deposit_density(grid, sigma)
    phi_P = solver.solve_poisson(rho)
    phi_curr = phi_P.copy()
    phi_prev = phi_P.copy() # Start static

    # Orbit Loop
    T_orbit = 2 * np.pi * np.sqrt(a**3 / 1.0)
    total_time = args.orbits * T_orbit
    dt = T_orbit / args.steps_per_orbit

    dt_wave = 0.2 * grid.dx / args.cs
    n_wave_substeps = int(np.ceil(dt / dt_wave))
    dt_wave_actual = dt / n_wave_substeps

    history_angle = []
    history_time = []

    # Initial Accel
    phi_val, grad_phi, acc_base = particles.sample_field_and_grad(phi_curr, grid)
    acc = compute_1pn_accel(acc_base, phi_val, grad_phi, particles.vel, args.cs, args.beta)

    t = 0
    import tqdm
    pbar = tqdm.tqdm(total=int(args.orbits * args.steps_per_orbit))

    while t < total_time:
        # Field
        for _ in range(n_wave_substeps):
            phi_next = solver.wave_step(phi_curr, phi_prev, rho, dt_wave_actual)
            phi_prev = phi_curr
            phi_curr = phi_next

        # Particle (Kick 1)
        particles.vel[1] += 0.5 * acc[1] * dt

        # Drift
        particles.pos[1] += particles.vel[1] * dt

        # Calc Accel
        phi_val, grad_phi, acc_base = particles.sample_field_and_grad(phi_curr, grid)
        # Note: We use the half-step velocity for the drag term v.grad_sigma
        acc = compute_1pn_accel(acc_base, phi_val, grad_phi, particles.vel, args.cs, args.beta)

        # Kick 2
        particles.vel[1] += 0.5 * acc[1] * dt

        # Measurement: Runge-Lenz Vector angle
        # e_vec = (v x h) / mu - r_hat
        if t > 0: # Skip t=0 to avoid div zero if initialized perfectly
            pos = particles.pos[1]
            vel = particles.vel[1]
            r_vec = pos
            v_vec = vel
            mu = 1.0

            h_vec = cp.cross(r_vec, v_vec)
            r_norm = cp.linalg.norm(r_vec)

            e_vec = cp.cross(v_vec, h_vec) / mu - (r_vec / r_norm)

            # Angle of perihelion in xy plane
            angle = cp.arctan2(e_vec[1], e_vec[0])

            # Unwrap to handle crossing 2pi
            if HAS_GPU:
                angle_cpu = float(angle.get())
                t_cpu = float(t)
            else:
                angle_cpu = float(angle)
                t_cpu = float(t)

            history_angle.append(angle_cpu)
            history_time.append(t_cpu)

        t += dt
        pbar.update(1)
    pbar.close()

    # Analysis
    times = np.array(history_time)
    angles = np.unwrap(np.array(history_angle))

    # --- SAVE RAW DATA ---
    np.savez("1pn_hero_run.npz", t=times, angle=angles, beta=args.beta, cs=args.cs, res=args.res)
    print("Raw data saved to 1pn_hero_run.npz")

    # Linear Regression to get slope (rad/time)
    res = linregress(times, angles)
    slope_rad_per_time = res.slope
    measured_precess_per_orbit = slope_rad_per_time * T_orbit

    print("\n" + "="*45)
    print("      1PN PRECESSION RESULTS")
    print("="*45)
    print(f"{'Sound Speed (cs)':<25} | {args.cs}")
    print(f"{'Beta Factor':<25} | {args.beta}")
    print("-" * 45)
    print(f"{'Theoretical (rad/orbit)':<25} | {precess_theory:.6f}")
    print(f"{'Measured (rad/orbit)':<25} | {measured_precess_per_orbit:.6f}")
    print("-" * 45)
    error_pct = abs(measured_precess_per_orbit - precess_theory) / precess_theory * 100
    print(f"{'Error (%)':<25} | {error_pct:.4f} %")
    print("="*45 + "\n")

    plt.figure()
    plt.plot(times, angles, label='Measured Angle')
    plt.plot(times, res.intercept + res.slope*times, 'r--', label='Fit')
    plt.xlabel('Time')
    plt.ylabel('Perihelion Angle (rad)')
    plt.title(f'1PN Precession (Beta={args.beta}, cs={args.cs})')
    plt.legend()
    plt.savefig('verify_1pn_results.png')
    print("Plot saved to verify_1pn_results.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--res", type=int, default=64)
    parser.add_argument("--box_size", type=float, default=20.0)
    parser.add_argument("--cs", type=float, default=20.0) # Low cs for visible effect
    parser.add_argument("--beta", type=float, default=3.0) # Target value
    parser.add_argument("--orbits", type=float, default=5.0)
    parser.add_argument("--steps_per_orbit", type=int, default=500)
    args = parser.parse_args()

    run_1pn_test(args)
