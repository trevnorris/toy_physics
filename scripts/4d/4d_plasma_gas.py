import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# --- GPU / CPU SWITCH ---
try:
    import cupy as xp
    GPU_AVAILABLE = True
    print(">> Using CuPy (GPU Acceleration Enabled)")
except ImportError:
    import numpy as xp
    GPU_AVAILABLE = False
    print(">> CuPy not found. Falling back to NumPy (CPU Mode)")

# ==============================================================================
# 1. PHYSICAL CONSTANTS (STIFF WATER MODEL)
# ==============================================================================
# Units: a0 = 1, rho0 = 1, cs = 1
RHO0 = 1.0
A0   = 1.0
L0   = 1.85 * A0
K    = 1.0
M_BARE = np.pi**2 * A0**3 * L0 * RHO0

# Force Calibrations (Corrected Signs)
# Stiff (1/a^13) and Quantum (1/a^3) are EXPANSIVE (+)
C_STIFF = (3 * K * M_BARE**5) / (25 * L0**4 * np.pi**8)
C_QUANT = C_STIFF * (A0**10) * 0.1

# Trap Force (Linear harmonic) needed to balance equilibrium at A0
F_PUSH_EQ = (C_STIFF / A0**13) + (C_QUANT / A0**3)
K_TRAP    = F_PUSH_EQ / A0

# Interaction Constants
# Short Range "Nuclear" (Gaussian Overlap)
V0_INT  = 5.0 * K * RHO0**5
# Long Range "Gravity" (Leakage/Poisson) -> F ~ G_eff / r^2
# We tune G_eff so gravity is weak compared to nuclear
G_EFF   = 0.05 * V0_INT

# ==============================================================================
# 2. THE SIMULATION CLASS
# ==============================================================================

class ThroatGas:
    def __init__(self, n_particles=50, box_size=40.0):
        self.N = n_particles
        self.L = box_size

        # --- Initialization ---
        # Positions: Random in box
        self.pos = xp.random.uniform(-self.L/2, self.L/2, size=(self.N, 2))

        # Velocities: Random thermal distribution
        v_thermal = 0.2 # fraction of cs
        self.vel = xp.random.normal(0, v_thermal, size=(self.N, 2))

        # Geometry: All start at equilibrium
        self.a = xp.ones(self.N) * A0
        self.va = xp.zeros(self.N)

        # Masses: Will be dynamic
        self.m_eff = xp.ones(self.N) * M_BARE * 1.5 # Initial guess

    def compute_forces(self):
        """
        Vectorized N-Body Force Calculation.
        Returns:
            F_pos: (N, 2) force on positions (x, y)
            F_geom: (N,) force on radius (a)
        """
        # 1. Pairwise Distances (Broadcasting)
        # r_vec[i, j, k] = pos[j, k] - pos[i, k]
        # Shape: (N, N, 2)
        r_vec = self.pos[None, :, :] - self.pos[:, None, :]

        # Distances: Shape (N, N)
        # Add epsilon to diagonal to avoid division by zero
        dist_sq = xp.sum(r_vec**2, axis=2)
        xp.fill_diagonal(dist_sq, 1.0) # Dummy value for self-interaction
        dist = xp.sqrt(dist_sq)

        # Unit Vectors: (N, N, 2)
        unit_vec = r_vec / dist[:, :, None]

        # --- A. Long Range Gravity (1/r^2) ---
        # F_grav = G * M_i * M_j / r^2
        # For simplicity, we assume 'mass charge' ~ 1 for gravity scaling here
        # (Refining this to use actual M_eff(t) is an option but unstable for start)
        f_grav_mag = G_EFF / dist_sq

        # --- B. Short Range Nuclear (Gaussian) ---
        # F_nuc ~ d * exp(-d^2/a^2)
        # Use mean radius for interaction scale
        a_mean_matrix = 0.5 * (self.a[:, None] + self.a[None, :])

        exp_term = xp.exp(-(4.0/5.0) * dist_sq / (a_mean_matrix**2))
        f_nuc_mag = (8.0 * dist * V0_INT) / (5.0 * a_mean_matrix**2) * exp_term

        # Total Pairwise Force Magnitude (Repulsive +, Attractive -)
        # Gravity is attractive (-), Nuclear is repulsive (+)
        f_total_mag = f_nuc_mag - f_grav_mag

        # Zero out self-interactions
        xp.fill_diagonal(f_total_mag, 0.0)

        # Sum forces on each particle i
        # F_pos[i] = Sum_j ( f_mag[i,j] * unit_vec[i,j] )
        # Note: We need a negative sign because r_vec was (j - i).
        # Force on i from j is repulsive if r_vec points to j? No.
        # If j is at +x, repulsive force on i should be -x.
        # r_vec points i -> j (+x). So we need -unit_vec for repulsion.
        F_pos = -xp.sum(f_total_mag[:, :, None] * unit_vec, axis=1)

        # --- C. Geometry Coupling Forces ---
        # The collision squeezes the throat: F_couple ~ -0.5 * F_nuc_mag
        # Sum of all squeezes from all neighbors
        squeeze_forces = -0.5 * xp.sum(f_nuc_mag, axis=1)

        # Internal Spring Forces
        # F_stiff (+) F_quant (+) F_trap (-)
        f_internal = (C_STIFF / self.a**13) + (C_QUANT / self.a**3) - (K_TRAP * self.a)

        F_geom = f_internal + squeeze_forces

        return F_pos, F_geom

    def step(self, dt=0.005):
        # 1. Update Geometry (Symplectic Euler)
        # Wall inertia
        m_geom = M_BARE * 0.5

        F_pos, F_geom = self.compute_forces()

        self.va += (F_geom / m_geom) * dt
        self.a  += self.va * dt

        # 2. Update Effective Mass (Variable Inertia!)
        # M_eff = M_bare * (a/A0)^3 * (1 + C_add)
        self.m_eff = M_BARE * (self.a / A0)**3 * 1.5

        # 3. Update Position (Symplectic Euler)
        # F = ma => a = F/m_eff
        acc_pos = F_pos / self.m_eff[:, None]

        self.vel += acc_pos * dt
        self.pos += self.vel * dt

        # 4. Periodic Boundary Conditions (Soft Box)
        # Instead of wrapping, let's put soft walls at box edge to keep gas contained
        # F_wall ~ -(pos - L/2)
        mask_hi = self.pos > self.L/2
        mask_lo = self.pos < -self.L/2

        self.vel[mask_hi] *= -0.9 # Bounce with damping
        self.pos[mask_hi] = self.L/2

        self.vel[mask_lo] *= -0.9
        self.pos[mask_lo] = -self.L/2

    def get_data_cpu(self):
        """Helper to get data back to CPU for plotting"""
        if GPU_AVAILABLE:
            return (xp.asnumpy(self.pos), xp.asnumpy(self.a), xp.asnumpy(self.vel))
        return (self.pos, self.a, self.vel)

# ==============================================================================
# 3. RUN & ANIMATE
# ==============================================================================

# Initialize
gas = ThroatGas(n_particles=60, box_size=30.0)

# Setup Plot
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(-16, 16)
ax.set_ylim(-16, 16)
ax.set_facecolor('black')
ax.set_title("4D Throat Gas (Variable Inertia & Geometry)")

# Scatter Plot:
# X, Y positions
# Size 's' will be mapped to Radius 'a'
scatter = ax.scatter([], [], c=[], cmap='plasma', s=[], edgecolors='white', alpha=0.9)

def init():
    scatter.set_offsets(np.empty((0, 2)))
    return scatter,

def update(frame):
    # Run a few physics steps per frame for smoothness
    for _ in range(4):
        gas.step(dt=0.005)

    pos, a, vel = gas.get_data_cpu()

    # Update Positions
    scatter.set_offsets(pos)

    # Update Visuals based on 4D State
    # Color mapped to Velocity Magnitude (Temperature)
    speed = np.linalg.norm(vel, axis=1)
    scatter.set_array(speed)

    # Size mapped to Throat Radius 'a'
    # Visual magnification for effect
    sizes = (a / A0) * 100
    scatter.set_sizes(sizes)

    return scatter,

# Create Animation
anim = FuncAnimation(fig, update, init_func=init, frames=200, interval=20, blit=True)

# To display in a notebook/lab, you might use:
# from IPython.display import HTML
# HTML(anim.to_jshtml())

plt.show()

# If running as script, this keeps window open
print(">> Simulation Running...")
print(">> Observe: When particles collide, they should 'blink' (shrink) and change color.")
print(">> This 'blink' is the 4D breathing mode absorbing energy.")
