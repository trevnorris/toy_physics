import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.collections import LineCollection

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
# 1. PHYSICAL CONSTANTS
# ==============================================================================
RHO0 = 1.0
A0   = 1.0
L0   = 1.85 * A0
K    = 1.0
M_BARE = np.pi**2 * A0**3 * L0 * RHO0

# Geometry & Trap
C_STIFF = (3 * K * M_BARE**5) / (25 * L0**4 * np.pi**8)
C_QUANT = C_STIFF * (A0**10) * 0.1
F_PUSH_EQ = (C_STIFF / A0**13) + (C_QUANT / A0**3)
K_TRAP    = F_PUSH_EQ / A0

# Forces
V0_INT  = 5.0 * K * RHO0**5
G_EFF   = 0.05 * V0_INT
ALPHA_B = 0.8 * V0_INT
ALPHA_C = 0.6 * V0_INT

# ==============================================================================
# 2. THE PLASMA CLASS
# ==============================================================================

class ThroatPlasma:
    def __init__(self, n_particles=50, box_size=40.0):
        self.N = n_particles
        self.L = box_size

        # --- SAFE SPAWN (Prevents relativistic explosions) ---
        # We try to place particles at least distance '2.0' apart.
        pos_cpu = []
        for i in range(self.N):
            while True:
                candidate = np.random.uniform(-self.L/2.5, self.L/2.5, size=2)
                if not pos_cpu: # First particle
                    pos_cpu.append(candidate)
                    break
                # Check distance to all existing
                dists = np.linalg.norm(np.array(pos_cpu) - candidate, axis=1)
                if np.all(dists > 2.0):
                    pos_cpu.append(candidate)
                    break
        self.pos = xp.array(pos_cpu)

        # Initial Velocities
        self.vel = xp.random.normal(0, 0.3, size=(self.N, 2))

        # Geometry
        self.a = xp.ones(self.N) * A0
        self.va = xp.zeros(self.N)
        self.m_eff = xp.ones(self.N) * M_BARE * 1.5

        # Spin
        self.spin = xp.random.choice(xp.array([-1.0, 0.0, 1.0]), size=self.N)

        # --- TRAIL HISTORY ---
        # We store the last 15 positions for the visual trail
        self.trail_len = 2000
        # Shape: (History, N, 2)
        self.history = xp.tile(self.pos[None, :, :], (self.trail_len, 1, 1))

    def compute_forces(self):
        # 1. Vectors
        r_vec = self.pos[None, :, :] - self.pos[:, None, :]
        dist_sq = xp.sum(r_vec**2, axis=2)
        xp.fill_diagonal(dist_sq, 1.0)
        dist = xp.sqrt(dist_sq)
        unit_vec = r_vec / dist[:, :, None]

        # 2. Forces (Scalar + Electric)
        f_grav = G_EFF / dist_sq

        a_mean = 0.5 * (self.a[:, None] + self.a[None, :])
        f_nuc = (8.0 * dist * V0_INT) / (5.0 * a_mean**2) * xp.exp(-(4.0/5.0) * dist_sq / (a_mean**2))

        q_prod = self.spin[:, None] * self.spin[None, :]
        f_coul = ALPHA_C * q_prod / dist_sq

        f_central = f_nuc + f_coul - f_grav
        xp.fill_diagonal(f_central, 0.0)
        F_central = -xp.sum(f_central[:, :, None] * unit_vec, axis=1)

        # 3. Magnetic (Lorentz)
        v_j = self.vel[None, :, :]
        r_source = -r_vec
        cross_z = v_j[:, :, 0] * r_source[:, :, 1] - v_j[:, :, 1] * r_source[:, :, 0]
        denom = (dist**2 + 0.5)**1.5
        b_contrib = ALPHA_B * self.spin[None, :] * cross_z / denom
        xp.fill_diagonal(b_contrib, 0.0)
        B_z = xp.sum(b_contrib, axis=1)

        F_lor_x =  self.spin * self.vel[:, 1] * B_z
        F_lor_y = -self.spin * self.vel[:, 0] * B_z
        F_vector = xp.stack([F_lor_x, F_lor_y], axis=1)

        # 4. Geometry
        squeeze = -0.5 * xp.sum(f_nuc, axis=1)
        f_int = (C_STIFF / self.a**13) + (C_QUANT / self.a**3) - (K_TRAP * self.a)
        F_geom = f_int + squeeze

        return F_central + F_vector, F_geom

    def step(self, dt=0.005):
        m_geom = M_BARE * 0.5
        F_pos, F_geom = self.compute_forces()

        # Geometry Update
        self.va += (F_geom / m_geom) * dt
        self.a  += self.va * dt
        self.m_eff = M_BARE * (self.a / A0)**3 * 1.5

        # Position Update
        acc_pos = F_pos / self.m_eff[:, None]
        self.vel += acc_pos * dt
        self.pos += self.vel * dt

        # Wall Bounce
        mask_hi = self.pos > self.L/2
        mask_lo = self.pos < -self.L/2
        self.vel[mask_hi] *= -0.9
        self.pos[mask_hi] = self.L/2
        self.vel[mask_lo] *= -0.9
        self.pos[mask_lo] = -self.L/2

        # Update History (Roll the array)
        # Shift everything back by 1
        self.history[:-1] = self.history[1:]
        # Set last element to current pos
        self.history[-1] = self.pos

    def get_data_cpu(self):
        if GPU_AVAILABLE:
            return (
                xp.asnumpy(self.pos),
                xp.asnumpy(self.a),
                xp.asnumpy(self.spin),
                xp.asnumpy(self.history)
            )
        return (self.pos, self.a, self.spin, self.history)

# ==============================================================================
# 3. VISUALIZATION
# ==============================================================================

sim = ThroatPlasma(n_particles=50, box_size=20.0)

fig, ax = plt.subplots(figsize=(9, 9))
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_facecolor('black')
ax.set_title("4D Plasma: Trails & Tunneling")

# 1. Trails (LineCollection)
# We need one collection per particle to color them differently,
# or one giant collection with segments colored by spin.
# Let's do one giant collection for speed.
lines = LineCollection([], linewidths=1.0, alpha=0.6)
ax.add_collection(lines)

# 2. Particles (Scatter)
cmap = plt.cm.coolwarm
scatter = ax.scatter([], [], c=[], cmap=cmap, s=[], edgecolors='white', alpha=1.0, vmin=-1, vmax=1)

def init():
    scatter.set_offsets(np.empty((0, 2)))
    lines.set_segments([])
    return scatter, lines

def update(frame):
    for _ in range(4):
        sim.step(dt=0.005)

    pos, a, spin, history = sim.get_data_cpu()

    # Update Particles
    scatter.set_offsets(pos)
    scatter.set_array(spin)
    scatter.set_sizes((a / A0) * 100)

    # Update Trails
    # history shape: (15, N, 2)
    # We need segments: (N, 15, 2) -> interpreted as lines
    # Transpose to (N, History, XY)
    trajectories = history.transpose(1, 0, 2)
    lines.set_segments(trajectories)

    # Color lines by spin (Red/Blue/White)
    # Map spin (-1, 0, 1) to colors manually or use cmap
    colors = cmap((spin + 1) / 2) # Norm -1..1 to 0..1
    lines.set_color(colors)

    return scatter, lines

anim = FuncAnimation(fig, update, init_func=init, frames=300, interval=20, blit=True)
plt.show()

print(">> System Status: Trails Active.")
print(">> Long Trails = High Speed.")
print(">> Observe the 'Delay' after a fast particle tunnels through a slow one.")
