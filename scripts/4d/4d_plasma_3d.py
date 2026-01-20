import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

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
# 2. THE 3D PLASMA CLASS
# ==============================================================================

class ThroatPlasma3D:
    def __init__(self, n_particles=40, box_size=30.0):
        self.N = n_particles
        self.L = box_size

        # --- SAFE SPAWN IN 3D ---
        pos_cpu = []
        for i in range(self.N):
            while True:
                # Spawn in a cube
                candidate = np.random.uniform(-self.L/2.5, self.L/2.5, size=3)
                if not pos_cpu:
                    pos_cpu.append(candidate)
                    break
                dists = np.linalg.norm(np.array(pos_cpu) - candidate, axis=1)
                if np.all(dists > 2.5): # Ensure good separation
                    pos_cpu.append(candidate)
                    break
        self.pos = xp.array(pos_cpu)

        # 3D Velocities
        self.vel = xp.random.normal(0, 0.3, size=(self.N, 3))

        # Geometry
        self.a = xp.ones(self.N) * A0
        self.va = xp.zeros(self.N)
        self.m_eff = xp.ones(self.N) * M_BARE * 1.5

        # Spin (Charge)
        self.spin = xp.random.choice(xp.array([-1.0, 0.0, 1.0]), size=self.N)

        # Trails (History)
        self.trail_len = 10
        self.history = xp.tile(self.pos[None, :, :], (self.trail_len, 1, 1))

    def compute_forces(self):
        # 1. Vectors (N, N, 3)
        # r_vec points FROM i TO j
        r_vec = self.pos[None, :, :] - self.pos[:, None, :]

        dist_sq = xp.sum(r_vec**2, axis=2)
        xp.fill_diagonal(dist_sq, 1.0)
        dist = xp.sqrt(dist_sq)
        unit_vec = r_vec / dist[:, :, None]

        # --- SCALAR + COULOMB SECTOR ---
        f_grav = G_EFF / dist_sq

        a_mean = 0.5 * (self.a[:, None] + self.a[None, :])
        f_nuc = (8.0 * dist * V0_INT) / (5.0 * a_mean**2) * xp.exp(-(4.0/5.0) * dist_sq / (a_mean**2))

        q_prod = self.spin[:, None] * self.spin[None, :]
        f_coul = ALPHA_C * q_prod / dist_sq

        # Net Central Magnitude (Repulsive > 0)
        f_central = f_nuc + f_coul - f_grav
        xp.fill_diagonal(f_central, 0.0)

        # Sum forces: F_i = Sum_j ( -mag * unit_vec_ij )
        F_central_vec = -xp.sum(f_central[:, :, None] * unit_vec, axis=1)

        # --- MAGNETIC SECTOR (3D CROSS PRODUCTS) ---
        # B-field at i comes from current j: B ~ v_j x r_ji
        # r_source (j->i) = -r_vec (i->j)
        r_source = -r_vec
        v_j = self.vel[None, :, :]

        # 3D Cross Product: v_j x r_source
        # Result shape (N, N, 3)
        cross_prod = xp.cross(v_j, r_source)

        # Softened Biot-Savart to prevent energy explosion
        # Note: '2.0' softening is higher here to fix the heating issue
        denom = (dist_sq + 2.0)**1.5

        b_contrib = ALPHA_B * self.spin[None, :, None] * cross_prod / denom[:, :, None]

        # Zero self-field
        # (CuPy/NumPy indexing trick for 3D diagonal)
        diag_idx = xp.arange(self.N)
        b_contrib[diag_idx, diag_idx, :] = 0.0

        # B_field at each particle (N, 3)
        B_field = xp.sum(b_contrib, axis=1)

        # Lorentz Force: F = q(v x B)
        F_lorentz = self.spin[:, None] * xp.cross(self.vel, B_field)

        # --- GEOMETRY ---
        squeeze = -0.5 * xp.sum(f_nuc, axis=1)
        f_int = (C_STIFF / self.a**13) + (C_QUANT / self.a**3) - (K_TRAP * self.a)
        F_geom = f_int + squeeze

        return F_central_vec + F_lorentz, F_geom

    def step(self, dt=0.005):
        m_geom = M_BARE * 0.5
        F_pos, F_geom = self.compute_forces()

        self.va += (F_geom / m_geom) * dt
        self.a  += self.va * dt
        self.m_eff = M_BARE * (self.a / A0)**3 * 1.5

        acc_pos = F_pos / self.m_eff[:, None]
        self.vel += acc_pos * dt
        self.pos += self.vel * dt

        # Soft Box (Reflective Cube)
        for dim in range(3):
            mask_hi = self.pos[:, dim] > self.L/2
            mask_lo = self.pos[:, dim] < -self.L/2
            self.vel[mask_hi, dim] *= -0.9
            self.pos[mask_hi, dim] = self.L/2
            self.vel[mask_lo, dim] *= -0.9
            self.pos[mask_lo, dim] = -self.L/2

        # Update History
        self.history[:-1] = self.history[1:]
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
# 3. 3D VISUALIZATION
# ==============================================================================

sim = ThroatPlasma3D(n_particles=500)

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_facecolor('black')
# Hide grid/panes for "Space" look
ax.grid(False)
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.set_title("4D Plasma in 3D Space (Spin & Breathing)")

ax.set_xlim(-15, 15)
ax.set_ylim(-15, 15)
ax.set_zlim(-15, 15)

# Particles
cmap = plt.cm.coolwarm
scatter = ax.scatter([], [], [], c=[], cmap=cmap, s=[], alpha=0.9, vmin=-1, vmax=1)

# Trails (Simple lines for 3D is tricky in FuncAnimation, using dots for trails is faster)
# We will just plot the heads for performance
# (Matplotlib 3D animation is slow with many lines)

def update(frame):
    for _ in range(4):
        sim.step(dt=0.005)

    pos, a, spin, hist = sim.get_data_cpu()

    # Update Scatter
    scatter._offsets3d = (pos[:,0], pos[:,1], pos[:,2])
    scatter.set_array(spin)

    # Map 4th Dimension (Radius a) to Size
    sizes = (a / A0) * 80
    scatter.set_sizes(sizes)

    return scatter,

anim = FuncAnimation(fig, update, frames=200, interval=30, blit=False)
plt.show()

print(">> 3D Simulation Running.")
print(">> Visuals: X,Y,Z is position.")
print(">> Dot Size is the 4th dimension (Throat Radius).")
print(">> Red/Blue are charged. White is neutral.")
