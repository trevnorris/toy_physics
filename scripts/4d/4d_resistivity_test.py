import numpy as np

# --- CONSTANTS ---
RHO0 = 1.0
A0   = 1.0
L0   = 1.85 * A0
K    = 1.0
M_BARE = np.pi**2 * A0**3 * L0 * RHO0

# Geometry & Forces
C_STIFF = (3 * K * M_BARE**5) / (25 * L0**4 * np.pi**8)
C_QUANT = C_STIFF * (A0**10) * 0.1
F_PUSH_EQ = (C_STIFF / A0**13) + (C_QUANT / A0**3)
K_TRAP    = F_PUSH_EQ / A0

V0_INT  = 5.0 * K * RHO0**5
G_EFF   = 0.05 * V0_INT

# Test Parameters
NUM_PARTICLES = 30
BOX_SIZE      = 20.0
DRIVE_FORCE   = 0.5   # The "Voltage" driving the current
DRAG_COEFF    = 0.1   # Base background friction
DT            = 0.005
STEPS         = 4000
PRINT_RATE    = 500

class ResistivityTest:
    def __init__(self):
        # Initialize positions same for both to be fair
        self.pos = np.random.uniform(-BOX_SIZE/3, BOX_SIZE/3, size=(NUM_PARTICLES, 2))

        # --- UNIVERSE A: STANDARD MODEL ---
        self.vel_std = np.zeros((NUM_PARTICLES, 2))
        self.m_std   = np.ones(NUM_PARTICLES) * M_BARE * 1.5 # Fixed Mass

        # --- UNIVERSE B: 4D MODEL ---
        self.vel_4d  = np.zeros((NUM_PARTICLES, 2))
        self.a_4d    = np.ones(NUM_PARTICLES) * A0
        self.va_4d   = np.zeros(NUM_PARTICLES)
        self.m_4d    = np.ones(NUM_PARTICLES) * M_BARE * 1.5 # Variable Mass

    def step(self):
        # --- 1. STANDARD PHYSICS STEP ---
        # Pairwise forces (Simple Repulsion)
        r_vec = self.pos[None, :, :] - self.pos[:, None, :]
        dist_sq = np.sum(r_vec**2, axis=2)
        np.fill_diagonal(dist_sq, 1.0)
        dist = np.sqrt(dist_sq)
        unit_vec = r_vec / dist[:, :, None]

        # Simple 1/r^2 repulsion for standard gas
        f_mag_std = V0_INT / (dist_sq + 1.0)
        np.fill_diagonal(f_mag_std, 0.0)
        F_coll_std = -np.sum(f_mag_std[:, :, None] * unit_vec, axis=1)

        # Drive + Drag
        F_drive = np.array([DRIVE_FORCE, 0.0]) # Push Right
        F_drag  = -DRAG_COEFF * self.vel_std
        F_net_std = F_coll_std + F_drive + F_drag

        # Update Standard
        self.vel_std += (F_net_std / self.m_std[:, None]) * DT
        # (We don't update pos here to keep interaction geometry roughly comparable for the test,
        # or we update pos but assume density is low enough that position divergence doesn't dominate yet)
        # Actually, let's let them drift.

        # --- 2. 4D PHYSICS STEP ---
        # Forces (Gaussian Overlap + Geometry)
        a_mean = 0.5 * (self.a_4d[:, None] + self.a_4d[None, :])
        # Gaussian Force
        f_nuc = (8.0 * dist * V0_INT) / (5.0 * a_mean**2) * np.exp(-(4.0/5.0) * dist_sq / (a_mean**2))
        np.fill_diagonal(f_nuc, 0.0)
        F_coll_4d = -np.sum(f_nuc[:, :, None] * unit_vec, axis=1)

        # Geometry Update (The Breathing)
        squeeze = -0.5 * np.sum(f_nuc, axis=1)
        f_int = (C_STIFF / self.a_4d**13) + (C_QUANT / self.a_4d**3) - (K_TRAP * self.a_4d)
        acc_a = (f_int + squeeze) / (M_BARE * 0.5)
        self.va_4d += acc_a * DT
        self.a_4d  += self.va_4d * DT

        # Variable Inertia
        self.m_4d = M_BARE * (self.a_4d / A0)**3 * 1.5

        # Drive + Drag
        F_drag_4d = -DRAG_COEFF * self.vel_4d
        F_net_4d  = F_coll_4d + F_drive + F_drag_4d

        # Update 4D
        self.vel_4d += (F_net_4d / self.m_4d[:, None]) * DT

        # Move both (simple updates)
        self.pos += (self.vel_4d + self.vel_std) * 0.5 * DT # Keep them roughly coupled in space for fair density comparison

        # Reset positions if they fly too far (Periodic-ish)
        self.pos = np.where(self.pos > BOX_SIZE, self.pos - 2*BOX_SIZE, self.pos)
        self.pos = np.where(self.pos < -BOX_SIZE, self.pos + 2*BOX_SIZE, self.pos)

    def measure_temperature(self, vel):
        # Temp = Variance of velocity (Random Kinetic Energy)
        # Subtract mean drift velocity first
        drift = np.mean(vel, axis=0)
        random_vel = vel - drift
        temp = np.sum(random_vel**2) / NUM_PARTICLES
        return drift[0], temp

    def run(self):
        print(f"{'STEP':<8} | {'DRIFT_STD':<10} | {'DRIFT_4D':<10} | {'TEMP_STD':<10} | {'TEMP_4D':<10} | {'RATIO(T)':<8}")
        print("-" * 75)

        for i in range(STEPS):
            self.step()

            if i % PRINT_RATE == 0:
                drift_std, temp_std = self.measure_temperature(self.vel_std)
                drift_4d, temp_4d   = self.measure_temperature(self.vel_4d)

                ratio = temp_4d / (temp_std + 1e-9)
                print(f"{i:<8} | {drift_std:<10.3f} | {drift_4d:<10.3f} | {temp_std:<10.3f} | {temp_4d:<10.3f} | {ratio:<8.2f}")

        print("-" * 75)
        print("INTERPRETATION:")
        print("DRIFT (Current):  If 4D is lower, the 'breathing' acts as extra resistance.")
        print("TEMP (Heat):      If 4D is higher, the model predicts 'Anomalous Heating'.")

if __name__ == "__main__":
    sim = ResistivityTest()
    sim.run()
