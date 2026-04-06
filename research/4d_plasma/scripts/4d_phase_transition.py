import numpy as np

# --- CONSTANTS ---
RHO0 = 1.0
A0   = 1.0
L0   = 1.85 * A0
K    = 1.0
M_BARE = np.pi**2 * A0**3 * L0 * RHO0

C_STIFF = (3 * K * M_BARE**5) / (25 * L0**4 * np.pi**8)
C_QUANT = C_STIFF * (A0**10) * 0.1
F_PUSH_EQ = (C_STIFF / A0**13) + (C_QUANT / A0**3)
K_TRAP    = F_PUSH_EQ / A0

V0_INT  = 5.0 * K * RHO0**5
G_EFF   = 0.05 * V0_INT

# Simulation Setup
NUM_PARTICLES = 40
START_BOX     = 40.0 # Low Density
END_BOX       = 10.0 # High Density (Compressed)
STEPS         = 8000
PRINT_RATE    = 800

class PhaseTransitionTest:
    def __init__(self):
        self.box_size = START_BOX
        self.pos_4d   = np.random.uniform(-self.box_size/2.2, self.box_size/2.2, size=(NUM_PARTICLES, 2))
        self.vel_4d   = np.random.normal(0, 0.2, size=(NUM_PARTICLES, 2))

        # 4D Internal State
        self.a_4d     = np.ones(NUM_PARTICLES) * A0
        self.va_4d    = np.zeros(NUM_PARTICLES)
        self.m_4d     = np.ones(NUM_PARTICLES) * M_BARE * 1.5

        # Standard Control (Fixed Mass, No Breathing)
        self.pos_std  = self.pos_4d.copy()
        self.vel_std  = self.vel_4d.copy()
        self.m_std    = np.ones(NUM_PARTICLES) * M_BARE * 1.5

    def step(self, i):
        # 1. COMPRESS THE BOX (Adiabatic Compression)
        progress = i / STEPS
        current_box = START_BOX - (START_BOX - END_BOX) * progress
        compression_rate = (START_BOX - END_BOX) / STEPS

        # Push particles inward if they hit the shrinking wall
        limit = current_box / 2.0

        # --- STANDARD PHYSICS ---
        # Simple Hard Sphere-ish Repulsion
        r_vec = self.pos_std[None, :, :] - self.pos_std[:, None, :]
        dist_sq = np.sum(r_vec**2, axis=2)
        np.fill_diagonal(dist_sq, 1.0)
        dist = np.sqrt(dist_sq)
        unit_vec = r_vec / dist[:, :, None]

        # Standard Force (1/r^6 Lennard-Jones style repulsion for gas)
        f_std = 50.0 / (dist_sq**3 + 0.1)
        np.fill_diagonal(f_std, 0.0)
        F_coll_std = -np.sum(f_std[:, :, None] * unit_vec, axis=1)

        self.vel_std += (F_coll_std / self.m_std[:, None]) * 0.005
        self.pos_std += self.vel_std * 0.005

        # Wall Bounce (Standard)
        # Add slight energy gain from compression (Piston effect)
        for dim in [0, 1]:
            mask_hi = self.pos_std[:, dim] > limit
            mask_lo = self.pos_std[:, dim] < -limit
            self.vel_std[mask_hi, dim] = -abs(self.vel_std[mask_hi, dim]) - compression_rate
            self.pos_std[mask_hi, dim] = limit
            self.vel_std[mask_lo, dim] = abs(self.vel_std[mask_lo, dim]) + compression_rate
            self.pos_std[mask_lo, dim] = -limit

        # --- 4D PHYSICS ---
        # Forces
        r_vec = self.pos_4d[None, :, :] - self.pos_4d[:, None, :]
        dist_sq = np.sum(r_vec**2, axis=2)
        np.fill_diagonal(dist_sq, 1.0)
        dist = np.sqrt(dist_sq)
        unit_vec = r_vec / dist[:, :, None]

        a_mean = 0.5 * (self.a_4d[:, None] + self.a_4d[None, :])
        f_nuc = (8.0 * dist * V0_INT) / (5.0 * a_mean**2) * np.exp(-(4.0/5.0) * dist_sq / (a_mean**2))
        np.fill_diagonal(f_nuc, 0.0)
        F_coll_4d = -np.sum(f_nuc[:, :, None] * unit_vec, axis=1)

        # Geometry
        squeeze = -0.5 * np.sum(f_nuc, axis=1)
        f_int = (C_STIFF / self.a_4d**13) + (C_QUANT / self.a_4d**3) - (K_TRAP * self.a_4d)
        acc_a = (f_int + squeeze) / (M_BARE * 0.5)
        self.va_4d += acc_a * 0.005
        self.a_4d  += self.va_4d * 0.005
        self.m_4d = M_BARE * (self.a_4d / A0)**3 * 1.5

        # Move
        self.vel_4d += (F_coll_4d / self.m_4d[:, None]) * 0.005
        self.pos_4d += self.vel_4d * 0.005

        # Wall Bounce (4D)
        for dim in [0, 1]:
            mask_hi = self.pos_4d[:, dim] > limit
            mask_lo = self.pos_4d[:, dim] < -limit
            # 4D particles absorb some wall impact into vibration (inelastic boundary)
            self.vel_4d[mask_hi, dim] = -abs(self.vel_4d[mask_hi, dim]) * 0.9 - compression_rate
            self.pos_4d[mask_hi, dim] = limit
            self.vel_4d[mask_lo, dim] = abs(self.vel_4d[mask_lo, dim]) * 0.9 + compression_rate
            self.pos_4d[mask_lo, dim] = -limit

        return current_box

    def measure_temp(self, vel):
        # Temperature = Mean Kinetic Energy
        speed_sq = np.sum(vel**2, axis=1)
        return np.mean(speed_sq)

    def run(self):
        print(f"{'STEP':<8} | {'BOX_SIZE':<8} | {'DENSITY':<8} | {'TEMP_STD':<10} | {'TEMP_4D':<10} | {'RATIO':<8}")
        print("-" * 75)

        for i in range(STEPS):
            box = self.step(i)

            if i % PRINT_RATE == 0:
                area = box * box
                density = NUM_PARTICLES / area

                t_std = self.measure_temp(self.vel_std)
                t_4d  = self.measure_temp(self.vel_4d)
                ratio = t_4d / (t_std + 1e-9)

                print(f"{i:<8} | {box:<8.1f} | {density:<8.3f} | {t_std:<10.3f} | {t_4d:<10.3f} | {ratio:<8.2f}")

        print("-" * 75)
        print("INTERPRETATION:")
        print("If RATIO drops as DENSITY rises, the 4D plasma is entering 'Super-Laminar' mode.")
        print("This mimics the H-Mode transition in fusion reactors.")

if __name__ == "__main__":
    sim = PhaseTransitionTest()
    sim.run()
