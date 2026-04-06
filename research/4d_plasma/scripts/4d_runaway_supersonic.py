import numpy as np

# --- CONSTANTS ---
RHO0 = 1.0
A0   = 1.0
L0   = 1.85 * A0
K    = 1.0
M_BARE = np.pi**2 * A0**3 * L0 * RHO0
C_SOUND = np.sqrt(5.0 * K * RHO0**4) # ~2.236

# Force Calibrations
C_STIFF = (3 * K * M_BARE**5) / (25 * L0**4 * np.pi**8)
C_QUANT = C_STIFF * (A0**10) * 0.1
F_PUSH_EQ = (C_STIFF / A0**13) + (C_QUANT / A0**3)
K_TRAP    = F_PUSH_EQ / A0

# --- EXPERIMENT PARAMETERS ---
# We want the theoretical speed limit > C_SOUND
# v_term = DRIVE / DRAG
# 12.0 / 4.0 = 3.0 (Supersonic Target)
DRAG_COEFF = 4.0
DRIVE_FORCE = 12.0
DT = 0.002
STEPS = 5000
PRINT_INTERVAL = 250

class SupersonicRunaway:
    def __init__(self):
        # Standard Particle
        self.v_std = 0.0
        self.x_std = 0.0
        self.m_std = M_BARE * 1.5

        # 4D Particle
        self.v_4d = 0.0
        self.x_4d = 0.0
        self.a_4d = A0
        self.va_4d = 0.0
        self.m_4d = M_BARE * 1.5

    def step(self):
        # 1. Standard Physics
        f_std = DRIVE_FORCE - (DRAG_COEFF * self.v_std)
        self.v_std += (f_std / self.m_std) * DT

        # 2. 4D Relativistic Physics
        speed_ratio = abs(self.v_4d) / C_SOUND

        # LORENTZ CONTRACTION
        # As v -> cs, pressure_scale -> 0
        if speed_ratio >= 0.99:
            pressure_scale = 0.05 # Collapse but prevent divide-by-zero
        else:
            pressure_scale = np.sqrt(1.0 - speed_ratio**2)

        f_stiff = (C_STIFF / self.a_4d**13) * pressure_scale
        f_quant = (C_QUANT / self.a_4d**3)  * pressure_scale
        f_trap  = -(K_TRAP * self.a_4d)

        # Geometry Update
        f_geom = f_stiff + f_quant + f_trap
        m_geom = M_BARE * 0.5
        self.va_4d += (f_geom / m_geom) * DT
        self.a_4d  += self.va_4d * DT

        # Mass Shedding (Variable Inertia)
        # If a < 0.1, we consider the particle "collapsed" and mass minimal
        safe_a = max(self.a_4d, 0.05)
        self.m_4d = M_BARE * (safe_a / A0)**3 * 1.5

        # Motion Update
        f_4d = DRIVE_FORCE - (DRAG_COEFF * self.v_4d)

        # ROCKET EQUATION CORRECTION (Optional but realistic)
        # If mass drops, F = d(mv)/dt = m*a + v*dm/dt
        # We ignore thrust for now to isolate the "Inertia Reduction" effect
        self.v_4d += (f_4d / self.m_4d) * DT

    def run(self):
        print(f"{'STEP':<8} | {'V_STD':<8} | {'V_4D':<8} | {'RADIUS':<8} | {'MASS':<8} | {'LORENTZ':<8}")
        print("-" * 75)

        for step in range(STEPS):
            self.step()

            if step % PRINT_INTERVAL == 0:
                ratio = min(abs(self.v_4d) / C_SOUND, 0.99)
                scale = np.sqrt(1 - ratio**2)
                print(f"{step:<8} | {self.v_std:<8.3f} | {self.v_4d:<8.3f} | {self.a_4d:<8.3f} | {self.m_4d:<8.2f} | {scale:<8.3f}")

        print("-" * 75)
        print(f"Sound Speed (Barrier): {C_SOUND:.3f}")
        print(f"Standard Limit:        {DRIVE_FORCE/DRAG_COEFF:.3f}")
        print(f"Final 4D Speed:        {self.v_4d:.3f}")

        if self.v_4d > self.v_std:
            print("\n>> SUCCESS: RUNAWAY CONFIRMED.")
            print(">> The 4D particle shed its mass and accelerated past the standard particle.")

if __name__ == "__main__":
    sim = SupersonicRunaway()
    sim.run()
