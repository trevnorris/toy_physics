import numpy as np

# --- CONSTANTS ---
RHO0 = 1.0
A0   = 1.0
L0   = 1.85 * A0
K    = 1.0
M_BARE = np.pi**2 * A0**3 * L0 * RHO0

# Speed of Sound in the Stiff Vacuum (P = K rho^5)
# cs = sqrt(dP/drho) = sqrt(5 * K * rho^4)
# At rho=1, K=1, cs = sqrt(5) approx 2.236
C_SOUND = np.sqrt(5.0 * K * RHO0**4)

# Force Calibrations
C_STIFF = (3 * K * M_BARE**5) / (25 * L0**4 * np.pi**8)
C_QUANT = C_STIFF * (A0**10) * 0.1
F_PUSH_EQ = (C_STIFF / A0**13) + (C_QUANT / A0**3)
K_TRAP    = F_PUSH_EQ / A0

# Simulation Controls
# We set drag/drive to push the standard particle to v ~ 1.0 (approx 50% cs)
# Standard Terminal Velocity = DRIVE / DRAG = 5.0 / 4.0 = 1.25
DRAG_COEFF = 4.0
DRIVE_FORCE = 5.0
DT = 0.002
STEPS = 5000
PRINT_INTERVAL = 500

class RelativisticRunaway:
    def __init__(self):
        # Standard
        self.v_std = 0.0
        self.x_std = 0.0
        self.m_std = M_BARE * 1.5

        # 4D Relativistic
        self.v_4d = 0.0
        self.x_4d = 0.0
        self.a_4d = A0
        self.va_4d = 0.0
        self.m_4d = M_BARE * 1.5

    def step(self):
        # 1. STANDARD PHYSICS (Baseline)
        f_net_std  = DRIVE_FORCE - (DRAG_COEFF * self.v_std)
        self.v_std += (f_net_std / self.m_std) * DT
        self.x_std += self.v_std * DT

        # 2. 4D RELATIVISTIC PHYSICS
        # Calculate Lorentz Factor for Pressure
        # Note: If v > cs, the pressure physics breaks (supersonic collapse), we clamp it.
        speed_ratio = abs(self.v_4d) / C_SOUND
        if speed_ratio >= 1.0:
            pressure_scale = 0.0 # Total collapse
        else:
            pressure_scale = np.sqrt(1.0 - speed_ratio**2)

        # Internal Forces (Scaled by Lorentz factor)
        # The standing wave holding the throat open gets weaker as we speed up.
        f_stiff = (C_STIFF / self.a_4d**13) * pressure_scale
        f_quant = (C_QUANT / self.a_4d**3)  * pressure_scale

        # Trap Force (External Harmonic Trap - usually static, so unscaled)
        f_trap  = -(K_TRAP * self.a_4d)

        # Net Geometry Force
        f_geom_net = f_stiff + f_quant + f_trap

        # Update Geometry
        m_geom = M_BARE * 0.5
        self.va_4d += (f_geom_net / m_geom) * DT
        self.a_4d  += self.va_4d * DT

        # Update Inertia (M ~ a^3)
        # Safety clamp to prevent negative mass if a < 0
        safe_a = max(self.a_4d, 0.01)
        self.m_4d = M_BARE * (safe_a / A0)**3 * 1.5

        # Update Motion
        f_net_4d = DRIVE_FORCE - (DRAG_COEFF * self.v_4d)
        self.v_4d += (f_net_4d / self.m_4d) * DT
        self.x_4d += self.v_4d * DT

    def run(self):
        print(f"{'STEP':<8} | {'V_STD':<8} | {'V_4D':<8} | {'RADIUS(a)':<10} | {'MASS_4D':<10} | {'SCALE_FAC':<10}")
        print("-" * 80)

        for step in range(STEPS):
            self.step()

            if step % PRINT_INTERVAL == 0:
                scale = np.sqrt(1 - (min(abs(self.v_4d), C_SOUND)/C_SOUND)**2)
                print(f"{step:<8} | {self.v_std:<8.4f} | {self.v_4d:<8.4f} | {self.a_4d:<10.4f} | {self.m_4d:<10.2f} | {scale:<10.4f}")

        print("-" * 80)
        print(f"Sound Speed (c_s): {C_SOUND:.4f}")
        print(f"Standard Limit:    {DRIVE_FORCE/DRAG_COEFF:.4f}")
        print(f"Final 4D Speed:    {self.v_4d:.4f}")

        if self.v_4d > self.v_std * 1.2:
            print("\n>> SUCCESS: RUNAWAY CONFIRMED.")
            print(">> The 'Photon Pressure' drop caused the throat to shrink, shedding mass.")
            print(">> The particle broke the standard speed limit.")
        else:
            print("\n>> RESULT: Still stable. Try pushing DRIVE_FORCE closer to c_s.")

if __name__ == "__main__":
    sim = RelativisticRunaway()
    sim.run()
