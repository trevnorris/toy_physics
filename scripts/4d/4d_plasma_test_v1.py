import numpy as np
import matplotlib.pyplot as plt

# --- GPU / CPU SWITCH ---
try:
    import cupy as xp
    print(">> Using CuPy (GPU Acceleration Enabled)")
except ImportError:
    import numpy as xp
    print(">> CuPy not found. Falling back to NumPy (CPU Mode)")

# ==============================================================================
# 1. PHYSICAL CONSTANTS (STIFF WATER MODEL)
# ==============================================================================
RHO0 = 1.0
A0   = 1.0
L0   = 1.85 * A0
K    = 1.0
M_BARE = np.pi**2 * A0**3 * L0 * RHO0

# Force Calibrations
# Force ~ 1/a^13 (Positive/Expansive)
C_STIFF = (3 * K * M_BARE**5) / (25 * L0**4 * np.pi**8)
C_QUANT = C_STIFF * (A0**10) * 0.1 # Reduced quantum vs stiff contribution

# Calculate Equilibrium Force to define the Trap Strength
# At a=A0, Internal Push = Trap Pull
F_PUSH_EQ = (C_STIFF / A0**13) + (C_QUANT / A0**3)
K_TRAP    = F_PUSH_EQ / A0  # Harmonic trap constant

V0_INT  = 5.0 * K * RHO0**5
WIDTH_SCALE = 0.8

# ==============================================================================
# 2. THE SIMULATION ENGINE
# ==============================================================================

def simulate_collision(v_impact=0.2, steps=1000, dt=0.001): # Smaller steps for stability
    print(f"\n--- Simulating Collision (v={v_impact} c_s) ---")

    x = xp.array([-3.0, 3.0])
    v = xp.array([v_impact, -v_impact])
    a = xp.array([A0, A0])
    va = xp.array([0.0, 0.0])

    history = {'t': [], 'x1': [], 'x2': [], 'a1': [], 'a2': [], 'ke': []}
    t = 0.0

    for step in range(steps):
        # 1. Forces
        dx = x[1] - x[0]
        dist = xp.abs(dx)

        # Interaction (Gaussian Repulsion)
        a_mean = 0.5 * (a[0] + a[1])
        exp_term = xp.exp(-(4.0/5.0) * (dist**2) / (a_mean**2))
        f_mag = (8.0 * dist * V0_INT) / (5.0 * a_mean**2) * exp_term
        f_int = xp.array([-f_mag, f_mag]) * xp.sign(dx)

        # Geometry Forces (Corrected Signs)
        # Internal Pressure (1/a^13) is EXPANSIVE (+)
        # Quantum Pressure (1/a^3) is EXPANSIVE (+)
        # External Trap (r) is COMPRESSIVE (-)
        f_internal = (C_STIFF / a**13) + (C_QUANT / a**3)
        f_external = -K_TRAP * a
        f_geom = f_internal + f_external

        # Coupling: Interaction squeezes the throat
        f_coupling = -0.5 * f_mag

        # 2. Updates
        m_geom = M_BARE * 0.5 # Wall inertia
        va += ((f_geom + f_coupling) / m_geom) * dt
        a  += va * dt

        # Variable Inertia
        m_eff = M_BARE * (a/A0)**3 * 1.5
        v += (f_int / m_eff) * dt
        x += v * dt

        if step % 5 == 0:
            history['t'].append(t)
            history['x1'].append(float(x[0]))
            history['x2'].append(float(x[1]))
            history['a1'].append(float(a[0]))
            history['a2'].append(float(a[1]))
            ke = 0.5 * float(xp.sum(m_eff * v**2))
            history['ke'].append(ke)
        t += dt

    return history

# ==============================================================================
# 3. RUN
# ==============================================================================
data = simulate_collision(v_impact=0.99, steps=25000)

fig, ax = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

ax[0].plot(data['t'], data['x1'], 'b')
ax[0].plot(data['t'], data['x2'], 'r')
ax[0].set_ylabel("Position")
ax[0].set_title("Trajectories")
ax[0].grid(True)

ax[1].plot(data['t'], data['a1'], 'b')
ax[1].plot(data['t'], data['a2'], 'r--')
ax[1].axhline(A0, color='k', linestyle=':')
ax[1].set_ylabel("Radius a")
ax[1].set_title("Breathing Mode (Look for Ringing)")
ax[1].grid(True)

ax[2].plot(data['t'], data['ke'], 'g')
ax[2].set_ylabel("Kinetic Energy")
ax[2].grid(True)

plt.tight_layout()
plt.show()
