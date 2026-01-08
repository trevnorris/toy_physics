import numpy as np
import matplotlib.pyplot as plt

# --- 1. SETUP PARAMETERS FROM TOY MODEL PAPERS ---
c_s = 343.0       # Speed of sound/light (m/s)
rho_0 = 1.2       # Background density (kg/m^3)
x_01 = 2.4048     # Bessel zero
aspect_ratio = 1.85

# Defect Parameters
a_eq = 1e-3       # 1 mm
L_eq = a_eq * aspect_ratio
mass_effective = rho_0 * np.pi * a_eq**2 * L_eq

# --- 2. RESONANT FREQUENCY ---
term_radial = (x_01 / a_eq)**2
term_axial = (np.pi / L_eq)**2
omega_natural = c_s * np.sqrt(term_radial + term_axial)
freq_natural = omega_natural / (2 * np.pi)

print(f"--- DEFECT PARAMETERS ---")
print(f"Throat Radius a: {a_eq*1000} mm")
print(f"Natural Resonant Frequency: {freq_natural/1000:.2f} kHz")

# --- 3. MODULE 1: JELLIFICATION (CORRECTED DYNAMICS) ---
def simulate_throat_dynamics(drive_freq_ratio, drive_force_amplitude, cycles=50):
    period = 1.0 / freq_natural
    dt = period / 100.0 # Higher resolution
    total_time = cycles * period
    t = np.arange(0, total_time, dt)

    omega_drive = drive_freq_ratio * omega_natural

    r = np.zeros_like(t)
    v = np.zeros_like(t)
    r[0] = a_eq

    # Damping: gamma has units of 1/s (frequency)
    # Equation: x'' + gamma*x' + w^2*x = F/m
    gamma = 0.05 * omega_natural

    # Spring constant k
    k_spring = mass_effective * omega_natural**2

    for i in range(1, len(t)):
        # Forces
        f_restore = -k_spring * (r[i-1] - a_eq)
        f_drive = drive_force_amplitude * np.cos(omega_drive * t[i-1])

        # Acceleration = Force/m - gamma*v
        accel = (f_restore + f_drive) / mass_effective - gamma * v[i-1]

        # Semi-implicit Euler
        v[i] = v[i-1] + accel * dt
        r[i] = r[i-1] + v[i] * dt

    return t, r

# Run Simulation
# Force: 1% of max restoring force
force_scale = (mass_effective * omega_natural**2) * (0.05 * a_eq)

# Scenario A: Normal (10% freq)
t_normal, r_normal = simulate_throat_dynamics(0.1, force_scale)

# Scenario B: Hutchison (100% freq)
t_jelly, r_jelly = simulate_throat_dynamics(1.0, force_scale)

# --- 4. MODULE 2: LEVITATION ---
g_gravity = 9.81
lambda_RF = c_s / freq_natural
k_RF = 2 * np.pi / lambda_RF

# Levitation Condition
# Force density f ~ k * E_density
E_density_required = (rho_0 * g_gravity) / k_RF
Delta_P_required = np.sqrt(2 * rho_0 * c_s**2 * E_density_required)
Delta_rho_fraction = Delta_P_required / (rho_0 * c_s**2)

print(f"\n--- LEVITATION REQUIREMENTS ---")
print(f"Wavelength: {lambda_RF*1000:.2f} mm")
print(f"Required Density Perturbation: {Delta_rho_fraction:.2e}")

# --- 5. FUSION CHECK ---
lattice_spacing = 3e-10
max_oscillation = np.max(np.abs(r_jelly - a_eq))
print(f"\n--- FUSION CRITERIA ---")
print(f"Max Oscillation: {max_oscillation:.2e} m")
print(f"Lattice Spacing: {lattice_spacing:.2e} m")
ratio = max_oscillation / lattice_spacing
print(f"Ratio: {ratio:.1f}")
if ratio > 1.0:
    print("RESULT: JELLIFICATION CONFIRMED")
else:
    print("RESULT: STABLE")

# --- PLOTTING ---
plt.figure(figsize=(10, 8))

# Jellification Plot
plt.subplot(2, 1, 1)
plt.plot(t_normal*1000, (r_normal - a_eq)*1e6, label='Off-Resonance', color='blue', alpha=0.5)
plt.plot(t_jelly*1000, (r_jelly - a_eq)*1e6, label='Resonance', color='red')
# Draw Lattice Spacing limit
plt.axhline(y=lattice_spacing*1e6, color='k', linestyle='--', label='Atomic Lattice Limit')
plt.axhline(y=-lattice_spacing*1e6, color='k', linestyle='--')
plt.title("Hutchison Effect: Defect Throat Destabilization")
plt.ylabel("Radius Deviation (microns)")
plt.xlabel("Time (ms)")
plt.legend()
plt.grid(True, alpha=0.3)

# Levitation Plot
z = np.linspace(0, lambda_RF, 200)
pot_g = rho_0 * g_gravity * z
pot_rad = - E_density_required * 5.0 * np.cos(k_RF * z)**2
# Align means for visual
pot_g -= np.mean(pot_g)
pot_rad -= np.mean(pot_rad)
pot_tot = pot_g + pot_rad

plt.subplot(2, 1, 2)
plt.plot(z*1000, pot_g, '--', color='gray', label='Gravity Potential')
plt.plot(z*1000, pot_rad, color='purple', label='RF Acoustic Potential')
plt.plot(z*1000, pot_tot, color='green', linewidth=2, label='Net Potential')
plt.fill_between(z*1000, pot_tot, max(pot_tot), where=(np.gradient(pot_tot)>0), color='green', alpha=0.1)
plt.title("Levitation Potential Wells")
plt.xlabel("Height (mm)")
plt.ylabel("Potential")
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('hutchison_results.png')
