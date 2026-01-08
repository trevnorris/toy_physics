import numpy as np
import matplotlib.pyplot as plt

# --- CONSTANTS & PARAMETERS ---
c_s = 343.0       # Speed of sound/light (m/s)
rho_0 = 1.2       # Background density (kg/m^3)
x_01 = 2.4048     # Bessel zero
aspect_ratio = 1.85
g_gravity = 9.81

# Material Parameters (Standard Defect)
a_eq = 1e-3       # 1 mm radius
L_eq = a_eq * aspect_ratio
mass_effective = rho_0 * np.pi * a_eq**2 * L_eq

# Natural Frequency Calculation
term_radial = (x_01 / a_eq)**2
term_axial = (np.pi / L_eq)**2
omega_natural = c_s * np.sqrt(term_radial + term_axial)
freq_natural = omega_natural / (2 * np.pi)

# --- SIMULATION 1: RESONANCE SPECTRUM (FREQUENCY SWEEP) ---
# We want to see how "tuned" the RF field needs to be.
def get_max_amplitude(drive_freq_ratio, drive_force_scale=0.01):
    # drive_force_scale is relative to restoring force at a_eq
    omega_drive = drive_freq_ratio * omega_natural
    k_spring = mass_effective * omega_natural**2
    F_drive = k_spring * a_eq * drive_force_scale

    # Analytical solution for steady-state amplitude of driven damped oscillator
    # A = F_0 / sqrt( (k - m*w^2)^2 + (gamma*w)^2 )
    gamma = 0.05 * omega_natural * mass_effective # Damping coefficient (force/velocity)

    # Denominator terms
    term1 = (k_spring - mass_effective * omega_drive**2)**2
    term2 = (gamma * omega_drive / mass_effective)**2 # Gamma here needs to be carefully handled.
    # Let's use standard form: x'' + gamma_damp*x' + w0^2*x = F/m
    # Amplitude = (F/m) / sqrt( (w0^2 - w^2)^2 + (gamma_damp*w)^2 )

    gamma_damp = 0.05 * omega_natural # Damping factor per unit mass
    force_per_mass = F_drive / mass_effective

    amplitude = force_per_mass / np.sqrt( (omega_natural**2 - omega_drive**2)**2 + (gamma_damp * omega_drive)**2 )
    return amplitude

freq_ratios = np.linspace(0.5, 1.5, 200)
amplitudes = [get_max_amplitude(r) for r in freq_ratios]

# --- SIMULATION 2: LEVITATION DYNAMICS (Particle in Trap) ---
# Simulate a mass dropped into the standing wave field.
def simulate_levitation_trajectory():
    # Standing Wave Parameters
    lambda_RF = c_s / freq_natural
    k_RF = 2 * np.pi / lambda_RF

    # Set power just above threshold (e.g., 1.2x gravity)
    # Force profile: F_rad(z) = F_max * sin(2*k*z)  (Gradient of cos^2)
    # We want max upward force to be 1.2 * mg
    F_max = 2.0 * mass_effective * g_gravity

    t_sim = np.linspace(0, 0.5, 1000)
    dt = t_sim[1] - t_sim[0]

    z = np.zeros_like(t_sim)
    v = np.zeros_like(t_sim)

    # Start slightly off the node to see oscillation
    z[0] = lambda_RF / 8.0 # Max force position for sin(2kz) is pi/4 / k ?
    # sin(2*k*z). Max at 2kz = pi/2 -> z = pi/4k = lambda/8.
    # Let's start it at the bottom of a potential well (node) with a kick, or drop it.
    # Let's drop it from a node position.
    z[0] = lambda_RF / 4.0 + 0.0001 # Start near a node

    for i in range(1, len(t_sim)):
        # Forces
        f_grav = -mass_effective * g_gravity
        # Radiation Force: Derived from potential U ~ cos^2(kz). Force ~ -dU/dz ~ sin(2kz)
        # We model it to push UP at certain zones.
        f_rad = F_max * np.sin(2 * k_RF * z[i-1])

        accel = (f_grav + f_rad) / mass_effective

        v[i] = v[i-1] + accel * dt
        z[i] = z[i-1] + v[i] * dt

    return t_sim, z, lambda_RF

t_lev, z_lev, wavelength = simulate_levitation_trajectory()


# --- SIMULATION 3: POWER SCALING ---
# Energy Density required vs Object Radius
radii = np.logspace(-4, -2, 50) # 0.1mm to 10cm
# For each radius, calculate natural freq, then wavelength, then required E_density
def get_power_req(r):
    L = r * aspect_ratio
    m = rho_0 * np.pi * r**2 * L

    # Freq
    w = c_s * np.sqrt( (x_01/r)**2 + (np.pi/L)**2 )
    freq = w / (2*np.pi)

    # Wavelength
    lam = c_s / freq
    k = 2*np.pi/lam

    # Required Energy Density to support mass
    # F_grad ~ k * E_dens > rho * g
    # Note: Force is per unit volume in the continuum limit.
    # If we treat the object as a point mass, Force = Volume * Force_density.
    # So the condition is actually independent of mass size if density is constant!
    # Force_density_rad > Force_density_grav
    # k * E_dens > rho_object * g
    # Wait, the object density is effectively rho_0 in this toy model (defect in fluid).
    # So E_dens > (rho_0 * g) / k

    E_req = (rho_0 * g_gravity) / k
    return E_req

power_reqs = [get_power_req(r) for r in radii]


# --- PLOTTING ---
plt.figure(figsize=(15, 5))

# Plot 1: Resonance Spectrum
plt.subplot(1, 3, 1)
plt.plot(freq_ratios, np.array(amplitudes)/a_eq, color='red', linewidth=2)
plt.axvline(x=1.0, linestyle='--', color='k', alpha=0.5)
plt.axhline(y=1.0, linestyle='--', color='blue', label='Disruption Threshold (a_eq)')
plt.title("Jellification Spectrum\n(Amplitude vs Freq Ratio)")
plt.xlabel("Frequency Ratio (f / f_0)")
plt.ylabel("Normalized Amplitude (A / a)")
plt.legend()
plt.grid(True, alpha=0.3)

# Plot 2: Levitation Dynamics
plt.subplot(1, 3, 2)
plt.plot(t_lev*1000, z_lev*1000, color='green')
# Draw "Shelves" (Nodes of standing wave)
node_spacing = wavelength / 2.0
for n in range(int(max(z_lev)/node_spacing) + 2):
    plt.axhline(y=n*node_spacing*1000, color='k', linestyle=':', alpha=0.3)
plt.title("Levitation Trajectory\n(Particle in Acoustic Trap)")
plt.xlabel("Time (ms)")
plt.ylabel("Height (mm)")
plt.grid(True, alpha=0.3)

# Plot 3: Scaling Law
plt.subplot(1, 3, 3)
plt.loglog(radii*1000, power_reqs, color='purple', linewidth=2)
plt.title("Power Scaling\n(Energy Density vs Object Size)")
plt.xlabel("Object Radius (mm)")
plt.ylabel("Required Vacuum Energy Density (J/m^3)")
plt.grid(True, which="both", alpha=0.3)

plt.tight_layout()
plt.savefig('hutchison_detailed_analysis.png')
