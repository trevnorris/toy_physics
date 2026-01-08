import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# --- PHYSICS SETUP ---
g = 9.81
mass = 1.0  # Normalized mass
freq_sim = 1.0 # Normalized frequency for visualization scale
wavelength = 1.0
k = 2 * np.pi / wavelength

# LEVITATION POWER
# We need enough force to counteract gravity (mg).
# The Ponderomotive Force F_p = -grad(Potential).
# If Potential U = A * cos^2(kz), then Max Force = A * k.
# We set A so that Max Force = 1.5 * mg (enough to lift).
potential_amplitude = (1.5 * mass * g) / k

# --- SIMULATION STATE ---
# We simulate the "Drift" (averaged movement), not the 40kHz vibration
particle_z = 2.5 * wavelength # Start high up
particle_v = 0.0
dt = 0.01 # Smooth time step (Real-time speed)

# --- VISUALIZATION SETUP ---
fig, ax = plt.subplots(figsize=(6, 8))
ax.set_xlim(-0.5, 0.5)
ax.set_ylim(0, 3 * wavelength)
ax.set_title("Hutchison Effect: Acoustic Trap (Time-Averaged)")
ax.set_ylabel("Height (Arbitrary Units)")
ax.get_xaxis().set_visible(False) # Hide X axis (1D physics)

# 1. Draw the "Invisible" Potential Wells (The Shelves)
z_background = np.linspace(0, 3*wavelength, 500)
# This curve shows where the vacuum pressure is strongest
pressure_field = potential_amplitude * np.cos(k * z_background)**2
# Normalize for display
pressure_display = (pressure_field / np.max(pressure_field)) * wavelength * 0.5
ax.plot(pressure_display, z_background, color='blue', alpha=0.3, label='Vacuum Pressure Field')
ax.fill_betweenx(z_background, 0, pressure_display, color='blue', alpha=0.1)

# 2. The Levitation Node Markers
for n in range(1, 7):
    node_height = (n * wavelength) / 2.0 - (wavelength/4.0) # Approx trap locations
    if node_height > 0:
        ax.axhline(y=node_height, color='green', linestyle=':', alpha=0.5)

# 3. The Particle
ball, = ax.plot([], [], 'ro', markersize=12, label='Matter Sample')
text_status = ax.text(0.05, 0.95, "", transform=ax.transAxes)

def update(frame):
    global particle_z, particle_v

    # --- PONDEROMOTIVE PHYSICS ---
    # 1. Gravity Force
    F_grav = -mass * g

    # 2. Levitation Force (Gradient of the Standing Wave Pressure)
    # F = -dU/dz.  U = A * cos^2(kz).  dU/dz = -A * k * sin(2kz)
    # So F_lev = + A * k * sin(2kz)
    F_lev = potential_amplitude * k * np.sin(2 * k * particle_z)

    # 3. Drag (Air Resistance) - stabilizes the bounce
    F_drag = -1.0 * particle_v

    # Newton's Law
    accel = (F_grav + F_lev + F_drag) / mass

    particle_v += accel * dt
    particle_z += particle_v * dt

    # Floor collision
    if particle_z < 0:
        particle_z = 0
        particle_v = -particle_v * 0.5

    # Update Graphic
    ball.set_data([0], [particle_z])

    # Status Text
    force_ratio = abs(F_lev / F_grav)
    status = "FALLING"
    if force_ratio > 1.0 and abs(particle_v) < 1.0:
        status = "LEVITATING"
    text_status.set_text(f"Status: {status}\nHeight: {particle_z:.2f}")

    return ball, text_status

ani = animation.FuncAnimation(fig, update, frames=300, interval=20, blit=True)
plt.show()
