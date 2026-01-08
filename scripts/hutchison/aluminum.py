import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# --- PHYSICS CONSTANTS ---
g = 9.81
rho_medium = 1.2    # Air density (kg/m^3)
c_s = 343.0         # Speed of sound (m/s)

# --- EXPERIMENTAL TARGET: 1cm ALUMINUM CUBE ---
mass = 0.0027       # 2.7 grams
side = 0.01         # 1 cm
area = side**2      # Bottom surface area

# Resonance Calculation (From Paper V)
# 28.89 kHz is the "Jellification" freq. We use this for the levitation carrier.
freq = 28890       
wavelength = c_s / freq
k = 2 * np.pi / wavelength

# --- POWER SETTINGS ---
# We need > 169.7 dB to lift. 
# Let's use 172 dB to give it a "Safety Factor" of ~1.7x Gravity.
input_db = 172.0 

# Calculate Forces
p_ref = 20e-6
pressure_rms = p_ref * (10**(input_db/20))
P_rad = (pressure_rms**2) / (rho_medium * c_s**2) # Radiation Pressure (Pa)
max_lift_force = P_rad * area
gravity_force = mass * g

print(f"--- PRE-FLIGHT CHECK ---")
print(f"Gravity Load: {gravity_force*1000:.2f} mN")
print(f"Max Lift Force: {max_lift_force*1000:.2f} mN")
print(f"Safety Factor: {max_lift_force/gravity_force:.2f}x")

# --- SIMULATION STATE ---
# Start AT a node height (approx lambda/4) to simulate stable capture
z = wavelength / 4.0 
v = 0.0
dt = 0.0001 # Fine time steps

# --- VISUALIZATION SETUP ---
fig, ax = plt.subplots(figsize=(6, 8))
ax.set_ylim(0, 0.05) # Zoom in to bottom 5cm
ax.set_xlim(-0.03, 0.03)
ax.set_title(f"Hutchison Experiment: 1cm Al Cube\nInput: {input_db} dB | Freq: {freq/1000:.1f} kHz")
ax.set_ylabel("Height (m)")
ax.set_xlabel("Vacuum Pressure Intensity")

# Draw the Levitation Shelves (Green Zones are Safe)
y_space = np.linspace(0, 0.05, 1000)
# Force profile is sin(2kz). Positive regions are "Push Up" zones.
force_profile = np.sin(2 * k * y_space)
# Mask only the upward pushing zones
ax.fill_betweenx(y_space, -0.03, -0.03 + 0.06*(force_profile > 0), color='green', alpha=0.1, label='Lift Zones')
ax.axhline(y=0, color='k', linewidth=3) # Floor

# The Cube
cube, = ax.plot([], [], 's', markersize=25, color='silver', markeredgecolor='black', label='Aluminum')
# Force Vector Arrow
arrow, = ax.plot([], [], 'r-', linewidth=2)
txt = ax.text(-0.025, 0.045, "", fontsize=9, bbox=dict(facecolor='white', alpha=0.9))

def update(frame):
    global z, v
    
    # 1. Forces
    f_g = -gravity_force
    
    # Levitation Force depends on height in the standing wave
    # F = F_max * sin(2kz)
    f_lev = max_lift_force * np.sin(2 * k * z)
    
    # Damping (Air Resistance) - CRITICAL for stability
    # F_drag = -b * v
    f_drag = -0.05 * v 
    
    # 2. Dynamics
    f_net = f_g + f_lev + f_drag
    accel = f_net / mass
    
    v += accel * dt
    z += v * dt
    
    # 3. Floor Boundary (Hard Stop)
    if z < 0: 
        z = 0
        v = 0 # Sticky floor if it crashes
    
    # 4. Update Graphics
    cube.set_data([0], [z])
    
    # Draw arrow representing Net Force
    arrow_scale = 0.005 # Visual scale
    arrow.set_data([0, 0], [z, z + (f_net/gravity_force)*arrow_scale])
    arrow.set_color('green' if f_net > 0 else 'red')
    
    status = "HOVERING" if z > 0.001 else "GROUNDED"
    txt.set_text(f"Height: {z*1000:.1f} mm\nNet Force: {f_net*1000:.1f} mN\nStatus: {status}")
    
    return cube, arrow, txt

ani = animation.FuncAnimation(fig, update, frames=200, interval=20, blit=True)
plt.show()
