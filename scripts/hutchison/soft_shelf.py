import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# --- PHYSICS ---
g = 9.81
mass = 0.0027 
c_s = 343.0
freq = 28890       
wavelength = c_s / freq
k = 2 * np.pi / wavelength

# Power Settings (172 dB)
input_db = 172.0 
p_ref = 20e-6
pressure_rms = p_ref * (10**(input_db/20))
rho_medium = 1.2
P_rad = (pressure_rms**2) / (rho_medium * c_s**2)
max_lift_force = P_rad * (0.01**2) # 1cm^2 area

# --- SCENARIO SETUP ---
# We simulate 4 Cubes with different starting conditions
# State format: [z_position, velocity, color_code]
cubes = [
    {'z': wavelength * 0.25, 'v': 0.0, 'c': 'green', 'label': 'Shelf 1 (Sit)'},    # Sitting in bottom trap
    {'z': wavelength * 1.25, 'v': 0.0, 'c': 'blue',  'label': 'Shelf 3 (Sit)'},    # Sitting in high trap
    {'z': wavelength * 0.85, 'v': 0.0, 'c': 'orange','label': 'Shelf 2 (Drop)'},   # Short drop into Shelf 2
    {'z': wavelength * 2.50, 'v': 0.0, 'c': 'red',   'label': 'Roof Drop'}         # Long drop (High Energy)
]

dt = 0.00005 # Ultra fine time step for collisions

# --- VISUALIZATION ---
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_ylim(0, 0.04) # Show bottom 4cm
ax.set_xlim(0, 5)    # 4 lanes
ax.set_title(f"Stability Test: The 'Soft Mattress' Effect\nPower: {input_db} dB")
ax.set_ylabel("Height (m)")
ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels([c['label'] for c in cubes])

# Draw the Potential Field (The "Mattresses")
z_space = np.linspace(0, 0.05, 1000)
force_profile = np.sin(2 * k * z_space)
# Draw Green zones where Force is UP
for i in range(1, 5):
    # Shift the background for each lane
    lane_center = i
    # Visualize the "Trap" zones
    ax.fill_betweenx(z_space, lane_center-0.4, lane_center+0.4, where=(force_profile>0), color='green', alpha=0.1)

# Create Plot Objects
plots = []
for i, c in enumerate(cubes):
    p, = ax.plot([], [], 's', markersize=20, color=c['c'], markeredgecolor='black')
    plots.append(p)

def update(frame):
    for i, cube in enumerate(cubes):
        # Physics Engine
        f_g = -mass * g
        f_lev = max_lift_force * np.sin(2 * k * cube['z'])
        f_drag = -0.1 * cube['v'] # Air resistance
        
        accel = (f_g + f_lev + f_drag) / mass
        
        cube['v'] += accel * dt
        cube['z'] += cube['v'] * dt
        
        # Floor
        if cube['z'] < 0:
            cube['z'] = 0
            cube['v'] = 0
            
        # Update Plot
        plots[i].set_data([i+1], [cube['z']])
        
    return plots

ani = animation.FuncAnimation(fig, update, frames=300, interval=20, blit=True)
plt.show()
