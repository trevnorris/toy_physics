import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# SUPERFLUID DEFECT SIMULATION: BINARY ORBIT (1PN) - UPDATED
# ==============================================================================
# Objective: Simulate the trajectory of a test mass using the FULL tensor
#            structure derived in Paper III.
#
# Updates:
# 1. Vector Sector now includes Longitudinal (Ram Pressure) term.
# 2. Coefficients match the derived alpha-tuned values (7/2 and 1/2).
# ==============================================================================

# --- Constants & Units ---
G = 1.0
M = 1.0
c = 100.0  # Speed of light/sound
beta = 1.5 # Inertial correction from Paper I

# Orbit Parameters
a = 1000.0
e = 0.5
r_peri = a * (1 - e)
v_peri = np.sqrt(G * M * (1 + e) / r_peri)

# Simulation Parameters
dt = 0.1
num_orbits = 5
steps = int(num_orbits * 2 * np.pi * np.sqrt(a**3 / (G*M)) / dt)

# --- State Initialization ---
pos = np.array([r_peri, 0.0])
vel = np.array([0.0, v_peri])

traj_x = []
traj_y = []

# --- Force Functions ---

def get_acceleration(pos, vel):
    r_sq = np.dot(pos, pos)
    r = np.sqrt(r_sq)
    n_vec = pos / r  # Radial unit vector (n)

    # Velocity components
    v_sq = np.dot(vel, vel)       # (v . v)
    v_dot_n = np.dot(vel, n_vec)  # (v . n) Longitudinal component

    # 1. Newtonian Acceleration
    acc_newton = - (G * M / r_sq) * n_vec

    # 2. Scalar Lag Acceleration (Paper I)
    # The lag field generates a 1/r^3 correction
    scalar_prefactor = 3.0 * (G * M) / (c**2)
    acc_scalar = - (scalar_prefactor * G * M / (r**3)) * n_vec

    # 3. Vector Magnetic Acceleration (Paper III - UPDATED)
    # Interaction Potential V ~ (G M / r) * [ 3.5 (v.v) + 0.5 (v.n)^2 ] / c^2
    # Force F = -Grad(V). The gradient of 1/r gives -1/r^2 * n_vec.
    # Note: We are simulating the force *derived* from the interaction energy.
    # (Assuming the velocity terms are effectively constant during the gradient op)

    coeff_parallel = 3.5      # From Transverse Vortex (tuned)
    coeff_longitudinal = 0.5  # From Scalar Offset + Longitudinal Dyon

    # The effective "Magnetic" Force Intensity
    # F_mag ~ (GM / r^2 c^2) * [ 3.5 v^2 + 0.5 (v.n)^2 ]
    # We apply this in the radial direction (approximate for 1PN central force)

    mag_intensity = (G * M / r_sq) * (1.0 / c**2) * (
        coeff_parallel * v_sq +
        coeff_longitudinal * (v_dot_n**2)
    )

    acc_vector = - mag_intensity * n_vec

    # Total Force
    acc_total_force = acc_newton + acc_scalar + acc_vector

    # 4. Inertial Correction (Paper I)
    # beta = 3/2 density-dependent mass
    sigma = beta * (G * M) / (c**2 * r)
    inertial_factor = 1.0 - sigma

    return acc_total_force * inertial_factor

# --- Integration Loop ---
print(f"Simulating {num_orbits} orbits with Tensor Force Law...")

acc = get_acceleration(pos, vel)

for i in range(steps):
    vel += 0.5 * acc * dt
    pos += vel * dt
    acc = get_acceleration(pos, vel)
    vel += 0.5 * acc * dt

    if i % 10 == 0:
        traj_x.append(pos[0])
        traj_y.append(pos[1])

# --- Visualization ---
plt.figure(figsize=(10, 10))
plt.plot(traj_x, traj_y, linewidth=0.5, label='Trajectory')
plt.scatter([0], [0], color='orange', s=100, label='Central Mass')
plt.axis('equal')
plt.grid(True, alpha=0.3)
plt.title(f"Superfluid Dyon Orbit (1PN)\nFull Tensor: 3.5(v.v) + 0.5(v.n)^2")
plt.legend()
plt.xlabel("x")
plt.ylabel("y")
plt.show()
