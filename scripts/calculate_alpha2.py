import numpy as np
from scipy import integrate

# ==========================================
# 1. SETUP & PARAMETERS
# ==========================================
a = 1.0       # Throat radius
L = 2.0       # Throat depth
rho0 = 1.0    # Reference density
epsilon = 0.1 # Deformation parameter

def legendre_p2(x):
    return 0.5 * (3 * x**2 - 1)

# ==========================================
# 2. DEFINE GEOMETRY (The Funnel)
# ==========================================
def radius_at_w(w):
    """
    Defines the radius of the throat on the brane at a given bulk depth w.
    Implements the 'Rounded/Flared' funnel shape.
    """
    # Example: A cylinder of radius 'a' that flares out near the mouth (w=0)
    # This matches the logic: current_radius = a + (0.5a) * exp(...)
    flared_radius = a + (0.5 * a) * np.exp(-w / (0.2 * a))
    return flared_radius

# ==========================================
# 3. DEFINE INTEGRANDS
# ==========================================
# We integrate in 'Cylindrical-like' coordinates:
# r     : Radius on the brane (0 to R(w))
# theta : Angle on the brane (0 to pi)
# w     : Depth in the bulk (0 to L)

def rho_val(r, theta, w):
    # The density includes the epsilon asymmetry term
    return rho0 * (1 + epsilon * legendre_p2(np.cos(theta)))

def mass_integrand(r, theta, w):
    # dV = r^2 * sin(theta) * dr * dtheta * dw
    return rho_val(r, theta, w) * r**2 * np.sin(theta)

def quad_integrand(r, theta, w):
    # Quadrupole moment includes an extra r^2 * P2(cos theta) term
    return rho_val(r, theta, w) * r**2 * legendre_p2(np.cos(theta)) * r**2 * np.sin(theta)

# ==========================================
# 4. DEFINE DYNAMIC LIMITS
# ==========================================
# nquad expects range functions to accept arguments for the outer variables.
# Order of integration: r, theta, w

def range_w():
    return [0, L]

def range_theta(w):
    return [0, np.pi]

def range_r(theta, w):
    # This is the Key Fix:
    # The limit for 'r' is determined by the geometry function R(w).
    # No if-statements needed!
    return [0, radius_at_w(w)]

# ==========================================
# 5. PERFORM INTEGRATION
# ==========================================
print(f"--- Starting Calculation (Finite Limits) ---")

# We pass the functions in the order of variables: [r, theta, w]
# ranges must match this order.
ranges = [range_r, range_theta, range_w]

print("Integrating Mass...")
# nquad handles the nested dependencies automatically
val_M, err_M = integrate.nquad(mass_integrand, ranges)
mass_total = val_M * 2 * np.pi  # Multiply by 2pi for phi integral

print("Integrating Quadrupole Q...")
val_Q, err_Q = integrate.nquad(quad_integrand, ranges)
quadrupole_total = val_Q * 2 * np.pi

print(f"Total Mass M = {mass_total:.6f}")
print(f"Quadrupole Q = {quadrupole_total:.6f}")

# ==========================================
# 6. COMPUTE ALPHA_2
# ==========================================
if mass_total != 0:
    alpha_2 = quadrupole_total / (mass_total * (a**2))
else:
    alpha_2 = 0

print("-" * 30)
print(f"CALCULATED RESULT:")
print(f"alpha_2 = {alpha_2:.6f}")
print("-" * 30)

"""
Output:

--- Starting Calculation (Finite Limits) ---
Integrating Mass...
Integrating Quadrupole Q...
Total Mass M = 9.983226
Quadrupole Q = 0.143266
------------------------------
CALCULATED RESULT:
alpha_2 = 0.014351
------------------------------
"""
