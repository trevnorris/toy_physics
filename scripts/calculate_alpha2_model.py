import numpy as np
from scipy import integrate

# ==========================================
# 1. SETUP & PARAMETERS
# ==========================================
# Geometric parameters
a = 1.0       # Throat radius (characteristic size)
L = 2.0       # Throat depth
rho0 = 1.0    # Reference density
epsilon = 0.1 # Deformation parameter (for the test model)

# ==========================================
# 2. DEFINE THE DENSITY MODEL
# ==========================================
def legendre_p2(x):
    """Returns the 2nd Legendre polynomial P2(x) = 0.5 * (3x^2 - 1)"""
    return 0.5 * (3 * x**2 - 1)

def rho_model(r, theta, w):
    """
    The 4D density profile rho(r, theta, w).
    
    Current: The Gaussian Toy Model from Paper Eq. (24) / Appendix A.
    rho ~ exp(-r^2/a^2) * exp(-w^2/L^2) * [1 + eps * P2(cos theta)]
    
    TODO: REPLACE this function with your physical Transition Region model.
    Example: A function that follows bent streamlines near the mouth.
    """
    
    # Radial localization (on brane)
    radial_part = np.exp(-(r**2) / (a**2))
    
    # Bulk localization (in throat)
    bulk_part = np.exp(-(w**2) / (L**2))
    
    # Angular structure (The source of alpha_2)
    # If this term is 0, alpha_2 will be 0.
    angular_part = 1.0 + epsilon * legendre_p2(np.cos(theta))
    
    return rho0 * radial_part * bulk_part * angular_part

# ==========================================
# 3. DEFINE INTEGRANDS
# ==========================================
# Volume element dV in effective 3D:
# We integrate out 'w' effectively, but we do full 4D integral for M and Q.
# Volume element = r^2 * sin(theta) * dr * dtheta * dw * dphi
# (The dphi integral is trivial = 2pi, so we factor it out)

def mass_integrand(w, theta, r):
    """Integrand for Total Mass M"""
    # Density * Jacobian (r^2 sin(theta))
    return rho_model(r, theta, w) * (r**2) * np.sin(theta)

def quadrupole_integrand(w, theta, r):
    """Integrand for Quadrupole Moment Q"""
    # Density * r^2 * P2(cos theta) * Jacobian
    # Total r power = r^4
    return rho_model(r, theta, w) * (r**2) * legendre_p2(np.cos(theta)) * (r**2) * np.sin(theta)

# ==========================================
# 4. PERFORM NUMERICAL INTEGRATION
# ==========================================
print(f"--- Starting Calculation ---")
print(f"Parameters: a={a}, L={L}, epsilon={epsilon}")

# Integration limits
# r: 0 to infinity
# theta: 0 to pi
# w: -infinity to infinity
limits_r = [0, np.inf]
limits_theta = [0, np.pi]
limits_w = [-np.inf, np.inf]

# NOTE: We use nquad for multi-dimensional integration.
# We multiply by 2*pi at the end for the phi integral.

# A. Calculate Mass M
print("Integrating Mass...")
val_M, err_M = integrate.nquad(
    mass_integrand, 
    [limits_w, limits_theta, limits_r]
)
mass_total = val_M * 2 * np.pi
print(f"Total Mass M = {mass_total:.6f}")

# B. Calculate Quadrupole Q
print("Integrating Quadrupole Q...")
val_Q, err_Q = integrate.nquad(
    quadrupole_integrand, 
    [limits_w, limits_theta, limits_r]
)
quadrupole_total = val_Q * 2 * np.pi
print(f"Quadrupole Q = {quadrupole_total:.6f}")

# ==========================================
# 5. COMPUTE ALPHA_2
# ==========================================
# Definition: Phi_quad = -G * Q / r^3
# Expansion:  Phi = -GM/r * [1 + alpha_2 * (a/r)^2 * P2 + ...]
#                   = -GM/r - GM * alpha_2 * a^2 / r^3 * P2
# Matching terms: Q = M * alpha_2 * a^2
# Therefore: alpha_2 = Q / (M * a^2)

if mass_total != 0:
    alpha_2 = quadrupole_total / (mass_total * (a**2))
else:
    alpha_2 = 0

print("-" * 30)
print(f"CALCULATED RESULT:")
print(f"alpha_2 = {alpha_2:.6f}")
print("-" * 30)

# ==========================================
# 6. VERIFICATION (Against Paper Analytical Result)
# ==========================================
# From Appendix A: Q/M = (3/10) * epsilon * a^2
# So theoretical alpha_2 = (3/10) * epsilon
theoretical_alpha2 = 0.3 * epsilon

print(f"VERIFICATION:")
print(f"Theoretical alpha_2 (Gaussian Model) = {theoretical_alpha2:.6f}")
print(f"Difference = {abs(alpha_2 - theoretical_alpha2):.6e}")

if abs(alpha_2 - theoretical_alpha2) < 1e-5:
    print(">> SUCCESS: Numerical result matches analytical theory.")
else:
    print(">> WARNING: Discrepancy detected.")
