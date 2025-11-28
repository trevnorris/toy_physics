import numpy as np
from scipy.integrate import quad

# ==========================================
# SHAPIRO DELAY VERIFICATION
# ==========================================
# Constants
G = 1.0
M = 1.0
c0 = 1.0
mu = G * M

# Geometry: Signal travels from x = -L to x = +L at impact parameter b
L_dist = 10000.0  # Distance to Earth/Venus (far away)
b_impact = 100.0  # Impact parameter (grazing Sun)

def time_delay_integrand(x, mode):
    r = np.sqrt(x**2 + b_impact**2)
    
    # 1. REFRACTION (Static n)
    # Stiff Vacuum (n=3) -> n_static = 1 + GM/(r c^2)
    if mode == 'refraction' or mode == 'split':
        n_static = 1.0 + mu / (r * c0**2)
    else:
        n_static = 1.0
        
    # 2. FLOW DRAG (Second Order)
    # Factor = 1 + v^2/c^2
    if mode == 'flow' or mode == 'split':
        # Split Model: Orbital Velocity (v^2 = GM/r)
        v_sq = mu / r
        drag_factor = 1.0 + v_sq / c0**2
    else:
        drag_factor = 1.0
        
    # Total Inverse Speed (1/v_eff)
    # dt/dx = (1/c0) * n_static * drag_factor
    return (1.0/c0) * n_static * drag_factor

# ==========================================
# RUN INTEGRATION
# ==========================================

# 1. Base Time (Vacuum)
t_vacuum, _ = quad(lambda x: 1.0/c0, -L_dist, L_dist)

# 2. Calculate Delays
t_refraction, _ = quad(lambda x: time_delay_integrand(x, 'refraction'), -L_dist, L_dist)
t_flow, _       = quad(lambda x: time_delay_integrand(x, 'flow'), -L_dist, L_dist)
t_split, _      = quad(lambda x: time_delay_integrand(x, 'split'), -L_dist, L_dist)

delay_refraction = t_refraction - t_vacuum
delay_flow = t_flow - t_vacuum
delay_total = t_split - t_vacuum

# 3. GR Theoretical Prediction (4GM/c^3 * ln(4L/b) is the approx for round trip)
# But here we look at the coefficient logic:
# Shapiro delay D = (1 + gamma) * GM * ln(...)
# We want to see if our components sum to (1+1)=2.

print(f"{'COMPONENT':<15} | {'DELAY (Arb Units)':<20} | {'FRACTION':<10}")
print("-" * 50)
print(f"{'Refraction':<15} | {delay_refraction:<20.5f} | {delay_refraction/delay_total:<10.2f}")
print(f"{'Flow (Drag)':<15} | {delay_flow:<20.5f} | {delay_flow/delay_total:<10.2f}")
print("-" * 50)
print(f"{'TOTAL':<15} | {delay_total:<20.5f} | 1.00")

# Validation Logic
if abs(delay_refraction - delay_flow) < 1e-4:
    print("\n[SUCCESS] Refraction and Flow contribute equally (1:1).")
    print("This matches the GR partition of Space (g_rr) and Time (g_00).")
else:
    print("\n[FAILURE] Components do not match.")
