"""Exploratory throat-spin geometry test.

This file is retained for note-level geometry exploration and is not part of the
paper's required verification bundle. It compares a long 4D cylindrical source
against the expected 3D dipole falloff on the brane.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# ==========================================
# 4D SPINNING THROAT (CYLINDER) TEST
# ==========================================

def test_cylindrical_spin():
    print("NOTE: exploratory geometry script; not a cited manuscript verification harness.")
    print("Computing 4D Cylinder -> 3D Projection...")
    
    # SETUP
    r_values = np.logspace(-1, 2, 20) # From r=0.1 to r=100
    L_throat = 1000.0  # The throat extends deep into the bulk
    
    # PHYSICS
    # We are integrating the 4D Green's function (1/R^2) along the w-axis.
    # Source: A cylinder of radius a=1, extending from w=-L to w=+L.
    # We observe at w=0 on the brane.
    
    def vector_potential_cylinder(r_obs):
        # We integrate the "Ring" result along the w-axis.
        # The distance R^2 becomes r_brane^2 + w^2.
        
        # Inner function: Potential of a slice at height w
        def slice_contribution(w):
            # The effective 3D distance from the observer to the ring slice
            # is sqrt(r_obs^2 + w^2) ???
            # Actually, let's go back to the raw integral.
            # Green's function is 1 / |X - X'|^2
            # |X - X'|^2 = (r_obs - a cos theta)^2 + (a sin theta)^2 + w^2
            #            = r_obs^2 - 2 r_obs a cos theta + a^2 + w^2
            
            # We integrate this 1/Dist^2 over theta (0..2pi) and w (-L..L)
            
            # Let's do the theta integral first (it's the same as before but with an extra w^2 term)
            def theta_integrand(theta):
                dist_sq = r_obs**2 - 2*r_obs*1.0*np.cos(theta) + 1.0**2 + w**2
                return np.cos(theta) / dist_sq # J_y component
            
            val_theta, _ = quad(theta_integrand, 0, 2*np.pi)
            return val_theta

        # Now integrate that slice contribution along w
        # Optimization: The integrand drops fast. We can integrate 0 to L and double it.
        val_w, _ = quad(slice_contribution, 0, L_throat)
        return 2 * val_w

    # COMPUTE PROFILE
    potentials = []
    print("Integrating... (this might take a moment due to double integration)")
    for r in r_values:
        # Optimization for speed: Only integrate deep if r is small? 
        # No, let's just do it.
        potentials.append(vector_potential_cylinder(r))
    
    potentials = np.array(potentials)
    
    # FIT POWER LAW
    # We expect Potential ~ 1/r^2 (because Velocity ~ 1/r^3 ?? Wait.)
    # Let's check the GR Target again.
    # GR Metric Shift (h_0i) is a POTENTIAL. It scales as 1/r^2 (Dipole).
    # Velocity of frame dragging is related to h_0i.
    
    # If the script outputs Slope = -2.0, we match GR Metric Potential.
    
    log_r = np.log10(r_values[-5:])
    log_A = np.log10(potentials[-5:])
    slope, intercept = np.polyfit(log_r, log_A, 1)
    
    print(f"\nRESULTS:")
    print(f"Calculated Power Law Slope: {slope:.4f}")
    
    print("\nINTERPRETATION:")
    print("GR Target for Dipole Potential: -2.0")
    if abs(slope - (-2.0)) < 0.2:
        print(">> MATCH! The Cylinder projects to a perfect 3D Dipole (1/r^2).")
        print(">> This implies the throat must be long (L >> r).")
    else:
        print(f">> MISMATCH. Slope is {slope:.2f}")

    # Plot
    plt.figure(figsize=(8,5))
    plt.loglog(r_values, potentials, 'o-', label='4D Cylinder Potential')
    plt.loglog(r_values, 10*r_values**-2, 'k--', label='GR Target (1/r^2)')
    plt.loglog(r_values, 100*r_values**-3, 'r:', label='Sphere Result (1/r^3)')
    
    plt.xlabel("Distance on Brane (r)")
    plt.ylabel("Gravitomagnetic Potential (h_0i)")
    plt.title(f"The Cylinder Spin Test: Slope = {slope:.2f}")
    plt.legend()
    plt.grid(True)
    backend = plt.get_backend().lower()
    if "agg" in backend:
        print("Plot display skipped on non-GUI backend.")
    else:
        plt.show()


if __name__ == "__main__":
    test_cylindrical_spin()
