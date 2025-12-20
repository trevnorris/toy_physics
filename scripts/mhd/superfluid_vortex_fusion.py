import numpy as np

def calculate_vortex_yield():
    """
    Compares the energy yield of Standard Nuclear Fusion vs.
    Superfluid Vortex Reconnection (Topological Fusion).

    Model: Two vortex rings (defects) merge into one larger ring.
    Constraint: Conservation of Impulse (Momentum) is standard for superfluid vortex reconnection.
    """

    # Constants
    # We work in ratios, so exact density (rho) and circulation (Gamma) cancel out
    # for efficiency calculations, but we track them for scaling.

    print("--- VORTEX FUSION EFFICIENCY AUDIT ---")

    # 1. SETUP: TWO PROTONS (Two Vortex Rings)
    # Assume the 'particle' is a vortex ring of radius R0
    # Energy of a vortex ring E ~ R * (ln(8R/a) - 1.61)
    # Impulse (Momentum) P ~ R^2

    R0 = 1.0  # Normalized radius of initial particles
    N = 2     # Number of particles

    # Initial Energy (System of 2 separated rings)
    # We ignore the log factor for the linear tension approximation (E propto Length)
    # Detailed form: E_initial = 2 * (2 * pi * R0) * Tension
    Length_initial = N * (2 * np.pi * R0)

    # 2. THE MERGER (Reconnection)
    # When vortex rings merge, they conserve IMPULSE (P ~ Area ~ R^2).
    # P_initial = 2 * (pi * R0^2)
    # P_final = pi * R_final^2
    # Conservation: 2 * pi * R0^2 = pi * R_final^2
    # => R_final = sqrt(2) * R0

    R_final = np.sqrt(2) * R0

    # Final Energy (One large ring)
    Length_final = 1 * (2 * np.pi * R_final)

    # 3. THE YIELD CALCULATION
    # Radiated Energy = E_initial - E_final
    # Since Energy is proportional to Length (String Tension approximation)

    delta_Length = Length_initial - Length_final
    efficiency = delta_Length / Length_initial

    print(f"Initial Defect Length: {Length_initial:.4f} units")
    print(f"Final Defect Length:   {Length_final:.4f} units")
    print(f"Length Reduced:        {delta_Length:.4f} units")
    print("-" * 30)
    print(f"Topological Mass Defect (Efficiency): {efficiency * 100:.2f}%")

    # 4. COMPARISON WITH NUCLEAR FUSION
    # Deuterium-Tritium Fusion releases ~17.6 MeV from a mass of ~4700 MeV.
    # Efficiency ~ 0.37% (approx 0.4 - 0.7% depending on the reaction chain).

    fusion_efficiency = 0.007 # 0.7% optimistic upper bound for stellar fusion

    print(f"Standard Nuclear Fusion Efficiency:     {fusion_efficiency * 100:.2f}%")

    ratio = efficiency / fusion_efficiency
    print("-" * 30)
    print(f"GAIN FACTOR: {ratio:.1f}x more energy per reaction than fusion.")

calculate_vortex_yield()

"""
Output:

--- VORTEX FUSION EFFICIENCY AUDIT ---
Initial Defect Length: 12.5664 units
Final Defect Length:   8.8858 units
Length Reduced:        3.6806 units
------------------------------
Topological Mass Defect (Efficiency): 29.29%
Standard Nuclear Fusion Efficiency:     0.70%
------------------------------
GAIN FACTOR: 41.8x more energy per reaction than fusion.
"""
