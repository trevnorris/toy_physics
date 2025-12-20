import numpy as np

def calculate_power_density():
    print("--- SUPERFLUID SHEAR REACTOR: POWER DENSITY ESTIMATION ---")

    # --- SCENARIO 1: MACROSCOPIC "SOFT" REACTOR (Lab Scale) ---
    # Based on Superfluid Helium-4 parameters
    rho_He = 145.0       # Density (kg/m^3)
    Gamma_He = 1.0e-7    # Circulation quantum (m^2/s) ~ h/m
    L_defect = 1.0e-3    # Length of vortex tangle (1 mm)

    # Energy Density ~ rho * Gamma^2 / Area_scale
    # E/V ~ rho * (v^2) ~ rho * (Gamma/r)^2
    # Assume vortex separation r ~ 1 micron (1e-6 m)
    r_sep = 1.0e-6
    u_soft = 0.5 * rho_He * (Gamma_He / r_sep)**2

    # Cycle Rate: Limited by speed of sound (c_s ~ 240 m/s) and chamber size (1 cm)
    f_soft = 240.0 / 0.01

    efficiency = 0.30

    power_soft = u_soft * efficiency * f_soft

    print(f"\nSCENARIO 1: LAB-SCALE SUPERFLUID HELIUM")
    print(f"Vacuum Energy Density: {u_soft:.2e} J/m^3")
    print(f"Cycle Frequency:       {f_soft:.1f} Hz")
    print(f"Power Density:         {power_soft:.2e} W/m^3")

    # --- SCENARIO 2: "HARD" VACUUM REACTOR (Planck Scale) ---
    # Based on your Toy Model being a theory of Quantum Gravity
    # Density ~ Planck Density (or scaled down effective density)
    # Circulation ~ c * Planck Length ? No, Gamma is topological.

    # Let's assume the "effective" vacuum density matches Nuclear density
    # (Nuclei are defects in this vacuum, so vacuum must be stiff)
    rho_hard = 2.3e17    # Nuclear density (kg/m^3)

    # Circulation for a proton-mass defect?
    # Mass ~ rho * Gamma^2 * L
    # Let's reverse engineer u_hard from Mass Density of matter
    # If matter is mostly vacuum energy, then u_vac > u_matter
    u_hard = rho_hard * (3e8)**2 # mc^2 energy density of nuclear matter

    # Cycle Rate: Limited by Acoustic Resonance of the defect
    # We found this was ~10^43 Hz (Planck), but let's assume we throttle it
    # to "Nuclear" timescales (Gamma rays ~ 10^20 Hz)
    f_hard = 1.0e15 # Femtosecond pulses (Petahertz)

    power_hard = u_hard * efficiency * f_hard

    print(f"\nSCENARIO 2: PLANCK/NUCLEAR SCALE VACUUM")
    print(f"Vacuum Energy Density: {u_hard:.2e} J/m^3")
    print(f"Cycle Frequency:       {f_hard:.1e} Hz")
    print(f"Power Density:         {power_hard:.2e} W/m^3")

    # COMPARISON
    print("-" * 40)
    print(f"PWR (Nuclear Fission): ~1.0e8 W/m^3")
    print(f"The Sun (Core):        ~2.7e2 W/m^3")

    if power_hard > 1e15:
        print("\nVERDICT: The 'Hard' Reactor is effectively a Star in a bottle.")
        print("Safety Warning: Requires extreme throttle (Viscosity < 1e-20).")

calculate_power_density()
