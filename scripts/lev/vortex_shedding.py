import numpy as np

def vortex_shedding_calculator():
    """
    Calculates the 'Detachment Threshold' for vortex rings in the proton throat.
    Determines the required pulse frequency and gradient strength to force
    vortex shedding into the bulk (Dark Thrust) rather than photon emission.
    """
    print("--- Vortex Shedding & Detachment Threshold Estimator ---")

    # 1. Physical Constants
    m_proton = 1.67e-27  # kg
    c = 3.0e8            # m/s
    a = 0.84e-15         # m (Proton radius)
    h_bar = 1.05e-34     # J*s

    # 2. Throat Geometry (from brane_bulk_ontology.tex)
    # Preferred aspect ratio Lambda_star = L/a ~ 1.85
    Lambda = 1.85
    L_throat = Lambda * a

    # 3. Fluid Parameters
    # Circulation Quantum Gamma = h/m (standard QM) or scaled?
    # From previous turn: Gamma ~ 2*pi*a*c is the MAX.
    # Standard QM circulation: Gamma_qm = h / m_proton
    # Let's check which one fits the 'Unit Charge' in the model.
    # em_fields.tex usually maps charge q to Gamma.
    # Let's use the QM spin circulation as the "flippable" unit.
    Gamma_spin = h_bar / m_proton # This is tiny.
    # Wait, the previous script used Gamma ~ a*c (~1e-6).
    # Gamma_spin ~ 1e-34 / 1e-27 ~ 1e-7.
    # Let's use the 'Unit Circulation' Gamma_0 from the model scaling
    # which makes the coupling dimensionless constants order 1.
    # If v_tip = c, then Gamma ~ a*c.
    Gamma_model = a * c

    print(f"Geometry & Fluid Parameters:")
    print(f"  Throat Radius (a): {a:.2e} m")
    print(f"  Throat Length (L): {L_throat:.2e} m (Lambda={Lambda})")
    print(f"  Model Circulation (Gamma ~ ac): {Gamma_model:.2e} m^2/s")

    # 4. Vortex Ring Self-Velocity (U)
    # Speed at which a generated ring travels down the throat.
    # Formula: U = (Gamma / 4*pi*R) * [ln(8R/core) - 0.5]
    # Assume R = a (fits in throat) and core ~ a/10 (stiff fluid core size?)
    R_ring = a
    r_core = a / 10.0 # Ansatz

    U_ring = (Gamma_model / (4 * np.pi * R_ring)) * (np.log(8 * R_ring / r_core) - 0.5)

    print(f"\nVortex Kinematics:")
    print(f"  Self-Induced Velocity (U): {U_ring:.2e} m/s")
    print(f"  Ratio U/c: {U_ring/c:.4f}")

    # 5. Shedding Frequency (Strouhal Limit)
    # To shed a ring, the drive frequency must match the 'clearing time' of the throat.
    # If you flip too slow, the fluid just stretches.
    # If you flip fast enough (Strouhal Number St ~ 0.2 - 0.4), it pulses packets.
    # f_shed = St * U / L

    St_crit = 0.3 # Typical jet shedding value
    f_shed = St_crit * U_ring / L_throat

    print(f"\nShedding Thresholds:")
    print(f"  Critical Frequency (f_shed): {f_shed:.2e} Hz")

    # Compare to known bands
    bands = {
        "NMR": 1e8,
        "EPR": 1e11,
        "Optical": 1e15,
        "X-Ray": 1e18,
        "Gamma": 1e22
    }

    print("  Comparison to Technology:")
    for name, f in bands.items():
        ratio = f / f_shed
        status = "Too Slow (Heat)" if ratio < 0.1 else "Too Fast (Blur)" if ratio > 10 else "RESONANT"
        print(f"    {name:<10} ({f:.0e} Hz): {ratio:.2e} x f_shed  -> {status}")

    # 6. The Gradient Trigger (Force Requirement)
    # If f_shed is too high (likely Gamma/X-ray), can we force it with a Gradient?
    # Force F = d(Momentum)/dt ~ P_ring / t_transit
    # t_transit = L / U
    # P_ring ~ rho * Gamma * pi * a^2 (from previous turn)
    # We need a Magnetic Field Gradient (dB/dz) that exerts this force on the magnetic moment mu.
    # F_mag = mu * (dB/dz)

    # Proton Magnetic Moment
    mu_p = 1.41e-26 # J/T

    # Estimate Fluid Force required to eject ring
    # Assume rho_0 ~ nuclear density (~1e17 kg/m^3)
    rho_0 = 1e17
    P_ring = rho_0 * Gamma_model * np.pi * a**2
    t_transit = L_throat / U_ring
    F_required = P_ring / t_transit

    # Required Gradient
    dB_dz_required = F_required / mu_p

    print(f"\nGradient Trigger Requirement:")
    print(f"  Momentum to Eject (P_ring): {P_ring:.2e} kg m/s")
    print(f"  Transit Time (t_transit): {t_transit:.2e} s")
    print(f"  Force Required (F): {F_required:.2e} N")
    print(f"  Required Magnetic Gradient (dB/dz): {dB_dz_required:.2e} T/m")

    # Is this gradient achievable?
    # Standard MRI gradient ~ 0.1 T/m
    # Lab High Gradient ~ 100 T/m
    # Laser Wakefield ~ 10^6 T/m?

    print("-" * 50)
    print("CONCLUSION:")
    if dB_dz_required > 1e9:
        print("  Standard Coils: IMPOSSIBLE.")
        print("  Requires: Laser Wakefield or Nuclear-Scale Field Gradients.")
    elif dB_dz_required > 100:
        print("  Standard Coils: HARD.")
        print("  Requires: Micro-coils or Nanostructured Magnetic Tips.")
    else:
        print("  Standard Coils: FEASIBLE.")

if __name__ == "__main__":
    vortex_shedding_calculator()
