# Unit/arith checks for the polarization survey r1 review (prints operands, then results)
import math
eV_J = 1.602176634e-19      # J per eV (exact SI)
c = 299792458.0             # m/s (exact SI)
eV_kg = eV_J / c**2         # kg per eV/c^2
print("kg per eV/c^2 =", eV_kg)
for label, m_eV in [("PDG adopted (Ryutov 07) <1e-18 eV", 1e-18),
                    ("Bonetti 16 FRB <1.8e-14 eV", 1.8e-14),
                    ("PDG-listed Retino 16 <1.9e-15 eV", 1.9e-15)]:
    print(f"{label}: {m_eV} eV -> {m_eV*eV_kg:.3e} kg")
for label, m_kg in [("Retino upper end 1.4e-49 kg", 1.4e-49), ("Retino lower end 3.4e-51 kg", 3.4e-51)]:
    print(f"{label}: {m_kg} kg -> {m_kg/eV_kg:.3e} eV")
print("ratio Retino tightest / PDG adopted (kg):", 3.4e-51 / (1e-18*eV_kg))
print("PVLAS sigma_dn / dn_QED = 17/2.5 =", 17/2.5)
print("SPTpol 0.10e-4 rad^2 in deg^2 =", 0.10e-4*(180/math.pi)**2)
print("Quinn&Martin rel. unc. 0.00076/5.66967 =", 0.00076/5.66967)
print("Planck normalization ratio g=3/g=2 =", 3/2, "; fractional excess =", 3/2-1)
print("excess / QM rel. unc. =", (3/2-1)/(0.00076/5.66967))
