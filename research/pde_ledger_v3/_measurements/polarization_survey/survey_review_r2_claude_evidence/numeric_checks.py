import math
eV_kg = 1.78266192e-36
print("PDG 1e-18 eV/c^2 in kg:", 1e-18*eV_kg)
print("SPTpol 0.10e-4 rad^2 in deg^2:", 0.10e-4*(180/math.pi)**2)
print("Quinn-Martin relative unc:", 0.00076/5.66967)
print("Quinn-Martin vs CODATA 5.670374419e-8 rel diff:", (5.670374419-5.66967)/5.670374419)
print("PVLAS sigma/QED:", 17/2.5)
print("Ran+ 3.5e-51 kg in eV:", 3.5e-51/eV_kg, " 6.5e-51:", 6.5e-51/eV_kg)
print("Lemos+ 29.4e-51 kg in eV:", 29.4e-51/eV_kg)
