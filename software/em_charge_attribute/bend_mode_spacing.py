"""Read the reflecting-end separation using the paper's dispersive support branch."""

import sympy as sp


def check(label, a, b):
    residual = sp.factor(a-b)
    print(f"{label}: measured={a}, reference={b}, residual={residual}", flush=True)
    if residual != 0:
        raise ValueError(label)


z = sp.symbols("z", real=True)
ell, speed, gap = [sp.Function(name, positive=True)(z) for name in ("ell", "c", "gap")]
k0 = sp.pi/(2*ell)
omega0_squared = gap**2 + speed**2*k0**2
omega1_squared = gap**2 + speed**2*(3*k0)**2
spacing = sp.factor(omega1_squared - omega0_squared)
print(f"squared-frequency spacing = {spacing}")
check("length change from two modes",
      sp.factor(ell*(speed.diff(z)/speed - spacing.diff(z)/(2*spacing))), ell.diff(z))

# Partial shape load. Also keep the distinct work terms when material
# dispersion varies independently during the bend.
beta = speed**2*k0**2/omega0_squared
log_drift = sp.diff(omega0_squared,z)/(2*omega0_squared)
check("frequency change per bend",
      log_drift, gap*gap.diff(z)/omega0_squared + beta*(speed.diff(z)/speed-ell.diff(z)/ell))
print("end force / energy at fixed material =",
      sp.factor(-sp.diff(sp.log(omega0_squared)/2,ell)))
check("end force / energy at fixed material",
      -sp.diff(sp.log(omega0_squared)/2,ell), beta/ell)
