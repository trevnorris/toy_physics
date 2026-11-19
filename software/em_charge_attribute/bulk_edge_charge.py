"""Translate a finite charge packet through a fixed cut in the Gaussian Maxwell benchmark."""
import sympy as sp


w, b, v, y = sp.symbols("w b v y", real=True)
width, scale, lam, k = sp.symbols("a sigma lambda k", positive=True)
charge, eps = sp.symbols("Q epsilon", real=True, nonzero=True)

# y is the offset from the packet centre in units of its half-width.
# The polynomial and its derivative join zero at either support endpoint.
shape = (1 - y**2)**2
shape /= sp.integrate(shape, (y, -1, 1))
h = shape.subs(y, (w-b)/width)/width
normal = sp.integrate(shape, (y, -1, 1))
norm_sq = sp.integrate(shape**2, (y, -1, 1))/width
z_int = sp.integrate(sp.exp(-(w/lam)**2), (w, -sp.oo, sp.oo))
projection = charge*normal/sp.sqrt(z_int)
residual = sp.factor(v*sp.diff(h,b) + sp.diff(v*h,w))
# The time-integrated flux at a fixed cut for a complete positive crossing.
# Changing from centre position to y reverses the integration limits.
crossing = sp.integrate(charge*shape, (y,-1,1))
print("depth normalization, squared norm, endpoint densities =", normal, norm_sq,
      [h.subs(w,b+s*width) for s in (-1,1)])
print("normalized massless-mode source =", projection)
print("continuity residual within support, integrated cut current =", residual, crossing)

# For a normalized 3D Gaussian with variance sigma^2 per coordinate, compute
# the integrated field energy for unit permittivity from its Fourier transform.
fourier_density = sp.exp(-scale**2*k**2/2)
phase_space = 4*sp.pi*k**2/(2*sp.pi)**3
coulomb_integral = sp.integrate(phase_space*fourier_density**2/k**2, (k,0,sp.oo))
energy_scale = charge**2*coulomb_integral/(2*eps)
far_depth = sp.Abs(b)+width
lower = sp.simplify(energy_scale/z_int)
upper = sp.simplify(energy_scale*norm_sq*sp.exp((far_depth/lam)**2))
print("Coulomb integral, energy lower bound, energy upper bound =",
      coulomb_integral, lower, upper)
assert normal == 1 and residual == 0 and crossing == charge
assert all(h.subs(w,b+s*width) == 0 for s in (-1,1))
