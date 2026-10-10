"""Stefan-Boltzmann constant vs the number g of thermally populated radiating modes per wavevector.

sigma(g, c_x) = g/2 * 2 pi^5 k_B^4 / (15 h^3 c^2) for modes at the light speed c.
A mode at a different speed c_x contributes emitted power per area ~ (1/2)(c/c_x)^2 of one
photon polarization (energy density ~ 1/c_x^3, flux ~ c_x * u / 4).
Measured: Quinn & Martin (NPL, 1985), as quoted in Gusev arXiv:1612.03199 eq. (11)-(12):
    sigma_exp = 5.66967(76)e-8 W m^-2 K^-4.
Prints operands, residuals and residual/uncertainty; states no conclusion.
"""
from mpmath import mp, mpf, pi
mp.dps = 30
h = mpf('6.62607015e-34'); kB = mpf('1.380649e-23'); c = mpf('299792458')   # SI exact (2019)
sigma_per_pol = 2 * pi**5 * kB**4 / (15 * h**3 * c**2) / 2
sig_exp = mpf('5.66967e-8'); u_exp = mpf('0.00076e-8')
print('sigma per polarization            =', mp.nstr(sigma_per_pol, 10))
cases = [('g=2 (two polarizations, no extra mode)', 2)]
cases.append(('third matter-coupled mode at c_x/c=1', 2 + 1))
for r in ['0.5', '2', '10', '100']:
    cases.append((f'third matter-coupled mode at c_x/c={r}', 2 + (1 / mpf(r))**2))
for label, g_eff in cases:
    sig_th = g_eff * sigma_per_pol
    resid = sig_exp - sig_th
    print(f'{label:28s} sigma_th = {mp.nstr(sig_th, 10)}  resid(exp-th) = {mp.nstr(resid, 6)}'
          f'  resid/u_exp = {mp.nstr(resid / u_exp, 6)}  relative = {mp.nstr(resid / sig_th, 6)}')
