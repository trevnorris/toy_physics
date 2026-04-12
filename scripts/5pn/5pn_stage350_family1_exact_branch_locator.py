#!/usr/bin/env python3
from __future__ import annotations
import mpmath as mp

"""
Stage 350 — exact Family-1 branch locator on the refreshed lambda_EM geometry.

This script numerically locates the explicit Family-1 support/source branch using
the exact refreshed lambda_EM geometry and the exact fixed-point equation

    Pe = Xi * Delta(Pe; kappa, eta)

for the two explicit wall-depth extractions carried in the notes.
"""

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line, flush=True)
    print(title, flush=True)
    print(line, flush=True)

def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line, flush=True)
    print(title, flush=True)
    print(line, flush=True)

mp.mp.dps = 50

x01 = mp.mpf("2.4048255576957727686216318793264546431242449091459671")
Lambda_EM = mp.sqrt(2) * mp.pi / x01
Lambda_ell = 20 * Lambda_EM
eta = Lambda_ell
chi_s = 10 * Lambda_EM
kappa = 4 * chi_s**2 + (mp.mpf(4) / 5) * Lambda_ell**2
alpha = mp.sqrt(kappa)

# exact Robin root in (0, pi/2)
f = lambda y: y * mp.tan(y) - eta
y = mp.findroot(f, (mp.mpf("1.52"), mp.mpf("1.54")))

den = alpha * mp.sinh(alpha) + eta * mp.cosh(alpha)

def K_kernel(x: mp.mpf) -> mp.mpf:
    return (
        mp.cosh(alpha * x)
        + (eta / alpha) * mp.sinh(alpha * x)
        - mp.cosh(alpha * (1 - x))
    ) / den

def Sigma_Pe(Pe: mp.mpf, x: mp.mpf) -> mp.mpf:
    if abs(Pe) < mp.mpf("1e-40"):
        return mp.mpf(1)
    return Pe * mp.e ** (Pe * (x - 1)) / (1 - mp.e ** (-Pe))

def Delta_of_Pe(Pe: mp.mpf) -> mp.mpf:
    return mp.quad(lambda xx: K_kernel(xx) * Sigma_Pe(Pe, xx), [0, 1])

Delta0 = eta * (mp.cosh(alpha) - 1) / (alpha**2 * den)
DeltaInf = (mp.cosh(alpha) + (eta / alpha) * mp.sinh(alpha) - 1) / den

A_K = (kappa + mp.pi**2 / 4) / (kappa + y**2)
zeta_max = A_K * (mp.pi**2 / 4)

def Omega_Pe(Pe: mp.mpf) -> mp.mpf:
    if abs(Pe) < mp.mpf("1e-40"):
        return mp.mpf(1)
    num = mp.pi * Pe * (2 * Pe + mp.pi * mp.e ** (-Pe))
    deno = (4 * Pe**2 + mp.pi**2) * (1 - mp.e ** (-Pe))
    return num / deno

def zeta_phys(Pe: mp.mpf) -> mp.mpf:
    return Omega_Pe(Pe) ** 2 * A_K

def solve_branch_root(Xi: mp.mpf) -> mp.mpf:
    fun = lambda Pe: Pe - Xi * Delta_of_Pe(Pe)
    # For the explicit Family-1 branches Xi is large and the root sits close to Xi*DeltaInf.
    guess_hi = Xi * DeltaInf
    guess_lo = guess_hi * mp.mpf("0.995")
    return mp.findroot(fun, (guess_lo, guess_hi))

Theta_chi = mp.mpf("4.06863235008162")
Theta_J = mp.mpf("0.927552032539308")

Xi_chi = 100 * Theta_chi * Lambda_ell**2
Xi_J = 100 * Theta_J * Lambda_ell**2

def fmt(x: mp.mpf, digits: int = 16) -> str:
    return mp.nstr(x, n=digits)

banner("STAGE 350 — FAMILY-1 EXACT BRANCH LOCATOR")

subbanner("I. Exact refreshed lambda_EM Family-1 geometry")
print("Lambda_EM =", fmt(Lambda_EM, 25))
print("Lambda_ell = L/ell =", fmt(Lambda_ell, 25))
print("chi_s =", fmt(chi_s, 25))
print("kappa =", fmt(kappa, 25))
print("alpha = sqrt(kappa) =", fmt(alpha, 25))
print("y solving y tan y = eta =", fmt(y, 25))
print("A_K =", fmt(A_K, 25))
print("Delta_0 =", fmt(Delta0, 25))
print("Delta_infty =", fmt(DeltaInf, 25))
print("zeta_max =", fmt(zeta_max, 25))

subbanner("II. Operator-selected fixed-point roots Pe_* from the two explicit wall-depth extractions")
for name, Theta, Xi in [
    ("chi-weighted", Theta_chi, Xi_chi),
    ("J-weighted", Theta_J, Xi_J),
]:
    root = solve_branch_root(Xi)
    Droot = Delta_of_Pe(root)
    Om = Omega_Pe(root)
    zeta = zeta_phys(root)
    rho_alpha_max = 1 + zeta
    print(f"{name} extraction:")
    print("  Theta_w =", fmt(Theta, 20))
    print("  Xi = W_wall =", fmt(Xi, 25))
    print("  Pe_* =", fmt(root, 25))
    print("  Delta(Pe_*) =", fmt(Droot, 25))
    print("  Xi * Delta(Pe_*) =", fmt(Xi * Droot, 25))
    print("  Omega_Pe(Pe_*) =", fmt(Om, 25))
    print("  zeta_phys(Pe_*) =", fmt(zeta, 25))
    print("  rho_alpha,max = 1 + zeta_phys =", fmt(rho_alpha_max, 25))
    print("  saturation gap zeta_max - zeta_phys =", fmt(zeta_max - zeta, 25))
    print()
