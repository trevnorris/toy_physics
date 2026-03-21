
#!/usr/bin/env python3
import math
import numpy as np
import sympy as sp

def section(title: str) -> None:
    print("\n" + "="*78)
    print(title)
    print("="*78)

# ---------------------------------------------------------------------
# Exact static support/source data and full 2PN cross-block reconstruction
# ---------------------------------------------------------------------
section("EXACT Family-1 static support/source law")

mu = sp.symbols('mu', real=True)
P2 = sp.Rational(1,2)*(3*mu**2 - 1)

z_base = sp.Rational(17,56) - sp.Rational(5,56)*mu**2
t_prof = sp.Rational(593,672) - sp.Rational(1553,672)*mu**2 + sp.Rational(7,8)*mu**4
S_prof = sp.simplify(10 - sp.Rational(63,2)*z_base)

print("z_base(mu) =", z_base)
print("t(mu)      =", t_prof)
print("S(mu)      =", S_prof)
print("Check S - (10 - 63/2 z_base) =", sp.simplify(S_prof - (10 - sp.Rational(63,2)*z_base)))

section("EXACT 2PN cross-block reconstruction from final throat module")

ax, ay, az, bx, by, bz, UA, UB = sp.symbols('ax ay az bx by bz UA UB', real=True)
uAB = ax*bx + ay*by
dA = az
dB = bz
vA2 = ax**2 + ay**2 + az**2
vB2 = bx**2 + by**2 + bz**2
vAB = uAB + dA*dB

# Odd carried-forward dipole wake + 2PN dressing
L1_wake = -sp.Rational(7,2)*uAB - 4*dA*dB
Lodd = sp.Rational(1,2)*(vA2 + vB2)*L1_wake - sp.Rational(15,4)*(UA + UB)*(uAB + dA*dB)

# Even P0 \oplus P2 layer
Pi0A = sp.sqrt(5)/2 * vA2
Pi0B = sp.sqrt(5)/2 * vB2
Pi20A = sp.Rational(1,2)*(3*dA**2 - vA2)
Pi20B = sp.Rational(1,2)*(3*dB**2 - vB2)
Pi21A = (sp.sqrt(2)*dA*ax, sp.sqrt(2)*dA*ay)
Pi21B = (sp.sqrt(2)*dB*bx, sp.sqrt(2)*dB*by)
Pi22A = (sp.Rational(1,1)/(2*sp.sqrt(2))*(ax**2-ay**2), sp.Rational(1,1)/(sp.sqrt(2))*ax*ay)
Pi22B = (sp.Rational(1,1)/(2*sp.sqrt(2))*(bx**2-by**2), sp.Rational(1,1)/(sp.sqrt(2))*bx*by)

dot21 = sp.simplify(Pi21A[0]*Pi21B[0] + Pi21A[1]*Pi21B[1])
dot22 = sp.simplify(Pi22A[0]*Pi22B[0] + Pi22A[1]*Pi22B[1])

J0 = 4/sp.sqrt(5)
J20 = sp.Rational(5,4)
Delta_geom = sp.Rational(281,80)

Leven = sp.simplify(
    Pi0A*Pi0B + Pi20A*Pi20B + dot21 + dot22
    + UA*(J0*Pi0B + J20*Pi20B)
    + UB*(J0*Pi0A + J20*Pi20A)
    + (J0**2 + J20**2 - Delta_geom)*UA*UB
)

Lmodule = sp.expand(Lodd + Leven)

Ltarget = sp.expand(
    -sp.Rational(7,4)*vAB*(vA2+vB2)
    -sp.Rational(1,4)*dA*dB*(vA2+vB2)
    +sp.Rational(11,8)*vA2*vB2
    +sp.Rational(1,4)*vAB**2
    -sp.Rational(5,8)*(vA2*dB**2 + vB2*dA**2)
    +sp.Rational(3,2)*vAB*dA*dB
    +sp.Rational(3,8)*dA**2*dB**2
    -sp.Rational(15,4)*(UA+UB)*vAB
    + UA*(sp.Rational(11,8)*vB2 + sp.Rational(15,8)*dB**2)
    + UB*(sp.Rational(11,8)*vA2 + sp.Rational(15,8)*dA**2)
    + sp.Rational(5,4)*UA*UB
)

residual = sp.expand(Lmodule - Ltarget)
print("Cross-block residual == 0 ?", residual == 0)
print("Residual =", residual)

section("DIRECT coupled sidewall + endcap full-profile wall completion")

# Balanced, thin-layer-consistent reference branch
epsr = 0.05
alphar = 10.0
pr = 2
epsz = 0.05          # d_z / L
chi_cap = 4.0
alphacap = 4.0 * epsz * chi_cap   # = 0.8
pz = 2

def smooth_step(x):
    return 0.5*(1.0 + np.tanh(x))

def gauss_interval(a, b, n):
    x, w = np.polynomial.legendre.leggauss(n)
    x = 0.5*(b-a)*x + 0.5*(a+b)
    w = 0.5*(b-a)*w
    return x, w

# composite Gauss-Legendre on [0,1] with extra boundary resolution
x1, w1 = gauss_interval(0.0, 0.93, 120)
x2, w2 = gauss_interval(0.93, 1.0, 260)
x = np.concatenate([x1, x2])
wx = np.concatenate([w1, w2])

y1, wy1 = gauss_interval(0.0, 0.93, 120)
y2, wy2 = gauss_interval(0.93, 1.0, 260)
y = np.concatenate([y1, y2])
wy = np.concatenate([wy1, wy2])

X, Y = np.meshgrid(x, y, indexing='xy')
WX, WY = np.meshgrid(wx, wy, indexing='xy')

rho_full_arg = (
    1.0
    - Y**2
    - alphar * smooth_step((X-1.0)/epsr)**pr
    - alphacap * smooth_step((Y-1.0)/(2.0*epsz))**pz
)
rho_full = np.maximum(rho_full_arg, 0.0)**0.25
weight = WX*WY

I2 = 2.0*np.sum(weight * X**2 * rho_full)
I4 = 2.0*np.sum(weight * X**4 * rho_full)
IW = 2.0*np.sum(weight * (Y**2/4.0) * X**2 * rho_full)

c0 = math.sqrt(math.pi)*math.gamma(1.25)/(2.0*math.gamma(1.75))
I2sharp = 2.0*c0/3.0
Rmass_full = I2 / I2sharp
Maa_full = I4 / I2
Mll_full = IW / I2

print(f"Rmass_full_exact = {Rmass_full:.15f}")
print(f"Maa_full_exact   = {Maa_full:.15f}")
print(f"Mll_full_exact   = {Mll_full:.15f}")

section("Compare direct full-profile values to separated-order carried-forward branch")

Rmass_sep = 0.8846236634
Maa_sep = 0.5623810783
Mll_sep = 0.0671965962

print(f"Rmass relative shift vs separated  = {(Rmass_full-Rmass_sep)/Rmass_sep:+.8%}")
print(f"Maa   relative shift vs separated  = {(Maa_full-Maa_sep)/Maa_sep:+.8%}")
print(f"Mll   relative shift vs separated  = {(Mll_full-Mll_sep)/Mll_sep:+.8%}")

section("Geometry breathing response from exact full-profile wall completion")

x01 = 2.40482555769577276862163187933
LamEM = math.sqrt(2.0)*math.pi/x01
rho_geom = 0.1
beta_geom = 12.0

H11 = 2.0*(4.0*math.pi*LamEM*(LamEM*rho_geom + LamEM + 2.0) + beta_geom)/LamEM
H12 = 4.0*math.pi*rho_geom + 8.0*math.pi - 2.0*beta_geom/LamEM**2
H22 = 2.0*beta_geom/LamEM**3
Hbar = np.array([[H11, H12], [H12, H22]], dtype=float)
gbar = np.array([3.0, 1.0/LamEM], dtype=float)

delta_unit = gbar @ np.linalg.inv(Hbar) @ gbar
Sigma_star = delta_unit / (109.0/280.0)

Mfull = np.array([[Maa_full, 0.0], [0.0, Mll_full]], dtype=float)
evals, evecs = np.linalg.eig(np.linalg.solve(Mfull, Hbar))
ordr = np.argsort(evals)
evals = evals[ordr]
evecs = evecs[:, ordr]

residues = []
for i in range(2):
    v = evecs[:, i]
    v = v / math.sqrt(v @ Mfull @ v)
    residues.append((gbar @ v)**2 / (Sigma_star * evals[i]))
residues = np.array(residues)

# sharp-wall TF baseline for physical ratio comparison
Mtf = np.array([[3.0/5.0, 0.0], [0.0, 1.0/14.0]], dtype=float)
evals_tf, _ = np.linalg.eig(np.linalg.solve(Mtf, Hbar))
evals_tf = np.sort(evals_tf)
omega_sq_ratios = (evals/evals_tf) / Rmass_full

lam_eff = residues.sum() / np.sum(residues/evals)
sgrid = np.linspace(0.0, 0.1*evals[0], 500)
Yexact = residues[0]/(1.0 - sgrid/evals[0]) + residues[1]/(1.0 - sgrid/evals[1])
Ypade = residues.sum()/(1.0 - sgrid/lam_eff)
max_rel_err = np.max(np.abs((Ypade - Yexact)/Yexact))

print(f"Lambda_EM  = {LamEM:.15f}")
print(f"Sigma_*    = {Sigma_star:.15f}")
print(f"lambda_-   = {evals[0]:.15f}")
print(f"lambda_+   = {evals[1]:.15f}")
print(f"R_-        = {residues[0]:.15f}")
print(f"R_+        = {residues[1]:.15f}")
print(f"R_- + R_+  = {residues.sum():.15f}")
print(f"lambda_eff = {lam_eff:.15f}")
print(f"max rel err on [0,0.1 lambda_-] = {max_rel_err:.15e}")
print(f"Omega^2/Omega_TF^2 (minus) = {omega_sq_ratios[0]:.15f}")
print(f"Omega^2/Omega_TF^2 (plus)  = {omega_sq_ratios[1]:.15f}")

K00raw = -757.0/2520.0
K00full0 = K00raw + residues.sum()
print(f"K00_raw(static local) = {K00raw:.15f}")
print(f"K00_full(0)           = {K00full0:.15f}")
print(f"4/45                  = {4.0/45.0:.15f}")

section("FINAL low-frequency Family-1 throat-response object")

print("Odd dipole residues:   R1_perp = 7/2, R10 = 4")
print("Odd dressings:         sigma = 1/2, eta_perp = 15/14, eta_parallel = 15/16")
print("Even local support:    K00_raw = -757/2520, K1_perp = 2/7, K10 = 1/4, K20 = 4/9, K21 = 2/3, K22 = 8/3")
print("Even source vector:    J = (4/sqrt(5), 5/4)")
print("Monopole closure:      K00(s) = -757/2520 + R_-/(1-s/lambda_-) + R_+/(1-s/lambda_+)")
print("Default balanced cap:  alpha_cap = 0.8, epsz = 0.05, chi_cap = 4")
print("Conservative 2PN status: the full added cross block is reproduced exactly at zero frequency;")
print("dynamic non-monopole pole scales remain genuine inner-throat observables for a beyond-2PN / Paper-7 extension.")
