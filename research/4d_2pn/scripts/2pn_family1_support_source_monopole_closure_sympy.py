
import sympy as sp


def section(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


mu = sp.symbols("mu", real=True)
P2 = sp.legendre(2, mu)

# -----------------------------------------------------------------------------
# 1) Carried-forward Family-1 local wall data
# -----------------------------------------------------------------------------
section("1) Carried-forward local Family-1 support/source data")

# Exact flare closure data from the passed support-sector reconstruction
xi = sp.Rational(1176, 1553)
chi = sp.expand(1 - xi * mu**2 + sp.Rational(1, 2) * xi**2 * mu**4)
Delta_chi = sp.expand(sp.diff((1 - mu**2) * sp.diff(chi, mu), mu))

tau0 = -sp.Rational(1714441, 790272)
tau1 = sp.Rational(2411809, 790272)

pi0 = sp.Rational(17671223, 241790304)
pi1 = sp.Rational(60295225, 241790304)
pi2 = sp.Rational(12059045, 967161216)

z_base = sp.expand(pi0 + pi1 * chi + pi2 * Delta_chi)
t_profile = sp.expand(tau0 + tau1 * chi)

z_target = sp.expand(sp.Rational(17, 56) - sp.Rational(5, 56) * mu**2)
t_target = sp.expand(
    sp.Rational(593, 672) - sp.Rational(1553, 672) * mu**2 + sp.Rational(7, 8) * mu**4
)

print(f"chi(mu)       = {chi}")
print(f"Delta chi(mu) = {Delta_chi}")
print("")
print(f"z_base(mu)    = {z_base}")
print(f"t(mu)         = {t_profile}")
print("")
print(f"z residual    = {sp.expand(z_base - z_target)}")
print(f"t residual    = {sp.expand(t_profile - t_target)}")

# -----------------------------------------------------------------------------
# 2) Exact source-from-support closure
# -----------------------------------------------------------------------------
section("2) Exact source-from-support closure")

source_target = sp.expand(sp.Rational(7, 16) + sp.Rational(45, 16) * mu**2)
source_p2 = sp.expand(sp.Rational(11, 8) + sp.Rational(15, 8) * P2)
source_from_z = sp.expand(10 - sp.Rational(63, 2) * z_base)

print(f"S_target(mu)      = {source_target}")
print(f"S_target(mu) [P2] = {source_p2}")
print(f"S_from_z(mu)      = {source_from_z}")
print("")
print(f"source residual   = {sp.expand(source_from_z - source_target)}")
print("")
print("So the wall/source sector is not independent:")
print("  S(mu) = 10 - (63/2) z_base(mu)")

# The same relation in the exact flare basis {1, chi, Delta chi}
s0 = sp.Rational(177262811, 23027648)
s1 = -sp.Rational(180885675, 23027648)
s2 = -sp.Rational(36177135, 92110592)

print("")
print("Exact flare-basis coefficient relation:")
print(f"  s0 - [10 - 63 pi0/2] = {sp.simplify(s0 - (10 - sp.Rational(63, 2) * pi0))}")
print(f"  s1 + 63 pi1/2        = {sp.simplify(s1 + sp.Rational(63, 2) * pi1)}")
print(f"  s2 + 63 pi2/2        = {sp.simplify(s2 + sp.Rational(63, 2) * pi2)}")

# -----------------------------------------------------------------------------
# 3) Dipole support determines the whole scalar source sector
# -----------------------------------------------------------------------------
section("3) Dipole support residues determine the scalar source sector")

K1perp, K10 = sp.symbols("K1perp K10", real=True)

# Any axisymmetric quadratic base profile is reconstructed from the l=1 data.
a0 = sp.simplify((K10 + 2 * K1perp) / 3)
a2 = sp.simplify(sp.Rational(5, 3) * (K10 - K1perp))
z_from_dipole = sp.expand(a0 + a2 * P2)

p_iso = sp.simplify(10 - sp.Rational(63, 2) * a0)
q_ax = sp.simplify(-sp.Rational(63, 2) * a2)

print("Reconstruct the base profile from dipole support data:")
print("  z_base(mu) = a0 + a2 P2(mu)")
print(f"  a0 = {a0}")
print(f"  a2 = {a2}")
print("")
print(f"  z_from_dipole(mu) = {z_from_dipole}")
print("")
print("Then the wall/source coefficients are fixed to")
print(f"  p_iso = {p_iso}")
print(f"  q_ax  = {q_ax}")

# Specialize to the passed dipole support targets
dipole_targets = {K1perp: sp.Rational(2, 7), K10: sp.Rational(1, 4)}

print("")
print("Specializing to the passed dipole support targets:")
print(f"  z_from_dipole - z_target = {sp.expand(z_from_dipole.subs(dipole_targets) - z_target)}")
print(f"  p_iso target residual    = {sp.simplify(p_iso.subs(dipole_targets) - sp.Rational(11, 8))}")
print(f"  q_ax target residual     = {sp.simplify(q_ax.subs(dipole_targets) - sp.Rational(15, 8))}")

# -----------------------------------------------------------------------------
# 4) Canonical scalar source vector from dipole support data
# -----------------------------------------------------------------------------
section("4) Canonical scalar source vector from the dipole support data")

# Mapping from the wall source profile S = p_iso + q_ax P2 to the canonical
# scalar source vector J=(J0,J20) used in the passed minimal irrep throat operator.
J0 = sp.simplify(sp.Rational(2, 1) / sp.sqrt(5) * (p_iso + q_ax / 3))
J20 = sp.simplify(sp.Rational(2, 3) * q_ax)

print("Canonical mapping:")
print("  J0  = (2/sqrt(5)) (p_iso + q_ax/3)")
print("  J20 = (2/3) q_ax")
print("")
print(f"  J0(K1perp,K10)  = {J0}")
print(f"  J20(K1perp,K10) = {J20}")

J_targets = {
    sp.Symbol("J0target"): sp.Rational(4, 1) / sp.sqrt(5),
    sp.Symbol("J20target"): sp.Rational(5, 4),
}

print("")
print("At the passed dipole support values:")
print(f"  J0 residual  = {sp.simplify(J0.subs(dipole_targets) - sp.Rational(4, 1) / sp.sqrt(5))}")
print(f"  J20 residual = {sp.simplify(J20.subs(dipole_targets) - sp.Rational(5, 4))}")

print("")
print("So the scalar source vector is fully determined by the dipole support residues:")
print("  no new source-sector fit is needed once the support base profile is fixed.")

# -----------------------------------------------------------------------------
# 5) Full static local wall operator + monopole auxiliary completion
# -----------------------------------------------------------------------------
section("5) Full static local wall operator + monopole auxiliary completion")

def density_lm(ell: int, m: int) -> sp.Expr:
    P = sp.assoc_legendre(ell, m, mu)
    norm = sp.Rational(2 * ell + 1, 2) * sp.factorial(ell - m) / sp.factorial(ell + m)
    return sp.simplify(norm * P**2)

def K_channel(ell: int, m: int) -> sp.Expr:
    w = density_lm(ell, m)
    return sp.simplify(
        sp.integrate(z_base * w, (mu, -1, 1))
        + (ell * (ell + 1) - 2) * sp.integrate(t_profile * w, (mu, -1, 1))
    )

K00_raw = sp.simplify(K_channel(0, 0))
K11 = sp.simplify(K_channel(1, 1))
K10_exact = sp.simplify(K_channel(1, 0))
K20 = sp.simplify(K_channel(2, 0))
K21 = sp.simplify(K_channel(2, 1))
K22 = sp.simplify(K_channel(2, 2))

DeltaK00 = sp.simplify(sp.Rational(4, 45) - K00_raw)

print("Local Family-1 wall operator")
print("  O_loc = z_base(mu) + 1/2 { -Delta_S - 2, t(mu) }")
print("")
print("Channel values from O_loc:")
print(f"  K00 raw = {K00_raw}")
print(f"  K11     = {K11}")
print(f"  K10     = {K10_exact}")
print(f"  K20     = {K20}")
print(f"  K21     = {K21}")
print(f"  K22     = {K22}")
print("")
print(f"Monopole add-on needed to reach the carried-forward target 4/45: DeltaK00 = {DeltaK00}")
print("")
print("Minimal full static even-wall completion:")
print("  O_full = O_loc + DeltaK00 * P_00")
print("with P_00 the monopole projector.")
print("")
print("Equivalent auxiliary breathing-channel realizations:")
print(f"  residue-normalized pole:   Y_mono(omega) = {DeltaK00}/(1 - omega^2/Omega_mono^2)")
print(f"  stiffness-normalized aux:  E[b,Psi00] = 1/2 b^2 - sqrt({DeltaK00}) b Psi00")

# -----------------------------------------------------------------------------
# 6) Compact closure statement
# -----------------------------------------------------------------------------
section("6) Compact closure statement")

print("Exact static Family-1 even-sector closure:")
print("  1) local support operator")
print("       O_loc = z_base + 1/2{ -Delta_S - 2, t }")
print("  2) source fixed by the base support profile")
print("       S = 10 - (63/2) z_base")
print("  3) one global monopole breathing add-on")
print(f"       DeltaK00 = {DeltaK00}")
print("")
print("So the entire static even wall sector is:")
print("  [local Family-1 support/source constitutive law] + [one monopole auxiliary channel]")
