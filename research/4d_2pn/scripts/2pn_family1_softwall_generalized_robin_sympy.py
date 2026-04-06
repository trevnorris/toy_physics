
from sympy import *
from sympy import Rational as Q

def line(msg=""):
    print(msg)

# ------------------------------------------------------------
# 1) Support-channel data from the passed 2PN Robin-wall step
# ------------------------------------------------------------
support_data = [
    {"name":"1perp", "l":1, "lam":Q(2), "p2":Q(-1,5), "p22":Q(1,7),  "Z":Q(2,7)},
    {"name":"10",    "l":1, "lam":Q(2), "p2":Q(2,5),  "p22":Q(11,35),"Z":Q(1,4)},
    {"name":"20",    "l":2, "lam":Q(6), "p2":Q(2,7),  "p22":Q(3,7),  "Z":Q(4,9)},
    {"name":"21",    "l":2, "lam":Q(6), "p2":Q(1,7),  "p22":Q(1,7),  "Z":Q(2,3)},
    {"name":"22",    "l":2, "lam":Q(6), "p2":Q(-2,7), "p22":Q(1,7),  "Z":Q(8,3)},
]

A0,A1,B0,B1,C = symbols('A0 A1 B0 B1 C')

eqs = []
for d in support_data:
    Lm2 = d["lam"] - 2
    pred = A0 + A1*Lm2 + (B0 + B1*Lm2)*d["p2"] + C*Lm2*d["p22"]
    eqs.append(Eq(pred, d["Z"]))

sol = solve(eqs, [A0,A1,B0,B1,C], dict=True)[0]

line("=== Minimal generalized Robin / wall-stress support operator ===")
line("Ansatz:")
line("  Z_lm = A0 + A1*(lambda_l-2) + (B0 + B1*(lambda_l-2))*<P2>_lm + C*(lambda_l-2)*<P2^2>_lm")
line("")
line("Solved coefficients:")
for k in [A0,A1,B0,B1,C]:
    line(f"  {k} = {sol[k]}    (~ {N(sol[k], 20)})")
line("")

line("Exact verification on the passed support data:")
for d in support_data:
    Lm2 = d["lam"] - 2
    pred = simplify(sol[A0] + sol[A1]*Lm2 + (sol[B0] + sol[B1]*Lm2)*d["p2"] + sol[C]*Lm2*d["p22"])
    line(f"  {d['name']:5s}: predicted Z = {pred}, target Z = {d['Z']}, residual = {simplify(pred - d['Z'])}")
line("")

# ------------------------------------------------------------
# 2) Operator profiles in P2 and in mu = cos(theta)
# ------------------------------------------------------------
mu = symbols('mu', real=True)
P2 = (3*mu**2 - 1)/2

zBaseP2 = simplify(sol[A0] + sol[B0]*P2)
zCurvP2 = simplify(sol[A1] + sol[B1]*P2 + sol[C]*P2**2)

zBaseMu = expand(zBaseP2)
zCurvMu = expand(zCurvP2)

line("Equivalent profile decomposition:")
line("  Z_lm = < zBase(mu) >_lm + (lambda_l-2) < zCurv(mu) >_lm")
line("")
line("Base profile:")
line(f"  zBase(P2) = {zBaseP2}")
line(f"  zBase(mu) = {zBaseMu}")
line("")
line("Curvature profile:")
line(f"  zCurv(P2) = {zCurvP2}")
line(f"  zCurv(mu) = {zCurvMu}")
line("")

# ------------------------------------------------------------
# 3) Rewrite in the explicit Family-1 flare basis {1, -mu^2, mu^4}
# ------------------------------------------------------------
c0, c2 = symbols('c0 c2')
solBaseMuBasis = solve([
    Eq(c0, zBaseMu.coeff(mu, 0)),
    Eq(-c2, zBaseMu.coeff(mu, 2))
], [c0, c2], dict=True)[0]

d0, d2, d4 = symbols('d0 d2 d4')
solCurvMuBasis = solve([
    Eq(d0, zCurvMu.coeff(mu, 0)),
    Eq(-d2, zCurvMu.coeff(mu, 2)),
    Eq(d4, zCurvMu.coeff(mu, 4))
], [d0, d2, d4], dict=True)[0]

line("Family-1 flare / soft-wall mouth basis:")
line("  h1(mu) = -mu^2")
line("  h2(mu) =  mu^4")
line("")
line("Base profile in {1, h1}:")
line(f"  zBase(mu) = {solBaseMuBasis[c0]} + {solBaseMuBasis[c2]} * h1(mu)")
line("")
line("Curvature profile in {1, h1, h2}:")
line(f"  zCurv(mu) = {solCurvMuBasis[d0]} + {solCurvMuBasis[d2]} * h1(mu) + {solCurvMuBasis[d4]} * h2(mu)")
line("")

# ------------------------------------------------------------
# 4) Family-1 flare expansion itself
# ------------------------------------------------------------
q, r = symbols('q r')
h = -q*mu**2 + r*mu**4
a,b,c = symbols('a b c')
solh = solve([
    Eq(expand(a + b*P2 + c*P2**2).coeff(mu,0), expand(h).coeff(mu,0)),
    Eq(expand(a + b*P2 + c*P2**2).coeff(mu,2), expand(h).coeff(mu,2)),
    Eq(expand(a + b*P2 + c*P2**2).coeff(mu,4), expand(h).coeff(mu,4)),
], [a,b,c], dict=True)[0]

line("Family-1 flare expansion near the mouth:")
line("  a(z)/a_m = 1 - q mu^2 + r mu^4 + O(mu^6),   mu = cos(theta)")
line("")
line("Rewritten in {1, P2, P2^2}:")
line(f"  -q mu^2 + r mu^4 = ({solh[a]}) + ({solh[b]}) P2 + ({solh[c]}) P2^2")
line("")

# ------------------------------------------------------------
# 5) Source profile in the same mu^2 basis
# ------------------------------------------------------------
Jiso = Q(11,8)
Jax  = Q(15,8)
sourceP2 = simplify(Jiso + Jax*P2)
sourceMu = expand(sourceP2)

s0, s2 = symbols('s0 s2')
solSource = solve([
    Eq(s0, sourceMu.coeff(mu,0)),
    Eq(s2, sourceMu.coeff(mu,2))
], [s0, s2], dict=True)[0]

line("Axisymmetric source profile from the earlier support/source solve:")
line(f"  S(P2) = {sourceP2}")
line(f"  S(mu) = {sourceMu} = {solSource[s0]} + {solSource[s2]} mu^2")
line("")

# ------------------------------------------------------------
# 6) Suggested generalized Robin operator
# ------------------------------------------------------------
line("Suggested minimal support-sector boundary operator (diagonal matrix elements on l=1,2):")
line("  B_eff = d_n + zBase(mu) + 1/2 * {(-Delta_S - 2), zCurv(mu)}")
line("")
line("with")
line(f"  zBase(mu) = {zBaseMu}")
line(f"  zCurv(mu) = {zCurvMu}")
line("")
line("This exactly reproduces the passed support impedances for the l=1,2 channels.")
