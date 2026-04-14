#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr) -> None:
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(sp.expand(expr))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


banner("STAGE 019 — RIGID-MOUTH DEPENDENT-PLANE PROJECTORS AND THE EQUAL-DRIFT DRESSING RAY")

subbanner("I. Exact rigid-mouth packet map and dependent microscopic section")

eps = sp.symbols("epsilon_eta_star", positive=True, real=True)
ceta = sp.simplify(eps / (1 - eps))

R1, E1 = sp.symbols("R1 E1", real=True)
xrm = sp.Matrix([R1, E1])
Mrm = sp.Matrix([[-1, -ceta], [0, 1]])
qnt, qeta = sp.symbols("q_nt q_eta", real=True)
qrm = sp.Matrix([qnt, qeta])
Sdep_rm = sp.Matrix([[0, 0], [0, -1], [1, -1]])

print("M_rm =")
sp.pprint(Mrm)
print("S_rm^(dep) =")
sp.pprint(Sdep_rm)

expect_zero("M_rm^2 - I_2", Mrm**2 - sp.eye(2))

yrm_from_q = sp.simplify(Sdep_rm * qrm)
print("Dependent correction S_rm^(dep) q_rm =")
sp.pprint(yrm_from_q)
expect_zero(
    "S_rm^(dep) q_rm - (0,-q_eta,q_nt-q_eta)",
    yrm_from_q - sp.Matrix([0, -qeta, qnt - qeta]),
)

Crm_dep = sp.simplify(Sdep_rm * Mrm)
print("C_rm^(dep) = S_rm^(dep) M_rm =")
sp.pprint(Crm_dep)
expect_zero(
    "C_rm^(dep) x_rm - (0,-E1,-R1-E1/(1-eps))",
    sp.simplify(Crm_dep * xrm - sp.Matrix([0, -E1, -R1 - E1 / (1 - eps)])),
)

subbanner("II. Exact left inverse and packet recovery on the dependent plane")

Ldep_rm = sp.Matrix([[0, -1, 1], [0, -1, 0]])
print("L_rm^(dep) =")
sp.pprint(Ldep_rm)
expect_zero("L_rm^(dep) S_rm^(dep) - I_2", sp.simplify(Ldep_rm * Sdep_rm - sp.eye(2)))

T, DKeta, Dmu = sp.symbols("Delta_T Delta_Keta Delta_mu", real=True)
ydep = sp.Matrix([T, DKeta, Dmu])
q_from_y = sp.simplify(Ldep_rm * ydep)
print("Packet recovery q_rm = L_rm^(dep) y =")
sp.pprint(q_from_y)
expect_zero(
    "q_rm recovery - (Delta_mu-Delta_Keta,-Delta_Keta)",
    q_from_y - sp.Matrix([Dmu - DKeta, -DKeta]),
)

subbanner("III. Exact dependent-plane packet projectors")

Qnt = sp.Matrix([[1, 0], [0, 0]])
Qeta = sp.Matrix([[0, 0], [0, 1]])
Pnt_dep = sp.simplify(Sdep_rm * Qnt * Ldep_rm)
Peta_dep = sp.simplify(Sdep_rm * Qeta * Ldep_rm)
Prm_dep = sp.Matrix([[0, 0, 0], [0, 1, 0], [0, 0, 1]])

print("P_nt^(dep) =")
sp.pprint(Pnt_dep)
print("P_eta^(dep) =")
sp.pprint(Peta_dep)

expect_zero("(P_nt^(dep))^2 - P_nt^(dep)", sp.simplify(Pnt_dep * Pnt_dep - Pnt_dep))
expect_zero("(P_eta^(dep))^2 - P_eta^(dep)", sp.simplify(Peta_dep * Peta_dep - Peta_dep))
expect_zero("P_nt^(dep) P_eta^(dep)", sp.simplify(Pnt_dep * Peta_dep))
expect_zero("P_eta^(dep) P_nt^(dep)", sp.simplify(Peta_dep * Pnt_dep))
expect_zero("P_nt^(dep)+P_eta^(dep)-P_rm^(dep)", sp.simplify(Pnt_dep + Peta_dep - Prm_dep))

# Restrict to the rigid-mouth dependent plane Delta_T = 0.
yrm_general = sp.Matrix([0, DKeta, Dmu])
ynt = sp.simplify(Pnt_dep * yrm_general)
yeta = sp.simplify(Peta_dep * yrm_general)
print("y_nt = P_nt^(dep) y_rm =")
sp.pprint(ynt)
print("y_eta = P_eta^(dep) y_rm =")
sp.pprint(yeta)
expect_zero("y_nt - (0,0,Delta_mu-Delta_Keta)", ynt - sp.Matrix([0, 0, Dmu - DKeta]))
expect_zero("y_eta - (0,Delta_Keta,Delta_Keta)", yeta - sp.Matrix([0, DKeta, DKeta]))
expect_zero("y_nt + y_eta - y_rm", ynt + yeta - yrm_general)

subbanner("IV. Static-blind line and equal-drift dressing ray")

Xi1 = sp.symbols("Xi1", real=True)
qrm_from_x = sp.simplify(Mrm * xrm)
print("q_rm = M_rm x_rm =")
sp.pprint(qrm_from_x)
expect_zero("q_nt - (-R1-ceta E1)", sp.simplify(qrm_from_x[0] - (-R1 - ceta * E1)))
expect_zero("q_eta - E1", sp.simplify(qrm_from_x[1] - E1))

# Apply projectors to the direct-image dependent correction.
yrm_from_x = sp.simplify(Crm_dep * xrm)
ynt_from_x = sp.simplify(Pnt_dep * yrm_from_x)
yeta_from_x = sp.simplify(Peta_dep * yrm_from_x)
print("y_rm(x_rm) =")
sp.pprint(yrm_from_x)
print("y_nt(x_rm) =")
sp.pprint(ynt_from_x)
print("y_eta(x_rm) =")
sp.pprint(yeta_from_x)
expect_zero("y_nt(x_rm) - (0,0,Xi1)", ynt_from_x - sp.Matrix([0, 0, -R1 - ceta * E1]))
expect_zero("y_eta(x_rm) - (0,-E1,-E1)", yeta_from_x - sp.Matrix([0, -E1, -E1]))

# Static strip Xi1 = 0 <=> R1 = -ceta E1
ystatic_strip = sp.simplify(yrm_from_x.subs(R1, -ceta * E1))
print("y_rm on the static strip Xi1=0 =")
sp.pprint(ystatic_strip)
expect_zero("static-strip dependent correction - (0,-E1,-E1)", ystatic_strip - sp.Matrix([0, -E1, -E1]))
static_norm_sq = sp.simplify((ystatic_strip.T * ystatic_strip)[0])
expect_zero("||y_eta||^2 - 2 E1^2", static_norm_sq - 2 * E1**2)

subbanner("V. Exact microscopic correction compilers")

# Static-only and full orbit corrections on the dependent plane
Dy_static = sp.simplify(-ynt_from_x)
Dy_orbit = sp.simplify(-yrm_from_x)
print("Delta y_static =")
sp.pprint(Dy_static)
print("Delta y_orbit =")
sp.pprint(Dy_orbit)
expect_zero("y_rm + Delta y_static - y_eta", sp.simplify(yrm_from_x + Dy_static - yeta_from_x))
expect_zero(
    "Delta y_orbit - (Delta y_static + q_eta (0,1,1))",
    sp.simplify(Dy_orbit - (Dy_static + qrm_from_x[1] * sp.Matrix([0, 1, 1]))),
)

banner("STAGE 019 LEDGER")
print("1. On the rigid-mouth slice q_tr=0, the exact microscopic dependent correction is")
print("      (Delta_T, Delta_Keta, Delta_mu) = (0, -q_eta, q_nt - q_eta).")
print("2. So the rigid-mouth quotient-failure image is the full plane Delta_T = 0.")
print("3. The exact packet recovery on that plane is")
print("      q_nt = Delta_mu - Delta_Keta,    q_eta = - Delta_Keta.")
print("4. The dependent-plane packet projectors split the correction into")
print("      y_nt  = (0, 0, q_nt),")
print("      y_eta = -q_eta (0,1,1).")
print("5. Therefore the Stage-018 static strip q_nt = 0 is exactly the diagonal line")
print("      Delta_mu = Delta_Keta")
print("   in the dependent microscopic plane.")
print("6. So after the first static ceiling is cleared, the remaining rigid-mouth")
print("   orbit defect is the equal-drift K_eta–mu dressing ray at fixed T_U.")
