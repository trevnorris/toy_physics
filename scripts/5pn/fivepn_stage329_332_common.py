from __future__ import annotations

import sympy as sp

from fivepn_stage323_325_common import banner, subbanner, expect_zero, split_epsilon, log_drift


def coherent_tracking_state(chi0, deltaU, ZW, epsW, epsEta, Lambda, zeta):
    eps = split_epsilon(epsW, deltaU)
    Rtr = sp.simplify((1 + chi0 / (1 + deltaU)) / (1 + chi0))
    delta = sp.symbols('delta', positive=True, real=True)
    Mmix = sp.simplify(8 * ZW * (1 + chi0) ** 2 / (sp.pi**2 * (1 - epsEta) * (1 - eps)))
    Msupp = sp.simplify(8 * zeta * ZW * (1 + chi0) ** 2 / (sp.pi**2 * (1 - epsEta) * (1 - zeta * eps)))
    S = sp.simplify(1 + zeta * (1 - eps) / (1 - zeta * eps))
    Mtr = sp.simplify(Mmix * S)
    Rtarget = sp.simplify(Lambda * (1 - epsEta) * (1 - eps) ** 2 / (ZW * (1 + chi0) ** 2))
    Cmix = sp.simplify(Rtarget * Mmix)
    Ctr = sp.simplify(8 * Lambda * (1 - eps) / sp.pi**2)
    return {
        'eps': eps,
        'Rtr': Rtr,
        'Mmix': Mmix,
        'Msupp': Msupp,
        'S': S,
        'Mtr': Mtr,
        'Rtarget': Rtarget,
        'Cmix': Cmix,
        'Ctr': Ctr,
    }


def direct_branch_packet(Rtr, Rtarget, epsEta, Rtr_ref, Rtarget_ref, epsEta_ref, chi0_star, deltaU_star):
    Cstar = sp.simplify((1 + chi0_star) * (1 + deltaU_star) * (1 + chi0_star + deltaU_star) / (chi0_star * deltaU_star))
    Bstar = sp.simplify(2 * (1 + chi0_star + deltaU_star) / deltaU_star)
    qtr = sp.simplify(-Cstar * sp.log(Rtr / Rtr_ref))
    qnt = sp.simplify(Bstar * sp.log(Rtr / Rtr_ref) + sp.log((1 - epsEta) / (1 - epsEta_ref)) - sp.log(Rtarget / Rtarget_ref))
    qeta = sp.simplify(sp.log(epsEta / epsEta_ref))
    return {'Cstar': Cstar, 'Bstar': Bstar, 'qtr': qtr, 'qnt': qnt, 'qeta': qeta}


def direct_branch_defects(Rtr, Rtarget, epsEta, dln_Rtr, dln_Rtarget, dln_epsEta, epsEta_star):
    Theta1 = sp.simplify(dln_Rtr)
    Xi1 = sp.simplify(-dln_Rtarget - epsEta_star * dln_epsEta / (1 - epsEta_star))
    R1 = sp.simplify(dln_Rtarget)
    return {'Theta1': Theta1, 'Xi1': Xi1, 'R1': R1}


def G_tr(xi, delta, R):
    return sp.simplify(9 * xi * (xi + delta) / (9 * delta + (9 + 2 * R**2) * xi))


def F_tr(xi, delta, R):
    num = (9 * delta + (9 + 2 * R**2) * xi) ** 2 * (9 * delta + (9 + 2 * R) * xi) ** 2
    den = 81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + (9 + 2 * R**2) * xi**2) ** 2
    return sp.simplify(num / den)


def Pi_tr(xi, delta, R):
    return sp.simplify(F_tr(xi, delta, R) * G_tr(xi, delta, R))


def S_req_from_branch(Pi, Cmix):
    return sp.simplify(Pi / Cmix)


def zeta_req_from_product(Pi, Cmix, eps):
    # Exact elimination of S_req using S_req = Pi / Cmix.
    return sp.simplify((Pi - Cmix) / (Cmix - eps * (2 * Cmix - Pi)))


def mixed_only_condition(Pi, Cmix):
    return sp.simplify(Pi - Cmix)


def lowest_twin_condition(Pi, Cmix):
    return sp.simplify(Pi - 2 * Cmix)


def support_regime_labels(Pi, Cmix):
    return {
        'mixed_only_boundary': sp.simplify(Cmix),
        'lowest_twin_boundary': sp.simplify(2 * Cmix),
    }


def zeta_phys_general(KWeff, Kphi_n_eff, I_n, I_0):
    return sp.simplify((KWeff / Kphi_n_eff) * (I_n / I_0) ** 2)


def I_uniform_ratio(n):
    return sp.simplify(1 / (2 * n + 1))


def zeta_twin_same_operator(n, x):
    return sp.simplify(1 / ((2 * n + 1) ** 2 * (1 + x * n * (n + 1))))


def S_of_zeta(zeta, eps):
    return sp.simplify(1 + zeta * (1 - eps) / (1 - zeta * eps))


def x_max_higher_harmonic(n, zeta_req):
    return sp.simplify((1 / (((2 * n + 1) ** 2) * zeta_req) - 1) / (n * (n + 1)))


def S_n_max(n, eps):
    return sp.simplify(1 + (1 - eps) / ((2 * n + 1) ** 2 - eps))
