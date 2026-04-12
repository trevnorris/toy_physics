from __future__ import annotations

import math
from dataclasses import dataclass

import mpmath as mp
import numpy as np
import sympy as sp

mp.mp.dps = 80


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


PI = mp.pi
LAMBDA_STAR = mp.mpf(37) / 20
rF1 = mp.sqrt((12 / PI**2) * LAMBDA_STAR**2 - 1)
g_minus_F1 = rF1 - mp.mpf("0.5") * mp.sqrt(1 + rF1**2)
g_plus_F1 = rF1 + mp.mpf("0.5") * mp.sqrt(1 + rF1**2)


def g_pi(Pi: mp.mpf) -> mp.mpf:
    Pi = mp.mpf(Pi)
    return 2 * Pi * (2 * Pi * mp.e**Pi + PI) / ((4 * Pi**2 + PI**2) * (mp.e**Pi - 1))



def S_q(Pi: mp.mpf) -> mp.mpf:
    Pi = mp.mpf(Pi)
    return Pi * ((PI / 2) * mp.tanh(PI / 2) + Pi * (mp.e ** (-Pi) * mp.sech(PI / 2) - 1)) / ((1 - mp.e ** (-Pi)) * ((PI**2) / 4 - Pi**2))


Pi_star = mp.findroot(lambda z: g_pi(z) - g_minus_F1, mp.mpf("1.5"))
S_star = S_q(Pi_star)
Sigma_m_star = Pi_star / (4 - S_star)
Sigma0_star = Pi_star / (1 - S_star / 4)
Tm_star = mp.sqrt(9 * Sigma0_star / 20)


def R_q(Pi: mp.mpf) -> mp.mpf:
    Pi = mp.mpf(Pi)
    return (g_pi(Pi) - rF1) ** 2 / (1 + rF1**2)



def Sigma0_of_Pi(Pi: mp.mpf) -> mp.mpf:
    Pi = mp.mpf(Pi)
    return Pi / (1 - R_q(Pi) * S_q(Pi))



def T_hat_of_Pi(Pi: mp.mpf) -> mp.mpf:
    Pi = mp.mpf(Pi)
    return mp.sqrt(9 * Sigma0_of_Pi(Pi) / 20)


Pi_match = mp.findroot(lambda z: g_pi(z) - PI / 4, mp.mpf("2.0"))

gprime_star = mp.diff(g_pi, Pi_star)
Sprime_star = mp.diff(S_q, Pi_star)
A_T = -(mp.mpf(9) / (40 * Tm_star)) * (
    1 / (gprime_star * (1 - S_star / 4))
    + (Pi_star * Sprime_star) / (4 * gprime_star * (1 - S_star / 4) ** 2)
)
B_T = (mp.mpf(9) / (40 * Tm_star)) * (Pi_star / (4 * (1 - S_star / 4) ** 2))



def c_kernel(x: mp.mpf) -> mp.mpf:
    x = mp.mpf(x)
    return mp.cos(PI * x / 2)



def Kq_kernel(x: mp.mpf) -> mp.mpf:
    x = mp.mpf(x)
    return mp.cosh(PI * (1 - x) / 2) / mp.cosh(PI / 2)



def Sigma_star(x: mp.mpf) -> mp.mpf:
    x = mp.mpf(x)
    return Pi_star * mp.e ** (-Pi_star * x) / (1 - mp.e ** (-Pi_star))


C_q_star = Pi_star / ((1 - mp.e ** (-Pi_star)) * ((PI**2) / 4 - Pi_star**2))
A_q_star = C_q_star * (((PI / 2) * mp.sinh(PI / 2)) + Pi_star * mp.e ** (-Pi_star)) / (((PI / 2) * mp.cosh(PI / 2)))



def T_s_star(x: mp.mpf) -> mp.mpf:
    x = mp.mpf(x)
    return (1 - mp.e ** (-Pi_star * x)) / (Pi_star * (1 - mp.e ** (-Pi_star))) - x * mp.e ** (-Pi_star) / (1 - mp.e ** (-Pi_star))



def T_q_star(x: mp.mpf) -> mp.mpf:
    x = mp.mpf(x)
    return A_q_star * mp.sinh(PI * x / 2) - C_q_star * mp.cosh(PI * x / 2) + C_q_star * mp.e ** (-Pi_star * x)



def R_star(x: mp.mpf) -> mp.mpf:
    x = mp.mpf(x)
    return Sigma_m_star * (4 * T_s_star(x) - T_q_star(x) - (4 - S_star) * x)



def weighted_average_sigma_star(func) -> mp.mpf:
    return mp.quad(lambda xx: Sigma_star(xx) * func(xx), [0, 1])



def covariance_sigma_star(func_f, func_h) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    mean_f = weighted_average_sigma_star(func_f)
    mean_h = weighted_average_sigma_star(func_h)
    mean_fh = mp.quad(lambda xx: Sigma_star(xx) * func_f(xx) * func_h(xx), [0, 1])
    return mean_fh - mean_f * mean_h, mean_f, mean_h


cov_c_R_star, _, mean_R_star = covariance_sigma_star(c_kernel, R_star)
cov_Kq_R_star, _, _ = covariance_sigma_star(Kq_kernel, R_star)
delta_g_act = -cov_c_R_star
delta_S_act = -cov_Kq_R_star
delta_Pi_act = -delta_g_act / gprime_star
delta_T_act = A_T * delta_g_act + B_T * delta_S_act
Pi_corr = Pi_star + delta_Pi_act
T_corr = Tm_star + delta_T_act


def Sigma1_unnormalized(x: mp.mpf) -> mp.mpf:
    x = mp.mpf(x)
    return mp.e ** (-Pi_star * x - R_star(x))


Z1 = mp.quad(lambda xx: Sigma1_unnormalized(xx), [0, 1])


def Sigma1(x: mp.mpf) -> mp.mpf:
    x = mp.mpf(x)
    return Sigma1_unnormalized(x) / Z1


g1 = mp.quad(lambda xx: Sigma1(xx) * c_kernel(xx), [0, 1])
S1 = mp.quad(lambda xx: Sigma1(xx) * Kq_kernel(xx), [0, 1])
Pi1 = Pi_star - (g1 - g_minus_F1) / gprime_star
T1 = Tm_star + A_T * (g1 - g_minus_F1) + B_T * (S1 - S_star)


def uniform_profile(x: mp.mpf) -> mp.mpf:
    return mp.mpf(1)



def derivative_profile(x: mp.mpf) -> mp.mpf:
    x = mp.mpf(x)
    return (PI / 2) * mp.cos(PI * x / 2)


g_u = mp.quad(lambda xx: uniform_profile(xx) * c_kernel(xx), [0, 1])
S_u = mp.quad(lambda xx: uniform_profile(xx) * Kq_kernel(xx), [0, 1])
g_d = mp.quad(lambda xx: derivative_profile(xx) * c_kernel(xx), [0, 1])
S_d = mp.quad(lambda xx: derivative_profile(xx) * Kq_kernel(xx), [0, 1])

delta_Pi_u = -(g_u - g_minus_F1) / gprime_star
delta_T_u = A_T * (g_u - g_minus_F1) + B_T * (S_u - S_star)
delta_Pi_d = -(g_d - g_minus_F1) / gprime_star
delta_T_d = A_T * (g_d - g_minus_F1) + B_T * (S_d - S_star)

lambda_Pi_zero = delta_Pi_u / (delta_Pi_u - delta_Pi_d)
lambda_T_zero = delta_T_u / (delta_T_u - delta_T_d)
lambda_eff_Pi = (delta_Pi_u - delta_Pi_act) / (delta_Pi_u - delta_Pi_d)
lambda_eff_T = (delta_T_u - delta_T_act) / (delta_T_u - delta_T_d)


class CoevolvingFamily1Solver:
    def __init__(self, N: int = 500):
        self.N = N
        self.x = np.linspace(0.0, 1.0, N)
        self.dx = self.x[1] - self.x[0]
        self.w = np.ones(N) * self.dx
        self.w[0] = self.w[-1] = self.dx / 2
        self.minxy = np.minimum.outer(self.x, self.x)
        kappa = math.pi / 2
        minm = np.minimum.outer(self.x, self.x)
        maxm = np.maximum.outer(self.x, self.x)
        self.Gq = np.sinh(kappa * minm) * np.cosh(kappa * (1 - maxm)) / (kappa * np.cosh(kappa))
        self.cvec = np.cos(math.pi * self.x / 2)
        self.Kqvec = np.cosh(math.pi * (1 - self.x) / 2) / np.cosh(math.pi / 2)
        self.seed = np.array([float(Pi_star * mp.e ** (-Pi_star * xx) / (1 - mp.e ** (-Pi_star))) for xx in self.x])
        self.seed /= np.sum(self.w * self.seed)

    def moments(self, sigma: np.ndarray) -> tuple[float, float]:
        return float(np.sum(self.w * sigma * self.cvec)), float(np.sum(self.w * sigma * self.Kqvec))

    def ratio(self, g: float) -> float:
        return float(((g - float(rF1)) ** 2) / (1 + float(rF1) ** 2))

    def T_s(self, sigma: np.ndarray) -> np.ndarray:
        return self.minxy.dot(self.w * sigma)

    def T_q(self, sigma: np.ndarray) -> np.ndarray:
        return self.Gq.dot(self.w * sigma)

    def fixed_point(self, Sigma0: float, init: np.ndarray | None = None, tol: float = 1e-12, maxit: int = 3000):
        if init is None:
            sigma = self.seed.copy()
        else:
            sigma = init.copy()
        sigma /= np.sum(self.w * sigma)
        for it in range(maxit):
            g, S = self.moments(sigma)
            R = self.ratio(g)
            Phi = Sigma0 * (self.T_s(sigma) - R * self.T_q(sigma))
            new = np.exp(-Phi)
            new /= np.sum(self.w * new)
            err = float(np.max(np.abs(new - sigma)))
            sigma = new
            if err < tol:
                break
        g, S = self.moments(sigma)
        R = self.ratio(g)
        Pi = Sigma0 * (1 - R * S)
        return {
            "sigma": sigma,
            "g": float(g),
            "S": float(S),
            "R": float(R),
            "Pi": float(Pi),
            "iterations": it + 1,
            "err": err,
            "x": self.x,
            "weights": self.w,
        }

    def canonical_sigma0(self, lo: float = 4.3, hi: float = 4.9, steps: int = 18):
        left = self.fixed_point(lo)
        right = self.fixed_point(hi, init=left["sigma"])
        if not (left["g"] < float(g_minus_F1) < right["g"]):
            raise ValueError("Bad bracket for canonical sigma0 root.")
        for _ in range(steps):
            mid = 0.5 * (lo + hi)
            init = left["sigma"] if abs(mid - lo) < abs(mid - hi) else right["sigma"]
            m = self.fixed_point(mid, init=init)
            if m["g"] < float(g_minus_F1):
                lo, left = mid, m
            else:
                hi, right = mid, m
        mid = 0.5 * (lo + hi)
        root = self.fixed_point(mid, init=left["sigma"])
        root["Sigma0"] = mid
        root["T_hat"] = math.sqrt(9 * mid / 20)
        return root


def mpf_str(x, digits: int = 15) -> str:
    return mp.nstr(x, digits)

