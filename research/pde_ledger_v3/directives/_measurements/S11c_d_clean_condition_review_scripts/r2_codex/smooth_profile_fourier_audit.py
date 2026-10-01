#!/usr/bin/env python3
"""Numerical counterexample to 'smooth profile implies exponential Fourier tail'.

The standard compact C-infinity bump exp[-1/(1-x^2)] for |x|<1 is smooth
but nonanalytic at its endpoints.  Its Fourier transform has a stretched-
exponential saddle tail, not a generic exp(-c |k|) tail.  We compute the
transform directly at high precision, splitting at half-periods.
"""

import mpmath as mp

mp.mp.dps = 80


def bump(x):
    if abs(x) >= 1:
        return mp.mpf("0")
    return mp.e ** (-1 / (1 - x*x))


def transform(k):
    step = mp.pi / k
    points = [mp.mpf("0")]
    x = step
    while x < 1:
        points.append(x)
        x += step
    points.append(mp.mpf("1"))
    return 2 * mp.fsum(mp.quad(lambda y: bump(y) * mp.cos(k*y), [a, b])
                       for a, b in zip(points[:-1], points[1:]))


print("PROFILE", "exp(-1/(1-x^2)) for |x|<1, else 0 (C-infinity, nonanalytic)")
for k in (20, 40, 80, 160, 320):
    value = transform(mp.mpf(k))
    decay = -mp.log(abs(value))
    print(
        "K", k,
        "ABS_F", mp.nstr(abs(value), 24),
        "NEG_LOG_OVER_K", mp.nstr(decay / k, 18),
        "NEG_LOG_OVER_SQRT_K", mp.nstr(decay / mp.sqrt(k), 18),
    )

