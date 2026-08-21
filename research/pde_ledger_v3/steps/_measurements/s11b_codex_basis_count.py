#!/usr/bin/env python3
"""Independent S11b quadratic-basis count modulo total divergences."""

import sympy as sp


k = sp.Matrix(sp.symbols("k1:4"))
u = sp.Matrix(sp.symbols("u1:4"))
G_fourier = k * u.T  # the common factor i is irrelevant to linear dependence

div_squared = sp.expand(sp.trace(G_fourier) ** 2)
trace_g_squared = sp.expand(sp.trace(G_fourier * G_fourier))
frobenius_squared = sp.expand(sp.trace(G_fourier.T * G_fourier))

assert sp.expand(div_squared - trace_g_squared) == 0
assert sp.expand(frobenius_squared - div_squared) != 0

sector_counts = {
    "grad_u modulo divergence": 2,
    "theta,e scalars": 3,
    "grad(theta),grad(e)": 3,
    "scalar-times-div(u)": 2,
}

print(f"Fourier image of (div u)^2 = {div_squared}")
print(f"Fourier image of tr((grad u)^2) = {trace_g_squared}")
print(f"difference = {sp.expand(div_squared - trace_g_squared)}")
for name, count in sector_counts.items():
    print(f"{name}: {count}")
print(f"basis count modulo total divergences = {sum(sector_counts.values())}")
