#!/usr/bin/env python3
"""Measure PROFILE_GRADE_SUBS vs directive Bridge D; protected-atom coexistence; HeldDiv drop."""
from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import sympy as sp

import S11c_b_brane_operator_sympy_audit as E
import S11c_b_cross_engine_comparator as C

subs = E.PROFILE_GRADE_SUBS
print(f"PROFILE_GRADE_SUBS_ENTRY_COUNT = {len(subs)}")
print(f"PROFILE_GRADE_SUBS_KEYS = {sorted(str(k) for k in subs)}")
print(f"SECOND_JET_W_bg_d1d1_IN_SUBS = {sp.Symbol('W_bg_d1d1') in subs}")
print(f"SIGMA_W_IDENTICAL_W0_ETA = {E.sigma_W == E.W0 * E.eta_bg}")

prot_names = [
    "gamma_s11cb_w_bg_07",
    "gamma_s11cb_w_bg_10",
    "gamma_s11cb_mu_r_bg_07",
    "gamma_s11cb_mu_r_bg_10",
    "gammaWidthDivGradTheta",
    "gammaWidthDivGradEw",
    "gammaModulusDivGradTheta",
    "gammaModulusDivGradEw",
]
prot = {sp.Symbol(n) for n in prot_names}
print(
    "PROTECTED_IN_PROFILE_GRADE_SUBS_KEYS = "
    f"{sorted(str(p) for p in prot & set(subs))}"
)

W_bg_d1 = E.grad_W[0]
g = sp.Symbol("gammaWidthDivGradTheta")
R = W_bg_d1 * sp.Symbol("u_1_t") + g * sp.Symbol("theta_probe_d1")
R_bridged = sp.expand(R.xreplace(subs))
print(f"TOY_R = {R}")
print(f"TOY_R_AFTER_BRIDGE_D = {R_bridged}")
print(
    "PROTECTED_COEFF_UNCHANGED = "
    f"{sp.expand(R_bridged.coeff(g) - R.coeff(g)) == 0}"
)

path = Path("/tmp/S11c_b_adjudicated_run.out")
prot_re = re.compile(
    r"^(gamma_s11cb_(?:w_bg|mu_r_bg)_(?:07|10)|"
    r"gamma(?:Width|Modulus)DivGrad(?:Theta|Ew))$"
)
stats: Counter[tuple[str, str]] = Counter()
cur = None
for line in path.read_text().splitlines():
    if line.startswith("ALGEBRAIC FLAG ") or line.startswith("CONTAINER FLAG "):
        fam = line.split()[2]
        cur = fam
    elif cur is not None and line.startswith("A_minus_B "):
        syms = set(re.findall(r"Symbol\('([^']+)'\)", line))
        hits = frozenset(s for s in syms if prot_re.match(s))
        stats[(cur, "HAS_PROTECTED" if hits else "NO_PROTECTED")] += 1
        stats[(cur, "HELDDIV" if "HeldDiv" in line else "NO_HELDDIV")] += 1
        cur = None
    elif line.startswith(
        (
            "ALGEBRAIC MATCH",
            "STRUCTURE",
            "COVERAGE",
            "NAMESPACE",
            "CONTAINER MATCH",
        )
    ):
        cur = None
for key in sorted(stats):
    print(f"BASELINE {key[0]} {key[1]} = {stats[key]}")

HeldDiv = sp.Function("HeldDiv")
expr = sp.Symbol("bulk") + HeldDiv(sp.Tuple(sp.Symbol("v1"), 0, 0))
dropped = C._drop_held_divergences(expr)
print(f"TOY_STRONG = {expr}")
print(f"AFTER_DROP_HELDDIV = {dropped}")

a, phi = sp.symbols("a phi")
a_d1, phi_d1 = sp.symbols("a_d1 phi_d1")
R_div = a_d1 * phi + a * phi_d1
bases = {"a", "phi"}


def split(name: str):
    match = re.fullmatch(r"^(.*)_(d[123](?:d[123])*)$", name)
    if match and match.group(1) in bases:
        return match.group(1), tuple(
            int(x) for x in re.findall(r"d([123])", match.group(2))
        )
    if name in bases:
        return name, ()
    return None


def jet(base: str, dirs: tuple[int, ...]) -> sp.Symbol:
    ordered = tuple(sorted(dirs))
    suffix = "" if not ordered else "_" + "".join(f"d{i}" for i in ordered)
    return sp.Symbol(base + suffix)


def formal_dx(expression: sp.Expr, direction: int) -> sp.Expr:
    terms = []
    for symbol in expression.free_symbols:
        decoded = split(symbol.name)
        if decoded is None:
            continue
        base, directions = decoded
        partial = sp.diff(expression, symbol)
        if partial != 0:
            terms.append(partial * jet(base, (*directions, direction)))
    return sp.Add(*terms)


V1 = a * phi
print(f"FIXTURE_R_DIV_MINUS_D1_V = {sp.expand(R_div - formal_dx(V1, 1))}")
print(f"FIXTURE_BULK_MINUS_D1_V = {sp.expand(a * phi_d1 - formal_dx(V1, 1))}")
print(f"EULER_PLACEHOLDER_BULK = {C.modulo_total_divergence(a * phi_d1)}")
