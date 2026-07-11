#!/usr/bin/env python3
"""MIDWAY knob-audit codimension dry-run: independent SymPy engine.

The algebraic dimension is obtained from the leading-monomial ideal of a
grevlex Groebner basis.  An exact positive smooth witness then certifies that
the complex Krull dimension is attained by the positive real locus required
by the audit contract.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from itertools import combinations

import sympy as sp


@dataclass(frozen=True)
class Block:
    name: str
    variables: tuple[sp.Symbol, ...]
    baseline: tuple[sp.Expr, ...]
    cases: tuple[tuple[str, tuple[sp.Expr, ...]], ...]
    witness: dict[sp.Symbol, sp.Expr]


def groebner_dimension(polynomials: tuple[sp.Expr, ...], variables: tuple[sp.Symbol, ...]) -> int:
    """Return Krull dimension using the initial-monomial ideal."""
    nonzero = tuple(sp.expand(p) for p in polynomials if sp.expand(p) != 0)
    if not nonzero:
        return len(variables)

    basis = sp.groebner(nonzero, *variables, order="grevlex")
    leading_supports: list[frozenset[int]] = []
    for polynomial in basis.polys:
        leading_monomial = tuple(polynomial.LM(order=basis.order))
        support = frozenset(i for i, exponent in enumerate(leading_monomial) if exponent)
        if not support:
            return -1
        leading_supports.append(support)

    # The dimension of a monomial quotient is the largest coordinate subset
    # that contains no generator support in full.
    for size in range(len(variables), -1, -1):
        for candidate in combinations(range(len(variables)), size):
            candidate_set = frozenset(candidate)
            if all(not support <= candidate_set for support in leading_supports):
                return size
    raise AssertionError("initial-ideal dimension search failed")


def certify_positive_real_dimension(
    equations: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
    witness: dict[sp.Symbol, sp.Expr],
    dimension: int,
) -> None:
    """Certify a positive smooth real point of the computed dimension."""
    for variable in variables:
        value = sp.simplify(witness[variable])
        assert value.is_real is True, f"{variable}: witness is not provably real"
        assert value.is_positive is True, f"{variable}: witness is not positive"
    for equation in equations:
        assert sp.simplify(equation.subs(witness)) == 0, f"witness violates {equation}"

    nonzero = tuple(equation for equation in equations if sp.expand(equation) != 0)
    if nonzero:
        jacobian = sp.Matrix(nonzero).jacobian(variables).subs(witness)
        rank = int(jacobian.rank())
    else:
        rank = 0
    assert len(variables) - rank == dimension, (
        f"positive witness is not smooth on a maximal-dimensional component: "
        f"ambient={len(variables)}, rank={rank}, dim={dimension}"
    )


def compute_case(block: Block, equations: tuple[sp.Expr, ...]) -> dict[str, int]:
    dim_before = groebner_dimension((), block.variables)
    dim_after = groebner_dimension(equations, block.variables)
    assert dim_before >= 0 and dim_after >= 0, "empty algebraic locus"
    certify_positive_real_dimension(equations, block.variables, block.witness, dim_after)
    return {
        "dim_before": dim_before,
        "dim_after": dim_after,
        "Delta": dim_before - dim_after,
    }


hbar, mass, big_k, rho0, cs0, xi_h, h0, scale_a, c_gamma, lambda_gamma = sp.symbols(
    "hbar mass K rho0 c_s0 xi_h h0 a c_gamma lambda_gamma", real=True
)
medium_variables = (hbar, mass, big_k, rho0, cs0, xi_h, h0, scale_a, c_gamma, lambda_gamma)
medium_baseline = (
    mass * cs0**2 - 5 * big_k * rho0**4,
    scale_a * mass * cs0 - hbar,
    xi_h**2 * mass**2 * cs0**2 - 2 * hbar**2,
    4 * h0 - mass * cs0**2,
    lambda_gamma * cs0 - c_gamma,
)
medium_witness = {
    hbar: sp.Integer(5),
    mass: sp.Integer(5),
    big_k: sp.Integer(1),
    rho0: sp.Integer(1),
    cs0: sp.Integer(1),
    xi_h: sp.sqrt(2),
    h0: sp.Rational(5, 4),
    scale_a: sp.Integer(1),
    c_gamma: sp.Integer(1),
    lambda_gamma: sp.Integer(1),
}
medium = Block(
    name="M",
    variables=medium_variables,
    baseline=medium_baseline,
    cases=(
        ("baseline", medium_baseline),
        ("C-M1", medium_baseline[:-1] + (lambda_gamma - lambda_gamma,)),
        ("C-M2", medium_baseline + (xi_h**2 - 2 * scale_a**2,)),
        ("C-M3", medium_baseline + (big_k - rho0,)),
    ),
    witness=medium_witness,
)

mu_eta, t_w, beta, k_eta, vp0_over_lc, t_omega, a_b, kappa_b, delta, sigma_wall = sp.symbols(
    "mu_eta T_w beta K_eta Vp0_over_lc T_Omega a_B kappa_B delta sigma_wall", real=True
)
wall_variables = (mu_eta, t_w, beta, k_eta, vp0_over_lc, t_omega, a_b, kappa_b, delta, sigma_wall)
wall_baseline = (
    k_eta - t_w * beta**2,
    2 * a_b * delta**2 - kappa_b,
    36 * sigma_wall**2 - 2 * a_b * kappa_b,
)
wall_witness = {
    mu_eta: sp.Integer(1),
    t_w: sp.Integer(1),
    beta: sp.Integer(1),
    k_eta: sp.Integer(1),
    vp0_over_lc: sp.Integer(1),
    t_omega: sp.Integer(1),
    a_b: sp.Integer(1),
    kappa_b: sp.Integer(1),
    delta: sp.sqrt(2) / 2,
    sigma_wall: sp.sqrt(2) / 6,
}
wall = Block(
    name="Wchi",
    variables=wall_variables,
    baseline=wall_baseline,
    cases=(
        ("baseline", wall_baseline),
        ("C-W1", (k_eta - k_eta,) + wall_baseline[1:]),
        ("C-W2", wall_baseline + (mu_eta - t_w,)),
    ),
    witness=wall_witness,
)


EXPECTED_PAYLOAD = {
    "M": {
        "baseline": {"dim_before": 10, "dim_after": 5, "Delta": 5},
        "C-M1": {"dim_before": 10, "dim_after": 6, "Delta": 4},
        "C-M2": {"dim_before": 10, "dim_after": 5, "Delta": 5},
        "C-M3": {"dim_before": 10, "dim_after": 4, "Delta": 6},
    },
    "Wchi": {
        "baseline": {"dim_before": 10, "dim_after": 7, "Delta": 3},
        "C-W1": {"dim_before": 10, "dim_after": 8, "Delta": 2},
        "C-W2": {"dim_before": 10, "dim_after": 6, "Delta": 4},
    },
}


def main() -> None:
    payload: dict[str, dict[str, dict[str, int]]] = {}
    for block in (medium, wall):
        payload[block.name] = {
            case_name: compute_case(block, equations) for case_name, equations in block.cases
        }

    # This is also the canonical cross-engine comparison structure printed by
    # the Mathematica engine.  Values above came only from the dimension routine.
    assert payload == EXPECTED_PAYLOAD, f"dimension payload mismatch: {payload!r}"

    print("BLOCK  CASE      dim_before  dim_after  Delta")
    for block_name, cases in payload.items():
        for case_name, record in cases.items():
            print(
                f"{block_name:<6} {case_name:<9} {record['dim_before']:>10}"
                f" {record['dim_after']:>10} {record['Delta']:>6}"
            )
    print(
        "NOTE: beta_2(w) is a PROFILE and the R35 radial-scalar reductions "
        "{Mtilde,Ktilde,Ttilde_Omega}=Integral[density*beta_2^2] are integral "
        "functionals, not polynomializable; they are outside this computed check "
        "and remain reasoned-tally only."
    )
    print("COMPARISON_PAYLOAD: " + json.dumps(payload, separators=(",", ":"), ensure_ascii=True))
    print("MIDWAY_CODIM_CHECK: CONSISTENT")


if __name__ == "__main__":
    main()
