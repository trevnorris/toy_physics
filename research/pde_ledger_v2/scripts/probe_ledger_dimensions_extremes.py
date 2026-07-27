#!/usr/bin/env python3
"""Exercise the stage038 and stage042 dimension requirements through the shared API."""

from __future__ import annotations

from fractions import Fraction
from functools import reduce
from operator import mul

import sympy as sp

from ledger_dimensions import Dimension, DimensionBasis


def dimension_sum(*dimensions: Dimension) -> Dimension:
    """Add exponent vectors using Dimension's public composition operation."""

    if not dimensions:
        raise ValueError("dimension_sum needs at least one dimension")
    return reduce(mul, dimensions)


def probe_stage038() -> None:
    basis = DimensionBasis("M", "L", "T", "E", render="tuple")
    zero = basis()
    a_e = basis(0, 1, 0, 1)
    q_t = basis(1, 0, -1, 0)
    c_gamma2 = basis(0, 2, -2, 0)
    mu_r = basis(2, 1, -4, -1)
    c_e2 = (mu_r * a_e) / (q_t**2)

    added = q_t * q_t
    ratio = (added * c_gamma2) / (mu_r * a_e)
    scaled = q_t**2

    assert basis.axes == ("M", "L", "T", "E")
    assert added == scaled
    assert ratio == zero
    assert c_e2 == c_gamma2

    actual_with_sentinel: list[Dimension | None] = [
        a_e,
        q_t,
        c_gamma2,
        mu_r,
        ratio,
        added,
        scaled,
        zero,
        None,
    ]
    expected_with_sentinel: list[Dimension | None] = [
        basis(0, 1, 0, 1),
        basis(1, 0, -1, 0),
        basis(0, 2, -2, 0),
        basis(2, 1, -4, -1),
        basis(),
        basis(2, 0, -2, 0),
        basis(2, 0, -2, 0),
        None,
        basis(),
    ]
    assert actual_with_sentinel != expected_with_sentinel

    print("stage038 axes: PASS ('M', 'L', 'T', 'E')")
    print("stage038 add/subtract/integer-scale/dimensionless equality: PASS")
    print("stage038 8-element heterogeneous list inequality with shifted None sentinel: PASS")


def probe_stage042() -> None:
    basis = DimensionBasis("stiffness", "L", "T", render="tuple")
    zero = basis()
    stiffness = basis(Fraction(1), Fraction(0), Fraction(0))
    length = basis(Fraction(0), Fraction(1), Fraction(0))
    speed = basis(Fraction(0), Fraction(1), Fraction(-1))
    charge = basis(Fraction(1, 2), Fraction(3, 2), Fraction(-1))

    assert basis.axes == ("stiffness", "L", "T")
    assert charge.exponents == {
        "stiffness": sp.Rational(1, 2),
        "L": sp.Rational(3, 2),
        "T": sp.Rational(-1),
    }
    assert not any(
        exponent.has(sp.Float) for exponent in charge.exponents.values()
    )

    n_way_added = dimension_sum(stiffness, charge, charge)
    scaled = charge**Fraction(2)
    subtracted = n_way_added / stiffness
    assert subtracted == scaled
    assert (charge, speed, length) == (
        basis(Fraction(1, 2), Fraction(3, 2), Fraction(-1)),
        basis(0, 1, -1),
        basis(0, 1, 0),
    )

    homogeneous_terms = (
        dimension_sum(stiffness, scaled),
        dimension_sum(stiffness, charge, charge),
        dimension_sum(charge, stiffness, charge),
    )
    assert len(set(homogeneous_terms)) == 1
    equal_left = dimension_sum(stiffness, scaled)
    equal_right = dimension_sum(stiffness, charge, charge)
    assert equal_left == equal_right
    assert hash(equal_left) == hash(equal_right)
    inverse_time = speed / length
    assert inverse_time == basis(0, 0, -1)

    print("stage042 axes: PASS ('stiffness', 'L', 'T')")
    print("stage042 Fraction exactness: PASS (1/2, 3/2, -1; no Float)")
    print("stage042 varargs add/scalar multiply/subtraction/tuple equality: PASS")
    print("stage042 set homogeneity/hashability/equal hashes: PASS")


def main() -> None:
    probe_stage038()
    probe_stage042()
    print("VERDICT: shared dimension module is not blocked by stage038 or stage042")


if __name__ == "__main__":
    main()
