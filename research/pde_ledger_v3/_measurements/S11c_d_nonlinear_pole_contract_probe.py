#!/usr/bin/env python3
"""Exact finite-dimensional witnesses for the pending nonlinear-pole contract issue.

This diagnostic neither imports the physics engine nor changes its authority.
All pencils and the coordinate z are dimensionless; no physical inputs are used.
"""
import hashlib
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
z = sp.Symbol('z')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix(value):
    return [[str(sp.simplify(v)) for v in row] for row in value.tolist()]


def vanishing_order(polynomial, point):
    shifted = sp.Poly(sp.expand(polynomial.subs(z, z + point)), z)
    return min(power[0] for power, coefficient in shifted.terms() if coefficient != 0)


def residue(value, point):
    return value.applyfunc(lambda v: sp.residue(v, z, point))


def compute(name, pencil, points):
    inverse = pencil.inv()
    derivative = pencil.diff(z)
    logarithmic_derivative = inverse * derivative
    integral = sum((residue(logarithmic_derivative, p) for p in points), sp.zeros(pencil.rows))
    poles = []
    for point in points:
        at_root = pencil.subs(z, point)
        right = sp.Matrix.hstack(*at_root.nullspace())
        left = sp.Matrix.hstack(*at_root.T.nullspace()).T
        pairing = left * derivative.subs(z, point) * right
        actual_residue = residue(inverse, point)
        orders = []
        for entry in inverse:
            if entry != 0:
                numerator, denominator = sp.fraction(sp.cancel(entry))
                orders.append(vanishing_order(denominator, point) - vanishing_order(numerator, point))
        record = {
            'point': str(point), 'nullity': len(at_root.nullspace()),
            'determinantZeroOrder': vanishing_order(pencil.det(), point),
            'inversePoleOrder': max(orders), 'residue': matrix(actual_residue),
            'derivativePairing': matrix(pairing), 'derivativePairingRank': pairing.rank(),
        }
        if pairing.rank() == pairing.rows:
            modal_residue = right * pairing.inv() * left
            local_projector = modal_residue * derivative.subs(z, point)
            record.update(
                modalResidueResidual=matrix(modal_residue - actual_residue),
                localProjector=matrix(local_projector),
                localIdempotencyResidual=matrix(local_projector**2 - local_projector),
                localContourResidual=matrix(residue(logarithmic_derivative, point) - local_projector),
            )
        poles.append(record)
    return {
        'name': name, 'pencil': matrix(pencil), 'inverse': matrix(inverse),
        'inverseResidual': matrix(pencil * inverse - sp.eye(pencil.rows)),
        'contourIntegral': matrix(integral), 'idempotencyResidual': matrix(integral**2 - integral),
        'trace': str(sp.trace(integral)), 'poles': poles,
    }


def main():
    cases = [
        ('scalarSimple', sp.Matrix([[z]]), (0,)),
        ('scalarDouble', sp.Matrix([[z**2]]), (0,)),
        ('semisimpleNullityTwo', sp.diag(z, 2*z, z + 1), (0,)),
        ('affineJordan', sp.Matrix([[z, -1], [0, z]]), (0,)),
        ('scalarOneOfTwoRoots', sp.Matrix([[z**2 - 1]]), (1,)),
        ('scalarBothRoots', sp.Matrix([[z**2 - 1]]), (-1, 1)),
    ]
    paths = [Path(__file__).resolve(), ROOT/'directives/S11c_d_SHARED_PHYSICS.md',
             ROOT/'scripts/S11c_d_mixing_scattering_sympy_audit.py',
             ROOT/'_measurements/S11c_d_wide_three_momentum.py']
    print(json.dumps({
        'scope': 'Exact rational matrix-pencil diagnostic; no S11c-d physical result or specification change.',
        'sympyVersion': sp.__version__,
        'metadata': {'unitLTM': [0, 0, 0], 'epsEtaSigmaWOrder': [0, 0, 0],
                     'lambdaOrder': 0, 'meaning': 'dimensionless synthetic pencils, exact ungraded calculation'},
        'sourceFiles': {str(p.relative_to(ROOT)): digest(p) for p in paths},
        'cases': [compute(*case) for case in cases],
    }, indent=2))


if __name__ == '__main__':
    main()
