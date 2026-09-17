#!/usr/bin/env python3
"""Exact controls for nonlinearPoleV2; no physical engine or numerical grids run."""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp
import S11c_d_nonlinear_pole_contract_probe as original

ROOT = original.ROOT
z = original.z
checks = []


def zero(name, operand):
    values = list(operand) if isinstance(operand, sp.MatrixBase) else [sp.sympify(operand)]
    reduced = [sp.simplify(v) for v in values]
    checks.append({'name': name, 'residuals': [str(v) for v in reduced],
                   'satisfied': all(v == 0 for v in reduced)})


def nonzero(name, operand):
    values = list(operand) if isinstance(operand, sp.MatrixBase) else [sp.sympify(operand)]
    reduced = [sp.simplify(v) for v in values]
    checks.append({'name': name, 'mutationResiduals': [str(v) for v in reduced],
                   'satisfied': any(v != 0 for v in reduced)})


def mat(value):
    return original.matrix(value)


def principal_part(name, pencil):
    inv = pencil.inv()
    orders = []
    for entry in inv:
        if entry != 0:
            num, den = sp.fraction(sp.cancel(entry))
            orders.append(original.vanishing_order(den, 0) - original.vanishing_order(num, 0))
    order = max(orders)
    coefficients = [original.residue(z**(j-1) * inv, 0) for j in range(1, order+1)]
    principal = sum((c/z**j for j, c in enumerate(coefficients, 1)), sp.zeros(pencil.rows))
    remainder = (inv-principal).applyfunc(sp.cancel)
    log_integral = original.residue(inv*pencil.diff(z), 0)
    log_from_coefficients = sum((c*pencil.diff(z, j).subs(z, 0)/sp.factorial(j-1)
                                 for j, c in enumerate(coefficients, 1)), sp.zeros(pencil.rows))
    zero(name+'LogarithmicLaurentConvolution', log_integral-log_from_coefficients)
    for j in range(1, order+1):
        zero(name+'RemainderCoefficient'+str(j), original.residue(z**(j-1)*remainder, 0))
        zero(name+'LeftLaurentEquation'+str(j), original.residue(z**(j-1)*pencil*principal, 0))
        zero(name+'RightLaurentEquation'+str(j), original.residue(z**(j-1)*principal*pencil, 0))
    nonzero(name+'HighestPrincipalCoefficient', coefficients[-1])
    return {'pencil': mat(pencil), 'inverse': mat(inv), 'order': order,
            'coefficients': [mat(c) for c in coefficients], 'holomorphicRemainder': mat(remainder),
            'logarithmicIntegral': mat(log_integral),
            'logarithmicIdempotencyResidual': mat(log_integral**2-log_integral)}


def chains(name, pencil, depth):
    size = pencil.rows
    records = []
    for length in range(1, depth+1):
        block = sp.zeros(size*length)
        for i in range(length):
            for j in range(i+1):
                block[i*size:(i+1)*size, j*size:(j+1)*size] = pencil.diff(z, i-j).subs(z, 0)/sp.factorial(i-j)
        basis = sp.Matrix.hstack(*block.nullspace())
        zero(name+'ChainResidual'+str(length), block*basis)
        zero(name+'ChainCompleteness'+str(length), basis.rank()+block.rank()-block.cols)
        records.append({'length': length, 'taylorSystem': mat(block), 'fullBasis': mat(basis),
                        'nullity': basis.cols, 'rank': block.rank()})
    return records


def semisimple():
    base = sp.diag(z, 2*z, z+1)
    e = sp.Matrix([[1, 1, 0], [0, 1, 1], [0, 0, 1]])
    f = sp.Matrix([[1, 0, 1], [1, 1, 0], [0, 0, 1]])
    pencil = e*base*f
    v = sp.Matrix.hstack(*pencil.subs(z, 0).nullspace())
    w = sp.Matrix.hstack(*pencil.subs(z, 0).T.nullspace()).T
    derivative = pencil.diff(z).subs(z, 0)
    pairing = w*derivative*v
    residue = v*pairing.inv()*w
    px, py = residue*derivative, derivative*residue
    s, t = sp.Matrix([[1, 1], [0, 1]]), sp.Matrix([[1, 0], [1, 1]])
    vb, wb = v*s, t*w
    rb = vb*(wb*derivative*vb).inv()*wb
    zero('SemisimpleDirectResidue', residue-original.residue(pencil.inv(), 0))
    zero('FullRightNullspace', pencil.subs(z, 0)*v)
    zero('FullLeftNullspace', w*pencil.subs(z, 0))
    zero('RightBasisCompleteness', v.rank()+pencil.subs(z, 0).rank()-pencil.cols)
    zero('LeftBasisCompleteness', w.rank()+pencil.subs(z, 0).rank()-pencil.rows)
    zero('FullPairingRank', pairing.rank()-v.cols)
    zero('FieldIdempotency', px**2-px)
    zero('SourceIdempotency', py**2-py)
    zero('FieldRangeIdentity', px*v-v)
    zero('ModalBasisCovariance', rb-residue)
    zero('FieldCoordinateCovariance', px-f.inv()*sp.diag(1, 1, 0)*f)
    zero('SourceCoordinateCovariance', py-e*sp.diag(1, 1, 0)*e.inv())
    nonzero('OmittedPairingBasisChangeRejected', vb*pairing.inv()*w-residue)
    injection = sp.Matrix([[1, 0], [0, 1], [1, 1]])
    observation = sp.Matrix([[1, 2, 0], [0, 1, 1]])
    direct = original.residue(observation*base.inv()*injection, 0)
    changed_injection, changed_observation = e*injection, observation*f
    response = changed_observation*residue*changed_injection
    zero('TypedOverlapCoordinateCovariance', response-direct)
    wrong = changed_observation*px*changed_injection
    nonzero('FieldProjectorOnSourceRejected', wrong-response)
    return {'pencil': mat(pencil), 'leftCoordinateMap': mat(e), 'rightCoordinateMap': mat(f),
            'rightBasis': mat(v), 'leftBasis': mat(w), 'pairing': mat(pairing),
            'rightBasisChange': mat(s), 'leftBasisChange': mat(t),
            'residue': mat(residue), 'fieldProjection': mat(px), 'sourceProjection': mat(py),
            'forcingInjection': mat(injection), 'observation': mat(observation),
            'responseResidue': mat(response), 'incorrectProjectorResponse': mat(wrong),
            'mapTypes': {'residue': ['Y', 'X'], 'fieldProjection': ['X', 'X'],
                         'sourceProjection': ['Y', 'Y'], 'forcingInjection': ['U', 'Y'],
                         'observation': ['X', 'Vo']}}


def realization():
    polynomial = z**2
    coefficients = sp.Poly(polynomial, z).all_coeffs()
    a = sp.Matrix([[0, 1], [-coefficients[-1]/coefficients[0], -coefficients[-2]/coefficients[0]]])
    b, c = sp.eye(2)[:, 1], sp.eye(2)[0, :]
    resolvent = (z*sp.eye(2)-a).inv()
    projector = original.residue(resolvent, 0)
    zero('RealizationTransfer', c*resolvent*b-sp.Matrix([[1/polynomial]]))
    zero('RealizedRieszIdempotency', projector**2-projector)
    zero('RealizationDeterminant', (z*sp.eye(2)-a).det()-polynomial)
    zero('StateNilpotentOrderTwo', a**2)
    zero('CompressedResidue', c*projector*b-original.residue(sp.Matrix([[1/polynomial]]), 0))
    nonzero('HigherOrderCouplingSurvivesZeroResidue', c*a*projector*b)
    source, observe = 1+3*z, 2+5*z
    response = sp.Matrix([[observe/polynomial*source]])
    residue = original.residue(response, 0)
    naive = observe.subs(z, 0)*original.residue(sp.Matrix([[1/polynomial]]), 0)*source.subs(z, 0)
    nonzero('FrozenForcingObservationAtDoublePoleRejected', naive-residue)
    degree = sp.degree(observe*source, z)
    product_route = sum(sp.diff(observe*source, z, n).subs(z, 0)*z**(n-2)/sp.factorial(n)
                        for n in range(degree+1))
    zero('AnalyticMapsFullLaurentReconstruction', response-sp.Matrix([[product_route]]))
    return {'scalarPencil': str(polynomial), 'stateOperator': mat(a), 'stateInjection': mat(b),
            'fieldRecovery': mat(c), 'stateResolvent': mat(resolvent), 'stateProjector': mat(projector),
            'forcing': str(source), 'observation': str(observe), 'physicalResponse': mat(response),
            'physicalResponseResidue': mat(residue), 'incorrectFrozenMapResidue': mat(naive)}


def perturbation():
    scale = sp.Symbol('perturbationScale', positive=True)
    angle = sp.Symbol('angle', real=True)
    original_pencil, perturbed = z**2, z**2-scale
    roots = sp.solve(perturbed, z)
    circle_norm = sp.simplify(sp.Abs((-scale/original_pencil).subs(z, sp.exp(sp.I*angle))))
    count = sum(sp.residue(sp.diff(perturbed, z)/perturbed, z, root) for root in roots)
    old_count = sp.residue(sp.diff(original_pencil, z)/original_pencil, z, 0)
    zero('CirclePerturbationNorm', circle_norm-scale)
    zero('HomotopyAlgebraicCount', count-old_count)
    for root in roots:
        zero('PerturbedRoot'+str(root), perturbed.subs(z, root))
    nonzero('DistinctRootCountPreservationRejected', len(roots)-len(sp.solve(original_pencil, z)))
    amplification = sp.limit(sp.sqrt(scale)/scale, scale, 0, dir='+')
    return {'original': str(original_pencil), 'perturbed': str(perturbed),
            'contour': 'unit circle, with 0 < perturbationScale < 1',
            'relativeErrorNorm': str(circle_norm), 'computedRoots': [str(r) for r in roots],
            'originalAlgebraicCount': str(old_count), 'perturbedAlgebraicCount': str(count),
            'displacementOverErrorLimit': str(amplification)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    cases = {'double': sp.Matrix([[z**2]]), 'mixed': sp.diag(z**2, z),
             'affineJordan': sp.Matrix([[z, -1], [0, z]]), 'semisimple': sp.diag(z, 2*z, z+1)}
    parts = {name: principal_part(name, pencil) for name, pencil in cases.items()}
    chain_records = {name: chains(name, pencil, 3) for name, pencil in cases.items()}
    semisimple_record = semisimple()
    realization_record = realization()
    perturbation_record = perturbation()
    for name in ('double', 'mixed'):
        pencil = cases[name]
        v = sp.Matrix.hstack(*pencil.subs(z, 0).nullspace())
        w = sp.Matrix.hstack(*pencil.subs(z, 0).T.nullspace()).T
        pairing = w*pencil.diff(z).subs(z, 0)*v
        nonzero(name+'SingularPairingFailsFullRank', pairing.rank()-v.cols)
        j = original.residue(pencil.inv()*pencil.diff(z), 0)
        nonzero(name+'UnrestrictedProjectorRejected', j**2-j)
    historical_path = ROOT/'_measurements/S11c_d_nonlinear_pole_contract_probe.json'
    historical = json.loads(historical_path.read_text())
    source_paths = [Path(__file__).resolve(), historical_path,
                    ROOT/'directives/S11c_d_NONLINEAR_POLE_CONTRACT.md',
                    ROOT/'directives/S11c_d_sympy_build_PROGRAM_BRIEF.md',
                    ROOT/'directives/S11c_d_sympy_build_directive.md',
                    ROOT/'directives/_measurements/S11c_d_sympy_build_completion_addendum.md']
    for name, expected in historical['sourceFiles'].items():
        actual = original.digest(ROOT/name)
        checks.append({'name': 'HistoricalSourceIdentity:'+name, 'expected': expected,
                       'actual': actual, 'satisfied': actual == expected})
    record = {'schema': 'nonlinearPoleV2ExactControls', 'scope': 'Synthetic dimensionless mathematical controls only.',
              'metadata': {'unitLTM': [0, 0, 0], 'epsEtaSigmaWOrder': [0, 0, 0], 'lambdaOrder': 0},
              'sourceFiles': {str(p.relative_to(ROOT)): original.digest(p) for p in source_paths},
              'principalParts': parts, 'fullChainSystems': chain_records,
              'semisimple': semisimple_record, 'realization': realization_record,
              'perturbation': perturbation_record, 'checks': checks}
    payload = json.dumps(record, indent=2)+'\n'
    args.output.write_text(payload)
    print(payload, end='')
    if not all(c['satisfied'] for c in checks):
        raise RuntimeError('Exact nonlinear-pole control failed; inspect saved operands and residuals')


if __name__ == '__main__':
    main()
