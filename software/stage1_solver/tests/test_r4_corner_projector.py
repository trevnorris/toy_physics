"""Captured four-corner projection regression, independent of saved PDE checkpoints."""
import argparse
import hashlib
from itertools import permutations
import json
from pathlib import Path
import sys
from unittest.mock import patch

import mpmath as mp
import numpy as np

TOOLS = Path(__file__).resolve().parents[1] / 'tools'
sys.path.insert(0, str(TOOLS))
import mouth_r4_corner_projection as projector

FIXTURE = Path(__file__).parent / 'fixtures' / 'r4_step288_corner_rows.json'


def reference(data):
    """Solve the binding equalities independently, then certify this active set."""
    with mp.workdps(90):
        rows = mp.matrix(data['normals'])
        n = mp.matrix([list(row / mp.norm(row)) for row in
                       (rows[i, :] for i in range(rows.rows))])
        diagonal = mp.diag(data['diagonal'])
        p = mp.matrix(data['raw'])
        binding = data['binding_rows_90digit']
        a = mp.matrix([[n[i, j] for j in range(n.cols)] for i in binding])
        # Independent projection in whitened primal coordinates using QR,
        # without solving the production dual Gram system.
        sqrt_d = mp.diag([mp.sqrt(d) for d in data['diagonal']])
        q, r = mp.qr((a * sqrt_d).T)
        free = len(binding)
        q_binding = q[:, :free]
        triangular = r[:free, :]
        correction = mp.lu_solve(triangular.T, -(a * p))
        d = p + sqrt_d * q_binding * correction
        multipliers = mp.matrix(rows.rows, 1)
        for i, value in zip(binding, mp.lu_solve(triangular, correction)):
            multipliers[i] = value
        slack = n * d
        assert min(multipliers) >= 0 and min(slack) >= -mp.mpf('1e-78')
        assert max(abs(multipliers[i] * slack[i]) for i in range(rows.rows)) <= mp.mpf('1e-78')
        assert mp.norm(d - p - diagonal * n.T * multipliers) <= mp.mpf('1e-78')
        return np.asarray(list(d), float), np.asarray(list(multipliers), float), np.asarray(list(slack), float)


def permutation_audit():
    data = json.loads(FIXTURE.read_bytes())
    expect, expect_lambda, expect_slack = reference(data)
    size = data['full_size']
    columns = data['indices']
    rows = np.zeros((4, size))
    rows[:, columns] = data['normals']
    # Untouched coordinates are zero-padded at their original indices so
    # production reductions see the original order of the nonzero terms.
    diagonal = np.ones(size)
    diagonal[columns] = data['diagonal']
    raw = np.zeros(size)
    raw[columns] = data['raw']
    expected = np.zeros(size)
    expected[columns] = expect
    # The original full energy gradient is validated on the PDE replay.
    # This separate test chooses the metric descent gradient for the same QP.
    gradient = -raw / diagonal
    unit = rows / np.linalg.norm(rows, axis=1)[:, None]
    rows_hp = unit[:, columns].astype(np.longdouble)
    result = []
    for order in permutations(range(4)):
        permutation = np.array(order)
        direction, info = projector.project_halfspaces(raw, rows[permutation], diagonal=diagonal,
                                                       gradient=gradient, return_info=True)
        multipliers = np.empty(4)
        multipliers[permutation] = info['multipliers']
        m = multipliers.astype(np.longdouble)
        local = direction[columns].astype(np.longdouble)
        slack = rows_hp @ local
        inverse_metric_residual = ((local - raw[columns]) / diagonal[columns]) - rows_hp.T @ m
        primal = float(max(0, -slack.min()))
        dual = float(max(0, -m.min()))
        complementarity = float(np.max(np.abs(m * slack)))
        stationarity = float(np.sqrt(inverse_metric_residual @ inverse_metric_residual))
        residual = float(np.sqrt(primal**2 + dual**2 + complementarity**2 + stationarity**2))
        direction_error = float(np.max(np.abs(direction - expected)))
        multiplier_error = float(np.max(np.abs(multipliers - expect_lambda)))
        slack_error = float(np.max(np.abs(slack - expect_slack)))
        slope = float(gradient @ direction)
        assert max(primal, complementarity, residual) <= 1e-10
        assert dual <= 1e-12 and slope < 0
        assert max(direction_error, multiplier_error, slack_error) <= 1e-10
        for k in ('primal_violation', 'complementarity', 'stationarity', 'kkt_norm'):
            assert 0 <= info[k] <= 1e-10, (order, k, info)
        assert 0 <= info['dual_violation'] <= 1e-12 and info['slope'] < 0
        result.append(dict(permutation=list(order), method=info['method'],
                           direction_error=direction_error, multiplier_error=multiplier_error,
                           slack_error=slack_error, primal_violation=primal, dual_violation=dual,
                           complementarity=complementarity, stationarity=stationarity,
                           kkt_norm=residual, slope=slope))
    return dict(precision_digits=90, compared_dofs=size, local_support_size=len(columns),
                reference_direction=expect.tolist(), reference_multipliers=expect_lambda.tolist(),
                reference_slacks=expect_slack.tolist(), permutations=result)


def test_all_step288_permutations():
    assert len(permutation_audit()['permutations']) == 24


def test_scaled_affine_redundant_rows():
    rows = np.array([[1., 0.], [20., 0.], [0., 2.]])
    rhs = np.array([1., 20., -2.])
    raw = np.array([-1., 2.])
    direction, info = projector.project_halfspaces(raw, rows, diagonal=np.array([1., 3.]), rhs=rhs,
                                                   gradient=np.array([-2., -1.]), return_info=True)
    np.testing.assert_allclose(direction, [1., 2.], rtol=0., atol=1e-10)
    assert max(info['primal_violation'], info['kkt_norm']) <= 1e-10


def test_corrupt_float64_claim_is_rejected():
    rows = np.array([[1., 0., 0.], [0., 1., 0.]])
    raw = np.array([-1., -2., 3.])
    with patch.object(projector, 'nnls', return_value=(np.zeros(2), 0.0)) as fake:
        direction, info = projector.project_halfspaces(raw, rows, diagonal=np.ones(3), gradient=-raw,
                                                       return_info=True)
    assert fake.called
    np.testing.assert_allclose(direction, [0., 0., 3.], rtol=0., atol=1e-10)
    assert max(info['primal_violation'], info['kkt_norm']) <= 1e-10
    return info['method']


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    path = parser.parse_args().output
    report = permutation_audit()
    test_scaled_affine_redundant_rows()
    method = test_corrupt_float64_claim_is_rejected()
    report.update(passed=True, affine_redundant=True, corrupt_nnls_method=method,
                  fixture_sha256=hashlib.sha256(FIXTURE.read_bytes()).hexdigest(),
                  production_sha256=hashlib.sha256(Path(projector.__file__).read_bytes()).hexdigest(),
                  suite_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    with path.open('x') as stream:
        json.dump(report, stream, allow_nan=False, indent=2)
        stream.write('\n')
    print(json.dumps({key: max(row[key] for row in report['permutations'])
                      for key in ('direction_error','multiplier_error','slack_error','kkt_norm')}))
