#!/usr/bin/env python3
"""Compact F1–F4 fidelity: selected original finite-current algebra only.

No production module is imported and no scientific packet is opened. These
small exact-in-binary NumPy fixtures identify conventions, not a physical solve.
Run inside scripts/s11c_guarded_run.py on this host.
"""
import ast
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np

BASE = Path(__file__).resolve().parents[1]
M = BASE / '_measurements'
PATHS = {
    'current': M / 'S11c_d_continuum_currents.py',
    'response': M / 'S11c_d_continuum_response.py',
    'engine': BASE / 'scripts/S11c_d_mixing_scattering_sympy_audit.py',
    'spec': BASE / 'directives/S11c_d_SHARED_PHYSICS.md',
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def selected(path, names, owner=None):
    body = ast.parse(path.read_text()).body
    if owner:
        body = next(n.body for n in body if isinstance(n, ast.ClassDef) and n.name == owner)
    nodes = [n for n in body if isinstance(n, ast.FunctionDef) and n.name in names]
    if {n.name for n in nodes} != set(names):
        raise ValueError('selected function missing')
    for n in nodes:
        n.decorator_list = []
    return ast.Module(body=nodes, type_ignores=[])


def main():
    before = {str(p.relative_to(BASE)): digest(p) for p in PATHS.values()}
    extracted = {}
    def functions(key, names, owner=None):
        tree = selected(PATHS[key], names, owner)
        extracted[key + ':' + ','.join(names)] = hashlib.sha256(ast.dump(tree).encode()).hexdigest()
        return tree
    jets = {'np': np}
    exec(compile(functions('engine', ['multiply'], 'RectangularModeJets'), '<selected native jets>', 'exec'), jets)
    response = {'J': SimpleNamespace(multiply=jets['multiply'])}
    response_tree = functions('response', ['adjoint', 'gram'])
    exec(compile(response_tree, '<selected native response>', 'exec'), response)
    current = {'np': np, 'response': SimpleNamespace(adjoint=response['adjoint'])}
    exec(compile(functions('current', ['multiply', 'quadratic']), '<selected native current>', 'exec'), current)
    checks = []
    def check(name, actual, expected):
        if not np.array_equal(actual, expected):
            raise AssertionError(name)
        checks.append({'name': name, 'passed': True})
    def reject(name, actual, wrong):
        if np.array_equal(actual, wrong):
            raise AssertionError('undetected control: ' + name)
        checks.append({'name': name, 'passed': True, 'kind': 'wrong-formula control'})
    def q(a, j):
        return current['quadratic']({(0, 0): a}, {(0, 0): j})[(0, 0)]
    # Fixed integer/Gaussian-integer values have exact binary arithmetic here.
    a = np.array([[1], [1j]], complex)
    j = np.array([[2, 1j], [-1j, 3]], complex)
    check('complex_current_with_cross_terms', q(a, j), np.array([[3]], complex))
    reject('transpose_is_not_adjoint', q(a, j), a.T @ j @ a)
    reject('diagonal_current_omits_interference', q(a, j), q(a, np.diag(np.diag(j))))
    reject('transposed_current_changes_form', q(a, j), q(a, j.T))
    x, y = np.array([[1], [0]], complex), np.array([[0], [1j]], complex)
    cross = x.conj().T @ j @ y + y.conj().T @ j @ x
    check('full_boundary_decomposition', q(x + y, j), q(x, j) + q(y, j) + cross)
    c = np.array([[2, 1j], [0, 2]], complex)
    pulled = response['gram']({(0, 0): c}, {(0, 0): j})[(0, 0)]
    check('native_gram_is_congruence', pulled, c.conj().T @ j @ c)
    check('basis_covariance', q(a, pulled), q(c @ a, j))
    reject('basis_metric_not_updated', q(c @ a, j), q(a, j))
    check('zero_amplitude', q(np.zeros((2, 1), complex), j), np.zeros((1, 1), complex))
    check('indefinite_null_nonzero', q(np.ones((2, 1), complex), np.diag([1, -1])), np.zeros((1, 1), complex))
    check('empty_amplitude_space', q(np.zeros((0, 1), complex), np.zeros((0, 0), complex)), np.zeros((1, 1), complex))
    # Mutate the original adjoint function, without changing canonical sources.
    class RemoveConjugation(ast.NodeTransformer):
        def visit_Call(self, node):
            node = self.generic_visit(node)
            if isinstance(node.func, ast.Attribute) and node.func.attr == 'conj' and not node.args:
                return node.func.value
            return node
    wrong = {'J': SimpleNamespace(multiply=jets['multiply'])}
    mutated = ast.fix_missing_locations(RemoveConjugation().visit(ast.parse(ast.unparse(response_tree))))
    exec(compile(mutated, '<isolated wrong native adjoint>', 'exec'), wrong)
    reject('mutated_original_adjoint_detected', pulled,
           wrong['gram']({(0, 0): c}, {(0, 0): j})[(0, 0)])
    # Read-only source anchors; these are not claimed as executed physical maps.
    open_metrics = next(n for n in ast.parse(PATHS['current'].read_text()).body
                        if isinstance(n, ast.FunctionDef) and n.name == 'open_metrics')
    text = ast.unparse(open_metrics)
    if 'g: s * a[' not in text or 'g: -s * a[' not in text:
        raise AssertionError('native outgoing/incident orientation anchor changed')
    if any(digest(BASE / p) != h for p, h in before.items()):
        raise RuntimeError('native source drift')
    result = {
        'status': 'PASS', 'source_sha256': before, 'instrument_sha256': digest(Path(__file__)),
        'selected_ast_sha256': extracted, 'checks': checks,
        'orientation_anchor': 'open_metrics: outgoing s*a, incoming -s*a; source inspection only',
        'scope': 'Grade-zero full complex current contractions and basis pullback, exact small fixtures.',
        'limits': ['Selected AST functions only; no production imports or saved physical operands.',
                   'No physical solve, current reality/positivity, numerical channel census or convergence certified.',
                   'Translation and NumPy arithmetic are outside the Lean kernel.',
                   'Full retained-grade bookkeeping belongs to the separate follow-on increment.'],
    }
    (M / 'S11_lean_flux_source_checks.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
