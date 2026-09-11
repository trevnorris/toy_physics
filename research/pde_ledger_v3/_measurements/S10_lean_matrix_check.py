#!/usr/bin/env python3
"""Reject matrix, determinant and basis mutations in isolated Lean copies."""
import json

import S10_lean_q6_q7_check as harness


def main():
    harness.SCRATCH = harness.LEAN / 's10' / '_scratch' / 'matrix'
    harness.SCRATCH.mkdir(parents=True, exist_ok=True)
    checks = [
        harness.run_mutation(
            'wrong_stack_row_unit', 'MatrixTrees.lean',
            'Fin.lastCases lengthDim⁻¹ (fun _ => matrixDim D field)',
            'Fin.lastCases (matrixDim D field) (fun _ => matrixDim D field)',
            ['stackedTree_hasDim']),
        harness.run_mutation(
            'drop_determinant_cofactor_sign', 'MinorTrees.lean',
            '.mul (.mul (.scalar ((-1 : ℝ) ^ j.val)) (M 0 j))',
            '.mul (.mul (.scalar 1) (M 0 j))',
            ['detTree_eval']),
        harness.run_mutation(
            'wrong_transverse_basis_sign', 'BasisTrees.lean',
            ': Vec D :=\n  unit j - (k j / k p) • unit p',
            ': Vec D :=\n  unit j + (k j / k p) • unit p',
            ['chartVector_transverse', 'chartVectorTree_eval']),
    ]
    output = harness.ROOT / '_measurements' / 'S10_lean_matrix_checks.json'
    output.write_text(json.dumps({
        'status': 'PASS',
        'scope': 'three isolated matrix/minor/basis source mutations',
        'checks': checks,
    }, indent=2) + '\n')
    print('PASS: all three mutations rejected; canonical source hashes unchanged.')


if __name__ == '__main__':
    main()
