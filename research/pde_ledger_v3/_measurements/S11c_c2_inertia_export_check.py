"""Compute saved/current c2 row differences and their zero-acceleration operands."""
from pathlib import Path
import hashlib
import json
import sys
import sympy as sp
from sympy.core.symbol import Str

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from ledger_fold import _restore
from S11c_c2_selfenergy_fold_sympy_audit import cas, grades
from S11c_inertia_artifact_audit import export_data

OLD = Path('/tmp/s11c-inertia-repair-baseline/scripts/S11c_c2_exports.py')
NEW = ROOT / 'scripts/S11c_c2_exports.py'
old_values, _, _ = export_data(OLD)
new_values, _, _ = export_data(NEW)


def named(obj, name):
    return next(value for key, value in obj if str(key) == name)


def leaves(value, units, path=()):
    if isinstance(units, sp.MatrixBase):
        yield path, value, units
        return
    association = all(isinstance(item, sp.Tuple) and len(item) == 2
                      and isinstance(item[0], Str) for item in value)
    for index, (item, unit) in enumerate(zip(value, units)):
        if association:
            yield from leaves(item[1], unit[1], path + (str(item[0]),))
        else:
            yield from leaves(item, unit, path + (index,))


def digest(value):
    return hashlib.sha256(sp.srepr(value).encode()).hexdigest()


print('INPUT_DIGESTS', json.dumps({str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in (OLD, NEW)}), flush=True)
summary = []
for root in ('s11cc2ClosedSlabOperator', 's11cc2ClosedCouplingKernel'):
    old, new = dict(_restore(old_values[root])), dict(_restore(new_values[root]))
    for case, payload in new.items():
        before_payload = old[case]
        before_leaves = {path: expr for path, expr, unit in
                         leaves(named(before_payload, 'VALUE'), named(before_payload, 'DIMENSION_L_T_M'))}
        body = named(payload, 'VALUE')
        atoms = body.atoms(sp.Symbol)
        knobs = tuple(next(a for a in atoms if a.name == name)
                      for name in ('epsilon_shape', 'eta_bg', 'sigma_W'))
        for path, after, unit in leaves(body, named(payload, 'DIMENSION_L_T_M')):
            before = before_leaves[path]
            delta = sp.expand(after - before, deep=False)
            accelerations = {a for a in (before.atoms(sp.Derivative) | after.atoms(sp.Derivative))
                             if sum(n for v, n in a.variable_count if str(v) == 's11cc2Time') >= 2}
            substitution = dict.fromkeys(accelerations, sp.S.Zero)
            before_static, after_static = before.xreplace(substitution), after.xreplace(substitution)
            residual = sp.expand(after_static - before_static, deep=False)

            def physical(value):
                return {'VALUE': value, 'MULTIGRADE': sorted(grades(value, *knobs)),
                        'DIMENSION_L_T_M': unit}

            def operand(value):
                return {'SHA256': digest(value), 'MULTIGRADE': sorted(grades(value, *knobs)),
                        'DIMENSION_L_T_M': unit}

            record = {'ROOT': root, 'CASE': case, 'PATH': path,
                      'BEFORE_OPERAND': operand(before), 'AFTER_OPERAND': operand(after),
                      'ROW_DELTA': physical(delta),
                      'ZERO_ACCELERATION_BEFORE_OPERAND': operand(before_static),
                      'ZERO_ACCELERATION_AFTER_OPERAND': operand(after_static),
                      'ZERO_ACCELERATION_RESIDUAL': physical(residual)}
            print('C2_INERTIA_RECORD', sp.srepr(cas(record)), flush=True)
            summary.append({'root': root, 'case': list(map(str, case)), 'path': list(path),
                            'delta_is_zero': delta == 0, 'zero_acceleration_residual': sp.sstr(residual)})
print('SUMMARY', json.dumps(summary), flush=True)
