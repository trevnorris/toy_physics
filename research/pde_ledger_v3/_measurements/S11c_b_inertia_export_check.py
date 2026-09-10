"""Compare regenerated slab rows with independent kinetic action and saved rows."""
from pathlib import Path
import hashlib
import sys
import sympy as sp
from sympy.core.symbol import Str

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import S11c_b_brane_operator_sympy_audit as b
from ledger_fold import _restore
from S11c_inertia_artifact_audit import ROOT, export_data

BASELINE = Path('/tmp/s11c-inertia-repair-baseline')
old_values, _, _ = export_data(BASELINE / 'scripts/S11c_b_exports.py')
new_values, _, _ = export_data(ROOT / 'scripts/S11c_b_exports.py')
old = dict(_restore(old_values['slab_operator']))
new = dict(_restore(new_values['slab_operator']))
t = sp.Symbol('inertiaExportTime', real=True)
fields = tuple(sp.Function('inertiaExportU' + str(i))(t) for i in range(3))
h = sp.Function('inertiaExportE')(t)
accelerations = (*b.u_tt, b.e_tt)
inverse = dict(zip((sp.diff(f, t, 2) for f in (*fields, h)), accelerations))

for case, payload in new.items():
    previous = b.named_tuple_row(old[case], 'VALUE')
    current = b.named_tuple_row(payload, 'VALUE')
    representative = str(case[1])
    density = b.density_pair(representative)[1]
    T = b.epsilon**2 * (density * sum(sp.diff(f, t)**2 for f in fields)
                        + b.mu_W * sp.diff(b.W_bg * h, t)**2) / 2
    action = sp.Tuple(*(b.first_shape_series_reference(
        (sp.diff(sp.diff(T, sp.diff(f, t)), t) / b.epsilon).xreplace(inverse))
        for f in (*fields, h)))

    def mechanical(body):
        return sp.Tuple(*b.named_tuple_row(b.named_tuple_row(body, 'U_BODY_BALANCE'), 'EXPANDED'),
                         b.named_tuple_row(b.named_tuple_row(body, 'E_W_BALANCE'), 'EXPANDED'))

    before, after = mechanical(previous), mechanical(current)
    delta = sp.Tuple(*(sp.expand(a - p) for a, p in zip(after, before)))
    inertia = sp.Tuple(*(sp.expand(sum(sp.diff(row, acc) * acc for acc in accelerations))
                         for row in after))
    action_residual = sp.Tuple(*(sp.expand(a - c) for a, c in zip(inertia, action)))
    zero_acceleration = dict.fromkeys(accelerations, sp.S.Zero)
    nonkinetic_before = before.xreplace(zero_acceleration)
    nonkinetic_after = after.xreplace(zero_acceleration)
    nonkinetic_residual = sp.Tuple(*(sp.expand(a - p)
                                  for a, p in zip(nonkinetic_after, nonkinetic_before)))
    objects = {'ROW_BEFORE_OPERAND': before, 'ROW_AFTER_OPERAND': after,
               'ROW_DELTA': delta, 'INERTIA_OPERAND': inertia, 'ACTION_OPERAND': action,
               'INERTIA_ACTION_RESIDUAL': action_residual,
               'NONKINETIC_BEFORE_OPERAND': nonkinetic_before,
               'NONKINETIC_AFTER_OPERAND': nonkinetic_after,
               'NONKINETIC_ROW_RESIDUAL': nonkinetic_residual}
    dimensions = b.dimension_object(action)
    print('CASE', sp.srepr(case), sp.srepr(b.case_payload(objects,
          {name: dimensions for name in objects})), flush=True)
    print('RESIDUALS', sp.srepr(b.case_payload(
        {case: sp.Tuple(action_residual, nonkinetic_residual)},
        {case: sp.Tuple(dimensions, dimensions)})), flush=True)
    untouched = []
    def identity_operands(left, right):
        return sp.Tuple(
            Str(hashlib.sha256(sp.srepr(left).encode()).hexdigest()),
            Str(hashlib.sha256(sp.srepr(right).encode()).hexdigest()),
            sp.sympify(left == right))
    for key, value in current:
        other = b.named_tuple_row(previous, str(key))
        if str(key) in ('U_BODY_BALANCE', 'E_W_BALANCE'):
            untouched.extend(sp.Tuple(key, slot, identity_operands(obj, b.named_tuple_row(other, str(slot))))
                             for slot, obj in value if str(slot) != 'EXPANDED')
        else:
            untouched.append(sp.Tuple(key, identity_operands(value, other)))
    print('OTHER_SLOTS_EQUAL', sp.srepr(sp.Tuple(case, sp.Tuple(*untouched))), flush=True)
