"""Independent time-variation check of the live raw mechanical assembly."""
from pathlib import Path
import sys
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import S11c_b_brane_operator_sympy_audit as b

t = sp.Symbol('inertiaCheckTime', real=True)
fields = tuple(sp.Function('inertiaCheckU' + str(i))(t) for i in range(3))
h = sp.Function('inertiaCheckE')(t)
jet_map = dict(zip((*b.u_tt, b.e_tt), (sp.diff(f, t, 2) for f in (*fields, h))))
for representative in b.DENSITY_REPS:
    density = b.density_pair(representative)[1]
    # Supplied T; independent functional time differentiation, not the helper.
    T = b.epsilon**2 * (density * sum(sp.diff(f, t)**2 for f in fields)
                        + b.mu_W * sp.diff(b.W_bg * h, t)**2) / 2
    action = sp.Tuple(*(sp.diff(sp.diff(T, sp.diff(f, t)), t) / b.epsilon
                        for f in (*fields, h)))
    raw, _ = b.operator_from_density(sp.S.Zero, representative)
    assembled = sp.Tuple(*b.named_tuple_row(raw['U_BODY_BALANCE'], 'EXPANDED'),
                         b.named_tuple_row(raw['E_W_BALANCE'], 'EXPANDED')).xreplace(jet_map)
    residual = sp.Tuple(*(sp.expand(a - c) for a, c in zip(assembled, action)))
    inverse = {v: k for k, v in jet_map.items()}
    objects = {name: obj.xreplace(inverse) for name, obj in
               (('ASSEMBLED', assembled), ('ACTION', action), ('RESIDUAL', residual))}
    dimensions = b.dimension_object(objects['ACTION'])
    print(representative, sp.srepr(b.case_payload(objects,
          {name: dimensions for name in objects})), flush=True)
