#!/usr/bin/env python3
"""Fixed reconcile manifest: K_normal, K_source_route, K_operand_swap.
The normal wrapper is confined to build_material_carrier and restores both
builder and cache in finally. Each FORM retains identity and same-site X2.
"""
from pathlib import Path
import sys
import shutil
sys.dont_write_bytecode = True
import S11c_c2_N6_covariance_ablation_harness as h


def keys(*names):
    return [('S11CC2_N6RC_' + s, None) for s in names]

NORMAL = "    rows, _, provenance = n.face_factory(a, b, inputs, alpha, rho, 'MATERIAL', mu_slot)"
WRAPPER = '''    original_builder = a.build_material_face_source
    original_cache = dict(a._FACE_CACHE)
    def mixed_material_face_source(*args, **kwargs):
        result = original_builder(*args, **kwargs)
        components = list(result.normal_exact)
        components[0] = components[0] + a.grad_W[1]
        length = sp.sqrt(a.dot(tuple(components), tuple(components)))
        return n.replace(result, normal_exact=tuple(value / length for value in components))
    a.build_material_face_source = mixed_material_face_source
    a._FACE_CACHE.clear()
    try:
        rows, _, provenance = n.face_factory(a, b, inputs, alpha, rho, 'MATERIAL', mu_slot)
    finally:
        a.build_material_face_source = original_builder
        a._FACE_CACHE.clear()
        a._FACE_CACHE.update(original_cache)'''
SOURCE = '    source = closed_response(comp, inputs, m_coeff, ds, kernels)'
DEAD_CARRIER = keys('CARRIER_EULERIAN', 'CARRIER_MATERIAL', 'CARRIER_BRIDGE_RESIDUAL')
CERTIFIED = keys('R_N6', 'SPLIT_CHECK', 'CARRIER_BRIDGE_RESIDUAL', 'SOURCE_BRIDGE_RESIDUAL',
                 'CARRIER_CHANNEL', 'SOURCE_CHANNEL', 'CROSS_CHANNEL', 'CARRIER_EULERIAN',
                 'CARRIER_MATERIAL', 'SOURCE_EULERIAN', 'SOURCE_MATERIAL', 'EULERIAN_OPERAND', 'MATERIAL_OPERAND')
KNIVES = [
    h.knife('K_normal', 'build_material_carrier', NORMAL, WRAPPER,
            WRAPPER.replace('components[0] = components[0] + a.grad_W[1]', 'components[0] = 2 * components[0]'),
            keys('SOURCE_EULERIAN', 'SOURCE_MATERIAL', 'SOURCE_BRIDGE_RESIDUAL', 'EULERIAN_OPERAND')),
    h.knife('K_source_route', 'run', SOURCE,
            '    source, _ = n.build_increment(comp, inputs, m_coeff, ds, kernels, slots)',
            SOURCE.replace(', ds,', ', {s: {w: 2 * value for w, value in terms.items()} for s, terms in ds.items()},'), DEAD_CARRIER),
    h.knife('K_operand_swap', 'run', 'ms[s].get(w, sp.S.Zero)', 'es[s].get(w, sp.S.Zero)',
            '2 * ms[s].get(w, sp.S.Zero)', DEAD_CARRIER),
]
CONFIG = dict(engine='S11c_c2_N6_reconcile_sympy.py', seed=110602,
              certified=CERTIFIED, knives=KNIVES, extra=[])
# run() hashes this declared provenance input in addition to the Python
# siblings and route specification. Copy it without changing its contents.
copy_python_tree = h.copy_tree


def copy_reconcile_tree(tree):
    copy_python_tree(tree)
    (tree / 'directives').mkdir()
    path = h.ROOT / 'directives/S11c_c2_SHARED_PHYSICS.md'
    shutil.copy2(path, tree / 'directives' / path.name)


if __name__ == '__main__':
    h.copy_tree = copy_reconcile_tree
    sys.exit(h.main(CONFIG, Path(__file__).resolve()))
