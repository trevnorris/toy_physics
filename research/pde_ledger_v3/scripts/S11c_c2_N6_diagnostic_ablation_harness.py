#!/usr/bin/env python3
"""Fixed diagnostic manifest: K_EW_rowdrop, K_slotdrop and K_ewsign.
FORM/IDENTITY/X2 variants use only the directive's literal construction sites.
Native TILT and N4_ADVECTION controls retain BASE, CORRUPTED, RESIDUAL order.
"""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
import S11c_c2_N6_covariance_ablation_harness as h


def keys(*names):
    return [('S11CC2_' + name, None) for name in names]

ROW = "    rows=flatten({'U':folded['U'],'E_W':folded['E_W'],'THETA':mass-correction})"
SLOTS = "    slots=tuple(inputs.a(prefix+label) for label in ('plus','minus') for prefix in ('delta_p_','d_w_delta_p_'))"
MU = keys('N6_MU_RECONSTRUCTION_IMPORTED', 'N6_MU_RECONSTRUCTION_NATIVE', 'N6_MU_RECONSTRUCTION_RESIDUAL')
CERTIFIED = keys('REP_INVARIANCE_EULERIAN_OPERAND', 'REP_INVARIANCE_MATERIAL_OPERAND', 'REP_INVARIANCE_RESIDUAL') + [
    ('S11CC2_N6_' + name + '_GUARD_RESIDUAL', route)
    for name in ('SLOT', 'CLOSURE') for route in ('EULERIAN', 'MATERIAL')] + [
    ('S11CC2_CONTROL_INDEPENDENCE_' + part, probe)
    for probe in ('TILT', 'N4_ADVECTION') for part in ('BASE', 'CORRUPTED', 'RESIDUAL')]
KNIVES = [
    h.knife('K_EW_rowdrop', 'face_factory', ROW, ROW + "\n    rows.pop('E_W')",
            ROW.replace("'E_W':folded['E_W']", "'E_W':2*folded['E_W']"),
            keys('REP_INVARIANCE_EULERIAN_OPERAND') + MU),
    h.knife('K_slotdrop', 'run', SLOTS,
            '''    slots=tuple(inputs.a(name) for name in
                ('d_w_delta_p_plus','delta_p_minus','d_w_delta_p_minus'))''',
            SLOTS + '''
    native_pressure_coefficients=globals()['pressure_coefficients']
    def pressure_coefficients(rows,active_slots):
        table=native_pressure_coefficients(rows,active_slots)
        return {key:(2*value if key[1].name=='delta_p_plus' else value)
                for key,value in table.items()}''',
            MU + [('S11CC2_N6_' + name, route) for name in ('MU_AMPLITUDE', 'FACE_VELOCITY')
                  for route in ('EULERIAN', 'MATERIAL')]),
]
CONFIG = dict(engine='S11c_c2_N6_diagnostic_sympy.py', seed=110602,
              certified=CERTIFIED, knives=KNIVES, extra=[h.knife('K_ewsign', 'face_factory', ROW,
              ROW.replace("'E_W':folded['E_W']", "'E_W':-folded['E_W']"), None, [])])
if __name__ == '__main__':
    sys.exit(h.main(CONFIG, Path(__file__).resolve()))
