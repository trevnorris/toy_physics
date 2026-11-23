"""Signed shape controls for the receiving-order stationarity at radius four."""
import argparse
from pathlib import Path

import numpy as np

from mouth_internal_collar import CollarPair
from mouth_joined_front_audit import Front
from mouth_poly_shallow_release import CUT, exterior
from mouth_receiver_order_equilibrium import receiver
from mouth_sealed_shape_relax import Meridian
from mouth_supercore_front import check_guard
from mouth_supercore_joined_front import emit, file_sha


def side_source(base, shape):
    pair = CollarPair(shape, base.co, base.na, base.nt)
    n = base.order_mesh
    t = np.linspace(0, 1, 2*n+1)
    x = np.linspace(0, 1, n+1)
    return pair.evaluate(base.labels, (t, x, base.c),
                         ((1-np.cos(np.pi*t))/2, .1*x+.9*x*x, base.outer))['outer_source']


def main(path):
    check_guard(path)
    base = Front(path, (13, 10, 14), 1)
    emit(kind='input', state=file_sha(path))
    for shift in (0., -.005, .005, -.05, .05):
        p = base.sh.params.copy()
        p[3] += shift
        shape = Meridian(4., p, extra=base.sh.extra)
        load = side_source(base, shape)
        source = lambda wa, xo: load
        profile, row = receiver(shape, base.order_mesh, CUT, source)
        try:
            old = exterior(shape, base.order_mesh, CUT, side_order_work=source)
            old_energy = old['energy']+old['static_work']
        except ArithmeticError:
            old_energy = None
        emit(kind='signed_control', shape_shift=shift, old_energy=old_energy,
             new_energy=row['energy']+row['static_work'], **row)
        if old_energy is not None and abs(old_energy-row['energy']-row['static_work']) > 1e-7:
            raise ArithmeticError("Receiver energy differs on matched chart")
        if shift == 0:
            delta = 2e-4
            def reduced(factor):
                _, result = receiver(shape, base.order_mesh, CUT, lambda wa, xo: (1+factor)*load)
                return result['energy'] + result['static_work']
            work = (reduced(delta)-reduced(-delta))/(2*delta)
            emit(kind='side_virtual_work', work=work, native_work=row['static_work'])
            if abs(work-row['static_work']) > 2e-7:
                raise ArithmeticError("Side order work fails load variation")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('state', type=Path)
    main(parser.parse_args().state)
