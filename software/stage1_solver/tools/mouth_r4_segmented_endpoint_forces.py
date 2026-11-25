"""Six signed forces on a held R4 state, with a local depth derivative on its material branch."""
import copy
import ctypes as ct
from pathlib import Path
import subprocess

import numpy as np

from mouth_cp10_label_gauge import depth_grid
from mouth_r4_aligned_crystal import AlignedGauge
from mouth_r4_aligned_derivatives import SignedGauge
from mouth_r4_eulerian_release import order_grid
from mouth_r4_full_label_precision import HPData, HPGrid, P
from mouth_r4_smooth_endpoint import LD, SHAPE_WORK_TOL, shape_work


def depth_stencil(directory, banks, labels, length, reference, active):
    """Select two consecutive widths by observed material branches before reading forces."""
    library = directory/'endpoint_branches.so'
    source = directory/'endpoint_branches.cpp'
    source.write_text('#include "audit.cpp"\n#include "mouth_r4_segmented_endpoint_forces.hpp"\n')
    subprocess.run(['g++','-std=gnu++17','-O3','-shared','-fPIC',
                    '-I',str(Path(__file__).resolve().parent),str(source),'-o',str(library)],check=True)
    count_changes = ct.CDLL(str(library)).endpoint_changed_rays
    count_changes.argtypes = [ct.POINTER(HPData), P, HPGrid, ct.c_void_p, P, ct.c_void_p]
    count_changes.restype = ct.c_int
    xi = np.ascontiguousarray(labels, LD)
    enabled = np.ascontiguousarray(active, np.int32)
    initial = depth_grid(length,reference)
    checks, consecutive = [], []
    for power in range(18,29):
        h = LD(2)**-power
        row = {'h':float(h)}
        for side, q in (('minus',-h), ('plus',h)):
            delta = np.ascontiguousarray(depth_grid(length*np.exp(q),reference)-initial,LD)
            first = np.full(5,-1,dtype=np.int32)
            n = count_changes(ct.byref(banks.data),xi.ctypes.data_as(P),banks.rule,enabled.ctypes.data,
                              delta.ctypes.data_as(P),first.ctypes.data)
            row[side] = dict(changed_rays=n,first_cell_and_ray=first.tolist() if n else None)
        checks.append(row)
        if row['minus']['changed_rays'] or row['plus']['changed_rays']:
            consecutive.clear()
        else:
            consecutive.append(float(h))
        if len(consecutive)==2:
            break
    return dict(selected_widths=consecutive,probes=checks,active_cells=int(enabled.sum()),
                selection='Two successive widths whose signed banks retain every active ray branch at its quadrature points.')


def endpoint_shape_work(driver, trial, solid, labels, fields, workdir):
    anchor = copy.copy(driver)
    anchor.folder, anchor.labels, anchor.fields = Path(workdir), labels, fields
    order = order_grid(trial,fields)
    energy, gradient, _, _ = solid.evaluate(labels,order)
    gauge = AlignedGauge(trial,driver.reference_depth,order)
    oracle = SignedGauge(trial,driver.reference_depth,order)
    rows = shape_work(anchor,trial,solid,order,energy,gradient,gauge=gauge,energy_gauge=oracle,
                      depth_stencil=depth_stencil)
    for i, row in enumerate(rows):
        widths = row['material_branches']['selected_widths'] if i==0 else [2.**-18,2.**-20]
        probes = [b for b in row['central'] if b['h'] in widths]
        row['passed'] = len(probes)==2 and all(abs(b['minus_analytic'])<=SHAPE_WORK_TOL for b in probes)
    return dict(six_total=[r['analytic'] for r in rows],rows=rows,tolerance=SHAPE_WORK_TOL,
                passed=len(rows)==6 and all(r['passed'] for r in rows))
