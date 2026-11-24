"""Check the scalar melted-cell barrier and its 24 nodal derivatives at the R4 bottlenecks."""
import argparse
import ctypes as ct
from itertools import product
import json
from pathlib import Path

import mpmath as mp
import numpy as np

from mouth_cp10_label_grid import get_cp10
from mouth_exact_checkpoint_control import digest, read


def unit_corner_work(values, widths):
    """Exact signed coordinate volumes from the physical cell's stored Q1 edge coefficients."""
    value = mp.mpf(0)
    for corner in product((0, 1),repeat=3):
        gradient = mp.matrix(3,3)
        for axis, h in enumerate(widths):
            low = list(corner)
            high = list(corner)
            low[axis] = 0
            high[axis] = 1
            for c in range(3):
                gradient[c,axis] = (values[tuple(high)+(c,)]-values[tuple(low)+(c,)])/h
        determinant = mp.det(gradient)
        value += (determinant - 1 - mp.log(determinant))/8
    return value


def compare(scratch, lib, source):
    cpp = Path(__file__).with_name('mouth_r4_melt_barrier.cpp')
    if digest(cpp) != source:
        raise ValueError('Barrier source changed during the check')
    fun = ct.CDLL(str(lib.resolve())).melt_barrier
    fun.argtypes = [ct.c_int] + [ct.c_void_p]*7
    fun.restype = ct.c_int
    _, _, axes = get_cp10(scratch)
    output = dict(source_sha256=source, library_sha256=digest(lib), cases={})
    for name in ('primary','secondary'):
        p = scratch/'r4_eulerian_stage2_01'/name/'label_stall.npz'
        readout = read(scratch/'r4_eulerian_corner01'/(name+'.json'))
        if digest(p) != readout['state_sha256']:
            raise ValueError('Unrecognized R4 label map')
        cell = readout['accepted_limiting_cell']
        coord = np.array(cell['cell'])
        with np.load(p,allow_pickle=False) as f:
            a = np.ascontiguousarray(f['labels'][tuple(slice(i,i+2) for i in coord)])
        widths = np.array([ax[i+1] - ax[i] for ax,i in zip(axes,coord)])
        ids = np.arange(8,dtype=np.int32)
        w = np.array([.125])
        g = np.zeros_like(a)
        det = np.zeros(8)
        energy = np.zeros(1,np.longdouble)
        rc = fun(1,*[v.ctypes.data for v in (ids,widths,w,a,g,det,energy)])
        if rc: raise ArithmeticError('Positive saved cell rejected')
        with mp.workdps(75):
            exact = np.array([mp.mpf(float(t)) for t in a.flat], dtype=object).reshape(a.shape)
            hs = [mp.mpf(float(t)) for t in widths]
            scalar = unit_corner_work(exact,hs)
            independent = np.empty_like(a)
            step = mp.mpf('1e-45')
            for idx in np.ndindex(a.shape):
                z = exact.copy()
                z[idx] += step*1j
                independent[idx] = float(mp.im(unit_corner_work(z,hs))/step)
        component_error = np.abs(g-independent)/np.maximum(1.,np.abs(independent))
        output['cases'][name] = dict(cell=cell['cell'],saved_state_sha256=digest(p),
                        production_barrier_per_unit_volume=float(energy[0]),
                        independent=float(scalar), energy_difference=float(energy[0])-float(scalar),
                        production_min_corner_det=float(min(det)), max_relative_derivative_error=float(component_error.max()),
                        max_gradient=float(np.abs(independent).max()), largest_error_component=np.unravel_index(component_error.argmax(),a.shape))
        output['cases'][name]['largest_error_component'] = list(map(int,output['cases'][name]['largest_error_component']))
    return output


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('scratch',type=Path)
    p.add_argument('library',type=Path)
    p.add_argument('source_sha256')
    p.add_argument('output',type=Path)
    v=p.parse_args()
    result = compare(v.scratch,v.library,v.source_sha256)
    with v.output.open('x') as out:
        json.dump(result,out,indent=2,allow_nan=False)
        out.write('\n')
    print(json.dumps(result),flush=True)
