"""Compare stored energy and five fixed-knot shape forces on one unfolded R4 control."""
import argparse
import ctypes as ct
import json
from pathlib import Path
import subprocess

import numpy as np

from mouth_cp10_nested_forces import subdivide
from mouth_cp10_spatial_crystal import MODES, signed_geometry
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_free_shape_solid import Data, Grid, Pointer, pointer


def run(root, area, output, degree):
    from mouth_r4_unfolded_control import load_control_context

    output.mkdir(parents=True, exist_ok=False)
    control = area / 'control.npz'
    fingerprint = digest(control)
    ctx = load_control_context(root, area)
    body = ctx[2]
    with np.load(control) as state:
        solid = state['solid']
        full_order = state['field'][2]
    q = body.pack(solid)
    native_order = full_order[body.ids].transpose(0, 2, 1)
    control_inputs = []
    arrays, spatial_axes = [], []
    for level in (1, 2, 4):
        path = area / str(level) / 'labels.npz'
        with np.load(path) as saved:
            if str(saved['source_sha256']) != fingerprint:
                raise ValueError('Label provenance differs from the unfolded control')
            order = saved['order']
            if not np.array_equal(order[::level, ::level, ::level], native_order):
                raise ValueError('The native order field has changed')
            axes = [np.ascontiguousarray(saved[c]) for c in 'xta']
            for coord, axis in zip('xta', axes):
                if not np.allclose(axis, subdivide(body.buffers['field_'+coord], level), atol=4e-16, rtol=0):
                    raise ValueError('Expected native chart knots and their subdivisions')
            spatial_axes.extend(axes)
            arrays.append(np.ascontiguousarray(saved['labels']))
        control_inputs.append(digest(path))

    src = Path(__file__).with_suffix('.cpp')
    library = output / 'shape_work.so'
    subprocess.run(['g++', '-O3', '-std=c++17', '-shared', '-fPIC', str(src), '-o', str(library)], check=True)
    fn = ct.CDLL(str(library)).unfolded_shape_work
    fn.argtypes = [ct.POINTER(Data), ct.c_void_p, ct.POINTER(Grid), ct.POINTER(Pointer),
                   ct.POINTER(ct.c_ubyte), ct.POINTER(Data), Pointer, Pointer, ct.POINTER(ct.c_longlong)]
    fn.restype = ct.c_int
    length = body.data.length
    knobs = [signed_geometry(ctx, mode, sign*1e-5) for mode in range(1, 6) for sign in (-1, 1)]
    for knob in knobs:
        if knob.data.length != length or not np.array_equal(knob.buffers['field_t'], body.buffers['field_t']):
            raise ValueError('An order knot moves in a nominally fixed-knot mode')
    probes = (Data*len(knobs))(*(knob.data for knob in knobs))
    axes_c = (Grid*9)(*(Grid(len(axis), pointer(axis)) for axis in spatial_axes))
    labels_c = (Pointer*3)(*(pointer(a) for a in arrays))
    # The mapped axial coordinate is barycentric on a triangle and linear axially
    # between its two polar faces. Its extrema on the prism occur at vertices.
    last = np.where(np.any(native_order > 0, axis=(0, 2)))[0][-1]
    t = body.buffers['field_t']
    landed = body.t[:,None,None]+solid[...,1]
    lo = np.minimum.reduce([landed[j:j+len(body.t)-1, i:i+len(body.x)-1, k:k+len(body.a)-1]
                            for j in (0,1) for i in (0,1) for k in (0,1)])
    mask = np.ascontiguousarray(lo <= t[min(last+1, len(t)-1)], dtype=np.uint8)
    points, weights = np.polynomial.legendre.leggauss(degree)
    points, weights = (points+1)/2, weights/2
    body.data.gauss = Grid(degree, pointer(points))
    body.data.weights = pointer(weights)
    rows = np.zeros((4, 6))
    count = ct.c_longlong()
    fd = np.zeros(1)
    body.work = lambda d, out: fn(d, out, axes_c, labels_c, mask.ctypes.data_as(ct.POINTER(ct.c_ubyte)),
                                  probes, pointer(rows), pointer(fd), ct.byref(count))
    body.evaluate(q, full_order, load=False, gradient=False)
    if not np.isfinite(rows).all() or np.any(rows[:,0] <= 0) or count.value <= 0 or fd[0] > 2e-7:
        raise ArithmeticError(f'Invalid force or signed-geometry derivative: {rows}, {fd}, {count.value}')

    prior_file = src.parents[3] / 'research/pde_ledger_v3/_measurements/mouth_r4_two_metric_labels.json'
    prior = read(prior_file)['data']['total_shape_forces']
    # Exterior material and bath fields are common to both representations. Their
    # absolute CP10 forces can be used only with identical geometry and bath fields.
    cp = root / 'r4_melted_survival01/coarse/state_00010.npz'
    if digest(cp) != CHECKPOINT_SHA:
        raise ValueError('Unexpected parent checkpoint')
    with np.load(cp) as parent, np.load(control) as state:
        for k in ('field', 'shape_params', 'shape_extra','x'):
            if not np.array_equal(parent[k], state[k]):
                raise ValueError('Common force changed for '+k)
    common = np.asarray(prior['field_analytic'])[1:]
    total = common + rows[:,1:]
    row = dict(parent_sha256=CHECKPOINT_SHA, control_sha256=fingerprint, label_sha256=control_inputs,
               cpp_sha256=digest(src), common_force_source_sha256=digest(prior_file), degree=degree,
               crystal_columns=['energy', *MODES[1:]], rows=['material','spatial 1x','spatial 2x','spatial 4x'],
               crystal_values=rows.tolist(), common_force=common.tolist(), total_force=total.tolist(),
               spatial_minus_material=(rows[1:] - rows[0]).tolist(), ordered_sample_count=count.value,
               sampled_signed_derivative_relative_error=float(fd[0]), signed_check_step=1e-5)
    with (output/'result.json').open('x') as stream:
        json.dump(row, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(row, allow_nan=False), flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root', type=Path)
    p.add_argument('area', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('degree', type=int, choices=(6,8))
    a=p.parse_args()
    run(a.root.resolve(),a.area.resolve(),a.output.resolve(),a.degree)
