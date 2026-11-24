"""Sample physical inverse singular values of the frozen checkpoint-10 two-metric label map."""
import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from mouth_cp10_label_grid import get_cp10
from mouth_cp10_label_orientation import basis_and_grad, native_state
from mouth_exact_checkpoint_control import digest
from mouth_supercore_shape import preparation


def physical_frame(body, points):
    x, t, a = np.moveaxis(points, -1, 0)
    length = body.sh.length
    radius, slope = body.sh.radius(t * length)
    y = (t - .06) / .56
    interior = (y > 0) & (y < 1)
    band = np.where(interior, 16 * body.epsilon * y**2 * (1-y)**2, 0)
    dband = np.where(interior, 32 * body.epsilon * y * (1-y) * (1-2*y) / .56, 0)
    p = (3*np.cos(a)**2 - 1) / 2
    r = x * radius + x**3 * band * p
    frame = np.zeros(points.shape + (3,))
    frame[...,0,0] = radius + 3*x**2*band*p
    frame[...,0,1] = x*length*slope + x**3*dband*p
    frame[...,0,2] = -3*x**3*band*np.sin(a)*np.cos(a)
    frame[...,1,1] = length
    frame[...,2,2] = r
    return frame, r * np.sin(a)


def reference_frame(body, ref, points):
    x, t, a = np.moveaxis(points, -1, 0)
    radius, slope = ref.radius(t * ref.length)
    scale = body.data.core_ratio
    frame = np.zeros(points.shape + (3,))
    frame[...,0,0] = radius * scale
    frame[...,0,1] = x * ref.length * slope * scale
    frame[...,1,1] = ref.length
    frame[...,2,2] = x * radius * scale
    collar = x >= body.cut
    xx, tt = x[collar], t[collar]
    tr, er = body.buffers['reg_t'], body.buffers['reg_e']
    e = -np.log(xx)/body.data.pressure
    j = np.clip(np.searchsorted(tr, tt, side='right') - 1, 0, len(tr)-2)
    i = np.clip(np.searchsorted(er, e, side='right') - 1, 0, len(er)-2)
    dt, de = np.diff(tr)[j], np.diff(er)[i]
    u, v = (e - er[i])/de, (tt - tr[j])/dt
    hi = (u + v > 1).astype(int)
    sc = 1 - 2*hi
    labels = body.buffers['labels'].reshape(len(tr), len(er), 2)
    c = labels[j+hi, i+hi]
    dc1 = labels[j+hi, i+1-hi] - c
    dc2 = labels[j+1-hi, i+hi] - c
    f = np.zeros((len(xx), 3, 3))
    f[:,:2,0] = dc1 * (-sc/(de*body.data.pressure*xx))[:,None]
    f[:,:2,1] = dc2 * (sc/dt)[:,None]
    f[:,2,2] = 4 + c[:,0] + np.where(hi, 1-u, u)*dc1[:,0] + np.where(hi, 1-v, v)*dc2[:,0]
    frame[collar] = f
    return frame, frame[...,2,2]*np.sin(a)


def readback(root, out):
    body, _, axes = get_cp10(root)
    ref = preparation(4)[0]
    path = out / 'labels.npz'
    with np.load(path, allow_pickle=False) as f:
        labels, order, melt = f['labels'], f['order'], f['melted_cells']
    regions = {key: {'samples': 0} for key in ('all_active', 'positive_chi_below_1e-6', 'exact_melt')}
    for gp in product((.5-.5/np.sqrt(3), .5+.5/np.sqrt(3)), repeat=3):
        val, grad = basis_and_grad(labels, axes, gp)
        chi = basis_and_grad(np.broadcast_to(order[...,None], labels.shape), axes, gp)[0][...,0]
        if np.any((chi > 0) != ~melt):
            raise ValueError('Positive quadrature order and active-cell mask disagree')
        p = np.stack(np.broadcast_arrays(*np.ix_(*[
            a[:-1]+g*np.diff(a) for a,g in zip(axes,gp)])), axis=-1)
        B, hoop_phys = physical_frame(body, p)
        C, hoop_ref = reference_frame(body, ref, val)
        inverse = C @ grad @ np.linalg.inv(B)
        hoop = hoop_ref/hoop_phys
        sv3 = np.linalg.svd(inverse, compute_uv=False)
        sv = np.concatenate((sv3, np.abs(hoop[...,None])), axis=-1)
        jinv = np.sign(np.linalg.det(grad))*np.prod(sv, axis=-1)
        condition = sv.max(axis=-1)/sv.min(axis=-1)
        for name, mask in [('all_active', ~melt), ('positive_chi_below_1e-6',(chi>0)&(chi<1e-6)),('exact_melt',melt)]:
            row = regions[name]
            row['samples'] += int(mask.sum())
            if not mask.any():
                continue
            for key, field, sign in [('min_Jinv', jinv, 1), ('max_Jinv', jinv, -1), ('max_condition', condition, -1)]:
                idx = np.unravel_index(np.argmin(np.where(mask, sign*field, np.inf)), mask.shape)
                if key in row and sign*row[key]['value'] <= sign*field[idx]:
                    continue
                independent = native_state(body, val[idx], p[idx], grad[idx])
                if not np.isclose(independent['J_inverse'], jinv[idx], rtol=1e-12, atol=0):
                    raise ValueError('Independent physical determinant mismatch at witness')
                J = 1/jinv[idx]
                W = (np.sum(1/sv[idx]**2)-4)/2 - np.log(J) + (J-1)**2/4
                if not np.isclose(independent['W'],W,rtol=1e-12,atol=1e-13):
                    raise ValueError('Independent physical strain energy mismatch at witness')
                row[key] = dict(value=float(field[idx]), chi=float(chi[idx]), Jinv=float(jinv[idx]),
                    condition=float(condition[idx]), principal_inverse_3d=sv3[idx].tolist(),
                    hoop_inverse=float(hoop[idx]), cell=list(map(int,idx)), chart=p[idx].tolist(),
                    gauss=gp, side='core' if val[idx][0]<body.cut else 'collar')
    return dict(label_sha256=digest(path), sampling='Two Gauss points per axis in every Q1 cell; excludes cell faces and axes. The extrema are sampled, not enclosures.',regions=regions)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('out', type=Path)
    args = parser.parse_args()
    result = readback(args.root, args.out)
    (args.out/'condition.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result,allow_nan=False))
