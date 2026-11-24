"""Independent physical readback of a point in the refined R4 collar.

Differentiate the original displacement field and its added vertex basis
directly. This reader uses neither the child determinant certificate nor the
native crystal integration kernel. No coordinate or field is advanced.
"""
from bisect import bisect_right
from fractions import Fraction
from itertools import product

import mpmath as mp
import numpy as np


def _order(body, order, point):
    if np.any(body.buffers['base']):
        raise ValueError('This reader assumes the zero-base continuous-order preparation')
    weights = []
    for key, xx in zip(('field_x', 'field_t', 'field_a'), point):
        nodes = [Fraction(float(v)) for v in body.buffers[key]]
        if not nodes[0] <= xx <= nodes[-1]:
            raise ValueError('Material readout outside the carrier grid')
        c = min(len(nodes)-2, bisect_right(nodes, xx)-1)
        h = (xx-nodes[c])/(nodes[c+1]-nodes[c])
        weights.append(((c,1-h), (c+1,h)))
    return sum(Fraction(float(order[body.ids[i,k,j]]))*wi*wj*wk
               for (i,wi),(j,wj),(k,wk) in product(*weights))


def _mp(v):
    return mp.mpf(v.numerator)/v.denominator if isinstance(v, Fraction) else mp.mpf(float(v))


def _collar_reference(body, x, t):
    """Flat registered crystal metric at the witness reference position."""
    pressure = _mp(body.data.pressure)
    e = -mp.log(x)/pressure
    er, tr = [list(map(_mp, body.buffers[k])) for k in ('reg_e','reg_t')]
    i = min(len(er)-2,max(0,bisect_right(er,e)-1))
    j = min(len(tr)-2,max(0,bisect_right(tr,t)-1))
    le,lt = er[i+1]-er[i], tr[j+1]-tr[j]
    u,v = (e-er[i])/le,(t-tr[j])/lt
    hi = int(u+v>1)
    sgn = 1-2*hi
    ij = ((j+hi,i+hi), (j+hi,i+1-hi), (j+1-hi,i+hi))
    labels = body.buffers['labels'].reshape(len(tr),len(er),2)
    a,b,c = [mp.matrix(list(map(_mp,labels[idx]))) for idx in ij]
    de,dt = (b-a)*sgn/le,(c-a)*sgn/lt
    w1,w2 = (u,v) if not hi else (1-u,1-v)
    radius = 4+a[0]+w1*(b[0]-a[0])+w2*(c[0]-a[0])
    return mp.matrix([[-de[0]/(pressure*x),dt[0],0],[-de[1]/(pressure*x),dt[1],0],
                      [0,0,radius]])


def measure(body, solid, surplus, order, witness, split_angle):
    """split_angle is the absolute material angle of the physical integration seam."""
    j,i,k = witness['cell_t_x_angle']
    tri = witness['triangle']
    if i < body.cut_node or k != 0:
        raise ValueError('Independent point readback covers the first collar angular prism')
    offsets = ((0,0),(1,0),(0,1)) if tri == 0 else ((1,1),(0,1),(1,0))
    sign = 1-2*tri
    q = lambda a: Fraction(float(a))
    axes = [[q(b) for b in arr[index:index+2]] for arr,index in
            zip((body.x,body.t,body.a),(i,j,k))]
    dx,dt,da = (a[1]-a[0] for a in axes)
    y = q(witness['angular_fraction'])
    s = (q(split_angle)-axes[2][0])/da
    if not 0 < s < 1 or y == s:
        raise ValueError('Specify a one-sided location off the angular derivative seam')
    hat, dz = (y/s, 1/(da*s)) if y < s else ((1-y)/(1-s), -1/(da*(1-s)))
    b1,b2 = (q(v) for v in witness['barycentric'][1:])
    bary = (1-b1-b2,b1,b2)
    bx,bt = (-sign/dx,sign/dx,0),(-sign/dt,0,sign/dt)
    pos = np.zeros(3,object)
    grad = np.zeros((3,3),object)
    ref = np.zeros(3,object)
    for vertex,(ix,jt) in enumerate(offsets):
        for side in (0,1):
            weight,slope = (y,1/da) if side else (1-y,-1/da)
            node = np.array([ax[n] for ax,n in zip(axes,(ix,jt,side))],dtype=object)
            mapped = node + [q(v) for v in solid[j+jt,i+ix,k+side]]
            ref += bary[vertex]*weight*node
            pos += bary[vertex]*weight*mapped
            grad += np.outer(mapped,[bx[vertex]*weight,bt[vertex]*weight,bary[vertex]*slope])
        if (j+jt,i+ix) == (5,18):
            extra = np.array([q(v) for v in surplus],dtype=object)
            pos += bary[vertex]*hat*extra
            grad += np.outer(extra,[bx[vertex]*hat,bt[vertex]*hat,bary[vertex]*dz])
    chi = _order(body,order,pos)
    with mp.workdps(70):
        x,t,a = map(_mp,ref)
        xx,tt,aa = map(_mp,pos)
        length = _mp(body.data.length)
        base, z_slope = [mp.mpf(float(v)) for v in body.sh.radius(float(length*tt))]
        h = (tt-_mp(.06))/_mp(.56)
        bump, derivative = mp.mpf(0), mp.mpf(0)
        if 0 < h < 1:
            bump = 16*_mp(body.epsilon)*h*h*(1-h)**2
            derivative = 32*_mp(body.epsilon)*h*(1-h)*(1-2*h)/_mp(.56)
        leg,dleg = (3*mp.cos(aa)**2-1)/2,-3*mp.cos(aa)*mp.sin(aa)
        radius = xx*base+xx**3*bump*leg
        curr = mp.matrix([[base+3*xx**2*bump*leg,xx*length*z_slope+xx**3*derivative*leg,
                           xx**3*bump*dleg],[0,length,0],[0,0,radius]])
        registered = _collar_reference(body,x,t)
        chart = mp.matrix([[_mp(v) for v in row] for row in grad])
        F = curr*chart*registered**-1
        hoop = (radius/registered[2,2])*(mp.sin(aa)/mp.sin(a) if a else chart[2,2])
        # Registered coordinates may reverse meridional orientation. The native
        # placement is positive when its chart is positive on this fixed branch.
        jac = mp.det(F)*hoop*mp.sign(mp.det(registered))
        energy = ((sum(f*f for f in F)+hoop**2-4)/2-mp.log(jac)+(jac-1)**2/4
                  if jac>0 else None)
        return dict(chart_J=float(mp.det(chart)),J=float(jac),chi=float(chi),
                    W=float(energy) if energy is not None else None,
                    chi_W=float(_mp(chi)*energy) if energy is not None else None,
                    placement_x_t_angle=[float(v) for v in pos],
                    upper_bound_ordered_J=float(jac) if chi>0 and jac>0 else None)


if __name__ == '__main__':
    import argparse
    import json
    from pathlib import Path
    from mouth_exact_checkpoint_control import acquired, digest

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_root',type=Path)
    parser.add_argument('snapshot',type=Path)
    parser.add_argument('witness_json',type=Path)
    parser.add_argument('witness_path',help='slash-separated object keys and list indices')
    args = parser.parse_args()
    root = args.run_root.resolve()
    build,*_ = acquired(root/'r4_melted_survival01',root/'r4_exact_support01/patch')
    witness = json.loads(args.witness_json.read_text())
    for key in args.witness_path.split('/'):
        witness = witness[int(key)] if isinstance(witness,list) else witness[key]
    with np.load(args.snapshot,allow_pickle=False) as saved:
        if int(saved['radius']) != 4 or int(saved['winding']) != 1:
            raise ValueError('Unexpected physical state')
        body = build(saved['x']).b
        material = saved['q']
        extra_dof = len(material)-body.dim
        if extra_dof not in (0,3):
            raise ValueError('Unexpected refined coordinate count')
        solid = body.unpack(material[:body.dim])
        if 'solid' in saved and not np.array_equal(saved['solid'],solid):
            raise ValueError('Inconsistent crystalline labels')
        surplus = material[body.dim:] if extra_dof else np.zeros(3)
        lo,hi = body.a[:2]
        edges = [a for a in body.buffers['angle_edges'] if lo < a < hi]
        seam = min(edges,key=lambda a:abs(a-0.5*(lo+hi)))
        row = measure(body,solid,surplus,saved['field'][2],witness,seam)
    print(json.dumps(dict(sha256=digest(args.snapshot),surplus=surplus.tolist(),seam=float(seam),
                          witness=witness,independent=row), allow_nan=False))
