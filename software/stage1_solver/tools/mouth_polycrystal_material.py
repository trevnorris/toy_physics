"""Shallow registered material on both sides, then a free internal receiver displacement."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import legval
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize

from mouth_finite_bore_support import Side, fields, graph_values
from mouth_sealed_shape_relax import Meridian
from mouth_unified_puncture import archived, RECORD, ROOT
from mouth_waist_outer_order import exterior


def stable_energy(trace, jac):
    return .5*(trace-4)-np.log(jac)+.25*(jac-1)**2


def bump(x, lower, upper):
    """C1 compact bump; derivative with respect to x."""
    t = (x-lower)/(upper-lower)
    on = (0 < t) & (t < 1)
    return (np.where(on, 16*(t*(1-t))**2, 0.),
            np.where(on, 32*t*(1-t)*(1-2*t)/(upper-lower), 0.))


def stationary_graph(side, shape):
    L = shape.length

    def objective(c):
        log_p, derivative = side.basis@c, side.deriv@c
        if np.max(log_p)>5 or np.min(log_p)<-32:
            return 1e70, np.zeros_like(c)
        r, rs, rz, h, jac, ap = graph_values(side.a, side.t, shape, log_p, derivative)
        trace = rs**2+rz**2+L**2+2*h**2
        q = -1+.5*jac*(jac-1)
        g = (rs**2*(1-ap) -(rz**2+2*h**2)*ap +rz*r*ap*derivative
             +q*(1-3*ap))
        gt = rz*r*ap
        return side.weight@stable_energy(trace,jac), (side.basis.T@(side.weight*g)
                                                       +side.deriv.T@(side.weight*gt))
    sol = minimize(objective, np.zeros(side.nt), jac=True, method='BFGS',
                   options={'gtol':2.e-8,'maxiter':3000})
    value, grad = objective(sol.x)
    if not np.isfinite(value) or np.max(abs(grad))>3.e-5:
        raise ArithmeticError(f'Unresolved normal-spacing amplitudes {value} {np.max(abs(grad))}')
    return sol.x, float(value), float(np.max(abs(grad)))


def free_receiver(shape, side, c, outer):
    """Geometry and all port traces remain fixed. Only the crystal label positions change."""
    w_nodes, x_nodes, outer_chi = outer
    current_order = RegularGridInterpolator((w_nodes,x_nodes),outer_chi)
    a,t = side.a,side.t
    log_p, dlog = side.basis@c, side.deriv@c
    r,rs,rz,h,j0,ap = graph_values(a,t,shape,log_p,dlog)
    z, zt = bump(t,.15,.85)
    b, ba = bump(-2*a,.04,.48)
    g,ga,gt = b*z,-2*ba*z,b*zt
    active = (g!=0)
    vol = side.weight[active]/side.chi[active]
    if not np.all(side.chi[active]>0):
        raise ArithmeticError('Missing reference volume at an active site')

    def values(e):
        ws,wz = e*ga,-shape.length-e*gt
        jac = (rs*wz-rz*ws)*h*h
        if min(jac)<=0:
            raise ArithmeticError('Orientation failure on receiver displacement')
        w = shape.length*t[active]+e*g[active]
        width = shape.radius(w)[0]
        xo = (r[active]-width)/(4.6-width)
        chi = current_order(np.c_[w,xo])
        trace = rs**2+rz**2+ws**2+wz**2+2*h*h
        part = stable_energy(trace[active],jac[active])
        return dict(step=e, reference_weight_energy=float(side.weight[active] @ part),
                    live_order_energy=float(vol@(chi*part)),
                    min_J=float(min(jac)), min_fractional_J=float(min(jac/j0)),
                    min_receiver_distance=float(min(r[active]-width)), max_r=float(max(r[active])))
    rows = [values(e) for e in (-.002,-.001,0.,.001,.002)]
    modes = {}
    for key in ('reference_weight_energy','live_order_energy'):
        modes[key] = dict(full_slope=(rows[4][key]-rows[0][key])/.004,
                          half_slope=(rows[3][key]-rows[1][key])/.002,
                          minus_work=rows[1][key]-rows[2][key],
                          plus_work=rows[3][key]-rows[2][key])
    return rows, modes


def run(radius,grade,dest):
    old = json.loads(RECORD.read_bytes())
    deep, source_chi, source = archived(radius,72,old)
    weights = fields(deep,source_chi)
    x=np.asarray(source['coordinates']).copy()
    x[0]-=np.log(deep.length)
    shape=Meridian(radius,x)
    _,outer=exterior(shape,48 if grade==1 else 70,return_side_field=True)
    coefficients={}
    for name,sign in [('ordered',1),('receiver',-1)]:
        side=Side(radius,weights,sign,grade)
        coeff,energy,residual=stationary_graph(side,shape)
        coefficients[name]=coeff
        r,rs,rz,h,jac,_=graph_values(side.a,side.t,shape,side.basis@coeff,side.deriv@coeff)
        print(json.dumps(dict(kind='graph', radius=radius, side=name, grade=grade, energy=energy,
                              max_normal_gradient=residual,min_J=float(min(jac)),max_J=float(max(jac)),
                              source_order_sha256=source['sha256']),allow_nan=False),flush=True)
        if sign==1:
            ts=np.linspace(0,1,1001)
            cut=shape.profile(ts)[0]*np.exp(-.5*np.exp(legval(2*ts-1,coeff)))
            print(json.dumps(dict(kind='material_cut',radius=radius,grade=grade,
                                  starting_material_radius=radius+.5,
                                  cut_radius_span=[float(min(cut)),float(max(cut))],
                                  entrance_cut_radius=float(cut[0]),end_cut_radius=float(cut[-1]),
                                  angular_stretch_span=(np.asarray([min(cut),max(cut)])/(radius+.5)).tolist(),
                                  slab_to_collar_lip_radial_tangent=float(-radius*np.exp(legval(-1,coeff)))),
                             allow_nan=False),flush=True)
        else:
            rows,slopes=free_receiver(shape,side,coeff,outer)
            print(json.dumps(dict(kind='released_receiver_mode',radius=radius,grade=grade,
                                  actual_order_side=True,rows=rows,slopes=slopes),allow_nan=False),flush=True)
    dest.parent.mkdir(parents=True,exist_ok=True)
    np.savez(dest,**coefficients,coordinates=x)
    print(json.dumps(dict(kind='snapshot',file=str(dest),sha256=hashlib.sha256(dest.read_bytes()).hexdigest())),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radius',type=int,choices=(1,2),required=True)
    parser.add_argument('--grade',type=int,choices=(1,2),required=True)
    parser.add_argument('--save',type=Path,required=True)
    args=parser.parse_args()
    run(args.radius,args.grade,args.save)
