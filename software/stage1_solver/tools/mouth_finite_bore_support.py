"""Finite elastic work of a registered collar formed from a straight flat-slab bore.

Coordinates S,Z are the radius and depth of a material tag in the unloaded
unit-thickness slab.  R is the half-order bore. Each side has a smooth positive radial stretch.
Axial registration and the radial-spacing form are imposed throughout the finite collar.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import Legendre, leggauss
from scipy.interpolate import RegularGridInterpolator
from scipy.linalg import solve_triangular
from scipy.optimize import minimize

from mouth_sealed_shape_relax import Meridian
from mouth_unified_puncture import archived, RECORD
from mouth_waist_outer_order import exterior


def fields(shape, source_order):
    # Native meridional material chart, sampled once on the long preparation.
    _, (w, xo, outside) = exterior(shape, 48, return_side_field=True)
    m = source_order.shape[1] - 1
    inside = RegularGridInterpolator((np.linspace(0, 1, 2*m+1), np.linspace(0, 1, m+1)),
                                     source_order)
    receiver = RegularGridInterpolator((w, xo), outside)

    def sample(a, t):
        b, _, _ = shape.profile(t)
        y = 1-a
        if np.all(a > 0):
            return inside(np.c_[t, y])
        return receiver(np.c_[shape.length*t, -a*b/(4.6-b)])

    return sample


def quadrature(side, nz, ns):
    z, wz = leggauss(nz)
    x, wx = leggauss(ns)
    a, t = np.meshgrid(side*(x+1)/4, (z+1)/2)
    return a.ravel(), t.ravel(), (wz[:, None]*wx[None, :]/8).ravel()


def graph(a,t,shape,c):
    """F in the S,Z basis at a material tag. w=L*t, r=b(t)*exp(-p(t)*a)."""
    basis = np.array([Legendre.basis(j)(2*t-1) for j in range(len(c))]).T
    deriv = np.array([2*Legendre.basis(j).deriv()(2*t-1) for j in range(len(c))]).T
    return graph_values(a,t,shape,basis@c,deriv@c)


def graph_values(a,t,shape,log_p,dlog_dt):
    b,slope,_ = shape.profile(t)
    L = shape.length
    p = np.exp(log_p)
    ap = a*p
    r = b*np.exp(-ap)
    rs = -p*r
    rz = -r*(slope*L/b - ap*dlog_dt)
    h = r/(shape.neck+a)
    J = -rs*L*h*h
    return r,rs,rz,h,J,ap


class Side:
    def __init__(self,radius,order,sign,grade):
        self.radius,self.sign = radius,sign
        self.nt,nz,ns = ((6,64,20),(12,100,32),(18,144,40))[grade]
        self.a,self.t,weight = quadrature(sign,nz,ns)
        self.S = radius+self.a
        self.chi = order(self.a,self.t)
        self.weight = 4*np.pi*self.S**2*weight*self.chi
        self.basis = np.array([Legendre.basis(j)(2*self.t-1) for j in range(self.nt)]).T
        self.deriv = np.array([2*Legendre.basis(j).deriv()(2*self.t-1) for j in range(self.nt)]).T

    def leg(self,shape,seed=None):
        # All labels keep affine axial registration and exponential radial spacing
        # on each side. The degree-dependent normal stretch relaxes freely.
        L=shape.length
        b,slope,_=shape.profile(self.t)
        def evaluate(c, full=False):
            logp=self.basis@c
            dlog=self.deriv@c
            if np.max(logp)>7 or np.min(logp)<-40:
                return 1e60,np.zeros_like(c)
            p=np.exp(logp)
            ap=self.a*p
            if max(abs(ap))>50:
                return 1e60,np.zeros_like(c)
            r=b*np.exp(-ap)
            rs=-p*r
            rz=-r*(slope*L/b-ap*dlog)
            h=r/self.S
            J=-rs*L*h*h
            ld=np.log(J)
            tr=rs**2+rz**2+L**2+2*h**2
            scale=1/np.sqrt(J)
            q=ld-scale*tr/4
            energy=self.weight@((scale*tr-4+ld**2)/2)
            g=self.weight*(scale*(rs**2*(1-ap) -(rz**2+2*h*h)*ap
                                  +rz*r*ap*dlog)+q*(1-3*ap))
            gd=self.weight*scale*rz*r*ap
            grad=self.basis.T@g+self.deriv.T@gd
            if not full:
                return energy,grad
            # Released axial-displacement diagnostic, dw=a*t*(1-t), zero
            # at the seam and on both flat ends. W would need zero work here
            # for an ungripped internal equilibrium in this direction.
            dws=self.t*(1-self.t)
            dwz=-self.a*(1-2*self.t)
            pws=-q*rz/(-rs*L)
            pwz=-scale*L-q/L
            free_w_work=self.weight@(pws*dws+pwz*dwz)
            # Envelope force from Piola virtual work on the very same shape
            # coordinates. Width derivative b_q=4t(1-t), and lengthening
            # keeps b(t) fixed, changing only wz by -1.
            bq=4*self.t*(1-self.t)
            dbq=4*(1-2*self.t)
            drq=np.exp(-ap)*bq
            drsq=-p*drq
            drzq=-np.exp(-ap)*(dbq-bq*ap*dlog)
            dhq=drq/self.S
            prs=scale*rs+q/rs
            prz=scale*rz
            ph=2*(scale*h+q/h)
            radial=-self.weight@(prs*drsq+prz*drzq+ph*dhq)
            axial=self.weight@pwz
            metrics=dict(energy=float(energy),grad_norm=float(np.linalg.norm(grad,np.inf)),
                         min_J=float(min(J)),max_J=float(max(J)),min_hoop=float(min(h)),
                         max_hoop=float(max(h)), min_p=float(min(p)),max_p=float(max(p)),
                         min_order=float(min(self.chi)),max_order=float(max(self.chi)),
                         signed_work_release_axial_mode=float(free_w_work),
                         center_force_from_stress=[float(radial),float(axial)])
            return metrics,dict(gradient=np.stack((rs,rz,0*rs,0*rs-L,h),axis=-1),J=J,
                                r=r,w=self.t*L,a=self.a,t=self.t,chi=self.chi,
                                weights_reference_with_order=self.weight,coeffs=c)
        if seed is None:
            seed=np.zeros(self.nt)
        sol=minimize(evaluate,seed,jac=True,method='BFGS',options=dict(gtol=3e-8,maxiter=3000))
        metrics,field=evaluate(sol.x,True)
        metrics.update(iterations=sol.nit,status=sol.message)
        if metrics['grad_norm']>2e-4 or not np.isfinite(metrics['energy']):
            raise ArithmeticError(f"Unresolved graph gradient {metrics}")
        return metrics,field,sol.x


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--radius', type=int, choices=(1,2), required=True)
    parser.add_argument('--grade', type=int, choices=(0,1,2), required=True)
    parser.add_argument('--save', type=Path, required=True)
    args = parser.parse_args()
    args.save.mkdir(parents=True,exist_ok=True)
    manifest = json.loads(RECORD.read_bytes())
    base, chi, source = archived(args.radius, 72, manifest)
    params = np.array(source['coordinates'])
    order = fields(base, chi)
    sides = [Side(args.radius, order, sign, args.grade) for sign in (1,-1)]
    inputs = dict(radius=args.radius, grade=args.grade, reference_width=1.,
                  depth_registration="Z=1-t", order_input=source['sha256'],
                  tool_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    for label, factor in (('long', 1.), ('unit', 1/base.length)):
        ref = params.copy()
        ref[0] += np.log(factor)
        center = Meridian(args.radius, ref)
        seeds = []
        for s in sides:
            metrics, field, seed = s.leg(center)
            seeds.append(seed)
            tag = 'ordered' if s.sign==1 else 'receiver'
            name = f'r{args.radius}_{label}_g{args.grade}_{tag}.npz'
            np.savez_compressed(args.save/name, **field, radius=args.radius,
                                side=s.sign,nt=s.nt,coordinates=ref)
            print(json.dumps(dict(**inputs, kind='central', depth=label, L=center.length,
                                  side=tag, field=name, **metrics)), flush=True)
        for axis in ('radial','axial'):
            rows = []
            for h in (-.002,-.001,.001,.002):
                pert = ref.copy()
                if axis=='radial':
                    pert[2] += h/args.radius
                else:
                    pert[0] += h
                shape = Meridian(args.radius, pert)
                legs=[]
                for s,seed in zip(sides,seeds):
                    metrics, field, _ = s.leg(shape,seed)
                    legs.append(metrics)
                    tag = 'ordered' if s.sign==1 else 'receiver'
                    name = f'r{args.radius}_{label}_g{args.grade}_{tag}_{axis}_{h}.npz'
                    np.savez_compressed(args.save/name, coeffs=field['coeffs'],
                                        radius=args.radius,side=s.sign,nt=s.nt,coordinates=pert)
                rows.append((h if axis=='radial' else shape.length,[r['energy'] for r in legs]))
            forces=[]
            for lo,hi in ((0,1),(2,3),(1,2)):
                forces.append(((np.array(rows[lo][1])-np.array(rows[hi][1]))/
                               (rows[hi][0]-rows[lo][0])).tolist())
            print(json.dumps(dict(**inputs,kind='forces',depth=label,axis=axis,rows=rows,
                                  forces_by_side=forces,total=np.sum(forces,axis=1).tolist())),flush=True)


if __name__=='__main__':
    main()
