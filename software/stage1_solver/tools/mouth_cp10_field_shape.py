"""Read-only six-force test on the CP10 field chart, excluding the interior crystal.

Evaluate the frozen Q1 fields directly on their moving quadrature grid. Complex
steps differentiate each smooth branch, including the distant flat trace, the
vortex subtraction surface and the existing exterior material integration.
"""
import argparse
import inspect
import json
from itertools import product
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import legval, legder, leggauss

from mouth_exact_checkpoint_control import CHECKPOINT_SHA, acquired


GAUSS = (.5 - .5 / np.sqrt(3), .5 + .5 / np.sqrt(3))
BITS = tuple(product((0, 1), repeat=3))
COMPONENTS = ("amplitude", "vortex_potential", "phase", "order_gradient",
              "double_well", "vortex_side", "exterior_crystal")


def spline_value(spline, position, order=0):
    """Evaluate real spline polynomial and its exact directional derivative."""
    result = spline(position.real, order)
    if np.iscomplexobj(position):
        result = result + 1j * position.imag * spline(position.real, order+1)
    return result


def bump(t, length):
    p = (t-.06)/.56
    active = (p.real > 0) & (p.real < 1)
    return (np.where(active, 16*(p*(1-p))**2, 0),
            np.where(active, 32*p*(1-p)*(1-2*p)/(.56*length), 0))


class Shape:
    def __init__(self, held, shift, offset):
        self.held = held
        self.par = held.params + offset[:4]*[1., .25, .25, .25]
        self.extra = held.extra+.25*offset[4]
        self.shift = shift+offset[5]
        self.length = held.length*np.exp(offset[0])

    def radius(self, z):
        t = z/self.length
        outside = (t.real < 0) | (t.real > 1)
        t = np.where(t.real < 0, 0, np.where(t.real > 1, 1, t))
        m, a = self.par[2:]
        exit_shift = 1.5*np.expm1(self.par[1])
        middle = t*(1-t)
        s = t*self.held.unit_length
        r = spline_value(self.held.reference, s) + exit_shift*t + 4*m*middle + 4*a*middle*(2*t-1) + 16*self.extra*middle**2
        slope = (spline_value(self.held.reference, s, 1)*self.held.unit_length + exit_shift
                 + 4*m*(1-2*t) + 4*a*(-6*t*t+6*t-1) + 32*self.extra*middle*(1-2*t))
        return 4*r, np.where(outside, 0, 4*slope/self.length)


def derivative_vortex(vortex, position):
    f, df = vortex(position.real)
    if np.iscomplexobj(position):
        env = inspect.getclosurevars(vortex).nonlocals
        lo, sol = env["lo"], env["sol"]
        f1, df1 = sol.sol(np.maximum(position.real, lo), 1)
        f = f + 1j*position.imag*np.where(position.real < lo, sol.y[0, 0]/lo, f1)
        df = df + 1j*position.imag*np.where(position.real < lo, 0, df1)
    return f, df


def axes(length, ref, model):
    nr = len(model.rr[0])-1
    nz = len(model.z[1])-nr-1
    s = np.linspace(0., 1., nr+1)
    origin, end = ref-.12, length-.12
    lower = origin*s/(1+np.expm1(24*np.log(origin/end))*s**24)**(1/24)
    near = np.linspace(end, length, nz+1)
    upper = np.linspace(length, length+.12, nz+1)
    far = length+.12+(22-length-.12)*(np.sinh(2*s)/np.sinh(2))**2
    return (model.z[0], np.r_[lower[:-1],near], np.r_[upper[:-1],far])


def corners(field):
    nx, na, nt = field.shape[-3:]
    return np.stack([field[...,i:i+nx-1,k:k+na-1,j:j+nt-1] for i,k,j in BITS])


def total(arr):
    return np.sum(arr, dtype=np.clongdouble if np.iscomplexobj(arr) else np.longdouble)


class Bath:
    def __init__(self, trial, base, u, reference_depth):
        from mouth_internal_collar import CollarPair
        from mouth_breathing_mobile_patch import WIDTH, WELL, ORDER_INERTIA, KAPPA
        self.model = m = trial.model
        self.width, self.well, self.inertia, self.kappa = WIDTH, WELL, ORDER_INERTIA, KAPPA
        self.ref = reference_depth
        self.u = u
        self.fixed = np.isin(trial.number,trial.fixed_order)
        self.blocks = {}
        for key, ids in m.ids.items():
            self.blocks[key] = u[:,ids]
        p = CollarPair(trial.shape,base.co,base.na,base.nt)
        ext = ~p.inside
        self.et, self.t = p.eta[ext], p.z[ext]
        self.reg = np.exp(legval(2*self.t-1, base.co))
        self.pt = 2*legval(2*self.t-1,legder(base.co))
        vals = base.labels.reshape(-1,2)[p.id[ext]]
        self.S = 4+np.einsum('qi,qi->q',p.B[ext],vals[:,:,0])
        inv = np.stack([np.einsum('qi,qik->qk',g[ext],vals) for g in (p.E,p.T)], axis=-1)
        self.D = np.linalg.det(inv)
        self.inv = np.linalg.inv(inv)
        self.weights = p.qw[ext]*self.S**2*self.D

    def field_block(self, key, zs):
        ids = self.model.ids[key]
        field = self.blocks[key].astype(np.result_type(zs,self.u))
        values = 1/(1+np.exp(zs/self.width))
        if key[0] == 2:
            values[-1] = 0
        field[2] = np.where(self.fixed[ids], values, field[2])
        return field

    def volume(self, sh, all_z):
        parts = np.zeros(5, dtype=np.result_type(sh.par, np.longdouble))
        for (slab, side), ids in self.model.ids.items():
            xs, aa, zs = self.model.rr[side], self.model.a, all_z[slab]
            dx, da, dz = np.diff(xs)[:, None, None], np.diff(aa)[None,:,None], np.diff(zs)[None,None,:]
            field = corners(self.field_block((slab, side),zs))
            for p, a, w in product(GAUSS, repeat=3):
                at = np.array([p,a,w])
                factors = np.where(BITS,at,1-at)
                basis = factors.prod(axis=1)
                raw = np.stack([(2*np.asarray(BITS)[:,i]-1)*factors[:,(i+1)%3]*factors[:,(i+2)%3]
                                 for i in range(3)])
                vals = np.einsum('i,ifxat->fxat',basis,field)
                along = np.einsum('bi,ifxat->bfxat',raw,field)
                x, theta, z = xs[:-1,None,None]+p*dx, aa[None,:-1,None]+a*da, zs[None,None,:-1]+w*dz
                b, slope = sh.radius(z)
                h, hz = bump(z/sh.length,sh.length)
                width = 28-b if side else b
                r = b+x*width if side else x*b
                k, kx = (1-x,-1) if side else (x**3,3*x**2)
                sn, co = np.sin(theta), np.cos(theta)
                polar = (3*co**2-1)/2
                rr = r+sh.shift*h*k*polar
                rx = width+sh.shift*h*kx*polar
                ra = -3*sh.shift*h*k*co*sn
                rz = slope*(1-x if side else x)+sh.shift*hz*k*polar
                wt = .5*np.pi*rr**2*sn*rx*dx*da*dz
                gr = along[0]/dx/rx
                gt = (along[1]/da-ra*gr)/rr
                gw = along[2]/dz-rz*gr
                squared = gr**2+gt**2+gw**2
                dot = gr[1]*gr[2]+gt[1]*gt[2]+gw[1]*gw[2]
                f,_ = derivative_vortex(self.model._vortex,rr*sn)
                amp,_,chi = vals
                n = (f*amp)**2
                y = amp**2
                phase = .5*n*(squared[1]-self.inertia*dot**2/(n+self.inertia*squared[2]))
                ed = (.5*f**2*squared[0], f**10*(y-1)**2*(y**3+2*y*y+3*y+4)/20,
                      phase, .5*self.kappa*squared[2], self.well*(chi*(1-chi))**2)
                parts += [total(wt*en) for en in ed]
        return parts

    def sides(self, sh, zs):
        a = self.model.a
        dt, dz = np.diff(a)[:,None], np.diff(zs)[None,:]
        out = 0.
        for p, t in product(GAUSS, repeat=2):
            theta = a[:-1,None]+p*dt
            sn, co = np.sin(theta), np.cos(theta)
            z = zs[None,:-1]+t*dz
            b, _ = sh.radius(z)
            h,_ = bump(z/sh.length,sh.length)
            r = b+sh.shift*h*(3*co**2-1)/2
            rt = -3*sh.shift*h*co*sn
            f,df = derivative_vortex(self.model._vortex,r*sn)
            area = np.pi*dt*dz*sn*r*(r*sn-rt*co)*f*df
            for side in (0,1):
                face = self.u[0,self.model.ids[1,side][0 if side else -1]]
                v = ((1-p)*((1-t)*face[:-1,:-1]+t*face[:-1,1:])
                     +p*((1-t)*face[1:,:-1]+t*face[1:,1:]))
                out += (1-2*side)*total(.5*area*(v*v-1))
        return out

    def exterior(self, sh, zs):
        L = sh.length
        z = self.t*L
        b, slope = sh.radius(z)
        r = b*np.exp(-self.reg*self.et)
        re = -self.reg*r
        rt = r*(L*slope/b-self.et*self.reg*self.pt)
        spatial = np.zeros(self.inv.shape,dtype=sh.par.dtype)
        spatial[:,0,0], spatial[:,0,1], spatial[:,1,1] = re,rt,L
        F = spatial@self.inv
        K = F@F.transpose(0,2,1)
        hoop = r/self.S
        jac = -re*L*hoop**2/self.D
        x = (r-b)/(28-b)
        def locate(nodes, values):
            i = np.clip(np.searchsorted(nodes.real, values.real)-1,0,len(nodes)-2)
            return i,(values-nodes[i])/(nodes[i+1]-nodes[i])
        ix, fx = locate(self.model.rr[1],x)
        iz, fz = locate(zs,z)
        field = self.field_block((1,1),zs)[2]
        h,hz = bump(self.t,L)
        k = (28-r)/(28-b)
        d,s,v = -h/(28-b), k*h/r, k*(hz+h*slope/(28-b))
        q, qw = leggauss(4)
        out = 0.
        for i,(angle,da) in enumerate(zip(self.model.a,np.diff(self.model.a))):
            for p,weight in zip((q+1)/2,qw/2):
                a = angle+p*da
                P, Pt = (3*np.cos(a)**2-1)/2, -3*np.sin(a)*np.cos(a)
                ar,rad,az = 1+sh.shift*d*P, 1+sh.shift*s*P, sh.shift*v*P
                J = jac*ar*rad**2
                tr = ar**2*K[:,0,0]+2*ar*az*K[:,0,1]+(1+az**2)*K[:,1,1]+hoop**2*((sh.shift*s*Pt)**2+2*rad**2)
                W = (tr-4)/2-np.log(J)+(J-1)**2/4
                chi = sum(wx*wz*wa*field[ix+di,i+dk,iz+dj]
                          for di,wx in ((0,1-fx),(1,fx)) for dj,wz in ((0,1-fz),(1,fz))
                          for dk,wa in ((0,1-p),(1,p)))
                out += total(self.weights*W*chi*(da*weight*np.sin(a)))
        return out

    def energy(self, offset):
        sh = Shape(self.model.shape,self.model.polar_shift,offset)
        zz = axes(sh.length,self.ref,self.model)
        return np.r_[self.volume(sh, zz),self.sides(sh,zz[1]),self.exterior(sh,zz[1])]


def run(folder):
    scratch = Path('_scratch/s11c').resolve()
    build,x,q,u,_,rec = acquired(scratch/'r4_melted_survival01',scratch/'r4_exact_support01/patch')
    trial = build(x)
    owner = inspect.getclosurevars(build).nonlocals["owner"]
    bath = Bath(trial,owner.base,u,owner.reference_depth)
    zero = np.zeros(6)
    neutral = bath.energy(zero)
    trial.model.material_load = trial.ext
    native = trial.model.evaluate(trial.transported(u),True)[0]
    report = dict(checkpoint=CHECKPOINT_SHA, components=COMPONENTS, fields=neutral.astype(float).tolist(),
                  native=float(native), static_delta=float(neutral.sum()-native), modes=[])
    if abs(report['static_delta']) > 2e-9:
        raise ValueError(f"Frozen field energy mismatch: {report}")
    for i in range(6):
        unit = np.eye(6)[i]
        cs = -bath.energy(unit*1e-24j).imag/1e-24
        case = dict(mode=i, analytic_components=cs.astype(float).tolist(), analytic=float(cs.sum()), central={})
        for step in (1e-4,5e-5,1e-5):
            minus, plus = (bath.energy(sign*step*unit) for sign in (-1,1))
            fd = (minus-plus)/(2*step)
            case['central'][str(step)] = dict(force=float(fd.sum()), delta=float((fd-cs).sum()),
                                               components=fd.astype(float).tolist())
        report['modes'].append(case)
        folder.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        print(i,case['analytic'],case['central']['1e-05']['delta'],flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("output",type=Path)
    run(ap.parse_args().output)
