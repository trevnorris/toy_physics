"""Variational attachment at fixed neck radius, prescribed tilt and wave action."""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial import Legendre, Polynomial
from numpy.polynomial.legendre import leggauss
from scipy.linalg import cholesky, eigh, solve_triangular
from scipy.optimize import minimize


RIM_BEND = 5.0
ACTION = 582.371869550
STEP = Polynomial([0, 0, 0, 10, -15, 6])
X = Polynomial([0, 1])


def jets(r, u, radial, angular):
    """Gradient and Hessian in the orthonormal spherical frame for f(r)P(u)."""
    s = np.sqrt(1-u*u)
    f, fr, frr = radial
    p, pu, puu = angular(u), angular.deriv()(u), angular.deriv(2)(u)
    g = fr*p/r - f*u*pu/r**2
    return np.array([fr*p, -s*f*pu/r, frr*p, -s*pu*(fr/r-f/r**2),
                     g+f*(1-u**2)*puu/r**2, g])


class Neck:
    def __init__(self, modes, angular, free, wave_r, wave_l, order):
        x, weight = leggauss(order)
        r = np.r_[1.5+x/2, 2.5+x/2]
        self.r, u = [v.ravel() for v in np.meshgrid(r, x, indexing="ij")]
        r = self.r
        self.wt = np.pi*r**2*np.kron(np.tile(weight, 2), weight)
        t = r-1
        inner = t < 1
        self.chi = STEP(np.clip(3-r, 0, 1))
        self.labels = [(l,m) for l in range(angular+1) for m in range(modes)]
        base = jets(r,u,((r*r-1)/10,r/5,r**0/5),Legendre.basis(0))
        taper = (1+X)*STEP
        tilt = tuple(-np.where(inner, taper.deriv(i)(t), [r,1,0][i])/5 for i in range(3))
        shapes = [base, jets(r,u,tilt,Legendre.basis(1))]
        for l,m in self.labels:
            poly = X**(1 if free else 2)*(1-X)**2*Legendre.basis(m).convert(kind=Polynomial)(2*X-1)
            radial = tuple(np.where(inner,poly.deriv(i)(t),0.) for i in range(3))
            shapes.append(jets(r,u,radial,Legendre.basis(l)))
        self.shapes = np.moveaxis(np.array(shapes),0,-1).copy()
        dp,dq,da,_,dd,df = self.shapes[:,:,2:]
        lap = da+dd+df
        gram = dp.T@(self.wt[:,None]*dp) + dq.T@(self.wt[:,None]*dq)
        gram += 2*RIM_BEND/9 * lap.T@(self.wt[:,None]*lap)
        self.pre = cholesky(gram,lower=True)

        # Independent endpoint hats and interval bubbles implement a node at
        # r=1, continuity at r=2 and a natural reflecting condition at r=3.
        local = np.where(inner,t,t-1)
        rv = [np.where(inner,t,2-t),np.maximum(0.,t-1)]
        dr = [np.where(inner,1.,-1.),np.where(inner,0.,1.)]
        for region in (inner,~inner):
            for m in range(wave_r-1):
                poly = X*(1-X)*Legendre.basis(m).convert(kind=Polynomial)(2*X-1)
                rv.append(poly(local)*region)
                dr.append(poly.deriv()(local)*region)
        ang = np.array([Legendre.basis(l)(u) for l in range(wave_l+1)]).T
        da = np.array([Legendre.basis(l).deriv()(u) for l in range(wave_l+1)]).T

        def tensor(rad, lat):
            return (np.array(rad).T[:,:,None]*lat[:,None,:]).reshape(len(r),-1).copy()
        self.val = tensor(rv,ang)
        self.dr = tensor(dr,ang)
        self.dt = tensor(rv,-np.sqrt(1-u**2)[:,None]*da/r[:,None])

    def evaluate(self,beta,coeff):
        p,q,a,b,d,f = self.shapes@np.r_[1.,beta,coeff]
        w = np.sqrt(1+p*p+q*q)
        ptp = a*p*p+2*b*p*q+d*q*q
        tr = a+d+f-ptp/w**2
        area = self.wt@w
        bend = self.wt@(tr**2/(9*w))
        sw = self.wt*self.chi/w

        def gram(mat,weights):
            return mat.T@(weights[:,None]*mat)
        d_cross = q[:,None]*self.dr-p[:,None]*self.dt
        stiff = gram(self.dr,sw)+gram(self.dt,sw)+gram(d_cross,sw)
        mass = gram(self.val,sw*w**2)
        lam,vec = eigh(stiff,mass,subset_by_index=(0,0),check_finite=False)
        lam,vec = lam[0],vec[:,0]
        omega = np.sqrt(lam)
        ur,ut,val = self.dr@vec,self.dt@vec,self.val@vec
        cross = q*ur-p*ut
        # Differentiate the wall integrand and the mass-normalized Rayleigh
        # quotient at its eigenfunction (fixed action, including the mass term).
        h = 2*RIM_BEND*tr/(9*w)
        dw = 1-RIM_BEND*tr**2/(9*w**2)+2*h*ptp/w**3
        v = ACTION*self.chi/(2*omega)
        dw -= v*((ur**2+ut**2+cross**2)/w**2+lam*val**2)
        gp = dw*p/w-2*h*(a*p+b*q)/w**2-2*v*cross*ut/w
        gq = dw*q/w-2*h*(b*p+d*q)/w**2+2*v*cross*ur/w
        local = np.array([gp,gq,h*(1-p*p/w**2),-2*h*p*q/w**2,h*(1-q*q/w**2),h])
        derivative = np.einsum("ni,nij->j",local*self.wt,self.shapes[:,:,1:])
        return np.array([area,bend,omega]),derivative

    def value(self,beta,coeff):
        parts,grad = self.evaluate(beta,coeff)
        return parts@np.array([1.,RIM_BEND,ACTION]),grad

    def relax(self,beta,initial):
        def objective(z):
            a=solve_triangular(self.pre.T,z,lower=False)
            energy,grad=self.value(beta,a)
            return energy,solve_triangular(self.pre,grad[1:],lower=True)
        fit=minimize(objective,self.pre.T@initial,jac=True,method="BFGS",
                     options={"gtol":3e-6})
        coeff=solve_triangular(self.pre.T,fit.x,lower=False)
        parts,grad=self.evaluate(beta,coeff)
        residual=np.linalg.norm(solve_triangular(self.pre,grad[1:],lower=True),np.inf)
        if residual>2e-5 or not np.isfinite(fit.fun):
            raise ArithmeticError(f"tilt {beta} unresolved: {fit.message}, gradient {residual}")
        return dict(beta=beta,energy=fit.fun,parts=parts.tolist(),slope=grad[0],
                    residual=residual,coeff=coeff.tolist())


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("output",type=Path)
    ap.add_argument("--shape",type=int,nargs=2,default=[3,2])
    ap.add_argument("--wave",type=int,nargs=2,default=[6,5])
    ap.add_argument("--quad",type=int,default=22)
    ap.add_argument("--clamped",action="store_true")
    ap.add_argument("--divisions",type=int,default=12)
    args=ap.parse_args()
    model=Neck(*args.shape,not args.clamped,*args.wave,args.quad)
    zero=np.zeros(len(model.labels))
    # Check the full shape gradient including the wave term away from the origin.
    test=0.03*np.cos(np.arange(len(zero)))
    direction=np.sin(1+np.arange(len(zero)))
    delta=1e-5
    e,grad=model.value(.7,test)
    fd=(model.value(.7+delta,test+delta*direction)[0]-
        model.value(.7-delta,test-delta*direction)[0])/(2*delta)
    print(f"gradient check fd={fd:.9f} analytic={grad[0]+grad[1:]@direction:.9f}",flush=True)
    if abs(fd-grad[0]-grad[1:]@direction)>2e-5:
        raise ArithmeticError("Wall and wave shape derivative failed")
    base=model.relax(0.,zero)
    initial=np.array(base["coeff"])
    # Check local stability in all admitted shape modes in the whitened coordinates.
    def hessian(beta,coeff):
        z=model.pre.T@coeff
        columns=[]
        for step in np.eye(len(z))*1e-3:
            def gradient(v):
                a=solve_triangular(model.pre.T,v,lower=False)
                return solve_triangular(model.pre,model.value(beta,a)[1][1:],lower=True)
            columns.append((gradient(z+step)-gradient(z-step))/0.002)
        matrix=np.array(columns).T
        return np.linalg.eigh((matrix+matrix.T)/2)
    curv,vec=hessian(0.,initial)
    print(f"symmetric start E={base['energy']:.9f} min_shape_curvature={curv[0]:+.8f}",flush=True)
    candidates=[base]
    for amp in (-1.,1.):
        a=initial+solve_triangular(model.pre.T,amp*vec[:,0],lower=False)
        candidates.append(model.relax(0.,a))
    lowest=min(item["energy"] for item in candidates)
    # Pick the reflected equivalent branch heading downhill with positive tilt.
    base=min((item for item in candidates if item["energy"]<=lowest+1e-7),key=lambda item:item["slope"])
    initial=np.array(base["coeff"])
    curv,_=hessian(0.,initial)
    states=[base]
    print(f"relaxed start E={lowest:.9f} min_shape_curvature={curv[0]:+.8f} "
          f"shape release lowers E by {model.value(0.,zero)[0]-lowest:.9f}",flush=True)
    if curv[0]<-2e-5:
        raise ArithmeticError("Starting attachment is a saddle")
    for beta in sorted(set(np.linspace(0,3,args.divisions+1)[1:])|{0.025,0.05,0.1,1.1,2.,2.5,3.}):
        state=model.relax(float(beta),initial)
        initial=np.array(state["coeff"])
        states.append(state)
        print(f"beta={beta:.6f} dE={state['energy']-base['energy']:+.9f} slope={state['slope']:+.8f} "
              f"residual={state['residual']:.2e}",flush=True)
    # A path connecting the actual states guards against downward jumps between
    # disconnected optimization minima. Interpolate the graph coefficients.
    largest=0.
    for left,right in zip(states,states[1:]):
        for frac in (0.25,0.5,0.75):
            beta=(1-frac)*left["beta"]+frac*right["beta"]
            coeff=(1-frac)*np.array(left["coeff"])+frac*np.array(right["coeff"])
            largest=max(largest,model.value(beta,coeff)[0]-base["energy"])
    print(f"joined path sampled excess maximum={max(largest,max(s['energy']-base['energy'] for s in states)):.9f}",flush=True)
    args.output.write_text(json.dumps(dict(shape=args.shape,wave=args.wave,quad=args.quad,
        free_angle=not args.clamped,stiffness=RIM_BEND,action=ACTION,states=states,
        min_initial_shape_curvature=curv[0],interpolated_excess_max=largest),indent=2)+"\n")


if __name__=="__main__":
    main()
