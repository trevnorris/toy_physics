"""Test the depths selected by a common pressure or a common area tension.

Calibrate just one axial stationary shape (R=H=1). Scan for minima in a
declared depth interval, then solve the force equation and check an independent
local stress integral at each interior minimum. Units for restoring coefficients
are divided by the squared common weak phase bias.
"""
import argparse
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.sparse.linalg import spsolve

from mouth_carrier_depth_work import axial_work,geo_work
from mouth_carrier_weak_radius import PhysicalDomain,translation
from mouth_flared_receiver import EXIT,emit


class DepthDomain(PhysicalDomain):
    def assemble(self,grids):
        super().assemble(grids)
        # Plane identification in the parent uses exact chart coordinates.
        # The true topological roof is known on these quadratic chart faces.
        for block in (3,4):
            pg,vg,indices=self.blocks[block]
            last=indices.reshape(len(pg)-1,len(vg)-1,3,3)[-1,:,2,:]
            self.points[last,1]=self.roof
            first=indices.reshape(len(pg)-1,len(vg)-1,3,3)[0,:,0,:]
            self.points[first,1]=round(EXIT,11)


def run(args):
    ref=Path('_scratch/s11c/carrier_depth_one_01')
    if json.loads((ref/'outcome.json').read_text())['exitCode'] or not json.loads((ref/'limit-validation.json').read_text())['verified']:
        raise ValueError('Unverified one-mouth calibration')
    reference=json.loads((ref/'stdout').read_text())
    def domain(h):
        roof=h*round(EXIT+16/h,9)
        return DepthDomain(args.radius,h,args.mesh,24.,roof,1.)
    def carrier(h):
        dom=domain(h)
        phi=dom.bulk.astype(float)
        free=dom.free
        phi[free]=spsolve(dom.stiff[free][:,free],-(dom.stiff@phi)[free])
        j=dom.stiff@phi
        energy=.5*phi@j
        if abs(2*energy+j[dom.feed].sum())>2e-7*(2*energy):
            raise ArithmeticError('Bath currents disagree')
        return dom,phi,energy
    @lru_cache(None)
    def energy(h):
        return carrier(h)[-1]
    # All area/volume expressions depend on the unit reference profile alone.
    geometry=domain(1.)
    @lru_cache(None)
    def terms(h):
        step=.003*h
        grad=(energy(h+step)-energy(h-step))/(2*step)
        dv,da=geo_work(geometry,args.radius,h)
        return grad,reference['required_pressure']*dv,reference['required_tension']*da
    depths=np.array([.25,.4,.6,.8,1.,1.25,1.6,2.,2.5])
    for fraction in (0.,.5,1.):
        def residual(h):
            grad,p,a=terms(h)
            return grad+fraction*p+(1-fraction)*a
        scan=np.array([residual(h) for h in depths])
        emit(test='restoring_scan',radius=args.radius,mesh=args.mesh,p_fraction=fraction,
             depths=depths.tolist(),gradients=scan.tolist())
        if fraction==0:
            g,p,a=terms(1.)
            emit(test='fixed_depth_range',radius=args.radius,required_outward=-g,
                 smallest_restoration=min(p,a),largest_restoration=max(p,a))
        for lo,hi,fl,fh in zip(depths[:-1],depths[1:],scan[:-1],scan[1:]):
            if fl>=0 or fh<0:
                continue
            root=brentq(residual,lo,hi,xtol=3e-5)
            dom,phi,e=carrier(root)
            m,work=translation(dom)
            grad,p,a=terms(root)
            local=axial_work(dom,phi,args.radius,root)
            step=.012*root
            curvature=(residual(root+step)-residual(root-step))/(2*step)
            emit(test='depth_minimum',p_fraction=fraction,radius=args.radius,mesh=args.mesh,
                 selected_depth=root,conductance=2*e,carrier_mass=m,ratio_per_bias=2*e/m,
                 relative_to_reference=(2*e/m)/reference['ratio_per_bias'],gradient=grad,
                 local_carrier_gradient=local,net_gradient=residual(root),curvature=curvature)
            if (not np.isclose(grad,local,rtol=.02,atol=2e-7) or curvature<=0
                    or abs(residual(root))>max(2e-6,1e-4*abs(grad)) or not np.isclose(m,work)):
                raise ArithmeticError('Interior minimum failed local independent check')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',type=float,required=True)
    p.add_argument('--mesh',type=int,default=22)
    run(p.parse_args())
