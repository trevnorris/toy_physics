"""Release a smooth internal waist bend on the jointly relaxed four-coordinate mouth."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from mouth_circulation_material import collar
from mouth_sealed_shape_relax import Meridian
from mouth_supercore_shape import Relax, emit


class WaistRelax(Relax):
    def shape(self,x):
        return Meridian(4.,self.ref.params+x[:4]*self.scale,extra=x[4]/4)

    def record(self,x):
        energy,gradient=self.eval(x)
        _,_,_,self.chi,info=self.cache
        sh=self.shape(x)
        self.q=sh.params.copy()
        path=self.folder/'state.npz'
        temp=path.with_suffix('.partial.npz')
        np.savez_compressed(temp,radius=4,shape_params=self.q,shape_extra=sh.extra,
                            order=self.chi,registration=self.co)
        temp.replace(path)
        emit(kind='accepted_waist',energy=energy,shape_params=self.q.tolist(),shape_extra=sh.extra,
             depth=sh.length,exit=sh.exit,waist=sh.wall()['minimum_radius'],
             shape_forces=(-gradient).tolist(),axial_outward=-gradient[0]/sh.length,
             exit_outward=-gradient[1]/(sh.exit*.25),localized_waist_outward=-gradient[4],
             pocket_mean_order=float(self.chi[self.seed<.01].mean()),
             pocket_min_order=float(self.chi[self.seed<.01].min()),
             shaft_mid_order=float(self.chi[self.n,0]),
             upper_collar=collar(sh,self.co,False),lower_collar=collar(sh,self.co,True),
             state_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),**info)

    def continue_shape(self,path):
        guard=json.loads((path.parent/'outcome.json').read_bytes())
        limits=json.loads((path.parent/'limit-validation.json').read_bytes())
        rows=[json.loads(l) for l in (path.parent/'stdout').read_bytes().splitlines()]
        if guard['exitCode'] or guard['childOutcome']['exitCode'] or not limits['verified'] or int(limits['actual']['memory.swap.max']):
            raise ValueError('Unverified four-coordinate source')
        state=[r for r in rows if r['kind']=='accepted'][-1]
        ident=hashlib.sha256(path.read_bytes()).hexdigest()
        if ident!=state['state_sha256']:
            raise ValueError('Material state mismatches four-coordinate endpoint')
        with np.load(path) as z:
            self.q=z['shape_params'].copy()
            self.chi=z['order'].copy()
            if self.chi.shape!=self.seed.shape or not np.array_equal(z['registration'],self.co):
                raise ValueError('Mismatched reference registration or material mesh')
        x=np.r_[self.pack(self.q),0.]
        emit(kind='input',four_coordinate_endpoint=ident,order_mesh=self.n,
             carrier_mesh=self.mesh,polar=self.polar,mode='epsilon 16 t²(1−t)² physical length')
        self.record(x)
        lim=[(-3.8,1.6),(-7.,2.7),(-2.,2.),(-2.,2.),(-2.,2.)]
        opt=minimize(self.eval,x,jac=True,method='L-BFGS-B',callback=self.record,bounds=lim,
                     options=dict(gtol=.0005,ftol=2e-15,maxiter=350,maxls=28))
        self.record(opt.x)
        lo,hi=np.array(lim).T
        active=np.any((opt.x-lo<1e-5)|(hi-opt.x<1e-5))
        emit(kind='optimization',opt_success=bool(opt.success),message=str(opt.message),
             shape_residual=float(np.max(abs(opt.jac))),active_numerical_limits=bool(active),
             stationary_in_five_coordinates=bool(np.max(abs(opt.jac))<.01 and not active))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('initial',type=Path)
    p.add_argument('folder',type=Path)
    p.add_argument('--order-mesh',type=int,required=True)
    p.add_argument('--carrier-mesh',type=int,required=True)
    p.add_argument('--polar',type=int,required=True)
    a=p.parse_args()
    WaistRelax(a.order_mesh,a.carrier_mesh,a.polar,a.folder).continue_shape(a.initial)
