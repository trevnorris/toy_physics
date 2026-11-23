"""Matched joined-to-sharp radius-four front and signed geometry-work diagnostic."""
import argparse
import contextlib
import gc
import io
from pathlib import Path

import numpy as np

from mouth_exchange_carrier import BIAS
from mouth_joint_collar_release import joined_collar_density
from mouth_internal_collar import CollarPair
from mouth_poly_shallow_release import CUT,source
from mouth_sealed_shape_relax import Meridian
from mouth_supercore_joined_front import Joined,emit,file_sha,registered_density,state_from
from mouth_supercore_shape import preparation
from mouth_supercore_terminal_order import terminal
from mouth_wound_carrier import Wound


class Front:
    def __init__(self,path,grid,winding):
        self.sh,self.c=state_from(path)
        self.grid=grid
        self.N=winding
        with np.load(path) as data:
            self.co=data['registration']
            self.labels=data['labels'] if 'labels' in data else None
            if self.labels is not None:
                self.na=(self.labels.shape[1]-1)//2
                self.nt=self.labels.shape[0]-1
                self.outer=data['outer_order']
        self.ref,_,self.co0=preparation(12)
        self.order_mesh=self.c.shape[1]-1

    def joined(self,sh,precondition=True):
        if self.labels is None:
            density=registered_density(sh,self.ref,self.co,self.co0)
            def outside(w,x):
                return source(sh,self.co,w/sh.length,x,True)
        else:
            pair=CollarPair(sh,self.co,self.na,self.nt)
            density=joined_collar_density(sh,self.co,self.ref,pair,self.labels)
            t=np.linspace(0,1,2*self.order_mesh+1)
            x=np.linspace(0,1,self.order_mesh+1)
            outside_work=pair.evaluate(self.labels,(t,x,self.c),((1-np.cos(np.pi*t))/2,.10*x+.90*x*x,self.outer))['outer_source']
            def outside(w,x):
                return outside_work
        profile=terminal(sh,self.order_mesh,CUT,outside)
        return Joined(sh,self.c,self.N,self.grid,elastic_density=density,precondition=precondition,terminal_profile=profile)

    def sharp(self,shape):
        return Wound(shape,self.c,self.N,self.grid[0],self.grid[1],28.,22.)


def run(args):
    base=Front(args.state,tuple(args.grid),args.winding)
    sh=base.sh
    emit(kind='preparation',source_sha256=file_sha(args.state),grid=args.grid,winding=args.winding,
         state=args.state.as_posix(),shape=sh.params.tolist(),depth=sh.length,exit_radius=sh.exit,
         released_archive_labels=base.labels is not None)
    jm=base.joined(sh)
    u,read=jm.run_state()
    u0,rd=jm.run_state(True)
    _,_,n,_,vel,t=jm.evaluate(u,True)
    t-=jm.evaluate(u0,True)[-1]
    _,r,a,z,_,_,s,cut,_,area,_=jm.D.T
    mach=np.sqrt(np.sum(vel**2,axis=1)+(base.N/s)**2)/n**2
    for sign in (-1,1):
        c=cut==sign; dense=c&(n>.5)&(mach<1)
        emit(kind='drive_minus_rest',sign=sign,axial=area[c]@t[c],dense=area[dense]@t[dense])
    sharp=base.sharp(sh)
    with contextlib.redirect_stdout(io.StringIO()):
        v,p,old=sharp.solve(BIAS)
    emit(kind='coupled_difference',joined_minus_sharp=read['excess_energy']-old['excess_energy'],
         joined_current=read['feed_current'],sharp_current=old['total_current'],sharp_energy=old['excess_energy'])
    # The constrained variables v, phase, delta-order are stationary at the center.
    # Match both energies at identical fixed pulls of their own material coordinates.
    del jm,sharp
    gc.collect()
    if args.save:
        np.savez_compressed(args.save,fields=u,rest=u0,sharp_amplitude=v,sharp_phase=p)
    if args.work:
        scales=np.array([1.,.25,.25,.25,.25])
        for axis in range(5):
            for h in ([.0004,.0002] if axis==1 else [.0004]):
                rows=[]
                for sign in (-1,1):
                    dq=sign*h*scales*np.eye(5)[axis]
                    new=Meridian(4.,sh.params+dq[:4],extra=sh.extra+dq[4])
                    m=base.joined(new,False)
                    live=m.evaluate(u)[0]
                    oldvar=base.sharp(new).amplitude(v,p)[0]
                    rows.append((live,oldvar))
                    del m
                forces=-(np.array(rows[1])-rows[0])/(2*h)
                emit(kind='shape_force',coordinate=axis,step=h,joined=forces[0],sharp=forces[1],
                     correction=forces[0]-forces[1],signed_energies=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path)
    p.add_argument('--grid',type=int,nargs=3,required=True)
    p.add_argument('--winding',type=int,choices=(0,1),required=True)
    p.add_argument('--work',action='store_true')
    p.add_argument('--save',type=Path)
    run(p.parse_args())
