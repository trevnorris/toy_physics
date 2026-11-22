"""Heat the same zero-load 4D material and probe shear and compression over one transit."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse.csgraph import connected_components

from mouth_bulk_isotropic_4d import strain_basis
from mouth_bulk_matched_pressure import diffraction
from mouth_bulk_relaxed_shape import PeriodicMaterial
from mouth_bulk_thermal_cell import ThermalCell, library
from neutral_patch_sites import emit

COLD = Path('_scratch/s11c/iso4_pressure_glass87_02/native/endpoint.npz')
TEMPERATURES = (.008, .064, .256, .512, 1.024)
TRANSIT = 48.88171243122866
FRAMES = 48
SHEAR_STEP = .004
PRESSURE_STEP = .25


def verified_state(path):
    area = path.parent.parent
    out,child,limits = (json.loads((area/name).read_text()) for name in
                        ('outcome.json','child-outcome.json','limit-validation.json'))
    if not (out['limitsVerified'] and limits['verified'] and out['exitCode']==child['exitCode']==0
            and out['childOutcome']==child and child['guardReason'] is None and int(limits['actual']['memory.swap.max'])==0):
        raise ValueError(f'Input state missing its successful guard: {path}')


def cold_data(path):
    verified_state(path)
    with np.load(path) as saved:
        y,box,rc,tangent = [saved[name].copy() for name in ('y','cell','rc','tangent')]
        assert int(saved['seed'])==87
    k=tangent[0,0]
    shear=tangent[1:,1:]-np.outer(tangent[1:,0],tangent[0,1:])/k
    vals,vectors=np.linalg.eigh(shear)
    axis=vectors[:,0]
    if axis[np.argmax(abs(axis))]<0:
        axis=-axis
    drive=np.einsum('q,qij->ij',axis,strain_basis()[1:])
    np.testing.assert_allclose([np.trace(drive),np.trace(drive@drive)],[0.,.5],atol=1e-12)
    return y@box,box,float(rc),drive,float(k),float(vals[0])


class Sample:
    def __init__(self,lib,pos,box,cutoff):
        self.lib,self.pos,self.box,self.cutoff=lib,pos.copy(),box.copy(),cutoff
        self.reset_cell()

    def reset_cell(self):
        self.volume=np.linalg.det(self.box)
        self.mat=ThermalCell(self.lib,self.box,len(self.pos),self.cutoff)

    def read(self,temp,axis):
        energy,grad,vir=self.mat.read(self.pos)
        density=len(self.pos)/self.volume
        pressure=density*temp-np.trace(vir)/(4*self.volume)
        tau=np.sum(axis*vir)/self.volume
        return dict(density=float(density),pressure=float(pressure),shear=float(tau),
                    energy_per_site=float(energy/len(self.pos)),ideal_pressure=float(density*temp))

    def move_volume(self,temp,axis,pext,mobility,interval,innovation):
        reading=self.read(temp,axis)
        # State measure in q=ln V at fixed fractional coordinates has factor
        # V^(N+1). The ideal part N*T/V is already in the internal pressure.
        step=mobility*interval*(self.volume*(reading['pressure']-pext)+temp)+innovation
        if not np.isfinite(step) or abs(step)>.15:
            raise ArithmeticError(f'Unresolved pressure integration increment {step}')
        scale=np.exp(step/4)
        self.pos*=scale
        self.box*=scale
        self.reset_cell()


class Contacts:
    def __init__(self,sample):
        self.metric=PeriodicMaterial(len(sample.pos),1.30)
        self.start=sample.pos@np.linalg.inv(sample.box)
        self.start_box=sample.box.copy()
        self.initial=self.edges(sample)[0]
        self.current=self.initial.copy()
        self.always=self.initial.copy()

    def edges(self,sample):
        y=sample.pos@np.linalg.inv(sample.box)
        i,j,offset=self.metric.pairs(y,sample.box)
        d=(y[j]-y[i]+offset)@sample.box
        r2=np.einsum('ij,ij->i',d,d)
        n=len(y)
        def graph(cut):
            keep=(i!=j)&(r2<cut*cut)
            ret=np.zeros((n,n),dtype=bool)
            ret[i[keep],j[keep]]=True
            return ret
        return graph(1.15),graph(1.30)

    def read(self,sample):
        new,keep=self.edges(sample)
        self.current=(self.current&keep)|new
        self.always &= self.current
        _,parts=connected_components(self.current,directed=False)
        delta=(sample.pos@np.linalg.inv(sample.box)-self.start)@self.start_box
        delta-=delta.mean(axis=0)
        return dict(initial=int(self.initial.sum()),contacts=int(self.current.sum()),
                    original_missing=int((self.initial&~self.current).sum()),
                    ever_missing=int((self.initial&~self.always).sum()),new=int((self.current&~self.initial).sum()),
                    largest_cluster=int(np.bincount(parts).max()),
                    rms_nonaffine=float(np.sqrt(np.mean(np.sum(delta*delta,axis=1)))))


def run(args):
    lib=library(args.work)
    x0,b0,cut,axis,k0,g0=cold_data(COLD)
    source=args.state
    if source:
        verified_state(source)
        with np.load(source) as saved:
            x,box=saved['pos'].copy(),saved['box'].copy()
            prev=int(saved['rung'])
            assert int(saved['history'])==args.history and str(saved['phase'])=='prepare'
            assert prev==args.rung-(args.phase=='prepare')
            np.testing.assert_allclose(saved['axis'],axis,atol=1e-12)
    else:
        assert args.rung==0 and args.phase in ('prepare','cold')
        x,box=x0,b0
    temp=0. if args.phase=='cold' else TEMPERATURES[args.rung]
    if args.phase=='shear':
        drives=[(s*SHEAR_STEP*(.5 if args.half else 1.),0.) for s in (1,-1)]
    elif args.phase=='pressure':
        drives=[(0.,s*PRESSURE_STEP*(.5 if args.half else 1.)) for s in (1,-1)]
    elif args.phase=='cold':
        drives=[(s*.0001,0.) for s in (1,-1)]+[(0.,s*.01) for s in (1,-1)]
    else:
        drives=[(0.,0.)]
    samples=[Sample(lib,x@expm(shift*axis),box@expm(shift*axis),cut) for shift,_ in drives]
    tracking=Contacts(samples[0]) if len(samples)==1 else None
    # The cell pressure response time is an imposed instrumental setting. Each
    # history shares the same noise across its undriven and signed branches.
    mobility=1/(np.linalg.det(b0)*k0*args.pressure_time)
    duration=8. if args.phase=='cold' else TRANSIT
    nframe=8 if args.phase=='cold' else FRAMES
    per_frame=16*args.refine
    interval=duration/nframe/per_frame
    base_steps=int(np.ceil((duration/nframe/16)/.00016))
    steps=base_steps  # halving the block time also halves the microscopic time step
    dt=interval/steps
    # Distinct independent draws across temperatures and independent histories.
    seed=[9017,args.history,args.rung,int(args.phase!='prepare')]
    rng=np.random.default_rng(seed)
    bridge=np.random.default_rng([*seed,2129])
    cell_rng=np.random.default_rng([*seed,2317])
    cell_bridge=np.random.default_rng([*seed,951])
    emit(test='setup',phase=args.phase,rung=args.rung,history=args.history,temperature=temp,refine=args.refine,
         applied_pressures=[p for _,p in drives],shears=[s for s,_ in drives],dt=dt,
         barostat_interval=interval,pressure_time=args.pressure_time,volume_mobility=mobility,
         transit=duration,sites=len(x),cutoff=cut,cold_shear=g0,cold_K_shape_fixed=k0,
         axis=axis.tolist(),source=str(source or COLD),source_sha256=hashlib.sha256((source or COLD).read_bytes()).hexdigest())
    if args.refine not in (1,2):
        raise ValueError('Only common-noise coarse and fine pairs supported')
    for frame in range(nframe+1):
        if frame:
            for block in range(16):
                coarse_dt=dt*args.refine
                w=rng.normal(size=(base_steps,len(x),4))*np.sqrt(2*temp*coarse_dt)
                q=cell_rng.normal()*np.sqrt(2*mobility*temp*interval*args.refine)
                if args.refine==2:
                    z=bridge.normal(size=w.shape)*np.sqrt(2*temp*coarse_dt)
                    ww=np.stack(((w+z)/2,(w-z)/2),axis=1).reshape(2*base_steps,len(x),4)
                    dq=cell_bridge.normal()*np.sqrt(2*mobility*temp*interval*2)
                    batches=[(ww[:base_steps],(q+dq)/2),(ww[base_steps:],(q-dq)/2)]
                else:
                    batches=[(w,q)]
                for noise,innovation in batches:
                    for sample,(_,pext) in zip(samples,drives):
                        sample.mat.advance(sample.pos,noise,dt)
                        sample.move_volume(temp,axis,pext,mobility,interval,innovation)
        readings=[sample.read(temp,axis) for sample in samples]
        row=dict(test='frame',frame=frame,elapsed=duration*frame/nframe,readings=readings)
        if tracking:
            row.update(tracking.read(samples[0]))
            if frame in (0,nframe):
                s=samples[0]
                row['diffraction']=diffraction(s.pos@np.linalg.inv(s.box),s.box)['peak']
        if args.phase in ('shear','cold'):
            row['shear_slope']=(readings[0]['shear']-readings[1]['shear'])/(drives[0][0]-drives[1][0])
        if args.phase in ('pressure','cold'):
            row['log_volume_compliance']=-np.log(samples[-2].volume/samples[-1].volume)/(drives[-2][1]-drives[-1][1])
        emit(**row)
    if args.phase=='cold':
        np.testing.assert_allclose([row['shear_slope'],row['log_volume_compliance']],[g0,1/k0],rtol=.005)
    if len(samples)==1:
        np.savez(args.work/'end.npz',pos=samples[0].pos,box=samples[0].box,axis=axis,phase=args.phase,
                 rung=args.rung,history=args.history,temperature=temp)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase',choices=('cold','prepare','undriven','shear','pressure'))
    p.add_argument('--rung',type=int,choices=range(len(TEMPERATURES)),default=0)
    p.add_argument('--history',type=int,default=0)
    p.add_argument('--state',type=Path)
    p.add_argument('--work',type=Path,required=True)
    p.add_argument('--half',action='store_true')
    p.add_argument('--refine',type=int,choices=(1,2),default=1)
    p.add_argument('--pressure-time',type=float,default=.5)
    run(p.parse_args())
