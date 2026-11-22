"""Finite-temperature common shear on the same converted patch with its brace removed."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_brace_grip import Material, clock, evolve_interval
from mouth_middle_tag_grip import cubic_trace
from neutral_patch_open_w import mix
from neutral_patch_sites import BETWEEN, TEMPERATURE, WITHIN, emit, morse


def shear_work(cell,pos,step):
    """Energy and deviatoric stress under a virtual area-preserving strain."""
    d=pos[cell.j]-pos[cell.i]
    d[:,:2]+=cell.shifts-cell.box*np.rint(d[:,:2]/cell.box)
    d[:,0]*=np.exp(step)
    d[:,1]*=np.exp(-step)
    r=np.linalg.norm(d,axis=1)
    use=r<cell.cutoff
    d,r=d[use],r[use]
    f=mix(d[:,2],cell.power)[0]
    ea,ga=morse(r,*WITHIN)
    eb,gb=morse(r,*BETWEEN)
    derivative=(ga+f*(gb-ga))*(d[:,0]**2-d[:,1]**2)/r
    return np.sum(ea+f*(eb-ea)),derivative.sum()/cell.area


def checked_stress(mat,pos):
    cell=mat.cell
    e,stress=shear_work(cell,pos,0.)
    np.testing.assert_allclose(cell.evaluate(pos)[0],e,atol=2e-9)
    delta=2.e-7
    numerical=(shear_work(cell,pos,delta)[0]-shear_work(cell,pos,-delta)[0])/(2*delta*cell.area)
    np.testing.assert_allclose(numerical,stress,atol=1.e-5,rtol=2e-6)
    return stress


class Strained:
    def __init__(self,mat,gamma,initial):
        self.mat=mat
        stretch=np.exp([gamma,-gamma])
        self.box=mat.cell.box*stretch
        self.image_shifts=mat.shifts*stretch
        self.energy_shifts=mat.cell.shifts*stretch
        self.pos=initial.copy()
        self.pos[:,:2]*=stretch

    def select(self):
        self.mat.cell.box=self.box
        self.mat.cell.shifts=self.energy_shifts
        self.mat.shifts=self.image_shifts


def run(args):
    args.work.mkdir(parents=True,exist_ok=True)
    depths,times,chi,_=clock(args.mesh)
    mat=Material(args)
    initial_grip=mat.validate()
    initial=mat.flat.copy()
    perturb=np.random.default_rng(args.seed+11).normal(size=initial.shape)
    initial+=.001*(perturb-perturb.mean(axis=0))
    rng=np.random.default_rng(args.seed)
    bridge=np.random.default_rng(args.seed+7013)
    for t0,t1 in zip(times[:-1],times[1:]):
        evolve_interval(mat,[initial],t0,t1,[chi],rng,bridge,args.refine)
    archive=Path(f'_scratch/s11c/brace_grip_p{args.power}s{args.seed}_01')
    outcome=json.loads((archive/'outcome.json').read_text())
    assert outcome['limitsVerified'] and outcome['exitCode']==outcome['childOutcome']['exitCode']==0
    reference,=[r for r in map(json.loads,(archive/'stdout').read_text().splitlines())
                if r['test']=='passage' and r['depth']==2. and r['history']=='converting']
    measured=mat.slope(initial)/initial_grip
    np.testing.assert_allclose(measured,reference['grip_tensor_over_ordered'],rtol=2e-10,atol=2e-10)
    trace,_=cubic_trace(args.mesh)
    transit=trace.transit(depths[-1],trace.edges[-1],8)
    amplitudes=np.array([0.,1.e-4,-1.e-4,.5e-4,-.5e-4])
    clones=[Strained(mat,a,initial) for a in amplitudes]
    # Each clone's native and Python image shifts use its own deformed box.
    emit(test='bulk_setup',seed=args.seed,power=args.power,temperature=mat.temperature,mobility=mat.mobility,
         enforced_order=0.,transit=transit,starting_depth=2.,restored_grip_error=np.max(abs(measured-reference['grip_tensor_over_ordered'])),
         imposed_log_strains=amplitudes.tolist(),stress='sigma_xx_minus_sigma_yy in patch force per projected area')
    rng=np.random.default_rng(args.seed+2103)
    nread=64
    for j in range(nread+1):
        if j:
            count=int(np.ceil(transit/nread/.001))
            dt=transit/(nread*count)
            for start in range(0,count,512):
                steps=min(512,count-start)
                noise=rng.normal(size=(steps,mat.n,3))*np.sqrt(2*mat.temperature*mat.mobility*dt)
                zero=np.zeros(steps)
                for clone in clones:
                    clone.select()
                    mat.step(clone.pos,dt,zero,noise)
        sample=[]
        for clone in clones:
            clone.select()
            sample.append(checked_stress(mat,clone.pos) if j in (0,nread) else shear_work(mat.cell,clone.pos,0.)[1])
        full=(sample[1]-sample[2])/(2*amplitudes[1])
        half=(sample[3]-sample[4])/(2*amplitudes[3])
        emit(test='bulk_hold',fraction=j/nread,elapsed=transit*j/nread,stress=sample,
             tangent_full=full,tangent_half=half)
        if not np.isfinite(sample).all():
            raise ArithmeticError('Nonfinite bulk shear stress')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--power',type=int,choices=(1,2),required=True)
    p.add_argument('--seed',type=int,required=True)
    p.add_argument('--work',type=Path,required=True)
    args=p.parse_args()
    args.mesh,args.side,args.cutoff,args.brace,args.temperature,args.mobility,args.refine=34,4,4.5,'elastic',TEMPERATURE,1.,1
    run(args)
