"""All free scalar-label derivatives of the frozen R4 pointwise smooth gauge and physical energy.

Graph colors permit simultaneous signed perturbations with disjoint element support.
The two sides of the weld retain their native integration cuts. All signed energies
are evaluated in extended precision; ambiguous narrow physical gradients are
resolved in binary128 on the affected stencil only.
"""
import argparse
import ctypes as ct
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from mouth_r4_smooth_static import Gauge, frozen, EXPECTED, digest


LD = np.longdouble
GAUGES = (1e-8,1e-10)
LIMIT = 1e-10


class CrystalBanks:
    def __init__(self, body, order, directory, quad=False):
        from mouth_r4_full_label_precision import HPData, HPGrid, P, promoted
        self.f = promoted(directory, quad)
        self.data = HPData()
        self.buffers = {k:np.ascontiguousarray(v,LD) for k,v in body.buffers.items()}
        self.pointer = P
        for name,dtype in self.data._fields_:
            if name in self.buffers:
                a = self.buffers[name]
                ptr = a.ctypes.data_as(P)
                setattr(self.data,name,HPGrid(len(a),ptr) if dtype==HPGrid else ptr)
            elif dtype != P:
                setattr(self.data,name,getattr(body.data,name))
        self.order = np.ascontiguousarray(order, LD)
        gp, gw = np.polynomial.legendre.leggauss(3)
        self.points = np.ascontiguousarray((1+gp)/2,LD)
        self.weights = np.ascontiguousarray(gw/2,LD)
        self.rule=HPGrid(3,self.points.ctypes.data_as(P))
        self.minimum=np.inf

    def energy(self,labels,selected):
        enabled=np.ascontiguousarray(selected,np.int32)
        value=np.zeros(len(selected),LD)
        if not np.any(enabled):
            return value
        floor=np.zeros(1,LD)
        P=self.pointer
        code=self.f(ct.byref(self.data),labels.ctypes.data_as(P),self.order.ctypes.data_as(P),
                    self.rule,self.weights.ctypes.data_as(P),enabled.ctypes.data,value.ctypes.data,floor.ctypes.data)
        if code or floor[0]<=0:
            raise ValueError(f'Invalid signed crystal energy {code}, orientation {floor[0]}')
        self.minimum=min(self.minimum,floor[0])
        return value


def checksum_array(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def run(area,out,c):
    out.mkdir(exist_ok=False,parents=True)
    driver,trial,_,labels,order=frozen(area,out)
    from mouth_r4_full_label_work import coloring
    from mouth_r4_eulerian_relax_labels import LabelRelaxation
    from mouth_r4_eulerian_solid import SpatialCrystal
    labels=np.ascontiguousarray(labels,LD)
    opt=LabelRelaxation(trial.b,np.asarray(labels,float),order,1e-8)
    cm=opt._components()[c]
    physical,pg,_,_=SpatialCrystal(trial.b).evaluate(labels,order)
    smooth=Gauge(trial,driver.reference_depth,order)
    zero=np.zeros(6,LD)
    unit,ug=smooth.evaluate(labels,zero,gradient=True)
    scalar=np.zeros((2,len(cm.held)),LD)
    for k,v in enumerate([pg,ug]):
        np.add.at(scalar[k], cm.groups,v[...,c].ravel())
    ids=smooth.ids
    group=cm.groups[ids]
    colors=coloring(cm,ids)
    active=opt.active.ravel()
    banks=CrystalBanks(trial.b,order,out)
    anchor=banks.energy(labels,active)
    delta=float(anchor.sum()-physical)
    if abs(delta)>1e-12:
        raise ValueError(f'Promoted physical anchor differs by {delta}')
    full_coverage=~cm.held
    def fd_engine(mask,steps,crystal):
        fd=np.zeros((len(steps),2,len(cm.held)),LD)
        for color in np.unique(colors[mask]):
            bank=mask & (colors==color)
            touched=bank[group]
            target=np.max(np.where(touched,group,-1),axis=1)
            selected=target>=0
            if np.any(touched & (group!=target[:,None])):
                raise ValueError('Signed stencils overlap a physical cell')
            moving=bank[cm.groups].reshape(order.shape)
            for k,h in enumerate(steps):
                physical_energy=[]; gauge_energy=[]
                for sign in (-1,1):
                    shifted=labels.copy()
                    shifted[...,c][moving]+=sign*h
                    if any(shifted[...,i].min()<ax[0]-2e-13 or shifted[...,i].max()>ax[-1]+2e-13
                           for i,ax in enumerate(opt.axes)):
                        raise ValueError('Signed bank left registered domain')
                    physical_energy.append(crystal.energy(shifted,selected & active)[selected])
                    gauge_energy.append(smooth.evaluate(shifted,zero,enabled=selected))
                for term,(minus,plus) in enumerate([physical_energy,gauge_energy]):
                    np.add.at(fd[k,term],target[selected],(plus-minus)/(2*h))
        return fd

    broad_steps=tuple(LD(2)**(-p) for p in (24,25))
    narrow_steps=tuple(LD(2)**(-p) for p in (36,37))
    broad=fd_engine(full_coverage,broad_steps,banks)
    def error(values,eps):
        return values[:,0]-scalar[0]+eps*(values[:,1]-scalar[1])
    unresolved=full_coverage & np.logical_or.reduce(
        [np.max(np.abs(error(broad,eps)),axis=0)>LIMIT/2 for eps in GAUGES])
    best=broad.copy()
    if np.any(unresolved):
        hp=out/'narrow';hp.mkdir()
        quad=CrystalBanks(trial.b,order,hp,quad=True)
        # Establish one binary128 baseline for all cells requested below.
        quad.energy(labels,active & np.any(unresolved[group],axis=1))
        narrow=fd_engine(unresolved,narrow_steps,quad)
        best[:,:,unresolved]=narrow[:,:,unresolved]
    else:
        narrow=None
    sectors={'active':cm.free,'fully_melted':cm.dead}
    if not np.array_equal(cm.free|cm.dead,full_coverage) or np.any(cm.free&cm.dead):
        raise ArithmeticError('Label partition is incomplete')
    info=dict(component=c, state_sha256=EXPECTED['held_stop'], script_sha256=digest(Path(__file__)),
              smooth_evaluator_sha256=digest(Path(__file__).with_name('mouth_r4_smooth_static.py')),
              command=sys.orig_argv, scalar_label_count=int(full_coverage.sum()),colors=int(colors.max()+1),
              active_count=int(cm.free.sum()),free_melt_count=int(cm.dead.sum()),narrow_refinements=int(unresolved.sum()),
              physical_energy_difference_from_native=delta, unscaled_smooth_gauge_energy=float(unit.sum()),
              native_precision_epsilon=float(np.finfo(LD).eps),broad_widths=[float(x) for x in broad_steps],
              narrow_binary128_physical_widths=[float(x) for x in narrow_steps], broad_minimum_ordered_chart=float(banks.minimum),
              status={}, tolerance=LIMIT)
    for eps in GAUGES:
        total=error(best,eps)
        entries={}
        for name,mask in sectors.items():
            groups=np.flatnonzero(mask)
            worst=groups[np.argmax(abs(total[-1,groups]))]
            entries[name]=dict(max_abs_full=float(np.abs(total[0,mask]).max(initial=0)),
                               max_abs_half=float(np.abs(total[1,mask]).max(initial=0)),
                               half_minus_full_peak=float(abs(total[1,mask]-total[0,mask]).max(initial=0)),
                               largest_error_group=int(worst),peak_node=list(map(int,np.unravel_index(cm.representative[worst],order.shape))),
                               analytic=float(scalar[0,worst]+eps*scalar[1,worst]),
                               direct_fd=float(best[-1,0,worst]+eps*best[-1,1,worst]),
                               max_abs_uncorrected_broad_error=float(abs(error(broad,eps)[:,mask]).max(initial=0)))
        info['status'][str(eps)]=entries
    info['passed']=all(row['max_abs_full']<=LIMIT and row['max_abs_half']<=LIMIT
                       for sector in info['status'].values() for row in sector.values())
    npz=out/'raw_signed.npz'
    np.savez_compressed(npz,analytic=scalar, free=cm.free, melted=cm.dead,
                        independently_narrowed=unresolved, broad_fd=broad,final_fd=best,
                        scalar_source=cm.representative,colors=colors)
    info['raw_sha256']=digest(npz)
    (out/'record.json').write_text(json.dumps(info,indent=2,allow_nan=False)+'\n')
    print(json.dumps(info),flush=True)


def collect(area,path):
    from mouth_r4_stage2_outcome import verify_guard
    rows=[]
    for c in range(3):
        verify_guard(area/f'guard_{c}',0)
        row=json.loads((area/f'component_{c}/record.json').read_bytes())
        if row['component']!=c or row['raw_sha256']!=digest(area/f'component_{c}/raw_signed.npz'):
            raise ValueError('Incomplete coordinate bank')
        rows.append(row)
    result=dict(state=EXPECTED['held_stop'], complete_free_coordinates=sum(r['scalar_label_count'] for r in rows),
                passed=all(r['passed'] for r in rows),components=rows,
                convention='Frozen Q1 smoothstep after interpolation. Energy signs evaluated per cell for joined free DOFs. Native physical derivative is compared to independently integrated long double/quad-precision physical energy; the gauge derivative is checked against long-double energy differences. Static result only.')
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('complete_free_coordinates','passed')}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('area',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--component',type=int,choices=(0,1,2))
    a=p.parse_args()
    if a.component is None:
        collect(a.area.resolve(),a.output.resolve())
    else:
        run(a.area.resolve(),a.output.resolve(),a.component)
