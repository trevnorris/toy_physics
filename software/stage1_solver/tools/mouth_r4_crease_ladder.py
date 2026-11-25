"""Read-only step-166 work and separate moving-cut directional derivatives on native quadrature rays."""
import argparse
import ctypes as ct
import json
from pathlib import Path
import subprocess

import numpy as np

from mouth_cp10_label_gauge import VERTICES
from mouth_r4_aligned_crystal import cell_call
from mouth_r4_aligned_derivatives import load
from mouth_r4_full_label_precision import HPData, HPGrid, P
from mouth_r4_smooth_label_derivatives import CrystalBanks, LD
from mouth_r4_smooth_static import RULE
from mouth_r4_stage2_outcome import digest, read


class GaugeRay:
    def __init__(self, geometry, labels, direction, order):
        ids=geometry.ids
        self.y=labels.reshape(-1,3)[ids]
        self.d=direction.reshape(-1,3)[ids]
        self.width=geometry.widths
        offsets=(labels-geometry.base).reshape(-1,3)[ids].transpose(0,2,1)
        self.weight=np.zeros(len(ids),LD)
        self.first=self.weight.copy()
        self.second=self.weight.copy()
        for q,vol,grad in zip(RULE,geometry.volume,geometry.physical_gradients):
            N=np.prod(np.where(VERTICES,q,1-np.array(q)),axis=1)
            chi=order.ravel()[ids]@N
            weight=vol*((1-chi)**2*(1+2*chi))
            du=self.d.transpose(0,2,1)@grad
            self.weight+=weight
            self.first+=weight*np.sum((offsets@grad)*du,axis=(1,2))
            self.second+=weight*np.sum(du*du,axis=(1,2))

    def evaluate(self,a):
        work,_=cell_call('signed_work',self.y,self.width,a*self.d)
        _,grad=cell_call('values',self.y+a*self.d,self.width,gradient=True)
        return (a*(self.first+a*self.second/2)+self.weight*work,
                self.first+a*self.second+self.weight*np.sum(grad*self.d,axis=(1,2)))


def main(root,out):
    out.mkdir(parents=True,exist_ok=False)
    dirname=root/'directions'
    rec=read(dirname/'result.json')
    cross=read(root/'crossings/result.json')
    if (digest(dirname/'directions.npz')!=rec['directions_sha256'] or
        rec['directions_sha256']!=cross['directions_sha256'] or rec['optimization_steps']!=0):
        raise ValueError('Frozen direction provenance mismatch')
    with np.load(dirname/'directions.npz') as source:
        y,d,chi=(np.ascontiguousarray(source[k],LD) for k in ('labels','pj','order'))
    _,trial,_,_=load(root.parent/'r4_rematch39_01','primary',out)
    geometry=trial.smooth_order_geometry
    active=np.max(chi.ravel()[geometry.ids],axis=1)>0
    dims=tuple(np.array(chi.shape)-1)
    banks=CrystalBanks(trial.b,chi,out)
    header=Path(__file__).with_suffix('.hpp').resolve()
    cpp=out/'ray.cpp'
    cpp.write_text('#include "audit.cpp"\n#include "'+str(header)+'"\n')
    library=out/'ray.so'
    subprocess.run(['g++','-std=gnu++17','-O3','-shared','-fPIC',str(cpp),'-o',str(library)],check=True)
    fn=ct.CDLL(str(library)).crease_ray_work
    fn.argtypes=[ct.POINTER(HPData),P,P,P,HPGrid,P,ct.c_longdouble,P]
    fn.restype=ct.c_int
    def ptr(array):
        return array.ctypes.data_as(P)
    gauge=GaugeRay(geometry,y,d,chi)
    def measure(a):
        ray=np.zeros((len(active),3,3,3),LD)
        rc=fn(ct.byref(banks.data),ptr(y),ptr(d),ptr(chi),banks.rule,ptr(banks.weights),a,ptr(ray))
        if rc:raise ArithmeticError(f'Invalid ray evaluation at {a}: {rc}')
        unit,slope=gauge.evaluate(LD(a))
        return ray,LD(1e-8)*unit,LD(1e-8)*slope
    first=LD(cross['alpha_cross'])
    cap=LD(2)**-13
    probe=[(r['multiple'],LD(r['alpha'])) for r in cross['ladder']]
    for v in (2.,4.):
        if not any(key==v for key,_ in probe):probe.append((v,first*v))
    for v in (0.99999,1.00001):
        probe.append((v,first*v))
    left,right=sorted((v*first for v in (0.99999,1.00001)))
    eps=LD(1e-8)
    zero,_,gz=measure(0.)
    E0=zero[...,0].sum((1,2))
    ref=banks.energy(y,active)
    native0=float(np.max(abs(ref-E0)))
    if native0>1e-14:raise ArithmeticError(f'Physical baseline mismatch {native0}')
    cells=[tuple(i) for i in cross['first']['cells']]
    dual=tuple(cross['first']['ray'])
    common=[np.ravel_multi_index(c,dims) for c in cells]
    def slope_ledger(ray,gs):
        ordinary,cut=ray[...,1].sum(),ray[...,2].sum()
        return dict(ordinary=float(ordinary),moving_cut=float(cut),gauge=float(gs.sum()),
                    total=float(ordinary+cut+gs.sum()))
    base=slope_ledger(zero,gz)
    if abs(base['total']-rec['slopes']['pj'])>1e-10:
        raise ArithmeticError(f'Base projected directional derivative {base} != {rec["slopes"]["pj"]}')
    native_g=LD(rec['slopes']['pj'])
    rows=[]
    raw={}
    perray={}
    for key,a in sorted(probe,key=lambda v:v[1]):
        ray,gw,gs=measure(a)
        # Subtract each paired ray before accumulating across cells.
        delta_rau=ray[...,0]-zero[...,0]
        physical=delta_rau.sum((1,2))
        dE=physical+gw
        dc=dE.sum()
        rows.append(dict(multiple=key,alpha=float(a),delta=float(dc),physical=float(physical.sum()),
                    gauge=float(gw.sum()),secant=float(dc/a),relative_error=float(abs(dc/a-native_g)/abs(native_g)),
                    directional=slope_ledger(ray,gs)))
        quantities=np.stack([dE,physical,gw,ray[...,1].sum((1,2)),ray[...,2].sum((1,2)),gs],axis=1)
        raw[float(a)]=quantities
        perray[float(a)]=[dict(cell=c,physical_energy_change=float(delta_rau[(idx,)+dual]),
                         ordinary=float(ray[(idx,)+dual+(1,)]),moving_cut=float(ray[(idx,)+dual+(2,)]))
                         for c,idx in zip(cells,common)]
        if key is None or key==.001:
            fd=banks.energy(np.ascontiguousarray(y+a*d),active)-ref
            rows[-1]['independent_signed_energy_minus_ray_sum']=float((fd-physical).sum())
        print(json.dumps(rows[-1]),flush=True)
    iwide=min(range(len(rows)),key=lambda i:abs(rows[i]['alpha']-float(cap)))
    if abs(rows[iwide]['delta']-rec['terminal_banks']['pj']['archived_delta'])>1e-13:
        raise ArithmeticError('Extended wide work differs from archive')
    wide=raw[float(cap)]
    peaks=np.argsort(-wide[:,0])[:10]
    ids=set(int(i) for i in peaks)|set(common)
    details=[]
    for idx in sorted(ids,key=lambda i:-wide[i,0]):
        def ds(col):
            return dict(ordinary=float(col[3]),moving_cut=float(col[4]),gauge=float(col[5]),
                        total=float(col[3:].sum()))
        basecell=np.array([0.,0.,0.,zero[idx,...,1].sum(),zero[idx,...,2].sum(),gz[idx]],LD)
        details.append(dict(cell=list(map(int,np.unravel_index(idx,dims))),wide_delta=float(wide[idx,0]),
                   rank=int(np.count_nonzero(wide[:,0]>wide[idx,0])+1),slope_at_zero=ds(basecell),
                   left=ds(raw[float(left)][idx]),right=ds(raw[float(right)][idx]),wide=ds(wide[idx])))
    npz=out/'cells.npz'
    alphas=sorted(raw)
    np.savez_compressed(npz,alpha=np.array(alphas),cells=np.stack([raw[a] for a in alphas]),shape=dims,
                    baseline=np.stack([E0,zero[...,1].sum((1,2)),zero[...,2].sum((1,2)),gz],axis=1))
    rr=dict(step=166,epsilon=float(eps),optimizations=0,state_sha256=cross['state_sha256'],
            directions_sha256=rec['directions_sha256'],cross_sha256=digest(root/'crossings/result.json'),
            alpha_cross=float(first),base=base,native_slope=float(native_g),base_cell_energy_max_error=native0,
            slope_probes=[float(left),float(right)],ray=dual,crossing_cells=cells,
            rows=rows,leading_cells=details,first_endpoint_rays={str(a):perray[a] for a in alphas},
            cells_archive_sha256=digest(npz),columns=['delta','physical_delta','gauge_delta',
            'physical_ordinary_slope','physical_moving_cut_slope','gauge_slope'])
    (out/'result.json').write_text(json.dumps(rr,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('output',type=Path)
    a=p.parse_args()
    main(a.root.resolve(),a.output.resolve())
