"""Promoted signed crystal quadrature and a narrower depth pair on the fixed R4 state."""
import argparse
import copy
import json
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from mouth_r4_smooth_label_derivatives import CrystalBanks
from mouth_r4_smooth_static import EXPECTED, Gauge, LD, digest, frozen


def run(area, first, output):
    native=json.loads(first.read_text())
    if native['state_sha256']!=EXPECTED['held_stop']:
        raise ValueError('Mismatched frozen shape-force measurements')
    with TemporaryDirectory() as directory:
        tmp=Path(directory)
        d,t,fields,labels,order=frozen(area,tmp)
        from mouth_cp10_field_shape import Bath
        from mouth_cp10_label_gauge import depth_grid
        from mouth_r4_depth_precision import extended
        banks=CrystalBanks(t.b,order,tmp)
        gauge=Gauge(t,d.reference_depth,order)
        bath=Bath(t,d.base,fields,d.reference_depth)
        field_t=banks.buffers['field_t'].copy()
        original=banks.buffers['shape'].copy()
        active=np.max(gauge.order,axis=1)>0
        # HPData's ctypes scalar accessor otherwise downcasts an 80-bit scalar
        # to Python float before use. Copy these two scalars by raw C ABI.
        import ctypes as ct
        def scalar(key,value):
            val=np.array(value,LD)
            ct.memmove(ct.addressof(banks.data)+getattr(type(banks.data),key).offset,val.ctypes.data,val.nbytes)
        length=LD(t.shape.length)
        base_depth=depth_grid(length,d.reference_depth)
        primary=LD(1e-8)
        secondary=LD(1e-10)
        def evaluate(dx):
            L=length*np.exp(dx[0])
            shift=depth_grid(L,d.reference_depth)-base_depth
            banks.buffers['field_t'][:]=field_t+shift
            scalar('length',L)
            scalar('epsilon',LD(t.polar_shift)+dx[5])
            p=banks.buffers['shape']
            p[:]=original
            p[0]+=(original[0]+LD(1.5))*np.expm1(dx[1]/4)
            p[1:4]+=dx[2:5]/4
            moved=labels.astype(LD)
            moved[...,1]+=shift[None,:,None]
            e=banks.energy(moved,active)
            # The native vortex table and exterior splines are tabulated in
            # double. Accumulate those terms per field; promote the two large
            # reservoir-order energies before any subtraction.
            be=bath.energy(np.asarray(dx,np.float64))
            shape=copy.copy(t.shape)
            shape.params=t.shape.params+dx[:4]*[1,.25,.25,.25]
            shape.extra=t.shape.extra+.25*dx[4]
            eb=copy.copy(bath)
            eb.model=copy.copy(bath.model)
            eb.model.polar_shift=t.polar_shift+dx[5]
            be[3:5]=extended(eb,shape,L,d.reference_depth)
            return e,be,gauge.evaluate(labels,dx)
        base=evaluate(np.zeros(6,LD))
        anch=float(base[0].sum()-native_crystal(area))
        if abs(anch)>1e-11:
            raise ArithmeticError(f'Promoted crystal energy differs by {anch}')
        rows=[]
        for i,unit in enumerate(np.eye(6,dtype=LD)):
            source=native['shapes'][i]
            widths=[LD(2)**-18,LD(2)**-20,LD(2)**-22,LD(2)**-24] if i==0 else [LD(2)**-18]
            samples=[]
            for h in widths:
                nn=evaluate(-h*unit)
                pp=evaluate(h*unit)
                crystal, carrier, reg=[np.sum(a-b,dtype=LD)/(2*h) for a,b in zip(nn,pp)]
                errs={str(float(e)):float((crystal-source['elastic'])+(carrier-source['bath'])
                       +e*(reg-source['unit_gauge'])) for e in (primary,secondary)}
                samples.append(dict(h=float(h),crystal_minus_analytic=float(crystal-source['elastic']),
                              bath_minus_analytic=float(carrier-source['bath']),unit_gauge_minus_analytic=float(reg-source['unit_gauge']),
                              smooth_total_central={str(float(e)):float(crystal+carrier+e*reg) for e in (primary,secondary)},
                              total_signed_errors=errs))
            rows.append(dict(mode=source['mode'],samples=samples))
        record=dict(command=sys.orig_argv,state_sha256=EXPECTED['held_stop'],source_sha256=digest(Path(__file__)),
            physical_bank_source_sha256=digest(Path(__file__).with_name('mouth_r4_smooth_label_derivatives.py')),
            analytic_bank_sha256=digest(first),promoted_crystal_static_delta=anch,rows=rows,
            accuracy='80-bit crystal quadrature with moving cuts, LD field and gauge sums, promoted extensive order terms; the existing carrier interpolation tables stay in native precision.')
        output.write_text(json.dumps(record,indent=2)+'\n')
        print([r['samples'][-1]['total_signed_errors'] for r in rows],flush=True)


def native_crystal(area):
    return json.loads((area/'primary/held_stop.json').read_bytes())['elastic_energy']


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('area','first','output'):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    run(a.area.resolve(),a.first.resolve(),a.output.resolve())
