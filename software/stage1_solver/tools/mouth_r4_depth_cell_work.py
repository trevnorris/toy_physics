"""Signed cell-by-cell log-depth work on radius four's fixed step-108 field rematches."""
import argparse
import ctypes as ct
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from mouth_cp10_label_gauge import depth_grid
from mouth_r4_aligned_derivatives import load
from mouth_r4_full_label_precision import HPData, HPGrid, P
from mouth_r4_smooth_label_derivatives import CrystalBanks
from mouth_r4_smooth_static import LD, digest


WIDTHS = 2.**-np.arange(18,21)
SOURCE = {'primary':'d21586ba828f891f67e4601e8a6145a0cf7f25f003273852b159c0094da63837',
          'secondary':'e2c2e2afb9af55e972b57a46b22527a0ba9983fbab173446b17164457564ed10'}


def gradient(out, banks, labels, rate, size):
    # Build an observables-only TU against the very implementation and promoted
    # math used for the scalar energy banks, without editing either source.
    extra = Path(__file__).resolve().with_suffix('.hpp')
    src = out/'depth_work.cpp'
    src.write_text('#include "audit.cpp"\n#include "'+extra.as_posix()+'"\n')
    lib = out/'depth_work.so'
    subprocess.run(['g++','-std=gnu++17','-O3','-shared','-fPIC',str(src),'-o',str(lib)],check=True)
    calc = ct.CDLL(str(lib)).depth_cell_gradient
    calc.argtypes = [ct.POINTER(HPData),P,P,HPGrid,P,P,P]
    calc.restype = ct.c_int
    terms = np.zeros((size,3),dtype=LD)
    points = (labels.ctypes.data_as(P),banks.order.ctypes.data_as(P),banks.rule,
              banks.weights.ctypes.data_as(P),rate.ctypes.data_as(P),terms.ctypes.data_as(P))
    if calc(ct.byref(banks.data),*points) or not np.isfinite(terms).all():
        raise ArithmeticError('Unresolved log-depth work on native cell/face branches')
    return terms


def summary(energy, analytic):
    rows = []
    for i,h in enumerate(WIDTHS):
        em,ep=energy[1+2*i:3+2*i]
        plus=(ep-energy[0])/h
        minus=(energy[0]-em)/h
        cent=(ep-em)/(2*h)
        error=cent-analytic
        rows.append(dict(h=float(h),energy_derivative_analytic=float(analytic),
                         energy_slope_plus=float(plus),energy_slope_minus=float(minus),
                         energy_slope_central=float(cent),energy_error_plus=float(plus-analytic),
                         energy_error_minus=float(minus-analytic),energy_error_central=float(error),
                         force_error_central=float(-error),effective_second_derivative=float((plus-minus)/h)))
    for a,b in zip(rows,rows[1:]):
        a['halving_ratio']=a['energy_error_central']/b['energy_error_central'] if b['energy_error_central'] else None
    rows[-1]['halving_ratio']=None
    return rows


def run(roots, case, folder):
    out = folder/case
    out.mkdir(parents=True,exist_ok=False)
    study = roots/'r4_rematch39_01'
    raw = study/'derivatives'/(case+'_shape')/'record.json'
    bank = json.loads(raw.read_bytes())
    if case not in SOURCE or bank['state_sha256']!=SOURCE[case] or digest(study/case/'rematched.npz')!=SOURCE[case]:
        raise ValueError('Expected the reviewed held state and native frozen shape record')
    dr,trial,solid,order = load(study,case,out)
    banks=CrystalBanks(trial.b,order,out)
    dims=tuple(n-1 for n in order.shape)
    vertices=[order[i:i+dims[0],j:j+dims[1],k:k+dims[2]] for i,j,k in np.ndindex(2,2,2)]
    active=np.max(vertices,axis=0).ravel()>0
    labels=np.array(dr.labels,LD)
    length=LD(trial.shape.length)
    orig=depth_grid(length,dr.reference_depth)
    tbase=banks.buffers['field_t'].copy()
    # Same native analytic rate as used by shape_work at the recorded shape.
    t=np.asarray(tbase,float)
    rate=np.array(np.where(t*length<length-.12,(t*length/(length-.12))**25,1.)-t,LD)
    terms=gradient(out,banks,labels,rate,len(active))
    analytic=terms.sum(axis=1)
    archived=bank['rows'][0]
    if archived['mode']!='log_depth' or abs(analytic.sum()+archived['elastic'])>5.e-11:
        raise ArithmeticError(f'Cell derivative sum differs from native analytical work: {analytic.sum()+archived["elastic"]}')
    states=[LD(0)]
    for h in WIDTHS:
        states.extend((LD(-h),LD(h)))
    energy=[]
    for probe in states:
        L=length*np.exp(probe)
        dt=depth_grid(L,dr.reference_depth)-orig
        banks.buffers['field_t'][:]=tbase+dt
        val=np.array(L,LD)
        ct.memmove(ct.addressof(banks.data)+HPData.length.offset,val.ctypes.data,val.nbytes)
        moved=labels.copy()
        moved[...,1]+=dt[None,:,None]
        energy.append(banks.energy(moved,active))
    energy=np.array(energy)
    h=WIDTHS[0]
    force_errors=(energy[1]-energy[2])/(2*h)+analytic
    selected=np.flatnonzero(active)[np.argsort(-abs(force_errors[active]),kind='stable')[:10]]
    levels=summary(np.sum(energy,axis=1,dtype=LD),analytic.sum(dtype=LD))
    # Sums of signed changes avoid losing digits before dividing by h.
    reduced=summary(np.sum(energy-energy[0],axis=1,dtype=LD),analytic.sum(dtype=LD))
    if any(abs(r['force_error_central']-a['force_error_central'])>1.e-12 for r,a in zip(reduced,levels)):
        raise ArithmeticError('Total-energy and cell-work paths disagree')
    for d,b in zip(reduced[::2],archived['central']):
        if abs(-d['energy_slope_central']-b['crystal'])>1e-11:
            raise ArithmeticError(f'Cellwise force differs from archived signed work at {d["h"]}')
    top=[dict(cell=list(map(int,np.unravel_index(ix,dims))),levels=summary(energy[:,ix],analytic[ix]),
              baseline_energy=float(energy[0,ix]),analytic_energy_terms=list(map(float,terms[ix]))) for ix in selected]
    ftotal=archived['central'][0]['minus_analytic']
    cumulative={}
    for number in (1,5,10):
        values=force_errors[selected[:number]]
        cumulative[str(number)]=dict(signed_force_error=float(values.sum()),relative_to_total=float(values.sum()/ftotal),
                                     absolute_concentration=float(np.abs(values).sum()/np.abs(force_errors[active]).sum()))
    path=out/'cells.npz'
    with path.open('xb') as f:
        np.savez_compressed(f,energy=energy,steps=np.array(states),active=active,dims=dims,
                            analytic_energy_gradient=analytic,analytic_energy_terms=terms,
                            ranked_top_indices=selected)
    result=dict(source_state_sha256=SOURCE[case],source_shape_sha256=digest(raw),raw_sha256=digest(path),
                source_files={p.name:digest(p) for p in [Path(__file__),Path(__file__).with_suffix('.hpp')]},
                command=sys.orig_argv,mantissa_bits=int(np.finfo(LD).nmant),dimensions=dims,active=int(active.sum()),
                depth=float(length),exit_radius=float(trial.shape.exit),force_sign='F = -dE/d(log depth)',
                shifts=list(map(float,states)),baseline_crystal_energy=float(energy[0].sum(dtype=LD)),
                analytic_energy_terms=['moving physical frame','transport of material labels','motion of radial integration cuts'],
                analytic_minus_native_energy_derivative=float(analytic.sum(dtype=LD)+archived['elastic']),
                levels=reduced,top=top,cumulative=cumulative,
                signed_sum_cell_force_error=float(force_errors.sum()),archived_total_force_error=ftotal,
                noncell_remainder=float(ftotal-force_errors.sum()))
    tmp=out/'cells.tmp'
    tmp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    tmp.replace(out/'cells.json')
    print(json.dumps(dict(case=case,source=SOURCE[case],active=int(active.sum()),levels=reduced,top=top[:1],
                          cumulative=cumulative),allow_nan=False),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('scratch',type=Path)
    p.add_argument('case',choices=SOURCE)
    p.add_argument('output',type=Path)
    a=p.parse_args()
    run(a.scratch.resolve(),a.case,a.output.resolve())
