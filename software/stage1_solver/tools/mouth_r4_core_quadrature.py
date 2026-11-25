"""Isolate the unchanged core's numerical energy shift when collar cuts are removed."""
import argparse
import ctypes as ct
import json
from pathlib import Path
import subprocess
import numpy as np
from mouth_cp10_label_grid import get_cp10
from mouth_cp10_label_condition import reference_frame
from mouth_exact_checkpoint_control import digest
from mouth_r4_full_label_precision import HPData,HPGrid,P
from mouth_r4_smooth_collar_reference import enable_smooth_collar
from mouth_r4_smooth_label_derivatives import CrystalBanks
from mouth_supercore_shape import preparation


def run(study,out):
    out.mkdir(exist_ok=False,parents=True)
    step=study.parent/'r4_segmented_held01/primary/step_00166.npz'
    rec=json.loads(step.with_suffix('.json').read_bytes())
    if rec['outer']!=166 or rec['sha256']!=digest(step):
        raise ValueError('Invalid step-166 anchor')
    b,cp,_=get_cp10(study.parent)
    with np.load(step,allow_pickle=False) as st:
        if not np.array_equal(st['x'],cp['x']):
            raise ValueError('Held shape mismatch')
        labels=np.ascontiguousarray(st['labels'],np.longdouble)
        order=np.ascontiguousarray(st['field'][2,b.ids].transpose(0,2,1),np.longdouble)
    cp.close()
    enable_smooth_collar(b)
    bank=CrystalBanks(b,order,out)
    code=out/'core.cpp'
    code.write_text('#include "audit.cpp"\n#include "'+str(Path(__file__).with_suffix('.hpp').resolve())+'"\n')
    so=out/'core.so'
    subprocess.run(['g++','-std=gnu++17','-O3','-shared','-fPIC',str(code),'-o',str(so)],check=True)
    f=ct.CDLL(str(so.resolve())).core_cutsum
    f.argtypes=[ct.POINTER(HPData),P,P,HPGrid,P,P];f.restype=ct.c_int
    values=np.zeros((2,7),np.longdouble)
    pt=lambda v:v.ctypes.data_as(P)
    if f(ct.byref(bank.data),pt(labels),pt(order),bank.rule,pt(bank.weights),pt(values))!=0:
        raise ArithmeticError('Core/collar cut attribution failed')
    axes=('<rays>','legacy_segments','smooth_segments','legacy_core_energy','smooth_core_energy',
          'smooth_minus_legacy','summed_absolute_ray_shift')
    probes=np.array([[1.,.184,.86],[np.nextafter(1.,np.inf),.184,.86]])
    frame,hoop=reference_frame(b,preparation(4)[0],probes)
    discrepancy=float(np.max(abs(frame[0]-frame[1])))
    result=dict(step166_sha256=rec['sha256'],scope='only positive-order physical crystal work inside the unchanged core',
                radial_endpoint_nextafter=dict(chart=probes[:,0].tolist(),frame_max_abs=discrepancy,
                                               hoop_max_abs=float(abs(hoop[0]-hoop[1]))),
                groups={g:dict(zip(axes,v.astype(float))) for g,v in zip(('core_only_rays','weld_crossing_ray_core_portion','both'),[*values,values.sum(axis=0)])})
    (out/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,allow_nan=False))
    if discrepancy>1.e-12 or np.min(values[:,3])<=0:raise ArithmeticError('Inconsistent edge or empty ray group')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('study',type=Path);p.add_argument('out',type=Path)
    v=p.parse_args();run(v.study.resolve(),v.out.resolve())
