"""Probe the two rays' incident cells under the bicubic registration at step 166."""
import argparse
import json
from pathlib import Path
import numpy as np

from mouth_cp10_label_grid import get_cp10
from mouth_r4_smooth_collar_reference import enable_smooth_collar
from mouth_r4_newton_audit import PromotedGradient
from mouth_r4_smooth_label_derivatives import CrystalBanks
from mouth_exact_checkpoint_control import digest


def run(study,output):
    output.mkdir(exist_ok=False,parents=True)
    cross=json.loads((study/'crossings/result.json').read_text())
    archive=study/'directions/directions.npz'
    if digest(archive)!=cross['directions_sha256']:
        raise ValueError('Displacement source mismatch')
    body,cp,_=get_cp10(study.parent)
    enable_smooth_collar(body)
    st=study.parent/'r4_segmented_held01/primary/step_00166.npz'
    if digest(st)!=cross['state_sha256']:
        raise ValueError('Wrong step 166')
    with np.load(st,allow_pickle=False) as step:
        if not np.array_equal(step['x'],cp['x']) or int(step['radius'])!=4:
            raise ValueError('Mouth shape changed')
        y=np.ascontiguousarray(step['labels'],np.longdouble)
        order=np.ascontiguousarray(step['field'][2,body.ids].transpose(0,2,1),np.longdouble)
    cp.close()
    with np.load(archive,allow_pickle=False) as data:
        d=data['pj'].astype(np.longdouble)
        if not np.array_equal(data['labels'],y) or not np.array_equal(data['order'],order):
            raise ValueError('Invalid held labels')
    dims=tuple(n-1 for n in y.shape[:3]); count=np.prod(dims)
    hp=PromotedGradient(body,output); fd=CrystalBanks(body,order,output)
    xc=np.longdouble(str(cross['alpha_cross'])); h=np.longdouble(2)**-22
    entries=[]
    for cell in ([23,2,17],[24,2,17]):
        idx=np.ravel_multi_index(cell,dims)
        mask=np.arange(count)==idx
        def energy(z):
            return fd.energy(z,mask)[idx]
        e0=energy(y)
        for fraction in [0.,.999,1.001, 1/(8192*xc)]:
            alpha=np.longdouble(fraction)*xc
            z=y+alpha*d
            grad=hp.gradient(z,order,mask)
            slope=np.sum(grad*d,dtype=np.longdouble)
            slope_fd=(energy(z+h*d)-energy(z-h*d))/(2*h)
            entries.append(dict(cell=cell,alpha=float(alpha),energy_delta=float(energy(z)-e0),
                                physical_directional_slope=float(slope),signed_central=float(slope_fd),
                                difference=float(slope_fd-slope)))
    record=dict(step_sha=cross['state_sha256'],source_pj=cross['directions_sha256'],
                reference='clamped weld-preserving bicubic',scope='entire cell, crystal term only',
                work=entries,max_slope_error=max(abs(r['difference']) for r in entries))
    (output/'result.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    print(json.dumps(record,allow_nan=False))
    if record['max_slope_error']>1.e-10:raise ArithmeticError('Promoted directional slope does not close')


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('study',type=Path);a.add_argument('out',type=Path)
    v=a.parse_args();run(v.study.resolve(),v.out.resolve())
