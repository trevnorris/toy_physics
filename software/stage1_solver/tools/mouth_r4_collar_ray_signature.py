"""Cross-check the line's physical triangle banks and first endpoint root from binary inputs."""
import argparse
import ctypes as ct
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from mouth_r4_collar_crossings import digest
from mouth_r4_aligned_derivatives import load
from mouth_r4_eulerian_release import order_grid
from mouth_r4_full_label_precision import HPData, HPGrid, P
from mouth_r4_smooth_label_derivatives import CrystalBanks

def run(area):
    original=area/'crossings/result.json'
    extent=area/'crossings_extension/result.json'
    r=json.loads(original.read_bytes()); extra=json.loads(extent.read_bytes())
    if extra['prior_sha256']!=digest(original):raise ValueError('Extension source mismatch')
    out=area/'ray_signatures'; out.mkdir(exist_ok=False)
    d,t,_,_=load(area.parent/'r4_rematch39_01','primary',out)
    st=area.parent/'r4_segmented_held01/primary/step_00166.npz'
    dp=area/'directions/directions.npz'
    if digest(st)!=r['state_sha256'] or digest(dp)!=r['directions_sha256']:raise ValueError('Unverified input')
    with np.load(st,allow_pickle=False) as saved,np.load(dp,allow_pickle=False) as bank:
        labels=np.asarray(saved['labels'],np.longdouble)
        direction=np.asarray(bank['pj'],np.longdouble)
        chi=order_grid(t,saved['field'])
        if not np.array_equal(saved['x'],d.x) or not np.array_equal(labels,bank['labels']):
            raise ValueError('Frozen shape or state differs')
    enabled=np.ascontiguousarray(np.max(np.lib.stride_tricks.sliding_window_view(chi,(2,2,2)),axis=(-3,-2,-1))>0,np.int32)
    banks=CrystalBanks(t.b,chi,out)
    src=out/'ray_signature.cpp'
    src.write_text('#include "audit.cpp"\n#include "mouth_r4_collar_ray_signature.hpp"\n')
    so=out/'ray_signature.so'
    subprocess.run(['g++','-std=gnu++17','-O3','-fPIC','-shared','-I',str(Path(__file__).resolve().parent),str(src),'-o',str(so)],check=True)
    fn=ct.CDLL(str(so)).direction_changed_rays
    fn.argtypes=[ct.POINTER(HPData),P,P,HPGrid,ct.c_void_p,ct.c_void_p]
    fn.restype=ct.c_int
    e=r['first']; ir,j,k=e['cell']; u,v,a=e['local']; pn=F(float(e['plane']))
    def affine(X,offset=F(0)):
        return sum((F(float(X[ir+int(u),j+b,k+c,1]))-offset)*(F(v) if b else 1-F(v))
                   *(F(a) if c else 1-F(a)) for b in (0,1) for c in (0,1))
    root=-affine(labels,pn)/affine(direction)
    if float(root)!=r['alpha_cross']:raise ValueError('Reported crossing differs from exact rational binary data')
    target=r['alpha_cross']
    sig=[]
    for alpha in sorted(set([w['alpha'] for w in r['ladder']+extra['ladder']]+[.999*target,1.001*target])):
        flags=np.zeros(enabled.shape+(banks.rule.n,banks.rule.n),np.int32)
        shift=np.ascontiguousarray(direction*np.longdouble(alpha))
        count=fn(ct.byref(banks.data),labels.ctypes.data_as(P),shift.ctypes.data_as(P),banks.rule,
                 enabled.ctypes.data,flags.ctypes.data)
        ids=np.argwhere(flags)
        if count!=len(ids):raise ArithmeticError('Invalid ray count')
        sig.append(dict(alpha=alpha,multiple=alpha/target,changed_cell_rays=count,all_cell_rays=ids.tolist()))
    return dict(command=sys.orig_argv,source_files={p.name:digest(p) for p in
                (Path(__file__),Path(__file__).with_suffix('.hpp'),Path(__file__).with_name('mouth_r4_segmented_endpoint_forces.hpp'))},
                exact_first_root_from_binary_inputs=dict(numerator=str(root.numerator),denominator=str(root.denominator),
                    depth_offset=float(affine(labels,pn)),depth_velocity=float(affine(direction)),root=float(root)),
                events_primary=digest(original),events_extension=digest(extent),active_cells=int(enabled.sum()),
                physical_triangle_signatures=sig,
                interpretation='At each tested alpha, compare actual sequences of collar reference triangles against zero, on every active production ray using extended-precision production cuts. One shared physical endpoint can change two incident cell rays. Counts can include changes of other radial or diagonal physical branches.')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('area',type=Path);a=p.parse_args()
    area=a.area.resolve(); result=run(area)
    (area/'ray_signatures'/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,allow_nan=False),flush=True)
