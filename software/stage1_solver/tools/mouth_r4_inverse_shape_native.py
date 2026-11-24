import argparse,json
from hashlib import sha256
from pathlib import Path

import numpy as np

from mouth_free_shape_solid import Solid
from mouth_r4_refined_solid import RefinedSolid
from mouth_r4_fold_order import axes_model
from mouth_sealed_shape_relax import Meridian

def read(p):return json.loads(p.read_bytes())
def digest(p):return sha256(p.read_bytes()).hexdigest()

def native(area,folder):
    cp=read(area/'r4_melted_survival01/coarse/latest.json')
    post=read(area/'r4_local_refinement01/first_step/result.json')
    old=area/'r4_feasible_descent01/coarse/stalled_00011.npz'
    prep=area/'supercore_joined_front01/adopt309/state.npz'
    depth_ref=read(area/'r4_whole_mouth01/coarse/endpoint.json')['depth']
    rows={}
    for name,p,expected in [('checkpoint10',area/'r4_melted_survival01/coarse'/cp['file'],cp['sha256']),
                             ('postshape',area/'r4_local_refinement01/first_step/material.npz',post['state_sha256'])]:
        if digest(p)!=expected or digest(old)!=post['source_sha256']:raise RuntimeError('Wrong material snapshot')
        with np.load(p,allow_pickle=False) as state,np.load(old,allow_pickle=False) as stalled:
            if str(state['state_sha'])!=digest(prep) or int(state['radius'])!=4 or int(state['winding'])!=1:
                raise RuntimeError('Wrong preparation')
            reference=state if name=='checkpoint10' else stalled
            params,extra,epsilon=reference['shape_params'],float(reference['shape_extra']),float(state['x'][5])
            if name=='postshape' and not np.array_equal(state['x'],stalled['x']):raise RuntimeError('Wrong geometry')
            table=[]
            for n in range(1,6):
                ee=[]
                for delta in [-1e-4,1e-4,-5e-5,5e-5]:
                    x=params.copy();v=epsilon;ve=extra
                    if n<4:x[n]+=.25*delta
                    elif n==4:ve+=.25*delta
                    else:v+=delta
                    sh=Meridian(4.,x,extra=ve)
                    b=Solid(prep,sh,v,8,24,6,axes_model(sh,depth_ref),degree=3,subdivisions=2)
                    b.melted_labels=True;b.data.exact_melt_support=1;b.buffers['base'].fill(0.)
                    obj=RefinedSolid(b,folder/'native_build') if name=='postshape' else b
                    ee.append(obj.evaluate(state['q'],state['field'][2],load=False,gradient=False)[0])
                table.append(ee)
            rows[name]=dict(state_sha256=expected,signed_energy=table,
                            forces_h=[(y[0]-y[1])/2e-4 for y in table],
                            forces_half=[(y[2]-y[3])/1e-4 for y in table])
            print(json.dumps({name:rows[name]}),flush=True)
    return rows
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('area',type=Path);p.add_argument('folder',type=Path)
    a=p.parse_args();out=a.folder/'native_work.json'
    if out.exists():raise FileExistsError(out)
    out.write_text(json.dumps(native(a.area,a.folder),indent=2,allow_nan=False)+'\n')
