"""Read the two disk traces and polar side work after returning relaxed strain to the front."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_joined_angular_work import deformed_work
from mouth_joined_front_audit import Front
from mouth_joined_front_trace import face
from mouth_supercore_front import check_guard


def run(a):
    check_guard(a.state);check_guard(a.fields)
    rows=list(map(json.loads,(a.fields.parent/'stdout').read_text().splitlines()))
    p=next(r for r in rows if r['kind']=='preparation')
    digest=hashlib.sha256(a.fields.read_bytes()).hexdigest()
    if (hashlib.sha256(a.state.read_bytes()).hexdigest()!=p['source_sha256'] or
            rows[-1]['field_sha256']!=digest):
        raise ValueError('Archive does not match the solved field')
    base=Front(a.state,tuple(p['grid']),p['winding'])
    model=base.joined(base.sh,False)
    # Both observables use the saved responding fields. Their geometry and
    # imposed order reference are identical before and after the elastic solve.
    # The polar probe is wholly outside the finite band where delta-chi moves.
    with np.load(a.fields) as f:
        u,rest=f['fields'],f['rest']
    print(json.dumps(dict(kind='source',grid=p['grid'],field_sha256=digest,
                          state_sha256=p['source_sha256'],winding=p['winding'])),flush=True)
    for side in (-1,1):
        vals,n,ma,j,T=face(model,base.c,u,side,16)
        T-=face(model,base.c,rest,side,16)[-1]
        choose=(n>.5)&(ma<1.);pos=np.maximum(j,0.)
        print(json.dumps(dict(kind='disk',side=side,signed=vals@j,dense_share=vals[choose]@pos[choose]/(vals@pos),
                             drive_stress=vals@T,dense_drive_stress=vals[choose]@T[choose]),allow_nan=False),flush=True)
    for step in (.0004,.0002):
        force=(deformed_work(model,u,-step)-deformed_work(model,u,step))/(2*step)
        print(json.dumps(dict(kind='polar_residual',step=step,depth_span=[.06,.62],force=force.sum()),allow_nan=False),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path);p.add_argument('fields',type=Path)
    run(p.parse_args())
