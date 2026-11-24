"""Frozen seven-component bath work at the two displacements used for smooth-map depth."""
import argparse
from copy import copy
import inspect
import json
from pathlib import Path

import numpy as np

from mouth_cp10_boundary_work import model_inputs
from mouth_cp10_field_shape import Bath
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read


def run(sc, dest):
    archive=Path(__file__).resolve().parents[3]/'research/pde_ledger_v3/_measurements/mouth_r4_two_metric_labels.json'
    previous=read(archive)['data']['held_fields']
    if previous['checkpoint']!=CHECKPOINT_SHA:raise ValueError('Unexpected held fields')
    build,x,_,u,_,_=model_inputs(sc)
    trial=build(x)
    owner=inspect.getclosurevars(build).nonlocals['owner']
    bath=Bath(trial,owner.base,u,owner.reference_depth)
    model=bath.model=copy(bath.model)
    cast=lambda a: np.asarray(a,dtype=np.longdouble)
    model.a=cast(model.a);model.rr=tuple(map(cast,model.rr));model.z=tuple(map(cast,model.z))
    bath.u=cast(bath.u)
    direction=cast([1,0,0,0,0,0])
    analytic=previous['modes'][0]['analytic']
    record=dict(checkpoint_sha256=CHECKPOINT_SHA,field_record_sha256=digest(archive),analytic=analytic,rows={})
    for h in (4e-7,2e-7):
        signed=np.asarray([bath.energy(s*np.longdouble(h)*direction) for s in (-1,1)])
        fd=(signed[0]-signed[1])/(2*np.longdouble(h))
        record['rows'][str(h)]=dict(components=fd.astype(float).tolist(),force=float(fd.sum()),
                                    analytic_minus_fd=float(analytic-fd.sum()))
    with dest.open('x') as f:
        json.dump(record,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(record,allow_nan=False),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('scratch',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();run(a.scratch.resolve(),a.output.resolve())
