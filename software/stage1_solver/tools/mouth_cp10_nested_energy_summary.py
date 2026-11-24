"""Reduce the saved pointwise CP10 refinement comparison across source prisms."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reduce(root, original):
    summary={}
    for factor in (2,4):
        location=root/str(factor)
        receipt=json.loads((root/('guard'+str(factor))/'outcome.json').read_bytes())
        row=json.loads((location/'result.json').read_bytes())
        if receipt['exitCode']!=0 or receipt['childOutcome']['exitCode']!=0 or not receipt['limitsVerified']:
            raise ValueError('A measured energy did not complete')
        if digest(location/'element_energies.npy')!=row['element_sha256']:
            raise ValueError('Per-element record changed')
        cells=np.load(location/'element_energies.npy')
        def energy_region(subset):
            val=subset.reshape(-1,5).sum(axis=0,dtype=np.longdouble)
            return dict(positive_points=int(val[0]),native=float(val[1]),projected=float(val[2]),
                        difference=float(val[2]-val[1]),sum_absolute_pointwise=float(val[4]))
        cut=cells.shape[1]*2//3
        bands={f'{start}/24 to {end}/24':energy_region(cells[start:end]) for start,end in ((0,2),(2,5),(5,9),(9,24))}
        sides={'core':energy_region(cells[:,:cut]),'collar':energy_region(cells[:,cut:])}
        summary[str(factor)]=dict(relative_deficit=row['signed_relative_deficit'],by_material_depth=bands,by_reference_side=sides)
    base=json.loads(original.read_bytes())
    d1=abs(base['energy_difference']/base['native'])
    d2,d4=(summary[str(i)]['relative_deficit'] for i in (2,4))
    summary.update(relative_deficit_1x=d1,successive_empirical_orders=[float(np.log2(d1/d2)),float(np.log2(d2/d4))])
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('root','original','output'):p.add_argument(key,type=Path)
    args=p.parse_args()
    r=reduce(args.root,args.original)
    with args.output.open('x') as f:json.dump(r,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(r,allow_nan=False))
