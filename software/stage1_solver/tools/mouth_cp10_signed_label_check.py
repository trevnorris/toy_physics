"""Compare signed depth-transport tangents with direct inversions of the saved crystal."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_cp10_label_grid import allnode_inverse, prisms
from mouth_cp10_nested_forces import subdivide
from mouth_cp10_spatial_crystal import load_context, signed_geometry
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest


def check(root, output):
    native = root/'r4_cp10_labels01/labels.npz'
    ctx = load_context(root,native)
    body = ctx[2]
    old_t = body.buffers['field_t']
    with np.load(root/'r4_melted_survival01/coarse/state_00010.npz') as saved:
        mesh = prisms(body,saved['solid'])
    result = dict(checkpoint_sha256=CHECKPOINT_SHA,grids={})
    for level in (1,2,4):
        path = native if level==1 else root/f'r4_nested_labels01/{level}/labels.npz'
        with np.load(path) as data:
            axes = [data[k] for k in 'xta']
            mask, xi = data['active_stencil'], data['labels']
            if level>1:
                minus,plus = (data['dxi_dt_'+s] for s in ('minus','plus'))
            else:
                with np.load(native.parent/'label_transport.npz') as rates:
                    if str(rates['label_sha256'])!=digest(native):
                        raise ValueError('Tangent input labels changed')
                    tangents = rates['dxi_dt']
                    minus,plus = tangents[...,0,:],tangents[...,1,:]
        row = dict(labels_sha256=digest(path),nodes=int(mask.sum()),probes=[])
        for h in (-4e-7,4e-7,-2e-7,2e-7):
            t = signed_geometry(ctx,0,h).buffers['field_t']
            signed_axes = [axes[0],subdivide(t,level),axes[2]]
            points = np.stack(np.meshgrid(*signed_axes,indexing='ij'),axis=-1)
            exact, n, duplicates, landing = allnode_inverse(signed_axes,points,mesh,mask)
            if duplicates or np.any(n[mask]!=1):
                raise ArithmeticError('Unresolved signed inverse on the closure of ordered cells')
            dt = subdivide(t-old_t,level)[None,:,None,None]
            rate = np.where(dt>=0,plus,minus)
            transported = xi+np.multiply(dt,rate,out=np.zeros_like(xi),where=dt!=0.)
            delta = np.max(abs(transported[mask]-exact[mask]),axis=-1)
            idx = np.argwhere(mask)[np.argmax(delta)]
            row['probes'].append(dict(h=h,max_abs_label_difference=float(delta.max()),
                                      max_difference_index=idx.tolist(),max_roundtrip=landing))
        result['grids'][str(level)] = row
    with output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
    print(json.dumps(result,allow_nan=False),flush=True)


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('root',type=Path)
    cli.add_argument('output',type=Path)
    args=cli.parse_args()
    check(args.root.resolve(),args.output.resolve())
