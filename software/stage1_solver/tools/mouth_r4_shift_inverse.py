"""Direct inverse labels throughout an unfolded R4 crystal, including its melted region."""
import argparse
from pathlib import Path

import numpy as np

from mouth_cp10_label_grid import allnode_inverse, prisms
from mouth_cp10_nested_forces import subdivide
from mouth_cp10_spatial_crystal import signed_geometry
from mouth_cp10_unfolded_depth import control_context, H, read, write_new
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest


def shifted(t):
    """One quarter of the downstream interval at every interior knot."""
    out=t.copy()
    out[1:-1]+=0.25*np.diff(t)[1:]
    if not np.all(np.diff(out)>0) or not np.array_equal(out[[0,-1]],t[[0,-1]]):
        raise ArithmeticError('Invalid interior shift')
    return out


def prepare(root,area,level):
    old=root/'r4_bucket_unfolded01'
    ctx,solid=control_context(root,old)
    dest=area/f'labels_{level}'
    dest.mkdir(exist_ok=False,parents=True)
    old_archive=read(old/'depth/prepared.json')
    old_static=old/str(level)/'labels.npz'
    old_signed=old/'depth'/f'signed_{level}.npz'
    auth=old_archive['grids'][str(level)]
    if auth['static_labels_sha256']!=digest(old_static) or auth['signed_sha256']!=digest(old_signed):
        raise ValueError('Old signed bank provenance')
    mesh=prisms(ctx[2],solid)
    report=dict(level=level,parent_sha256=CHECKPOINT_SHA,control_sha256=digest(old/'control.npz'),
                h=H,static_reference_sha256=digest(old_static),signed_reference_sha256=digest(old_signed),cases={})
    with np.load(old_static) as base,np.load(old_signed) as archive:
        x,a=base['x'],base['a']
        order=base['order']; chi_mask=base['active_stencil']
        shape=order.shape
        for kind in ('unshifted','shifted'):
            content=dict(x=x,a=a,order=order)
            tests={}
            for sign,tag in ((0,'zero'),(-1,'minus'),(1,'plus')):
                body=signed_geometry(ctx,0,sign*H)
                knots=subdivide(body.buffers['field_t'],level)
                t=shifted(knots) if kind=='shifted' else knots
                axes=(x,t,a)
                pos=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1)
                labels,count,duplicates,landing=allnode_inverse(axes,pos,mesh,np.ones(shape,bool))
                if not (count==1).all() or duplicates or not np.isfinite(labels).all():
                    raise ArithmeticError(f'Direct inverse failed at {level} {kind} {tag}')
                # At the registered component traces the coordinates are prescribed exactly.
                fixed=np.zeros_like(labels,dtype=bool)
                fixed[:,[0,-1]]=True; fixed[-1]=True; fixed[0,:,:,0]=True; fixed[:,:,[0,-1],2]=True
                trace=np.max(abs(labels[fixed]-pos[fixed]))
                if trace>2e-12:raise ArithmeticError(f'Boundary inverse failed {trace}')
                ref=base['labels'] if sign==0 else archive[tag+'_labels']
                # Signed banks in 344030b4 used raw direct roots, with no trace snapping;
                # reproduce that convention to compare the actual energy function.
                if sign==0:labels[fixed]=pos[fixed]
                info=dict(maximum_landing_error=landing,registration_error=float(trace),n_unique=int(count.size),
                          max_depth_shift=float((t-knots).max()),closed_faces_displacement=0.)
                if kind=='unshifted':
                    if not np.array_equal(knots,base['t'] if sign==0 else archive[tag+'_knots']):
                        raise ValueError('Changed order grid')
                    info['active_max_abs_difference']=float(abs(labels[chi_mask]-ref[chi_mask]).max())
                    info['active_equal']=bool(np.array_equal(labels[chi_mask],ref[chi_mask]))
                    info['inactive_max_abs_difference']=float(abs(labels[~chi_mask]-ref[~chi_mask]).max())
                    if info['active_max_abs_difference']>1e-13:
                        raise ArithmeticError('Old active spline and direct inverse disagree')
                content[tag+'_labels']=labels
                content[tag+'_order_t']=knots
                content[tag+'_label_t']=t
                tests[tag]=info
            path=dest/f'{kind}.npz'
            np.savez_compressed(path,**content)
            report['cases'][kind]=dict(tests=tests,sha256=digest(path))
    write_new(dest/'result.json',report)
    print(report,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path);p.add_argument('output',type=Path);p.add_argument('level',type=int,choices=(1,2,4))
    a=p.parse_args();prepare(a.root.resolve(),a.output.resolve(),a.level)
