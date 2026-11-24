"""Locate the synthetic weld in the held mouth and check normalized and elastic volume ratios."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import binary_dilation
from scipy.optimize import brentq

from mouth_cp10_label_grid import get_cp10
from mouth_cp10_label_orientation import native_state
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest
from mouth_supercore_shape import preparation


AMPLITUDE=.12


def ranges(det, selected):
    return [float(det[selected].min()),float(det[selected].max())]


def audit(root, destination):
    destination.mkdir(parents=True,exist_ok=True)
    body,state,axes=get_cp10(root)
    amplitude=AMPLITUDE
    cut=body.cut
    x,t,a=axes
    old_path=root/'r4_bucket_geometry01/masks.npz'
    with np.load(old_path,allow_pickle=False) as old:
        if not all(np.array_equal(old[k],axis) for k,axis in zip('xta',axes)):
            raise ValueError('The CP10 comparison buckets belong to a different chart')
        seam=old['seam'].astype(bool)
        old_c=(old['bucket']==2)
        if not np.array_equal(old_c,seam):
            raise ValueError('CP10 legacy fold neighborhood extends outside the physical seam mask')
        old_weld=old['weld'].astype(bool)
        old_bucket=old['bucket'].copy()
    s=np.sin(np.pi*t)
    low=np.minimum(s[:-1],s[1:])
    high=np.maximum(s[:-1],s[1:])
    high[(t[:-1]<=.5) & (.5<=t[1:])]=1.
    # xi_r strictly increases in x; its minimum/maximum at fixed x occur
    # at the extremes of sin(pi*t). Its closed range therefore determines
    # whether the analytic synthetic weld intersects the entire native cell.
    xlo=x[:-1,None]; xhi=x[1:,None]
    lo=xlo+amplitude*xlo*(1-xlo*xlo)*low
    hi=xhi+amplitude*xhi*(1-xhi*xhi)*high
    meet=(lo<=cut)&(cut<=hi)
    if not np.array_equal(meet,(lo<=cut+1e-13)&(cut-1e-13<=hi)):
        raise ArithmeticError('Synthetic weld mask is sensitive to endpoint roundoff')
    hit=np.broadcast_to(meet[:,:,None],(len(x)-1,len(t)-1,len(a)-1)).copy()
    weld=binary_dilation(hit,structure=np.ones((3,3,3),bool))
    bucket=np.where(seam,2,np.where(weld,1,0)).astype(np.uint8)
    order=state['field'][2,body.ids].transpose(0,2,1)
    active=np.lib.stride_tricks.sliding_window_view(order,(2,2,2)).max(axis=(3,4,5))>0
    # Independent 1-D weld radius by root inversion at both endpoints and
    # the possible interior extremum in t. There is no angular dependence.
    for j in range(len(t)-1):
        roots=[brentq(lambda u:u+amplitude*u*(1-u*u)*ss-cut,0.,1.,xtol=5e-16)
               for ss in (low[j],high[j])]
        check=(x[1:]>=min(roots)-2e-14)&(x[:-1]<=max(roots)+2e-14)
        if not np.array_equal(check,meet[:,j]):
            raise ArithmeticError('Monotone inverse and range-intersection masks differ')
    cache=destination/'masks.npz'
    if cache.exists():
        raise FileExistsError(cache)
    np.savez_compressed(cache,x=x,t=t,a=a,bucket=bucket,hit=hit,weld=weld,seam=seam,
                        legacy_bucket=old_bucket,legacy_weld=old_weld)
    # The chart's inverse determinant is xi_x. In a Cartesian 3-brane
    # radial reconstruction the two angular factors contribute (xi_r/x)^2.
    determinant_bounds=[1-2*amplitude,1+amplitude]
    cartesian_bounds=[determinant_bounds[0],(1+amplitude)**3]
    extrema={k:[float('inf'),-float('inf')] for k in ('ordered','melted','core','collar')}
    witness={}
    ref=preparation(4)[0]
    grid_shape=(len(x)-1,len(t)-1,len(a)-1)
    for frac in itertools.product((.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)),repeat=3):
        X,T,A=np.meshgrid(*[axis[:-1]+p*np.diff(axis) for p,axis in zip(frac,axes)],indexing='ij')
        b=amplitude*np.sin(np.pi*T)
        label=X*(1+b*(1-X**2))
        g=1+b*(1-3*X**2)
        shape=body.sh.radius(T*body.sh.length)[0]
        y=(T-.06)/.56
        band=np.where((y>0)&(y<1),16*y*y*(1-y)**2*body.epsilon,0.)
        p2=(3*np.cos(A)**2-1)/2
        R=X*shape+X**3*band*p2
        Rx=shape+3*X**2*band*p2
        reference_radius=ref.radius(T*ref.length)[0]
        S=label*reference_radius*body.data.core_ratio
        area=reference_radius*body.data.core_ratio*ref.length
        collar=label>cut
        er,tr=body.buffers['reg_e'],body.buffers['reg_t']
        l=label[collar]
        e=-np.log(l)/body.data.pressure
        i=(np.searchsorted(er,e,side='right')-1).clip(0,len(er)-2)
        j=(np.searchsorted(tr,T[collar],side='right')-1).clip(0,len(tr)-2)
        u=(e-er[i])/np.diff(er)[i];v=(T[collar]-tr[j])/np.diff(tr)[j]
        upper=(u+v>1).astype(int)
        q=body.buffers['labels'].reshape(len(tr),len(er),2)
        c=q[j+upper,i+upper]
        du=q[j+upper,i+1-upper]-c;dv=q[j+1-upper,i+upper]-c
        S[collar]=4+c[:,0]+np.where(upper,1-u,u)*du[:,0]+np.where(upper,1-v,v)*dv[:,0]
        area[collar]=abs((du[:,0]*dv[:,1]-du[:,1]*dv[:,0])/(np.diff(er)[i]*np.diff(tr)[j]*l*body.data.pressure))
        forward=body.sh.length*Rx*R*R/(area*S*S*g)
        if not np.isfinite(forward).all() or (forward<=0).any():
            raise ArithmeticError('Nonpositive sampled physical deformation')
        for key,mask in (('ordered',active),('melted',~active),('collar',collar),('core',~collar)):
            lo,hi=ranges(forward,mask)
            first=lo<extrema[key][0]
            extrema[key]=[min(extrema[key][0],lo),max(extrema[key][1],hi)]
            if first:
                idx=np.unravel_index(np.where(mask,forward,np.inf).argmin(),grid_shape)
                xx,tt,aa=(float(ar[idx]) for ar in (X,T,A))
                G=np.eye(3);G[0]=[g[idx],amplitude*np.pi*xx*(1-xx**2)*np.cos(np.pi*tt),0.]
                reference=native_state(body,(label[idx],tt,aa),(xx,tt,aa),G)
                if not np.isclose(reference['J'],lo,atol=1e-13,rtol=1e-12):
                    raise ArithmeticError('Independent pullback with both hoop factors differs')
                witness[key]=dict(chart=[xx,tt,aa],J=lo,side='collar' if collar[idx] else 'core',
                                  independent_SVD=reference)
    report=dict(checkpoint_sha256=CHECKPOINT_SHA,amplitude=amplitude,normalized_chart_cut=cut,
                legacy_masks_sha256=digest(old_path),masks_sha256=digest(cache),changed_native_cells=int((old_bucket!=bucket).sum()),
                changed_active_native_cells=int(((old_bucket!=bucket)&active).sum()),
                weld_changed_including_C=int(np.sum(weld!=old_weld)),
                normalized_inverse_chart_bounds=determinant_bounds,normalized_cartesian_inverse_bounds=cartesian_bounds,
                normalized_cartesian_forward_bounds=[1/cartesian_bounds[1],1/cartesian_bounds[0]],
                normalized_cartesian_scope='Radial 3D vector scaling plus normalized depth; independent of the two stored stress-free reference metrics.',
                synthetic_region_scope='Exact image of the synthetic weld in held native chart, one cell halo; C uses CP10 order-front neighborhood and gets overlap priority.',
                active_cell_counts=[int(np.sum(active&(bucket==k))) for k in range(3)],
                physical_forward_J_sample=extrema,sample_points=8*int(np.prod(grid_shape)),minimum_witnesses=witness)
    with (destination/'result.json').open('x') as f:
        json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(report,allow_nan=False),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path);p.add_argument('destination',type=Path)
    a=p.parse_args();audit(a.root.resolve(),a.destination.resolve())
