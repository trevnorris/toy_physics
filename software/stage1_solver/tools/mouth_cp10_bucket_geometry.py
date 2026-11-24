"""Classify native physical cells by the actual saved weld, folds and inverse/gauge seam."""
import argparse
import ctypes as ct
from itertools import product
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.ndimage import binary_dilation

from mouth_cp10_label_grid import get_cp10, prisms
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest


def folded_prisms(corners):
    """The prism's determinant is affine in barycentric coordinates and quadratic in angle."""
    def at(alpha):
        verts=corners[:,0]*(1-alpha)+corners[:,1]*alpha
        normal=np.cross(verts[:,1]-verts[:,0],verts[:,2]-verts[:,0])
        return np.einsum('ei,evi->ev',normal,corners[:,1]-corners[:,0])
    a,m,b=at(0),at(.5),at(1)
    quadratic=2*((b-m)-(m-a))
    linear=b-a-quadratic
    turn=np.divide(-linear,2*quadratic,out=np.zeros_like(a),where=quadratic>0)
    valley=(quadratic>0)&(0<turn)&(turn<1)
    minimum=np.minimum(a,b)
    minimum[valley]=np.minimum(minimum[valley],(a+turn*(linear+turn*quadratic))[valley])
    return minimum.min(axis=1)<=0.


def raster(axes,patches,library,pad):
    """AABB enumerates candidates. Swept triangle/box feasibility certifies contact."""
    dims=tuple(len(a)-1 for a in axes)
    hits=np.zeros(dims,bool)
    candidates=0
    for patch in patches:
        patch=np.ascontiguousarray(patch,dtype=float)
        low=patch.min(axis=(0,1))
        high=patch.max(axis=(0,1))
        ranges=[]
        for grid,l,h in zip(axes,low,high):
            gap=pad*np.max(np.diff(grid))
            ranges.append(range(max(0,np.searchsorted(grid,l-gap,side='left')-1),
                                min(len(grid)-1,np.searchsorted(grid,h+gap,side='right'))))
        for idx in product(*ranges):
            if hits[idx]:
                continue
            candidates+=1
            l=np.array([a[i] for a,i in zip(axes,idx)])
            h=np.array([a[i+1] for a,i in zip(axes,idx)])
            result=library(patch,l,h,pad)
            if result<0:raise ArithmeticError('Degenerate Eulerian cell')
            hits[idx]=bool(result)
    return hits,candidates


def geometry(root,out):
    out.mkdir(parents=True,exist_ok=True)
    if (out/'masks.npz').exists():
        raise FileExistsError(out/'masks.npz')
    src=Path(__file__).with_name('mouth_cp10_bucket_intersection.cpp')
    lib=out/'intersections.so'
    subprocess.run(['g++','-O3','-std=c++17','-fPIC','-shared',str(src),'-o',str(lib)],check=True)
    intersects=ct.CDLL(str(lib)).intersects
    arr=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    intersects.argtypes=[arr,arr,arr,ct.c_double];intersects.restype=ct.c_int
    body,state,axes=get_cp10(root)
    labels=root/'r4_cp10_labels01/labels.npz'
    earlier=Path(__file__).resolve().parents[3]/'research/pde_ledger_v3/_measurements/mouth_r4_checkpoint10_boundary.json'
    old=json.loads(earlier.read_bytes())
    if old['checkpoint_sha256']!=CHECKPOINT_SHA:
        raise ValueError('Wrong saved fold control')
    with np.load(labels,allow_pickle=False) as data:
        if str(data['source_sha256'])!=CHECKPOINT_SHA or not all(np.array_equal(data[n],a) for n,a in zip('xta',axes)):
            raise ValueError('Inconsistent native chart')
        stencil=data['active_stencil']
        order=data['order']
        melt=data['melted_cells']
    chi_vertices=np.lib.stride_tricks.sliding_window_view(order,(2,2,2))
    if not np.array_equal(melt,chi_vertices.max(axis=(-1,-2,-3))==0):
        raise ValueError('Saved order extent disagrees with supplied harmonic/spatial partition')
    all_native=np.lib.stride_tricks.sliding_window_view(stencil,(2,2,2)).all(axis=(-1,-2,-3))
    active_with_harmonic_only=int(np.sum(~melt&~all_native))
    if active_with_harmonic_only:
        raise ValueError('Positive-order cell uses harmonic-only labels')
    keys,_,landed,_,_=prisms(body,state['solid'])
    folded=folded_prisms(landed)
    folds=[keys[i] for i in np.where(folded)[0]]
    expected=[tuple(f['cell_t_x_angle_triangle']) for f in old['folded_prisms']]
    if sorted(folds)!=sorted(expected) or len(folds)!=8:
        raise ValueError('Saved folds differ from the determinant minimum at checkpoint 10')
    T,X,A=np.meshgrid(body.t,body.x,body.a,indexing='ij')
    coords=np.stack((X,T,A),axis=-1)+state['solid']
    ring=coords[:,body.cut_node]
    weld_patches=[]
    for j,k in np.ndindex(len(body.t)-1,len(body.a)-1):
        ends=ring[j:j+2,k:k+2].transpose(1,0,2)
        weld_patches.append(ends[:,[0,1,1]])
    # The seam is the topological common closure of cut cells with positive
    # Q1 order in their interior and closed, completely melted cells.
    neighbor=np.ones((3,3,3),bool)
    dilate=lambda m: binary_dilation(m,structure=neighbor)
    seam_hit=(melt&dilate(~melt))|(~melt&dilate(melt))
    seam=dilate(seam_hit)
    data=dict(seam_hit=seam_hit,seam=seam,active=~melt)
    report=dict(checkpoint_sha256=CHECKPOINT_SHA,labels_sha256=digest(labels),fold_reference_sha256=digest(earlier),
                cpp_sha256=digest(src),native_cell_shape=list(melt.shape),folded_material_prisms=folds,
                positive_order_cells_touching_harmonic_only_nodes=active_with_harmonic_only,
                active_rule='Q1 cell with any chi > 0 vertex; its entire nodal closure is material-inverted.',
                seam_rule='Closed active and fully melted native Q1 cells meeting across a face, edge or vertex.',
                adjacency='Closed cell hit plus one full native Chebyshev layer (sharing a face, edge or vertex).',
                precedence='Primary C=(fold OR seam) before B=weld. Individual masks retain all overlaps.',
                fold_scope='Cells intersecting the whole mapped image of a prism with nonpositive determinant somewhere.',
                weld_scope='Image of the original shared core/collar material surface under the saved physical placement.',
                cut_node=int(body.cut_node),material_weld_chart_radius=float(body.x[body.cut_node]),checks={})
    for name,patches in [('weld',weld_patches),('fold',landed[folded])]:
        mask,n=raster(axes,patches,intersects,2e-12)
        # Alternate exact subdivision does not change the underlying prism.
        # A second pad in units of local cell size checks roundoff at contacts.
        split=[half for patch in patches for half in
               (np.stack((patch[0],(patch[0]+patch[1])*.5)),
                np.stack(((patch[0]+patch[1])*.5,patch[1])))]
        other,_=raster(axes,split,intersects,4e-12)
        difference=int(np.sum(mask!=other))
        if difference:raise ArithmeticError(f'{name} swept-cell map changed under exact angular subdivision: {difference}')
        report['checks'][name]=dict(unique_direct_hits=int(mask.sum()),tested_candidates=n,
                                    subdivision_and_pad_changed_cells=difference)
        data[name+'_hit']=mask
        data[name]=dilate(mask)
    bucket=np.where(data['fold']|data['seam'],2,np.where(data['weld'],1,0)).astype(np.uint8)
    bits=data['weld'].astype(np.uint8)+2*data['fold']+4*data['seam']
    bits=bits.astype(np.uint8)
    np.savez_compressed(out/'masks.npz',**{k:v.astype(np.uint8) for k,v in data.items()},
                        bits=bits,bucket=bucket,x=axes[0],t=axes[1],a=axes[2])
    report.update(cell_counts={k:int(v.sum()) for k,v in data.items()},
                  active_counts={k:int(np.sum(v&~melt)) for k,v in data.items()},
                  bit_code='weld + 2*fold + 4*seam, all with one neighbor layer',
                  independent_counts_all=np.bincount(bits.ravel(),minlength=8).tolist(),
                  independent_counts_active=np.bincount(bits[~melt],minlength=8).tolist(),
                  masks_sha256=digest(out/'masks.npz'))
    (out/'geometry.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,allow_nan=False),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    geometry(args.root.resolve(),args.output.resolve())
