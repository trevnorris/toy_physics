"""Unfitted two-sided and bounded order-force checks on every free R4 order DOF."""
import argparse
from itertools import product
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np
from scipy import sparse

from mouth_r4_projected_fields import EXPECTED, MANIFEST, digest, read
from mouth_r4_smooth_order_boundary import LD, PhysicalOrderEnergy, gauge_coefficients


class GaugeOrderBanks:
    """Independent two-sided gauge work with other coordinates held fixed."""
    def __init__(self, trial, labels, fields, *, barrier=None):
        from mouth_cp10_label_condition import physical_frame
        from mouth_cp10_label_gauge import VERTICES, gradient_basis
        from mouth_r4_eulerian_release import block_axes
        axes=block_axes(trial.b)
        sites=trial.number[trial.b.ids.transpose(0,2,1)]
        cells=np.indices(np.subtract(sites.shape,1)).reshape(3,-1).T
        gidx=np.stack([sites[tuple((cells+p).T)] for p in VERTICES],axis=1)
        state=np.stack([labels[tuple((cells+p).T)] for p in VERTICES],axis=1).astype(LD)
        lo=np.column_stack([a[cells[:,i]] for i,a in enumerate(axes)])
        widths=np.column_stack([np.diff(a)[cells[:,i]] for i,a in enumerate(axes)])
        corner_chart=lo[:,None,:]+VERTICES[None,:,:]*widths[:,None,:]
        offset=state-corner_chart.astype(LD)
        if barrier is None:
            density=np.zeros(len(cells),dtype=LD)
            for k,p in enumerate(VERTICES):
                edges=[]
                for i in range(3):
                    v=p.copy()
                    v[i]=1-v[i]
                    other=int(v@np.array([4,2,1]))
                    edges.append((1-2*p[i])*(state[:,other]-state[:,k])/widths[:,i,None])
                d=np.sum(edges[0]*np.cross(edges[1],edges[2]),axis=1)
                if np.any(d<=0):
                    raise ValueError('Negative corner determinant')
                density+=((d-1)-np.log1p(d-1))/8
        else:
            density=barrier(state,widths)
        matrices=[]
        energy=[]
        points=fields[2,trial.first]
        rows=np.repeat(np.arange(len(cells)),8)
        chis=[]
        ones=np.all(points[gidx]==1,axis=1)
        for q in product((.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)),repeat=3):
            frame,hoop=physical_frame(trial.b,lo+q*widths)
            weight=(np.prod(widths,axis=1)*4*np.pi/8*hoop*np.linalg.det(frame)).astype(LD)
            g=np.einsum('nvc,vd->ncd',offset,gradient_basis(q).astype(LD))/widths[:,None,:]
            physical=g@np.linalg.inv(frame).astype(LD)
            energy.append(weight*(density+np.sum(physical**2,axis=(1,2))/2))
            phi=np.where(VERTICES,q,1-np.asarray(q)).astype(LD).prod(axis=1)
            phi/=np.sum(phi)
            mat=sparse.coo_matrix((np.tile(phi,len(cells)),(rows,gidx.ravel())),shape=(len(cells),len(points))).tocsr()
            c=mat@points.astype(LD)
            c[ones]=1 # exact partition at a constant unit trace
            chis.append(np.clip(c,0,1))
            matrices.append(mat)
        self.weight=np.concatenate(energy)
        self.c=np.concatenate(chis)
        self.N=sparse.vstack(matrices).tocsc()

    def at(self,node):
        start,end=self.N.indptr[node:node+2]
        rows=self.N.indices[start:end]
        basis=self.N.data[start:end]
        a,w=self.c[rows],self.weight[rows]
        def quotient(h):
            # Exact integral of -6t(1-t) between the clipped endpoints a,b.
            # Compute effective length before adding it to a to retain small
            # representable steps at either bound. This is a direct signed
            # energy difference of the clamped spline, without using G'(a).
            inc=np.minimum(np.maximum(h*basis,-a),1-a)
            b=a+inc
            sdiff_per_h=-(inc/h)*(2*a*(1-a)+2*b*(1-b)+a*(1-b)+b*(1-a))
            return np.sum(w*sdiff_per_h)
        return quotient


def worst(error,selected,ids):
    indexes=np.flatnonzero(selected)
    if not len(indexes):
        return None
    at=indexes[np.argmax(abs(error[indexes]))]
    return dict(signed=float(error[at]),absolute=float(abs(error[at])),group=int(ids[at]))


def run(area,output):
    src=area/'pin/software/stage1_solver/tools'
    if digest(area/'sources.json')!=MANIFEST or any(digest(src/p)!=sha for p,sha in read(area/'sources.json').items()):
        raise ValueError('Different held sources')
    sys.path.insert(0,str(src))
    from mouth_continuous_depth_mesh import axes_for_depth
    from mouth_continuous_shape_relax import ShapeTrial
    from mouth_r4_eulerian_release import Driver,order_grid
    from mouth_r4_eulerian_solid import SpatialCrystal
    from mouth_r4_stage2_outcome import verify_guard
    verify_guard(area/'primary_guard',1)
    state=area/'primary/held_stop.npz'
    if digest(state)!=EXPECTED['held_stop'] or np.finfo(LD).nmant<63:
        raise ValueError('Wrong state or unavailable extended precision')
    with np.load(state,allow_pickle=False) as st:
        labels,fields,x=(st[k] for k in ('labels','field','x'))
        if int(st['radius'])!=4 or int(st['winding'])!=1:raise ValueError('Wrong radius or winding')
    with TemporaryDirectory(dir=area.parent) as temp:
        d=Driver(area.parent,Path(temp),1e-8)
        if not np.array_equal(d.x,x):raise ValueError('Changed geometry')
        shape=d.shape(x)
        trial=ShapeTrial(d.base,d.prep,x,8,24,6,3,2,precondition=False,axes=axes_for_depth(shape.length,d.reference_depth))
        order=order_grid(trial,fields)
        _,_,source,_=SpatialCrystal(trial.b).evaluate(labels,order,False)
        trial.model.material_load=trial.ext+np.bincount(trial.b.ids.ravel(),source.transpose(0,2,1).ravel(),minlength=trial.model.size)
        _,gradient=trial.model.evaluate(fields)
        raw=trial.pull.T@gradient[2]
        analytic_unit=gauge_coefficients(trial,labels,fields)[0][0]
        physical=PhysicalOrderEnergy(trial,fields)
        gauge=GaugeOrderBanks(trial,labels,fields)
        ids=np.flatnonzero(trial.free)
        c=fields[2,trial.first][ids]
        widths=np.array([1e-8,1e-9,1e-10],LD)
        # Last axis holds negative and positive signed energy quotients.
        pfd=np.empty((len(ids),len(widths),2),LD)
        gfd=np.empty_like(pfd)
        for j,ix in enumerate(ids):
            p=physical.at(ix)
            gu=gauge.at(ix)
            for k,h in enumerate(widths):
                for m,dh in enumerate((-h,h)):
                    pfd[j,k,m]=p(dh)
                    gfd[j,k,m]=gu(dh)
        whole=np.ones_like(c,dtype=bool)
        groups={'all_free':whole,'lower_bound':c==0,'interior':(c>0)&(c<1),'upper_bound':c==1}
        rows=[]
        for k,h in enumerate(widths):
            phy_c=pfd[:,k].mean(axis=1)
            gau_c=gfd[:,k].mean(axis=1)
            row=dict(h=float(h),physical=worst(phy_c-raw[ids],whole,ids),
                     unscaled_gauge=worst(gau_c-analytic_unit[ids],whole,ids),both_gauges={})
            for eps in (1e-8,1e-10):
                analytic=raw[ids].astype(LD)+LD(eps)*analytic_unit[ids]
                central=phy_c+LD(eps)*gau_c
                entry={name:worst(central-analytic,select,ids) for name,select in groups.items()}
                for name,select,side in [('lower_bound',c==0,1),('upper_bound',c==1,0)]:
                    one=pfd[:,k,side]+LD(eps)*gfd[:,k,side]
                    entry[name+'_feasible_one_sided']=worst(one-analytic,select,ids)
                    if k:
                        coarse=pfd[:,k-1,side]+LD(eps)*gfd[:,k-1,side]
                        ratio=widths[k-1]/h
                        entry[name+'_feasible_extrapolated']=worst((ratio*one-coarse)/(ratio-1)-analytic,select,ids)
                row['both_gauges'][str(eps)]=entry
            rows.append(row)
        result=dict(state_sha256=digest(state),source_sha256=digest(Path(__file__)),
                    boundary_energy_evaluator_sha256=digest(Path(__file__).with_name('mouth_r4_smooth_order_boundary.py')),
                    pinned_sources_sha256=MANIFEST,precision_mantissa=np.finfo(LD).nmant,
                    counts={name:int(mask.sum()) for name,mask in groups.items()},widths=rows,
                    analytic='Production complete physical residual plus analytic derivative of the quadrature-weighted cubic.',
                    differences='Per-node long-double signed changes of the local discrete energies. Harmonic and determinant density recomputed for the gauge. All joined order DOFs, all six blocks, every other coordinate held fixed.',
                    boundary_extension='Physical polynomial and kinetic terms continue algebraically for central differences. The cubic gauge holds its value constant outside [0,1] at each quadrature point. One-sided differences use only the admissible nodal direction.',
                    target=1e-7)
        result['all_order_central_within_target']=all(rows[-1]['both_gauges'][str(e)]['all_free']['absolute']<=1e-7 for e in (1e-8,1e-10))
    with output.open('x') as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write('\n')
    print(json.dumps({'counts':result['counts'],'smallest_width':rows[-1]},allow_nan=False),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('area',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    run(args.area.resolve(),args.output.resolve())
