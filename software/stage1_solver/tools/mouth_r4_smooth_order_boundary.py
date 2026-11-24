"""Frozen R4 one-sided nodal order work with a quadrature-weighted cubic gauge."""
import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np
from scipy import sparse

from mouth_r4_projected_fields import EXPECTED, MANIFEST, digest, read

LD = np.longdouble
STEPS = (1e-9, 1e-10, 1e-12)
PROBES = {'previous_blocking_node': (16, 10, 1), 'gauge_witness': (24, 9, 17)}


def gauge_coefficients(trial, labels, fields):
    """Coefficients of the exact cubic work for a single free zero node.

    Preserve the two-point physical-volume harmonic integration and original
    eight one-sided corner determinant evaluations. The mean of the latter is
    constant within each cell, multiplied by s(chi_q) in the volume integral.
    """
    from mouth_cp10_label_condition import physical_frame
    from mouth_cp10_label_gauge import VERTICES, gradient_basis
    from mouth_r4_eulerian_release import block_axes
    axes = block_axes(trial.b)
    sites = trial.number[trial.b.ids.transpose(0,2,1)]
    cells = np.indices(tuple(n-1 for n in sites.shape)).reshape(3,-1).T
    nodes = np.stack([sites[tuple((cells+v).T)] for v in VERTICES], axis=1)
    points = np.stack([labels[tuple((cells+v).T)] for v in VERTICES], axis=1)
    base = np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1)
    base_corners = np.stack([base[tuple((cells+v).T)] for v in VERTICES],axis=1)
    delta = points-base_corners
    widths = np.column_stack([np.diff(ax)[cells[:,i]] for i,ax in enumerate(axes)])
    barrier = np.zeros(len(cells))
    for v in VERTICES:
        grad = np.einsum('nvc,vd->ncd',points-points[:,:1],gradient_basis(v))/widths[:,None,:]
        det = np.linalg.det(grad)
        if np.any(det <= 0):
            raise ValueError('Nonpositive determinant in saved map')
        barrier += ((det-1)-np.log1p(det-1))/8
    chis = fields[2,trial.first]
    pure = np.max(chis[nodes],axis=1)==0
    coeff = np.zeros((3,len(chis)),dtype=LD)
    total_old = LD(0)
    total_new = LD(0)
    vol = np.prod(widths,axis=1)
    rows=np.repeat(np.arange(len(cells)),8)
    for q in product((.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)),repeat=3):
        coord = base_corners[:,0]+np.asarray(q)*widths
        frame,hoop = physical_frame(trial.b,coord)
        dg = np.einsum('nvc,vd->ncd',delta,gradient_basis(q))/widths[:,None,:]
        grad_x = dg@np.linalg.inv(frame)
        weight = 4*np.pi/8*vol*hoop*np.linalg.det(frame)
        work = weight*(barrier+.5*np.sum(grad_x**2,axis=(1,2)))
        basis = np.where(VERTICES,q,1-np.asarray(q)).prod(axis=1)
        interpolator = sparse.coo_matrix((np.tile(basis,len(cells)),(rows,nodes.ravel())),
                                        shape=(len(cells),len(chis))).tocsr()
        # Coalescing repeated IDs before raising basis values to powers is
        # essential on the polar axis, which has one shared order node.
        chi = np.clip(interpolator@chis,0.,1.)
        total_old += np.sum(work[pure],dtype=LD)
        total_new += np.sum(work*(1-chi)**2*(1+2*chi),dtype=LD)
        for m,coef in enumerate((-6*chi*(1-chi),-3+6*chi,np.full_like(chi,2.)),start=1):
            coeff[m-1] += interpolator.power(m).astype(LD).T@(work.astype(LD)*coef)
    return coeff,sites,dict(old_unit=float(total_old),new_unit=float(total_new))


class PhysicalOrderEnergy:
    """Exact local energy differences for one joined Q1 order coefficient."""
    def __init__(self, trial, fields):
        from mouth_continuous_order import I,KAPPA,WELL
        m = trial.model
        self.I, self.k, self.well = map(LD,(I,KAPPA,WELL))
        self.c = m.V@fields[2]
        self.n = (m.D[:,4]*(m.V@fields[0]))**2
        self.grad_c = np.column_stack([g@fields[2] for g in m.G])
        self.grad_p = np.column_stack([g@fields[1] for g in m.G])
        self.volume = m.hq
        self.elastic = m.D[:,-1]
        self.load = trial.pull.T@m.material_load
        self.V = trial.V.tocsc()
        self.G = [g.tocsc() for g in trial.G]
        for matrix in (self.V,*self.G):
            matrix.sort_indices()

    def at(self,node):
        sl = slice(*self.V.indptr[node:node+2])
        loc = self.V.indices[sl]
        phi = self.V.data[sl].astype(LD)
        cols = []
        for G in self.G:
            gs = slice(*G.indptr[node:node+2])
            ids = G.indices[gs]
            at = np.searchsorted(loc, ids)
            if not np.array_equal(loc[at],ids):
                raise ArithmeticError('Gradient has support outside the Q1 element')
            v = np.zeros(len(loc),dtype=LD)
            v[at] = G.data[gs]
            cols.append(v)
        basisgrad = np.stack(cols,axis=1)
        gc = self.grad_c[loc].astype(LD)
        gp = self.grad_p[loc].astype(LD)
        gcb = np.sum(gc*basisgrad,axis=1)
        b2 = np.sum(basisgrad*basisgrad,axis=1)
        pgb = np.sum(gp*basisgrad,axis=1)
        dot = np.sum(gp*gc,axis=1)
        n = self.n[loc].astype(LD)
        den = n+self.I*np.sum(gc*gc,axis=1)
        w = self.volume[loc].astype(LD)
        chi = self.c[loc].astype(LD)
        potential = chi*(1-chi)
        linear_elastic = LD(self.load[node])+np.sum(w*self.elastic[loc]*phi)
        def increment(step):
            # Divide the direct energy difference by h algebraically before
            # accumulation. The h=0 limit is an independent derivative.
            h=LD(step)
            dd = self.I*(2*gcb+h*b2)
            dpot = phi*(1-2*chi-h*phi)
            flow = -.5*self.I*n*(2*dot*pgb+h*pgb**2-dot**2*dd/den)/(den+h*dd)
            local = self.k*(gcb+.5*h*b2)+self.well*dpot*(2*potential+h*dpot)+flow
            return linear_elastic+np.sum(w*local)
        return increment


def run(area, output):
    pin = area/'pin/software/stage1_solver/tools'
    if digest(area/'sources.json') != MANIFEST or any(digest(pin/name)!=v for name,v in read(area/'sources.json').items()):
        raise ValueError('Unverified held sources')
    sys.path.insert(0,str(pin))
    from mouth_continuous_depth_mesh import axes_for_depth
    from mouth_continuous_shape_relax import ShapeTrial
    from mouth_r4_eulerian_release import Driver, order_grid
    from mouth_r4_eulerian_solid import SpatialCrystal
    from mouth_r4_stage2_outcome import verify_guard
    verify_guard(area/'primary_guard',1)
    state = area/'primary/held_stop.npz'
    if digest(state)!=EXPECTED['held_stop']:
        raise ValueError('Wrong retained state')
    with np.load(state,allow_pickle=False) as st:
        fields, labels, x = (st[k] for k in ('field','labels','x'))
        if (int(st['radius']),int(st['winding']))!=(4,1):
            raise ValueError('Wrong winding or radius')
    with TemporaryDirectory(dir=area.parent) as tmp:
        d=Driver(area.parent,Path(tmp),1e-8)
        if not np.array_equal(x,d.x):
            raise ValueError('Held geometry changed')
        shape=d.shape(x)
        trial=ShapeTrial(d.base,d.prep,x,8,24,6,3,2,precondition=False,
                         axes=axes_for_depth(shape.length,d.reference_depth))
        order=order_grid(trial,fields)
        _,_,source,_=SpatialCrystal(trial.b).evaluate(labels,order,False)
        material=np.bincount(trial.b.ids.ravel(),source.transpose(0,2,1).ravel(),minlength=trial.model.size)
        trial.model.material_load=trial.ext+material
        _,g=trial.model.evaluate(fields)
        joined=trial.pull.T@g[2]
        c=fields[2,trial.first]
        if not np.array_equal(c[trial.number],fields[2]):
            raise ValueError('Split physical order')
        physical=PhysicalOrderEnergy(trial,fields)
        coeff,sites,gauges=gauge_coefficients(trial,labels,fields)
        positive_q=(physical.c>0)&(trial.model.hq>0)
        front=((trial.V.T@positive_q.astype(float))>0)&(c==0)&trial.free
        checks={name:int(sites[idx]) for name,idx in PROBES.items()}
        # Include distant pure-melt coordinates to catch dipoles of the gauge
        # with no physical first derivative, plus every front group in all baths.
        ids=np.flatnonzero((c==0)&trial.free)
        signed=np.empty((len(ids),len(STEPS)),dtype=LD)
        deriv=np.empty(len(ids),dtype=LD)
        for j,ix in enumerate(ids):
            variation=physical.at(ix)
            deriv[j]=variation(0)
            signed[j]=[LD(h)*variation(h) for h in STEPS]
        maxdg=np.max(np.abs(deriv-joined[ids]))
        if maxdg>5e-10:
            raise ArithmeticError(f'Independent physical derivative disagrees by {maxdg}')
        summary=[]
        target={}
        for j,h in enumerate(STEPS):
            H=LD(h)
            unit=H*(coeff[0]+H*(coeff[1]+H*coeff[2]))
            for eps in (1e-8,1e-10):
                gauge=LD(eps)*unit[ids]
                allwork=signed[:,j]+gauge
                def batch(select):
                    jj=np.flatnonzero(select)
                    pos=jj[signed[jj,j]>0]
                    gain=jj[(signed[jj,j]>0)&(allwork[jj]<=0)]
                    descent=jj[allwork[jj]<0]
                    at=jj[np.argmin(allwork[jj])]
                    maxratio=np.max(np.abs(gauge[pos])/signed[pos,j],initial=0)
                    return dict(count=len(jj),positive_physical=len(pos),physical_positive_total_nonpositive=len(gain),
                                total_strict_descents=len(descent),exact_zero=int(np.sum(allwork[jj]==0)),
                                sign_failure_groups=ids[gain].tolist(),min_group=int(ids[at]),
                                min_total=float(allwork[at]),physical_at_min=float(signed[at,j]),
                                gauge_at_min=float(gauge[at]),max_abs_gauge_over_positive_physical=float(maxratio))
                summary.append(dict(h=h,epsilon=eps,whole_free_zero=batch(np.ones(len(ids),bool)),front=batch(front[ids])))
                for name,ix in checks.items():
                    row=int(np.searchsorted(ids,ix))
                    if ids[row]!=ix:
                        raise ValueError('Witness left free zero set')
                    target.setdefault(name,dict(node=PROBES[name],group=ix,joins_positive_q=bool(front[ix]),
                                                slope_physical=float(deriv[row]),unit_gauge_linear=float(coeff[0,ix]),
                                                unit_gauge_quadratic=float(coeff[1,ix]),unit_gauge_cubic=float(coeff[2,ix]),samples=[]))['samples'].append(
                               dict(h=h,epsilon=eps,physical=float(signed[row,j]),gauge=float(gauge[row]),total=float(allwork[row])))
        ref=read(Path(__file__).resolve().parents[3]/'research/pde_ledger_v3/_measurements/mouth_r4_projected_fields.json')['held_stop']
        if abs(gauges['old_unit']-ref['gauge_energy']/1e-8)>5e-10:
            raise ArithmeticError('Old gauge reconstruction disagrees')
        m=front&(coeff[0]!=0)
        location=np.unravel_index(np.argmax(np.abs(np.where(m,coeff[0],0))[sites]),sites.shape)
        i=int(sites[location])
        if np.count_nonzero(m):
            max_node=dict(node=list(map(int,location)),group=i,slope=float(coeff[0,i]))
        else:
            max_node=None
        result=dict(state_sha256=EXPECTED['held_stop'],sources_sha256=MANIFEST,script_sha256=digest(Path(__file__)),
                    frozen_shape=x.tolist(),chi_gauge='1-3 chi^2+2 chi^3 on [0,1], constant outside; applied at Gauss points',
                    unit_gauge=gauges, physical_derivative_max_error=float(maxdg),
                    free_zero_count=len(ids),free_zero_front=int(front.sum()),
                    front_with_linear_gauge_term=int(np.count_nonzero(m)),greatest_absolute_front_gauge_slope=max_node,
                    quadrature='2x2x2 physical Gauss weights; harmonic density at Gauss points and original mean of eight corner barrier densities per cell.',
                    physical_work='Factored long-double direct cell energy changes with full joined Q1 support. Directed differences keep shape, other order coefficients, labels, carriers fixed.',
                    scope='Every free zero order group in all six carrier blocks including all zero groups incident to positive order at a physical quadrature point. Numerical gauge acts only on the original spatial-label grid.',
                    witnesses=target,summary=summary)
    with output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps({k:result[k] for k in ('unit_gauge','free_zero_count','free_zero_front','front_with_linear_gauge_term','summary')},allow_nan=False),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('area',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    run(args.area.resolve(),args.output.resolve())
