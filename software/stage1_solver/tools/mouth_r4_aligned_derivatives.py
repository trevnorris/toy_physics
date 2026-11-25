"""Signed first derivatives on the two individually rematched step-108 R4 states.

Each invocation consumes the passing Gate 1 readbacks from both guarded rematches.
Raw positive/negative quotients are saved with every free scalar group in the bank.
"""
import argparse
import json
from itertools import product
from pathlib import Path
import sys

import numpy as np

from mouth_cp10_label_gauge import VERTICES
from mouth_r4_aligned_stencil import BASIS, EDGES, EDGE_WEIGHT, WEIGHTS
from mouth_r4_eulerian_relax_labels import chart_orientation
from mouth_r4_newton_stop import STATES
from mouth_r4_smooth_label_derivatives import CrystalBanks, LIMIT
from mouth_r4_smooth_static import Gauge, LD
from mouth_r4_stage2_outcome import digest, read, verify_guard


EPS = {'primary': 1e-8, 'secondary': 1e-10}
TAIL = 1.e-3
BROAD = (LD(2)**-24, LD(2)**-25)
NARROW = (LD(2)**-36, LD(2)**-37)
FINE = (LD(2)**-39, LD(2)**-40)
ORDER_WIDTHS = np.array([1e-8, 1e-9, 1e-10], LD)


def direct_rule(points, widths, *, minima=False):
    """Evaluate the Q1 gradient at each moving location directly.

    Each column of the gradient is affine along an edge. Interpolate the
    endpoint *columns*, then contract them at the new point. The production
    rule instead interpolates the three determinants of the edge quadratic.
    """
    points = np.asarray(points,dtype=LD)
    offsets = (points-points[:,:1]).transpose(0,2,1)
    def det(matrix):
        a,b,c = matrix.transpose(2,0,1)
        return np.sum(a*np.cross(b,c),axis=1)
    def density(j):
        if np.any(j<=0) or not np.all(np.isfinite(j)):
            raise ValueError('Infeasible signed determinant sample')
        return j-1-np.log(j)
    corners = {}
    fixed = []
    accumulated = np.zeros(len(points), LD)
    smallest = np.full(len(points), np.inf, LD)
    for q,(basis,weight) in enumerate(zip(BASIS,WEIGHTS)):
        jac = offsets@basis / widths[:,None,:]
        j = det(jac)
        fixed.append(j)
        accumulated += weight*density(j)
        np.minimum(smallest,j,out=smallest)
        if weight==1:
            corners[q] = jac
    for left,_,right in EDGES:
        a,c = np.sqrt(fixed[left]),np.sqrt(fixed[right])
        t = (a/(a+c))[:,None,None]
        j = det(corners[left]+t*(corners[right]-corners[left]))
        accumulated += EDGE_WEIGHT*density(j)
        np.minimum(smallest,j,out=smallest)
    return (accumulated/216,smallest) if minima else accumulated/216


class SignedGauge(Gauge):
    """Transport quadrature and harmonic work as in the held shape bank."""
    def cell_barrier(self, offsets, widths, gradient):
        if gradient:
            raise ValueError('The signed-energy oracle has no analytic gradient')
        return direct_rule(offsets+VERTICES*widths[:,None,:],widths),None


def load(study, case, output):
    from mouth_continuous_depth_mesh import axes_for_depth
    from mouth_continuous_shape_relax import ShapeTrial
    from mouth_exact_checkpoint_control import CHECKPOINT_SHA
    from mouth_r4_eulerian_release import order_grid
    from mouth_r4_eulerian_solid import SpatialCrystal
    from mouth_r4_smooth_fields import GaugeGeometry
    from mouth_r4_smooth_held import SmoothHeld
    if np.finfo(LD).nmant<63:
        raise ValueError('Signed banks require extended precision')
    for name,epsilon in EPS.items():
        root = study/name
        verify_guard(study/(name+'_field_guard'),0)
        verify_guard(study/(name+'_sweep_guard'),0)
        verify_guard(study/(name+'_labels_guard'),0)
        record = read(root/'rematched.json')
        gate = read(root/'gate1.json')
        label_read = read(root/'labels.json')
        if not (record['field_matched'] and record['epsilon']==epsilon and
                record['source_sha256']==STATES['last_active'][1] and
                gate['passed'] and gate['epsilon']==epsilon
                and gate['archive_sha256']==digest(root/'order_work.npz')
                and label_read['passed'] and label_read['state_sha256']==record['sha256']
                and gate['rematched_sha256']==record['sha256']==digest(root/'rematched.npz')):
            raise ValueError(f'Uncertified rematch in {root}')
    source = study/case/'rematched.npz'
    d = SmoothHeld(study.parent,output,EPS[case])
    original = study.parent/'r4_newton_held01/primary/last_active.npz'
    if digest(original)!=STATES['last_active'][1]:
        raise ValueError('Unverified step-108 labels')
    with np.load(source,allow_pickle=False) as state, np.load(original,allow_pickle=False) as start:
        if not (int(state['radius'])==4 and int(state['winding'])==1 and int(state['outer'])==108
                and int(state['stencil_points'])==39 and float(state['epsilon'])==EPS[case]
                and str(state['parent_sha256'])==CHECKPOINT_SHA and str(state['source_sha256'])==STATES['last_active'][1]
                and np.array_equal(state['x'],d.x) and np.array_equal(state['labels'],start['labels'])):
            raise ValueError('Changed held shape, ancestry, labels or preparation')
        d.labels,d.fields = state['labels'],state['field']
    t = ShapeTrial(d.base,d.prep,d.x,8,24,6,3,2,precondition=False,
                   axes=axes_for_depth(d.shape(d.x).length,d.reference_depth))
    t.model.extended_energy_sum = True
    t.smooth_order_geometry = GaugeGeometry(t,aligned=True)
    if not np.array_equal(d.fields[2,t.first][t.number],d.fields[2]):
        raise ValueError('Order copies disagree')
    return d,t,SpatialCrystal(t.b),order_grid(t,d.fields)


def save_raw(out, **arrays):
    raw = out/'raw_signed.npz'
    with raw.open('xb') as stream:
        np.savez_compressed(stream,**arrays)
    return digest(raw)


def label_bank(d,t,solid,order,out,component):
    from mouth_r4_aligned_crystal import AlignedGauge
    from mouth_r4_eulerian_relax_labels import LabelRelaxation
    from mouth_r4_full_label_work import coloring
    label = np.asarray(d.labels,LD)
    reg = LabelRelaxation(t.b,d.labels,order,d.epsilon)
    cm = reg._components()[component]
    if np.any(cm.free&cm.dead) or not np.array_equal(cm.free|cm.dead,~cm.held):
        raise ArithmeticError('Incomplete free label partition')
    gauge = SignedGauge(t,d.reference_depth,order)
    native = AlignedGauge(t,d.reference_depth,order)
    zero = np.zeros(6,LD)
    elastic,ge,_,_ = solid.evaluate(d.labels,order)
    gauge_cell,gg = native.evaluate(d.labels,zero,gradient=True)
    anchor_g = gauge.evaluate(label,zero)
    gauge_anchor_error = np.max(abs(gauge_cell-anchor_g))
    if gauge_anchor_error>1e-12:
        raise ArithmeticError(f'Independent cell gauge energies differ by {gauge_anchor_error}')
    analytic = np.zeros((2,len(cm.held)), LD)
    for dest,grad in zip(analytic,(ge,gg)):
        np.add.at(dest,cm.groups,grad[...,component].ravel())
    ids = gauge.ids
    group = cm.groups[ids]
    colors = coloring(cm,ids)
    active = reg.active.ravel()
    banks = CrystalBanks(t.b,order,out)
    anchor_p = banks.energy(label,active)
    physical_anchor_error = anchor_p.sum()-elastic
    if abs(physical_anchor_error)>1e-12:
        raise ArithmeticError(f'Independent crystal energy differs by {physical_anchor_error}')
    _,smallest = direct_rule(label.reshape(-1,3)[ids], t.smooth_order_geometry.widths, minima=True)
    tail = np.zeros(len(cm.held),bool)
    tail[group[smallest<TAIL].ravel()] = True
    tail &= ~cm.held
    min_chart = np.inf

    def sweep(mask,steps,crystal,origin):
        nonlocal min_chart
        # width, physical/unit gauge, negative/positive, joined scalar
        slopes = np.zeros((len(steps),2,2,len(cm.held)),LD)
        for color in np.unique(colors[mask]):
            chosen = mask & (colors==color)
            touched = chosen[group]
            target = np.where(touched,group,-1).max(axis=1)
            selected = target>=0
            if np.any(touched & (group!=target[:,None])):
                raise ArithmeticError('Two signed DOFs share a physical cell')
            nodes = chosen[cm.groups].reshape(order.shape)
            for k,h in enumerate(steps):
                for side,dh in enumerate((-h,h)):
                    moved = label.copy()
                    moved[...,component][nodes] += dh
                    cert = chart_orientation(reg.axes,moved,reg.active)
                    min_chart = min(min_chart,cert['chart_det_global_lower'])
                    crystal_work = crystal.energy(moved,selected & active)[selected]-origin[selected]
                    gauge_work = gauge.evaluate(moved,zero,enabled=selected)-anchor_g[selected]
                    for term,work in enumerate((crystal_work,gauge_work)):
                        np.add.at(slopes[k,term,side],target[selected],work/dh)
        return slopes

    # The same pair of broad widths is used outside every group incident on
    # the small-J tail. Large signed steps can invert that tail. Those groups
    # receive the independent narrow pair from the start.
    other = ~cm.held & ~tail
    broad = sweep(other,BROAD,banks,anchor_p)
    def errors(signed):
        central = signed.mean(axis=2)
        return (central[:,0]-analytic[0])+d.epsilon*(central[:,1]-analytic[1])
    narrow_mask = tail | (other & (np.max(abs(errors(broad)),axis=0)>LIMIT/2))
    final = broad.copy()
    fine_mask = np.zeros_like(narrow_mask)
    refined = np.empty((2,2,2,0),LD)
    widths = np.where(narrow_mask[None,:],np.array(NARROW)[:,None],np.array(BROAD)[:,None])
    if np.any(narrow_mask):
        narrow_dir = out/'narrow'
        narrow_dir.mkdir()
        quad = CrystalBanks(t.b,order,narrow_dir,quad=True)
        # The promoted bridge subtracts this same baseline before rounding.
        anchor = quad.energy(label,active & np.any(narrow_mask[group],axis=1))
        refined = sweep(narrow_mask,NARROW,quad,anchor)
        final[:,:,:,narrow_mask] = refined[:,:,:,narrow_mask]
        fine_mask = narrow_mask & (np.max(abs(errors(final)),axis=0)>LIMIT)
        if fine_mask.any():
            fine = sweep(fine_mask,FINE,quad,anchor)
            final[:,:,:,fine_mask] = fine[:,:,:,fine_mask]
            widths[:,fine_mask] = np.array(FINE)[:,None]
    err = errors(final)
    def witness(i):
        signs = final[:,0,:,i] + d.epsilon*final[:,1,:,i]
        return dict(group=int(i),node=list(map(int,np.unravel_index(cm.representative[i],order.shape))),
                    analytic=float(analytic[0,i]+d.epsilon*analytic[1,i]),
                    physical=float(analytic[0,i]),unit_gauge=float(analytic[1,i]),
                    widths=widths[:,i].astype(float).tolist(),negative_positive=signs.astype(float).tolist(),
                    central_minus_analytic=err[:,i].astype(float).tolist(),tail=bool(tail[i]))
    summary = {}
    for name,mask in [('all',~cm.held),('active',cm.free),('fully_melted',cm.dead),('small_J_tail',tail)]:
        indexes = np.flatnonzero(mask)
        summary[name] = dict(count=len(indexes),max_abs_by_width=np.max(abs(err[:,mask]),axis=1,initial=0).astype(float).tolist(),
                             worst=witness(indexes[np.argmax(abs(err[-1,mask]))]) if len(indexes) else None)
    failures = np.flatnonzero((np.max(abs(err),axis=0)>LIMIT)&~cm.held)
    return dict(component=component,checks=summary,threshold=LIMIT,passed=bool(not len(failures) and np.isfinite(err[:,~cm.held]).all()),
                failures=[witness(i) for i in failures],colors=int(colors.max()+1),narrow=int(narrow_mask.sum()),
                additional_fine=int(fine_mask.sum()),narrow_widths=list(map(float,NARROW)),
                additional_widths=list(map(float,FINE)),
                tail_cells=int((smallest<TAIL).sum()),tail_threshold=TAIL,minimum_sampled_determinant=float(smallest.min()),
                minimum_signed_whole_chart_bound=float(min_chart),gauge_cell_anchor_max=float(gauge_anchor_error),
                physical_anchor_difference=float(physical_anchor_error),
                raw_sha256=save_raw(out,analytic=analytic,broad_signed=broad,final_signed=final,negative_positive=[-1,1],
                    broad_widths=BROAD,final_widths=widths,tail=tail,free=cm.free,melted=cm.dead,narrow=narrow_mask,
                    narrow_signed=refined[:,:,:,narrow_mask] if narrow_mask.any() else refined,
                    narrow_widths=NARROW,additional_fine=fine_mask,
                    representative=cm.representative,colors=colors,cell_min_det=smallest))


def independent_linear_load(t,labels,order,out):
    """Integrate each native crystal order basis, including fully melted cells."""
    banks = CrystalBanks(t.b,order,out)
    shape = order.shape
    sites = t.number[t.b.ids.transpose(0,2,1)]
    cells = np.indices(np.subtract(shape,1)).reshape(3,-1).T
    enabled = np.ones(len(cells),bool)
    joined = np.zeros(len(t.first),LD)
    coords = np.indices(shape)
    for parity in product((0,1),repeat=3):
        indicator = np.all(coords%2 == np.array(parity)[:,None,None,None],axis=0)
        banks.order[:] = indicator
        values = banks.energy(labels,enabled)
        node = cells + ((parity-cells)%2)
        np.add.at(joined,sites[tuple(node.T)],values)
    return joined + np.asarray(t.pull.T@t.ext,LD)


def order_bank(d,t,solid,order,out):
    from mouth_r4_smooth_fields import SmoothOrder,material_load
    from mouth_r4_smooth_order_boundary import PhysicalOrderEnergy
    from mouth_r4_smooth_order_derivatives import GaugeOrderBanks
    labels,fields = d.labels,d.fields
    t.model.material_load,_ = material_load(t,solid,labels,fields)
    _,g = t.model.evaluate(fields)
    c = fields[2,t.first]
    ids = np.flatnonzero(t.free)
    solver = SmoothOrder(t,labels,d.epsilon)
    energy,unit = solver.gauge(c)
    analytic = np.stack(((t.pull.T@g[2])[ids],unit[ids])).astype(LD)
    physics = PhysicalOrderEnergy(t,fields)
    physics.load = independent_linear_load(t,np.ascontiguousarray(labels,LD),order,out)
    gauge = GaugeOrderBanks(t,labels,fields,barrier=direct_rule)
    anchor = np.sum(gauge.weight*(1-gauge.c)**2*(1+2*gauge.c))-energy
    if abs(anchor)>1e-11:
        raise ArithmeticError(f'Independent order-gauge anchor differs by {anchor}')
    signed = np.empty((len(ORDER_WIDTHS),2,2,len(ids)),LD)
    for j,node in enumerate(ids):
        p,gu = physics.at(node),gauge.at(node)
        for k,h in enumerate(ORDER_WIDTHS):
            for side,dh in enumerate((-h,h)):
                signed[k,0,side,j],signed[k,1,side,j] = p(dh),gu(dh)
    exact = analytic[0]+d.epsilon*analytic[1]
    both = signed[:,0]+d.epsilon*signed[:,1]
    central = both.mean(axis=1)
    sites = t.number[t.b.ids.transpose(0,2,1)]
    geom = t.smooth_order_geometry
    _,dets = direct_rule(labels.reshape(-1,3)[geom.ids],geom.widths,minima=True)
    incident = np.unique(sites.ravel()[geom.ids[dets<TAIL]])
    groups = {'all':np.ones(len(ids),bool),'lower':c[ids]==0,'interior':(c[ids]>0)&(c[ids]<1),
              'upper':c[ids]==1,'small_J_tail':np.isin(ids,incident)}
    summary = {}
    def compare(values,mask):
        indices = np.flatnonzero(mask)
        if not len(indices):
            return None
        err = values-exact
        ix = indices[np.argmax(abs(err[indices]))]
        return dict(max_abs=float(abs(err[ix])),signed_error=float(err[ix]),direct=float(values[ix]),
                    analytic=float(exact[ix]),physical=float(analytic[0,ix]),unit_gauge=float(analytic[1,ix]),
                    node=int(ids[ix]))
    for name,mask in groups.items():
        rows = []
        for k,h in enumerate(ORDER_WIDTHS):
            row = dict(h=float(h),central=compare(central[k],mask),
                       negative=compare(both[k,0],mask),positive=compare(both[k,1],mask))
            if k and name in ('lower','upper'):
                side = 1 if name=='lower' else 0
                ratio = ORDER_WIDTHS[k-1]/h
                row['feasible_extrapolated'] = compare((ratio*both[k,side]-both[k-1,side])/(ratio-1),mask)
            rows.append(row)
        summary[name] = dict(count=int(mask.sum()),widths=rows)
    checks = [summary['all']['widths'][-1]['central']]
    checks += [summary[name]['widths'][-1]['feasible_extrapolated'] for name in ('lower','upper') if groups[name].any()]
    return dict(checks=summary,threshold=1e-7,passed=all(q['max_abs']<=1e-7 for q in checks) and bool(np.isfinite(both).all()),
                unit_gauge_energy_difference=float(anchor),tail_threshold=TAIL,
                raw_sha256=save_raw(out,ids=ids,order=c[ids],analytic=analytic,signed=signed,widths=ORDER_WIDTHS,
                                    negative_positive=[-1,1],small_det=groups['small_J_tail']))


def shape_bank(d,t,solid,order):
    from mouth_r4_aligned_crystal import AlignedGauge
    from mouth_r4_smooth_endpoint import shape_work,SHAPE_WORK_TOL
    energy,grad,_,_ = solid.evaluate(d.labels,order)
    rows = shape_work(d,t,solid,order,energy,grad,gauge=AlignedGauge(t,d.reference_depth,order),
                      energy_gauge=SignedGauge(t,d.reference_depth,order))
    return dict(rows=rows,threshold=SHAPE_WORK_TOL,
                passed=len(rows)==6 and all(abs(b['minus_analytic'])<=SHAPE_WORK_TOL for r in rows for b in r['central']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study',type=Path)
    parser.add_argument('case',choices=EPS)
    parser.add_argument('sector',choices=('label0','label1','label2','order','shape'))
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    study = args.study.resolve()
    d,t,solid,order = load(study,args.case,out)
    if args.sector.startswith('label'):
        result = label_bank(d,t,solid,order,out,int(args.sector[-1]))
    elif args.sector=='order':
        result = order_bank(d,t,solid,order,out)
    else:
        result = shape_bank(d,t,solid,order)
    src = Path(__file__).resolve().parent
    modules = (Path(m.__file__).resolve() for m in sys.modules.values() if getattr(m,'__file__',None))
    used = {p for p in modules if p.parent==src}
    used.update(p for p in src.iterdir() if p.suffix in ('.cpp','.hpp'))
    result.update(state_sha256=digest(study/args.case/'rematched.npz'),case=args.case,epsilon=d.epsilon,
                  sector=args.sector,command=sys.orig_argv,mantissa_bits=int(np.finfo(LD).nmant),
                  source_sha256=digest(Path(__file__)),sources={p.name:digest(p) for p in sorted(used)})
    with (out/'record.json').open('x') as stream:
        json.dump(result,stream,allow_nan=False,indent=2)
        stream.write('\n')
    print(json.dumps(result,allow_nan=False),flush=True)
    if not result['passed']:
        raise SystemExit(1)


if __name__=='__main__':
    main()
