"""Fixed-order, fixed-mouth tests at the last accepted radius-four barrier state."""
import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from mouth_cp10_label_condition import physical_frame
from mouth_cp10_label_grid import get_cp10
from mouth_exact_checkpoint_control import digest
from mouth_r4_corner_projection import BITS, point_normal, packed_lift, project_halfspaces
from mouth_r4_eulerian_relax_labels import InfeasibleLabels, LabelRelaxation, chart_orientation
from mouth_r4_eulerian_solid import SpatialCrystal
from mouth_r4_melt_barrier import MeltBarrier
from mouth_r4_stage2_outcome import verify_guard
from mouth_supercore_shape import preparation

CELL = (15, 9, 1)
LOCAL = (1, 1, 0)
SHA = '418376c3eb1922f6335543f43dfee11e531cebbb6958d37afd73e844571b653a'


def first_root(coef):
    roots = np.polynomial.polynomial.polyroots(coef)
    positive = [z.real for z in roots if abs(z.imag)<1e-10*max(1.,abs(z.real)) and z.real>0]
    return float(min(positive)) if positive else None


def local_polynomial(opt, labels, delta):
    vertex = np.array(CELL) + LOCAL
    def matrix(f):
        edges = []
        for d, sign in enumerate(1-2*np.array(LOCAL)):
            other = vertex.copy()
            other[d] += sign
            h = opt.axes[d][other[d]]-opt.axes[d][vertex[d]]
            edges.append((f[tuple(other)]-f[tuple(vertex)])/h)
        return np.array(edges)
    M, V = matrix(labels), matrix(delta)
    co = np.zeros(4)
    for bits in BITS:
        rows = np.where(bits[:, None], V, M)
        co[bits.sum()] += rows[0] @ np.cross(rows[1],rows[2])
    return co


def determinant_controls(axes, labels):
    """Whole-mesh Bernstein coefficients and orientation roundoff scale."""
    dims = tuple(len(a)-1 for a in axes)
    slices = lambda b: tuple(slice(k,k+n) for k,n in zip(b,dims))
    cols = []
    norms = np.ones(dims)
    for d, ax in enumerate(axes):
        tangents = [i for i in range(3) if i!=d]
        w = np.diff(ax).reshape(tuple(dims[i] if i==d else 1 for i in range(3))+(1,))
        corners = []
        for a,b in product(range(2), repeat=2):
            base = [0]*3
            base[tangents[0]],base[tangents[1]] = a,b
            left = labels[slices(base)]
            base[d] = 1
            corners.append((labels[slices(base)]-left)/w)
        cols.append([(1-u)*((1-v)*corners[0]+v*corners[1])+u*((1-v)*corners[2]+v*corners[3])
                     for u,v in product((0.,.5,1.), repeat=2)])
        norms *= np.max([np.linalg.norm(c,axis=-1) for c in corners],axis=0)
    out = np.empty(dims+(3,3,3))
    for i,j,k in product(range(3), repeat=3):
        out[...,i,j,k] = np.einsum('...i,...i->...',cols[0][3*j+k], np.cross(cols[1][3*i+k],cols[2][3*i+j]))
    for k in [3,4,5]:
        view = np.moveaxis(out,k,-1)
        view[...,1] = 2*view[...,1]-.5*(view[...,0]+view[...,2])
    return out,norms


def certify_path(opt, labels, delta, alpha):
    """Sufficient joint space/step determinant bound for the entire continuous prefix."""
    from mouth_r4_eulerian_relax_labels import certify_subcell
    samples = [determinant_controls(opt.axes, labels + (a*alpha)*delta) for a in [0.,1/3,2/3,1.]]
    p0,p1,p2,p3 = [s[0] for s in samples]
    c1 = 3*p1 - 1.5*p2 - (5/6)*p0 + p3/3
    c2 = 3*p2 - 1.5*p1 - (5/6)*p3 + p0/3
    margin = np.maximum.reduce([s[1] for s in samples])*4096*np.finfo(float).eps
    bound = np.inf
    for sheet in (p0,c1,c2,p3):
        lower = sheet.min(axis=(-3,-2,-1))-margin
        for row in np.argwhere(lower<=0):
            cell = tuple(row)
            lower[cell] = certify_subcell(sheet[cell],margin[cell],list(map(int,cell)))
        bound = min(bound,float(lower.min()))
    return bound


def ray_limit(opt, labels, delta):
    """Bracket first invalid fraction and certify the whole continuous prefix below it."""
    def test(alpha):
        try:
            return True,chart_orientation(opt.axes,labels+alpha*delta,opt.active)
        except InfeasibleLabels as exc:
            return False,dict(exc.witness)
    lo,hi = 0.,1.
    ok,issue = test(hi)
    while ok and hi < 1e4:
        lo,hi = hi,hi*2
        ok,issue = test(hi)
    if ok:
        raise ArithmeticError('No upper bound in measured ray')
    for _ in range(34):
        mid = (lo+hi)/2
        ok,item = test(mid)
        if ok:
            lo = mid
        else:
            hi,issue = mid,item
    # The witness is an upper bound on first loss. Certify that there are no
    # hidden earlier inversions anywhere in space or along the trial ray.
    for factor in (1-1e-8,1-1e-6,.9999,.99,.95,.9,.5):
        try:
            lower = lo*factor
            value = certify_path(opt,labels,delta,lower)
            if value>0:
                return dict(positive_through=lower,whole_path_chart_lower=value,
                            last_feasible_endpoint=lo,first_rejection_upper=hi,upper_issue=issue)
        except InfeasibleLabels:
            pass
    raise ArithmeticError('Endpoint brackets did not certify a positive continuous ray')


def coordinate(opt,node):
    r = np.ravel_multi_index(node,opt.stencil.shape)
    comp = opt._components()[0]
    return list(np.flatnonzero(comp.free)).index(comp.groups[r])


def angle(normal, direction):
    cos = normal@direction/(np.linalg.norm(normal)*np.linalg.norm(direction))
    return float(np.rad2deg(np.arccos(np.clip(cos,-1.,1.))))


def nearest_points(opt, solid, labels, chi, ref, witness):
    """Native cut radial quadrature in the lone ordered cell incident to the witness."""
    ids = BITS+CELL
    nodes = labels[tuple(ids.T)]
    if nodes[:,0].max() >= opt.body.cut:
        raise ValueError('This core-cell reconstruction would require material-registry collar cuts')
    def position(row):
        chart = np.array(row['point'])
        frame,_ = physical_frame(opt.body,chart[None])
        r = frame[0,2,2]
        return np.array([r*np.sin(chart[2]),r*np.cos(chart[2]),chart[1]*opt.body.sh.length])
    origin = position(witness)
    rows = []
    for t,a in product(solid.points,repeat=2):
        borders = []
        for r in (0,1):
            l = np.array([r,t,a]); b = np.where(BITS,l,1-l).prod(axis=1)
            borders.append(b@nodes)
        left,right = borders
        cuts = [0.,1.]
        for k,key in [(0,'x'),(1,'t'),(2,'a'),(1,'reg_t')]:
            if left[k]!=right[k]:
                ax = opt.body.buffers[key]
                cuts.extend(((ax[(ax>min(left[k],right[k]))&(ax<max(left[k],right[k]))]-left[k])/(right[k]-left[k])).tolist())
        cuts = np.unique(cuts)
        for l,h in zip(cuts[:-1],cuts[1:]):
            for x in l+(h-l)*solid.points:
                _,row = point_normal(opt,labels,CELL,(x,t,a),ref)
                row['chi'] = float(np.array(row['basis'])@chi[tuple(ids.T)])
                row['distance'] = float(np.linalg.norm(position(row)-origin))
                if row['chi']>0: rows.append(row)
    return sorted(rows,key=lambda row:row['distance'])[:3]


def run(scratch,out):
    out.mkdir(exist_ok=True,parents=True)
    path = scratch/'r4_barrier_held01/primary/held_stop.npz'
    guard = scratch/'r4_barrier_held01/primary_guard'
    verify_guard(guard,1)
    if digest(path)!=SHA:
        raise ValueError('Wrong saved active crystal')
    saved = json.loads(path.with_suffix('.json').read_text())
    body,cp,_ = get_cp10(scratch)
    with np.load(path) as f:
        if not np.array_equal(f['x'],cp['x']) or f['radius']!=4 or f['winding']!=1:
            raise ValueError('Geometry or winding differs')
        labels = f['labels']
        order = np.ascontiguousarray(f['field'][2,body.ids].transpose(0,2,1))
        eps = float(f['epsilon'])
    opt = LabelRelaxation(body,labels,order,eps)
    solid = SpatialCrystal(body)
    barrier = MeltBarrier(body,opt,out)
    ref = preparation(4)[0]
    elastic,ge,_,diag = solid.evaluate(labels,order)
    unit,gg,_,_ = barrier.sample(labels)
    g = barrier.extract(ge+eps*gg,'free')
    metric = opt._precondition(order)
    d = -metric(g)
    basis,witness = point_normal(opt,labels,CELL,LOCAL,ref)
    chart,_ = point_normal(opt,labels,CELL,LOCAL)
    if abs(witness['determinant']-.0005014362395531604)>3e-12:
        raise ArithmeticError('Different physical inverse volume')
    witness['chi'] = float(order[tuple(np.array(CELL)+LOCAL)])
    motion = packed_lift(opt,d)
    poly = local_polynomial(opt,labels,motion)
    section = np.flatnonzero(basis) # exclude zero components at the other four vertices
    readback = out/'corner.json'
    if readback.exists():
        results = json.loads(readback.read_text())
        if results['source_sha256']!=digest(path) or results['metric_restart']['actual_chart_polynomial']!=poly.tolist():
            raise ValueError('Retained angular readout changed source or direction')
    else:
        corner_rec = dict(accepted=witness, closest_ordered_points_in_incident_cell=nearest_points(opt,solid,labels,order,ref,witness))
        results = dict(source_sha256=SHA, epsilon=eps,shape_sha256=digest(scratch/'r4_melted_survival01/coarse/state_00010.npz'),
                       initial_material_force=float(np.max(np.abs(g))),corner=corner_rec)
        all_probes=saved['stop']['probes']
        def history(entries):
            hrows = [r for r in entries if r['predictor']=='harmonic_increment']
            data = [r for r in hrows if r.get('witness',{}).get('cell')==list(CELL)
                    and r['witness'].get('chart_fractions')==list(LOCAL)]
            a=np.array([r['alpha'] for r in data]); y=np.array([r['witness']['chart_det'] for r in data])
            norm=a.max(); coeff=np.polynomial.polynomial.polyvander(a/norm,3)
            fit=np.linalg.lstsq(coeff[:,1:],y-poly[0],rcond=None)[0]
            co=np.r_[poly[0],fit/norm**np.arange(1,4)]
            # A sampled amount passed domain checks even when its first det witness was elsewhere.
            safe=max(r['alpha'] for r in hrows if 'cell' in r.get('witness',{}))
            z=np.array(witness['label']); beta=np.array([2/z[0], 3*ref.length*ref.radius(z[1]*ref.length)[1]/ref.radius(z[1]*ref.length)[0], 1/np.tan(z[2])])
            bounds=np.array([[(ax[0]-v-2e-13)/safe,(ax[-1]-v+2e-13)/safe] for ax,v in zip(opt.axes,z)])
            metric_work=np.sort(witness['determinant']*bounds*beta[:,None],axis=1).sum(axis=0)
            return dict(chart_polynomial=co.tolist(), fractions=a.tolist(), saved_chart_determinants=y.tolist(),
                        max_chart_fit_error=float(np.max(abs(np.polynomial.polynomial.polyval(a,co)-y))),
                        chart_slope=float(co[1]),physical_slope_bounds=(witness['metric_factor']*co[1]+metric_work).tolist(),
                        first_corner_zero=first_root(co),largest_domain_passing_alpha=safe)
        split=len(all_probes)//2
        retained,restarted=history(all_probes[:split]),history(all_probes[split:])
        retained['angle_degrees']=None
        retained['identifiability']='Retained L-BFGS direction and secants were not saved; one scalar determinant polynomial cannot identify the angle or full direction norm.'
        restarted.update(actual_chart_polynomial=poly.tolist(),actual_physical_slope=float(basis@d),
                         actual_chart_slope=float(chart@d),global_packed_angle_deg=angle(basis,d),
                         local_angles_definition='Euclidean packed label components on the four supporting vertices',
                         local_angle_deg=angle(basis[section],d[section]),g_dot_d=float(g@d))
        restarted['saved_probe_max_error']=float(np.max(np.abs(np.polynomial.polynomial.polyval(restarted['fractions'],poly)-restarted['saved_chart_determinants'])))
        if restarted['saved_probe_max_error']>1e-10:raise ArithmeticError('Restart direction does not reproduce rejected probes')
        results.update(retained_history=retained,metric_restart=restarted)
        out.joinpath('corner.json').write_text(json.dumps(results,indent=2)+'\n')
    ep0=elastic+eps*unit
    def trial(direction,fraction,tag):
        dest=labels+fraction*packed_lift(opt,direction)
        before=chart_orientation(opt.axes,dest,opt.active)
        inner,meta=barrier.solve(dest,tag,out)
        if not meta['stationary']:raise ArithmeticError('Nonstationary melted extension')
        e,_,_,diagnostics=solid.evaluate(inner,order,False)
        u,_,_,_=barrier.sample(inner)
        after=chart_orientation(opt.axes,inner,opt.active)
        return dict(fraction=fraction,elastic_change=e-elastic, reduced_energy_change=float(e+eps*u-ep0),
                    inner_residual=meta['residual'],min_Jinv_active=diagnostics['min_Jinv_active'],
                    initial_whole_cell_bounds=before,settled_whole_cell_bounds=after)
    for name,direction in [('metric_restart',d),('projected_restart',project_halfspaces(d,[basis],metric))]:
        if name not in results:results[name]={}
        rr=results[name]; motion=packed_lift(opt,direction)
        rr.update(g_dot_d=float(g@direction),physical_J_slope=float(basis@direction),
                  chart_J_slope=float(chart@direction),four_node_first_zero=first_root(local_polynomial(opt,labels,motion)),
                  global_predictor=ray_limit(opt,labels,motion))
        if name=='projected_restart':
            for n,factor in enumerate((.5,.125,.03125,.0078125)):
                alpha=min(1.,rr['global_predictor']['positive_through']*factor)
                data=trial(direction,alpha,'projected_'+str(n))
                if data['reduced_energy_change']<0:
                    rr['settled_descent']=data
                    break
            else:raise ArithmeticError('Projected direction lacks a measured downhill step')
        out.joinpath('partial.json').write_text(json.dumps(results,indent=2)+'\n')
    index=int(np.argmax(abs(g)))
    node=(17,1,19)
    if coordinate(opt,node)!=index:raise ArithmeticError('Entrance force moved')
    unitdir=np.zeros_like(g); unitdir[index]=-np.sign(g[index])
    slope=abs(float(d[index]))
    entrance={'node':list(node),'force':float(-g[index]),'restart_displacement_per_fraction':float(d[index]),
              'unit_label_direction':float(unitdir[index]),'corner_slope_for_entrance_only':float(basis@unitdir)}
    move=packed_lift(opt,unitdir)
    if np.any(move[opt.stencil].reshape(-1,3)[:,1:]) or np.count_nonzero(move)>1:
        raise ArithmeticError('Chosen entrance coordinate touches additional groups/melt')
    entrance['unit_label_global_limit']=ray_limit(opt,labels,move)
    def residual(t):
        q=labels+t*move
        _,a,_,_=solid.evaluate(q,order)
        _,b,_,_=barrier.sample(q)
        return float(barrier.extract(a+eps*b,'free')[index]*unitdir[index])
    stop_fraction=.9*entrance['unit_label_global_limit']['positive_through']
    r=residual(stop_fraction)
    if r>0:
        minimum=brentq(residual,0.,stop_fraction,xtol=5e-14)
        entrance.update(local_force_zero_label_displacement=minimum,local_force_zero_metric_fraction=minimum/slope,
                        residual_at_local_force_zero=residual(minimum),settled_descent=trial(unitdir,minimum,'entrance_only'))
    else:
        entrance.update(probed_label_displacement=stop_fraction,remaining_signed_gradient=r,
                        settled_descent=trial(unitdir,stop_fraction,'entrance_only'))
    results['entrance_only']=entrance
    results['conditioning_note']='Ray upper/lower limits concern the specified line and include label-domain bounds. A melt re-equilibration can change purely melted constraints. The four active corner nodes do not change in that subsolve.'
    out.joinpath('measurement.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('scratch',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    run(args.scratch,args.output)
