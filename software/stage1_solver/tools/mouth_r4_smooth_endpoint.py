"""Read and compare the two stationary smooth-gauge R4 states with shape fixed."""
import argparse
import copy
import ctypes as ct
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from mouth_cp10_field_shape import Bath
from mouth_cp10_label_gauge import depth_grid
from mouth_cp10_nested_forces import integral, kernel
from mouth_cp10_spatial_crystal import MODES
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_exchange_carrier import BIAS
from mouth_r4_depth_precision import extended
from mouth_r4_eulerian_relax_labels import FORCE_TOL, LabelRelaxation, absmax, chart_orientation
from mouth_r4_eulerian_release import order_grid
from mouth_r4_fold_order import all_low_after
from mouth_r4_slaved_melt import FLOW
from mouth_r4_smooth_label_derivatives import CrystalBanks
from mouth_r4_smooth_order_boundary import PhysicalOrderEnergy, gauge_coefficients
from mouth_r4_smooth_static import EXPECTED, Gauge, LD
from mouth_r4_stage2_outcome import verify_guard
from mouth_whole_r4_relax import atomic_json, save_fields

INNER_TARGET = 3e-8
SHAPE_WORK_TOL = 1e-7
ORDER_STEP = LD(1e-9)
MAX_RELATIVE_BOUNDARY_GAUGE_WORK = 1e-4
BASELINE = dict(mean_width=.705499, scaled_log_exit=-.641236)


def promoted_bath(bath, trial, reference, dx):
    """Retain the tabulated vortex and promote the two extensive bath order terms."""
    offset = np.array(dx, dtype=np.result_type(dx, LD))
    native_offset = np.asarray(dx, dtype=complex if np.iscomplexobj(dx) else float)
    parts = bath.energy(native_offset)
    shape = copy.copy(trial.shape)
    shape.params = shape.params + offset[:4]*[1.,.25,.25,.25]
    shape.extra += .25*offset[4]
    bank = copy.copy(bath)
    bank.model = copy.copy(bath.model)
    bank.model.polar_shift += offset[5]
    parts[3:5] = extended(bank,shape,LD(trial.shape.length)*np.exp(offset[0]),reference)
    return parts


def shape_work(driver, t, solid, order, energy, label_gradient):
    labels = driver.labels
    epsilon = driver.epsilon
    gauge = Gauge(t,driver.reference_depth,order)
    bath = Bath(t,driver.base,driver.fields,driver.reference_depth)
    length = LD(t.shape.length)
    depth = gauge.axes[1]
    dz = np.where(depth*length < length-.12, (depth*length/(length-.12))**25, 1.)-depth
    with TemporaryDirectory(dir=driver.folder) as directory:
        tmp = Path(directory)
        native = integral(kernel(tmp),t.b,gauge.axes,np.asarray(dz,float),labels,order,
                          (solid.points,solid.weights),True)
        if abs(native['energy']-energy)>1e-11:
            raise ArithmeticError('Shape work and physical elastic anchor disagree')
        elastic = np.array(native['analytic'],dtype=LD)
        elastic[0] -= np.sum(label_gradient[...,1]*dz[None,:,None],dtype=LD)
        banks = CrystalBanks(t.b,order,tmp)
        active = np.max(gauge.order,axis=1)>0
        tbuf = banks.buffers['field_t'].copy()
        sbuf = banks.buffers['shape'].copy()
        zero = np.zeros(6,dtype=LD)
        origin = depth_grid(length,driver.reference_depth)

        def set_scalar(name,value):
            val = np.array(value,LD)
            ct.memmove(ct.addressof(banks.data)+getattr(type(banks.data),name).offset,val.ctypes.data,val.nbytes)

        def energies(dx):
            L = length*np.exp(dx[0])
            dt = depth_grid(L,driver.reference_depth)-origin
            banks.buffers['field_t'][:] = tbuf+dt
            set_scalar('length',L)
            set_scalar('epsilon',LD(t.polar_shift)+dx[5])
            shape = banks.buffers['shape']
            shape[:] = sbuf
            shape[0] += (sbuf[0]+LD(1.5))*np.expm1(dx[1]/4)
            shape[1:4] += dx[2:5]/4
            xi = labels.astype(LD)
            xi[...,1] += dt[None,:,None]
            return (banks.energy(xi,active),promoted_bath(bath,t,driver.reference_depth,dx),
                    gauge.evaluate(labels,dx))

        base = energies(zero)
        if abs(np.sum(base[0])-energy)>1e-11:
            raise ArithmeticError('Promoted finite-difference crystal anchor differs')
        rows = []
        for i,unit in enumerate(np.eye(6,dtype=LD)):
            tangent = unit.astype(np.clongdouble)*1.e-24j
            ext = -np.sum(promoted_bath(bath,t,driver.reference_depth,tangent)).imag/LD(1.e-24)
            reg = -np.sum(gauge.evaluate(labels,tangent)).imag/LD(1.e-24)
            force = elastic[i]+ext+epsilon*reg
            signed = []
            for h in (LD(2)**-18, LD(2)**-20):
                minus,plus = energies(-h*unit),energies(h*unit)
                cs,bs,gs = [np.sum(m-p,dtype=LD)/(2*h) for m,p in zip(minus,plus)]
                central = cs+bs+epsilon*gs
                signed.append(dict(h=float(h),central=float(central),minus_analytic=float(central-force),
                                   crystal=float(cs),bath=float(bs),unit_gauge=float(gs)))
            rows.append(dict(mode=MODES[i],analytic=float(force),elastic=float(elastic[i]),
                             bath=float(ext),unit_gauge=float(reg),central=signed))
    return rows


def boundary_sweep(t, fields, coefficient, physical_gradient, epsilon):
    physical = PhysicalOrderEnergy(t,fields)
    chi = fields[2,t.first]
    free_zero = np.flatnonzero((chi==0.) & t.free)
    physical_work = np.empty(len(free_zero),dtype=LD)
    max_error = 0.
    for j,ix in enumerate(free_zero):
        increment = physical.at(ix)
        max_error = max(max_error, abs(float(increment(0)-physical_gradient[ix])))
        physical_work[j] = ORDER_STEP*increment(ORDER_STEP)
    if max_error>5e-10:
        raise ArithmeticError('Boundary energy path and physical order derivative disagree')
    a,b,c = coefficient[:,free_zero]
    numerical_work = LD(epsilon)*ORDER_STEP*(a+ORDER_STEP*(b+ORDER_STEP*c))
    total = physical_work+numerical_work
    positive = physical_work>0
    positive_rows = np.flatnonzero(positive)
    fraction = np.abs(numerical_work[positive])/physical_work[positive]
    largest_ratio = float(np.max(fraction,initial=0.))
    negligible = largest_ratio <= MAX_RELATIVE_BOUNDARY_GAUGE_WORK
    bad = (total<0) | (positive & (total<=0))
    incident = (t.V.T @ ((physical.c>0)&(t.model.hq>0)).astype(float))[free_zero]>0
    def witness(k):
        return dict(joined_node=int(free_zero[k]),physical=float(physical_work[k]),
                    gauge=float(numerical_work[k]),total=float(total[k]),incident_positive=bool(incident[k]))
    return dict(step=float(ORDER_STEP),free_zero=len(free_zero),front=int(incident.sum()),
                near_front_with_linear_term=int(np.sum(incident & (a!=0))),
                max_physical_gradient_error=max_error,physical_positive=int(positive.sum()),
                max_abs_gauge_fraction=largest_ratio,relative_limit=MAX_RELATIVE_BOUNDARY_GAUGE_WORK,
                worst_fraction_node=witness(positive_rows[np.argmax(fraction)]) if len(positive_rows) else None,
                negligible_gauge_work=negligible,
                minimum_work=witness(np.argmin(total)) if len(total) else None,
                physical_positive_total_nonpositive=int(np.sum(positive & (total<=0))),
                negative_total=int(np.sum(total<0)),nonfinite=int(np.sum(~np.isfinite(total))),
                failures=[witness(k) for k in np.flatnonzero(bad)],
                passed=bool(np.all(np.isfinite(total)) and not np.any(bad) and negligible))


def write_endpoint(driver,t,solid,details):
    labels, fields, eps = driver.labels,driver.fields,driver.epsilon
    if eps not in (1e-8,1e-10) or not details['material']['stationary']:
        raise ValueError('Stationary shape forces require the native target gauge and stationary internal labels')
    case = driver.folder
    state = case/'stationary_state.npz'
    if state.exists():
        raise FileExistsError(state)
    save_fields(state,labels=labels,field=fields,x=driver.x,radius=4,winding=1,epsilon=eps,
                iteration=driver.iteration,parent_sha256=CHECKPOINT_SHA,initial_sha256=details['starting_sha256'])
    sha = digest(state)
    order = order_grid(t,fields)
    chart = Gauge(t,driver.reference_depth,order)
    opt = LabelRelaxation(t.b,labels,order,eps)
    e_s,derivative,src,phys = solid.evaluate(labels,order)
    load = t.ext + np.bincount(t.b.ids.ravel(),src.transpose(0,2,1).ravel(),minlength=t.model.size)
    t.model.material_load=load
    energy, gr=t.model.evaluate(fields)
    currents=t.measure(fields,load)
    ec,g = chart.evaluate(labels,np.zeros(6,dtype=LD),gradient=True)
    unit = float(ec.sum())
    coeff,_,units = gauge_coefficients(t,labels,fields)
    if abs(unit-units['new_unit']) > 1e-11*max(unit,1):
        raise ArithmeticError('Independent order and label gauge energies disagree')
    chi=fields[2,t.first]
    if (not np.array_equal(chi[t.number],fields[2]) or np.any((chi<0)|(chi>1)) or
            np.any(fields[1,t.model.feed]!=0.) or np.any(fields[1,t.model.bulk]!=BIAS)):
        raise ValueError('Unmatched order, distant phase or inadmissible order coefficients')
    bare_order = t.pull.T@gr[2]
    joined = bare_order+np.asarray(eps*coeff[0],dtype=float)
    feasible = np.where(chi==0., np.minimum(joined,0.), np.where(chi==1.,np.maximum(joined,0.),joined))
    active = melt = 0.
    for c,cm in enumerate(opt._components()):
        total_component = cm.sum(derivative[...,c]+eps*np.asarray(g[...,c],float))
        active = max(active,absmax(total_component[cm.free]))
        melt = max(melt,absmax(cm.sum(np.asarray(g[...,c],float))[cm.dead]))
    cert=chart_orientation(opt.axes,labels,opt.active)
    boundary=boundary_sweep(t,fields,coeff,bare_order,eps)
    force_rows=shape_work(driver,t,solid,order,e_s,derivative)
    forces=[row['analytic'] for row in force_rows]
    front=all_low_after(t.model.z[1],order.transpose(0,2,1),0.)
    residuals=dict(label=active,unit_melt=melt,order_kkt=absmax(feasible[t.free]),
                   amplitude=absmax(gr[0,t.model.masks[0]]),phase=absmax(gr[1,t.model.masks[1]]))
    stepped = details['smallest_accepted_active_fraction']
    checks=dict(internal=max(active,residuals['order_kkt'],residuals['phase'],residuals['amplitude'])<=FORCE_TOL,
                unit_melt=melt<=INNER_TARGET,whole_cell_orientation=cert['chart_det_global_lower']>0,
                global_sampled_orientation=phys['min_Jinv_global']>0,ordered_sampled_orientation=phys['min_Jinv_active']>0,
                conservative_current=(currents['feed_current']>0 and
                    max(abs(currents[k]) for k in ('port_error','disk_minus_feed_current'))<=FORCE_TOL),
                negligible_boundary_gauge_work=boundary['negligible_gauge_work'],
                finite_steps=stepped is None or stepped>1e-4,boundary_work=boundary['passed'],
                energy=abs(float(energy)+eps*unit-details['jointly_matched_energy'])<=1e-7,
                all_six_signed=max(abs(b['minus_analytic']) for row in force_rows for b in row['central'])<=SHAPE_WORK_TOL)
    result=dict(state_sha256=sha,initial_sha256=details['starting_sha256'],case=details['case'],epsilon=eps,shape=driver.x.tolist(),
        depth=float(t.shape.length),exit_radius=float(t.shape.exit),phase_drop=BIAS,total_energy=float(energy)+eps*unit,
        physical_energy=float(energy),crystal_energy=e_s,gauge_energy=eps*unit,gauge_to_crystal=eps*unit/e_s,
        current={k:float(currents[k]) for k in FLOW},sampled_melted_front=front,
        order_max_by_depth=np.max(order,axis=(0,2)).tolist(),depth_knots=t.model.z[1].tolist(),native_residuals=residuals,
        min_sampled_physical_Jinv=phys['min_Jinv_global'],min_sampled_ordered_Jinv=phys['min_Jinv_active'],orientation=cert,
        smallest_accepted_active_fraction=stepped,order_boundary=boundary,shape_forces=forces,shape_audit=force_rows,
        cp10_comparison={m:dict(cp10=v,now=forces[MODES.index(m)],delta=forces[MODES.index(m)]-v) for m,v in BASELINE.items()},
        checks={k:bool(v) for k,v in checks.items()},passed=all(checks.values()))
    atomic_json(case/'endpoint.json',result)
    if not result['passed']:
        raise ArithmeticError(f'Fresh stationary endpoint readback failed: {result["checks"]}')
    atomic_json(case/'HELD.json',dict(saved=sha,endpoint=digest(case/'endpoint.json'),epsilon=eps,passed=True))
    return result


def certified_state(study, name, epsilon, initial):
    folder = study/name
    verify_guard(study/(name+'_guard'), 0)
    cert = read(folder/'HELD.json')
    observed = read(folder/'endpoint.json')
    state = folder/'stationary_state.npz'
    if not (cert['passed'] and observed['passed'] and all(observed['checks'].values())
            and cert['saved']==observed['state_sha256']==digest(state)
            and cert['endpoint']==digest(folder/'endpoint.json')
            and cert['epsilon']==observed['epsilon']==epsilon
            and observed['initial_sha256']==initial and observed['case']==name):
        raise ValueError(f'Invalid stationary endpoint provenance in {folder}')
    with np.load(state,allow_pickle=False) as data:
        if not (int(data['radius'])==4 and int(data['winding'])==1 and float(data['epsilon'])==epsilon
                and str(data['parent_sha256'])==CHECKPOINT_SHA and str(data['initial_sha256'])==initial
                and np.array_equal(data['x'],observed['shape']) and observed['phase_drop']==BIAS):
            raise ValueError(f'Fixed-shape state mismatch in {folder}')
        return observed, data['field'], data['labels']


def difference_profile(study, first, second):
    from mouth_cp10_label_condition import physical_frame, reference_frame
    from mouth_cp10_label_grid import get_cp10
    from mouth_cp10_label_orientation import basis_and_grad
    from mouth_r4_smooth_static import BITS, RULE
    from mouth_supercore_shape import preparation
    body, cp, axes=get_cp10(study.parent)
    cp.close()
    ref=preparation(4)[0]

    def gather(fields,labels):
        order=np.broadcast_to(fields[2,body.ids].transpose(0,2,1)[...,None],labels.shape)
        # Pointwise gauge density carries the eight-corner mean barrier per cell.
        barrier=0.
        for vertex in BITS:
            _, grad=basis_and_grad(labels,axes,vertex)
            det=np.linalg.det(grad)
            if np.any(det<=0):
                raise ValueError('Inverted corner on independently certified pair')
            barrier=barrier+(det-1-np.log(det))/8
        records=[]
        for pt in RULE:
            label,grad=basis_and_grad(labels,axes,pt)
            chi=basis_and_grad(order,axes,pt)[0][...,0]
            pos=np.stack(np.broadcast_arrays(*np.ix_(*[a[:-1]+p*np.diff(a) for a,p in zip(axes,pt)])),axis=-1)
            b,hoop=physical_frame(body,pos)
            c,ref_hoop=reference_frame(body,ref,label)
            metric=np.linalg.inv(b)
            inverse=c@grad@metric
            ratio=ref_hoop/hoop
            jinv=np.sign(np.linalg.det(grad))*np.abs(np.linalg.det(inverse)*ratio)
            if np.any(jinv<=0):
                raise ValueError('Inverted quadrature point on independently certified pair')
            W=.5*(np.sum(np.linalg.inv(inverse)**2,axis=(-1,-2))+1/ratio**2-4)+np.log(jinv)+(1/jinv-1)**2/4
            reg=barrier+.5*np.sum(((grad-np.eye(3))@metric)**2,axis=(-1,-2))
            records.append(np.stack((chi,jinv,W,(1-chi)**2*(1+2*chi)*reg),axis=-1))
        return np.stack(records)

    primary=gather(*first)
    secondary=gather(*second)
    output=study/'gauge_split_profile.npz'
    names=['chi','Jinv','W_crystal','s_chi_W_reg']
    save_fields(output,primary=primary,secondary=secondary,components=names,
                quadrature=np.array(RULE),chart_x=axes[0],chart_t=axes[1],chart_angle=axes[2])
    witnesses={}
    for i,name in enumerate(names):
        delta=primary[...,i]-secondary[...,i]
        idx=np.unravel_index(np.argmax(abs(delta)),delta.shape)
        witnesses[name]=dict(quadrature=int(idx[0]),cell=list(map(int,idx[1:])),primary=float(primary[(*idx,i)]),
                             secondary=float(secondary[(*idx,i)]),primary_minus_secondary=float(delta[idx]))
    return dict(file=output.name,sha256=digest(output),maximum_absolute_difference=witnesses)


def compare(study, *, primary_initial=None):
    if primary_initial is None:
        primary_initial=EXPECTED['held_stop']
    primary,first,labels_a=certified_state(study,'primary',1e-8,primary_initial)
    secondary,second,labels_b=certified_state(study,'secondary',1e-10,primary['state_sha256'])
    if primary['shape']!=secondary['shape'] or primary['depth_knots']!=secondary['depth_knots']:
        raise ValueError('Held shapes differ between the two gauge equilibria')
    a,b=first[2],second[2]
    delta=a-b
    mismatch=(a==0.)!=(b==0.)
    # Count distinct global degree-of-freedom copies as stored and record whether
    # any support change exceeds the precision of the held equilibrium.
    mismatch_count=int(np.sum(mismatch))
    mismatch_size=float(np.max(np.abs(delta[mismatch]),initial=0.))
    diff={key:primary[key]-secondary[key] for key in
          ('total_energy','physical_energy','crystal_energy','gauge_energy','gauge_to_crystal',
           'min_sampled_physical_Jinv','min_sampled_ordered_Jinv')}
    current={key:primary['current'][key]-secondary['current'][key] for key in FLOW}
    forces=np.asarray(primary['shape_forces'])-secondary['shape_forces']
    # Sampled fully melted cross-sections use identical physical knots.
    front_equal=primary['sampled_melted_front']==secondary['sampled_melted_front']
    chi_max=absmax(delta)
    checks=dict(energy=max(abs(diff[key]) for key in ('total_energy','physical_energy','crystal_energy'))<=1e-6,
                current=max(abs(current[k]) for k in FLOW if k!='total_axial_stress')<=FORCE_TOL,
                order=chi_max<=FORCE_TOL,front=front_equal,support=mismatch_size<=FORCE_TOL,
                stationary_shape_forces=absmax(forces)<=FORCE_TOL)
    result=dict(input_states={k:r['state_sha256'] for k,r in [('primary',primary),('secondary',secondary)]},
                energies_and_orientation=diff,current=current,shape_force_delta=dict(zip(MODES,forces.tolist())),
                max_absolute_order_difference=chi_max,order_difference_node=int(np.argmax(abs(delta))),
                support_difference_full_array_count=mismatch_count,maximum_order_in_support_difference=mismatch_size,
                sampled_front=[primary['sampled_melted_front'],secondary['sampled_melted_front']],
                native_residuals=[primary['native_residuals'],secondary['native_residuals']],
                limits=dict(energy=1e-6,force_order_current=FORCE_TOL,front='same physical cross-section'),
                signed_differences='primary (1e-8) minus secondary (1e-10)',checks=checks,
                status='Held for independent review' if all(checks.values()) else 'Gauge split at held shape',
                passed=all(checks.values()))
    if not result['passed']:
        result['local_profile']=difference_profile(study,(first,labels_a),(second,labels_b))
    atomic_json(study/'paired_held.json',result)
    if not result['passed']:
        raise ArithmeticError('Held gauge endpoints differ; see paired_held.json')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation',choices=('pair',))
    parser.add_argument('study',type=Path)
    args=parser.parse_args()
    print(json.dumps(compare(args.study.resolve()),allow_nan=False),flush=True)
