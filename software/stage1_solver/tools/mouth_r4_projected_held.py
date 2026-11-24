"""Relax the R4 active labels at fixed shape with a slaved melted gauge.

Projection changes the search metric only. Every finite proposal uses the
original volume certificate, reduced energy and full active gauge gradient.
The order and carrier are rematched after each accepted active step.
"""
import argparse
from pathlib import Path

import numpy as np

from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_cp10_field_shape import Bath
from mouth_cp10_nested_forces import integral, kernel
from mouth_r4_barrier_held import BarrierHeldDriver, HeldStop, fractions
from mouth_r4_barrier_shape import unit_work
from mouth_r4_corner_projection import packed_lift, point_normal, project_halfspaces
from mouth_r4_eulerian_relax_labels import FORCE_TOL, InfeasibleLabels, absmax, chart_orientation
from mouth_r4_eulerian_release import block_axes, melted_cells, order_grid
from mouth_r4_eulerian_solid import InvalidSpatialCrystal, SpatialCrystal
from mouth_supercore_shape import preparation
from mouth_whole_r4_relax import FIELD_TOL, atomic_json, log

OLD_STATE = '418376c3eb1922f6335543f43dfee11e531cebbb6958d37afd73e844571b653a'
CORNER = ((15,9,1), (1.,1.,0.))


class ActiveStepLimit(HeldStop):
    """Only a healthy-fraction search failure is eligible for stronger-gauge continuation."""


def analytic_forces(driver,t,solid):
    axes=block_axes(t.b)
    order=order_grid(t,driver.fields)
    depth=t.shape.length
    z=axes[1]
    dz=np.where(z*depth<depth-.12,(z*depth/(depth-.12))**25,1.)-z
    energy,gradient,_,_=solid.evaluate(driver.labels,order)
    shape_work=integral(kernel(driver.folder),t.b,axes,dz,driver.labels,order,
                        (solid.points,solid.weights),True)
    if abs(shape_work['energy']-energy)>1.e-12:
        raise ArithmeticError('Independent fixed-label energy differs on the held bank')
    crystal=np.array(shape_work['analytic'],dtype=float)
    crystal[0]-=np.sum(gradient[...,1]*dz[None,:,None],dtype=np.longdouble)
    bath=Bath(t,driver.base,driver.fields,driver.reference_depth)
    for i,unit in enumerate(np.eye(6)):
        dq=1.e-24j*unit
        exterior=bath.energy(dq).sum()
        gauge=unit_work(t.shape,t.polar_shift,driver.reference_depth,axes,driver.labels,melted_cells(order),dq)
        crystal[i]-=(exterior+driver.epsilon*gauge).imag/1.e-24
    return crystal


def local_inverse_diagonal(opt, order, grad):
    """Positive Jacobi metric with a quarter-cell trust scale per independent label."""
    # _precondition constructs the same integrated active-cell elliptic matrix.
    grids = np.moveaxis(opt.base, -1, 0)
    radius = opt.body.sh.radius(grids[1]*opt.body.sh.length)[0]
    lengths = (radius, np.full_like(radius,opt.body.sh.length),
               np.maximum(grids[0]*radius,.5*np.diff(opt.axes[0])[0]*radius))
    cell_chi = np.lib.stride_tricks.sliding_window_view(order,(2,2,2)).mean(axis=(-1,-2,-3))
    weight = np.zeros_like(order)
    count = np.zeros_like(order)
    for offset in np.ndindex(2,2,2):
        ids = tuple(slice(b,b+n) for b,n in zip(offset,opt.active.shape))
        weight[ids] += cell_chi
        count[ids] += opt.active
    mean_order = np.maximum(weight/np.maximum(count,1.),1.e-4)
    diagonal, widths = [], []
    for c, cm in enumerate(opt._components()):
        src = cm.source[cm.free]
        k = (cm.C.T @ opt.Kactive @ cm.C).diagonal()[cm.free]
        diagonal.append(np.maximum(k,1.e-12)*(mean_order*lengths[c]**2).ravel()[src])
        dx = np.diff(opt.axes[c])
        node_width = np.r_[dx[0],np.minimum(dx[1:],dx[:-1]),dx[-1]]
        widths.append(node_width[np.unravel_index(src,order.shape)[c]])
    inv = 1./np.maximum(np.concatenate(diagonal), abs(grad)/(.25*np.concatenate(widths)))
    return lambda v: inv*v


class ProjectedHeld(BarrierHeldDriver):
    def __init__(self,scratch,folder,epsilon):
        super().__init__(scratch,folder,epsilon)
        self.smallest_alpha = None
        self.corners = [CORNER]
        self.reference = preparation(4)[0]

    def point(self,labels,order,solid,barrier):
        value,gradient,info = super().point(labels,order,solid,barrier)
        if self.epsilon > 1.e-8:
            info['stationary'] = (absmax(gradient)<=FORCE_TOL and info['max_unit_melt_residual']<=3.e-8)
        return value,gradient,info

    def active_step(self,t,solid,order,barrier,current,energy,grad,state):
        opt = barrier.opt
        elliptic = opt._precondition(order)
        diagonal = local_inverse_diagonal(opt,order,grad)
        tests = []

        def probe(direction, alpha, step, family):
            proposed = current + alpha*step
            chart_orientation(opt.axes,proposed,opt.active)
            slope = float(grad@direction)
            elastic = solid.evaluate(proposed,order,False)[0]
            if elastic >= energy+1.e-4*alpha*slope:
                raise InfeasibleLabels('Elastic lower bound misses sufficient decrease')
            try:
                labels,inner = self.extend(barrier,proposed)
                value,_,after = self.point(labels,order,solid,barrier)
            except HeldStop as exc:
                if exc.details['reason']=='Changed registry on an active trace':
                    raise
                raise InfeasibleLabels('Unaccepted proposal failed the inner or orientation gate',
                                       rejected_inner=exc.details) from exc
            if not (value<energy and value<=energy+1.e-4*alpha*slope):
                raise InfeasibleLabels('Reduced energy misses sufficient decrease',delta=value-energy)
            if absmax(opt._pack(labels)-(opt._pack(current)+alpha*direction))>5.e-13:
                raise ArithmeticError('Slaved extension changed a proposed active label')
            evidence = dict(family=family,alpha=alpha,slope=slope,energy_change=value-energy,
                            physical_delta=after['elastic_energy']-state['elastic_energy'],
                            gauge_delta=after['gauge_energy']-state['gauge_energy'],
                            min_Jinv=after['min_Jinv_global'],min_ordered_Jinv=after['min_Jinv_active'],
                            whole_cell_bound=after['chart']['chart_det_global_lower'],
                            active_force=after['max_free_label_force'],inner=inner,rejections=tests)
            np.savez_compressed(self.folder/f'direction_{self.outer_sweep:06d}.npz',
                                direction=direction,gradient=grad,active_mask=opt.active)
            return labels,evidence

        for name, metric in [('elliptic',elliptic),('local_diagonal',diagonal)]:
            # One correction per new witness. Repeated witnesses go through the
            # finite line search. A separate maximum-force direction remains
            # available if the local tangents block the coupled proposal.
            for refinement in range(4):
                normals = []
                for cell,local in self.corners:
                    physical = self.reference if (cell,local)==CORNER else None
                    normal,site = point_normal(opt,current,cell,local,reference=physical)
                    # Only near-inverting cells need local one-sided constraints.
                    if site['chart_det']<.10 and np.any(normal):
                        normals.append(normal)
                raw = -metric(grad)
                direction = project_halfspaces(raw,normals,metric) if normals else raw
                if grad@direction >= 0:
                    break
                new_corner = None
                for harmonic in (True,False):
                    step = packed_lift(opt,direction,harmonic=harmonic)
                    for alpha in fractions():
                        if alpha<=1.e-4:
                            break
                        try:
                            return probe(direction,alpha,step,name+('_harmonic' if harmonic else '_fixed_melt'))
                        except (InvalidSpatialCrystal, InfeasibleLabels) as exc:
                            issue = getattr(exc,'witness',getattr(exc,'diagnostics',{}))
                            tests.append(dict(alpha=alpha,family=name,harmonic=harmonic,reason=str(exc),site=issue))
                            frac = issue.get('chart_fractions',issue.get('local'))
                            if 'cell' in issue and frac is not None:
                                key = (tuple(issue['cell']),tuple(frac))
                                if key not in self.corners and opt.active[key[0]]:
                                    _,pt = point_normal(opt,current,*key)
                                    if pt['chart_det']<.10:
                                        new_corner=key
                if new_corner is None:
                    break
                self.corners.append(new_corner)
        # A single far-field load may be geometrically independent of a tight
        # corner. Its natural Jacobi unit is reported, alongside the fraction.
        idx = int(np.argmax(np.abs(grad)))
        d = np.zeros_like(grad)
        d[idx] = -diagonal(grad)[idx]
        if d[idx] != 0.:
            step = packed_lift(opt,d,harmonic=False)
            for alpha in fractions():
                if alpha<=1.e-4:
                    break
                try:
                    return probe(d,alpha,step,'largest_residual_coordinate')
                except (InvalidSpatialCrystal, InfeasibleLabels) as exc:
                    tests.append(dict(alpha=alpha,family='largest_residual_coordinate',reason=str(exc),
                                 site=getattr(exc,'witness',getattr(exc,'diagnostics',{}))))
        raise ActiveStepLimit('Projected, diagonal and isolated active steps failed at healthy fractions',
                              current,dict(material=state,tests=tests,corners=self.corners))

    def held(self,t,solid):
        previous = None
        while True:
            self.fields,physical,load,field = self.coherent(t,solid,self.labels,self.fields)
            order = order_grid(t,self.fields)
            barrier = self.context(t,self.labels,order,melted_cells(order))
            try:
                # Extension cannot alter physical crystal work with respect to
                # active labels, but it can change the order source at χ=0.
                labels,inner = self.extend(barrier,self.labels)
                extension_changed = not np.array_equal(labels,self.labels)
                self.labels = labels
                value,grad,material = self.point(labels,order,solid,barrier)
                total = physical+material['gauge_energy']
                if previous is not None and total>previous+64*abs(np.spacing(previous)):
                    raise ArithmeticError('A coupled held sweep increased the actual energy')
                previous=total
                details=dict(material=material,field=field,jointly_matched_energy=total,outer=self.outer_sweep,
                             smallest_accepted_active_fraction=self.smallest_alpha)
                self.save_internal('matched',t,labels,self.fields,details)
                if material['stationary']:
                    # Recheck fields against the FINAL extended labels. The
                    # finite mask change is accounted for by coherent().
                    if extension_changed:
                        # Reassemble the order loads on the final labels. Check
                        # physical obstacle stationarity even on the χ=0 side.
                        _,_,nodes,_=solid.evaluate(labels,order,False)
                        load=t.ext+np.bincount(t.b.ids.ravel(),weights=nodes.transpose(0,2,1).ravel(),
                                                minlength=t.model.size)
                        t.model.material_load=load
                        physical,gr=t.model.evaluate(self.fields)
                        joined=t.pull.T@gr[2]
                        chi=self.fields[2,t.first]
                        kkt=np.where(chi==0.,np.minimum(joined,0.),np.where(chi==1.,np.maximum(joined,0.),joined))
                        field.update(order_raw=absmax(joined[t.free]),order_kkt=absmax(kkt[t.free]),
                                     order_projected=absmax(t.projected(chi,joined)),
                                     order_knot_transport=self.forces_at_order_knots(t,self.fields,gr))
                        total=float(physical)+material['gauge_energy']
                        details['jointly_matched_energy']=total
                        if max(field['order_kkt'],field['order_projected'])>=FIELD_TOL:
                            continue
                    if max(field['amplitude_raw'],field['phase_raw'])>FORCE_TOL:
                        count,_=t.carrier(self.fields,False,1.e-10)
                        if count==0:
                            raise HeldStop('Tighter carrier solve leaves an unresolved raw force',labels,details)
                        continue
                    details['one_sided_order_work']=self.mask_work(t,barrier,order)
                    return total,load,details
                changed,search = self.active_step(t,solid,order,barrier,labels,value,grad,material)
                self.labels=changed
                alpha=search['alpha']
                self.smallest_alpha=min(alpha,self.smallest_alpha if self.smallest_alpha is not None else 1.)
                self.save_internal('last_active',t,changed,self.fields,dict(outer=self.outer_sweep,search=search,
                                   smallest_accepted_active_fraction=self.smallest_alpha))
                log('projected_active_step',outer=self.outer_sweep, family=search['family'],alpha=alpha,
                    energy_change=search['energy_change'],material_force=search['active_force'],
                    chart_lower=search['whole_cell_bound'],order_kkt_before=field['order_kkt'])
                self.outer_sweep+=1
            except HeldStop as exc:
                # A failed inner TRIAL is not the latest accepted material state.
                if not isinstance(exc,ActiveStepLimit):
                    self.save_internal('failed_trial',t,exc.labels,self.fields,exc.details)
                self.write_stop(t,solid,exc,field)
                raise

    def mask_work(self,t,barrier,order):
        """A fixed-mask KKT can coexist with an admissible negative order variation."""
        _,_,_,det=barrier.sample(self.labels)
        # The nonnegative corner energy selects a witness. Integrate BOTH
        # gauge terms on all disappearing cells of that joined order DOF.
        a=det-1.
        per_cell=barrier.weight*np.maximum(a-np.log1p(a),0.).sum(axis=1)
        groups=t.number[t.b.ids.transpose(0,2,1)].ravel()
        cell_dofs=groups[barrier.ids]
        free=t.free[cell_dofs]
        if not np.any(free):
            return dict(free_melt_site=False)
        i,k=np.unravel_index(np.argmax(np.where(free,per_cell[:,None],-1.)),free.shape)
        ix=cell_dofs[i,k]
        if self.fields[2,t.first[ix]]!=0.:
            raise ArithmeticError('Mask probe node has positive order')
        mask=melted_cells(order)
        selected=mask.copy()
        selected[mask]=np.any(cell_dofs==ix,axis=1)
        work=-self.epsilon*float(unit_work(t.shape,t.polar_shift,self.reference_depth,
            block_axes(t.b),self.labels,selected,np.zeros(6)))
        _,g=t.model.evaluate(self.fields)
        return dict(free_melt_site=True,joined_node=int(ix),cells=int(selected.sum()),
                    physical_right_slope=float((t.pull.T@g[2])[ix]),gauge_right_jump=work,
                    two_sided_shape_held=True)


def run(study,case,start):
    eligibility=read(study/'eligibility.json')
    if not eligibility['passed'] or digest(start)!=eligibility['starting_sha256']:
        raise ValueError('Missing Part 1 or derivative clearance')
    epsilon=1.e-8 if case=='primary' else 1.e-10
    folder=study/case
    folder.mkdir()
    try:
        with np.load(start) as data:
            if (str(data['source_sha256']),int(data['radius']),int(data['winding']),
                    str(data['parent_sha256']))!=(OLD_STATE,4,1,CHECKPOINT_SHA):
                raise ValueError('Wrong held starting state')
            fields,labels,x = (data[key] for key in ('field','labels','x'))
        driver=ProjectedHeld(study.parent,folder,epsilon)
        if not np.array_equal(driver.x,x):
            raise ValueError('Held shape differs')
        driver.fields,driver.labels=fields,labels
        t=driver.trial()
        solid=SpatialCrystal(t.b)
        try:
            energy,load,details=driver.held(t,solid)
        except ActiveStepLimit:
            # Restart from the LAST ACCEPTED fields only. Each gauge gets its
            # own inner solve, metric and energies. No shape coordinate changes.
            chain=[1.e-4,1.e-6,1.e-8]+([1.e-10] if case=='secondary' else [])
            for number,strength in enumerate(chain):
                stage=folder/f'continuation_{number}'
                stage.mkdir()
                next_driver=ProjectedHeld(study.parent,stage,strength)
                next_driver.fields,next_driver.labels=driver.fields,driver.labels
                driver=next_driver
                energy,load,details=driver.held(t,solid)
                atomic_json(stage/'stationary.json',details)
            driver.folder=folder
        if driver.epsilon!=epsilon:
            raise ValueError('Intermediate gauge offered as a final solution')
        analytic=analytic_forces(driver,t,solid)
        force,curvature,base,_=driver.signed_work(t,solid)
        half,_,_,_=driver.signed_work(t,solid,half=True)
        atomic_json(folder/'signed_forces.json',dict(analytic=analytic.tolist(),full=force.tolist(),half=half.tolist(),
                    analytic_minus_half=(analytic-half).tolist(),anchor_energy_difference=base-energy))
        if abs(base-energy)>2.e-8 or max(absmax(force-half),absmax(analytic-half))>1.e-4:
            raise ArithmeticError('Signed held-shape work fails energy or half-step consistency')
        details['half_shape_forces']=half.tolist()
        details['analytic_shape_forces']=analytic.tolist()
        driver.record(t,energy,load,details,force,curvature)
        atomic_json(folder/'HELD.json',dict(epsilon=epsilon,saved=digest(folder/'checkpoint_00010.npz'),
                   force=force.tolist(),initial=digest(start),clearance=digest(study/'eligibility.json')))
    except Exception as exc:
        atomic_json(folder/'ISSUE.json',dict(error=type(exc).__name__,message=str(exc)))
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('study',type=Path)
    p.add_argument('case',choices=('primary','secondary'))
    p.add_argument('start',type=Path)
    a=p.parse_args()
    run(a.study.resolve(),a.case,a.start.resolve())
