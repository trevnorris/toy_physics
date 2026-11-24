"""Held R4 coupled equilibration with the ordered trace and a slaved melt-only barrier.

Each accepted active-label trial lowers elastic plus minimized barrier energy
at fixed physical fields. The carrier/order rematch accounts separately for
the gauge energy of cells entering and leaving exact melt. A paired external
gate is required before any later release of the six shape coordinates.
"""
import argparse
from pathlib import Path

import numpy as np

from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_exchange_carrier import BIAS
from mouth_r4_barrier_shape import BarrierShapeDriver
from mouth_r4_eulerian_relax_labels import FORCE_TOL, InfeasibleLabels, LabelRelaxation, absmax, chart_orientation
from mouth_r4_eulerian_release import melted_cells, order_grid
from mouth_r4_eulerian_solid import InvalidSpatialCrystal, SpatialCrystal
from mouth_r4_fold_order import all_low_after
from mouth_r4_melt_barrier import MeltBarrier
from mouth_whole_r4_relax import FIELD_TOL, atomic_json, log


class HeldStop(RuntimeError):
    def __init__(self, reason, labels, details):
        super().__init__(reason)
        self.labels=labels.copy()
        self.details=dict(reason=reason,**details)


def initial_metric_direction(gradient,metric,history):
    """Apply an inverse L-BFGS operator using the elastic elliptic metric."""
    direction=gradient.copy()
    work=[]
    for step,change,inv in reversed(history):
        dot=inv*(step@direction)
        direction-=dot*change
        work.append(dot)
    direction=metric(direction)
    for (step,change,inv),dot in zip(history,reversed(work)):
        direction+=step*(dot-inv*(change@direction))
    return -direction


def fractions():
    value=1.
    while value>1.e-4:
        yield value
        value*=.5
    yield 1.e-4


class BarrierHeldDriver(BarrierShapeDriver):
    def __init__(self,scratch,folder,epsilon):
        super().__init__(scratch,folder,epsilon)
        self.contexts={}
        self.gauge_number=0
        self.operation=0
        self.outer_sweep=0

    def context(self,t,labels,order,mask):
        key=(t.b,np.packbits(mask).tobytes())
        if key not in self.contexts:
            # A rematch needs only the entering and leaving masks. Old geometry
            # contexts otherwise retain their factorizations and dense buffers.
            if len(self.contexts)>=2:
                self.contexts.pop(next(iter(self.contexts)))
            directory=self.folder/f'melt_mask_{self.gauge_number:05d}'
            self.gauge_number+=1
            directory.mkdir()
            opt=LabelRelaxation(t.b,labels,order,self.epsilon)
            if not np.array_equal(mask,~opt.active):raise ValueError('Melt gauge and physical order disagree')
            self.contexts[key]=MeltBarrier(t.b,opt,directory)
        return self.contexts[key]

    def order_mask_gauge(self,t,labels,order,mask):
        barrier=self.context(t,labels,order,mask)
        return self.epsilon*float(barrier.sample(labels)[0])

    def coherent(self,t,solid,labels,fields):
        """Rematch the physical fields and measure the actual finite mask-exchange work."""
        model=t.model
        u=t.transported(fields)
        u[0,~model.masks[0]]=1.
        u[1,model.feed],u[1,model.bulk]=0.,BIAS
        old_order=order_grid(t,u)
        old_mask=melted_cells(old_order)
        old_gauge=self.order_mask_gauge(t,labels,old_order,old_mask)
        elastic,_,nodes,_=solid.evaluate(labels,old_order,False)
        distributed=np.bincount(t.b.ids.ravel(),weights=nodes.transpose(0,2,1).ravel(),minlength=model.size)
        exchange=elastic-distributed@u[2]
        if abs(exchange)>2e-10*max(1.,abs(elastic)):
            raise ArithmeticError(f'Order deposition differs from crystal energy by {exchange}')
        load=model.material_load=t.ext+distributed
        joined=t.pull.T@load
        before=model.evaluate(u)[0]
        while True:
            os,_=t.solve_order(u,joined,False,2e-9)
            cs,carrier=t.carrier(u,False,8e-8)
            physical,g=model.evaluate(u)
            chi=u[2,t.first]
            grad=t.pull.T@g[2]
            kkt=np.where(chi==0,np.minimum(grad,0),np.where(chi==1,np.maximum(grad,0),grad))
            projected=absmax(t.projected(chi,grad))
            if max(absmax(kkt[t.free]),projected,carrier)<FIELD_TOL:break
            if os==cs==0:raise ArithmeticError('Carrier and bounded order cannot reduce their residual')
        order=order_grid(t,u)
        mask=melted_cells(order)
        gauge=self.order_mask_gauge(t,labels,order,mask) if not np.array_equal(mask,old_mask) else old_gauge
        change=physical-before+gauge-old_gauge
        info=dict(phase_raw=absmax(g[1,model.masks[1]]),amplitude_raw=absmax(g[0,model.masks[0]]),
            order_raw=absmax(grad[t.free]),order_kkt=absmax(kkt[t.free]),order_projected=projected,
            carrier_metric=carrier,shared_work_error=float(exchange),
            order_knot_transport=self.forces_at_order_knots(t,u,g),
            order_step_energy=dict(physical_change=physical-before,old_gauge=old_gauge,new_gauge=gauge,
               delta=change,entered_melt=int((mask&~old_mask).sum()),left_melt=int((old_mask&~mask).sum())))
        if change>64*np.spacing(before):
            raise ArithmeticError(f'Order/carrier rematch raised full energy: {info["order_step_energy"]}')
        return u,float(physical),load,info

    def extend(self,barrier,labels):
        opt=barrier.opt
        corrected=labels.copy()
        # A changed order mask can register a formerly independent dead-axis
        # unknown with its neighbor; preserve each active/held source exactly.
        for component,cm in enumerate(barrier.components):
            old=(labels-opt.base)[...,component].ravel()
            values=old[cm.representative].copy()
            values[cm.free]=old[cm.source[cm.free]]
            values[cm.held]=0.
            dead=~opt.stencil
            physical=opt.base[...,component]+values[cm.groups].reshape(opt.stencil.shape)
            if absmax((physical-labels[...,component])[~dead])>2.e-13:
                raise HeldStop('Changed registry on an active trace',labels,dict(component=component))
            corrected[...,component][dead]=physical[dead]
        chart_orientation(opt.axes,corrected,opt.active)
        op=f'inner_{self.operation:07d}'
        self.operation+=1
        result,inner=barrier.solve(corrected,op,self.folder)
        if not inner['stationary']:
            raise HeldStop('Melted inner solve failed stationarity',result,dict(inner=inner))
        if not np.array_equal(result[opt.stencil],labels[opt.stencil]):
            raise ArithmeticError('Melt-only solution modified the active trace')
        return result,inner

    def point(self,labels,order,solid,barrier):
        opt=barrier.opt
        chart=chart_orientation(opt.axes,labels,opt.active)
        elastic,ge,_,physical=solid.evaluate(labels,order)
        unit,gg,_,_=barrier.sample(labels)
        grad=barrier.extract(ge+self.epsilon*gg,'free')
        active=absmax(barrier.extract(ge,'free'))
        melted=absmax(barrier.extract(ge,'dead'))
        inner_res=absmax(barrier.extract(gg,'dead'))
        if melted>1e-12:raise ArithmeticError('Physical elastic load on a strictly melted label')
        if not np.isfinite(unit) or not np.isfinite(grad).all():raise ArithmeticError('Nonfinite gauge or active force')
        if physical['min_Jinv_active']<=np.sqrt(np.finfo(float).eps):
            raise HeldStop('Ordered physical inverse volume approached numerical zero',labels,dict(physical=physical))
        # Report a material point directly if the whole-cell certificate can no
        # longer resolve active orientations on the scale of working precision.
        if chart['chart_det_active_lower']<=128*np.finfo(float).eps:
            raise HeldStop('Active whole-cell orientation bound approached roundoff',labels,dict(chart=chart,physical=physical))
        gauge=self.epsilon*float(unit)
        info=dict(**physical,chart=chart,elastic_energy=elastic,gauge_energy=gauge,
            gauge_to_elastic=gauge/elastic,max_free_label_force=absmax(grad),
            max_free_elastic_force=active,max_unit_melt_residual=inner_res,max_dead_elastic_force=melted,
            stationary=max(active,absmax(grad))<=FORCE_TOL and inner_res<=3.e-8)
        return elastic+gauge,grad,info

    def relax_active(self,t,solid,order,barrier,labels):
        opt=barrier.opt
        try:
            current,inner=self.extend(barrier,labels)
        except HeldStop as exc:
            self.save_internal('held_unaccepted_extension',t,exc.labels,self.fields,exc.details)
            raise HeldStop('Initial extension did not clear stationarity',labels,
                dict(probe='held_unaccepted_extension',inner=exc.details)) from exc
        x=opt._pack(current)
        energy,grad,state=self.point(current,order,solid,barrier)
        if state['stationary']:
            return current,dict(accepted_steps=0,material=state)
        metric=opt._precondition(order)
        history=[]; accepted=0; min_alpha=1.;last_alpha=None
        while not state['stationary']:
            rejections=[];ok=False
            for restart in (False,True):
                if restart:
                    if not history:break
                    history.clear()
                d=initial_metric_direction(grad,metric,history)
                slope=grad@d
                if slope>=0:continue
                harmonic=opt._lift(x)
                ok=False
                for alpha in fractions():
                    trial_x=x+alpha*d
                    moving=opt._lift(trial_x)
                    for predictor in ('harmonic_increment','fixed_melt'):
                        proposed=current+(moving-harmonic)
                        if predictor=='fixed_melt':
                            # Keep each dead group fixed while moving active and tied groups.
                            for c,cm in enumerate(barrier.components):
                                dead=cm.dead[cm.groups].reshape(opt.stencil.shape)
                                proposed[...,c][dead]=current[...,c][dead]
                        try:
                            chart_orientation(opt.axes,proposed,opt.active)
                            if not np.allclose(opt._pack(proposed),trial_x,rtol=0.,atol=4e-15):
                                raise ArithmeticError('Active trace changed by an extension predictor')
                            # Unit barrier is nonnegative for any extension. A
                            # physically excessive elastic energy cannot meet Armijo.
                            ep=solid.evaluate(proposed,order,False)[0]
                            if ep>=energy+1.e-4*alpha*slope:
                                rejections.append(dict(alpha=alpha,predictor=predictor,reason='elastic lower bound'))
                                break
                            probe,probe_inner=self.extend(barrier,proposed)
                            new_energy,new_grad,new_state=self.point(probe,order,solid,barrier)
                        except HeldStop as exc:
                            self.save_internal('held_unaccepted_probe',t,exc.labels,self.fields,
                                dict(alpha=alpha,predictor=predictor,probe=exc.details))
                            raise HeldStop('Unaccepted active trial hit a material or gauge stop',current,
                                dict(accepted_steps=accepted,min_alpha=min_alpha if accepted else None,
                                last_accepted_alpha=last_alpha,material=state,trial_alpha=alpha,
                                rejected_file='held_unaccepted_probe.npz',probe=exc.details)) from exc
                        except (InfeasibleLabels,InvalidSpatialCrystal) as exc:
                            rejections.append(dict(alpha=alpha,predictor=predictor,reason=str(exc),
                                witness=getattr(exc,'witness',getattr(exc,'diagnostics',{}))))
                            continue
                        if new_energy<energy and new_energy<=energy+1.e-4*alpha*slope:
                            ok=True
                        else:
                            rejections.append(dict(alpha=alpha,predictor=predictor,reason='reduced Armijo',
                                delta=new_energy-energy))
                        break
                    if ok:break
                if ok:break
            else:
                ok=False
            if not ok:
                raise HeldStop('No certified reducing active step at or above 0.0001',current,
                    dict(accepted_steps=accepted,last_accepted_alpha=last_alpha,min_alpha=min_alpha if accepted else None,
                    material=state,probes=rejections))
            step=opt._pack(probe)-x;change=new_grad-grad;sy=step@change
            if sy>1e-12*np.linalg.norm(step)*np.linalg.norm(change):
                history.append((step,change,1./sy));history=history[-9:]
            accepted+=1;min_alpha=min(min_alpha,alpha);last_alpha=alpha
            change_E=new_energy-energy
            x=opt._pack(probe);current=probe;energy=new_energy;grad=new_grad;state=new_state
            log('barrier_active_step',sweep=self.outer_sweep,step=accepted,alpha=alpha,
                material_force=state['max_free_label_force'],energy_change=change_E,
                inner_residual=probe_inner['residual'],chart_lower=state['chart']['chart_det_global_lower'],
                physical_min_active_Jinv=state['min_Jinv_active'])
            self.save_internal('last_active',t,current,self.fields,dict(outer=self.outer_sweep,
                material=state,active_iteration=accepted,alpha=alpha,min_alpha=min_alpha))
        return current,dict(accepted_steps=accepted,min_alpha=min_alpha,last_accepted_alpha=last_alpha,material=state)

    def settle(self,t,solid):
        previous=None; self.outer_sweep=0
        while True:
            self.fields,_,_,before=self.coherent(t,solid,self.labels,self.fields)
            order=order_grid(t,self.fields)
            barrier=self.context(t,self.labels,order,melted_cells(order))
            try:
                self.labels,search=self.relax_active(t,solid,order,barrier,self.labels)
            except HeldStop as exc:
                self.labels=exc.labels
                self.write_stop(t,solid,exc,before)
                raise
            self.fields,physical,load,field=self.coherent(t,solid,self.labels,self.fields)
            order=order_grid(t,self.fields)
            updated=self.context(t,self.labels,order,melted_cells(order))
            try:
                _,_,material=self.point(self.labels,order,solid,updated)
            except HeldStop as exc:
                self.write_stop(t,solid,exc,field)
                raise
            total=physical+material['gauge_energy']
            details=dict(material=material,field=field,label_search=search,
                jointly_matched_energy=total,outer=self.outer_sweep)
            self.save_internal('internal',t,self.labels,self.fields,details)
            if previous is not None and total>previous+64*np.spacing(previous):
                raise ArithmeticError('Joint field/label rematch raised the reduced objective')
            if (material['stationary'] and field['carrier_metric']<FIELD_TOL
                    and field['order_projected']<FIELD_TOL and field['order_kkt']<FIELD_TOL):
                return total,load,details
            if previous is not None and previous==total and search['accepted_steps']==0:
                stop=HeldStop('Joint rematch stopped changing at unresolved physical KKT',self.labels,details)
                self.write_stop(t,solid,stop,field)
                raise stop
            previous=total;self.outer_sweep+=1

    def write_stop(self,t,solid,exc,field):
        order=order_grid(t,self.fields)
        energy,_,nodes,physical=solid.evaluate(self.labels,order)
        load=t.ext+np.bincount(t.b.ids.ravel(),weights=nodes.transpose(0,2,1).ravel(),minlength=t.model.size)
        self.save_internal('held_stop',t,self.labels,self.fields,dict(stop=exc.details,field=field,
            physical=physical,elastic_energy=energy,current=t.measure(self.fields,load),
            front=all_low_after(t.model.z[1],order.transpose(0,2,1),0.)))


def run(scratch,case,folder,gate):
    if not read(gate)['passed']:raise ValueError('Missing paired reduced-energy audit clearance')
    barrier_data=scratch/'r4_slaved_melt_barrier01'
    measured=read(barrier_data/'measured.json')['cases'][case]
    if not measured['inner']['stationary'] or not measured['outer_probe']['inner']['stationary']:
        raise ValueError('Fixed-trace or active probe unfinished')
    if measured['outer_probe']['alpha']<1e-4 or measured['outer_probe']['reduced_energy_change']>=0:
        raise ValueError('Active descent probe did not pass')
    parent=scratch/'r4_eulerian_stage2_01'/case/'label_stall.npz'
    init=barrier_data/case/'fixed_active.npz'
    if digest(parent)!=measured['sha256'] or digest(init)!=measured['fixed_trace_state_sha256']:
        raise ValueError('Changed held starting state')
    folder.mkdir(parents=True,exist_ok=False)
    d=BarrierHeldDriver(scratch,folder,measured['epsilon'])
    with np.load(parent) as p,np.load(init) as s:
        if str(p['parent_sha256'])!=CHECKPOINT_SHA or int(p['radius'])!=4 or not np.array_equal(p['x'],d.x):
            raise ValueError('Wrong held geometry or checkpoint')
        d.labels=s['labels'];d.fields=p['field']
    t=d.trial();solid=SpatialCrystal(t.b)
    try:
        energy,load,details=d.settle(t,solid)
        force,curvature,base,_=d.signed_work(t,solid)
        half,_,_,_=d.signed_work(t,solid,half=True)
        if abs(base-energy)>2e-8 or absmax(force-half)>1e-4:
            raise ArithmeticError(f'Held signed work disagreement {base-energy}, {absmax(force-half)}')
        details['half_shape_forces']=half.tolist()
        info=d.record(t,energy,load,details,force,curvature)
        atomic_json(folder/'HELD.json',dict(outcome='held_stationary',epsilon=d.epsilon,iteration=d.iteration,
            gate_sha256=digest(gate),input_sha256=digest(parent),saved=digest(folder/'checkpoint_00010.npz'),
            signed_shape_forces=info['shape_force'],energy=energy))
    except Exception as exc:
        atomic_json(folder/'ISSUE.json',dict(error=type(exc).__name__,message=str(exc)))
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('scratch',type=Path);p.add_argument('case',choices=('primary','secondary'))
    p.add_argument('folder',type=Path);p.add_argument('gate',type=Path)
    a=p.parse_args();run(a.scratch.resolve(),a.case,a.folder.resolve(),a.gate.resolve())
