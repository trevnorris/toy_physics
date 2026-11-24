"""Advance the held R4 two-metric Eulerian crystal, continuous order and wound carrier.

Each shape has a bonded spatial label field. During shape work, nodal label
offsets from the normalized identity chart stay fixed. This preserves all
registered surface labels when the physical depth knots move. The order and
carrier coefficients follow those knots, with the analytic far trace restored.
"""
import argparse
import gc
import json
from pathlib import Path

import numpy as np

from mouth_breathing_mobile_patch import WELL
from mouth_continuous_depth_mesh import axes_for_depth
from mouth_continuous_order import norm
from mouth_continuous_shape_relax import ShapeTrial
from mouth_cp10_field_shape import Bath
from mouth_cp10_label_gauge import work as gauge_work
from mouth_cp10_spatial_crystal import MODES
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_exchange_carrier import BIAS
from mouth_free_shape_solid import Solid
from mouth_joined_front_audit import Front
from mouth_joint_endpoint_audit import geometry
from mouth_last_rings_fluid import GRIDS
from mouth_r4_eulerian_relax_labels import LabelRelaxation, LabelStall
from mouth_r4_eulerian_solid import InvalidSpatialCrystal, SpatialCrystal
from mouth_r4_fold_order import all_low_after, axes_model
from mouth_sealed_shape_relax import Meridian
from mouth_whole_r4_relax import FIELD_TOL, SHAPE_TOL, atomic_json, log, save_fields


def order_grid(trial, fields):
    return np.ascontiguousarray(fields[2, trial.b.ids].transpose(0, 2, 1))


def block_axes(body):
    return tuple(body.buffers['field_' + s] for s in ('x','t','a'))


def melted_cells(c):
    windows = np.lib.stride_tricks.sliding_window_view(c, (2,2,2))
    return ~np.any(windows > 0., axis=(-3,-2,-1))


def emit(data):
    print(json.dumps(data, allow_nan=False), flush=True)


class Driver:
    def __init__(self, scratch, folder, epsilon):
        self.scratch, self.folder, self.epsilon = scratch, folder, epsilon
        self.prep = scratch/'supercore_joined_front01/adopt309/state.npz'
        endpoint = read(scratch/'r4_whole_mouth01/coarse/endpoint.json')
        self.reference_depth = endpoint['depth']
        self.base = Front(self.prep, GRIDS['double'], 1)
        archive = scratch/'r4_melted_survival01/coarse/state_00010.npz'
        file = scratch/'r4_cp10_labels01/labels.npz'
        self.metric = None
        self.last_x = self.last_gradient = None
        if (digest(archive) != CHECKPOINT_SHA or BIAS != .4510084282857874
                or digest(file) != 'd40b285faa21ce8c6a6e173ce5f8f3e237d78d1df83807cb82b9711ab74146f5'):
            raise ValueError('Frozen checkpoint, labels or reservoir changed')
        with np.load(archive) as cp, np.load(file) as m:
            if int(cp['radius']) != 4 or int(cp['winding']) != 1 or str(cp['state_sha']) != digest(self.prep):
                raise ValueError('Wrong radius, winding or source')
            self.x, self.fields = cp['x'].copy(), cp['field'].copy()
            self.labels = m['labels'].copy()
            self.input_order = m['order'].copy()
        self.iteration = 10

    def shape(self, x):
        return Meridian(4., self.base.sh.params + x[:4]*[1.,.25,.25,.25],
                        extra=self.base.sh.extra+.25*x[4])

    def trial(self):
        shape = self.shape(self.x)
        t = ShapeTrial(self.base, self.prep, self.x, 8,24,6,3,2, precondition=True,
                       axes=axes_for_depth(shape.length, self.reference_depth))
        t.model.extended_energy_sum = True
        t.f = None  # continuous fields and solid are addressed directly
        return t

    def forces_at_order_knots(self, t, u, gradient):
        # Difference in nodal work between advecting the coefficients and
        # resampling their piecewise-linear depth profile on moving knots.
        # Use both one-sided slopes at kinks; sample in each slab's own
        # normalized coordinate, which keeps its two physical endpoints fixed.
        full = t.pull.T @ gradient[2]
        outputs = []
        for sign in (-1,1):
            test = np.zeros(t.model.size)
            covered = np.zeros(t.model.size, bool)
            for (s, side), ids in t.model.ids.items():
                zz = t.model.z[s]
                step = sign*1e-7
                new = axes_for_depth(t.shape.length*np.exp(step),self.reference_depth)[1][s]
                z = (zz-zz[0])/(zz[-1]-zz[0])
                moved = (new-new[0])/(new[-1]-new[0])
                v = (moved-z)/step
                idx = np.clip(np.searchsorted(z,moved,side='right')-1,0,len(z)-2)
                slope = np.diff(u[2,ids],axis=-1)/np.diff(z)
                change = slope[...,idx]*v
                if np.any(covered[ids]) and norm((test[ids]-change)[covered[ids]]) > 1e-7:
                    raise ArithmeticError('Different order transports on the same node')
                test[ids] = change
                covered[ids] = True
            joined = test[t.first]
            if norm(joined[t.number]-test)>1e-7:
                raise ArithmeticError('Split spatial order transport')
            c = u[2,t.first]
            groups = {'free_interior': t.free & (c>0.) & (c<1.),
                      'lower_obstacle': t.free & (c==0.), 'upper_obstacle': t.free & (c==1.)}
            parts = {k: -float(full[sel]@joined[sel]) for k,sel in groups.items()}
            outputs.append(dict(parts, total=sum(parts.values())))
        return outputs

    def coherent(self, t, spatial, labels, fields):
        m = t.model
        u = t.transported(fields)
        u[0,~m.masks[0]] = 1.
        u[1,m.feed], u[1,m.bulk] = 0., BIAS
        old_order = order_grid(t,u)
        old_melt = melted_cells(old_order)
        energy, _, nodes, _ = spatial.evaluate(labels,old_order, gradient=False)
        inner = np.bincount(t.b.ids.ravel(), weights=nodes.transpose(0,2,1).ravel(), minlength=m.size)
        exchange = float(energy - inner@u[2])
        if abs(exchange) > 2e-10*max(1.,abs(energy)):
            raise ArithmeticError(f'Unbalanced order/crystal exchange {exchange}')
        load = t.ext + inner
        joined = t.pull.T @ load
        m.material_load = load
        before = m.evaluate(u)[0]
        while True:
            os, _ = t.solve_order(u, joined, False, 2e-9)
            cs, carrier = t.carrier(u, False, 8e-8)
            total, g = m.evaluate(u)
            c = u[2,t.first]
            gc_order = t.pull.T @ g[2]
            projected = norm(t.projected(c,gc_order))
            if max(carrier,projected) < FIELD_TOL:
                break
            if os == 0 and cs == 0:
                raise ArithmeticError('Coherent fields cannot reduce their residual')
        kkt = np.where(c==0,np.minimum(gc_order,0),np.where(c==1,np.maximum(gc_order,0),gc_order))
        info = dict(phase_raw=norm(g[1,m.masks[1]]), amplitude_raw=norm(g[0,m.masks[0]]),
                    order_raw=norm(gc_order[t.free]), order_kkt=norm(kkt[t.free]),
                    order_projected=projected, carrier_metric=carrier, shared_work_error=exchange,
                    order_knot_transport=self.forces_at_order_knots(t,u,g))
        new_melt = melted_cells(order_grid(t,u))
        def gauge(mask):
            return float(self.epsilon*gauge_work(t.shape,t.polar_shift,self.reference_depth,
                         block_axes(t.b),labels,mask,np.zeros(6)))
        old_gauge = gauge(old_melt)
        new_gauge = gauge(new_melt) if np.any(old_melt!=new_melt) else old_gauge
        change = (float(total)-before)+(new_gauge-old_gauge)
        info['order_step_energy'] = dict(physical_change=total-before,old_gauge=old_gauge,
                    new_gauge=new_gauge,delta=change,entered_melt=int((new_melt & ~old_melt).sum()),
                    left_melt=int((old_melt & ~new_melt).sum()))
        if change>64*np.spacing(before):
            raise ArithmeticError(f'Order and carrier rematch raised same-objective energy: {info["order_step_energy"]}')
        return u, float(total), load, info

    def save_internal(self, name, t, labels, fields, details):
        path = self.folder/(name+'.npz')
        save_fields(path, labels=labels, field=fields, x=self.x, epsilon=self.epsilon,
                    iteration=self.iteration, radius=4, winding=1, parent_sha256=CHECKPOINT_SHA)
        atomic_json(self.folder/(name+'.json'), dict(file=path.name, sha256=digest(path),
                    iteration=self.iteration, **geometry(t.shape,t.polar_shift), **details))

    def settle(self, t, spatial):
        previous = None
        outer = 0
        while True:
            self.fields, e0, first_load, first_fields = self.coherent(t,spatial,self.labels,self.fields)
            c = order_grid(t,self.fields)
            opt = LabelRelaxation(t.b,self.labels,c,self.epsilon)
            try:
                self.labels, search = opt.solve(spatial,self.labels,c,emit)
            except LabelStall as exc:
                # Fields are held during the failed label step. Evaluate currents
                # on those saved fields and identify them as fixed-field diagnostics.
                cur=t.measure(self.fields,first_load)
                self.save_internal('label_stall',t,exc.labels,self.fields,
                    dict(label=exc.diagnostics,field_before_labels=first_fields,held_field_current=cur,
                         frozen_field_front=all_low_after(t.model.z[1],c.transpose(0,2,1),0.)))
                raise
            self.fields, energy, load, field = self.coherent(t,spatial,self.labels,self.fields)
            c = order_grid(t,self.fields)
            opt = LabelRelaxation(t.b,self.labels,c,self.epsilon)
            material = opt.measure(spatial,self.labels,c)
            energy += material['gauge_energy']
            details = dict(field=field,material=material,label_search=search,
                           jointly_matched_energy=energy,outer=outer)
            self.save_internal('internal',t,self.labels,self.fields,details)
            log('coupled_sweep',iteration=self.iteration,outer=outer,energy=energy,
                stationary=material['stationary'],knot_work=field['order_knot_transport'])
            if previous is not None and energy > previous + 64*np.spacing(previous):
                raise ArithmeticError(f'Joint field/label sweep raised energy by {energy-previous}')
            if material['stationary']:
                return energy, load, details
            if previous is not None and previous == energy and e0 == energy-material['gauge_energy']:
                raise ArithmeticError('Joint sweep has no resolved descent at nonstationary labels')
            previous = energy
            outer += 1

    def parts_at(self,t,spatial):
        # Evaluate the common scalar objective at the six actual signed shapes.
        # The carrier and continuous order on both sides use the matched bath.
        # The lighter geometry constructor makes no field solve on a signed probe.
        bath = Bath(t,self.base,self.fields,self.reference_depth)
        c = order_grid(t,self.fields)
        old = block_axes(t.b)
        melt = melted_cells(c)
        def parts(dx):
            shape = self.shape(self.x+dx)
            angle = t.polar_shift+dx[5]
            b = t.b if not np.any(dx) else Solid(self.prep,shape,angle,8,24,6,
                        axes_model(shape,self.reference_depth),degree=3,subdivisions=2)
            labels = self.labels.copy()
            labels[...,1] += (b.buffers['field_t']-old[1])[None,:,None]
            solid = spatial if b is t.b else SpatialCrystal(b)
            stored = solid.evaluate(labels,c,gradient=False)[0]
            gauge = self.epsilon*gauge_work(t.shape,t.polar_shift,self.reference_depth,old,labels,melt,dx)
            if not np.isfinite(gauge):
                raise ArithmeticError('Nonfinite gauge work')
            return np.r_[bath.energy(dx),stored,gauge].astype(np.longdouble)
        return parts

    def signed_work(self, t, spatial, half=False):
        parts = self.parts_at(t,spatial)
        zero = parts(np.zeros(6))
        force, curvature = [], []
        for i,unit in enumerate(np.eye(6)):
            h = (1e-5 if i == 0 else 1e-4)*(.5 if half else 1.)
            minus, plus = parts(-h*unit), parts(h*unit)
            force.append(float(np.sum(minus-plus)/(2*h)))
            curvature.append(float(np.sum((minus-zero)+(plus-zero))/h**2))
        return np.array(force), np.array(curvature), float(zero.sum()), parts

    def record(self,t,energy,load,details,force,curvature):
        m = t.model
        current = t.measure(self.fields,load)
        c = order_grid(t,self.fields)
        front = all_low_after(m.z[1],c.transpose(0,2,1),0.)
        info = dict(**details, shape_force=force.tolist(),shape_curvature=curvature.tolist(),
                    sampled_melted_front=front,current=current,epsilon=self.epsilon,
                    order_max_by_depth=c.max(axis=(0,2)).tolist(),depth_knots=m.z[1].tolist(),
                    energy=energy,phase_drop=BIAS,input_checkpoint_sha256=CHECKPOINT_SHA,
                    shape_alpha=getattr(self,'shape_alpha',None))
        self.save_internal(f'checkpoint_{self.iteration:05d}',t,self.labels,self.fields,info)
        log('checkpoint',iteration=self.iteration,epsilon=self.epsilon,force=force.tolist(),
            exit=t.shape.exit,depth=t.shape.length,current=current['feed_current'],front=front)
        if (current['port_error'] > 3e-5 or abs(current['disk_minus_feed_current'])>3e-5):
            raise ArithmeticError('Conservative inlet, disk and bulk currents disagree')
        return info

    def step(self,force,curvature,parts,energy,t):
        gradient = -force
        fallback = np.diag(1./np.maximum(abs(curvature),1.))
        if self.metric is None:
            self.metric = fallback
        elif self.last_x is not None:
            dx, dg = self.x-self.last_x, gradient-self.last_gradient
            product = dx@dg
            if product>1e-10*np.linalg.norm(dx)*np.linalg.norm(dg):
                b = np.eye(6)-np.outer(dx,dg)/product
                self.metric = b@self.metric@b.T+np.outer(dx,dx)/product
        for metric in (self.metric,fallback):
            direction = metric@force
            direction *= min(1., .035/max(1e-20,norm(direction)))
            slope = gradient@direction
            if slope>=0:
                continue
            alpha = 1.
            while alpha >= 2.**-25:
                delta = alpha*direction
                # Only geometrically invalid trials are backtracked here.
                # Field, quadrature or work-accounting errors must surface.
                try:
                    trial_energy = float(parts(delta).sum())
                except InvalidSpatialCrystal:
                    trial_energy = np.inf
                if trial_energy < energy + 1e-4*alpha*slope:
                    self.last_x, self.last_gradient = self.x.copy(), gradient.copy()
                    self.x = self.x + delta
                    self.metric = metric
                    self.shape_alpha = alpha
                    old = t.b.buffers['field_t']
                    new_shape = self.shape(self.x)
                    new = axes_for_depth(new_shape.length,self.reference_depth)[1][1]/new_shape.length
                    self.labels[...,1] += (new-old)[None,:,None]
                    return trial_energy
                alpha *= .5
        raise ArithmeticError(f'No resolved open-shape descent with forces {force.tolist()}')

    def run(self, mode, prior=None):
        if mode == 'finish':
            if prior is None:
                raise ValueError('Verified checkpoint 11 required')
            gate = read(prior.parent/'gauge_11.json')
            state = prior/'checkpoint_00011.npz'
            if not gate['passed'] or digest(state)!=gate['states'][str(self.epsilon)]:
                raise ValueError('Unverified checkpoint 11 restart')
            rec = read(prior/'checkpoint_00011.json')
            with np.load(state) as a:
                self.labels,self.fields,self.x = (a[k].copy() for k in ('labels','field','x'))
            history=prior/'secant_00011.npz'
            with np.load(history) as a:
                if digest(history)!=rec['secant_sha256']:
                    raise ValueError('Changed shape secant')
                self.last_x,self.last_gradient,self.metric=(a[k].copy() for k in ('last_x','last_gradient','metric'))
            self.iteration=11
            t=self.trial()
            spatial=SpatialCrystal(t.b)
            parts=self.parts_at(t,spatial)
            energy_bound=self.step(np.array(rec['shape_force']),np.array(rec['shape_curvature']),
                                   parts,rec['energy'],t)
            self.iteration=12
            del t,spatial,parts
            gc.collect()
        else:
            log('start',checkpoint=CHECKPOINT_SHA,epsilon=self.epsilon,coordinates=MODES,
                radius=4,winding=1,phase_drop=BIAS,material_mesh=[24,24],label_nodes=list(self.labels.shape[:3]))
            energy_bound = None
        while True:
            t = self.trial()
            spatial = SpatialCrystal(t.b)
            if self.iteration == 10 and not np.array_equal(order_grid(t,self.fields),self.input_order):
                raise ValueError('Order differs from checkpoint 10')
            energy,load,details = self.settle(t,spatial)
            if energy_bound is not None and energy>energy_bound+64*np.spacing(energy):
                raise ArithmeticError('Coupled rematch increased accepted shape energy')
            force,curvature,base,parts = self.signed_work(t,spatial)
            if abs(base-energy)>2e-8:
                raise ArithmeticError(f'Shape-work energy differs from coupled energy by {base-energy}')
            info=self.record(t,energy,load,details,force,curvature)
            front=info['sampled_melted_front']
            if front == 0. or t.shape.exit <= 0.:
                intact = [v for v in details['material']['intact_by_depth'] if v is not None]
                closest = next((v for v in reversed(intact)
                                if front is None or v['physical_depth']<=front),None)
                atomic_json(self.folder/'physical_stop.json', dict(iteration=self.iteration,
                        front=front,exit=t.shape.exit,well_threshold=WELL/3,
                        last_intact_band_before_front=closest,material=details['material']))
                return
            if self.iteration>=11 and max(abs(z['total']) for z in details['field']['order_knot_transport'])>=1e-4:
                raise ArithmeticError('Relaxed-state one-sided order transport exceeds envelope tolerance')
            if norm(force)<=SHAPE_TOL:
                half,_,_,_=self.signed_work(t,spatial,half=True)
                if norm(half-force)>1e-4 or norm(half)>SHAPE_TOL:
                    raise ArithmeticError('Coarse stationarity fails half-step shape-force check')
                cur = info['current']
                if cur['feed_current']<=0 or not (cur['dense_total_subsonic_share'] and cur['dense_total_subsonic_share']>0):
                    raise ArithmeticError('Stationary geometry lacks verified inward subsonic drainage')
                atomic_json(self.folder/'equilibrium.json',dict(iteration=self.iteration,
                            force=force.tolist(),half_step_force=half.tolist(),energy=energy))
                return
            if mode=='cp11' and self.iteration==11:
                path=self.folder/'secant_00011.npz'
                save_fields(path,metric=self.metric,last_x=self.last_x,last_gradient=self.last_gradient)
                record=self.folder/'checkpoint_00011.json'
                atomic_json(record,dict(read(record),secant_sha256=digest(path)))
                return
            energy_bound=self.step(force,curvature,parts,base,t)
            self.iteration+=1
            del parts,spatial,t
            gc.collect()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('scratch',type=Path)
    p.add_argument('folder',type=Path)
    p.add_argument('epsilon',type=float,choices=(1e-8,1e-10))
    p.add_argument('mode',choices=('cp11','finish'))
    p.add_argument('--prior',type=Path)
    a=p.parse_args()
    a.folder.mkdir(exist_ok=False,parents=True)
    try:
        Driver(a.scratch.resolve(),a.folder.resolve(),a.epsilon).run(a.mode,a.prior)
    except Exception as exc:
        atomic_json(a.folder/'issue.json',dict(error=type(exc).__name__,message=str(exc)))
        raise


if __name__=='__main__':
    main()
