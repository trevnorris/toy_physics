"""Two fixed-shape R4 equilibrations with the smooth quadrature-point gauge.

Only the second run may start from a certified primary endpoint. Interior
melt-only labels affect identically melted cells, where the smooth weight is
one. Their existing precision solver is consequently the exact restriction
of the full smooth functional, including on a changing order support.
"""
import argparse
from pathlib import Path

import numpy as np

from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_r4_barrier_held import HeldStop
from mouth_r4_eulerian_relax_labels import FORCE_TOL, absmax, chart_orientation
from mouth_r4_eulerian_release import melted_cells, order_grid
from mouth_r4_eulerian_solid import SpatialCrystal
from mouth_r4_projected_fields import EXPECTED, location
from mouth_r4_projected_held import ActiveStepLimit, ProjectedHeld
from mouth_r4_smooth_fields import coherent, read_fields
from mouth_r4_smooth_static import Gauge
from mouth_whole_r4_relax import atomic_json, log


class SmoothHeld(ProjectedHeld):
    def __init__(self, scratch, folder, epsilon):
        super().__init__(scratch, folder, epsilon)
        self.outer_force = FORCE_TOL
        self.point_cache = None

    def extend(self, barrier, labels):
        return super().extend(barrier, labels, outer_force=self.outer_force)

    def point(self, labels, order, solid, barrier):
        if self.point_cache is not None:
            at, on_order, result = self.point_cache
            if at is labels and on_order is order:
                return result
        opt = barrier.opt
        certificate = chart_orientation(opt.axes, labels, opt.active)
        elastic, ge, _, physical = solid.evaluate(labels, order)
        if (physical['min_Jinv_active'] <= np.sqrt(np.finfo(float).eps)
                or certificate['chart_det_active_lower'] <= 128*np.finfo(float).eps):
            raise HeldStop('Positive-order crystal orientation is unresolved', labels,
                           dict(physical=physical, chart=certificate))
        cells, unit_grad = self.full_gauge.evaluate(labels, np.zeros(6), gradient=True)
        unit_energy = float(cells.sum())
        unit_grad = np.asarray(unit_grad, float)
        _, restricted, _, _ = barrier.sample(labels)
        melt_grad = barrier.extract(unit_grad, 'dead')
        restriction_error = absmax(melt_grad - barrier.extract(restricted, 'dead'))
        if restriction_error > 1.e-9 * max(1., absmax(melt_grad)):
            raise ArithmeticError('Restricted melted derivative differs from the pointwise gauge')
        gradient = barrier.extract(ge + self.epsilon*unit_grad, 'free')
        magnitude = absmax(gradient)
        dead_elastic = absmax(barrier.extract(ge,'dead'))
        if dead_elastic > 1.e-12 or not np.all(np.isfinite(gradient)) or not np.isfinite(unit_energy):
            raise ArithmeticError('Unresolved physical dead-label load or nonfinite smooth energy')
        gauge = self.epsilon*unit_energy
        info = dict(**physical, chart=certificate, elastic_energy=elastic,
                    unit_gauge_energy=unit_energy, gauge_energy=gauge, gauge_to_elastic=gauge/elastic,
                    max_free_label_force=magnitude, max_free_elastic_force=absmax(barrier.extract(ge,'free')),
                    max_unit_melt_residual=absmax(melt_grad), max_dead_elastic_force=dead_elastic,
                    max_melt_restriction_error=restriction_error,
                    stationary=magnitude<=FORCE_TOL and absmax(melt_grad)<=3.e-8)
        k = int(np.argmax(np.abs(gradient)))
        for c, cm in enumerate(barrier.components):
            free = np.flatnonzero(cm.free)
            if k < len(free):
                group = free[k]
                node = tuple(np.unravel_index(cm.source[group], order.shape))
                info['worst_free_active'] = dict(component=c, group=int(group), signed_force=float(
                    -cm.sum((ge+self.epsilon*unit_grad)[...,c])[group]), chi=float(order[node]),
                    **location(self.current_trial,opt.axes,node))
                break
            k -= len(free)
        result = elastic+gauge, gradient, info
        self.point_cache = labels, order, result
        return result

    def held(self, t, solid):
        self.current_trial = t
        previous = None
        last_field = {}
        rematches = 0
        try:
            while True:
                self.fields, before_extension, _, last_field = coherent(self,t,solid,self.labels,self.fields)
                order = order_grid(t,self.fields)
                barrier = self.context(t,self.labels,order,melted_cells(order))
                self.full_gauge = Gauge(t,self.reference_depth,order)
                _, before_gradient, before_state = self.point(self.labels,order,solid,barrier)
                total_before = before_extension + before_state['gauge_energy']
                self.outer_force = max(FORCE_TOL,absmax(before_gradient))
                labels,inner = self.extend(barrier,self.labels)
                value,gradient,material = self.point(labels,order,solid,barrier)
                if (material['elastic_energy'] != before_state['elastic_energy']
                        or not np.array_equal(labels[barrier.opt.stencil],self.labels[barrier.opt.stencil])):
                    raise ArithmeticError('Melt extension changed physically loaded crystal')
                self.labels = labels
                physical,load,last_field = read_fields(self,t,solid,labels,self.fields)
                total = physical + material['gauge_energy']
                roundoff = 64*abs(np.spacing(float(total)))
                if (total>total_before+roundoff
                        or (previous is not None and total>previous+roundoff)):
                    raise HeldStop('A complete held sweep increased the smooth total energy',labels,
                                   dict(previous=previous,before_extension=total_before,total=total,material=material))
                if abs(last_field['gauge_energy']-material['gauge_energy'])>1.e-12:
                    raise ArithmeticError('Independent smooth gauge energy disagrees between fields and labels')
                previous = total
                rematches += 1
                details = dict(material=material,field=last_field,jointly_matched_energy=total,
                               outer=self.outer_sweep,field_rematches=rematches,inner=inner,
                               smallest_accepted_active_fraction=self.smallest_alpha,
                               starting_sha256=self.start_sha,case=self.case,fixed_shape=self.x.tolist())
                self.save_internal('matched',t,labels,self.fields,details)
                field_max=max(last_field[key] for key in
                              ('order_kkt','order_projected','carrier_metric','amplitude_raw','phase_raw'))
                log('smooth_rematched',outer=self.outer_sweep,energy=total,active_force=absmax(gradient),
                    field_force=field_max,unit_inner=material['max_unit_melt_residual'],
                    min_Jinv=material['min_Jinv_global'])
                if material['stationary'] and field_max<=FORCE_TOL:
                    return details
                if absmax(gradient)<=FORCE_TOL:
                    # The fresh extension may have changed the order source or
                    # arrived with a loose inner tolerance from its starting force.
                    if absmax(before_gradient)<=FORCE_TOL and field_max<=FORCE_TOL:
                        raise HeldStop('Final inner residual exceeds target after a strict extension',
                                       labels,details)
                    continue

                # Reuse the fixed physical order/phase and full gauge for up
                # to four accepted moves, then re-equilibrate all fields.
                for burst in range(4):
                    self.outer_force = max(FORCE_TOL,absmax(gradient))
                    new_labels,search = self.active_step(t,solid,order,barrier,labels,value,gradient,
                                                         material,diagonal_first=True)
                    next_value,gradient,material = self.point(new_labels,order,solid,barrier)
                    if not (search['alpha']>1.e-4 and next_value<value):
                        raise ArithmeticError('Accepted active step missed the energy or fraction gate')
                    self.labels = labels = new_labels
                    previous += next_value-value
                    value = next_value
                    alpha = search['alpha']
                    self.smallest_alpha = min(alpha, self.smallest_alpha or 1.)
                    self.outer_sweep += 1
                    self.save_internal('last_active',t,labels,self.fields,
                        dict(outer=self.outer_sweep,burst=burst+1,search=search,material=material,
                             smallest_accepted_active_fraction=self.smallest_alpha))
                    log('smooth_active_step',outer=self.outer_sweep,burst=burst+1,alpha=alpha,
                        family=search['family'],energy_change=search['energy_change'],
                        active_force=absmax(gradient),min_Jinv=material['min_Jinv_global'],
                        min_ordered_Jinv=material['min_Jinv_active'])
                    if absmax(gradient)<=FORCE_TOL:
                        break
        except Exception as error:
            failure = error if isinstance(error,HeldStop) else HeldStop(
                type(error).__name__+': '+str(error),self.labels,{})
            if isinstance(error,HeldStop) and not isinstance(error,ActiveStepLimit):
                self.save_internal('failed_trial',t,error.labels,self.fields,error.details)
            self.write_stop(t,solid,failure,last_field)
            raise


def run(study, case, start):
    gate = read(study/'eligibility.json')
    if not gate['passed'] or gate['starting_sha256']!=EXPECTED['held_stop']:
        raise ValueError('Missing frozen-state clearance')
    epsilon = 1e-8 if case=='primary' else 1e-10
    expected = EXPECTED['held_stop']
    if case=='secondary':
        held = read(study/'primary/HELD.json')
        expected = held['saved']
        if start != study/'primary/stationary_state.npz' or held['epsilon']!=1e-8:
            raise ValueError('Expected the signed primary stationary state for the secondary warm start')
    if digest(start)!=expected:
        raise ValueError('Held starting state digest differs')
    folder=study/case
    folder.mkdir()
    try:
        driver=SmoothHeld(study.parent,folder,epsilon)
        driver.start_sha,driver.case=expected,case
        with np.load(start,allow_pickle=False) as data:
            if ((str(data['parent_sha256']),int(data['radius']),int(data['winding']))!=(CHECKPOINT_SHA,4,1)
                    or not np.array_equal(data['x'],driver.x)):
                raise ValueError('Unexpected radius, winding, parent or fixed shape')
            driver.fields,driver.labels=data['field'].copy(),data['labels'].copy()
        t=driver.trial()
        solid=SpatialCrystal(t.b)
        log('smooth_held_start',case=case,epsilon=epsilon,starting_sha256=expected,
            depth=t.shape.length,exit_radius=t.shape.exit)
        details=driver.held(t,solid)
        from mouth_r4_smooth_endpoint import write_endpoint
        write_endpoint(driver,t,solid,details)
    except Exception as error:
        atomic_json(folder/'ISSUE.json',dict(error=type(error).__name__,message=str(error)))
        raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study',type=Path)
    parser.add_argument('case',choices=('primary','secondary'))
    parser.add_argument('start',type=Path)
    args=parser.parse_args()
    run(args.study.resolve(),args.case,args.start.resolve())
