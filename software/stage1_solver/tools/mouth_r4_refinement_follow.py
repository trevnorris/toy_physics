"""Read the first conforming R4 step; continue only while the physical descent gate holds."""
from pathlib import Path
import argparse

import numpy as np

from mouth_backtrack_r4_release import BacktrackRelease
from mouth_continuous_shape_relax import ContinuousRelease
from mouth_feasible_crystal_step import feasible_crystal_block, TINY_ALPHA
from mouth_joint_endpoint_audit import geometry
from mouth_melted_collar_survival import Survival
from mouth_r4_refined_bounds import refined_bounds
from mouth_r4_refinement_audit import inputs
from mouth_r4_refinement_descent import RefinedTrial, hard_gate, paired_point, source_gate
from mouth_rim_pinned_comparison import digest, read
from mouth_whole_r4_relax import SHAPE_TOL, atomic_json, log, save_fields


class ScientificStop(Exception):
    """Requested gate reached; record already written."""


class RefinedRelease(BacktrackRelease):
    trial_type=RefinedTrial

    def __init__(self,area,owner,record,material):
        self.__dict__.update(vars(owner))
        self.args.folder=area/'coarse'
        self.args.folder.mkdir(exist_ok=False)
        self.iteration=11
        self.q,self.u=material['q'].copy(),material['field'].copy()
        self.history=list(zip(material['history_step'],material['history_gradient'],material['history_inverse']))
        self.last=self.cp
        self.parent_x=np.array(self.cp['x'])
        self.parent_g=-np.array(self.cp['shape_force'])
        self.previous_force=-self.parent_g.copy()
        self.shape_metric=np.asarray(owner.proposal['metric'])
        self.last_accepted=dict(iteration=10,energy=self.cp['energy'],sha256=self.cp['sha256'])
        self.energy=float(material['energy'])
        self.trial=self.build(self.x,True)
        self.original=owner.stall['material_blocks'][-1]['last_overlap']['witness']
        self.first_report=record['signed_material_search']
        self.locations=[self.original]
        self.purpose='refined shape and ordered crystal'
        self.prior_overlap=None
        self.previous_pin=None
        self.interleaved=False

    def build(self,x,solve_fields=False):
        return ContinuousRelease.build(self,x,solve_fields)

    def bounds(self,t,q,u):
        return refined_bounds(t.b.parent,t.b.unpack(q),q[-3:],u[2],self.locations)

    def stop(self,t,q,u,rep,reason,bounds):
        path=self.args.folder/'material_stop.npz'
        save_fields(path,x=t.x,q=q,surplus=q[-3:],solid=t.b.unpack(q),field=u,
                    energy=rep['energy_after'],state_sha=digest(self.args.state),radius=4,winding=1)
        flow=t.measure(u,t.crystal_load_for_match(q,u)+t.ext)
        row=dict(decision=reason, iteration=self.iteration,last_completed=self.last_accepted,
            stopped_state_sha256=digest(path), search=rep, bounds=bounds,
            physical_check=paired_point(t,q,u[2],bounds['best_J_witness']),
            flow=flow,shape=geometry(t.shape,t.polar_shift),trial_state='fixed-field material update')
        atomic_json(self.args.folder/'scientific_stop.json',row)
        raise ScientificStop(reason)

    def material_update(self,t,q,u):
        self.block_reports=[]
        total=0.
        for i in range(36):
            if self.first_report is not None:
                rep=self.first_report
                self.first_report=None
            else:
                q,rep=feasible_crystal_block(t,q,u,self.history,steps=1,
                    on_event=lambda row:log('refined_crystal_search',checkpoint=self.iteration,**row))
            total+=rep['delta_total_energy']
            self.block_reports.append(rep)
            if rep['last_overlap']:
                self.locations.append(rep['last_overlap']['witness'])
            if rep['status']=='material_tolerance':
                break
            if rep['min_alpha'] is None or rep['min_alpha']<TINY_ALPHA:
                bounds=self.bounds(t,q,u)
                if hard_gate(rep,[dict(primary=bounds['best_J_witness'])]):
                    self.stop(t,q,u,rep,'refined_child_barrier',bounds)
                if not rep['steps'] or -rep['delta_total_energy']<=rep['negligible_drop']:
                    self.stop(t,q,u,rep,'healthy_descent_unresolved',bounds)
        field=self.rematch(t,q,u,'after_refined_crystal',total)
        return q,sum(b['steps'] for b in self.block_reports),rep['status'],field

    def begin_material_record(self,t,q,u,energy,internal,info):
        save_fields(self.args.folder/'pending.npz',x=t.x,q=q,surplus=q[-3:],solid=t.b.unpack(q),field=u,
            energy=energy,state_sha=digest(self.args.state),radius=4,winding=1)

    def complete_material_record(self,t,q,u,energy,internal,info,values):
        cert=t.check_material_domain(q,u)
        forces=np.asarray(values['shape_force'])
        if cert['status']!='certified' or forces.shape!=(6,) or not np.isfinite(forces).all():
            raise ArithmeticError('Incomplete shape forces or ordered-domain evidence')
        if energy>=self.last_accepted['energy']:
            raise ArithmeticError('Checkpoint energy descent unresolved')
        flow=values['current']
        front=flow['shaft_order']['terminal_from_fraction']['disk']['0.0']
        bounds=self.bounds(t,q,u)
        path=self.args.folder/f'state_{self.iteration:05d}.npz'
        if path.exists():
            raise FileExistsError(path)
        (self.args.folder/'pending.npz').rename(path)
        fractions=[b['min_alpha'] for b in self.block_reports if b['min_alpha'] is not None]
        trend=(abs(forces)<abs(self.previous_force)).tolist()
        self.last=dict(**values,**geometry(t.shape,t.polar_shift),x=t.x.tolist(),energy=energy,
            file=path.name,sha256=digest(path),material=internal,orientation=cert,region_bounds=bounds,
            refined_displacement=q[-3:].tolist(),min_alpha=min(fractions,default=None),
            melted_front_depth=front*t.shape.length if front is not None else None,
            shape_magnitudes_decreased=trend,both_exit_and_middle_decrease=trend[1] and trend[2])
        atomic_json(path.with_suffix('.json'),self.last)
        atomic_json(self.args.folder/'latest.json',self.last)
        self.check_progress()
        log('checkpoint',iteration=self.iteration,forces=forces.tolist(),front=self.last['melted_front_depth'])
        self.previous_force=forces
        self.last_accepted=dict(iteration=self.iteration,energy=energy,sha256=self.last['sha256'])

    def finish(self,t,q,u,load,measured):
        if (max(map(abs,self.last['half_shape_force']))>=SHAPE_TOL or
                self.last['field'].get('excluded_support',{}).get('negative_count',0)):
            raise ArithmeticError('Half-step or one-sided endpoint force unresolved')
        Survival.finish(self,t,q,u,load,measured)


def follow(area):
    source_gate(area)
    g=read(area/'guard_step/outcome.json')
    if not(g['limitsVerified'] and g['exitCode']==g['childOutcome']['exitCode']==0
           and g['childOutcome']['guardReason'] is None):
        raise ValueError('First refinement step lacks a verified endpoint')
    step=area/'first_step'
    record=read(step/'result.json')
    if digest(step/'material.npz')!=record['state_sha256']:
        raise ValueError('First crystal step hash changed')
    owner=inputs(area.parent)
    owner.trial_type=RefinedTrial
    with np.load(step/'material.npz',allow_pickle=False) as material:
        t=ContinuousRelease.build(owner,owner.x,False)
        q,u=material['q'],material['field']
        if (not np.array_equal(t.b.unpack(q),material['solid']) or not np.array_equal(material['surplus'],q[-3:])
                or not np.array_equal(owner.u,u) or not np.array_equal(owner.x,material['x'])):
            raise ValueError('Intermediate material state inconsistent')
        locations=[p['primary'] for p in record['physical_points_on_feasible_state']]
        bounds=refined_bounds(t.b.parent,t.b.unpack(q),q[-3:],u[2],locations)
        verified=paired_point(t,q,u[2],bounds['best_J_witness'])
        rep=record['signed_material_search']
        blocked=record['gate']=='refined_child_barrier' or hard_gate(rep,[verified])
        decision='refined_child_barrier' if blocked else 'continue'
        atomic_json(area/'first_step_bounds.json',dict(bounds=bounds,independent=verified,decision=decision))
        del t
        if blocked:
            return decision
        run=RefinedRelease(area,owner,record,material)
    try:
        run.run(matched_energy=run.energy,last_x=run.parent_x,last_g=run.parent_g)
    except ScientificStop as stop:
        return str(stop)
    return 'coarse_equilibrium'


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('area',type=Path)
    area=p.parse_args().area.resolve()
    outcome=follow(area)
    atomic_json(area/'decision.json',dict(decision=outcome))
