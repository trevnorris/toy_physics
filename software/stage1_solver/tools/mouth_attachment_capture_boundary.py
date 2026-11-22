"""Integrate early joint capture over clearance for fixed local incident rays.

At the positive pole, reversing the drive gives the other pole's launch
condition exactly. The detector is positive time spent with 0<=z<1.5 and
negative tag energy before first exit from the fixed mouth-depth gate.
Both endpoints and the whole clearance range are tested, allowing multiple
capture intervals. A finite grid can miss an interval between equal labels.
"""
import argparse
from pathlib import Path
import json
import subprocess

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import ndtri
from scipy.stats import qmc

from mouth_attachment_impact_trial import sha
from mouth_attachment_incidence import HalfOrderBand, thermal_at_fixed_axial
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import TEMPERATURE, emit


DRIVES = (0., .002, -.002, .001, -.001)
BRACKET_WIDTH = 1e-8


class Trials:
    def __init__(self, args, stream):
        self.band = HalfOrderBand.load(args.section)
        self.probes = qmc.Sobol(4,scramble=True,seed=args.seed).random_base2(args.power)
        self.stream = stream
        self.drives = (0.,.002,-.002) if args.full_strength_only else DRIVES
        self.control = args.control
        self.command = [str(args.engine.resolve()), '--side',str(args.side),'--pitch',str(self.band.pitch),
                        '--dt',str(args.dt),'--steps',str(round(8/args.dt)),'--host','mobile',
                        '--mixed','1','--gate-slope',str(self.band.wall_slope),'--input','-']
        self.cache = {}
        self.min_energy = np.inf
        self.max_defect = 0.
        emit(test='boundary_setup', seed=args.seed, probe_count=len(self.probes),
             offset=1.5, depth_gate=.01, native_time=8., dt=args.dt, cosine=1.,
             source_measure='uniform registry, isotropic yaw, Gaussian radial kick at fixed axial velocity',
             coarse_cells=args.cells//2,fine_cells=args.cells,bracket_width=BRACKET_WIDTH,
             drives=list(self.drives),engine=str(args.engine),engine_sha256=sha(args.engine),
             section=str(args.section),section_sha256=sha(args.section),command=self.command)
        self.reused = 0
        if args.reference:
            saved = args.reference
            status = json.loads((saved/'outcome.json').read_text())
            if status['exitCode'] or status['childOutcome']['exitCode'] or not status['limitsVerified']:
                raise ValueError('Unverified reference')
            rows = [json.loads(l) for l in (saved/'stdout').read_text().splitlines()]
            setup, = [r for r in rows if r['test']=='boundary_setup']
            audit, = [r for r in rows if r['test']=='boundary_audit']
            if (setup['seed']!=args.seed or setup['probe_count']!=len(self.probes)
                    or setup['command']!=self.command or setup['section_sha256']!=sha(args.section)
                    or setup['engine_sha256']!=sha(args.engine)
                    or sha(Path(audit['path']))!=audit['sha256']):
                raise ValueError('Reference has a different trajectory law or preparation')
            with Path(audit['path']).open() as stream:
                for line in stream:
                    r = json.loads(line)
                    if r['drive']==0.:
                        key = (r['probe'],0.,r['clearance'])
                        self.cache[key] = r['hit']
                        self.stream.write(line)
                        self.min_energy = min(self.min_energy,r['energy_tag_initial'])
                        self.max_defect = max(self.max_defect,r['max_energy_error'])
            self.reused = len(self.cache)
        if args.reference or args.control!='full' or args.side!=5 or args.full_strength_only:
            emit(test='boundary_control',mode=args.control,side=args.side,reused_undriven_samples=self.reused,
                 reference=None if args.reference is None else str(args.reference),
                 flux_weights_use_signed_drive=True)

    def evaluate(self, keys):
        keys = list(dict.fromkeys(k for k in keys if k not in self.cache))
        if not keys:
            return
        rows = []
        band = self.band
        for i,(probe,u,s) in enumerate(keys):
            x,y,yaw,q = self.probes[probe]
            cs,sn = np.cos(2*np.pi*yaw),np.sin(2*np.pi*yaw)
            velocity_drive = 0. if self.control=='order_only' else u
            order_drive = 0. if self.control=='kinematics_only' else u
            a = thermal_at_fixed_axial(band,band.launch(velocity_drive,1.,s),np.sqrt(TEMPERATURE)*ndtri(q))
            vx = a['tangent']*cs-a['polar']*sn
            vy = a['tangent']*sn+a['polar']*cs
            init = (i,x,y,1.5+a['z_offset'],vx,vy,a['normal'],band.scalar_order+order_drive*band.order_dipole,
                    cs,sn,.01*np.hypot(1,band.wall_slope)/PITCH)
            rows.append(' '.join(map(str,init)))
        result = subprocess.run(self.command,input='\n'.join(rows)+'\n',text=True,
                                check=True,stdout=subprocess.PIPE)
        for i,(key,line) in enumerate(zip(keys,result.stdout.splitlines(),strict=True)):
            r = json.loads(line)
            if r['id']!=i or r['energy_tag_initial']<=0:
                raise ArithmeticError('Unordered trajectory or initially bound participant')
            self.min_energy = min(self.min_energy,r['energy_tag_initial'])
            self.max_defect = max(self.max_defect,r['max_energy_error'])
            hit = r['gated_joint_residence']>0
            self.cache[key] = hit
            record = dict(probe=key[0],drive=key[1],clearance=key[2],hit=hit)
            for k in ('gated_joint_residence','first_joint_time','gate_exit_time','energy_tag_initial','max_energy_error'):
                record[k] = r[k]
            self.stream.write(json.dumps(record,allow_nan=False)+'\n')
        self.stream.flush()


def flux_integral(band, u, intervals):
    # A quadratic value times a quadratic w-gradient times r² is degree six
    # on each radial element. Four-point Gauss integrates it exactly.
    nodes, weights = leggauss(4)
    total = 0.
    for a,b in intervals:
        edges = np.r_[a, band.edges[(a<band.edges)&(band.edges<b)], b]
        half = np.diff(edges)/2
        s = edges[:-1,None]+half[:,None]*(1+nodes)
        total += (half[:,None]*weights*band.launch(u,1.,s)['weight']).sum()
    return float(total)


def boundaries(trials, cells):
    points = np.linspace(0.,1.,cells+1)**2
    ids = [(k,u) for k in range(len(trials.probes)) for u in trials.drives]
    changes = {}
    for k,u in ids:
        hits = [trials.cache[k,u,s] for s in points]
        changes[k,u] = [[a,b,ha] for a,b,ha,hb in zip(points[:-1],points[1:],hits[:-1],hits[1:],strict=True)
                         if ha!=hb]
    while True:
        keys = [(k,u,(a+b)/2) for (k,u),rs in changes.items() for a,b,_ in rs if b-a>BRACKET_WIDTH]
        if not keys:
            break
        trials.evaluate(keys)
        for (k,u),rs in changes.items():
            for bracket in rs:
                a,b,ha = bracket
                if b-a>BRACKET_WIDTH:
                    mid = (a+b)/2
                    bracket[0 if trials.cache[k,u,mid]==ha else 1] = mid
    for k,u in ids:
        roots = changes[k,u]
        edges = [0.,*[(a+b)/2 for a,b,_ in roots],1.]
        first = 0 if trials.cache[k,u,0.] else 1
        intervals = [(edges[i],edges[i+1]) for i in range(first,len(edges)-1,2)]
        mass = flux_integral(trials.band,u,intervals)
        uncertainty = flux_integral(trials.band,u,[(a,b) for a,b,_ in roots])
        emit(test='clearance_integral',cells=cells,probe=k,drive=u,draw=trials.probes[k].tolist(),
             intervals=intervals,mass_over_undriven_incident_rate=mass,
             localized_boundary_mass_uncertainty=uncertainty,
             start_label=trials.cache[k,u,0.],end_label=trials.cache[k,u,1.])


def run(args):
    if args.cells<4 or args.cells%2 or not 0<=args.power<=8 or abs(8/args.dt-round(8/args.dt))>1e-10:
        raise ValueError('Require even mesh, valid probe power, and dt dividing the observation')
    with args.output.open('x') as audit:
        trial = Trials(args,audit)
        points = np.linspace(0.,1.,args.cells+1)**2
        # Batch probes separately to bound captured stdout while retaining
        # every queried label for the nested and signed comparisons.
        for k in range(len(trial.probes)):
            trial.evaluate((k,u,s) for u in trial.drives for s in points)
        for cells in (args.cells//2, args.cells):
            boundaries(trial,cells)
        emit(test='boundary_audit',trajectories=len(trial.cache),seed=args.seed,
             min_incoming_tag_energy=trial.min_energy,max_total_energy_error=trial.max_defect,
             path=str(args.output),sha256=sha(args.output))


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--engine',type=Path,required=True)
    p.add_argument('--section',type=Path,required=True)
    p.add_argument('--seed',type=int,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--power',type=int,default=3)
    p.add_argument('--cells',type=int,default=32)
    p.add_argument('--dt',type=float,default=.002)
    p.add_argument('--side',type=int,choices=(5,8),default=5)
    p.add_argument('--full-strength-only',action='store_true',help='Use only the undriven case and ±.002')
    p.add_argument('--control',choices=('full','kinematics_only','order_only'),default='full')
    p.add_argument('--reference',type=Path,help='Reuse identical undriven trajectories from this verified block')
    run(p.parse_args())
