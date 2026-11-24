"""Compare independently settled R4 outputs on the two numerical gauge strengths."""
import argparse
from pathlib import Path

import numpy as np

from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_whole_r4_relax import SHAPE_TOL, atomic_json


def compare(area, name):
    rows = {}
    states = {}
    ends = []
    for epsilon, case in ((1.e-8,'primary'), (1.e-10,'secondary')):
        case += '' if name=='cp11' else '_final'
        path = area/case
        guard = read(area/(case+'_guard')/'outcome.json')
        if guard['exitCode'] or not guard['limitsVerified'] or guard['childOutcome']['exitCode']:
            raise ValueError(f'Failed run {case}')
        ending = 'physical_stop' if (path/'physical_stop.json').exists() else (
                  'equilibrium' if (path/'equilibrium.json').exists() else 'cp11')
        ends.append(ending)
        if ending=='cp11' and name!='cp11':
            raise ValueError('Missing terminal state')
        i = 11 if ending=='cp11' else read(path/(ending+'.json'))['iteration']
        r = read(path/f'checkpoint_{i:05d}.json')
        if (r['input_checkpoint_sha256']!=CHECKPOINT_SHA or r['epsilon']!=epsilon or
            r['phase_drop']!=.4510084282857874 or not r['material']['stationary'] or
            r['material']['min_chart_Jinv_bound']<=0 or r['material']['min_Jinv_global']<=0 or
            digest(path/r['file'])!=r['sha256']):
            raise ValueError(f'Unverified gauge preparation {case}')
        cur=r['current']
        rows[epsilon] = dict(iteration=i, energy=r['energy'],physical_energy=r['energy']-r['material']['gauge_energy'],
            gauge_energy=r['material']['gauge_energy'],forces=r['shape_force'],x=r['x'] if 'x' in r else None,
            exit=r['exit_radius'],depth=r['depth'],front=r['sampled_melted_front'],
            current=cur['feed_current'],subsonic=cur['dense_total_subsonic_share'],
            middle_change=abs(r['shape_force'][2])-.70549897,exit_force_change=abs(r['shape_force'][1])-.64123594)
        with np.load(path/r['file']) as saved:
            rows[epsilon]['x'] = saved['x'].tolist()
        states[str(epsilon)] = r['sha256']
    a,b=rows.values()
    # Absolute tolerances in solver units. The energy is ~9.2e4 due to the
    # common far bath, so a relative tolerance there would mask real differences.
    limits=dict(energy=1.e-6,physical_energy=1.e-6,forces=1.e-4,x=1.e-4,exit=1.e-4,
                depth=1.e-4,front=1.e-4,current=3.e-5,subsonic=1.e-4)
    differences={}
    for key in limits:
        if a[key] is None or b[key] is None:
            differences[key]=0. if a[key] is None and b[key] is None else None
        else:
            differences[key]=float(np.max(np.abs(np.asarray(a[key])-b[key])))
    passed = ends[0]==ends[1] and all(v is not None and v<=limits[k] for k,v in differences.items())
    if ends==['equilibrium','equilibrium']:
        passed = passed and all(max(abs(x) for x in r['forces'])<=SHAPE_TOL for r in rows.values())
    result=dict(passed=passed,phase=name,outcomes=ends,states=states,rows=rows,
                absolute_differences=differences,limits=limits)
    atomic_json(area/('gauge_11.json' if name=='cp11' else 'gauge_final.json'),result)
    if not passed:
        raise ArithmeticError(f'Gauge or physical outcome differs: {differences}, {ends}')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('area',type=Path)
    p.add_argument('phase',choices=('cp11','final'))
    a=p.parse_args()
    compare(a.area,a.phase)
