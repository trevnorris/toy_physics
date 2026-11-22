"""Finite-sample kinetic plus potential energy and the conditional inherited time conversion."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from statistics import mean


def read(path): return json.loads(path.read_bytes())


def clock():
    rec=read(Path('research/pde_ledger_v3/_measurements/mouth_attachment_impact.json'))['guarded_runs']['ax15_fine']
    saved=Path(rec['stdout_path']).read_bytes()
    assert rec['guard_exit']==rec['child_exit']==0 and rec['limits_verified'] and sha256(saved).hexdigest()==rec['stdout_sha256']
    setting=json.loads(saved.splitlines()[0]);assert setting['test']=='preparation'
    return setting['model_time']/setting['native_time']


def energy(r):
    # Four physical velocity components with unit constituent mass. Kinetic T uses all velocities.
    last=r['blocks'][-4:]; count=sum(b['steps'] for b in last)
    return sum((b['mean']['potential_per_site']+2*b['mean']['temperature'])*b['steps'] for b in last)/count


def cases(root):
    for name, signature in read(root/'DONE.json')['results'].items():
        content=(root/'runs'/(name+'.json')).read_bytes()
        assert sha256(content).hexdigest()==signature
        r=json.loads(content)
        if r['pressure_external']==r['strain']==0.:
            yield r


def main(primary,screen,out):
    rate=clock();ref=list(cases(primary)); times={r['measure'] for r in ref};assert len(times)==1 and len(ref)==6
    transit=times.pop(); baselines={h:mean(energy(r) for r in ref if r['history']==h) for h in range(3)}
    point=[]
    for r in sorted(cases(screen),key=lambda r:(r['temperature'],r['history'])):
        v=energy(r);point.append(dict(temperature=r['temperature'],history=r['history'],energy_per_site=v,excess_over_same_history_at_1_024=v-baselines[r['history']]))
    result=dict(conditional_mouth_time_per_native_time=rate,nominal_transit_mouth_time=transit,native_time_for_transit_if_transferred=transit/rate,
                fraction_of_transit_covered_if_transferred=rate,primary_late_energy_per_site=baselines,screen=point,
                scope='Conditional impact-calibrated clock if constituent mass, pitch and well energy transfer. The gamma=1 thermal friction is a model assumption. Energies are endpoint-window comparisons at externally imposed zero pressure on finite prepared samples; no heat of conversion inferred.')
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,allow_nan=False))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('primary',type=Path);p.add_argument('screen',type=Path);p.add_argument('out',type=Path)
    a=p.parse_args();main(a.primary,a.screen,a.out)
