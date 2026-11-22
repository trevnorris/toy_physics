"""Validate two-sided shape work against signed energy and a distinct drive."""
import hashlib
import json
from pathlib import Path

ROOT = Path('_scratch/s11c/exchange_carrier01')


def digest(file):
    return hashlib.sha256(file.read_bytes()).hexdigest()


def gather(name, worker_code=0):
    folder=ROOT/name
    def read(key):
        return json.loads((folder/key).read_text())
    out=read('outcome.json')
    inv=read('invocation.json')
    lim=read('limit-validation.json')
    dur=read('duration-validation.json')
    if (not out['limitsVerified'] or out['exitCode']!=worker_code or out['childOutcome']['exitCode']!=worker_code
        or not lim['verified'] or lim['actual']['memory.swap.max']!='0' or not dur['verified']
        or dur['wallDeadlineSeconds'] is not None):
        raise ValueError(f'Invalid guarded carrier readout {name}')
    txt=(folder/'stdout').read_text()
    info=dict(worker=inv['command'],guard_exit=out['exitCode'],child_exit=worker_code,
              mem_bytes=inv['memoryMax'],swap_bytes=0,wall_deadline=None,stdout_sha256=digest(folder/'stdout'),
              logdir=str(folder),pool=inv['pool'])
    return [json.loads(s) for s in txt.splitlines()],info


def signed(m,z,p,h):
    energies=(p['energy']-m['energy'])/(2*h)
    K=(p['energy']-z['energy']+m['energy']-z['energy'])/h**2
    virtual=(p['work']-m['work'])/(2*h)
    by_region={k:(p['sector_work'][k]-m['sector_work'][k])/(2*h) for k in z['sector_work']}
    if abs(energies-z['work'])>2e-5 or abs(virtual-K)>1e-4:
        raise ArithmeticError('Integrated field stress and varied stationary energy disagree')
    return dict(half_width=h,energy_slope=energies,virtual_slope=z['work'],
                energy_stiffness=K,virtual_stiffness=virtual,sector_stiffness=by_region)


def main():
    history={}
    weak={}
    for name in ('n40','n64','n80','widefront_b','plain','box','interior','interior_plain','n24'):
        entries,history[name]=gather(name,1 if name=='n24' else 0)
        weak[name]=entries
    # The coarse job produced all states, then stopped on an unextrapolated
    # central derivative at h=.01. Its results are not used as a passed run.
    coarse=weak.pop('n24')
    history['n24']['failed_check']=(ROOT/'n24/stderr').read_text().splitlines()[-1]
    history['n24']['raw_full_step']=coarse[1]
    history['n24']['source_sha256']=digest(ROOT/'first.py')
    for name in ('n40','n64','n80','widefront_b','plain','box'):
        history[name]['source_sha256']=digest(ROOT/'reference_solver.py')
    history['widefront']={'outcome':'guard refused admission during cgroup teardown; no worker or scientific output'}
    finite={}
    rows,history['f24_0']=gather('f24_0')
    finite['coarse_center']=rows[0]
    for mesh,drive,steps in ((40,1,(0,-.005,.005)),(64,1,(0,-.005,.005,-.01,.01)),(64,.5,(0,-.005,.005))):
        group={}
        for s in steps:
            name=f'finite_{mesh}_{s}_0_{drive}'
            rows,history[name]=gather(name)
            group[s]=rows[0]
        pairs=[signed(group[-h],group[0],group[h],h) for h in steps if h>0]
        finite[f'mesh{mesh}_drive{drive}']=dict(states=list(group.values()),signed=pairs)
    interior=[]
    for s in (-.005,.005):
        name=f'finite_64_{s}_1_1'
        rows,history[name]=gather(name)
        interior.append(rows[0])
    slope=(interior[1]['energy']-interior[0]['energy'])/.01
    mean=(interior[1]['work']+interior[0]['work'])/2
    assert abs(slope-mean)<2e-6
    finite['held_ends']=dict(states=interior,slope=slope,mean_work=mean)
    tool=Path(__file__).parent
    paths=[tool/name for name in ('mouth_exchange_carrier.py','mouth_exchange_carrier_density.py',
           'mouth_exchange_carrier_interior.py','mouth_exchange_geometry.py','mouth_exchange_interior.py',
           'mouth_conversion_shape_mode.py','mouth_exchange_carrier_readout.py')]
    print(json.dumps(dict(base='017b95e2',scope='R1 turning flare, stationary two-sided thin-layer Route A',
                          history=history,weak=weak,finite=finite,
                          source_sha256={str(p):digest(p) for p in paths}),indent=2,allow_nan=False))


if __name__=='__main__':
    main()
