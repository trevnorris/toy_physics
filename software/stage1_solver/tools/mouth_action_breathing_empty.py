"""Remove the archived wave's quasistatic work from the breathing-aperture curvature."""
import hashlib
import json
import numpy as np

from mouth_action_breathing import coordinates, exterior_probe
from mouth_exit_witness import inputs
from mouth_poly_shallow_release import ROOT, Trial, emit
from mouth_trapping_modes import Sector, pole


def run():
    parent=ROOT/'_scratch/s11c/mouth_action_breathing01/run'
    meta=json.loads((parent/'outcome.json').read_text())
    if not meta['limitsVerified'] or meta['exitCode'] or meta['childOutcome']['exitCode']:
        raise ValueError('Unvalidated breathing trial')
    raw=(parent/'stdout').read_bytes()
    emit(kind='origin',breathing_sha256=hashlib.sha256(raw).hexdigest())
    allrows=list(map(json.loads,raw.splitlines()))
    for r in (1,2):
        for n in (48,72):
            trial=Trial(r,n,bonded_core=True)
            path=ROOT/f'_scratch/s11c/mouth_exit_release01/r{r}n{n}/endpoint.npz'
            sh,chi,arch,digest,status,*_=inputs(path,r,n,trial)
            rows=[x for x in allrows if x['radius']==r and x['n']==n]
            setup,=[x for x in rows if x['kind']=='setup']
            if status!='balanced' or setup['endpoint_sha256']!=digest:raise ValueError('Input mismatch')
            b=setup['physical_exit']
            orig=json.loads((path.parent/'stdout').read_text().splitlines()[0])['exit_radius']
            ee={}; masses={}
            for x in rows:
                if x['kind']!='basis_state':continue
                u=x['relative_exit']
                shift,shape=coordinates(trial,sh,b*(1+u))
                op=Sector(trial,shift,chi,4)
                # The archived wave omits axis nodes; Sector preserves an all-point map.
                full=np.zeros(len(op.live),complex)
                full[op.allpoints[:,0]>0]=arch[1]
                seed,data=pole(op,(arch[0],full[op.live]))
                if data['seed_overlap']<.97 or (u==0 and abs(seed[0]-arch[0])>1e-7):
                    raise ArithmeticError('Occupied profile changed along empty-state subtraction')
                masses[u]=x['mass'][1][1]/orig**2
                ee[u]=x['energy']-seed[0].real
            for h in (.012,.006):
                freq=[]
                for i in (-1,0,1):
                    curv=(ee[(i+1)*h]+ee[(i-1)*h]-2*ee[i*h])/(b*h)**2
                    freq.append(np.sqrt(curv/masses[i*h]))
                emit(kind='shape_held_empty_exit',radius=r,n=n,step=h,
                     empty_slope=(ee[h]-ee[-h])/(2*h*b),
                     trial_frequencies=freq,domega_dexit=(freq[2]-freq[0])/(2*h*b),
                     **exterior_probe(freq[1],b))


if __name__=='__main__':
    run()
