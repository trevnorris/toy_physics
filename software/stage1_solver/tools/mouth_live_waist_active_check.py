"""Open-field continuation into zero-shear cells and a completely cut modal cavity."""
from pathlib import Path
import json
import numpy as np
from mouth_sealed_shape_relax import Meridian
from mouth_live_waist_fields import CoupledLoads
from mouth_live_waist_active_wave import ActiveWave


def main():
    path=Path('research/pde_ledger_v3/_measurements/mouth_sealed_shape_relax.json')
    state=json.loads(path.read_bytes())['training']['2']['states'][-2]
    shape=Meridian(2,state['coordinates'])
    n=48
    x,t=np.meshgrid(np.linspace(0,1,n+1),np.linspace(0,1,2*n+1))
    chi=.5+.49*np.tanh((1-x)*shape.profile(t)[0]/.08)*np.tanh((1-t)*shape.length/.08)
    loads=ActiveWave.install(CoupledLoads(shape,n,2*n,12,target_frequency=state['wave_frequency']))
    start,_,baseline=loads.wave.response(chi)
    seed=loads.wave.seed
    for name,mask in [('local', (t>.28)&(t<.48)&(x>.3)&(x<.7)),('cut',(t>.28)&(t<.48))]:
        loads.wave.seed=seed; loads.wave.target=start
        last_positive=None
        for weight in (.25,.001,0.):
            order=chi.copy();order[mask]*=weight
            w,g,d=loads.wave.response(order)
            print(json.dumps(dict(case=name,remaining_order_fraction=weight,frequency=w,**d)),flush=True)
            if weight: last_positive=w
            else:
                assert d['inactive_nodes']>0 and d['mode_overlap']>.97 and abs(w-last_positive)<.0001
                assert bool(d['detached_mode'])==(name=='cut')
                assert np.isfinite(g).all()


if __name__=='__main__':
    main()
