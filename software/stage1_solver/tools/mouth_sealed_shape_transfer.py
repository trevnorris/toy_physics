"""Continue the occupied wave through the accepted radius-two shapes on a finer grid."""
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_sealed_shape_relax import Experiment, emit


def states(path):
    out=json.loads((path/'outcome.json').read_text())
    if out['exitCode'] or out['childOutcome']['exitCode'] or not out['limitsVerified']:
        raise ValueError('Unverified source geometry')
    data=(path/'stdout').read_bytes()
    emit(test='geometry_input',source=str(path),sha256=hashlib.sha256(data).hexdigest())
    return [r for r in map(json.loads,data.splitlines()) if r['test']=='state']


def main():
    base=Path('_scratch/s11c/route_a_duct_relax01')
    path=states(base/'r2m10s')[-2:-1]+states(base/'r2continuation')
    exp=Experiment(2,16,5.2)
    for i,selected in enumerate(path):
        x=np.array(selected['coordinates'])
        if i:
            distance=float(np.linalg.norm(x-old))
            steps=max(1,int(np.ceil(distance/.045)))
            xs=[old+factor*(x-old) for factor in np.linspace(0.,1.,steps+1)[1:]]
        else:
            xs=[x]
        for xx in xs:
            row,v,mass=exp.fields(tuple(xx))
            _,last=exp.previous_mode
            overlap=1. if last is None else abs(np.vdot(last,mass@v))**2/np.vdot(last,mass@last).real
            if overlap<.9:
                raise ArithmeticError('Finer-mesh occupied mode switched')
            exp.previous_mode=complex(row['wave_frequency'],-row['wave_rate']/2),v
        emit(test='independent_mesh', index=i,coordinates=x.tolist(),overlap=overlap,**row)
        old=x


if __name__=='__main__':
    main()
