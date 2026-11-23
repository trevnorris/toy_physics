"""Directional derivative checks of the reciprocal mixed-order energy at the R4 rim."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_last_rings_fluid import GRIDS, read_source
from mouth_material_order_transfer import MixedOrder
from mouth_shared_rim_shift import MovedCarrier, SharedRimChart
from mouth_supercore_joined_front import file_sha


def run(state,wall,axes):
    x, disp, u, _ = read_source(state,wall)
    base = Front(state,GRIDS['native'],1)
    m = Trial(base,state,x,11,32,False).model
    if u.shape != m.seed.shape:
        raise ValueError('Native field topology does not match material checkpoint')
    # Work with the physical block plus a smooth nonstationary order test,
    # supported on free nodes in the terminal duct and receiving region.
    rng = np.random.default_rng(8279)
    field = u.copy()
    perturb = np.zeros_like(u)
    for k in (0,1):
        ids = m.ids[1,k][:,:,-5:]
        perturb[:,ids] = rng.normal(size=(3,*ids.shape))
        ids = m.ids[2,k][:,:,:5]
        perturb[:,ids] = rng.normal(size=(3,*ids.shape))
    perturb *= m.masks
    field += 3.e-4*perturb
    chart = SharedRimChart(m,m,axes)
    old = MixedOrder(m,base.c)
    rows=[]
    for n in (0,4,8):
        carrier = chart.stationary if n==0 else (lambda side,r,z:chart.carrier(side,r,z,n))
        material = chart.still if n==0 else chart.physical
        moved = MovedCarrier(m,chart,carrier,material)
        for step in (-2e-5,2e-5):
            moved.apply(step)
            old.install(step,material,carrier,m.material_load)
            energy, force = m.evaluate(field)
            native_energy, native_force = old.evaluate_native(field)
            # The common displacement also strains the material. The native
            # difference isolates the part implemented by the mixed operator.
            metrics = []
            for k in range(3):
                delta=np.zeros_like(field)
                delta[k]=perturb[k]
                h=3.e-6
                xp,xm=field+h*delta,field-h*delta
                ep,em=(m.evaluate(v)[0] for v in (xp,xm))
                enp,enm=(old.evaluate_native(v)[0] for v in (xp,xm))
                fd = ((ep-enp)-(em-enm))/(2*h)
                analytic = (force[k]-native_force[k])@delta[k]
                metrics.append(dict(field=k,finite=fd,analytic=analytic,error=fd-analytic))
                if abs(fd-analytic)>1.5e-6:
                    raise ArithmeticError(f'Mixed energy/gradient check failed: {k} {fd} {analytic}')
            rows.append(dict(width=n,step=step, correction=energy-native_energy,
                             reconstructed_blocks=old.active,probes=metrics))
            old.reset()
            moved.reset()
    return dict(wall_sha=file_sha(wall),state_sha=file_sha(state),axes_sha=file_sha(axes),checks=rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('state','wall','axes'):
        parser.add_argument(name,type=Path)
    a=parser.parse_args()
    print(json.dumps(run(a.state,a.wall,a.axes),allow_nan=False))
