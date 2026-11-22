"""Reconstruct 4D slab budgets and sampled stress from saved mobile-coupon data."""
import argparse
import json
from pathlib import Path

import numpy as np


def verify(r):
    prev=r['prepared']; max_storage=max_profile=max_endpoint=0.; all_heat=all_cell=0.;attempts=accepted=0
    sources=np.array(r['origin']); n=r['n']; refrows=np.array(prev['box'])[:3]; assert np.max(abs(refrows[:,3]))<1.e-12
    for b in r['blocks']:
        faces=b['balance'];nb=len(faces);assert nb==r['bins'] and len(b['profile'])==nb
        forward=np.array([f['p_plus'] for f in faces]);reverse=np.array([f['p_minus'] for f in faces]);flux=forward-reverse
        nc=np.array([np.array(f['cross_plus'])-f['cross_minus'] for f in faces]);assert nc.shape==(nb,2)
        diff=np.array([np.array(f['p1'])-f['p0'] for f in faces]);external=np.array([np.array(f['pair'])+f['bath'] for f in faces]);res=diff-external-flux+np.roll(flux,-1,axis=0)
        max_storage=max(max_storage,float(np.abs(res).max()));np.testing.assert_allclose(res,[f['residual'] for f in faces],atol=1.e-10)
        end=b['state']; hist=[]
        for state,key in ((prev,'0'),(end,'1')):
            box=np.array(state['box']);np.testing.assert_array_equal(box[:3],refrows)
            x=np.array(state['pos']);v=np.array(state['vel']);labels=np.floor(nb*x[:,3]/box[3,3]).astype(int)%nb
            count=np.array([np.bincount(labels[sources==parent],minlength=nb) for parent in range(2)]).T;hist.append(count)
            np.testing.assert_array_equal(count.sum(axis=1),[f['n'+key] for f in faces]);assert count.sum()==n
            ep=np.array([np.bincount(labels,weights=v[:,k],minlength=nb) for k in range(4)]).T
            max_endpoint=max(max_endpoint,float(np.abs(ep-np.array([f['p'+key] for f in faces])).max()))
        np.testing.assert_array_equal(hist[1]-hist[0],nc-np.roll(nc,-1,axis=0))
        local=np.array([k['pressure_w'] for k in b['profile']]);global_w=-np.array(b['mean']['stress_tensile'])[:,3];max_profile=max(max_profile,float(np.abs(local.mean(axis=0)-global_w).max()))
        np.testing.assert_allclose(np.mean([u['rho'] for u in b['profile']]),b['mean']['density'],atol=1.e-10)
        all_heat+=b['bath_heat'];all_cell+=b['cell_dU'];attempts+=b['attempts'];accepted+=b['accepted'];prev=end
    assert 0<accepted<attempts and max_storage<1e-6 and max_profile<1e-8 and max_endpoint<1e-8
    def energy(t):return n*(t['potential_per_site']+2*t['temperature'])
    residual=energy(r['blocks'][-1]['end_reading'])-energy(r['prepared_reading'])-all_cell-all_heat
    return dict(sites=n,blocks=len(r['blocks']),bins_per_block=nb,measured_steps=sum(b['steps'] for b in r['blocks']),
        normal_cell_attempts=attempts,accepted=accepted,rejected=attempts-accepted,
        max_abs_four_momentum_residual=max_storage,max_abs_profile_to_whole_normal_traction=max_profile,
        max_abs_saved_endpoint_momentum_error=max_endpoint,bath_heat=all_heat,accepted_cell_potential_change=all_cell,
        integrated_energy_residual_per_site=residual/n)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('raw',type=Path);args=p.parse_args()
    print(json.dumps(verify(json.loads(args.raw.read_bytes())),allow_nan=False))
