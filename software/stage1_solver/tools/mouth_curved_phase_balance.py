"""Read both sides of the held curved graph and check the proposed normal load."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import simpson
from mouth_curved_phase import Domain


def balance(folder):
    outcome=json.loads((folder/'outcome.json').read_text())
    assert outcome['limitsVerified'] and outcome['exitCode']==outcome['childOutcome']['exitCode']==0
    state=json.loads((folder/'stdout').read_text().splitlines()[0])
    domain=Domain(state['n'],state['far'])
    with np.load(folder/'field.npz') as packet:
        np.testing.assert_array_equal(domain.points,packet['points'])
        field=packet['psi']; r=packet['r']; h=packet['h_0.4']; q=packet['q_0.4']
    sigma,K=1/3,.4
    surface=CubicSpline(r,h)
    rad=.5
    z=surface(rad); slope=surface(rad,1); normal=np.array([-slope,1])/np.hypot(1,slope)
    for distance in (.02,.01,.005):
        samples=np.array([[-3,-2,-1],[3,2,1]])*distance
        rs=rad+samples*normal[0]; zs=z+samples*normal[1]
        assert np.min(np.hypot(domain.rc-rs,domain.rounding-zs))>domain.rounding
        val,*_=domain.read(field,rs,zs); phase=np.angle(val); pressure=abs(val)**10/5
        # Separately extrapolate each side from three physical positions.
        weights=np.array([1,-3,3]); p_limits=pressure@weights; phi_limits=phase@weights
        print(json.dumps(dict(step=distance,feed_phase=phi_limits[0],bulk_phase=phi_limits[1],
            feed_pressure=p_limits[0],bulk_pressure=p_limits[1],pressure_difference_across_samples=pressure[0,-1]-pressure[1,-1],
            extrapolated_feed_minus_bulk=p_limits[0]-p_limits[1])),flush=True)
    projected=4*np.pi/3
    actual_paired_trace=0.0  # Same continuous n and EOS on this transparent internal graph.
    required=sigma/K*simpson(4*np.pi*r**2*q,x=r)
    exact=sigma*state['rate']/K
    print(json.dumps(dict(integrated_required=required,required_from_bath_rate=exact,mean_required=exact/projected,
        required_at_midradius=sigma/K*np.interp(rad,r,q),same_material_jump=actual_paired_trace,
        min_sampled_before_rim=q[:-1].min(),all_interior_points_inward=bool(np.all(q[:-1]>0)),
        interior_sample_count=int(len(q)-1))),flush=True)
    assert np.all(q[:-1]>0) and abs(required-exact)<2e-5 and abs(p_limits[0]-p_limits[1])<1e-5
    assert abs(phi_limits[0]-phi_limits[1])<1e-5


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('folder',type=Path)
    balance(parser.parse_args().folder)
