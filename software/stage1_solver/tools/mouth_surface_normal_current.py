"""Test steady cross-sheet current under a phase bias in the same flat two-bath candidate."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_bvp

from mouth_curved_site_front import sha
from mouth_surface_carrier_reflection import PITCH, sheet_mean, stationary
from neutral_patch_sites import emit


def run(folder):
    out = json.loads((folder/'outcome.json').read_bytes())
    assert out['exitCode'] == out['childOutcome']['exitCode'] == 0 and out['limitsVerified']
    lines = [json.loads(line) for line in (folder/'stdout').read_text().splitlines()]
    setup, = [r for r in lines if r['test']=='setup']
    equilibrium, = [r for r in lines if r['test']=='base' and r['strength']==1]
    height = setup['far']
    nb = setup['background_density']
    def potential(x):
        return setup['energy_scale']*sheet_mean(x/PITCH,setup['pitch'])[0]

    mesh = np.unique(np.r_[np.linspace(0,.06,151),np.linspace(.06,1,111),np.linspace(1,height,500)])
    base = stationary(mesh,potential,nb,7e-8)
    x = base.x
    f, df = base.y
    dist = cumulative_trapezoid(1/f**2,x,initial=0.)
    def ode(z, values):
        a, b, da, db = values
        twice_h = 2*potential(z)+((a*a+b*b)**4-nb**4)/2
        return np.array([da,db,twice_h*a,twice_h*b])

    probe = np.r_[np.linspace(0,.03,30),np.linspace(.1,height,80)]
    readings = {}
    for delta in (.002, -.002, .001, -.001, 0.):
        phase = (delta/2)*dist/dist[-1]
        deriv = (delta/2)/(dist[-1]*f*f)
        trial = np.array([f*np.cos(phase),f*np.sin(phase),
                          df*np.cos(phase)-f*np.sin(phase)*deriv,
                          df*np.sin(phase)+f*np.cos(phase)*deriv])
        def bc(left,right):
            return np.array([left[2],left[1],right[0]-np.sqrt(nb)*np.cos(delta/2),
                             right[1]-np.sqrt(nb)*np.sin(delta/2)])
        wave = solve_bvp(ode,bc,x,trial,tol=7e-8,max_nodes=50000)
        if not wave.success:
            raise ArithmeticError(wave.message)
        a,b,da,db = wave.sol(probe)
        currents = a*db-b*da
        flux = float(currents[0])
        readings[delta] = flux
        emit(test='steady_normal', height=height, total_phase_bias=delta,
             current=flux, max_current_difference=np.max(abs(currents-flux)),
             density_at_sheet=a[0]**2+b[0]**2, both_bath_densities=nb, nodes=len(wave.x))
        assert np.max(abs(currents-flux)) < 5e-10*max(abs(flux),1e-2)
    excess = equilibrium['excess_phase_length']
    predicted = nb/(2*height+excess)
    emit(test='linear_check', height=height, source_stdout_sha256=sha(folder/'stdout'),
         predicted_current_per_total_phase=predicted, zero_bias_current=readings[0.],
         signed_full=(readings[.002]-readings[-.002])/.004,
         signed_half=(readings[.001]-readings[-.001])/.002,
         current_per_excess_phase=nb/excess)
    for amp in (.002,.001):
        assert np.isclose((readings[amp]-readings[-amp])/(2*amp),predicted,rtol=5e-5,atol=0)
    assert abs(readings[0.])<1e-12


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('collective_run',type=Path)
    run(parser.parse_args().collective_run)
