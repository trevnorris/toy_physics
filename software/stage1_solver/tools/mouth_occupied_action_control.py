"""Test successive doublings of the wave action on the same draining geometry."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh

from mouth_flared_receiver import emit


def complex_row(row):
    return np.asarray(row['real'])+1j*np.asarray(row['imag'])


def run(folder):
    outcome = json.loads((folder/'outcome.json').read_text())
    child = json.loads((folder/'child-outcome.json').read_text())
    if (not outcome['limitsVerified'] or outcome['exitCode'] or child['exitCode']
            or child != outcome['childOutcome'] or (folder/'stderr').stat().st_size):
        raise ValueError('Unverified primary field calculation')
    rows = [json.loads(s) for s in (folder/'stdout').read_text().splitlines()]
    wave, = [r for r in rows if r['test']=='wave']
    static, = [r for r in rows if r['test']=='static']
    field, = [r for r in rows if r['test']=='preparation']
    kc = np.array(static['carrier'])
    kw = np.array(wave['static'])
    baseline = wave['action']
    # The source rate is a specified control. The same unit-density,
    # unit-thickness slab zero mode gives its feeding-brane incident speed.
    source_rate = .25
    separation = 64.
    u = source_rate/(4*np.pi*separation**2)
    emit(test='control', mesh=field['mesh'], target_intake=field['rate'], source_intake=source_rate,
         separation=separation, leading_incident_speed=u, source='feeding slab zero mode')
    for action in (2*baseline,4*baseline):
        k = kc+action*kw
        krel = k[1:,1:]-np.outer(k[1:,0],k[0,1:])/k[0,0]
        least = min(eigvalsh(krel))
        emit(test='occupation', action=action, added_rigid_inertia=action*wave['mass'][0][0],
             least_static_bend=least, independent_wall_tension=0.)
        if least<=0 or k[0,0]<=0:
            continue
        for row in rows:
            if row['test']!='operators':
                continue
            nu = row['frequency']
            z = complex_row(row['carrier'])+action*np.array(row['wave'])
            drive = complex_row(row['drive'])
            x = np.linalg.solve(z,drive)
            rigid = -nu**2*drive[0]/z[0,0]
            accel = -nu**2*x[0]
            # Eliminate both bends separately to verify the acceleration.
            abr = z[0,1:]@np.linalg.solve(z[1:,1:],np.column_stack((drive[1:],z[1:,0])))
            by_elimination = -nu**2*(drive[0]-abr[0])/(z[0,0]-abr[1])
            np.testing.assert_allclose(by_elimination,accel,rtol=2e-10)
            # Independently solve the physically scaled forcing in station
            # coordinates (entrance, middle, flare) instead of bend offsets.
            transform = np.array([[1.,0.,0.],[-1.,1.,0.],[0.,-1.,1.]])
            at_source = np.linalg.solve(transform.T@z@transform,transform.T@drive*u)
            np.testing.assert_allclose(at_source,u*np.cumsum(x),rtol=2e-9)
            emit(test='occupied_wave_only', action=action, frequency=nu,
                 entrance_a_per_u=[accel.real,accel.imag], rigid_a_per_u=[rigid.real,rigid.imag],
                 in_phase_departure=(accel.real-rigid.real)/rigid.real,
                 complex_rigid_relative_error=abs(accel-rigid)/abs(rigid),
                 flare_mismatch=abs(x[1:].sum()/x[0]),
                 entrance_acceleration_at_separation=[-nu**2*at_source[0].real,-nu**2*at_source[0].imag],
                 peak_entrance_displacement_at_separation=abs(at_source[0]))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('guarded_input',type=Path)
    run(parser.parse_args().guarded_input)
