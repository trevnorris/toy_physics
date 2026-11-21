"""Independently calibrate the carrier on one mesh and read hidden loads on another.

Balanced auxiliary actuators calibrate the force readout. Their forces are
external known standards and are never used to supply a material attachment.
"""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.sparse.linalg import splu

from mouth_carrier_lateral_response import prepare
from mouth_flared_receiver import emit
from mouth_free_bend_drive import Response
from mouth_free_bend_release import pencil
from mouth_two_depth_load import PEAK_FORCE


def setup(mesh, extension):
    path = Path(f'_scratch/s11c/mouth_flared_cubic_n{mesh}')
    dom, field = prepare(SimpleNamespace(source=path))
    return Response(dom, field, extension)


def gauges_to_coordinates(gauges):
    return np.vstack((gauges[0],np.diff(gauges,axis=0)))


def retained_motion(label):
    path=Path(f'_scratch/s11c/two_depth_observation_{label}_01')
    outcome=json.loads((path/'outcome.json').read_text())
    guard=json.loads((path/'limit-validation.json').read_text())
    if (outcome['exitCode'] or outcome['childOutcome']['exitCode'] or not outcome['limitsVerified']
            or not guard['verified'] or guard['actual']['memory.swap.max']!='0'):
        raise ValueError('Unverified time history')
    records=[json.loads(line) for line in (path/'stdout').read_text().splitlines()]
    row, = [r for r in records if r['test']=='motion' and r['elapsed']==96.]
    return np.array(row['per_peak_force_real'])+1j*np.array(row['per_peak_force_imag'])


def main():
    rate=.30+.40j
    coarse=setup(26,'compact').force_matrix(rate)
    fine=setup(34,'wide')
    read=np.tri(3)
    ports=read.T
    # Published drives are at middle and flare. Each also receives a separate
    # equal/opposite force pair with known magnitude and phase for verification.
    primary=ports[:,1:].copy()
    auxiliary=np.array([[-.08j,0.], [.08j,-.05], [0.,.05]],complex)
    energy,time=pencil(fine,0.)
    rhs=np.zeros((energy.shape[0],2),complex)
    rhs[-3:]=primary+ports@auxiliary
    state=splu((energy-rate*time).tocsc()).solve(rhs)
    gauges=read@state[-3:]
    fitted=np.linalg.solve(ports,coarse@gauges_to_coordinates(gauges)-primary)
    emit(test='auxiliary_force', peak_scale=PEAK_FORCE,
         applied_real=auxiliary.real.tolist(), applied_imag=auxiliary.imag.tolist(),
         inferred_real=fitted.real.tolist(), inferred_imag=fitted.imag.tolist(),
         max_error_per_primary=float(np.max(abs(fitted-auxiliary))),
         peak_relative_error_by_drive=(np.max(abs(fitted-auxiliary),axis=0)/np.max(abs(auxiliary),axis=0)).tolist())
    if not np.isfinite(fitted).all() or np.max(abs(fitted-auxiliary))>.01:
        raise ArithmeticError('Hidden forces missed one percent of the primary load')
    measured=retained_motion('n34')
    ordered=measured[:,[2,0,1]].copy()
    ordered[:,0]/=-.5j
    q=gauges_to_coordinates(ordered)
    difference=np.diff(q,axis=1)
    residual=ports-coarse@q
    emit(test='force_resolution', singular_values=np.linalg.svd(q[1:,1:],compute_uv=False).tolist(),
         condition_single=float(np.linalg.cond(q[1:,1:])),
         condition_balanced=float(np.linalg.cond(difference[1:])),
         null_error_per_primary=np.max(abs(np.linalg.solve(ports,residual)),axis=0).tolist())
    # Peak times or flight times cannot be read out of wrapped phase.
    for row, name in enumerate(('entrance','middle','flare')):
        emit(test='station_summary', station=name,
             amplitude_per_peak=abs(measured[row,:2]).tolist(),
             signed_lag_degrees=(-np.angle(measured[row,:2],deg=True)).tolist())


if __name__=='__main__':
    main()
