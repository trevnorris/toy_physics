"""Independent thickness thresholds and Wronskian power from the archived toroidal port traces."""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import spherical_jn, spherical_yn

from mouth_exchange_shear import assemble
from mouth_poly_shallow_release import Trial
from mouth_sealed_shape_relax import Meridian
from mouth_shallow_puncture_fields import PORT


SCRATCH = Path('_scratch/s11c')


def read_rows(folder):
    guard = json.loads((folder/'outcome.json').read_bytes())
    if guard['exitCode'] or guard['childOutcome']['exitCode'] or not guard['limitsVerified']:
        raise ValueError(f'Unverified input {folder}')
    return list(map(json.loads,(folder/'stdout').read_bytes().splitlines()))


def inner(a,b,h):
    # Exact cell integral for the product of complex, piecewise-linear functions.
    ac = a.conj()
    return h/6*np.sum((2*ac[:-1]+ac[1:])*b[:-1] + (ac[:-1]+2*ac[1:])*b[1:])


def run(radius, n, detached):
    key = f'r{radius}n{n}'
    family = 'material' if detached else 'release'
    parent = SCRATCH/f'mouth_exit_{family}01'/key
    raw = (parent/'endpoint.npz').read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    if sha not in [p.get('snapshot_sha256') for p in read_rows(parent)]:
        raise ValueError('Unrecognized saved vector')
    scan = SCRATCH/'mouth_trapping_modes01'/(key+'l4'+('_detached' if detached else ''))
    setup, row, *_ = read_rows(scan)
    if (setup['sha256']!=sha or row['label']!='archived' or row['detached']!=detached):
        raise ValueError('Different archived reference for the angular sector')
    trial = Trial(radius,n,bonded_core=True)
    with np.load(parent/'endpoint.npz') as f:
        shape = Meridian(radius,trial.x+f['shifts'][:4])
        vector = f['vector']
        omega = f['pole'].item()
    N = trial.mesh
    mesh, *_ = assemble(shape,N,PORT,0.,int(np.ceil(2.77*radius*N)))
    mask = abs(mesh[:,0]-PORT)<1e-9
    ix = np.argsort(mesh[mask,1])
    z = mesh[mask,1][ix]
    trace = vector[mask][ix]
    H = z[-1]-z[0]
    h = H/N
    assert len(trace)==N+1 and max(abs(np.diff(z)-h))<2e-10
    # Neumann endpoint rows of the consistent Q1 pencil give the same
    # discrete cosine eigenvectors as the bulk three-point recurrence.
    indices = np.arange(N+1)
    angles = indices*np.pi/N
    cut = np.sqrt(12/h**2/(2+np.cos(angles)))*np.sin(angles/2)
    coeffs = []
    powers = []
    for j,cutoff in enumerate(cut):
        basis = np.cos(angles[j]*indices)
        basis /= np.sqrt(inner(basis,basis,h))
        amp = inner(basis,trace,h)
        coeffs.append(amp)
        if omega.real>cutoff:
            k = np.sqrt(omega.real**2-cutoff**2)
            arg = k*PORT
            norm2 = spherical_jn(4,arg)**2 + spherical_yn(4,arg)**2
            imag_traction = 1/(k*PORT**2*norm2)
            powers.append(float(.5*omega.real*PORT**4*abs(amp)**2*imag_traction))
    error_cut = max(abs(cut[:4]-row['depth_thresholds']))
    parseval = abs(inner(trace,trace,h)-np.vdot(coeffs,coeffs))
    if detached:
        scaled = powers
    else:
        old = read_rows(SCRATCH/'mouth_exit_referee01/leakage')
        i, = [i for i,p in enumerate(old) if p.get('endpoint_sha256')==sha]
        scaled = (np.array(powers)/old[i+2]['stored_energy_to_port']).tolist()
    diff = np.max(abs(np.array(scaled)-row['channel_power_per_energy']))
    print(json.dumps(dict(radius=radius,mesh=n,port_radius=PORT,height=H,wave_cells_in_height=N,
                          detached=detached,sha256=sha,frequency=omega.real,
                          continuum_first_odd_cutoff=np.pi/H,discrete_depth_cutoffs=cut[:4].tolist(),
                          depth_channels=len(powers),port_norm=inner(trace,trace,h).real,
                          power_per_stored_by_depth=scaled,cutoff_error=error_cut,
                          parseval_error=parseval,power_error=diff)),flush=True)
    if error_cut>1e-9 or parseval>1e-12 or diff>1e-9:
        raise ArithmeticError('Independent thickness-channel check failed')


if __name__=='__main__':
    for rad in (1,2):
        for mesh in (48,72):
            run(rad,mesh,False)
    run(2,48,True)
