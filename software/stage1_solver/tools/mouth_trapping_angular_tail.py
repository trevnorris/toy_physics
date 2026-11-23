"""Continue the degree-six cavity pole in angular order on each saved bonded shape."""
import argparse
import copy
import hashlib
import json

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.special import hankel1e

from mouth_exit_witness import inputs
from mouth_poly_shallow_release import ROOT, Trial, emit
from mouth_shallow_puncture_fields import PORT
from mouth_trapping_modes import Sector, pole


class AngularPort:
    """Curl-form outgoing traction at real angular order, including evanescent channels."""
    def __init__(self, port, degree):
        self.ids, self.cutoff, self.project = port.ids, port.cutoff, port.project
        self.degree = degree

    def traction(self, omega):
        k = np.sqrt(omega**2-self.cutoff**2+0j)
        k = np.where((omega.real < self.cutoff) & (k.imag < 0), -k, k)
        nu = self.degree+.5
        z = PORT*k
        return k*hankel1e(nu-1, z)/hankel1e(nu, z)-self.degree/PORT

    def block(self, omega):
        return PORT**4*(self.project.T*self.traction(omega))@self.project


def recorded_seed(op, row):
    """Reconstruct the archived pole vector at its measured frequency by a bordered linear solve."""
    omega = complex(row['frequency'], -row['loss']/2)
    use = op.labels == op.connected[0]
    k, m = op.K[use][:, use], op.M[use][:, use]
    ids = np.searchsorted(np.flatnonzero(use), op.port.ids)
    coords = np.repeat(ids, len(ids)), np.tile(ids, len(ids))
    a = k-omega**2*m-sparse.coo_matrix((op.port.block(omega).ravel(), coords), shape=m.shape)
    s = sparse.diags(1/np.sqrt(m.diagonal()))
    z = np.random.default_rng(740).standard_normal(m.shape[0])
    extended = sparse.bmat([[s@a@s, sparse.csc_matrix(z[:, None])],
                            [sparse.csc_matrix(z[None, :]), None]], format='csc')
    part = s@spsolve(extended, np.r_[np.zeros_like(z), 1.])[:-1]
    v = np.zeros(op.K.shape[0], complex)
    v[use] = part/np.sqrt(np.vdot(part, m@part).real)
    seed, data = pole(op, (omega, v))
    if (abs(seed[0]-omega)>1e-8 or abs(data['cavity_fraction']-row['cavity_fraction'])>1e-7
            or data['seed_overlap']<.99999):
        raise ArithmeticError('Failed to recover the archived degree-six mode')
    return seed, data


def integer_readout(op, seed, data):
    o, v = seed
    native = op.port.traction(o)
    angular_error = np.max(abs(native-AngularPort(op.port,op.degree).traction(o)))/np.max(abs(native))
    if not np.isfinite(angular_error) or angular_error>1e-10:
        raise ArithmeticError('Integer and continuous angular traction disagree')
    face = op.port.project@v[op.port.ids]
    old = .25*(abs(o)**2*np.vdot(v, op.M@v).real+np.vdot(v, op.K@v).real)
    stored = old-.5*PORT**3*np.vdot(face, face).real
    flux = .5*o.real*PORT**4*abs(face)**2*op.port.traction(o.real).imag
    channels = op.port.cutoff<o.real
    if stored<=0 or np.sum(flux[channels])<=0 or data['loss']<=0:
        raise ArithmeticError('Unresolved outgoing power or linewidth')
    return dict(**data, integer_traction_relative_error=angular_error,
                physical_power_per_energy=np.sum(flux[channels])/stored,
                physical_channel_power=(flux[channels]/stored).tolist(), open_channels=int(sum(channels)),
                nearest_depth_cutoff=float(op.port.cutoff[np.argmin(abs(op.port.cutoff-o.real))]))


def follow(left, right, seed):
    # Both material forms are affine functions of lambda=ell(ell+1) at fixed order
    # and shape. Only the interior angular term, the flat angular term and DtN vary.
    lo = left.degree
    delta_k = right.K-left.K
    if not np.array_equal(left.live, right.live) or (left.M-right.M).nnz:
        raise ArithmeticError('Angular continuation changed material coordinates or inertia')
    def at(d):
        if d==lo:
            return left
        if d==lo+1:
            return right
        op = copy.copy(left)
        fraction = (d*(d+1)-lo*(lo+1))/(2*lo+2)
        op.K = left.K+fraction*delta_k
        op.port = AngularPort(left.port, d)
        return op

    start = seed
    accepted = [float(lo)]
    h = .25
    overlap = 1.
    while accepted[-1] < lo+1:
        d = min(lo+1, accepted[-1]+h)
        try:
            new, read = pole(at(d), seed)
            if read['seed_overlap']<.985 or abs(new[0].real/seed[0].real-1)>.05:
                raise ArithmeticError('Angular shape overlap or frequency step unresolved')
        except (ArithmeticError, ValueError, RuntimeError) as err:
            h *= .5
            if h<2**-14:
                raise ArithmeticError(f'Forward correspondence unresolved at ell={d}: {err}') from err
            continue
        accepted.append(d)
        seed = new
        overlap = min(overlap, read['seed_overlap'])
        h = min(.25, h*1.5)
    reverse = seed
    for d in reversed(accepted[:-1]):
        reverse, control = pole(at(d), reverse)
        if control['seed_overlap']<.985:
            raise ArithmeticError(f'Reverse step lost correspondence at ell={d}')
        overlap = min(overlap, control['seed_overlap'])
    check = abs(np.vdot(start[1], left.M@reverse[1]))**2
    if abs(reverse[0]-start[0])>1e-8 or check<.99999:
        raise ArithmeticError('Reverse continuation missed preceding pole')
    emit(kind='angular_correspondence',degree=lo+1,accepted_degrees=accepted,minimum_step_overlap=overlap,
         return_overlap=check, return_pole_error=abs(reverse[0]-start[0]))
    return seed, read


def exit_work(trial, sh, c, op, center):
    for h in (.005, .0025):
        vals = []
        for sign in (-1, 1):
            shifts = sh.copy()
            shifts[1] += sign*h
            at = Sector(trial, shifts, c, op.degree)
            if not np.array_equal(at.live, op.live):
                raise ArithmeticError('Signed geometry displaced a wave label')
            sol, data = pole(at, center)
            if data['seed_overlap']<.98:
                raise ArithmeticError('Signed exit work crossed branches')
            vals.append((at.shape.exit, sol[0], data['seed_overlap']))
        (bm, m, a), (bp, p, b) = vals
        emit(kind='signed_exit', degree=op.degree, log_step=h, minus_work=m.real-center[0].real,
             plus_work=p.real-center[0].real, domega_dexit=(p.real-m.real)/(bp-bm),
             dloss_dexit=-2*(p.imag-m.imag)/(bp-bm), minimum_overlap=min(a,b))


def run(radius, n):
    trial = Trial(radius, n, bonded_core=True)
    key = f'r{radius}n{n}'
    base = ROOT/'_scratch/s11c/mouth_exit_release01'/key/'endpoint.npz'
    sh,c,_,digest,status,*_ = inputs(base, radius, n, trial)
    folder = ROOT/'_scratch/s11c/mouth_trapping_modes01'/(key+'l6')
    guard = json.loads((folder/'outcome.json').read_bytes())
    source = (folder/'stdout').read_bytes()
    rows = list(map(json.loads, source.splitlines()))
    if (status!='balanced' or rows[0]['sha256']!=digest or not guard['limitsVerified']
            or guard['exitCode'] or guard['childOutcome']['exitCode']):
        raise ValueError('Archived spectrum or shape failed verification')
    row, = [r for r in rows if r['kind']=='mode' and r['label']=='clamped_1']
    op = Sector(trial, sh, c, 6)
    seed, data = recorded_seed(op, row)
    emit(kind='start',radius=radius,mesh=n,degree=6,endpoint_sha256=digest,
         spectrum_stdout_sha256=hashlib.sha256(source).hexdigest(),**integer_readout(op,seed,data))
    for degree in range(7, 13):
        right = Sector(trial, sh, c, degree)
        try:
            seed, data = follow(op, right, seed)
        except (ArithmeticError, ValueError, RuntimeError) as exc:
            emit(kind='unresolved_angular_step',degree=degree,reason=str(exc))
            return
        op = right
        emit(kind='angular_mode',degree=degree,**integer_readout(op,seed,data))
        try:
            exit_work(trial,sh,c,op,seed)
        except (ArithmeticError, ValueError, RuntimeError) as exc:
            emit(kind='unresolved_exit_work',degree=degree,reason=str(exc))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',choices=(1,2),type=int,required=True)
    p.add_argument('--mesh',choices=(48,72),type=int,required=True)
    a=p.parse_args()
    run(a.radius,a.mesh)
