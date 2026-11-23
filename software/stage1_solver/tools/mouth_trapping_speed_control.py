"""Follow the longest-lived sampled connected cavity resonance to unit local shear speed."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_exit_witness import inputs
from mouth_poly_shallow_release import Trial, emit
from mouth_strained_wave import cell_tangents
from mouth_trapping_modes import Sector, closed_seeds, form, pole
from types import SimpleNamespace


def read_rows(path):
    guard = json.loads((path/'outcome.json').read_bytes())
    if guard['exitCode'] or guard['childOutcome']['exitCode'] or not guard['limitsVerified']:
        raise ValueError(f'Input guard failed at {path}')
    return [json.loads(line) for line in (path/'stdout').read_text().splitlines()]


def stiffness_difference(op, data):
    mesh = op.mesh
    nt = len(op.ids)//mesh
    cell = data['cell']
    i, j = np.divmod(cell, nt)
    r, x, t = data['radius'], data['x'], data['t']
    width, slope, _ = op.shape.profile(t)
    u, z = x*mesh-i, t*nt-j
    n = np.stack(((1-u)*(1-z), u*(1-z), (1-u)*z, u*z), axis=1)
    nr = mesh*np.stack((z-1, 1-z, -z, z), axis=1)/width[:, None]
    nw = nt*np.stack((u-1, -u, 1-u, u), axis=1)/op.shape.length - (x*slope)[:, None]*nr
    measure = data['current_volume_per_solid_angle']/data['jacobian']
    if not np.array_equal(data['wave_cell_chi'], op.weights[cell]):
        raise ArithmeticError('The speed map and spectrum use different shear order')
    def outer(a,b):
        return a[:, :, None]*b[:, None, :]
    def integrate(q):
        result = np.zeros_like(op.mc)
        np.add.at(result, cell, measure[:, None, None]*q)
        return result
    grad_r, grad_w, amp = outer(nr,nr), outer(nw,nw), outer(n,n)
    angular = op.degree*(op.degree+1)-2
    r2 = r[:, None, None]**2
    flat = integrate(r2*(grad_r+grad_w)+angular*amp)
    delta = data['c_meridional_min']**2 - data['c_meridional_max']**2
    brw = delta*data['slow_meridional_nr']*data['slow_meridional_nw']
    measured = integrate(r2*(data['c_radial'][:, None, None]**2*grad_r
                             + data['c_depth'][:, None, None]**2*grad_w
                             + brw[:, None, None]*(outer(nr,nw)+outer(nw,nr)))
                         + angular*data['c_hoop'][:, None, None]**2*amp)
    mass = integrate(r2*amp)
    dummy = SimpleNamespace(mesh=mesh, ids=op.ids, kc=op.mc, mc=op.mc)
    native, _ = cell_tangents(dummy, op.shape, op.coeff, degree=op.degree, core_reference=op.reference)
    err_k = np.max(abs(measured-native))/np.max(abs(native))
    err_m = np.max(abs(mass-op.mc))/np.max(abs(op.mc))
    if max(err_k,err_m)>2e-12:
        raise ArithmeticError(f'Change of tangent also changed its native controls: {err_k}, {err_m}')
    delta = form(op.ids, op.weights[:, None, None]*(flat-native), len(op.live))
    return delta[op.live][:, op.live].tocsc(), float(err_k), float(err_m)


def run(radius):
    case = f'r{radius}n72'
    scratch = Path('_scratch/s11c')
    spectrum = scratch/'mouth_trapping_modes01'
    modes = []
    def scan(degree):
        # One guarded scan was admitted in the replacement directory.
        candidates = [p for p in (spectrum/f'{case}l{degree}', spectrum/f'{case}l{degree}_adopt')
                      if (p/'outcome.json').is_file()]
        if len(candidates)!=1:
            raise ValueError(f'Expected one completed scan for {case}, degree {degree}')
        return read_rows(candidates[0])
    for degree in range(1,7):
        rows = scan(degree)
        modes.extend((s['loss'], degree, s) for s in rows if s['kind']=='mode' and
                     not s['detached'] and s['cavity_fraction']>.5 and s['loss']>0)
    _, degree, longest = min(modes, key=lambda v:v[0])
    tr = Trial(radius, 72, bonded_core=True)
    path = scratch/'mouth_exit_release01'/case/'endpoint.npz'
    shifts, chi, _, digest, status, *_ = inputs(path, radius, 72, tr)
    if status != 'balanced':
        raise ValueError('Different material state')
    from_modes = scan(degree)[0]
    if from_modes['sha256'] != digest:
        raise ValueError('Spectrum fingerprint differs from the target material state')
    maps = scratch/'mouth_trapping_speeds01'
    row, = [s for s in read_rows(maps/'measured') if s['kind']=='map' and s['radius']==radius
            and s['order_mesh']==72]
    grid = Path(row['path'])
    if row['endpoint_sha256']!=digest or row['sha256']!=hashlib.sha256(grid.read_bytes()).hexdigest():
        raise ValueError('Speed map fingerprint differs from the target material state')
    op = Sector(tr, shifts, chi, degree)
    initial = closed_seeds(op,3)[int(longest['label'].split('_')[1])-1]
    origin, info = pole(op, initial)
    if (abs(info['loss']-longest['loss']) > 1e-8*longest['loss']
            or abs(info['frequency']-longest['frequency']) > 1e-8):
        raise ArithmeticError('Could not reconstruct the selected center pole')
    with np.load(grid,allow_pickle=False) as data:
        correction, err_k, err_m = stiffness_difference(op,data)
    original = op.K
    emit(kind='material',radius=radius,degree=degree,center_sha256=digest,
         speed_map_sha256=row['sha256'],mode=longest['label'],native_K_error=err_k,
         native_M_error=err_m,unchanged_mass=True,**info)
    seed = origin
    for leg, target, initial_width in (('flatten',1.,.05),('restore',0.,.025)):
        s = 1.-target
        width = initial_width
        count=0
        minimum = 1.
        while s != target:
            nxt = s+np.sign(target-s)*min(width,abs(target-s))
            if abs(nxt-target)<1e-12:
                nxt = target
            op.K = original+nxt*correction
            try:
                trial, stats = pole(op,seed)
                if stats['seed_overlap']<.985 or abs(trial[0].real/seed[0].real-1.)>.025:
                    raise ArithmeticError('Mode correspondence unresolved on proposed step')
            except (ArithmeticError,RuntimeError,ValueError) as exc:
                width *= .5
                if width<1e-6:
                    emit(kind='unresolved',radius=radius,leg=leg,last_accepted=s,
                         proposed=nxt,reason=str(exc))
                    return
                continue
            s, seed = nxt, trial
            minimum = min(minimum,stats['seed_overlap'])
            count += 1
            emit(kind='continued',radius=radius,leg=leg,fraction_flat=s,**stats)
            width = min(initial_width, width*1.4)
        emit(kind='endpoint',radius=radius,leg=leg,steps=count,min_step_overlap=minimum,
             original_mass_overlap=float(abs(np.vdot(origin[1],op.M@seed[1]))**2),
             freq_change=seed[0].real-origin[0].real,**stats)
        if leg=='restore' and (abs(seed[0]-origin[0])>1e-8 or
                               abs(np.vdot(origin[1],op.M@seed[1]))**2<.99999):
            raise ArithmeticError('Reverse continuation returned a different resonance')


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',type=int,choices=(1,2),required=True)
    run(p.parse_args().radius)
