"""Absolute outgoing sound from the measured breathing sweep treated as a specified 4D mass source."""
import hashlib
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import spherical_jn

from mouth_exit_witness import inputs
from mouth_poly_shallow_release import ROOT, Trial, emit
from mouth_sealed_shape_relax import Meridian


def ball(x):
    return 3*spherical_jn(1, x)/x


def calibration(shape, basis, exit_scale, omega, quadrature):
    """Fourier surface sweep per unit mass-normalized generalized displacement."""
    points, weights = leggauss(quadrature)
    t = (points+1)/2
    wt = weights/2
    angle = (points+1)*np.pi/2
    # On the isotropic four-dimensional shell, dOmega/S3=(2/pi)sin(theta)^2 dtheta.
    angular = weights*np.sin(angle)**2
    k = np.sqrt(2*omega**2/(1+np.sqrt(1+omega**2)))
    kt = k*np.sin(angle)[:, None]
    kw = k*np.cos(angle)[:, None]
    r, drdw, _ = shape.profile(t)
    depth = shape.length
    exit_radius = shape.exit
    dl = depth*basis[0]
    db = basis[1]*exit_scale*t+4*shape.neck*basis[2]*t*(1-t)
    side_measure = 4*np.pi*r*r*depth*(db-drdw*t*dl)*wt
    phase = np.exp(-1j*kw*(depth*t))
    side = (phase*spherical_jn(0, kt*r)) @ side_measure
    disk = np.exp(-1j*kw[:,0]*depth)*ball(kt[:,0]*exit_radius)
    cap_area = 4*np.pi*exit_radius**3/3
    locked = side+dl*cap_area*disk
    s0 = float(side_measure.sum()+dl*cap_area)

    # Differentiate the Fourier transform of the four-volume as a separate
    # computation of the boundary transport formula over curved wall and end.
    def swept_volume(h):
        shift = shape.params.copy()
        shift[0] += h*basis[0]
        shift[1] += np.log1p(h*basis[1]*exit_scale/exit_radius)
        shift[2] += h*basis[2]
        moved = Meridian(shape.neck, shift)
        radii = moved.profile(t)[0]
        weight = (4*np.pi/3)*radii**3*moved.length*wt
        return (np.exp(-1j*kw*(moved.length*t))*ball(kt*radii)) @ weight

    errors = []
    for h in (1.e-4, 5.e-5):
        derivative = (swept_volume(h)-swept_volume(-h))/(2*h)
        errors.append(float(np.max(abs(derivative-locked))))

    group_speed = k*(1+k*k/2)/omega
    # From d<H>/dt = integral chemical_potential * s for a prescribed source
    # in the continuity equation. Fourier convention d^4k/(2pi)^4.
    rate = k**3*(1+k*k/4)/(32*np.pi*group_speed)
    # Alternate residue from (k^2-k0^2)(k^2+k0^2+4)/4 in density response.
    alternate_rate = omega*k**2*(4+k*k)/(64*np.pi*(2+k*k))
    # Mean cavity-reference energy=omega^2 |q|^2/2 at unit modal mass.
    # Net loading is evaluated on the prescribed displacement, with no backreaction.
    def source(name, strength):
        src = locked+strength*disk
        moment = float(angular @ abs(src)**2)
        per_velocity = rate*moment
        return dict(source=name, added_disk_volume_per_q=[strength.real,strength.imag],
                    source_squared_angular_mean=moment, power_per_squared_modal_speed=per_velocity,
                    work_per_cycle_over_reference_energy=4*np.pi*per_velocity/omega,
                    reference_omega_energy_over_power=omega/(2*per_velocity))

    optimum = -(angular @ (np.conj(disk)*locked))/(angular @ abs(disk)**2)
    cases = [source('full_geometric_sweep',0.+0j), source('cancel_integrated_volume',-s0+0j),
             source('least_power_uniform_disk', optimum)]
    # A point-monopole check from a 4D outgoing Helmholtz Green function.
    # G=i k H_1^(1)(kr)/(8pi r), |G|^2 r^3 -> k/(32 pi^3).
    # The density wave from a unit continuity source is
    # -i omega G/(1+k^2/2); multiply its mean energy density by c_g S3.
    radial_flux = (2*np.pi**2)*group_speed*.5*(1+k*k/4)*omega**2*k / (
        32*np.pi**3*(1+k*k/2)**2)
    return dict(k=k,group_speed=group_speed,point_power=rate,
                alternate_point_power=alternate_rate,point_power_by_far_flux=radial_flux,
                unsigned_volume_norm=abs(s0),side_volume=side_measure.sum(),cap_volume=cap_area*dl,
                transformed_volume_errors=errors,
                disk_squared_angular_mean=float(angular @ abs(disk)**2),cases=cases)


def main():
    record = ROOT/'research/pde_ledger_v3/_measurements/mouth_gallery_action.json'
    data = json.loads(record.read_bytes())
    run = data['runs']['mouth_action_breathing01/run']
    raw = run['literal_stdout']
    digest = hashlib.sha256(raw.encode()).hexdigest()
    if run['guard_exit'] or run['child_exit'] or digest != run['stdout_sha256']:
        raise ValueError('Invalid archived mode source')
    rows = list(map(json.loads,raw.splitlines()))
    for radius in (1,2):
        n = 72
        trial = Trial(radius,n,bonded_core=True)
        path = ROOT/f'_scratch/s11c/mouth_exit_release01/r{radius}n{n}/endpoint.npz'
        shifts, _, _, background_hash, status, *_ = inputs(path, radius, n, trial)
        origin = json.loads((path.parent/'stdout').read_bytes().splitlines()[0])
        selected = [r for r in rows if r['radius']==radius and r['n']==n]
        setup, = [r for r in selected if r['kind']=='setup']
        if setup['endpoint_sha256']!=background_hash or status!='balanced':
            raise ValueError('Different mouth preparation')
        shape = Meridian(radius, trial.x+shifts[:4])
        M = np.array(next(r['mass'] for r in selected if r['kind']=='basis_state' and r['relative_exit']==0))
        emit(kind='inputs',radius=radius,material_mesh=n,endpoint_sha256=background_hash,
             archived_modes_sha256=digest,depth=shape.length,exit=shape.exit,phasor='exp(-i omega t)')
        for row in selected:
            if row['kind']!='three_shape_trial':
                continue
            v = np.array(row['basis'])
            mass = v@M@v
            freq = row['freq']
            readings = [calibration(shape,v,origin['exit_radius'],freq,res) for res in (64,128)]
            lower, measured = readings
            step_error = max(abs(a['power_per_squared_modal_speed']-b['power_per_squared_modal_speed'])
                             for a,b in zip(lower['cases'],measured['cases']))
            emit(kind='radiation_readout',radius=radius,index=row['index'],frequency=freq,modal_mass=mass,
                 power_quadrature_error=step_error,**measured)
            if (abs(mass-1)>1e-8 or abs(measured['unsigned_volume_norm']-row['abs_dVdq_per_sqrt_mass'])>1e-7
                or abs(measured['disk_squared_angular_mean']-row['bulk_unit_disk_monopole_relative_power'])>1e-10
                or measured['transformed_volume_errors'][-1]>5e-7
                or step_error>1e-8
                or max(abs(measured[s]/measured['point_power']-1) for s in (
                    'point_power_by_far_flux','alternate_point_power'))>1e-12):
                raise ArithmeticError('Outgoing sound or shape calibration mismatch')


if __name__=='__main__':
    main()
