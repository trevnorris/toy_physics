"""Measure the apparent normal reflection plane of a held half-order sheet.

The atom and triangle-center rays are invariant under 120-degree rotations.
An incoming point tag therefore has a specular return on either ray. Compare
the flight time from the conserved energy with the existing Cartesian
trajectory engine, without imposing either a hard wall or a contact detector.
Heights are in native lengths above the plane of the host nuclei.
"""
import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from mouth_attachment_impact_trial import sha
from mouth_attachment_incidence import HalfOrderBand, incident_normal_speed
from mouth_surface_site_transfer import sheet_mean
from neutral_patch_periodic_compression import BASIS
from neutral_patch_sites import BETWEEN, TEMPERATURE, emit, morse


QUANTILES = (.1, .5, .9)
HEIGHTS = (5.5, 7.)
REGISTRIES = {'atom': (0., 0.), 'hollow': (1/3, 1/3)}
SWITCH_START, CUTOFF = 4.5, 5.


class Ray:
    def __init__(self, band, registry):
        self.pitch = band.pitch
        self.chi = band.scalar_order
        # Each reciprocal basis coordinate varies by <= cutoff/(pitch sin60)
        # within the contact disk. Include every site able to enter that disk.
        extent = CUTOFF/(self.pitch*BASIS[1, 1])
        uv = np.asarray(registry)
        axes = [np.arange(np.floor(c-extent), np.ceil(c+extent)+1) for c in uv]
        tiles = np.array(np.meshgrid(*axes, indexing='ij')).reshape(2, -1).T
        xy = self.pitch*((tiles-uv) @ BASIS)
        r2 = np.einsum('ij,ij->i', xy, xy)
        self.r2 = r2[r2 < CUTOFF**2]

    def energy(self, z):
        r = np.sqrt(self.r2 + z*z)
        r = r[r < CUTOFF]
        s = np.clip((r-SWITCH_START)/(CUTOFF-SWITCH_START), 0., 1.)
        switch = (1-s)**3*(1+3*s+6*s*s)
        discrete = switch @ morse(r, *BETWEEN)[0]
        uniform = sheet_mean(z, self.pitch)[0]
        return self.chi*discrete+(1-self.chi)*uniform

    def flight(self, energy, height):
        speed = np.sqrt(2*(energy-self.energy(height)))
        turn = brentq(lambda z: self.energy(z)-energy, .3, BETWEEN[2], xtol=1.e-15)
        # z = turn + s² removes the square-root singularity in the flight
        # time at the turning point; the outward leg retraces the inward leg.
        def time_element(s):
            return 2*s/np.sqrt(2*(energy-self.energy(turn+s*s)))

        half_time, error = quad(time_element, 0., np.sqrt(height-turn),
                                epsabs=1.e-9, epsrel=1.e-11)
        return dict(launch_energy=energy, speed=speed, read_height=height,
                    energy_at_height=self.energy(height), z_turn=turn,
                    return_time=2*half_time, z_tof=height-speed*half_time,
                    z_tof_quad_error_estimate=speed*error)


def cartesian(args, band, energy, height, samples):
    # Integrate the original unconstrained three-dimensional rigid-host engine
    # through the return crossing. The observation span is in physical time.
    observation = max(r['return_time'] for r in samples)+2.
    dt = args.dt
    cmd = [str(args.engine.resolve()), '--host', 'rigid', '--mixed', '1', '--direct', '1',
           '--pitch', str(band.pitch), '--side', '5', '--dt', str(dt),
           '--steps', str(int(np.ceil(observation/dt))), '--input', '-']
    rows = []
    for i, (registry, r) in enumerate(zip(REGISTRIES.values(), samples, strict=True)):
        rows.append(' '.join(map(str, (i, *registry, height, 0., 0., -r['speed'], band.scalar_order))))
    run = subprocess.run(cmd, input='\n'.join(rows)+'\n', text=True, check=True,
                         stdout=subprocess.PIPE)
    trajectories = [json.loads(line) for line in run.stdout.splitlines()]
    for i, (r, t) in enumerate(zip(samples, trajectories, strict=True)):
        if t['id'] != i or t['escape_time'] <= 0. or t['plane_crossing_time'] >= 0.:
            raise ArithmeticError('Unresolved return to read plane')
        native_height = height-r['speed']*t['escape_time']/2
        r.update(cartesian_z_tof=native_height, cartesian_z_min=t['min_z'],
                 cartesian_min_energy=t['min_energy_tag'],
                 cartesian_max_energy_defect=t['max_energy_error'],
                 cartesian_lateral_end_displacement=np.hypot(t['dx_end'], t['dy_end']),
                 cartesian_outward_speed_defect=t['vz_end']-r['speed'])
        if abs(t['energy_tag_initial']-energy) > 1.e-14:
            raise ArithmeticError('Launch total energy differs from requested energy')


def run(args):
    band = HalfOrderBand.load(args.section)
    rays = {name: Ray(band, uv) for name, uv in REGISTRIES.items()}
    emit(test='reflection_setup', section=str(args.section), section_sha256=sha(args.section),
         engine=str(args.engine), engine_sha256=sha(args.engine), model='held sheet, retained mean',
         host_site_plane_height=0., pitch=band.pitch, chi=band.scalar_order, well=BETWEEN,
         taper=[SWITCH_START, CUTOFF], temperature=TEMPERATURE, read_heights=HEIGHTS,
         normal_beam_flux_quantiles=QUANTILES, primitive_registries=REGISTRIES,
         cartesian_dt=args.dt, native_mass=1.)
    records = []
    for q in QUANTILES:
        energy = incident_normal_speed(q)**2/2
        for height in HEIGHTS:
            rows = [dict(ray.flight(energy, height), registry=name, quantile=q)
                    for name, ray in rays.items()]
            cartesian(args, band, energy, height, rows)
            for row in rows:
                emit(test='reflection', **row)
            records.extend(rows)
    read_height_delta = max(abs(a['z_tof']-b['z_tof']) for a in records for b in records
                            if a['quantile']==b['quantile'] and a['registry']==b['registry'])
    at_first_height = [r for r in records if r['read_height']==HEIGHTS[0]]
    summary = dict(
        z_tof_atom_minus_hollow=[a['z_tof']-b['z_tof'] for a, b in
                                 zip(at_first_height[::2], at_first_height[1::2], strict=True)],
        z_tof_range=[min(r['z_tof'] for r in records), max(r['z_tof'] for r in records)],
        z_turn_range=[min(r['z_turn'] for r in records), max(r['z_turn'] for r in records)],
        max_abs_read_height_difference=read_height_delta,
        max_abs_cartesian_tof_difference=max(abs(r['cartesian_z_tof']-r['z_tof']) for r in records),
        max_abs_cartesian_turn_difference=max(abs(r['cartesian_z_min']-r['z_turn']) for r in records),
        max_absolute_tail_energy=max(abs(r['energy_at_height']) for r in records),
        max_lateral_displacement=max(r['cartesian_lateral_end_displacement'] for r in records),
        max_abs_outward_speed_defect=max(abs(r['cartesian_outward_speed_defect']) for r in records),
        min_cartesian_tag_energy=min(r['cartesian_min_energy'] for r in records))
    emit(test='reflection_summary', **summary)
    if (read_height_delta > 1.e-8 or summary['max_abs_cartesian_tof_difference'] > 5.e-5
            or summary['max_abs_cartesian_turn_difference'] > 5.e-5
            or summary['max_lateral_displacement'] > 1.e-8
            or summary['max_abs_outward_speed_defect'] > 1.e-7
            or summary['min_cartesian_tag_energy'] <= 0.):
        raise ArithmeticError('Unresolved reflection geometry or Cartesian accuracy')


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--section', type=Path, required=True)
    p.add_argument('--engine', type=Path, required=True)
    p.add_argument('--dt', type=float, default=.001)
    args = p.parse_args()
    if not 0 < args.dt < np.inf:
        p.error('dt must be finite and positive')
    run(args)
