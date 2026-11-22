"""Release two warm patches into a common pair law and test shear and projected compression."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_brace_grip import Material, clock
from mouth_bulk_contact import Tracking, converted, material_pair, observations
from mouth_middle_tag_grip import cubic_trace
from neutral_patch_open_w import mix
from neutral_patch_sites import BETWEEN, TEMPERATURE, WITHIN, emit, morse
from neutral_patch_structural_endpoints import PITCH


REFERENCE = Path('research/pde_ledger_v3/_measurements/mouth_bulk_contact.json')
FRAMES = 48
STEP = .0001


def reference_rows():
    run = json.loads(REFERENCE.read_text())['runs']['p1']
    if not run['limits_verified'] or run['exit_code'] or run['child_exit']:
        raise ValueError('Missing verified reference trajectory')
    return {(row['test'], round(FRAMES*row['fraction'])): row
            for row in map(json.loads, run['literal_stdout'].splitlines())
            if row.get('contact') == 'joined'}


def coordinates(material, pos):
    cell = material.cell
    d = pos[cell.j] - pos[cell.i]
    d[:, :2] += cell.shifts - cell.box*np.rint(d[:, :2]/cell.box)
    return d


def pair_values(material, displacements):
    radii = np.linalg.norm(displacements, axis=1)
    keep = radii < material.cutoff
    d, r = displacements[keep], radii[keep]
    ea, da = morse(r, *material.wells[0])
    eb, db = morse(r, *material.wells[1])
    f, df, _ = mix(d[:, 2], material.power)
    energy = ea + f*(eb-ea)
    gradient = d*((da+f*(db-da))/r)[:, None]
    gradient[:, 2] += df*(eb-ea)
    return keep, energy, gradient


def measure(material, pos, check=False):
    cell = material.cell
    displacements = coordinates(material, pos)
    keep, energy, gradient = pair_values(material, displacements)
    virial = np.sum(gradient[:, :2]*displacements[keep, :2], axis=0)
    area = np.prod(cell.box)
    stress = (virial[0]-virial[1])/area
    pressure = (material.n*material.temperature-.5*sum(virial))/area
    if check:
        directions = np.array([[1., -1., 0.], [.5, .5, 0.]])
        for direction, expected in zip(directions, [stress, material.n*material.temperature/area-pressure]):
            delta = 2.e-7
            ep, em = [sum(pair_values(material, displacements*np.exp(sign*delta*direction))[1])
                      for sign in (1, -1)]
            np.testing.assert_allclose((ep-em)/(2*delta*area), expected, rtol=4e-6, atol=5e-6)
        # The native integrator omits self-image forces (equal and opposite).
        i, j = cell.i[keep], cell.j[keep]
        force = np.array([np.bincount(i, g, minlength=cell.n)-np.bincount(j, g, minlength=cell.n)
                          for g in gradient.T]).T
        probe = pos.copy()
        dt = 1e-6
        material.step(probe, dt, np.zeros(1), np.zeros((1, material.n, 3)))
        np.testing.assert_allclose((probe-pos)/dt, material.mobility*force, rtol=1e-8, atol=4e-9)
    return float(energy.sum()), float(stress), float(pressure)


def prepare(args, base):
    depths, times, order, intake = clock(34)
    trace, _ = cubic_trace(34)
    transit = trace.transit(depths[-1], trace.edges[-1], 8)
    bottom = converted(base, 319, times, order)
    top = converted(base, 817, times, order)
    top[:, 2] += bottom[:, 2].max()-top[:, 2].min()+PITCH
    top[:, :2] += PITCH*np.array([.5, np.sqrt(3)/6])
    pos = np.vstack([bottom, top])
    pos[:, 2] -= pos[:, 2].mean()
    mat = material_pair(base, True)
    shear, normal = observations(mat.cell, pos)
    ref = reference_rows()['preparation', 0]
    restored = [*shear, normal, np.ptp(pos[:, 2])]
    np.testing.assert_allclose(restored, [*ref['rest_stress'], ref['bond_normal_traction'],
                                          ref['height_span']], rtol=0., atol=1e-10)
    original_energy = measure(mat, pos, check=True)[0]
    mat.wells = [np.array(WITHIN), np.array(WITHIN)]
    within_energy = measure(mat, pos, check=True)[0]
    matched_depth = WITHIN[0]*original_energy/within_energy
    if matched_depth <= 0:
        raise ArithmeticError('A positive Morse depth cannot match the starting energy')
    with args.state.open('xb') as file:
        np.savez(file, pos=pos, transit=transit, matched_depth=matched_depth)
    emit(test='prepare', seeds=[319, 817], registration=1, mix_power=1,
         sites_per_patch=base.n, area=mat.cell.area, T=base.temperature,
         external_normal_traction=0., conversion_intake=intake, transit=transit,
         starting_stresses_normal_traction_height=restored,
         directional_energy=original_energy, within_energy=within_energy,
         matched_isotropic_depth=matched_depth, state_sha256=hashlib.sha256(args.state.read_bytes()).hexdigest())


class Deformed:
    def __init__(self, mat, pos, shear=0., log_area=0.):
        self.mat = mat
        stretch = np.exp(.5*log_area + np.array([shear, -shear]))
        self.box = mat.cell.box*stretch
        self.images = mat.shifts*stretch
        self.energy_images = mat.cell.shifts*stretch
        self.pos = pos.copy()
        self.pos[:, :2] *= stretch

    def select(self):
        self.mat.cell.box = self.box
        self.mat.cell.shifts = self.energy_images
        self.mat.shifts = self.images


def run_law(args, base):
    with np.load(args.state) as start:
        initial = start['pos']
        transit = float(start['transit'])
        depth = float(start['matched_depth'])
    mat = material_pair(base, True)
    if args.mode != 'directional':
        well = np.array(WITHIN)
        if args.mode == 'energy-matched':
            well[0] = depth
        mat.wells = [well, well]
    count = int(np.ceil(transit/(FRAMES*.001)))
    dt = transit/(FRAMES*count)
    rng = np.random.default_rng(319+23*817)
    ref = reference_rows()

    def interval(states, tracking):
        for start in range(0, count, 128):
            steps = min(128, count-start)
            noise = rng.normal(size=(steps, mat.n, 3))*np.sqrt(2*mat.temperature*mat.mobility*dt)
            order = np.zeros(steps)
            for state in states:
                state.select()
                mat.step(state.pos, dt, order, noise)
            tracking.update(states[0].pos, states[0].box)

    start = Deformed(mat, initial)
    tracker = Tracking(start.pos, start.box, True)
    emit(test='law', law=args.mode, initial_sha256=hashlib.sha256(args.state.read_bytes()).hexdigest(),
         wells=[w.tolist() for w in mat.wells], initial_observables=measure(mat, initial, True),
         dt=dt, noise_seed=319+23*817, transit=transit, site_mobility=mat.mobility,
         temperature=mat.temperature, area=float(np.prod(start.box)), normal_load=0., brace_order=0.,
         shear_source='archived original at identical state and noise' if args.mode=='directional' else 'measured')
    for _ in range(FRAMES):
        interval([start], tracker)
    start.select()
    emit(test='unstrained_contact', law=args.mode, observables=measure(mat, start.pos, True),
         **tracker.update(start.pos, start.box))
    if args.mode == 'directional':
        whole, normal = observations(mat.cell, start.pos)
        old = ref['preparation', FRAMES]
        np.testing.assert_allclose([*whole, normal, np.ptp(start.pos[:, 2])],
                                   [*old['rest_stress'], old['bond_normal_traction'], old['height_span']],
                                   atol=1e-10, rtol=0.)

    amps = [STEP, -STEP, STEP/2, -STEP/2]
    area_states = [Deformed(mat, start.pos, log_area=a) for a in amps]
    shear_states = [] if args.mode=='directional' else [Deformed(mat, start.pos, shear=a) for a in amps]
    states = [start, *area_states, *shear_states]
    tracker = Tracking(start.pos, start.box, True)
    samples = []
    for frame in range(FRAMES+1):
        if frame:
            interval(states, tracker)
        readings = []
        for state in states:
            state.select()
            readings.append(measure(mat, state.pos, frame in (0, FRAMES)))
        pressure = np.array([r[2] for r in readings[1:5]])
        ka = (pressure[1]-pressure[0])/(2*STEP)
        ka_half = (pressure[3]-pressure[2])/STEP
        if args.mode == 'directional':
            old = ref['hold', frame]
            np.testing.assert_allclose([readings[0][1], np.ptp(start.pos[:, 2])],
                                       [old['rest_stress'][0], old['height_span']], atol=1e-10, rtol=0.)
            g, g_half = old['full_derivative'][0], old['half_derivative'][0]
        else:
            sigma = np.array([r[1] for r in readings[5:]])
            g = (sigma[0]-sigma[1])/(2*STEP)
            g_half = (sigma[2]-sigma[3])/STEP
        diag = tracker.update(start.pos, start.box)
        emit(test='hold', law=args.mode, frame=frame, unstrained_energy_stress_pressure=readings[0],
             area_modulus=ka, area_modulus_half=ka_half, shear=g, shear_half=g_half, **diag)
        samples.append([ka, ka_half, g, g_half])
    values = np.array(samples)
    late = values[-12:].mean(axis=0)
    emit(test='late', law=args.mode, frames=[37, 48], area_start=values[0,0], area_late=late[0],
         area_late_half=late[1], shear_start=values[0,2], shear_late=late[2], shear_late_half=late[3],
         area_retained=late[0]/values[0,0], shear_retained=late[2]/values[0,2],
         max_full_half_difference=np.abs(values[:, ::2]-values[:, 1::2]).max(axis=0).tolist(),
         **diag)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('prepare', 'directional', 'within', 'energy-matched'))
    parser.add_argument('--work', required=True, type=Path)
    parser.add_argument('--state', required=True, type=Path)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    args.side, args.cutoff, args.power, args.brace = 4, 4.5, 1, 'elastic'
    args.mobility, args.temperature = 1., TEMPERATURE
    base = Material(args)
    if args.mode == 'prepare':
        prepare(args, base)
    else:
        run_law(args, base)
