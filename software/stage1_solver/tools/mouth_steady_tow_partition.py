"""Compare separate stationary-stream drives to archived physical-chart towing.

Use the same rigid, sealed-face/open-exit coherent carrier as the drag audit.
Report medium-on-mouth forces for positive stationary reservoir velocities and
the external actuator force for positive lab tow velocity. All readings are
linear-response coefficients at zero frequency and a prescribed intake.
"""
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_carrier_lateral_inertia import WINDOWS, volume
from mouth_carrier_lateral_response import Dipole, surfaces, trace
from mouth_flared_receiver import EXIT, EXIT_RADIUS
from mouth_occupied_push_baths import digest, load_bath


DATA = Path('research/pde_ledger_v3/_measurements')
BASE = Path('_scratch/s11c/occupied_push_bath18_n34_01')


def bulk_shell(dom, fields, radius):
    """Outgoing x momentum on a hemisphere about the exit and its flat shelf.

    Split the quadrature at every meridional mesh edge. The open outlet disk
    closes the receiving control volume. The shelf outside it is impermeable.
    """
    angles = [0., np.pi/2]
    for block in (3, 4):
        pg, vg, _ = dom.blocks[block]
        radial = dom.chart(block, np.full_like(vg, EXIT), vg)[0]
        angles.extend(np.arcsin(r/radius) for r in radial if 0 < r < radius)
        angles.extend(np.arccos((w-EXIT)/radius) for w in pg if EXIT < w < EXIT+radius)
    edges = np.unique(angles)
    root, weights = leggauss(8)
    half = np.diff(edges)[:, None]/2
    phi = (edges[:-1, None]+half*(root+1)).ravel()
    dphi = (half*weights).ravel()
    r, w = radius*np.sin(phi), EXIT+radius*np.cos(phi)
    value = 0.
    for block in (3, 4):
        mask = r < EXIT_RADIUS if block == 3 else r >= EXIT_RADIUS
        rr, ww = r[mask], w[mask]
        v = 1-rr/EXIT_RADIUS if block == 3 else (rr-EXIT_RADIUS)/(dom.far-EXIT_RADIUS)
        *_, f, fr, fw = trace(dom, block, ww, v, fields)
        lag = (fr[:, 0].conj()*fr[:, 1]+fw[:, 0].conj()*fw[:, 1]).real
        lag += (abs(f[:, 0])**8-1)/2*(f[:, 0].conj()*f[:, 1]).real
        dx = fr[:, 1]+2*f[:, 1]/rr
        rad = (fr[:, 0].conj()*fr[:, 1]+dx.conj()*fr[:, 0]).real-lag
        axial = (fr[:, 0].conj()*fw[:, 1]+dx.conj()*fw[:, 0]).real
        value += (4*np.pi/3)*np.sum(rr**2*radius*dphi[mask]*
                                     (rad*np.sin(phi[mask])+axial*np.cos(phi[mask])))
    from mouth_carrier_lateral_response import surface
    shelf, _ = surface(dom, fields, 4, 'v', EXIT, (0., -1.),
                       (0., (radius-EXIT_RADIUS)/(dom.far-EXIT_RADIUS)))
    disk, _ = surface(dom, fields, 3, 'v', EXIT, (0., -1.))
    return dict(hemisphere=value, shelf_raw=shelf, disk_outward=disk,
                closure=value+shelf+disk)


def measure(dom, f, h):
    fields = np.column_stack([f, h])
    _, ports, neck = surfaces(dom, fields)
    windows = [float(volume(dom, fields, win, False)['virtual_wall_force']) for win in WINDOWS]
    return dict(medium_on_rigid=-sum(ports.values()), outward=ports, neck_plus_w=neck,
                interior=windows)


def main():
    dom, f = load_bath(BASE)
    archive = json.loads((DATA/'mouth_rigid_drag_checks.json').read_text())['runs']['steady']
    log = Path(archive['command'][archive['command'].index('--log-directory')+1])
    outcome = json.loads((log/'outcome.json').read_text())
    assert outcome['limitsVerified'] and outcome['exitCode'] == outcome['childOutcome']['exitCode'] == 0
    assert digest(log/'stdout') == archive['stdout_sha256']
    meta, steady = map(json.loads, archive['literal_stdout'].splitlines())
    assert meta['sha256'] == digest(BASE/'field.npz') and steady['omega'] == 0
    assert (log/'stdout').read_text() == archive['literal_stdout']
    model = Dipole(dom, f)
    amp = .001
    result = {}
    for label, drives in [('brane', (amp, 0.)), ('bulk', (0., amp)), ('common_minus', (-amp, -amp))]:
        h, err = model.steady(*drives)
        result[label] = measure(dom, f, h/amp)
        result[label]['residual'] = err
        if label == 'brane':
            fields = np.column_stack([f, h/amp])
            result[label]['receiving_shells'] = {r: bulk_shell(dom, fields, r) for r in (3., 6., 12.)}
    # The phase boost tests the reference-frame ledger. Its node interpolation
    # is an additional field readout, not an independent force measurement.
    result['nodal_lab_boost'] = measure(dom, f, 1j*dom.points[:, 0]*f)
    p = sum(result[k]['medium_on_rigid'] for k in ('brane', 'bulk'))
    intp = np.sum([result[k]['interior'] for k in ('brane','bulk')], axis=0)
    actual = steady['FoverV_wall_via_domain'][0]
    result['comparison'] = dict(positive_lab_tow_external=actual, common_flow_prediction=p,
                                common_dirichlet_linear_error=p+result['common_minus']['medium_on_rigid'],
                                separate_frame_difference=actual-p, window_prediction=intp.tolist(),
                                archived_tow_windows=[row[0] for row in steady['FoverV_two_windows']],
                                archived_tow_ports=steady['FoverV_ports'])
    if abs(result['comparison']['common_dirichlet_linear_error']) > 1e-9:
        raise ValueError('The common-stream and two separate-stream solves differ')
    print(json.dumps(dict(model='prescribed cubic rigid R=1; native box (18,20); intake 0.25',
                          base_sha256=meta['sha256'], tow_stdout_sha256=archive['stdout_sha256'],
                          drive_amplitude=amp, results=result), allow_nan=False))


if __name__ == '__main__':
    main()
