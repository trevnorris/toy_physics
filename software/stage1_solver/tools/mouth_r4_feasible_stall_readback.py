"""Measure the saved R4 stall state without rematching fields or advancing a coordinate."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_exact_checkpoint_control import acquired, digest, read


def flow_readings(row, depth):
    fractions = row['shaft_order']['terminal_from_fraction']
    return dict(feed=row['feed_current'], bulk=row['bulk_current'], weak_disk=row['weak_disk_current'],
                incident_cut=row['incident_slice_current'], receiving_cut=row['receiving_slice_current'],
                feed_bulk_mismatch=row['port_error'], disk_minus_feed=row['disk_minus_feed_current'],
                dense_subsonic_share=row['dense_total_subsonic_share'],
                full_disk_zero_from=depth*fractions['disk']['0.0'],
                axis_zero_from=depth*fractions['axis']['0.0'], shaft_order=row['shaft_order'])


def measure(area):
    source = area.parent/'r4_melted_survival01'
    previous = area.parent/'r4_exact_support01'
    replay = area.parent/'r4_exact_fulltrial01'
    case = area/'coarse'
    st = read(case/'stall.json')
    guard = read(area/'guard/outcome.json')
    error = read(case/'issue.json')
    assert guard['limitsVerified'] and guard['childOutcome']['guardReason'] is None
    assert guard['exitCode'] == guard['childOutcome']['exitCode'] == 1
    assert error['error'] == 'MaterialStall' and not (area/'DONE.json').exists()
    assert st['last_accepted']['iteration'] == 10 and st['attempted_checkpoint'] == 11
    assert not list(case.glob('state_*.json')) and not (case/'endpoint.json').exists()
    assert all(digest(area/key)==sha for key,sha in read(area/'sources.json').items())
    assert all(digest(previous/key)==sha for key,sha in read(previous/'inputs.json').items())
    gate = read(case/'full_trial_gate.json')
    path = case/st['stall_state']
    assert digest(path)==st['sha256'] and digest(replay/'trial.npz')==gate['replayed_sha']
    build, _, _, _, _, cp = acquired(source, previous/'patch')
    with np.load(path,allow_pickle=False) as f, np.load(replay/'trial.npz',allow_pickle=False) as start:
        assert int(f['radius'])==4 and int(f['winding'])==1
        assert str(f['state_sha'])==str(start['state_sha']) and np.array_equal(f['x'],start['x'])
        x, q, u = f['x'], f['q'], f['field']
        trial = build(x)
        assert np.array_equal(trial.b.unpack(q),f['solid'])
        assert np.array_equal(trial.transported(u),u)
    cert = trial.check_material_domain(q,u)
    if cert['status']!='certified':
        raise ArithmeticError(f'Stalled placement failed exact support check: {cert}')
    load = trial.crystal_load_for_match(q,u)+trial.ext
    flow = trial.measure(u,load)
    trial.model.material_load = load
    energy, _ = trial.model.evaluate(u)
    solid_energy, gq, _, _ = trial.b.evaluate(q,u[2],load=False)
    residual = trial.b.residual(q,gq)
    authority = st['material_blocks'][-1]['field_update']
    assert energy==authority['energy_before'] and authority['delta_total_energy']==0.
    assert float(u[2,trial.first[51603]])==authority['order51603_after']
    assert flow['port_error']<3.e-5 and abs(flow['disk_minus_feed_current'])<3.e-5
    return dict(stalled_snapshot_sha256=st['sha256'],last_accepted_snapshot_sha256=cp['sha256'],
                status='Fixed-shape diagnostic after hard stop, preceding checkpoint 11 signed probes.',
                radius=4,winding=1,phase_drop=0.4510084282857874,exit=trial.shape.exit,
                depth=trial.shape.length,energy=float(energy),solid_energy=float(solid_energy),
                material_residual=residual,certificate=cert,order51603=float(u[2,trial.first[51603]]),
                flow=flow_readings(flow,trial.shape.length),
                cp10=dict(exit=cp['exit_radius'],depth=cp['depth'],shape_force=cp['shape_force'],
                          flow=flow_readings(cp['current'],cp['depth'])),
                cleared_shape=dict(exit=gate['exit'],depth=gate['depth'],flow=flow_readings(gate['current'],gate['depth'])))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('run',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    with args.output.open('x') as f:
        json.dump(measure(args.run.resolve()),f,indent=2,allow_nan=False)
        f.write('\n')
