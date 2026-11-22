"""Test a necessary stationary power balance before fitting free mouth shapes.

Use the two calibration radii only. Number-current units equal mass-current
units here. A phase derivative of stationary energy measures the current; a
time derivative of the port phase determines its power.
"""
import hashlib
import json
from pathlib import Path

from mouth_radius_pivot_summary import read_run

RECORDS = Path('research/pde_ledger_v3/_measurements')
BIAS = 0.4510084282857874


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def archived_conversion(manifest, label):
    record = manifest['runs'][label]
    data = record['literal_stdout'].encode()
    if (hashlib.sha256(data).hexdigest() != record['stdout_sha256']
            or record['guard_exit'] or record['worker_exit'] or not record['limits_verified']):
        raise ValueError(f'Archived conversion control failed verification: {label}')
    return [json.loads(line) for line in data.splitlines()]


def main():
    radius_file = RECORDS / 'mouth_radius_pivot.json'
    conversion_file = RECORDS / 'mouth_coupled_conversion.json'
    cases = json.loads(radius_file.read_text())['cases']
    conversion = json.loads(conversion_file.read_text())
    measured = []
    for radius in (1, 2):
        case, = (c for c in cases if c['name'] == f'r{radius}_n30_box28')
        path = Path(case['stdout_path']).parent
        guard, rows = read_run(path, case['stdout_sha256'])
        if any(guard['exit_codes']) or guard['wall_deadline'] is not None:
            raise ValueError(f'Invalid steady intake record at radius {radius}')
        base, = (r for r in rows if r['test'] == 'stationary_action' and r['bias'] == BIAS)
        feed, bulk = base['feed'], base['bulk']  # positive into the device and out, respectively
        rate = (feed + bulk) / 2
        if rate != case['D'] or min(feed, bulk) <= 0 or abs(feed-bulk) > 2e-8*rate:
            raise ArithmeticError('Current mismatch in two-port calibration')
        # For the native common frequency mu=0 every boundary work integrand
        # vanishes. Choosing a common unit energy origin tests conservation;
        # its tiny remainder is the port current residual of the same field.
        measured.append(dict(radius=radius, spatial_phase_bias=BIAS, feed=feed, bulk=bulk,
                             number_current=rate, native_net_port_power=0.*(feed-bulk),
                             common_unit_energy_offset_current_residual=feed-bulk,
                             minimum_power_times_conversion_capacity_at_this_rate=rate**2))

    native = archived_conversion(conversion, 'native_power_48')
    prep, = (r for r in native if r['test'] == 'preparation')
    front, = (r for r in native if r['test'] == 'conversion' and r['q'] == 0)
    lower = prep['feed']**2 / front['marker_norm']
    if not 0 < lower < front['loss'] or abs(lower-front['bound_from_bath_rate']) > 1e-13:
        raise ArithmeticError('Flux lower bound disagrees with the native conversion reading')

    controls = archived_conversion(conversion, 'curved_matched_spatial_phase')
    comparisons = []
    for state in controls:
        if state['test'] != 'free_spherical_front':
            continue
        head = state['head_brane_minus_bulk']
        supplied = state['current_per_inner_area'] * head
        loss = state['conversion_loss']
        if head <= 0 or loss <= 0 or abs(supplied/loss-1) > 1e-7:
            raise ArithmeticError('Powered curved control fails the open energy budget')
        comparisons.append(dict(head=head, loss_per_inner_area=loss, power_per_inner_area=supplied,
                                relative_balance_error=supplied/loss-1))
    if len(comparisons) != 2:
        raise ValueError('Both independently meshed powered controls are required')

    print(json.dumps(dict(
        input_sha256={str(p): digest(p) for p in (radius_file, conversion_file)},
        rigid_reference_ports=measured,
        known_order_field=dict(radius=1, intake=prep['feed'], phase_bias=prep['phase_bias'],
                               measured_loss=front['loss'], capacity=front['marker_norm'],
                               lower_bound=lower, required_energy_per_constituent=front['loss']/prep['feed']),
        powered_annulus_controls=comparisons), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
