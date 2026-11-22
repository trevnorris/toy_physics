"""Required residence response for the earlier, local, stationary order-force calibration."""
import argparse
import hashlib
import json
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(lines):
    return [json.loads(line) for line in lines.splitlines() if line]


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('order_record', type=Path)
parser.add_argument('capture_run', type=Path)
args = parser.parse_args()
record = json.loads(args.order_record.read_text())
calibration, = [r for r in rows(record['guarded_runs']['mesh48']['literal_stdout'])
                if r.get('test') == 'binding' and r['location'] == 'half'
                and r['law'] == 'retain_mean_attraction' and r['width'] == 0.01]
outcome = json.loads((args.capture_run / 'outcome.json').read_text())
assert outcome['limitsVerified'] and outcome['exitCode'] == outcome['childOutcome']['exitCode'] == 0
capture_path = args.capture_run / 'stdout'
results = [r for r in rows(capture_path.read_text()) if r.get('test') == 'capture_response'
           and r['scope'] == 'whole' and r['cells'] == 32]
assert sorted(r['signed_magnitude'] for r in results) == [0.001, 0.002]
fixed_force = calibration['fixed_count_medium_force_per_mass_U']
count_force = calibration['count_dipole_force_coefficient']
assert count_force < 0
critical_count_slope = -fixed_force / count_force
for r in results:
    capture_slope = r['fractional_odd_supply']
    print(json.dumps(dict(
        test='local_stationary_residence_condition', drive=r['signed_magnitude'],
        order_record_sha256=digest(args.order_record), capture_stdout_sha256=digest(capture_path),
        captured_supply_slope=capture_slope, critical_resident_count_slope=critical_count_slope,
        lifetime_slope_must_be_below=critical_count_slope - capture_slope,
        assumptions='dilute stationary renewal, capture and residence at the same depth and angle, linear pole response, relaxed-host order-force calibration')))
