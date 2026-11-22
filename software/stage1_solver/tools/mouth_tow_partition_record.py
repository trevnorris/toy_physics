"""Collect the guarded steady stream, frame and tow readouts without replay."""
import hashlib
import json
from pathlib import Path

ROOT = Path('_scratch/s11c/tow_pivot01')
TOOLS = Path('software/stage1_solver/tools')
MEASURE = Path('research/pde_ledger_v3/_measurements')


def load(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def guarded(name, expected):
    folder = ROOT/name
    call = load(folder/'invocation.json')
    limits = load(folder/'effective-limits.json')
    result = load(folder/'outcome.json')
    if (call['command'] != expected or call['wallDeadlineSeconds'] is not None
            or int(limits['memory.swap.max']) or int(limits['memory.max']) != call['memoryMax']
            or not result['limitsVerified'] or result['exitCode'] != 0
            or result['childOutcome']['exitCode'] != 0 or result['childOutcome']['guardReason']
            or (folder/'stderr').stat().st_size):
        raise ValueError(f'Guard or command mismatch for {name}')
    log = folder/'stdout'
    return dict(command=expected, cap_bytes=call['memoryMax'], swap_bytes=0, limits_verified=True,
                exit_code=0, stdout_sha256=sha(log), stdout_file=str(log),
                rows=[json.loads(line) for line in log.read_text().splitlines()])


def main():
    partition = TOOLS/'mouth_steady_tow_partition.py'
    frame = TOOLS/'mouth_tow_frame_check.py'
    runs = {'pair18':guarded('pair18', ['python3', str(partition)])}
    for far in (9,18):
        runs[f'frame{far}'] = guarded(f'frame{far}', ['python3',str(frame),str(far)])
    small = runs['frame9']['rows'][-1]
    large = runs['pair18']['rows'][0]['results']
    small_rows = {r['drive']:r for r in small['old_stream_rows']}
    small_base = load(MEASURE/'mouth_order_and_flare.json')['output_field_sha256']
    if small['field_sha256'] != small_base['_scratch/s11c/mouth_flared_cubic_n34/field.npz']:
        raise ValueError('Small-box base changed')
    if runs['pair18']['rows'][0]['base_sha256'] != runs['frame18']['rows'][0]['field_sha256']:
        raise ValueError('Large-box tests used different backgrounds')
    table = []
    for far in (9,18):
        if far == 9:
            components = [small_rows[k]['outer_cut'] for k in ('feed','bulk')]
            windows = [small_rows[k]['interior'] for k in ('feed','bulk')]
            direct = -small['lab_chart_tow']['medium_on_rigid']
            direct_windows = [-v for v in small['lab_chart_tow']['interior']]
            boost = small['uniform_nodal_boost']
            common = small_rows['tow']['outer_cut']
        else:
            components = [large[k]['medium_on_rigid'] for k in ('brane','bulk')]
            windows = [large[k]['interior'] for k in ('brane','bulk')]
            direct = large['comparison']['positive_lab_tow_external']
            direct_windows = large['comparison']['archived_tow_windows']
            boost = large['nodal_lab_boost']
            common = large['common_minus']['medium_on_rigid']
        summed = sum(components)
        native_mass = -runs[f'frame{far}']['rows'][-1]['analytic_product_boost_medium']
        pred = summed-boost['medium_on_rigid']
        pred_windows = [sum(pair)-term for pair,term in zip(zip(*windows),boost['interior'])]
        if (abs(summed+common) > 1e-10 or abs(pred-direct)>1e-9
                or abs(summed+native_mass-direct)>1e-9
                or max(abs(a-b) for a,b in zip(pred_windows,direct_windows)) > 1e-9):
            raise ValueError(f'Frame comparison failed in box {far}')
        table.append(dict(far=far, stream_medium_coefficients=components,
                          raw_predicted_tow_external=summed, stationary_common_linear_residual=summed+common,
                          base_trace_outward_mass=native_mass, prediction_analytic_product=summed+native_mass,
                          prediction_nodal_boost=pred, measured_tow_external=direct,
                          nodal_prediction_minus_tow=pred-direct, window_predictions=pred_windows,
                          measured_tow_windows=direct_windows))
    sources = [partition, frame, Path(__file__), TOOLS/'mouth_carrier_lateral_response.py',
               TOOLS/'mouth_carrier_lateral_inertia.py', TOOLS/'mouth_carrier_free_center.py',
               TOOLS/'mouth_occupied_push_baths.py', TOOLS/'mouth_flared_receiver.py',
               MEASURE/'mouth_carrier_0pn.json', MEASURE/'mouth_rigid_drag_checks.json']
    record = dict(scope='Steady linear coherent carrier on rigid prescribed cubic reference, sealed faces, open exit',
                  source_sha256={str(p.relative_to(Path.cwd()) if p.is_absolute() else p):sha(p) for p in sources},
                  runs=runs, frame_comparison=table)
    destination = MEASURE/'mouth_steady_tow_partition.json'
    destination.write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
    print(json.dumps(table, allow_nan=False))


if __name__ == '__main__':
    main()
