#!/usr/bin/env python3
"""Read complete LAB operator inputs and prior returns without assembly/science."""
import argparse
import ast
import json
from pathlib import Path
import resource
import signal
import time
import S11c_d_remaining_case_frequency_remainder_inputs_finish as metadata

saved = metadata.original.saved
io = metadata.original.io
np, sp = saved.np, saved.sp
M, F, require, same = saved.M, saved.F, saved.require, saved.same
NAME = 'S11c_d_remaining_case_frequency_lab_operator_inputs'
PINS = {
    'remainder_rows': 'd6c9ec1e6492a87782caeebfb36ce60134c682dad2d5967662cd773b0b619e3a',
    'row_pilot': '0c2448600fbd33b82b5e9c1a0d6b16ddbf5dbb423decff540ceec7d221098e15',
    'rows_1d': '2907d7da91fb048b345e0c4e4aa238a182b90e021672f81ecb453f7cf1fac2d5',
    'row_2d': 'a4e1a9cc50f126a40730dac2d5f8886ad527f8d5b667be5ee3be83ccd8706229',
    'row_47': '2d2c98959f64f2edf48c10decd7f25051f8ac660bcd96ff4855966987eb7ad18',
    'scalar': 'd1961c25f725959e896e62e49c45d389fdd503575d27803a89055175f6308b50'}
ROW_INVENTORY_SHA = '2d025399874fe772b8e3199916ccd97978db0ca33f8ef1cb0bbc73a2f1340a4f'


def view(value, depth=0):
    """Bounded metadata only; exact physical values stay at their packet keys."""
    if isinstance(value, np.ndarray):
        return {'type': type(value).__name__, 'shape': list(value.shape), 'dtype': value.dtype.str}
    if isinstance(value, dict):
        return {'type': 'dict', 'count': len(value), 'keys': [repr(k) for k in value],
                'first': {repr(k): view(v, depth+1) for k, v in list(value.items())[:2]} if depth < 3 else None}
    if isinstance(value, (tuple, list)):
        return {'type': type(value).__name__, 'count': len(value), 'first': [view(v, depth+1) for v in value[:2]] if depth < 3 else None}
    if isinstance(value, sp.Basic):
        return {'type': type(value).__name__, 'representation': repr(value)[:3000], 'freeSymbols': sorted(map(str, value.free_symbols))}
    return {'type': type(value).__name__, 'representation': repr(value)[:3000]}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--run-directory', type=Path, required=True)
    base = ap.parse_args().run_directory.resolve(); base.relative_to(io.REPO/'_scratch/s11c'); base.mkdir(parents=True, exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3)); signal.alarm(900); started = time.monotonic()
    io.digest = saved.digest
    reader, journal = saved.Reader(), metadata.MetadataJournal(base)
    cache = {}
    def packet(rec):
        path = rec.get('logical', rec.get('path')); r = reader.retain(path, rec['sha256'])
        if r['canonical'] not in cache: cache[r['canonical']] = reader.packet(path, rec['sha256'])
        return cache[r['canonical']]
    def address(route):
        v = packet(route['packet'])
        for key in route['keys']: v = v[tuple(key) if isinstance(key, list) else key]
        return v
    def accepted(cp, name):
        rec = reader.retain(Path(cp['runDirectory'])/name, cp['artifacts'][name]['sha256'])
        return packet(rec), rec
    cps = {}
    for short, h in PINS.items():
        filename = 'S11c_d_remaining_case_frequency_'+('scalar_bindings' if short == 'scalar' else short)+'_checkpoint.json'
        cps[short] = reader.json(M/filename, h)
        reader.retain(Path(cps[short]['runDirectory'])/'checks.json', cps[short]['checksSha256'])
    bcp = reader.json(M/'S11c_d_frequency_matrix_checkpoint.json', '0c350c5f9c554c2e585920745157be826daf51775e05676a0554d8670a750663')
    reader.retain(Path(bcp['runDirectory'])/'checks.json', bcp['checksSha256'])
    caller = {}
    for name, names in (
        ('S11c_d_frequency_matrix.py', ('bind', 'assemble', 'system_and_solve', 'end_maps', 'complex_case')),
        ('S11c_d_continuum_matrices.py', ('values', 'assemble', 'recombine'))):
        h = bcp['sourceFiles']['_measurements/'+name]
        current = reader.retain(M/name, h); frozen = reader.retain(Path(bcp['runDirectory'])/'source/_measurements'/name, h)
        tree = ast.parse(Path(current['canonical']).read_text())
        caller[name] = {'current': current, 'frozen': frozen, 'bodies': {n.name: ast.unparse(n) for n in tree.body if getattr(n, 'name', None) in names}}
    journal.json('native-operator-callers.json', caller)
    baseline, baseline_route = accepted(bcp, 'complex/frequency-binding.pickle')
    baseline_rows, baseline_rows_route = accepted(bcp, 'complex/frequency-rows.pickle')
    old_returns = {}
    for name in ('complex/frequency-interior.pickle', 'complex/frequency-system.pickle',
                 'complex/frequency-solution.pickle', 'complex/frequency-end-maps.pickle', 'seed/frequency-system.pickle'):
        val, rec = accepted(bcp, name); old_returns[name] = {'packet': rec, 'savedSchema': view(val)}
    journal.json('saved-baseline-whole-operator-returns.json', old_returns)
    prior = reader.json(F/'row-inputs/completed-input-artifact-inventory.json', ROW_INVENTORY_SHA)
    # This is an input consumer, not a repeat of the380-row/full-source reader.
    def old_json(name):
        return reader.json(F/'row-inputs/complete'/name, prior[name]['sha256'])
    news = {28: {'checkpoint': 'row_pilot', 'input': cps['row_pilot']['rowInput'], 'value': cps['row_pilot']['rowReturn']}}
    for row in cps['rows_1d']['rows']: news[row['rowIndex']] = {'checkpoint': 'rows_1d', 'input': row['rowInput'], 'value': row['rowReturn']}
    for short, index in (('row_2d', 46), ('row_47', 47)):
        news[index] = {'checkpoint': short, 'input': cps[short]['rowInput'], 'value': {'packet': cps[short]['rowReturn'], 'keys': []}}
    for row in cps['remainder_rows']['rowResults']:
        news[row['rowIndex']] = {'checkpoint': 'remainder_rows', 'input': row['input'], 'value': {'packet': row['fineValue'], 'keys': []}}
    require(set(news) == {28,29,30,46,47,51,52,53,54}, 'actual nine accepted LAB RHOBR row owners')
    for route in news.values():
        if 'packet' not in route['value']: route['value'] = {'packet': route['value'], 'keys': []}
    scalar_cp = cps['scalar']; case_results = {}
    baseline_scalars = {tuple(v['address']): v for v in baseline['sourceJoins']}
    def forbidden(*args, **kwargs): raise RuntimeError('operator input inspection disables new science')
    journal.write = forbidden
    io.native.f.source_jets = io.native.f.polynomial_basis = io.native.f.BasisMomentum.prepare_basis = forbidden
    io.native.Pair.__init__ = io.native.maps = io.native.continue_pair = forbidden
    for name in ('diff','lambdify','cancel','expand','factor','solve','gcd','resultant','integrate'): setattr(sp, name, forbidden)
    for label in ('LAB_HELD__RHO4_CONSTANT', 'LAB_HELD__RHOBR_CONSTANT'):
        sr = Path(scalar_cp['runDirectory']); prefix = 'cases/'+label
        scalar, scalar_route = accepted(scalar_cp, prefix+'/scalar-bindings.pickle')
        routes = reader.json(sr/prefix/'input-routes.json', scalar_cp['artifacts'][prefix+'/input-routes.json']['sha256'])
        raw = packet(routes['originalNativeRowsSourcesProfilesSettings']); physical = packet(routes['physicalRowsTermsCharacters'])
        source = packet(routes['source']); context = packet(routes['context']); analytic = packet(routes['analytic'])
        journal.json(prefix+'/physical-operator-inputs.json', {'scalar': scalar_route, 'routes': routes, 'sourceSchema': view(source),
            'nativeRowsSchema': view(raw), 'physicalTermsSchema': view(physical), 'contextSchema': view(context),
            'analyticSchema': view(analytic), 'frequency': repr(scalar['frequency'])})
        summary = old_json(prefix+'/summary.json'); selected_rows = []; same_baseline_rows = []
        for row in summary['rows']:
            index = row['rowIndex']; old = old_json(prefix+'/row-'+str(index)+'.json')
            if old['status'] == 'SAVED_COMPLETE_NUMERICAL_ROW':
                chosen = old['savedCompleteMatches'][0]; route = chosen['result']; row_input = chosen['binding']; disposition = 'accepted complete row input match'
            else:
                require(label == 'LAB_HELD__RHOBR_CONSTANT' and index in news, 'all actual missing LAB row inputs now have accepted results')
                chosen = news[index]; route = chosen['value']; row_input = chosen['input']; disposition = 'accepted new own numerical row'
                packet(row_input)
            array = address(route)
            require(isinstance(array,np.ndarray) and array.shape == (129,129) and array.dtype == np.dtype(complex) and np.isfinite(array).all(), 'full actual consumer row array')
            exact = index in baseline_rows['rows'] and same(array, baseline_rows['rows'][index])
            same_baseline_rows.append(exact)
            selected_rows.append({'rowIndex': index, 'input': row_input, 'value': route, 'disposition': disposition,
                                  'savedFullRowInputJoin': reader.retain(F/'row-inputs/complete'/prefix/('row-'+str(index)+'.json')),
                                  'equalsSameIndexBaselineArray': exact})
        journal.json(prefix+'/complete-row-routes.json', selected_rows)
        scalar_views = []
        for adr, expression in scalar['actual'].items():
            if adr[0] not in ('local','cell'): continue
            prior_scalar = baseline_scalars.get(adr)
            scalar_views.append({'address': list(adr), 'actualBound': view(expression), 'baselineSameAddress': prior_scalar is not None,
                'sameSavedBoundValue': prior_scalar is not None and same(expression, prior_scalar['bound']),
                'scalarPacket': scalar_route, 'actualKey': ['actual', list(adr)],
                'keyEncoding': 'The address is an actual tuple dictionary key; JSON lists are metadata only.'})
        journal.json(prefix+'/local-and-cell-scalar-inputs.json', scalar_views)
        end_routes = reader.json(sr/prefix/'end-map-routes.json', scalar_cp['artifacts'][prefix+'/end-map-routes.json']['sha256'])
        end_views = {}
        for side, route in end_routes.items():
            actual = packet(route['map']['file'])
            for key in route['map']['keys']: actual = actual[key]
            reader.retain(route['wholeInputJoin']['logical'], route['wholeInputJoin']['sha256'])
            for rec in route['physicalSources']: reader.retain(rec['logical'], rec['sha256'])
            reader.retain(route['selection']['logical'], route['selection']['sha256'])
            end_views[side] = {'route': route, 'actualSavedMapSchema': view(actual)}
        journal.json(prefix+'/saved-forcing-observation-end-inputs.json', end_views)
        case_results[label] = {'rows': len(selected_rows), 'localCellScalars': len(scalar_views),
            'scalarSameAddressBaselineMatches': sum(v['sameSavedBoundValue'] for v in scalar_views),
            'rowSameIndexBaselineMatches': sum(same_baseline_rows), 'sourceSummary': source.get('summary'),
            'sourceRecordCount': len(source['records']), 'physicalTermView': view(physical['terms']),
            'settings': raw['settings'], 'wholeOperatorReuseAccepted': False,
            'scope': 'Saved component comparisons only. Whole actual grade/cell/native caller and both end/forcing/observation inputs must join before any complete operator or response reuse.'}
        journal.json(prefix+'/summary.json', case_results[label])
        # Drop bulky case scientific data; retained input references remain exact.
        del raw, physical, source, context, analytic
    for path in (Path(__file__).resolve(), M/(NAME+'_plan.md'), Path(metadata.__file__),
                 M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):
        rec = reader.retain(path); target = base/'source'/path.name; target.parent.mkdir(exist_ok=True)
        with target.open('xb') as out: out.write(Path(rec['canonical']).read_bytes())
        reader.retain(target,rec['sha256'])
    reader.postcheck(); journal.json('inputs.json', {'consumedRoutes': reader.routes, 'checkpoints': PINS})
    checks = {'status': 'COMPLETED_SAVED_LAB_OPERATOR_INPUT_INSPECTION', 'cases': case_results, 'newScientificCalls': 0,
              'newMatrixMapSolveCalls': 0, 'allConsumedHashesUnchanged': True, 'consumedLogicalPaths': len(reader.routes),
              'artifacts': journal.artifacts, 'wallSeconds': time.monotonic()-started,
              'scope': 'Practical analog toy-model input inspection for actual LAB frequency operators. No assembly/response/reuse certification or numerical operation.'}
    journal.json('checks.json', checks); signal.alarm(0); print(json.dumps(checks,indent=2))


if __name__ == '__main__': main()
