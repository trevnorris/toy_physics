#!/usr/bin/env python3
"""Read saved whole-cluster inputs and completed numerical end results.

This is a bounded preparation for the toy-model continuation pilot. It writes
only metadata addresses into immutable accepted packets; it performs no
scientific construction, numerical continuation, or mathematical validation.
"""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import pickle
import resource
import signal
import time

import numpy as np
import sympy as sp
import S11c_d_frequency_end as native

M = Path(__file__).resolve().parent
LEDGER = M.parent
REPO = LEDGER.parent.parent
CP = M / 'S11c_d_remaining_case_frequency_end_complex_threshold_checkpoint.json'
CP_SHA = '38c47b224f9c192aed3cec1d677bbb2e64b1ce4ac7a9efacb975dbb26a6029e7'
PLAN = M / 'S11c_d_remaining_case_frequency_end_continuation_inputs_plan.md'


def require(value, message):
    if not value:
        raise AssertionError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def same(a, b, memo=None):
    """Exact saved identity, without normalization, coercion, or tolerance."""
    if type(a) is not type(b):
        return False
    if a is b:
        return True
    memo = set() if memo is None else memo
    pair = (id(a), id(b))
    if pair in memo:
        return True
    memo.add(pair)
    if isinstance(a, np.ndarray):
        return a.dtype == b.dtype and a.shape == b.shape and np.array_equal(a, b)
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k], memo) for k in a)
    if isinstance(a, (tuple, list)):
        return len(a) == len(b) and all(same(x, y, memo) for x, y in zip(a, b))
    return bool(a == b)


def save(base, name, value):
    path = base / name
    path.relative_to(base)
    require('..' not in Path(name).parts, 'no output escape')
    parent = path.parent
    while parent != base:
        require(not parent.is_symlink(), 'no output through reference parent')
        parent = parent.parent
    require(not path.exists() and not path.is_symlink(), 'new metadata only')
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as out:
        out.write(json.dumps(value, indent=2) + '\n')


class Reader:
    def __init__(self):
        self.routes = {}
        self.physical = {}

    def retain(self, path, expected=None):
        path = Path(os.path.abspath(path))
        canonical = path.resolve(strict=True)
        route = {'logical': str(path), 'canonical': str(canonical),
                 'directLink': str(path.readlink()) if path.is_symlink() else None,
                 'bytes': canonical.stat().st_size}
        if str(canonical) not in self.physical:
            self.physical[str(canonical)] = digest(canonical)
        route['sha256'] = self.physical[str(canonical)]
        require(expected is None or route['sha256'] == expected, ('accepted input hash', str(path)))
        if str(path) in self.routes:
            require(self.routes[str(path)] == route, 'stable logical/canonical route')
        self.routes[str(path)] = route
        return route

    def json(self, path, expected=None):
        route = self.retain(path, expected)
        return json.loads(Path(route['canonical']).read_text())

    def packet(self, path, expected):
        route = self.retain(path, expected)
        with Path(route['canonical']).open('rb') as stream:
            return pickle.load(stream)

    def postcheck(self):
        for path, expected in self.physical.items():
            require(digest(Path(path)) == expected, ('unchanged consumed bytes', path))
        for path, prior in list(self.routes.items()):
            require(self.retain(path) == prior, ('unchanged consumed route', path))


def source_join(reader, cp, origin):
    # The accepted checkpoint retains all ancestral inventories. Check consumed
    # current/frozen sources directly, without copying or nesting those files.
    for name, expected in cp['sourceFiles'].items():
        reader.retain(LEDGER / name, expected)
        reader.retain(origin / 'source' / name, expected)
    sources = {}
    for path, names in (
        (M / 'S11c_d_frequency_end.py', ('Pair', 'seeds', 'maps', 'continue_pair')),
        (M / 'S11c_d_frequency_matrix.py', ('end_maps',)),
        (M / 'S11c_d_frequency_contour.py', ('adapters', 'close_end_paths')),
        (M / 'S11c_d_uniform_response.py', ('match',)),
    ):
        text = path.read_text()
        module = ast.parse(text)
        bodies = {}
        for name in names:
            node = next(n for n in module.body if getattr(n, 'name', None) == name)
            bodies[name] = {'source': ast.get_source_segment(text, node),
                            'wholeBodyAST': hashlib.sha256(ast.dump(node).encode()).hexdigest()}
        sources[str(path)] = {'file': reader.retain(path), 'bodies': bodies}
    for path in (Path(__file__).resolve(), PLAN, LEDGER / 'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):
        reader.retain(path)
    return sources


def prohibit():
    def forbidden(*args, **kwargs):
        raise RuntimeError('saved continuation input reader cannot execute science')
    for name in ('load', 'seeds', 'maps', 'continue_pair', 'focused', 'polynomial', 'main'):
        setattr(native, name, forbidden)
    for name in ('__init__', 'initial', 'coefficients', 'equation', 'jacobian', 'solve'):
        setattr(native.Pair, name, forbidden)
    for name in ('diff', 'lambdify', 'cancel', 'factor', 'expand', 'simplify', 'solve', 'resultant', 'gcd', 'gcdex', 'integrate'):
        setattr(sp, name, forbidden)
    for name in ('diff', 'subs', 'xreplace'):
        setattr(sp.Basic, name, forbidden)
    for lib in (np.linalg, native.la):
        for name in ('solve', 'inv', 'svd', 'eig', 'eigh', 'lstsq', 'pinv', 'expm'):
            if hasattr(lib, name):
                setattr(lib, name, forbidden)
    native.f.atomic_pickle = forbidden
    pickle.dump = forbidden
    pickle.dumps = forbidden


def address(route, *keys):
    return {'file': route, 'keys': list(keys)}


def array_view(value):
    return {'type': type(value).__module__ + '.' + type(value).__name__,
            'shape': list(value.shape), 'dtype': str(value.dtype)}


def numerical_atlas(reader, base):
    result = {}
    for tag, filename in (
        ('focused', 'S11c_d_frequency_end_focused.json'),
        ('matrix', 'S11c_d_frequency_matrix_checkpoint.json'),
        ('contour16', 'S11c_d_frequency_contour_checkpoint.json'),
        ('contour32', 'S11c_d_frequency_contour_refine_checkpoint.json'),
    ):
        cp = reader.json(M / filename)
        root = Path(cp['runDirectory'])
        if 'artifactInventory' in cp:
            reference = cp['artifactInventory']
            checks = reader.json(reference['path'], reference['sha256'])
            inventory = checks[reference['key']]
        else:
            reader.retain(root / 'checks.json', cp['checksSha256'])
            inventory = cp['artifacts']
        for name, expected in cp['sourceFiles'].items():
            reader.retain(LEDGER / name, expected)
        selected = {}
        for name, item in inventory.items():
            if tag == 'focused' or name.endswith('/frequency-end-maps.pickle') or name == 'accepted-end-seeds.pickle':
                selected[name] = reader.retain(root / name, item['sha256'])
        maps = []
        for name, route in selected.items():
            if not name.endswith('/frequency-end-maps.pickle'):
                continue
            value = reader.packet(root / name, route['sha256'])
            ends = {}
            for end, packet in value.items():
                clusters = packet['clusters']
                ends[end] = {'frequencyRepresentations': [repr(v['state']['frequency']) for v in clusters],
                             'clusterIndices': [int(v['seed']['index']) for v in clusters],
                             'right': array_view(packet['right']),
                             'incomingRight': array_view(packet['incomingRight']),
                             'completeSavedMap': address(route, end)}
            maps.append({'name': name, 'ends': ends})
            del value
        result[tag] = {'checkpoint': reader.retain(M / filename), 'status': cp['status'],
                       'allAcceptedArtifactNames': list(inventory), 'selectedSavedResults': selected,
                       'savedMaps': maps, 'newNumericalCalls': 0}
        save(base, 'baseline/' + tag + '.json', result[tag])
    return result


def table_values(table):
    # Literal fields read by the whole native Pair constructor. Keep the full
    # entry records as well as all source fields actually used there.
    return {'entries': table['entries'], 'source': {k: table['source'][k] for k in
            ('frequency', 'momentum', 'radical', 'wave', 'livePencil')}}


def inspect_inputs(reader, base, cp, origin):
    def route(name):
        return reader.retain(origin / name, cp['artifacts'][name]['sha256'])

    def packet(name):
        return reader.packet(origin / name, cp['artifacts'][name]['sha256'])

    saved = packet('end-accepted/continuation/end-seeds.pickle')
    saved_route = route('end-accepted/continuation/end-seeds.pickle')
    baseline_tables = {end: packet('end-accepted/chart/' + end.lower() + '-rational-end.pickle')
                       for end in ('LEFT', 'RIGHT')}
    binding_routes = {}
    for name in ('PENCIL_PLUS', 'FREQUENCY_PENCIL_PLUS', 'NORMAL_PENCIL_PLUS',
                 'FREQUENCY_TRANSPORT', 'NORMAL_TRANSPORT', 'wave'):
        stem = 'binding-operations/' + name + '/'
        receipt = reader.json(origin / (stem + 'completed.json'), cp['artifacts'][stem + 'completed.json']['sha256'])
        ip, vp = route(stem + 'input.pickle'), route(stem + 'value.pickle')
        require(ip['sha256'] == receipt['inputSha256'] and vp['sha256'] == receipt['valueSha256'], 'actual saved binding receipt')
        value = packet(stem + 'value.pickle')
        binding_routes[name] = {'input': ip, 'value': vp, 'completed': route(stem + 'completed.json'),
                                'savedValueType': str(type(value)),
                                'savedValueRepresentation': repr(value)}
    save(base, 'saved-binding-and-tangent-routes.json', binding_routes)
    families = []
    results = []
    for case in cp['cases']:
        context_name = 'end-input-cases/' + case + '/context-pairs.pickle'
        context_route = route(context_name)
        context = packet(context_name)
        for end in ('LEFT', 'RIGHT'):
            directory = 'end-input-cases/' + case + '/' + end.lower() + '/'
            full_name = directory + 'full-inputs.pickle'
            selection_name = directory + 'whole-cluster-inputs.pickle'
            full = packet(full_name)
            selection = packet(selection_name)
            rational_name = 'rational-cases/' + case + '/' + end.lower() + '/route.json'
            rational = reader.json(origin / rational_name, cp['artifacts'][rational_name]['sha256'])
            table_name = ('new-rational-end/right-rational-end.pickle' if rational['mode'] == 'shared-new-rational-table'
                          else 'end-accepted/chart/' + end.lower() + '-rational-end.pickle')
            table = packet(table_name)
            require(Path(rational['tablePath']).resolve() == (origin / table_name).resolve(), 'actual accepted rational source route')
            require(same(full['address'], (case, end)), 'own physical source address')
            modes = {v['info']['INDEX']: v for v in selection['allCandidateModes']}
            incoming = {v['RECORD_INDEX'] for v in selection['channel']['incoming']}
            rows = []
            for index, items in selection['groups'].items():
                mode = modes[index]
                info = mode['info']
                # Locate an actual saved complete matrix whose columns are the
                # native selected vectors in order. Never assemble a new R.
                bases = []
                for key in ('right', 'rawRight', 'fluxRight'):
                    value = mode.get(key)
                    if isinstance(value, np.ndarray) and value.shape == (5, len(items)):
                        if all(same(value[:, j], item['vector']) for j, item in enumerate(items)):
                            bases.append(key)
                raw_basis = None if not bases else mode[bases[0]]
                candidates = []
                direction = 'incoming' if index in incoming else 'outgoing'
                for position, old in enumerate(saved[end]['clusters']):
                    fields = {'fullTable': same(table, baseline_tables[end]),
                              'consumedTableFields': same(table_values(table), table_values(baseline_tables[end])),
                              'fullItems': same(items, old['items']), 'fullOriginalMode': same(mode, old['originalMode']),
                              'wholeBasis': raw_basis is not None and same(raw_basis, old['R']),
                              'rawK': same(info['K'], old['originalMode']['info']['K']),
                              'rawQ': same(info['Q'], old['originalMode']['info']['Q']),
                              'kind': same(items[0]['kind'], old['kind']), 'direction': direction == old['direction']}
                    candidates.append({'index': int(old['index']), 'fields': fields,
                                       'completeSavedCallerMatch': all(fields.values()),
                                       'savedInputAndSeedState': address(saved_route, end, 'clusters', position)})
                request = {'table': table_values(table), 'items': items, 'mode': mode,
                           'basis': raw_basis, 'kind': items[0]['kind'], 'direction': direction,
                           'context': context['context'], 'seedInput': context['seedInput'],
                           'fieldUnits': full['signature']['fieldUnits'], 'equationUnits': full['signature']['equationUnits']}
                matches = [i for i, (_, prior) in enumerate(families) if same(request, prior)]
                family = matches[0] if matches else len(families)
                if not matches:
                    families.append(((case, end, int(index)), request))
                record = {'case': case, 'end': end, 'index': int(index), 'family': family,
                          'firstOwner': families[family][0], 'nullity': int(info['NULLITY']),
                          'basisColumns': [int(v['BASIS_COLUMN']) for v in items], 'kind': items[0]['kind'],
                          'direction': direction, 'savedWholeBasisKeys': bases,
                          'rawK': {'representation': repr(info['K']), 'type': str(type(info['K']))},
                          'rawQ': {'representation': repr(info['Q']), 'type': str(type(info['Q']))},
                          'fullNativeSelection': address(route(selection_name), 'groups', int(index)),
                          'fullOriginalMode': address(route(selection_name), 'allCandidateModes',
                                                     next(i for i, v in enumerate(selection['allCandidateModes']) if v is mode)),
                          'context': context_route, 'fullOwnInput': route(full_name),
                          'rationalTable': route(table_name), 'acceptedSourceRoute': route(rational_name),
                          'baselineCandidates': candidates,
                          'numericalReuseAccepted': False, 'newScientificCalls': 0}
                save(base, 'cases/' + case + '/' + end.lower() + '/cluster-' + str(index) + '.json', record)
                rows.append(record)
            result = {'case': case, 'end': end, 'candidateCount': len(modes),
                      'selectedClusters': len(rows), 'selectedDirections': len(selection['selectedItems']),
                      'clusters': rows, 'physicalSourceKeptSeparate': True}
            save(base, 'cases/' + case + '/' + end.lower() + '/summary.json', result)
            results.append(result)
    return results, [owner for owner, _ in families]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-directory', type=Path, required=True)
    args = parser.parse_args()
    base = args.run_directory.resolve()
    base.relative_to(REPO / '_scratch/s11c')
    base.mkdir(parents=True, exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS, (2 * 1024**3, 2 * 1024**3))
    signal.alarm(900)
    started = time.monotonic()
    reader = Reader()
    cp = reader.json(CP, CP_SHA)
    require(cp['status'] == 'ACCEPTED_CASE_FREQUENCY_END_COMPLEX_THRESHOLD_CANDIDATES', 'accepted upstream scope')
    origin = Path(cp['runDirectory'])
    reader.retain(origin / 'checks.json', cp['checksSha256'])
    validation = cp['validation']
    reader.retain(Path(validation['runDirectory']) / 'checks.json', validation['checksSha256'])
    joins = source_join(reader, cp, origin)
    save(base, 'native-callers.json', joins)
    prohibit()
    atlas = numerical_atlas(reader, base)
    ends, families = inspect_inputs(reader, base, cp, origin)
    clusters = [v for end in ends for v in end['clusters']]
    reader.postcheck()
    save(base, 'inputs.json', {'acceptedCheckpoint': reader.routes[str(CP)],
         'ancestralInventoriesRemainInCheckpoint': True, 'consumedRoutes': reader.routes,
         'consumedPhysicalHashes': reader.physical, 'noPhysicalFilesCopiedOrWritten': True})
    checks = {'status': 'COMPLETED_SAVED_END_CONTINUATION_INPUT_INSPECTION',
              'physicalEnds': len(ends), 'clusterUses': len(clusters), 'literalFullInputFamilies': len(families),
              'familyOwners': families,
              'wholeBasisRoutesFound': sum(bool(v['savedWholeBasisKeys']) for v in clusters),
              'completeBaselineCallerCandidates': sum(any(x['completeSavedCallerMatch'] for x in v['baselineCandidates']) for v in clusters),
              'savedMapFiles': {k: len(v['savedMaps']) for k, v in atlas.items()},
              'consumedLogicalPaths': len(reader.routes), 'uniquePhysicalFiles': len(reader.physical),
              'allConsumedSourcesAndInputsUnchanged': True, 'newScientificCalls': 0,
              'newPhysicalPickleWrites': 0, 'numericalReuseAccepted': False,
              'next': 'Use these actual saved full input routes for one bounded numerical end-continuation pilot near 1-0.01i. No global-domain survey or repeated baseline continuation.',
              'scope': 'Analog toy-model preparation only. Algebraic candidates are not physical thresholds or poles. No new continuation result is claimed.',
              'wallSeconds': time.monotonic() - started,
              'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    save(base, 'checks.json', checks)
    signal.alarm(0)
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
