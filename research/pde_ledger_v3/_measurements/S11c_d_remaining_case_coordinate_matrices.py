#!/usr/bin/env python3
"""Material row reuse and common-Eulerian interior assembly from saved bindings."""
import argparse
import ast
import copy
import gc
import hashlib
import inspect
import json
from pathlib import Path
import resource
import shutil
import signal
import textwrap
import time
from types import SimpleNamespace

import numpy as np
import S11c_d_remaining_case_coordinate_bindings as h
import S11c_d_coordinate_response_production as material_native

f, engine, sp, native, modes = h.f, h.engine, h.sp, h.native, h.modes
matrices, interior, c = h.matrices, material_native.interior, h.c
BASELINE = h.BASELINE
CP = f.M/'S11c_d_remaining_case_coordinate_bindings_checkpoint.json'
ECP = f.M/'S11c_d_remaining_case_matrices_checkpoint.json'
FCP = f.M/'S11c_d_remaining_case_coordinate_matrices_focused.json'
PLAN = f.M/'S11c_d_remaining_case_coordinate_matrices_plan.md'
SCOPE = ('Material interior operators in the transported common Eulerian trial basis. '
         'Only new numerical rows; no boundary replacement, mode/current or response solve.')


def body(function):
    return hashlib.sha256(ast.dump(ast.parse(textwrap.dedent(inspect.getsource(function)))).encode()).hexdigest()


def native_joins():
    function, join = material_native.native_accumulator()
    f.require(function.__code__.co_code == material_native.MaterialProduction.matrix_group.__code__.co_code,
              'actual native material accumulator bytecode')
    return {'accumulator': join, 'bodies': {name: body(fun) for name, fun in (
        ('MaterialProduction', material_native.MaterialProduction), ('MaterialMomentum', c.MaterialMomentum),
        ('Chart', c.Chart), ('BasisMomentum', f.BasisMomentum), ('assemble', interior.assemble),
        ('direct_cells', matrices.direct_cells), ('values', interior.values), ('recombine', interior.recombine),
        ('mass_controls', material_native.mass_controls), ('restore', h.restore), ('restore_chart', h.restore_chart))}}


def accepted(path, status):
    checkpoint = json.loads(path.read_text()); root = Path(checkpoint['runDirectory'])
    f.require(checkpoint['status'] == status and f.digest(root/'checks.json') == checkpoint['checksSha256'],
              ('accepted complete checkpoint', str(path)))
    checks = json.loads((root/'checks.json').read_text())
    f.require(all(checkpoint[key] == checks[key] for key in ('sourceFiles', 'inputPackets', 'artifacts')),
              'checkpoint and actual producer identities')
    return root, checks, checkpoint


def hash_check(base, manifest):
    for n, sha in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/n) == f.digest(base/'source'/n) == sha, ('current/frozen source', n))
    for n, sha in manifest['inputPackets'].items(): f.require(f.digest(Path(n)) == sha, ('original input', n))
    for n, sha in manifest['copiedInputs'].items(): f.require(f.digest(base/n) == sha, ('unchanged copy', n))


def load(base, resume):
    br, bc, bcp = accepted(CP, 'ACCEPTED_CASE_MATERIAL_BINDINGS')
    er, ec, ecp = accepted(ECP, 'ACCEPTED_FOUR_CASE_INTERIOR_MATRICES')
    manifest = {'runDirectory': str(base), 'sourceFiles': {}, 'inputPackets': {}, 'copiedInputs': {},
                'input': bc['input'], 'settings': bc['settings'], 'scope': SCOPE,
                'bindingChecksSha256': bcp['checksSha256'], 'eulerianChecksSha256': ecp['checksSha256']}
    for root, checks in ((br, bc), (er, ec)):
        for n, sha in checks['sourceFiles'].items():
            f.require(n not in manifest['sourceFiles'] or manifest['sourceFiles'][n] == sha, 'common native source')
            f.require(f.digest(f.ROOT/n) == f.digest(root/'source'/n) == sha, 'accepted source/current join')
            manifest['sourceFiles'][n] = sha
        for n, sha in checks['inputPackets'].items():
            f.require(n not in manifest['inputPackets'] or manifest['inputPackets'][n] == sha, 'common original input')
            manifest['inputPackets'][n] = sha
        manifest['inputPackets'][str(root/'checks.json')] = f.digest(root/'checks.json')
    for path in (Path(__file__).resolve(), PLAN, CP, ECP, Path(h.__file__).resolve(),
                 Path(material_native.__file__).resolve(), Path(matrices.__file__).resolve()):
        manifest['sourceFiles'][str(path.relative_to(f.ROOT))] = f.digest(path)
    if resume:
        old = json.loads((resume/'checks.json').read_text()); focus = json.loads(FCP.read_text())
        f.require(focus['status'] == 'ACCEPTED_CASE_MATERIAL_MATRIX_INPUTS' and
                  focus['runDirectory'] == str(resume) and focus['checksSha256'] == f.digest(resume/'checks.json'),
                  'accepted material matrix input focus')
        f.require(old['sourceFiles'] == manifest['sourceFiles'], 'unchanged complete focus implementation')
        for n, item in old['artifacts'].items(): modes.retain(resume/n, base/n, manifest, item['sha256'])
        modes.retain(resume/'checks.json', base/'accepted-focus-checks.json', manifest)
        modes.retain(resume/'inputs.json', base/'accepted-focus-inputs.json', manifest)
        manifest['completedFocusReuse'] = {'directory': str(resume), 'artifacts': len(old['artifacts']),
                                         'checksSha256': f.digest(resume/'checks.json')}
    else:
        for n, item in bc['artifacts'].items(): modes.retain(br/n, base/'bindings'/n, manifest, item['sha256'])
        modes.retain(br/'checks.json', base/'accepted-binding-checks.json', manifest)
        modes.retain(br/'inputs.json', base/'accepted-binding-inputs.json', manifest)
        vr = Path(bcp['validator']['runDirectory'])
        f.require(f.digest(vr/'checks.json') == bcp['validator']['checksSha256'], 'accepted independent binding validator')
        modes.retain(vr/'checks.json', base/'accepted-binding-validation.json', manifest)
        names = ['accepted-finite-system.pickle']
        for label in ec['cases']:
            names += ['cases/'+label+'/'+kind+'.pickle' for kind in
                      ('row-matrices', 'interior-matrices', 'direct-independent', 'direct-unsplit')]
        for n in names: modes.retain(er/n, base/'eulerian'/n, manifest, ec['artifacts'][n]['sha256'])
        modes.retain(er/'checks.json', base/'accepted-eulerian-checks.json', manifest)
    for n, sha in manifest['sourceFiles'].items():
        path = base/'source'/n; path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(f.ROOT/n, path); f.require(f.digest(path) == sha, 'frozen matrix helper source')
    f.save(base/'inputs.json', manifest)
    hash_check(base, manifest)
    return manifest, tuple(bc['cases'])


class SavedCoefficients:
    """Address/unit-checked saved images behind the unchanged assembler API."""
    def __init__(self, case, material, origin):
        self.input = SimpleNamespace(origin=origin)
        self.entries = []; self.buckets = {}; self.pending = None; self.position = 0
        for item in case['grades']['records'].values():
            address = tuple(item['address']); record = item['record']
            if address[0] not in ('local', 'cell'): continue
            f.require(tuple(material['bindings'][address]['unit']) == tuple(record['UNIT']), 'full coefficient unit')
            self.add(record['ORIGINAL'], material['bindings'][address]['material'], record['UNIT'], address, None)
            f.require(material['coefficients'][address].keys() == record['COEFFICIENTS'].keys(), 'all saved grades')
            gu = [tuple(engine.PHYSICAL_METADATA.dimensions.measure(v)) for v in case['grades']['generators']]
            for grade, expression in record['COEFFICIENTS'].items():
                unit = tuple(v-sum(n*u[j] for n, u in zip(grade, gu)) for j, v in enumerate(record['UNIT']))
                self.add(expression, material['coefficients'][address][grade], unit, address, grade)

    def add(self, expression, value, unit, address, grade):
        entry = {'expression': expression, 'value': value, 'unit': tuple(unit), 'address': address, 'grade': grade}
        bucket = self.buckets.setdefault(hash(expression), [])
        bucket.append(entry); self.entries.append(entry)

    def bind(self, expression):
        f.require(self.pending is not None and self.position < len(self.pending), 'declared native physical coefficient call sequence')
        entry = self.pending[self.position]
        f.require(native.same(entry['expression'], expression), 'exact next native source expression at declared physical address')
        self.position += 1
        return entry['value']

    def begin(self, kind, direct_addresses=None):
        f.require(self.pending is None, 'previous native coefficient sequence consumed')
        if kind in ('graded', 'unsplit'):
            self.pending = [v for v in self.entries if (v['grade'] is not None) == (kind == 'graded')]
        else:
            f.require(kind == 'direct' and direct_addresses is not None, 'actual independent native cell sequence')
            by_address = {v['address']: v for v in self.entries if v['grade'] is None}
            self.pending = [by_address[address] for address in direct_addresses]
        self.position = 0

    def finish(self):
        f.require(self.pending is not None and self.position == len(self.pending), 'every native coefficient call consumed exactly once')
        self.pending = None

    def addressed(self, expression, unit, address, grade):
        matches = [v for v in self.buckets.get(hash(expression), ()) if native.same(v['expression'], expression)
                   and v['unit'] == tuple(unit) and v['address'] == address and v['grade'] == grade]
        f.require(len(matches) == 1, 'exact coefficient source/address/unit/grade route')
        return matches[0]['value']


def case_inputs(base, label, manifest):
    data, case, packets = h.restore(base/'bindings', label, manifest, saved=True)
    if label == BASELINE:
        material = f.unpickle(base/'bindings/baseline-material/material-binding.pickle')
        fresh, reused = [], [{'row': i, 'fromCase': label, 'fromRow': i} for i in range(len(material['rows']))]
    else:
        saved = f.unpickle(base/'bindings/material-cases'/label/'case-material-binding.pickle')
        f.require(modes.same(case, saved['originalCase']) and native.same(data['coordinate'], saved['coordinate']),
                  'actual full own-case material source and original operator')
        material, fresh, reused = saved['material'], saved['newRows'], saved['reusedRows']
    state = f.unpickle(base/'bindings/chart-state.pickle'); chart = h.restore_chart(data, state)
    system = f.unpickle(base/'eulerian/accepted-finite-system.pickle')
    f.require(modes.same(system, f.unpickle(base/'bindings/accepted-finite-system.pickle')), 'exact original trial system')
    data['system'] = system
    return data, case, packets, material, chart, fresh, reused


def saved_layouts(base, system):
    rows = f.unpickle(base/'bindings/baseline-material-row-matrices.pickle'); covered = {}; evidence = []
    for count in (1, 2, 3):
        path = base/'bindings/baseline-material'/f'layout-{count}.pickle'; group = f.unpickle(path)
        f.require(native.same(group['setting'], system['settings']), 'complete saved material quadrature settings')
        f.require(np.isfinite(group['matrices']).all() and group['matrices'].shape[1:] == (129, 129), 'full material row arrays')
        f.require(abs(group['massResidual']) < 1e-9*(1+abs(group['mass'])), 'saved actual material volume')
        f.require(np.max(abs(group['actionResidual']))/(1+np.max(abs(group['direct']))) < 1e-10, 'saved independent direct action')
        for index, array in zip(group['rowIndices'], group['matrices']):
            f.require(index not in covered and np.array_equal(array, rows[index]), 'unique exact baseline material array')
            covered[index] = True
        changed = group['direct']*1.001
        f.require(not np.array_equal(changed, group['direct']), 'actual saved direct-action weight control responds')
        evidence.append({'path': str(path), 'sha256': f.digest(path), 'rowIndices': group['rowIndices'],
                         'massResidual': group['massResidual'], 'actionMaximum': float(np.max(abs(group['actionResidual']))),
                         'changedWeightMaximum': float(np.max(abs(changed-group['direct'])))})
    f.require(set(covered) == set(rows), 'all accepted baseline rows addressed')
    return evidence


def focus(base, manifest, labels, joins):
    state = f.unpickle(base/'bindings/chart-state.pickle'); common = f.unpickle(base/'bindings/contexts'/BASELINE/'frame.pickle')
    system = f.unpickle(base/'eulerian/accepted-finite-system.pickle')
    f.require(len(system['nodes']) == 129 and json.loads(json.dumps(system['settings'])) == manifest['settings'], 'approved complete trial basis and exact JSON settings view')
    layout_evidence = saved_layouts(base, system); f.atomic_pickle(base/'accepted-layout-evidence.pickle', layout_evidence)
    local_witness = f.unpickle(base/'bindings/baseline-material/material-local.pickle')
    # The common chart and exact saved field-jet compositions certify the reused
    # Eulerian derivative basis; no baseline local assembly is repeated.
    f.require(max(local_witness['derivativeResiduals'].values()) < 1e-9, 'accepted actual transported derivative witness')
    for (order, q), value in state['values']['composition'].items():
        f.require(value == (sp.eye(5) if order == q else sp.zeros(5)), 'saved complete field-jet composition')
    trial = {'nodes': system['nodes'], 'derivativeMatrices': system['derivativeMatrices'],
             'settings': system['settings'], 'chartState': state, 'acceptedLocalWitness': local_witness,
             'scope': 'Exact common field-jet composition plus accepted native transported derivative residuals; no new basis or local matrix.'}
    f.atomic_pickle(base/'transported-trial-reuse.pickle', trial)
    atlas = {}; summaries = {}; total_controls = 0
    for label in labels:
        data, case, packets, material, chart, fresh, reused = case_inputs(base, label, manifest)
        folder = base/'preparation'/label; folder.mkdir(parents=True)
        memo = {}
        same = lambda a, b: native.same(a, b, memo)
        f.require(same(h.frame(data, case), common) and same(material['geometry'], state['values']['g'])
                  and modes.same(material['chartComposition'], state['values']['composition']), 'complete common transported frame')
        f.require(same(material['settings'], system['settings']), 'full case material settings')
        old = f.unpickle(base/'eulerian/cases'/label/'interior-matrices.pickle')
        f.require(modes.same(old['gradeOrigin'], data['adapter'].input.origin) and
                  same(old['fieldUnits'], case['binding']['fieldUnits']) and same(old['equationUnits'], case['binding']['equationUnits'])
                  and same(old['generators'], case['grades']['generators']) and same(old['settings'], material['settings']),
                  'actual own-case matrix units and independent grade origin')
        eulerian_rows_path = base/'eulerian/cases'/label/'row-matrices.pickle'
        eulerian_rows = f.unpickle(eulerian_rows_path)
        f.require(set(eulerian_rows['rows']) == set(range(len(material['rows']))) and
                  same(eulerian_rows['settings'], system['settings']), 'complete own-case saved Eulerian row census/settings')
        f.require(all(a.shape == (129, 129) and np.isfinite(a).all() for a in eulerian_rows['rows'].values()),
                  'actual full own-case Eulerian row arrays')
        array_units = {row['index']: {'rowUnit': row['unit'],
                      'columns': tuple(sorted({material['jets'][v['sourceIndex']]['column'] for v in row['factors']}))}
                       for row in material['rows']}
        f.require(all(len(v['columns']) == 1 for v in array_units.values()), 'one actual source column per saved integral row')
        f.atomic_pickle(folder/'eulerian-row-evidence.pickle', {'path': str(eulerian_rows_path),
                        'sha256': f.digest(eulerian_rows_path), 'arrayUnits': array_units,
                        'fieldUnits': old['fieldUnits'], 'equationUnits': old['equationUnits'],
                        'nodes': system['nodes'], 'derivativeMatrices': system['derivativeMatrices']})
        del eulerian_rows
        adapter = SavedCoefficients(case, material, old['gradeOrigin'])
        for entry in adapter.entries:
            f.require(same(adapter.addressed(entry['expression'], entry['unit'], entry['address'], entry['grade']), entry['value']),
                      'saved typed native assembly argument route')
        for kind in ('graded', 'unsplit'):
            adapter.begin(kind)
            for entry in adapter.pending:
                f.require(same(adapter.bind(entry['expression']), entry['value']), 'saved exact native call sequence')
            adapter.finish()
        # Preserve native full-cell expressions and addresses independently of
        # the coefficient-record traversal used by assemble().
        cell_routes = []; lookup = matrices.binding.integral_addresses(packets['assembly']['result'], packets['factorization']['result'])
        records = {tuple(v['address']): v['record'] for v in case['grades']['records'].values()}
        for cell in packets['assembly']['result']['ROWS']:
            i, j = cell['ROW'], cell['COLUMN']
            for ti, (integral, coefficient) in enumerate(cell['NONLOCAL']):
                index = lookup[id(integral)]; address = ('cell', i, j, ti)
                f.require(same(records[address]['ORIGINAL'], coefficient) and same(integral, material['rows'][index]['original']),
                          'whole original native cell/row identity')
                f.require(all(case['binding']['jets'][v['sourceIndex']]['column'] == j for v in material['rows'][index]['factors']),
                          'all actual source fields in native cell')
                value = adapter.addressed(coefficient, records[address]['UNIT'], address, None)
                cell_routes.append({'address': address, 'row': index, 'originalIntegral': integral,
                                    'coefficient': coefficient, 'materialValue': value, 'unit': records[address]['UNIT']})
        adapter.begin('direct', [v['address'] for v in cell_routes])
        for entry in cell_routes:
            f.require(same(adapter.bind(entry['coefficient']), entry['materialValue']), 'actual independent native cell call sequence')
        adapter.finish()
        for item in case['grades']['records'].values():
            kind, *address = item['address']; record = item['record']
            if kind in ('factor', 'source'): f.require(set(record['COEFFICIENTS']) <= {(0, 0, 0)}, 'actual grade-free quadrature inputs')
            if kind == 'factor': f.require(same(record['ORIGINAL'], case['binding']['bound']['rows'][address[0]]['factors'][address[1]]['symbolicCoefficient']), 'factor source join')
            if kind == 'source': f.require(same(record['ORIGINAL'], case['binding']['bound']['sources'][0, address[0]]['symbolicAmplitude']), 'source amplitude join')
        for term in case['grades']['termJoins']:
            f.require(same(term['originalIntegral'], material['rows'][term['integralIndex']]['original']), 'full term integral')
            f.require(all(v['factorGrades'][1:] == ((0, 0, 0), (0, 0, 0)) for factor in term['factors'] for v in factor['combinations']), 'complete independent grade convolution')
        f.require(len(cell_routes) == len(case['grades']['termJoins']), 'all native terms included')
        routes = []
        by_reuse = {v['row']: v for v in reused}
        f.require(set(fresh).isdisjoint(by_reuse) and set(fresh) | set(by_reuse) == set(range(len(material['rows']))), 'complete disjoint numerical partition')
        for row in material['rows']:
            sig = h.signature(row, material, case, common)
            if label == BASELINE or row['index'] in fresh:
                owner = (label, row['index']); atlas[owner] = sig
            else:
                reuse = by_reuse[row['index']]; owner = (reuse['fromCase'], reuse['fromRow'])
                f.require(owner in atlas and same(sig, atlas[owner]), 'actual owner full material numerical signature')
            routes.append({'row': row['index'], 'owner': owner, 'signature': sig})
        for i in fresh: f.require(len(material['rows'][i]['limits']) <= 2, 'no new triple momentum integration')
        native_worker = material_native.MaterialProduction(material, data, chart)
        measures = material_native.mass_controls(native_worker, material)
        controls = []
        entry = next(v for v in adapter.entries if v['expression'] != 0 and v['value'] != 0)
        for name, changed in [('unit', dict(entry, unit=(entry['unit'][0]+1, *entry['unit'][1:]))),
                              ('address', dict(entry, address=(*entry['address'][:-1], entry['address'][-1]+100))),
                              ('coefficient', dict(entry, expression=2*entry['expression']))]:
            rejected = False
            try: adapter.addressed(changed['expression'], changed['unit'], changed['address'], changed['grade'])
            except (ValueError, AssertionError, RuntimeError): rejected = True
            f.require(rejected, ('actual saved coefficient route mutation', name))
            controls.append({'kind': name, 'original': entry, 'changed': changed, 'rejected': rejected})
        row = material['rows'][0]; changed = copy.deepcopy(row); v, a, b = changed['limits'][0]
        changed['limits'] = ((v, a, b+1), *changed['limits'][1:])
        f.require(not native.same(h.signature(row, material, case, common), h.signature(changed, material, case, common)), 'changed internal limit rejects')
        controls.append({'kind': 'limit', 'original': row, 'changed': changed, 'rejected': True})
        changed_frame = copy.deepcopy(common); changed_frame['input']['parameters']['omega'] *= 2
        f.require(not native.same(changed_frame, common), 'changed actual frequency rejects common trial reuse')
        controls.append({'kind': 'frequency', 'original': common['input']['parameters'], 'changed': changed_frame['input']['parameters'], 'rejected': True})
        changed_material = dict(material, settings=dict(material['settings'], sourceNodes=material['settings']['sourceNodes']+1))
        f.require(not native.same(h.signature(row, material, case, common), h.signature(row, changed_material, case, common)), 'changed source quadrature order rejects row reuse')
        controls.append({'kind': 'source-order', 'original': material['settings'], 'changed': changed_material['settings'], 'rejected': True})
        f.atomic_pickle(folder/'coefficient-routes.pickle', adapter.entries)
        f.atomic_pickle(folder/'native-cell-routes.pickle', cell_routes)
        f.atomic_pickle(folder/'row-routes.pickle', routes)
        f.atomic_pickle(folder/'measure-controls.pickle', measures)
        f.atomic_pickle(folder/'mutation-controls.pickle', controls)
        layouts = sorted({tuple(v[0] for v in material['rows'][i]['limits']) for i in fresh}, key=lambda v: (len(v), str(v)))
        f.atomic_pickle(folder/'integration-plan.pickle', {'newRows': fresh, 'reusedRows': reused, 'layouts': layouts,
                        'settings': system['settings'], 'nodes': system['nodes'], 'gradeOrigin': old['gradeOrigin'],
                        'fieldUnits': old['fieldUnits'], 'equationUnits': old['equationUnits']})
        summaries[label] = {'rows': len(material['rows']), 'sources': len(material['jets']), 'terms': len(cell_routes),
                            'records': len(case['grades']['records']), 'newRows': fresh, 'reusedRows': len(reused),
                            'layouts': len(layouts), 'coefficientRoutes': len(adapter.entries), 'controls': len(controls)}
        total_controls += len(controls); f.save(base/'case-inventory.json', summaries)
        memo.clear(); del data, case, packets, material, chart, old, adapter, native_worker, routes, cell_routes, records, controls
        gc.collect()
    f.require(sum(v['rows'] for v in summaries.values()) == 300 and sum(v['terms'] for v in summaries.values()) == 647,
              'complete accepted case census')
    result = {'cases': summaries, 'newUnionRows': sum(len(v['newRows']) for v in summaries.values()),
              'reusedCaseRows': sum(v['reusedRows'] for v in summaries.values()), 'controls': total_controls,
              'nativeJoins': joins, 'baselineLayouts': layout_evidence,
              'newQuadratureNodes': 0, 'newMatrices': 0, 'newBasisConstructions': 0, 'newSolves': 0}
    f.save(base/'preflight.json', result)
    return result


def assemble_material(target, case, packets, material, adapter, r, system, chart, rows, old, independent):
    # Coefficients use X=x/d; derivative matrices already act on the common
    # Eulerian trial coefficients, with complete saved field-jet witnesses.
    view = dict(system, nodes=system['nodes']/float(chart.d))
    full_local = {order: sp.ImmutableMatrix(5, 5, lambda i, j: material['bindings']['local', order, i, j]['material'])
                  for order in case['binding']['local']}
    selected = dict(case, binding=dict(case['binding'], local=full_local))
    adapter.begin('graded')
    graded = interior.assemble(case['grades']['records'], case['grades']['termJoins'], adapter, r, view, rows, True)
    adapter.finish()
    f.atomic_pickle(target/'coefficient-matrices.pickle', graded)
    direct_addresses = [('cell', cell['ROW'], cell['COLUMN'], ti)
                        for cell in packets['assembly']['result']['ROWS'] for ti, _ in enumerate(cell['NONLOCAL'])]
    adapter.begin('direct', direct_addresses)
    direct = matrices.direct_cells(selected, packets, adapter, r, view, rows)
    adapter.finish()
    f.atomic_pickle(target/'direct-native-cells.pickle', direct)
    adapter.begin('unsplit')
    unsplit = interior.assemble(case['grades']['records'], case['grades']['termJoins'], adapter, r, view, rows, False)
    adapter.finish()
    f.atomic_pickle(target/'direct-unsplit.pickle', unsplit)
    size = len(system['nodes']); comparisons = {}
    for kind in ('local', 'nonlocal', 'total'):
        comparisons['native_'+kind] = interior.differences(unsplit[kind][0, 0, 0], direct[kind], size)
        comparisons['approved_'+kind] = interior.differences(interior.recombine(graded[kind], adapter.input.origin, case['grades']['generators']), direct[kind], size)
        for grade, array in graded[kind].items():
            comparisons['eulerian_'+kind+str(grade)] = interior.differences(array, old['matrices'][kind][grade], size)
    origin = dict(adapter.input.origin)
    origin.update({case['grades']['generators'][1]: sp.Rational(1, 137), case['grades']['generators'][2]: sp.Rational(1, 911)})
    actual = interior.recombine(graded['total'], origin, case['grades']['generators'])
    comparisons['independent_saved_eulerian'] = interior.differences(actual, independent['total'][0, 0, 0], size)
    full = interior.recombine(graded['total'], adapter.input.origin, case['grades']['generators'])
    omitted = interior.recombine({g: a for g, a in graded['total'].items() if g != (0, 1, 1)}, adapter.input.origin, case['grades']['generators'])
    mutation = interior.differences(omitted, full, size)
    f.atomic_pickle(target/'comparisons.pickle', {'comparisons': comparisons, 'mixedOmission': mutation,
                    'independentFormalPoint': origin, 'independentMaterialMatrix': actual,
                    'scope': 'Saved own-case Eulerian arithmetic controls in the identical common trial frame; not new physical cases.'})
    f.require(max(v['maximumScaledReferenceFrame'] for v in comparisons.values()) < 1e-9, 'complete material operator/cell/common-frame joins')
    if (0, 1, 1) in graded['total'] and np.any(graded['total'][0, 1, 1]):
        f.require(mutation['maximumScaledReferenceFrame'] > 0, 'actual mixed omission responds')
    result = {k: old[k] for k in ('size', 'generators', 'gradeOrigin', 'blockUnits', 'fieldUnits', 'equationUnits', 'dimensionState', 'settings')}
    result.update(matrices=graded, direct=direct, comparisons=comparisons, mixedOmission=mutation,
                  scope=SCOPE, materialCoordinates=view['nodes'])
    f.atomic_pickle(target/'interior-matrices.pickle', result)
    return {'maximumComparison': max(v['maximumScaledReferenceFrame'] for v in comparisons.values()),
            'mixedOmission': mutation['maximumScaledReferenceFrame'], 'nativeTerms': direct['terms']}


def construct(base, manifest, labels, joins):
    plan = json.loads((base/'preflight.json').read_text())
    all_rows = {BASELINE: f.unpickle(base/'bindings/baseline-material-row-matrices.pickle')}
    summaries = {BASELINE: {'rows': len(all_rows[BASELINE]), 'newRows': [], 'newNodes': 0, 'wholeMaterialOperatorReused': True}}
    for label in labels:
        if label == BASELINE: continue
        target = base/'cases'/label; target.mkdir(parents=True)
        data, case, packets, material, chart, fresh, reused = case_inputs(base, label, manifest)
        system = data['system']; settings = system['settings']; size = len(system['nodes'])
        old = f.unpickle(base/'eulerian/cases'/label/'interior-matrices.pickle')
        adapter = SavedCoefficients(case, material, old['gradeOrigin'])
        f.require(modes.same(adapter.entries, f.unpickle(base/'preparation'/label/'coefficient-routes.pickle')), 'accepted exact assembly adapter')
        f.require(fresh == plan['cases'][label]['newRows'], 'accepted planned new rows only')
        rows = {v['row']: all_rows[v['fromCase']][v['fromRow']] for v in reused}; groups = []
        worker = material_native.MaterialProduction(material, data, chart)
        worker.rows = [material['rows'][i] for i in fresh]
        worker.prepare_basis(material['jets'], system['nodes'], settings['sourceBound'], size, settings['sourceNodes'])
        prepared = {k: getattr(worker, k) for k in ('source_nodes', 'source_weights', 'size', 'jet_data', 'amplitudes', 'basis_checks')}
        prepared['original'] = {k: getattr(worker.original, k) for k in ('source_nodes', 'source_weights', 'size', 'jet_data', 'amplitudes')}
        f.atomic_pickle(target/'prepared-source-basis.pickle', prepared)
        width = float(case['binding']['bound']['abel']['width'].subs(data['r'].regulator, settings['regulator']))
        layouts = sorted({tuple(v[0] for v in row['limits']) for row in worker.rows}, key=lambda v: (len(v), str(v)))
        for variables in layouts:
            folder = target/('layout-'+str(len(groups))); folder.mkdir()
            group = worker.matrix_group(variables, settings, case['binding']['bound']['pairs'], width, system['nodes'], folder)
            groups.append(group); rows.update(zip(group['rowIndices'], group['matrices']))
            f.save(target/'layout-inventory.json', {str(i): {'rowIndices': v['rowIndices'], 'nodes': v['nodes'],
                   'sha256': f.digest(target/('layout-'+str(i))/f'layout-{len(v["variables"])}.pickle')} for i, v in enumerate(groups)})
        f.require(set(rows) == set(range(len(material['rows']))), 'complete actual material row arrays')
        f.atomic_pickle(target/'row-matrices.pickle', {'rows': rows, 'newGroups': groups, 'reusedRows': reused, 'settings': settings})
        eulerian = f.unpickle(base/'eulerian/cases'/label/'row-matrices.pickle')['rows']
        differences = {i: {'residual': a-eulerian[i], 'scaled': float(np.max(abs(a-eulerian[i])))/(1+float(np.max(abs(eulerian[i]))))} for i, a in rows.items()}
        f.atomic_pickle(target/'row-comparisons.pickle', differences)
        f.require(max(v['scaled'] for v in differences.values()) < 1e-9, 'all actual material/own-Eulerian rows')
        for v in reused: f.require(np.array_equal(rows[v['row']], all_rows[v['fromCase']][v['fromRow']]), 'every reused full array exact')
        independent = f.unpickle(base/'eulerian/cases'/label/'direct-independent.pickle')
        summary = assemble_material(target, case, packets, material, adapter, data['r'], system, chart, rows, old, independent)
        summary.update(rows=len(rows), newRows=fresh, reusedRows=len(reused), layouts=len(groups),
                       newNodes=sum(v['nodes'] for v in groups), rowMaximum=max(v['scaled'] for v in differences.values()))
        summaries[label] = summary; all_rows[label] = rows; f.save(base/'case-inventory.json', summaries)
        del data, case, packets, material, chart, old, adapter, worker, prepared, eulerian, differences, independent, groups
        gc.collect()
    return {'cases': summaries, 'newNodes': sum(v['newNodes'] for v in summaries.values()),
            'newRows': sum(len(v['newRows']) for v in summaries.values()), 'nativeJoins': joins}


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--run-directory', required=True, type=Path)
    parser.add_argument('--focused', action='store_true'); parser.add_argument('--resume-focused', type=Path)
    args = parser.parse_args(); f.require(args.focused != bool(args.resume_focused), 'focused preparation or accepted-focus production')
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3)); signal.alarm(900); start = time.monotonic()
    base = args.run_directory.resolve(); base.relative_to(f.STORE); base.mkdir(parents=True, exist_ok=False)
    manifest, labels = load(base, args.resume_focused)
    joins = native_joins()
    if not (base/'native-matrix-joins.json').exists(): f.save(base/'native-matrix-joins.json', joins)
    else: f.require(json.loads((base/'native-matrix-joins.json').read_text()) == joins, 'complete accepted native code joins')
    def forbidden(*args, **kwargs): raise RuntimeError('completed scalar/coordinate/current/mode/solve construction disabled')
    engine.NumericalReducedAction.bind = forbidden
    c.bind = forbidden; c.Chart.__init__ = forbidden; c.Chart.bind_image = forbidden; c.Chart.image = forbidden
    c.source.coordinate_change = forbidden
    for key in ('solve', 'inv', 'pinv', 'svd', 'eig', 'eigh', 'eigvals', 'eigvalsh'): setattr(np.linalg, key, forbidden)
    if args.focused:
        material_native.MaterialProduction.prepare_basis = forbidden
        material_native.MaterialProduction.matrix_group = forbidden
        interior.assemble = forbidden; matrices.direct_cells = forbidden
        result = focus(base, manifest, labels, joins)
    else: result = construct(base, manifest, labels, joins)
    hash_check(base, manifest)
    artifacts = {str(p.relative_to(base)): {'sha256': f.digest(p), 'bytes': p.stat().st_size}
                 for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts
                 and p not in (base/'inputs.json', base/'checks.json')}
    checks = {**manifest, 'status': 'COMPLETED_CASE_MATERIAL_MATRIX_INPUTS' if args.focused else 'COMPLETED_CASE_MATERIAL_INTERIORS',
              'result': result, 'nativeJoins': joins, 'artifacts': artifacts, 'wallSeconds': time.monotonic()-start}
    f.save(base/'checks.json', checks); signal.alarm(0); print(json.dumps(checks, indent=2))


if __name__ == '__main__': main()
