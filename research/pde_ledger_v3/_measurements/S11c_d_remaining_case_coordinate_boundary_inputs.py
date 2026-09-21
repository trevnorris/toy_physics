#!/usr/bin/env python3
"""Join saved material charts to each actual finite/continuum end input."""
import argparse
import gc
import json
from pathlib import Path
import resource
import shutil
import signal
import time

import numpy as np
import scipy.linalg as la
import S11c_d_remaining_case_coordinate_matrices as matrices
import S11c_d_remaining_case_boundary as ends_native

f, h, c, engine, modes = matrices.f, matrices.h, matrices.c, matrices.engine, matrices.modes
BASELINE = matrices.BASELINE
ICP = f.M/'S11c_d_remaining_case_coordinate_matrices_checkpoint.json'
BCP = f.M/'S11c_d_remaining_case_boundary_checkpoint.json'
MCP = f.M/'S11c_d_coordinate_response_checkpoint.json'
PLAN = f.M/'S11c_d_remaining_case_coordinate_boundary_plan.md'
SCOPE = ('Saved full material chart and own-case end/current/phase input joins. '
         'No boundary construction, current contraction, quadrature or response solve. '
         'The historical material response is continuum only; finite controls remain new work.')


def load(base):
    accepted = [matrices.accepted(path, status) for path, status in (
        (ICP, 'ACCEPTED_CASE_MATERIAL_INTERIOR_MATRICES'),
        (BCP, 'ACCEPTED_FOUR_CASE_BOUNDARY_MAPS'), (MCP, 'PUBLISHED_ANNEX_VERIFIED'))]
    (ir, ic, ip), (br, bc, bp), (mr, mc, mp) = accepted
    manifest = {'runDirectory': str(base), 'sourceFiles': {}, 'inputPackets': {},
                'copiedInputs': {}, 'input': ic['input'], 'settings': ic['settings'], 'scope': SCOPE}
    for root, checks, checkpoint in accepted:
        for n, sha in checks['sourceFiles'].items():
            f.require(n not in manifest['sourceFiles'] or manifest['sourceFiles'][n] == sha,
                      'same consumed native/source implementation')
            f.require(f.digest(f.ROOT/n) == f.digest(root/'source'/n) == sha, 'accepted current/frozen source')
            manifest['sourceFiles'][n] = sha
        for n, sha in checks['inputPackets'].items():
            f.require(n not in manifest['inputPackets'] or manifest['inputPackets'][n] == sha,
                      'same complete original input identity')
            manifest['inputPackets'][n] = sha
        modes.retain(root/'checks.json', base/('accepted-'+root.parent.parent.name+'-checks.json'),
                     manifest, checkpoint['checksSha256'])
    labels = tuple(ic['result']['cases'])
    def retain(root, checks, name, destination=None):
        modes.retain(root/name, base/(destination or name), manifest, checks['artifacts'][name]['sha256'])
    for name in ('bindings/chart-state.pickle', 'eulerian/accepted-finite-system.pickle'):
        retain(ir, ic, name)
    for label in labels:
        for name in ('input-state', 'frame', 'profile-cache-pair'):
            retain(ir, ic, 'bindings/contexts/'+label+'/'+name+'.pickle')
        retain(ir, ic, 'bindings/sources/accepted-bindings/'+label+'/case-binding.pickle')
        retain(ir, ic, 'bindings/sources/coordinate-cases/'+label+'/coordinate-source.pickle')
        for name in ('reduced-action', 'actions', 'assembly', 'factorization'):
            retain(ir, ic, 'bindings/sources/accepted-cases/'+label+'/'+name+'.pickle')
        if label != BASELINE:
            retain(ir, ic, 'cases/'+label+'/interior-matrices.pickle', 'interiors/'+label+'/interior-matrices.pickle')
        retain(br, bc, 'boundary-cases/'+label+'/case-boundary.pickle')
        for end in ('left', 'right'):
            for name in ('continuum-boundary', 'finite-boundary', 'channels', 'array-units'):
                retain(br, bc, 'boundary-cases/'+label+'/'+end+'/'+name+'.pickle')
    retain(br, bc, 'accepted-continuum-boundary.pickle')
    for name in ('material-boundary.pickle', 'material-binding.pickle'):
        retain(mr, mc, name, 'baseline-material/'+name)
    for end in ('left', 'right'):
        for name in ('slab', 'bulk'):
            filename = end+'-'+name+'-material-current-tables.pickle'
            retain(mr, mc, filename, 'baseline-material/'+filename)
        original = [Path(p) for p in mc['inputPackets'] if Path(p).name == end+'-material-boundary.pickle']
        f.require(len(original) == 1, 'exact original per-end material producer')
        modes.retain(original[0], base/'baseline-material'/original[0].name,
                     manifest, mc['inputPackets'][str(original[0])])
    original = [p for p in mc['inputPackets'] if Path(p).name == 'continuum-boundary.pickle']
    f.require(len(original) == 1 and mc['inputPackets'][original[0]] ==
              bc['artifacts']['accepted-continuum-boundary.pickle']['sha256'],
              'original native material constructor consumed this exact full boundary source')
    manifest['acceptedInputs'] = {
        'materialInteriors': {'checkpoint': str(ICP), 'checksSha256': ip['checksSha256']},
        'caseBoundaries': {'checkpoint': str(BCP), 'checksSha256': bp['checksSha256']},
        'baselineMaterial': {'checkpoint': str(MCP), 'checksSha256': mp['checksSha256']},
        'baselineOriginalBoundary': {'path': original[0], 'sha256': mc['inputPackets'][original[0]]}}
    for path in (Path(__file__).resolve(), PLAN, ICP, BCP, MCP, Path(ends_native.__file__),
                 Path(matrices.__file__), Path(c.__file__), Path(h.__file__)):
        manifest['sourceFiles'][str(path.relative_to(f.ROOT))] = f.digest(path)
    for n, sha in manifest['sourceFiles'].items():
        target = base/'source'/n; target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(f.ROOT/n, target); f.require(f.digest(target) == sha, 'frozen full input implementation')
    f.save(base/'inputs.json', manifest)
    matrices.hash_check(base, manifest)
    return manifest, labels


def signature(continuum, finite, state, settings, units):
    # Full packets, including candidate dispositions, complete bases, currents,
    # phase maps and diagnostics. No representation or physical field is waived.
    return {'continuum': continuum, 'finite': finite, 'chart': state,
            'settings': settings, 'units': units}


def actual_address(base, label, end, labels):
    f.require(label in labels and end in ('LEFT', 'RIGHT'), 'declared physical case/end address')
    case = f.unpickle(base/'boundary-cases'/label/'case-boundary.pickle')
    f.require(case['case'] == label, 'case packet address identity')
    packet = case['ends'][end]
    f.require(packet['orientation'] == {'LEFT': -1, 'RIGHT': 1}[end], 'physical end orientation')
    return case


def controls(base, label, end, original, labels):
    records = []
    folder = base/'preparation'/label/end.lower()/'controls'; folder.mkdir()
    def check(path, before, after, changed):
        result = modes.same(original, changed)
        records.append({'path': path, 'original': before, 'changed': after, 'rejected': not result})
        f.atomic_pickle(folder/(str(len(records))+'.pickle'), records[-1])
        f.require(not result, ('actual full end/chart input mutation', path))
    def change_array(kind, key):
        old = original[kind][key][(0, 0)]
        changed = old.copy(); changed.flat[0] += 1
        target = dict(original[kind]); target[key] = dict(target[key])
        target[key][(0, 0)] = changed
        check((kind, key, (0, 0), 0), old.flat[0], changed.flat[0], dict(original, **{kind: target}))
    change_array('continuum', 'incomingOriginPhase')
    change_array('continuum', 'insertion')
    clusters = list(original['continuum']['clusters']); cluster = dict(clusters[0])
    old = cluster['K'][(0, 0)]; changed = old.copy(); changed.flat[0] += 1
    cluster['K'] = dict(cluster['K']); cluster['K'][(0, 0)] = changed; clusters[0] = cluster
    check(('continuum', 'clusters', 0, 'K', (0, 0), 0), old.flat[0], changed.flat[0],
          dict(original, continuum=dict(original['continuum'], clusters=clusters)))
    tables = original['continuum']['currentTables']; name = next(iter(tables))
    key = next(k for k, v in tables[name].items() if v != matrices.sp.zeros(*v.shape))
    before = tables[name][key]; after = 2*before
    altered = dict(tables); altered[name] = dict(tables[name]); altered[name][key] = after
    check(('continuum', 'currentTables', name, key), before, after,
          dict(original, continuum=dict(original['continuum'], currentTables=altered)))
    before = original['finite']['incomingBoundaryData']; after = before.copy(); after.flat[0] += 1
    check(('finite', 'incomingBoundaryData', 0), before.flat[0], after.flat[0],
          dict(original, finite=dict(original['finite'], incomingBoundaryData=after)))
    chart = dict(original['chart']); values = dict(chart['values']); before = values['kappa']
    values['kappa'] = before+1; chart['values'] = values
    check(('chart', 'values', 'kappa'), before, values['kappa'], dict(original, chart=chart))
    units = dict(original['units']); fields = list(units['fieldUnits']); old = tuple(fields[0])
    fields[0] = (old[0]+1, *old[1:]); units['fieldUnits'] = fields
    check(('units', 'fieldUnits', 0), old, fields[0], dict(original, units=units))
    rejected = False
    try: actual_address(base, label, end+'_WRONG', labels)
    except ValueError: rejected = True
    records.append({'path': ('physicalAddress', label, end+'_WRONG'), 'rejected': rejected})
    f.atomic_pickle(base/'preparation'/label/end.lower()/'mutation-controls.pickle', records)
    f.require(rejected, 'actual wrong physical end address')
    return len(records)


def prepare(base, manifest, labels):
    state = f.unpickle(base/'bindings/chart-state.pickle')
    system = f.unpickle(base/'eulerian/accepted-finite-system.pickle')
    settings = system['settings']
    f.require(json.loads(json.dumps(settings)) == manifest['settings'] and len(system['nodes']) == 129,
              'approved actual basis and complete settings')
    accepted = f.unpickle(base/'accepted-continuum-boundary.pickle')
    prior = f.unpickle(base/'baseline-material/material-boundary.pickle')
    baseline = actual_address(base, BASELINE, 'LEFT', labels)
    unit_names = ('fieldUnits', 'rowUnits', 'currentUnit')
    units = {n: baseline[n] for n in unit_names}
    f.require(all(modes.same(accepted[n], baseline[n]) and modes.same(prior['ends'][n], baseline[n])
                  for n in unit_names), 'original material and actual case field/equation/current units')
    baseline_proofs = {}
    for end in ('LEFT', 'RIGHT'):
        f.require(modes.same(baseline['ends'][end], accepted['ends'][end]), 'full original end consumed by material baseline')
        saved = f.unpickle(base/'baseline-material'/(end.lower()+'-material-boundary.pickle'))
        f.require(modes.same(saved['material'], prior['materialEnds'][end]) and
                  modes.same(saved['commonEulerian'], prior['ends']['ends'][end]) and
                  modes.same(saved['proofs'], prior['proofs'][end]), 'entire original end packet producer join')
        material = saved['material']; common = saved['commonEulerian']; g = state['values']
        f.require(material['position'] == accepted['ends'][end]['orientation']*settings['sourceBound']/float(g['d'])
                  and material['tangentialMeasure'] == float(g['A']), 'actual chart endpoint and Fourier measure')
        f.require(modes.same(common['currentTables'], accepted['ends'][end]['currentTables']),
                  'all original source current tables retained with common maps')
        table_joins = {}
        for name in ('slab', 'bulk'):
            path = base/'baseline-material'/(end.lower()+'-'+name+'-material-current-tables.pickle')
            table = f.unpickle(path)
            f.require(modes.same(table['original'], accepted['ends'][end]['currentTables'][name]),
                      'entire actual baseline current-table source')
            f.require(table['piolaNormal'] == table['tangentialMeasure'] == g['A'] and
                      modes.same(table['fieldComponents'], g['B']), 'actual saved Piola/field map')
            kl, kr, ql, qr = table['originalArguments']; ml, mr, ql1, qr1 = table['materialArguments']
            f.require(ql == ql1 and qr == qr1 and table['covectors'] ==
                      {kl: (ml-g['kappa'])/g['d'], kr: (mr-g['kappa'])/g['d']},
                      'complete actual normal-covector chart and unaltered depth arguments')
            table_joins[name] = {'path': str(path), 'sha256': f.digest(path), 'entries': len(table['original'])}
        maximum = c.boundary.norm(saved['proofs']['differences'])
        mutation = c.boundary.norm(saved['proofs']['omittedTangentialMeasure'])
        f.require(maximum < 1e-8 and mutation > 0, 'accepted full material end/current/phase and omitted measure evidence')
        baseline_proofs[end] = {'currentTables': table_joins, 'maximumResidual': maximum,
                                'omittedMeasureMaximum': mutation}
    f.save(base/'baseline-material-reuse.json', baseline_proofs)
    owners = {end: signature(baseline['ends'][end], baseline['finite'][end], state, settings, units)
              for end in ('LEFT', 'RIGHT')}
    inventory = {}; case_inventory = {}; new_families = {}
    for label in labels:
        target = base/'preparation'/label; target.mkdir(parents=True)
        case = actual_address(base, label, 'LEFT', labels)
        frame = f.unpickle(base/'bindings/contexts'/label/'frame.pickle')
        binding = f.unpickle(base/'bindings/sources/accepted-bindings'/label/'case-binding.pickle')
        coordinate = f.unpickle(base/'bindings/sources/coordinate-cases'/label/'coordinate-source.pickle')
        f.require(modes.same(frame['geometry'], state['values']['g']) and
                  modes.same(frame['fieldJets'], state['values']['jets']) and
                  modes.same(frame['geometry'], coordinate['chart']) and
                  modes.same(frame['fieldJets'], coordinate['fieldJets']), 'own source chart and full field-jet identity')
        f.require(modes.same(frame['settings'], settings) and modes.same(binding['binding']['settings'], settings),
                  'own native source/binding settings')
        f.require(modes.same(frame['fieldUnits'], case['fieldUnits']) and
                  modes.same(frame['equationUnits'], case['rowUnits']) and
                  all(modes.same(case[k], units[k]) for k in unit_names), 'complete own field/equation/current coordinates')
        if label != BASELINE:
            interior = f.unpickle(base/'interiors'/label/'interior-matrices.pickle')
            f.require(modes.same(interior['fieldUnits'], case['fieldUnits']) and
                      modes.same(interior['equationUnits'], case['rowUnits']) and
                      modes.same(interior['settings'], settings) and interior['size'] == 129,
                      'actual material interior and own boundary coordinate inputs')
            del interior
        raw = {'frame': frame, 'chart': state, 'sourceChart': coordinate['chart'],
               'sourceFieldJets': coordinate['fieldJets'], 'units': {n: case[n] for n in unit_names}}
        f.atomic_pickle(target/'chart-input-pairs.pickle', raw)
        for end in ('LEFT', 'RIGHT'):
            folder = target/end.lower(); folder.mkdir()
            source = base/'boundary-cases'/label/end.lower()
            current = signature(case['ends'][end], case['finite'][end], state, settings, units)
            for name, expected in (('continuum-boundary', case['ends'][end]), ('finite-boundary', case['finite'][end])):
                f.require(modes.same(f.unpickle(source/(name+'.pickle')), expected), 'entire case/end packet address join')
            channels = f.unpickle(source/'channels.pickle')
            f.require(np.array_equal(channels['OUTWARD_CURRENT'], current['finite']['current']),
                      'actual full native finite open-current form')
            f.require(len(channels['CANDIDATES']) == len(current['continuum']['census']), 'complete root candidate census')
            f.atomic_pickle(folder/'end-inputs.pickle', current)
            flags = {name: modes.same(current, operand) for name, operand in owners.items()}
            f.atomic_pickle(folder/'source-reuse-pairs.pickle', {'actual': current, 'owners': owners, 'equal': flags})
            matches = [name for name, value in flags.items() if value]
            f.require(len(matches) <= 1, 'unique full physical material-end family')
            if matches:
                owner = matches[0]; disposition = 'accepted-baseline-material-end' if owner in ('LEFT', 'RIGHT') else 'same-new-full-input-family'
            else:
                owner = label+'__'+end; owners[owner] = current; new_families[owner] = str(folder/'end-inputs.pickle')
                disposition = 'new-material-end-needed'
            controls_count = controls(base, label, end, current, labels)
            inventory[label+'__'+end] = {'owner': owner, 'disposition': disposition,
                'candidates': len(channels['CANDIDATES']), 'continuumDirections': int(current['continuum']['offsets'][-1]),
                'finiteIncoming': len(current['finite']['incoming']), 'finiteOutgoing': len(current['finite']['outgoing']),
                'finiteCurrentShape': list(current['finite']['current'].shape), 'controls': controls_count,
                'inputSha256': f.digest(folder/'end-inputs.pickle'), 'unitPacketSha256': f.digest(source/'array-units.pickle')}
            f.save(base/'end-inventory.json', inventory)
        case_inventory[label] = {'rows': len(binding['binding']['rows']),
            'terms': len(binding['grades']['termJoins']), 'sourceAmplitudes': len(binding['binding']['jets']),
            'endAddresses': [label+'__'+end for end in ('LEFT', 'RIGHT')]}
        f.save(base/'case-inventory.json', case_inventory)
        del case, frame, binding, coordinate, raw, channels, current
        gc.collect()
    f.require(sum(v['rows'] for v in case_inventory.values()) == 300 and
              sum(v['terms'] for v in case_inventory.values()) == 647, 'all actual physical row and native-cell addresses')
    result = {'cases': case_inventory, 'ends': inventory, 'newMaterialFamilies': new_families,
              'baselineReuse': baseline_proofs, 'newBoundaryConstructions': 0, 'newCurrentContractions': 0,
              'newFiniteResponses': 0, 'newContinuumResponses': 0, 'baselineMaterialControlIsContinuumOnly': True}
    f.save(base/'preflight.json', result)
    return result


def disable_constructors():
    def forbidden(*args, **kwargs): raise RuntimeError('scientific construction prohibited in saved material-boundary input focus')
    for cls in (engine.NumericalReducedAction, engine.ModalCurrentSubspaces, engine.TwoEndedMatchingChannels):
        cls.__init__ = forbidden
    engine.NumericalReducedAction.bind = forbidden
    c.Chart.__init__ = forbidden; c.Chart.maps = forbidden; c.Chart.image = forbidden; c.Chart.basis = forbidden
    c.material_ends = forbidden; c.source.coordinate_change = forbidden; c.bind = forbidden
    c.MaterialMomentum.prepare_basis = forbidden; matrices.interior.assemble = forbidden
    matrices.matrices.direct_cells = forbidden; f.boundary_map = forbidden; f.source_jets = forbidden
    c.boundary.current_pair = forbidden; c.boundary.construct_end = forbidden
    c.response.systems = forbidden; c.response.solve = forbidden; c.response.channels = forbidden; c.response.open_flux = forbidden
    for module in (np.linalg, la):
        for name in ('solve', 'inv', 'pinv', 'svd', 'eig', 'eigh', 'eigvals', 'eigvalsh', 'lu_factor', 'lu_solve'):
            if hasattr(module, name): setattr(module, name, forbidden)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--run-directory', type=Path, required=True); args = ap.parse_args()
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3)); signal.alarm(900); start = time.monotonic()
    base = args.run_directory.resolve(); base.relative_to(f.STORE); base.mkdir(parents=True, exist_ok=False)
    joins = {name: matrices.body(fun) for name, fun in (
        ('material_ends', c.material_ends), ('Chart', c.Chart), ('boundary_map', f.boundary_map),
        ('current_pair', c.boundary.current_pair), ('array_units', ends_native.array_units),
        ('restore', h.restore), ('restore_chart', h.restore_chart))}
    disable_constructors(); manifest, labels = load(base)
    f.save(base/'native-boundary-joins.json', joins)
    result = prepare(base, manifest, labels); matrices.hash_check(base, manifest)
    artifacts = {str(p.relative_to(base)): {'sha256': f.digest(p), 'bytes': p.stat().st_size}
        for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts
        and p not in (base/'inputs.json', base/'checks.json')}
    checks = {**manifest, 'status': 'COMPLETED_CASE_MATERIAL_BOUNDARY_INPUTS', 'result': result,
              'nativeJoins': joins, 'artifacts': artifacts, 'wallSeconds': time.monotonic()-start}
    f.save(base/'checks.json', checks); signal.alarm(0); print(json.dumps(checks, indent=2))


if __name__ == '__main__': main()
