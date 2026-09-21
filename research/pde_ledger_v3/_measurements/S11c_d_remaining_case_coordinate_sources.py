#!/usr/bin/env python3
"""Native material-source images for only the missing saved case operands."""
import argparse
import ast
import copy
import hashlib
import inspect
import json
from pathlib import Path
import resource
import shutil
import signal
import time

import S11c_d_remaining_case_coordinate_inputs as h

f, engine, sp, native, modes = h.f, h.engine, h.sp, h.native, h.modes
coordinate, source = h.coordinate, h.source
CP = f.M/'S11c_d_remaining_case_coordinate_inputs_checkpoint.json'
JCP = f.M/'S11c_d_remaining_case_first_jet_sources_checkpoint.json'
PLAN = f.M/'S11c_d_remaining_case_coordinate_sources_plan.md'
WIRING = f.M/'S11c_d_remaining_case_coordinate_sources_wiring.json'
SCOPE = ('Only missing union coordinate-source images, with saved chart/field jets '
         'and first-jet/density controls reused. No new numerical or physical response.')


def observe_record(record, path):
    strings = {name: sp.srepr(record[name]) for name in
               ('original', 'encoded', 'coordinateImage', 'coordinateReplay', 'coordinateResidual')}
    f.atomic_pickle(path.with_name(path.stem+'-live.pickle'), {
        'strings': strings, 'sha256': {k: hashlib.sha256(v.encode()).hexdigest() for k, v in strings.items()},
        'address': record['address'], 'unit': record['unit']})


def native_record_body():
    tree = ast.parse(inspect.getsource(coordinate.construct)).body[0]
    loop = next(n for n in tree.body if isinstance(n, ast.For)
                and ast.unparse(n.target) == '(key, item)')
    original = copy.deepcopy(loop.body); body = copy.deepcopy(original); changed = []
    for n in body:
        if isinstance(n, ast.Assign) and ast.unparse(n.targets[0]) == 'mutation':
            expected = ast.parse('first_jet_mutation(original, r)', mode='eval').body
            f.require(ast.dump(n.value) == ast.dump(expected), 'actual native mutation routing call')
            n.value = ast.Name(id='saved_mutation', ctx=ast.Load()); changed.append(1)
    observer = ast.parse('observe_record(record, path)').body[0]
    index = next(i for i, n in enumerate(body) if isinstance(n, ast.Expr)
                 and isinstance(n.value, ast.Call) and ast.unparse(n.value.func) == 'f.require')
    body.insert(index, observer)
    reverse = copy.deepcopy(body); reverse.pop(index)
    for n in reverse:
        if isinstance(n, ast.Assign) and ast.unparse(n.targets[0]) == 'mutation':
            n.value = ast.parse('first_jet_mutation(original, r)', mode='eval').body
    f.require(len(changed) == 1 and ast.dump(ast.Module(body=reverse, type_ignores=[]))
              == ast.dump(ast.Module(body=original, type_ignores=[])), 'whole native record-loop reverse AST')
    wrapper = ast.parse('''def create_record(base, r, geometry, jet, material_jets, source, key, item, saved_mutation, progress):
    target = base
    inventory = {}
    records = {}
''').body[0]
    wrapper.body += body + ast.parse('return record').body
    namespace = dict(vars(coordinate), observe_record=observe_record)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[wrapper], type_ignores=[])), __file__, 'exec'), namespace)
    return namespace['create_record'], {'wholeNativeRecordLoopReverseAst': True,
        'savedMutationRoutingEdits': 1, 'liveOperandObserverEdits': 1,
        'nativeConstructAstSha256': source.body(coordinate.construct),
        'nativeCoordinateChangeAstSha256': source.body(coordinate.coordinate_change),
        'nativeSourceJetsAstSha256': source.body(coordinate.linear_source)}


def native_limit_body():
    tree = ast.parse(inspect.getsource(coordinate.construct)).body[0]
    loop = next(n for n in tree.body if isinstance(n, ast.For) and ast.unparse(n.target) == 'term')
    wrapper = ast.parse('def create_limit(term, geometry):\n    limits = []\n').body[0]
    wrapper.body += copy.deepcopy(loop.body) + ast.parse('return limits[0]').body
    f.require(ast.dump(ast.Module(body=wrapper.body[1:-1], type_ignores=[]))
              == ast.dump(ast.Module(body=loop.body, type_ignores=[])), 'whole unchanged native ordered-limit body')
    namespace = dict(vars(coordinate))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[wrapper], type_ignores=[])), __file__, 'exec'), namespace)
    return namespace['create_limit'], {'wholeNativeLimitBodyUnchanged': True}


def load(base):
    old, checks, cp = source.provenance.accepted(CP, 'ACCEPTED_CASE_COORDINATE_SOURCE_INPUTS')
    jr, jc, _ = source.provenance.accepted(JCP, 'ACCEPTED_CASE_FIRST_JET_SOURCE_CONTROLS')
    manifest = {'runDirectory': str(base), 'sourceFiles': dict(checks['sourceFiles']),
                'inputPackets': dict(checks['inputPackets']), 'copiedInputs': {},
                'scope': SCOPE, 'input': checks['input'], 'settings': checks['settings'],
                'acceptedInputs': {'runDirectory': str(old), 'checksSha256': cp['checksSha256'],
                                   'artifacts': checks['artifacts']}}
    for name, digest in jc['sourceFiles'].items():
        f.require(name not in manifest['sourceFiles'] or manifest['sourceFiles'][name] == digest,
                  'same complete native source before first-jet reuse')
        manifest['sourceFiles'][name] = digest
    for name, digest in jc['inputPackets'].items():
        f.require(name not in manifest['inputPackets'] or manifest['inputPackets'][name] == digest,
                  'same original first-jet input')
        manifest['inputPackets'][name] = digest
    for name, item in checks['artifacts'].items():
        modes.retain(old/name, base/name, manifest, item['sha256'])
    modes.retain(old/'inputs.json', base/'accepted-input-manifest.json', manifest)
    modes.retain(old/'checks.json', base/'accepted-input-checks.json', manifest)
    modes.retain(jr/'checks.json', base/'accepted-first-jet-checks.json', manifest)
    for label in checks['cases']:
        if label == native.BASELINE: continue
        for key in f.unpickle(base/'accepted-bindings'/label/'case-binding.pickle')['grades']['records']:
            for kind in ('raw-records', 'mutation-proofs'):
                relative = f'cases/{label}/{kind}/{key}.pickle'
                modes.retain(jr/relative, base/'accepted-first-jet'/relative, manifest,
                             jc['artifacts'][relative]['sha256'])
    for path in (Path(__file__).resolve(), PLAN, WIRING, CP, JCP, Path(h.__file__).resolve()):
        manifest['sourceFiles'][str(path.relative_to(f.ROOT))] = f.digest(path)
    for name, digest in manifest['sourceFiles'].items():
        target = base/'source'/name; target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists(): f.require(f.digest(target) == digest, 'reused source snapshot')
        else: shutil.copyfile(f.ROOT/name, target)
        f.require(f.digest(target) == digest, 'actual frozen current source identity')
    f.save(base/'inputs.json', manifest)
    return manifest, checks


def construct(base, manifest, inputs, progress):
    create_record, record_join = native_record_body(); create_limit, limit_join = native_limit_body()
    wiring = json.loads(WIRING.read_text())
    f.require(wiring['helperSha256'] == f.digest(Path(__file__))
              and wiring['recordJoin'] == record_join and wiring['limitJoin'] == limit_join
              and wiring['actualToolOutcome']['exitCode'] == 0, 'actual completed constructor wiring')
    f.save(base/'constructor-joins.json', {**record_join, **limit_join})
    accepted = f.unpickle(base/'accepted-coordinate/coordinate-source.pickle')
    geometry, jet = accepted['chart'], accepted['fieldJets']
    known = dict(accepted['dimensionState']['known']); cache = {}; limit_cache = []
    for key, value in accepted['records'].items(): cache[native.BASELINE, key] = value
    for value in accepted['orderedLimits']:
        limit_cache.append((value['originalSourceLimit'], value['originalLimits'], value,
                            {'kind': 'accepted-coordinate', 'address': (value['row'], value['column'], value['term'])}))
    summaries = {}; paths = {}; source_serial = 1000; new_total = new_sources = new_limits = 0
    for label, census in inputs['cases'].items():
        target = base/'coordinate-cases'/label; target.mkdir(parents=True)
        if label == native.BASELINE:
            modes.retain(base/'accepted-coordinate/coordinate-source.pickle', target/'coordinate-source.pickle', manifest)
            paths[label] = str(target/'coordinate-source.pickle')
            summaries[label] = {'wholeCoordinatePacketReused': True, 'records': len(accepted['records']),
                                'terms': len(accepted['orderedLimits']), 'newImages': 0, 'newSourceJets': 0}
            continue
        case = f.unpickle(base/'accepted-bindings'/label/'case-binding.pickle'); grade = case['grades']
        reduction = f.unpickle(base/'accepted-cases'/label/'reduced-action.pickle')
        r, dims = f.prior.domain.momentum.source.native.source.restore_context(reduction)
        dims.__dict__.update(f.unpickle(base/'cases'/label/'dimension-state.pickle'))
        for atom, unit in known.items():
            f.require(atom not in dims.known or dims.known[atom] == unit, 'actual reused source-jet unit')
            dims.known[atom] = unit
        material_jets = {(column, order): sp.Symbol(f's11cdMaterialField{column}Jet{order}')
                         for column in range(5) for order in range(jet['maximumOrder']+1)}
        for (column, order), atom in material_jets.items():
            f.require(dims.known[atom] == tuple(v-order*w for v, w in zip(grade['fieldUnits'][column], (1,0,0))),
                      'actual saved full material derivative symbol unit')
        routes = f.unpickle(base/'cases'/label/'source-routes.pickle'); records = {}; aliases = {}
        new_count = source_count = limit_count = 0; controls = []; owners = []
        for key, route in routes.items():
            item = grade['records'][key]; address = item['address']; original = item['record']['ORIGINAL']; unit = item['record']['UNIT']
            owner = route['owner']; owner_key = (owner['case'], owner['key'])
            raw = f.unpickle(base/'accepted-first-jet/cases'/label/'raw-records'/(key+'.pickle'))
            proof = f.unpickle(base/'accepted-first-jet/cases'/label/'mutation-proofs'/(key+'.pickle'))
            f.require(native.same((raw['address'], raw['original'], raw['unit']), (address, original, unit))
                      and proof['involutionResidual'] == 0 and native.same(proof['restored'], original),
                      'actual completed first-jet input and involution reuse')
            if owner_key not in cache:
                f.require(owner_key == (label, key) and owner['kind'] == 'unmapped-source', 'exact new union owner')
                folder = base/'new-records'/label/key; folder.mkdir(parents=True)
                source_jets = {}
                if address[0] == 'source':
                    source_jets[address[1]] = coordinate.linear_source(original, r, source_serial, grade['fieldUnits'])
                    f.atomic_pickle(folder/'source-jets.pickle', source_jets[address[1]])
                    f.save(folder/'source-namespace.json', {'case': label, 'key': key, 'address': address,
                           'internalJetIndex': source_serial, 'physicalSourceIndex': address[1]})
                    f.require(len(source_jets[address[1]]['fields'])-1 <= jet['maximumOrder'], 'complete accepted field-jet order')
                    source_serial += 1; source_count += 1
                elif address[0] == 'local':
                    f.require(address[1] <= jet['maximumOrder'], 'actual local derivative order within saved map')
                value = create_record(folder, r, geometry, jet, material_jets, source_jets, key, item, raw['mutation'], progress)
                cache[owner_key] = value; new_count += 1
                f.require(value['coordinateResidual'] == 0, 'actual new coordinate replay')
                if value['encoded'] != 0:
                    changed = value['coordinateReplay']+value['encoded']
                    controls.append({'key': key, 'original': value['encoded'], 'changed': changed,
                                     'residual': changed-value['encoded'], 'unit': unit})
                    f.atomic_pickle(folder/'changed-coefficient.pickle', controls[-1])
                    f.require(controls[-1]['residual'] != 0, 'actual same-unit coefficient mutation responds')
                known.update(dims.known)
            value = cache[owner_key]
            f.require(native.same((original, unit), (value['original'], value['unit']))
                      and address[0] == value['address'][0]
                      and native.same(raw['mutation'], value['shape']), 'full scalar image/shape/unit owner join')
            records[key] = dict(value, address=address)
            aliases[key] = {'address': address, 'owner': owner, 'ownerAddress': value['address'],
                            'originalUnitPair': (original, unit, value['original'], value['unit'])}
            owners.append(owner_key)
        f.atomic_pickle(target/'record-aliases.pickle', aliases)
        f.atomic_pickle(target/'records.pickle', records)
        limits = []; limit_aliases = []
        for term in grade['termJoins']:
            old = term['remainingLimits']; sl = term['sourceLimit']
            match = next((v for v in limit_cache if native.same((sl, old), v[:2])), None)
            if match is None:
                value = create_limit(term, geometry); limit_count += 1
                path = base/'new-limits'/f'{new_limits+limit_count}.pickle'; path.parent.mkdir(exist_ok=True)
                f.atomic_pickle(path, value)
                match = (sl, old, value, {'kind': 'new-limit', 'case': label,
                                         'address': (term['row'], term['column'], term['term'])})
                limit_cache.append(match)
            value = dict(match[2], row=term['row'], column=term['column'], term=term['term'], integralIndex=term['integralIndex'])
            limits.append(value); limit_aliases.append({'original': term, 'owner': match[3], 'material': value})
        f.atomic_pickle(target/'ordered-limit-aliases.pickle', limit_aliases)
        result = {'chart': geometry, 'fieldJets': jet, 'records': records, 'orderedLimits': limits,
                  'densityAdvection': f.unpickle(base/'cases'/label/'density-operands.pickle'),
                  'fieldUnits': grade['fieldUnits'], 'equationUnits': grade['equationUnits'],
                  'dimensionState': dict(vars(dims)), 'sourceFiles': manifest['sourceFiles'],
                  'inputPackets': manifest['inputPackets'], 'scope': SCOPE}
        f.atomic_pickle(target/'coordinate-source.pickle', result)
        f.require(len(records) == census['records'] and len(limits) == census['terms']
                  and new_count == census['newUnionOperands'], 'complete native case source/term/union census')
        paths[label] = str(target/'coordinate-source.pickle')
        summaries[label] = {'wholeCoordinatePacketReused': False, 'records': len(records),
            'terms': len(limits), 'newImages': new_count, 'newSourceJets': source_count,
            'newLimitImages': limit_count, 'changedCoefficientControls': len(controls),
            'recordReuseUses': len(records)-new_count, 'sourceNamespaceEnd': source_serial}
        new_total += new_count; new_sources += source_count; new_limits += limit_count
        f.save(base/'coordinate-case-inventory.json', summaries)
    f.save(base/'inputs.json', manifest)
    f.require(new_total == sum(x['newUnionOperands'] for x in inputs['cases'].values()), 'only missing union source images')
    f.atomic_pickle(base/'remaining-case-coordinate-sources.pickle', {'cases': paths, 'inventory': summaries,
                    'sourceFiles': manifest['sourceFiles'], 'inputPackets': manifest['inputPackets'], 'scope': SCOPE})
    return summaries, paths, {'newImages': new_total, 'newSourceJets': new_sources, 'newLimitImages': new_limits}


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--run-directory', type=Path, required=True)
    args = parser.parse_args(); base = args.run_directory.resolve(); base.relative_to(f.STORE)
    base.mkdir(parents=True, exist_ok=False); started = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3)); signal.alarm(900)
    def progress(stage):
        with (base/'progress.jsonl').open('a') as stream:
            stream.write(json.dumps({'stage': stage, 'wallSeconds': time.monotonic()-started})+'\n')
    manifest, inputs = load(base); cases, paths, counts = construct(base, manifest, inputs, progress)
    for name, digest in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name) == f.digest(base/'source'/name) == digest, 'actual current/frozen source pre/post join')
    for name, digest in manifest['inputPackets'].items(): f.require(f.digest(Path(name)) == digest, 'original input pre/post identity')
    for name, digest in manifest['copiedInputs'].items(): f.require(f.digest(base/name) == digest, 'copied completed operand identity')
    artifacts = {str(p.relative_to(base)): {'sha256': f.digest(p), 'bytes': p.stat().st_size}
                 for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts
                 and p not in (base/'inputs.json', base/'checks.json')}
    checks = {**manifest, 'status': 'COMPLETED_CASE_COORDINATE_SOURCES', 'cases': cases, 'casePackets': paths,
              **counts, 'newFieldJetMaps': 0, 'newFirstJetMutations': 0, 'newDensityDerivatives': 0,
              'newNumericalWork': 0, 'artifacts': artifacts, 'wallSeconds': time.monotonic()-started}
    f.save(base/'checks.json', checks); signal.alarm(0); print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
