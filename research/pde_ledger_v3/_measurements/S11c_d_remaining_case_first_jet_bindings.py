#!/usr/bin/env python3
"""Bind only new literal first-derivative operands; identify exact row reuse."""
import argparse
import ast
import hashlib
import inspect
import json
from pathlib import Path
import resource
import shutil
import signal
import textwrap
import time

import sympy as sp
from sympy.core.function import AppliedUndef
import S11c_d_remaining_case_first_jet_sources as source

f, engine, native, modes = source.f, source.engine, source.native, source.modes
CP = f.M/'S11c_d_remaining_case_first_jet_sources_checkpoint.json'
PLAN = f.M/'S11c_d_remaining_case_first_jet_bindings_plan.md'
BASELINE = source.BASELINE
SCOPE = ('Bindings of the three missing literal first-w-derivative controls. '
         'All baseline control results and unchanged native operands are reused. '
         'Row signatures establish eligibility; saved row arrays and trial bases '
         'must still join before numerical reuse. No new quadrature or response.')


def native_body(function):
    return hashlib.sha256(ast.dump(ast.parse(textwrap.dedent(inspect.getsource(function)))).encode()).hexdigest()


def load(base):
    origin, checks, _ = source.provenance.accepted(CP, 'ACCEPTED_CASE_FIRST_JET_SOURCE_CONTROLS')
    pins = dict(checks['sourceFiles'])
    for path in (Path(__file__).resolve(), PLAN, CP):
        pins[str(path.relative_to(f.ROOT))] = f.digest(path)
    manifest = {'runDirectory': str(base), 'sourceFiles': pins, 'inputPackets': dict(checks['inputPackets']),
                'copiedInputs': {}, 'input': checks['input'], 'settings': checks['settings'], 'scope': SCOPE}
    for name, item in checks['artifacts'].items():
        modes.retain(origin/name, base/name, manifest, item['sha256'])
    modes.retain(origin/'inputs.json', base/'source-control-inputs.json', manifest)
    modes.retain(origin/'checks.json', base/'source-control-checks.json', manifest)
    manifest['settings'] = f.unpickle(base/'accepted-first-jet/first-jet-binding.pickle')['settings']
    source.provenance.restore_settings(checks['settings'], manifest['settings'])
    for name, digest in pins.items():
        p = base/'source'/name; p.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(f.ROOT/name, p); f.require(f.digest(p) == digest, 'binding frozen source identity')
    manifest['nativeBodies'] = {'bind': native_body(engine.NumericalReducedAction.bind),
                               'sourceJets': native_body(f.source_jets),
                               'sourceSubstitution': native_body(engine.dag_substitute)}
    f.save(base/'inputs.json', manifest)
    return manifest, tuple(checks['cases'])


def binding_context(r, adapter, binding):
    # These are every mutable/contextual input consumed by native bind/source_jets
    # and the subsequent cutoff substitution, plus the actual inherited unit frame.
    return {'limits': adapter.input.limits, 'origin': adapter.input.origin,
            'parameters': adapter.input.parameters, 'profiles': adapter.input.profiles,
            'profileFunctions': r.profiles, 'xi': r.xi, 'z': r.z, 'zp': r.zp,
            'regulator': r.regulator, 'probes': tuple(adapter.pencil.probes),
            'cutoffs': binding['bound']['cutoffBindings'], 'tests': binding['bound']['tests'],
            'fieldUnits': binding['fieldUnits'], 'equationUnits': binding['equationUnits'],
            'settings': binding['settings']}


def row_signature(row, binding):
    bound = binding['bound']
    used_profiles = set().union(*(v['coefficient'].atoms(sp.Integral) for v in row['factors']))
    f.require(used_profiles <= set(bound['profileUnits']), 'all consumed row profile units present')
    # BasisMomentum consumes the row coefficients and their actual nested
    # integrals. Other rows' unused cache entries do not enter this quadrature.
    return (native.signature(row, binding),
            tuple(bound['sources'][0, v['sourceIndex']]['boundCharacter'] for v in row['factors']),
            binding['fieldUnits'], binding['equationUnits'],
            {v: bound['profileUnits'][v] for v in used_profiles}, bound['abel'], binding['settings'])


def scalar_entries(binding, records):
    cells = {(v['row'], v['column']): v for v in binding['bound']['cells'] if v['test'] == 0}
    rows = {v['index']: v for v in binding['bound']['rows']}
    for key, item in records.items():
        kind, *address = item['address']; record = item['record']
        if kind == 'local':
            order, i, j = address; value = binding['local'][order][i, j]
        elif kind == 'cell':
            i, j, index = address; value = cells[i, j]['terms'][index][1]
        elif kind == 'factor':
            index, fi = address; value = rows[index]['factors'][fi]['coefficient']
            f.require(native.same(rows[index]['factors'][fi]['symbolicCoefficient'], record['ORIGINAL']),
                      'actual accepted scalar factor source')
        else:
            continue
        yield key, record['ORIGINAL'], record['UNIT'], kind == 'factor', value


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--run-directory', type=Path, required=True)
    args = ap.parse_args(); start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    def timeout(*_):
        raise TimeoutError('first-derivative binding budget; retain all completed pairs and operands')
    signal.signal(signal.SIGALRM, timeout); signal.alarm(900)
    base = args.run_directory.resolve(); base.relative_to(f.STORE); base.mkdir(parents=True, exist_ok=False)
    manifest, labels = load(base)
    old = {label: f.unpickle(base/'accepted-bindings'/label/'case-binding.pickle') for label in labels}
    baseline = f.unpickle(base/'accepted-first-jet/first-jet-binding.pickle')
    contexts = {}; context_pairs = {}; reference = None
    for label, case in old.items():
        r, adapter, packets = source.matrices.context(base, label, case, manifest['input'])
        signature = binding_context(r, adapter, case['binding'])
        context_pairs[label] = signature
        f.atomic_pickle(base/('binding-context-'+label+'.pickle'), signature)
        if reference is None:
            reference = signature
        f.require(native.same(signature, reference), 'actual complete native binder/source-jet context identity')
        contexts[label] = (r, adapter, packets)
    wrong = dict(reference, origin=dict(reference['origin']))
    changed_key = next(k for k, v in wrong['origin'].items() if v != 0)
    wrong['origin'][changed_key] *= 2
    f.require(not native.same(wrong, reference), 'actual material-origin change rejects binding reuse')
    f.atomic_pickle(base/'binding-context-control.pickle', {'actual': reference, 'mutated': wrong})

    baseline_view = {k: baseline[k] for k in ('bound', 'jets', 'local', 'settings')}
    f.require(baseline['changedRecordCounts']['cell'] == 0, 'actual original baseline unchanged cell coefficients')
    baseline_view.update(fieldUnits=old[BASELINE]['binding']['fieldUnits'],
                         equationUnits=old[BASELINE]['binding']['equationUnits'],
                         dimensionState=old[BASELINE]['grades']['dimensionState'])
    # The historical continuum-only baseline control did not consume Gaussian
    # cell-local witnesses. Keep that packet unchanged; the view derives those
    # witnesses from its already accepted changed local matrix without rebinding.
    r = contexts[BASELINE][0]; cells = []
    for cell in baseline_view['bound']['cells']:
        ti, i, j = cell['test'], cell['row'], cell['column']; width, k = baseline_view['bound']['tests'][ti]
        field = sp.exp(-(r.z/width)**2+sp.I*k*r.z)
        cells.append(dict(cell, local=sum(v[i, j]*sp.diff(field, r.z, n) for n, v in baseline_view['local'].items())))
    baseline_view['bound'] = dict(baseline_view['bound'], cells=cells)
    f.atomic_pickle(base/'baseline-first-jet-binding-view.pickle', baseline_view)
    all_views = [('original', label, value['binding'], value['grades']['records'],
                  str(base/'accepted-bindings'/label/'case-binding.pickle')) for label, value in old.items()]
    all_views.append(('first-jet', BASELINE, baseline_view, baseline['records'],
                      str(base/'accepted-first-jet/first-jet-binding.pickle')))
    scalars = {}; source_cache = {}; profile_units = {}; row_cache = {}
    for kind, label, binding, records, path in all_views:
        for key, expression, unit, cut, value in scalar_entries(binding, records):
            scalars.setdefault((hash(expression), tuple(unit), cut), []).append(
                {'expression': expression, 'value': value, 'packet': path, 'case': label, 'key': key, 'kind': kind})
        for si, jet in binding['jets'].items():
            s = binding['bound']['sources'][0, si]; expression = s['symbolicAmplitude']
            proofs = {v['test']: v for v in binding.get('gaussianComparisons', ()) if v['sourceIndex'] == si}
            source_cache.setdefault((hash(expression), tuple(s['amplitudeUnit']), tuple(s['integralUnit'])), []).append(
                {'expression': expression, 'jet': jet, 'proofs': proofs,
                 'origin': {'packet': path, 'case': label, 'sourceIndex': si, 'kind': kind}})
        for integral, unit in binding['bound']['profileUnits'].items():
            f.require(integral not in profile_units or profile_units[integral] == unit, 'accepted nested profile unit joins')
            profile_units[integral] = unit
        for row in binding['bound']['rows']:
            row_cache.setdefault(hash(row['original']), []).append(
                (row_signature(row, binding), {'kind': kind, 'fromCase': label, 'fromRow': row['index'], 'packet': path}))

    summaries = {BASELINE: {'reusedWholeControl': True, 'rows': len(baseline_view['bound']['rows']),
                           'sources': len(baseline_view['jets']), 'newRows': [], 'newBindings': 0}}
    outputs = {BASELINE: str(base/'baseline-first-jet-binding-view.pickle')}
    for label in labels:
        if label == BASELINE:
            continue
        target = base/'bindings'/label; target.mkdir(parents=True, exist_ok=False)
        selected = f.unpickle(base/'cases'/label/'first-jet-sources.pickle')
        r, adapter, packets = contexts[label]
        engine.PHYSICAL_METADATA.dimensions.__dict__.update(selected['dimensionState'])
        prior = old[label]['binding']; cuts = prior['bound']['cutoffBindings']
        bound_values = {}; operand_joins = {}; fresh_count = 0; reused_count = 0; new_jets = 0; new_gaussians = 0
        def scalar(key, expression, unit, cut):
            nonlocal fresh_count, reused_count
            candidates = scalars.get((hash(expression), tuple(unit), cut), ())
            cached = next((v for v in candidates if native.same(v['expression'], expression)), None)
            if cached:
                value = cached['value']; route = {k: v for k, v in cached.items() if k not in ('value', 'expression')}; reused_count += 1
            else:
                value = adapter.bind(expression)
                if cut:
                    value = engine.memo_xreplace(value, cuts)
                route = None; fresh_count += 1
            result = {'expression': expression, 'unit': unit, 'cutoffsApplied': cut, 'value': value,
                      'reusedFrom': route, 'bindingContext': str(base/('binding-context-'+label+'.pickle'))}
            path = target/'bound-operands'/(key+'.pickle'); path.parent.mkdir(exist_ok=True); f.atomic_pickle(path, result)
            operand_joins[key] = {'path': str(path.relative_to(base)), 'sha256': f.digest(path), 'reusedFrom': route}
            if not cached:
                scalars.setdefault((hash(expression), tuple(unit), cut), []).append(
                    {'expression': expression, 'value': value, 'packet': str(path), 'case': label, 'key': key, 'kind': 'completed-new'})
            return value
        addresses = {tuple(v['address']): v['record'] for v in selected['records'].values()}
        for key, item in selected['records'].items():
            kind = item['address'][0]
            if kind not in ('local', 'cell', 'factor'):
                continue
            rec = item['record']; bound_values[tuple(item['address'])] = scalar(key, rec['ORIGINAL'], rec['UNIT'], kind == 'factor')
            if kind == 'factor':
                f.require(set(rec['COEFFICIENTS']) <= {(0, 0, 0)}, 'actual grade-free momentum factor')
        rows = []; profiles = {}; source_routes = []; gaussian = []; jets = {}; sources = {}
        for row in prior['bound']['rows']:
            factors = []
            for fi, factor in enumerate(row['factors']):
                value = bound_values['factor', row['index'], fi]; rec = addresses['factor', row['index'], fi]
                f.require(not value.has(AppliedUndef, sp.Derivative, sp.Subs)
                          and not (engine.dag_free_symbols(value)-{r.z, r.regulator, *(v[0] for v in row['limits'])}),
                          'complete bound actual momentum/source coefficient')
                factors.append(dict(factor, coefficient=value, symbolicCoefficient=rec['ORIGINAL']))
                unknown = value.atoms(sp.Integral)-set(profile_units)
                if unknown:
                    for number, integral in enumerate(sorted(rec['ORIGINAL'].atoms(sp.Integral), key=sp.default_sort_key)):
                        unit = engine.PHYSICAL_METADATA.dimensions.measure(integral)
                        actual = scalar(f'profile{row["index"]}Factor{fi}Integral{number}', integral, unit, True)
                        f.require(actual not in profile_units or profile_units[actual] == unit, 'new actual nested integral unit')
                        profile_units[actual] = unit
                for integral in value.atoms(sp.Integral):
                    f.require(integral in profile_units, 'no missing nested profile unit')
                    profiles[integral] = profile_units[integral]
            actual = dict(row, factors=factors); rows.append(actual)
            path = target/'rows'/f'{row["index"]:03}.pickle'; path.parent.mkdir(exist_ok=True); f.atomic_pickle(path, actual)
        for si, old_jet in prior['jets'].items():
            rec = addresses['source', si]; old_source = prior['bound']['sources'][0, si]
            f.require(set(rec['COEFFICIENTS']) <= {(0, 0, 0)}, 'actual grade-free source amplitude')
            new_source = dict(old_source, symbolicAmplitude=rec['ORIGINAL'])
            lookup = (hash(rec['ORIGINAL']), tuple(new_source['amplitudeUnit']), tuple(new_source['integralUnit']))
            cached = next((v for v in source_cache.get(lookup, ()) if native.same(v['expression'], rec['ORIGINAL'])), None)
            if cached is None:
                jet = f.source_jets(new_source, adapter, r); new_jets += 1
                cached = {'expression': rec['ORIGINAL'], 'jet': jet, 'proofs': {},
                          'origin': {'packet': str(target/'sources'/f'{si:03}.pickle'), 'case': label, 'sourceIndex': si, 'kind': 'completed-new'}}
                source_cache.setdefault(lookup, []).append(cached)
            jet = cached['jet']; jets[si] = jet
            f.require(jet['column'] == old_jet['column'] and native.same(jet['probe'], old_jet['probe']),
                      'same actual field/probe/source position after reversal')
            f.require(jet['residual'] == 0 and jet['amplitudeUnit'] == new_source['amplitudeUnit']
                      and jet['integralUnit'] == new_source['integralUnit'], 'actual complete jet proof and units')
            request = {'sourceIndex': si, 'selected': new_source, 'jet': jet, 'reusedFrom': cached['origin']}
            path = target/'sources'/f'{si:03}.pickle'; path.parent.mkdir(exist_ok=True); f.atomic_pickle(path, request)
            source_routes.append(request)
            for ti, (width, k) in enumerate(prior['bound']['tests']):
                previous_proof = cached['proofs'].get(ti)
                if previous_proof is None:
                    field = lambda z: sp.exp(-(z/width)**2+sp.I*k*z)
                    mapping = {probe: field for probe in adapter.pencil.probes}
                    direct = adapter.bind(engine.dag_substitute(rec['ORIGINAL'], mapping))
                    reconstructed = sum(value*sp.diff(field(r.zp), r.zp, n) for n, value in enumerate(jet['coefficients']))
                    previous_proof = {'sourceIndex': si, 'test': ti, 'direct': direct, 'jet': reconstructed,
                                      'rawResidual': direct-reconstructed, 'expandedResidual': sp.expand(direct-reconstructed)}
                    new_gaussians += 1
                    cached['proofs'][ti] = previous_proof
                proof = {'requestCase': label, 'requestSource': si, 'sourceExpression': rec['ORIGINAL'],
                         'originalSourceRoute': cached['origin'], 'actualProof': previous_proof}
                path = target/'source-pairs'/f'{si:03}-{ti}.pickle'; path.parent.mkdir(exist_ok=True); f.atomic_pickle(path, proof)
                f.require(previous_proof['expandedResidual'] == 0, 'actual direct Gaussian/generic-source residual')
                gaussian.append(dict(previous_proof, sourceIndex=si, originalSourceRoute=cached['origin']))
                sources[ti, si] = dict(prior['bound']['sources'][ti, si], symbolicAmplitude=rec['ORIGINAL'],
                                      boundAmplitude=previous_proof['direct'])
        orders = sorted({a[1] for a in addresses if a[0] == 'local'})
        local = {order: sp.ImmutableMatrix(5, 5, [bound_values['local', order, i, j] for i in range(5) for j in range(5)]) for order in orders}
        cells = []
        for cell in prior['bound']['cells']:
            ti, i, j = cell['test'], cell['row'], cell['column']; width, k = prior['bound']['tests'][ti]
            field = sp.exp(-(r.z/width)**2+sp.I*k*r.z)
            terms = [(index, bound_values['cell', i, j, term]) for term, (index, _) in enumerate(cell['terms'])]
            cells.append(dict(cell, terms=terms, local=sum(value[i, j]*sp.diff(field, r.z, order) for order, value in local.items())))
        rebound = dict(prior['bound'], rows=rows, sources=sources, cells=cells,
                       profiles=[{'bound': p, 'unit': u} for p, u in profiles.items()], profileUnits=profiles)
        actual = dict(prior, bound=rebound, jets=jets, local=local, gaussianComparisons=gaussian,
                      dimensionState=dict(vars(engine.PHYSICAL_METADATA.dimensions)))
        f.atomic_pickle(target/'bound-operands-complete.pickle', actual)
        f.atomic_pickle(target/'scalar-binding-joins.pickle', operand_joins)
        f.atomic_pickle(target/'source-binding-joins.pickle', source_routes)
        # The even Abel source/density and its accepted half-height proof have
        # unchanged complete binder inputs. Check actual new denominators against
        # those retained sampling panels before admitting numerical row reuse.
        moments = {r.normal_map[g[2]] for g in r.momentum_groups}; denominators = {}; seen = set()
        def visit(node):
            if id(node) in seen:
                return
            seen.add(id(node))
            if node.is_Pow and node.exp.is_negative and node.base.has(r.regulator):
                pair = tuple(sorted(engine.dag_free_symbols(node.base) & moments, key=sp.default_sort_key))
                if pair:
                    denominators[node.base] = pair
            for child in node.args:
                visit(child)
        for row in rows:
            for factor in row['factors']:
                visit(factor['coefficient'])
        f.require(set(denominators.values()) <= set(prior['bound']['pairs']), 'actual transfer denominators covered by accepted panels')
        abel = {'acceptedSourceJoin': prior['abelSourceJoin'], 'actualDenominatorPairs': denominators,
                'actualPairs': tuple(sorted(set(denominators.values()), key=str)),
                'bindingContextJoin': str(base/('binding-context-'+label+'.pickle'))}
        f.atomic_pickle(target/'abel-source-join.pickle', abel)
        f.require(prior['abelSourceJoin']['halfHeightResidual'] == 0, 'saved actual Abel half-height proof')
        actual['abelSourceJoin'] = abel
        f.atomic_pickle(target/'binding.pickle', actual)
        reused = []; fresh = []; row_pairs = []
        for row in rows:
            signature = row_signature(row, actual); bucket = row_cache.setdefault(hash(row['original']), [])
            matched = next((route for value, route in bucket if native.same(value, signature)), None)
            if matched:
                reused.append(dict(matched, row=row['index']))
            else:
                fresh.append(row['index']); bucket.append((signature, {'kind': 'first-jet', 'fromCase': label,
                    'fromRow': row['index'], 'packet': str(target/'case-binding.pickle')}))
            row_pairs.append({'row': row['index'], 'signature': signature, 'reusedFrom': matched})
            f.require(len({jets[v['sourceIndex']]['column'] for v in row['factors']}) == 1, 'full actual row input field')
        wrong = dict(rows[0]); limit = wrong['limits'][0]; wrong['limits'] = (sp.Tuple(limit[0], limit[1], limit[2]+1), *wrong['limits'][1:])
        f.require(not native.same(row_signature(wrong, actual), row_signature(rows[0], actual)), 'actual changed ordered limit rejects row reuse')
        wrong = dict(rows[0], factors=[dict(v) for v in rows[0]['factors']]); fi = next(i for i, v in enumerate(wrong['factors']) if v['coefficient'] != 0)
        wrong['factors'][fi]['coefficient'] *= 2
        f.require(not native.same(row_signature(wrong, actual), row_signature(rows[0], actual)), 'actual changed coefficient rejects row reuse')
        f.require(set(fresh).isdisjoint(v['row'] for v in reused) and set(fresh)|{v['row'] for v in reused} == {r['index'] for r in rows}, 'complete actual new/reused partition')
        f.atomic_pickle(target/'row-reuse-pairs.pickle', row_pairs)
        result = {'binding': actual, 'grades': selected, 'newRows': fresh, 'reusedRows': reused,
                  'endpoints': selected['endpointOperands'], 'sourceFiles': manifest['sourceFiles'],
                  'inputPackets': manifest['inputPackets'], 'scope': SCOPE}
        f.atomic_pickle(target/'case-binding.pickle', result); outputs[label] = str(target/'case-binding.pickle')
        summaries[label] = {'rows': len(rows), 'sources': len(jets), 'terms': len(selected['termJoins']),
                            'changedSourceRecordCounts': selected['changedRecordCounts'],
                            'newScalarBindings': fresh_count, 'reusedScalarBindings': reused_count,
                            'newSourceJets': new_jets, 'newGaussianChecks': new_gaussians,
                            'gaussianChecks': len(gaussian), 'nestedProfiles': len(profiles),
                            'newRows': fresh, 'reusedRows': reused, 'rowMutationsRejected': 2}
        f.save(base/'binding-case-inventory.json', summaries)
    f.atomic_pickle(base/'remaining-case-first-jet-bindings.pickle', {'cases': outputs, 'inventory': summaries,
                    'sourceFiles': manifest['sourceFiles'], 'inputPackets': manifest['inputPackets'], 'scope': SCOPE})
    for name, digest in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name) == f.digest(base/'source'/name) == digest, 'binding source/frozen pre/post')
    for name, digest in manifest['inputPackets'].items():
        f.require(f.digest(Path(name)) == digest, 'binding original input pre/post')
    for name, digest in manifest['copiedInputs'].items():
        f.require(f.digest(base/name) == digest, 'binding copied input pre/post')
    artifacts = {str(p.relative_to(base)): {'sha256': f.digest(p), 'bytes': p.stat().st_size}
                 for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts
                 and p not in (base/'inputs.json', base/'checks.json')}
    checks = {**manifest, 'status': 'COMPLETED_CASE_FIRST_JET_BINDINGS', 'cases': summaries, 'casePackets': outputs,
              'newUnionRows': sum(len(v['newRows']) for v in summaries.values()),
              'newQuadratureNodes': 0, 'newModes': 0, 'newResponseSolves': 0,
              'artifacts': artifacts, 'wallSeconds': time.monotonic()-start,
              'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json', checks); signal.alarm(0); print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
