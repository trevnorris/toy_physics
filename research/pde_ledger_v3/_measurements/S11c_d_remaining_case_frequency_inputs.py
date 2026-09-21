#!/usr/bin/env python3
"""Saved-source catalogue for the three remaining bounded frequency searches."""
import argparse
import ast
import collections
import gc
import hashlib
import inspect
import json
from pathlib import Path
import resource
import shutil
import signal
import time

import sympy as sp
import S11c_d_remaining_case_bindings as native
import S11c_d_remaining_case_first_jet_sources as source
import S11c_d_frequency_source as frequency
import S11c_d_frequency_chart as chart

f, engine = native.f, native.engine
PLAN = f.M / 'S11c_d_remaining_case_frequency_inputs_plan.md'
SOURCE_CP = f.M / 'S11c_d_remaining_case_first_jet_sources_preflight.json'
FREQUENCY_CP = f.M / 'S11c_d_frequency_source_checkpoint.json'
CHART_CP = f.M / 'S11c_d_frequency_chart_checkpoint.json'
SCOPE = ('Exact saved expression/unit/source-address catalogue for remaining-case frequency work. '
         'Candidate symbolic reuse only: no new frequency binding, derivative, analytic lift, '
         'end continuation, quadrature, operator, inverse, contour point or physical output. '
         'The accepted baseline 32-point search remains closed without a resolved candidate, '
         'not a certified empty spectrum.')


def body(fn):
    return hashlib.sha256(ast.dump(ast.parse(inspect.getsource(fn))).encode()).hexdigest()


def prohibit():
    def forbidden(*args, **kwargs): raise RuntimeError('saved frequency input focus cannot construct scientific operands')
    for module, names in (
        (native, ('load', 'bind_case', 'operator_grades')),
        (native.grades, ('load', 'split', 'check', 'specs', 'term_joins')),
        (frequency, ('load', 'sources', 'end_sources', 'binding_comparison')),
        (chart, ('load', 'root_chart', 'source_chart', 'denominator_chart', 'path_checks', 'end_tables', 'lift'))):
        for name in names:
            if hasattr(module, name): setattr(module, name, forbidden)
    engine.NumericalReducedAction.__init__ = forbidden
    engine.NumericalReducedAction.bind = forbidden
    for name in ('solve', 'inv', 'pinv', 'lstsq', 'svd', 'eig', 'eigh', 'eigvals', 'eigvalsh', 'matrix_rank'):
        if hasattr(native.np.linalg, name): setattr(native.np.linalg, name, forbidden)


def reference(base, manifest, origin, relative, target, expected):
    original = origin / relative; dest = base / target
    f.require(original.is_file() and f.digest(original) == expected, ('accepted reference hash', str(original)))
    f.require(not dest.exists() and not dest.is_symlink(), 'fresh reference path')
    dest.parent.mkdir(parents=True, exist_ok=True); dest.symlink_to(original)
    f.require(dest.resolve() == original.resolve() and f.digest(dest) == expected, 'exact referenced producer bytes')
    manifest['inputPackets'][str(original)] = expected
    manifest['referencedInputs'][target] = {'original':str(original), 'resolvedOriginal':str(original.resolve()), 'sha256':expected, 'bytes':original.stat().st_size}


def load(base):
    manifest = {'runDirectory':str(base), 'sourceFiles':{}, 'inputPackets':{}, 'referencedInputs':{}, 'acceptedOrigins':{}, 'scope':SCOPE}
    checkpoints = {}
    for label, path, status in (
        ('cases', SOURCE_CP, 'ACCEPTED_CASE_FIRST_JET_SOURCE_INPUTS'),
        ('frequency', FREQUENCY_CP, 'PUBLISHED_ANNEX_VERIFIED'),
        ('chart', CHART_CP, 'PUBLISHED_ANNEX_VERIFIED')):
        cp = json.loads(path.read_text()); origin = Path(cp['runDirectory'])
        f.require(cp['status'] == status and f.digest(origin/'checks.json') == cp['checksSha256'], 'accepted producer checkpoint')
        checks = json.loads((origin/'checks.json').read_text())
        f.require(cp['sourceFiles'] == checks['sourceFiles'], 'whole original source inventory')
        for name, digest in checks['sourceFiles'].items():
            f.require(f.digest(f.ROOT/name) == f.digest(origin/'source'/name) == digest, ('current/frozen original helper', name))
            f.require(name not in manifest['sourceFiles'] or manifest['sourceFiles'][name] == digest, 'shared unchanged source')
            manifest['sourceFiles'][name] = digest
        for name, digest in checks['inputPackets'].items():
            f.require(name not in manifest['inputPackets'] or manifest['inputPackets'][name] == digest, 'shared original input')
            manifest['inputPackets'][name] = digest
        for p in (path, origin/'inputs.json', origin/'checks.json'): manifest['inputPackets'][str(p)] = f.digest(p)
        if 'publication' in cp:
            p = f.ROOT / cp['publication']['path']; f.require(p.is_symlink() and f.digest(p) == cp['publication']['sha256'], 'actual accepted published source evidence')
            manifest['inputPackets'][str(p)] = cp['publication']['sha256']
        manifest['acceptedOrigins'][label] = {'checkpoint':str(path), 'checkpointSha256':f.digest(path), 'runDirectory':str(origin), 'checksSha256':cp['checksSha256']}
        checkpoints[label] = (cp, checks, origin)
        reference(base, manifest, origin, 'inputs.json', 'accepted-'+label+'/inputs.json', f.digest(origin/'inputs.json'))
        reference(base, manifest, origin, 'checks.json', 'accepted-'+label+'/checks.json', cp['checksSha256'])
    cp, checks, origin = checkpoints['cases']; labels = tuple(checks['cases'])
    f.require(labels[0] == native.BASELINE and len(labels) == 4, 'all four actual case inputs')
    manifest.update(input=checks['input'], settings=checks['settings'])
    for label in labels:
        for name in (f'accepted-bindings/{label}/case-binding.pickle',
                     *(f'accepted-cases/{label}/{kind}.pickle' for kind in ('reduced-action','actions','assembly','factorization')),
                     f'input-joins/{label}/source-record-pairs.pickle'):
            reference(base, manifest, origin, name, name, checks['artifacts'][name]['sha256'])
    for label, names in (
        ('frequency', ('frequency-sources.pickle','baseline-binding.pickle','frequency-source.pickle',
                       'reference-frequency-pencil.pickle','left-frequency-pencil.pickle','right-frequency-pencil.pickle')),
        ('chart', ('analytic-sources.pickle','root-chart.pickle','denominator-chart.pickle','frequency-chart.pickle'))):
        cp, checks, origin = checkpoints[label]
        for name in names: reference(base, manifest, origin, name, 'accepted-'+label+'/'+name, cp['artifacts'][name]['sha256'])
    for label in ('frequency','chart'):
        old = json.loads((base/('accepted-'+label)/'inputs.json').read_text())
        f.require(old['input'] == manifest['input'] and old['settings'] == manifest['settings'], 'actual unchanged physical input/profile/settings JSON')
    controls_path = f.M/'S11c_d_remaining_case_coordinate_output_checkpoint.json'
    controls = json.loads(controls_path.read_text())
    f.require(controls['status'] == 'PUBLISHED_ANNEX_VERIFIED' and controls['publication']['independentSha256Verified'], 'preceding all-case coordinate output published')
    published = f.ROOT/controls['publication']['path']
    f.require(published.is_symlink() and f.digest(published) == controls['publication']['sha256'], 'actual preceding control transcript')
    manifest['inputPackets'][str(published)] = controls['publication']['sha256']
    for p in (controls_path, f.M/'S11c_d_frequency_contour_refine_checkpoint.json', f.M/'S11c_d_frequency_contour_refine_report.md'):
        reference(base, manifest, p.parent, p.name, 'accepted-decisions/'+p.name, f.digest(p))
    manifest['baselineSearchDecision'] = 'Closed 32-point bounded baseline search without a resolved candidate; no baseline doubling or empty-spectrum assertion.'
    for path in (Path(__file__).resolve(), PLAN, SOURCE_CP, FREQUENCY_CP, CHART_CP,
                 Path(native.__file__), Path(source.__file__), Path(frequency.__file__), Path(chart.__file__),
                 f.ROOT/'directives/S11c_d_NONLINEAR_POLE_CONTRACT.md', f.ROOT/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):
        name = str(path.relative_to(f.ROOT)); digest = f.digest(path)
        f.require(name not in manifest['sourceFiles'] or manifest['sourceFiles'][name] == digest, 'helper pin remains unchanged')
        manifest['sourceFiles'][name] = digest
    manifest['nativeBodies'] = {'literalEquality':body(native.same), 'actualRecordJoin':body(source.require_record),
        'frequencySourceConstructor':body(frequency.sources), 'frequencyCensus':body(frequency.census),
        'analyticSourceConstructor':body(chart.source_chart), 'analyticLift':body(chart.lift)}
    for name, digest in manifest['sourceFiles'].items():
        dest = base/'source'/name; dest.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(f.ROOT/name,dest)
        f.require(f.digest(dest) == digest, 'frozen helper copy')
    f.save(base/'reference-inventory.json', manifest['referencedInputs']); f.save(base/'inputs.json', manifest)
    return manifest, labels


def signature(kind, expression, unit):
    # This is an index, never an equality certificate or a numerical reuse rule.
    return kind, hash(expression), tuple(unit)


def prepare(base, manifest, labels):
    old = f.unpickle(base/'accepted-frequency/frequency-sources.pickle')
    analytic = f.unpickle(base/'accepted-chart/analytic-sources.pickle')
    f.require(set(old['records']) == set(analytic), 'all baseline frequency and analytic source records')
    atlas = {}; baseline_pairs = []
    for key, record in old['records'].items():
        a = analytic[key]
        f.require(native.same(a['address'],record['address']) and native.same(a['originalLive'],record['liveFrequencyAndGrades'])
                  and native.same(tuple(a['unit']),tuple(record['unit'])), 'full saved frequency-to-analytic source/unit join')
        info = {'kind':'accepted-baseline-frequency','case':native.BASELINE,'key':key,'address':record['address']}
        original, unit = record['original'], tuple(record['unit'])
        atlas.setdefault(signature(record['address'][0],original,unit), []).append((info,original,unit))
        baseline_pairs.append((key,record['address'],original,unit,record['liveFrequencyAndGrades'],a['originalLive']))
    f.atomic_pickle(base/'baseline-frequency-chart-pairs.pickle',baseline_pairs)
    f.atomic_pickle(base/'baseline-frequency-variables.pickle', {k:old[k] for k in ('frequency','referenceFrequency','origin')})
    cases = {}; first_basis = None; new_operands = []; original_scalar_records = 0
    for label in labels:
        folder = base/'cases'/label; folder.mkdir(parents=True)
        case = f.unpickle(base/'accepted-bindings'/label/'case-binding.pickle')
        grade, binding = case['grades'], case['binding']
        pairs = f.unpickle(base/'input-joins'/label/'source-record-pairs.pickle')
        f.require(len(pairs) == len(grade['records']) and {v[0] for v in pairs} == set(grade['records']), 'complete accepted physical source-record pairs')
        basis = {'generators':grade['generators'],'fieldUnits':binding['fieldUnits'],'equationUnits':binding['equationUnits'],
                 'settings':binding['settings'],'cutoffBindings':binding['bound']['cutoffBindings']}
        if first_basis is None: first_basis = basis
        f.atomic_pickle(folder/'common-input-pairs.pickle',(basis,first_basis))
        f.require(native.same(basis,first_basis) and json.loads(json.dumps(binding['settings'])) == manifest['settings'],
                  'full actual grade/field/equation/cutoff/settings join; explicit JSON settings encoding only')
        aliases = {}; raw_pairs = []; counts = collections.Counter(); kinds = collections.Counter()
        for key,address,expression,unit,accepted in pairs:
            source.require_record(grade['records'][key],address,expression,unit)
            source.require_record(accepted,address,expression,unit)
            unit = tuple(unit); kind = address[0]; kinds[kind] += 1
            if label == native.BASELINE:
                record = old['records'][key]
                f.require(native.same((address,expression,unit),(record['address'],record['original'],tuple(record['unit']))), 'all baseline actual original addresses/expressions/units')
            bucket = atlas.setdefault(signature(kind,expression,unit),[])
            owner = next((v for v in bucket if native.same(expression,v[1]) and native.same(unit,v[2])),None)
            if owner is None:
                info = {'kind':'uncomputed-frequency-source','case':label,'key':key,'address':address}
                owner = (info,expression,unit); bucket.append(owner); counts['newUnionOperands'] += 1
                new_operands.append(owner)
            elif owner[0]['kind'] == 'accepted-baseline-frequency': counts['baselineFrequencyCandidates'] += 1
            else: counts['additionalNewOperandUses'] += 1
            raw_pairs.append((key,address,expression,unit,owner))
            aliases[key] = {'address':address,'unit':unit,'owner':owner[0],
                            'bindingReuseAccepted':False,'numericalRowReuseAccepted':False}
        f.atomic_pickle(folder/'frequency-source-candidate-pairs.pickle',raw_pairs)
        f.atomic_pickle(folder/'source-routes.pickle',aliases)
        f.atomic_pickle(folder/'term-inputs.pickle',grade['termJoins'])
        f.atomic_pickle(folder/'row-source-profile-inputs.pickle',{'rows':binding['bound']['rows'],'sources':binding['bound']['sources'],
            'profiles':binding['bound']['profiles'],'profileUnits':binding['bound']['profileUnits'],'jets':binding['jets'],
            'abelSourceJoin':binding['abelSourceJoin'],'settings':binding['settings']})
        # Controls alter real saved operands and are checkpointed before guards.
        nonzero = next(v for v in pairs if v[2] != 0); key,address,expr,unit,accepted = nonzero
        row = binding['bound']['rows'][0]; limit = row['limits'][0]
        wrong_address = (address[0],999,*address[2:]); wrong_unit = (unit[0]+1,*unit[1:])
        wrong_limit = sp.Tuple(limit[0],limit[1],limit[2]+1)
        controls_input = {'record':nonzero,'changedExpression':2*expr,'changedAddress':wrong_address,'changedUnit':wrong_unit,
                          'limit':limit,'changedLimit':wrong_limit,'frequency':old['frequency'],'changedFrequency':old['frequency']+1}
        f.atomic_pickle(folder/'mutation-operands.pickle',controls_input)
        controls = {'coefficient':source.rejects(lambda:source.require_record(accepted,address,2*expr,unit)),
            'address':source.rejects(lambda:source.require_record(accepted,wrong_address,expr,unit)),
            'unit':source.rejects(lambda:source.require_record(accepted,address,expr,wrong_unit)),
            'orderedLimit':not native.same(limit,wrong_limit),'frequency':not native.same(old['frequency'],old['frequency']+1)}
        f.require(all(controls.values()), 'actual source/address/unit/ordered-limit/frequency controls respond')
        f.save(folder/'mutation-controls.json',controls)
        cases[label] = {'records':len(pairs),'recordKinds':dict(kinds),'rows':len(binding['bound']['rows']),
            'terms':len(grade['termJoins']),'sources':len(binding['jets']),**{k:counts[k] for k in ('newUnionOperands','baselineFrequencyCandidates','additionalNewOperandUses')},
            'controls':controls,'newNumericalRows':'not classified; requires full live frequency signatures'}
        original_scalar_records += len(pairs); f.save(base/'case-inventory.json',cases)
        del case,grade,binding,pairs,raw_pairs;gc.collect()
    f.atomic_pickle(base/'uncomputed-frequency-operands.pickle',new_operands)
    f.require(sum(v['records'] for v in cases.values()) == original_scalar_records and
              sum(v['newUnionOperands'] for v in cases.values()) == len(new_operands), 'complete actual union operand census')
    return cases,len(new_operands),len(old['records'])


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True);args=parser.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic()
    manifest,labels=load(base);prohibit();cases,new_count,baseline_count=prepare(base,manifest,labels)
    for name,digest in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/name)==f.digest(base/'source'/name)==digest,'all current/frozen helper hashes')
    for name,digest in manifest['inputPackets'].items():f.require(f.digest(Path(name))==digest,'all original input pre/post hashes')
    for name,v in manifest['referencedInputs'].items():
        p=base/name;f.require(p.is_symlink() and str(p.readlink())==v['original'] and str(p.resolve())==v['resolvedOriginal'] and p.stat().st_size==v['bytes'] and f.digest(p)==v['sha256'],'all original reference address/hash/size joins')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*')
               if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    checks={**manifest,'status':'COMPLETED_CASE_FREQUENCY_SOURCE_INPUTS','cases':cases,'acceptedBaselineFrequencyRecords':baseline_count,
        'newUnionSourceOperands':new_count,'newBindings':0,'newFrequencyDerivatives':0,'newAnalyticLifts':0,'newNumericalWork':0,
        'artifacts':artifacts,'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
