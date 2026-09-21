#!/usr/bin/env python3
"""Join saved live sources to the actual accepted scalar analytic chart."""
import argparse
import ast
import collections
import copy
import gc
import hashlib
import inspect
import json
from pathlib import Path
import resource
import shutil
import signal
import time

import S11c_d_remaining_case_frequency_characters as characters

source = characters.source
f, native, engine, sp, q = source.f, source.native, source.engine, source.sp, source.q
chart = source.inputs.chart
CP = f.M/'S11c_d_remaining_case_frequency_characters_checkpoint.json'
PLAN = f.M/'S11c_d_remaining_case_frequency_analytic_inputs_plan.md'
SCOPE = ('Saved live-source/analytic-image input and scalar-domain routes only. '
         'No lift, differentiation, certificate reconstruction, chart/end continuation, '
         'binding, quadrature, matrix, current, mode, solve or emission. '
         'Numerical row reuse and own-case outgoing end domains remain pending.')


def protect_references(base):
    old_pickle, old_json = f.atomic_pickle, f.save
    def checked(path):
        path = Path(path)
        relative = path.relative_to(base)
        current = base
        for part in relative.parts:
            current = current/part
            f.require(not current.is_symlink(), ('write through reference prohibited',str(current)))
        return path
    def save_pickle(path,value): return old_pickle(checked(path),value)
    def save_json(path,value): return old_json(checked(path),value)
    f.atomic_pickle, f.save = save_pickle, save_json


def load(base):
    cp = json.loads(CP.read_text()); origin = Path(cp['runDirectory'])
    vr = Path(cp['validation']['runDirectory'])
    f.require(cp['status'] == 'ACCEPTED_CASE_FREQUENCY_CHARACTER_INPUTS'
              and f.digest(origin/'checks.json') == cp['checksSha256'], 'accepted character inputs')
    source.receipts.inspect_guard(vr, 'validate')
    f.require(f.digest(vr/'checks.json') == cp['validation']['checksSha256']
              and (vr/'checks.json').read_bytes() == (vr/'validate.stdout').read_bytes(), 'accepted final validator')
    manifest = {'runDirectory':str(base), 'sourceFiles':dict(cp['sourceFiles']),
                'inputPackets':dict(cp['inputPackets']), 'referencedInputs':{},
                'input':cp['input'], 'settings':cp['settings'], 'scope':SCOPE,
                'acceptedCharacters':{'checkpoint':str(CP), 'runDirectory':str(origin),
                    'checksSha256':cp['checksSha256'], 'validatorChecksSha256':cp['validation']['checksSha256']}}
    for name, item in cp['artifacts'].items():
        source.reference(base, manifest, origin/name, name, item['sha256'])
    for path, name in ((CP,'accepted-character-checkpoint.json'),
                       (origin/'checks.json','accepted-character-checks.json'),
                       (origin/'inputs.json','accepted-character-inputs.json'),
                       (vr/'checks.json','accepted-character-validation.json')):
        source.reference(base, manifest, path, name, f.digest(path))
    for path in (Path(__file__).resolve(), PLAN, CP):
        name = str(path.relative_to(f.ROOT)); value = f.digest(path)
        f.require(name not in manifest['sourceFiles'] or manifest['sourceFiles'][name] == value, 'unchanged source pin')
        manifest['sourceFiles'][name] = value
    for name, value in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name) == value, ('current source',name))
        if name in cp['sourceFiles']:
            f.require(f.digest(origin/'source'/name) == value, 'accepted frozen source')
        dest = base/'source'/name; dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(f.ROOT/name,dest); f.require(f.digest(dest) == value, 'frozen source')
    for name, value in manifest['inputPackets'].items():
        f.require(f.digest(Path(name)) == value, ('original input prehash',name))
    f.save(base/'inputs.json',manifest)
    return manifest, tuple(cp['cases'])


def native_join():
    tree = ast.parse(inspect.getsource(chart.source_chart))
    fields = set()
    for node in ast.walk(tree):
        if (isinstance(node,ast.Subscript) and isinstance(node.value,ast.Name)
                and node.value.id == 'old' and isinstance(node.slice,ast.Constant)):
            fields.add(node.slice.value)
    f.require(fields == {'address','liveFrequencyAndGrades','unit','acceptedBinding',
                         'controlFrequency','controlExpected','census'}, 'complete native per-record input fields')
    return {'wholeNativeSourceChartSha256':source.inputs.body(chart.source_chart),
            'wholeNativeLiftSha256':source.inputs.body(chart.lift),
            'wholeNativeRootChartSha256':source.inputs.body(chart.root_chart),
            'wholeNativeDenominatorChartSha256':source.inputs.body(chart.denominator_chart),
            'wholeNativeSeedComparisonSha256':source.inputs.body(chart.seed_comparison),
            'wholeNativeSeedRootIdentitySha256':source.inputs.body(chart.seed_root_identity),
            'consumedRecordFields':sorted(fields), 'nativeBodiesUnchanged':True,
            'nativeConstructorsCalled':False,
            'addressRouting':'Native record address is retained separately for every case and owner; no physical address is discarded.'}


def prohibit():
    characters.prohibit()
    def forbidden(*args,**kwargs):
        raise RuntimeError('saved analytic inputs cannot construct scientific operands')
    for name in ('seed_comparison','seed_root_identity','path_checks','emit_result','main'):
        setattr(chart,name,forbidden)
    engine.BoundedSourceFourierAssembly.reconstruction_certificate = forbidden
    q.binding_comparison = q.census = forbidden
    characters.load = characters.prepare = characters.row_constructor = forbidden


def image_input(record, variables):
    # These are all physical operands read by the native source_chart loop.
    # Its output address is an explicit independent alias, never an equality waiver.
    return {'liveFrequencyAndGrades':record['liveFrequencyAndGrades'], 'unit':record['unit'],
            'acceptedBinding':record['acceptedBinding'], 'controlFrequency':record['controlFrequency'],
            'controlExpected':record['controlExpected'], 'integralLimits':record['census']['integralLimits'],
            'variables':variables}


def index_key(value):
    return (hash(value['liveFrequencyAndGrades']),tuple(value['unit']),
            hash(value['acceptedBinding']),hash(value['controlExpected']))


def require_same(actual, expected):
    f.require(native.same(actual,expected), 'literal full analytic input/domain/address route')


def prepare(base, manifest, labels):
    accepted = base/'accepted'
    old = f.unpickle(accepted/'accepted-frequency/frequency-sources.pickle')
    analytic = f.unpickle(accepted/'accepted-chart/analytic-sources.pickle')
    roots = f.unpickle(accepted/'accepted-chart/root-chart.pickle')
    denominators = f.unpickle(accepted/'accepted-chart/denominator-chart.pickle')
    original_chart = f.unpickle(accepted/'accepted-chart/frequency-chart.pickle')
    context = f.unpickle(accepted/'contexts'/native.BASELINE/'native-input.pickle')
    variables = {k:old[k] for k in ('frequency','referenceFrequency','origin')}
    raw_chart = {'rootPacket':roots, 'savedChartRoots':original_chart['chart'],
                 'denominatorPacket':denominators, 'savedChartDenominators':original_chart['denominators'],
                 'variables':variables, 'dimensionState':original_chart['dimensionState'],
                 'context':context, 'scope':roots['scope']}
    f.atomic_pickle(base/'analytic-chart-input-pairs.pickle',raw_chart)
    require_same(roots,original_chart['chart']); require_same(denominators,original_chart['denominators'])
    f.require(roots['frequency'] == variables['frequency'] and roots['center'] == variables['referenceFrequency'],
              'actual common frequency and seed chart')
    known = original_chart['dimensionState']['known']
    f.require(tuple(known[variables['frequency']]) == (0,-1,0), 'actual formal frequency unit')
    for radicand, root in roots['roots'].items():
        f.require(native.same(radicand,root['radicand']) and root['momentum'].is_real is True
                  and tuple(known[root['carrier']]) == tuple(root['unit']) == (0,-1,0), 'actual chart radical/unit/domain address')
    f.require(set(old['records']) == set(analytic), 'all accepted analytic source records')
    atlas = {}; saved_pairs = []
    for key, record in old['records'].items():
        a = analytic[key]; value = image_input(record,variables)
        pair = {'key':key, 'input':value, 'source':record, 'analytic':a}
        saved_pairs.append(pair)
        f.require(native.same((a['address'],a['originalLive'],a['unit'],a['limits']),
                             (record['address'],value['liveFrequencyAndGrades'],value['unit'],value['integralLimits'])),
                  'actual baseline analytic record address/live/unit/limits')
        for name, expected in (('seed',record['acceptedBinding']),('control',record['controlExpected'])):
            f.require(native.same(a['bindingComparisons'][name]['operands']['right'],expected),
                      'actual saved analytic comparison input')
        owner = {'kind':'accepted-baseline-analytic', 'case':native.BASELINE, 'key':key, 'address':a['address']}
        atlas.setdefault(index_key(value),[]).append({'owner':owner,'input':value})
    f.atomic_pickle(base/'baseline-analytic-input-pairs.pickle',saved_pairs)
    del saved_pairs, original_chart
    pending = []; pending_roots = {}; pending_denominators = {}; cases = {}; owner_records = {}
    first_basis = None
    for label in labels:
        target = base/'analytic-input-cases'/label; target.mkdir(parents=True)
        packet = f.unpickle(accepted/'frequency-cases'/label/'frequency-source.pickle')
        chars = f.unpickle(base/'character-cases'/label/'frequency-characters.pickle')
        own_context = f.unpickle(accepted/'contexts'/label/'native-input.pickle')
        basis = {k:packet[k] for k in ('generators','fieldUnits','equationUnits')}
        if first_basis is None: first_basis = basis
        common = {'contextPair':(own_context,context), 'basisPair':(basis,first_basis),
                  'variablesPair':({k:packet[k] for k in variables},variables),
                  'characterVariablesPair':({k:chars[k] for k in variables},variables),
                  'chartSource':str(accepted/'accepted-chart/root-chart.pickle'),
                  'chartSha256':f.digest(accepted/'accepted-chart/root-chart.pickle')}
        f.atomic_pickle(target/'common-input-pairs.pickle',common)
        for name in ('contextPair','basisPair','variablesPair','characterVariablesPair'):
            require_same(*common[name])
        f.require(chars['summary']['uncomputedCharacterBindings'] == 0
                  and chars['summary']['rowCensusComplete'], 'actual complete character/row inputs')
        rows = f.unpickle(base/'character-cases'/label/'full-row-input-pair.pickle')
        terms = f.unpickle(base/'character-cases'/label/'term-inputs.pickle')
        f.atomic_pickle(target/'physical-row-term-character-inputs.pickle',
                        {'rows':rows, 'terms':terms, 'characters':chars,
                         'source':str(accepted/'frequency-cases'/label/'frequency-source.pickle')})
        routes = {}; domain_routes = {}; counts = collections.Counter()
        for key, record in packet['records'].items():
            own = record['owner']; where = (own['kind'],own['case'],own['key'])
            f.require(own['kind'] in ('accepted-baseline-frequency','uncomputed-frequency-source'), 'actual source owner kind')
            if where not in owner_records:
                owner_records[where] = (old['records'][own['key']] if own['kind'] == 'accepted-baseline-frequency'
                    else f.unpickle(accepted/'new-records'/own['case']/'records'/(own['key']+'.pickle')))
            producer = owner_records[where]
            expected = dict(producer,address=record['address'],ownerAddress=producer['address'],owner=own)
            value = image_input(record,variables)
            bucket = atlas.setdefault(index_key(value),[])
            match = next((v for v in bucket if native.same(v['input'],value)),None)
            if match is None:
                owner = {'kind':'uncomputed-analytic-image','case':label,'key':key,'address':record['address']}
                match = {'owner':owner,'input':value}; bucket.append(match)
                pending.append({'owner':owner,'record':record,'input':value,'context':own_context,
                                'basis':basis,'chartPath':common['chartSource'],'chartSha256':common['chartSha256']})
                counts['newUnionImages'] += 1
            elif match['owner']['kind'] == 'accepted-baseline-analytic': counts['acceptedAnalyticUses'] += 1
            else: counts['additionalNewImageUses'] += 1
            route = {'case':label,'key':key,'address':record['address'],'sourceOwner':own,
                     'sourceOwnerAddress':record['ownerAddress'],'input':value,'analyticOwner':match['owner'],
                     'chartSha256':common['chartSha256'],'numericalRowReuse':False}
            radical_pairs = []; denominator_pairs = []
            for power in record['census']['fractionalPowers']:
                saved = roots['roots'].get(power.base)
                radical_pairs.append({'power':power,'radicand':power.base,'savedRoot':saved})
                if saved is None: pending_roots.setdefault(power.base,[]).append((label,key,record['address']))
            for denominator in record['census']['literalDenominatorBases']:
                saved = next((v for v in denominators['records'] if native.same(v['original'],denominator)),None)
                denominator_pairs.append({'original':denominator,'savedDenominator':saved})
                if saved is None: pending_denominators.setdefault(denominator,[]).append((label,key,record['address']))
            domain = {'case':label,'key':key,'address':record['address'], 'census':record['census'],
                      'radicals':radical_pairs,'denominators':denominator_pairs,
                      'scope':roots['scope'],'chartSha256':common['chartSha256']}
            folder = target/'records'; folder.mkdir(exist_ok=True)
            f.atomic_pickle(folder/(key+'.pickle'),{'source':record,'ownerSource':producer,'route':route,'domain':domain})
            require_same(record,expected)
            f.require(native.same(own['address'],producer['address']), 'actual source owner address')
            for pair in radical_pairs:
                if pair['savedRoot'] is not None: require_same(pair['radicand'],pair['savedRoot']['radicand'])
            for pair in denominator_pairs:
                if pair['savedDenominator'] is not None: require_same(pair['original'],pair['savedDenominator']['original'])
            routes[key] = route; domain_routes[key] = domain
        f.atomic_pickle(target/'analytic-routes.pickle',routes)
        f.atomic_pickle(target/'domain-routes.pickle',domain_routes)
        # Change real operands passed through the literal join, not stand-in booleans.
        key = next(k for k,v in routes.items() if v['input']['liveFrequencyAndGrades'] != 0)
        good = routes[key]; changed = []
        for name in ('liveFrequencyAndGrades','acceptedBinding','controlExpected','controlFrequency'):
            bad = copy.deepcopy(good); bad['input'][name] += 1; changed.append((name,bad))
        bad = copy.deepcopy(good); unit = bad['input']['unit']; bad['input']['unit'] = (unit[0]+1,*unit[1:]); changed.append(('unit',bad))
        bad = copy.deepcopy(good); bad['address'] = ('wrong-address',*bad['address']); changed.append(('address',bad))
        bad = copy.deepcopy(good); bad['chartSha256'] = '0'*64; changed.append(('chart',bad))
        key_with_limit = next(k for k,v in routes.items() if v['input']['integralLimits'])
        limit_good = routes[key_with_limit]; limit_bad = copy.deepcopy(limit_good)
        limits = limit_bad['input']['integralLimits']; limit_bad['input']['integralLimits'] = tuple(reversed(limits))+limits
        root = next(iter(roots['roots'].values())); denominator = denominators['records'][0]
        domain_changes = [('radicand',root['radicand'],root['radicand']+1),
                          ('chartRadius',roots['radius'],roots['radius']+1),
                          ('denominator',denominator['original'],denominator['original']+1)]
        f.atomic_pickle(target/'mutation-operands.pickle',{'original':good,'changed':changed,
                        'orderedLimitOriginal':limit_good,'orderedLimitChanged':limit_bad,'domainChanges':domain_changes})
        controls = {name:source.inputs.source.rejects(lambda bad=bad:require_same(bad,good)) for name,bad in changed}
        controls['orderedLimit'] = source.inputs.source.rejects(lambda:require_same(limit_bad,limit_good))
        controls.update({name:source.inputs.source.rejects(lambda a=a,b=b:require_same(a,b)) for name,a,b in domain_changes})
        f.save(target/'mutation-controls.json',controls); f.require(all(controls.values()), 'actual analytic source/unit/domain controls')
        summary = {'records':len(routes),'rows':len(rows[0]['rows']),'terms':len(terms),'sources':len(rows[0]['jets']),
                   **{k:counts[k] for k in ('newUnionImages','acceptedAnalyticUses','additionalNewImageUses')},
                   'radicalUses':sum(len(v['radicals']) for v in domain_routes.values()),
                   'denominatorUses':sum(len(v['denominators']) for v in domain_routes.values()),'controls':controls}
        cases[label] = summary; f.save(base/'analytic-input-case-inventory.json',cases)
        del packet, chars, rows, terms, routes, domain_routes
        gc.collect()
    f.atomic_pickle(base/'uncomputed-analytic-image-inputs.pickle',pending)
    f.atomic_pickle(base/'uncomputed-radical-inputs.pickle',pending_roots)
    f.atomic_pickle(base/'uncomputed-denominator-inputs.pickle',pending_denominators)
    f.require(tuple(sum(v[k] for v in cases.values()) for k in ('records','rows','terms','sources')) == (1467,300,647,120),
              'all own-case source/row/term/column identities')
    return {'cases':cases,'acceptedBaselineAnalyticRecords':len(analytic),'newUnionImages':len(pending),
            'uncomputedRadicals':len(pending_roots),'uncomputedDenominators':len(pending_denominators),
            'newAnalyticLifts':0,'newFrequencyDerivatives':0,'newNumericalWork':0,
            'scalarDomainScope':roots['scope'],'endFamilyContinuationPending':True}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--run-directory',type=Path,required=True); args = ap.parse_args()
    base = args.run_directory.resolve(); base.relative_to(f.STORE); base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3)); signal.alarm(900); started = time.monotonic()
    protect_references(base)
    manifest, labels = load(base); join = native_join(); f.save(base/'native-analytic-input-joins.json',join)
    prohibit(); result = prepare(base,manifest,labels)
    for name,value in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name) == f.digest(base/'source'/name) == value, 'current/frozen posthash')
    for name,value in manifest['inputPackets'].items(): f.require(f.digest(Path(name)) == value, 'original input posthash')
    for name,item in manifest['referencedInputs'].items():
        path = base/name
        f.require(path.is_symlink() and str(path.readlink()) == item['original'] and str(path.resolve()) == item['resolvedOriginal']
                  and path.stat().st_size == item['bytes'] and f.digest(path) == item['sha256'], 'complete reference identity')
    artifacts = {str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*')
                 if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    checks = {**manifest,**result,'status':'COMPLETED_CASE_FREQUENCY_ANALYTIC_INPUTS','nativeJoins':join,'artifacts':artifacts,
              'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json',checks); signal.alarm(0); print(json.dumps(checks,indent=2))


if __name__ == '__main__': main()
