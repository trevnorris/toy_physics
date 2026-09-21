#!/usr/bin/env python3
"""Saved-source material bindings with exact scalar and baseline reuse."""
import argparse
import ast
import copy
import inspect
import hashlib
import json
from pathlib import Path
import resource
import shutil
import signal
import time
import textwrap
from types import SimpleNamespace

import numpy as np
import S11c_d_coordinate_response as c
import S11c_d_remaining_case_coordinate_sources as h
import S11c_d_remaining_case_profile_response_recover as profile_cache

f, engine, sp, native, modes = h.f, h.engine, h.sp, h.native, h.modes
matrices = h.source.matrices
BASELINE = native.BASELINE
CP = f.M/'S11c_d_remaining_case_coordinate_sources_checkpoint.json'
MCP = f.M/'S11c_d_coordinate_response_checkpoint.json'
BCP = f.M/'S11c_d_remaining_case_bindings_checkpoint.json'
FCP = f.M/'S11c_d_remaining_case_coordinate_bindings_focused.json'
PLAN = f.M/'S11c_d_remaining_case_coordinate_bindings_plan.md'
SCOPE = ('Actual material scalar/source/grade bindings and numerical row reuse routes. '
         'No quadrature, mode/current construction, matrix assembly or response solve.')


def native_binder():
    """Generalize only source count and route numerical calls through exact reuse."""
    original = ast.parse(inspect.getsource(c.bind)).body[0]
    changes = [
        ('range(35)', "sorted(data['old']['source-binding.pickle']['jets'])"),
        ("chart.bind_image(rec['coordinateImage'])", "cache.get('bound-image', rec['coordinateImage'], rec['unit'], ('record', rec['address']))"),
        ("[chart.bind_image(v) for v in rec['materialSourceCoefficients']]",
         "[cache.get('bound-image', v, cache.jet_unit(rec, order), ('source-jet', si, order)) for order, v in enumerate(rec['materialSourceCoefficients'])]"),
        ("chart.image(old['symbolicFrequency'])", "cache.get('image', old['symbolicFrequency'], (-1, 0, 0), ('frequency', ti, si))"),
        ("chart.image(v)", "cache.get('image', v, cache.grade_unit(original['UNIT'], g), ('grade', address, g))"),
        ("engine.memo_xreplace(data['adapter'].bind(rec['original']),chart.cuts)",
         "cache.get('original', rec['original'], rec['unit'], ('original', address))"),
    ]
    def replace(tree, a, b):
        class Edit(ast.NodeTransformer):
            count = 0
            def visit(self, node):
                if ast.dump(node) == ast.dump(a):
                    self.count += 1
                    return copy.deepcopy(b)
                return super().visit(node)
        edit = Edit(); result = edit.visit(tree)
        return result, edit.count
    node = copy.deepcopy(original); counts = []
    for a, b in changes:
        node, count = replace(node, ast.parse(a, mode='eval').body, ast.parse(b, mode='eval').body)
        counts.append(count)
    f.require(counts == [1, 2, 1, 1, 1, 1], 'actual complete native binding call census')
    reverse = copy.deepcopy(node)
    for a, b in reversed(changes):
        reverse, _ = replace(reverse, ast.parse(b, mode='eval').body, ast.parse(a, mode='eval').body)
    f.require(ast.dump(reverse) == ast.dump(original), 'whole native material binder reverse AST')
    node.args.args.append(ast.arg(arg='cache'))
    namespace = dict(vars(c))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[node], type_ignores=[])), __file__, 'exec'), namespace)
    return namespace['bind'], {'wholeNativeBinderReverseAst': True, 'edits': counts,
        'nativeBinderAstSha256': h.source.body(c.bind),
        'nativeChartAstSha256': h.source.body(c.Chart),
        'nativeNumericalBindAstSha256': hashlib.sha256(ast.dump(ast.parse(textwrap.dedent(inspect.getsource(engine.NumericalReducedAction.bind)))).encode()).hexdigest(),
        'cacheComparisonAstSha256': h.source.body(native.same),
        'profileCacheJoinAstSha256': h.source.body(profile_cache.cache_evidence)}


def load(base, resume):
    sr, sc, source_cp = h.source.provenance.accepted(CP, 'ACCEPTED_CASE_COORDINATE_SOURCES')
    mr, mc, material_cp = h.source.provenance.accepted(MCP, 'PUBLISHED_ANNEX_VERIFIED')
    pins = dict(sc['sourceFiles'])
    for n, sha in mc['sourceFiles'].items():
        f.require(n not in pins or pins[n] == sha, ('exact common material source', n)); pins[n] = sha
    for p in (Path(__file__).resolve(), PLAN, CP, MCP, BCP, Path(c.__file__).resolve(),
              Path(profile_cache.__file__).resolve()):
        pins[str(p.relative_to(f.ROOT))] = f.digest(p)
    manifest = {'runDirectory': str(base), 'sourceFiles': pins, 'inputPackets': {},
                'copiedInputs': {}, 'scope': SCOPE, 'input': sc['input'],
                'settings': sc['settings'], 'sourceAcceptance': source_cp['checksSha256'],
                'baselineMaterialAcceptance': material_cp['checksSha256']}
    for root, checks in ((sr, sc), (mr, mc)):
        manifest['inputPackets'][str(root/'checks.json')] = f.digest(root/'checks.json')
        for n, sha in checks['inputPackets'].items():
            f.require(n not in manifest['inputPackets'] or manifest['inputPackets'][n] == sha,
                      'exact common original material input'); manifest['inputPackets'][n] = sha
    if resume:
        old = json.loads((resume/'checks.json').read_text()); accepted = json.loads(FCP.read_text())
        f.require(accepted['status'] == 'ACCEPTED_CASE_MATERIAL_BINDING_INPUTS'
                  and accepted['checksSha256'] == f.digest(resume/'checks.json'), 'accepted material binding focus')
        f.require(old['sourceFiles'] == pins, 'unchanged complete focus implementation')
        for n, item in old['artifacts'].items():
            modes.retain(resume/n, base/n, manifest, item['sha256'])
        modes.retain(resume/'checks.json', base/'accepted-focus-checks.json', manifest)
        modes.retain(resume/'inputs.json', base/'accepted-focus-inputs.json', manifest)
        manifest['completedFocusReuse'] = {'directory': str(resume), 'artifacts': len(old['artifacts']),
                                         'checksSha256': f.digest(resume/'checks.json')}
    else:
        for n, item in sc['artifacts'].items(): modes.retain(sr/n, base/'sources'/n, manifest, item['sha256'])
        modes.retain(sr/'checks.json', base/'accepted-source-checks.json', manifest)
        modes.retain(sr/'inputs.json', base/'accepted-source-inputs.json', manifest)
        for n, item in mc['artifacts'].items(): modes.retain(mr/n, base/'baseline-material'/n, manifest, item['sha256'])
        modes.retain(mr/'inputs.json', base/'baseline-material/inputs.json', manifest)
        modes.retain(mr/'checks.json', base/'baseline-material/checks.json', manifest)
        bc = json.loads(BCP.read_text()); br = Path(bc['runDirectory'])
        f.require(bc['status'] == 'ACCEPTED_BINDINGS_AND_GRADES' and f.digest(br/'checks.json') == bc['checksSha256'], 'actual original case trial-system producer')
        modes.retain(br/'accepted-finite-system.pickle', base/'accepted-finite-system.pickle', manifest,
                     bc['artifacts']['accepted-finite-system.pickle']['sha256'])
    for n, sha in pins.items():
        dst = base/'source'/n; dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(f.ROOT/n, dst); f.require(f.digest(dst) == sha, 'current/frozen binding source')
    f.save(base/'inputs.json', manifest)
    return manifest, tuple(sc['cases'])


def input_state(adapter):
    return {k: v for k, v in vars(adapter.input).items() if k != 'r'}


def restore(base, label, manifest, saved=False):
    source = base/'sources'; case = f.unpickle(source/'accepted-bindings'/label/'case-binding.pickle')
    coordinate = f.unpickle(source/'coordinate-cases'/label/'coordinate-source.pickle')
    packets = {kind: f.unpickle(source/'accepted-cases'/label/(kind+'.pickle'))
               for kind in ('reduced-action', 'actions', 'assembly', 'factorization')}
    r, dims, pencil = native.factors.context(packets['reduced-action'], packets['actions'], packets['assembly'])
    dims.__dict__.update(coordinate['dimensionState'])
    if saved:
        adapter = engine.NumericalReducedAction.__new__(engine.NumericalReducedAction)
        adapter.r, adapter.pencil, adapter.assembly = r, pencil, packets['assembly']['result']
        adapter.input = engine.ChannelInput.__new__(engine.ChannelInput)
        adapter.input.__dict__.update(f.unpickle(base/'contexts'/label/'input-state.pickle')); adapter.input.r = r
    else:
        adapter = engine.NumericalReducedAction(pencil, packets['assembly']['result'], manifest['input'])
        path = base/'contexts'/label; path.mkdir(parents=True)
        f.atomic_pickle(path/'input-state.pickle', input_state(adapter))
    data = {'r': r, 'adapter': adapter, 'coordinate': coordinate, 'packet': case['grades'],
            'old': {'domain-binding.pickle': {'bound': case['binding']['bound']},
                    'source-binding.pickle': {'jets': case['binding']['jets']}},
            'system': {'settings': case['binding']['settings']}}
    return data, case, packets


def frame(data, case):
    r = data['r']; b = case['binding']
    return {'input': input_state(data['adapter']), 'geometry': data['coordinate']['chart'],
            'fieldJets': data['coordinate']['fieldJets'], 'fieldUnits': b['fieldUnits'],
            'equationUnits': b['equationUnits'], 'settings': b['settings'],
            'cuts': b['bound']['cutoffBindings'], 'abel': b['bound']['abel'],
            'profileUnits': b['bound']['profileUnits'], 'samplingPairs': b['bound']['pairs'],
            'coordinates': (r.z, r.zp, r.xi, r.ell, r.omega, r.tangents),
            'profiles': r.profiles, 'generators': case['grades']['generators']}


def chart_state(chart):
    names = ('g', 'jets', 'd', 'A', 'B', 'kappa', 'rename', 'cuts', 'forward', 'inverse', 'composition')
    return {'values': {k: getattr(chart, k) for k in names}, 'materialInput': input_state(chart.material)}


def restore_chart(data, state):
    chart = c.Chart.__new__(c.Chart)
    chart.__dict__.update(state['values']); chart.r = data['r']; chart.saved = data['coordinate']; chart.adapter = data['adapter']
    chart.material = engine.NumericalReducedAction.__new__(engine.NumericalReducedAction)
    chart.material.r = chart.r; chart.material.pencil = SimpleNamespace(r=chart.r); chart.material.assembly = {}
    chart.material.input = engine.ChannelInput.__new__(engine.ChannelInput)
    chart.material.input.__dict__.update(state['materialInput']); chart.material.input.r = chart.r
    return chart


class ScalarCache:
    def __init__(self, folder, entries=()):
        self.folder = folder; self.entries = list(entries); self.buckets = {}; self.routes = []; self.new = 0
        for i, v in enumerate(self.entries): self.buckets.setdefault((v['kind'], hash(v['expression']), tuple(v['unit'])), []).append(i)

    def select(self, kind, expression, unit):
        return [i for i in self.buckets.get((kind, hash(expression), tuple(unit)), ())
                if native.same(self.entries[i]['expression'], expression)]

    def add(self, kind, expression, unit, value, owner):
        found = self.select(kind, expression, unit)
        if found:
            f.require(modes.same(self.entries[found[0]]['value'], value), 'same literal binding input has same saved value')
            return found[0]
        i = len(self.entries); record = {'kind': kind, 'expression': expression,
             'unit': tuple(unit), 'value': value, 'owner': owner}
        self.entries.append(record); self.buckets.setdefault((kind, hash(expression), tuple(unit)), []).append(i)
        return i

    def configure(self, data, case, chart, label, allow_new):
        self.data, self.case, self.chart, self.label, self.allow_new = data, case, chart, label, allow_new
        dims = engine.PHYSICAL_METADATA.dimensions
        self.generator_units = [tuple(dims.measure(g)) for g in case['grades']['generators']]

    def grade_unit(self, unit, grade):
        return tuple(v-sum(n*u[j] for n, u in zip(grade, self.generator_units)) for j, v in enumerate(unit))

    def jet_unit(self, record, order):
        column = record['sourceJets']['column']; field = self.case['binding']['fieldUnits'][column]
        return tuple(v-w+order*length for v, w, length in zip(record['unit'], field, (1,0,0)))

    def get(self, kind, expression, unit, address):
        matches = self.select(kind, expression, unit)
        if matches: index = matches[0]
        else:
            f.require(self.allow_new and kind != 'original', ('missing saved scalar binding', self.label, kind, address))
            index = len(self.entries); target = self.folder/f'{index:05d}'; target.mkdir(parents=True)
            raw = {'case': self.label, 'kind': kind, 'expression': expression, 'unit': tuple(unit), 'address': address}
            f.atomic_pickle(target/'input.pickle', raw)
            image = c.source.coordinate_change(expression, self.chart.g['definitions'], self.chart.g['forward']) if kind == 'image' else expression
            f.atomic_pickle(target/'material-image.pickle', image)
            value = self.chart.bind_image(image)
            f.atomic_pickle(target/'value.pickle', value)
            index = self.add(kind, expression, unit, value, {'case': self.label, 'address': address, 'path': str(target)})
            self.new += 1
            f.atomic_pickle(target/'record.pickle', self.entries[index])
            f.save(self.folder/'inventory.json', {'newScalars': self.new, 'entries': len(self.entries)})
        route = {'case': self.label, 'address': address, 'kind': kind,
                 'expression': expression, 'unit': tuple(unit), 'entry': index}
        self.routes.append(route)
        return self.entries[index]['value']


def seed_original(cache, data, case, packets, label):
    b = case['binding']; cells = {(v['row'], v['column']): v for v in b['bound']['cells'] if v['test'] == 0}
    count = 0
    for key, item in case['grades']['records'].items():
        address = item['address']; record = item['record']; kind = address[0]
        h.source.require_record(item, data['coordinate']['records'][key]['address'],
                                data['coordinate']['records'][key]['original'], data['coordinate']['records'][key]['unit'])
        if kind == 'local':
            _, order, i, j = address; value = b['local'][order][i,j]
        elif kind == 'cell':
            _, i, j, term = address; native_cell = next(v for v in packets['assembly']['result']['ROWS'] if (v['ROW'],v['COLUMN']) == (i,j))
            integral, coefficient = native_cell['NONLOCAL'][term]
            index, value = cells[i,j]['terms'][term]
            f.require(native.same((record['ORIGINAL'], integral), (coefficient, b['bound']['rows'][index]['original'])), 'full original cell coefficient/integral binding join')
        else: continue
        cache.add('original', record['ORIGINAL'], record['UNIT'], value, {'case': label, 'key': key, 'address': address})
        count += 1
    return count


def seed_material(cache, data, case, accepted):
    by_address = {tuple(v['address']): v for v in data['coordinate']['records'].values()}
    for address, value in accepted['bindings'].items():
        rec = by_address[address]
        cache.add('bound-image', rec['coordinateImage'], rec['unit'], value['material'], {'case': BASELINE, 'address': address})
    for si, jet in accepted['jets'].items():
        rec = by_address['source', si]
        for order, (expr, value) in enumerate(zip(rec['materialSourceCoefficients'], jet['coefficients'])):
            cache.add('bound-image', expr, cache.jet_unit(rec, order), value, {'case': BASELINE, 'source': si, 'order': order})
        f.require(len(rec['materialSourceCoefficients']) == len(jet['coefficients']), 'complete saved source derivative coefficient count')
    for (ti, si), value in accepted['sources'].items():
        cache.add('image', value['symbolicFrequency'], (-1,0,0), value['frequency']/cache.chart.d, {'case': BASELINE, 'source': si, 'test': ti})
    for item in case['grades']['records'].values():
        address = item['address']; record = item['record']
        if address[0] not in ('local', 'cell'): continue
        for grade, expr in record['COEFFICIENTS'].items():
            cache.add('image', expr, cache.grade_unit(record['UNIT'], grade), accepted['coefficients'][address][grade],
                      {'case': BASELINE, 'address': address, 'grade': grade})


def signature(row, material, case, shared_frame):
    original = case['binding']['bound']; factors = []
    for value in row['factors']:
        si = value['sourceIndex']; src = material['sources'][0,si]; jet = material['jets'][si]
        factors.append((value['symbolicCoefficient'], value['coefficient'], value['unit'],
            src['originalSourceIntegral'], src['symbolicAmplitude'], src['symbolicFrequency'],
            src['frequency'], src['boundCharacter'], jet['column'], jet['probe'],
            tuple(jet['coefficients']), jet['amplitudeUnit'], jet['integralUnit']))
    return (row['original'], row['symbolicLimits'], row['limits'], row['sourceLimit'], row['unit'],
            tuple(factors), shared_frame, original['profileUnits'], original['abel'], material['settings'])


def focused(base, manifest, labels, bind, joins):
    cache = ScalarCache(base/'new-scalars'); contexts = {}; comparisons = {}; original_counts = {}
    for label in labels:
        data, case, packets = restore(base, label, manifest)
        contexts[label] = (data, case, packets)
        actual = frame(data, case)
        if label == BASELINE: common = actual
        else:
            f.require(native.same(actual, common), 'complete common material coefficient frame')
            proof = profile_cache.cache_evidence(case['binding']['bound'], contexts[BASELINE][1]['binding']['bound'])
            f.atomic_pickle(base/'contexts'/label/'profile-cache-pair.pickle', proof)
        f.atomic_pickle(base/'contexts'/label/'frame.pickle', actual)
        original_counts[label] = seed_original(cache, data, case, packets, label)
        comparisons[label] = {'records': len(case['grades']['records']), 'rows': len(case['binding']['bound']['rows']),
                              'sources': len(case['binding']['jets']), 'terms': len(case['grades']['termJoins'])}
    data, case, _ = contexts[BASELINE]
    chart = c.Chart(data); state = chart_state(chart); f.atomic_pickle(base/'chart-state.pickle', state)
    cache.configure(data, case, chart, BASELINE, False)
    accepted = f.unpickle(base/'baseline-material/material-binding.pickle')
    f.require(native.same(accepted['geometry'], chart.g) and modes.same(accepted['chartComposition'], chart.composition), 'actual accepted full chart composition')
    seed_material(cache, data, case, accepted)
    target = base/'baseline-replay'; target.mkdir()
    replay = bind(target, data, chart, cache)
    f.require(cache.new == 0 and modes.same(replay, accepted), 'full baseline native binding replay from saved scalars, no rebinding')
    rows = {}; layouts = []; system = f.unpickle(base/'accepted-finite-system.pickle')
    f.require(native.same(system['settings'], accepted['settings']) and len(system['nodes']) == 129,
              'actual accepted material/Eulerian trial nodes and full settings')
    for count in (1,2,3):
        group = f.unpickle(base/'baseline-material'/f'layout-{count}.pickle')
        f.require(native.same(group['setting'], accepted['settings']) and np.isfinite(group['matrices']).all(), 'actual baseline material layout settings/arrays')
        f.require(abs(group['massResidual']) < 1e-9*(1+abs(group['mass'])) and np.max(abs(group['actionResidual']))/(1+np.max(abs(group['direct']))) < 1e-10, 'saved complete material measure/direct-action evidence')
        f.require(group['matrices'].shape[1:] == (len(system['nodes']), len(system['nodes'])), 'actual full trial-matrix dimensions')
        rows.update(zip(group['rowIndices'], group['matrices'])); layouts.append({'layout': count, 'rows': group['rowIndices'], 'nodes': group['nodes']})
    f.require(set(rows) == set(range(len(accepted['rows']))), 'complete actual material row matrices')
    f.atomic_pickle(base/'baseline-material-row-matrices.pickle', rows)
    controls = []
    item = next(v for v in cache.entries if v['expression'] != 0 and v['value'] != 0)
    for changed in (dict(item, expression=2*item['expression']), dict(item, unit=(item['unit'][0]+1,*item['unit'][1:]))):
        rejected = not native.same((item['expression'],item['unit']), (changed['expression'],changed['unit']))
        f.require(rejected, 'actual changed scalar/unit rejects reuse'); controls.append({'original': item, 'changed': changed, 'rejected': rejected})
    changed = copy.deepcopy(common); changed['input']['parameters']['omega'] *= 2
    f.require(not native.same(changed, common), 'actual changed physical frequency frame rejects')
    controls.append({'originalFrame': common, 'changedFrame': changed, 'rejected': True})
    actual = next(row for row in accepted['rows'] if any(v['coefficient'] != 0 for v in row['factors']))
    for kind in ('limit', 'coefficient'):
        changed = copy.deepcopy(actual)
        if kind == 'limit':
            v, a, b = changed['limits'][0]; changed['limits'] = ((v,a,b+1),*changed['limits'][1:])
        else:
            factor = next(v for v in changed['factors'] if v['coefficient'] != 0); factor['coefficient'] *= 2
        f.require(not native.same(signature(actual,accepted,case,common), signature(changed,accepted,case,common)), 'actual changed material row reuse control')
        controls.append({'kind': kind, 'original': actual, 'changed': changed, 'rejected': True})
    item = next(iter(case['grades']['records'].values())); wrong = (*item['address'][:-1],item['address'][-1]+1)
    rejected = False
    try: h.source.require_record(item,wrong,item['record']['ORIGINAL'],item['record']['UNIT'])
    except (ValueError, AssertionError, RuntimeError): rejected = True
    f.require(rejected,'actual wrong physical record address rejects')
    controls.append({'kind':'address','original':item['address'],'changed':wrong,'rejected':rejected})
    f.atomic_pickle(base/'reuse-controls.pickle', controls)
    f.atomic_pickle(base/'binding-cache.pickle', cache.entries)
    f.atomic_pickle(base/'baseline-binding-routes.pickle', cache.routes)
    f.save(base/'preflight.json', {'cases': comparisons, 'originalScalarUses': original_counts,
        'savedMaterialRows': len(rows), 'layouts': layouts, 'baselineBindingCallsReused': len(cache.routes),
        'cachedScalarEntries': len(cache.entries), 'newScalarBindings': 0, 'controls': len(controls),
        'nativeJoins': joins, 'commonChartContextOnly': True})
    return comparisons


def construct(base, manifest, labels, bind):
    state = f.unpickle(base/'chart-state.pickle'); cache = ScalarCache(base/'new-scalars', f.unpickle(base/'binding-cache.pickle'))
    common = f.unpickle(base/'contexts'/BASELINE/'frame.pickle')
    baseline = f.unpickle(base/'baseline-material/material-binding.pickle')
    basecase = f.unpickle(base/'sources/accepted-bindings'/BASELINE/'case-binding.pickle')
    atlas = {}; results = {}; inventory = {BASELINE: {'rows': len(baseline['rows']), 'newRows': [], 'wholeBindingReused': True}}
    for row in baseline['rows']:
        atlas.setdefault(hash(row['original']), []).append((BASELINE,row['index'],signature(row,baseline,basecase,common)))
    for label in labels:
        if label == BASELINE: continue
        data, case, _ = restore(base, label, manifest, saved=True)
        f.require(native.same(frame(data,case),common), 'complete saved/current material binding frame')
        chart = restore_chart(data,state); cache.configure(data,case,chart,label,True)
        folder = base/'material-cases'/label; folder.mkdir(parents=True)
        start = len(cache.routes); before = cache.new; material = bind(folder,data,chart,cache)
        f.atomic_pickle(folder/'binding-routes.pickle',cache.routes[start:])
        reused = []; fresh = []; row_pairs = []
        for row in material['rows']:
            sig = signature(row,material,case,common); bucket = atlas.setdefault(hash(row['original']),[])
            matches = [(a,i) for a,i,v in bucket if native.same(sig,v)]
            if matches: reused.append({'row':row['index'],'fromCase':matches[0][0],'fromRow':matches[0][1]})
            else: fresh.append(row['index']);bucket.append((label,row['index'],sig))
            row_pairs.append({'row':row['index'],'signature':sig,'owner':matches[0] if matches else (label,row['index'])})
        f.atomic_pickle(folder/'row-reuse-pairs.pickle',row_pairs)
        actual=material['rows'][0];changed=copy.deepcopy(actual);v,a,b=changed['limits'][0];changed['limits']=((v,a,b+1),*changed['limits'][1:])
        f.require(not native.same(signature(actual,material,case,common),signature(changed,material,case,common)), 'actual changed ordered numerical limit rejects')
        f.atomic_pickle(folder/'row-mutation.pickle',{'original':actual,'changed':changed})
        result={'material':material,'originalCase':case,'coordinate':data['coordinate'],'newRows':fresh,'reusedRows':reused,
                'frame':common,'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets']}
        f.atomic_pickle(folder/'case-material-binding.pickle',result);results[label]=str(folder/'case-material-binding.pickle')
        inventory[label]={'rows':len(material['rows']),'sources':len(material['jets']), 'terms':len(case['grades']['termJoins']),
                          'newRows':fresh,'reusedRows':reused,'newScalarBindings':cache.new-before}
        f.atomic_pickle(folder/'completed-binding-cache.pickle',cache.entries)
        f.save(base/'case-inventory.json',inventory)
    f.atomic_pickle(base/'final-binding-cache.pickle',cache.entries)
    f.atomic_pickle(base/'remaining-case-material-bindings.pickle',{'cases':results,'inventory':inventory,'baseline':str(base/'baseline-material/material-binding.pickle'),'scope':SCOPE})
    return inventory


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);ap.add_argument('--focused',action='store_true');ap.add_argument('--resume-focused',type=Path);args=ap.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);start=time.monotonic()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    f.require(args.focused != bool(args.resume_focused), 'focus or accepted-focus production required')
    bind,joins=native_binder();manifest,labels=load(base,args.resume_focused)
    f.save(base/'native-binding-joins.json',joins) if not (base/'native-binding-joins.json').exists() else f.require(json.loads((base/'native-binding-joins.json').read_text())==joins,'copied native binding joins')
    cases=focused(base,manifest,labels,bind,joins) if args.focused else construct(base,manifest,labels,bind)
    for n,v in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(base/'source'/n)==v,'current/frozen source pre/post')
    for n,v in manifest['inputPackets'].items():f.require(f.digest(Path(n))==v,'original input pre/post')
    for n,v in manifest['copiedInputs'].items():f.require(f.digest(base/n)==v,'copied input pre/post')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    checks={**manifest,'status':'COMPLETED_MATERIAL_BINDING_FOCUS' if args.focused else 'COMPLETED_CASE_MATERIAL_BINDINGS','cases':cases,'artifacts':artifacts,'nativeJoins':joins,'newQuadratureNodes':0,'newModeCurrentConstructions':0,'newScatteringSolves':0,'wallSeconds':time.monotonic()-start}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
