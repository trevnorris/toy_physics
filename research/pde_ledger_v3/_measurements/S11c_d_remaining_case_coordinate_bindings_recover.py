#!/usr/bin/env python3
"""Finish saved material-binding focus; keep the original native binder untouched."""
import ast
import copy
import gc
import inspect
import json
from pathlib import Path
import resource
import signal
import sys
import time

import S11c_d_remaining_case_coordinate_bindings as h

f, sp, native, modes = h.f, h.sp, h.native, h.modes
ORIGIN = f.STORE/'s11c-remaining-case-coordinate-20260921/bindings/focused/complete'
PLAN = f.M/'S11c_d_remaining_case_coordinate_bindings_recovery_plan.md'
REPAIR = f.M/'S11c_d_remaining_case_coordinate_bindings_schema_repair.json'
DIAGNOSTICS = frozenset(('frequencies', 'boundSource', 'assignmentResidual', 'assignments'))
ORIGINAL_LOAD = h.load
ORIGINAL_FOCUS = h.focused


def baseline_source_join(replay, accepted, historical):
    """Compare complete packets with an explicit view of historical diagnostics."""
    f.require(replay.keys() == accepted.keys(), 'same complete material packet schema')
    f.require(modes.same({k:v for k,v in replay.items() if k != 'sources'},
                         {k:v for k,v in accepted.items() if k != 'sources'}),
              'all nonsource material payloads exactly identical')
    current, prior = replay['sources'], accepted['sources']
    f.require(current.keys() == prior.keys() == historical.keys(), 'all actual source/test addresses')
    proofs = []
    for address in sorted(current):
        left, right, original = current[address], prior[address], historical[address]
        f.require(set(right)-set(left) == DIAGNOSTICS and not set(left)-set(right),
                  'only four recorded historical source diagnostic keys')
        f.require(right.keys() == original.keys(), 'complete original diagnostic schema')
        f.require(modes.same({k:v for k,v in right.items() if k != 'frequency'},
                             {k:v for k,v in original.items() if k != 'frequency'}),
                  'material diagnostics and inherited fields exactly join original producer')
        # This conversion is confined to this recorded numeric width, not expressions.
        a, b = left['profileWidth'], right['profileWidth']
        f.require(isinstance(a, sp.Integer) and type(b) is float and a.is_positive is True,
                  'actual saved Integer/float positive width types')
        f.require(h.np.isfinite(b) and sp.Rational(b) == a,
                  'exact width rational identity, no tolerance or normalization')
        f.require(modes.same({k:v for k,v in left.items() if k != 'profileWidth'},
                             {k:v for k,v in right.items() if k not in DIAGNOSTICS | {'profileWidth'}}),
                  'every shared physical source field exactly identical')
        f.require(left['sourceIndex'] == address[1] and left['test'] == address[0],
                  'actual source/test physical address')
        proofs.append({'address': address, 'current': left, 'acceptedMaterial': right,
                       'originalProducer': original, 'widthResidual': a-sp.Rational(b),
                       'retainedDiagnosticKeys': tuple(sorted(DIAGNOSTICS))})
    return proofs


def original_domain(root):
    manifest = json.loads((root/'baseline-material/inputs.json').read_text())
    paths = [Path(n) for n in manifest['inputPackets'] if Path(n).name == 'domain-binding.pickle']
    f.require(len(paths) == 1, 'one exact historical diagnostic producer')
    path = paths[0]
    f.require(f.digest(path) == manifest['inputPackets'][str(path)], 'original source diagnostic packet hash')
    return path, f.unpickle(path)['bound']['sources']


def regression(base):
    base.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    replay = f.unpickle(ORIGIN/'baseline-replay/material-binding.pickle')
    accepted = f.unpickle(ORIGIN/'baseline-material/material-binding.pickle')
    historical_path, historical = original_domain(ORIGIN)
    proofs = baseline_source_join(replay, accepted, historical)
    f.atomic_pickle(base/'source-view-pairs.pickle', proofs)
    controls = []
    address = next(a for a, v in replay['sources'].items() if v['boundAmplitude'] != 0 and v['frequency'] != 0)
    source = replay['sources'][address]
    changes = [('amplitude', 'boundAmplitude', 2*source['boundAmplitude']),
               ('frequency', 'frequency', 2*source['frequency']),
               ('width-value', 'profileWidth', source['profileWidth']+1),
               ('width-type', 'profileWidth', float(source['profileWidth'])),
               ('width-bool', 'profileWidth', True),
               ('source-address', 'sourceIndex', source['sourceIndex']+1),
               ('unit', 'amplitudeUnit', (source['amplitudeUnit'][0]+1,*source['amplitudeUnit'][1:]))]
    for label, key, value in changes:
        changed = dict(replay, sources=dict(replay['sources']))
        changed['sources'][address] = dict(source, **{key:value})
        rejected = False
        try: baseline_source_join(changed, accepted, historical)
        except (ValueError, AssertionError, RuntimeError): rejected = True
        f.require(rejected, ('actual changed source view rejected', label))
        controls.append({'kind':label,'address':address,'original':source,'changed':changed['sources'][address],'rejected':rejected})
    for key in sorted(DIAGNOSTICS):
        changed = dict(accepted, sources=dict(accepted['sources']))
        changed['sources'][address] = dict(accepted['sources'][address])
        del changed['sources'][address][key]
        rejected = False
        try: baseline_source_join(replay, changed, historical)
        except (ValueError, AssertionError, RuntimeError): rejected = True
        f.require(rejected, ('missing historical diagnostic rejected',key))
        controls.append({'kind':'missing-diagnostic','address':address,'key':key,'rejected':rejected})
    changed = dict(accepted, sources=dict(accepted['sources']))
    changed['sources'][address] = dict(accepted['sources'][address])
    original = changed['sources'][address]['boundSource']
    changed['sources'][address]['boundSource'] = 2*original
    f.require(not modes.same(original,2*original), 'actual diagnostic coefficient responds')
    rejected = False
    try: baseline_source_join(replay, changed, historical)
    except (ValueError, AssertionError, RuntimeError): rejected = True
    f.require(rejected, 'changed historical diagnostic rejects')
    controls.append({'kind':'diagnostic-coefficient','address':address,'original':original,'changed':2*original,'rejected':rejected})
    f.atomic_pickle(base/'mutation-controls.pickle', controls)
    _, wiring = recovered_focus()
    paths = [Path(__file__).resolve(),Path(h.__file__).resolve(),PLAN,
             ORIGIN/'inputs.json',ORIGIN/'baseline-replay/material-binding.pickle',
             ORIGIN/'baseline-material/material-binding.pickle',historical_path]
    checks = {'status':'PASSED_SAVED_MATERIAL_BINDING_SCHEMA_REPAIR','sourcePairs':len(proofs),
              'widthResidualsZero':all(v['widthResidual'] == 0 for v in proofs),
              'retainedHistoricalDiagnostics':len(proofs)*len(DIAGNOSTICS),
              'respondingControls':len(controls),'wiring':wiring,
              'inputs':{str(p):f.digest(p) for p in paths},
              'artifacts':{str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size}
                           for p in base.iterdir() if p.is_file()},
              'newBindings':0,'newScientificConstruction':0,'wallSeconds':time.monotonic()-start}
    f.save(base/'checks.json',checks); print(json.dumps(checks,indent=2))


def load(base, resume):
    origin = resume if resume else ORIGIN
    old = json.loads((origin/'inputs.json').read_text())
    if resume:
        cp = json.loads(h.FCP.read_text()); checks = json.loads((origin/'checks.json').read_text())
        f.require(cp['status'] == 'ACCEPTED_CASE_MATERIAL_BINDING_INPUTS'
                  and cp['checksSha256'] == f.digest(origin/'checks.json'), 'accepted recovered focus')
        old = checks
    # The original helper and each frozen source remain exact. The separate wrapper
    # is added explicitly, never substituted for a historical source fingerprint.
    for n, sha in old['sourceFiles'].items():
        f.require(f.digest(f.ROOT/n) == f.digest(origin/'source'/n) == sha, 'original/current/frozen source join')
    for n, sha in old['inputPackets'].items(): f.require(f.digest(Path(n)) == sha, 'original input hash')
    for n, sha in old['copiedInputs'].items(): f.require(f.digest(origin/n) == sha, 'original completed input copy')
    manifest = {k:copy.deepcopy(old[k]) for k in ('sourceFiles','inputPackets','scope','input','settings','sourceAcceptance','baselineMaterialAcceptance')}
    manifest.update(runDirectory=str(base),copiedInputs={})
    if resume:
        for n, item in checks['artifacts'].items(): modes.retain(origin/n,base/n,manifest,item['sha256'])
        modes.retain(origin/'checks.json',base/'accepted-focus-checks.json',manifest)
        modes.retain(origin/'inputs.json',base/'accepted-focus-inputs.json',manifest)
        manifest['completedFocusReuse']={'directory':str(origin),'artifacts':len(checks['artifacts']),
                                        'checksSha256':f.digest(origin/'checks.json')}
    else:
        # Copy every completed packet, including frames and the baseline replay.
        for p in sorted(origin.rglob('*')):
            if p.is_file() and p != origin/'inputs.json': modes.retain(p,base/p.relative_to(origin),manifest)
        modes.retain(origin/'inputs.json',base/'original-focus-inputs.json',manifest)
        for p in sorted(origin.parent.rglob('*')):
            if p.is_file() and 'complete' not in p.relative_to(origin.parent).parts and 'completion-watcher' not in p.relative_to(origin.parent).parts:
                modes.retain(p,base/'original-focus-logs'/p.relative_to(origin.parent),manifest)
        manifest['completedInputReuse']={'directory':str(origin),'copiedFiles':len(manifest['copiedInputs']),
            'savedBaselineReplaySha256':f.digest(origin/'baseline-replay/material-binding.pickle'),
            'nativeBinderReexecuted':False,'chartReconstructed':False}
    for p in (Path(__file__).resolve(),PLAN,REPAIR):
        manifest['sourceFiles'][str(p.relative_to(f.ROOT))]=f.digest(p)
    repair = json.loads(REPAIR.read_text())
    f.require(repair['status'] == 'ACCEPTED_SAVED_SOURCE_VIEW_REPAIR', 'focused actual schema regression accepted')
    rr = Path(repair['runDirectory']); rc = json.loads((rr/'checks.json').read_text())
    f.require(f.digest(rr/'checks.json') == repair['checksSha256'], 'actual schema regression hash')
    for n,sha in rc['inputs'].items(): f.require(f.digest(Path(n)) == sha, 'exact repair source/operand join')
    manifest['inputPackets'][str(rr/'checks.json')]=repair['checksSha256']
    for n,item in rc['artifacts'].items(): manifest['inputPackets'][str(rr/n)]=item['sha256']
    for n,sha in manifest['sourceFiles'].items():
        dst=base/'source'/n
        if dst.exists(): f.require(f.digest(dst)==sha,'preserved original source snapshot')
        else: modes.retain(f.ROOT/n,dst,manifest,sha)
    f.save(base/'inputs.json',manifest)
    labels=tuple(json.loads((base/'accepted-source-checks.json').read_text())['cases'])
    return manifest,labels


def restore_routes(cache, data, case, replay):
    """Recover cache addresses from saved operands, without calling the binder."""
    by_address={tuple(v['address']):v for v in data['coordinate']['records'].values()}
    proofs=[]
    def lookup(kind,expr,unit,address,expected):
        value=cache.get(kind,expr,unit,address)
        f.require(modes.same(value,expected),('saved scalar route/value identity',address))
        proofs.append({'route':cache.routes[-1],'savedValue':expected})
    for address,value in replay['bindings'].items():
        rec=by_address[address]
        lookup('bound-image',rec['coordinateImage'],rec['unit'],('record',rec['address']),value['material'])
        if address[0] in ('local','cell'):
            lookup('original',rec['original'],rec['unit'],('original',address),value['original'])
    for si,jet in replay['jets'].items():
        rec=by_address['source',si]
        f.require(len(rec['materialSourceCoefficients'])==len(jet['coefficients']),'all saved jet orders')
        for order,(expr,value) in enumerate(zip(rec['materialSourceCoefficients'],jet['coefficients'])):
            lookup('bound-image',expr,cache.jet_unit(rec,order),('source-jet',si,order),value)
    for (ti,si),value in replay['sources'].items():
        original=data['old']['domain-binding.pickle']['bound']['sources'][ti,si]
        lookup('image',original['symbolicFrequency'],(-1,0,0),('frequency',ti,si),value['frequency']/cache.chart.d)
    for item in case['grades']['records'].values():
        address=item['address'];rec=item['record']
        if address[0] not in ('local','cell'):continue
        for grade,expr in rec['COEFFICIENTS'].items():
            lookup('image',expr,cache.grade_unit(rec['UNIT'],grade),('grade',address,grade),replay['coefficients'][address][grade])
    f.require(cache.new==0,'no new scalar binding during route restoration')
    return proofs


def restore_completed_focus(base,manifest,labels,joins):
    def forbidden(*args,**kwargs): raise RuntimeError('scientific construction disabled in saved focus recovery')
    # These operations completed in the original focus. Recovery only reads them.
    h.c.Chart.__init__=forbidden;h.c.Chart.bind_image=forbidden;h.c.Chart.image=forbidden
    h.engine.NumericalReducedAction.bind=forbidden;h.c.source.coordinate_change=forbidden
    cache=h.ScalarCache(base/'new-scalars');comparisons={};original_counts={}
    common=f.unpickle(base/'contexts'/h.BASELINE/'frame.pickle')
    # Release each large restored context before the next one. Cache entries own
    # only their consumed expressions and bound values, not all four case graphs.
    for label in labels:
        data,case,packets=h.restore(base,label,manifest,saved=True)
        f.require(native.same(h.frame(data,case),common),'saved complete common material frame')
        f.require(native.same(f.unpickle(base/'contexts'/label/'frame.pickle'),common),'original completed frame join')
        original_counts[label]=h.seed_original(cache,data,case,packets,label)
        comparisons[label]={'records':len(case['grades']['records']),'rows':len(case['binding']['bound']['rows']),
                            'sources':len(case['binding']['jets']),'terms':len(case['grades']['termJoins'])}
        del data,case,packets;gc.collect()
    data,case,_=h.restore(base,h.BASELINE,manifest,saved=True)
    chart=h.restore_chart(data,f.unpickle(base/'chart-state.pickle'))
    cache.configure(data,case,chart,h.BASELINE,False)
    accepted=f.unpickle(base/'baseline-material/material-binding.pickle')
    replay=f.unpickle(base/'baseline-replay/material-binding.pickle')
    f.require(native.same(accepted['geometry'],chart.g) and modes.same(accepted['chartComposition'],chart.composition),
              'original complete chart composition')
    path,historical=original_domain(base)
    pairs=baseline_source_join(replay,accepted,historical)
    f.atomic_pickle(base/'baseline-source-view-pairs.pickle',pairs)
    manifest['inputPackets'][str(path)]=f.digest(path)
    h.seed_material(cache,data,case,accepted)
    routes=restore_routes(cache,data,case,replay)
    f.atomic_pickle(base/'baseline-route-value-joins.pickle',routes)
    f.save(base/'completed-binding-reuse.json',{'sourcePairs':len(pairs),'historicalDiagnosticFields':len(pairs)*4,
        'scalarRoutes':len(routes),'routeOrder':'reconstructed from saved physical addresses; not claimed original call order',
        'newScalarBindings':cache.new,'binderReexecuted':False,'chartReconstructed':False,
        'savedReplaySha256':f.digest(base/'baseline-replay/material-binding.pickle')})
    f.save(base/'inputs.json',manifest)
    return cache,common,comparisons,original_counts,case,accepted


def recovered_focus():
    original=ast.parse(inspect.getsource(ORIGINAL_FOCUS)).body[0]
    index=next(i for i,n in enumerate(original.body) if isinstance(n,ast.Assign)
               and any(isinstance(t,ast.Name) and t.id=='rows' for t in n.targets))
    replacement=ast.parse('cache, common, comparisons, original_counts, case, accepted = restore_completed_focus(base, manifest, labels, joins)').body[0]
    node=copy.deepcopy(original);node.body=[replacement]+node.body[index:]
    reverse=copy.deepcopy(node);reverse.body=copy.deepcopy(original.body[:index])+reverse.body[1:]
    f.require(ast.dump(reverse)==ast.dump(original),'whole focus reverse AST with saved-prefix replacement')
    f.require(ast.dump(ast.Module(body=node.body[1:],type_ignores=[]))==ast.dump(ast.Module(body=original.body[index:],type_ignores=[])),
              'entire original unfinished layout/cache/control/output tail unchanged')
    ns=dict(vars(h));ns['restore_completed_focus']=restore_completed_focus
    exec(compile(ast.fix_missing_locations(ast.Module(body=[node],type_ignores=[])),__file__,'exec'),ns)
    return ns['focused'],{'wholeFocusReverseAst':True,'prefixStatementsReplaced':index,
        'wholeUnfinishedTailIdentical':True,'originalMainUnchanged':True,
        'originalNativeBinderUnchanged':True,'nativeBinderCallsInRecovery':0,
        'focusAstSha256':h.h.source.body(ORIGINAL_FOCUS),'mainAstSha256':h.h.source.body(h.main),
        'constructAstSha256':h.h.source.body(h.construct)}


def main():
    if '--regression-directory' in sys.argv:
        index=sys.argv.index('--regression-directory');base=Path(sys.argv[index+1]).resolve();base.relative_to(f.STORE)
        resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(120)
        regression(base);return
    h.focused,wiring=recovered_focus();h.load=load
    h.main()


if __name__=='__main__':main()
