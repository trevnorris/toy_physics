#!/usr/bin/env python3
"""Only missing live-frequency source records, with native saved-record reuse."""
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
import textwrap
import time
from types import SimpleNamespace

import S11c_d_remaining_case_frequency_inputs as inputs
import S11c_d_remaining_case_first_jet_bindings as binding_routes
import S11c_d_remaining_case_coordinate_output_staged as receipts

f, engine, native, q, sp = inputs.f, inputs.engine, inputs.native, inputs.frequency, inputs.sp
PLAN = f.M/'S11c_d_remaining_case_frequency_sources_plan.md'
CP = f.M/'S11c_d_remaining_case_frequency_inputs_checkpoint.json'
CONTEXT_CP = f.M/'S11c_d_remaining_case_coordinate_bindings_focused.json'
SCOPE = ('Only missing frequency-live source bindings and first/second derivatives. '
         'Baseline records and all original scalar bindings remain saved inputs. '
         'No analytic lift, threshold/end construction, numerical row reuse, quadrature, '
         'matrix, mode/current, inverse, scattering solve, contour or physical transcript.')


def reference(base, manifest, path, name, expected):
    path = Path(path); dest = base/name
    f.require(path.is_file() and f.digest(path) == expected, ('accepted source bytes',str(path)))
    f.require(not dest.exists() and not dest.is_symlink(), ('fresh reference',name))
    dest.parent.mkdir(parents=True,exist_ok=True); dest.symlink_to(path)
    f.require(dest.resolve() == path.resolve() and f.digest(dest) == expected,'reference target identity')
    manifest['inputPackets'][str(path)] = expected
    manifest['referencedInputs'][name] = {'original':str(path),'resolvedOriginal':str(path.resolve()),
                                        'sha256':expected,'bytes':path.stat().st_size}


def load(base):
    cp = json.loads(CP.read_text()); origin = Path(cp['runDirectory'])
    f.require(cp['status'] == 'ACCEPTED_CASE_FREQUENCY_SOURCE_INPUTS' and
              f.digest(origin/'checks.json') == cp['checksSha256'],'accepted saved frequency input stage')
    vr = Path(cp['validation']['runDirectory'])
    f.require(f.digest(vr/'checks.json') == cp['validation']['checksSha256'],'accepted input validator')
    receipts.inspect_guard(vr,'validate')
    f.require((vr/'checks.json').read_bytes() == (vr/'validate.stdout').read_bytes(),'validator final stdout identity')
    manifest = {'runDirectory':str(base),'sourceFiles':dict(cp['sourceFiles']),
                'inputPackets':dict(cp['inputPackets']),'referencedInputs':{},
                'input':cp['input'],'settings':cp['settings'],'scope':SCOPE,
                'acceptedInput':{'checkpoint':str(CP),'checksSha256':cp['checksSha256'],
                                 'runDirectory':str(origin),'validatorChecksSha256':cp['validation']['checksSha256']}}
    for name, value in cp['artifacts'].items():
        reference(base,manifest,origin/name,name,value['sha256'])
    for path,name in ((origin/'inputs.json','accepted-inputs.json'),(origin/'checks.json','accepted-input-checks.json'),
                      (CP,'accepted-input-checkpoint.json'),(vr/'checks.json','accepted-input-validation.json')):
        reference(base,manifest,path,name,f.digest(path))
    context = json.loads(CONTEXT_CP.read_text()); cr = Path(context['runDirectory'])
    f.require(context['status'] == 'ACCEPTED_CASE_MATERIAL_BINDING_INPUTS' and
              f.digest(cr/'checks.json') == context['checksSha256'],'accepted original binder contexts')
    for name, value in context['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name) == f.digest(cr/'source'/name) == value,'current/frozen context helper')
        f.require(name not in manifest['sourceFiles'] or manifest['sourceFiles'][name] == value,'shared context helper')
        manifest['sourceFiles'][name] = value
    reference(base,manifest,CONTEXT_CP,'accepted-context-checkpoint.json',f.digest(CONTEXT_CP))
    reference(base,manifest,cr/'checks.json','accepted-context-checks.json',context['checksSha256'])
    for label in cp['cases']:
        for filename in ('input-state.pickle','frame.pickle'):
            name = 'contexts/'+label+'/'+filename
            reference(base,manifest,cr/name,'accepted-'+name,context['artifacts'][name]['sha256'])
    for path in (Path(__file__).resolve(),PLAN,CP,CONTEXT_CP,Path(binding_routes.__file__),Path(receipts.__file__)):
        name = str(path.relative_to(f.ROOT)); value = f.digest(path)
        f.require(name not in manifest['sourceFiles'] or manifest['sourceFiles'][name] == value,'new source pin consistency')
        manifest['sourceFiles'][name] = value
    for name,value in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name) == value,('current source',name))
        dest = base/'source'/name; dest.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(f.ROOT/name,dest)
        f.require(f.digest(dest) == value,'frozen source')
    for name,value in manifest['inputPackets'].items(): f.require(f.digest(Path(name)) == value,('input prehash',name))
    f.save(base/'inputs.json',manifest)
    return manifest,tuple(cp['cases'])


def constructor():
    """Literal complete native record loop; add persistence observations only."""
    tree = ast.parse(inspect.getsource(q.sources)); fn = tree.body[0]
    loop = next(n for n in fn.body if isinstance(n,ast.For) and isinstance(n.target,ast.Tuple)
                and [v.id for v in n.target.elts] == ['key','item'])
    changed = copy.deepcopy(loop); conditional = changed.body[0]
    f.require(isinstance(conditional,ast.If) and ast.unparse(conditional.test) == 'completed is not None and key in completed',
              'actual native completed-record branch')
    observed = []
    for statement in conditional.orelse:
        observed.append(statement)
        if isinstance(statement,ast.Assign) and len(statement.targets)==1 and isinstance(statement.targets[0],ast.Name):
            name = statement.targets[0].id
            if name in ('original','live','expected','actual','control_expected','control_actual','first','second'):
                observed.append(ast.parse(f"checkpoint_operand(base,key,{name!r},{name},item,data)").body[0])
    conditional.orelse = observed
    reverse = copy.deepcopy(changed)
    reverse.body[0].orelse = [n for n in reverse.body[0].orelse if not (
        isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Name)
        and n.value.func.id == 'checkpoint_operand')]
    f.require(ast.dump(reverse) == ast.dump(loop),'whole native source-record loop reverse AST')
    first_native = next(i for i,n in enumerate(fn.body) if isinstance(n,ast.Assign)
                        and any(isinstance(t,ast.Name) and t.id=='native' for t in n.targets))
    start = next(i for i,n in enumerate(fn.body) if isinstance(n,ast.Assign)
                 and any(isinstance(t,ast.Name) and t.id=='records' for t in n.targets))
    stop = fn.body.index(loop)
    body = copy.deepcopy(fn.body[:first_native] + fn.body[start:stop]) + [changed]
    body += ast.parse("return {'records':records,'inventory':inventory,'derivativeCounts':derivative_counts,'respondingFrequencyControls':frequency_controls}").body
    wrapper = ast.FunctionDef(name='construct_records',args=ast.arguments(posonlyargs=[],
        args=[ast.arg(arg=v) for v in ('base','data','completed')],kwonlyargs=[],kw_defaults=[],defaults=[ast.Constant(None)]),
        body=body,decorator_list=[])
    module=ast.fix_missing_locations(ast.Module(body=[wrapper],type_ignores=[]))
    namespace=dict(vars(q),checkpoint_operand=checkpoint_operand)
    exec(compile(module,'<unchanged-frequency-record-loop>','exec'),namespace)
    return namespace['construct_records'],{'wholeNativeRecordLoopReverseAST':True,'addedOperandCheckpointSlots':8,
        'nativeLoopAstSha256':hashlib.sha256(ast.dump(loop).encode()).hexdigest(),
        'wholeNativeSourcesSha256':inputs.body(q.sources),'nativeBindingComparisonSha256':inputs.body(q.binding_comparison),
        'nativeBinderSha256':hashlib.sha256(ast.dump(ast.parse(textwrap.dedent(inspect.getsource(engine.NumericalReducedAction.bind)))).encode()).hexdigest(),
        'nativeCensusSha256':inputs.body(q.census),
        'excludedCompletedPrefix':'Baseline rebinding/local matrix replay and source-character construction are not called.'}


def checkpoint_operand(base,key,slot,value,item,data):
    path=base/'operand-checkpoints'/key;path.mkdir(parents=True,exist_ok=True)
    if slot == 'original':
        data['adapter'].select(key,item)
        f.atomic_pickle(path/'input.pickle',{'key':key,'item':item,'bindingContext':data['contextPath']})
    f.atomic_pickle(path/(slot+'.pickle'),value)


def prohibit():
    def forbidden(*args,**kwargs): raise RuntimeError('completed scientific constructor is disabled in frequency source stage')
    engine.NumericalReducedAction.__init__=engine.ChannelInput.__init__=engine.ReducedPencil.__init__=forbidden
    for module,names in ((q,('load','sources','resume_sources','end_sources','emit_result','main')),
                         (q.p,('load','rebind','baseline_check')),
                         (native,('load','bind_case','operator_grades')),
                         (native.grades,('load','specs','split','check','term_joins')),
                         (inputs.chart,('load','root_chart','source_chart','denominator_chart','end_tables','lift'))):
        for name in names:
            if hasattr(module,name):setattr(module,name,forbidden)
    for name in ('solve','inv','pinv','lstsq','svd','eig','eigh','eigvals','eigvalsh','matrix_rank'):
        setattr(native.np.linalg,name,forbidden)
    native.np.polynomial.legendre.leggauss=forbidden


def adapter(r,pencil,assembly,state):
    a=engine.NumericalReducedAction.__new__(engine.NumericalReducedAction)
    a.r,a.pencil,a.assembly=r,pencil,assembly
    a.input=engine.ChannelInput.__new__(engine.ChannelInput)
    a.input.__dict__.update(copy.deepcopy(state));a.input.r=r
    return a


def restore(base,label,manifest,old):
    case=f.unpickle(base/'accepted-bindings'/label/'case-binding.pickle')
    packets={kind:f.unpickle(base/'accepted-cases'/label/(kind+'.pickle'))
             for kind in ('reduced-action','actions','assembly','factorization')}
    r,dimensions=native.factors.cases.source.restore_context(packets['reduced-action'])
    dimensions.__dict__.update(case['grades']['dimensionState'])
    pencil=SimpleNamespace(r=r,fields=packets['actions']['fields'],probes=packets['actions']['probes'])
    state=f.unpickle(base/'accepted-contexts'/label/'input-state.pickle')
    frame=f.unpickle(base/'accepted-contexts'/label/'frame.pickle')
    f.require(native.same(state,frame['input']) and state['specification']==manifest['input'],'actual saved original input state')
    f.require(native.same(frame['coordinates'],(r.z,r.zp,r.xi,r.ell,r.omega,r.tangents))
              and native.same(frame['profiles'],r.profiles),'actual saved coordinates/profile functions')
    b=case['binding'];g=case['grades']
    for name,value in (('fieldUnits',b['fieldUnits']),('equationUnits',b['equationUnits']),
                       ('settings',b['settings']),('cuts',b['bound']['cutoffBindings']),
                       ('abel',b['bound']['abel']),('profileUnits',b['bound']['profileUnits']),('generators',g['generators'])):
        f.require(native.same(frame[name],value),('saved current binder frame',name))
    seed=adapter(r,pencil,packets['assembly']['result'],state)
    context=binding_routes.binding_context(r,seed,b)
    context.update(specification=state['specification'],unitFrame=state['frame'],fields=pencil.fields,
                   generators=g['generators'],frequency=old['frequency'],referenceFrequency=old['referenceFrequency'])
    f.require(native.same(state['origin'],old['origin']) and state['parameters']['omega']==old['referenceFrequency'],
              'same completed baseline independent grades and frequency seed')
    engine.PHYSICAL_METADATA.dimensions.known[old['frequency']]=(0,-1,0)
    live_state=copy.deepcopy(state);live_state['parameters']['omega']=old['frequency']
    live_state['origin']={v:v for v in state['origin']}
    control_state=copy.deepcopy(state);control_state['parameters']['omega']=sp.Rational(11,10)
    control_state['specification']=copy.deepcopy(state['specification']);control_state['specification']['parameters']['omega']='11/10'
    path=base/'contexts'/label;path.mkdir(parents=True)
    for name,value in (('seed',state),('live',live_state),('control',control_state),('native-input',context)):
        f.atomic_pickle(path/(name+'.pickle'),value)
    live=adapter(r,pencil,packets['assembly']['result'],live_state)
    control=adapter(r,pencil,packets['assembly']['result'],control_state)
    f.require(native.same({k:v for k,v in live_state.items() if k not in ('parameters','origin')},
                          {k:v for k,v in state.items() if k not in ('parameters','origin')}),'only actual formal frequency/grade changes')
    f.require(native.same({k:v for k,v in control_state.items() if k not in ('parameters','specification')},
                          {k:v for k,v in state.items() if k not in ('parameters','specification')}),'only actual independent control frequency changes')
    return case,packets,r,seed,live,control,context


class SeedRouter:
    """Supply completed scalar values to the literal native record loop."""
    def __init__(self,real,case,directory,context_path):
        self.real=real;self.input=real.input;self.records=case['grades']['records'];self.folder=directory
        self.context_path=context_path;self.current=None;self.routes={}
        self.values={key:(expr,unit,value) for key,expr,unit,cut,value in
                     binding_routes.scalar_entries(case['binding'],self.records)}
        for key,item in self.records.items():
            if item['address'][0]!='source':continue
            si=item['address'][1];original=case['binding']['bound']['sources'][0,si];jet=case['binding']['jets'][si]
            f.require(native.same(original['symbolicAmplitude'],item['record']['ORIGINAL']) and
                      native.same(tuple(jet['amplitudeUnit']),tuple(item['record']['UNIT'])),
                      'complete saved generic source amplitude and unit route')
            self.values[key]=(item['record']['ORIGINAL'],item['record']['UNIT'],jet['originalBoundAmplitude'])
        f.require(set(self.values)==set(self.records),'every native seed operand already saved')

    def select(self,key,item):
        inputs.source.require_record(self.records[key],item['address'],item['record']['ORIGINAL'],item['record']['UNIT'])
        self.current=(key,item)

    def bind(self,expression):
        f.require(self.current is not None,'actual native top-level seed call selected')
        key,item=self.current;self.current=None
        f.require(native.same(expression,item['record']['ORIGINAL']),'actual native seed expression call')
        path=self.folder/key;path.mkdir(parents=True,exist_ok=True)
        f.require(key in self.values,'actual native seed call has complete saved operand')
        original,unit,value=self.values[key]
        f.require(native.same((expression,tuple(item['record']['UNIT'])),(original,tuple(unit))),
                  'actual completed scalar expression/unit route')
        kind='accepted-original-source-amplitude' if item['address'][0]=='source' else 'accepted-original-scalar'
        f.atomic_pickle(path/'seed-value.pickle',{'address':item['address'],'expression':expression,
                        'unit':item['record']['UNIT'],'value':value,'route':kind,'context':self.context_path})
        self.routes[key]={'address':item['address'],'kind':kind,'path':str(path/'seed-value.pickle')}
        return value


def construct(base,manifest,labels,run):
    old=f.unpickle(base/'accepted-frequency/frequency-sources.pickle')
    baseline_binding=f.unpickle(base/'accepted-frequency/baseline-binding.pickle')
    owners={};all_cases={};case_inventory={};first_context=None;new_count=0;seed_counts=collections.Counter()
    for key,record in old['records'].items():owners[('accepted-baseline-frequency',native.BASELINE,key)]=record
    for label in labels:
        case,packets,r,seed,live,control,context=restore(base,label,manifest,old)
        if first_context is None:first_context=context
        f.atomic_pickle(base/'contexts'/label/'common-context-pair.pickle',(context,first_context))
        f.require(native.same(context,first_context),'complete actual original/live/control binder context identity')
        if label==native.BASELINE:
            pairs={'local':(case['binding']['local'],baseline_binding['local']),
                   'jets':(case['binding']['jets'],baseline_binding['jets']),
                   'cuts':(case['binding']['bound']['cutoffBindings'],baseline_binding['bound']['cutoffBindings'])}
            f.atomic_pickle(base/'baseline-binding-context-pairs.pickle',pairs)
            for name,(a,b) in pairs.items():f.require(native.same(a,b),('baseline completed binder operand',name))
        candidates=f.unpickle(base/'cases'/label/'frequency-source-candidate-pairs.pickle')
        new={key:case['grades']['records'][key] for key,address,expr,unit,owner in candidates
             if owner[0]['kind']=='uncomputed-frequency-source' and owner[0]['case']==label and owner[0]['key']==key}
        folder=base/'new-records'/label;folder.mkdir(parents=True)
        context_path=str(base/'contexts'/label/'native-input.pickle')
        seed_route=SeedRouter(seed,case,folder/'seed-routes',context_path)
        data={'r':r,'frequency':old['frequency'],'adapter':seed_route,'live':live,'control':control,
              'old':{'domain-binding.pickle':{'bound':case['binding']['bound']}},'packet':{'records':new},'contextPath':context_path}
        f.atomic_pickle(folder/'requested-records.pickle',{'records':new,'context':context,'scope':SCOPE})
        if new:result=run(folder,data)
        else:result={'records':{},'inventory':{},'derivativeCounts':{},'respondingFrequencyControls':0}
        f.save(folder/'seed-routes.json',seed_route.routes)
        for v in seed_route.routes.values():seed_counts[v['kind']]+=1
        for key,record in result['records'].items():owners[('uncomputed-frequency-source',label,key)]=record
        new_count+=len(result['records']);records={};aliases={}
        for key,address,expr,unit,owner in candidates:
            info,owner_expr,owner_unit=owner;where=(info['kind'],info['case'],info['key']);record=owners[where]
            f.require(native.same((expr,tuple(unit)),(record['original'],tuple(record['unit']))),'complete owner source/unit identity')
            f.require(native.same(info['address'],record['address']) and native.same(expr,owner_expr)
                      and native.same(tuple(unit),owner_unit),'physical owner address identity')
            for comparison in record['bindingComparisons'].values():
                f.require(comparison['normalizedResidual']==0 and all(v==0 for v in comparison['proofResiduals']),
                          'completed native source binding proof')
            records[key]=dict(record,address=address,ownerAddress=record['address'],owner=info)
            aliases[key]={'address':address,'owner':info,'context':context_path,'unit':unit}
        target=base/'frequency-cases'/label;target.mkdir(parents=True)
        f.atomic_pickle(target/'records.pickle',records);f.atomic_pickle(target/'record-aliases.pickle',aliases)
        f.atomic_pickle(target/'term-inputs.pickle',case['grades']['termJoins'])
        f.atomic_pickle(target/'original-row-inputs.pickle',case['binding'])
        counts=collections.Counter();dependent=0
        for record in records.values():
            counts[record['address'][0]]+=int(record['firstFrequencyDerivative']!=0)
            dependent+=int(record['frozenFrequencyDifference']!=0)
        summary={'records':len(records),'newRecords':len(new),'reusedRecords':len(records)-len(new),
                 'rows':len(case['binding']['bound']['rows']),'terms':len(case['grades']['termJoins']),
                 'sources':len(case['binding']['jets']),'firstDerivativeCounts':dict(counts),
                 'respondingFrequencyControls':dependent,
                 'nonanalyticRecords':sum(bool(v['census']['nonanalyticNodes']) for v in records.values())}
        f.require(len(new)==len(result['records']) and dependent>0,'complete new record and actual frequency controls')
        packet={'records':records,'aliases':aliases,'frequency':old['frequency'],'referenceFrequency':old['referenceFrequency'],
                'origin':old['origin'],'summary':summary,'generators':case['grades']['generators'],
                'fieldUnits':case['binding']['fieldUnits'],'equationUnits':case['binding']['equationUnits'],
                'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions)),
                'scope':SCOPE,'numericalRowReuseAccepted':False,'sourceFrequencyCharactersPendingSeparateJoin':True}
        f.atomic_pickle(target/'frequency-source.pickle',packet)
        all_cases[label]=str(target/'frequency-source.pickle');case_inventory[label]=summary
        f.save(base/'frequency-case-inventory.json',case_inventory)
        engine.NumericalReducedAction.bind.cache_clear()
        del case,packets,data,seed_route,seed,live,control,records,result
        gc.collect()
    expected=f.unpickle(base/'uncomputed-frequency-operands.pickle')
    f.require(new_count==len(expected),'all and only actual uncomputed source operands')
    result={'cases':all_cases,'summaries':case_inventory,'newSourceRecords':new_count,'seedBindingRoutes':dict(seed_counts),
            'baselineRecords':len(old['records']),'scope':SCOPE,'sourceFiles':manifest['sourceFiles'],
            'inputPackets':manifest['inputPackets']}
    f.atomic_pickle(base/'remaining-case-frequency-sources.pickle',result)
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True);args=parser.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic()
    manifest,labels=load(base);run,join=constructor();f.save(base/'native-source-joins.json',join)
    prohibit();result=construct(base,manifest,labels,run)
    for name,value in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/name)==f.digest(base/'source'/name)==value,'current/frozen posthash')
    for name,value in manifest['inputPackets'].items():f.require(f.digest(Path(name))==value,'original input posthash')
    for name,item in manifest['referencedInputs'].items():
        path=base/name;f.require(path.is_symlink() and str(path.readlink())==item['original'] and str(path.resolve())==item['resolvedOriginal']
            and path.stat().st_size==item['bytes'] and f.digest(path)==item['sha256'],'all original reference pre/post joins')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*')
               if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    checks={**manifest,'status':'COMPLETED_CASE_FREQUENCY_SOURCE_RECORDS','cases':result['summaries'],
            'newSourceRecords':result['newSourceRecords'],'seedBindingRoutes':result['seedBindingRoutes'],
            'baselineRecords':result['baselineRecords'],'nativeJoins':join,'newNumericalWork':0,
            'artifacts':artifacts,'wallSeconds':time.monotonic()-started,
            'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
