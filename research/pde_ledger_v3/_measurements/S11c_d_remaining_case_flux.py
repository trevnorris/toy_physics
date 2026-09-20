#!/usr/bin/env python3
"""Complete case current bookkeeping from accepted response and end operands."""
import argparse,ast,copy,json,resource,shutil,signal,time
from pathlib import Path
import S11c_d_continuum_currents as native
import S11c_d_remaining_case_response as h
import S11c_d_remaining_case_response_finish as output

f,b,m=h.f,h.b,h.m
PLAN=f.M/'S11c_d_remaining_case_flux_plan.md'
RCP=f.M/'S11c_d_remaining_case_response_checkpoint.json'
CCP=f.M/'S11c_d_continuum_currents_checkpoint.json'
BCP=f.M/'S11c_d_continuum_boundary_checkpoint.json'
FCP=f.M/'S11c_d_remaining_case_flux_focused.json'


def emitter(label):
    tree=ast.parse(Path(native.__file__).read_text());emit=copy.deepcopy(h.function(tree,'emit_result'));main=h.function(tree,'main')
    first=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value.func)=='engine.EMISSION_LINES.clear')
    stop=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='summary' for t in n.targets))
    body=copy.deepcopy(main.body[first:stop]);original=copy.deepcopy(body);prefix='s11cd'+label.replace('__','_')+'ContinuumCurrents';changes=[]
    class Keys(ast.NodeTransformer):
        def visit_Constant(self,n):
            if n.value=='s11cdContinuumCurrents':changes.append(True);return ast.Constant(prefix)
            return n
    body=[Keys().visit(n) for n in body];f.require(len(changes)==1,'one actual current export namespace')
    class Undo(ast.NodeTransformer):
        def visit_Constant(self,n):return ast.Constant('s11cdContinuumCurrents') if n.value==prefix else n
    restored=[Undo().visit(copy.deepcopy(n)) for n in body]
    f.require(ast.dump(ast.Module(body=restored,type_ignores=[]))==ast.dump(ast.Module(body=original,type_ignores=[])),'whole original current output and validation tail')
    node=ast.FunctionDef(name='emit_case',args=ast.arguments(posonlyargs=[],args=[ast.arg(n) for n in ('base','result','r','pins','operands','before')],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[]),body=body+[ast.Return(ast.Dict(keys=[ast.Constant(k) for k in ('tags','keys','metadataPaths')],values=[ast.Call(ast.Name('len',ast.Load()),[ast.Name('entries',ast.Load())],[]),ast.Name('keys',ast.Load()),ast.Name('paths',ast.Load())]))],decorator_list=[])
    namespace=dict(vars(native),PREFIX='CONTINUUM_CURRENTS_'+label.replace('__','_'));namespace['emit_result']=h.compile_function(emit,namespace)
    return h.compile_function(node,namespace)


def load(base,resume=None):
    cp,origin=h.checked_checkpoint(RCP,'PUBLISHED_ANNEX_VERIFIED');cc,cr=h.checked_checkpoint(CCP,'PUBLISHED_ANNEX_VERIFIED');bc,br=h.checked_checkpoint(BCP,'PUBLISHED_ANNEX_VERIFIED')
    pins=dict(cp['sourceFiles'])
    for name,sha in cc['sourceFiles'].items():
        f.require(name not in pins or pins[name]==sha,'common current source version');pins[name]=sha
    for p in (Path(__file__).resolve(),PLAN,RCP,CCP,BCP,Path(native.__file__).resolve()):pins[str(p.relative_to(f.ROOT))]=f.digest(p)
    manifest={'runDirectory':str(base),'sourceFiles':pins,'inputPackets':{str(p/'checks.json'):v['checksSha256'] for v,p in ((cp,origin),(cc,cr),(bc,br))},'copiedInputs':{},'settings':cp['settings'],'input':cp['input'],'scope':'Actual four-case finite continuum current/amplitude bookkeeping, positive regulator and approximate boundaries.'}
    labels=tuple(cp['validation'])
    if resume:
        focused=json.loads((resume/'checks.json').read_text());accepted=json.loads(FCP.read_text())
        f.require(accepted['status']=='ACCEPTED_FOUR_CASE_CURRENT_INPUTS' and accepted['checksSha256']==f.digest(resume/'checks.json'),'accepted focused current input routes')
        f.require(focused['mode']=='focused' and focused['status']=='PASSED_FOUR_CASE_CURRENT_INPUTS' and focused['sourceFiles']==pins,'exact focused constructor/sources')
        for name,v in focused['artifacts'].items():m.retain(resume/name,base/name,manifest,v['sha256'])
        manifest['completedFocusedReuse']={'directory':str(resume),'checksSha256':f.digest(resume/'checks.json'),'artifacts':len(focused['artifacts'])}
        manifest['inputPackets'][str(resume/'checks.json')]=f.digest(resume/'checks.json')
    else:
        wanted=[n for n in cp['artifacts'] if n.startswith(('interiors/accepted-cases/','interiors/accepted-bindings/')) or n=='accepted-modes/reference/modal.pickle' or n.startswith('boundary-cases/')]
        for name in wanted:m.retain(origin/name,base/name,manifest,cp['artifacts'][name]['sha256'])
        for label in labels:
            name='cases/'+label+'/continuum/continuum-response.pickle';m.retain(origin/name,base/'case-response'/label/'continuum-response.pickle',manifest,cp['artifacts'][name]['sha256'])
        for name,v in cc['artifacts'].items():m.retain(cr/name,base/'cases'/h.BASELINE/'continuum'/name,manifest,v['sha256'])
        m.retain(cr/'checks.json',base/'cases'/h.BASELINE/'continuum/checks.json',manifest,cc['checksSha256'])
        for end in ('left','right'):
            name=end+'-pencil.pickle';m.retain(br/name,base/'original-pencils'/name,manifest,bc['artifacts'][name]['sha256'])
        m.retain(br/'continuum-boundary.pickle',base/'original-boundary.pickle',manifest,bc['artifacts']['continuum-boundary.pickle']['sha256'])
    for n,v in pins.items():
        path=base/'source'/n;path.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,path);f.require(f.digest(path)==v,'frozen current source')
    f.save(base/'inputs.json',manifest)
    return manifest,labels,cc


def input_routes(base,labels,accepted):
    accepted=json.loads((Path(accepted['runDirectory'])/'checks.json').read_text())
    reference=f.unpickle(base/'accepted-modes/reference/modal.pickle')[0]
    baseline=f.unpickle(base/'original-boundary.pickle');domains=f.unpickle(base/'cases'/h.BASELINE/'continuum/bulk-kinematics.pickle')
    old_source=next(Path(n) for n in accepted['inputPackets'] if Path(n).name=='continuum-response.pickle')
    f.require(f.digest(base/'case-response'/h.BASELINE/'continuum-response.pickle')==accepted['inputPackets'][str(old_source)],'exact original baseline response consumed by current bookkeeping')
    new=f.unpickle(base/'boundary-cases/LAB_HELD__RHOBR_CONSTANT/case-boundary.pickle')['ends']['RIGHT'];routes={}
    for label in labels:
        source=f.unpickle(base/'case-response'/label/'continuum-response.pickle');ends=f.unpickle(base/'boundary-cases'/label/'case-boundary.pickle')
        f.require(source['fieldUnits']==ends['fieldUnits'] and source['rowUnits']==ends['rowUnits'] and source['currentUnit']==ends['currentUnit'],'actual case response/current unit maps')
        selection=native.selectors(source,ends,reference);case_domains={};sources={}
        for end in ('LEFT','RIGHT'):
            if m.same(ends['ends'][end],baseline['ends'][end]):path=base/'original-pencils'/(end.lower()+'-pencil.pickle')
            else:
                f.require(end=='RIGHT' and m.same(ends['ends'][end],new),'actual new RHOBR-right map family');path=base/'boundary-cases/LAB_HELD__RHOBR_CONSTANT/right/right-pencil.pickle'
            pencil=f.unpickle(path);original_waves=tuple(row['wave'] for row in domains[end])
            f.require(m.same(tuple(pencil['waves']),original_waves),'actual acoustic waves before kinematic-domain reuse')
            f.require(not m.same(tuple(2*v for v in pencil['waves']),original_waves),'changed acoustic coefficient rejection')
            case_domains[end]=domains[end];sources[end]={'pencil':str(path.relative_to(base)),'sha256':f.digest(path),'literalWaveJoin':True,'waveCoefficientMutationRejected':True}
        if label==h.BASELINE:
            f.require(m.same(selection,f.unpickle(base/'cases'/label/'continuum/channel-selectors.pickle')),'complete original classifier reuse')
            f.require(m.same(ends['ends'],baseline['ends']),'complete original baseline end input')
        else:
            folder=base/'cases'/label/'continuum';folder.mkdir(parents=True,exist_ok=True)
            f.atomic_pickle(folder/'channel-selectors.pickle',selection);f.atomic_pickle(folder/'bulk-kinematics.pickle',case_domains)
            emitter(label)
        routes[label]={'openRanks':selection['ranks'],'candidateCounts':{e:len(ends['ends'][e]['census']) for e in ('LEFT','RIGHT')},'waveSources':sources,'currentUnit':list(map(str,source['currentUnit'])),'fullSelectedDirections':sum(sum(len(v['columns']) for v in rows) for rows in selection['census'].values())}
    f.save(base/'input-routes.json',routes);return routes


def construct(base,label,manifest):
    source=f.unpickle(base/'case-response'/label/'continuum-response.pickle');ends=f.unpickle(base/'boundary-cases'/label/'case-boundary.pickle');folder=base/'cases'/label/'continuum'
    selection=f.unpickle(folder/'channel-selectors.pickle');bulk=f.unpickle(folder/'bulk-kinematics.pickle')
    metrics,mres=native.open_metrics(source,ends,selection);opened,mutation=native.construct_open(source,metrics,selection)
    f.atomic_pickle(folder/'open-channel-currents.pickle',{'metrics':metrics,'records':opened,'residuals':mres,'mutation':mutation})
    closed={}
    for end in ('LEFT','RIGHT'):
        indices=[i for row in selection['census'][end] if row['kind']=='evanescent' and row['direction']=='outgoing' for i in row['columns']]
        closed[end]=native.amplitude_bookkeeping({g:a[indices] for g,a in source['response']['outgoingFieldAmplitudes'][end].items()})
    finite=native.end_currents(source,ends);f.atomic_pickle(folder/'finite-end-currents.pickle',finite)
    result={'selection':selection,'metrics':metrics,'metricResiduals':mres,'open':opened,'closedAmplitudes':closed,'endCurrents':finite,'bulkDomains':bulk,
            'currentMutation':mutation,'currentUnit':source['currentUnit'],'ratio':source['ratio'],'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],
            'settings':source['settings'],'dimensionState':source['dimensionState'],
            'scope':'Case finite retained-rectangle currents and matching amplitudes. Bulk normal current integrated over depth is not depth-boundary escape; no physical tiny gain/loss or pole result.'}
    f.atomic_pickle(folder/'continuum-currents.pickle',result)
    summary={'openRanks':selection['ranks'],'bulkRealPropagationSets':{end:[str(row['realPropagationSet']) for row in values] for end,values in bulk.items()},
             'openResidualMaximum':max(b.norm(v['residuals']) for v in opened.values()),'finiteEndResidualMaximum':max(b.norm(v['residuals']) for v in finite.values()),
             'currentMutationMaximum':b.norm(mutation),'newQuadratureNodes':0,'newScatteringSolves':0}
    f.save(folder/'current-checks.json',summary);return summary


def emit_saved(base,label,manifest):
    directory=base/'cases'/label/'continuum';result=f.unpickle(directory/'continuum-currents.pickle');case=f.unpickle(base/'interiors/accepted-bindings'/label/'case-binding.pickle')
    r,adapter,packets=h.matrices.context(base/'interiors',label,case,manifest['input']);f.engine.PHYSICAL_METADATA.dimensions.__dict__.update(result['dimensionState'])
    emitted=emitter(label)(directory,result,r,manifest['sourceFiles'],manifest['inputPackets'],f.digest(directory/'continuum-currents.pickle'))
    f.save(directory/'emission-checks.json',emitted);return emitted


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);ap.add_argument('--focused',action='store_true');ap.add_argument('--resume-focused',type=Path);args=ap.parse_args()
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900)
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    manifest,labels,accepted=load(base,args.resume_focused)
    routes=json.loads((args.resume_focused/'input-routes.json').read_text()) if args.resume_focused else input_routes(base,labels,accepted)
    counts={};emissions={};combined=None
    if not args.focused:
        for label in labels:
            if label==h.BASELINE:continue
            counts[label]=construct(base,label,manifest);f.save(base/'case-inventory.json',counts)
        for label in labels:
            if label!=h.BASELINE:emissions[label]=emit_saved(base,label,manifest)
        with (base/'original-combined.out').open('xb') as stream:
            for label in labels:stream.write((base/'cases'/label/'continuum/full.out').read_bytes())
        combined=output.aggregate(base,labels);f.save(base/'aggregation-checks.json',combined)
        f.atomic_pickle(base/'remaining-case-flux.pickle',{'cases':counts,'inputRoutes':routes,'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'scope':manifest['scope']})
    for n,v in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(base/'source'/n)==v,'current source pre/post identity')
    for n,v in manifest['inputPackets'].items():f.require(f.digest(Path(n))==v,'original current input pre/post identity')
    for n,v in manifest['copiedInputs'].items():f.require(f.digest(base/n)==v,'copied current input pre/post identity')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts and p.name not in ('inputs.json','checks.json')}
    # Baseline checks are an explicit consumed artifact, not the final job checks.
    p=base/'cases'/h.BASELINE/'continuum/checks.json';artifacts[str(p.relative_to(base))]={'sha256':f.digest(p),'bytes':p.stat().st_size}
    checks={**manifest,'status':'PASSED_FOUR_CASE_CURRENT_INPUTS' if args.focused else 'COMPLETED_FOUR_CASE_CURRENT_BOOKKEEPING','mode':'focused' if args.focused else 'construct',
            'inputRoutes':routes,'cases':counts,'emissions':emissions,'combinedOutput':combined,'artifacts':artifacts,'newCaseContractions':len(counts),'baselineContractionsRepeated':0,'wallSeconds':time.monotonic()-start}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
