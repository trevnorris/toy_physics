#!/usr/bin/env python3
"""Catalogue actual saved end-family inputs before frequency continuation."""
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

import S11c_d_remaining_case_frequency_analytic_inputs as inputs
import S11c_d_frequency_end as end_native

source=inputs.source
f,native,engine,sp=inputs.f,inputs.native,inputs.engine,inputs.sp
CP=f.M/'S11c_d_remaining_case_frequency_analytic_checkpoint.json'
PLAN=f.M/'S11c_d_remaining_case_frequency_end_inputs_plan.md'
BASELINE=native.BASELINE
SCOPE=('Saved end-family source/context/unit/current/mode/continuation input catalogue. '
       'Candidate identity only; no live end binding, analytic lift, derivative, rational table, '
       'seed/current/mode construction, continuation point, matrix, solve or physical output. '
       'Numerical end and row reuse require subsequent complete frequency input joins.')


def same(a,b):f.require(native.same(a,b),'full typed saved end input identity')


def checkpoint(base,manifest,key,name,status,names):
    path=f.M/name;cp=json.loads(path.read_text());origin=Path(cp['runDirectory'])
    f.require(cp['status']==status and f.digest(origin/'checks.json')==cp['checksSha256'],('accepted end input checkpoint',name))
    for name,value in cp['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==f.digest(origin/'source'/name)==value,('accepted current/frozen source',name))
        f.require(name not in manifest['sourceFiles'] or manifest['sourceFiles'][name]==value,'shared source versions')
        manifest['sourceFiles'][name]=value
    source.reference(base,manifest,path,'end-origins/'+key+'-checkpoint.json',f.digest(path))
    source.reference(base,manifest,origin/'checks.json','end-origins/'+key+'-checks.json',cp['checksSha256'])
    if (origin/'inputs.json').is_file():source.reference(base,manifest,origin/'inputs.json','end-origins/'+key+'-inputs.json',f.digest(origin/'inputs.json'))
    for name in names(cp):
        item=cp['artifacts'][name]
        source.reference(base,manifest,origin/name,'end-accepted/'+key+'/'+name,item['sha256'])
    return cp


def load(base):
    cp=json.loads(CP.read_text());origin=Path(cp['runDirectory']);vr=Path(cp['validation']['runDirectory'])
    f.require(cp['status']=='ACCEPTED_CASE_FREQUENCY_ANALYTIC_SOURCES' and f.digest(origin/'checks.json')==cp['checksSha256'],'accepted full analytic sources')
    source.receipts.inspect_guard(vr,'validate')
    f.require(f.digest(vr/'checks.json')==cp['validation']['checksSha256'] and
              (vr/'checks.json').read_bytes()==(vr/'validate.stdout').read_bytes(),'accepted final analytic validator')
    manifest={'runDirectory':str(base),'sourceFiles':dict(cp['sourceFiles']),'inputPackets':dict(cp['inputPackets']),
              'referencedInputs':{},'input':cp['input'],'settings':cp['settings'],'scope':SCOPE,
              'acceptedAnalyticSources':{'checkpoint':str(CP),'checksSha256':cp['checksSha256'],
                  'validatorChecksSha256':cp['validation']['checksSha256'],'runDirectory':str(origin)}}
    for name,item in cp['artifacts'].items():source.reference(base,manifest,origin/name,name,item['sha256'])
    for name,item in cp['referencedInputs'].items():
        path=origin/name
        f.require(path.is_symlink() and str(path.readlink())==item['original'] and
                  str(path.resolve())==item['resolvedOriginal'] and path.stat().st_size==item['bytes']
                  and f.digest(path)==item['sha256'],'accepted analytic reference address identity')
    for path,name in ((CP,'accepted-analytic-checkpoint.json'),(origin/'checks.json','accepted-analytic-checks.json'),
                      (origin/'inputs.json','accepted-analytic-inputs.json'),(vr/'checks.json','accepted-analytic-validation.json')):
        source.reference(base,manifest,path,name,f.digest(path))
    labels=tuple(cp['cases'])
    checkpoint(base,manifest,'sources','S11c_d_remaining_case_end_sources_checkpoint.json','ACCEPTED_CASE_END_AND_CURRENT_INPUT_SOURCES',
        lambda c:['remaining-case-end-sources.pickle']+['cases/'+label+'/case-end-sources.pickle' for label in labels])
    checkpoint(base,manifest,'uniform','S11c_d_remaining_case_uniform_checkpoint.json','PUBLISHED_ANNEX_VERIFIED',
        lambda c:['remaining-case-currents.pickle','input-routes.json','accepted-uniform/uniform-response.pickle']+
        ['cases/'+label+'/continuum/uniform-response.pickle' for label in labels]+
        ['cases/'+label+'/'+end.lower()+'/'+name for label in labels for end in ('REFERENCE','LEFT','RIGHT') for name in ('mode-inputs.pickle','modal.pickle')])
    checkpoint(base,manifest,'frequency','S11c_d_frequency_source_checkpoint.json','PUBLISHED_ANNEX_VERIFIED',
        lambda c:[end.lower()+'-'+name for end in ('REFERENCE','LEFT','RIGHT') for name in ('frequency-pencil.pickle','threshold-candidates.pickle')])
    checkpoint(base,manifest,'chart','S11c_d_frequency_chart_checkpoint.json','PUBLISHED_ANNEX_VERIFIED',
        lambda c:['left-rational-end.pickle','right-rational-end.pickle'])
    checkpoint(base,manifest,'continuation','S11c_d_frequency_end_focused.json','ACCEPTED_FOCUSED_END_CONTINUATION',lambda c:list(c['artifacts']))
    for path in (Path(__file__).resolve(),PLAN,CP,Path(end_native.__file__)):
        name=str(path.relative_to(f.ROOT));value=f.digest(path)
        f.require(name not in manifest['sourceFiles'] or manifest['sourceFiles'][name]==value,'new end input source pin')
        manifest['sourceFiles'][name]=value
    for name,value in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==value,'current input helper')
        if name in cp['sourceFiles']:f.require(f.digest(origin/'source'/name)==value,'accepted frozen analytic helper')
        path=base/'source'/name;path.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/name,path)
        f.require(f.digest(path)==value,'frozen end input helper')
    for name,value in manifest['inputPackets'].items():f.require(f.digest(Path(name))==value,('input prehash',name))
    f.save(base/'inputs.json',manifest)
    return manifest,labels


def native_join():
    names={'endSources':source.q.end_sources,'rationalEndTables':inputs.chart.end_tables,
           'endSeeds':end_native.seeds,'invariantPair':end_native.Pair,'endMaps':end_native.maps,
           'continuePair':end_native.continue_pair,'endFocused':end_native.focused,
           'analyticPencil':engine.FullPencilModes.analytic,'channelMapping':engine.ChannelInput.mapping}
    result={}
    for name,value in names.items():
        text=textwrap.dedent(inspect.getsource(value));path=Path(inspect.getsourcefile(value))
        result[name]={'sourceFile':str(path),'sourceSha256':f.digest(path),
                      'wholeBodyAST':hashlib.sha256(ast.dump(ast.parse(text)).encode()).hexdigest(),'source':text}
    return {'nativeBodies':result,'nativeConstructorsCalled':False,
            'sourceSignatureFields':['strong','weak','curl','weakUnits','profileBindings','sourceEnergy','constraint','mass',
                                     'context','seedInput','frame','fieldUnits','equationUnits','end'],
            'interpretation':'Literal full saved input candidate catalogue; no live end operator or numerical reuse conclusion.'}


def prohibit():
    inputs.prohibit()
    def forbidden(*args,**kwargs):raise RuntimeError('end input catalogue cannot reconstruct accepted or new science')
    for name in ('load','seeds','maps','continue_pair','focused','main','polynomial'):
        setattr(end_native,name,forbidden)
    for name in ('__init__','initial','coefficients','equation','jacobian','solve'):setattr(end_native.Pair,name,forbidden)
    for cls,names in ((engine.FullPencilModes,('__init__','analytic','solve','rational_determinant')),
                      (engine.ChannelInput,('__init__','mapping'))):
        for name in names:
            if hasattr(cls,name):setattr(cls,name,forbidden)
    engine.NumericalReducedAction.bind=forbidden
    for name in ('diff','integrate','lambdify','simplify','cancel','expand','factor','solve','resultant'):
        setattr(sp,name,forbidden)
    sp.Basic.diff=forbidden
    for name in ('solve','inv','svd','eig','eigh','expm'):setattr(end_native.la,name,forbidden)


def prepare(base,manifest,labels):
    root=base/'end-accepted'
    sources=f.unpickle(root/'sources/remaining-case-end-sources.pickle')
    currents=f.unpickle(root/'uniform/remaining-case-currents.pickle')
    baseline_response=f.unpickle(root/'uniform/accepted-uniform/uniform-response.pickle')
    saved_seeds=f.unpickle(root/'continuation/end-seeds.pickle')
    routes=json.loads((root/'uniform/input-routes.json').read_text())
    baseline_context=f.unpickle(base/'accepted/contexts'/BASELINE/'native-input.pickle')
    signatures={};groups=[];inventory={};pending=[]
    for label in labels:
        directory=base/'end-input-cases'/label;directory.mkdir(parents=True)
        context=f.unpickle(base/'accepted/contexts'/label/'native-input.pickle')
        state=f.unpickle(base/'accepted/contexts'/label/'seed.pickle')
        frame=f.unpickle(base/'accepted/accepted-contexts'/label/'frame.pickle')
        analytic=f.unpickle(base/'analytic-cases'/label/'frequency-analytic.pickle')
        case=f.unpickle(root/'sources/cases'/label/'case-end-sources.pickle')
        uniform=sources['sources'][label]
        response=f.unpickle(root/'uniform/cases'/label/'continuum/uniform-response.pickle')
        f.atomic_pickle(directory/'context-pairs.pickle',{'context':context,'baselineContext':baseline_context,
            'seedInput':state,'frame':frame,'analyticVariables':analytic['variables'],
            'ownSourcePath':str(root/'sources/cases'/label/'case-end-sources.pickle'),
            'analyticPacketPath':str(base/'analytic-cases'/label/'frequency-analytic.pickle')})
        same(context,baseline_context);same(state,frame['input']);same(state['specification'],manifest['input'])
        same(case['source'],uniform);same(case['currentInputs'],sources['currentInputs'][label]);f.require(case['case']==label,'own physical end source case')
        same(analytic['variables']['frequency'],context['frequency']);same(analytic['variables']['referenceFrequency'],state['parameters']['omega'])
        same(analytic['variables']['origin'],state['origin'])
        same(tuple(response['fieldUnits']),tuple(frame['fieldUnits']))
        summary={}
        for end in ('LEFT','RIGHT'):
            target=directory/end.lower();target.mkdir()
            record=uniform['records'][end];ci=case['currentInputs'][end];cv=currents['cases'][label][end]
            inp=f.unpickle(root/'uniform/cases'/label/end.lower()/'mode-inputs.pickle')
            modal,known=f.unpickle(root/'uniform/cases'/label/end.lower()/'modal.pickle')
            background=response['backgrounds'][end];channel=background['response']['channels'][end]
            signature={'strong':record['strong'],'weak':record['weak'],'curl':uniform['curl'],'weakUnits':uniform['units']['weak'],
                'profileBindings':ci['profileBindings'],'sourceEnergy':ci['sourceEnergy'],
                'constraint':ci['MATERIAL_CONSTRAINT_COEFFICIENTS'],'mass':ci['ZERO_TRANSFER_MASS_ROW'],
                'context':context,'seedInput':state,'frame':frame,'fieldUnits':frame['fieldUnits'],
                'equationUnits':frame['equationUnits'],'end':end}
            raw={'address':(label,end),'signature':signature,'sourceRecord':record,'currentInput':ci,
                 'modeInput':inp,'knownUnits':known,'uniformBackground':background,
                 'fullCurrentPacketPath':str(root/'uniform/remaining-case-currents.pickle'),
                 'fullModalPacketPath':str(root/'uniform/cases'/label/end.lower()/'modal.pickle'),
                 'acceptedUniformInputRoute':routes['routes'][label+'__'+end]}
            f.atomic_pickle(target/'full-inputs.pickle',raw)
            same(ci['strong'],record['strong']);same(cv['strong'],record['strong']);same(ci['profileBindings'],state['limits'])
            same((inp['case'],inp['end']),(label,end));same(inp['pairing'],cv['pairing']);same(inp['acoustic'],cv['acoustic'])
            same(tuple(inp['fieldUnits']),tuple(frame['fieldUnits']));same(inp['frequency'],state['parameters']['omega'])
            same(modal['NATIVE_RECORDS'],inp['rootPacket']['records'])
            f.require(len(background['modes'])==len(modal['RECORDS'])==18,'all accepted root/lift dispositions')
            selected=channel['outgoing']+channel['incoming'];by_index={}
            for item in selected:by_index.setdefault(item['RECORD_INDEX'],[]).append(item)
            modes={item['info']['INDEX']:item for item in background['modes']}
            selection={'selectedItems':selected,'groups':by_index,'allCandidateModes':background['modes'],
                       'channel':channel,'allRootRecords':inp['rootPacket']['records']}
            f.atomic_pickle(target/'whole-cluster-inputs.pickle',selection)
            for index,items in by_index.items():
                info=modes[index]['info']
                f.require(len(items)==info['NULLITY'] and sorted(x['BASIS_COLUMN'] for x in items)==list(range(info['NULLITY'])),'full selected degenerate cluster')
            f.require(len(by_index)==5 and len(selected)==7,'full five clusters/seven directions per end')
            base_signature=signature if label==BASELINE else signatures[BASELINE,end]
            candidate=native.same(signature,base_signature)
            pair={'actual':signature,'baseline':base_signature,'sameFullSourceInput':candidate,
                  'actualAddress':(label,end),'baselineAddress':(BASELINE,end)}
            f.atomic_pickle(target/'baseline-source-input-pair.pickle',pair)
            matches=[i for i,(_,other) in enumerate(groups) if native.same(signature,other)]
            if matches:group=matches[0]
            else:group=len(groups);groups.append(((label,end),signature))
            if label==BASELINE:
                table=f.unpickle(root/'chart'/(end.lower()+'-rational-end.pickle'))
                pencil=f.unpickle(root/'frequency'/(end.lower()+'-frequency-pencil.pickle'))
                seeds=saved_seeds[end]
                f.atomic_pickle(target/'baseline-accepted-proof-pair.pickle',{'tableSource':table['source'],'frequencyPencil':pencil,
                    'actualBackground':background,'acceptedBackground':baseline_response['backgrounds'][end],
                    'actualChannel':channel,'acceptedSeedChannel':seeds['acceptedChannel'],
                    'seedClusters':seeds['clusters'],'actualSourceInput':signature})
                same(table['source'],pencil);same(background,baseline_response['backgrounds'][end]);same(channel,seeds['acceptedChannel'])
                same(modes,{x['info']['INDEX']:x for x in seeds['allCandidates']})
                for seed in seeds['clusters']:same(seed['items'],by_index[seed['index']]);same(seed['originalMode'],modes[seed['index']])
                same(pencil['frequency'],context['frequency']);same(pencil['origin'],state['origin'])
            if not candidate:pending.append({'address':(label,end),'family':group,'owner':groups[group][0],
                'fullInputPath':str(target/'full-inputs.pickle'),'reason':'Full literal source/context input differs from this end baseline; no reuse is assumed.'})
            changed=copy.deepcopy(signature);mat=sp.MutableDenseMatrix(changed['strong']);mat[0,0]+=1;changed['strong']=sp.ImmutableMatrix(mat)
            unit=copy.deepcopy(signature);unit['fieldUnits']=tuple(unit['fieldUnits']);u=tuple(unit['fieldUnits'][0]);unit['fieldUnits']=((u[0]+1,u[1],u[2]),)+unit['fieldUnits'][1:]
            frequency=copy.deepcopy(signature);frequency['context']['referenceFrequency']+=1
            address=(label,'WRONG_END')
            incoming=copy.deepcopy(selection);incoming['selectedItems'][0]['vector']=-incoming['selectedItems'][0]['vector']
            f.atomic_pickle(target/'mutation-operands.pickle',{'original':signature,'changedCoefficient':changed,'changedUnit':unit,
                'changedFrequency':frequency,'originalAddress':(label,end),'changedAddress':address,
                'originalSelection':selection,'changedSelection':incoming})
            controls={'coefficient':not native.same(changed,signature),'unit':not native.same(unit,signature),
                'frequency':not native.same(frequency,signature),'address':address!=(label,end),'selectedVector':not native.same(incoming,selection)}
            f.save(target/'mutation-controls.json',controls);f.require(all(controls.values()),'actual saved end input mutations respond')
            summary[end]={'sourceInputFamily':group,'familyOwner':groups[group][0],'baselineSourceCandidate':candidate,
                'candidateDispositions':len(modes),'selectedClusters':len(by_index),'selectedDirections':len(selected),
                'actualInputSha256':f.digest(target/'full-inputs.pickle'),'controls':controls,
                'newEndConstruction':False,'completeFrequencyReuseAccepted':False}
            signatures[label,end]=signature
        inventory[label]=summary;f.save(directory/'case-inventory.json',summary);f.save(base/'end-input-case-inventory.json',inventory)
        del analytic,case,response;gc.collect()
    f.atomic_pickle(base/'pending-end-frequency-inputs.pickle',pending)
    f.atomic_pickle(base/'end-source-input-families.pickle',groups)
    return {'cases':inventory,'actualEnds':sum(len(v) for v in inventory.values()),'sourceInputFamilies':len(groups),
        'baselineSourceCandidates':sum(v['baselineSourceCandidate'] for c in inventory.values() for v in c.values()),
        'unmatchedEndUses':len(pending),'newEndConstructions':0,'newContinuationPoints':0,
        'frequencyPencilOrNumericalReuseAccepted':False,'endFamilyConstructionPending':True}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);args=ap.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic()
    inputs.protect_references(base);manifest,labels=load(base)
    join=native_join();f.save(base/'native-end-input-joins.json',join);prohibit();result=prepare(base,manifest,labels)
    for name,value in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/name)==f.digest(base/'source'/name)==value,'current/frozen posthash')
    for name,value in manifest['inputPackets'].items():f.require(f.digest(Path(name))==value,'original input posthash')
    for name,item in manifest['referencedInputs'].items():
        path=base/name;f.require(path.is_symlink() and str(path.readlink())==item['original'] and str(path.resolve())==item['resolvedOriginal']
            and path.stat().st_size==item['bytes'] and f.digest(path)==item['sha256'],'full reference identity')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*')
               if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    checks={**manifest,**result,'status':'COMPLETED_CASE_FREQUENCY_END_INPUTS','nativeJoins':join,'artifacts':artifacts,
            'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
