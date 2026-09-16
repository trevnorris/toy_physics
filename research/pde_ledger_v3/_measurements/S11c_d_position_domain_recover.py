#!/usr/bin/env python3
"""Finish position-domain preflight validation from immutable saved operands."""
import argparse,ast,copy,hashlib,inspect,json,resource,shutil,textwrap,time
from pathlib import Path
import numpy as np
import sympy as sp
import S11c_d_position_domain as domain
from S11c_d_position_domain import ROOT,STORE,engine,require,digest,save,atomic_pickle,unpickle
from S11c_d_output_codec import decoded_lines
from ledger_fold import _restore


def repair_join(origin):
    path=Path(domain.__file__).resolve();old=ast.parse((origin/'source'/path.relative_to(ROOT)).read_text());new=ast.parse(path.read_text())
    names={'cache_for_tag','replay_adapter'}
    additions=[n for n in new.body if isinstance(n,ast.FunctionDef) and n.name in names]
    require({n.name for n in additions}==names,'repair helper census')
    for node in additions:new.body.remove(node)
    a=next(n for n in old.body if isinstance(n,ast.Import));b=next(n for n in new.body if isinstance(n,ast.Import))
    require({n.name for n in b.names}-{n.name for n in a.names}=={'ast','hashlib','inspect','textwrap'},'repair imports')
    b.names=copy.deepcopy(a.names)
    prior=next(n for n in old.body if isinstance(n,ast.FunctionDef) and n.name=='finish')
    current=next(n for n in new.body if isinstance(n,ast.FunctionDef) and n.name=='finish')
    start=next(i for i,n in enumerate(current.body) if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='(replay, method_join)')
    require(ast.unparse(current.body[start].value)=='replay_adapter()','repair replay call')
    stop=next(i for i,n in enumerate(current.body) if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='(entries, keys, paths)')
    require(stop==start+3 and ast.unparse(current.body[stop].value)=="replay(base, result, data['r'], data['bound'], data['bound']['rows'], data['finest'], data['provenance'])",'repair replay invocation')
    prior_start=next(i for i,n in enumerate(prior.body) if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='(previous_emit, previous_prefix)')
    prior_stop=next(i for i,n in enumerate(prior.body) if isinstance(n,ast.Try))
    require(ast.dump(current.body[start+1])==ast.dump(prior.body[prior_start+2]) and ast.dump(current.body[start+2])==ast.dump(prior.body[prior_start+3]),'unchanged encoder resets')
    current.body[start:stop+1]=copy.deepcopy(prior.body[prior_start:prior_stop+1])
    summary=next(n.value for n in current.body if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='summary')
    indexes=[i for i,k in enumerate(summary.keys) if isinstance(k,ast.Constant) and k.value=='replayMethodJoin'];require(len(indexes)==1,'repair summary field')
    i=indexes[0];require(ast.unparse(summary.values[i])=='method_join','repair summary operand');summary.keys.pop(i);summary.values.pop(i)
    require(ast.dump(new)==ast.dump(old),'whole checker unchanged outside cache adapter wiring')
    _,replay=domain.replay_adapter()
    return {'originalCheckerSha256':digest(origin/'source'/path.relative_to(ROOT)),'repairedCheckerSha256':digest(path),
        'restoredWholeCheckerAstSha256':hashlib.sha256(ast.dump(new).encode()).hexdigest(),'originalWholeCheckerAstSha256':hashlib.sha256(ast.dump(old).encode()).hexdigest(),'nativeReplayJoin':replay}


def metadata_fields(body):
    for record in body:
        if len(record)==2 and isinstance(record[0],sp.Tuple):yield {str(k):v for k,v in record[1]}
        else:yield {str(k):v for k,v in record}


def metadata_guard():
    # Exercise the actual restored native guard, including its cache lookup,
    # against altered metadata. No parallel handwritten unit validator.
    source=ast.parse(textwrap.dedent(inspect.getsource(domain.native.emit_and_replay))).body[0]
    start=next(i for i,n in enumerate(source.body) if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='metadata_paths')
    function=ast.parse('def validate(entries,result,r,bound):\n    pass').body[0]
    function.body=copy.deepcopy(source.body[start:-1])+[ast.Return(value=ast.Name(id='metadata_paths',ctx=ast.Load()))]
    lookup=next(n for n in ast.walk(function) if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='cache')
    require(ast.unparse(lookup.value)=="result['cacheChecks'][int(tag.rsplit('_', 1)[1])]",'native metadata guard lookup')
    lookup.value=ast.parse('cache_for_tag(result,tag)',mode='eval').body
    namespace=dict(vars(domain.native),cache_for_tag=domain.cache_for_tag)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[function],type_ignores=[])),'<native-metadata-negative-controls>','exec'),namespace)
    return namespace['validate']


def validate_saved(origin,base):
    preflight=json.loads((origin/'preflight.json').read_text());bound_packet=unpickle(origin/'bound-position-domains.pickle');packet=unpickle(origin/'position-domain.pickle')
    require(packet['sourceFiles']==bound_packet['sourceFiles']==preflight['sourceFiles'] and packet['provenance']==bound_packet['provenance']==preflight['provenance'],'saved packet/source joins')
    changed=str(Path(domain.__file__).resolve().relative_to(ROOT))
    for n,h in preflight['sourceFiles'].items():
        require(digest(origin/'source'/n)==h,'original frozen source')
        if n!=changed:require(digest(ROOT/n)==h,'unchanged current source')
    join=repair_join(origin)
    artifacts={str(p.relative_to(origin)):{'sha256':digest(p),'bytes':p.stat().st_size} for p in origin.rglob('*') if p.is_file()}
    save(base/'original-artifacts.json',artifacts)
    # Preserve every original payload byte-for-byte. Re-emission has its own file.
    for relative in artifacts:
        source=origin/relative;target=base/relative
        if relative=='full.out':target=base/'original-full.out'
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
    cp,source_base=domain.momentum.source.accepted(domain.momentum.source.SOURCE)
    r,dimensions=domain.momentum.source.native.source.restore_context(unpickle(source_base/'reduced-action.pickle'))
    dimensions.__dict__.update(packet['dimensionState'])
    domains=bound_packet['domains'];bound=domains[0];result=packet['result'];references=bound_packet['references']
    require(result['joins']==bound_packet['joins'],'domain join packets')
    for index,choice in enumerate(domain.CHOICES):
        rebound,proof=domain.rebind(r,bound,choice)
        require(rebound['rows']==domains[index]['rows'] and rebound['profileUnits']==domains[index]['profileUnits'] and proof==result['joins'][index],'all cutoff/integrand/reverse joins')
    focused=json.loads((origin/'focused-before-smoke.json').read_text());require(len(focused['transformArtifacts'])==82 and len(focused['prefixArtifacts'])==6,'focused census')
    transform_evidence=[]
    for a in focused['transformArtifacts']:
        require(digest(origin/a['path'])==a['sha256'],'transform artifact')
        v=unpickle(origin/a['path']);require(np.array_equal(v['sourceResidual'],v['values']-v['direct']) and np.array_equal(v['orderChange'],v['refined']-v['values']) and np.array_equal(v['mutationResidual'],v['mutationValues']-v['values']),'transform residual operands')
        require(np.all(np.isfinite(v['values'])) and np.all(np.isfinite(v['refined'])) and np.max(abs(v['mutationResidual']))>0,'transform finiteness/sensitivity')
        scale=1+abs(v['refined']);require(np.max(abs(v['sourceResidual'])/scale)<1e-9 and np.max(abs(v['orderChange'])/scale)<1e-9,'transform agreement')
        if 'sourceIndex' in v:
            original=bound['sources'][(v['test'],v['sourceIndex'])]
            require(v['original']==original['boundSource'] and v['unit']==original['integralUnit'] and np.array_equal(v['frequencies'],original['frequencies']) and np.array_equal(v['assignments'],original['assignments']),'native source/field/frequency joins')
        else:
            proof=result['joins'][v['index']]['profiles'][v['profile']]
            require(v['integral']==proof['new'] and v['unit']==proof['unit'] and v['range']['residual']==0,'profile unit/operand/range joins')
        transform_evidence.append({'path':a['path'],'maxScaledSourceResidual':float(np.max(abs(v['sourceResidual'])/scale)),'maxScaledOrderChange':float(np.max(abs(v['orderChange'])/scale))})
    for a in focused['prefixArtifacts']:
        require(digest(origin/a['path'])==a['sha256'],'prefix hash');v=unpickle(origin/a['path']);x,y=v['original'],v['rebound']
        require(np.array_equal(v['valueResidual'],x['values']-y['values']) and not np.any(v['valueResidual']) and np.array_equal(v['mutationResidual'],x['measureMutationValues']-y['measureMutationValues']) and not np.any(v['mutationResidual']) and v['massResidual']==x['quadratureMass']-y['quadratureMass']==0,'prefix residuals')
        require(x['rowIndices']==y['rowIndices'] and x['variables']==y['variables'] and x['nodeCount']==y['nodeCount'] and x['setting']==y['setting'],'prefix assignments')
    manifest=json.loads((origin/'workers.json').read_text());require(manifest['status']=='completed' and len(manifest['outcomes'])==4 and all(v['exitCode']==v['stderrBytes']==0 for v in manifest['outcomes']),'all worker outcomes')
    worker_evidence=[];coarse_artifacts=[]
    for key,h in manifest['resultFiles'].items():
        directory=origin/f'worker-{key}';w=unpickle(directory/'result.pickle');ti,index=w['task']
        require(digest(directory/'result.pickle')==h==json.loads((directory/'checks.json').read_text())['resultSha256'] and (directory/'stdout').stat().st_size==(directory/'stderr').stat().st_size==0,'worker hash/logs')
        require(w['sourceFiles']==preflight['sourceFiles'] and w['provenance']==preflight['provenance'] and w['settings']==domain.settings(references[ti],domain.CHOICES[index],True),'worker provenance/settings')
        indices=[i for g in w['groups'] for i in g['rowIndices']];require(len(indices)==len(set(indices))==80 and set(indices)==set(range(80)),'worker complete row census')
        require(set(w['telemetry']['sourceFrequencyCensus'])=={(ti,si) for t,si in bound['sources'] if t==ti} and w['telemetry']['evaluatedProfileIntegrals']==set(domains[index]['profileUnits']),'complete worker source/profile census')
        for cache in w['telemetry']['cacheChecks']:require(not cache['nodesResidual'].any() and not cache['weightsResidual'].any() and cache['readOnly'],'exact cache residual')
        for g in w['groups']:
            require(abs(g['volumeResidual'])<1e-10*(1+abs(g['boxVolume'])) and np.all(np.isfinite(g['values'])) and not np.any((abs(g['values'])>1e-9)&(abs(g['measureMutationResidual'])<1e-12)),'mass/actual measure')
            require(g['test']==ti and g['rowIndices']==[row['index'] for row in domains[index]['rows'] if tuple(l[0] for l in row['limits'])==tuple(g['variables'])],'row/limit/field join')
        file=origin/f'coarse-{ti}-{index}.pickle';coarse=unpickle(file)
        require(np.array_equal(coarse['workerAction'],w['action']) and np.array_equal(coarse['residual'],coarse['directAction']-w['action']) and not np.any(coarse['residual']),'independent complete action residual')
        for g,hg in zip(coarse['groups'],w['groups']):require(np.array_equal(g['values'],hg['values']) and np.array_equal(g['measureMutationValues'],hg['measureMutationValues']) and g['quadratureMass']==hg['quadratureMass'],'serial/worker raw operands')
        coarse_artifacts.append(domain.artifact(origin,file));worker_evidence.append({'task':w['task'],'rowCount':len(indices),'sourceCount':len(w['telemetry']['sourceFrequencyCensus']),'profileCount':len(w['telemetry']['evaluatedProfileIntegrals']),'maxMassResidual':max(abs(g['volumeResidual']) for g in w['groups'])})
    for a in (*result['recordArtifacts'],*result['workerArtifacts']):require(digest(origin/a['path'])==a['sha256'],'record/worker artifacts')
    for item in result['records']:
        replay=domain.assemble(r,domains[item['index']],references[item['test']],item['groups'])
        require(np.array_equal(replay['action'],item['action']),'all-row saved action reconstruction')
        for a,b in zip(replay['terms'],item['terms']):require(a['index']==b['index'] and a['rows']==b['rows'] and np.array_equal(a['values'],b['values']) and all(b['changed']),'all native term joins')
        if item['index']:
            previous=next(v for v in result['records'] if v['test']==item['test'] and v['index']==item['index']-1)
            require(np.array_equal(item['precedingActionDifference'],item['action']-previous['action']),'saved finite action change')
            for a,b,c in zip(item['integralChanges'],item['groups'],previous['groups']):require(np.array_equal(a,b['values']-c['values']),'saved raw integral change')
            for a,b,c in zip(item['termChanges'],item['terms'],previous['terms']):require(np.array_equal(a,np.asarray(b['values'])-np.asarray(c['values'])),'saved raw term change')
    # Recheck only the literal domain guard, with no quadrature execution.
    worker=engine.BoundedSourceFourierQuadrature.ThreeMomentum(bound['rows'],bound['sources'],r)
    wrong=next(iter(bound['profileUnits']));rejected=False
    try:worker.profile_value(wrong,{},domain.settings(references[0],domain.CHOICES[2])[3],{})
    except ValueError as e:rejected=str(e)=='nested profile cutoff mismatch'
    require(rejected,'wrong profile limit rejected')
    return r,bound,result,references,preflight,artifacts,join,focused,transform_evidence,worker_evidence,coarse_artifacts


def main():
    p=argparse.ArgumentParser();p.add_argument('--origin',type=Path,required=True);p.add_argument('--run-directory',type=Path,required=True);args=p.parse_args()
    origin=args.origin.resolve();base=args.run_directory.resolve();origin.relative_to(STORE);base.relative_to(STORE);base.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    r,bound,result,references,preflight,original_artifacts,join,focused,transforms,workers,coarse=validate_saved(origin,base)
    before={p.name:digest(p) for p in base.glob('*.pickle')}
    replay,method=domain.replay_adapter();engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=engine.PayloadEncoder()
    entries,keys,paths=replay(base,result,r,bound,bound['rows'],references,preflight['provenance'])
    original={}
    for line in decoded_lines(origin/'full.out'):
        tag,_,body=line.rstrip('\n').partition(': ');require(tag not in original,'original duplicate tag');original[tag]=_restore(body)
    differences=[tag for tag in set(original)|set(entries) if original.get(tag)!=entries.get(tag)]
    save(base/'emission-differences.json',differences);require(not differences,'original/recovery payload identity')
    controls=[];guard=metadata_guard();require(guard(entries,result,r,bound)==paths,'native metadata guard path census')
    for tag,body in entries.items():
        if tag.startswith('PY_S11CD_METADATA_'+domain.PREFIX+'_CACHE_') and any('_'+name+'_' in tag for name in ('NODESRESIDUAL','WEIGHTSRESIDUAL')):
            cache=domain.cache_for_tag(result,tag);expected=tuple(engine.PHYSICAL_METADATA.dimensions.measure(r.xi if cache['profileRule'] else r.zp))
            for fields in metadata_fields(body):
                actual=tuple(fields['DIMENSION_L_T_M']);require(actual==expected,'exact cache unit')
            replacements={v:sp.Tuple(v[0],sp.Tuple(v[1][0]+1,*v[1][1:])) for v in sp.preorder_traversal(body) if isinstance(v,sp.Tuple) and len(v)==2 and str(v[0])=='DIMENSION_L_T_M'}
            require(bool(replacements),'metadata mutation sites');mutated=body.xreplace(replacements);rejected=False
            try:guard({tag:mutated},result,r,bound)
            except ValueError as e:rejected=e.args[0][0]=='restored quadrature unit'
            require(rejected,'native guard rejects altered cache unit')
            controls.append({'tag':tag,'unit':list(map(int,expected)),'changedUnitRejected':True})
    for suffix in ('0_0_0','2_1_0','0_1_999','0_1_-1','0_1'):
        rejected=False
        try:domain.cache_for_tag(result,'PY_S11CD_METADATA_'+domain.PREFIX+'_CACHE_NODESRESIDUAL_'+suffix)
        except ValueError:rejected=True
        require(rejected,'invalid cache address rejected');controls.append({'invalidAddress':suffix,'rejected':True})
    save(base/'metadata-controls.json',controls)
    require(before=={n:digest(base/n) for n in before} and not engine.PHYSICAL_METADATA.dimensions.constraints,'unchanged packet/dimension state')
    for n,a in original_artifacts.items():
        require(digest(origin/n)==a['sha256'],'original post hash')
        copied=base/('original-full.out' if n=='full.out' else n);require(digest(copied)==a['sha256'],'byte-identical copied operand')
    current=dict(preflight['sourceFiles']);current[str(Path(domain.__file__).resolve().relative_to(ROOT))]=digest(Path(domain.__file__))
    current[str(Path(__file__).resolve().relative_to(ROOT))]=digest(Path(__file__))
    for n,h in current.items():
        require(digest(ROOT/n)==h,'current repair source');target=base/'repair-source'/n;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/n,target)
    summary={'runDirectory':str(base),'sourceRunDirectory':str(origin),'status':'VALIDATED_INSTRUMENT_PREFLIGHT','sourceFiles':current,'originalSourceFiles':preflight['sourceFiles'],'provenance':preflight['provenance'],
        'repairJoin':join,'replayMethodJoin':method,'rows':80,'sources':70,'profiles':6,'records':6,'tagCount':len(entries),'writeKeyCount':len(keys),'metadataPaths':paths,
        'transformChecks':focused['transformChecks'],'transformArtifacts':focused['transformArtifacts'],'prefixArtifacts':focused['prefixArtifacts'],'coarseArtifacts':coarse,
        'recordArtifacts':result['recordArtifacts'],'workerArtifacts':result['workerArtifacts'],'workerManifest':json.loads((base/'workers.json').read_text()),'workerEvidence':workers,
        'metadataControls':len(controls),'originalRecoveryPayloadDifferences':len(differences),'profileLimitMutationRejected':True,'packetHashesBeforeEmission':before,'packetHashesAfterEmission':{n:digest(base/n) for n in before},
        'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in base.iterdir() if p.suffix in ('.pickle','.out')},
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'scope':'Saved-operand instrument preflight and cache metadata repair; no numerical quadrature repeated, no production position-domain result or physical tail/Abel/scattering claim.'}
    save(base/'checks.json',summary);save(base/'recovery.json',{'origin':str(origin),'repairJoin':join,'originalTranscriptSha256':digest(origin/'full.out'),'recoveryTranscriptSha256':digest(base/'full.out'),'originalRecoveryPayloadDifferences':len(differences),'numericalIntegrationRepeated':False})
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
