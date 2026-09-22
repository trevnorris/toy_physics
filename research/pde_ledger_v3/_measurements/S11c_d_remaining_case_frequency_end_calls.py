#!/usr/bin/env python3
"""Join saved native end-call operands without evaluating a scientific call."""
import argparse
import ast
import copy
import gc
import hashlib
import inspect
import json
from pathlib import Path
import resource
import signal
import textwrap
import time

import S11c_d_remaining_case_frequency_end_inputs as h
import S11c_d_remaining_case_modes as modes_reader

f,sp,engine=h.f,h.sp,h.engine
CP=f.M/'S11c_d_remaining_case_frequency_end_inputs_checkpoint.json'
PLAN=f.M/'S11c_d_remaining_case_frequency_end_calls_plan.md'
SCOPE=('Saved actual analytic/mapping and full frequency-end call inputs. '
       'No analytic, binding, derivative, rational-table, domain, seed, mode, '
       'current, continuation, matrix or solve operation is performed. '
       'A saved-call candidate is not acceptance of a frequency operator or numerical reuse.')


def same(a,b):
    f.require(modes_reader.same(a,b),'exact complete typed end-call input')


def protect(base):
    h.inputs.protect_references(base)
    old_pickle,old_save=f.atomic_pickle,f.save
    def fresh(path):
        path=Path(path);path.relative_to(base)
        f.require(not path.exists() and not path.is_symlink(),('fresh saved call evidence',str(path)))
        return path
    f.atomic_pickle=lambda path,value:old_pickle(fresh(path),value)
    f.save=lambda path,value:old_save(fresh(path),value)


def body(value):
    text=textwrap.dedent(inspect.getsource(value));path=Path(inspect.getsourcefile(value))
    return {'path':str(path),'sha256':f.digest(path),'source':text,
            'wholeBodyAST':hashlib.sha256(ast.dump(ast.parse(text)).encode()).hexdigest()}


def native_joins():
    values={'analytic':engine.FullPencilModes.analytic,'mapping':engine.ChannelInput.mapping,
            'endSources':h.source.q.end_sources,'endTables':h.inputs.chart.end_tables,
            'pairing':engine.ClosedCurrentPairing.construct,
            'modalPrepare':engine.ModalCurrentSubspaces.prepare,
            'wholePair':h.end_native.Pair,'seeds':h.end_native.seeds,'maps':h.end_native.maps,
            'continuation':h.end_native.continue_pair,'exactArrayReader':modes_reader.same}
    joins={name:body(value) for name,value in values.items()}
    # Restrict the consumed-state claim to the complete actual native bodies.
    for name,expected in (('analytic',{'r','k','q'}),('mapping',{'parameters','limits'})):
        tree=ast.parse(joins[name]['source'])
        used={node.attr for node in ast.walk(tree) if isinstance(node,ast.Attribute)
              and isinstance(node.value,ast.Name) and node.value.id=='self'}
        f.require(used==expected,('whole native consumed-state fields',name,used))
    joins['interpretation']='Whole bodies are pinned and inspected only; none is called.'
    return joins


def load(base):
    cp=json.loads(CP.read_text());origin=Path(cp['runDirectory'])
    f.require(cp['status']=='ACCEPTED_CASE_FREQUENCY_END_INPUTS'
              and f.digest(origin/'checks.json')==cp['checksSha256'],'accepted complete end inputs')
    vr=Path(cp['validation']['runDirectory']);h.source.receipts.inspect_guard(vr,'validate')
    f.require(f.digest(vr/'checks.json')==cp['validation']['checksSha256']
              and (vr/'checks.json').read_bytes()==(vr/'validate.stdout').read_bytes(),'final independent end input acceptance')
    manifest={'runDirectory':str(base),'sourceFiles':dict(cp['sourceFiles']),
              'inputPackets':dict(cp['inputPackets']),'referencedInputs':{},
              'input':cp['input'],'settings':cp['settings'],'scope':SCOPE,
              'acceptedEndInputs':{'checkpoint':str(CP),'checkpointSha256':f.digest(CP),
                  'runDirectory':str(origin),'checksSha256':cp['checksSha256'],
                  'validatorDirectory':str(vr),'validatorChecksSha256':cp['validation']['checksSha256']}}
    for name,item in cp['artifacts'].items():
        h.source.reference(base,manifest,origin/name,name,item['sha256'])
    for name,item in cp['referencedInputs'].items():
        path=origin/name
        f.require(path.is_symlink() and str(path.readlink())==item['original']
                  and str(path.resolve())==item['resolvedOriginal'] and path.stat().st_size==item['bytes']
                  and f.digest(path)==item['sha256'],'unchanged accepted original reference')
    for path,name in ((CP,'accepted-end-input-checkpoint.json'),(origin/'checks.json','accepted-end-input-checks.json'),
                      (origin/'inputs.json','accepted-end-input-manifest.json'),(vr/'checks.json','accepted-end-input-validation.json')):
        h.source.reference(base,manifest,path,name,f.digest(path))
    # The native baseline frequency producer actually read this source. Join
    # its accepted producer address instead of assuming the case alias is it.
    path=f.M/'S11c_d_uniform_source_checkpoint.json';uc=json.loads(path.read_text());uo=Path(uc['runDirectory'])
    f.require(uc['status']=='PUBLISHED_ANNEX_VERIFIED' and f.digest(uo/'checks.json')==uc['checksSha256'],
              'accepted original uniform source')
    h.source.reference(base,manifest,path,'end-call-origins/uniform-source-checkpoint.json',f.digest(path))
    item=uc['artifacts']['uniform-source.pickle']
    h.source.reference(base,manifest,uo/'uniform-source.pickle','end-call-origins/uniform-source.pickle',item['sha256'])
    fc=json.loads((base/'end-origins/frequency-checkpoint.json').read_text())
    f.require(fc['inputPackets'][str(uo/'uniform-source.pickle')]==item['sha256'],
              'actual native frequency caller consumed this uniform source')
    for name,value in uc['checks']['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==f.digest(uo/'source'/name)==value,'native baseline uniform source pin')
        f.require(name not in manifest['sourceFiles'] or manifest['sourceFiles'][name]==value,'shared actual source pin')
        if name not in manifest['sourceFiles']:
            manifest['sourceFiles'][name]=value
            h.source.reference(base,manifest,uo/'source'/name,'source/'+name,value)
    # Retain frozen sources by explicit immutable references. In particular,
    # do not copy the 146 MB inherited generated source/checkpoint inventory.
    for name,value in cp['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==f.digest(origin/'source'/name)==value,'current and frozen source')
        h.source.reference(base,manifest,origin/'source'/name,'source/'+name,value)
    for path in (Path(__file__).resolve(),PLAN,Path(modes_reader.__file__)):
        name=str(path.relative_to(f.ROOT));value=f.digest(path)
        f.require(name not in manifest['sourceFiles'] or manifest['sourceFiles'][name]==value,'new source version')
        if name not in manifest['sourceFiles']:
            manifest['sourceFiles'][name]=value
            dst=base/'source'/name;dst.parent.mkdir(parents=True,exist_ok=True);dst.write_bytes(path.read_bytes())
    for name,value in manifest['inputPackets'].items():f.require(f.digest(Path(name))==value,('original input hash',name))
    f.save(base/'inputs.json',manifest)
    return manifest,tuple(cp['cases'])


def route(path):
    return {'path':str(path),'resolved':str(path.resolve()),'rawLink':str(path.readlink()) if path.is_symlink() else None,
            'sha256':f.digest(path),'bytes':path.stat().st_size}


def unit(known,atom):
    f.require(atom in known,('actual saved unit declaration',str(atom)))
    return known[atom]


def analytic_input(strong,omega,k,q,units):
    return {'quotient':strong,'omega':omega,'momentum':k,'radical':q,'units':units}


def prepare(base,labels):
    baseline=f.unpickle(base/'end-call-origins/uniform-source.pickle')
    atlas=[];summaries={};pending=[];end_families=[]
    # Actual saved returns are never reconstructed from a fixed-frequency
    # pencil. The native producer stores the entire original algebraic return.
    ctxpath=base/'end-input-cases'/h.BASELINE/'context-pairs.pickle'
    ctx=f.unpickle(ctxpath);omega=ctx['frame']['coordinates'][4];seed=ctx['seedInput']
    for end in ('LEFT','RIGHT'):
        rawpath=base/'end-input-cases'/h.BASELINE/end.lower()/'full-inputs.pickle';raw=f.unpickle(rawpath)
        pp=base/'end-accepted/frequency'/(end.lower()+'-frequency-pencil.pickle');p=f.unpickle(pp)
        same(raw['sourceRecord'],baseline['records'][end]);same(raw['signature']['curl'],baseline['curl'])
        same(raw['signature']['weakUnits'],baseline['units']['weak'])
        k,q=p['momentum'],p['radical'];known=raw['knownUnits']
        units={'omega':unit(known,omega),'momentum':unit(known,k),'radical':unit(known,q),
               'fieldUnits':ctx['frame']['fieldUnits'],'equationUnits':ctx['frame']['equationUnits']}
        call=analytic_input(raw['signature']['strong'],omega,k,q,units)
        mapping={'algebraic':p['originalAlgebraic'],'relation':p['originalRelation'],
                 'live':(omega,k,q,*seed['origin']),'parameters':seed['parameters'],'limits':seed['limits']}
        result={'analytic':(p['originalAlgebraic'],p['originalRelation'],p['branchResiduals']),
                'mapping':p['mapping'],'frequencyPencil':p}
        item={'owner':(h.BASELINE,end),'analyticInput':call,'mappingInput':mapping,'result':result,
              'sourceRoute':route(rawpath),'contextRoute':route(ctxpath),'resultRoute':route(pp)}
        target=base/'end-call-baseline'/end.lower();target.mkdir(parents=True)
        f.atomic_pickle(target/'input-value.pickle',item)
        same(p['origin'],seed['origin']);same(p['frequency'],ctx['analyticVariables']['frequency'])
        same(p['originalRelation'],raw['modeInput']['relation'])
        f.require(all(v==0 for v in p['branchResiduals']),'accepted baseline branch receipt')
        atlas.append(item)
    for label in labels:
        cp=base/'end-input-cases'/label/'context-pairs.pickle';context=f.unpickle(cp)
        state,frame=context['seedInput'],context['frame'];omega=frame['coordinates'][4]
        case_summary={}
        for end in ('LEFT','RIGHT'):
            target=base/'end-call-cases'/label/end.lower();target.mkdir(parents=True)
            rawpath=base/'end-input-cases'/label/end.lower()/'full-inputs.pickle';raw=f.unpickle(rawpath)
            sig=raw['signature'];k,q=atlas[0]['analyticInput']['momentum'],atlas[0]['analyticInput']['radical']
            known=raw['knownUnits'];units={'omega':unit(known,omega),'momentum':unit(known,k),'radical':unit(known,q),
                  'fieldUnits':frame['fieldUnits'],'equationUnits':frame['equationUnits']}
            call=analytic_input(sig['strong'],omega,k,q,units)
            owners=[i for i,item in enumerate(atlas) if modes_reader.same(call,item['analyticInput'])]
            mapping=None;mapping_owners=[]
            if owners:
                chosen=atlas[owners[0]];a,relation,branch=chosen['result']['analytic']
                mapping={'algebraic':a,'relation':relation,'live':(omega,k,q,*state['origin']),
                         'parameters':state['parameters'],'limits':state['limits']}
                mapping_owners=[i for i,item in enumerate(atlas) if modes_reader.same(mapping,item['mappingInput'])]
            full={'analyticInput':call,'mappingInput':mapping,'origin':state['origin'],
                  'frequency':context['analyticVariables']['frequency'],'seedInput':state,
                  'frame':frame,'context':context['context'],'sourceSignature':sig}
            matches=[i for i,item in enumerate(end_families) if modes_reader.same(full,item['input'])]
            family=matches[0] if matches else len(end_families)
            if not matches:end_families.append({'owner':(label,end),'input':full})
            packet={'address':(label,end),'input':full,'analyticCandidates':[atlas[i]['owner'] for i in owners],
                    'mappingCandidates':[atlas[i]['owner'] for i in mapping_owners],
                    'sourceRoute':route(rawpath),'contextRoute':route(cp),'family':family,
                    'familyOwner':end_families[family]['owner'],
                    'baselinePairs':[(item['owner'],item['analyticInput'],item['mappingInput']) for item in atlas]}
            f.atomic_pickle(target/'call-input-pairs.pickle',packet)
            same(state,frame['input']);same(sig['seedInput'],state);same(sig['frame'],frame)
            same(sig['context'],context['context']);same(raw['address'],(label,end))
            same(tuple(units['omega']),(0,-1,0));same(tuple(units['momentum']),(-1,0,0));same(tuple(units['radical']),(0,-1,0))
            if owners:
                same(raw['modeInput']['relation'],relation)
                same(raw['modeInput']['acoustic']['SOURCE_BRANCH_JOINS'],branch)
            # Preserve already-completed frequency-dependent two-leg operands
            # for a subsequent exact native caller/coordinate join. No inverse
            # substitution or derivative is used to manufacture an analytic return.
            pairing=raw['modeInput']['pairing']
            modalpath=base/'end-accepted/uniform/cases'/label/end.lower()/'modal.pickle'
            modal,modal_known=f.unpickle(modalpath)
            fields=('FREQUENCY_LEGS','NORMAL_LEGS','BULK_LEGS','CLOSED_PENCIL_LEGS','ACOUSTIC_WAVE_ROWS','SOURCE_BRANCH_JOINS')
            two_leg={key:pairing[key] for key in fields}
            symbolic=modal['SYMBOLIC_OPERANDS'];scalar=modal['SCALAR_OPERANDS']
            saved={'pairing':two_leg,'symbolic':{key:symbolic[key] for key in
                   ('PENCIL_PLUS','PENCIL_MINUS','NORMAL_PENCIL_PLUS','FREQUENCY_PENCIL_PLUS')},
                   'scalar':{key:scalar[key] for key in ('RADICAL_SCALE','NORMAL_TRANSPORT','FREQUENCY_TRANSPORT')},
                   'binding':raw['modeInput']['binding'],'knownUnits':modal_known,
                   'modalRoute':route(modalpath),'sourceRoute':route(rawpath),
                   'interpretation':'Original full frequency-leg values, not yet routed to a native end call.'}
            f.atomic_pickle(target/'saved-frequency-leg-operands.pickle',saved)
            same(saved['symbolic']['PENCIL_PLUS'],two_leg['CLOSED_PENCIL_LEGS'][0])
            same(saved['symbolic']['PENCIL_MINUS'],two_leg['CLOSED_PENCIL_LEGS'][1])
            # Mutate actual guard operands, retaining all unrelated operands.
            changed_coefficient=copy.deepcopy(call)
            matrix=sp.MutableDenseMatrix(call['quotient']);matrix[0,0]=matrix[0,0]+1
            changed_coefficient['quotient']=sp.ImmutableMatrix(matrix)
            changed_unit=copy.deepcopy(call);changed_unit['units']['omega']=(0,0,0)
            changed_parameter=copy.deepcopy(state);changed_parameter['parameters']['omega']=state['parameters']['omega']+1
            actual_address=(label,end);changed_address=(label,end+'_WRONG')
            operands={'actualCall':call,'changedCoefficient':changed_coefficient,'changedUnit':changed_unit,
                      'actualState':state,'changedParameter':changed_parameter,
                      'actualAddress':actual_address,'changedAddress':changed_address}
            f.atomic_pickle(target/'mutation-operands.pickle',operands)
            controls={'coefficient':not modes_reader.same(call,changed_coefficient),
                      'unit':not modes_reader.same(call,changed_unit),
                      'parameter':not modes_reader.same(state,changed_parameter),
                      'address':not modes_reader.same(actual_address,changed_address)}
            f.require(all(controls.values()),'actual saved call-input mutations respond')
            f.save(target/'mutation-controls.json',controls)
            if not owners or not mapping_owners:
                pending.append({'address':(label,end),'inputPath':str(target/'call-input-pairs.pickle'),
                                'missingSavedAnalyticReturn':not owners,'missingSavedMappingReturn':not mapping_owners,
                                'completedFrequencyLegOperands':str(target/'saved-frequency-leg-operands.pickle')})
            case_summary[end]={'analyticCandidateCount':len(owners),'mappingCandidateCount':len(mapping_owners),
                               'family':family,'familyOwner':end_families[family]['owner'],
                               'fullSavedFrequencyLegOperands':True,'actualMutations':len(controls),
                               'newScientificCalls':0,'frequencyEndReuseAccepted':False}
            f.save(target/'call-summary.json',case_summary[end]);del modal,modal_known;gc.collect()
        summaries[label]=case_summary
        f.save(base/'end-call-cases'/label/'case-summary.json',case_summary)
    f.atomic_pickle(base/'pending-end-call-inputs.pickle',pending)
    f.atomic_pickle(base/'end-call-input-families.pickle',end_families)
    return {'cases':summaries,'physicalEnds':sum(len(v) for v in summaries.values()),
            'callInputFamilies':len(end_families),'pendingSavedReturnUses':len(pending),
            'analyticCandidateUses':sum(v['analyticCandidateCount']>0 for c in summaries.values() for v in c.values()),
            'mappingCandidateUses':sum(v['mappingCandidateCount']>0 for c in summaries.values() for v in c.values()),
            'newScientificCalls':0,'frequencyEndReuseAccepted':False}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);args=ap.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic()
    protect(base);manifest,labels=load(base)
    joins=native_joins();f.save(base/'native-end-call-joins.json',joins);h.prohibit()
    result=prepare(base,labels)
    for name,value in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==f.digest(base/'source'/name)==value,'current/frozen end call source posthash')
    for name,value in manifest['inputPackets'].items():f.require(f.digest(Path(name))==value,('original input posthash',name))
    for name,item in manifest['referencedInputs'].items():
        path=base/name
        f.require(path.is_symlink() and str(path.readlink())==item['original']
                  and str(path.resolve())==item['resolvedOriginal'] and path.stat().st_size==item['bytes']
                  and f.digest(path)==item['sha256'],'full end-call reference postidentity')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*')
               if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    checks={**manifest,**result,'status':'COMPLETED_SAVED_END_CALL_INPUTS','nativeJoins':joins,'artifacts':artifacts,
            'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
