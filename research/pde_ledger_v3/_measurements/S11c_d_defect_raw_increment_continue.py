#!/usr/bin/env python3
"""Restore saved raw construction and execute only its unfinished controls.

No scientific imports or payload restoration before the pinned pooled guard.
The original source, failed tree and all review verdicts remain unchanged.
"""
import argparse
import ast
import copy
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import sys
import time
import traceback

ROOT=Path('/var/projects/toy_physics')
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS',
         'NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
HELPERS=('require','sha','save','replace_json','posthash_records','containment',
         'Journal','decode','one_symbol')


def require(value,message):
    if value is not True:raise ValueError(message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()


def save(path,value):
    with Path(path).open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())


def string_keys(mapping):
    """Only evidence labels change; reject a collision rather than discard data."""
    result={str(k):v for k,v in mapping.items()}
    require(len(result)==len(mapping),'distinct evidence key labels')
    return result


def source_parts(source):
    tree=ast.parse(source)
    outer=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='scientific_work')
    start=[i for i,n in enumerate(outer.body) if isinstance(n,ast.Assign)
           and any(isinstance(t,ast.Name) and t.id=='control_k' for t in n.targets)]
    require(len(start)==1,'unique unfinished control boundary')
    tail=copy.deepcopy(outer.body[start[0]:])
    changes=0
    for n in ast.walk(ast.Module(body=tail,type_ignores=[])):
        if isinstance(n,ast.Dict):
            for i,(k,v) in enumerate(zip(n.keys,n.values)):
                if isinstance(k,ast.Constant) and k.value=='depths' and isinstance(v,ast.Name) and v.id=='control_depths':
                    n.values[i]=ast.Call(func=ast.Name(id='string_keys',ctx=ast.Load()),args=[v],keywords=[]);changes+=1
    require(changes==1,'exactly one evidence serialization repair')
    fn=ast.FunctionDef(name='resumed_controls',args=ast.arguments(posonlyargs=[],args=[],kwonlyargs=[],kw_defaults=[],defaults=[]),body=tail,decorator_list=[])
    q=next(n for n in outer.body if isinstance(n,ast.FunctionDef) and n.name=='q')
    helpers=[n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in HELPERS]
    require({n.name for n in helpers}==set(HELPERS),'complete helper definitions')
    return ast.fix_missing_locations(ast.Module(body=helpers,type_ignores=[])),ast.fix_missing_locations(ast.Module(body=[copy.deepcopy(q),fn],type_ignores=[]))


def check_file(path,record):
    require(Path(path).is_file() and Path(path).stat().st_size==record['bytes'] and sha(path)==record['sha256'],'immutable file '+str(path))


def copy_prior(manifest,out):
    prior=Path(manifest['priorComplete']);inventory=json.loads(Path(manifest['priorInventory']).read_text())
    require(inventory['root']==str(prior),'prior inventory root')
    files=inventory['files']
    actual={str(p.relative_to(prior)) for p in prior.rglob('*') if p.is_file()}
    require(actual==set(files),'entire prior result tree')
    target=out/'prior';target.mkdir();receipts={}
    for name,record in files.items():
        source=prior/name;destination=target/name
        check_file(source,record);destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(source,destination);check_file(destination,record)
        receipts[name]={'source':str(source),'copy':str(destination.relative_to(out)),**record}
    save(out/'prior-copy-index.json',receipts)
    return target,files


def verify_gate(path,manifest_path,manifest):
    gate=json.loads(Path(path).read_text())
    require(gate['workerSha256']==sha(__file__) and gate['manifestSha256']==sha(manifest_path),'worker and manifest gate')
    require(gate['sourcePins']==manifest['sourcePins'],'gate pins identical')
    for p,h in gate['sourcePins'].items():require(sha(p)==h,'gate pin '+p)
    for key in ('sharedGuard','supervisor','authority','repairRecord','priorCompletion','priorGate','buildReviewRecord','methodRecord'):
        require(sha(gate[key])==gate[key+'Sha256'],'gate identity '+key)
    require(gate['independentBuildClearance'] is False and gate['localToolingContinuation'] is True,'literal build status')
    authority=json.loads(Path(gate['authority']).read_text())
    repair=json.loads(Path(gate['repairRecord']).read_text())
    prior=json.loads(Path(gate['priorCompletion']).read_text())
    oldgate=json.loads(Path(gate['priorGate']).read_text())
    review=json.loads(Path(gate['buildReviewRecord']).read_text())
    method=json.loads(Path(gate['methodRecord']).read_text())
    require(authority['userDirection']=='Continue' and authority['executionsAuthorized']==1 and authority['automaticRetry'] is False,'renewed continuation authority')
    require(repair['workerSha256']==gate['workerSha256'] and repair['toolingOnly'] is True and repair['testsPassed'] is True,'tested serialization continuation')
    require(prior['status']=='FAILED_FORMATTER_PARTIAL_CONSTRUCTION_PRESERVED' and prior['allInspectionChecksPassed'] is True,'actual preserved failure')
    require(prior['workerSha256']==oldgate['workerSha256']==sha(manifest['originalWorker']),'original worker join')
    require(oldgate['buildReviewRecordSha256']==gate['buildReviewRecordSha256'] and oldgate['guardSha256']==gate['sharedGuardSha256'] and oldgate['supervisorSha256']==gate['supervisorSha256'],'prior helpers and review join')
    require(review['reports']['claude']['literalVerdict']=='CLEAR FOR THIS BOUNDED RAW-INCREMENT BUILD' and review['reports']['grok']['literalVerdict']=='NEEDS REVISION','literal reports unchanged')
    require(method['jointIndependentMethodClearance'] is True,'scoped method record')
    require(gate['scope']==manifest['scope'] and gate['scientificRunsAuthorized']==1 and gate['durationLimits'] is None,'bounded no-deadline continuation')
    return gate


def run_controls(manifest,J,helpers,parts):
    prior,files=copy_prior(manifest,J.out)
    original_manifest=json.loads(Path(manifest['originalManifest']).read_text())
    artifacts=json.loads((prior/'artifact-index.json').read_text())
    operations=json.loads((prior/'operation-index.json').read_text())
    reused=json.loads((prior/'saved-return-reuse.json').read_text())
    require(len(operations)==4 and len(reused)==17,'prior completed-work census')
    for alias,record in reused.items():
        expected=original_manifest['savedOperands'][alias]
        require(record['source']==expected['path'] and record['sha256']==expected['sha256'] and record['functionCalled'] is False,'original saved route '+alias)
        require(sha(prior/record['copy'])==sha(record['source'])==record['sha256'],'original saved bytes '+alias)
    observed={}
    def load(name,field=None):
        path=prior/name;check_file(path,files[name]);raw=json.loads(path.read_text())
        observed.setdefault(name,{'sha256':files[name]['sha256'],'fields':[]})['fields'].append(field)
        return helpers['decode'](raw if field is None else raw[field])
    mixed=load('saved-operands/bare.json','mixed')
    modes=load('saved-operands/modes.json')
    point=load('saved-operands/point.json')
    binding=load('saved-operands/binding.json')
    physical=json.loads(Path(original_manifest['physicalInput']).read_text())
    require(binding['physicalInput']==physical and point['originalParameters']==physical['parameters'],'physical input ancestry')
    qi,qh,qs,qo=(helpers['one_symbol']([mixed],n) for n in ('q_i','q_h','q_s','q_o'))
    k,H,om,rho=(helpers['one_symbol']([mixed],n) for n in ('k','H','omega','rho_m'))
    factor_input=load('physical-factorization-input.json')
    Q=helpers['one_symbol'](factor_input['momenta'],'increment_difference')
    cs=load('physical-sheet.json','cs');freq=point['omega'];edge=tuple(point['edge'])
    require(freq==binding['frequency']==load('physical-sheet.json','frequency'),'frequency ancestry')
    require(edge==tuple(load('physical-sheet.json','edge')),'edge ancestry')
    values={'c_s0':sp.sympify(physical['parameters']['c_s0'])};mass=sp.sympify(physical['parameters']['rho_m'])
    eta,sigma=binding['independentGrades'];eps=binding['epsilon'];h,s=modes['h'],modes['s']
    qpoint={qi:point['inputDepth'],qo:point['outputDepth']}
    # Only tiny argument/context expressions are rebuilt, never completed functions.
    lower_input=load('lower-boundary-input.json')
    native_normal=load('lower-boundary-operands.json','normalSource')
    expected_inputs={'lower-boundary':{'upperSaved':mixed,'nativeNormal':native_normal,'face':-1},
       'physical-factorization':{'coefficient':mixed,'momenta':[k,k+H,k+Q-H,k+Q]},
       'lower-wrong-height-control':{'labHeightSign':1,'required':-1},
       'lower-slope-omission-control':{'slopeMultiplier':0}}
    restored={}
    for op in operations:
        name=op['name'];require(name in expected_inputs,'known completed stage')
        for key in ('input','result'):
            record=op[key];require(record==artifacts[record['path']],'original journal receipt')
            check_file(prior/record['path'],record)
        actual=load(op['input']['path'])
        J.emit('restore-'+name+'-argument-join',{'actual':actual,'expected':expected_inputs[name],'sourceReceipt':op['input'],'sourceResult':op['result']})
        require(actual==expected_inputs[name],'actual prior arguments '+name)
        restored[name]=load(op['result']['path'])
        receipt={'name':name,'status':'RESTORED_PRIOR_COMPLETE_RETURN','functionCalled':False,
                 'input':{**op['input'],'path':'prior/'+op['input']['path']},
                 'result':{**op['result'],'path':'prior/'+op['result']['path']}}
        J.completed.append(receipt)
        helpers['replace_json'](J.out/'operation-index.json',J.completed)
    require(set(restored)==set(expected_inputs),'all four complete returns restored')
    lower=restored['lower-boundary'];wrongheight=restored['lower-wrong-height-control'];omitted=restored['lower-slope-omission-control']
    increments={};consumer_bases={};slot_dicts=[];factors={}
    for row in original_manifest['rowNames']:
        increments[row]=load(row+'-retained-increment.json','mixedCoefficient')
        require(load(row+'-retained-increment.json','grades')==[1,1] and load(row+'-retained-increment.json','epsilon')==eps,'saved grade/epsilon '+row)
        if row.startswith('U'):require(increments[row]==0,'saved zero U increment')
        consumer_bases[row]=load(row+'-full-increment-input.json','completePublishedBoundPressureSlotSum')
        require(load(row+'-full-increment-input.json','completeCensus') is True,'saved complete pressure census '+row)
        slot_dicts.append(load(row+'-full-increment-input.json','slotInsertion'))
    require(all(x==slot_dicts[0] for x in slot_dicts),'identical saved slot insertions across rows')
    replacements={helpers['one_symbol'](list(consumer_bases.values()),name):value for name,value in slot_dicts[0].items()}
    Dplus,Dminus=(helpers['one_symbol'](list(increments.values()),name) for name in ('increment_raw_plus','increment_raw_minus'))
    for label in ('plus','minus'):
        factors[label]=(load(label+'-closed-before-cancel.json','reference'),load(label+'-closed-before-cancel.json','jet'))
    J.emit('restored-context',{'symbols':{'qi':qi,'qh':qh,'qs':qs,'qo':qo,'k':k,'H':H,'Q':Q,'cs':cs,'omega':om,'rho':rho,'eta':eta,'sigma':sigma,'epsilon':eps,'Dplus':Dplus,'Dminus':Dminus},
        'frequency':freq,'edge':edge,'mass':mass,'values':values,'qpoint':[[a,b] for a,b in qpoint.items()],
        'replacements':[[a,b] for a,b in replacements.items()],'factors':factors,'observedSourceFields':observed,
        'reconstructedOnly':['q callable from exact original source','four small expected stage argument expressions','exact numeric cs/rho from pinned physical inputs','original unfinished control point declarations'],
        'completedFunctionsCalled':False})
    context={'sp':sp,'string_keys':string_keys,'J':J,'mixed':mixed,'lower':lower,'wrongheight':wrongheight,'omitted':omitted,
        'qi':qi,'qh':qh,'qs':qs,'qo':qo,'k':k,'H':H,'Q':Q,'cs':cs,'freq':freq,'edge':edge,'om':om,'rho':rho,'mass':mass,'values':values,
        'increments':increments,'consumer_bases':consumer_bases,'Dplus':Dplus,'Dminus':Dminus,'eps':eps,'eta':eta,'sigma':sigma,
        'replacements':replacements,'factors':factors,'qpoint':qpoint,'h':h,'s':s}
    J.emit('unfinished-control-source',{'originalWorker':manifest['originalWorker'],'originalSha256':sha(manifest['originalWorker']),
        'executedFragment':ast.unparse(parts),'serializationOnlyChange':"depths: control_depths -> depths: string_keys(control_depths)"})
    exec(compile(parts,manifest['originalWorker']+'#unfinished-controls','exec'),context)
    result=context['resumed_controls']()
    result.update(restoredCompleteStages=4,priorFilesCopied=len(files),priorFunctionsCalled=False,
                  continuationOnly=True,scientificAcceptance=False)
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--gate',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();manifest=json.loads(args.inputs.read_text());gate=verify_gate(args.gate,args.inputs,manifest)
    pins={**manifest['sourcePins'],str(args.inputs):sha(args.inputs),str(args.gate):sha(args.gate)}
    args.out.resolve().relative_to(ROOT/'_scratch/s11c');args.out.mkdir(exist_ok=False)
    J=None;result={};code=1;start=time.monotonic()
    try:
        source=Path(manifest['originalWorker']).read_text();helper_ast,parts=source_parts(source)
        helpers={'ast':ast,'hashlib':hashlib,'json':json,'os':os,'Path':Path,'resource':resource,'THREADS':THREADS}
        exec(compile(helper_ast,manifest['originalWorker']+'#helpers','exec'),helpers)
        save(args.out/'containment.json',helpers['containment']())
        global sp
        import sympy as sp
        from sympy.core.symbol import Str
        helpers.update(sp=sp,Str=Str)
        J=helpers['Journal'](args.out)
        result=J.stage('unfinished-controls',{'priorCompletionSha256':gate['priorCompletionSha256'],'priorInventorySha256':sha(manifest['priorInventory']),
            'originalWorkerSha256':sha(manifest['originalWorker'])},lambda:run_controls(manifest,J,helpers,parts))
        result['executionStatus']='COMPLETED_SAVED_RAW_INCREMENT_CONTROLS';code=0
    except BaseException:
        result={'executionStatus':'FAILED_PRESERVED','traceback':traceback.format_exc(),'incompleteOperation':None if J is None else J.active,'automaticRetry':False}
        save(args.out/'failure.json',result)
    finally:
        prior=args.out/'prior'
        if prior.exists():
            inventory=json.loads(Path(manifest['priorInventory']).read_text())['files']
            pins.update({str(prior/name):r['sha256'] for name,r in inventory.items() if (prior/name).exists()})
        records={}
        for path,expected in pins.items():
            try:records[path]={'expected':expected,'actual':sha(path),'error':None}
            except OSError as error:records[path]={'expected':expected,'actual':None,'error':str(error)}
        save(args.out/'posthashes.json',records)
        if any(x['expected']!=x['actual'] for x in records.values()):result['integrityFailure']=True;code=1
        if J is not None:
            helpers['replace_json'](args.out/'operation-index.json',J.completed)
        result.update(wallSeconds=time.monotonic()-start,scientificAcceptance=False)
        save(args.out/'checks.json',result);sys.stdout.write((args.out/'checks.json').read_text())
    return code


if __name__=='__main__':sys.exit(main())
