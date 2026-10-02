#!/usr/bin/env python3
"""Certify a saved control's exact finite domain; finish only untouched controls."""
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
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
HELPERS=('require','sha','save','replace_json','posthash_records','containment','Journal','decode','one_symbol')


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

def source_parts(source):
    tree=ast.parse(source)
    outer=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='scientific_work')
    loops=[(i,n) for i,n in enumerate(outer.body) if isinstance(n,ast.For)
           and isinstance(n.target,ast.Name) and n.target.id=='row'
           and ast.unparse(n.iter)=="('THETA_BALANCE', 'E_W_BALANCE')"]
    require(len(loops)==1,'unique unfinished row loop');position,loop=loops[0]
    theta=copy.deepcopy(loop);theta.iter=ast.Tuple(elts=[ast.Constant('THETA_BALANCE')],ctx=ast.Load())
    omitted=0;restored=0;body=[]
    for n in theta.body:
        if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='without' for t in n.targets):
            omitted+=1;continue
        if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value.func)=='J.nonzero' and ast.unparse(n.value.args[0])=="row + '-actual-lower-row-control'":
            omitted+=1;continue
        if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='good' for t in n.targets):
            n.value=ast.Name(id='saved_good',ctx=ast.Load());restored+=1
        body.append(n)
    require(omitted==2 and restored==1,'restore completed THETA baseline/omission only');theta.body=body
    ew=copy.deepcopy(loop);ew.iter=ast.Tuple(elts=[ast.Constant('E_W_BALANCE')],ctx=ast.Load())
    fn=ast.FunctionDef(name='unfinished_rows',args=ast.arguments(posonlyargs=[],args=[],kwonlyargs=[],kw_defaults=[],defaults=[]),
                      body=[theta,ew,*copy.deepcopy(outer.body[position+1:])],decorator_list=[])
    helpers=[n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in HELPERS]
    require({n.name for n in helpers}==set(HELPERS),'original helpers complete')
    return ast.fix_missing_locations(ast.Module(body=helpers,type_ignores=[])),ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[]))


def exact_finite_certificate(value,J,name):
    """Same finite predicate, proved from exact constant components if undecided."""
    direct=value.is_finite
    record={'value':value,'directFinite':direct,'accepted':direct is True,'route':'direct-flag'}
    if direct is None:
        record.update(route='exact-constant-fraction',freeSymbols=list(value.free_symbols),isNumber=value.is_number)
        if not value.free_symbols and value.is_number is True:
            numerator,denominator=sp.fraction(value)
            real,imaginary=denominator.as_real_imag()
            residual=sp.cancel(denominator-(real+sp.I*imaginary))
            finite=[numerator.is_finite,real.is_finite,imaginary.is_finite]
            nonzero=[real.is_zero,imaginary.is_zero]
            record.update(numerator=numerator,denominator=denominator,denominatorReal=real,denominatorImaginary=imaginary,
                          reconstructionResidual=residual,componentFinite=finite,componentZero=nonzero,
                          realComponents=[real.is_real,imaginary.is_real])
            record['accepted']=(residual==0 and all(x is True for x in finite)
                                and real.is_real is True and imaginary.is_real is True
                                and any(x is False for x in nonzero))
    J.emit(name+'-finite-certificate',record)
    return record


def certificate_decision(record):
    """Exact facts only; unknown/zero/nonfinite cannot be treated as nonzero."""
    if record['directFinite'] is True:return True
    if record['directFinite'] is not None:return False
    return (record.get('isNumber') is True and record.get('freeSymbols')==[]
            and record.get('reconstructionResidual')==0
            and all(x is True for x in record.get('componentFinite',[]))
            and len(record.get('componentFinite',[]))==3
            and record.get('realComponents')==[True,True]
            and any(x is False for x in record.get('componentZero',[])))


def make_journal(parent):
    class CertifiedJournal(parent):
        def nonzero(self,name,baseline,corrupt):
            previous=self.active;self.active=name
            self.emit(name+'-input',{'baseline':baseline,'corrupt':corrupt})
            movement=sp.simplify(corrupt-baseline)
            self.emit(name+'-return',{'movement':movement,'zero':movement.is_zero,'finite':movement.is_finite})
            certificate=exact_finite_certificate(movement,self,name)
            require(certificate['accepted']==certificate_decision(certificate),'finite certificate decision joins facts')
            require(movement.is_zero is False and certificate['accepted'] is True,name)
            self.active=previous
    return CertifiedJournal


def verify_gate(path,manifest_path,manifest):
    g=json.loads(Path(path).read_text())
    require(g['workerSha256']==sha(__file__) and g['manifestSha256']==sha(manifest_path),'worker/manifest identity')
    require(g['sourcePins']==manifest['sourcePins'],'same source pins')
    for p,h in g['sourcePins'].items():require(sha(p)==h,'source pin '+p)
    for k in ('sharedGuard','supervisor','authority','repairRecord','priorCompletion','priorGate','buildReviewRecord','methodRecord'):
        require(sha(g[k])==g[k+'Sha256'],'gate identity '+k)
    authority=json.loads(Path(g['authority']).read_text());repair=json.loads(Path(g['repairRecord']).read_text())
    prior=json.loads(Path(g['priorCompletion']).read_text());oldgate=json.loads(Path(g['priorGate']).read_text())
    review=json.loads(Path(g['buildReviewRecord']).read_text());method=json.loads(Path(g['methodRecord']).read_text())
    require(authority['continuedBoundedWorkAuthorized'] is True and authority['automaticScientificRetry'] is False,'standing continued-work authority')
    require(repair['toolingOnly'] is True and repair['testsPassed'] is True and repair['workerSha256']==g['workerSha256'],'tested exact predicate implementation')
    require(prior['allInspectionChecksPassed'] is True and prior['status']=='FAILED_INDETERMINATE_FINITE_FLAG_SAVED_CONTROLS_PRESERVED','actual preserved finite-flag failure')
    require(prior['workerSha256']==oldgate['workerSha256']==sha(manifest['priorWorker']),'prior worker ancestry')
    require(oldgate['sharedGuardSha256']==g['sharedGuardSha256'] and oldgate['supervisorSha256']==g['supervisorSha256'] and oldgate['buildReviewRecordSha256']==g['buildReviewRecordSha256'],'helper/review ancestry')
    require(review['reports']['claude']['literalVerdict']=='CLEAR FOR THIS BOUNDED RAW-INCREMENT BUILD' and review['reports']['grok']['literalVerdict']=='NEEDS REVISION','literal reports')
    require(g['independentBuildClearance'] is False and method['jointIndependentMethodClearance'] is True,'scoped method versus build status')
    require(g['scope']==manifest['scope'] and g['scientificRunsAuthorized']==1 and g['durationLimits'] is None,'one launch/no deadline')
    return g


def run_finish(manifest,J,helpers,parts):
    prior,files=copy_prior(manifest,J.out);original=prior/'prior'
    observations={}
    def load(name,field=None):
        check_file(prior/name,files[name]);value=json.loads((prior/name).read_text())
        observations.setdefault(name,{'sha256':files[name]['sha256'],'fields':[]})['fields'].append(field)
        return helpers['decode'](value if field is None else value[field])
    context=load('restored-context.json');symbols=context['symbols']
    require(context['completedFunctionsCalled'] is False,'saved context provenance')
    oldops=json.loads((prior/'operation-index.json').read_text())
    require(len(oldops)==4,'four original nested returns')
    for op in oldops:
        require(op['status']=='RESTORED_PRIOR_COMPLETE_RETURN' and op['functionCalled'] is False,'prior restoration status')
        for key in ('input','result'):check_file(prior/op[key]['path'],op[key])
        value=load(op['result']['path'])
        J.emit('restore-'+op['name'],{'value':value,'originalInput':op['input'],'originalResult':op['result'],'functionCalled':False})
        J.completed.append({**op,'input':{**op['input'],'path':'prior/'+op['input']['path']},'result':{**op['result'],'path':'prior/'+op['result']['path']}})
        helpers['replace_json'](J.out/'operation-index.json',J.completed)
    lower=load('prior/lower-boundary-return.json');mixed=load('prior/saved-operands/bare.json','mixed');modes=load('prior/saved-operands/modes.json');binding=load('prior/saved-operands/binding.json');point=load('prior/saved-operands/point.json')
    original_manifest=json.loads(Path(manifest['originalManifest']).read_text());physical=json.loads(Path(original_manifest['physicalInput']).read_text())
    require(binding['physicalInput']==physical and point['originalParameters']==physical['parameters'],'original physical input')
    require(context['frequency']==point['omega']==binding['frequency'] and context['edge']==point['edge'],'frequency/edge context')
    for key,name in [('qi','q_i'),('qh','q_h'),('qs','q_s'),('qo','q_o'),('k','k'),('H','H'),('omega','omega'),('rho','rho_m')]:
        require(symbols[key]==helpers['one_symbol']([mixed],name),'actual source symbol '+name)
    require([symbols['eta'],symbols['sigma']]==binding['independentGrades'] and symbols['epsilon']==binding['epsilon'],'independent grade/epsilon ancestry')
    increments={};consumer_bases={}
    replacements=dict(context['replacements']);factors=context['factors'];qpoint=dict(context['qpoint'])
    for row in original_manifest['rowNames']:
        increments[row]=load('prior/'+row+'-retained-increment.json','mixedCoefficient')
        require(load('prior/'+row+'-retained-increment.json','grades')==[1,1],'grade '+row)
        consumer_bases[row]=load('prior/'+row+'-full-increment-input.json','completePublishedBoundPressureSlotSum')
        slots=load('prior/'+row+'-full-increment-input.json','slotInsertion')
        require({str(k):v for k,v in replacements.items()}==slots,'actual restored slots '+row)
    for side in ('plus','minus'):
        require(factors[side]==[load('prior/'+side+'-closed-before-cancel.json','reference'),load('prior/'+side+'-closed-before-cancel.json','jet')],'actual external factors '+side)
    savedpoint=load('nonzero-input-control-point.json')
    require(savedpoint['cs']==context['values']['c_s0'],'saved physical speed context')
    require(qpoint[symbols['qi']]==point['inputDepth'] and qpoint[symbols['qo']]==point['outputDepth'],'saved reference depth context')
    all_symbols={v.name:v for v in symbols.values()}
    require(len(all_symbols)==len(symbols),'unique context symbol names')
    require(set(savedpoint['sample'])<=set(all_symbols),'saved sample symbol coverage')
    sample={all_symbols[name]:value for name,value in savedpoint['sample'].items()}
    completed_checks=['control-height-dispersion','control-slope-dispersion','wrong-lab-face-response','native-slope-omission-response','wrong-sheet-response']
    artifact_index=json.loads((prior/'artifact-index.json').read_text())
    for name in completed_checks:
        operand=load(name+'-input.json');value=load(name+'-return.json')
        if name.startswith('control-'):require(value['cancelled']==0,'saved completed dispersion')
        else:require(value['zero'] is False and value['finite'] is True,'saved completed response')
        J.emit('reuse-'+name,{'input':operand,'return':value,'functionCalled':False})
        J.completed.append({'name':name,'status':'RESTORED_PRIOR_COMPLETE_CHECK','functionCalled':False,
            'input':{**artifact_index[name+'-input.json'],'path':'prior/'+name+'-input.json'},
            'result':{**artifact_index[name+'-return.json'],'path':'prior/'+name+'-return.json'}})
        helpers['replace_json'](J.out/'operation-index.json',J.completed)
    name='THETA_BALANCE-actual-lower-row-control';saved_input=load(name+'-input.json');saved_return=load(name+'-return.json')
    J.active=name+'-saved-finite-certificate'
    J.emit(name+'-saved-operand-return',{'input':saved_input,'return':saved_return,'movementRecomputed':False,'argumentRoute':artifact_index[name+'-input.json'],'returnRoute':artifact_index[name+'-return.json']})
    require(saved_return['zero'] is False and saved_return['finite'] is None,'actual incomplete finite predicate')
    movement=saved_return['movement'];certificate=exact_finite_certificate(movement,J,name)
    require(certificate['accepted']==certificate_decision(certificate),'saved certificate decision joins facts')
    require(certificate['accepted'] is True and movement.is_zero is False,'saved movement finite and nonzero')
    J.emit(name+'-certificate-accepted',{'zeroFlag':saved_return['zero'],'originalFiniteFlag':None,'certifiedFinite':True,'movementRecomputed':False})
    J.active='remaining-row-controls'
    env={'sp':sp,'J':J,'increments':increments,'consumer_bases':consumer_bases,'replacements':replacements,
         'factors':factors,'qpoint':qpoint,'sample':sample,'saved_good':saved_input['baseline'],'mixed':mixed,'lower':lower,
         'h':modes['h'],'s':modes['s'],'freq':context['frequency'],'edge':tuple(context['edge'])}
    mapping={'qi':'qi','qh':'qh','qs':'qs','qo':'qo','k':'k','H':'H','Q':'Q','cs':'cs','eta':'eta','sigma':'sigma','epsilon':'eps','Dplus':'Dplus','Dminus':'Dminus'}
    env.update({target:symbols[key] for key,target in mapping.items()})
    J.emit('unfinished-source-and-context',{'source':ast.unparse(parts),'observedFields':observations,'symbols':symbols,'sample':[[a,b] for a,b in sample.items()],
        'savedThetaBaseline':saved_input['baseline'],'reconstructedContextOnly':'Small unsaved THETA probe dictionary and the original changed_insertion callable; no completed good/without/response recomputation.'})
    exec(compile(parts,manifest['originalWorker']+'#remaining-rows','exec'),env)
    result=env['unfinished_rows']();result.update(restoredNestedReturns=4,restoredPassedChecks=5,savedMovementCertifiedWithoutRecompute=True,
        priorFilesCopied=len(files),completedFunctionsReplayed=False,scientificAcceptance=False)
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--gate',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();manifest=json.loads(args.inputs.read_text());gate=verify_gate(args.gate,args.inputs,manifest)
    pins={**manifest['sourcePins'],str(args.inputs):sha(args.inputs),str(args.gate):sha(args.gate)}
    args.out.resolve().relative_to(ROOT/'_scratch/s11c');args.out.mkdir(exist_ok=False)
    J=None;result={};code=1;start=time.monotonic()
    try:
        helper_ast,parts=source_parts(Path(manifest['originalWorker']).read_text())
        helpers={'ast':ast,'hashlib':hashlib,'json':json,'os':os,'Path':Path,'resource':resource,'THREADS':THREADS}
        exec(compile(helper_ast,manifest['originalWorker']+'#helpers','exec'),helpers)
        save(args.out/'containment.json',helpers['containment']())
        global sp
        import sympy as sp
        from sympy.core.symbol import Str
        helpers.update(sp=sp,Str=Str)
        J=make_journal(helpers['Journal'])(args.out)
        result=J.stage('remaining-controls',{'priorCompletionSha256':gate['priorCompletionSha256'],'priorInventorySha256':sha(manifest['priorInventory'])},lambda:run_finish(manifest,J,helpers,parts))
        result['executionStatus']='COMPLETED_SAVED_RAW_INCREMENT_CONTROL_FINISH';code=0
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
            except OSError as e:records[path]={'expected':expected,'actual':None,'error':str(e)}
        save(args.out/'posthashes.json',records)
        if any(x['expected']!=x['actual'] for x in records.values()):result['integrityFailure']=True;code=1
        if J is not None:helpers['replace_json'](args.out/'operation-index.json',J.completed)
        result.update(wallSeconds=time.monotonic()-start,scientificAcceptance=False)
        save(args.out/'checks.json',result);sys.stdout.write((args.out/'checks.json').read_text())
    return code

if __name__=='__main__':sys.exit(main())
