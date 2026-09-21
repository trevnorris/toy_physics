#!/usr/bin/env python3
"""Finish saved character inputs using their accepted amplitude certificates."""
import ast
import copy
import hashlib
import inspect
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import S11c_d_remaining_case_frequency_characters as h

f,native=h.f,h.native
PREVIOUS=f.STORE/'s11c-remaining-case-frequency-20260921/characters'
DIAGNOSTIC=PREVIOUS.parent/'characters-proof-diagnostic'
PLAN=f.M/'S11c_d_remaining_case_frequency_characters_recovery_plan.md'
REPAIR=f.M/'S11c_d_remaining_case_frequency_characters_proof_repair.json'
FACTOR_CP=f.M/'S11c_d_remaining_case_factors_checkpoint.json'
original_prepare,original_main=h.prepare,h.main
REUSED={};WRITES=[];PROOFS={};MEMO={};CURRENT_CASE=None


def adapter():
    original=ast.parse(inspect.getsource(original_prepare));changed=copy.deepcopy(original)
    targets=[]
    for parent in ast.walk(changed):
        for field,value in ast.iter_fields(parent):
            if isinstance(value,list):
                for i,node in enumerate(value):
                    if isinstance(node,ast.For) and ast.unparse(node.target)=='(_, _, v)' and ast.unparse(node.iter)=='uses':
                        targets.append((value,i))
    f.require(len(targets)==1,'one native saved amplitude assertion reader')
    sequence,index=targets[0];old=copy.deepcopy(sequence[index])
    f.require('AMPLITUDE_RECONSTRUCTION_RESIDUAL' in ast.unparse(old),'exact failed raw amplitude assertion')
    sequence[index]=ast.parse('check_source_proofs(base,label,si,physical,factors_packet)').body[0]
    baseline_guard=next(n for n in ast.walk(changed) if isinstance(n,ast.Call) and ast.unparse(n.func)=='f.require'
        and len(n.args)==2 and isinstance(n.args[1],ast.Constant) and n.args[1].value=='full baseline completed row/jet/character reuse')
    old_baseline=copy.deepcopy(baseline_guard)
    baseline_guard.func=ast.Name(id='reuse_baseline_rows',ctx=ast.Load())
    baseline_guard.args=[ast.Name(id=k,ctx=ast.Load()) for k in ('base','source_values','old')]
    mkdir=next(n for n in ast.walk(changed) if isinstance(n,ast.Call) and ast.unparse(n.func)=='target.mkdir')
    f.require([k.arg for k in mkdir.keywords]==['parents'],'exact precreated case directory call')
    mkdir.keywords.append(ast.keyword(arg='exist_ok',value=ast.Constant(True)))
    reverse=copy.deepcopy(changed);restored=0
    call=next(n for n in ast.walk(reverse) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='reuse_baseline_rows')
    call.func=old_baseline.func;call.args=old_baseline.args;call.keywords=old_baseline.keywords
    call=next(n for n in ast.walk(reverse) if isinstance(n,ast.Call) and ast.unparse(n.func)=='target.mkdir')
    call.keywords=[k for k in call.keywords if k.arg!='exist_ok']
    for parent in ast.walk(reverse):
        for field,value in ast.iter_fields(parent):
            if isinstance(value,list):
                for i,node in enumerate(value):
                    if isinstance(node,ast.Expr) and isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Name) and node.value.func.id=='check_source_proofs':
                        value[i]=copy.deepcopy(old);restored+=1
    f.require(restored==1 and ast.dump(reverse)==ast.dump(original),'whole original prepare reverse AST: proof reader, saved baseline guard and directory reuse')
    namespace=dict(vars(h),check_source_proofs=check_source_proofs,reuse_baseline_rows=reuse_baseline_rows)
    exec(compile(ast.fix_missing_locations(changed),str(Path(h.__file__))+':saved-amplitude-certificates','exec'),namespace)
    return namespace['prepare'],{'wholePrepareReverseAST':True,'proofReaderReplacements':1,'completedBaselineGuardReuses':1,'directoryKeywordEdits':1,
        'originalPrepareSha256':hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'originalMainSha256':h.source.inputs.body(original_main),'originalRowConstructorSha256':h.source.inputs.body(h.row_constructor),
        'nativeFactorPairsSha256':h.source.inputs.body(native.factors.proof_output.pairs),
        'nativeProofScalarsSha256':h.source.inputs.body(native.factors.proof_scalars)}


def safe_pickle(path,value):
    path=Path(path);name=str(path.relative_to(STATE['base']))
    if name in REUSED:
        f.require(path.is_symlink() and f.digest(path)==REUSED[name]['sha256'],'preserved completed input bytes')
        f.require(native.same(f.unpickle(path),value),'complete typed requested/saved input identity')
        WRITES.append({'path':name,'source':REUSED[name]['original'],'sha256':REUSED[name]['sha256'],'typedIdentity':True})
        f.save(STATE['base']/'completed-input-write-joins.json',WRITES)
        return
    f.require(not path.exists() and not path.is_symlink(),('new packet creates once',name))
    f.atomic_pickle(path,value)


def safe_save(path,value):
    path=Path(path);name=str(path.relative_to(STATE['base']))
    if name in REUSED:
        f.require(path.is_symlink() and f.digest(path)==REUSED[name]['sha256'],'preserved completed JSON bytes')
        f.require(json.loads(path.read_text())==json.loads(json.dumps(value)),'complete requested/saved JSON identity')
        WRITES.append({'path':name,'source':REUSED[name]['original'],'sha256':REUSED[name]['sha256'],'jsonIdentity':True})
        f.save(STATE['base']/'completed-input-write-joins.json',WRITES);return
    f.require(not path.is_symlink(),'no write through preserved reference')
    f.save(path,value)


def reuse_baseline_rows(base,source_values,old):
    checks=STATE['diagnosticChecks']
    paths={
        'accepted/accepted-bindings/'+native.BASELINE+'/case-binding.pickle':
            str(PREVIOUS/'complete/accepted/accepted-bindings'/native.BASELINE/'case-binding.pickle'),
        'accepted/accepted-frequency/baseline-binding.pickle':
            str(PREVIOUS/'complete/accepted/accepted-frequency/baseline-binding.pickle')}
    joined={}
    for name,original in paths.items():
        expected=checks['paths'][original]
        f.require(f.digest(base/name)==f.digest(Path(original))==expected,'actual unchanged baseline input producer')
        joined[name]={'originalDiagnosticInput':original,'currentReference':str(base/name),'sha256':expected}
    f.require(checks['baselineRows']==80 and checks['baselineRowsExact'] and checks['baselineJetsExact'],
              'completed full baseline row/jet equality evidence')
    f.require(native.same(source_values,old['sourceFrequencies']),'all actual saved baseline characters')
    f.save(base/'baseline-census-input-reuse.json',{'inputs':joined,'diagnosticChecksSha256':f.digest(DIAGNOSTIC/'checks.json'),
        'originalPrepareAST':STATE['adapterJoin']['originalPrepareSha256'],'baselineRows':80,'completedComparisonNotRepeated':True})


def proof_join(factor,source,record,owner,fi,memo=None):
    same=lambda a,b:native.same(a,b,memo=memo)
    f.require(same((factor['AMPLITUDE'],factor['FREQUENCY']),(source['symbolicAmplitude'],source['symbolicFrequency'])),
              'exact physical amplitude/frequency source join')
    f.require(factor['CHARACTER_EQUATION_RESIDUAL']==factor['CHARACTER_NORMALIZATION_RESIDUAL']==0,
              'unchanged native actual character equations')
    f.require(record['row']==owner['index'] and record['kind']=='AMPLITUDE_'+str(fi),'actual saved certificate owner and factor address')
    certificate=record['certificate']
    f.require(same(record['raw'],factor['AMPLITUDE_RECONSTRUCTION_RESIDUAL']) and same(certificate['LEFT'],factor['SOURCE'])
              and same(certificate['RIGHT'],factor['CHARACTER']*factor['AMPLITUDE']),'complete saved certificate/raw/source pair')
    f.require(certificate['RESIDUAL']==0 and all(v==0 for v in native.factors.proof_scalars(certificate)),
              'accepted actual amplitude reconstruction proof')


def check_source_proofs(base,label,si,physical,factors_packet):
    global CURRENT_CASE,MEMO
    if CURRENT_CASE!=label:MEMO={};CURRENT_CASE=label
    f.require(f.digest(base/'accepted/accepted-cases'/label/'factorization.pickle')==STATE['factorCases'][label],
              'whole actual case factor packet from accepted proof producer')
    owners={v['caseIndex']:v for v in factors_packet['addresses']};joined=[]
    directory=base/'character-proof-joins'/label;directory.mkdir(parents=True,exist_ok=True)
    for ri,fi,factor in physical['uses']:
        owner=owners[ri];key=(owner['kind'],owner['index'],fi)
        if owner['kind']=='accepted':
            location=base/'factor-proofs/accepted-factor-certificates.pickle'
            if 'baseline' not in PROOFS:PROOFS['baseline']=f.unpickle(location)
            record=next(v for v in PROOFS['baseline']['records'] if v['row']==owner['index'] and v['kind']=='AMPLITUDE_'+str(fi))
        else:
            f.require(owner['kind']=='new','actual recorded factor owner disposition')
            location=base/'factor-proofs/rows'/str(owner['index']).zfill(3)/('amplitude_'+str(fi)+'-certificate.pickle')
            if key not in PROOFS:PROOFS[key]=f.unpickle(location)
            record=PROOFS[key]
        saved={'case':label,'sourceIndex':si,'row':ri,'factorIndex':fi,'owner':owner,'factor':factor,
               'source':physical['source'],'certificateRecord':record,'certificatePath':str(location),
               'certificateSha256':f.digest(location),'sourceInput':str(base/'character-cases'/label/'source-inputs'/f'{si:03}.pickle')}
        path=directory/f'{si:03}-{ri:03}-{fi}.pickle';safe_pickle(path,saved)
        if label==native.BASELINE and si==0:
            # The diagnostic already checked these six exact producer pairs.
            old=next(v for v in STATE['diagnosticPairs'] if v['owner']==owner)
            f.require(native.same(old['factor'],factor) and native.same(old['proof'],record),'exact completed diagnostic certificate pair reuse')
            evidence=next(v for v in STATE['diagnosticChecks']['failedSource'] if v['row']==ri and v['factor']==fi)
            f.require(all(evidence[n] for n in ('amplitudeSourceSame','frequencySourceSame','characterEquationZero','normalizationZero',
                'certificateRawSame','certificateLeftSame','certificateRightSame','certificateResidualZero','proofScalarsZero')),
                'all completed source/certificate diagnostic guards')
            route='completed-diagnostic'
        else:
            proof_join(factor,physical['source'],record,owner,fi,MEMO);route='accepted-saved-certificate'
        joined.append({'row':ri,'factor':fi,'owner':owner,'certificatePath':str(location),'route':route,
                       'rawAmplitudeZero':factor['AMPLITUDE_RECONSTRUCTION_RESIDUAL']==0,'certifiedResidualZero':True})
    if si==0:
        ri,fi,factor=physical['uses'][0];owner=owners[ri]
        saved=f.unpickle(directory/f'{si:03}-{ri:03}-{fi}.pickle');record=saved['certificateRecord']
        changed_owner=dict(owner,index=owner['index']+10000)
        wrong_left=copy.deepcopy(record);wrong_left['certificate']['LEFT']+=1
        wrong_right=copy.deepcopy(record);wrong_right['certificate']['RIGHT']+=1
        wrong_proof=copy.deepcopy(record);wrong_proof['certificate']['RESIDUAL']=1
        operands={'owner':owner,'changedOwner':changed_owner,'record':record,'wrongLeft':wrong_left,'wrongRight':wrong_right,'wrongProof':wrong_proof}
        safe_pickle(directory/'mutation-operands.pickle',operands)
        rejects=h.source.inputs.source.rejects
        controls={'owner':rejects(lambda:proof_join(factor,physical['source'],record,changed_owner,fi)),
          'left':rejects(lambda:proof_join(factor,physical['source'],wrong_left,owner,fi)),
          'right':rejects(lambda:proof_join(factor,physical['source'],wrong_right,owner,fi)),
          'proof':rejects(lambda:proof_join(factor,physical['source'],wrong_proof,owner,fi))}
        f.save(directory/'mutation-controls.json',controls);f.require(all(controls.values()),'actual certificate owner/source/proof mutations reject')
    f.save(directory/(f'{si:03}-checks.json'),{'sourceIndex':si,'joins':joined})


def native_proof_joins(base,origin):
    joins=[]
    for path,cls,name in ((Path(native.factors.__file__),None,'proof_scalars'),
                          (Path(native.factors.proof_output.__file__),None,'pairs'),
                          (h.engine.HERE,'BoundedSourceFourierAssembly','reconstruction_certificate')):
        frozen=origin/'source'/path.relative_to(f.ROOT)
        def definition(p):
            nodes=ast.parse(p.read_text()).body
            if cls:nodes=next(n for n in nodes if isinstance(n,ast.ClassDef) and n.name==cls).body
            return next(n for n in nodes if isinstance(n,ast.FunctionDef) and n.name==name)
        left,right=definition(path),definition(frozen)
        f.require(ast.dump(left)==ast.dump(right),'unchanged actual native factor proof/pair reader body')
        joins.append({'current':str(path),'frozen':str(frozen),'class':cls,'function':name,
                      'currentSha256':f.digest(path),'frozenSha256':f.digest(frozen),
                      'bodySha256':hashlib.sha256(ast.dump(left).encode()).hexdigest()})
    f.save(base/'native-factor-proof-joins.json',joins)


def load(base):
    repair=json.loads(REPAIR.read_text());oldbase=PREVIOUS/'complete';old=json.loads((oldbase/'inputs.json').read_text())
    f.require(f.digest(Path(h.__file__))==repair['originalHelperSha256']==f.digest(PREVIOUS/'helper-source.py'),
              'original character helper/current/frozen source unchanged')
    f.require(f.digest(oldbase/'inputs.json')==repair['originalInputsSha256'],'original completed loader manifest')
    diagnostic=h.source.receipts.inspect_guard(DIAGNOSTIC,'diagnose')
    f.require((DIAGNOSTIC/'checks.json').read_bytes()==(DIAGNOSTIC/'diagnose.stdout').read_bytes()
              and f.digest(DIAGNOSTIC/'checks.json')==repair['diagnosticChecksSha256'],'clean saved-only diagnostic')
    manifest=dict(old,runDirectory=str(base),sourceFiles=dict(old['sourceFiles']),inputPackets=dict(old['inputPackets']),referencedInputs={})
    actual={str(p.relative_to(oldbase)) for p in oldbase.rglob('*') if p.is_file() and 'source' not in p.relative_to(oldbase).parts and p!=oldbase/'inputs.json'}
    f.require(actual==set(repair['completedFiles']),'all and only original completed files retained')
    for name,item in repair['completedFiles'].items():
        oldpath=oldbase/name
        f.require(f.digest(oldpath)==item['sha256'] and oldpath.stat().st_size==item['bytes'],'completed original bytes')
        f.require((str(oldpath.readlink()) if oldpath.is_symlink() else None)==item['rawLink'] and str(oldpath.resolve())==item['resolved'],
                  'completed original raw/resolved source address')
        h.source.reference(base,manifest,oldpath,name,item['sha256'])
        if not oldpath.is_symlink():REUSED[name]=manifest['referencedInputs'][name]
    h.source.reference(base,manifest,oldbase/'inputs.json','original-character-inputs.json',repair['originalInputsSha256'])
    for root,prefix,inventory in ((PREVIOUS,'original-character-logs',repair['originalLogs']),
                                 (DIAGNOSTIC,'accepted-proof-diagnostic',repair['diagnosticFiles'])):
        for name,item in inventory.items():h.source.reference(base,manifest,root/name,prefix+'/'+name,item['sha256'])
    cp=json.loads(FACTOR_CP.read_text());origin=Path(cp['runDirectory'])
    f.require(cp['status']=='PUBLISHED_ANNEX_VERIFIED' and f.digest(origin/'checks.json')==cp['checksSha256'],'accepted complete factor/certificate producer')
    h.source.reference(base,manifest,FACTOR_CP,'accepted-factor-proof-checkpoint.json',f.digest(FACTOR_CP))
    h.source.reference(base,manifest,origin/'checks.json','accepted-factor-proof-checks.json',cp['checksSha256'])
    factor_cases={}
    for label in repair['cases']:
        name='cases/'+label+'/factorization.pickle';factor_cases[label]=cp['artifacts'][name]['sha256']
        f.require(f.digest(base/'accepted/accepted-cases'/label/'factorization.pickle')==factor_cases[label],
                  'actual full source case packet equals accepted factor proof owner')
    for name,item in cp['artifacts'].items():
        if name=='accepted-factor-certificates.pickle' or (name.startswith('rows/') and '/amplitude_' in name and name.endswith('-certificate.pickle')):
            h.source.reference(base,manifest,origin/name,'factor-proofs/'+name,item['sha256'])
    for name,value in cp['sourceFiles'].items():
        f.require(f.digest(origin/'source'/name)==value,'accepted factor source snapshot')
        manifest['inputPackets'][str(origin/'source'/name)]=value
    native_proof_joins(base,origin)
    # Current native pair/proof reader definitions are already pinned in the
    # accepted remaining-case factor and source stages; retain their full files.
    for p in (Path(__file__).resolve(),PLAN,REPAIR,FACTOR_CP,Path(native.factors.proof_output.__file__)):
        name=str(p.relative_to(f.ROOT));value=f.digest(p)
        f.require(name not in manifest['sourceFiles'] or manifest['sourceFiles'][name]==value,'current proof helper joins')
        manifest['sourceFiles'][name]=value
    for name,value in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==value,'current source pin')
        if name in old['sourceFiles']:f.require(f.digest(oldbase/'source'/name)==value,'original source snapshot')
        dst=base/'source'/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/name,dst)
        f.require(f.digest(dst)==value,'frozen current proof helper')
    for name,value in manifest['inputPackets'].items():f.require(f.digest(Path(name))==value,'all input prehashes')
    STATE.update(base=base,factorCases=factor_cases,diagnosticChecks=json.loads((DIAGNOSTIC/'checks.json').read_text()),
                 diagnosticPairs=f.unpickle(base/'accepted-proof-diagnostic/failed-guard-source-certificate-pairs.pickle'))
    manifest['savedCharacterInputRecovery']={'original':str(PREVIOUS),'diagnostic':str(DIAGNOSTIC),
        'originalFiles':len(repair['completedFiles']),'savedPacketWriteIdentities':len(REUSED),
        'originalManifestSha256':repair['originalInputsSha256'],'diagnosticGuard':diagnostic,
        'scope':'Saved amplitude certificate reader only; original physical operands and native row census unchanged.'}
    f.save(base/'inputs.json',manifest);f.save(base/'completed-character-input-reuse.json',manifest['referencedInputs'])
    f.save(base/'character-proof-reader-join.json',STATE['adapterJoin'])
    return manifest,tuple(repair['cases'])


STATE={}
if __name__=='__main__':
    proxy=SimpleNamespace(**vars(f));proxy.atomic_pickle=safe_pickle;proxy.save=safe_save
    h.f=proxy
    fn,join=adapter();repair=json.loads(REPAIR.read_text())
    f.require(join==repair['staticJoins'],'actual proof reader adapter equals reviewed whole-prepare join')
    STATE['adapterJoin']=join
    h.load=load;h.prepare=fn
    original_main()
