#!/usr/bin/env python3
"""Resume only unfinished near-unity point checks using exact domain certificates.

Inert at import. All scientific restoration and certification occur inside the
unchanged pooled guard and native containment. No producer/source/limit replay.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import pickle
import shutil
import sys
import traceback
import S11c_d_near_unity_uniform as original
from S11c_d_numerical_radiating_blob_store import BlobStore

ROOT=original.ROOT
STAGE='near_unity_uniform_continue'
require=original.require
digest=original.digest
save=original.save


def exact_structure(a,b):
    """Exact saved structure/array bytes; never printed-summary equivalence."""
    if type(a) is not type(b):return False
    if isinstance(a,dict):return a.keys()==b.keys() and all(exact_structure(a[k],b[k]) for k in a)
    if isinstance(a,(tuple,list)):return len(a)==len(b) and all(exact_structure(x,y) for x,y in zip(a,b))
    if hasattr(a,'dtype') and hasattr(a,'shape'):
        return a.dtype==b.dtype and a.shape==b.shape and a.tobytes()==b.tobytes()
    return bool(a==b)


def certificate_acceptance(finite,free_symbols,reconstruction_zero,real_is_real,imag_is_real,real_sign,imag_sign):
    return bool(finite and not free_symbols and reconstruction_zero and real_is_real is True
        and imag_is_real is True and (real_sign in (-1,1) or imag_sign in (-1,1)))


def exact_nonzero(value):
    """Exact real/imaginary witness; unknown is never promoted by a decimal."""
    finite=original.finite(value);flag=value.is_zero
    result=dict(value=value,finite=finite,originalZeroFlag=flag,freeSymbols=tuple(sorted(value.free_symbols,key=str)))
    if not finite or value.free_symbols:
        return result|dict(certified=False,route='NONFINITE_OR_UNBOUND')
    if flag is False:
        return result|dict(certified=True,route='ORIGINAL_EXACT_NONZERO_PROPERTY')
    real,imag=value.as_real_imag(deep=True)
    real=sp.simplify(real);imag=sp.simplify(imag)
    residual=sp.cancel(value-real-sp.I*imag)
    def sign(x):
        if x.is_positive is True:return 1
        if x.is_negative is True:return -1
        if x.is_zero is True:return 0
        return None
    rs,ims=sign(real),sign(imag)
    accepted=certificate_acceptance(finite,result['freeSymbols'],residual==0,
        real.is_real,imag.is_real,rs,ims)
    return result|dict(certified=accepted,route='EXACT_REAL_IMAGINARY_SIGN_WITNESS' if accepted else 'UNRESOLVED_EXACT_NONZERO',
        real=real,imaginary=imag,reconstructionResidual=residual,realIsReal=real.is_real,
        imaginaryIsReal=imag.is_real,realSign=rs,imaginarySign=ims)


class ContinuationJournal(original.Journal):
    def __init__(self,out,manifest):
        super().__init__(out);self.manifest=manifest;self.restored={};self.prior_inputs={};self.evidence_cache={}
        self.restored_count=0;self.certificates={};self.certificate_uses=0
        prior=manifest['priorComplete'];dest=out/'prior-complete';dest.mkdir()
        copies={}
        for relative,pin in prior['files'].items():
            source=Path(prior['path'])/relative;target=dest/relative;target.parent.mkdir(parents=True,exist_ok=True)
            require(digest(source)==pin['sha256'] and source.stat().st_size==pin['bytes'],'prior file integrity '+relative)
            shutil.copyfile(source,target)
            require(digest(target)==pin['sha256'],'byte-identical prior copy '+relative)
            copies[relative]=dict(source=str(source),copy=str(target),**pin)
        save(out/'prior-copy-index.json',copies)
        self.prior_store=BlobStore(dest/'operations.sqlite');self.prior_store.integrity_check()
        self.blob_refs={}
        for name,h,n in self.prior_store.connection.execute('select name,sha256,bytes from blobs'):
            ref=dict(storage='sqlite',member=name,sha256=h,bytes=n)
            self.prior_store.get(ref);self.blob_refs[name]=ref
        require(len(self.blob_refs)==5383,'complete prior blob inventory')
        events=[json.loads(line) for line in (dest/'operation-index.jsonl').read_text().splitlines()]
        self.terminals={e['name']:e for e in events if e['status'] in ('COMPLETE','UNRESOLVED','FAILED_CHECK')}
        require(len(self.terminals)==71 and sum(e['status']=='COMPLETE' for e in self.terminals.values())==27,
            'exact prior operation inventory')
        require(all(e['status']!='FAILED_CHECK' for e in self.terminals.values()),'no unidentified prior scientific failure')

    def reference(self,ref):return dict(ref,storage='prior-sqlite',database='prior-complete/operations.sqlite')

    def read_prior(self,ref):return original.decode(self.prior_store.get(ref))

    def evidence(self,name):
        if name not in self.evidence_cache:
            self.evidence_cache[name]=self.read_prior(self.blob_refs[name])
            self.event(dict(name=name,status='REUSED_PRIOR_COMPLETE_EVIDENCE',value=self.reference(self.blob_refs[name]),functionCalled=False))
        return self.evidence_cache[name]

    def restore_completed(self):
        for name,e in self.terminals.items():
            if e['status']!='COMPLETE':continue
            args=self.read_prior(e['input']);value=self.read_prior(e['result'])
            self.prior_inputs[name]=args;self.restored[name]=value;self.restored_count+=1
            self.event(dict(name=name,status='RESTORED_PRIOR_COMPLETE_RETURN',input=self.reference(e['input']),
                result=self.reference(e['result']),functionCalled=False,priorStatus='COMPLETE'))
        require(self.restored_count==27,'restore every completed return, including prior aggregate as history only')

    def join(self,name,actual):
        old=self.terminals[name];old_args=self.prior_inputs.get(name)
        if old_args is None:old_args=self.read_prior(old['input']);self.prior_inputs[name]=old_args
        raw=pickle.dumps(actual,protocol=5)
        byte_join=raw==self.prior_store.get(old['input'])
        structure_join=byte_join or exact_structure(old_args,actual)
        ref=self.emit('argument-joins/'+name,dict(prior=self.reference(old['input']),actual=actual,
            serializedBytesEqual=byte_join,exactStructuralEqual=structure_join))
        if not structure_join:
            save(self.out/'argument-mismatch.json',dict(operation=name,prior=self.reference(old['input']),current=ref,
                classification='INTEGRITY_FAILURE_NOT_PHYSICS'))
            raise original.IntegrityError('actual prior argument mismatch: '+name)
        return 'EXACT_SERIALIZED_BYTES' if byte_join else 'EXACT_STRUCTURAL_IDENTITY'


def certified_domain(N,mapping,J,prefix,grazing=False):
    require(not grazing,'completed grazing points must be restored, not recalculated')
    old=J.evidence(prefix+'/raw-source-domain');rows=old['factors'];bad=[];certificates={}
    require(rows.keys()==N['domains'].keys(),'saved domain object names')
    for name,factors in N['domains'].items():
        require(len(rows[name])==len(factors),'saved denominator inventory length '+name)
        certificates[name]=[]
        for i,(factor,row) in enumerate(zip(factors,rows[name])):
            require(exact_structure(factor,row['factor']),'saved actual denominator factor '+name)
            value=row['value'];J.certificate_uses+=1
            if value not in J.certificates:
                cert=exact_nonzero(value)
                ref=J.emit('nonzero-certificates/%05d'%len(J.certificates),cert)
                J.certificates[value]=dict(certificate=cert,receipt=ref)
            saved=J.certificates[value];cert=saved['certificate']
            certificates[name].append(dict(index=i,factor=factor,value=value,certificate=saved['receipt'],certified=cert['certified']))
            if not cert['certified']:bad.append(dict(object=name,index=i,factor=factor,certificate=cert))
    J.emit(prefix+'/certified-source-domain',dict(originalDomain=J.reference(J.blob_refs[prefix+'/raw-source-domain']),
        originalExclusions=old['excludedAtDirectSubstitution'],mapping=mapping,factors=certificates,unresolved=bad,
        argumentJoin='Actual pending schedule/sign/source/limit/probe arguments joined before this call; unchanged physical substitution recipe.',
        minimumNonzeroMagnitude=old['minimumNonzeroMagnitude'],magnitudeScope='Prior certified subset only; not a new global margin'))
    if bad:raise original.ScientificIssue('exact nongrazing denominator certificate unresolved',bad)
    return dict(minimumNonzeroMagnitude=old['minimumNonzeroMagnitude'],directSubstitutionExclusions=[],
        certificateRoute='SAVED_DENOMINATOR_EXACT_NONZERO',priorExclusions=old['excludedAtDirectSubstitution'],
        magnitudeScope='Prior certified subset only')


def certificate_self_checks(J):
    sample=sp.Integer(2250)+sp.sqrt(266)*(-1365+900*sp.I)/4+7500*sp.I
    examples=[('saved-shape',sample,True),('positive',sp.Integer(2),True),('negative',sp.Integer(-2),True),
        ('pure-imaginary',3*sp.I,True),('zero',sp.Integer(0),False),('infinity',sp.oo,False),
        ('unbound',sp.Symbol('certificate_unknown'),False)]
    data=[dict(name=name,expected=expected,certificate=exact_nonzero(value)) for name,value,expected in examples]
    J.check('certificate-self-checks',all(d['certificate']['certified']==d['expected'] for d in data),data)
    J.report('certificate-self-checks',data)


def continue_science(manifest,J):
    J.restore_completed();certificate_self_checks(J)
    restored={name:J.restored['restore/'+name] for name in original.PACKETS}
    for name in original.PACKETS:
        J.join('restore/'+name,dict(pin=manifest['packets'][name],originalBytes=J.blob_refs['original-inputs/'+name+'.pickle']))
    params=dict(J.prior_inputs['selected-lift']['parameters'])
    require(params['omega']==3 and params['s11cdTangentialMomentum1']==sp.Rational(1,5)
        and params['s11cdTangentialMomentum2']==sp.Rational(1,10),'unchanged fixed physical point')
    w,cs=sp.symbols('uniformFrequency uniformSoundSpeed',positive=True)
    k,km=sp.symbols('uniformNormal uniformMinusNormal',real=True)
    q,qm,Q=sp.symbols('uniformPhysicalDepth uniformMinusPhysicalDepth uniformAlgebraicRadical')
    lift=J.restored['selected-lift'];units=J.evidence('physical-units')
    require(exact_structure(units['field'],restored['uniformResponse']['fieldUnits'])
        and exact_structure(units['current'],restored['uniformResponse']['currentUnit']),'saved unit joins')
    engine=next(path for path in manifest['sourcePins'] if path.endswith('/S11c_d_mixing_scattering_sympy_audit.py'))
    J.join('selected-lift',dict(uniform=restored['uniformSource']['curl'],common=restored['uniformCommon'],
        parameters=params,source=manifest['sourcePins'][engine]))
    lift=dict(lift,fieldUnits=units['field'],currentUnit=units['current'])
    J.evidence('profiles/constant-end-and-zero-jets')
    ends={}
    for label in ('LEFT','RIGHT'):
        B=original.EndBinding(label,restored,params,w,cs,k,km,q,qm,Q,J);B.fixed_frequency=sp.Integer(3)
        R=J.restored[label+'/restriction'];N=J.restored[label+'/native-reconstruction']
        J.join(label+'/restriction',dict(actualSource=B.r['CLOSED_PENCIL_LEGS'],frequency=B.frequency,lift=lift,parameters=params))
        J.join(label+'/applicability',dict(source=B.uniform['records'][label],native=B.native['profileBindings']))
        J.join(label+'/native-reconstruction',dict(pairing=B.packet,native=B.native,parameters=params,frequency=3,cs=cs))
        J.join(label+'/native-source-controls',dict(nativeFaces=B.r['FACE_LEG_OBJECTS']))
        source_args=J.prior_inputs[label+'/native-reconstruction']
        require(exact_structure(source_args['pairing'],B.packet) and exact_structure(source_args['native'],B.native)
            and exact_structure(source_args['parameters'],params) and source_args['frequency']==3 and source_args['cs']==cs,
            'restored source/binding context identity '+label)
        require(exact_structure(R['lift'],lift['lift']),'selected lift source identity '+label)
        ends[label]=dict(B=B,R=R,N=N,omissions=J.restored[label+'/native-source-controls'],
            applicability=J.restored[label+'/applicability'],
            limits={sign:J.restored[label+'/exact-limits/'+str(sign)] for sign in (-1,1)})
        for sign in (-1,1):
            J.join(label+'/exact-limits/'+str(sign),dict(restriction=R,native=N,sign=sign))
    saved_schedule=J.evidence('exact-speed-schedule');schedule=saved_schedule['schedule']
    probes=saved_schedule['declaredSheetProbes'];require(len(schedule)==12 and probes['normalSign']==1,'saved declared schedule')
    J.emit('recovered-context',dict(parameters=params,symbols=(w,cs,k,km,q,qm,Q),schedule=saved_schedule,
        reconstruction='Only EndBinding attribute context and exact symbols; no bind/source reconstruction, derivative, limit or producer invoked.',
        physicalSubstitution='Unchanged small fixed-point algebraic substitution is reconstituted because it was not separately returned before the prior domain refusal. Saved denominator values are reused, not rebound.'))
    original.domain_at_point=certified_domain
    rows=[];probed=set()
    for entry in schedule:
        for label,item in ends.items():
            results=[]
            for sign in (-1,1):
                prefix='points/%02d/%s/%s'%(entry['index'],label,sign)
                limit=item['limits'][sign];probe=entry['index']==probes['indices'][label] and sign==probes['normalSign']
                args=dict(schedule=entry,sign=sign,limit=limit,actualSource=item['N'],probeSheet=probe)
                J.join(prefix,args)
                if J.terminals[prefix]['status']=='COMPLETE':
                    result=J.restored[prefix]
                else:
                    require(J.terminals[prefix]['status']=='UNRESOLVED','only declared unresolved points resumed')
                    result=J.call(prefix,args,lambda:original.evaluate_sign(item['B'],item['R'],item['N'],entry['speed'],sign,
                        limit,item['omissions'],J,prefix,probe),soft=True)
                if probe and result.get('status')=='SUPPORTED_ON_SELECTED_UNIFORM_SUBSPACE':probed.add(label)
                results.append(result)
            statuses=[r['status'] for r in results]
            status=('FAILED_CHECK' if 'FAILED_CHECK' in statuses else 'SUPPORTED_ON_SELECTED_UNIFORM_SUBSPACE'
                if statuses==['SUPPORTED_ON_SELECTED_UNIFORM_SUBSPACE']*2 else 'UNRESOLVED')
            row=dict(schedule=entry,end=label,status=status,momentumSigns=results,
                applicability=item['applicability']['applicability'],dependencySnapshot=manifest['dependencySnapshot'])
            rows.append(row);J.report('point-%02d-%s'%(entry['index'],label),row)
    statuses=[r['status'] for r in rows]
    overall=('FAILED_CHECK' if 'FAILED_CHECK' in statuses else 'SUPPORTED_ON_SELECTED_UNIFORM_SUBSPACE'
        if all(s=='SUPPORTED_ON_SELECTED_UNIFORM_SUBSPACE' for s in statuses) and probed=={'LEFT','RIGHT'} else 'UNRESOLVED')
    result=dict(J.restored['uniform-science'],status=overall,rows=rows,nativeSheetProbeEnds=sorted(probed))
    J.emit('continued-uniform-result',result)
    J.report('nonzero-certificate-summary',dict(uniqueValues=len(J.certificates),uses=J.certificate_uses,
        certified=sum(v['certificate']['certified'] for v in J.certificates.values()),
        rows=[dict(receipt=v['receipt'],certificate=v['certificate']) for v in J.certificates.values()]))
    return result


def gate_check(args):
    manifest=json.loads(args.inputs.read_text());gate=json.loads(args.gate.read_text())
    require(gate['status']=='READY_FOR_SAVED_UNIFORM_DOMAIN_CONTINUATION' and gate['pooledExecution'] is True,'fresh continuation gate')
    require(gate['workerSha256']==digest(__file__) and gate['manifestSha256']==digest(args.inputs),'worker/manifest pins')
    require(gate['sourcePins']==manifest['sourcePins'] and gate['outputDirectory']==str(args.out.resolve()),'gate scope/source joins')
    require(manifest['plan']['sha256']==original.PLAN_SHA and manifest['physicalInput']['sha256']==original.INPUT_SHA,'unchanged method/input')
    for path,h in manifest['sourcePins'].items():require(digest(path)==h,'changed source '+path)
    for pin in (*manifest['packets'].values(),manifest['plan'],manifest['physicalInput'],manifest['reviewRecord'],manifest['priorManifest'],manifest['authorityRecord']):original.verify(pin)
    require(gate['authorityRecord']==manifest['authorityRecord'],'explicit continuation authority pin')
    review=json.loads(Path(manifest['reviewRecord']['path']).read_text())
    require(review['status']=='INDEPENDENT_METHOD_CLEARANCE_BOTH_REVIEWERS_NO_RUNTIME_RESULT' and review['allChecksPassed'],'unchanged cleared method record')
    require(len(manifest['priorComplete']['files'])==32,'exact preserved prior output inventory')
    args.out.resolve().relative_to(ROOT/'_scratch/s11c')
    require(not args.out.exists(),'fresh one-run output; no automatic retry')
    return manifest,gate


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--inputs',type=Path,required=True);parser.add_argument('--gate',type=Path,required=True)
    args=parser.parse_args();manifest,gate=gate_check(args);gate_hash=digest(args.gate)
    limits=original.containment(manifest)
    global sp
    import sympy as sp
    import numpy as np
    original.sp=sp;original.np=np
    args.out.mkdir(parents=True,exist_ok=False);J=None;result=None;exit_code=2
    save(args.out/'invocation.json',dict(workerSha256=digest(__file__),manifestSha256=digest(args.inputs),gateSha256=gate_hash,
        containment=limits,restoration='All scientific imports and saved-object restoration follow containment.'))
    try:
        J=ContinuationJournal(args.out,manifest)
        result=J.call('uniform-domain-continuation',dict(manifest=manifest,gate=gate,containment=limits),lambda:continue_science(manifest,J),soft=True)
        exit_code=0 if result.get('status')=='SUPPORTED_ON_SELECTED_UNIFORM_SUBSPACE' else 2
    except BaseException as error:
        result=dict(status='FAILED_CHECK',reason=str(error),exceptionType=type(error).__name__,traceback=traceback.format_exc(),
            activeOperations=[] if J is None else list(J.active))
        save(args.out/'failure.json',original.readable(result))
        if J is not None:J.emit('fatal-failure',result)
    finally:
        post={path:digest(path)==h for path,h in manifest['sourcePins'].items()}
        post.update({pin['path']:digest(pin['path'])==pin['sha256'] for pin in manifest['packets'].values()})
        post.update({str(Path(manifest['priorComplete']['path'])/rel):digest(Path(manifest['priorComplete']['path'])/rel)==pin['sha256'] for rel,pin in manifest['priorComplete']['files'].items()})
        post[str(args.inputs.resolve())]=digest(args.inputs)==gate['manifestSha256'];post[str(args.gate.resolve())]=digest(args.gate)==gate_hash
        save(args.out/'posthashes.json',post)
        count=0
        if J is not None:
            J.emit('posthashes',post);J.store.integrity_check();count=J.store.count();J.store.close();J.prior_store.close()
        if not all(post.values()):result=dict(status='FAILED_CHECK',reason='source/input/prior posthash changed',computed=result);exit_code=2
        checks=dict(status=result['status'],result=result,posthashes=post,blobCount=count,
            completedOperations=0 if J is None else J.count,restoredCompleteOperations=0 if J is None else J.restored_count,
            scientificFailures=[] if J is None else J.failures,containment=limits,
            manifestSha256=digest(args.inputs),workerSha256=digest(__file__))
        rendered=json.dumps(original.readable(checks),indent=2,allow_nan=False)+'\n'
        with (args.out/'checks.json').open('x') as stream:stream.write(rendered);stream.flush();original.os.fsync(stream.fileno())
        sys.stdout.write(rendered);sys.stdout.flush()
    return exit_code


if __name__=='__main__':raise SystemExit(main())
