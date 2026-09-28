#!/usr/bin/env python3
"""Corrected saved-output validation; restore all 217 prior complete returns."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import pickle
import resource
import signal
import time
import traceback

STORE=Path('/var/projects/toy_physics/_scratch/s11c')
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS',
    'NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')


def require(value,message):
    if not value:raise ValueError(message)


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')
        stream.flush();os.fsync(stream.fileno())


def route(path):return {'path':str(path),'bytes':path.stat().st_size,'sha256':sha(path)}


def containment():
    group=next(v[3:] for v in Path('/proc/self/cgroup').read_text().splitlines() if v.startswith('0::'))
    root=Path('/sys/fs/cgroup')/group.lstrip('/')
    actual={key:(root/key).read_text().strip() for key in ('memory.max','memory.swap.max','pids.max')}
    actual.update(nice=os.getpriority(os.PRIO_PROCESS,0),affinity=sorted(os.sched_getaffinity(0)),
        threads={key:os.environ.get(key) for key in THREADS})
    require(actual['memory.max']==str(2*1024**3) and actual['memory.swap.max']=='0'
        and actual['pids.max']=='32' and actual['nice']>=15 and len(actual['affinity'])==1
        and all(v=='1' for v in actual['threads'].values()),'required validation containment missing')
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    def timeout(*_):raise TimeoutError('840-second saved-output validator limit')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(840)
    actual['nativeWallSeconds']=840
    return actual


class SavedCodec(pickle.Unpickler):
    def find_class(self,module,name):
        if not (module.startswith(('sympy.','numpy.')) or module in ('sympy','numpy','builtins','collections')):
            raise pickle.UnpicklingError((module,name))
        return super().find_class(module,name)


class Journal:
    def __init__(self,out,reuse):
        self.out,self.records=out,[]
        self.reuse=list(reuse);self.reused=0

    def value(self,name,value):
        p=self.out/name;p.parent.mkdir(parents=True,exist_ok=True)
        with p.open('xb') as f:
            pickle.dump(value,f,protocol=4);f.flush();os.fsync(f.fileno())
        return route(p)

    def op(self,name,function,*args):
        folder='operations/%03d-%s'%(len(self.records),name)
        r={'name':name,'input':self.value(folder+'/input.pickle',args),
           'startedUtc':datetime.now(timezone.utc).isoformat()}
        save(self.out/folder/'started.json',r)
        if self.reused<len(self.reuse):
            previous=self.reuse[self.reused]
            require(previous['name']==name,'ordered saved-return continuation boundary')
            path=Path(previous['value']['path'])
            require(route(path)==previous['value'],'prior saved return hash changed')
            with path.open('rb') as f:
                result=SavedCodec(f).load()
            r.update(execution='RESTORED_PRIOR_VALIDATION_RETURN',priorOperation=previous)
            self.reused+=1
        else:
            r['execution']='NEW_UNFINISHED_OPERATION'
            result=function(*args)
        r.update(value=self.value(folder+'/value.pickle',result),finishedUtc=datetime.now(timezone.utc).isoformat())
        save(self.out/folder/'completed.json',r);self.records.append(r)
        return result


def validate(source,manifest,out):
    import sympy as sp
    import numpy as np
    started=time.monotonic();pins=json.loads(manifest.read_text());journal=Journal(out,pins['resumeOperations'])
    require(pins['sourceRoot']==str(source) and pins['validatorSha256']==sha(Path(__file__)),'validator source/root pins')
    hashes={str(p.relative_to(source)):sha(p) for p in sorted(source.rglob('*')) if p.is_file()}
    save(out/'input-manifest.json',pins);save(out/'prehashes.json',hashes)
    require(hashes==pins['files'],'all saved result bytes must match inventory')
    resume_root=Path(pins['resumeRoot']).resolve();resume_root.relative_to(STORE)
    resumed_hashes={str(p.relative_to(resume_root)):sha(p) for p in sorted(resume_root.rglob('*')) if p.is_file()}
    save(out/'resume-prehashes.json',resumed_hashes)
    require(resumed_hashes==pins['resumeFiles'],'prior validation artifacts changed')
    operations=json.loads((source/'operation-index.json').read_text())
    named={r['name']:r for r in operations};cache={}
    def read_saved(name):
        require(name in hashes,'unmanifested saved artifact')
        path=source/name;path.resolve().relative_to(source)
        with path.open('rb') as stream:return SavedCodec(stream).load()
    def load(name):
        if name not in cache:cache[name]=journal.op('load-'+str(len(cache)),read_saved,name)
        return cache[name]
    def operation(name,part='value'):
        return load(str(Path(named[name][part]['path']).relative_to(source)))
    candidate=load('outgoing-prescription-candidate.pickle')
    domain=load('prescription-domain.pickle')
    context=operation('restore-kernel-context');units=operation('restore-kernel-units')
    binding=load('reference-physical-binding.pickle')
    branch=load('branch-and-coordinate-evidence.pickle')
    field_units,row_units,saved_units=load('row-field-unit-join.pickle')
    kn=context['normalMomentum'];fourier=context['fourierMass'];phase=context['normalCharacter']
    fixed_source=operation('bind-accepted-symbol');inverse=operation('bind-saved-inverse')
    determinant=operation('bind-saved-determinant');adjugate=operation('bind-saved-cofactors')
    density=operation('bind-saved-spectral-density')
    source_join=operation('actual-symbol-join');inverse_join=operation('inverse-determinant-adjugate-join')
    density_join=operation('spectral-density-phase-join')
    inv_args=operation('inverse-determinant-adjugate-join','input')
    density_args=operation('spectral-density-phase-join','input')
    structural={
        'candidateSourceJoin':candidate['fixedSymbol']==fixed_source,
        'candidateInverseJoin':candidate['fixedInverse']==inverse,
        'candidateDensityJoin':candidate['regularDensity']==density,
        'fullFiveByFiveShapes':all(a.shape==(5,5) for a in (fixed_source,inverse,adjugate,density)),
        'savedSourceResidualZero':source_join==sp.zeros(5),
        'savedInverseAdjugateResidualZero':inverse_join==sp.zeros(5),
        'inverseAdjugateOperands':inv_args==(inverse,determinant,adjugate),
        'savedDensityResidualZero':density_join==sp.zeros(5),
        'densityResidualOperands':density_args==(density,inverse,phase,fourier),
        'physicalInputPreserved':candidate['fixedPhysicalInput']==binding['suppliedPhysicalInput'],
        'referenceBindingPreserved':candidate['referenceBinding']==binding['actualSourceBinding'],
        'actualReferenceGradesZero':candidate['effectiveReferenceGrades']==binding['effectiveReferenceGrades']=={'eta_bg':0,'sigma_W':0},
        'bulkWaveRelationJoin':candidate['bulkWaveRelation']==operation('bind-physical-wave'),
        'evanescentBranchJoin':candidate['physicalBulkBranch']==sp.I*sp.sqrt(-operation('physical-q-square')),
        'strictEvanescentRealAxis':branch['degree']==2 and branch['linear']==0 and branch['constant']>0 and branch['quadratic']>0,
        'normalBulkMomentumUnits':tuple(branch['normalUnits'])==tuple(branch['physicalBulkUnits'])==tuple(units['spectralMeasureUnit']),
        'rowFieldUnitsPreserved':saved_units==units and field_units==list(map(tuple,units['fieldUnits'])) and row_units==list(map(tuple,units['rowUnits'])),
        'kernelAndResidueUnits':candidate['kernelEntryUnits']==candidate['residueEntryUnits']==units['kernelEntryUnits'],
        'domainEvidenceJoin':candidate['domainEvidence']==domain,
        'pointwiseDomainJoin':candidate['pointwiseDomain']==domain['pointwiseDomain'],
        'separatedPointDomain':isinstance(domain['pointwiseDomain'],sp.Unequality) and len(domain['pointwiseDomain'].free_symbols)==2,
        'normalAxisReal':domain['normalMomentumDomain']==sp.S.Reals,
        'noPolynomialContactTerm':domain['polynomialContactPart']==sp.zeros(5),
        'noDiagonalExtension':domain['diagonalExcluded'] is True and domain['diagonalDistributionExtensionClaimed'] is False,
        'noRetardedOrFormClaim':candidate['complexFrequencyRetardedEquivalenceClaimed'] is False and candidate['fullFORM'] is False,
    }
    # Saved exact branch-exclusion certificates, not new root/sheet searches.
    for index in (14,15):
        prefix='excluded-%d-'%index
        args=operation(prefix+'opposite-sheet-join','input')
        inherited=load('candidate-%d-inherited-operands.pickle'%index)['record']
        structural['excluded%dSavedZero'%index]=operation(prefix+'opposite-sheet-join')==0
        structural['excluded%dOperandJoin'%index]=args==(operation(prefix+'physical-q'),operation(prefix+'radical-scale'),inherited['EXACT_LOW_DEGREE_LIFT']['Q'])
    coverage=load('inherited-real-axis-coverage.pickle')
    structural['inheritedCoverageCertificates']=coverage['nativeCoverage']['FINITE_POLYNOMIAL_ROOT_COVERAGE'] in (True,sp.true) and all(v==0 for v in coverage['normalReality']['CHECKS'].values()) and all(v==0 for v in coverage['exactSpectrumResiduals'].values())
    save(out/'structural-checks.json',{k:bool(v) for k,v in structural.items()})
    journal.value('structural-operands.pickle',{'candidate':candidate,'domain':domain,'context':context,'units':units,'branch':branch,'binding':binding,'checks':structural})
    require(all(structural.values()),'saved structural/source/domain checks')

    certificates=[]
    for name,value in [('determinant',determinant)]+[('adjugate-%d-%d'%(i,j),adjugate[i,j]) for i in range(5) for j in range(5)]:
        chart=operation(name+'-algebraic-chart');chart_args=operation(name+'-algebraic-chart','input')
        cert=operation(name+'-denominator-certificate');cert_args=operation(name+'-denominator-certificate','input')
        proof=operation(name+'-positive-radical-denominator-proof');proof_args=operation(name+'-positive-radical-denominator-proof','input')
        def check_denominator(chart,chart_args,cert,cert_args,proof,proof_args,value):
            radical=chart_args[3];component=proof['component'];coefficients=proof['coefficients']
            evidence=proof['evidence'][component]['polynomial']
            actual_component=sp.re(sp.expand_complex(cert['reducedOnBranch'])) if component=='real' else sp.im(sp.expand_complex(cert['reducedOnBranch']))
            reconstructed=sp.Add(*(v*radical**(len(coefficients)-1-i) for i,v in enumerate(coefficients)))
            return {'sourceEntry':chart_args[0]==value,'chartResidualZero':chart['backJoin']==0,
                'rationalChart':chart['rationalInCoordinates'] is True,
                'denominatorInput':cert_args[0]==chart['rational'],
                'proofInput':proof_args[0]==cert['reducedOnBranch'],
                'positiveRadical':radical.is_positive is True and proof['positiveRadical']==(radical>0),
                'componentCoefficientJoin':sp.expand(reconstructed-evidence)==0 and sp.expand(actual_component-evidence)==0,
                'strictOneSidedCoefficientSign':all((proof['sign']*c).is_nonnegative is True for c in coefficients) and any((proof['sign']*c).is_positive is True for c in coefficients),
                'savedNonzeroCertificate':proof['nonzero'] is True}
        result=journal.op(name+'-saved-certificate-check',check_denominator,chart,chart_args,cert,cert_args,proof,proof_args,value)
        summary={'name':name,'checks':{k:bool(v) for k,v in result.items()}}
        save(out/(name+'-certificate.json'),summary);certificates.append(summary)
        require(all(result.values()),'saved denominator certificate mismatch: '+name)
    tails=[]
    for i in range(5):
        for j in range(5):
            label='inverse-%d-%d'%(i,j);chart=operation(label+'-algebraic-chart')
            chart_args=operation(label+'-algebraic-chart','input')
            require(chart_args[0]==inverse[i,j] and chart['backJoin']==0 and chart['rationalInCoordinates'],'saved inverse chart join')
            for side in (-1,1):
                name=label+('-tail-minus' if side<0 else '-tail-plus')
                tail=operation(name);args=operation(name,'input')
                def check_tail(tail,args,rational,side):
                    t=args[3];power=tail['power'];coefficient=tail['coefficient']
                    series_polynomial=tail['series'].removeO()
                    return {'rationalOperand':args[0]==rational,'side':args[4]==side,
                        'positiveReciprocalCoordinate':t.is_positive is True,
                        'positiveIntegerPower':tail['identicallyZero'] or (power.is_Integer is True and power>0 and not coefficient.has(t)),
                        'savedLeadingTerm':sp.expand(tail['leading']-coefficient*t**power)==0,
                        'savedSeriesLeadingTerm':sp.expand(series_polynomial-tail['leading'])==0,
                        'savedDecay':tail['decays'] in (True,sp.true)}
                checked=journal.op(name+'-saved-check',check_tail,tail,args,chart['rational'],side)
                if name=='inverse-0-0-tail-minus':
                    require(journal.reused==217 and len(journal.records)==217,'failed validation resume boundary')
                    def correct_literal_flag(flag,prior):
                        require(prior['savedDecay'] is False and all(v for k,v in prior.items() if k!='savedDecay'),
                            'only the saved Boolean identity comparison may be corrected')
                        corrected=dict(prior);corrected['savedDecay']=flag in (True,sp.true)
                        return {'checks':corrected,'savedFlagType':type(flag).__module__+'.'+type(flag).__name__,
                            'priorSavedDecay':prior['savedDecay'],'otherChecksUnchanged':all(corrected[k]==v for k,v in prior.items() if k!='savedDecay')}
                    correction=journal.op(name+'-literal-boolean-correction',correct_literal_flag,tail['decays'],checked)
                    save(out/'first-tail-boolean-correction.json',correction)
                    checked=correction['checks']
                result={'row':i,'column':j,'side':side,'power':str(tail['power']),'checks':{k:bool(v) for k,v in checked.items()}}
                save(out/(name+'.json'),result);tails.append(result)
                require(all(checked.values()),'saved tail certificate mismatch: '+name)
    require(len(domain['largeMomentumTails'])==len(tails)==50 and len(domain['denominatorCertificates'])==len(certificates)==26,'complete tail/denominator inventory')

    def numeric(matrix):return np.asarray(sp.ImmutableMatrix(matrix).evalf(50).tolist(),dtype=complex)
    def norm(value):return float(np.linalg.norm(value))
    block_results=[]
    for block in candidate['blocks']:
        index=block['index'];label='block-%d-'%index
        record=load(label+'accepted-operands.pickle');projected=load(label+'projected-operands.pickle')
        bundle=load(label+'checks-and-mutation.pickle')
        point=operation(label+'source-at-pole');q0=operation(label+'physical-q-at-pole')
        residue=operation(label+'exact-residue');d2=operation(label+'det-value-2');a1=operation(label+'adj-value-1')
        require(block['k']==record['EXACT_LOW_DEGREE_LIFT']['K'] and block['residue']==residue==bundle['exactResidue'],'saved exact block identities')
        require(operation(label+'det-value-0')==operation(label+'det-value-1')==0 and d2!=0 and operation(label+'adj-value-0')==sp.zeros(5),'saved multiplicity-two zero orders')
        require(operation(label+'native-physical-q-join')==0,'saved block branch join')
        def check_block(block,record,projected,bundle,point,q0,residue,d2,a1,omega):
            right,left,dk,dw,c,saved_projected,sv=projected
            n=int(record['NULLITY']);direction=block['direction'];a=numeric(residue);p=numeric(point)
            forms=record['FORMS'];current=np.asarray(forms['SIGNED_CURRENT'],complex)
            projector=a@dk;independent=bundle['independentBlockResidue']
            residuals={'normalCoordinateJoin':complex(sp.N(block['k']-record['K'],50)),
                'physicalBulkCoordinateJoin':complex(sp.N(q0-record['PHYSICAL_Q'],50)),
                'frequencyCoordinateJoin':complex(sp.N(omega-record['OMEGA'],50)),
                'actualPencilJoin':p-np.asarray(record['OPERANDS']['PENCIL_PLUS'],complex),
                'projectedDerivativeJoin':left.conj().T@dk@right-c,
                'cofactorVsBlockResidue':a-independent,'leftLaurent':p@a,'rightLaurent':a@p,
                'derivativeProjector':projector@projector-projector,
                'fluxFrequencyNormalization':np.asarray(forms['FLUX_LEFT'],complex).conj().T@dw@np.asarray(forms['FLUX_RIGHT'],complex)-np.eye(n),
                'outgoingCurrentSign':direction*current-np.eye(n)}
            norms={k:norm(v) for k,v in residuals.items()};norms['cofactorVsBlockResidue']/=max(1.,norm(independent))
            roundtrip=max(norm(v-bundle['residuals'][k]) for k,v in residuals.items())
            mutation=(-direction)*current-np.eye(n)
            delta=sp.I*sp.pi*direction*residue
            delta_join=block['deltaCoefficient']==bundle['deltaCoefficient']==delta
            changed_delta=-delta;movement=norm(numeric(changed_delta-delta))
            exact_residual=(residue*d2-2*a1).applyfunc(lambda v:sp.cancel(sp.expand(v)))
            summary={'index':index,'nullity':n,'direction':direction,'norms':norms,
                'residualRoundTripNorm':roundtrip,'savedNormDifference':max(abs(norms[k]-bundle['norms'][k]) for k in norms),
                'savedProjectedDifference':norm(left.conj().T@dk@right-saved_projected),
                'currentOperandDifference':norm(current-bundle['signedCurrent']),
                'mutationRoundTripNorm':norm(mutation-bundle['mutatedCurrentResidual']),
                'mutationNorm':norm(mutation),'deltaMovementNorm':movement,
                'savedDeltaMovementDifference':abs(movement-bundle['deltaMovementNorm']),
                'deltaCoefficientJoin':delta_join,'mutatedDeltaJoin':changed_delta==bundle['mutatedDeltaCoefficient'],
                'exactCofactorResidueIdentity':exact_residual==sp.zeros(5),
                'fullBlockShapes':right.shape==left.shape==(5,n) and c.shape==(n,n),
                'savedProjectedRegular':bool(sv[-1]>1e-10*max(1.,sv[0]))}
            return {'residuals':residuals,'mutationResidual':mutation,'exactResidueResidual':exact_residual,'summary':summary}
        result=journal.op(label+'saved-residual-validation',check_block,block,record,projected,bundle,point,q0,residue,d2,a1,
            sp.Rational(candidate['fixedPhysicalInput']['parameters']['omega']))
        s=result['summary'];save(out/(label+'validation.json'),s);block_results.append(s)
        require(all(v<1e-8 for v in s['norms'].values()),'full saved block residuals')
        require(all(s[k]<1e-13 for k in ['residualRoundTripNorm','savedNormDifference','savedProjectedDifference','currentOperandDifference','mutationRoundTripNorm','savedDeltaMovementDifference']),'saved block roundtrip')
        require(s['mutationNorm']>1 and s['deltaMovementNorm']>1e-8,'responsive saved direction mutation')
        require(all(s[k] for k in ['deltaCoefficientJoin','mutatedDeltaJoin','exactCofactorResidueIdentity','fullBlockShapes','savedProjectedRegular']),'exact block and signed residue structure')
    require([(s['index'],s['direction']) for s in block_results]==[(17,-1),(16,1)],'complete opposite outgoing blocks')

    def check_assembly(candidate,context,domain):
        kn=context['normalMomentum'];fourier=context['fourierMass'];phase=context['normalCharacter']
        pv=candidate['principalValueKernel'];intervals=candidate['principalValueExclusionIntervals']
        exclusion=pv[0,0].args[1];points=[b['k'] for b in candidate['blocks']]
        expected_intervals=tuple(zip([-sp.oo]+[k+exclusion for k in points],[k-exclusion for k in points]+[sp.oo]))
        integral_checks=[]
        for i in range(5):
            for j in range(5):
                entry=pv[i,j]
                terms=sp.Add.make_args(entry.args[0])
                integral_checks.append(isinstance(entry,sp.Limit) and entry.args[1]==exclusion
                    and entry.args[2]==0 and str(entry.args[3])=='+'
                    and len(terms)==len(intervals)
                    and all(isinstance(term,sp.Integral) and term.function==candidate['regularDensity'][i,j] for term in terms)
                    and {term.limits for term in terms}=={((kn,a,b),) for a,b in intervals})
        correction=sp.zeros(5);distribution=sp.zeros(5)
        for block in candidate['blocks']:
            correction+=phase.subs(kn,block['k'])*block['deltaCoefficient']/fourier
            distribution+=block['deltaCoefficient']*sp.DiracDelta(kn-block['k'])
        return {'checks':{'symmetricSavedExclusions':intervals==expected_intervals,
            'positiveExclusion':exclusion.is_positive is True,'all25ExplicitPVEntries':all(integral_checks),
            'deltaCorrectionAssembly':candidate['singularKernelCorrection']==correction,
            'momentumDeltaAssembly':candidate['momentumDeltaCorrection']==distribution,
            'totalKernelAssembly':candidate['candidateKernel']==pv+correction,
            'exclusionUnits':candidate['exclusionMomentumUnit']==domain['radicalUnit']},
            'entryChecks':integral_checks,'expectedCorrection':correction,'expectedMomentumDelta':distribution,
            'expectedIntervals':expected_intervals}
    assembly=journal.op('saved-integral-assembly-validation',check_assembly,candidate,context,domain)
    save(out/'assembly-checks.json',{k:bool(v) for k,v in assembly['checks'].items()})
    require(all(assembly['checks'].values()),'explicit saved PV/delta assembly')
    support=load('spatial-sign-identity.pickle')
    require(len(support)==4 and all(r['residual']==0 for r in support),'saved four contour wiring identities')
    post={str(p.relative_to(source)):sha(p) for p in sorted(source.rglob('*')) if p.is_file()}
    save(out/'posthashes.json',post);require(post==hashes,'saved candidate artifacts unchanged')
    resume_post={str(p.relative_to(resume_root)):sha(p) for p in sorted(resume_root.rglob('*')) if p.is_file()}
    save(out/'resume-posthashes.json',resume_post);require(resume_post==resumed_hashes,'prior validation bytes changed')
    save(out/'operation-index.json',journal.records)
    result={'status':'SAVED_OUTGOING_PRESCRIPTION_CANDIDATE_VALIDATED','sourceRoot':str(source),
        'savedFiles':len(hashes),'validationOperations':len(journal.records),
        'restoredPriorValidationOperations':journal.reused,'newValidationOperations':len(journal.records)-journal.reused,
        'priorValidationArtifactsUnchanged':True,
        'structuralChecks':{k:bool(v) for k,v in structural.items()},'blockChecks':block_results,
        'denominatorCertificates':len(certificates),'tailChecks':len(tails),
        'assemblyChecks':{k:bool(v) for k,v in assembly['checks'].items()},'allSourceHashesUnchanged':True,
        'newConstructors':0,'newModeSolves':0,'newLUSolves':0,'newRootSearches':0,
        'integralsEvaluated':False,'contourIdentityIsWiringOnly':True,
        'pointwiseDomain':str(candidate['pointwiseDomain']),'diagonalDistributionExtensionClaimed':False,
        'outgoingGreenOperatorAccepted':False,'formCompleted':False,'a11Cleared':False,'a12Cleared':False,
        'freshIndependentMethodClearClaimed':False,'wallSeconds':time.monotonic()-started}
    save(out/'checks.json',result);print(json.dumps(result,indent=2,allow_nan=False))


def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--run-directory',type=Path,required=True)
    args=parser.parse_args();observed=containment()
    source=args.source.resolve();source.relative_to(STORE)
    out=args.run_directory.resolve();out.relative_to(STORE)
    require(out!=source and source not in out.parents,'separate saved-validation output required')
    out.mkdir(parents=True,exist_ok=False);save(out/'native-containment.json',observed)
    try:validate(source,args.manifest,out)
    except BaseException:
        save(out/'failure.json',{'traceback':traceback.format_exc(),'automaticRetry':False});raise


if __name__=='__main__':main()
