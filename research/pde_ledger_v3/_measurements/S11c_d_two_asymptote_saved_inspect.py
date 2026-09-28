#!/usr/bin/env python3
"""Resume only unfinished two-asymptote end-lift/forcing work.

No producer imports, root/mode solver, inverse reconstruction, numerical
response, integral evaluation or automatic continuation. A fresh method gate
and the shared whole-job guard are required before importing SymPy.
"""
import argparse
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import pickle
import re
import resource
import signal
import shutil
import time
import traceback

ROOT = Path('/var/projects/toy_physics')
STORE = ROOT / '_scratch/s11c'
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')
GRADES = ((0, 0), (1, 0), (0, 1))


def require(test, message):
    if not test:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def route(path):
    path = Path(path)
    return {'path': str(path), 'canonicalPath': str(path.resolve(strict=True)),
            'bytes': path.stat().st_size, 'sha256': digest(path)}


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n'); stream.flush(); os.fsync(stream.fileno())


def containment():
    group = next(s[3:] for s in Path('/proc/self/cgroup').read_text().splitlines()
                 if s.startswith('0::'))
    root = Path('/sys/fs/cgroup') / group.lstrip('/')
    result = {k: (root/k).read_text().strip()
              for k in ('memory.max', 'memory.swap.max', 'pids.max')}
    result.update(nice=os.getpriority(os.PRIO_PROCESS, 0),
                  affinity=sorted(os.sched_getaffinity(0)),
                  threads={k: os.environ.get(k) for k in THREADS})
    require(result['memory.max'] == str(2*1024**3) and result['memory.swap.max'] == '0'
            and result['pids.max'] == '32' and result['nice'] >= 15
            and len(result['affinity']) == 1
            and all(v == '1' for v in result['threads'].values()),
            'required whole-job containment absent; no science imported')
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    def timeout(*_):
        raise TimeoutError('840-second native limit; preserve incomplete input')
    signal.signal(signal.SIGALRM, timeout); signal.alarm(840)
    return {**result, 'nativeWallSeconds': 840}


class SavedCodec(pickle.Unpickler):
    def find_class(self, module, name):
        if not (module.startswith(('sympy.', 'numpy.'))
                or module in ('sympy', 'numpy', 'builtins', 'collections')):
            raise pickle.UnpicklingError(('unapproved saved class', module, name))
        return super().find_class(module, name)




class Journal:
    """Restore prior returns without calling their functions; journal new work."""
    def __init__(self, out):
        self.out, self.records, self.stack = out, [], []
        self.next_number = 0
        self.restored_parent = self.restored_diagnostic = 0
        self.restored_artifacts = []

    @property
    def active(self):
        return self.stack[-1] if self.stack else None

    def value(self, name, value):
        path = self.out/name; path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as stream:
            pickle.dump(value, stream, protocol=4)
            stream.flush(); os.fsync(stream.fileno())
        return route(path)

    def copy(self, name, record):
        path = Path(record['path'])
        require(route(path) == record, ('saved bytes changed', path))
        target = self.out/name; target.parent.mkdir(parents=True, exist_ok=True)
        with path.open('rb') as src, target.open('xb') as dest:
            shutil.copyfileobj(src, dest); dest.flush(); os.fsync(dest.fileno())
        result = route(target)
        require(result['sha256'] == record['sha256'], 'copied return identity')
        return result

    def begin(self, name):
        folder = 'operations/%04d-%s' % (self.next_number, name)
        self.next_number += 1
        return folder, {'name': name, 'startedUtc': datetime.now(timezone.utc).isoformat()}

    def restore(self, previous, diagnostic=False):
        name = ('diagnostic-' if diagnostic else '')+previous['name']
        folder, record = self.begin(name)
        record.update(input=self.copy(folder+'/input.pickle', previous['input']),
                      execution=('RESTORED_DIAGNOSTIC_COMPLETE_RETURN' if diagnostic
                                 else 'RESTORED_PRIOR_COMPLETE_RETURN'), priorOperation=previous)
        save(self.out/folder/'started.json', record); self.stack.append(record)
        record['value'] = self.copy(folder+'/value.pickle', previous['value'])
        with Path(record['value']['path']).open('rb') as stream:
            result = SavedCodec(stream).load()
        record['finishedUtc'] = datetime.now(timezone.utc).isoformat()
        save(self.out/folder/'completed.json', record)
        self.records.append(record); self.stack.pop()
        if diagnostic: self.restored_diagnostic += 1
        else: self.restored_parent += 1
        return result

    def artifact(self, name, record, restore=False):
        result = self.copy(name, record)
        self.restored_artifacts.append({'source':record,'copy':result,
                                       'execution':'COPIED_PRIOR_COMPLETE_ARTIFACT'})
        if restore:
            with Path(result['path']).open('rb') as stream:
                return SavedCodec(stream).load()
        return result

    def op(self, name, function, *args):
        folder, record = self.begin(name)
        record.update(input=self.value(folder+'/input.pickle', args),
                      execution='NEW_UNFINISHED_OPERATION')
        save(self.out/folder/'started.json', record); self.stack.append(record)
        result = function(*args)
        record.update(value=self.value(folder+'/value.pickle', result),
                      finishedUtc=datetime.now(timezone.utc).isoformat())
        save(self.out/folder/'completed.json', record)
        self.records.append(record); self.stack.pop()
        return result

    def reuse(self, name, previous, result, args):
        folder, record = self.begin(name)
        record.update(input=self.value(folder+'/input.pickle',args),
                      execution='REUSED_COMPLETE_RETURN_IDENTICAL_OPERANDS',priorOperation=previous)
        save(self.out/folder/'started.json',record);self.stack.append(record)
        record.update(value=self.copy(folder+'/value.pickle',previous['value']),
                      finishedUtc=datetime.now(timezone.utc).isoformat())
        save(self.out/folder/'completed.json',record)
        self.records.append(record);self.stack.pop()
        return result

    def resume_input(self, name, function, previous_input):
        # The failed call never returned. Copy its exact saved input before
        # executing only that unfinished function under the existing guard.
        folder, record = self.begin(name)
        record.update(input=self.copy(folder+'/input.pickle',previous_input),
                      execution='NEW_UNFINISHED_OPERATION_FROM_SAVED_INPUT',
                      priorIncompleteInput=previous_input)
        save(self.out/folder/'started.json',record);self.stack.append(record)
        with Path(record['input']['path']).open('rb') as stream:
            args = SavedCodec(stream).load()
        result = function(*args)
        record.update(value=self.value(folder+'/value.pickle',result),
                      finishedUtc=datetime.now(timezone.utc).isoformat())
        save(self.out/folder/'completed.json',record)
        self.records.append(record);self.stack.pop()
        return result


def inspect_saved(spec, out, journal):
    import sympy as sp
    from collections import Counter
    j = journal
    source = Path(spec['sourceRoot']).resolve()
    source.relative_to(STORE)
    actual = {name:route(r['path']) for name,r in spec['inputs'].items()}
    save(out/'prehashes.json',actual)
    require(actual == spec['inputs'],'saved-output/source hashes changed')
    def metadata(name):
        require(name in spec['sourceFiles'],'unmanifested source JSON')
        return json.loads((source/name).read_text())
    operations = {r['name']:r for r in metadata('operation-index.json')}
    cache = {}
    def load(relative):
        require(relative in spec['sourceFiles'],'unmanifested saved object')
        if relative not in cache:
            def read_saved(name):
                path = source/name
                path.resolve(strict=True).relative_to(source)
                require(route(path) == actual[spec['sourceFiles'][name]],'saved object hash mismatch')
                with path.open('rb') as stream:return SavedCodec(stream).load()
            cache[relative] = j.op('read-saved-%03d'%len(cache),read_saved,relative)
        return cache[relative]
    def operation(name,part='value'):
        require(name in operations,('missing completed operation',name))
        relative = str(Path(operations[name][part]['path']).relative_to(source))
        return load(relative)
    def text_matrix(matrix):
        return {'shape':list(matrix.shape),'entries':[
            {'row':i,'column':k,'expression':sp.sstr(matrix[i,k])}
            for i in range(matrix.rows) for k in range(matrix.cols) if matrix[i,k] != 0]}
    def integral_text(value):
        return [{'expression':sp.sstr(v),'orderedLimits':[
            [sp.sstr(w) for w in lim] for lim in v.limits]}
            for v in sorted(value.atoms(sp.Integral),key=sp.default_sort_key)]
    def json_value(value):
        if value is None or isinstance(value,(str,bool,int,float)):return value
        if isinstance(value,dict):return {str(k):json_value(v) for k,v in value.items()}
        if isinstance(value,(tuple,list)):return [json_value(v) for v in value]
        if isinstance(value,sp.Integer):return int(value)
        if value is sp.true or value is sp.false:return bool(value)
        return sp.sstr(value)
    def emit(name,result):
        save(out/(name+'.json'),json_value(result))
        require(all(result['checks'].values()),('saved payload join/check failed',name))

    binding = load('units-and-binding.pickle')
    partition = load('auxiliary-partition.pickle')
    premise = load('plane-jet-premise.pickle')
    prepared = operation('bind-native-action')
    endpoint = operation('endpoint-binding-source-evidence')
    def inspect_context(binding,partition,premise,prepared,endpoint):
        cells=[]
        for (grade,i,k),cell in sorted(prepared.items()):
            local=[{'order':n,'coefficient':sp.sstr(c)} for n,c in cell['local'] if c != 0]
            nonlocal_terms=[{'term':t,'coefficient':sp.sstr(c),'integral':sp.sstr(v),
                            'orderedIntegrals':integral_text(v)}
                           for t,c,v in cell['nonlocal'] if c != 0 and v != 0]
            cells.append({'grade':list(grade),'row':i,'column':k,
                          'local':local,'nonlocal':nonlocal_terms})
        return {'checks':{
            'physicalBindingPresent':bool(binding['physicalInput']),
            'partitionDoesNotChangePhysics':not bool(partition['physicalInputChanged']),
            'partitionHeldFixed':bool(partition['heldFixedUnderBackgroundGrades']),
            'nativeSubtractionUnchanged':bool(partition['sourceCanonicalHeavisideSubtractionUnchanged']),
            'sameRegulator':partition['profileAbelRegulator']==binding['profileRegulator'],
            'actualEndpointMapReused':endpoint['sourceEndpointMap']==binding['profileEndpoints'],
            'allActualEndpointLimitsMatched':not endpoint['unknown'] and len(endpoint['matched'])==len(endpoint['actualNativeLimits']),
            'nonlocalPlaneJetRemainsUnverified':premise['status']=='UNVERIFIED_NONLOCAL_PLANE_JET_EXTENSION'
                and not bool(premise['translationInvarianceEstablished'])
                and not bool(premise['momentumDerivativePassesOrderedIntegralsAndAbelLimitEstablished'])},
            'physicalInput':binding['physicalInput'],'profileRegulator':sp.sstr(binding['profileRegulator']),
            'profileFormulas':{k:sp.sstr(v) for k,v in binding['profileFormulas'].items()},
            'partition':{k:sp.sstr(v) for k,v in partition['chi'].items()},
            'endpointBindings':[{'limit':sp.sstr(k),'value':sp.sstr(v)} for k,v in endpoint['matched']],
            'preparedNativeCells':cells,'cellCount':len(cells),
            'zeroGradeRegulatorPresent':bool(premise['zeroGradeContainsAbelRegulator']),
            'zeroGradeProfileOccurrences':[sp.sstr(v) for v in premise['zeroGradeProfileOccurrences']],
            'forcingUnits':binding['forcingUnits'],'fieldEntryUnits':binding['endFieldUnits'],
            'integralsOrLimitsEvaluated':False,'newPhysicsMethod':False}
    context = j.op('inspect-saved-context-and-native-cells',inspect_context,
                   binding,partition,premise,prepared,endpoint)
    emit('saved-context',context)

    cases={};case_results=[];end_results=[]
    for block in (17,16):
        for grade in ((1,0),(0,1)):
            tag='block%d-grade%d%d'%(block,*grade)
            data={'forcing':load(tag+'-forcing-operands.pickle'),
                  'plane':load(tag+'-plane-jet-reduction-operands.pickle'),
                  'domain':load(tag+'-domain.pickle'),
                  'localTail':load(tag+'-local-tail-input.pickle'),
                  'carriers':operation(tag+'-reconstruction-carriers'),
                  'carrierInput':operation(tag+'-reconstruction-carriers','input'),
                  'reduced':operation(tag+'-reconstruction-coefficients-zero'),
                  'reductionInput':operation(tag+'-reconstruction-coefficients-zero','input'),
                  'forcingReturn':operation(tag+'-native-forcing'),
                  'liftReturn':operation(tag+'-native-lift-action'),
                  'assembly':operation(tag+'-assemble-decomposition'),
                  'sideReturns':{e:operation(tag+'-'+e+'-cutoff-action') for e in ('LEFT','RIGHT')},
                  'planeReturns':{e:operation(tag+'-'+e+'-native-plane-action') for e in ('LEFT','RIGHT')},
                  'summary':metadata(tag+'-summary.json')}
            cases[tag]=data
            def inspect_case(tag,block,grade,data,regulator):
                f,p,d,tail,carriers = (data[k] for k in ('forcing','plane','domain','localTail','carriers'))
                obligations=d['nonlocalObligations']
                source_terms=[]
                for side,record in [('PROFILE_FORCING',f['forcing']),*f['sideActions'].items()]:
                    for term in record['terms']:
                        if term['value'] != 0:
                            source_terms.append((side,{k:term[k] for k in ('row','sourceColumn','frameColumn','term')},term['value']))
                obligation_join=len(obligations)==len(source_terms)
                for item,(side,address,value) in zip(obligations,source_terms):
                    obligation_join = obligation_join and item['side']==side and item['address']==address and item['operand']==value
                    obligation_join = obligation_join and item['orderedIntegralLimits']==[
                        v.limits for v in sorted(value.atoms(sp.Integral),key=sp.default_sort_key)]
                    obligation_join = obligation_join and bool(item['hasProfileAbelRegulator'])==bool(value.has(regulator))
                conditions={
                    'fiveByFiveDirectAndDecomposed':f['direct'].shape==f['decomposed'].shape==(5,5),
                    'forcingOperationJoin':f['forcing']==data['forcingReturn'],
                    'liftOperationJoin':f['liftedAction']==data['liftReturn'],
                    'sideActionOperationJoins':f['sideActions']==data['sideReturns'],
                    'assemblyReturnJoin':f['decomposed']==data['assembly']['decomposed']
                        and f['conditionalCommutators']==data['assembly']['conditionalCommutators'],
                    'directFromSavedTotals':f['direct']==-f['forcing']['total']-f['liftedAction']['total'],
                    'carrierInputJoin':data['carrierInput']==(f['direct']-f['decomposed'],),
                    'reductionInputJoin':data['reductionInput']==(carriers['coefficients'],),
                    'savedReducedCoefficientsZero':data['reduced']==sp.zeros(*data['reduced'].shape),
                    'rawPlaneOperationJoins':p['nativePlaneActions']==data['planeReturns'],
                    'sourcePlaneOperandJoin':p['sourceSymbolPlaneActions']==f['referencePlaneJetActions'],
                    'planePremiseJoin':p['premise']==f['planeJetReductionPremise'],
                    'directDomainOperandJoin':d['directForcing']==f['direct'],
                    'profileForcingTermJoin':d['profileForcingTerms']==f['forcing']['terms'],
                    'nativeObligationsAndOrderedLimitsJoin':obligation_join,
                    'domainCaseIdentity':d['block']==block and tuple(d['grade'])==grade,
                    'nonlocalObligationCount':len(obligations)==data['summary']['nonlocalObligations'],
                    'localLimitCount':len(d['localWeightedTailLimits'])==data['summary']['localWeightedLimits']==50,
                    'localLimitsUnevaluated':all(isinstance(v['operand'],sp.Limit) for v in d['localWeightedTailLimits']),
                    'fullDefectLimitsPreserved':set(d['fullDefectEndLimitOperands'])=={'-oo','oo'}
                        and all(v.shape==(5,5) and all(isinstance(x,sp.Limit) for x in v) for v in d['fullDefectEndLimitOperands'].values()),
                    'localTailOperandJoin':tail['localPart']==-f['forcing']['local']-f['liftedAction']['local'],
                    'pairingStillUnresolved':not bool(f['wholeLinePairingEstablished'])
                        and not bool(d['schwartzMembershipEstablished']) and not bool(d['wholeLineOutgoingActionEstablished'])
                        and d['status']=='FORCING_DOMAIN_UNRESOLVED',
                    'independentPoleExclusionRetained':bool(d['realPoleExclusionIndependentOfProfileAbel'])}
                terms=[{'side':v['side'],'address':v['address'],'operand':sp.sstr(v['operand']),
                        'orderedIntegrals':integral_text(v['operand']),
                        'hasProfileAbelRegulator':bool(v['hasProfileAbelRegulator']),'status':v['status']}
                       for v in obligations]
                return {'case':tag,'checks':conditions,'directForcing':text_matrix(f['direct']),
                    'localForcingPart':text_matrix(tail['localPart']),
                    'lift':text_matrix(f['lift']),'seed':text_matrix(f['seed']),
                    'carrierCount':len(carriers['carriers']),
                    'unreducedCoefficientZeroFlags':[str(v.is_zero) for v in carriers['coefficients']],
                    'reducedShape':list(data['reduced'].shape),
                    'nativeTermsBySide':dict(Counter(v['side'] for v in obligations)),
                    'regulatorPresentCount':sum(bool(v['hasProfileAbelRegulator']) for v in obligations),
                    'nonlocalDomainOperands':terms,'profileAbelLimitOrder':d['profileAbelLimitOrder'],
                    'planeJetStructuralEquality':{e:p['nativePlaneActions'][e]['total']==p['sourceSymbolPlaneActions'][e] for e in ('LEFT','RIGHT')},
                    'planeJetEstablished':False,'forcingPairingEstablished':False,
                    'localLimitOperands':[{'row':v['row'],'column':v['column'],'end':v['end'],
                                          'operand':sp.sstr(v['operand'])} for v in d['localWeightedTailLimits']],
                    'integralsOrLimitsEvaluated':False}
            result=j.op(tag+'-inspect-saved-forcing',inspect_case,tag,block,grade,data,binding['profileRegulator'])
            emit(tag+'-inspection',result)
            case_results.append({'case':tag,'checks':result['checks'],'carrierCount':result['carrierCount'],
                                 'nativeTermsBySide':result['nativeTermsBySide'],
                                 'regulatorPresentCount':result['regulatorPresentCount'],
                                 'planeJetStructuralEquality':result['planeJetStructuralEquality']})
            for end in ('LEFT','RIGHT'):
                end_tag=end+'-'+tag
                field=load(end_tag+'-field.pickle')
                residuals=[operation(end_tag+'-end-action-power'+str(n)) for n in (0,1)]
                principal=operation(end_tag+'-principal-part')
                def inspect_end(name,field,residuals,principal,units):
                    return {'endField':name,'checks':{
                        'fieldShape':field['field'].shape==field['polynomial'].shape==(5,5),
                        'principalPartReturnJoin':field['principalPart']==principal,
                        'bothSavedEquationResidualsZero':all(v==sp.zeros(5) for v in residuals),
                        'entryUnits':field['fieldEntryUnits']==units['endFieldUnits'],
                        'columnUnits':field['columnUnits']==units['residueColumnUnits']},
                        'field':text_matrix(field['field']),'polynomial':text_matrix(field['polynomial']),
                        'integralsOrLimitsEvaluated':False}
                observed=j.op(end_tag+'-inspect-saved-field',inspect_end,end_tag,field,residuals,principal,binding)
                emit(end_tag+'-inspection',observed)
                end_results.append({'endField':end_tag,'checks':observed['checks']})

    tag='block17-grade10'
    mutation=load(tag+'-nonlocal-mutation-0.pickle')
    end_mutation=load(tag+'-end-mutation.pickle')
    mutation_returns={'carriers':operation(tag+'-mutated-carriers-0'),
        'reduced':operation(tag+'-mutated-coefficients-0'),
        'omission':operation(tag+'-nonlocal-omission-0'),
        'carrierInput':operation(tag+'-mutated-carriers-0','input'),
        'reductionInput':operation(tag+'-mutated-coefficients-0','input'),
        'endReturns':{e:operation(tag+'-'+e+'-end-mutation') for e in ('LEFT','RIGHT')},
        'endFields':{e:load(e+'-'+tag+'-field.pickle') for e in ('LEFT','RIGHT')}}
    def inspect_mutations(m,em,returns,baseline):
        f=baseline['forcing'];side,row,column,frame,term_id=m['address']
        actual=[v for v in f['sideActions'][side]['terms']
                if (v['row'],v['sourceColumn'],v['frameColumn'],v['term'])==(row,column,frame,term_id)]
        carrier=m['carrierReduction'];reduced=m['reducedCoefficients']
        responsive=[i for i,(pos,v) in enumerate(zip(carrier['positions'],reduced))
                    if any(pos['carrierPowers']) and v.is_zero is False]
        end_matches=[e for e in ('LEFT','RIGHT') if em['removedEndLift']==returns['endFields'][e]['field']
                     and em['heldFixedEndForcing']==f['endPlaneActions'][e]
                     and em['residualWithoutLift']==returns['endReturns'][e]]
        end_nonzero=[(i,k,v) for i in range(5) for k in range(5)
                     if (v:=em['residualWithoutLift'][i,k]).is_zero is False]
        conditions={'removedOneActualNativeTerm':len(actual)==1 and m['removedActualTerms']==actual,
            'omissionReturnJoin':m['mutatedDecomposition']==returns['omission']['decomposed']
                and returns['omission']['omission']==m['address'],
            'baselineHeldFixed':m['baselineDirect']==f['direct'] and m['baselineDecomposition']==f['decomposed'],
            'savedMutationResidualJoin':m['actualResidual']==m['baselineDirect']-m['mutatedDecomposition'],
            'carrierReturnJoin':carrier==returns['carriers'] and returns['carrierInput']==(m['actualResidual'],),
            'reductionReturnJoin':reduced==returns['reduced'] and returns['reductionInput']==(carrier['coefficients'],),
            'responsivePositionsMatchActualCoefficients':bool(responsive) and responsive==m['exactNonzeroCarrierCoefficientPositions'],
            'formalControlFlag':bool(m['formalControlResponsive']),
            'noIntegratedMutationClaim':not bool(m['integratedResponseNonzeroEstablished']),
            'endControlUsesActualFieldAndReturn':bool(end_matches) and bool(end_nonzero),
            'endControlSeparate':bool(em['responsiveExactCoefficient']) and not bool(em['nonlocalControlSubstitute'])}
        return {'checks':conditions,'actualOmissionAddress':list(m['address']),
                'removedNativeTerm':sp.sstr(actual[0]['value']) if actual else None,
                'responsiveCarrierCoefficients':[{'index':i,'position':carrier['positions'][i],
                    'coefficient':sp.sstr(reduced[i]),'isZero':str(reduced[i].is_zero)} for i in responsive],
                'matchingEndControlSides':end_matches,
                'nonzeroEndResidualEntries':[{'row':i,'column':k,'coefficient':sp.sstr(v),'isZero':str(v.is_zero)}
                                            for i,k,v in end_nonzero],
                'integratedMutationEstablished':False,'integralsOrLimitsEvaluated':False}
    controls=j.op('inspect-saved-responsive-controls',inspect_mutations,
                  mutation,end_mutation,mutation_returns,cases[tag])
    emit('saved-controls',controls)
    summary={'status':'SAVED_END_LIFT_FORCING_PAYLOADS_INSPECTED_PAIRING_UNRESOLVED',
        'contextChecks':context['checks'],'caseChecks':case_results,'endFieldChecks':end_results,
        'controlChecks':controls['checks'],'actualOmissionAddress':controls['actualOmissionAddress'],
        'responsiveCarrierCoefficientCount':len(controls['responsiveCarrierCoefficients']),
        'nonzeroEndResidualCount':len(controls['nonzeroEndResidualEntries']),
        'savedPayloadsRead':len(cache),'completedOperations':len(j.records),
        'sourceCases':4,'savedEndFields':8,'integralsEvaluated':False,'limitsEvaluated':False,
        'producerCalls':0,'newRoots':0,'newModeSolves':0,'newLUSolves':0,
        'derivativesReconstructed':False,'constructorReplayed':False,
        'forcingPairingEstablished':False,'planeJetExtensionEstablished':False,
        'responseConstructed':False,'formCompleted':False,'a11Cleared':False,'a12Cleared':False,
        'freshIndependentClearClaimed':False,'automaticRetry':False,
        'scope':'Saved-object structural joins, returned zero matrices and actual responsive coefficients; readable native/domain operands. No integral/limit evaluation or replay of constructor functions.'}
    summary=json_value(summary)
    save(out/'payload-inspection-summary.json',summary)
    return summary


def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument('--input-manifest',type=Path,required=True)
    parser.add_argument('--gate-receipt',type=Path,required=True)
    parser.add_argument('--run-directory',type=Path,required=True)
    args=parser.parse_args()
    spec=json.loads(args.input_manifest.read_text());gate=json.loads(args.gate_receipt.read_text())
    require(gate['status']=='READY_FOR_ONE_GUARDED_TWO_ASYMPTOTE_SAVED_INSPECTION','method gate incomplete')
    require(gate['workerSha256']==digest(Path(__file__)) and gate['inputManifestSha256']==digest(args.input_manifest),
            'method gate worker/input identity')
    require(gate['scopeExplicitlyApproved'] is True and gate['substantiveReviewFindingsClosed'] is True,
            'scope or independent method-review disposition missing')
    require(gate['seconds']==900 and gate['nativeSeconds']==840 and gate['automaticRetry'] is False,
            'bounded stage duration/retry contract')
    require(gate['sharedGuardSha256']==digest(ROOT/'scripts/s11c_guarded_run.py'), 'shared guard identity')
    require(gate['supervisorSha256']==digest(ROOT/'research/pde_ledger_v3/_measurements/S11c_d_end_normalization_run.py'), 'supervisor identity')
    approval_path=Path(gate['scopeApprovalPath'])
    require(gate['scopeApprovalSha256']==digest(approval_path),'continuation approval identity')
    approval=json.loads(approval_path.read_text())
    require(approval['status']=='AUTHORIZED_ONE_GUARDED_TWO_ASYMPTOTE_SAVED_INSPECTION'
            and approval['workerSha256']==digest(Path(__file__))
            and approval['inputManifestSha256']==digest(args.input_manifest), 'continuation approval scope')
    require(gate['correctionDispositionSha256']==digest(Path(gate['correctionDispositionPath'])),
            'bounded correction disposition identity')
    require(str(args.run_directory.resolve())==gate['resultDirectory'], 'continuation result directory')
    require(gate['savedOutputInspectionOnly'] is True,'saved-output scope required')
    observed=containment()
    out=args.run_directory.resolve();out.relative_to(STORE);out.mkdir(parents=True,exist_ok=False)
    save(out/'native-containment.json',observed);save(out/'input-manifest.json',spec);save(out/'gate-receipt.json',gate)
    journal=Journal(out);started=time.monotonic()
    try:
        result=inspect_saved(spec,out,journal)
        posthashes={name:route(v['path']) for name,v in spec['inputs'].items()}
        save(out/'posthashes.json',posthashes)
        require(posthashes==spec['inputs'],'input posthash mismatch')
        result.update(wallSeconds=time.monotonic()-started,allSourceHashesUnchanged=True,
                      restoredParentOperations=journal.restored_parent,
                      restoredDiagnosticOperations=journal.restored_diagnostic,
                      priorFunctionsExecuted=0)
        save(out/'operation-index.json',journal.records)
        save(out/'checks.json',result)
        save(out/'artifact-index.json',{str(p.relative_to(out)):route(p) for p in sorted(out.rglob('*')) if p.is_file()})
        print(json.dumps(result,indent=2,allow_nan=False))
    except BaseException:
        failure_traceback = traceback.format_exc()
        signal.alarm(0)  # bookkeeping only; no further science after a failure
        post = {}
        for name, record in spec['inputs'].items():
            try: post[name] = route(record['path'])
            except OSError as error: post[name] = {'error':str(error),'path':record['path']}
        save(out/'failure-posthashes.json',post)
        save(out/'partial-operation-index.json',journal.records)
        save(out/'failure.json',{'traceback':failure_traceback,'completedOperations':len(journal.records),
                                'incompleteOperation':journal.active,'incompleteOperationStack':journal.stack,
                                'restoredParentOperations':journal.restored_parent,
                                'restoredDiagnosticOperations':journal.restored_diagnostic,
                                'restoredArtifacts':journal.restored_artifacts,
                                'wallSeconds':time.monotonic()-started,'automaticRetry':False})
        save(out/'failure-artifact-index.json',{str(p.relative_to(out)):route(p) for p in sorted(out.rglob('*')) if p.is_file()})
        raise


if __name__=='__main__':
    main()
