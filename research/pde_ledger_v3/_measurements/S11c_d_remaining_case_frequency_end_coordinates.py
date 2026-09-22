#!/usr/bin/env python3
"""Missing coordinate images of accepted live-frequency current operands.

No analytic, mapping, binding, derivative, current or mode call is repeated.
The new images are explicitly coordinate conversions, not reconstructed
serialized returns from an old process or accepted numerical end operators.
"""
import argparse
import ast
import copy
import hashlib
import inspect
import json
from pathlib import Path
import resource
import signal
import textwrap
import time
import S11c_d_remaining_case_frequency_end_calls as previous

h,f,sp,engine=previous.h,previous.f,previous.sp,previous.engine
CP=f.M/'S11c_d_remaining_case_frequency_end_calls_checkpoint.json'
MCP=f.M/'S11c_d_remaining_case_modes_checkpoint.json'
PLAN=f.M/'S11c_d_remaining_case_frequency_end_coordinates_plan.md'
CP_SHA='5b2b392b0d7d4637aef02f7f604f85a7d1af31cb2f4df1f06a1472c99bb9a869'
OWNER=('LAB_HELD__RHOBR_CONSTANT','RIGHT')


def same(a,b):f.require(previous.modes_reader.same(a,b),'exact full saved coordinate input')


def source_body(value):
    source=textwrap.dedent(inspect.getsource(value));path=Path(inspect.getsourcefile(inspect.unwrap(value)))
    return {'path':str(path),'sha256':f.digest(path),'source':source,
            'wholeBodyAST':hashlib.sha256(ast.dump(ast.parse(source)).encode()).hexdigest()}


def load(base):
    cp=json.loads(CP.read_text());origin=Path(cp['runDirectory']);validation=cp['validation'];vr=Path(validation['runDirectory'])
    f.require(f.digest(CP)==CP_SHA and cp['status']=='ACCEPTED_CASE_FREQUENCY_END_CALL_INPUTS','accepted exact full call inputs')
    h.source.receipts.inspect_guard(origin.parent,'frequency_end_calls');h.source.receipts.inspect_guard(vr,'validate')
    f.require(f.digest(origin/'checks.json')==cp['checksSha256'] and f.digest(vr/'checks.json')==validation['checksSha256'],'accepted checks hashes')
    f.require((vr/'checks.json').read_bytes()==(vr/'validate.stdout').read_bytes(),'independent checks/stdout identity')
    manifest={'runDirectory':str(base),'sourceFiles':dict(cp['sourceFiles']),'inputPackets':dict(cp['inputPackets']),
              'referencedInputs':{},'input':cp['input'],'settings':cp['settings'],
              'acceptedCallInputs':{'checkpoint':str(CP),'checkpointSha256':CP_SHA,'runDirectory':str(origin),
                  'checksSha256':cp['checksSha256'],'validationDirectory':str(vr),'validatorChecksSha256':validation['checksSha256']}}
    reference=lambda path,name,sha:h.source.reference(base,manifest,path,name,sha)
    for name,item in cp['artifacts'].items():reference(origin/name,name,item['sha256'])
    for name,value in cp['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==f.digest(origin/'source'/name)==value,'accepted current/frozen source')
        reference(origin/'source'/name,'source/'+name,value)
    for path,name in ((CP,'accepted-end-call-checkpoint.json'),(origin/'checks.json','accepted-end-call-checks.json'),
                      (origin/'inputs.json','accepted-end-call-manifest.json'),(vr/'checks.json','accepted-end-call-validation.json')):
        reference(path,name,f.digest(path))
    for name,item in validation['artifacts'].items():reference(vr/name,'accepted-call-validation/'+name,item['sha256'])
    modes=json.loads(MCP.read_text());mo=Path(modes['runDirectory'])
    f.require(modes['status']=='ACCEPTED_FOUR_CASE_MODE_SUBSPACES' and f.digest(mo/'checks.json')==modes['checksSha256'],'accepted full current/modal operations')
    reference(MCP,'coordinate-origins/modes-checkpoint.json',f.digest(MCP))
    reference(mo/'checks.json','coordinate-origins/modes-checks.json',modes['checksSha256'])
    for name,value in modes['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==f.digest(mo/'source'/name)==value,'native modal source pin')
        if name in manifest['sourceFiles']:f.require(manifest['sourceFiles'][name]==value,'shared exact modal source')
        else:
            manifest['sourceFiles'][name]=value;reference(mo/'source'/name,'source/'+name,value)
    for name in ('cases/LAB_HELD__RHOBR_CONSTANT/right/remainder-context/complete.pickle',
                 'cases/LAB_HELD__RHOBR_CONSTANT/right/pairing-derivatives.pickle',
                 'cases/LAB_HELD__RHOBR_CONSTANT/right/prepared-modal.pickle'):
        reference(mo/name,'accepted-coordinate-modal/'+name,modes['artifacts'][name]['sha256'])
    for path in (Path(__file__).resolve(),PLAN):
        name=str(path.relative_to(f.ROOT));value=f.digest(path);f.require(name not in manifest['sourceFiles'],'new coordinate source pin')
        manifest['sourceFiles'][name]=value;target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(path.read_bytes())
    for name,value in manifest['inputPackets'].items():f.require(f.digest(Path(name))==value,('original input prehash',name))
    f.save(base/'inputs.json',manifest)
    return cp,manifest


def joined_sources(base):
    native={'pairing':engine.ClosedCurrentPairing.construct,'modalPrepare':engine.ModalCurrentSubspaces.prepare,
            'acoustic':engine.ClosedAcousticEnergy.construct,'analytic':engine.FullPencilModes.analytic,
            'mapping':engine.ChannelInput.mapping,'endSources':h.source.q.end_sources}
    result={name:source_body(value) for name,value in native.items()}
    old=json.loads((base/'native-end-call-joins.json').read_text())
    for name in ('pairing','modalPrepare','analytic','mapping','endSources'):same(result[name],old[name])
    # Require the literal entire caller and its exact forward-coordinate
    # expressions, not merely a method name or inferred equivalence.
    pairing_tree=ast.parse(result['pairing']['source'])
    expressions={ast.unparse(n.value) for n in ast.walk(pairing_tree) if isinstance(n,ast.Assign)}
    f.require("algebraic.subs(self.modes.q, a.qlegs[1] / acoustic['ACOUSTIC_RADICAL_SCALE'])" in expressions,'actual native radical coordinate substitution')
    f.require('tuple((physical_q.xreplace(mapping) for mapping in leg_maps))' in expressions,'actual whole native forward pencil leg maps')
    result['coordinateRule']={'input':'saved PENCIL_PLUS / saved full tangent derivatives',
        'inverse':{'omegaRight':'omega','kRight':'native normal momentum','qRight':'actual acoustic scale * native radical'},
        'transportRule':'transform saved physical-q transport and divide by actual scale',
        'scope':'New coordinate images only; original analytic/mapping/derivative bodies are not called.'}
    f.save(base/'native-end-coordinate-joins.json',result)
    return result


def prohibit():
    h.prohibit()
    for module,names in ((previous,('load','prepare','main')),(h,('load','prepare','main')),
                         (h.source.q,('end_sources','sources','main')),(h.inputs.chart,('end_tables','source_chart','lift','main'))):
        for name in names:
            if hasattr(module,name):setattr(module,name,forbidden)


def forbidden(*args,**kwargs):raise RuntimeError('coordinate conversion cannot replay accepted scientific operations')


def atomic_operation(base,name,inputs,operation):
    path=base/'operations'/name;path.mkdir(parents=True)
    f.atomic_pickle(path/'input.pickle',inputs)
    value=operation()
    f.atomic_pickle(path/'value.pickle',value)
    f.save(path/'completed.json',{'inputSha256':f.digest(path/'input.pickle'),'valueSha256':f.digest(path/'value.pickle'),
           'operation':'one new exact frequency-leg coordinate image','nativeDerivativeOrBindingCalled':False})
    return value


def construct(base,cp):
    proof_summary=json.loads((base/'accepted-call-validation/accepted-source-proof-routing-summary.json').read_text())
    proof_by_address={tuple(v['address']):v for v in proof_summary['records']}
    inputs={};aliases={};baseline={}
    for end in ('LEFT','RIGHT'):baseline[end]=f.unpickle(base/'end-call-baseline'/end.lower()/'input-value.pickle')
    for label in cp['cases']:
        aliases[label]={}
        for end in ('LEFT','RIGHT'):
            address=(label,end);folder=base/'end-call-cases'/label/end.lower()
            raw=f.unpickle(base/'end-input-cases'/label/end.lower()/'full-inputs.pickle')
            ctx=f.unpickle(base/'end-input-cases'/label/'context-pairs.pickle')
            call=f.unpickle(folder/'call-input-pairs.pickle');legs=f.unpickle(folder/'saved-frequency-leg-operands.pickle')
            proof=proof_by_address[address]
            same(call['input']['sourceSignature'],raw['signature']);same(raw['address'],address)
            same(call['input']['seedInput'],ctx['seedInput']);same(call['input']['frame'],ctx['frame'])
            f.require(f.digest(base/'end-input-cases'/label/end.lower()/'full-inputs.pickle')==proof['actualSourceSha256'],'accepted actual source proof consumer hash')
            f.require(proof['otherNativeAnalyticInputsIdentical'] and proof['baselineMappingInputsIdentical'] and proof['relationMatchesBaseline'] and proof['branchReceiptMatchesBaseline'],'complete saved native source-call joins')
            for key in ('FREQUENCY_LEGS','NORMAL_LEGS','BULK_LEGS','CLOSED_PENCIL_LEGS','ACOUSTIC_WAVE_ROWS','SOURCE_BRANCH_JOINS'):
                same(legs['pairing'][key],raw['modeInput']['pairing'][key])
            route={'address':address,'sourceInput':call,'sourceProofReceipt':proof,'rawInput':raw,'context':ctx,'frequencyLegOperands':legs}
            inputs[address]=route
            if proof['acceptedSemanticCurrentInputReuse']:
                selected=baseline[end]
                same(raw['modeInput']['relation'],selected['result']['analytic'][1]);same(raw['modeInput']['acoustic']['SOURCE_BRANCH_JOINS'],selected['result']['analytic'][2])
                alias={'address':address,'mode':'accepted-baseline-source-proof-route','owner':selected['owner'],
                    'baselineResultPath':str(base/'end-call-baseline'/end.lower()/'input-value.pickle'),
                    'fullInputPath':str(base/'end-call-cases'/label/end.lower()/'call-input-pairs.pickle'),
                    'proofReceipt':proof,'numericalEndReuseAccepted':False}
            else:
                f.require(end=='RIGHT' and label in ('LAB_HELD__RHOBR_CONSTANT','MATERIAL_ADVECTED__RHOBR_CONSTANT'),'actual distinct saved current family')
                alias={'address':address,'mode':'missing-coordinate-image','owner':OWNER,'proofReceipt':proof,'numericalEndReuseAccepted':False}
            directory=base/'coordinate-cases'/label/end.lower();directory.mkdir(parents=True)
            f.atomic_pickle(directory/'full-input-route.pickle',route)
            f.save(directory/'source-route.json',alias);aliases[label][end]=alias
    own=inputs[OWNER];shared=inputs[('MATERIAL_ADVECTED__RHOBR_CONSTANT','RIGHT')]
    # Addresses differ. Each complete physical value and consumed binder/context
    # field must match independently of the metadata provenance strings.
    for key in ('pairing','symbolic','scalar','binding','knownUnits'):same(own['frequencyLegOperands'][key],shared['frequencyLegOperands'][key])
    for key in ('analyticInput','origin','frequency','seedInput','frame','context'):same(own['sourceInput']['input'][key],shared['sourceInput']['input'][key])
    legs=own['frequencyLegOperands'];raw=own['rawInput'];ctx=own['context'];call=own['sourceInput']['input']['analyticInput']
    omega,k,q=(call[n] for n in ('omega','momentum','radical'))
    wr,kr,qr=(legs['pairing'][name][1] for name in ('FREQUENCY_LEGS','NORMAL_LEGS','BULK_LEGS'))
    scale=legs['scalar']['RADICAL_SCALE'];known=legs['knownUnits']
    f.require(scale.is_Pow and scale.exp==-1 and scale.base in known,'actual saved inverse sound-speed scale')
    f.require(tuple(known[scale.base])==(1,-1,0) and tuple(known[qr])==(-1,0,0),'actual sound-speed and bulk-momentum dimensions')
    f.require(tuple(call['units']['radical'])==(0,-1,0) and tuple(known[wr])==(0,-1,0) and tuple(known[kr])==(-1,0,0),'actual native coordinate units')
    f.require(not scale.has(wr,kr,qr,omega,k,q),'actual coordinate scale independent of differentiated variables')
    same(scale,raw['modeInput']['acoustic']['ACOUSTIC_RADICAL_SCALE'])
    f.require(raw['modeInput']['acoustic']['ACOUSTIC_RADICAL_JOIN_RESIDUAL']==0,'saved native acoustic radical proof')
    inverse={wr:omega,kr:k,qr:scale*q};forward={omega:wr,k:kr,q:qr/scale}
    coordinate_input={'owner':OWNER,'inverse':inverse,'forward':forward,'scale':scale,'knownUnits':known,
        'rawInput':raw,'context':ctx,'callInput':own['sourceInput'],'operands':legs,
        'acceptedAcousticRelation':raw['modeInput']['relation'],'acceptedAcousticResidual':raw['modeInput']['acoustic']['ACOUSTIC_RADICAL_JOIN_RESIDUAL']}
    f.atomic_pickle(base/'coordinate-map-input.pickle',coordinate_input)
    results={};roundtrips={}
    for name in ('PENCIL_PLUS','FREQUENCY_PENCIL_PLUS','NORMAL_PENCIL_PLUS','FREQUENCY_TRANSPORT','NORMAL_TRANSPORT'):
        transport=name.endswith('_TRANSPORT');source=legs['scalar' if transport else 'symbolic'][name]
        operation_input={'owner':OWNER,'name':name,'source':source,'inverse':inverse,'forward':forward,'scale':scale,'transport':transport}
        value=atomic_operation(base,name,operation_input,lambda source=source,transport=transport:source.xreplace(inverse)/scale if transport else source.xreplace(inverse))
        # Save the actual reversed operand and raw equation before its guard.
        reverse=(value*scale if transport else value).xreplace(forward)
        residual=reverse-source
        proof={'source':source,'image':value,'reverse':reverse,'rawResidual':residual,'inverse':inverse,'forward':forward}
        f.atomic_pickle(base/'operations'/name/'roundtrip.pickle',proof)
        values=list(residual) if isinstance(residual,sp.MatrixBase) else [residual]
        f.require(all(v==0 for v in values),('exact new coordinate roundtrip',name))
        f.require(not value.has(wr,kr,qr),'new native-coordinate image contains no old plus-leg variables')
        f.save(base/'operations'/name/'checks.json',{'roundtripResidualScalars':len(values),'allZero':True,'newDerivativeOrBindingCalls':0})
        results[name]=value;roundtrips[name]=len(values)
    # The original native mapping call already completed in derivative_packet.
    # Read its complete return, full current source and preserved known units.
    sourcepath=base/'accepted-coordinate-modal/cases/LAB_HELD__RHOBR_CONSTANT/right/remainder-context/complete.pickle'
    remainder,remknown=f.unpickle(sourcepath);same(remainder['result'],raw['modeInput']['pairing'])
    mapping=remainder['bindings'];relation=raw['modeInput']['relation'];live=(omega,k,q,*ctx['seedInput']['origin'])
    atoms=(results['PENCIL_PLUS'].free_symbols|relation.free_symbols)-set(live)
    limits=results['PENCIL_PLUS'].atoms(sp.Limit)
    requested={'algebraic':results['PENCIL_PLUS'],'relation':relation,'live':live,'parameters':ctx['seedInput']['parameters'],
               'limits':ctx['seedInput']['limits'],'savedMapping':mapping,'actualAtoms':atoms,'actualLimits':limits,
               'originalMappingPacketPath':str(sourcepath),'originalMappingSourceSha256':f.digest(sourcepath)}
    f.atomic_pickle(base/'saved-native-mapping-input-return.pickle',requested)
    f.require(set(mapping)==atoms|limits,'complete native consumed mapping return fieldset')
    for atom in atoms:f.require(atom.name in requested['parameters'] and mapping[atom]==requested['parameters'][atom.name],'actual original parameter value')
    for limit in limits:f.require(limit in requested['limits'] and mapping[limit]==requested['limits'][limit],'actual original ordered-limit value')
    f.require(not any(atom in mapping for atom in live),'native frequency/momentum/independent grades remain live')
    # New route controls alter actual coordinate/source/owner/unit/mapping inputs,
    # never a scientific result. Preserve all changed operands before checks.
    changed_inverse=dict(inverse);changed_inverse[qr]=2*inverse[qr]
    changed_owner=(OWNER[0],'LEFT');changed_unit=dict(call['units']);changed_unit['radical']=(0,0,0)
    atom=next(iter(mapping));changed_mapping=dict(mapping);changed_mapping[atom]=mapping[atom]+1
    controls={'actualInverse':inverse,'changedInverse':changed_inverse,'actualOwner':OWNER,'changedOwner':changed_owner,
              'actualUnits':call['units'],'changedUnits':changed_unit,'actualMapping':mapping,'changedMapping':changed_mapping}
    f.atomic_pickle(base/'route-mutation-operands.pickle',controls)
    receipt={'coordinateCoefficient':not previous.modes_reader.same(inverse,changed_inverse),'owner':OWNER!=changed_owner,
             'unit':not previous.modes_reader.same(call['units'],changed_unit),'mappingCoefficient':not previous.modes_reader.same(mapping,changed_mapping)}
    f.require(all(receipt.values()),'actual coordinate/caller route controls respond');f.save(base/'route-mutation-controls.json',receipt)
    output={'owner':OWNER,'coordinateImages':results,'mapping':mapping,'originalRelation':relation,
            'branchResiduals':raw['modeInput']['acoustic']['SOURCE_BRANCH_JOINS'],'sourceInput':own,
            'sharedInput':shared,'roundtripScalarCounts':roundtrips,'aliases':aliases,
            'originalAnalyticReturnReconstructed':False,'nativeCoordinateImageOnly':True,
            'scope':'Missing exact coordinate images from accepted live-frequency current values; binding to the live end parameter and native end domains/tables remain unfinished.'}
    f.atomic_pickle(base/'remaining-case-end-coordinate-images.pickle',output)
    f.save(base/'coordinate-case-inventory.json',aliases)
    return {'physicalEnds':8,'baselineSourceAliases':sum(v['mode']=='accepted-baseline-source-proof-route' for c in aliases.values() for v in c.values()),
            'newCoordinateOwners':1,'sharedNewCoordinateUses':2,'newCoordinateImages':len(results),'roundtripResidualScalars':sum(roundtrips.values()),
            'savedNativeMappingReused':True,'actualMutationControls':len(receipt),'newDerivativeBindingAnalyticCurrentModeCalls':0,
            'frequencyEndOperatorOrDomainAccepted':False,'cases':aliases}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);args=ap.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);start=time.monotonic()
    previous.protect(base);cp,manifest=load(base);joins=joined_sources(base);prohibit()
    result=construct(base,cp)
    for name,value in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/name)==f.digest(base/'source'/name)==value,'current/frozen coordinate source posthash')
    for name,value in manifest['inputPackets'].items():f.require(f.digest(Path(name))==value,'original coordinate input posthash')
    for name,item in manifest['referencedInputs'].items():
        path=base/name;f.require(path.is_symlink() and str(path.readlink())==item['original'] and str(path.resolve())==item['resolvedOriginal'] and path.stat().st_size==item['bytes'] and f.digest(path)==item['sha256'],'coordinate reference postidentity')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    checks={**manifest,**result,'status':'COMPLETED_CASE_FREQUENCY_END_COORDINATE_IMAGES','nativeJoins':joins,'artifacts':artifacts,
            'wallSeconds':time.monotonic()-start,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
