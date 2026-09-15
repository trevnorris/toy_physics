#!/usr/bin/env python3
"""Finish the saved transcript comparison without repeating validated work."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import shutil
import time

import S11c_d_source_fourier_quadrature_recover as recovery
from S11c_d_source_fourier_quadrature_recover import ROOT,STORE,engine,digest,save,unpickle,equal

CHECKPOINT=ROOT/'_measurements/S11c_d_source_fourier_quadrature_metadata_repair.json'
PLAN=ROOT/'_measurements/S11c_d_source_fourier_quadrature_finish_plan.md'


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--source-run-directory',type=Path,required=True)
    parser.add_argument('--run-directory',type=Path,required=True);args=parser.parse_args()
    previous=args.source_run_directory.resolve();previous.relative_to(STORE)
    base=args.run_directory.resolve();base.relative_to(STORE);base.mkdir(parents=True,exist_ok=False)
    begin=time.monotonic();checkpoint=json.loads(CHECKPOINT.read_text())
    equal(str(previous),checkpoint['sourceRunDirectory'])
    for name,item in checkpoint['sourceArtifacts'].items():equal(digest(previous/name),item['sha256'],('pinned-artifact',name))
    checks=json.loads((previous/'checks.json').read_text());preflight=json.loads((previous/'preflight.json').read_text())
    original=Path(preflight['sourceRunDirectory']);original.relative_to(STORE)
    invocation=json.loads((previous.parent/'source_recover.invocation.json').read_text())
    equal(invocation,checkpoint['sourceInvocation'])
    # The frozen program reached only its final difference-census guard. Its
    # original guards, full emission replay and inventories completed first.
    error=(previous.parent/'source_recover.stderr').read_text()
    if invocation['exitCode']!=1 or 'unexpected recovery emission differences' not in error:
        raise ValueError('unexpected predecessor outcome')
    repair_name=str(Path(recovery.__file__).relative_to(ROOT))
    for name,sha in preflight['sourceFiles'].items():
        equal(digest(previous/'source'/name),sha,('frozen-source',name))
        if name!=repair_name:equal(digest(ROOT/name),sha,('unchanged-helper',name))
    old=ast.parse((previous/'source'/repair_name).read_text());new=ast.parse(Path(recovery.__file__).read_text())
    old_members={n.name:ast.dump(n) for n in old.body if isinstance(n,ast.FunctionDef)}
    new_members={n.name:ast.dump(n) for n in new.body if isinstance(n,ast.FunctionDef)}
    for name in old_members:
        if name!='compare_transcripts':equal(old_members[name],new_members[name],('unchanged-recovery-function',name))
    equal(digest(Path(recovery.__file__)),checkpoint['repairedRecoverySha256'])
    equal(hashlib.sha256(old_members['compare_transcripts'].encode()).hexdigest(),checkpoint['oldComparisonAstSha256'])
    equal(hashlib.sha256(new_members['compare_transcripts'].encode()).hexdigest(),checkpoint['newComparisonAstSha256'])
    pins={name:digest(ROOT/name) for name in preflight['sourceFiles']}
    for path in (Path(__file__).resolve(),CHECKPOINT,PLAN):pins[str(path.relative_to(ROOT))]=digest(path)
    for name in pins:
        destination=base/'source'/name;destination.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/name,destination)
    save(base/'preflight.json',{'sourceFiles':pins,'sourceRunDirectory':str(previous),'sourceChecksSha256':digest(previous/'checks.json'),
        'originalNumericalRunDirectory':str(original),'sourceArtifacts':checkpoint['sourceArtifacts'],
        'unchangedPrecomparisonGuards':True,'provenance':checks['provenance']})
    for name,item in checks['artifacts'].items():
        equal(digest(previous/name),item['sha256'],('saved-output',name));shutil.copyfile(previous/name,base/name)
    for collection in ('integralArtifacts','boundRecordArtifacts'):
        for item in checks[collection]:
            equal(digest(previous/item['path']),item['sha256'],('saved-record',item['path']))
            destination=base/item['path'];destination.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(previous/item['path'],destination)
    replay_inventory=json.loads((previous/'source-binding-replay-inventory.json').read_text())
    expected_keys={(test,source,kind) for test in range(2) for source in range(35) for kind in ('boundAmplitude','boundSource')}
    if len(replay_inventory)!=len(expected_keys) or {tuple(i['key']) for i in replay_inventory}!=expected_keys:
        raise ValueError('source replay pair coverage')
    proofs=[]
    for item in replay_inventory:
        directory=previous/item['path'];target=base/item['path'];target.mkdir(parents=True)
        for name,artifact in item['artifacts'].items():
            equal(digest(directory/name),artifact['sha256'],('replay-artifact',item['key'],name));shutil.copyfile(directory/name,target/name)
        operands=unpickle(directory/'operands.pickle')
        equal(operands['sourceBoundPacketSha256'],checks['boundPacketSha256BeforeEmission'])
        equal(digest(directory/'operands.pickle'),item['operandSha256'])
        for value,text in zip((operands['left'],operands['right']),operands['representations']):equal(value,recovery._restore(text))
        if item['liveExactEqual']:
            equal(operands['left'],operands['right']);equal(operands['rawResidual'],0)
            if not item['normalizedResidualZero']:raise ValueError('exact replay guard')
        else:
            worker=json.loads((directory/'checks.json').read_text())
            equal(worker,item['checks']);equal(json.loads((directory/'stdout').read_text()),worker)
            if item['exitCode']!=0 or item['stderrBytes'] or (directory/'stderr').stat().st_size:
                raise ValueError('source replay worker outcome')
            equal(worker['representationSha256'],[hashlib.sha256(s.encode()).hexdigest() for s in operands['representations']])
            cert=unpickle(directory/'certificate.pickle');mutation=unpickle(directory/'mutation.pickle')
            equal(cert['LEFT'],operands['left']);equal(cert['RIGHT'],operands['right']);equal(cert['RESIDUAL'],0)
            equal(mutation['LEFT'],cert['LEFT']);equal(mutation['RIGHT'],2*cert['RIGHT'])
            if mutation['RESIDUAL']==0 or not worker['coefficientMutationNonzero']:raise ValueError('coefficient mutation')
            values=tuple(cert['REPLAY_RESIDUALS'])+tuple(v[1] for v in cert['PHASE_SPLITS'].values())+tuple(v[2] for v in cert['RADICAL_POWERS'].values())
            equal(len(values),worker['proofResidualCount'])
            if any(v!=0 for v in values) or worker['nonzeroProofResiduals']:raise ValueError('source replay proof')
            proofs.extend(values)
    for name in ('source-binding-replay-inventory.json','integral-inventory.json','bound-record-inventory.json','representation-joins.json'):
        shutil.copyfile(previous/name,base/name)
    bound_packet=unpickle(base/'bound-sources.pickle');quadrature=unpickle(base/'quadrature.pickle')
    bound=bound_packet['result'];representations=json.loads((base/'representation-joins.json').read_text())
    equal(bound_packet['provenance'],checks['provenance']);equal(quadrature['provenance'],checks['provenance'])
    equal(quadrature['boundPacketSha256'],digest(base/'bound-sources.pickle'))
    source,source_base=recovery.original.accepted(recovery.original.SOURCE)
    r,dimensions=recovery.original.native.source.restore_context(unpickle(source_base/'reduced-action.pickle'))
    dimensions.__dict__.update(unpickle(base/'dimensions-after-emission.pickle'))
    comparison=recovery.compare_transcripts(original,base,bound,representations)
    # Test that the narrow zero-support proof rejects changes of units or
    # nonzero support, rather than merely permitting these twenty tag names.
    old_entries=recovery.entries(original/'full.out');new_entries=recovery.entries(base/'full.out')
    controls=[]
    for record in comparison['certifiedRawZeroMetadataTransitions']:
        item=bound['nativeTestIntegralComparisons'][record['occurrence']]
        old_body=recovery._restore(old_entries[record['tag']]);new_body=recovery._restore(new_entries[record['tag']])
        for field in ('DIMENSION_L_T_M','MULTIGRADE'):
            changed={str(k):v for k,v in old_body[0]}
            v=changed[field]
            changed[field]=(v[0]+1,*v[1:]) if field=='DIMENSION_L_T_M' else ()
            try:recovery.zero_metadata_transition(engine.cas([changed]),new_body,item)
            except ValueError:controls.append({'occurrence':record['occurrence'],'field':field,'rejected':True})
            else:raise ValueError(('metadata mutation accepted',field))
    save(base/'metadata-mutation-controls.json',controls)
    if dimensions.constraints:raise ValueError('unresolved dimensions')
    if not all(checks['originalLimitJoins']) or not all(checks['nativeTestIntegralJoins']) or checks['nonzeroNormalizedBindingResiduals'] or checks['nonzeroBindingProofResiduals']:
        raise ValueError('prior completed exact guards')
    if (checks['maxAssignmentResidual']>1e-12 or checks['nonzeroAffineResiduals'] or
            any(n['maxScaledSourceResidual']>1e-10 or n['maxScaledAdaptiveResidual']>1e-9 or n['maxMutationDifference']<=1e-12 for n in checks['norms'])):
        raise ValueError('prior completed numerical guards')
    for name in ('bound-sources.pickle','quadrature.pickle','full.out','dimensions-after-emission.pickle'):
        equal(digest(base/name),checks['artifacts'][name]['sha256'],('unchanged-accepted-artifact',name))
    for name,item in checkpoint['sourceArtifacts'].items():equal(digest(previous/name),item['sha256'],('unchanged-predecessor',name))
    equal(pins,{name:digest(ROOT/name) for name in pins},('post-completion-sources',))
    result=checks|{'runDirectory':str(base),'sourceFiles':pins,'sourceRunDirectory':str(previous),
        'originalNumericalRunDirectory':str(original),'sourceChecksSha256':digest(previous/'checks.json'),
        'sourceTranscriptSha256':digest(original/'full.out'),'sourceRecoveryTranscriptSha256':digest(previous/'full.out'),
        'unchangedPrecomparisonGuards':True,'savedOperandAndInventoryJoins':True,'quadratureRecomputed':False,
        'sourceReplayOrEmissionRepeated':False,'emissionComparison':comparison,
        'restoredExactBindingOccurrences':sum(v['restoredExactEqual'] for v in representations),
        'nonzeroLiveRawBindingRepresentations':sum(v['comparison']['representationStrings'][2]!=recovery.sp.srepr(recovery.sp.S.Zero) for v in bound['nativeTestIntegralComparisons']),
        'sourceBindingReplayArtifacts':replay_inventory,'sourceBindingReplayProofCount':len(proofs),
        'sourceBindingReplayNonzeroProofCount':sum(v!=0 for v in proofs),'metadataMutationControls':controls,
        'sourceBindingReplayInventorySha256':digest(base/'source-binding-replay-inventory.json'),
        'representationJoinsSha256':digest(base/'representation-joins.json'),'emissionDifferencesSha256':digest(base/'emission-differences.json'),
        'finishWallSeconds':time.monotonic()-begin,'status':'VALIDATED_SAVED_SOURCE_QUADRATURE'}
    save(base/'checks.json',result);save(base/'recovery.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
