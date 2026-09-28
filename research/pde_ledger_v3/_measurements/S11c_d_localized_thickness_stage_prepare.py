#!/usr/bin/env python3
"""Pin bounded review corrections and a pending guarded stage; never run science."""
import ast
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import pickle
import shutil
import subprocess

ROOT=Path('/var/projects/toy_physics');M=ROOT/'research/pde_ledger_v3/_measurements'
BASE=ROOT/'_scratch/s11c/s11c-d-localized-thickness-20260928'
REVIEW=BASE/'build-review';P='S11c_d_localized_thickness_'
WORKER=M/(P+'response.py');LAUNCHER=M/(P+'stage_launch.py')


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def route(p):
 p=Path(p);return {'path':str(p),'canonicalPath':str(p.resolve(strict=True)),'bytes':p.stat().st_size,'sha256':sha(p)}
def save(p,x):
 with p.open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')
def named(tree,name):
 found=[n for n in ast.walk(tree) if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name==name]
 assert len(found)==1,(name,len(found));return found[0]


def main():
 out=BASE/'runtime-preparation';out.mkdir(exist_ok=False)
 state=json.loads((REVIEW/'state.json').read_text());packet=Path(state['packetDirectory'])
 assert state['status']=='AWAITING_ADJUDICATION'
 for name,h in state['fileHashes'].items():assert sha(packet/name)==h,name
 assert sha(REVIEW/'packet.tar.gz')==state['archiveSha256']
 reviews={}
 for engine,key in [('claude','result'),('grok','text')]:
  receipt=json.loads((REVIEW/(engine+'-run.json')).read_text());result=json.loads((REVIEW/(engine+'.json')).read_text())
  assert receipt['exitStatus']==0 and receipt['finishedUtc'] and receipt['outputSha256']==sha(REVIEW/(engine+'.json'))
  assert receipt['stderrSha256']==sha(REVIEW/(engine+'.stderr'))
  assert (M/(P+engine+'_review.md')).read_text()==result[key]
  assert 'NEEDS REVISION' in result[key]
  reviews[engine]={'literalVerdict':'NEEDS REVISION','literalReport':route(M/(P+engine+'_review.md')),
    'rawOutput':route(REVIEW/(engine+'.json')),'runReceipt':route(REVIEW/(engine+'-run.json')),
    'stderr':route(REVIEW/(engine+'.stderr')),'finishedUtc':receipt['finishedUtc'],'wallSeconds':receipt['wallSeconds']}
 assert reviews['claude']['finishedUtc']<=json.loads((REVIEW/'grok-run.json').read_text())['startedUtc']
 old=ast.parse((REVIEW/'reviewed-source'/WORKER.name).read_text());new=ast.parse(WORKER.read_text())
 compile(new,str(WORKER),'exec');compile(LAUNCHER.read_text(),str(LAUNCHER),'exec')
 unchanged=('require','digest','route','save','containment','SavedCodec','Journal',
            'json_form','basis_polynomial','base_transform','profile_transform','amplitude_domain',
            'certificate','numerical_matrix','norm','equation_probe','profile_probe','main')
 for name in unchanged:assert ast.dump(named(old,name))==ast.dump(named(new,name)),name
 forbidden={'lu_solve','LUsolve','inv','inverse','solve','nsolve','roots','nroots','series','diff','integrate','doit','limit'}
 calls={n.func.attr if isinstance(n.func,ast.Attribute) else n.func.id for n in ast.walk(new)
        if isinstance(n,ast.Call) and isinstance(n.func,(ast.Attribute,ast.Name))}
 assert not calls&forbidden,calls&forbidden
 assert 'sp.exp(-sp.expand(phase))' in ast.unparse(named(new,'remove_phase'))
 assert not any(isinstance(n,ast.ImportFrom) and n.module and 'S11c' in n.module for n in ast.walk(new))
 # Regression fixtures use Python integers/matrices, not SymPy or actual source data.
 class Matrix:
  def __init__(self,n=5,data=None):self.data=[r[:] for r in data] if data is not None else [[0]*n for _ in range(n)]
  def __getitem__(self,k):return self.data[k[0]][k[1]]
  def __setitem__(self,k,v):self.data[k[0]][k[1]]=v
  def __mul__(self,s):return Matrix(data=[[v*s for v in row] for row in self.data])
  def __eq__(self,other):return isinstance(other,Matrix) and self.data==other.data
 class Stub:
  @staticmethod
  def zeros(n):return Matrix(n)
  @staticmethod
  def ImmutableMatrix(v):return Matrix(data=v.data)
 ns={'sp':Stub,'base_transform':lambda _:1}
 for name in ('assemble_transformed_source','unit_contractions'):
  exec(compile(ast.Module(body=[named(new,name)],type_ignores=[]),'<stdlib-fixture>','exec'),ns)
 local=[{'row':0,'column':0,'amplitudeWithoutB0':7}]
 address={'row':0,'sourceColumn':2,'frameColumn':0,'term':4}
 native=[{'address':address,'amplitudeWithoutB0':11}]
 complete=ns['assemble_transformed_source'](local,native,0)
 omitted=ns['assemble_transformed_source'](local,native,0,address)
 assert complete['forcing'][0,0]==18 and omitted['forcing'][0,0]==7
 broken_baseline=Matrix();broken_baseline[0,0]=7
 assert complete['forcing']!=broken_baseline,'full rebuild must reject a dropped baseline contribution'
 duplicate=ns['assemble_transformed_source'](local,native*2,0,address)
 assert not duplicate['checks']['exactAddressOmissionCount']
 units=[[(0,0,0) for _ in range(5)] for _ in range(5)]
 assert ns['unit_contractions'](units,units,units)['checks']['allUnitContractionsMatch']
 bad=[r[:] for r in units];bad[0][0]=(1,0,0)
 fixture_result=ns['unit_contractions'](units,units,bad)
 assert len(fixture_result['contractions'])==125 and not fixture_result['checks']['allUnitContractionsMatch']
 # Check that this failure's return and JSON persist before checked raises.
 fixture=out/'unit-persistence-fixture';fixture.mkdir()
 harness={'Path':Path,'hashlib':hashlib,'json':json,'os':os,'pickle':pickle,'shutil':shutil,
          'datetime':datetime,'timezone':timezone,'out':fixture,'json_form':lambda v:v}
 for name in ('require','digest','route','save','SavedCodec','Journal','checked'):
  exec(compile(ast.Module(body=[named(new,name)],type_ignores=[]),'<stdlib-persistence>','exec'),harness)
 harness['j']=harness['Journal'](fixture)
 try:harness['checked']('bad-synthetic-units',lambda value:value,fixture_result)
 except ValueError:pass
 else:raise AssertionError('bad synthetic units accepted')
 assert (fixture/'bad-synthetic-units.json').exists()
 assert (fixture/'operations/0000-bad-synthetic-units/value.pickle').exists()
 assert len(harness['j'].records)==1 and not harness['j'].stack
 # No fixture pickle was restored; data above are deliberately non-scientific.
 manifest=M/(P+'corrected_inputs.json');spec=json.loads((M/(P+'inputs.json')).read_text())
 spec['status']='PINNED_CORRECTED_INPUTS_PENDING_ONE_STAGE_APPROVAL'
 spec['inputs']['methodDescription']=route(ROOT/'research/pde_ledger_v3/directives/S11c_d_localized_thickness_implementation.md')
 for name,p in [('reviewDisposition',M/(P+'review_disposition.md')),('phaseSourceEvidence',M/(P+'phase_source_evidence.json')),
  ('claudeLiteralReview',M/(P+'claude_review.md')),('grokLiteralReview',M/(P+'grok_review.md')),
  ('completionHook',ROOT/'scripts/codex_job_watch.py'),('stageLauncher',LAUNCHER),
  ('completionMessage',M/(P+'completion_message.txt'))]:spec['inputs'][name]=route(p)
 for name in spec['scientificPayloads']:
  assert route(spec['inputs'][name]['path'])==spec['inputs'][name],name
 for name,r in spec['inputs'].items():assert route(r['path'])==r,name
 save(manifest,spec)
 static={'status':'PASS_SOURCE_AND_STDLIB_INSTRUMENT_CHECKS_NO_SCIENCE',
  'unchangedHelperASTs':list(unchanged),'phaseExpansionPresent':True,
  'actualAddressReassemblyFixture':True,'droppedBaselineContributionDetected':True,
  'duplicateOmissionRefused':True,'all125UnitFixtureRecordsReturned':True,
  'unitFailureReturnAndJsonPersistBeforeGuard':True,'scientificImportsOrRestorations':0,
  'workerSha256':sha(WORKER),'inputManifestSha256':sha(manifest),
  'launcherSha256':sha(LAUNCHER),'sourceInputs':len(spec['inputs'])}
 adjudication=M/(P+'review_adjudication.json')
 save(adjudication,{'status':'BOUNDED_FINDINGS_LOCALLY_DISPOSED_PENDING_STAGE_APPROVAL',
  'recordedUtc':datetime.now(timezone.utc).isoformat(),'packetSha256':state['packetSha256'],
  'archiveSha256':state['archiveSha256'],'reviewers':reviews,
  'originalReviewedWorkerSha256':sha(REVIEW/'reviewed-source'/WORKER.name),
  'correctedWorkerSha256':sha(WORKER),'correctedInputManifestSha256':sha(manifest),
  'reviewDisposition':route(M/(P+'review_disposition.md')),
  'phaseCitationEvidence':route(M/(P+'phase_source_evidence.json')),
  'findings':{'Claude S1':'ADOPTED_EXPLICIT_PHASE_EXPANSION',
   'Claude S2':'ADOPTED_FULL_BASELINE_REBUILD_AND_ACTUAL_ADDRESS_OMISSION',
   'Claude S3':'ADOPTED_COMPLETE_UNIT_RECORDS_PERSISTED_BEFORE_GUARD',
   'Grok G1':'CITATION_DISPROVED_BY_OWN_CASE_PACKET_EVIDENCE_EXPLICIT_RUNTIME_PHASE_JOIN_ADOPTED'},
  'substantiveFindingsClosed':True,'freshIndependentClearClaimed':False,
  'reviewRerunOrPeerSharing':False,'scientificStageAuthorized':False,
  'scientificObjectsRestored':False,'staticCorrectionChecks':static,
  'reviewStderrDisposition':'Claude empty; Grok retained configuration warnings and one read_file tool error with completed literal report.',
  'completionEventHandled':True,'protectedBuilderSuffixSha256':'f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2'})
 proposal=M/(P+'proposal.json')
 save(proposal,{'status':'PENDING_EXPLICIT_APPROVAL_FOR_ONE_GUARDED_LOCALIZED_THICKNESS_RESPONSE',
  'workerSha256':sha(WORKER),'inputManifestSha256':sha(manifest),
  'launcherSha256':sha(LAUNCHER),'completionMessageSha256':sha(M/(P+'completion_message.txt')),
  'sharedGuardSha256':sha(ROOT/'scripts/s11c_guarded_run.py'),
  'supervisorSha256':sha(M/'S11c_d_end_normalization_run.py'),
  'completionHookSha256':sha(ROOT/'scripts/codex_job_watch.py'),
  'adjudicationPath':str(adjudication),'adjudicationSha256':sha(adjudication),
  'reviewDispositionPath':str(M/(P+'review_disposition.md')),'reviewDispositionSha256':sha(M/(P+'review_disposition.md')),
  'scopeExplicitlyApproved':False,'substantiveReviewFindingsClosed':True,
  'scope':'One fixed-input two-case thickness forcing transform and full saved coupled PV/delta action candidate; original operands/physical input unchanged.',
  'seconds':900,'nativeSeconds':840,'memoryGiB':2,'swapMax':0,'cpuCount':1,'nice':15,'tasksMax':32,'nativeThreads':1,
  'automaticRetry':False,'resultDirectory':str(BASE/'production/complete'),
  'expectedLimitedStop':'FIXED_INPUT_THICKNESS_OUTGOING_ACTION_CANDIDATE_BUILT',
  'resultAcceptancePendingSavedOutputInspection':True,'reviewVerdicts':{'claude':'NEEDS REVISION','grok':'NEEDS REVISION'},
  'freshIndependentClearClaimed':False,'currentThread':'01a0e01b-ef84-7192-817f-584cda5d339b'})
 result=subprocess.run(['python3',str(WORKER),'--input-manifest',str(manifest),'--gate-receipt',str(proposal),
    '--run-directory',str(BASE/'production/complete')],capture_output=True,text=True)
 assert result.returncode and 'method gate incomplete' in result.stderr and not (BASE/'production').exists()
 (out/'pending-gate.stdout').write_text(result.stdout);(out/'pending-gate.stderr').write_text(result.stderr)
 result=subprocess.run(['python3',str(LAUNCHER)],capture_output=True,text=True)
 assert result.returncode and 'S11c_d_localized_thickness_authorization.json' in result.stderr
 assert not (M/(P+'gate.json')).exists() and not (BASE/'production').exists()
 (out/'missing-stage-approval.stdout').write_text(result.stdout);(out/'missing-stage-approval.stderr').write_text(result.stderr)
 static.update(pendingWorkerGateRefused=True,missingStageApprovalRefused=True,noRuntimeDirectoryOrReadyGateCreated=True)
 save(M/(P+'corrected_static_checks.json'),static)
 save(M/(P+'stage_preparation.json'),{'status':'CORRECTED_STAGE_PREPARED_PENDING_USER_APPROVAL',
  'recordedUtc':datetime.now(timezone.utc).isoformat(),'worker':route(WORKER),'inputManifest':route(manifest),
  'proposal':route(proposal),'launcher':route(LAUNCHER),'adjudication':route(adjudication),
  'staticChecks':route(M/(P+'corrected_static_checks.json')),'sourceInputs':len(spec['inputs']),
  'scientificPayloads':len(spec['scientificPayloads']),'productionRoot':str(BASE/'production'),
  'scientificLaunchAuthorized':False,'scientificOperationsDuringPreparation':0,'newExternalReviews':0,
  'currentThread':'01a0e01b-ef84-7192-817f-584cda5d339b'})
 print(json.dumps({'status':'CORRECTED_STAGE_PREPARED_PENDING_USER_APPROVAL','sourceInputs':len(spec['inputs']),
 'scientificPayloads':len(spec['scientificPayloads']),'correctedWorkerSha256':sha(WORKER),'readyGateExists':False,
 'scienceRun':False,'reviewRerun':False}))


if __name__=='__main__':main()
