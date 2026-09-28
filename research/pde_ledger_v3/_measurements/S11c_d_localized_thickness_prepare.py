#!/usr/bin/env python3
"""Build the fixed review packet using source/JSON/hash inspection only.

No scientific imports or execution. No external submission. Refuse to replace
an existing review root. Runtime copies and archives remain ignored scratch;
canonical workers, launcher, prompt, manifest and concise receipts live here.
"""
import ast
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

ROOT=Path('/var/projects/toy_physics')
M=ROOT/'research/pde_ledger_v3/_measurements'
D=ROOT/'research/pde_ledger_v3/directives'
RUN=ROOT/'_scratch/s11c/s11c-d-localized-thickness-20260928/build-review'
THREAD='01a0e01b-ef84-7192-817f-584cda5d339b'
PREFIX='S11c_d_localized_thickness_'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def route(p):
 p=Path(p)
 return {'path':str(p),'canonicalPath':str(p.resolve(strict=True)),
         'bytes':p.stat().st_size,'sha256':sha(p)}

def save(p,x):
 with p.open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')


def main():
 RUN.mkdir(parents=True,exist_ok=False)
 packet=RUN/'packet';packet.mkdir()
 parent=ROOT/'_scratch/s11c/s11c-d-two-asymptote-20260927/forcing-continuation/complete'
 inspection=ROOT/'_scratch/s11c/s11c-d-two-asymptote-20260927/saved-inspection/complete'
 reference=ROOT/'_scratch/s11c/s11c-d-reference-kernel-20260926/production/complete'
 outgoing=ROOT/'_scratch/s11c/s11c-d-outgoing-prescription-20260927/unlimited/complete'
 parent_index=json.loads((parent/'artifact-index.json').read_text())
 inspection_cp=json.loads((M/'S11c_d_two_asymptote_saved_inspection_checkpoint.json').read_text())
 kernel_cp=json.loads((M/'S11c_d_reference_kernel_checkpoint.json').read_text())
 outgoing_cp=json.loads((M/'S11c_d_outgoing_prescription_checkpoint.json').read_text())
 scientific={'candidate':outgoing/'outgoing-prescription-candidate.pickle',
             'context':reference/'outgoing-source-context.pickle',
             'units':reference/'physical-units.pickle','binding':parent/'units-and-binding.pickle'}
 for block in (17,16):
  tag='block%d-grade01'%block
  scientific[tag+'-forcing']=parent/(tag+'-forcing-operands.pickle')
  scientific[tag+'-domain']=parent/(tag+'-domain.pickle')
 for name,p in scientific.items():
  if p.parent==parent:expected=parent_index[p.name]['sha256']
  elif p.parent==reference:expected=kernel_cp['artifactSha256'][p.name]
  else:expected=outgoing_cp['primaryArtifacts'][p.name]['sha256']
  assert sha(p)==expected,(name,'accepted saved artifact identity')
 for name,record in [('saved-context.json',inspection_cp['context'])]+[
   (k+'-inspection.json',v) for k,v in inspection_cp['caseReports'].items() if k.endswith('grade01')]:
  assert sha(inspection/name)==record['sha256']
 worker=M/(PREFIX+'response.py');runner=M/(PREFIX+'review_run.py');launcher=M/(PREFIX+'review_launch.py')
 manifest=M/(PREFIX+'inputs.json');pending=M/(PREFIX+'pending_gate.json')
 inputs={name:route(p) for name,p in scientific.items()}
 for name,p in [('kernelCheckpoint',M/'S11c_d_reference_kernel_checkpoint.json'),
   ('outgoingCheckpoint',M/'S11c_d_outgoing_prescription_checkpoint.json'),
   ('twoAsymptoteCheckpoint',M/'S11c_d_two_asymptote_checkpoint.json'),
   ('inspectionCheckpoint',M/'S11c_d_two_asymptote_saved_inspection_checkpoint.json'),
   ('physicalInput',M/'S11c_d_variable_profile_development_input.json'),
   ('guard',ROOT/'scripts/s11c_guarded_run.py'),('supervisor',M/'S11c_d_end_normalization_run.py'),
   ('methodDescription',D/'S11c_d_localized_thickness_implementation.md')]:
  inputs[name]=route(p)
 save(manifest,{'status':'PINNED_INPUTS_FOR_UNSUBMITTED_BUILD_REVIEW',
   'scope':'Two existing thickness-grade sources; no new physical input or producer replay.',
   'case':'LAB_HELD_RHO4_CONSTANT','blocks':[17,16],'grade':[0,1],
   'inputs':inputs,'scientificPayloads':list(scientific),'automaticRetry':False})
 save(pending,{'status':'PENDING_INDEPENDENT_METHOD_REVIEW_AND_STAGE_APPROVAL',
   'workerSha256':sha(worker),'inputManifestSha256':sha(manifest),
   'scopeExplicitlyApproved':False,'substantiveReviewFindingsClosed':False,
   'seconds':900,'nativeSeconds':840,'automaticRetry':False,
   'sharedGuardSha256':sha(ROOT/'scripts/s11c_guarded_run.py'),
   'supervisorSha256':sha(M/'S11c_d_end_normalization_run.py')})
 context=json.loads((inspection/'saved-context.json').read_text())
 evidence={'provenance':{'inspectionCheckpoint':route(M/'S11c_d_two_asymptote_saved_inspection_checkpoint.json'),
   'savedContext':route(inspection/'saved-context.json')},
   'physicalInput':context['physicalInput'],'forcingUnits':context['forcingUnits'],
   'fieldEntryUnits':context['fieldEntryUnits'],
   'nativeThicknessCells':[v for v in context['preparedNativeCells'] if v['grade']==[0,1]],'cases':[]}
 for block in (17,16):
  p=inspection/('block%d-grade01-inspection.json'%block);data=json.loads(p.read_text())
  selected={k:data[k] for k in ('case','checks','directForcing','localForcingPart','lift','seed',
    'carrierCount','unreducedCoefficientZeroFlags','nativeTermsBySide','regulatorPresentCount',
    'planeJetStructuralEquality','planeJetEstablished','forcingPairingEstablished')}
  selected['sourceReport']=route(p)
  selected['nonlocalDomainOperands']=[{**{k:v[k] for k in ('side','address','operand','hasProfileAbelRegulator','status')},
    'orderedLimits':[x['orderedLimits'] for x in v['orderedIntegrals']]} for v in data['nonlocalDomainOperands']]
  evidence['cases'].append(selected)
 save(packet/'saved-source-evidence.json',evidence)
 (packet/'source-schema.md').write_text('''# Saved schemas used by this proposal

All scientific payloads are pinned in the input manifest. No pickle was
restored to prepare this packet. saved-source-evidence.json is selected JSON
from the completed guarded reader, with original report routes/hashes.

The forcing bundle has seed, lift, forcing, liftedAction, direct, decomposed,
forcingUnits and the unresolved plane-jet premise. forcing has local/nonlocal/
total matrices and addressed terms with nativeIntegral, coefficient and value.
The domain bundle retains directForcing and profileForcingTerms, joined by the
saved reader. Native term values already contain the selected plane field.

The outgoing candidate contains fixedPhysicalInput, fixedSymbol, fixedInverse,
blocks (index/k/residue/direction/deltaCoefficient), symmetric exclusion
intervals, branch/unit information and the explicitly limited scope.
The reference context supplies normalMomentum/fourierMass/profileRegulator;
physical-units supplies inverseEntryUnits and spectralMeasureUnit.
units-and-binding supplies physicalInput/forcingUnits/endFieldUnits.

The two source grades have zero liftedAction and no profile Abel regulator.
Their force is the negative saved native total. The worker preserves each
actual term and does not import its producer. Included previous source files
show serialization schemas; they are read-only evidence, never a run queue.
''')
 source_paths=[worker,manifest,pending,runner,launcher,Path(__file__),
  D/'S11c_d_localized_thickness_implementation.md',D/'S11c_d_localized_thickness_response_plan.md',
  D/'S11c_d_SCATTERING_FORM_AMENDMENT.md',D/'S11c_d_FORM_build_directive.md',
  M/'S11c_d_two_asymptote_saved_inspection_report.md',M/'S11c_d_two_asymptote_checkpoint.json',
  M/'S11c_d_outgoing_prescription_checkpoint.json',M/'S11c_d_outgoing_prescription_candidate_report.md',
  M/'S11c_d_reference_kernel_report.md',M/'S11c_d_reference_kernel.py',
  M/'S11c_d_two_asymptote_forcing_continue.py',M/'S11c_d_outgoing_prescription_unlimited.py',
  M/'S11c_d_variable_profile_development_input.json',ROOT/'scripts/s11c_guarded_run.py',
  M/'S11c_d_end_normalization_run.py']
 source_records=[]
 for p in source_paths:
  relative=str(p.relative_to(ROOT));dest=packet/relative;dest.parent.mkdir(parents=True,exist_ok=True)
  shutil.copyfile(p,dest);source_records.append({'path':relative,'sourceSha256':sha(p)})
 prompt=M/(PREFIX+'review_prompt.md');shutil.copyfile(prompt,packet/'review-prompt.md')
 source_records.append({'path':str(prompt.relative_to(ROOT)),'sourceSha256':sha(prompt)})
 save(packet/'packet-index.json',{'files':sorted(str(p.relative_to(packet)) for p in packet.rglob('*') if p.is_file()),
   'readFirst':['review-prompt.md','research/pde_ledger_v3/directives/S11c_d_localized_thickness_implementation.md',
                'research/pde_ledger_v3/_measurements/S11c_d_localized_thickness_response.py','saved-source-evidence.json'],
   'note':'Read only indexed packet files; absolute scientific artifact routes are metadata, not external read/restore authorization.'})
 hashes={str(p.relative_to(packet)):sha(p) for p in sorted(packet.rglob('*')) if p.is_file()}
 packet_hash=hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest()
 with tarfile.open(RUN/'packet.tar.gz','w:gz') as tar:tar.add(packet,arcname='packet')
 state={'status':'PREPARED_NOT_LAUNCHED','currentThread':THREAD,'reviewers':['claude','grok'],
        'packetDirectory':str(packet),'fileHashes':hashes,'packetSha256':packet_hash,
        'archiveSha256':sha(RUN/'packet.tar.gz'),'sourceRecords':source_records,
        'fileCount':len(hashes),'bytes':sum(p.stat().st_size for p in packet.rglob('*') if p.is_file()),
        'newExternalSubmissionAuthorized':False,'productionAuthorized':False,'scienceCalls':0}
 save(RUN/'state.json',state)
 (RUN/'completion_message.txt').write_text('''User-authorized local script watcher: completion.
The explicitly approved localized-thickness method/build reviews have finished or reported an error.
Root: /var/projects/toy_physics/_scratch/s11c/s11c-d-localized-thickness-20260928/build-review.
Read actual state/launch/static/coordinator outcomes, both run receipts, literal Claude result/Grok text,
stderr, fixed packet/archive/prompt/index and explicit-packet-approval. Both reports must finish before
adjudication or editing. No peer sharing, automatic rerun or optional-wording review cycle. Preserve
literal verdicts and exact reviewed bytes; no fresh independent CLEAR from local corrections.
No science launch is authorized by this event. The worker proposes transforms and the full saved
PV/delta action on only two fixed thickness-grade forcing cases. Assess source-specific distributional
reading, Fourier phases/mass, profile removable values, branch/domain checks and actual controls.
No new physical input, mode/root/LU/producer/end-field replay/current normalization or profile solve.
Any later approved stage needs a fresh pinned gate and the normal shared guard around the normalization
supervisor: 900s outer/840s native,2GiB,zero swap,one CPU,nice15,32 tasks,one native thread. Hook first
for session 01a0e01b-ef84-7192-817f-584cda5d339b; no unlimited inheritance, overlap, fallback or retry.
Retain practical toy-model scope, all earlier science/failures/review literals and incident history.
Full response/Green/FORM/A11/A12, density/mixed response, current normalization, general diagonal
extension, retarded equivalence and radiating coverage remain open. Leave Lean/S11_lean/shared guard
and protected builder suffix f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2 untouched.
Read JSON/source/hash metadata only outside containment. Update durable execution/recovery status.
No model polling or recurring task. Scratch copies/archives are runtime artifacts, never committed.
''')
 # Syntax and existing instrument identities, without importing scientific code.
 tree=ast.parse(worker.read_text());compile(tree,str(worker),'exec')
 old=ast.parse((M/'S11c_d_two_asymptote_saved_inspect.py').read_text())
 names=('require','digest','route','save','containment','SavedCodec','Journal')
 for name in names:
  left=next(n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name==name)
  right=next(n for n in old.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name==name)
  assert ast.dump(left,include_attributes=False)==ast.dump(right,include_attributes=False),name
 for p in (runner,launcher,Path(__file__)):compile(p.read_text(),str(p),'exec')
 assert not any(isinstance(n,ast.ImportFrom) and n.module and 'S11c' in n.module for n in ast.walk(tree))
 forbidden={'lu_solve','LUsolve','inv','inverse','solve','nsolve','roots','nroots','series','diff','integrate','doit','limit'}
 called={n.func.attr if isinstance(n.func,ast.Attribute) else n.func.id for n in ast.walk(tree)
         if isinstance(n,ast.Call) and isinstance(n.func,(ast.Attribute,ast.Name))}
 assert not forbidden&called,forbidden&called
 # A pending gate must refuse before containment/imports/output creation.
 result=subprocess.run(['python3',str(worker),'--input-manifest',str(manifest),'--gate-receipt',str(pending),
                        '--run-directory',str(RUN/'must-not-exist')],capture_output=True,text=True)
 assert result.returncode!=0 and 'method gate incomplete' in result.stderr and not (RUN/'must-not-exist').exists()
 (RUN/'pending-gate-refusal.stdout').write_text(result.stdout)
 (RUN/'pending-gate-refusal.stderr').write_text(result.stderr)
 # Verify archive members byte-for-byte; no external tool or scientific payload is run.
 with tarfile.open(RUN/'packet.tar.gz') as tar:
  archived={str(Path(v.name).relative_to('packet')):hashlib.sha256(tar.extractfile(v).read()).hexdigest()
            for v in tar.getmembers() if v.isfile()}
 assert archived==hashes
 # The verifier's missing explicit approval is an expected refusal, not a reviewer run.
 namespace={'__name__':'preparation_verify_only','__file__':str(runner)}
 exec(compile(runner.read_text(),str(runner),'exec'),namespace)
 try:namespace['verify'](state)
 except FileNotFoundError as error:assert error.filename==str(RUN/'explicit-packet-approval.json')
 else:raise AssertionError('missing external approval was accepted')
 static={'status':'PASS_SOURCE_AND_INSTRUMENT_CHECKS_NO_SCIENCE',
   'workerSha256':sha(worker),'inputManifestSha256':sha(manifest),'runnerSha256':sha(runner),
   'launcherSha256':sha(launcher),'completionMessageSha256':sha(RUN/'completion_message.txt'),
   'unchangedInstrumentASTs':list(names),'scientificImportsOrRestorations':0,
   'forbiddenProducerSolveReplayCalls':0,'pendingGateRefused':True,'missingReviewApprovalRefused':True,
   'archiveMatchesPacket':True,'allSelectedAcceptedArtifactHashesMatch':True}
 save(RUN/'static-review.json',static)
 save(M/(PREFIX+'preparation.json'),{'status':'PREPARED_NOT_SUBMITTED_NOT_RUN',
   'recordedUtc':datetime.now(timezone.utc).isoformat(),'worker':route(worker),'inputManifest':route(manifest),
   'pendingGate':route(pending),'reviewRoot':str(RUN),'packetSha256':packet_hash,
   'archiveSha256':state['archiveSha256'],'packetFiles':state['fileCount'],'packetBytes':state['bytes'],
   'sourceEvidence':route(packet/'saved-source-evidence.json'),'staticChecks':static,
   'readOnlyReviewers':['claude','grok'],'newExternalSubmissionAuthorized':False,
   'scientificLaunchAuthorized':False,'automaticRetry':False,'currentThread':THREAD,
   'scratchPolicy':'Runtime packet copies/archive/logs stay ignored. Canonical implementation and preparation files are under research/pde_ledger_v3.'})
 print(json.dumps({'status':'PREPARED_NOT_SUBMITTED_NOT_RUN','packetFiles':state['fileCount'],
                   'packetBytes':state['bytes'],'packetSha256':packet_hash,'sourceFiles':len(inputs)}))


if __name__=='__main__':main()
