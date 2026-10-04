from pathlib import Path
from datetime import datetime,timezone
from decimal import Decimal
from collections import Counter
import json,hashlib,sqlite3
R=Path('/var/projects/toy_physics');M=R/'research/pde_ledger_v3/_measurements';B=R/'_scratch/s11c/s11c-defect-packet-action-20261003/local-01';C=B/'complete';P='S11c_d_defect_packet_local'
read=lambda p:json.loads(Path(p).read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
counts=Counter()
def ck(v,label):
 if v is not True:raise AssertionError(label)
 counts[label.split(':',1)[0]]+=1
ZERO={'text':'0','srepr':'Integer(0)'}
# Decimal is used only to inspect/report persisted scalar decimal fields. No MP,
# symbolic or numerical-library objects are restored and no quadrature is rerun.
def dec(v):return Decimal(v['decimal'])
def same_numeric_encoding(a,b):
 if isinstance(a,dict) and 'mpf' in a:return isinstance(b,dict) and a['mpf']==b.get('mpf')
 if isinstance(a,dict) and 'mpc' in a:return isinstance(b,dict) and len(b.get('mpc',[]))==2 and all(same_numeric_encoding(x,y) for x,y in zip(a['mpc'],b['mpc']))
 return a==b

g=read(M/(P+'_gate.json'));m=read(M/(P+'_inputs.json'));launch=read(B/'launch.json');stage=read(B/'stage-start.json');active=read(B/'active.json');guard=read(B/'resource-guard/outcome.json');duration=read(B/'resource-guard/duration-validation.json');limit=read(B/'resource-guard/limit-validation.json');inv=read(B/'resource-guard/invocation.json');native=read(C/'containment.json');watch=read(B/'completion-watcher/state.json')
ck(launch['command']==stage['command']==g['command'] and launch['gateSha256']==stage['gateSha256']==sha(M/(P+'_gate.json')),'launch commands')
ck(active['command']==g['command'][g['command'].index('-u')-1:] and inv['command']==g['command'][g['command'].index('--')+1:],'supervisor/guard commands')
ck(active['status']=='completed' and active['exitCode']==guard['exitCode']==guard['childOutcome']['exitCode']==read(B/'coordinator-outcome.json')['exitCode']==0 and guard['childOutcome']['guardReason'] is None,'completed no failure')
ck(duration['verified'] and duration['actual']=={'RuntimeMaxUSec':'infinity','Restart':'no'} and limit['verified'] and guard['limitsVerified'],'no-deadline enforced')
a=limit['actual'];ck(a==read(B/'resource-guard/effective-limits.json') and a['memory.max']=='4294967296' and a['memory.swap.max']=='0' and a['pids.max']=='32' and a['nativeAddressSpaceBytes']==[4294967296]*2 and len(a['affinity'])==1 and set(a['threads'].values())=={'1'},'resource limits')
ck(native['nativeAddressSpace']==4294967296 and native['durationLimits'] is None and native['affinity']==a['affinity'] and native['cgroup']==a['group'],'native limits')
ck(inv['pool']=='s11c-near-unity' and inv['poolMemoryMax']==17179869184 and inv['priorityPolicy']=='desktop-managed','pool')
ck(watch['status']=='finished' and watch['thread']=='01a0e01b-ef84-7192-817f-584cda5d339b' and watch['startedUtc']<stage['utc'] and watch['notifications']['completion']['status']=='queued','hook first/completed')
for n in ['defect_packet_local.stderr','coordinator.stderr','guard-launch.stderr','resource-guard/stderr','watcher.stderr']:ck((B/n).stat().st_size==0,'empty stderr:'+n)
ck((B/'defect_packet_local.stdout').read_bytes()==(C/'checks.json').read_bytes(),'stdout/checks exact bytes');ck(not (C/'failure.json').exists(),'no failure file')
samples=[json.loads(x) for x in (B/'resource-guard/resource-samples.jsonl').read_text().splitlines()];events=Counter()
for s in samples:
 ev={k:int(v) for k,v in (x.split() for x in s['memory.events'].splitlines())}
 for k,v in ev.items():events[k]=max(events[k],v)
 ck(int(s['memory.swap.current'])==0 and all(ev[k]==0 for k in ['oom','oom_kill','oom_group_kill']) and int(s['pids.current'])<=32 and s['hostAvailableBytes']>=4*1024**3,'resource sample')
post=read(C/'posthashes.json');snap=read(B/'source-snapshots.json');copies=read(C/'saved-copy-index.json')
ck(set(snap)==set(g['sourcePins'])|{str(M/(P+'_gate.json')),str(M/(P+'_inputs.json'))},'source snapshot census')
for p,r in post.items():ck(r['error'] is None and r['expected']==r['actual']==sha(p),'posthash:'+p)
for p,r in snap.items():ck(sha(p)==r['sha256']==sha(r['path']),'snapshot:'+p)
ck(set(copies)==set(m['savedInputs']) and len(copies)==69,'copy census')
for n,r in copies.items():ck(r['source']==m['savedInputs'][n]['path'] and r['sha256']==m['savedInputs'][n]['sha256']==sha(r['source'])==sha(C/r['path']) and (C/r['path']).stat().st_size==r['bytes'],'copy:'+n)
extra=read(g['additionalIdentityCopyIndex']);ck(sha(g['additionalIdentityCopyIndex'])==g['additionalIdentityCopyIndexSha256'],'additional index')
for p,r in extra.items():ck(sha(p)==r['sha256']==sha(r['path'])==g['additionalIdentityPins'][p],'additional identity:'+p)
arts=read(C/'artifact-index.json');ops=read(C/'operation-index.json')
for n,r in arts.items():ck(n==r['path'] and sha(C/n)==r['sha256'] and (C/n).stat().st_size==r['bytes'],'artifact:'+n)
ck(len(ops)==1 and ops[0]['name']=='packet-local' and all(ops[0][k]==arts[ops[0][k]['path']] for k in ['input','result']),'complete top operation')
ck(read(C/'packet-local-input.json')=={'manifestSha256':g['manifestSha256'],'buildReviewSha256':g['buildReviewRecordSha256']},'stage actual inputs')
result=read(C/'checks.json');ret=read(C/'packet-local-return.json');ck({k:v for k,v in result.items() if k!='wallSeconds'}==ret and ret['status']=='BOUNDED_LOCAL_PACKET_ACTIONS_COMPLETE_PRESSURE_PENDING' and ret['cells']==16 and ret['packets']==2 and ret['completePacketAction'] is False and ret['pressureActionEvaluated'] is False,'complete bounded return')
zeros=[]
for name in arts:
 if name.endswith('-return.json') and (C/(name[:-12]+'-raw.json')).exists():
  stem=name[:-12];r=read(C/name);ck(r['cancelled']==ZERO and r['raw']==read(C/(stem+'-raw.json'))['residual'] and set(read(C/(stem+'-input.json')))=={'left','right'},'complete new zero:'+stem);zeros.append(stem)
saved=lambda n:read(C/'saved'/n)
def eq(a,b):
 if isinstance(a,dict):
  if 'mpf' in a:return isinstance(b,dict) and a['mpf']==b.get('mpf')
  if 'mpc' in a:return isinstance(b,dict) and len(b.get('mpc',[]))==2 and all(eq(x,y) for x,y in zip(a['mpc'],b['mpc']))
  return isinstance(b,dict) and set(a)==set(b) and all(eq(v,b[k]) for k,v in a.items())
 if isinstance(a,list):return isinstance(b,list) and len(a)==len(b) and all(eq(x,y) for x,y in zip(a,b))
 return a==b

def finite(v):
 if isinstance(v,dict):
  if 'mpf' in v:return set(v)=={'mpf','decimal'} and len(v['mpf'])==4 and dec(v).is_finite()
  if 'mpc' in v:return set(v)=={'mpc'} and len(v['mpc'])==2 and all(finite(x) for x in v['mpc'])
  return all(finite(x) for x in v.values())
 if isinstance(v,list):return all(finite(x) for x in v)
 return True

def root(v,nonzero=False):
 if 'mpf' in v:return dec(v)>0 if nonzero else dec(v)>=0
 return dec(v['mpc'][0])==0 and (dec(v['mpc'][1])>0 if nonzero else dec(v['mpc'][1])>=0)

def receipt(name,h,p):return {'record':name,'sha256':h,'bytes':len(p)}

def get(name):
 row=db.execute('SELECT sha256,payload FROM records WHERE name=?',(name,)).fetchone();ck(row is not None,'referenced record exists')
 h,p=row;ck(hashlib.sha256(p).hexdigest()==h,'referenced record hash');return json.loads(p),receipt(name,h,p)

from fractions import Fraction
cells=saved('selected/local-cells.json')['selected'];partition=saved('native/THETA-partition.json');children={};inherited=0
ck(cells==[c for c in saved('local/all-local-cells.json') if c['row']=='THETA_BALANCE' and c['field']=='e_W'] and len(cells)==16,'complete native selection')
for alias in m['localBatches']:
 batch=saved(alias);inp=saved(alias.replace('-return','-input'))
 ck(inp['children']==[next(v for v in partition['children'] if v['childIndex']==e['childIndex']) for e in batch],'native batch inputs')
 for e in batch:
  orig=next(v for v in partition['children'] if v['childIndex']==e['childIndex'])
  ck(e['completed'] is True and e['sourceConstructor']==orig['constructorText'] and e['sourceSha256']==orig['sha256']==hashlib.sha256(e['sourceConstructor'].encode()).hexdigest(),'native child actual constructor')
  for proof in e['identities']:ck(proof['cancelled']==ZERO,'inherited child zero');inherited+=1
  ck(e['childIndex'] not in children,'unique child');children[e['childIndex']]=e
ck(set(children)==set(partition['localChildIndices']),'complete native local partition')
source=Path(partition['source']['source'])
with source.open('rb') as f:line=next(line for i,line in enumerate(f,1) if i==partition['source']['valueLine'])
ck(hashlib.sha256(line).hexdigest()==partition['source']['sourceLineSha256'],'original native source line')
native_row=saved('native/THETA-row.json');unit_joins=saved('units/source-text-joins.json')
ck(partition['source']==native_row['source']==unit_joins['provenance'] and hashlib.sha256(native_row['fullConstructor'].encode()).hexdigest()==unit_joins['rowHashes']['THETA_BALANCE'],'actual native unit source join')
specs=[];used=[];matches_record={}
for i,cell in enumerate(cells):
 matches=[v for v in children.values() if v['fieldColumn']==4 and v['xOrder']==cell['xOrder']]
 summands=[next(r['value'] for r in v['mappedGrades'] if r['grade']==cell['grade']) for v in matches]
 ancestry=read(C/f'cell-{i}-ancestry.json')
 ck(ancestry=={'cell':cell,'childIds':[v['childIndex'] for v in matches],'summands':summands} and cell['sourceChildren']==ancestry['childIds'] and cell['summands']==summands,'complete native ancestry')
 for proof in cell['identities']:ck(proof['cancelled']==ZERO,'inherited cell zero');inherited+=1
 poly=cell['polynomial'];orig=read(C/f'cell-{i}-certificate-originals.json');match=read(C/f'cell-{i}-certificate-match.json')
 ck(orig['coefficients']==poly['coefficients'] and orig['certificates']==poly['coefficientCertificates'] and match['originalCertificatesUnchanged'] is True,'unaltered certificate originals')
 indices=match['certificateIndicesInCoefficientOrder'];matches_record[str(i)]=indices
 ck(sorted(indices)==list(range(len(orig['certificates']))) and len(indices)==len(orig['coefficients']),'certificate exact multiplicity')
 for j,k in enumerate(indices):
  cert=orig['certificates'][k];ck(cert['finite'] is True and cert['value']==orig['coefficients'][j]['value'],'actual exact certificate match')
  op=read(C/f'new-cell-{i}-certificate-components-{j}-input.json');ck(op['left']==cert['value'],'actual certificate proof operand')
 adapter=read(C/f'cell-{i}-adapter.json')
 ck(adapter['cell']==cell and adapter['nativeSourceChildren']==cell['sourceChildren'] and adapter['epsilonAlreadyExtracted'] is True and adapter['timeAndTangentsAlreadyAbsorbed'] is True,'actual adapter source argument')
 ck(adapter['pairingUnit']=='U_dual_THETA * M_ref * L_ref^-2 * T_ref^-1' and adapter['dxUnit']==[1,0,0],'local integration unit')
 for pair in adapter['coefficientsAscending']:
  ck(len(pair)==2 and all(str(Fraction(v))==v for v in pair),'literal rational transport')
 ck(cell['epsilonPower']==(0 if poly['zero'] else 1),'epsilon once')
 specs.append({'cellIndex':i,'xOrder':cell['xOrder'],'grade':cell['grade'],'coefficientsAscending':adapter['coefficientsAscending'],'zero':poly['zero']})
 for e in matches:
  idx=e['childIndex']
  if idx in used:continue
  used.append(idx);u=read(C/f'new-unit-child-{idx}.json')
  ck(u['originalChild']==e and u['noCoefficientRecalculation'] is True and u['mappedCoefficientUnit']==u['expected'] and u['rowUnit']==adapter['rowUnit'],'actual unbound unit operands/result')
unit=read(C/'local-unit-conclusion.json');ck(unit['checkedNativeChildren']==used and unit['pressureSummandUnitsEstablished'] is False and unit['powerInterpretation'] is False,'complete narrow unit conclusion')
ck(len(used)==len(set(used)) and len(specs)==16 and sum(s['zero'] for s in specs)==8,'native unit/cell census')
physical=saved('preflight/physical-plan.json');context=saved('local/context.json')
ck(context['physical']==saved('physical-input.json') and context['frequencyOverride']=={'old':'1','actual':3} and context['physical']['parameters']['L_W']=='10' and physical['frequency']==3 and read(C/'new-local-tail-arguments.json')['actualPacket']==physical,'actual physical/profile inputs')
DB=C/'local-evidence.sqlite';dbreceipt=read(C/'numerical-journal-receipt.json');ck(DB.stat().st_size==dbreceipt['bytes'] and sha(DB)==dbreceipt['sha256'],'SQLite full byte receipt')
db=sqlite3.connect('file:'+str(DB)+'?mode=ro',uri=True);db.execute('PRAGMA query_only=ON');ck(db.execute('PRAGMA integrity_check').fetchone()==('ok',),'SQLite integrity')
states={};returns={};requests={};completed={};grades={};controls=None;panels_count=Counter();total_rows=0;payload_bytes=0;maxima={'difference':Decimal(0),'toleranceFraction':Decimal(0),'BempiricalError':Decimal(0),'envelope':Decimal(0)};max_at={};comparison_count=0;analytic_count=0
for name,h,payload in db.execute('SELECT name,sha256,payload FROM records ORDER BY rowid'):
 total_rows+=1;payload_bytes+=len(payload);ck(hashlib.sha256(payload).hexdigest()==h,'SQLite row bytes');d=json.loads(payload);r=receipt(name,h,payload)
 ck('failed' not in name and finite(d),'finite complete saved record')
 if name=='restored-rules':
  ck(d['constructorsCalled'] is False and d['momentReturnsInherited'] is True and d['completeOperands']=={n:saved('rules/'+n+'.json') for n in ['A-GL24','A-GL48','B-G7-K15']},'complete rules restored no constructor');continue
 if name=='controls/complete':controls={'data':d,'receipt':r};continue
 if name=='complete-local-actions':
  ck(d=={'packets':{p:v['receipt'] for p,v in grades.items()},'controls':controls['receipt'],'completePacketAction':False,'pressurePending':True},'complete local receipt joins');finalreceipt=r;continue
 path=name.split('/');ck(path[0]=='packet' and path[1] in ['0','kappa'],'fixed packet family');carrier=path[1]
 if path[2]=='local-grade-totals':
  ck(d['pressureIncluded'] is False and d['finiteGradeWeightsApplied'] is False and len(d['grades'])==4,'local grade totals scope')
  for v,g0 in zip(d['grades'],[[0,0],[1,0],[0,1],[1,1]]):
   ck(v['grade']==g0 and v['completeCellReceipts']==[completed[f'packet/{carrier}/cell/{s["cellIndex"]}']['receipt'] for s in specs if s['grade']==g0],'actual complete grade cell receipts')
  grades[carrier]={'data':d,'receipt':r};continue
 ck(path[2]=='cell','native cell request');i=int(path[3]);ck(0<=i<16,'fixed cell index');prefix='/'.join(path[:4]);at=4;mutant=None
 if path[at].startswith('mutant-'):
  mutant=path[at][7:];ck((mutant,i,carrier) in [('Leibniz',6,'kappa'),('conjugation',0,'kappa')],'predeclared mutant address');prefix+='/'+path[at];at+=1
 kind=path[at]
 if kind=='input':
  ck(d['spec']==specs[i] and d['carrier']==carrier and d['mutant']==mutant and d['centers']==['-5/2','5/2'] and d['width']==8 and d['length']==10 and d['nativeTimeTangentsAlreadyAbsorbed'] is True,'native actual local request')
  trials=d['radiusTrials'];ck([dec(v['R']) for v in trials]==list(range(8,int(dec(d['R']))+1,8)) and all(dec(v['tail'])>Decimal('1e-14') for v in trials[:-1]) and Decimal(0)<=dec(trials[-1]['tail'])<=Decimal('1e-14'),'actual full tail prefix')
  requests[prefix]=d;continue
 if kind=='complete':
  ck(d['spec']==specs[i] and d['carrier']==carrier,'complete native cell arguments')
  if specs[i]['zero']:
   ck(d['explicitZero'] is True and d['quadratureCalls']==0 and dec(d['value'])==0 and prefix not in requests and all(a==b=='0' for a,b in specs[i]['coefficientsAscending']),'actual explicit zero no quadrature')
  else:
   ck(d['mutant']==mutant and d['quadratureErrorProof'] is False and set(d['routes'])=={'A24','A48','A48Rplus8','B'},'complete four routes')
   for route,value in d['routes'].items():ck(eq(value,returns[prefix+'/'+route]['value'][0]),'actual numerical completed return join')
   ck(eq(d['value'],d['routes']['A48']) and eq(d['BActualEmpiricalError'],returns[prefix+'/B']['summedEmpiricalErrors'][0]),'actual reference and empirical error')
   expected=['A24','A48Rplus8','B']+(['analyticConstant'] if len(specs[i]['coefficientsAscending'])==1 and mutant is None else [])
   ck([v['route'] for v in d['comparisons']]==expected,'complete comparison census')
   for v in d['comparisons']:
    comparison_count+=1;analytic_count+=v['route']=='analyticConstant'
    ck(v['passed'] is True and eq(v['value'],d['analyticConstant'] if v['route']=='analyticConstant' else d['routes'][v['route']]) and Decimal(0)<=dec(v['difference'])<=dec(v['tolerance']),'actual saved absolute comparison')
    for key,value in [('difference',dec(v['difference'])),('toleranceFraction',dec(v['difference'])/dec(v['tolerance']))]:
     if value>maxima[key]:maxima[key]=value;max_at[key]={'record':name,'comparison':v}
   ck(all(Decimal(0)<=dec(v)<=Decimal('1e-14') for v in d['tails'].values()) and dec(d['empiricalEnvelope'])>=max(dec(v['difference']) for v in d['comparisons']),'saved tails and empirical envelope')
   maxima['envelope']=max(maxima['envelope'],dec(d['empiricalEnvelope']))
  completed[prefix]={'data':d,'receipt':r};continue
 route=kind;key=prefix+'/'+route;sub=path[at+1];ck(route in ['A24','A48','A48Rplus8','B'],'known route')
 if sub=='input':
  inp=requests[prefix];rad=int(dec(inp['R']))+(8 if route=='A48Rplus8' else 0)
  ck(d['spec']==specs[i] and d['carrier']==carrier and d['mutant']==mutant and d['R']==rad and d['precision']==(50 if route=='B' else 30) and d['order']==(None if route=='B' else int(route[1:3])),'actual independent route arguments')
  ps=d['panels'];ck(len(ps)==(3 if route=='B' else 4)*rad and dec(ps[0][0])==-rad and dec(ps[-1][1])==rad and all(dec(x)<dec(y) for x,y in ps) and all(eq(ps[j][1],ps[j+1][0]) for j in range(len(ps)-1)),'exact route finite coverage')
  states[key]={'panels':ps,'seen':set(),'leaves':{},'count':0};continue
 st=states[key]
 if sub=='panel':
  st['count']+=1;panels_count[route]+=1;n=15 if route=='B' else int(route[1:3])
  ck(len(d['points'])==len(d['values'])==len(d['kernelOperands'])==n and all(len(v)==1 for v in d['values']),'full panel nodes/values/operands')
  lo,hi=map(dec,d['bounds']);ck(all(lo<dec(v)<hi for v in d['points']),'open physical nodes')
  for x,e in zip(d['points'],d['kernelOperands']):
   ck(set(e)=={'x','T','coefficient','waveDerivative','packetProduct','coefficientDerivativeCorrection','mutant'} and eq(x,e['x']) and e['mutant']==mutant and Decimal(-1)<dec(e['T'])<Decimal(1),'actual native integrand operand record')
   if mutant!='Leibniz':ck(dec(e['coefficientDerivativeCorrection'])==0,'native derivative before coefficient')
  if route!='B':
   panel,side=int(path[at+2]),int(path[at+3]);ck(d['bounds']==st['panels'][panel] and d['side']==side and (panel,side) not in st['seen'],'actual squared route panel')
   st['seen'].add((panel,side));ck(len(d['jacobians'])==n and all(dec(v)>0 for v in d['jacobians']) and len(d['return'])==1,'positive Jacobians and return')
  else:
   label=path[at+2];ck(label not in st['leaves'] and len(d['K'])==len(d['G'])==len(d['error'])==1 and dec(d['error'][0])>=0,'adaptive full operands')
   if label[1:].isdigit():ck(d['bounds']==st['panels'][int(label[1:])],'adaptive initial bounds')
   st['leaves'][label]={'bounds':d['bounds'],'K':d['K'],'error':d['error']}
  continue
 if sub=='refinement':
  parent=d['replacedLeaf'];left,right=d['children'];old=st['leaves'][parent];ll=st['leaves'][left];rr=st['leaves'][right]
  ck(left==parent+'L' and right==parent+'R' and eq(ll['bounds'][0],old['bounds'][0]) and eq(rr['bounds'][1],old['bounds'][1]) and eq(ll['bounds'][1],rr['bounds'][0]) and d['parentError']==old['error'] and d['childErrors']==[ll['error'],rr['error']],'actual adaptive replacement');del st['leaves'][parent];continue
 if sub in ['global-initial','global-sum']:
  ck(set(d['activeLeaves'])==set(st['leaves']) and len(d['summedEmpiricalErrors'])==1,'actual global active partition');continue
 if sub=='return':
  ck(len(d['value'])==1,'route complete scalar')
  if route=='B':
   ck(set(d['activeLeaves'])==set(st['leaves']) and d['panelsEvaluated']==st['count'] and d['physicalVariable'] is True and all(Decimal(0)<=dec(v)<=dec(d['budgetPerComponent']) for v in d['summedEmpiricalErrors']),'global empirical criterion')
   maxima['BempiricalError']=max(maxima['BempiricalError'],dec(d['summedEmpiricalErrors'][0]))
  else:ck(st['seen']=={(i,j) for i in range(len(st['panels'])) for j in [0,1]} and d['panels']==len(st['panels']) and d['order']==int(route[1:3]) and d['squareSubstitution'] is True,'complete A panel census')
  returns[key]=d;del states[key];continue
 raise AssertionError('unknown record '+name)
ck(not states and len(completed)==34 and len(requests)==18 and len(returns)==72 and len(grades)==2 and len(controls['data'])==2 and analytic_count==4,'all bounded local work complete')
for v in controls['data']:
 base=completed[v['baseline']['record'][:-9]];changed=completed[v['changed']['record'][:-9]]
 ck(v['baseline']==base['receipt'] and v['changed']==changed['receipt'] and v['cell']==base['data']['spec']==changed['data']['spec'] and v['responsive'] is True,'actual control completed operands')
 ck(max(abs(dec(x)) for x in v['movement']['mpc'])>max(Decimal('1e-12'),dec(v['tenTimesEmpiricalEnvelope'])),'saved control movement exceeds full envelope')
ck(ret['localResult']['completeReceipt']==finalreceipt and ret['localResult']['controlReceipt']==controls['receipt'] and ret['localResult']['packetGradeReceipts']=={p:v['receipt'] for p,v in grades.items()},'top scientific return journal joins')
db.close();ck(sha(DB)==dbreceipt['sha256'],'SQLite postinspection immutable')
record={'status':'ACCEPTED_TWO_LOCAL_GAUSSIAN_ACTIONS_PRESSURE_PENDING','utc':datetime.now(timezone.utc).isoformat(),'root':str(B),'reviewCommit':'abc6c012','readyCommit':'f262339d','launchCommit':'f2981980','allChecksPassed':True,'inspectionCounts':dict(counts),'inspectionChecks':sum(counts.values()),'inspectionKind':'Read-only stdlib source/JSON/hash/opaque SQLite traversal and persisted scalar display comparison. No scientific restoration or calculation replay.','result':result,'newLiteralZeroReturns':len(zeros),'inheritedLiteralZeroReturns':inherited,'nativeUnitChildren':len(used),'coefficientCertificateMatching':matches_record,'jsonArtifacts':len(arts),'savedCopies':len(copies),'sourcePosthashes':len(post),'sourceSnapshots':len(snap),'additionalIdentityCopies':len(extra),'sqlite':{'receipt':dbreceipt,'records':total_rows,'payloadBytes':payload_bytes,'panelCounts':dict(panels_count),'readonlyIntegrity':'ok'},'comparisons':comparison_count,'analyticConstantComparisons':analytic_count,'savedScalarMaxima':{k:str(v) for k,v in maxima.items()},'maximumRecords':max_at,'gradeResults':grades,'controls':controls,'resources':{'samples':len(samples),'memoryPeakBytes':max(int(s['memory.peak']) for s in samples),'hostAvailableMinimumBytes':min(s['hostAvailableBytes'] for s in samples),'eventsMax':dict(events),'swapBytes':0,'affinity':a['affinity'],'duration':duration['actual'],'interpretation':'Actual cgroup counters, no process/cache split inferred.'},'strictStderrBytes':0,'stdoutChecksByteIdentity':True,'independentBuildClearance':False,'literalBuildHistory':'Claude NEEDS REVISION / Grok CLEAR; tested exact certificate matching and evidence persistence repair under standing tooling authority, not independent CLEAR.','limitations':['Local THETA/e_W part only:16cells,2carriers,2sensitivity controls; pressure and complete packet action pending.','Adaptive errors/cross-route differences empirical, not rigorous quadrature error bounds.','Full pressure summand units, each actual outer Fourier/inner request, collision-aware independent full actions and full addressed controls remain.','No scattering/current/loss, inverse, drain, calibration or sweep.'],'next':'Continue bounded source-bound pressure summand unit and outer-action readiness; reuse all completed inputs/returns. No routine permission question.'}
with (M/(P+'_result_record.json')).open('x') as f:json.dump(record,f,indent=2);f.write('\n')
print(json.dumps({k:record[k] for k in ['status','inspectionChecks','newLiteralZeroReturns','inheritedLiteralZeroReturns','nativeUnitChildren','sqlite','savedScalarMaxima','resources']},indent=2))
