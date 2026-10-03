from pathlib import Path
from datetime import datetime,timezone
from decimal import Decimal
from collections import Counter
import json,hashlib,sqlite3
R=Path('/var/projects/toy_physics');M=R/'research/pde_ledger_v3/_measurements';B=R/'_scratch/s11c/s11c-defect-packet-action-20261003/inner-01';C=B/'complete';P='S11c_d_defect_packet_inner'
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
for n in ['defect_packet_inner.stderr','coordinator.stderr','guard-launch.stderr','resource-guard/stderr','watcher.stderr']:ck((B/n).stat().st_size==0,'empty stderr:'+n)
ck((B/'defect_packet_inner.stdout').read_bytes()==(C/'checks.json').read_bytes(),'stdout/checks exact bytes');ck(not (C/'failure.json').exists(),'no failure file')
samples=[json.loads(x) for x in (B/'resource-guard/resource-samples.jsonl').read_text().splitlines()];events=Counter()
for s in samples:
 ev={k:int(v) for k,v in (x.split() for x in s['memory.events'].splitlines())}
 for k,v in ev.items():events[k]=max(events[k],v)
 ck(int(s['memory.swap.current'])==0 and all(ev[k]==0 for k in ['oom','oom_kill','oom_group_kill']) and int(s['pids.current'])<=32 and s['hostAvailableBytes']>=4*1024**3,'resource sample')
post=read(C/'posthashes.json');snap=read(B/'source-snapshots.json');copies=read(C/'saved-copy-index.json')
ck(set(snap)==set(g['sourcePins'])|{str(M/(P+'_gate.json')),str(M/(P+'_inputs.json'))},'source snapshot census')
for p,r in post.items():ck(r['error'] is None and r['expected']==r['actual']==sha(p),'posthash:'+p)
for p,r in snap.items():ck(sha(p)==r['sha256']==sha(r['path']),'snapshot:'+p)
ck(set(copies)==set(m['savedInputs']) and len(copies)==312,'copy census')
for n,r in copies.items():ck(r['source']==m['savedInputs'][n]['path'] and r['sha256']==m['savedInputs'][n]['sha256']==sha(r['source'])==sha(C/r['path']) and (C/r['path']).stat().st_size==r['bytes'],'copy:'+n)
extra=read(g['additionalIdentityCopyIndex']);ck(sha(g['additionalIdentityCopyIndex'])==g['additionalIdentityCopyIndexSha256'],'additional index')
for p,r in extra.items():ck(sha(p)==r['sha256']==sha(r['path'])==g['additionalIdentityPins'][p],'additional identity:'+p)
arts=read(C/'artifact-index.json');ops=read(C/'operation-index.json')
for n,r in arts.items():ck(n==r['path'] and sha(C/n)==r['sha256'] and (C/n).stat().st_size==r['bytes'],'artifact:'+n)
ck(len(ops)==1 and ops[0]['name']=='packet-inner' and all(ops[0][k]==arts[ops[0][k]['path']] for k in ['input','result']),'complete top operation')
ck(read(C/'packet-inner-input.json')=={'manifestSha256':g['manifestSha256'],'buildReviewSha256':g['buildReviewRecordSha256']},'stage actual inputs')
result=read(C/'checks.json');ret=read(C/'packet-inner-return.json');ck({k:v for k,v in result.items() if k!='wallSeconds'}==ret and ret['status']=='BOUNDED_INNER_KERNEL_BANK_COMPLETE_NO_PACKET_ACTION' and ret['points']==38 and ret['distinctHRequests']==37 and ret['controls']==3 and ret['exactHReuses']==1,'complete bounded return')
zeros=[]
for name in arts:
 if name.endswith('-return.json') and (C/(name[:-12]+'-raw.json')).exists():
  stem=name[:-12];r=read(C/name);ck(r['cancelled']==ZERO and r['raw']==read(C/(stem+'-raw.json'))['residual'] and set(read(C/(stem+'-input.json')))=={'left','right'},'complete new zero:'+stem);zeros.append(stem)
saved=lambda n:read(C/'saved'/n)
selected=saved('selected/pressure-addresses.json')['selected'];adapters=saved('preflight/numeric-factor-adapters.json')
ck(selected==[v for v in saved('inventory/THETA_BALANCE-ordered-addresses.json') if v['jet']['channel']=='e_W'] and len(selected)==544,'complete native selection')
for address,route in zip(selected,adapters['addressJoins']):
 label='-'.join([address['face'],address['slot'],address['component']]);entry=adapters['definitions'][label];proof=saved('factors/'+address['fullFactorProof']['proof']+'-operands.json')
 ck(read(C/('address-'+str(address['addressId'])+'.json'))=={'address':address,'route':route,'completeAdapter':entry},'native address full operands')
 ck(route=={'addressId':address['addressId'],'adapter':label} and entry['original']==proof['mappedAddressFactor'] and proof['addressNormalOriginal']==address['normalOriginal'] and proof['requiredMap']==address['fullFactorProof']['completeNormalMap'],'native factor/normal argument joins')
prefixes=['preflight/numeric-factor-'+label for label in adapters['definitions']]+['preflight/new-J-numerical-adapter','preflight/new-D-numerical-adapter','native/plus-bound-native-height','native/minus-bound-native-height','preflight/new-H-contact-adapter','preflight/new-H-subtracted-adapter']
for i,p in enumerate(prefixes):
 r=read(C/('inherited-'+str(i)+'.json'));ck(r=={'input':saved(p+'-input.json'),'return':saved(p+'-return.json'),'functionCalled':False} and r['return']['cancelled']==ZERO,'inherited exact operands/return')
for label,entry in adapters['definitions'].items():
 op=saved('preflight/numeric-factor-'+label+'-input.json');args=saved('preflight/numeric-factor-'+label+'-arguments.json')
 ck(op['left']==entry['mapped'] and op['right']==entry['template'] and args['actualCompleteFactor']==entry['original'] and args['address']['addressId']==entry['firstAddressId'],'native template argument join')
for tag,alias in m['wholeDefinitionInputs'].items():ck(saved('pressure/whole-tags.json')[tag]['savedDefinition']==saved(alias) and saved('pressure/whole-tags.json')[tag]['sha256']==m['savedInputs'][alias]['sha256'],'whole definition actual arguments')
lemma=read(C/'physical-H-and-PV-lemma.json');ck(lemma['nativeBinding']==saved('native/binding.json') and lemma['nativeBinding']['physicalInput']==saved('physical-input.json') and lemma['nativeGeometry']==saved('native/lower-geometry.json') and lemma['nativeFourierContract']==saved('native/fourier-contract.json') and lemma['savedProfileLemma']==saved('preflight/Fourier-profile-lemma.json'),'native physical profile operands')
scale=read(C/'new-native-slope-scale-operands.json');ck(scale['savedScale']==saved('inventory/native-profile-scale-join.json') and scale['physicalLength']==saved('physical-input.json')['parameters']['L_W']=='10','actual native length')
tails=read(C/'new-inner-tail-certificates.json');htail=read(C/'new-H-tail-domination.json');ck(tails['actualBounds']==saved('preflight/saved-global-bound-inputs.json') and htail['actualProfileBound']==tails['actualBounds']['profile'] and htail['T']==tails['T']=={'text':'122','srepr':'Integer(122)'} and htail['integratedCoefficientUpper']=={'text':'605/122','srepr':'Rational(605, 122)'},'actual saved tail operands')
plan=read(C/'fixed-point-plan.json')['points'];ck(len(plan)==38 and len({(v['k'],v['l']) for v in plan})==38,'complete fixed plan')
DB=C/'inner-evidence.sqlite';dbreceipt=read(C/'numerical-journal-receipt.json');ck(DB.stat().st_size==dbreceipt['bytes'] and sha(DB)==dbreceipt['sha256'],'SQLite byte receipt')
db=sqlite3.connect('file:'+str(DB)+'?mode=ro',uri=True);db.execute('PRAGMA query_only=ON');ck(db.execute('PRAGMA integrity_check').fetchone()==('ok',),'SQLite integrity')
# Only literal JSON records and scalar display fields are inspected. No scientific
# objects, kernels, quadrature, MP arithmetic or symbolic residuals are restored.
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

states={};returns={};comps={};comparisons={};points={};Hresults={};Hinputs={};momenta={};controls={};panels_count=Counter();maxima={'difference':Decimal(0),'toleranceFraction':Decimal(0),'BempiricalError':Decimal(0)};max_at={};total_rows=0;payload_bytes=0;reuses=[]
restored_rules=None
for name,h,payload in db.execute('SELECT name,sha256,payload FROM records ORDER BY rowid'):
 total_rows+=1;payload_bytes+=len(payload);ck(hashlib.sha256(payload).hexdigest()==h,'SQLite row bytes');d=json.loads(payload);r=receipt(name,h,payload)
 ck('failed' not in name,'no failed numerical record')
 if name=='restored-rules':
  ck(d['constructorsCalled'] is False and d['momentReturnsInherited'] is True and d['completeOperands']=={n:saved('rules/'+n+'.json') for n in ['A-GL24','A-GL48','B-G7-K15']},'actual restored complete rules');restored_rules=d;continue
 if name=='complete-point-receipts':ck(d==[points[i]['receipt'] for i in range(38)],'complete point receipt list');continue
 if name.startswith('H-reuse/'):
  ck(d['settingsIdentical'] is True and d['completedReturn']==Hresults[d['exactQCoefficient'].replace('/','_')]['receipt'],'exact H completed reuse');reuses.append(d);continue
 if name.startswith('controls/'):
  ck(d['point']==plan[-1] and finite(d),'actual finite control operands');movement=d['movement'];parts=movement.get('mpc',[movement]);ck(any(abs(dec(v))>Decimal('1e-12') for v in parts),'actual responsive control')
  controls[name]=d;continue
 path=name.split('/');prefix='/'.join(path[:2]);kind=path[2]
 ck(path[0] in ['H','point'],'known request family')
 if path[0]=='point':
  index=int(path[1]);ck(0<=index<38,'fixed point index')
  if kind=='request':ck(d['point']==plan[index] and d['routeSettings']==[['A',30,24,122],['A',30,48,122],['A',30,48,124],['B',50,'G7/K15',122]],'actual point request');continue
 if kind=='input' and path[0]=='H':
  ck(d['T']==122 and d['QCoefficient'].replace('/','_')==path[1] and finite(d),'actual H argument/settings');Hinputs[prefix]=d
  ps=d['panels'];ck(all(dec(a)<dec(b) for a,b in ps) and all(eq(ps[i][1],ps[i+1][0]) for i in range(len(ps)-1)),'H input coverage')
  for route in ['A24','A48']:states[prefix+'/'+route]={'panels':ps,'seen':set(),'leaves':{},'count':0,'components':1}
  continue
 if kind=='momentum-assembly':
  ck(eq(d['integral'],returns[prefix+'/A48']['value'][0]) and d['integralRecord']==prefix+'/A48/return' and d['integralSelector']=='value[0]' and d['panelInputRecord']==prefix+'/input' and eq(d['contact'],Hinputs[prefix]['contact']) and d['precision']==30 and d['order']==48 and d['includeContact'] is True,'H actual completed integral assembly');momenta[prefix]={'data':d,'receipt':r};continue
 if kind=='comparison':
  if len(path)==4:
   ck(path[3]=='operands' and d['referenceComponents']==(5 if path[0]=='point' else 1) and all(v==d['referenceComponents'] for v in d['candidateComponents'].values()),'comparison vector census');comps[prefix]=d
  else:
   original=comps[prefix];ck(len(d)==3*original['referenceComponents'],'all route components compared')
   for v in d:
    ck(v['passed'] is True and eq(v['reference'],original['reference'][v['component']]) and eq(v['other'],original['candidates'][v['route']][v['component']]) and Decimal(0)<=dec(v['difference'])<=dec(v['tolerance']),'actual saved comparison')
    for label,value in [('difference',dec(v['difference'])),('toleranceFraction',dec(v['difference'])/dec(v['tolerance']))]:
     if value>maxima[label]:maxima[label]=value;max_at[label]={'record':name,'route':v['route'],'component':v['component'],'difference':v['difference'],'tolerance':v['tolerance']}
   comparisons[prefix]={'data':d,'receipt':r}
  continue
 if kind=='complete':
  ck(d['comparison']==comparisons[prefix]['receipt'],'completed comparison receipt')
  if path[0]=='point':
   ck(d['point']==plan[index] and len(d['pressureMixedAndDirect'])==2 and set(d['normalMixedAndDirect'])=={'plus','minus'} and all(len(v)==2 for v in d['normalMixedAndDirect'].values()) and finite(d),'complete both-face actual point operands')
   for route,values in d['routes'].items():ck(len(values)==5 and eq(values[:4],returns[prefix+'/'+route]['value']),'point complete primitive return join')
   hd,hr=get(d['H']['record']);ck(d['H']==hr and hd['comparison']==comparisons['/'.join(d['H']['record'].split('/')[:2])]['receipt'],'point actual H receipt')
   points[index]={'data':d,'receipt':r}
  else:
   ck(d['momentumAssembly']==momenta[prefix]['receipt'] and eq(d['A48'][0],momenta[prefix]['data']['baseline']) and eq(d['B'],returns[prefix+'/B']['value']) and eq(d['A24'],comps[prefix]['candidates']['A24']) and dec(d['tTailBound'])>0,'H complete input/return joins');Hresults[path[1]]={'data':d,'receipt':r}
  continue
 # Route rows
 route=kind;key=prefix+'/'+route;subkind=path[3]
 ck(route in ['A24','A48','A48T124','B'],'known numerical route')
 if subkind=='cut-operands':
  ck(d['point']==plan[index] and d['T']==(124 if route=='A48T124' else 122),'cut actual point/settings')
  cuts=d['resolvedCuts'];ck(all(dec(cuts[i][0])<dec(cuts[i+1][0]) for i in range(len(cuts)-1)),'strict resolved cut order');continue
 if subkind=='input':
  ps=d['panels'];ck(all(dec(a)<dec(b) for a,b in ps) and all(eq(ps[i][1],ps[i+1][0]) for i in range(len(ps)-1)),'route initial panel coverage')
  if path[0]=='point':
   ck(d['point']==plan[index] and d['T']==(124 if route=='A48T124' else 122) and d['componentOrder']==['J','Dreflected','Dheight','Dquadratic'],'primitive order/point join')
   ck(dec(ps[0][0])==-d['T'] and dec(ps[-1][1])==d['T'],'actual route radius endpoints')
   resolved=d['cuts']['resolvedCuts'];boundary={json.dumps(x['mpf']) for pair in ps for x in pair};ck(all(json.dumps(v['mpf']) in boundary for v,_ in resolved),'all exact labelled cuts retained')
  else:
   ck(route=='B' and d['QCoefficient']==Hinputs[prefix]['QCoefficient'] and d['profileScale']==10 and d['FourierNormalization']=='1/(2pi)','H physical independent inputs')
   trials=d['radiusTrials'];ck(dec(trials[-1]['tail'])<=Decimal('1e-14') and all(dec(x['tail'])>Decimal('1e-14') for x in trials[:-1]),'H physical tail radius selection')
  states[key]={'panels':ps,'seen':set(),'leaves':{},'count':0,'components':4 if path[0]=='point' else 1};continue
 st=states[key];components=st['components']
 if subkind=='panel':
  ck(finite(d),'complete finite panel encoding');st['count']+=1;panels_count[route]+=1
  n=15 if route=='B' else int(route[1:3]);ck(len(d['points'])==len(d['values'])==len(d['kernelOperands'])==n and all(len(v)==components for v in d['values']),'full panel operand dimensions')
  lo,hi=map(dec,d['bounds']);ck(all(lo<dec(v)<hi for v in d['points']),'open physical nodes')
  if route!='B':
   panel,side=int(path[4]),int(path[5]);ck(d['bounds']==st['panels'][panel] and d['side']==side and (panel,side) not in st['seen'],'actual squared panel arguments');st['seen'].add((panel,side))
   ck(len(d['jacobians'])==n and all(dec(j)>0 for j in d['jacobians']) and len(d['return'])==components,'positive squared Jacobians/returns')
  else:
   label=path[4];ck(label not in st['leaves'] and len(d['K'])==len(d['G'])==len(d['error'])==components and all(dec(e)>=0 for e in d['error']),'actual adaptive panel vector')
   if label[1:].isdigit():ck(d['bounds']==st['panels'][int(label[1:])],'initial adaptive panel bound join')
   st['leaves'][label]={'bounds':d['bounds'],'K':d['K'],'error':d['error']}
  if path[0]=='point':
   for e in d['kernelOperands']:
    ck(set(e)=={'qi','qo','qh','qs','used_qs','A1','A2'} and eq(e['qs'],e['used_qs']) and root(e['qi'],True) and root(e['qo'],True) and root(e['qh'],True) and root(e['qs']),'actual outgoing/reflected depth record')
  continue
 if subkind=='refinement':
  parent=d['replacedLeaf'];left,right=d['children'];ck(left==parent+'L' and right==parent+'R','adaptive child labels')
  old=st['leaves'][parent];ll=st['leaves'][left];rr=st['leaves'][right]
  ck(eq(ll['bounds'][0],old['bounds'][0]) and eq(rr['bounds'][1],old['bounds'][1]) and eq(ll['bounds'][1],rr['bounds'][0]) and d['parentError']==old['error'] and d['childErrors']==[ll['error'],rr['error']],'actual adaptive replacement operands');del st['leaves'][parent];continue
 if subkind in ['global-initial','global-sum']:
  ck(set(d['activeLeaves'])==set(st['leaves']) and len(d['summedEmpiricalErrors'])==components,'actual active global partition');continue
 if subkind=='return':
  ck(len(d['value'])==components,'complete route return components')
  if route=='B':
   ck(set(d['activeLeaves'])==set(st['leaves']) and d['panelsEvaluated']==st['count'] and all(Decimal(0)<=dec(e)<=dec(d['budgetPerComponent']) for e in d['summedEmpiricalErrors']) and d['physicalVariable'] is True,'completed global empirical criterion')
   for e in d['summedEmpiricalErrors']:
    if dec(e)>maxima['BempiricalError']:maxima['BempiricalError']=dec(e);max_at['BempiricalError']={'record':name,'savedError':e,'budget':d['budgetPerComponent']}
  else:ck(st['seen']=={(i,j) for i in range(len(st['panels'])) for j in [0,1]} and d['panels']==len(st['panels']) and d['order']==int(route[1:3]) and d['squareSubstitution'] is True,'completed squared panel census')
  returns[key]=d;del states[key];continue
 raise AssertionError('unknown record '+name)
db.close()
ck(len(points)==38 and len(Hresults)==37 and len(reuses)==1 and len(comparisons)==75 and len(returns)==263 and not states and len(controls)==3,'all bounded points/routes/comparisons complete')
hc=controls['controls/H-contact'];hp='H/'+hc['exactQCoefficient'].replace('/','_');ck(eq(hc['savedMomentum'],{**momenta[hp]['data'],'assemblyReceipt':momenta[hp]['receipt']}) and eq(hc['baseline'],momenta[hp]['data']['baseline']) and eq(hc['mutated'],momenta[hp]['data']['integral']) and hc['newQuadratureCalls']==0 and hc['completedIntegralReused'] is True,'actual H ablation operand/return reuse')
lc=controls['controls/lower-normal'];ck(eq(lc['baseline'],points[37]['data']['normalMixedAndDirect']['minus'][1]) and eq(lc['wholeDirectValue'],points[37]['data']['routes']['B'][-1]),'actual lower-normal whole value')
rc=controls['controls/reflected-root'];ck(eq(rc['operands']['qs'],rc['operands']['used_qs']) and eq(rc['mutatedOperands']['used_qs'],rc['operands']['qh']) and not eq(rc['operands']['used_qs'],rc['mutatedOperands']['used_qs']),'actual reflected root corruption')
ck(sha(DB)==dbreceipt['sha256'],'SQLite postinspection identity')
record={'status':'ACCEPTED_38_FIXED_INNER_KERNEL_POINTS_NO_PACKET_ACTION','utc':datetime.now(timezone.utc).isoformat(),'root':str(B),'sourceCommit':'f509899b','reviewCommit':'ef63ee63','readyCommit':'ca0969f6','launchCommit':'8cadda1e','allChecksPassed':True,'inspectionCounts':dict(counts),'inspectionChecks':sum(counts.values()),'inspectionKind':'stdlib JSON/hash/opaque SQLite traversal; saved scalar decimal displays inspected only. No scientific payload restoration, symbolic work or quadrature replay.','result':result,'newLiteralZeroReturns':len(zeros),'inheritedLiteralZeroReturns':len(prefixes),'jsonArtifacts':len(arts),'savedCopies':len(copies),'sourcePosthashes':len(post),'sourceSnapshots':len(snap),'additionalIdentityCopies':len(extra),'sqlite':{'receipt':dbreceipt,'records':total_rows,'payloadBytes':payload_bytes,'panelCounts':dict(panels_count),'readonlyIntegrity':'ok'},'savedScalarMaxima':{k:str(v) for k,v in maxima.items()},'maximumRecords':max_at,'controls':controls,'pointComparisons':[{'point':plan[i],'comparisons':comparisons['point/'+str(i)]} for i in range(38)],'HComparisons':[{'key':p,'comparisons':v} for p,v in comparisons.items() if p.startswith('H/')],'resources':{'samples':len(samples),'memoryPeakBytes':max(int(s['memory.peak']) for s in samples),'hostAvailableMinimumBytes':min(s['hostAvailableBytes'] for s in samples),'eventsMax':dict(events),'swapBytes':0,'affinity':a['affinity'],'duration':duration['actual'],'interpretation':'Actual cgroup event counts retained. No process/cache split inferred.'},'strictStderrBytes':0,'stdoutChecksByteIdentity':True,'limitations':['Only38fixed inner points,37distinct H requests and3narrow sensitivity controls. No uniform quadrature error or outer packet value.','Embedded adaptive errors/cross-route differences are empirical; shared analytic conventions not independently validated by route agreement.','Full summand units, local actions/tails, actual outer Fourier requests, collision-aware full actions and full addressed controls remain.','No finite inverse, scattering/current/loss, drain, calibration or sweep.'],'next':'Continue bounded packet preparation from accepted immutable inner/Fourier/preflight results. Derive and assess only missing full-summand units/local action and outer composition prerequisites; no completed integral replay or routine permission question.'}
with (M/(P+'_result_record.json')).open('x') as f:json.dump(record,f,indent=2);f.write('\n')
print(json.dumps({k:record[k] for k in ['status','inspectionChecks','newLiteralZeroReturns','sqlite','savedScalarMaxima','resources']},indent=2))
