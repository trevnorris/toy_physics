"""Static pickle opcode census: no imports of scientific packages or unpickling."""
import collections,hashlib,json,pickletools,signal,resource,time
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path('/var/projects/toy_physics');M=ROOT/'research/pde_ledger_v3/_measurements'
manifest=json.loads((M/'S11c_d_transverse_face_saved_reader_continue_inputs.json').read_text())
resource.setrlimit(resource.RLIMIT_AS,(512*1024**2,512*1024**2))
signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('metadata census budget')));signal.alarm(55)
start=time.monotonic();results=[];seen=set();MARK=object();UNKNOWN=object()
for key in manifest['priorRestoreKeys']:
 rec=manifest['inputs'][key];path=Path(rec['path'])
 if rec['sha256'] in seen:continue
 seen.add(rec['sha256']);raw=path.read_bytes();assert len(raw)==rec['bytes'] and hashlib.sha256(raw).hexdigest()==rec['sha256']
 stack=[];memo={};globals_found=set();counts=collections.Counter();unsupported=[]
 for op,arg,pos in pickletools.genops(raw):
  n=op.name;counts[n]+=1
  if n in ('PROTO','FRAME'):continue
  if n=='MARK':stack.append(MARK)
  elif n in ('SHORT_BINUNICODE','BINUNICODE','BINUNICODE8','UNICODE','STRING','BINSTRING','SHORT_BINSTRING'):stack.append(arg)
  elif n in ('PUT','BINPUT','LONG_BINPUT'):memo[int(arg)]=stack[-1]
  elif n=='MEMOIZE':memo[len(memo)]=stack[-1]
  elif n in ('GET','BINGET','LONG_BINGET'):stack.append(memo[int(arg)])
  elif n=='GLOBAL':
   module,name=arg.split(' ',1);globals_found.add((module,name));stack.append(UNKNOWN)
  elif n=='STACK_GLOBAL':
   name,module=stack.pop(),stack.pop()
   assert isinstance(name,str) and isinstance(module,str),(path,pos,'unresolved GLOBAL')
   globals_found.add((module,name));stack.append(UNKNOWN)
  elif n=='POP':stack.pop()
  elif n=='DUP':stack.append(stack[-1])
  elif n=='POP_MARK':
   while stack.pop() is not MARK:pass
  elif n in ('TUPLE','LIST','DICT','FROZENSET'):
   while stack.pop() is not MARK:pass
   stack.append(UNKNOWN)
  elif n in ('APPENDS','SETITEMS','ADDITEMS'):
   while stack.pop() is not MARK:pass
  elif n in ('APPEND','BUILD'):stack.pop()
  elif n=='SETITEM':stack.pop();stack.pop()
  elif n in ('TUPLE1','TUPLE2','TUPLE3'):
   for _ in range(int(n[-1])):stack.pop()
   stack.append(UNKNOWN)
  elif n in ('REDUCE','NEWOBJ'):
   stack.pop();stack.pop();stack.append(UNKNOWN)
  elif n=='NEWOBJ_EX':
   stack.pop();stack.pop();stack.pop();stack.append(UNKNOWN)
  elif n=='STOP':assert len(stack)==1;break
  elif n in ('NONE','NEWTRUE','NEWFALSE','INT','BININT','BININT1','BININT2','LONG','LONG1','LONG4','FLOAT','BINFLOAT','BINBYTES','SHORT_BINBYTES','BINBYTES8','BYTEARRAY8','EMPTY_DICT','EMPTY_LIST','EMPTY_TUPLE','EMPTY_SET'):stack.append(UNKNOWN)
  else:raise ValueError(('unsupported opcode',n,str(path),pos))
 results.append(dict(keys=[k for k in manifest['priorRestoreKeys'] if manifest['inputs'][k]['sha256']==rec['sha256']],route=rec,globals=[list(x) for x in sorted(globals_found)],setOpcodes={n:counts[n] for n in ('EMPTY_SET','ADDITEMS','FROZENSET')},setSliceGlobals=[list(x) for x in sorted(globals_found) if x in [('builtins','set'),('builtins','frozenset'),('builtins','slice')]],opcodeCount=sum(counts.values())))
signal.alarm(0)
record=dict(status='STATIC_SERIALIZATION_METADATA_ONLY',recordedUtc=datetime.now(timezone.utc).isoformat(),manifestSha256=hashlib.sha256((M/'S11c_d_transverse_face_saved_reader_continue_inputs.json').read_bytes()).hexdigest(),uniquePickles=len(results),wallSeconds=time.monotonic()-start,scientificImports=0,scientificRestorations=0,reducersExecuted=0,results=results)
p=M/'S11c_d_transverse_face_saved_reader_continue_codec_census.json'
with p.open('x') as f:json.dump(record,f,indent=2);f.write('\n')
print(json.dumps({'uniquePickles':len(results),'wallSeconds':record['wallSeconds'],'setOrSliceFiles':[{k:v for k,v in r.items() if k in ['keys','setOpcodes','setSliceGlobals']} for r in results if any(r['setOpcodes'].values()) or r['setSliceGlobals']],'allGlobals':sorted({'.'.join(g) for r in results for g in r['globals']})},indent=2))
