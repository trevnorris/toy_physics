"""Source/JSON/byte audits plus synthetic exact arithmetic. No scientific imports or pickle loads."""
from pathlib import Path
from datetime import datetime, timezone
from fractions import Fraction
from math import comb
import ast, hashlib, json, tempfile, types, traceback
R=Path('/var/projects/toy_physics');M=R/'research/pde_ledger_v3/_measurements';P='S11c_d_transverse_face_right_final';OLD='S11c_d_transverse_face_right_finish'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def require(x,message):
 if not x:raise AssertionError(message)
checks={}
def check(name,value):require(value,name);checks[name]=bool(value)
source=(M/(P+'.py')).read_text();tree=ast.parse(source);oldtree=ast.parse((M/(OLD+'.py')).read_text());gtree=ast.parse((M/(P+'_guard.py')).read_text());spec=read(M/(P+'_inputs.json'))
for key,pin in spec['inputs'].items():check('input-'+key,sha(pin['path'])==pin['sha256'] and Path(pin['path']).stat().st_size==pin['bytes'])
functions=lambda t:{n.name:ast.dump(n,include_attributes=False) for n in ast.walk(t) if isinstance(n,ast.FunctionDef)}
a,b=functions(oldtree),functions(tree)
changed=sorted(name for name in b if b[name]!=a.get(name))
expected={'json_bytes','save','strict_json','json_differences','reject_constant','epsilon_coefficients','multiply','visit','reuse_summary','extract','grade','construct','main','verify_failed_attempt'}
check('bounded-function-delta',set(changed)==expected)
for name in ('faces','saved_seed','bound_end','show','same','bind','loaded_control','numeric','numerical_face_legs','face_address_census','op','restore_suboperation','configure_prior','configure_resume','containment'):
 check('unchanged-'+name,a[name]==b[name])
check('no-executable-Poly-call',not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='Poly' for n in ast.walk(tree)))
# Exact stand-in expression trees: no SymPy or source payload. Evaluate independently
# at rational values after the actual new algorithm returns its coefficients.
class E:
 def __init__(self,kind,*args):self.kind,self.args=kind,args
 def __hash__(self):return hash((self.kind,self.args))
 def __eq__(self,x):return isinstance(x,E) and (self.kind,self.args)==(x.kind,x.args)
 def __add__(self,x):return E('add',self,wrap(x))
 def __radd__(self,x):return wrap(x)+self
 def __mul__(self,x):return E('mul',self,wrap(x))
 def __rmul__(self,x):return wrap(x)*self
 def __pow__(self,x):return E('pow',self,wrap(x))
 def has(self,e):return self==e or any(isinstance(x,E) and x.has(e) for x in self.args)
 @property
 def is_Add(self):return self.kind=='add'
 @property
 def is_Mul(self):return self.kind=='mul'
 @property
 def is_Pow(self):return self.kind=='pow'
 @property
 def is_Integer(self):return self.kind=='number' and self.args[0].denominator==1
 @property
 def is_nonnegative(self):return self.kind=='number' and self.args[0]>=0
 @property
 def exp(self):return self.args[1]
 @property
 def base(self):return self.args[0]
 @property
 def func(self):return self.kind
 def __int__(self):return int(self.args[0])
def wrap(x):return x if isinstance(x,E) else E('number',Fraction(x))
def evaluate(x,values):
 if x.kind=='number':return x.args[0]
 if x.kind=='symbol':return values[x.args[0]]
 if x.kind=='add':return sum(evaluate(y,values) for y in x.args)
 if x.kind=='mul':
  out=Fraction(1)
  for y in x.args:out*=evaluate(y,values)
  return out
 if x.kind=='pow':return evaluate(x.base,values)**int(x.exp)
 raise AssertionError('unsupported independent synthetic evaluation')
selected=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ('epsilon_coefficients','json_bytes','strict_json','json_differences')]
ns=dict(require=require,json=json)
exec(compile(ast.Module(body=selected,type_ignores=[]),'<actual-helper-synthetic-tests>','exec'),ns)
e=E('symbol','epsilon'); y=E('symbol','eta'); z=E('symbol','sigma'); sp=types.SimpleNamespace(S=types.SimpleNamespace(Zero=wrap(0),One=wrap(1)))
cases=[('zero',wrap(0),(0,0,0)),('constant',wrap(7),(7,0,0)),('linear',e,(0,1,0)),('mixed-degrees',3+2*e+7*e**2+11*e**3,(3,2,7)),('cross-terms',(2+3*e+5*e**2)*(7+11*e+13*e**2),(14,43,94)),('high-power',(2+3*e)**9,(2**9,9*2**8*3,comb(9,2)*2**7*9)),('epsilon-high-only',e**20,(0,0,0)),('opaque-grades',(y+z*e+y*z*e**2)**2,(Fraction(4,9),Fraction(20,21),Fraction(25,49)+Fraction(40,63)))]
for name,expression,expected_values in cases:
 result=ns['epsilon_coefficients'](expression,e,sp)
 actual=tuple(evaluate(x,{'eta':Fraction(2,3),'sigma':Fraction(5,7)}) for x in result['coefficients'])
 check('exact-synthetic-'+name,actual==expected_values and all(not x.has(e) for x in result['coefficients']))
for name,expression in [('denominator',(1+e)**-1),('function',E('sin',e))]:
 try:ns['epsilon_coefficients'](expression,e,sp)
 except ValueError:checks['reject-'+name]=True
 else:raise AssertionError('unsupported dependence was accepted')
# Exercise actual summary method with normal JSON only, checking persistence on failure.
method=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='reuse_summary')
class IntegrityError(RuntimeError):pass
def route(p):p=Path(p);return dict(path=str(p),canonicalPath=str(p.resolve()),bytes=p.stat().st_size,sha256=sha(p))
def publish(p,b):p.parent.mkdir(parents=True,exist_ok=True);p.open('xb').write(b)
def save(p,value):publish(p,ns['json_bytes'](value))
ns.update(Path=Path,route=route,publish=publish,save=save,IntegrityError=IntegrityError,traceback=traceback)
exec(compile(ast.Module(body=[method],type_ignores=[]),'<actual-summary-method-synthetic-tests>','exec'),ns)
with tempfile.TemporaryDirectory(prefix='s11cd-final-synthetic-') as tmp:
 root=Path(tmp)
 def trial(name,old,fresh,expect_error=False):
  prior=root/(name+'-prior.json');save(prior,old);out=root/name
  j=types.SimpleNamespace(out=out,spec=dict(priorSummaryFiles=['summary.json'],priorFiles={'summary.json':'prior'},inputs={'prior':route(prior)}),summary_reuse=[],blob=lambda x:dict(synthetic=True))
  try:ns['reuse_summary'](j,'summary.json',fresh)
  except IntegrityError:require(expect_error,'unexpected strict summary failure')
  else:require(not expect_error,'missed strict summary failure')
  return out
 out=trial('tuple-keys',{'3':[1,2]}, {3:(1,2)})
 check('same-encoder-roundtrip',read(out/'summary-joins/summary.json.join.json')['actualArgumentJoin'] is True)
 out=trial('different',{'x':[1,2],'a/b~':0},{'x':[1,3],'a/b~':1},True)
 differences=read(out/'summary-comparisons/summary.json/differences.json')
 check('persisted-differing-paths',{d['path'] for d in differences}=={'/x/1','/a~1b~0'})
 check('both-mismatch-summaries-exist',(out/'summary-comparisons/summary.json/saved.json').is_file() and (out/'summary-comparisons/summary.json/fresh.json').is_file())
 check('no-published-success-summary-after-mismatch',not (out/'summary.json').exists())
 for name,value in [('nan',float('nan')),('infinity',float('inf'))]:
  try:ns['json_bytes']({'bad':value})
  except ValueError:checks['reject-fresh-'+name]=True
  else:raise AssertionError('nonfinite accepted')
 for raw in ('{"bad":NaN}','{"bad":Infinity}'):
  try:ns['strict_json'](raw)
  except ValueError:checks['reject-saved-'+raw]=True
  else:raise AssertionError('nonfinite accepted')
# The real wrapper must restore 52..54 and compare input55 before extracting55.
extract=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='extract')
def scalar_trial(mutate=False):
 class J:
  carrier_entry_number=52
  same=staticmethod(lambda x,y:x==y)
  def __init__(self):self.values={};self.restored=[]
  def resume_has(self,n):return n in self.values
  def restore_suboperation(self,n):self.restored.append(n);return self.values[n]
 j=J();calls=[];published=[]
 for i in range(52,56):
  j.values['carrier-entry-%06d-input'%i]=dict(entry=i,epsilon='epsilon')
  if i<55:j.values['carrier-entry-%06d-return'%i]=dict(entry=i,epsilon='epsilon',coefficient=100+i)
 if mutate:j.values['carrier-entry-000055-input']['entry']='wrong'
 def exact(entry,epsilon,sp):calls.append(entry);return dict(coefficients=(0,0,100+entry),method='SYNTHETIC')
 env=dict(journal=j,e='epsilon',sp=None,require=require,epsilon_coefficients=exact,evidence=lambda n,v:published.append((n,v)))
 exec(compile(ast.Module(body=[extract],type_ignores=[]),'<actual-wrapper-synthetic-tests>','exec'),env)
 try:values=[env['extract'](i) for i in range(52,57)]
 except AssertionError:require(mutate and calls==[],'mutation did not stop before extraction');return dict(refused=True,calls=calls)
 require(not mutate and values==[152,153,154,155,156] and calls==[55,56],'completed coefficient replay')
 return dict(calls=calls,restored=j.restored)
wrapper=dict(success=scalar_trial(),mutation=scalar_trial(True))
# Guard changes are path/ordinal/final-attempt authorization only; core controls unchanged.
oldguard=ast.parse((M/(OLD+'_guard.py')).read_text());ga,gb=functions(oldguard),functions(gtree)
check('only-duration-identity-guard-change',[k for k in gb if gb[k]!=ga[k]]==['validate_duration'])
check('shared-guard-unchanged',sha(R/'scripts/s11c_guarded_run.py')=='9194b3e04b43d345828082e6312b0a2d7a9293f6296153f9c6ccea2defb931d2')
raw=(M/'S11c_d_sympy_builder_report.md').read_bytes();check('protected-builder-suffix',hashlib.sha256(raw[raw.index(b'## Retained user-approved solver/export contract'):]).hexdigest()=='f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2')
# Every execution6 file is checked as bytes; zero objects restored.
verify=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='verify_failed_attempt')
ns.update(route=route)
exec(compile(ast.Module(body=[verify],type_ignores=[]),'<actual-metadata-preservation-check>','exec'),ns)
failed=ns['verify_failed_attempt'](spec)
result=dict(status='SOURCE_JSON_HASH_AND_SYNTHETIC_ARITHMETIC_ONLY',recordedUtc=datetime.now(timezone.utc).isoformat(),workerSha256=sha(M/(P+'.py')),guardSha256=sha(M/(P+'_guard.py')),manifestSha256=sha(M/(P+'_inputs.json')),testSourceSha256=sha(__file__),allChecksPassed=all(checks.values()),inputPins=len(spec['inputs']),checks=len(checks),changedFunctions=changed,nonInputChecks={k:v for k,v in checks.items() if not k.startswith('input-')},syntheticWrapper=wrapper,failedAttemptPreservation=failed,scientificImports=0,payloadRestorations=0,workerImportedOrLaunched=False,independentClearance=False)
with (M/(P+'_static_checks.json')).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({k:result[k] for k in ('status','allChecksPassed','inputPins','checks','syntheticWrapper','failedAttemptPreservation')},indent=2))
