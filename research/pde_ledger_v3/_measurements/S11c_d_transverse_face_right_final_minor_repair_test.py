"""Stdlib synthetic NEWOBJ test and source identity audit; never load science."""
from pathlib import Path
from datetime import datetime,timezone
import ast,hashlib,io,json,pickle,pickletools,types
R=Path('/var/projects/toy_physics');M=R/'research/pde_ledger_v3/_measurements';P='S11c_d_transverse_face_right_final';D=R/'_scratch/s11c/s11c-d-transverse-face-fixed-point-20260929/right-final-build-review'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def pin(p):p=Path(p);return dict(path=str(p),canonicalPath=str(p.resolve()),bytes=p.stat().st_size,sha256=sha(p))
def read(p):return json.loads(Path(p).read_text())
review=read(M/(P+'_review_record.json'));worker=M/(P+'.py');baseline=D/'packet'/worker.relative_to(R)
assert sha(baseline)==review['workerSha256'];before=ast.parse(baseline.read_text());after=ast.parse(worker.read_text())
functions=lambda t:{n.name:ast.dump(n,include_attributes=False) for n in ast.walk(t) if isinstance(n,ast.FunctionDef)}
a,b=functions(before),functions(after);changed=sorted(n for n in b if b[n]!=a.get(n))
assert set(changed)=={'stored_inequality_class','__new__','find_class','main'},changed
adapter=next(n for n in after.body if isinstance(n,ast.FunctionDef) and n.name=='stored_inequality_class');ns={}
exec(compile(ast.Module(body=[adapter],type_ignores=[]),'<actual-codec-adapter-synthetic>','exec'),ns)
class SyntheticRelation:
 def __new__(cls,lhs,rhs,evaluate=True):
  if evaluate:return lhs>rhs
  obj=object.__new__(cls);obj.args=(lhs,rhs);return obj
 def __getnewargs__(self):return self.args
 def __getstate__(self):return None
 def __str__(self):return str(self.args[0])+' > '+str(self.args[1])
class SyntheticCodec(pickle.Unpickler):
 def find_class(self,module,name):
  assert module==__name__ and name=='SyntheticRelation'
  return ns['stored_inequality_class'](SyntheticRelation)
def restore_synthetic(lhs,rhs):
 original=SyntheticRelation(lhs,rhs,evaluate=False);payload=pickle.dumps(original,protocol=4)
 assert any(op.name=='NEWOBJ' for op,_,_ in pickletools.genops(payload))
 default=pickle.loads(payload);restored=SyntheticCodec(io.BytesIO(payload)).load()
 assert type(default) is bool and default==(lhs>rhs)
 assert type(restored) is SyntheticRelation and restored.args==(lhs,rhs)
 assert str(restored)==str(original)
 return dict(arguments=[lhs,rhs],defaultLoaderValue=default,preservedText=str(restored),originalClassReturned=True,newobjPresent=True)
results=[restore_synthetic(a,b) for a,b in [(1,0),(0,1),(0,0),(7,-3)]]
# A changed argument stays changed; the strict JSON guard still observes it.
assert str(SyntheticCodec(io.BytesIO(pickle.dumps(SyntheticRelation(2,0,evaluate=False),protocol=4))).load())!='1 > 0'
for name in ('reuse_summary','json_bytes','strict_json','json_differences','epsilon_coefficients','grade','extract','faces','bound_end','saved_seed','same','op','configure_prior','restore_suboperation','containment','verify_failed_attempt'):
 assert a[name]==b[name],name
assert sha(M/(P+'_guard.py'))==review['guardSha256']
assert sha(M/(P+'_inputs.json'))==review['inputManifestSha256']
assert sha(M/(P+'_scope.md'))==review['scopeSha256']
assert sha(M/(P+'_authorization.json'))==review['authorizationSha256']
assert sha(R/'scripts/s11c_guarded_run.py')==review['sharedGuardSha256']
assert sha(M/'S11c_d_end_normalization_run.py')==review['supervisorSha256']
find=next(n for n in ast.walk(after) if isinstance(n,ast.FunctionDef) and n.name=='find_class')
branches=[n for n in ast.walk(find) if isinstance(n,ast.Return) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Name) and n.value.func.id=='stored_inequality_class'];assert len(branches)==1
result=dict(status='SYNTHETIC_NEWOBJ_AND_SOURCE_IDENTITY_PASS',recordedUtc=datetime.now(timezone.utc).isoformat(),workerSha256=sha(worker),reviewedWorkerSha256=sha(baseline),testSource=pin(__file__),changedFunctions=changed,onlyCodecAndAuthorizationChanged=True,syntheticCases=results,changedArgumentsRemainUnequal=True,strictSummaryAndAllScientificFunctionsUnchanged=True,guardManifestScopeOriginalAuthorizationUnchanged=True,scientificImports=0,scientificPayloadRestorations=0,independentBuildClearance=False)
with (M/(P+'_minor_repair_checks.json')).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
