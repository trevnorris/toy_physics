#!/usr/bin/env python3
"""Read-only Q9 dependency trace; standard-library AST/hash inspection only.

Does not import a CAS engine, reconstruct symbolic exports, calculate a nullspace,
change active inputs, or certify the proposed Q9 correction.
"""
import ast,hashlib,json,pickletools,re,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
STORE=REPO/'_scratch/s11c/s11c-q9-impact-20260916';STORE.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(ROOT/'scripts'))
from S11c_d_output_codec import decoded_lines
PRODUCTION=REPO/'_scratch/s11c/s11c-momentum-domain-action-20260916/production'

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def tree(path):return ast.parse(path.read_text())
def assigned(module,name):return next(n.value for n in module.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id==name for t in n.targets))
def function(module,name):return next(n for n in module.body if isinstance(n,ast.FunctionDef) and n.name==name)
def strings(node):return [n.value for n in ast.walk(node) if isinstance(n,ast.Constant) and isinstance(n.value,str)]
def rows(path):
 value=assigned(tree(path),'_LEDGER');return {ast.literal_eval(k):{ast.literal_eval(a):b for a,b in zip(v.keys,v.values)} for k,v in zip(value.keys,value.values)}
def serial(node):return ast.dump(node,include_attributes=False)
def write(path,value):path.write_text(json.dumps(value,indent=2)+'\n')

producer=ROOT/'scripts/S11_stray_longitudinal_sympy_audit.py';current=tree(producer)
head_text=subprocess.check_output(['git','show','HEAD:'+str(producer.relative_to(REPO))],cwd=REPO,text=True);head=ast.parse(head_text)
old_functions={n.name:serial(n) for n in head.body if isinstance(n,(ast.FunctionDef,ast.ClassDef))};new_functions={n.name:serial(n) for n in current.body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
changed=[k for k in old_functions.keys()|new_functions.keys() if old_functions.get(k)!=new_functions.get(k)]
assert changed==['compute_q9']
qnames=ast.literal_eval(next(n.iter for n in ast.walk(function(current,'emit_q9')) if isinstance(n,ast.For)))+('PD_TERM',)
qkeys={name.lower()+'_d'+str(n) for name in qnames for n in (2,3,4,5)}
exports={};allrows={}
for name in ['S11','S11b','S11c_a','S11c_b','S11c_c1','S11c_c2']:
 path=ROOT/f'scripts/{name}_exports.py';values=rows(path);allrows[name]=values
 present=sorted(qkeys&values.keys());baseline=allrows['S11']
 identical=[k for k in present if serial(values[k]['value'])==serial(baseline[k]['value'])]
 exports[name]={'path':str(path.relative_to(ROOT)),'sha256':digest(path),'rowCount':len(values),'q9FamilyKeys':present,'q9ValueIdentityToS11':len(identical),'rawCopiedQ9Family':len(present)==len(identical)}

manifest={}
for label,file in [('c1','S11c_c1_bulk_closure_sympy_audit.py'),('c2','S11c_c2_selfenergy_fold_sympy_audit.py'),('d','S11c_d_mixing_scattering_sympy_audit.py')]:
 path=ROOT/'scripts'/file;module=tree(path)
 if label!='d':keys=ast.literal_eval(assigned(module,'IMPORT_KEYS'))
 else:
  assignments=[n for n in module.body if isinstance(n,(ast.Assign,ast.AugAssign)) and ((isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id in ('CLOSED_KEYS','IMPORT_KEYS','DIMENSION_CARRIERS') for t in n.targets)) or (isinstance(n,ast.AugAssign) and isinstance(n.target,ast.Name) and n.target.id=='IMPORT_KEYS'))]
  namespace={};exec(compile(ast.Module(body=assignments,type_ignores=[]),'<literal-import-manifest>','exec'),namespace);keys=namespace['IMPORT_KEYS']
 manifest[label]={'path':str(path.relative_to(ROOT)),'sha256':digest(path),'keys':list(keys),'q9Intersection':sorted(qkeys&set(keys))}
 assert not manifest[label]['q9Intersection']

# Execute only the actual pure S11c-d binder against opaque, unevaluated AST rows.
# Delete every carried Q9 family row and separately mutate an actual closed row.
fold={}
for name in ('S11c_b','S11c_c1','S11c_c2'):
 for key,value in allrows[name].items():
  fold[key]={k:(ast.literal_eval(v) if isinstance(v,ast.Constant) else serial(v)) for k,v in value.items()}
dmodule=tree(ROOT/'scripts/S11c_d_mixing_scattering_sympy_audit.py');binder=function(dmodule,'bind');namespace={'IMPORT_KEYS':tuple(manifest['d']['keys']),'CLOSED_KEYS':ast.literal_eval(assigned(dmodule,'CLOSED_KEYS'))};exec(compile(ast.Module(body=[binder],type_ignores=[]),'<actual-S11c-d-bind>','exec'),namespace)
selected=namespace['bind'](fold);removed=namespace['bind']({k:v for k,v in fold.items() if k not in qkeys});assert selected==removed
mutated=dict(fold);key=namespace['CLOSED_KEYS'][0];mutated[key]=dict(fold[key],value='Q9_DEPENDENCY_NEGATIVE_CONTROL');assert namespace['bind'](mutated)!=selected

# Direct Q9 action edge is confined to the XFORM_EXTRA branch of package_build.
build=function(current,'package_build');mapping=next(n.value for n in ast.walk(build) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='stiffness_map' for t in n.targets))
pd_packages=[ast.literal_eval(k) for k,v in zip(mapping.keys,mapping.values) if any(isinstance(n,ast.Name) and n.id=='pd_density' for n in ast.walk(v))];assert pd_packages==['XFORM_EXTRA']
q7_uses=sum(isinstance(n,ast.Name) and n.id=='q9_data' for n in ast.walk(function(current,'q7_objects')));assert q7_uses==0

# Current raw-row caches contain only the selected import rows, not carried Q9
# fields. Parse pickle string opcodes without unpickling or importing SymPy.
cache_checks=[]
for checkpoint,filenames in [('S11c_d_reduced_action_source_checkpoint.json',['reduced-action.pickle','actions.pickle']),('S11c_d_momentum_action_checkpoint.json',['bound-momentum.pickle'])]:
 cp=ROOT/'_measurements'/checkpoint;c=json.loads(cp.read_text());directory=Path(c['runDirectory'])
 for name in filenames:
  p=directory/name;sha=digest(p);assert sha==c['artifacts'][name]['sha256']
  literal_strings={v for op,v,_ in pickletools.genops(p.read_bytes()) if isinstance(v,str)}
  hits=sorted(qkeys&literal_strings);assert not hits
  cache_checks.append({'checkpoint':str(cp.relative_to(ROOT)),'checkpointSha256':digest(cp),'path':str(p),'sha256':sha,'q9RowNameOpcodes':hits,'note':'String-opcode census supports the source-level bind trace; it is not a standalone algebraic independence proof.'})
source_cp=json.loads((ROOT/'_measurements/S11c_d_reduced_action_source_checkpoint.json').read_text());directory=Path(source_cp['runDirectory']);imports={}
for line in decoded_lines(directory/'reduction.out'):
 tag,sep,body=line.partition(': ')
 if tag in ('PY_S11CD_IMPORT_LOOKUPS','PY_S11CD_IMPORT_CLOSURE'):
  literals=strings(ast.parse(body,mode='eval'));imports[tag]={'sha256':hashlib.sha256(line.encode()).hexdigest(),'q9Intersection':sorted(qkeys&set(literals)),'strings':literals if tag.endswith('LOOKUPS') else None}
assert len(imports)==2 and all(not v['q9Intersection'] for v in imports.values())
assert set(imports['PY_S11CD_IMPORT_LOOKUPS']['strings'])==set(manifest['d']['keys'])

# Capture the active stage and its full source-hash boundary once; no polling.
active=json.loads((PRODUCTION/'active.json').read_text());workers=json.loads((PRODUCTION/'complete/workers.json').read_text());preflight=json.loads((PRODUCTION/'complete/preflight.json').read_text());pins=preflight['sourceFiles']
for n,h in pins.items():assert digest(ROOT/n)==h==digest(PRODUCTION/'complete/source'/n),n
assert str(producer.relative_to(ROOT)) not in pins
pinned_exports={n:h for n,h in pins.items() if n.endswith('_exports.py')}

# Collect every literal inherited-key and export-merge access for human review.
accesses={}
for name in ['S11b_interface_coupling_law_sympy_audit.py','S11c_a_interface_geometry_sympy_audit.py','S11c_b_brane_operator_sympy_audit.py']:
 path=ROOT/'scripts'/name;module=tree(path);records=[]
 for n in ast.walk(module):
  if isinstance(n,ast.Subscript) and isinstance(n.value,ast.Name) and n.value.id=='INCOMING_LEDGER':records.append({'line':n.lineno,'operation':ast.unparse(n)})
  if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and isinstance(n.func.value,ast.Name) and n.func.value.id=='INCOMING_LEDGER':records.append({'line':n.lineno,'operation':ast.unparse(n)})
  if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id in ('inherited','inherited_symbol','bind_additional_inherited') and n.args:records.append({'line':n.lineno,'operation':ast.unparse(n.args[0])})
 accesses[name]={'sha256':digest(path),'accesses':sorted(records,key=lambda x:x['line'])}

result={'status':'DEPENDENCY_TRACE_NO_S11CD_NUMERICAL_RERUN_REQUIRED_FOR_Q9_ONLY_PATCH','checkedUtc':datetime.now(timezone.utc).isoformat(),'scope':'Read-only standard-library AST, hash, raw-key binder and pickle-opcode trace; no CAS or Q9 D3-D5 validation.','instrumentSha256':digest(Path(__file__)),'producer':{'path':str(producer.relative_to(ROOT)),'workingSha256':digest(producer),'headSha256':hashlib.sha256(head_text.encode()).hexdigest(),'changedDefinitions':changed},'q9FamilyQuantities':list(qnames),'q9FamilyKeyCount':len(qkeys),'q9ActionPackages':pd_packages,'q7DirectQ9ReadCount':q7_uses,'exports':exports,'manifests':manifest,'actualDBinderControl':{'removedQ9Keys':len(qkeys),'selectedRows':len(selected),'sameBindingsAfterRemoval':True,'changedConsumedRootDetected':True},'upstreamLedgerAccesses':accesses,'reducedCacheImportEvidence':imports,'cacheChecks':cache_checks,'activeJob':{'directory':str(PRODUCTION),'active':active,'workers':workers,'pinnedSourceCount':len(pins),'allCurrentFrozenHashesMatch':True,'producerIsPinned':False,'pinnedExports':pinned_exports,'stderrBytes':(PRODUCTION/'momentum_construct.stderr').stat().st_size},'disposition':{'continueActiveS11cdJob':True,'recomputeQ9AndCheckD3D5':'Owned by upstream repair session after its review; counts alone are insufficient.','recheckExtraAction':'Revalidate PD_TERM and XFORM_EXTRA for each dimension; do not presume the parity-odd density changes before comparison.','refreshAccumulatedExports':'Q9 family values are carried unchanged through S11/S11b/S11c-a/S11c-b; these containers require a controlled row/provenance refresh after repair acceptance.','retainS11cdOperands':'No identified Q9 value edge into the closed slab/kernel, energy-current input or bound numerical operators. Preserve caches and compare actually consumed roots/closure before any later rebase.','activeInputImmutability':'Do not replace pinned S11c-b/c1/c2 exports while production runs; even an unused-row change alters their whole-file hashes and causes a provenance guard failure.','notEstablished':'No D3-D5 Q9 correctness claim, no blanket downstream correctness claim, no accepted repair or new physical result.'}}
write(STORE/'dependency-trace.json',result);write(ROOT/'_measurements/S11c_d_q9_dependency_checkpoint.json',result)
print(json.dumps({'status':result['status'],'q9FamilyKeys':len(qkeys),'exportRowsCarryingFamily':{n:len(v['q9FamilyKeys']) for n,v in exports.items()},'manifestCounts':{n:len(v['keys']) for n,v in manifest.items()},'binder':result['actualDBinderControl'],'pinnedSourceCount':len(pins),'activeStage':active['stage'],'workerCount':len(workers.get('workers',[]))},indent=2))
