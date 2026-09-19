#!/usr/bin/env python3
"""One midpoint refinement of the actual saved finite-frequency contour."""
import argparse,ast,copy,json,os,resource,shutil,signal,subprocess,sys,time
from pathlib import Path
import numpy as np
import S11c_d_frequency_contour as c
f=c.f;end=c.end;COUNT=32;WORKERS=4
PLAN=f.M/'S11c_d_frequency_contour_refine_plan.md'

def frequencies():return [c.CENTER+c.RADIUS*np.exp(2j*np.pi*i/COUNT) for i in range(COUNT)]

def load(base):
 data=c.load(base);path=f.M/'S11c_d_frequency_contour_checkpoint.json';cp=json.loads(path.read_text())
 f.require(cp['status']=='VALIDATED_NUMERICAL_CONTOUR_DIAGNOSTIC','accepted sixteen-point contour')
 for name,h in cp['sourceFiles'].items():f.require(data['pins'][name]==h==f.digest(f.ROOT/name),'unchanged accepted contour sources')
 old=Path(cp['runDirectory']);f.require(f.digest(old/'checks.json')==cp['checksSha256'],'accepted complete contour inventory')
 checks=json.loads((old/'checks.json').read_text());inventory=json.loads((old/'point-inventory.json').read_text())
 f.require([v['index'] for v in inventory]==list(range(16)),'all sixteen accepted points')
 data.update(oldRun=old,oldChecks=checks,oldInventory=inventory)
 data['operands'][str(old/'checks.json')]=cp['checksSha256']
 for p in (Path(__file__),PLAN,path):data['pins'][str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
 for name,h in data['pins'].items():
  target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/name,target)
 data['manifest'].update(sourceFiles=data['pins'],inputPackets=data['operands'],scope='One angular refinement to32 actual finite contour points, reusing every accepted even point and computing only16 midpoints. Numerical spectral evidence, not certified completeness or a physical bound-pole set.')
 f.save(base/'inputs.json',data['manifest']);f.save(base/'accepted-point-inventory.json',inventory);return data

def adapters():
 original=c.tree(c.worker);node=copy.deepcopy(original);node.name='midpoint_worker'
 loop=next(x for x in ast.walk(node) if isinstance(x,ast.For) and ast.unparse(x.iter)=='range(index, COUNT, WORKERS)');old=copy.deepcopy(loop.iter);loop.iter=ast.parse('range(2*index+1,COUNT,2*WORKERS)',mode='eval').body
 reverse=copy.deepcopy(node);reverse.name=original.name;next(x for x in ast.walk(reverse) if isinstance(x,ast.For) and ast.unparse(x.iter)=='range(2 * index + 1, COUNT, 2 * WORKERS)').iter=old
 f.require(ast.dump(reverse)==ast.dump(original),'whole midpoint-worker scheduling reverse join')
 namespace=dict(vars(c),load=load,frequencies=frequencies,COUNT=COUNT,WORKERS=WORKERS)
 worker=c.compile_function(node,namespace)
 original=c.tree(c.summarize);node=copy.deepcopy(original);node.name='summarize_refinement'
 class Orders(ast.NodeTransformer):
  def visit_Constant(self,n):
   if type(n.value) is int and n.value in (8,16):return ast.copy_location(ast.Constant(value={8:16,16:32}[n.value]),n)
   return n
 node=Orders().visit(node)
 address=next(x for x in ast.walk(node) if isinstance(x,ast.Assign) and len(x.targets)==1 and ast.unparse(x.targets[0])=='i');old_address=copy.deepcopy(address.value);address.value=ast.parse("record['index']",mode='eval').body
 meta=next(x for x in ast.walk(node) if isinstance(x,ast.Call) and ast.unparse(x.func)=='metadata.append');old_meta=copy.deepcopy(meta.args[0]);meta.args[0]=ast.Call(ast.Name('dict',ast.Load()),[meta.args[0]],[ast.keyword('index',ast.Name('i',ast.Load())),ast.keyword('originalIndex',ast.parse("x['index']",mode='eval').body)])
 reverse=copy.deepcopy(node);reverse.name=original.name
 class ReverseOrders(ast.NodeTransformer):
  def visit_Constant(self,n):
   if type(n.value) is int and n.value in (16,32):return ast.copy_location(ast.Constant(value={16:8,32:16}[n.value]),n)
   return n
 reverse=ReverseOrders().visit(reverse)
 next(x for x in ast.walk(reverse) if isinstance(x,ast.Assign) and len(x.targets)==1 and ast.unparse(x.targets[0])=='i').value=old_address
 next(x for x in ast.walk(reverse) if isinstance(x,ast.Call) and ast.unparse(x.func)=='metadata.append').args[0]=old_meta
 f.require(ast.dump(reverse)==ast.dump(original),'whole moment-summary reverse join after order and saved-index addressing')
 summarize=c.compile_function(node,namespace)
 loop=c.compile_function(c.tree(c.close_end_paths),namespace)
 return worker,summarize,loop,{'wholeMidpointWorkerJoin':True,'wholeMomentSummaryJoin':True,'unchangedEndLoopBody':True,'unchangedPointConstructorInverseAndMaps':True}

def original_hashes(data):
 for name,v in data['oldChecks']['artifacts'].items():f.require(f.digest(data['oldRun']/name)==v['sha256'] and (data['oldRun']/name).stat().st_size==v['bytes'],'every accepted artifact unchanged')

def reused_records(data):
 records=[]
 for old in data['oldInventory']:
  p=Path(old['directory'])/'contour-point.pickle';f.require(f.digest(p)==old['sha256'],'unchanged reused point packet');x=f.unpickle(p);i=2*old['index']
  f.require(x['index']==old['index'] and x['frequency']==frequencies()[i],'literal reused index and frequency address')
  records.append(dict(old,index=i,originalIndex=old['index'],reused=True))
 return records

def focused(base):
 data=load(base);original_hashes(data);_,_,_,joins=adapters();records=reused_records(data)
 moments=[np.zeros((645,645),complex) for _ in range(4)];responses=[np.zeros((4,4),complex) for _ in range(4)]
 for record in records:
  x=f.unpickle(Path(record['directory'])/'contour-point.pickle');i=record['index'];delta=frequencies()[i]-c.CENTER;weight=c.RADIUS*np.exp(2j*np.pi*i/COUNT)/16
  for order in range(4):moments[order]+=weight*delta**order*x['balancedInverse'];responses[order]+=weight*delta**order*x['openResponse']
 old=f.unpickle(data['oldRun']/'contour-16.pickle')
 for value,saved in zip(moments+responses,old['inverseMoments']+old['responseMoments']):f.require(np.array_equal(value,saved),'exact saved sixteen-point moments under refined addressing')
 f.require(all(frequencies()[r['index']+1]!=f.unpickle(Path(r['directory'])/'contour-point.pickle')['frequency'] for r in records),'actual wrong-address response')
 owned=[list(range(2*i+1,COUNT,2*WORKERS)) for i in range(WORKERS)];f.require(sorted(sum(owned,[]))==list(range(1,COUNT,2)),'exact once-only midpoint ownership')
 result={'status':'ACCEPTED_FOCUSED_MIDPOINT_ADAPTER','joins':joins,'sourceFiles':data['pins'],'inputPackets':data['operands'],'reusedPoints':16,'newPoints':16,'workerIndices':owned,'exactSavedInverseAndResponseMoments':True,'wrongAddressesRejected':16,'scope':'Source/hash/address and exact saved-moment checks only; no numerical integration or solve repeated.'}
 f.save(base/'checks.json',result);return result

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--run-directory',required=True,type=Path);ap.add_argument('--worker',type=int);ap.add_argument('--focused',action='store_true');args=ap.parse_args();base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
 resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));children=[];logs=[]
 def timeout(*_):raise TimeoutError('bounded midpoint contour budget; retain all completed points')
 signal.signal(signal.SIGALRM,timeout);signal.alarm(900);start=time.monotonic()
 if args.focused:result=focused(base);signal.alarm(0);print(json.dumps(result,indent=2));return
 worker,summarize,loop,joins=adapters()
 if args.worker is not None:result=worker(base,args.worker);signal.alarm(0);print(json.dumps(result,indent=2));return
 data=load(base);original_hashes(data);inventory=reused_records(data);outcomes=[]
 oldcontrol=data['oldRun']/'contour-controls.pickle';shutil.copyfile(oldcontrol,base/'accepted-contour-controls.pickle');f.require(f.digest(base/'accepted-contour-controls.pickle')==data['oldChecks']['artifacts']['contour-controls.pickle']['sha256'],'unchanged accepted nonlinear-pole controls')
 f.save(base/'preflight.json',{'joins':joins,'totalPoints':32,'reusedPoints':16,'newPoints':16,'center':str(c.CENTER),'radius':c.RADIUS,'workers':WORKERS})
 try:
  for i in range(WORKERS):
   directory=base/f'worker-{i}';out=(base/f'worker-{i}.stdout').open('xb');err=(base/f'worker-{i}.stderr').open('xb');logs.extend((out,err));command=[sys.executable,str(Path(__file__).resolve()),'--run-directory',str(directory),'--worker',str(i)]
   child=subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=out,stderr=err);children.append((i,child,directory))
  for i,child,directory in children:
   code=child.wait();outcomes.append({'index':i,'pid':child.pid,'exitCode':code,'stderrBytes':(base/f'worker-{i}.stderr').stat().st_size});f.save(base/'workers.json',outcomes)
  f.require(all(v['exitCode']==0 and v['stderrBytes']==0 for v in outcomes),'all midpoint workers clean')
  for i,child,directory in children:
   checks=json.loads((directory/'worker-checks.json').read_text());f.require(checks==json.loads((base/f'worker-{i}.stdout').read_text()),'worker checks/stdout identity')
   f.require(checks['sourceFiles']==data['pins'] and checks['inputPackets']==data['operands'],'worker/coordinator source and input identity')
   f.require([x['index'] for x in checks['records']]==list(range(2*i+1,COUNT,2*WORKERS)),'actual worker midpoint ownership')
   for record in checks['records']:
    for name,v in record['artifacts'].items():f.require(f.digest(Path(record['directory'])/name)==v['sha256'],'every completed midpoint artifact hash')
    x=f.unpickle(Path(record['directory'])/'contour-point.pickle');f.require(x['index']==record['index'] and x['frequency']==frequencies()[record['index']],'actual new midpoint address')
    inventory.append(dict(record,originalIndex=record['index'],reused=False))
  inventory.sort(key=lambda x:x['index']);f.require([x['index'] for x in inventory]==list(range(COUNT)),'all32 actual/reused points');f.save(base/'point-inventory.json',inventory)
  loopdir=base/'end-loops';loopdir.mkdir();_,continued,_=c.adapters();closures=loop(loopdir,data,inventory,continued)
  results,comparison,metadata=summarize(base,data,inventory)
  old=f.unpickle(data['oldRun']/'contour-16.pickle')
  for key in ('inverseMoments','responseMoments'):
   for x,y in zip(results[16][key],old[key]):f.require(np.array_equal(x,y),'exact retained16 moment comparison')
  f.require(np.array_equal(results[16]['winding']['increments'],old['winding']['increments']),'exact retained16 winding increments')
  original_hashes(data)
  f.require(all(f.digest(f.ROOT/n)==h for n,h in data['pins'].items()) and all(f.digest(Path(n))==h for n,h in data['operands'].items()),'all source/input post hashes')
  result={'status':'COMPLETED_NUMERICAL_CONTOUR_MIDPOINT_REFINEMENT','sourceFiles':data['pins'],'inputPackets':data['operands'],'workers':outcomes,'joins':joins,'points':32,'reusedPoints':16,'newPoints':16,'center':str(c.CENTER),'radius':c.RADIUS,'winding':{str(n):{k:v for k,v in results[n]['winding'].items() if k!='increments'} for n in results},'inverseMomentNorms':{str(n):[end.norm(v) for v in results[n]['inverseMoments']] for n in results},'responseMomentNorms':{str(n):[end.norm(v) for v in results[n]['responseMoments']] for n in results},'inverseMomentChangeNorms':[end.norm(v) for v in comparison['inverseMomentDifferences']],'responseMomentChangeNorms':[end.norm(v) for v in comparison['responseMomentDifferences']],'maximumPointCondition':max(v['condition'] for v in metadata),'minimumPointSingularValue':min(v['smallestSingularValue'] for v in metadata),'maximumEndClosureDifference':max(end.norm(v['closureDifference']) for v in closures),'newMomentumNodes':sum(x['newMomentumNodes'] for x in inventory if not x['reused']),'unchangedAcceptedArtifacts':len(data['oldChecks']['artifacts']),'artifacts':{str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*.pickle') if 'source' not in p.relative_to(base).parts},'wallSeconds':time.monotonic()-start,'peakCoordinatorRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':data['manifest']['scope'],'limitations':['No certified contour-interior exceptional-locus or spectral completeness claim.','Inverse moments are not assumed projectors or simple residues.','No physical bound-pole set or global empty spectrum is supplied.']}
  f.save(base/'checks.json',result);signal.alarm(0);print(json.dumps(result,indent=2))
 finally:
  for _,child,_ in children:
   if child.poll() is None:child.terminate()
  for _,child,_ in children:
   if child.poll() is None:
    try:child.wait(timeout=5)
    except subprocess.TimeoutExpired:child.kill();child.wait()
  for handle in logs:handle.close()

if __name__=='__main__':main()
