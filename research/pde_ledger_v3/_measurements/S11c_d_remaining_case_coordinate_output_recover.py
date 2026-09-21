#!/usr/bin/env python3
"""Continue saved output preparation using actual chart-owned unit declarations."""
import ast, copy, json, shutil
from pathlib import Path
import S11c_d_remaining_case_coordinate_output as h

f,m=h.f,h.m
ORIGIN=f.STORE/'s11c-remaining-case-coordinate-20260921/response/output/focused/complete'
PLAN=f.M/'S11c_d_remaining_case_coordinate_output_recovery_plan.md'
REPAIR=f.M/'S11c_d_remaining_case_coordinate_output_unit_reader_repair.json'
ORIGINAL_MAIN=h.main
BASELINE=h.h.BASELINE


def unit_join(chart, frame, coordinates, known, raw, address, expected):
    f.require(address==expected, 'actual chart unit producer case/address')
    f.require(m.same(chart,frame) and m.same(chart['coordinates'],coordinates), 'full saved chart/frame/coordinate identity')
    f.require(len(coordinates)==3 and len(set(coordinates))==3, 'three distinct actual material coordinates')
    joined={}
    for atom in coordinates:
        f.require(atom in known and tuple(known[atom])==(1,0,0), 'actual chart-owned coordinate declaration')
        if atom in raw:f.require(m.same(raw[atom],known[atom]), 'existing response declaration agrees')
        joined[atom]=known[atom]
    return joined


def coordinate_units(base,label,value,folder):
    manifest=json.loads((base/'inputs.json').read_text())
    path=base/'numerical/boundary-inputs/bindings/sources/coordinate-cases'/label/'coordinate-source.pickle'
    frame_path=base/'numerical/boundary-inputs/bindings/contexts'/label/'frame.pickle'
    for p in (path,frame_path):f.require(f.digest(p)==manifest['copiedInputs'][str(p.relative_to(base))], 'unchanged actual chart/frame input')
    producer=f.unpickle(path);frame=f.unpickle(frame_path)
    chart=value['chartState']['values']['g'];coords=chart['coordinates'];known=producer['dimensionState']['known'];raw=value['dimensionState']['known']
    address=str(path.relative_to(base));expected='numerical/boundary-inputs/bindings/sources/coordinate-cases/'+label+'/coordinate-source.pickle'
    native=f.M/'S11c_d_coordinate_source.py';tree=ast.parse(native.read_text());fn=h.h.native.function(tree,'chart')
    loop=next(n for n in fn.body if isinstance(n,ast.For) and ast.unparse(n.target)=='v' and ast.unparse(n.iter)=='x')
    f.require(ast.dump(loop)==ast.dump(ast.parse('for v in x: dimensions.known[v]=(1,0,0)').body[0]), 'unchanged native chart unit declaration body')
    f.require(f.digest(native)==manifest['sourceFiles'][str(native.relative_to(f.ROOT))], 'accepted native chart source hash')
    # Preserve the full producer/frame and raw response history before new guards.
    pairs={'case':label,'address':address,'producerSha256':f.digest(path),'frameSha256':f.digest(frame_path),
        'chart':chart,'producerChart':producer['chart'],'frameGeometry':frame['geometry'],'coordinates':coords,
        'producerKnown':{a:known[a] for a in coords if a in known},
        'rawResponseKnown':{a:raw[a] for a in coords if a in raw},'rawAbsent':tuple(a for a in coords if a not in raw),
        'nativeSourceSha256':f.digest(native),'nativeUnitLoop':ast.dump(loop),'responsePacketChanged':False}
    f.atomic_pickle(folder/'coordinate-unit-source-pairs.pickle',pairs)
    f.require(m.same(producer['chart'],chart), 'complete actual coordinate producer geometry')
    joined=unit_join(chart,frame['geometry'],coords,known,raw,address,expected)
    mutations=[]
    for atom in coords:
        changed=dict(known);changed[atom]=(2,0,0)
        mutations.append(('producer-unit',coords,changed,raw,address))
        changed=dict(raw);changed[atom]=(0,1,0)
        mutations.append(('response-unit',coords,known,changed,address))
        changed=dict(known);del changed[atom]
        mutations.append(('missing-producer-unit',coords,changed,raw,address))
    mutations.extend([('coordinate-order',tuple(reversed(coords)),known,raw,address),
        ('duplicate-coordinate',(coords[0],coords[0],coords[2]),known,raw,address),
        ('wrong-producer-address',coords,known,raw,address+'.wrong')])
    records=[]
    for name,actual,k,r,a in mutations:
        rejected=False
        try:unit_join(chart,frame['geometry'],actual,k,r,a,expected)
        except (ValueError,AssertionError,RuntimeError):rejected=True
        records.append({'kind':name,'coordinates':actual,'known':{x:k[x] for x in coords if x in k},
            'raw':{x:r[x] for x in coords if x in r},'address':a,'rejected':rejected})
    f.atomic_pickle(folder/'coordinate-unit-mutation-controls.pickle',records)
    f.require(all(x['rejected'] for x in records), 'actual changed coordinate unit/address controls reject')
    f.save(folder/'coordinate-unit-checks.json',{'coordinates':len(joined),'rawAbsent':len(pairs['rawAbsent']),
        'respondingMutations':len(records),'nativeUnitDeclarationJoined':True,'rawResponseStateUnchanged':True,
        'producerSha256':pairs['producerSha256'],'frameSha256':pairs['frameSha256']})


def saved_prefix(base,manifest,label):
    f.require(label==BASELINE, 'only the actually completed baseline bundle prefix is reused')
    folder=base/'bundles'/label;packet=folder/'coordinate-output.pickle';origin=ORIGIN/'bundles'/label/'coordinate-output.pickle'
    f.require(f.digest(packet)==f.digest(origin)==manifest['copiedInputs'][str(packet.relative_to(base))], 'complete original bundle byte identity')
    value=f.unpickle(packet)
    f.require(value['case']==label and value['sourceFiles']==json.loads((ORIGIN/'inputs.json').read_text())['sourceFiles'], 'original bundle physical address and source pins')
    # The immutable traceback reaches the next statement only after these
    # original loops finished. Restore their counter; do not rebuild the bundle.
    count=sum(1 for _ in h.boundary_arrays(value))
    f.save(folder/'completed-focus-prefix-reuse.json',{'originalPacket':str(origin),'packetSha256':f.digest(packet),
        'bytes':packet.stat().st_size,'boundaryArrayViews':count,'completedBoundaryAndFinitePrefixReused':True,
        'tracebackSha256':f.digest(ORIGIN.parent/'coordinate_output_inputs.stderr'),
        'originalHelperSha256':f.digest(Path(h.__file__)),'newBundleConstruction':False,'scientificWork':0})
    return folder,value,packet,count


def focus_adapter():
    original=h.h.native.function(ast.parse(Path(h.__file__).read_text()),'focused');node=copy.deepcopy(original)
    loop=next(n for n in node.body if isinstance(n,ast.For) and ast.unparse(n.iter)=='labels')
    unit_index=next(i for i,n in enumerate(loop.body) if isinstance(n,ast.For) and ast.unparse(n.iter)=="value['chartState']['values']['g']['coordinates']")
    prefix=copy.deepcopy(loop.body[:unit_index]);unit_loop=copy.deepcopy(loop.body[unit_index])
    replacement=ast.parse('if label==BASELINE:\n folder,value,packet,count=saved_prefix(base,manifest,label)\nelse:\n pass').body[0]
    replacement.orelse=copy.deepcopy(prefix)
    loop.body=[replacement,ast.parse('coordinate_units(base,label,value,folder)').body[0],*loop.body[unit_index+1:]]
    back=copy.deepcopy(node);back_loop=next(n for n in back.body if isinstance(n,ast.For) and ast.unparse(n.iter)=='labels')
    back_loop.body=[*back_loop.body[0].orelse,unit_loop,*back_loop.body[2:]]
    f.require(ast.dump(back)==ast.dump(original),'whole original focus reverse AST, saved prefix and unit reader only')
    namespace=dict(vars(h),BASELINE=BASELINE,saved_prefix=saved_prefix,coordinate_units=coordinate_units)
    return h.h.native.compile_function(node,namespace),{'wholeFocusReverseAST':True,'savedBaselinePrefixOnly':True,
        'coordinateUnitReaderOnly':True,'allUnfinishedEmitterSampleAndHashGuardsUnchanged':True,
        'originalHelperSha256':f.digest(Path(h.__file__)),'originalMainUnchanged':h.main is ORIGINAL_MAIN}


def load(base):
    old=json.loads((ORIGIN/'inputs.json').read_text());repair=json.loads(REPAIR.read_text())
    f.require(repair['originalHelperSha256']==f.digest(Path(h.__file__)) and repair['originalInputsSha256']==f.digest(ORIGIN/'inputs.json'), 'explicit original focus/helper join')
    f.require(repair['recoveryHelperSha256']==f.digest(Path(__file__)) and repair['planSha256']==f.digest(PLAN), 'reviewed recovery sources')
    inv=json.loads((ORIGIN.parent/'coordinate_output_inputs.invocation.json').read_text());guard=json.loads((ORIGIN.parent/'resource-guard/outcome.json').read_text())
    f.require(inv['exitCode']==guard['exitCode']==guard['childOutcome']['exitCode']==1 and guard['limitsVerified'] and guard['childOutcome']['guardReason'] is None and not (ORIGIN/'checks.json').exists(), 'actual stopped focus, not acceptance')
    f.require(f.digest(ORIGIN.parent/'coordinate_output_inputs.stderr')==repair['originalStderrSha256'], 'exact original failed statement')
    for n,sha in old['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(ORIGIN/'source'/n)==sha, 'original/current/frozen output sources')
    for n,sha in old['inputPackets'].items():f.require(f.digest(Path(n))==sha, 'original accepted input identity')
    for n,sha in old['copiedInputs'].items():f.require(f.digest(ORIGIN/n)==sha, 'every completed original copy')
    bundles=list((ORIGIN/'bundles').glob('*/coordinate-output.pickle'))
    f.require(bundles==[ORIGIN/'bundles'/BASELINE/'coordinate-output.pickle'], 'exact actually completed bundle census')
    manifest={k:copy.deepcopy(v) for k,v in old.items() if k not in ('runDirectory','copiedInputs')}
    manifest.update(runDirectory=str(base),copiedInputs={})
    for p in (Path(__file__),PLAN,REPAIR):manifest['sourceFiles'][str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    inventory={}
    for p in sorted(ORIGIN.rglob('*')):
        if not p.is_file():continue
        n=str(p.relative_to(ORIGIN));target='original-output-inputs.json' if n=='inputs.json' else 'original-output-source/'+n[len('source/'):] if n.startswith('source/') else n
        sha=f.digest(p);m.retain(p,base/target,manifest,sha);inventory[n]={'sha256':sha,'bytes':p.stat().st_size,'copy':target}
    for name in ('active.json','coordinate_output_inputs.invocation.json','coordinate_output_inputs.stderr','coordinate_output_inputs.stdout','guard.stdout','guard.stderr','launch.json','static-review.json'):
        m.retain(ORIGIN.parent/name,base/'original-output-logs'/name,manifest)
    for p in sorted((ORIGIN.parent/'resource-guard').glob('*')):
        if p.is_file():m.retain(p,base/'original-output-logs/resource-guard'/p.name,manifest)
    _,join=focus_adapter();f.save(base/'recovery-focus-join.json',join)
    reuse={'originalDirectory':str(ORIGIN),'originalOutcome':inv,'files':inventory,'completedBundles':1,'completedNewPhysicalTranscripts':0,
        'newScientificOperations':0,'originalFailedGuardPreserved':True,'originalMainUnchanged':True}
    f.save(base/'completed-output-input-reuse.json',reuse);manifest['completedOutputInputReuse']=reuse
    for n,sha in manifest['sourceFiles'].items():
        p=base/'source'/n;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,p);f.require(f.digest(p)==sha,'frozen recovery source')
    f.save(base/'inputs.json',manifest)
    cp=json.loads(h.CP.read_text());f.require(cp['status']=='ACCEPTED_CASE_MATERIAL_NUMERICAL_RESPONSES' and cp['checksSha256']==old['numericalChecksSha256'], 'unchanged numerical acceptance')
    return manifest,tuple(cp['checks']['cases']),json.loads((base/'part-inventory.json').read_text())


if __name__=='__main__':
    adapter,join=focus_adapter();h.focused=adapter;h.load=load
    f.require(h.main is ORIGINAL_MAIN,'complete original main function unchanged')
    h.main()
