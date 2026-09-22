#!/usr/bin/env python3
"""Route saved threshold inputs directly to their unchanged resolved files."""
import ast
import copy
import hashlib
import inspect
import json
from pathlib import Path
import sys
import textwrap

sys.dont_write_bytecode=True
import S11c_d_remaining_case_frequency_end_complex_threshold_finish as finish
original=finish.original
f=finish.f
PLAN=f.M/'S11c_d_remaining_case_frequency_end_complex_threshold_link_finish_plan.md'
REPAIR=f.M/'S11c_d_remaining_case_frequency_end_complex_threshold_link_repair.json'
ROUTES={}
STATE={}


def logical_route(path):
    path=Path(path);resolved=path.resolve();current=path;seen=set();chain=[]
    while current.is_symlink():
        f.require(str(current) not in seen,'acyclic original reference route')
        seen.add(str(current));target=current.readlink()
        chain.append({'path':str(current),'target':str(target)})
        current=target if target.is_absolute() else current.parent/target
    f.require(current.resolve()==resolved,'actual complete resolved reference route')
    return {'path':str(path),'rawLink':str(path.readlink()) if path.is_symlink() else None,
        'resolved':str(resolved),'bytes':resolved.stat().st_size,'symlinkChain':chain}


def canonical_reference(base,manifest,path,name,expected):
    route=logical_route(path);route.update(sha256=expected,consumer=name)
    f.require(name not in ROUTES,'unique logical reference consumer')
    ROUTES[name]=route
    # Preserve every original logical hop above; only the new link's target
    # changes to the exact canonical file. The native hash/size writer remains.
    original.h.source.reference(base,manifest,Path(route['resolved']),name,expected)


def retain_boundary_evidence(base,manifest):
    repair=STATE['repair'];failed=Path(repair['failedDirectory'])
    inventory=json.loads(Path(repair['fileInventory']['path']).read_text())
    for name,item in inventory.items():
        route=logical_route(failed/'complete'/name)
        for key in ('path','rawLink','resolved','bytes','symlinkChain'):
            f.require(route[key]==item[key],('preserved failed logical reference',name,key))
        canonical_reference(base,manifest,failed/'complete'/name,'failed-link-workspace/'+name,item['sha256'])
    for name,item in repair['failedLogs'].items():
        route=logical_route(failed/name)
        for key in ('path','rawLink','resolved','bytes','symlinkChain'):
            f.require(route[key]==item[key],('preserved failed log route',name,key))
        canonical_reference(base,manifest,failed/name,'failed-link-logs/'+name,item['sha256'])
    for key,name in (('fileInventory','failed-link-file-inventory.json'),('outcomeInspection','failed-link-outcome.json')):
        item=repair[key];canonical_reference(base,manifest,Path(item['path']),name,item['sha256'])
    for path in (Path(__file__).resolve(),PLAN,REPAIR):
        name=str(path.relative_to(f.ROOT));f.require(name not in manifest['sourceFiles'],'fresh link-boundary source pin')
        value=f.digest(path);manifest['sourceFiles'][name]=value;target=base/'source'/name
        target.parent.mkdir(parents=True,exist_ok=True)
        f.require(not target.exists() and not target.is_symlink(),'fresh link-boundary frozen source')
        target.write_bytes(path.read_bytes())
    manifest['referenceLinkBoundary']={'failedDirectory':str(failed),'failedOutcome':repair['failedOutcome'],
        'wholeFinishHelperSha256':repair['originalFinishHelperSha256'],'wholeLoadJoin':STATE['loadJoin'],
        'wholePacketReaderJoin':STATE['packetJoin'],'failedPartialReferencesRetained':len(inventory),
        'logicalRoutes':len(ROUTES),'newLinksUseExactCanonicalFiles':True,'logicalSourceRoutesPreserved':True,
        'originalManifestAndOperandLoadingDidNotComplete':True,'completedScienceRepeated':False}
    f.save(base/'reference-link-boundary-joins.json',{'routes':ROUTES,'continuation':manifest['referenceLinkBoundary']})


def load_adapter():
    before=ast.parse(textwrap.dedent(inspect.getsource(finish.load))).body[0]
    index=next(i for i,n in enumerate(before.body) if isinstance(n,ast.Assign)
        and any(isinstance(t,ast.Name) and t.id=='ref' for t in n.targets))
    saved=copy.deepcopy(before.body[index]);changed=copy.deepcopy(before)
    f.require(ast.unparse(saved.value.body.func)=='original.h.source.reference','exact native reference boundary')
    changed.body[index].value.body.func=ast.Name(id='canonical_reference',ctx=ast.Load())
    position=next(i for i,n in enumerate(changed.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call)
        and ast.unparse(n.value.func)=='STATE.update')
    changed.body.insert(position,ast.parse('retain_boundary_evidence(base,manifest)').body[0])
    reverse=copy.deepcopy(changed);reverse.body.pop(position);reverse.body[index]=saved
    f.require(ast.dump(reverse)==ast.dump(before),'whole finish load: canonical reference and extra preserved evidence only')
    env=dict(vars(finish),canonical_reference=canonical_reference,retain_boundary_evidence=retain_boundary_evidence)
    module=ast.fix_missing_locations(ast.Module(body=[changed],type_ignores=[]))
    exec(compile(module,'<same saved loader with direct reference boundary>','exec'),env)
    return env['load'],{'entireFinishLoadReverseAST':True,'originalLoadAST':hashlib.sha256(ast.dump(before).encode()).hexdigest(),
        'originalLoaderNotCalled':True,'originalDataChecksAndPacketOutputRouteUnchanged':True,'adapterSource':ast.unparse(module)}


def packet_adapter():
    before=ast.parse(textwrap.dedent(inspect.getsource(finish.packet))).body[0]
    changed=copy.deepcopy(before);hits=0
    class Route(ast.NodeTransformer):
        def visit_Compare(self,node):
            nonlocal hits
            if ast.unparse(node.left)=='path.readlink()':
                f.require(len(node.comparators)==1 and ast.unparse(node.comparators[0])=="PREVIOUS / 'complete' / name",'exact saved packet address boundary')
                node.comparators[0]=ast.Call(func=ast.Attribute(value=node.comparators[0],attr='resolve',ctx=ast.Load()),args=[],keywords=[])
                hits+=1
            return self.generic_visit(node)
    changed=Route().visit(changed);f.require(hits==1,'one canonical packet link expectation')
    class Reverse(ast.NodeTransformer):
        def visit_Compare(self,node):
            if ast.unparse(node.left)=='path.readlink()':node.comparators[0]=node.comparators[0].func.value
            return self.generic_visit(node)
    f.require(ast.dump(Reverse().visit(copy.deepcopy(changed)))==ast.dump(before),'whole packet reader: canonical address only')
    module=ast.fix_missing_locations(ast.Module(body=[changed],type_ignores=[]));env=dict(vars(finish))
    exec(compile(module,'<same typed packet reader with canonical link address>','exec'),env)
    return env['packet'],{'entireFinishPacketReaderReverseAST':True,'originalPacketAST':hashlib.sha256(ast.dump(before).encode()).hexdigest(),
        'oneResolvedPathExpectation':True,'allSizeHashValueAndSourceGuardsUnchanged':True,'adapterSource':ast.unparse(module)}


def load(base):
    repair=json.loads(REPAIR.read_text());STATE['repair']=repair
    f.require(f.digest(Path(finish.__file__))==repair['originalFinishHelperSha256'],'immutable failed packet-finish helper')
    f.require(f.digest(finish.PLAN)==repair['originalFinishPlanSha256'] and f.digest(finish.REPAIR)==repair['originalFinishRepairSha256'],'immutable failed plan/repair')
    for key in ('fileInventory','outcomeInspection'):
        item=repair[key];f.require(f.digest(Path(item['path']))==item['sha256'],'failed reference boundary inspection identity')
    tree=ast.parse(Path(finish.__file__).read_text())
    for name,value in repair['wholeFinishBodies'].items():
        node=next(n for n in tree.body if getattr(n,'name',None)==name)
        f.require(hashlib.sha256(ast.dump(node).encode()).hexdigest()==value,('whole failed finish body',name))
    reader,reader_join=packet_adapter();loader,loader_join=load_adapter()
    STATE.update(packetJoin=reader_join,loadJoin=loader_join)
    finish.packet=reader
    return loader(base)


def joined_sources(base,manifest):
    functions,joins,dispatch=finish.joined_sources(base,manifest)
    joins['referenceLinkBoundary']={'wholeFinishHelperSha256':STATE['repair']['originalFinishHelperSha256'],
        'loadJoin':STATE['loadJoin'],'packetJoin':STATE['packetJoin'],
        'referenceWriter':original.source_body(original.h.source.reference),
        'canonicalReference':original.source_body(canonical_reference),'logicalRoute':original.source_body(logical_route),
        'failedOutcome':STATE['repair']['failedOutcome'],'physicalOperandOrNativeAlgorithmChanged':False}
    f.save(base/'link-boundary-native-end-complex-threshold-joins.json',joins)
    return functions,joins,dispatch


def construct(base,cp,functions,joins,dispatch):
    result,router=finish.construct(base,cp,functions,joins,dispatch)
    edges={}
    for item in ROUTES.values():
        for link in item['symlinkChain']:
            if link['path'] in edges:f.require(edges[link['path']]==link['target'],'consistent shared original reference hop')
            edges[link['path']]=link['target']
        path=Path(item['path']);canonical=Path(item['resolved'])
        f.require(str(path.resolve())==item['resolved'] and canonical.stat().st_size==item['bytes'],'logical source postidentity')
        f.require(finish.STATE['manifest']['inputPackets'][str(canonical)]==item['sha256'],'canonical source remains in mandatory final hash inventory')
    for path,target in edges.items():
        p=Path(path);f.require(p.is_symlink() and str(p.readlink())==target,'unchanged original raw reference hop')
    f.save(base/'validated-reference-link-boundary.json',{'logicalRoutes':len(ROUTES),'uniqueOriginalLinkHops':len(edges),
        'allOriginalRawLinksAndResolvedSizesUnchanged':True,'allCanonicalFileHashesRequiredByOriginalFinalGuards':True,
        'newScientificOperations':0})
    result['referenceLinkBoundary']=joins['referenceLinkBoundary']
    return result,router


def main():
    before=ast.parse(textwrap.dedent(inspect.getsource(original.main))).body[0]
    repair=json.loads(finish.REPAIR.read_text())
    f.require(hashlib.sha256(ast.dump(before).encode()).hexdigest()==repair['wholeOriginalBodies']['main'],'whole unchanged original main/final hashes')
    env=dict(vars(original),load=load,joined_sources=joined_sources,construct=construct,prohibit=finish.prohibit)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[before],type_ignores=[])),'<whole original main with explicit reference boundary slots>','exec'),env)
    env['main']()

if __name__=='__main__':main()
