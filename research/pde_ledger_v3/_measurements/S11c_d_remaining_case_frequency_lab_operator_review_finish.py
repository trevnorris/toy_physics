#!/usr/bin/env python3
"""Finish the saved response review using the actual native SVDResult class."""
import ast
import copy
import importlib
import json
from pathlib import Path
import S11c_d_remaining_case_frequency_lab_operator_review as original

NAME='S11c_d_remaining_case_frequency_lab_operator_review_finish'
ORIGINAL_SHA='6f37c1fe553132e17784857ec558c4ba863f5814fd8b88ed650da37ba0a1667b'
FAILED=original.F/'lab-operator-review-recovery-01/failed-review-file-inventory.json'
FAILED_SHA='822d4e6096a1358d40b7b316f298660ba488d93db82bfbb9be93b86e0f96d26a'
SVD_SOURCE_SHA='90354ad460716d48e5460c570a811f911266f32acdaf522eee1327dab28af111'


def svd_matches(value,local):
    cls=importlib.import_module('numpy.linalg.linalg').SVDResult
    return (type(value) is cls and cls._fields==('U','S','Vh') and
            original.same(value.U,local['u']) and original.same(value.S,local['s']) and original.same(value.Vh,local['vh']))


def restore(reader,journal,checks,arts,packet,read,info,read_packets):
    require=original.require
    reader.retain(original.__file__,ORIGINAL_SHA)
    reader.retain(original.M/'S11c_d_remaining_case_frequency_lab_operator_review_plan.md','55458f082bde33b8ea276ccd74e8ad1d1a7b8ee4887c8cbe397a1571118cb557')
    inventory=reader.json(FAILED,FAILED_SHA);oldbase=original.F/'lab-operator-review/complete';count=0
    for path,rec in inventory.items():
        actual=reader.retain(path,rec['sha256']);require(actual['canonical']==rec['canonical'] and actual['directLink']==rec['directLink'],'preserved failed raw and canonical route')
        p=Path(path)
        if p.is_relative_to(oldbase):
            name=str(p.relative_to(oldbase));target=journal.base/name
            require(not target.exists() and not target.is_symlink(),'new saved-summary reference only');target.parent.mkdir(parents=True,exist_ok=True)
            target.symlink_to(actual['canonical']);journal.artifacts[name]={'path':str(target),'sha256':rec['sha256'],'bytes':rec['bytes']};reader.retain(target,rec['sha256']);count+=1
    require(count==191,'all actual completed review summaries retained without rerunning prefix')
    caller=info('native-and-numerical-callers.json');rec=caller['linearAlgebraSources']['np.linalg.svd'];reader.retain(rec['logical'],SVD_SOURCE_SHA)
    tree=ast.parse(Path(rec['canonical']).read_text());cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='SVDResult')
    require([n.target.id for n in cls.body if isinstance(n,ast.AnnAssign)]==['U','S','Vh'],'actual pinned native named result fields')
    returns=[ast.unparse(n) for n in ast.walk(tree) if isinstance(n,ast.Return) and isinstance(n.value,ast.Call) and ast.unparse(n.value.func)=='SVDResult']
    require('return SVDResult(wrap(u), s, wrap(vh))' in returns,'actual unmodified native SVD return')
    baseline='LAB_HELD__RHO4_CONSTANT';own='LAB_HELD__RHOBR_CONSTANT';prefix='cases/'+own;solve=prefix+'/solve'
    interior=json.loads((oldbase/'validated-interior.json').read_text());new_blocks=interior['newBlocks'];saved_blocks=interior['savedBlocks'];products=interior['products']
    sy=read(solve+'/frequency-system.pickle');result=read(solve+'/frequency-solution.pickle');summary=info(prefix+'/response-summary.json')
    local={site:read(solve+'/locals/'+site+'.pickle') for site in ('10','25')}
    # These packet reads are proved complete by191 emitted summaries plus the literal failure trace.
    # Mark their exact routes so the unchanged trailing coverage loop reads only unfinished packets.
    covered=[]
    for name,record in arts.items():
        complete=(name.endswith('/operator-input.pickle') or
                  name.startswith('coefficient-values/') and name.endswith(('/input.pickle','/value.pickle')) or
                  '/blocks/' in name and name.endswith(('/input.pickle','/value.pickle','/accumulated.pickle','/selected-direct-comparison.pickle')) or
                  name in (prefix+'/frequency-interior.pickle',solve+'/frequency-system.pickle',solve+'/frequency-solution.pickle') or
                  name.startswith(solve+'/locals/'))
        if complete:
            route=reader.retain(record['path'],record['sha256']);read_packets.add(route['canonical']);covered.append(record)
    journal.json('completed-review-prefix-reuse.json',{'source':reader.retain(original.__file__,ORIGINAL_SHA),'failedFiles':reader.retain(FAILED,FAILED_SHA),
        'completedSummaries':191,'completedPacketReads':covered,'coverageBasis':'Exact full original prefix loops and191 summaries; post-local guards precede literal SVDResult failure trace. No invented independent old pass receipt.',
        'newScientificCalls':0,'nativeSVDReturnClass':{'file':rec,'definition':ast.unparse(cls),'returns':returns},
        'typedFieldsOnly':'Actual SVDResult object retained; compare U/S/Vh individually to actual unpacked local arrays, with exact class and no tuple coercion.',
        'wholeOriginalMainReverseAST':True,'oldNumericalResponseRecomputed':False})
    return solve,sy,result,summary,local,prefix,baseline,own,new_blocks,saved_blocks,products


def adapted(source):
    fn=next(n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef) and n.name=='main');before=copy.deepcopy(fn)
    start=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='caller')
    end=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.For) and ast.unparse(n.iter)=="enumerate(('np.linalg.svd', 'la.lu_factor', 'la.lu_solve'))")
    removed=copy.deepcopy(fn.body[start:end]);replacement=ast.parse('solve,sy,result,summary,local,prefix,baseline,own,new_blocks,saved_blocks,products=restore(reader,journal,checks,arts,packet,read,info,read_packets)').body
    fn.body[start:end]=replacement
    count=0
    class Repair(ast.NodeTransformer):
        def visit_Call(self,node):
            nonlocal count
            if ast.unparse(node)=="same(val, (local['10']['u'], local['10']['s'], local['10']['vh']))":
                count+=1;return ast.copy_location(ast.parse("svd_matches(val,local['10'])",mode='eval').body,node)
            return self.generic_visit(node)
    fn=Repair().visit(fn);original.require(count==1,'one actual named-return comparison repair')
    restored=copy.deepcopy(fn);restored.body[start:start+len(replacement)]=removed
    class Reverse(ast.NodeTransformer):
        def visit_Call(self,node):
            if ast.unparse(node)=="svd_matches(val, local['10'])":return ast.copy_location(ast.parse("same(val,(local['10']['u'],local['10']['s'],local['10']['vh']))",mode='eval').body,node)
            return self.generic_visit(node)
    original.require(ast.dump(Reverse().visit(restored))==ast.dump(before),'whole original saved review prefix and result-type reverse AST')
    return ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[]))


def main():
    original.require(original.saved.digest(original.__file__)==ORIGINAL_SHA,'immutable original review')
    namespace=dict(vars(original),__file__=str(Path(__file__).resolve()),NAME=NAME,restore=restore,svd_matches=svd_matches)
    module=adapted(Path(original.__file__).read_text());exec(compile(module,'<unfinished saved response checks>','exec'),namespace);namespace['main']()


if __name__=='__main__':main()
