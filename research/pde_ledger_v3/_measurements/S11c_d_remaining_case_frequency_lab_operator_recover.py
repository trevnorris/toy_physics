#!/usr/bin/env python3
"""Fresh continuation after source-dispatch inspection failed before any science."""
import ast
import copy
import inspect
from pathlib import Path
import S11c_d_remaining_case_frequency_lab_operator as original

NAME='S11c_d_remaining_case_frequency_lab_operator_recover'
ORIGINAL_SHA='06aa307c7b43c97db370cbebbb7a179042d41a374203ca16ee76bc31bf5869db'
FAILED=original.F/'lab-operator-recovery-01/failed-file-inventory.json'
FAILED_SHA='6d2251257066609e5f81e23ba4eb950c29a1925c8514f3e8e312b5012594cf66'


def preserve_failure(reader,journal):
    reader.retain(original.__file__,ORIGINAL_SHA)
    reader.retain(original.M/'S11c_d_remaining_case_frequency_lab_operator_plan.md','e8238c0a48ebfc130fe2679e0db6a079f56b8a3a587b9f2c60ba670f3f8f25b9')
    inventory=reader.json(FAILED,FAILED_SHA)
    for path,rec in inventory.items():reader.retain(path,rec['sha256'])
    journal.json('source-dispatch-recovery.json',{'originalHelper':reader.retain(original.__file__,ORIGINAL_SHA),
        'failedFileInventory':reader.retain(FAILED,FAILED_SHA),'wholeMainReverseAST':True,
        'changedExpression':'inspect.getsourcefile(fn) -> inspect.getsourcefile(inspect.unwrap(fn))',
        'numericalCallableUnchanged':True,'oldCompletedScientificPackets':0,'oldCompletedScientificCalls':0,
        'oldFailure':'NumPy _ArrayFunctionDispatcher source lookup; no scientific input restoration or operation had started.'})


def adapted(source):
    fn=next(n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef) and n.name=='main');before=copy.deepcopy(fn)
    count=0
    class Lookup(ast.NodeTransformer):
        def visit_Call(self,node):
            nonlocal count
            if ast.unparse(node)=='inspect.getsourcefile(fn)':
                count+=1;return ast.copy_location(ast.parse('inspect.getsourcefile(inspect.unwrap(fn))',mode='eval').body,node)
            return self.generic_visit(node)
    fn=Lookup().visit(fn)
    original.require(count==1,'one source-only NumPy dispatcher lookup')
    at=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='cache')
    fn.body[at:at]=ast.parse('preserve_failure(reader,journal)').body
    class Reverse(ast.NodeTransformer):
        def visit_Expr(self,node):
            if isinstance(node.value,ast.Call) and ast.unparse(node.value.func)=='preserve_failure':return None
            return self.generic_visit(node)
        def visit_Call(self,node):
            if ast.unparse(node)=='inspect.getsourcefile(inspect.unwrap(fn))':return ast.copy_location(ast.parse('inspect.getsourcefile(fn)',mode='eval').body,node)
            return self.generic_visit(node)
    original.require(ast.dump(Reverse().visit(copy.deepcopy(fn)))==ast.dump(before),'whole original main with source-only recovery')
    return ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[]))


def main():
    original.require(original.saved.digest(original.__file__)==ORIGINAL_SHA,'immutable original operator helper')
    module=adapted(Path(original.__file__).read_text())
    namespace=dict(vars(original),__file__=str(Path(__file__).resolve()),NAME=NAME,preserve_failure=preserve_failure)
    exec(compile(module,'<source-dispatch-metadata-recovery>','exec'),namespace);namespace['main']()


if __name__=='__main__':main()
