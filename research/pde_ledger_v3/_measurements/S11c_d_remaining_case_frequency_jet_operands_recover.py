#!/usr/bin/env python3
"""Fresh saved-reader continuation: one derivative-order JSON metadata cast."""
import ast,copy,json
from pathlib import Path
import S11c_d_remaining_case_frequency_jet_operands as original

M=original.M
REPAIR=M/'S11c_d_remaining_case_frequency_jet_operands_reader_repair.json'
PLAN=M/'S11c_d_remaining_case_frequency_jet_operands_recovery_plan.md'


def main():
    repair=json.loads(REPAIR.read_text());path=Path(original.__file__).resolve()
    original.io.require(original.io.digest(path)==repair['originalHelperSha256'],'immutable failed helper')
    tree=ast.parse(path.read_text());fn=next(v for v in tree.body if isinstance(v,ast.FunctionDef) and v.name=='main')
    changed=copy.deepcopy(fn);slots=[]
    for node in ast.walk(changed):
        if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name) and node.targets[0].id=='degree':
            slots.append(node);node.value=ast.Call(func=ast.Name(id='int',ctx=ast.Load()),args=[node.value],keywords=[])
    original.io.require(len(slots)==1,'only derivative-order JSON metadata slot')
    reverse=copy.deepcopy(changed)
    for node in ast.walk(reverse):
        if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name) and node.targets[0].id=='degree':node.value=node.value.args[0]
    original.io.require(ast.dump(reverse)==ast.dump(fn),'whole saved-reader main reverse AST')
    old_reader=original.io.Reader
    class ProvenanceReader(old_reader):
        def __init__(self):
            super().__init__()
            for item in (path,Path(__file__).resolve(),PLAN,REPAIR):self.retain(item)
            for name,record in repair['originalFiles'].items():
                actual=self.retain(name,record['sha256'])
                original.io.require(actual['bytes']==record['bytes'],'preserved failed reader metadata/log bytes')
    original.io.Reader=ProvenanceReader
    env=dict(vars(original))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[changed],type_ignores=[])),'<unchanged saved reader with one metadata int conversion>','exec'),env)
    env['main']()


if __name__=='__main__':main()
