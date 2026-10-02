#!/usr/bin/env python3
"""Stdlib-only persistence/source-routing tests; no science imported or restored."""
import ast
from fractions import Fraction
import json
from pathlib import Path
import runpy
import tempfile
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
WORKER=HERE/'S11c_d_defect_raw_increment.py'


def run():
    source=WORKER.read_text();tree=ast.parse(source)
    namespace=runpy.run_path(str(WORKER),run_name='storage_test_only')
    journal_type=namespace['Journal']
    # Journal test payloads are ordinary Python containers, never symbolic data.
    class NotScientific: pass
    journal_type.encode.__globals__['sp']=SimpleNamespace(Basic=NotScientific,MatrixBase=NotScientific)
    tests=[]
    def check(name,condition):
        if not condition:raise AssertionError(name)
        tests.append(name)
    top_imports=[n for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom))]
    check('scientific imports not at top level',not any('sympy' in ast.unparse(n) or 'numpy' in ast.unparse(n) for n in top_imports))
    main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
    text=ast.get_source_segment(source,main)
    check('containment before scientific import',text.index('containment()')<text.index('import sympy'))
    check('gate before scientific import',text.index('verify_gate(')<text.index('import sympy'))
    check('no scientific timer calls',not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr in ('alarm','setitimer') for n in ast.walk(tree)))
    with tempfile.TemporaryDirectory() as directory:
        out=Path(directory);J=journal_type(out)
        value=J.stage('complete',{'value':[1,2]},lambda:{'return':3})
        check('complete return recorded',value=={'return':3} and json.loads((out/'operation-index.json').read_text())[0]['name']=='complete')
        check('input and return receipts hash bytes',all(namespace['sha'](out/r['path'])==r['sha256'] for r in J.completed[0].values() if isinstance(r,dict)))
        def fail():
            J.emit('partial-result',{'value':4});raise ValueError('intentional')
        try:J.stage('unfinished',{'value':5},fail)
        except ValueError:pass
        else:raise AssertionError('failure swallowed')
        check('failed input and partial preserved',(out/'unfinished-input.json').exists() and (out/'partial-result.json').exists())
        check('failed stage not falsely complete',len(json.loads((out/'operation-index.json').read_text()))==1 and J.active=='unfinished')
        check('partial artifact indexed','partial-result.json' in json.loads((out/'artifact-index.json').read_text()))
        try:J.emit('partial-result',{'value':99})
        except FileExistsError:pass
        else:raise AssertionError('overwrite allowed')
        check('immutable result not overwritten',json.loads((out/'partial-result.json').read_text())=={'value':4})
    geom=(HERE.parent/'scripts/S11c_a_interface_geometry_sympy_audit.py').read_text()
    nodes=[n for n in ast.walk(ast.parse(geom)) if isinstance(n,ast.Assign) and any(isinstance(a,ast.Name) and a.id=='normal_exact' for a in n.targets) and 'grad_h' in ast.get_source_segment(geom,n)]
    check('unique native normal source',len(nodes)==1)
    snippet=ast.get_source_segment(geom,nodes[0])
    # Synthetic rational operands only; catches Python comprehension namespace lookup.
    local={'__builtins__':{},'tuple':tuple,'face':-1,'grad_h':(Fraction(-2),0,0),'denominator':Fraction(3)}
    exec(snippet,local)
    check('native snippet same-namespace execution',local['normal_exact']==(Fraction(-2,3),0,0,Fraction(-1,3)))
    for name,target in [('kernel_bridge','three_inverse')]:
        c2=(HERE.parent/'scripts/S11c_c2_selfenergy_fold_sympy_audit.py').read_text()
        fragment=namespace['assignment_source'](c2,name,target)
        compile(fragment,'native_assignment','exec');check('unique native assignment '+target,target in fragment)
    check('missing native function refused',_raises(lambda:namespace['function_source']('def present(): pass','absent')))
    return {'status':'PASS_STDLIB_TOOLING_ONLY','tests':tests,'count':len(tests),'scientificImport':False,'scientificPayloadRestoration':False,'physicsValidated':False}


def _raises(fn):
    try:fn()
    except ValueError:return True
    return False


if __name__=='__main__':print(json.dumps(run(),indent=2))
