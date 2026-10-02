#!/usr/bin/env python3
"""Standard-library continuation/storage/source tests; no scientific imports."""
import ast
import json
from pathlib import Path
import runpy
import symtable
import tempfile
from types import SimpleNamespace

M=Path(__file__).resolve().parent
WORKER=M/'S11c_d_defect_raw_increment_continue.py'
ORIGINAL=M/'S11c_d_defect_raw_increment.py'

def run():
    tests=[]
    def check(name,condition):
        if not condition:raise AssertionError(name)
        tests.append(name)
    def raises(fn,exception=ValueError):
        try:fn()
        except exception:return True
        return False
    source=WORKER.read_text();original=ORIGINAL.read_text();ns=runpy.run_path(str(WORKER),run_name='stdlib_tests_only')
    helpers,parts=ns['source_parts'](original)
    compile(helpers,'helpers','exec');compile(parts,'controls','exec');compile(source,str(WORKER),'exec')
    helper_names={n.name for n in helpers.body}
    check('only selected helper definitions',helper_names==set(ns['HELPERS']) and 'scientific_work' not in helper_names and 'main' not in helper_names)
    check('only q and unfinished controls executable',{n.name for n in parts.body}=={'q','resumed_controls'})
    outer=next(n for n in ast.parse(original).body if isinstance(n,ast.FunctionDef) and n.name=='scientific_work')
    first=next(i for i,n in enumerate(outer.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='control_k' for t in n.targets))
    old_tail=ast.Module(body=outer.body[first:],type_ignores=[]);new_tail=ast.Module(body=parts.body[1].body,type_ignores=[])
    assignments=lambda tree:[ast.dump(n,include_attributes=False) for n in ast.walk(tree) if isinstance(n,(ast.Assign,ast.AugAssign,ast.AnnAssign))]
    check('all unfinished mathematical assignment ASTs identical',assignments(old_tail)==assignments(new_tail))
    old=ast.dump(old_tail,include_attributes=False);new=ast.dump(new_tail,include_attributes=False)
    old_node="Name(id='control_depths', ctx=Load())"
    repair="Call(func=Name(id='string_keys', ctx=Load()), args=[Name(id='control_depths', ctx=Load())], keywords=[])"
    check('exactly one evidence-only AST change',new.count(repair)==1 and new.replace(repair,old_node)==old)
    check('q definition unchanged',ast.dump(parts.body[0],include_attributes=False)==ast.dump(next(n for n in outer.body if isinstance(n,ast.FunctionDef) and n.name=='q'),include_attributes=False))
    check('duplicate or absent repair refused',raises(lambda:ns['source_parts'](original.replace("'depths':control_depths","'other':control_depths"))))
    tree=ast.parse(source);main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main');main_text=ast.get_source_segment(source,main)
    check('gate and containment before science',main_text.index('verify_gate(')<main_text.index("helpers['containment']()")<main_text.index('import sympy'))
    check('no top-level science import',not any('sympy' in ast.unparse(n) or 'numpy' in ast.unparse(n) for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom))))
    check('no computation timers',not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr in ['alarm','setitimer'] for n in ast.walk(tree)))
    def find_context():
        fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run_controls')
        return next(n.value for n in fn.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='context' for t in n.targets))
    context={ast.literal_eval(k) for k in find_context().keys}
    symbols=symtable.symtable(ast.unparse(parts),'tail','exec').get_children()
    import builtins
    globals_needed={s.get_name() for block in symbols for s in block.get_symbols() if s.is_global() and s.is_referenced()}
    check('unfinished fragment namespace complete',not (globals_needed-context-{'q'}-set(vars(builtins))))
    class Inert:
        def __init__(self,name):self.name=name
        def __str__(self):return self.name
    depths={Inert(n):object() for n in ['q_i','q_h','q_s','q_o']};fixed=ns['string_keys'](depths)
    check('serialization key conversion retains all values',len(fixed)==4 and all(fixed[str(k)] is v for k,v in depths.items()))
    check('duplicate labels refused',raises(lambda:ns['string_keys']({Inert('q'):1,Inert('q'):2})))
    class Basic:pass
    class MatrixBase:pass
    env={'ast':ast,'hashlib':__import__('hashlib'),'json':json,'os':__import__('os'),'Path':Path,'resource':__import__('resource'),'THREADS':ns['THREADS'],'sp':SimpleNamespace(Basic=Basic,MatrixBase=MatrixBase)}
    exec(compile(helpers,'helpers','exec'),env)
    with tempfile.TemporaryDirectory() as directory:
        p=Path(directory);J=env['Journal'](p)
        check('original encoder failure reproduced',raises(lambda:J.encode({'depths':{Inert('q'):1}})))
        check('fixed encoder JSON roundtrip',json.loads(json.dumps(J.encode({'depths':ns['string_keys']({Inert('q'):1})})))=={'depths':{'q':1}})
        J.completed.append({'name':'prior','status':'RESTORED_PRIOR_COMPLETE_RETURN','functionCalled':False})
        def fail():J.emit('partial',{'unchanged':[1,2]});raise ValueError('intentional')
        check('new failure propagates',raises(lambda:J.stage('unfinished',{'a':1},fail)))
        check('partial evidence survives failure',(p/'unfinished-input.json').exists() and (p/'partial.json').exists() and J.active=='unfinished')
        check('prior complete remains distinct',len(J.completed)==1 and J.completed[0]['status']=='RESTORED_PRIOR_COMPLETE_RETURN')
        check('prior files immutable',raises(lambda:J.emit('partial',{'changed':2}),FileExistsError))
    with tempfile.TemporaryDirectory() as directory:
        p=Path(directory);prior=p/'old';prior.mkdir();(prior/'a.json').write_text('{"srepr":"inert, never decoded"}\n');(prior/'nested').mkdir();(prior/'nested/b').write_bytes(b'opaque\x00bytes')
        inventory={'root':str(prior),'files':{str(f.relative_to(prior)):{'bytes':f.stat().st_size,'sha256':ns['sha'](f)} for f in prior.rglob('*') if f.is_file()}}
        ip=p/'inventory.json';ip.write_text(json.dumps(inventory));out=p/'new';out.mkdir()
        target,files=ns['copy_prior']({'priorComplete':str(prior),'priorInventory':str(ip)},out)
        check('byte-identical entire prior copy',all((prior/n).read_bytes()==(target/n).read_bytes() for n in files))
        (prior/'a.json').write_text('changed');out2=p/'new2';out2.mkdir()
        check('changed prior input rejected',raises(lambda:ns['copy_prior']({'priorComplete':str(prior),'priorInventory':str(ip)},out2)))
    check('no completed construction in tail',not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id in ['scientific_work','lower_boundary','factorization'] for n in ast.walk(parts)))
    # Runtime argument equality is enforced on actual restored data, not text rendering.
    check('actual prior argument guard present',"require(actual==expected_inputs[name]" in source and 'RESTORED_PRIOR_COMPLETE_RETURN' in source)
    return {'status':'PASS_STDLIB_TOOLING_ONLY','count':len(tests),'tests':tests,'mathematicalAssignmentCount':len(assignments(old_tail)),
            'scienceImported':False,'scientificPayloadRestored':False,'physicsValidated':False}

if __name__=='__main__':print(json.dumps(run(),indent=2))
