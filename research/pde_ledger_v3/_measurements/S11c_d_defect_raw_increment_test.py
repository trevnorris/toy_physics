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
    # Full-row inspection is syntax only. No symbolic constructor is evaluated.
    census=namespace['pressure_census']
    synthetic="Add(Symbol('unrelated'), Mul(Integer(2), Symbol('delta_p_plus')), Symbol('d_w_delta_p_minus'))"
    observed=census(synthetic)
    check('complete synthetic pressure census',observed['totalChildren']==3 and len(observed['selected'])==2
          and observed['occurrenceCounts']['delta_p_plus']==1 and observed['maximumPressureDegree']==1)
    check('nonaffine pressure term detected',census("Pow(Symbol('delta_p_plus'), Integer(2))")['maximumPressureDegree']==2)
    check('pressure reciprocal refused',_raises(lambda:census("Pow(Symbol('delta_p_plus'), Integer(-1))")))
    check('unknown pressure function refused',_raises(lambda:census("Function('delta_p_extra')(Symbol('x'))")))
    source_contract=namespace['fourier_contract']
    contract=source_contract(c2)
    check('native executable Fourier factors',contract['profileForwardPower']==-3 and contract['sourceInversePower']==-3
          and contract['coordinates']==3 and contract['nativeProfileDimension']==[3,0,0])
    changed=c2.replace('phase*local_field/(2*sp.pi)**3','phase*local_field/(2*sp.pi)**2')
    check('altered profile measure refused',changed!=c2 and _raises(lambda:source_contract(changed)))
    changed=c2.replace('phase1 * second * local_source / (2*sp.pi)**3','phase1 * second * local_source / (2*sp.pi)**2')
    check('altered second-slot measure refused',changed!=c2 and _raises(lambda:source_contract(changed)))
    changed=c2.replace('NEW_DIMENSIONS[function]=(3,0,0)','NEW_DIMENSIONS[function]=(2,0,0)')
    check('altered profile dimension refused',changed!=c2 and _raises(lambda:source_contract(changed)))
    with tempfile.TemporaryDirectory() as directory:
        present=Path(directory)/'present';present.write_text('ordinary fixture')
        missing=Path(directory)/'missing'
        pins={str(missing):'unavailable',str(present):namespace['sha'](present)}
        records=namespace['posthash_records'](pins)
        check('missing pin preserves subsequent hash',records[str(missing)]['actual'] is None
              and records[str(missing)]['error'] and records[str(present)]['actual']==pins[str(present)])
    # Stand-ins test exact expansion routing and rejection, not SymPy algebra.
    class Atom:
        def __init__(self,arg):self.args=(arg,)
        def __hash__(self):return hash(self.args)
        def __eq__(self,other):return isinstance(other,Atom) and self.args==other.args
    class Expression:
        def __init__(self,atom,scale=1):self.atom=atom;self.scale=scale
        def atoms(self,kind):return {self.atom}
        def xreplace(self,mapping):return Expression(mapping.get(self.atom,self.atom),self.scale)
        def __eq__(self,other):return isinstance(other,Expression) and self.atom==other.atom and self.scale==other.scale
        def __sub__(self,other):return ('unreduced difference',self,other)
    expansion=lambda arg:('expanded',Fraction(1,2),Fraction(-5)) if arg=='factored' else arg
    fake=SimpleNamespace(sinh=Atom,expand=expansion)
    namespace['expanded_sinh_arguments'].__globals__['sp']=fake
    class Trace:
        sinh_zero=journal_type.sinh_zero
        def __init__(self):self.active='parent';self.events=[]
        def emit(self,name,value):self.events.append(name)
        def zero(self,name,left,right):
            self.events.append(name)
            if left!=right:raise ValueError('remaining nonzero')
    left=Expression(Atom('factored'));right=Expression(Atom(expansion('factored')))
    trace=Trace();trace.sinh_zero('identity',left,right)
    check('sinh comparison preserves original operands',left.atom.args==('factored',) and trace.events[:2]==['identity-original-input','identity-original-raw'])
    check('sinh rewrite saved before exact acceptance',trace.events[2:] == ['identity-argument-expansion','identity-canonical'] and trace.active=='parent')
    trace=Trace()
    check('sinh rewrite does not hide unequal scale',_raises(lambda:trace.sinh_zero('unequal',left,Expression(right.atom,2)))
          and trace.active=='unequal' and 'unequal-argument-expansion' in trace.events)
    return {'status':'PASS_STDLIB_TOOLING_ONLY','tests':tests,'count':len(tests),'scientificImport':False,'scientificPayloadRestoration':False,'physicsValidated':False}


def _raises(fn):
    try:fn()
    except ValueError:return True
    return False


if __name__=='__main__':print(json.dumps(run(),indent=2))
