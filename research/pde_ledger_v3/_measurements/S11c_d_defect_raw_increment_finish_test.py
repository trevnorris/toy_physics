#!/usr/bin/env python3
"""Stdlib-only source, exact-predicate decision and restoration tests."""
import ast,json,runpy,symtable,tempfile
from pathlib import Path
from types import SimpleNamespace
M=Path(__file__).resolve().parent

def run():
    p=M/'S11c_d_defect_raw_increment_finish.py';source=p.read_text();ns=runpy.run_path(str(p),run_name='tooling_test_only');old=(M/'S11c_d_defect_raw_increment.py').read_text();helper,parts=ns['source_parts'](old)
    tests=[]
    def check(n,v):
        if not v:raise AssertionError(n)
        tests.append(n)
    def raises(f):
        try:f()
        except ValueError:return True
        return False
    compile(source,str(p),'exec');compile(helper,'helper','exec');compile(parts,'tail','exec')
    outer=next(n for n in ast.parse(old).body if isinstance(n,ast.FunctionDef) and n.name=='scientific_work')
    i,loop=next((i,n) for i,n in enumerate(outer.body) if isinstance(n,ast.For) and ast.unparse(n.iter)=="('THETA_BALANCE', 'E_W_BALANCE')")
    theta,ew=parts.body[0].body[:2]
    dump=lambda n:ast.dump(n,include_attributes=False)
    check('EW body exactly original',list(map(dump,ew.body))==list(map(dump,loop.body)))
    check('late jet/domain/end body exactly original',list(map(dump,parts.body[0].body[2:]))==list(map(dump,outer.body[i+1:])))
    check('only unfinished row selections',ast.literal_eval(theta.iter)==('THETA_BALANCE',) and ast.literal_eval(ew.iter)==('E_W_BALANCE',))
    original_without=[n for n in loop.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='without' for t in n.targets)]
    check('completed THETA omission not replayed',len(original_without)==1 and not any(isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='without' for t in n.targets) for n in theta.body))
    good=next(n for n in theta.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='good' for t in n.targets))
    check('completed THETA baseline restored',ast.unparse(good.value)=='saved_good')
    check('completed pre-row responses not called',not any(isinstance(n,ast.Constant) and n.value in ['wrong-lab-face-response','native-slope-omission-response','wrong-sheet-response'] for n in ast.walk(parts)))
    check('no prior construction callable',not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id in ['scientific_work','factorization','lower_boundary','q'] for n in ast.walk(parts)))
    tree=ast.parse(source);main=ast.get_source_segment(source,next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main'))
    check('gate and containment before science',main.index('verify_gate(')<main.index("helpers['containment']()")<main.index('import sympy'))
    check('no top-level science',not any('sympy' in ast.unparse(n) or 'numpy' in ast.unparse(n) for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom))))
    check('no numerical fallback or timer',not any(isinstance(n,ast.Attribute) and n.attr in ['evalf','N','alarm','setitimer'] for n in ast.walk(tree)))
    class Component:
        def __init__(self,label,finite=True,zero=False,real=True):self.label=label;self.is_finite=finite;self.is_zero=zero;self.is_real=real
        def __add__(self,v):return ('sum',self.label,v)
    class I:
        def __mul__(self,v):return ('i',v.label)
    class Den:
        def __init__(self,r,im,residual):self.r=r;self.im=im;self.residual=residual
        def as_real_imag(self):return self.r,self.im
        def __sub__(self,v):assert v==('sum',self.r.label,('i',self.im.label));return self.residual
    class Value:
        def __init__(self,direct=None,free=None,number=True):self.is_finite=direct;self.free_symbols=set(free or []);self.is_number=number
    class Log:
        def __init__(self):self.saved=[]
        def emit(self,name,value):self.saved.append((name,value))
    def fixture(direct=None,free=None,number=True,finite=(True,True,True),zeros=(False,False),residual=0,real=(True,True)):
        n=Component('n',finite[0]);r=Component('r',finite[1],zeros[0],real[0]);im=Component('im',finite[2],zeros[1],real[1]);den=Den(r,im,residual)
        fake=SimpleNamespace(fraction=lambda value:(n,den),cancel=lambda value:value,I=I())
        ns['exact_finite_certificate'].__globals__['sp']=fake;log=Log();record=ns['exact_finite_certificate'](Value(direct,free,number),log,'test')
        assert record['accepted']==ns['certificate_decision'](record)
        assert len(log.saved)==1 and log.saved[0][0]=='test-finite-certificate'
        return record
    check('finite direct accepted',fixture(direct=True)['accepted'])
    check('known nonfinite never overridden',not fixture(direct=False)['accepted'])
    check('exact nonzero real denominator accepted',fixture(zeros=(False,True))['accepted'])
    check('exact nonzero imaginary denominator accepted',fixture(zeros=(True,False))['accepted'])
    check('zero denominator refused',not fixture(zeros=(True,True))['accepted'])
    check('unknown denominator nonzero refused',not fixture(zeros=(None,None))['accepted'])
    check('unknown component finiteness refused',not fixture(finite=(True,None,True))['accepted'])
    check('nonfinite numerator refused',not fixture(finite=(False,True,True))['accepted'])
    check('unbound expression refused',not fixture(free=['x'])['accepted'])
    check('nonnumber refused',not fixture(number=False)['accepted'])
    check('failed reconstruction refused',not fixture(residual=1)['accepted'])
    check('nonreal denominator components refused',not fixture(real=(None,True))['accepted'])
    # Namespace inventory without executing any scientific expression.
    fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run_finish')
    env=next(n.value for n in fn.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='env' for t in n.targets));keys={ast.literal_eval(k) for k in env.keys}
    mapping=next(n.value for n in fn.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='mapping' for t in n.targets));keys.update(ast.literal_eval(v) for v in mapping.values)
    import builtins
    st=symtable.symtable(ast.unparse(parts),'remaining','exec').get_children()[0];needed={s.get_name() for s in st.get_symbols() if s.is_global() and s.is_referenced()}
    check('remaining-control namespace complete',not (needed-keys-set(vars(builtins))))
    with tempfile.TemporaryDirectory() as d:
        p=Path(d);oldp=p/'old';oldp.mkdir();(oldp/'x').write_bytes(b'opaque\x00saved');inventory={'root':str(oldp),'files':{'x':{'bytes':12,'sha256':ns['sha'](oldp/'x')}}};ip=p/'inventory';ip.write_text(json.dumps(inventory));out=p/'new';out.mkdir()
        target,files=ns['copy_prior']({'priorComplete':str(oldp),'priorInventory':str(ip)},out)
        check('entire prior copied byte-identically',(target/'x').read_bytes()==(oldp/'x').read_bytes())
        (oldp/'extra').write_text('not pinned');out2=p/'new2';out2.mkdir()
        check('prior census mismatch refused',raises(lambda:ns['copy_prior']({'priorComplete':str(oldp),'priorInventory':str(ip)},out2)))
    return {'status':'PASS_STDLIB_TOOLING_ONLY','count':len(tests),'tests':tests,'scientificImports':False,'scientificRestoration':False,'independentPhysicsValidation':False}

if __name__=='__main__':print(json.dumps(run(),indent=2))
