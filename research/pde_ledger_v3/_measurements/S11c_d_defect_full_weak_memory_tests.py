#!/usr/bin/env python3
"""Exact stdlib stand-ins and source metadata only; never restores scientific payloads."""
import ast, copy, hashlib, json, tempfile, unittest
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace
M=Path('/var/projects/toy_physics/research/pde_ledger_v3/_measurements')
W=M/'S11c_d_defect_full_weak.py';source=W.read_text();tree=ast.parse(source)
names={'require','sha','memory_domain_certificate','verify_build_assessment'}
ns={'json':json,'Path':Path,'hashlib':hashlib}
exec(compile(ast.Module([n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names],type_ignores=[]),'actual-exact-predicate-functions','exec'),ns)
class P:
    def __init__(self,x=0):
        self.d=x.copy() if isinstance(x,dict) else {(x,):Fraction(1)} if isinstance(x,str) else {():Fraction(x)}
        self.d={k:v for k,v in self.d.items() if v}
    def __add__(self,o):
        o=as_p(o);d=self.d.copy()
        for k,v in o.d.items():d[k]=d.get(k,0)+v
        return P(d)
    __radd__=__add__
    def __neg__(self):return P({k:-v for k,v in self.d.items()})
    def __sub__(self,o):return self+-as_p(o)
    def __rsub__(self,o):return as_p(o)+-self
    def __mul__(self,o):
        d={}
        for a,x in self.d.items():
            for b,y in as_p(o).d.items():
                k=a+b;n=k.count('I');k=tuple(sorted(v for v in k if v!='I'))+('I',)*(n%2)
                k=tuple(sorted(k));d[k]=d.get(k,0)+x*y*(-1)**(n//2)
        return P(d)
    __rmul__=__mul__
    def __pow__(self,n):
        assert type(n) is int and n>=0
        v=P(1)
        for _ in range(n):v=v*self
        return v
    def __eq__(self,o):return isinstance(o,(P,int,Fraction)) and self.d==as_p(o).d
    @property
    def name(self):
        assert len(self.d)==1
        k,v=next(iter(self.d.items()));assert len(k)==1 and v==1
        return k[0]
    @property
    def free_symbols(self):return [P(n) for n in sorted({n for k in self.d for n in k if n!='I'})]
    def subs(self,s,v):
        out=P()
        for k,c in self.d.items():out+=P({tuple(n for n in k if n!=s.name):c})*as_p(v)**k.count(s.name)
        return out

def as_p(x):return x if isinstance(x,P) else P(x)
def terms(p,om,unused):
    out={}
    for k,v in p.d.items():
        key=(k.count(om.name),k.count(unused.name));a=P({tuple(n for n in k if n not in (om.name,unused.name)):v});out[key]=out.get(key,P())+a
    return out
sp=SimpleNamespace(I=P('I'),Symbol=P,cancel=lambda p:p)
def bind(p):
    for s in p.free_symbols:
        if s.name in {'omega','tau_A','tau_X','tau_V','W_0','L_W','rho_br'}:p=p.subs(s,{'omega':3,'tau_A':Fraction(1,5),'tau_X':Fraction(1,4),'tau_V':Fraction(1,3),'W_0':2,'L_W':10,'rho_br':3}[s.name])
    return p

def constant(v,allow_zero):
    if not isinstance(v,P) or v.free_symbols or (not allow_zero and v==0):raise ValueError('nonfinite/unbound/zero')
    return {'finite':True,'nonzero':v!=0}
def zero(rec,label,a,b):
    rec['residual']=a-b
    if a!=b:raise ValueError(label)
def call(p,binder=bind):
    events=[]
    result=ns['memory_domain_certificate'](p,binder,constant,terms,zero,lambda k,v:events.append((k,copy.deepcopy(v))),sp)
    return result,events
class Tests(unittest.TestCase):
    def test_actual_forms(self):
        o=P('omega');i=sp.I
        for tau in ('tau_A','tau_X','tau_V'):
            for pref in (P(1),P('rho_br'),P('W_0')*P('rho_br'),P('L_W')**2*P('W_0')**4*P('rho_br')**2):
                for power in (1,2):
                    r,e=call(pref*(o*P(tau)+i)**power)
                    self.assertEqual(r['power'],power);self.assertEqual(r['residual'],0)
                    self.assertEqual([x[0] for x in e],['memory-domain-input','memory-sign-operands','memory-domain-return'])
    def test_original_predicate_rejects_known_prefactor_and_square(self):
        names={'omega','rho_br','tau_A','W_0'}
        self.assertFalse(names<={'omega','rho_br','tau_A'})
        d=P('rho_br')**2*(P('omega')*P('tau_A')+sp.I)**2
        a0=d.subs(P('omega'),0)
        self.assertNotEqual(d,a0*(1-sp.I*P('omega')*P('tau_A')))
    def test_nonmemory_constants(self):
        for p in (P('W_0'),P('L_W')):self.assertNotIn('power',call(p)[0])
    def test_wrong_time_sign_refused(self):
        with self.assertRaises(ValueError):call(P('omega')*P('tau_A')-sp.I)
    def test_mixture_not_single_memory_power(self):
        with self.assertRaises(ValueError):call((P('omega')*P('tau_A')+sp.I)**2+P('omega')*P('tau_A'))
    def test_two_taus_refused(self):
        with self.assertRaises(ValueError):call(P('omega')*P('tau_A')+P('tau_X')+sp.I)
    def test_unknown_prefactor_refused(self):
        with self.assertRaises(ValueError):call(P('other')*(P('omega')*P('tau_A')+sp.I))
    def test_unbound_or_nonfinite_refused(self):
        for b in (lambda p:p,lambda p:float('inf')):
            with self.assertRaises(ValueError):call(P('omega')*P('tau_A')+sp.I,b)
    def test_zero_bound_refused(self):
        with self.assertRaises(ValueError):call(P('omega')*P('tau_A')+sp.I,lambda p:P())
    def test_zero_prefactor_refused(self):
        with self.assertRaises(ValueError):call(P('omega')*P('tau_A'))
    def test_uncensused_power_refused(self):
        with self.assertRaises(ValueError):call((P('omega')*P('tau_A')+sp.I)**3)
    def test_operands_saved_before_wrong_sign_failure(self):
        events=[]
        with self.assertRaises(ValueError):ns['memory_domain_certificate'](P('omega')*P('tau_A')-sp.I,bind,constant,terms,zero,lambda k,v:events.append(k),sp)
        self.assertEqual(events,['memory-domain-input','memory-sign-operands'])
    def test_scientific_body_unchanged_except_predicate_call(self):
        old=ast.parse((Path('/var/projects/toy_physics/_scratch/s11c/s11c-defect-full-weak-20261002/build-review-r2/packet/worker.py')).read_text())
        a=next(n for n in old.body if isinstance(n,ast.FunctionDef) and n.name=='run_science');b=copy.deepcopy(next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run_science'))
        def neutralize(f):
            f.body=[n for n in f.body if not (isinstance(n,ast.For) and ast.unparse(n.target)=='(i, den)' and 'raw_denominators.values()' in ast.unparse(n.iter))]
            return ast.dump(f)
        self.assertEqual(neutralize(a),neutralize(b))
    def test_other_definitions_unchanged(self):
        old=ast.parse((Path('/var/projects/toy_physics/_scratch/s11c/s11c-defect-full-weak-20261002/build-review-r2/packet/worker.py')).read_text())
        before={n.name:ast.dump(n) for n in old.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name not in ('run_science','verify_gate')}
        after={n.name:ast.dump(n) for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in before}
        self.assertEqual(before,after)
    def test_review_repair_identity_guards(self):
        review={'methodAssessed':True,'allChecksPassed':True,'independentBuildClearance':False,'reports':{'claude':{'literalVerdict':'CLEAR FOR THIS BOUNDED FULL-WEAK BUILD'},'grok':{'literalVerdict':'NEEDS REVISION'}},'methodSha256':'method','sharedGuardSha256':'guard','supervisorSha256':'supervisor','workerSha256':'oldworker','manifestSha256':'oldmanifest','launcherSha256':'oldlauncher'}
        repair={'toolingOnly':True,'testsPassed':True,'noScientificPayloadRestored':True,'reviewRecordSha256':'review','reviewed':review,'workerSha256':'worker','manifestSha256':'manifest','launcherSha256':'launcher','methodSha256':'method','evidencePins':{}}
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'repair';p.write_text(json.dumps(repair))
            gate={'repairRecord':str(p),'repairRecordSha256':ns['sha'](p),'buildReviewRecordSha256':'review','workerSha256':'worker','manifestSha256':'manifest','launcherSha256':'launcher','independentBuildClearance':False,'localToolingRepairAccepted':True,'sharedGuardSha256':'guard','supervisorSha256':'supervisor'}
            ns['verify_build_assessment'](gate,review)
            for key in ('repairRecordSha256','buildReviewRecordSha256','workerSha256','manifestSha256','launcherSha256','sharedGuardSha256','supervisorSha256','localToolingRepairAccepted','independentBuildClearance'):
                g=copy.deepcopy(gate);g[key]='wrong'
                with self.assertRaises(ValueError,msg=key):ns['verify_build_assessment'](g,review)
if __name__=='__main__':unittest.main()
