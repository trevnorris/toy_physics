import sympy as sp, random, sys
from sympy import Add,Mul,Pow,Integer,Rational,Symbol,I,Function,Tuple
HFT=Function('HeldInactiveFourierTransform');qOut=Function('qOut');w1Profile=Function('w1Profile');rhoField=Function('rhoBrBgRho4Constant')
for nm in ('chiOne','chiTwo','chiThree'): globals()[nm]=Function(nm)
ns=dict(Add=Add,Mul=Mul,Pow=Pow,Integer=Integer,Rational=Rational,Symbol=Symbol,I=I,Function=Function,Tuple=Tuple,DiracDelta=sp.DiracDelta,HeldInactiveFourierTransform=HFT,qOut=qOut,w1Profile=w1Profile,rhoBrBgRho4Constant=rhoField,chiOne=chiOne,chiTwo=chiTwo,chiThree=chiThree,oo=sp.oo,conjugate=sp.conjugate)
d=sys.argv[1]; lines=open(f"{d}/light_run.out").read().split("\n")
def plon(e): return e.xreplace({s:Symbol(s.name) for s in e.free_symbols})
kOsym=[Symbol('kOne'),Symbol('kTwo'),Symbol('kThree')];kIsym=[Symbol('kPrimeOne'),Symbol('kPrimeTwo'),Symbol('kPrimeThree')]
QOUT=Symbol('QOUT');QIN=Symbol('QIN');HAT=Symbol('HAT');RHO=Symbol('RHO')
def canon(e):
    e=plon(e)
    e=e.replace(lambda x:x.func==qOut and tuple(x.args[1])==tuple(kOsym),lambda x:QOUT)
    e=e.replace(lambda x:x.func==qOut and tuple(x.args[1])==tuple(kIsym),lambda x:QIN)
    e=e.replace(lambda x:x.func==HFT,lambda x:HAT);e=e.replace(lambda x:x.func==rhoField,lambda x:RHO)
    r={Symbol('s11cc1_q_out_output'):QOUT,Symbol('s11cc1_q_out_input'):QIN,Symbol('s11cc1_w1_profile_hat_transfer'):HAT,Symbol('rho_br_bg_rho4_constant'):RHO}
    for i in range(3):
        r[Symbol(f's11cc1_k_input_{i+1}')]=kIsym[i];r[Symbol(f's11cc1_w1_profile_jet_hat_{i+1}')]=I*Symbol('L_W')*(kOsym[i]-kIsym[i])*HAT
    return e.xreplace(r)
cases=[]
for i,l in enumerate(lines):
    if l.startswith("CASE family=DTN_KERNEL ") and "LEAF=KERNEL_EXPRESSION" in l:
        a=lines[i+1][12:];b=lines[i+2][12:]
        if "<MISSING>" not in a and "<MISSING>" not in b: cases.append((l.split("key=")[1],a,b))
for key,a,b in cases:
    A=canon(eval(a,ns));B=canon(eval(b,ns)); worst=0
    for t in range(8):
        R=lambda:sp.Rational(random.randint(3,40),random.randint(2,15))
        env={Symbol(k):R() for k in['L_W','W_0','c_s0','rho_m','eta_bg']}
        env[Symbol('sigma_W')]=env[Symbol('eta_bg')]*env[Symbol('W_0')]/env[Symbol('L_W')];env[Symbol('omega')]=sp.Integer(80)
        for s in kOsym+kIsym: env[s]=R()/6
        cs0=env[Symbol('c_s0')];om=env[Symbol('omega')]
        env[QOUT]=sp.sqrt(om**2/cs0**2-sum(env[s]**2 for s in kOsym));env[QIN]=sp.sqrt(om**2/cs0**2-sum(env[s]**2 for s in kIsym));env[HAT]=R()+I*R();env[RHO]=R()
        aa=complex(sp.N(A.xreplace(env)));bb=complex(sp.N(B.xreplace(env)));worst=max(worst,abs(aa-bb)/(abs(aa)+abs(bb)+1e-30))
    print(f"  {key[:60]:60} worst_rel={worst:.1e} {'AGREE' if worst<1e-9 else 'MISMATCH'}")
