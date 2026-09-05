import sympy as sp, random, sys
from sympy import Add,Mul,Pow,Integer,Rational,Symbol,I,Function,Tuple
HFT=Function('HeldInactiveFourierTransform');qOut=Function('qOut');w1Profile=Function('w1Profile');rhoField=Function('rhoBrBgRho4Constant')
for nm in ('chiOne','chiTwo','chiThree'): globals()[nm]=Function(nm)
ns=dict(Add=Add,Mul=Mul,Pow=Pow,Integer=Integer,Rational=Rational,Symbol=Symbol,I=I,Function=Function,Tuple=Tuple,DiracDelta=sp.DiracDelta,HeldInactiveFourierTransform=HFT,qOut=qOut,w1Profile=w1Profile,rhoBrBgRho4Constant=rhoField,chiOne=chiOne,chiTwo=chiTwo,chiThree=chiThree,oo=sp.oo)
d=sys.argv[1]
def plon(e): return e.xreplace({s:Symbol(s.name) for s in e.free_symbols})
kOsym=[Symbol('kOne'),Symbol('kTwo'),Symbol('kThree')]; kIsym=[Symbol('kPrimeOne'),Symbol('kPrimeTwo'),Symbol('kPrimeThree')]
QOUT=Symbol('QOUT');QIN=Symbol('QIN');HAT=Symbol('HAT');RHO=Symbol('RHO'); EPS=Symbol('epsilon_shape')
def canon(expr):
    e=plon(expr)
    e=e.replace(lambda x:x.func==qOut and tuple(x.args[1])==tuple(kOsym),lambda x:QOUT)
    e=e.replace(lambda x:x.func==qOut and tuple(x.args[1])==tuple(kIsym),lambda x:QIN)
    e=e.replace(lambda x:x.func==HFT,lambda x:HAT); e=e.replace(lambda x:x.func==rhoField,lambda x:RHO)
    r={Symbol('s11cc1_q_out_output'):QOUT,Symbol('s11cc1_q_out_input'):QIN,Symbol('s11cc1_w1_profile_hat_transfer'):HAT,Symbol('rho_br_bg_rho4_constant'):RHO}
    for i in range(3):
        r[Symbol(f's11cc1_k_input_{i+1}')]=kIsym[i]; r[Symbol(f's11cc1_w1_profile_jet_hat_{i+1}')]=I*Symbol('L_W')*(kOsym[i]-kIsym[i])*HAT
    return e.xreplace(r)
def test(fA,fB,label,epsA=True,trials=10):
    A=canon(eval(open(f"{d}/{fA}").read(),ns)); B=canon(eval(open(f"{d}/{fB}").read(),ns))
    if epsA: A=EPS*A     # PY coeff * epsilon  (PY baked eps into mu_theta; WL keeps it in coeff)
    worst=0.0; ok=True
    for t in range(trials):
        R=lambda: sp.Rational(random.randint(3,40),random.randint(2,15))
        env={Symbol(k):R() for k in ['L_W','W_0','c_s0','rho_m','eta_bg','Lambda_A_0','Lambda_V_0','Lambda_X_0','tau_A','tau_V','tau_X','epsilon_shape']}
        env[Symbol('sigma_W')]=env[Symbol('eta_bg')]*env[Symbol('W_0')]/env[Symbol('L_W')]; env[Symbol('omega')]=sp.Integer(80)
        for s in kOsym+kIsym: env[s]=R()/6
        cs0=env[Symbol('c_s0')]; om=env[Symbol('omega')]
        env[QOUT]=sp.sqrt(om**2/cs0**2-sum(env[s]**2 for s in kOsym)); env[QIN]=sp.sqrt(om**2/cs0**2-sum(env[s]**2 for s in kIsym))
        env[HAT]=R(); env[RHO]=R()
        a=complex(sp.N(A.xreplace(env))); b=complex(sp.N(B.xreplace(env)))
        rel=abs(a-b)/(abs(a)+abs(b)+1e-30); worst=max(worst,rel); ok=ok and rel<1e-9
    print(f"[{label}] (eps*A) worst_rel={worst:.2e}  {'AGREE' if ok else 'MISMATCH'}")
test("mtc_A.txt","mtc_B.txt","MU_THETA_COEFF")
test("fvc_A.txt","fvc_B.txt","FACE_VELOCITY_COEFF")
