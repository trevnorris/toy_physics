import sympy as sp, random, sys
from sympy import Add,Mul,Pow,Integer,Rational,Symbol,I,Function,Tuple
HFT=Function('HeldInactiveFourierTransform'); qOut=Function('qOut'); w1Profile=Function('w1Profile'); rhoField=Function('rhoBrBgRho4Constant')
for nm in ('chiOne','chiTwo','chiThree'): globals()[nm]=Function(nm)
ns=dict(Add=Add,Mul=Mul,Pow=Pow,Integer=Integer,Rational=Rational,Symbol=Symbol,I=I,Function=Function,Tuple=Tuple,DiracDelta=sp.DiracDelta,HeldInactiveFourierTransform=HFT,qOut=qOut,w1Profile=w1Profile,rhoBrBgRho4Constant=rhoField,chiOne=chiOne,chiTwo=chiTwo,chiThree=chiThree,oo=sp.oo)
d=sys.argv[1]
def plon(e): return e.xreplace({s:Symbol(s.name) for s in e.free_symbols})
kOsym=[Symbol('kOne'),Symbol('kTwo'),Symbol('kThree')]; kIsym=[Symbol('kPrimeOne'),Symbol('kPrimeTwo'),Symbol('kPrimeThree')]
QOUT=Symbol('QOUT'); QIN=Symbol('QIN'); HAT=Symbol('HAT'); RHO=Symbol('RHO')
def canon(expr):
    e=plon(expr)
    # WL function apps -> placeholders (match on symbolic args)
    e=e.replace(lambda x:x.func==qOut and tuple(x.args[1])==tuple(kOsym), lambda x:QOUT)
    e=e.replace(lambda x:x.func==qOut and tuple(x.args[1])==tuple(kIsym), lambda x:QIN)
    e=e.replace(lambda x:x.func==HFT, lambda x:HAT)
    e=e.replace(lambda x:x.func==rhoField, lambda x:RHO)
    # PY opaque -> placeholders + jet identity
    r={Symbol('s11cc1_q_out_output'):QOUT, Symbol('s11cc1_q_out_input'):QIN,
       Symbol('s11cc1_w1_profile_hat_transfer'):HAT, Symbol('rho_br_bg_rho4_constant'):RHO}
    for i in range(3):
        r[Symbol(f's11cc1_k_input_{i+1}')]=kIsym[i]
        r[Symbol(f's11cc1_w1_profile_jet_hat_{i+1}')]=I*Symbol('L_W')*(kOsym[i]-kIsym[i])*HAT
    return e.xreplace(r)
def test(fA,fB,label,trials=10):
    A=canon(eval(open(f"{d}/{fA}").read(),ns)); B=canon(eval(open(f"{d}/{fB}").read(),ns))
    freeA=A.free_symbols; freeB=B.free_symbols
    worst=0.0; ok=True
    for t in range(trials):
        R=lambda: sp.Rational(random.randint(3,40), random.randint(2,15))
        env={s:R() for s in ['L_W','W_0','c_s0','rho_m','eta_bg','Lambda_A_0','Lambda_V_0','Lambda_X_0','tau_A','tau_V','tau_X','epsilon_shape']}
        env=dict((Symbol(k),v) for k,v in env.items())
        env[Symbol('sigma_W')]=env[Symbol('eta_bg')]*env[Symbol('W_0')]/env[Symbol('L_W')]  # binding
        env[Symbol('omega')]=sp.Integer(80)
        for s in kOsym+kIsym: env[s]=R()/6
        cs0=env[Symbol('c_s0')]; om=env[Symbol('omega')]
        env[QOUT]=sp.sqrt(om**2/cs0**2 - sum(env[s]**2 for s in kOsym))
        env[QIN]=sp.sqrt(om**2/cs0**2 - sum(env[s]**2 for s in kIsym))
        env[HAT]=R(); env[RHO]=R()
        miss=(freeA|freeB)-set(env)
        if miss: print(f"  [{label}] UNSUBSTITUTED: {sorted(str(m) for m in miss)}"); ok=False; break
        a=complex(sp.N(A.xreplace(env))); b=complex(sp.N(B.xreplace(env)))
        rel=abs(a-b)/(abs(a)+abs(b)+1e-30); worst=max(worst,rel)
        if rel>1e-9: ok=False
    print(f"[{label}] trials={trials} worst_rel={worst:.2e}  {'AGREE (residual==0 under identifications)' if ok else 'MISMATCH/incomplete'}")
test("opA.txt","opB.txt","DTN_KERNEL LH F1 (seals 2,3)")
test("mtc_A.txt","mtc_B.txt","FACE_RESP_COEFFS MU_THETA_COEFF (seals 1,5)")
