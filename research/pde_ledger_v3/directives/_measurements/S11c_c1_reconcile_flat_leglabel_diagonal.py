import sympy as sp
from sympy import Add,Mul,Pow,Integer,Rational,Symbol,I,Function,Tuple
HFT=Function('HeldInactiveFourierTransform');qOut=Function('qOut');w1Profile=Function('w1Profile');rhoField=Function('rhoBrBgRho4Constant')
for nm in ('chiOne','chiTwo','chiThree'): globals()[nm]=Function(nm)
ns=dict(Add=Add,Mul=Mul,Pow=Pow,Integer=Integer,Rational=Rational,Symbol=Symbol,I=I,Function=Function,Tuple=Tuple,DiracDelta=sp.DiracDelta,HeldInactiveFourierTransform=HFT,qOut=qOut,w1Profile=w1Profile,rhoBrBgRho4Constant=rhoField,chiOne=chiOne,chiTwo=chiTwo,chiThree=chiThree,oo=sp.oo,conjugate=sp.conjugate)
d="/tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/reconcile"
A=eval(open(f"{d}/ma_A.txt").read(),ns);B=eval(open(f"{d}/ma_B.txt").read(),ns)
def plon(e): return e.xreplace({s:Symbol(s.name) for s in e.free_symbols})
kOsym=[Symbol('kOne'),Symbol('kTwo'),Symbol('kThree')];kIsym=[Symbol('kPrimeOne'),Symbol('kPrimeTwo'),Symbol('kPrimeThree')]
QOUT=Symbol('QOUT');QIN=Symbol('QIN');HAT=Symbol('HAT');RHO=Symbol('RHO');EPS=Symbol('epsilon_shape')
def canon(e):
    e=plon(e)
    e=e.replace(lambda x:x.func==qOut and tuple(x.args[1])==tuple(kOsym),lambda x:QOUT)
    e=e.replace(lambda x:x.func==qOut and tuple(x.args[1])==tuple(kIsym),lambda x:QIN)
    e=e.replace(lambda x:x.func==HFT,lambda x:HAT);e=e.replace(lambda x:x.func==rhoField,lambda x:RHO)
    r={Symbol('s11cc1_q_out_output'):QOUT,Symbol('s11cc1_q_out_input'):QIN,Symbol('s11cc1_w1_profile_hat_transfer'):HAT,Symbol('rho_br_bg_rho4_constant'):RHO}
    for i in range(3):
        r[Symbol(f's11cc1_k_input_{i+1}')]=kIsym[i];r[Symbol(f's11cc1_w1_profile_jet_hat_{i+1}')]=I*Symbol('L_W')*(kOsym[i]-kIsym[i])*HAT
    for s in e.free_symbols:
        if 'mu_theta' in s.name: r[s]=Symbol('MUTH')
    return e.xreplace(r)
A=EPS*canon(A);B=canon(B)
print("A free:",sorted(str(s) for s in A.free_symbols))
print("B free:",sorted(str(s) for s in B.free_symbols))
env={Symbol(k):Rational(p,q) for k,(p,q) in {'L_W':(7,3),'W_0':(5,2),'c_s0':(9,4),'rho_m':(11,5),'eta_bg':(3,2),'Lambda_A_0':(4,3),'Lambda_V_0':(5,7),'Lambda_X_0':(6,5),'tau_A':(2,9),'tau_V':(3,8),'tau_X':(4,7),'epsilon_shape':(5,3)}.items()}
env[Symbol('sigma_W')]=env[Symbol('eta_bg')]*env[Symbol('W_0')]/env[Symbol('L_W')];env[Symbol('omega')]=Integer(80)
for i,s in enumerate(kOsym): env[s]=Rational(i+2,6)
for i,s in enumerate(kIsym): env[s]=Rational(i+3,7)
cs0=env[Symbol('c_s0')];om=env[Symbol('omega')]
env[QOUT]=sp.sqrt(om**2/cs0**2-sum(env[s]**2 for s in kOsym));env[QIN]=sp.sqrt(om**2/cs0**2-sum(env[s]**2 for s in kIsym))
env[HAT]=Rational(7,5)+I*Rational(2,3);env[RHO]=Rational(8,5);env[Symbol('MUTH')]=Rational(3,2)+I*Rational(5,4)
miss=(A.free_symbols|B.free_symbols)-set(env)
if miss: print("MISS",sorted(str(m) for m in miss))
else:
    a=sp.N(A.xreplace(env),50);b=sp.N(B.xreplace(env),50)
    print("a=",a); print("b=",b); print("rel=",sp.N(abs(a-b)/(abs(a)+abs(b)),8))
