#!/usr/bin/env python3
"""Comprehensive S11c-c1 FACE_RESPONSE_COEFFS cross-engine sweep — CORRECT per-leaf kinematics.

Each joined coefficient leaf is evaluated at its PHYSICAL kinematics: FLAT leaves are the diagonal
half-space object (k=k'), FIRST_SHAPE leaves are the off-diagonal two-momentum kernel (k!=k').  WL emits
the traction as a zero-padded 4-vector (0,0,0,scalar); PY emits the scalar — compared by dropping WL's
identically-zero components.  Prints AGREE / NONZERO / and the traction-padding count.  (An earlier bridge
that used uniform off-diagonal kinematics wrongly flagged the FLAT diagonal objects as mismatches — the
'anchoring split' was a test artifact, corrected here.)

  python3 S11c_c1_reconcile_coeff_sweep.py --runlog <comparator run over FACE_RESPONSE_COEFFS>
"""
import argparse, random, re, sys
import sympy as sp
from sympy import Add, Mul, Pow, Integer, Rational, Symbol, I, Function, Tuple
HFT=Function('HeldInactiveFourierTransform'); qOut=Function('qOut'); w1Profile=Function('w1Profile'); rhoField=Function('rhoBrBgRho4Constant')
for _n in ('chiOne','chiTwo','chiThree'): globals()[_n]=Function(_n)
NS=dict(Add=Add,Mul=Mul,Pow=Pow,Integer=Integer,Rational=Rational,Symbol=Symbol,I=I,Function=Function,Tuple=Tuple,DiracDelta=sp.DiracDelta,HeldInactiveFourierTransform=HFT,qOut=qOut,w1Profile=w1Profile,rhoBrBgRho4Constant=rhoField,chiOne=chiOne,chiTwo=chiTwo,chiThree=chiThree,oo=sp.oo,conjugate=sp.conjugate)
kO=[Symbol('kOne'),Symbol('kTwo'),Symbol('kThree')]; kI=[Symbol('kPrimeOne'),Symbol('kPrimeTwo'),Symbol('kPrimeThree')]
QOUT,QIN,HAT,RHO,MUTH,EPS=(Symbol('QOUT'),Symbol('QIN'),Symbol('HAT'),Symbol('RHO'),Symbol('MUTH'),Symbol('epsilon_shape'))
def plon(e): return e.xreplace({s:Symbol(s.name) for s in e.free_symbols})
def canon(e):
    e=plon(e)
    e=e.replace(lambda x:x.func==qOut and tuple(x.args[1])==tuple(kO),lambda x:QOUT)
    e=e.replace(lambda x:x.func==qOut and tuple(x.args[1])==tuple(kI),lambda x:QIN)
    e=e.replace(lambda x:x.func==HFT,lambda x:HAT); e=e.replace(lambda x:x.func==rhoField,lambda x:RHO)
    r={Symbol('s11cc1_q_out_output'):QOUT,Symbol('s11cc1_q_out_input'):QIN,Symbol('s11cc1_w1_profile_hat_transfer'):HAT,Symbol('rho_br_bg_rho4_constant'):RHO}
    for i in range(3):
        r[Symbol(f's11cc1_k_input_{i+1}')]=kI[i]; r[Symbol(f's11cc1_w1_profile_jet_hat_{i+1}')]=I*Symbol('L_W')*(kO[i]-kI[i])*HAT
    for s in e.free_symbols:
        if 'mu_theta' in s.name: r[s]=MUTH
    return e.xreplace(r)
def comps(x):  # flatten python tuple / sympy Tuple to scalar list
    if isinstance(x,tuple): out=[]; [out.extend(comps(e)) for e in x]; return out
    if getattr(x,'func',None)==Tuple: out=[]; [out.extend(comps(e)) for e in x.args]; return out
    return [x]
def env_for(diagonal):
    R=lambda:Rational(random.randint(3,40),random.randint(2,15))
    env={Symbol(k):R() for k in['L_W','W_0','c_s0','rho_m','eta_bg','Lambda_A_0','Lambda_V_0','Lambda_X_0','tau_A','tau_V','tau_X','epsilon_shape']}
    env[Symbol('sigma_W')]=env[Symbol('eta_bg')]*env[Symbol('W_0')]/env[Symbol('L_W')]; env[Symbol('omega')]=Integer(80)
    kv=[R()/6 for _ in range(3)]
    for i in range(3): env[kO[i]]=kv[i]; env[kI[i]]=kv[i] if diagonal else R()/7
    cs0,om=env[Symbol('c_s0')],env[Symbol('omega')]
    env[QOUT]=sp.sqrt(om**2/cs0**2-sum(env[s]**2 for s in kO)); env[QIN]=sp.sqrt(om**2/cs0**2-sum(env[s]**2 for s in kI))
    env[HAT]=R()+I*R(); env[RHO]=R(); env[MUTH]=R()+I*R()
    return env
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--runlog",required=True); a=ap.parse_args()
    lines=open(a.runlog).read().split("\n"); cases=[]
    for i,l in enumerate(lines):
        if l.startswith("CASE family=FACE_RESPONSE_COEFFS "):
            pa=lines[i+1][12:]; pb=lines[i+2][12:]
            if "<MISSING>" not in pa and "<MISSING>" not in pb: cases.append((l.split("key=")[1],pa,pb))
    agree=nz=incomp=padded=0; nzk=[]
    for key,pa,pb in cases:
        diag=".FLAT." in key
        try:
            A=[EPS*canon(x) for x in comps(eval(pa,NS))]
            B=[canon(x) for x in comps(eval(pb,NS))]
        except Exception: incomp+=1; continue
        Bnz=[x for x in B if x!=0]
        if len(A)==1 and len(Bnz)==1 and len(B)>1: padded+=1  # WL zero-padded vector vs PY scalar
        Buse=Bnz if (len(A)==len(Bnz)) else B
        if len(A)!=len(Buse): incomp+=1; continue
        worst=0; ok=True
        for _ in range(4):
            env=env_for(diag)
            try:
                for aa,bb in zip(A,Buse):
                    miss=(aa.free_symbols|bb.free_symbols)-set(env)
                    if miss: raise KeyError(miss)
                    va=complex(sp.N(aa.xreplace(env))); vb=complex(sp.N(bb.xreplace(env)))
                    worst=max(worst,abs(va-vb)/(abs(va)+abs(vb)+1e-30))
            except Exception: ok=False; break
        if not ok: incomp+=1
        elif worst<1e-9: agree+=1
        else: nz+=1; nzk.append(re.search(r"LEAF=([^)]+)\)",key).group(1))
    print(f"COEFF SWEEP (FLAT=diagonal, FIRST_SHAPE=off-diagonal, traction zero-padding removed):")
    print(f"  AGREE={agree}  NONZERO={nz}  INCOMPLETE={incomp}  of {len(cases)}  (traction zero-padded WL-vec vs PY-scalar: {padded})")
    if nz: 
        from collections import Counter
        print("  NONZERO leaves:", dict(Counter(nzk)))
if __name__=="__main__": main()
