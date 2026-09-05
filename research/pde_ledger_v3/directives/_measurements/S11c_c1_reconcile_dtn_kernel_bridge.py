import sympy as sp
from sympy import Add, Mul, Pow, Integer, Rational, Symbol, I, Function, Tuple
HeldInactiveFourierTransform=Function('HeldInactiveFourierTransform'); qOut=Function('qOut'); w1Profile=Function('w1Profile')
ns=dict(Add=Add,Mul=Mul,Pow=Pow,Integer=Integer,Rational=Rational,Symbol=Symbol,I=I,Function=Function,Tuple=Tuple,
        DiracDelta=sp.DiracDelta,HeldInactiveFourierTransform=HeldInactiveFourierTransform,qOut=qOut,w1Profile=w1Profile,oo=sp.oo)
SCR="/tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/reconcile"
A=eval(open(f"{SCR}/opA.txt").read(),ns); B=eval(open(f"{SCR}/opB.txt").read(),ns)
def plainize(e): return e.xreplace({s:Symbol(s.name) for s in e.free_symbols})
A=plainize(A)               # <-- FIX: strip assumptions on PY side FIRST
kO=[Symbol('kOne'),Symbol('kTwo'),Symbol('kThree')]; kI=[Symbol('kPrimeOne'),Symbol('kPrimeTwo'),Symbol('kPrimeThree')]
LW=Symbol('L_W'); W0=Symbol('W_0'); cs0=Symbol('c_s0'); om=Symbol('omega')
Qout=Symbol('Q_OUT'); Qin=Symbol('Q_IN'); H=Symbol('HAT')
# WL live -> common opaque
B=B.replace(lambda e:e.func==qOut and e.args[1]==Tuple(*kO),lambda e:Qout)
B=B.replace(lambda e:e.func==qOut and e.args[1]==Tuple(*kI),lambda e:Qin)
B=B.replace(lambda e:e.func==HeldInactiveFourierTransform,lambda e:H)
B=plainize(B)
# PY opaque -> common opaque  (+ jet identity + k_input->kPrime)
A_repl={Symbol('s11cc1_q_out_output'):Qout,Symbol('s11cc1_q_out_input'):Qin,Symbol('s11cc1_w1_profile_hat_transfer'):H,
        Symbol('s11cc1_k_input_1'):kI[0],Symbol('s11cc1_k_input_2'):kI[1],Symbol('s11cc1_k_input_3'):kI[2]}
jet_repl={Symbol(f's11cc1_w1_profile_jet_hat_{i+1}'):I*LW*(kO[i]-kI[i])*H for i in range(3)}
A2=A.xreplace({**A_repl,**jet_repl})
r2=sp.simplify(A2-B)
print("=== STAGE 2 (maps applied, assumptions fixed) ===")
print("free:",sorted(str(s) for s in r2.free_symbols)); print("r2 =",r2); print("ZERO?",r2==0)
# STAGE 3: on-shell dispersion  Q_IN^2 = om^2/cs0^2 - |k'|^2 ; Q_OUT^2 = om^2/cs0^2 - |k|^2
kin2=sum(k**2 for k in kI); kout2=sum(k**2 for k in kO)
onshell={Qin**2: om**2/cs0**2 - kin2, Qout**2: om**2/cs0**2 - kout2}
# apply to numerator by substituting squares (use .subs on expanded form)
r3=sp.simplify(sp.together(r2).rewrite(sp.Pow))
num,den=sp.fraction(sp.together(r2))
num_os=sp.expand(num).subs(onshell)
r3=sp.simplify(num_os/den)
print("\n=== STAGE 3 (+ on-shell Q^2 = om^2/cs0^2 - |k|^2) ===")
print("r3 =",r3); print("ZERO?",sp.simplify(r3)==0)
