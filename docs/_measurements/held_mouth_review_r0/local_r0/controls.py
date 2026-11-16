from pathlib import Path
from contextlib import redirect_stdout
from io import StringIO
import numpy as np
import sympy as s
import sys
base=Path(__file__).resolve().parent
cand=base.parent/'r0/candidate/software/stage1_solver/tools'
sys.path.insert(0,str(cand));sys.path.insert(1,str(Path('software/stage1_solver/tools').resolve()))
import held_crossflow_budget as B
import unequal_mouth_boundary_work as U
from mouth_feed_weak_check import check_file
# Unequal mouth normalization defect is resolvable using an allowed asymmetric perturbation.
x=np.arange(-1.1,1.11,.1);r=np.arange(.05,1.15,.1)
cut=np.maximum(1-x[:,None]**2-r**2,0)
j=(1+0.9*x[:,None])*cut**3
a=1+.02*cut
f,df,analytic,valid=U.probe(x,r,1.,a,a,j,True)
src=(base/'mutated/unequal_mouth_boundary_work_fixed_rate_fd_omitted.py').read_text()
ns={'__name__':'formwitness'}
exec(compile(src,'unequal_FD_normalization','exec'),ns)
_,_,ref,mutant=ns['probe'](x,r,1.,a,a,j,True)
print("asymmetric trace contraction: normalized analytic/difference, unnormalized difference:",ref,valid,mutant)
assert np.isclose(analytic,valid,rtol=3e-8) and not np.isclose(ref,mutant,rtol=3e-8)
# Local independent derivatives using explicitly varying fields/energy, dimensional c and hbar/m.
y,z=s.symbols('y z',real=True);c=s.Rational(7,11); nh=s.Rational(3,7)
u=1+y/7+(y**2)/20
n=u**2
phase=(1+y+y**3)/5
p=s.Rational(2,3)*(1+y/3)
left=(1+y)/17
N=s.symbols('N')
W=(1-y**2)**3*(N**2/3+N**3/7+y*N**2/9)
Wn=W.diff(N).subs(N,n)
right=left+Wn*u/c
def tensor(v):
    nz=2*u*v
    return nh*n*p*phase.diff(y)+c*(n.diff(y)*nz/n - nz.diff(y))/2
drive=s.cancel(tensor(left)-tensor(right))
explicit=W.diff(y).subs(N,n)
res=s.cancel(drive-s.diff(n*Wn-W.subs(N,n),y)-explicit)
ints=[s.integrate(s.expand(expr),(y,-1,1)) for expr in [drive,explicit]]
flat=s.diff(W,N,y).subs(N,n)*n
frozen=s.cancel(drive-flat)
print("explicit scalar surface polynomial: local chain residual, independent integrated stress, profile derivative",res,ints)
print("freezing n(y) inside d_y e_n at y=0 loses amount",frozen.subs(y,0))
assert res==0 and ints[0]==ints[1]!=0 and frozen.subs(y,0)!=0
for lam in (0,1,2):
    g=(1-y*y)**3; el=lam*c*g*(s.sqrt(N)-1)**2
    en=el.diff(N).subs(N,u**2).subs(s.Abs(u),u)
    # Positive u on support -1<=y<=1.
    lim=(en*u/c).subs(y,0.5)
    u0=float(u.subs(y,.5));g0=float(g.subs(y,.5))
    print('lambda',lam,'normal-slope receiving, direct formula',float(-g0*(u0-1)+lim),(lam-1)*g0*(u0-1))
# Exact first derivative about shape held everywhere vs just around center;
# r^2 integrates with physical measure. Smooth example for literal R(q) level set.
r=s.symbols("r",positive=True);a=s.Integer(1)
n0=s.exp(-r);V=(r-a)**2
slope= -4*s.pi/3*s.integrate(n0*r**2*(V.diff(r,2)+2*V.diff(r)/a),(r,0,s.oo))
print("held mouth q(0)=0 but q_U(w) nonzero inside; level-set geometry slope per q_U",slope)
assert slope!=0
# Test imported on-disk guard readouts one operand at a time on real archived fields.
path=B.archived_packet('mouth_intake_redistribution.json','n32_restoring_slope_u0.01.npz')
fn=Path(check_file.__code__.co_filename)
for label,old,new in [
    ("rate","abs(current.sum()*dx**3 - 0.25)","abs(2*current.sum()*dx**3 - 0.25)"),
    ("aperture","np.array_equal(aperture, mouth)","np.array_equal(aperture, ~mouth)"),
    ("bulk","residual = -lap/2","residual = -lap/3"),
    ("normal_flux","np.max(abs(link_rates-0.25)) >","np.max(abs(2*link_rates-0.25)) >"),
    ("feed_phase","abs(incident-speed) >","abs(-incident-speed) >"),
]:
    source=fn.read_text();assert source.count(old)==1
    out=base/'mutated'/('mouth_feed_weak_check_'+label+'.py')
    out.write_text(source.replace(old,new))
    ms={'__name__': 'field_witness'};exec(compile(out.read_text(),str(out),'exec'),ms)
    try:
        ms['check_file'](path,'restoring_slope')
    except ArithmeticError as e:
        print('imported readout',label,type(e).__name__,str(e));continue
    raise AssertionError(label)
# Test the imported unequal archive packet gating on deliberate fake input, no original file modification.
from unequal_mouth_holding_readout import packet
dummy=base/'forged.npz'
dummy.write_bytes(b"review checksum negative control")
row=dict(guard_exit=0,worker_exit=0,limits_verified=True,fields={str(dummy):"0"*64})
for label,delta in [('guard',{'guard_exit':1}),('worker',{'worker_exit':1}),('limits',{'limits_verified':False}),('checksum',{})]:
    try:packet({**row,**delta},dummy.name)
    except ValueError as e:print('unequal packet',label,type(e).__name__,str(e));continue
    raise AssertionError(label)
