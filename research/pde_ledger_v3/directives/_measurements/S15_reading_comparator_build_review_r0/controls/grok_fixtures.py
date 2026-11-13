"""Run the remote reviewer's synthetic controls with the frozen native readers locally."""
import json
from pathlib import Path
import sys
import sympy as s
base=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(base/'packet/research/pde_ledger_v3/scripts'))
import S15_reading_math as m
ans=json.loads((base/'grok_list/answer.txt').read_text())['exact_results']
inner,outer=ans['radii']['inner'],ans['radii']['outer']
assert all(type(r) is int and r>0 for r in (inner,outer)) and inner<outer

def surface(matrix,r,engine):
    n=('sin(theta)*cos(phi)','sin(theta)*sin(phi)','cos(theta)') if engine=='py' else ('Sin[theta]*Cos[phi]','Sin[theta]*Sin[phi]','Cos[theta]')
    assert all(type(a) is int for v in matrix for a in v)
    expression='+'.join(f'({a})*({n[i]})*({n[j]})' for i,v in enumerate(matrix) for j,a in enumerate(v) if a)
    measure='sin(theta)' if engine=='py' else 'Sin[theta]'
    body=f'{r**3}*{measure}*({expression})'
    return (f'Integral({body},(theta,0,pi),(phi,0,2*pi))' if engine=='py'
            else f'Inactive[Integrate][{body},{{theta,0,Pi}},{{phi,0,2*Pi}}]')

def exact(x):
    # These known trigonometric-polynomial fixtures admit direct exact integration.
    return s.simplify(x.doit()/s.pi)
for i,mat in enumerate(ans['matrices']):
    pieces={e:[surface(mat,r,e) for r in (outer,inner)] for e in ('py','wl')}
    x,y=[m.parse(f'{hi}-({lo})',e) for e,(hi,lo) in pieces.items()]
    bad=m.parse('+'.join(pieces['wl']),'wl')
    left, right, flipped=map(exact,(x,y,bad))
    delta=exact(m.exact_residual(x,bad))
    print('outward flux / pi, component',i+1,'SymPy',left,'Wolfram',right,'Wolfram inner sign flipped',flipped,'signed difference',delta,flush=True)
    assert left==right==s.Rational(ans['reference_outer_plus_inner_in_units_of_pi'][i])
    assert flipped==s.Rational(ans['one_engine_inner_sign_flip_outer_plus_inner_in_units_of_pi'][i])
    assert delta==left-flipped and delta!=0
mats=[row[:] for row in ans['matrices'][0]]
mats[0][0]=2
a=m.parse(surface(ans['matrices'][0],outer,'py'),'py')
b=m.parse(surface(mats,outer,'wl'),'wl')
delta=exact(m.exact_residual(a,b))
print('outer component 1 diagonal change, reference minus changed / pi',delta,flush=True)
assert delta==s.Rational(ans['structure_change']['reference_minus_changed'])
py=m.parse('Piecewise(((a_1**2-b_1**2)/(a_1-b_1), Ne(a_1,b_1)), (2*a_1,True))','py')
wl=m.parse('Piecewise[{{aa[1],aa[1]==bb[1]}},(aa[1]^2-bb[1]^2)/(aa[1]-bb[1])]','wl')
for bval in (s.Rational(1,3),s.Rational(1,2)):
    bind={s.Symbol('a_1'):s.Rational(1,2),s.Symbol('b_1'):bval}
    x,y=py.subs(bind,simultaneous=True),wl.subs(bind,simultaneous=True)
    print('exponents',str(bind),'SymPy',x,'Wolfram changed',y,'signed difference',x-y,flush=True)
    key='generic_oriented_residual' if bval!=s.Rational(1,2) else 'joint_oriented_residual'
    assert x-y==s.Rational(ans['piecewise'][key])
print('The exponent fixture uses one equality hyperplane. Its codimension is one.',flush=True)
