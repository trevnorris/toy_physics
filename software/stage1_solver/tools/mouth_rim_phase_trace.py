"""Local carrier phase and number conservation on a tilted, permeable graph."""
import sympy as sp


x,y,z,t,b = sp.symbols('x y z t b', real=True)  # t=w-h; b=local second depth derivative
coords=(x,y,z)
r2=x*x+y*y+z*z
height=(1-r2)*(1+x/10)/15
surface_phase=(x*x+y*y+z*z)/9+x/30+x*y/12
number=1+x/50+y*y/40
depth_gradient=(1+y)/30
q=(1-r2)**3*(1+y/9)/3

def grad(f):
    return sp.Matrix([sp.diff(f,i) for i in coords])

def lap(f):
    return sum(sp.diff(f,i,2) for i in coords)

slope,g=grad(height),grad(surface_phase)
metric=1+slope.dot(slope)
a=(q/number+slope.dot(g))/metric
parallel=g-a*slope
b_required=(2*slope.dot(grad(a))+a*lap(height)-lap(surface_phase)
            -(grad(number).dot(parallel)+depth_gradient*q/number)/number)/metric
# Independently take laboratory derivatives of a local phase expansion to
# evaluate its normal current and divergence at the physical point w=h.
phase=surface_phase+t*a+t*t*b/2
n=number+t*depth_gradient

def at_fixed_w(f,i):
    return f.diff(i)-height.diff(i)*f.diff(t)

jx=sp.Matrix([n*at_fixed_w(phase,i) for i in coords]); jw=n*phase.diff(t)
flux=(jw-jx.dot(slope)).subs(t,0)
divergence=(sum(at_fixed_w(jx[j],i) for j,i in enumerate(coords))+jw.diff(t)).subs(t,0)
for point in ((0,0,0),(sp.Rational(1,3),sp.Rational(1,4),sp.Rational(-1,5)),(1,0,0),(-1,0,0)):
    p=dict(zip(coords,point,strict=True)); target=q.subs(p)
    normal=flux.subs(p); coefficient=b_required.subs(p)
    residual=divergence.subs(p).subs(b,coefficient)
    uncorrected=divergence.subs(p).subs(b,0)
    print(f'point={point} prescribed={float(target):.12g} crossing={float(normal):.12g} '
          f'j_w={float(jw.subs(t,0).subs(p)):.12g} '
          f'horizontal_projection={float(jx.dot(slope).subs(t,0).subs(p)):.12g} '
          f'divergence_before={float(uncorrected):.12g} after={float(residual):.12g}')
    assert normal==target and residual==0 and uncorrected!=0
