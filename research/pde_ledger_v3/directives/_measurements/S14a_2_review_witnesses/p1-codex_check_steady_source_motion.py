#!/usr/bin/env python3
"""Small admissible shell witness, not a construction of the proposed build.

All rates below are divided by the common angular factor 4*pi. Distribution
moments use the stipulated full weight at each closed-slab endpoint.
"""
import resource
resource.setrlimit(resource.RLIMIT_AS, (1024**3, 1024**3))
import sympy as s

r, t, w = s.symbols('r t w', real=True)
lo, hi = s.Integer(0), 1 + t/4
n, ur, uw = r**-2, s.Integer(-2), s.Integer(-1)
radial_div = lambda f: s.diff(r*r*f, r)/(r*r)
shell = lambda f: s.integrate(r*r*f, (r, 1, 2))
material = s.diff(n,t) + radial_div(n*ur) + s.diff(n*uw,w)
a_lo, a_hi = uw-s.diff(lo,t), uw-s.diff(hi,t)
# Directly differentiate the ordered current. Exterior J_chi is zero here.
chi = s.Heaviside(w-lo)-s.Heaviside(w-hi)
order_density, order_radial, order_w = n*chi, n*chi*ur, n*chi*uw
order_source = s.expand(s.diff(order_density,t)+radial_div(order_radial)+s.diff(order_w,w))
low_atom = s.simplify(order_source.coeff(s.DiracDelta(w)))
high_atom = s.simplify(order_source.coeff(s.DiracDelta(w-hi)))
ret, drain = high_atom, -low_atom
ret_total, drain_total = shell(ret), shell(drain)
ret_first_w_moment = shell(ret*hi)
sigma, flux = n*(hi-lo), n*ur*(hi-lo)
storage = s.diff(shell(sigma),t)
face = shell(-n*a_lo+n*a_hi)
core = -(r*r*flux).subs(r,1)
outer = (r*r*flux).subs(r,2)

print('domain: 1 <= r <= 2; t >= 0; 0 <= w <= 1+t/4; core r<1 excluded')
print('fields: n=1/r^2 in all regions, u_r=-2, u_w=-1, interior chi=1, J_chi=0')
print('material regular residual =', material)
print('material face jump residuals =', (s.simplify(n*a_lo-n*a_lo),s.simplify(n*a_hi-n*a_hi)))
print('order-source lower/upper delta coefficients =', (low_atom,high_atom))
print('return total/(4*pi) =', ret_total)
print('drain total/(4*pi) =', drain_total)
print('time derivatives of those totals =', (s.diff(ret_total,t),s.diff(drain_total,t)))
print('return first w moment/(4*pi) =', ret_first_w_moment)
print('time derivative of return first w moment/(4*pi) =', s.diff(ret_first_w_moment,t))
print('projected material storage/(4*pi) =', storage)
print('material shell budget residual/(4*pi) =', s.simplify(storage+face+core+outer))
print('order shell budget residual/(4*pi) =', s.simplify(storage+core+outer-ret_total+drain_total))

# Converse witness: stationary ordinary source densities over a moving slab
# need not yield constant integrated rates, even without any source atoms.
n2, ur2 = s.Integer(1), -2/r**2
f2 = w*(hi-w)/hi**2
primitive2 = s.integrate(f2,w)
j2w = -s.diff(primitive2,t)
j2top = s.simplify(j2w.subs(w,hi))
regular_material2 = s.diff(n2,t)+radial_div(n2*ur2)
regular_order2 = s.simplify(s.diff(f2,t)+radial_div(f2*ur2)+s.diff(j2w,w))
source_density2 = s.Integer(1)  # Both independent inputs n*Gamma_return and n*Gamma_drain.
each_total2 = shell(source_density2*(hi-lo))
atom2_lower = s.simplify(f2.subs(w,lo)*(-s.diff(lo,t))+j2w.subs(w,lo))
atom2_upper = s.simplify(-f2.subs(w,hi)*(-s.diff(hi,t))+j2top-j2w.subs(w,hi))
print('converse witness: n=1, u_r=-2/r^2, u_w=0, chi_inside=w*(h-w)/h^2, h=1+t/4')
print('converse J_chi^w inside =',s.factor(j2w))
print('converse J_chi^w below/above =',(0,j2top))
print('converse material/order regular residuals =',(s.simplify(regular_material2),regular_order2))
print('converse order face coefficients =',(atom2_lower,atom2_upper))
print('converse each source density/time derivative =',(source_density2,s.diff(source_density2,t)))
print('converse each return/drain total/(4*pi) =',each_total2)
print('converse time derivative of each total/(4*pi) =',s.diff(each_total2,t))
