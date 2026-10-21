#!/usr/bin/env python3
"""Exact spherical witness to the information absent from constant sources.

This only checks a single family against the supplied laws. It neither solves
for source forms nor implements the requested S14a build.
"""
import resource
resource.setrlimit(resource.RLIMIT_AS, (1024**3, 1024**3))
import sympy as s

r,t,w,eps,tau = s.symbols('r t w epsilon tau', real=True)
n = 1+s.sin(t)/4
ur = (-r*s.cos(t)/12-r**-2)/n
gamma_ret = gamma_drain = 1/n
radial_div = lambda f: s.diff(r*r*f,r)/(r*r)
shell = lambda f: s.integrate(r*r*f,(r,1,2))
material = s.simplify(s.diff(n,t)+radial_div(n*ur))
order = s.simplify(s.diff(n/2,t)+radial_div(n*ur/2)-n*(gamma_ret-gamma_drain))
chi=s.Rational(1,2)*(s.Heaviside(w)-s.Heaviside(w-1))
global_order=s.expand(s.diff(n*chi,t)+radial_div(n*ur*chi))
face_coefficients=tuple(s.simplify(global_order.coeff(s.DiracDelta(w-face))) for face in (0,1))
storage = s.simplify(s.diff(shell(n),t))
core = s.simplify(-(r*r*n*ur).subs(r,1))
outer = s.simplify((r*r*n*ur).subs(r,2))

print('domain: 1 <= r <= 2; fixed faces w=0,1; all real t; core r<1 excluded')
print('fields: n=1+sin(t)/4; u_r=(-r*cos(t)/12-1/r^2)/n; u_w=0')
print('interior chi=1/2; exterior chi=0; J_chi=0; Gamma_return=Gamma_drain=1/n')
print('density range =',(1-abs(n.coeff(s.sin(t))),1+abs(n.coeff(s.sin(t)))))
print('material regular residual =',material)
print('order regular residual =',order)
print('order face coefficients =',face_coefficients)
print('material face jump coefficients =',tuple(s.simplify(n*0-n*0) for face in (0,1)))
print('source densities =',(s.simplify(n*gamma_ret),s.simplify(n*gamma_drain)))
print('source-density time derivatives =',(s.diff(s.simplify(n*gamma_ret),t),s.diff(s.simplify(n*gamma_drain),t)))
print('return/drain shell totals/(4*pi) =',(shell(n*gamma_ret),shell(n*gamma_drain)))
print('projected density time derivative =',s.diff(n,t))
print('material storage/(4*pi) =',storage)
print('core channel J/(4*pi) =',core)
print('material shell budget residual/(4*pi) =',s.simplify(storage+core+outer))
print('density at t=0, pi/2 =',(n.subs(t,0),n.subs(t,s.pi/2)))

# One possible slow family: the source densities remain exactly the same.
n_eps = 1+s.sin(eps*t)/4
ur_eps = (-eps*r*s.cos(eps*t)/12-r**-2)/n_eps
K_eps = s.simplify(s.diff(shell(n_eps),t))
print('possible slow-family material residual =',s.simplify(s.diff(n_eps,t)+radial_div(n_eps*ur_eps)))
print('fixed-t epsilon->0 density/storage =',(s.limit(n_eps,eps,0),s.limit(K_eps,eps,0)))
print('at t=tau/epsilon density/storage divided by epsilon =',(n_eps.subs(t,tau/eps),s.simplify((K_eps/eps).subs(t,tau/eps))))
