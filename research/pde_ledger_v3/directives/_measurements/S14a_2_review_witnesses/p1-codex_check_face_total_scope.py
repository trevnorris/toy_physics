#!/usr/bin/env python3
"""A fixed admissible field configuration distinguishes shell total meanings.

Full face weights and rates divided by 4*pi, as in the other diagnostics.
"""
import resource
resource.setrlimit(resource.RLIMIT_AS, (1024**3, 1024**3))
import sympy as s

r,w,t=s.symbols('r w t', real=True)
n=r**-2
ur,uw=s.Integer(-2),s.Integer(-1)
chi=s.Heaviside(w)-s.Heaviside(w-1)
# Extra order flux is zero below/inside, 1/r^2 above: live regional fields,
# with a permitted jump and no concentrated input part.
jw=s.Heaviside(w-1)/r**2
divr=lambda f:s.diff(r*r*f,r)/(r*r)
source=s.expand(s.diff(n*chi,t)+divr(n*chi*ur)+s.diff(n*chi*uw+jw,w))
c_lo=s.simplify(source.coeff(s.DiracDelta(w)))
c_hi=s.simplify(source.coeff(s.DiracDelta(w-1)))
shell=lambda f:s.integrate(r*r*f,(r,1,2))
face_source=shell(c_lo+c_hi)
# An explicit independent source split for this witness, rather than a rule
# assigning positive/negative parts in the general problem.
return_measure=2*s.DiracDelta(w-1)/r**2
drain_measure=s.DiracDelta(w)/r**2
full_return=shell(return_measure.coeff(s.DiracDelta(w-1)))
full_drain=shell(drain_measure.coeff(s.DiracDelta(w)))
# Exterior order flux supplies the closed-projection face term.
face_order_flux=shell(n)
regular_source=s.simplify(source-c_lo*s.DiracDelta(w)-c_hi*s.DiracDelta(w-1))
regular_return=regular_drain=s.Integer(0)

print('domain: 1 <= r <= 2; fixed faces w=0,1; core r<1 excluded')
print('fields: n=1/r^2, u_r=-2, u_w=-1 in all regions; interior chi=1')
print('J_chi^w: below=0, inside=0, above=1/r^2; all other J_chi components=0')
print('material regular residual =',s.simplify(s.diff(n,t)+divr(n*ur)+s.diff(n*uw,w)))
print('material face jump residuals =',tuple(s.simplify(n*uw-n*uw) for face in (0,1)))
print('regular order source =',regular_source)
print('computed lower/upper face source coefficients =',(c_lo,c_hi))
print('chosen full source split minus computed order source =',s.simplify(return_measure-drain_measure-source))
print('regular-only return/drain totals/(4*pi) =',(shell(regular_return),shell(regular_drain)))
print('full-source return/drain totals/(4*pi) =',(full_return,full_drain))
print('separately printed face-source total/(4*pi) =',face_source)
print('net in-plane order flux through both spheres/(4*pi) =',s.simplify((-r*r*n*ur).subs(r,1)+(r*r*n*ur).subs(r,2)))
print('closed-projection order flux through faces/(4*pi) =',face_order_flux)
print('order budget residual with regular totals plus face source =',s.simplify(face_order_flux-face_source))
print('order budget residual with full totals, faces counted within them =',s.simplify(face_order_flux-full_return+full_drain))
print('order budget residual if full totals plus face source are added =',s.simplify(face_order_flux-full_return+full_drain-face_source))
