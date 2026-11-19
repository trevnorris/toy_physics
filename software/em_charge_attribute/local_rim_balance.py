"""Conservative lateral junction balance, incremental about a matched intake."""

import sympy as sp


def emit(label, value):
    print(f"{label} = {value}", flush=True)


# Configurational stress for the diffuse order energy on a full 4D cut.
coords = sp.symbols("x y z w", real=True)
x = coords[0]
chi = sp.Function("chi")(*coords)
pot = sp.Function("V")(chi)
kap = sp.symbols("kappa", positive=True)
grad = [sp.diff(chi, t) for t in coords]
energy = pot + kap*sum(g**2 for g in grad)/2
row = [int(j == 0)*energy - kap*grad[0]*g for j,g in enumerate(grad)]
div_stress = sum(sp.diff(s,t) for s,t in zip(row,coords))
chem = sp.diff(pot,chi) - kap*sum(sp.diff(g,t) for g,t in zip(grad,coords))
identity = sp.simplify(div_stress-chem*grad[0])
emit("order-stress divergence minus chemical potential times chi_x", identity)
nu, cr, cw, qp = sp.symbols("nu cr cw qp", real=True)
section = sp.factor(-2*sp.pi*kap*sp.integrate(nu*cr*(cw-qp*nu*cr),(nu,-1,1)))
emit("order load across positive-depth cut per r^2 dr", section)

# Axisymmetric l=1 vector motion with polar axis +x. The measured traces
# are f(a)=normal at the pole, g(a)=tangential along +x at the equator.
r,a,mu,bulk = sp.symbols("r a mu B", positive=True)
f = sp.Function("f")(r)
g = sp.Function("g")(r)
d = sp.diff(f,r) + 2*(f-g)/r
h = (f-g)/r-sp.diff(g,r)
lag = 2*sp.pi*r**2*(bulk*d**2+2*mu*h**2)/3
equations = [sp.diff(sp.diff(lag,sp.diff(v,r)),r)-sp.diff(lag,v) for v in (f,g)]
def insert(expr, fr, gr):
    return sp.factor(expr.subs({f:fr,g:gr}).doit())

u,v,p = sp.symbols("u v p", real=True)
euler = sp.Matrix([insert(eq,u*r**p,v*r**p).subs(r,1) for eq in equations])
det = sp.factor(euler.jacobian([u,v]).det())
emit("radial characteristic determinant",det)
roots = sp.solve(det,p)
emit("radial exponents", sorted(roots))
zn,zt = sp.symbols("Zn Zt", real=True)
fr=gr=0
amps=[]
for j,root in enumerate(sorted(q for q in roots if q < 0)):
    direction, = euler.jacobian([u,v]).subs(p,root).nullspace()
    c=sp.Symbol(f"coef{j}")
    amps.append(c)
    fr+=c*direction[0]*(r/a)**root
    gr+=c*direction[1]*(r/a)**root
constants=sp.solve([fr.subs(r,a)-zn,gr.subs(r,a)-zt],amps)
fr,gr=[sp.factor(val.subs(constants)) for val in (fr,gr)]
traces=[sp.simplify(fr.subs(r,a)-zn),sp.simplify(gr.subs(r,a)-zt)]
emit("wall trace and bulk equation residuals", traces+[insert(eq,fr,gr) for eq in equations])
# At the inner edge the exterior domain normal is -e_r. On any concentric
# cut use the same orientation to compare the force passed out of the rim.
loads=[insert(-sp.diff(lag,sp.diff(shape,r)),fr,gr) for shape in (f,g)]
load=sp.factor(sum(loads))
emit("normal and tangential virtual-work coefficients",loads)
emit("integrated exterior load on any sphere",load)
Kb=sp.factor(load.subs({zn:1,zt:1}))
emit("exterior stiffness for a carried full trace",Kb)
E=sp.factor(sp.integrate(insert(lag,fr,gr),(r,a,sp.oo)))
emit("exterior energy",E)
rank=sp.hessian(E,(zn,zt)).rank()
free,=sp.solve(sp.diff(E,zt),zt)
emit("trace Hessian rank, free tangential trace, exterior load there",
     [rank,free,sp.simplify(load.subs(zt,free))])
emit("zero-load normal-follow exterior field and div/curl",
     [sp.factor(value.subs(zt,free)) for value in (fr,gr,insert(d,fr,gr),insert(h,fr,gr))])

# Joining previous layer stiffness to the exterior; the held geometric rim
# has X=0. Local contact and continuation act on the same near trace Z.
kl,ke,c=sp.symbols("Kl Ke c",positive=True)
z,y,X=sp.symbols("Z Y X",real=True)
Lenergy=(kl*(y-z)**2+ke*z**2+c*(z-X)**2)/2
z0=sp.solve(sp.diff(Lenergy,z).subs(X,0),z)[0]
P,S,H=[sp.factor(expr.subs(z,z0)) for expr in
       (kl*(y-z),ke*z,sp.diff(Lenergy,X).subs(X,0))]
emit("incoming, exterior, external hold", [P,S,H])
emit("local force residual H+P-S",sp.factor(H+P-S))
emit("no contact: transmitted and held", [sp.factor(expr.subs(c,0)) for expr in (S,H)])
emit("pinned: transmitted and held", [sp.limit(expr,c,sp.oo) for expr in (S,H)])
emit("free normal-only exterior and no contact: incoming, transmitted, held",
     [sp.factor(expr.subs(c,0).subs(ke,0)) for expr in (P,S,H)])
assert identity==0 and all(vv==0 for vv in traces+[insert(eq,fr,gr) for eq in equations])
assert all(sp.simplify(sp.diff(E,b)-s)==0 for b,s in zip((zn,zt),loads))
assert sp.factor(H+P-S)==0
