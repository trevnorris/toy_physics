"""Signed tilt drive from round drain stresses and a displaced wave confinement."""

import sympy as sp


def emit(name, value):
    print(f"{name} = {value}", flush=True)


def round_mouth():
    mu = sp.symbols("mu", real=True)  # projection of radial unit vector on lab x
    r, w, nr, nw, speed = sp.symbols("r w_relative n_r n_w v_cut", real=True)
    rr, rw, ww, pr, pw = sp.symbols("T_rr T_rw T_ww p_r p_w", real=True)
    tr = nr*rr + nw*rw - speed*pr
    tw = nr*rw + nw*ww - speed*pw
    for_direction = r*mu*tw - w*mu*tr
    average = lambda expr: sp.integrate(expr, (mu, -1, 1))/2
    emit("angular average of moving-cut mixed angular momentum flux", average(for_direction))
    dipole = average(mu*for_direction)
    emit("dipole-weighted sample divided by local radial torque factor",
         sp.factor(dipole/(r*tw-w*tr)))
    # Under a lateral wall displacement f(w) the derivative of any radial
    # confinement is proportional to mu. Density is radial and P^2 ~ mu^2.
    emit("off-axis wall force angular factors, carrier and unseeded wave",
         (average(mu), average(mu**3)))
    assert average(for_direction)==0 and sp.simplify(3*dipole-(r*tw-w*tr))==0


def splitting_wave():
    x, w, d, b0, b1, b2, k, angle = sp.symbols("x w d b0 b1 b2 k beta", real=True)
    h = w**2-sp.Rational(1,2)
    shift = b0 + b1*w + b2*h
    potential = (x-shift)**2 + (w**2-d**2)**2/4
    # Actual modes are sqrt(2/pi) times these polynomials times exp(-r^2/2).
    p, t = x, h
    def moment(n):
        return 0 if n%2 else sp.gamma(sp.Rational(n+1,2))
    def inner_integral(poly):
        return sp.factor(2*sum(coef*moment(i)*moment(j)
                               for (i,j),coef in sp.Poly(poly,x,w).terms())/sp.pi)
    def grad(poly, axis):
        return sp.diff(poly, axis) - axis*poly
    def operator(poly):
        return -sum(grad(grad(poly,z),z) for z in (x,w)) + potential*poly

    assert [inner_integral(v) for v in (p*p,t*t,p*t)] == [1,1,0]
    pp, pt, tt = [inner_integral(u*operator(v)) for u,v in ((p,p),(t,p),(t,t))]
    emit("instantaneous target acceleration from a unit parent", -pt)
    emit("fixed-wall difference of target and parent Rayleigh quotients",
         sp.factor((tt-pp).subs({b0:0,b1:0,b2:0})))
    energy = (pp*sp.cos(angle)**2 + 2*pt*sp.sin(angle)*sp.cos(angle)
              + tt*sp.sin(angle)**2)/2 + k*b2**2/2
    energy = energy.subs({b0:0,b1:0})
    origin = {b2:0,angle:0}
    slope = [energy.diff(z).subs(origin) for z in (b2,angle)]
    curvature = sp.simplify(sp.hessian(energy,(b2,angle)).subs(origin))
    relaxed = sp.factor(curvature[1,1]-curvature[1,0]**2/curvature[0,0])
    emit("wall and tilt slopes", slope)
    emit("wall and tilt curvature matrix", curvature)
    emit("tilt curvature with the lateral wall relaxed", relaxed)
    emit("fixed-wall tilt curvatures at d=0 and d=2",
         tuple(curvature[1,1].subs(d,depth) for depth in (0,2)))
    reflection = sp.trigsimp(energy-energy.xreplace({angle:-angle,b2:-b2}))
    emit("energy difference on reversing both lateral bow and wave tilt", reflection)
    assert -pt == b2 and slope==[0,0] and reflection==0


if __name__=="__main__":
    round_mouth()
    splitting_wave()
