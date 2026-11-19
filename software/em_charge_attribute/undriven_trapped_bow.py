"""Undriven scalar cavity with translated spherical sections and fixed wave action."""

from functools import lru_cache
from pathlib import Path

import numpy as np
import sympy as sp
import yaml
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh
from scipy.optimize import brentq
from scipy.special import eval_legendre, spherical_jn


def out(label, value):
    print(f"{label} = {value}", flush=True)


@lru_cache(None)
def roots(order, count):
    stop = np.pi*(count+order+2)
    grid = np.arange(0.01, stop, 0.25)
    values = spherical_jn(order, grid)
    return [brentq(lambda v: spherical_jn(order, v), x, y, xtol=1e-14)
            for x, y, a, b in zip(grid, grid[1:], values, values[1:]) if a*b < 0][:count]


def straight(R, L, c, gap):
    dipole = roots(1,1)[0]
    wave_numbers = [(np.pi,0), (dipole,0), (np.pi,1), (2*np.pi,0)]
    return [np.sqrt(gap**2+c**2*((kr/R)**2+((j+.5)*np.pi/L)**2))
            for kr,j in wave_numbers]


def gram(left, right, weights):
    return (left*weights) @ right.T


def operators(R, L, c, gap, radial_count, max_order, axial_count):
    # Pull the bowed domain back to the straight ball. The axial derivative is
    # d_z - b f'(z) d_x. Use normalized m=0 harmonics about the bow direction.
    r, wr = leggauss(56)
    r, wr = (r+1)*R/2, wr*R/2
    u, wu = leggauss(24)
    weights = (2*np.pi*wr[:,None]*r[:,None]**2*wu[None,:]).ravel()
    fields, gradients, kr2 = [], [], []
    for order in range(max_order+1):
        leg = eval_legendre(order,u)
        angular = np.sqrt((2*order+1)/(4*np.pi))
        weighted_derivative = (order*(eval_legendre(order-1,u)-u*leg)
                               if order else np.zeros_like(u))
        for zero in roots(order,radial_count):
            k = zero/R
            norm = np.sqrt(2/R**3)/abs(spherical_jn(order+1,zero))
            val = norm*spherical_jn(order,k*r)[:,None]
            dr = norm*k*spherical_jn(order,k*r,derivative=True)[:,None]
            fields.append((angular*val*leg).ravel())
            gradients.append((angular*(dr*u*leg+val*weighted_derivative/r[:,None])).ravel())
            kr2.append(k*k)
    fields, gradients = np.array(fields), np.array(gradients)
    v, wv = leggauss(48)
    z, wz = (v+1)*L/2, wv*L/2
    kz = (np.arange(axial_count)+.5)*np.pi/L
    axial = np.sqrt(2/L)*np.sin(kz[:,None]*z)
    dz = np.sqrt(2/L)*kz[:,None]*np.cos(kz[:,None]*z)
    fp = -np.pi/L*np.sin(np.pi*z/L)
    mixed = np.kron(gram(gradients,fields,weights),gram(axial,dz,wz*fp))
    linear = -c*c*(mixed+mixed.T)
    quadratic = c*c*np.kron(gram(gradients,gradients,weights),gram(axial,axial,wz*fp**2))
    diag = gap**2+c*c*(np.asarray(kr2)[:,None]+kz**2)
    return np.diag(diag.ravel()), linear, quadratic


def check_bow(R, L, c, gap, dims, stiffness, actions):
    base, linear, quadratic = operators(R,L,c,gap,*dims)
    # First excited branch in the sector invariant around the chosen bow axis.
    # The two spectator p harmonics have zero first-order splitting.
    def mode(b, index):
        value, vec = eigh(base+b*linear+b*b*quadratic, subset_by_index=[index,index])
        omega = np.sqrt(value[0])
        dw = vec[:,0] @ (linear+2*b*quadratic) @ vec[:,0]/(2*omega)
        return omega,dw
    p0, s1 = dims[0]*dims[2], 1
    slope = abs(linear[p0,s1])/(2*np.sqrt(base[s1,s1]))
    out(f"volume projection dims {dims}, slope per action", slope)
    ground_change = mode(.08,0)[0]-np.sqrt(base[0,0])
    out("ground frequency change at bow .08", ground_change)
    rows = []
    for action in actions:
        def gradient(b):
            return stiffness*b+action*mode(b,1)[1]
        b = brentq(gradient,1e-7,.15,xtol=1e-13)
        omega, dw = mode(b,1)
        mirror_omega, mirror_dw = mode(-b,1)
        curvature = (gradient(b+1e-5)-gradient(b-1e-5))/2e-5
        energy_change = stiffness*b*b/2+action*(omega-np.sqrt(base[s1,s1]))
        out("action, positive bow, internal load, energy change, curvature",
            [float(x) for x in (action,b,-action*dw,energy_change,curvature)])
        assert abs(gradient(b)) < 1e-10 and abs(omega-mirror_omega)<1e-10
        assert abs(dw+mirror_dw)<1e-8 and energy_change<0 and curvature>0
        rows.append(b)
    assert ground_change>0
    return slope,np.array(rows)


def main():
    archive = Path(__file__).resolve().parents[1]/"stage1_solver/reports/pathA_26_derrick_results.yaml"
    data = yaml.safe_load(archive.read_text())
    R,L = data["phase_A"]["q_star"]
    cfg = data["parameters"]
    c,gap,cut = [cfg[key] for key in ("c_w","mu_in","mu_out")]
    freq = straight(R,L,c,gap)
    out("archived s0 p0 s1 second-radial s0 frequencies", freq)
    zeta = roots(1,1)[0]
    aspect = np.pi*np.sqrt(2/(zeta*zeta-np.pi**2))
    min_length = c*np.sqrt((aspect*zeta)**2+(np.pi/2)**2)/np.sqrt(cut**2-gap**2)
    out("crossing aspect, threshold length, optional slab cap", [aspect,min_length,cfg["d_slab"]])
    # Axial and angular projections in the side-wall shape derivative.
    v, u, phase = sp.symbols("v u phase", real=True)
    axial = 2*sp.sin(sp.pi*v/2)*sp.sin(3*sp.pi*v/2)
    overlap = sp.integrate(sp.expand_trig(axial*sp.cos(sp.pi*v)), (v,0,1))
    rigid = sp.integrate(axial,(v,0,1))
    angular = sp.sqrt(3)*sp.integrate(u**2,(u,-1,1))/2
    mix = sp.Matrix([1,sp.exp(sp.I*phase)])/sp.sqrt(2)
    pressure = sp.simplify(sp.expand_complex((mix.H*sp.Matrix([[0,1],[1,0]])*mix)[0]))
    out("axial bow, rigid translation, angular projection, phase factor", [overlap,rigid,angular,pressure])
    # Separate scale equilibrium test for geometrically similar shapes.
    rad, beta, speed, mass, action = sp.symbols("R beta c m J", positive=True)
    area, volume, curvature = sp.symbols("A V W", nonnegative=True)
    eig = sp.Function("Lambda")(beta)
    spectrum = sp.sqrt(mass**2+speed**2*eig/rad**2)
    geom = area*rad**3+volume*rad**4+curvature*rad
    scale_action = sp.solve(sp.diff(action*spectrum+geom,rad),action)[0]
    load = sp.factor((-action*sp.diff(spectrum,beta)/rad).subs(action,scale_action))
    out("bow load at scale equilibrium for a fixed shape", load)
    num = [*freq, aspect,min_length]
    for radius in (2.,4.):
        length = aspect*radius
        omega = straight(radius,length,c,gap)[1]
        force = c*c*np.pi*zeta/(2*np.sqrt(3)*radius**3*omega)
        # At b -> 0+ the mixed eigenstate has half of each axial occupation.
        axial_force = c*c*5*np.pi**2/(4*omega*length**3)
        out("R, tuned frequency, bow force/J, anchored axial load magnitude/J (ell=1)",
            [radius,omega,force,axial_force])
        num.extend((omega,force,axial_force))
    coarse = check_bow(2.,2*aspect,c,gap,(2,2,4),1.,[.05,.025])
    fine = check_bow(2.,2*aspect,c,gap,(4,3,6),1.,[.05,.025])
    out("largest refinement change in the two equilibrium bows", np.max(abs(coarse[1]-fine[1])))
    print("comparison",*(f"{x:.12g}" for x in num),flush=True)
    assert abs(freq[0]-data["phase_A"]["omega"])<1e-9
    assert overlap==sp.Rational(1,2) and rigid==0 and pressure==sp.cos(phase)
    assert abs(fine[0]-num[7])<1e-10 and np.max(abs(coarse[1]-fine[1]))<1e-4


if __name__ == "__main__":
    main()
