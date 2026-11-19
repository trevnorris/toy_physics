"""Weak, cyclic constituent detector reading the outgoing pair's signed twist.

Use a square in the equatorial xy plane, one constituent on each face. An
individual bond can see more than the net flux. The cyclically averaged
circulating current and translation-invariant preparation are essential.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_pair_translating_cuts import AMPLITUDE, Pair, hankel, polynomial


def fields(pair, coeff, k, x, y):
    # Exact background and linear harmonic amplitudes (peak phasors) on z=0.
    cyl = np.hypot(x,y)
    dc = np.zeros_like(cyl)
    dz = np.zeros_like(cyl,dtype=complex)
    angle = np.zeros_like(dz)
    for i,z in enumerate(pair.centers):
        r = np.hypot(cyl,z)
        co, si = -z/r, cyl/r
        dc += pair.monopoles[i]*co/r**2
        for j in range(pair.lmax):
            l = j+1
            p,t = polynomial(l,co)
            dc += pair.static[i*pair.lmax+j]*((l+1)*p*co-t*si)/r**(l+2)
            h = hankel(l,k)
            f = hankel(l,k*r)/h
            df = k*hankel(l,k*r,True)/h
            amplitude = coeff[i*pair.lmax+j]
            dz += amplitude*(co*l*(l+1)*p*f/r + si*(df+f/r)*t)
            angle += amplitude*f*t
    return dc, dz, -angle*y/cyl, angle*x/cyl


def loop_flux(pair, coeff, k, r, side):
    grid, weights = leggauss(6)
    weights *= side/2
    grid *= side/2
    x,y = np.meshgrid(r+grid,grid,indexing="ij")
    dc,dz,_,_ = fields(pair,coeff,k,x,y)
    flux = weights@dz@weights
    low_y = fields(pair,coeff,k,r+grid,-side/2*np.ones_like(grid))[2]
    high_y = fields(pair,coeff,k,r+grid,side/2*np.ones_like(grid))[2]
    low_x = fields(pair,coeff,k,r-side/2+0*grid,grid)[3]
    high_x = fields(pair,coeff,k,r+side/2+0*grid,grid)[3]
    circulation = weights@(low_y+high_x-high_y-low_x)
    return weights@dc@weights,flux,circulation


def ring_response(dc, omega):
    # Unit t, retaining each local bond and site port. The scalar readout is
    # -dH/dF with F split uniformly over the links for its definition only.
    size = 4
    eye = np.eye(size)
    h = np.diag((16.*(1-eye)).ravel()).astype(complex)
    first, second = [], []
    for i in range(size):
        j = (i+1)%size
        cols = [np.zeros_like(h) for _ in range(3)]
        for face in (1,-1):
            bond = np.zeros_like(eye,dtype=complex)
            bond[j,i] = -np.exp(1j*face*dc/size)
            operators = (bond+bond.conj().T, 1j*face*(bond-bond.conj().T), -(bond+bond.conj().T))
            for col,op in zip(cols,operators):
                col += np.kron(op,eye) if face==1 else np.kron(eye,op)
        h += cols[0]
        first.append(cols[1])
        second.append(cols[2])
    e,v = np.linalg.eigh(h)
    ground = v[:,0]
    gap = e[1:]-e[0]
    j_op = sum(first)/size
    j_mg = v[:,1:].conj().T@(j_op@ground)
    def response(port,direct):
        port_mg = v[:,1:].conj().T@(port@ground)
        mix = j_mg.conj()*port_mg
        return (ground.conj()@direct@ground - sum(mix/(gap-omega)+mix.conj()/(gap+omega)))
    local = np.array([response(op,h2/size) for op,h2 in zip(first,second)])
    site = [response(np.kron(np.diag(row),eye),0*j_op) for row in eye]
    # K=d^2 E/dF^2 in the slow limit, opposite sign for current -dH/dF.
    return float(local.mean().real),max(abs(local-local.mean())),max(abs(np.array(site)))


def main():
    k=.12
    pair = Pair(np.array([-2.,2.]),10,8.)
    relative = pair.emit(k,[AMPLITUDE,-AMPLITUDE])
    common = pair.emit(k,[AMPLITUDE,AMPLITUDE])
    previous = None
    for kr in (8.,16.,32.,64.):
        r = kr/k
        dc,flux,line = loop_flux(pair,relative,k,r,.1)
        _,control,_ = loop_flux(pair,common,k,r,.1)
        kernel,phase_spread,clock_response = ring_response(dc,.02)
        current = -kernel*flux
        err = abs(line-flux)/abs(flux)
        print(f"kr={kr:g} Fdc={dc:.10g} F_relative_abs={abs(flux):.12g} "
              f"F_common_abs={abs(control):.3g} K_over_t={kernel:.12f} I_ac_abs_over_t={abs(current):.12g} "
              f"stokes_rel={err:.3g} local_phase_spread={phase_spread:.3g} site_response={clock_response:.3g}")
        assert err<1e-8 and abs(control)<1e-10*abs(flux)
        assert kernel>.001 and phase_spread<1e-12 and clock_response<1e-12
        if previous is not None:
            ratio = current/previous*np.exp(-1j*kr/2)
            print(f"kr={kr:g} double_distance_amplitude={abs(ratio):.9f} "
                  f"phase_after_travel_radians={np.angle(ratio):.9f}")
            if kr==64:
                assert abs(abs(ratio)-.5)<.002 and abs(np.angle(ratio))<.025
        previous = current


if __name__ == "__main__":
    main()
