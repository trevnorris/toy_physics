"""Material labels and transverse stability of a cold connected strip on the cubic wall.

The crystalline tangent plane and affine axial/hoop transport are the hypotheses.
A perturbation shared by sites at the same meridional position is an invariant
subspace of the triangular lattice. Within that subspace all sites can move.
"""
import argparse

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq

from mouth_middle_tag_grip import cubic_trace
from mouth_order_wave_pulse import emit
from mouth_segment_contact_hold import PITCH
from mouth_two_depth_load import DEPTHS
from neutral_patch_sites import WITHIN, morse

WAIST = 0.5
MIDDLE = float(DEPTHS[1])
# Wavelengths in unstretched meridional lattice spacings.
WAVES = np.array([4., 8., 16., 32.])


def bond_derivatives(r):
    depth, alpha, rest = WITHIN
    w = np.exp(-alpha*(r-rest))
    tension = morse(r, *WITHIN)[1]
    return tension, 2*depth*alpha**2*w*(2*w-1)


def lattice(n):
    m, j = np.mgrid[-n:n+1, -n:n+1]
    x = m + j/2
    y = j*np.sqrt(3.)/2
    mask = (m != 0) | (j != 0)
    return np.column_stack([x[mask], y[mask]])


def natural_pitch():
    length = np.linalg.norm(lattice(12), axis=1)
    def pressure(p):
        return np.dot(length, bond_derivatives(p*length)[0])
    return brentq(pressure, .90, 1.02, xtol=1e-14)


def arc_speed(trace, w):
    return trace.speed(w)*np.hypot(1., trace.geometry.radius(w)[1])


def path_grid(trace, subdivisions):
    edges = np.unique(np.r_[WAIST, trace.edges[(trace.edges>WAIST)&(trace.edges<MIDDLE)], MIDDLE])
    return np.r_[np.concatenate([np.linspace(a,b,subdivisions,endpoint=False)
                                 for a,b in zip(edges[:-1],edges[1:])]), MIDDLE]


def integrals(fun, edges):
    x, w = leggauss(6)
    h = np.diff(edges)
    return np.cumsum(h/2 * (fun(edges[:-1,None]+(1+x)*h[:,None]/2) @ w))


def rates(bonds, phases, axial, hoop):
    x, y = bonds.T
    xx, yy = x*axial, y*hoop
    r = np.hypot(xx, yy)
    tension, curvature = bond_derivatives(r)
    side = tension/r
    yy_hessian = side + (curvature-side)*(yy/r)**2
    xy_hessian = (curvature-side)*xx*yy/r**2
    return yy_hessian @ phases, xy_hessian @ phases


def run(args):
    trace, intake = cubic_trace(args.mesh)
    w = path_grid(trace, args.sub)
    reference = WAIST if args.feed=='waist' else MIDDLE
    p = natural_pitch()
    # One native Morse length = 0.01 mouth radii, chosen as in the preceding model.
    spacing = PITCH*p
    u_ref = arc_speed(trace, reference)
    period = spacing/u_ref
    arc = np.r_[0., integrals(lambda v: np.hypot(1.,trace.geometry.radius(v)[1]),w)]
    time = np.r_[0., integrals(lambda v: 1/trace.speed(v),w)]
    u = arc_speed(trace,w)
    axial = u/u_ref
    hoop = trace.geometry.radius(w)[0]/trace.geometry.radius(reference)[0]
    # At the material crossing event: distances to the preceding and following tags,
    # using exact time-of-flight coordinates. Uniform leads outside the sampled band.
    def at(t):
        return (np.interp(t,time,arc)+np.minimum(t,0)*u[0]
                +np.maximum(t-time[-1],0)*u[-1])
    read = []
    for t in (0.,time[-1]):
        s = at(t)
        read.append([float(s-at(t-period)),float(at(t+period)-s)])
    # Independently count a common integer label sequence crossing both gates.
    ids = np.arange(32)
    gates = [time[0], time[-1]]
    measured=[]
    for gate in gates:
        hits=[brentq(lambda t: at(t-i*period)-at(gate), gate+(i-.5)*period,
                     gate+(i+.5)*period) for i in ids]
        measured.append(float(np.mean(np.diff(hits))))
    np.testing.assert_allclose(measured,period,rtol=2e-11)
    np.testing.assert_allclose(time[-1],trace.transit(WAIST,MIDDLE,8),rtol=1e-10)
    emit(test='arrivals',intake=intake,mesh=args.mesh,reference=args.feed,
         natural_pitch_in_patch=p,native_length_in_R=PITCH,local_arc_speed=[u[0],u[-1]],
         period_from_crossings=measured,rate=1/period,
         local_pitch_in_R=(period*u[[0,-1]]).tolist(),
         simultaneous_tag_gaps_in_R=read,waist_middle_time=time[-1],
         hoop_middle_over_waist=hoop[-1]/hoop[0],axial_middle_over_waist=axial[-1]/axial[0])

    # Retain a fixed symmetric list containing every neighbor that can approach
    # within the nominal range anywhere along this transport. Extra corners are kept.
    ix=int(np.ceil(args.range/(p*axial.min())))+2
    iy=int(np.ceil(args.range/(p*hoop.min()*np.sqrt(3)/2)))+2
    j, doubled_x = np.mgrid[-iy:iy+1,-2*ix:2*ix+1]
    keep=((j-doubled_x)%2==0)&((j!=0)|(doubled_x!=0))
    reduced_x=doubled_x[keep]/2
    bonds=p*np.column_stack([reduced_x,j[keep]*np.sqrt(3)/2])
    phases=1-np.cos(2*np.pi*reduced_x[:,None]/WAVES)
    modes=np.array([rates(bonds,phases,a,b) for a,b in zip(axial,hoop)])
    lam=modes[:,0,:]
    if np.max(abs(modes[:,1,:]))>2e-8*max(1.,np.max(abs(lam))):
        raise ArithmeticError('Mode mixes axial and transverse response')
    gain=-cumulative_trapezoid(lam,time,axis=0,initial=0)
    emit(test='stability',reference=args.feed,range=args.range,sub=args.sub,neighbors=len(bonds),
         waves=WAVES.tolist(),rates_at_waist=lam[0].tolist(),rates_at_middle=lam[-1].tolist(),
         min_rate=lam.min(axis=0).tolist(),max_rate=lam.max(axis=0).tolist(),
         log_gain_end=gain[-1].tolist())
    for col,wave in enumerate(WAVES):
        neg=np.flatnonzero(lam[:,col]<0)
        # Observable is the crest-to-trough relative offset, divided by hoop stretch.
        # A fresh infinitesimal seed may enter at any point. Record its largest
        # log growth and first tenfold growth, before extrapolation to nonlinear slip.
        floor=np.minimum.accumulate(gain[:,col])
        best=gain[:,col]-floor
        ten=np.flatnonzero(best>=np.log(10))
        out=dict(test='tracked_mode',wave=wave,max_log_amplification=float(max(best)))
        if len(neg):
            i=neg[0]
            out.update(first_negative_depth_bracket=[w[max(0,i-1)],w[i]],first_negative_time=time[i])
        if len(ten):
            k=ten[0]; start=int(np.argmin(gain[:k+1,col]))
            out.update(tenfold_seed_depth=w[start],tenfold_depth_bracket=[w[k-1],w[k]],
                       elapsed_to_tenfold=time[k]-time[start])
        emit(**out)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',type=int,choices=[26,34],default=34)
    parser.add_argument('--feed',choices=['waist','middle'],required=True)
    parser.add_argument('--sub',type=int,default=32)
    parser.add_argument('--range',type=float,default=5.)
    run(parser.parse_args())
