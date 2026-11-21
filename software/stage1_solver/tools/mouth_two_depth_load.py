"""Two independent force drives and three material-shape readings at fixed baseline intake.

No added wall law. The growing sinusoid is switched on from rest. Evolving
its complex envelope gives the cosine and sine quadratures in one solve.
"""
import argparse
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse.linalg import splu

from mouth_carrier_lateral_response import prepare, surface, trace
from mouth_empty_bend_stiffness import shapes
from mouth_flared_receiver import emit
from mouth_free_bend_drive import Response
from mouth_free_bend_release import pencil


DEPTHS = np.array([.525, 1.175, 2.25])
PEAK_FORCE = 1e-9


def polar(values):
    """Gain and signed lag for the exp(+i omega t) convention."""
    return dict(gain=abs(values).tolist(), lag_degrees=(-np.angle(values, deg=True)).tolist())


def cut_flux(dom, field, perturbation, depth):
    # A constant-depth cut intersects a ball in the three brane coordinates.
    # Each chosen cut slides rigidly in x. Its normal speed in w is zero.
    vgrid = dom.blocks[2][1]
    nodes, weights = leggauss(8)
    half = np.diff(vgrid)[:, None]/2
    v = (vgrid[:-1,None]+half*(1+nodes)).ravel()
    w = np.full_like(v, depth)
    fields = np.column_stack([field, perturbation])
    r, _, _, rv, _, _, psi, radial, axial = trace(dom, 2, w, v, fields)
    f, h = psi.T
    fr, hr = radial.T
    fw, hw = axial.T
    # Decompose the off-diagonal gradient stress into flow and amplitude terms.
    a = abs(f)
    ar, aw = a*(fr/f).real, a*(fw/f).real
    ur, uw = (fr/f).imag, (fw/f).imag
    g = h/f
    gr, gw = (hr-g*fr)/f, (hw-g*fw)/f
    b = a*g.real
    br, bw = ar*g.real+a*gr.real, aw*g.real+a*gw.real
    advective = a*a*((gr.imag+2*g.imag/r)*uw+ur*gw.imag)+2*a*b*ur*uw
    amplitude = (br+2*b/r)*aw+ar*bw
    measure = (4*np.pi/3)*r*r*abs(rv)*(half*weights).ravel()
    split = np.array([measure@advective, measure@amplitude])
    direct, isotropic = surface(dom, fields, 2, 'v', depth, (0.,1.))
    np.testing.assert_allclose(split.sum(), direct, rtol=1e-9, atol=1e-11)
    if isotropic != 0:
        raise ArithmeticError('Expected only off-diagonal stress at a constant-depth cut')
    return np.r_[direct, split]


def observe_cuts(dom, field, response, state, columns):
    parts = np.zeros((2*len(field), columns),complex)
    parts[response.dipole.free] = state[:-3,:columns]
    real, imag = np.split(parts, 2)
    for depth in DEPTHS:
        readings = []
        for i in range(columns):
            cosine = cut_flux(dom, field, real[:,i].real+1j*imag[:,i].real, depth)
            sine = cut_flux(dom, field, real[:,i].imag+1j*imag[:,i].imag, depth)
            readings.append(cosine+1j*sine)
        answer = np.array(readings)
        emit(test='cut', depth=float(depth), channels=['total','advective','amplitude_gradient'],
             per_peak_force_real=answer.real.tolist(), per_peak_force_imag=answer.imag.tolist(),
             **polar(answer[:,0]))


def run(args):
    dom, field = prepare(args)
    response = Response(dom, field, args.extension)
    energy, time = pencil(response, 0.)
    radius, _ = dom.radius(DEPTHS)
    bends, sr, sw = shapes(radius, DEPTHS, args.extension, dom.far, dom.roof)
    read = np.column_stack((np.ones(3),bends))
    np.testing.assert_array_equal(read,np.tri(3))
    np.testing.assert_array_equal(sr, np.zeros_like(sr))
    np.testing.assert_array_equal(sw, np.zeros_like(sw))
    # Columns are unit pushes on the middle and flare. The third is a
    # withheld, phase-shifted half push on the entrance, different from both.
    ports = np.column_stack([read[1], read[2], -.5j*read[0]])
    rhs = np.zeros((energy.shape[0],3),complex)
    rhs[-3:] = ports
    rate = .30+.40j
    operator = energy-rate*time
    particular = splu(operator.tocsc()).solve(rhs)
    steady_error = np.max(abs(operator@particular-rhs))
    emit(test='instrument', depths=DEPTHS.tolist(), forcing_depths=DEPTHS[1:].tolist(),
         envelope_growth=rate.real, omega=rate.imag, peak_force=PEAK_FORCE, dt=args.dt,
         third_drive='entrance, -i/2', field_balance=steady_error)
    # Physical load for -T<t<=0 is Re(F exp(s t)). Normalizing by its
    # current exponential removes the small starting amplitude. State starts
    # from zero. No unknown tail or periodically balanced initial state is used.
    dt = args.dt
    advance = splu((operator*.5-time/dt).tocsc())
    state = np.zeros_like(rhs)
    for i in range(round(96/dt)):
        state += advance.solve(rhs-operator@state)
        elapsed = (i+1)*dt
        if elapsed in (64.,96.):
            measured = read@state[-3:]
            exact = read@particular[-3:]
            relative = np.linalg.norm(measured-exact,axis=0)/np.linalg.norm(exact,axis=0)
            emit(test='motion', elapsed=elapsed, per_peak_force_real=measured.real.tolist(),
                 per_peak_force_imag=measured.imag.tolist(), **polar(measured),
                 switch_on_relative_error=relative.tolist())
    observe_cuts(dom, field, response, state, 2)
    # At this complex frequency the exact transformed connector load is
    # F-Z_car Q. Keep the residual from finite switch-on as a control.
    known = response.force_matrix(rate)
    for name, experiment in [('finite_record',state),('particular',particular)]:
        coordinates = experiment[-3:]
        residual = ports-known@coordinates
        offsets = coordinates[1:,:2]
        condition = np.linalg.cond(offsets)
        connection = np.linalg.solve(offsets.T,residual[1:,:2].T).T
        predicted = known@coordinates[:,2]+np.r_[0.,connection@coordinates[1:,2]]
        error = np.max(abs(predicted-ports[:,2]))
        emit(test='inverse',data=name,relative_offset_condition=condition,
             common_force_residual=float(max(abs(residual[0]))),
             force_per_peak_residual=float(np.max(abs(residual))),
             extra_force_matrix_real=connection.real.tolist(),extra_force_matrix_imag=connection.imag.tolist(),
             held_out_force_error=float(error))
        if name=='particular' and max(error,np.max(abs(residual)))>1e-7:
            raise ArithmeticError('Independent reduced and full field equations disagree')
    if max(relative)>1e-5 or steady_error>1e-7 or not np.isfinite(condition) or condition>1e4:
        raise ArithmeticError('Switch-on, field residual or independent depth excitation failed')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--extension',choices=['compact','wide'],required=True)
    p.add_argument('--dt',type=float,choices=[.25,.125],required=True)
    run(p.parse_args())
