"""Elastic holding-force limits and a separate conserved-stream bending control."""
import argparse

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse.linalg import splu

from mouth_calibrated_flare_coating import read_inputs
from mouth_carrier_lateral_inertia import taper
from mouth_empty_bend_stiffness import INTERVALS
from mouth_flared_receiver import emit
from mouth_middle_tag_grip import cubic_trace
from mouth_replacing_flare_slip import assemble


def minimum(k, free):
    """Impose unit flare displacement and relax the selected tangential slips."""
    y = np.zeros(k.shape[0])
    y[-1] = 1.
    if len(free):
        block = k[free][:,free].tocsc()
        y[free] = splu(block).solve(-k[free,-1].toarray().ravel())
    forces = k @ y
    residual = np.max(abs(forces[free]), initial=0.)
    return y, float(forces[-1]), float(y @ forces), float(residual)


def elastic(mesh, refine, trace):
    patch, base, areas = read_inputs(mesh)
    points, registry, mats = assemble(trace, patch, refine)
    membrane, between, contact = mats[:3]
    slips = np.arange(membrane.shape[0]-1)
    neighbor = slips[slips%4 >= 2]
    bare = float(base[2,2])
    targets = [float(bare+t*areas[2]) for t in (.02,.04)]
    emit(test='fixture', mesh=mesh, refine=refine, stations=len(points), stiffness_order='force per flare displacement',
         carrier_and_wave=bare, area_per_unit_tension=float(areas[2]), candidate_total=targets,
         wall_facing_registry=registry, zero_first_and_middle_displacement=True)
    # This family is positive semidefinite in the wall-facing tangential slips.
    # Restricting that sheet to zero slip supplies an upper comparator for ANY
    # passive quadratic grip on those slips that vanishes at perfect registration.
    modes = [('smooth', 0., slips), ('neutral_registry', 1., slips), ('tenfold_registry', 10., slips),
             ('face_locked', 0., neighbor), ('both_locked', 0., np.array([],dtype=int))]
    previous = None
    for name, factor, free in modes:
        k = membrane + between + factor*contact
        y, modulus, curvature, residual = minimum(k, free)
        steps = np.array([1.e-5,-1.e-5,.5e-5,-.5e-5])
        signed = (k @ (y[:,None]*steps))[-1] + bare*steps
        total = bare+modulus
        emit(test='elastic', attachment=name, coating=modulus, combined_total=total,
             quadratic_curvature=curvature, free_balance=residual, imposed=steps.tolist(),
             holding_forces=signed.tolist(), deficit_below_trial=[v-total for v in targets])
        np.testing.assert_allclose(signed/steps,total,rtol=2e-9)
        np.testing.assert_allclose(modulus,curvature,rtol=2e-8,atol=1e-11)
        if residual>1e-8 or (previous is not None and modulus < previous-1e-8):
            raise ArithmeticError('Elastic relaxation or grip ordering failed')
        previous=modulus


def transport(trace):
    # Keep every passing point at fixed angular label and section depth relative
    # to the given geometry. The bend translates its entire spherical section.
    # This is a separate fully guided material current, with no local drag law.
    lo,hi = INTERVALS[1]
    edges = np.r_[lo,trace.edges[(trace.edges>lo)&(trace.edges<hi)],hi]
    g,weights = leggauss(10)
    half = np.diff(edges)[:,None]/2
    z = edges[:-1,None]+half*(1+g)
    velocity = trace.speed(z)
    dh = -taper(z,lo,hi)[1]
    weak = -np.sum(half*weights*velocity*dh**2)
    # Independent momentum balance over adjacent small slabs, including changes
    # in the discrete carrier tangent across element boundaries.
    for divisions in (16,64):
        stations = np.r_[np.concatenate([np.linspace(a,b,divisions,endpoint=False)
                         for a,b in zip(edges[:-1],edges[1:])]), hi]
        gate, derivative = taper(stations,lo,hi)
        h = 1-gate
        vx = -trace.speed(stations)*derivative   # change in vx per bend displacement
        dv = np.diff(vx)
        generalized = float(((h[:-1]+h[1:])/2) @ dv)
        net = float(dv.sum())
        emit(test='transport',subdivisions=divisions,depth=[lo,hi],
             speed_range=[float(velocity.min()),float(velocity.max())],
             holding_slope_per_mass_current=weak, from_weighted_momentum_changes=generalized,
             net_transverse_momentum_change_per_mass_current=net,
             inward_end_velocity_change_per_displacement=float(vx[0]),
             outward_end_velocity_change_per_displacement=float(vx[-1]))
        if abs(generalized/weak-1) > .002 or abs(net)>1.e-12:
            raise ArithmeticError('Guided stream failed section momentum balance')


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=(26,34),required=True)
    p.add_argument('--refine',type=int,required=True,choices=(1,2,4))
    args=p.parse_args()
    trace,intake=cubic_trace(args.mesh)
    emit(test='intake',intake=intake)
    np.testing.assert_allclose(intake,.25,atol=1.e-8,rtol=0.)
    elastic(args.mesh,args.refine,trace)
    transport(trace)
