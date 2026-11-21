"""Bound the differential shift of any two orthogonal material rows in 4D."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from mouth_contact_age import PRIOR, emit, trace_path
from mouth_free_contact import Carrier


def cell(fields, spacing, resolution):
    t,r,w,q = trace_path(fields, resolution)
    carrier = Carrier(Path(fields["carrier_file"]).parent)
    offsets = np.array([[0,0],[-.5,0],[.5,0],[0,-.5],[0,.5]])
    initial = np.array([r[0],w[0]])+spacing*offsets
    def flow(time, points):
        rr,ww = points.reshape(-1,2).T
        radius = carrier.radius(ww)[0]
        _,_,_,vr,vw = carrier.eval(rr/radius, ww)
        return np.column_stack((vr,vw)).ravel()
    sol = solve_ivp(flow,(0.,t[-1]),initial.ravel(),t_eval=t,rtol=2e-9,atol=2e-11,max_step=t[-1]/100)
    if not sol.success:
        raise ValueError(sol.message)
    positions = sol.y.T.reshape(-1,5,2)
    F = np.zeros((len(t),4,4))
    F[:,:2,0] = (positions[:,2]-positions[:,1])/spacing
    F[:,:2,1] = (positions[:,4]-positions[:,3])/spacing
    F[:,2,2] = F[:,3,3] = positions[:,0,0]/positions[0,0,0]
    return F,q


def examine(fields, spacing, resolution, fraction_only=False):
    F,q = cell(fields, spacing, resolution)
    # In the deformed cell at time s, split the neighbor vector into its
    # projection on the row and the perpendicular remainder. For the relative
    # deformation A from s to exit, the extra projection is bounded by
    # (cond(A)-1/cond(A))/2 times the remainder-to-row length ratio.
    condition = np.linalg.cond(F[:-1])
    rel_cond = np.linalg.cond(F[-1] @ np.linalg.inv(F[:-1]))
    shift_bound = .5*(rel_cond-1/rel_cond)*condition
    upper = 2*shift_bound/(q[:-1]-q[-1])
    if fraction_only:
        # Interpolate the positive metric C and q linearly between samples.
        # For any fixed row and neighbor, d=(n C h)/(n C n) then lies
        # between its endpoint values. Allow both edges of the exit interval
        # to lose the largest possible overhang observed on any segment.
        bounds = np.r_[shift_bound,0.]
        overhang = np.maximum(bounds[:-1],bounds[1:])-(q[1:]-q[-1])/2
        lower_fraction = max(0.,1-2*max(0.,overhang.max())/q[-1])
        sampled_overhang = shift_bound-(q[:-1]-q[-1])/2
        lower_at_samples = max(0.,1-2*max(0.,sampled_overhang.max())/q[-1])
        emit(profile=fields["profile"],state=fields["state"],mesh=fields["mesh"],
             spacing=spacing,steps_per_radius=resolution,minimum_retention_bound_samples=lower_at_samples,
             minimum_retention_bound_interpolated=lower_fraction, max_overhang=overhang.max(),
             exit_order=q[-1])
        assert 0. <= lower_fraction <= lower_at_samples+1e-12 <= 1.+1e-12
        return
    # Lower bound on the maximum: for each initial row n, choose the
    # maximizing unit neighbor perpendicular to n, with no loss of its sign.
    # Both angular eigenvalues agree, so one angular component suffices here.
    meridian, polar = np.meshgrid(np.linspace(0,np.pi,101),np.linspace(0,np.pi,51))
    axes = np.vstack([np.sin(polar).ravel()*np.cos(meridian).ravel(),
                      np.sin(polar).ravel()*np.sin(meridian).ravel(),np.cos(polar).ravel()])
    metric = F.transpose(0,2,1) @ F
    # Also sample reflected radial/axial orientations, consequential when
    # the meridional off-diagonal of C is nonzero.
    axes = np.hstack([axes,axes*np.array([[-1.],[1.],[1.]])])
    def projected(C):
        v = C[:3,:3] @ axes
        return v/np.sum(v*axes, axis=0)
    final = projected(metric[-1])
    largest = 0.
    smallest_survival_example = 1.
    for M,order in zip(metric[:-1],q[:-1],strict=True):
        difference = projected(M)-final
        magnitudes = np.linalg.norm(difference,axis=0)
        magnitude = magnitudes.max()
        largest = max(largest, 2*magnitude/(order-q[-1]))
        # Common length of just these two intervals upper-bounds continuous
        # survival for the corresponding fixed row/neighbor orientation.
        survival_upper = min(1.,max(0.,((order+q[-1])/2-magnitude)/q[-1]))
        smallest_survival_example = min(smallest_survival_example,survival_upper)
    emit(profile=fields["profile"],state=fields["state"],mesh=fields["mesh"],
         spacing=spacing,steps_per_radius=resolution,lower_on_max_shift_over_margin=largest,
         upper_on_max_shift_over_margin=upper.max(),
         smallest_two_section_survival_upper=smallest_survival_example,
         min_determinant=np.linalg.det(F).min(), transit_points=len(q))
    assert np.linalg.det(F).min()>0 and largest <= upper.max()+1e-7


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retained-fraction", action="store_true")
    args = parser.parse_args()
    archive = json.loads(PRIOR.read_text())["runs"]["mouth_mixed_order_ripple02"]
    for line in archive["stdout_lines"]:
        p=json.loads(line)
        if p["test"]=="archive_patch":
            examine(p,.01,1000,args.retained_fraction)
            if p["mesh"]==48:
                examine(p,.005,2000,args.retained_fraction)
