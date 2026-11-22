"""Same clamped flare with only cross-layer attractions weakened by conversion."""
import argparse

import numpy as np

from mouth_calibrated_flare_coating import read_inputs, surface_rule
from mouth_empty_bend_stiffness import INTERVALS
from mouth_middle_tag_grip import cubic_trace
from mouth_segment_contact_hold import PITCH
from neutral_patch_periodic_compression import Bilayer
from neutral_patch_sites import GAP, emit, morse


def material_tangents(patch, trace):
    cell = Bilayer(1,1,True,patch['shell'])
    a, b = (cell.vectors*patch['pitch']).T
    r = np.sqrt(a*a + b*b + GAP**2*cell.cross)
    depth, alpha, length = cell.params
    _, first = morse(r,depth,alpha,length)
    y = np.exp(-alpha*(r-length))
    second = 2*depth*alpha**2*y*(2*y-1)
    scale = patch['shared_scale']*2*GAP*PITCH
    readings = {name: [] for name in ('persistent','ordered_fraction','cross_only_weakens','circumferential')}
    for i, span in enumerate(INTERVALS, start=1):
        area, h_s, radial, order = surface_rule(trace,span)
        # Average (a*T_x + b*t_x)^2 on the round section. T is the
        # meridional tangent, t is a uniformly oriented angular tangent.
        projection2 = (np.outer(radial**2, a*a) + b*b)/3
        # Second derivative of a single Morse energy under the actual displacement.
        # Retain first/r: separate bond families carry stress at the common pitch.
        per_bond = scale*np.outer(area*h_s**2,a*a)*(
            first/r + (second-first/r)*projection2/r**2)
        same = per_bond[:,~cell.cross].sum(axis=1)
        cross = per_bond[:,cell.cross].sum(axis=1)
        readings['persistent'].append(float((same+cross).sum()))
        readings['ordered_fraction'].append(float(order@(same+cross)))
        readings['cross_only_weakens'].append(float(same.sum()+order@cross))
        readings['circumferential'].append(0.)
        emit(test='selective_bonds',bend=i,same_layer=same.sum(),cross_when_ordered=cross.sum(),
             cross_when_weakened=order@cross,material_tangent=same.sum()+order@cross)
    return readings


def compare(mesh):
    patch, baseline, area = read_inputs(mesh)
    trace, _ = cubic_trace(mesh)
    for case, values in material_tangents(patch,trace).items():
        # Bend gradients have disjoint support. Common translation leaves all
        # local lengths unchanged. Only the carrier and wave couple the bends.
        stiffness = baseline + np.diag([0.,*values])
        if stiffness[0,0] <= 0:
            raise ArithmeticError('Finite-bath static elimination requires positive spring')
        internal = stiffness[1:,1:] - np.outer(stiffness[1:,0],stiffness[0,1:])/stiffness[0,0]
        sigma = values[1]/area[2]
        emit(test='two_bends',mesh=mesh,case=case,first=float(stiffness[1,1]),
             flare=float(stiffness[2,2]), clamped_cross=float(stiffness[1,2]),
             effective_tension_from_flare=sigma, first_from_flare_tension=baseline[1,1]+sigma*area[1],
             common_relaxed_bend_eigenvalues=np.linalg.eigvalsh(internal).tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=(26,34),required=True)
    compare(p.parse_args().mesh)
