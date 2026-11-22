"""Read the tangential wall traction directly on a steadily displaced, replenished flare."""
import argparse

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse.linalg import splu

from mouth_calibrated_flare_coating import read_inputs
from mouth_carrier_lateral_inertia import taper
from mouth_empty_bend_stiffness import INTERVALS
from mouth_middle_tag_grip import cubic_trace
from mouth_replacing_flare_slip import assemble
from neutral_patch_sites import emit


def measure(mesh, refine):
    trace, _ = cubic_trace(mesh)
    patch, _, _ = read_inputs(mesh)
    points,grip,(membrane,between,lining,_,transport) = assemble(trace,patch,refine)
    time = trace.transit(points[0],points[-1],8)
    roots,weights = leggauss(8)
    half = np.diff(points)[:,None]/2
    w = points[:-1,None]+half*(1+roots)
    r, slope = trace.geometry.radius(w)
    arc = np.hypot(1.,slope)
    h = 1-taper(w,*INTERVALS[1])[0]
    jac = 4*np.pi*r*r/3*arc*half*weights
    for tau in (0.,1/32,1/8):
        for registry in (False,True):
            zeta = grip*tau*time
            elastic = membrane+between+registry*lining
            op = elastic+zeta*transport
            slip = splu(op[:-1,:-1].tocsc()).solve(-op[:-1,-1].toarray().ravel())
            fields = np.r_[np.zeros(4),slip].reshape(-1,4)
            a, b = [v[:-1,None]*(1-roots)/2+v[1:,None]*(1+roots)/2 for v in fields[:,:2].T]
            tangent = -registry*grip*np.sum(jac*h*(slope/arc*a+2*b))
            state = np.r_[slip,1.]
            force = (op@state)[-1]
            emit(test='steady_wall_split',tau=tau,wall='lattice' if registry else 'smooth',total=force,
                 tangent_on_facing_sheet=tangent,remaining_normal_constraint=force-tangent,
                 steady_equation_error=float(np.max(abs((op@state)[:-1]))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,default=34)
    p.add_argument('--refine',type=int,default=2)
    args=p.parse_args()
    measure(args.mesh,args.refine)
