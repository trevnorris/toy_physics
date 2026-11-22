"""Upper limit on static grip when a bond samples the wall's order within its own length.

The relative slip probe freezes each bond's strength and measures its incremental
Morse force. The comparison includes unfavorable signs of individual curvatures.
"""
import argparse

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_calibrated_flare_coating import read_inputs
from mouth_flare_grip_memory import PassingGrip
from mouth_flare_support_bounds import minimum
from mouth_flared_receiver import emit
from mouth_middle_tag_grip import cubic_trace
from mouth_replacing_flare_slip import assemble
from mouth_segment_contact_hold import PITCH
from neutral_patch_periodic_compression import Bilayer
from neutral_patch_sites import GAP, morse


def run(mesh,refine):
    trace,intake=cubic_trace(mesh)
    patch,base,area=read_inputs(mesh)
    points,grip,(membrane,_,wall,_,_)=assemble(trace,patch,refine)
    process=PassingGrip(trace,points,grip,8)
    edges=trace.edges
    left,mid,right=trace.chi.T
    endpoint_slope=np.r_[(-3*left+4*mid-right)/np.diff(edges),
                         (left-4*mid+3*right)/np.diff(edges)]
    jumps=left[1:]-right[:-1]
    if max(endpoint_slope)>0 or np.any(jumps>0):
        raise ArithmeticError('Use interval extrema for a nonmonotone order trace')
    cell=Bilayer(1,1,True,patch['shell'])
    vectors=patch['pitch']*cell.vectors[cell.cross]
    length=np.sqrt(np.einsum('ij,ij->i',vectors,vectors)+GAP**2)
    depth,alpha,rest=[p[cell.cross] for p in cell.params]
    _,first=morse(length,depth,alpha,rest)
    decay=np.exp(-alpha*(length-rest))
    second=2*depth*alpha**2*decay*(2*decay-1)
    perp=first/length
    parallel=perp+(second-perp)*np.sum(vectors**2,axis=1)/length**2
    # In a complete sixfold orbit the positive and negative spectral parts each
    # sum to half their summed traces times the identity in the sheet plane.
    positive=(np.maximum(parallel,0)+np.maximum(perp,0))/2
    negative=(np.minimum(parallel,0)+np.minimum(perp,0))/2
    normalizer=(positive+negative).sum()
    scale=patch['shared_scale']*2*GAP/PITCH
    np.testing.assert_allclose(normalizer*scale,grip,rtol=2.e-12)
    roots,_=leggauss(6)
    z=(points[:-1,None]+np.diff(points)[:,None]*(1+roots)/2).ravel()
    location=trace.crossing(.01)
    samples=np.r_[z,location]
    reach=PITCH*length
    highest=trace.sample(trace.chi,np.maximum(edges[0],samples[:,None]-reach))[0]
    lowest=trace.sample(trace.chi,np.minimum(edges[-1],samples[:,None]+reach))[0]
    highest,lowest=np.clip(highest,0,1),np.clip(lowest,0,1)
    upper=(highest@positive+lowest@negative)/normalizer
    # Independent two-dimensional construction of the matrix envelope. This also
    # checks isotropy if birth order or floating bond lengths split nominal shells.
    normalized=np.divide(vectors,np.linalg.norm(vectors,axis=1)[:,None],out=np.zeros_like(vectors),
                         where=np.linalg.norm(vectors,axis=1)[:,None]>0)
    orient=normalized[:,:,None]*normalized[:,None,:]
    def spectral(fn):
        return fn(perp)[:,None,None]*np.eye(2)+(fn(parallel)-fn(perp))[:,None,None]*orient
    pos=spectral(lambda x:np.maximum(x,0))
    neg=spectral(lambda x:np.minimum(x,0))
    upper_matrix=(highest@pos.reshape(-1,4)+lowest@neg.reshape(-1,4)).reshape(-1,2,2)/normalizer
    discrepancy=np.max(abs(upper_matrix-upper[:,None,None]*np.eye(2)))
    if discrepancy>1.e-10:
        raise ArithmeticError('Scalar upper comparison lost the material orientation')
    profile=np.zeros_like(process.order)
    # Use the largest eigenvalue at each point to retain the full matrix envelope.
    profile[process.quad_idx]=np.linalg.eigvalsh(upper_matrix[:-1])[:,1]
    envelope=membrane+process.coupling(profile)
    all_slips=np.arange(membrane.shape[0]-1)
    emit(test='reach',mesh=mesh,refine=refine,intake=intake,bond_count=len(reach),
         largest_bond_length=float(reach.max()),gauge_depth=location,gauge_order=trace.sample(trace.chi,location)[0],
         upper_grip_fraction_at_gauge=upper[-1],signed_negative_curvature_over_neutral=negative.sum()/normalizer,
         order_derivative_max=endpoint_slope.max(),max_order_trace_jump=np.max(abs(jumps)),
         matrix_isotropy_error=discrepancy)
    for name,free,k in [('locked',all_slips[all_slips%4>=2],envelope),('neutral',all_slips,envelope+wall)]:
        _,force,energy,residual=minimum(k,free)
        total=force+base[2,2]
        targets=base[2,2]+np.array([.02,.04])*area[2]
        emit(test='bond_upper',wall=name,total_max=total,candidate_total=targets.tolist(),
             target_minus_upper=(targets-total).tolist(),coating=force,energy=energy,residual=residual)
        np.testing.assert_allclose(energy,force,rtol=2.e-8)
        if residual>1.e-8: raise ArithmeticError('Unbalanced upper comparison')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=(26,34),required=True)
    p.add_argument('--refine',type=int,choices=(1,2,4),required=True)
    args=p.parse_args()
    run(args.mesh,args.refine)
