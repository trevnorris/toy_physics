"""Ordinary zero-drift carrier current close to the surface from injection to discharge.

For the regular, stationary, axisymmetric continuum with a reflecting side and
positive depth current, the outer cumulative current is a streamline label.
The largest/smallest near-wall label at a list of sections describes union/
intersection at those sections. Neither observable identifies mechanical bonds.
"""
import argparse
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_carrier_lateral_response import prepare, trace
from mouth_flared_receiver import EXIT
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import BETWEEN, emit


ENTRY = 1.175
# A common physical location for all three meshes; near the fine 1%-order cut.
ORDER_CUT = 1.42748


def section(dom, field, w, widths):
    radius = float(dom.radius(w)[0])
    fractions = widths/radius
    edges = np.unique(np.r_[dom.blocks[2][1],fractions])
    q, wt = leggauss(6)
    h = np.diff(edges)[:,None]/2
    v = (edges[:-1,None]+h*(1+q)).ravel()
    r, _, _, rv, _, _, f, _, fw = trace(dom,2,np.full_like(v,w),v,field)
    jw = (f[:,0].conj()*fw[:,0]).imag
    integrals = (4*np.pi*r*r*abs(rv)*jw).reshape(h.size,-1)*h*wt
    outer = np.r_[0.,np.cumsum(integrals.sum(axis=1))]
    return outer[np.searchsorted(edges,fractions)],float(outer[-1]),float(jw.min())


def side_flux(dom, field, z, widths):
    """Outward current on surfaces r=R(w)−width, per unit w."""
    radius, derivative = dom.radius(z[:,None])
    r, _, _, _, _, _, psi, fr, fw = trace(dom,2,np.broadcast_to(z[:,None],(len(z),len(widths))),
                                            widths[None,:]/radius,field)
    current = (psi[...,0].conj()*(fr[...,0]-derivative*fw[...,0])).imag
    return 4*np.pi*r*r*current


def main(args):
    dom,psi=prepare(SimpleNamespace(source=Path(f'_scratch/s11c/mouth_flared_cubic_n{args.mesh}')))
    field=psi[:,None]
    radius0=float(dom.radius(ENTRY)[0])
    widths=np.array([PITCH*np.log1p(1/np.sqrt(2))/BETWEEN[1],PITCH,2*PITCH])
    targets=np.array([ENTRY,ORDER_CUT,2.5,3.,EXIT])
    pg=dom.blocks[2][0]
    breaks=np.unique(np.r_[targets,pg[(pg>ENTRY)&(pg<EXIT)]])
    # Include both one-sided values at native Q2 derivative jumps.
    samples=[]
    for a,b in zip(breaks[:-1],breaks[1:],strict=True):
        nodes=np.linspace(a,b,args.splits+1)
        nodes[0],nodes[-1]=np.nextafter(a,b),np.nextafter(b,a)
        samples.extend(nodes)
    samples=np.array(samples)
    flux,whole,jmin=zip(*(section(dom,field,z,widths) for z in samples),strict=True)
    flux=np.array(flux)
    whole=np.array(whole)
    fractions=flux/whole[:,None]
    emit(test='ordinary_flow_setup',mesh=args.mesh,splits=args.splits,widths=widths.tolist(),
         minimum_sampled_jw=min(jmin),outer_first_radial_cell_at_injection=radius0*dom.blocks[2][1][1],
         outer_first_radial_cell_at_exit=float(dom.radius(EXIT)[0])*dom.blocks[2][1][1],
         nonnegative_native_separation_requires_inner_depth_below=BETWEEN[2]*PITCH,
         tag_law_exponential_log_at_axis=BETWEEN[1]*radius0/PITCH)
    if min(jmin)<=0: raise ValueError('Nonpositive section flux prevents flux rank assignment')
    # j_r−R'j_w checks the weakly imposed impermeable wall and independent
    # radial exchange through the inner edge. End sections are taken in the
    # subinterval. Native trace errors, especially at the outlet corner, remain.
    q,wt=leggauss(12)
    for low,high in zip(targets[:-1],targets[1:],strict=True):
        inside=(samples>low)&(samples<high)
        zi=samples[inside]
        measured=flux[inside]
        contact_fractions=fractions[inside]
        stops=np.r_[low,pg[(pg>low)&(pg<high)],high]
        h=np.diff(stops)[:,None]/2
        z=(stops[:-1,None]+h*(1+q)).ravel()
        radial_flux=(h*wt).ravel()@side_flux(dom,field,z,np.r_[0.,widths])
        predicted=radial_flux[1:]-radial_flux[0]
        for k,width in enumerate(widths):
            emit(test='carrier_annulus',interval=[low,high],width=width,
                 start_fraction=contact_fractions[0,k],end_fraction=contact_fractions[-1,k],
                 current_gain=measured[-1,k]-measured[0,k],
                 inner_to_band=radial_flux[k+1],outer_wall_leak=radial_flux[0],
                 band_continuity_defect=measured[-1,k]-measured[0,k]-predicted[k],
                 min_section_current=whole[inside].min(),max_section_current=whole[inside].max())
    for high in [ORDER_CUT,2.5,3.,EXIT]:
        keep=samples<=high
        depths=samples[keep]
        for k,width in enumerate(widths):
            fraction=fractions[keep,k]
            min_i,max_i=np.argmin(fraction),np.argmax(fraction)
            emit(test='conservative_flux_rank',sampled_interval=[ENTRY,high],width=width,
                 section_count=len(depths),min_fraction=fraction[min_i],min_at=depths[min_i],
                 max_fraction=fraction[max_i],max_at=depths[max_i],
                 min_raw_current=flux[keep,k].min(),max_raw_current=flux[keep,k].max(),
                 start_fraction=fraction[0],end_fraction=fraction[-1])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',type=int,choices=[26,34,48],required=True)
    parser.add_argument('--splits',type=int,default=12)
    main(parser.parse_args())
