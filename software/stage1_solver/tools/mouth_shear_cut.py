"""Replace the unresolved tail by a free or locally absorbing physical cross section."""
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import leggauss
from mouth_order_wave_pulse import linear, quadratic


def truncate(dom,chi,desired):
    old=list(dom.order_blocks)
    pg,vg,ids=old[2]
    i=int(np.argmin(abs(pg-desired)))
    if not 0<i<len(pg)-1:
        raise ValueError('Cut must lie inside the flared conduit')
    keep_ids=ids.copy().reshape(len(pg)-1,len(vg)-1,9)
    keep_ids[i:]=len(chi)
    old[2]=(pg,vg,keep_ids.reshape(-1,9))
    for block in (3,4):
        p,v,origin=old[block]
        old[block]=(p,v,np.full_like(origin,len(chi)))
    return SimpleNamespace(**dict(vars(dom),order_blocks=old)),np.r_[chi,0.],float(pg[i])


def impedance(dom,blocks,pts,chi,cut,floor):
    # Integral of sqrt(rho*chi) against the torsional area measure on the cut.
    ref_p,ref_v,old=dom.order_blocks[2]
    p,v,ids=blocks[2]
    row=np.searchsorted(p,cut)-1
    old_row=np.searchsorted(ref_p,cut)-1
    if max(abs(p[row+1]-cut),abs(ref_p[old_row+1]-cut))>1e-10:
        raise ValueError('Cut not aligned to original and wave grids')
    nodes,weight=leggauss(5)
    center=(v[:-1]+v[1:])/2
    half=np.diff(v)/2
    pos=center[:,None]+half[:,None]*nodes
    old_j=np.searchsorted(ref_v,center)-1
    arg=2*(pos-ref_v[old_j,None])/np.diff(ref_v)[old_j,None]-1.
    old_ids=old.reshape(len(ref_p)-1,len(ref_v)-1,3,3)[old_row,old_j,-1,:]
    order=np.clip(np.einsum('jqi,ji->jq',quadratic(arg)[0],chi[old_ids]),0,1)
    r,_,_,rv,_,_=dom.chart(2,np.full_like(pos,cut),pos)
    w=(8*np.pi/3)*r**4*abs(rv)*weight*half[:,None]*np.sqrt(order*(floor+(1-floor)*order))
    edge=ids.reshape(len(p)-1,len(v)-1,2,2)[row,:,-1,:]
    b=np.zeros(len(pts))
    np.add.at(b,edge.ravel(),(w@linear(nodes)).ravel())
    return b,float(order.min()),float(order.max())
