"""Power reflection from either neck junction under calibrated incident modes.

The source is in a uniform straight extension of the neck. The remote slab
or half-space radiates through the same stretched-coordinate domain as the
bend experiment. Report total power returned into all guided polarizations.
"""
import argparse

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from mouth_open_reservoirs import ReservoirDomain, settings
from mouth_carrier_bend_work import IDENTITY, J, lead, radial_basis, rotate


def truncate(domain, side):
    """Replace the other junction by a flux-controlled uniform launch lead."""
    d=domain
    z=d.z_index
    keep=z<d.neck_end if side=="brane" else z>=d.nfeed
    left,right=d.left,d.right
    active=keep[left]&keep[right]
    ids=np.full(len(z),-1)
    ids[keep]=np.arange(sum(keep))
    d.left,d.right=ids[left[active]],ids[right[active]]
    d.r_index,d.z_index=d.r_index[keep],z[keep]
    d.w_index=d.w_index[keep]
    launch=np.nonzero(d.z_index==(d.neck_end-1 if side=="brane" else d.nfeed))[0]
    shell,_,_=radial_basis(len(launch),d.dr)
    end=np.zeros(sum(keep))
    end[launch]=6*shell/d.dw
    def clip(metric):
        volume,links,inlet,outlet,angular,_=metric
        inlet=inlet[keep] if side=="brane" else end
        outlet=end if side=="brane" else outlet[keep]
        links=links[active]
        adj=sparse.coo_matrix((-links,(d.left,d.right)),shape=(len(inlet),len(inlet))).tocsr()
        adj=adj+adj.T
        op=d.kinetic*(adj+sparse.diags(inlet+outlet-np.asarray(adj.sum(axis=1)).ravel()))
        return volume[keep],links,inlet,outlet,angular[keep],op
    real=clip((d.volume,d.conductance,d.inlet,d.outlet,d.ac[4],d.operator))
    d.volume,d.conductance,d.inlet,d.outlet,_,d.operator=real
    d.ac=clip(d.ac)
    return launch


def incoming(lam,omega,d,delta):
    hop=d.kinetic/d.dw**2
    r=rotate(delta)
    restoring=5*d.eos_k*d.bath_density**4
    central=(2*hop*np.cos(delta)+d.kinetic*lam)*IDENTITY+np.diag([2*restoring,0])+1j*d.hbar*omega*J
    comp=np.block([[np.zeros((2,2)),IDENTITY],[-r.T@r.T,r.T@central/hop]])
    roots,vectors=np.linalg.eig(comp)
    for z,v in zip(roots,vectors[:2].T):
        flux=(omega*d.kinetic/d.dw)*(z*np.vdot(v,r@v)).imag
        if abs(abs(z)-1)<1e-7 and flux < -1e-10:
            return z,v/np.sqrt(-flux)
    return None


def scatter(d,ids,omega,phase,outflow,side):
    n=len(d.psi)
    vol,_,_,_,angular,lap=d.ac
    weight=vol/3
    dens=abs(d.psi)**2
    hc=1.25*d.eos_k
    h=lap/3+sparse.diags(d.kinetic*angular+weight*hc*(dens**4-d.bath_density**4))
    v=8*hc*dens**3*weight
    xy=v*d.psi.real*d.psi.imag
    iz=1j*omega*d.hbar*weight
    A=sparse.bmat([[h+sparse.diags(v*d.psi.real**2),sparse.diags(xy-iz)],
                   [sparse.diags(xy+iz),h+sparse.diags(v*d.psi.imag**2)]],format="csc")
    area,lam,basis=radial_basis(len(ids),d.dr)
    delta=np.arcsin(outflow*d.constituent_mass*d.dw/(3*sum(area)*d.bath_density*d.hbar))
    port=lead(d,ids,omega,phase,outflow)
    t=d.kinetic*np.tile(area,2)/d.dw
    idx=np.r_[ids,ids+n]
    A+=sparse.coo_matrix(((port-np.diag(2*t)).ravel(),(np.repeat(idx,len(idx)),np.tile(idx,len(idx)))),shape=A.shape).tocsc()
    factor=splu(A)
    reflected=[]
    transmissions=[]
    for j,value in enumerate(lam):
        mode=incoming(value,omega,d,delta)
        if mode is None:
            continue
        z,v=mode
        profile=basis[:,j]/np.sqrt(area)
        inc=np.kron(rotate(phase-delta/2)@v,profile)
        ghost=z*np.kron(rotate(phase+delta/2)@v,profile)
        rhs=np.zeros(2*n,complex)
        rhs[idx]=t*(ghost-inc)+port@inc
        u=factor.solve(rhs)
        reflected.append(u[idx]-inc)
        emit=d.cuts(omega,u,round(1.5*d.mouth_cells),round(d.mouth_cells))
        transmissions.append(emit[0] if side=="brane" else emit[1])
    reflected=np.asarray(reflected).T
    W=-(port-port.T.conj())*omega/(2j)
    power=reflected.T.conj()@W@reflected
    vals=np.linalg.eigvalsh(power)
    residual=max(abs(power.diagonal().real+transmissions-1))
    assert min(vals)>-1e-8 and max(vals)<1+1e-7 and residual<1e-7
    print(f"junction={side} omega={omega:.12g} propagating_modes={len(vals)} "
          f"fundamental_return={power[0,0].real:.12g} fundamental_transmission={transmissions[0]:.12g} "
          f"min_coherent_return={min(vals):.12g} max_coherent_return={max(vals):.12g} "
          f"channel_power_residual={residual:.3g}",flush=True)


def run(args):
    _,_,bend,flow,archive=settings()
    p=archive["parameters"]
    R=archive["phase_A"]["q_star"][0]
    rate=0. if args.rest else float(flow["reference_inward_current"])/p["m"]
    print(f"settings cells={args.cells} stub_over_R={args.stub:g} slab_width_over_R=1 DC={rate:.12g}",flush=True)
    for side in ("brane","bulk"):
        d=ReservoirDomain(R,R*args.stub,R,args.cells,2.,2.,3.,p)
        ids=truncate(d,side)
        d.psi,phi,error=d.stationary(rate)
        feed=(d.hbar/d.constituent_mass)*sum(d.inlet*(np.sqrt(d.bath_density)*np.exp(-1j*phi)*d.psi).imag)
        out=(d.hbar/d.constituent_mass)*sum(d.outlet*(np.sqrt(d.bath_density)*d.psi.conj()).imag)
        assert abs(feed-rate)<1e-8 and abs(out-rate)<1e-8
        for tone in ("difference","sum"):
            scatter(d,ids,float(bend[tone]),0. if side=="brane" else phi,
                    rate if side=="brane" else -rate,side)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells",type=int,required=True)
    parser.add_argument("--stub",type=float,default=1.)
    parser.add_argument("--rest",action="store_true")
    run(parser.parse_args())
