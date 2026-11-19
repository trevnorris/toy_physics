"""Angular carrier radiation from a neck joining a brane slab to a 4D half-space.

The nonlinear stationary solve sets the intake at the far radial edge of
the slab. Both untouched slab faces and the bulk plane outside the aperture
reflect the scalar carrier. Outgoing waves are continued in complex radial
and axial coordinates outside the monitoring cuts.
"""
import argparse

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from mouth_carrier_bend_budget import settings
from mouth_carrier_bend_work import OpenSleeve, wall_source
from mouth_carrier_bend_spectrum import dipole_widths


class ReservoirDomain(OpenSleeve):
    def __init__(self,radius,half_length,width,cells,clear,layer,stretch,params):
        self.mouth_cells=self.sleeve_cells=cells
        self.dr=radius/cells
        self.dw=half_length/cells
        self.radius=radius
        self.length=2*half_length
        self.stretch=stretch
        self.layer=layer*radius
        self.clear=clear*radius
        self.hbar=params["hbar"]
        self.constituent_mass=params["m"]
        self.eos_k=params["K_eos"]
        self.bath_density=params["rho0"]/self.constituent_mass
        self.kinetic=self.hbar**2/(2*self.constituent_mass)
        self.nr=round((clear+layer)*cells)
        self.nfeed=round(width/self.dr)
        nb=round((clear+layer)*cells)
        self.re=self.dr*np.arange(self.nr+1)
        self.ze=np.r_[np.linspace(-width,0,self.nfeed+1),self.dw*np.arange(1,cells+1),
                      half_length+self.dr*np.arange(1,nb+1)]
        self.neck_end=self.nfeed+cells
        self.nw=len(self.ze)-1
        rr,zz=np.indices((self.nr,self.nw))
        mask=(zz<self.nfeed)|(zz>=self.neck_end)|(rr<cells)
        ids=np.full(mask.shape,-1)
        ids[mask]=np.arange(mask.sum())
        self.r_index,self.z_index=rr[mask],zz[mask]
        # The common sidewall formula uses neck-local axial cell indices.
        self.w_index=np.where(self.z_index<self.nfeed,self.nw+self.z_index,self.z_index-self.nfeed)
        self.lr,self.lz=np.nonzero(mask[:-1,:]&mask[1:,:])
        self.ar,self.az=np.nonzero(mask[:,:-1]&mask[:,1:])
        self.left=np.r_[ids[self.lr,self.lz],ids[self.ar,self.az]]
        self.right=np.r_[ids[self.lr+1,self.lz],ids[self.ar,self.az+1]]
        self.volume,self.conductance,self.inlet,self.outlet,_,self.operator=self.metric(self.re,self.ze)
        def stretch_axis(x,start):
            u=np.maximum(0.,x-start)
            return x+1j*self.stretch*u**3/(3*self.layer**2)
        self.ac=self.metric(stretch_axis(self.re,self.clear),
                            stretch_axis(self.ze,half_length+self.clear))

    def metric(self,r,z):
        dr,dz=np.diff(r),np.diff(z)
        shell=4*np.pi*np.diff(r**3)/3
        ir,iz=self.r_index,self.z_index
        volume=shell[ir]*dz[iz]
        links=np.r_[4*np.pi*r[self.lr+1]**2*dz[self.lz]/((dr[self.lr]+dr[self.lr+1])/2),
                    shell[self.ar]/((dz[self.az]+dz[self.az+1])/2)]
        radial_end=np.where(ir==self.nr-1,8*np.pi*r[-1]**2*dz[iz]/dr[-1],0.)
        inlet=np.where(iz<self.nfeed,radial_end,0.)
        outlet=np.where(iz>=self.neck_end,radial_end,0.)+np.where(iz==self.nw-1,2*shell[ir]/dz[-1],0.)
        adj=sparse.coo_matrix((-links,(self.left,self.right)),shape=(len(ir),len(ir))).tocsr()
        adj=adj+adj.T
        operator=self.kinetic*(adj+sparse.diags(inlet+outlet-np.asarray(adj.sum(axis=1)).ravel()))
        angular=8*np.pi*dr[ir]*dz[iz]/3
        return volume,links,inlet,outlet,angular,operator

    def stationary(self,rate):
        # A unit feed phase solves the linearized current problem. Starting
        # on that phase profile keeps the nonlinear iteration well scaled even
        # when the feed is far from the aperture.
        bath=np.sqrt(self.bath_density)
        hc=1.25*self.eos_k
        scale=self.hbar/self.constituent_mass
        if rate==0:
            return np.full(len(self.volume),bath,dtype=complex),0.,0.
        shape=spsolve(self.operator,self.kinetic*self.inlet)
        phi=rate/(scale*bath**2*np.dot(self.inlet,shape-1))
        psi=bath*np.exp(1j*phi*shape)
        def evaluate(field,angle):
            feed=bath*np.exp(1j*angle)
            n=abs(field)**2
            local=self.volume*hc*(n**4-self.bath_density**4)
            edge=self.conductance*(field[self.left]-field[self.right])
            residual=self.inlet*(field-feed)+self.outlet*(field-bath)
            np.add.at(residual,self.left,edge)
            np.add.at(residual,self.right,-edge)
            residual=self.kinetic*residual+local*field
            flux=scale*np.dot(self.inlet,(feed.conjugate()*field).imag)-rate
            return np.r_[residual.real,residual.imag,flux],local,feed
        for _ in range(20):
            residual,local,feed=evaluate(psi,phi)
            err=max(abs(residual))
            if err<3e-11:
                # This cut stays inside the neck even in isolated-junction tests.
                mid=self.nfeed+self.sleeve_cells//2
                internal=(self.z_index[self.left]==mid-1)&(self.z_index[self.right]==mid)
                through=scale*np.dot(self.conductance[internal],
                    (psi[self.left[internal]].conj()*psi[self.right[internal]]).imag)
                exterior=scale*np.dot(self.outlet,(psi.conjugate()*bath).imag)
                if max(abs(through-rate),abs(exterior-rate),abs(residual[-1]))<1e-9:
                    return psi,phi,err
            x,y=psi.real,psi.imag
            slope=self.volume*8*hc*abs(psi)**6
            base=self.operator+sparse.diags(local)
            xy=sparse.diags(slope*x*y)
            phase_col=-1j*self.kinetic*self.inlet*feed
            grad=-scale*np.dot(self.inlet,feed.real*x+feed.imag*y)
            jac=sparse.bmat([[base+sparse.diags(slope*x*x),xy,sparse.csr_matrix(phase_col.real[:,None])],
                             [xy,base+sparse.diags(slope*y*y),sparse.csr_matrix(phase_col.imag[:,None])],
                             [sparse.csr_matrix(-scale*feed.imag*self.inlet),sparse.csr_matrix(scale*feed.real*self.inlet),
                              sparse.csr_matrix([[grad]])]],format="csc")
            delta=spsolve(jac,-residual)
            a,b=np.split(delta[:-1],2)
            for power in range(15):
                factor=.5**power
                candidate=psi+factor*(a+1j*b)
                angle=phi+factor*delta[-1]
                if max(abs(evaluate(candidate,angle)[0])) < err:
                    psi,phi=candidate,angle
                    break
            else:
                raise ArithmeticError(f"Stationary residual stopped decreasing at {err:g}")
        raise ArithmeticError(f"Stationary residual {err:g}")

    def set_current(self,rate):
        self.psi,phase,res=self.stationary(rate)
        s=self.hbar/self.constituent_mass
        self.rates=(s*sum(self.inlet*(np.sqrt(self.bath_density)*np.exp(-1j*phase)*self.psi).imag),
                    s*sum(self.outlet*(self.psi.conj()*np.sqrt(self.bath_density)).imag))
        j=s*self.conductance*(self.psi[self.left].conj()*self.psi[self.right]).imag
        neck=sum(j[(self.z_index[self.left]==self.nfeed-1)&(self.z_index[self.right]==self.nfeed)])
        if max(abs(v-rate) for v in (*self.rates,neck))>1e-8:
            raise ArithmeticError("Background flow is not conserved")
        print(f"DC rate={rate:.12g} inlet={self.rates[0]:.12g} neck={neck:.12g} bulk={self.rates[1]:.12g} residual={res:.3g}",flush=True)
        self.spatial,self.temporal=wall_source(self,self.psi,self.length,2)

    def cuts(self,omega,u,cut_r,cut_w):
        # Integrate signed power out of a box entirely within real coordinates.
        # Feed and bulk are disconnected across the reflecting shoulders outside R.
        r=self.r_index
        z=self.z_index
        l,h=self.left,self.right
        inside=(r<cut_r)&(z<self.neck_end+cut_w)
        cross=inside[l]&~inside[h]
        a,b=np.split(u,2)
        link=omega*self.kinetic*self.conductance/3*(a[l].conj()*a[h]+b[l].conj()*b[h]).imag
        return np.array([sum(link[cross&(z[l]<self.nfeed)]),sum(link[cross&(z[l]>=self.neck_end)])])

    def response(self,omega):
        volume,_,_,_,ang,lap=self.ac
        weight=volume/3
        p=self.psi
        n=abs(p)**2
        hc=1.25*self.eos_k
        h=lap/3+sparse.diags(self.kinetic*ang+hc*(n**4-self.bath_density**4)*weight)
        coupling=8*hc*n**3*weight
        xy=coupling*p.real*p.imag
        symp=1j*omega*self.hbar*weight
        a=sparse.bmat([[h+sparse.diags(coupling*p.real**2),sparse.diags(xy-symp)],
                       [sparse.diags(xy+symp),h+sparse.diags(coupling*p.imag**2)]],format="csc")
        src=self.spatial+1j*omega*self.temporal
        u=spsolve(a,src/2)
        reaction=np.vdot(src,u)
        power=.5*omega*reaction.imag
        # Second energy ledger on two distinct interior monitoring surfaces.
        cuts=[self.cuts(omega,u,round(f*self.mouth_cells),round((f-.5)*self.mouth_cells))
              for f in (1.25,self.clear/self.radius-.125)]
        if not np.allclose(cuts[0],cuts[1],atol=4e-8,rtol=1e-8) or abs(sum(cuts[0])-power)>4e-8*max(1,abs(power)):
            raise ArithmeticError(f"Wall and remote cut powers disagree: {power}, {cuts}")
        return reaction,power,*cuts[0]


def calculate(args):
    point,energy,bend,current,archive=settings()
    rate=float(current["reference_inward_current"])/archive["parameters"]["m"]
    r=archive["phase_A"]["q_star"][0]
    domain=ReservoirDomain(r,float(point["L"])/2,args.width*r,args.cells,args.clear,args.layer,args.stretch,archive["parameters"])
    domain.set_current(0. if args.rest else rate)
    print(f"geometry cells={args.cells} slab_width_over_R={args.width:g} PML_start_over_R={args.clear:g} "
          f"PML_width_over_R={args.layer:g} stretch={args.stretch:g} active_cells={len(domain.volume)}",flush=True)
    if args.spectrum:
        static=domain.response(0.)[0].real
        offsets=np.array([-.016,-.012,-.009,-.0075,-.006,-.004,-.002,0.,.004])
    else:
        static=0.
        offsets=np.array([0.])
    forces=[]
    for tone in ("difference","sum"):
        row=[]
        for offset in offsets:
            omega=float(bend[tone])+offset
            force,power,feed,bulk=domain.response(omega)
            if min(feed,bulk)<=0:
                raise ArithmeticError("Nonpositive emitted power")
            print(f"tone={tone} offset={offset:.8g} omega={omega:.12g} power_per_b2={power:.12g} "
                  f"feed_per_b2={feed:.12g} bulk_per_b2={bulk:.12g} feed_fraction={feed/power:.10g} "
                  f"effective_C={2*power/omega**2:.12g} cycle_b1em5={2*np.pi*1e-10*power/omega:.12g} "
                  f"induced_real={force.real:.12g}",flush=True)
            row.append(force-static)
        forces.append(np.array(row))
    if args.spectrum:
        ws,wd=float(point["omega_s"]),float(point["omega_d"])
        g=float(bend["overlap"])*np.sqrt(2*archive["parameters"]["I_wave"]/ws)
        gap=abs(float(energy["U_over_action"]))
        center,width,asym,count,err=dipole_widths(offsets,forces,ws,wd,g,gap,1)
        print(f"spectrum arms=1 static_induced={static:.12g} shift={center:.12g} FWHM={width:.12g} "
              f"width_over_gap={width/gap:.12g} asymmetry={asym:.8g} frequency_sampling_change={err:.4g} peaks={count}",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells",required=True,type=int)
    parser.add_argument("--width",default=1.,type=float,help="slab thickness in mouth radii")
    parser.add_argument("--clear",default=2.,type=float,help="radius and bulk depth where absorption starts, in R")
    parser.add_argument("--layer",default=2.,type=float,help="absorbing layer thickness in R")
    parser.add_argument("--stretch",default=3.,type=float)
    parser.add_argument("--rest",action="store_true")
    parser.add_argument("--spectrum",action="store_true")
    calculate(parser.parse_args())
