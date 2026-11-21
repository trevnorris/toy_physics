"""Hold the first bend of the 0PN cubic mouth with moving microscopic partners.

Each family has wall tips tethered to fixed material stations and brane tips
carried past them by the measured flow. Incremental slip couples the tips
through short overlapping bindings. Both tips can deflect against their own
backbone. All wall stations share the four-dimensional exterior shear field.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import linalg, sparse

from mouth_bend_renewal import volume_cost
from mouth_bend_static_shear import ANGULAR
from mouth_carrier_lateral_inertia import taper
from mouth_empty_bend_stiffness import INTERVALS
from mouth_middle_tag_grip import cubic_trace
from mouth_order_wave_pulse import ROOT, emit
from mouth_two_depth_load import DEPTHS
from neutral_patch_conversion import NEUTRAL
from neutral_patch_sites import BETWEEN, COUNT, GAP, morse
from neutral_patch_turning_neighbor import AREA_PER_SITE, ENTRY, LAYERS, sheet_modulus

PITCH = .01
RECORD = Path('research/pde_ledger_v3/_measurements/neutral_patch_conversion.json')
AMPLITUDE = 1e-5


def coupon():
    ids = np.array(json.loads(RECORD.read_text())['independent_id_extraction']['laws']['central']['strict_surviving_ids'])
    coords = np.array(json.loads(NEUTRAL.read_text())['saved_measurements']['frames1321.json']['contents']['central']['rested_xy'])
    offsets = coords[ids[:,1]]-coords[ids[:,0]]
    lateral2 = np.sum(offsets**2,axis=1)
    length = np.sqrt(GAP**2+lateral2)
    depth,alpha,rest = BETWEEN
    _,tension = morse(length,*BETWEEN)
    exp = np.exp(-alpha*(length-rest))
    second = 2*depth*alpha**2*exp*(2*exp-1)
    # Isotropic directions of the unperturbed offset in the three wall tangents.
    tangent = tension/length+(second-tension/length)*lateral2/(3*length**2)
    # Independent bond-energy curvature check for the orientation average.
    h=1e-4
    trace=np.zeros(len(ids))
    for proj in (np.sqrt(lateral2),np.zeros(len(ids)),np.zeros(len(ids))):
        e0=morse(length,*BETWEEN)[0]
        ep=morse(np.sqrt(length**2+2*h*proj+h*h),*BETWEEN)[0]
        em=morse(np.sqrt(length**2-2*h*proj+h*h),*BETWEEN)[0]
        trace+=(ep+em-2*e0)/(3*h*h)
    np.testing.assert_allclose(trace,tangent,rtol=2e-6)
    if min(tangent)<=0 or len(ids)!=9:
        raise ValueError('Neutral calibrator did not yield stable registered pairs')
    per_site_shear=sheet_modulus(ENTRY)*AREA_PER_SITE/LAYERS
    return ids,tangent,per_site_shear


class Segment:
    def __init__(self,args):
        trace,intake=cubic_trace(args.mesh)
        geo=trace.geometry
        path=ROOT/f'mouth_free_order_cubic_n{args.mesh}_d04/order_seed1.55.npz'
        with np.load(path) as data:
            edges,live,bulk=volume_cost(geo,data['chi'],args.sub,5,1e-8)
        self.transit=trace.transit(edges[0],2.,8)
        roots,weights=leggauss(5)
        waist,middle=edges[0],DEPTHS[1]
        last=edges[np.flatnonzero(live)[-1]+1]
        # Split at material edges, readout, and ends of the held shape.
        cuts=np.unique(np.r_[edges[edges<=last],middle,INTERVALS[0]])
        left=np.repeat(cuts[:-1],len(roots))
        length=np.repeat(np.diff(cuts),len(roots))
        w=left+length*np.tile((1+roots)/2,len(cuts)-1)
        dw=length*np.tile(weights/2,len(cuts)-1)
        cell=np.searchsorted(edges,w)-1
        use=live[cell]|live[cell+1]
        w,dw,cell=w[use],dw[use],cell[use]
        width=np.diff(edges)[cell]
        f=(w-edges[cell])/width
        a,slope=geo.radius(w)
        arc=np.hypot(1,slope)
        chi=np.clip(trace.sample(trace.chi,w)[0],0,1)
        if np.min(chi)<=0: raise ValueError('Non-positive order on active trace')
        row=np.arange(len(w))
        def basis(a0,a1):
            return sparse.coo_matrix((np.r_[a0,a1],(np.r_[row,row],np.r_[cell,cell+1])),
                                     shape=(len(w),len(edges))).tocsc()[:,live].tocsr()
        shape=basis(1-f,f)
        mer=sparse.diags(1/a**2)@shape
        ang=sparse.diags(1/(2*a*arc))@basis(-1/width,1/width)
        self.J=sparse.vstack((mer,ang),format='csr')
        area=a*a*arc*dw
        W=np.r_[area,2*area]
        self.scatter=(self.J.T@sparse.diags(W)).tocsr()
        load=shape.T@(dw*a*taper(w,*INTERVALS[0])[1])
        ids,k,per_site=coupon()
        self.families=len(k)
        self.points=len(w)
        self.grip=args.brane_grip
        self.density=args.coverage*(len(ids)/COUNT)/(PITCH*per_site)
        self.k=chi[:,None]*k[None,:]
        self.mean_k=chi*np.mean(k)
        self.share=k/k.sum()
        inv=np.tile(1/(self.density*self.mean_k),2)
        S=bulk+(self.scatter@sparse.diags(inv)@self.J).toarray()
        norm=1/np.sqrt(S.diagonal())
        factor=linalg.cho_factor(norm[:,None]*S*norm)
        self.solve=lambda rhs: norm*linalg.cho_solve(factor,norm*rhs,check_finite=False)
        self.force=load*args.amplitude
        self.amplitude=args.amplitude
        i=np.searchsorted(edges,middle)-1
        part=(middle-edges[i])/(edges[i+1]-edges[i])
        cut=np.zeros(len(edges));cut[i:i+2]=[1-part,part]
        self.cut=-ANGULAR*geo.radius(np.array([middle]))[0][0]*cut[live]
        self.band=w<middle
        self.xarea=-ANGULAR*np.array([area*slope/arc,2*area])
        r0,s0=geo.radius(np.array([waist]))
        self.period=PITCH/(trace.speed(waist)*np.hypot(1,s0[0]))
        self.age=np.array([trace.transit(waist,wi,6) for wi in w])/self.period
        self.multiplicity=area/(r0[0]**2*PITCH)
        self.phase=(np.arange(len(k))/len(k)+args.phase)%1
        self.offset=int(np.ceil(max(self.age)))+3
        self.parcels=self.offset+int(np.ceil(self.transit/self.period))+5
        self.labels=np.arange(len(k))*self.parcels
        emit(test='setup',mesh=args.mesh,sub=args.sub,dt=args.dt,coverage=args.coverage,
             brane_grip=self.grip,wall_grip=1.,tip_drags=[1.,1.],pitch=PITCH,phase=args.phase,
             intake=intake,transit=self.transit,inlet_duration_per_pitch=self.period,
             families=ids.tolist(),neutral_tangents=k.tolist(),shear_energy_curvature_per_bulk_site=per_site,
             contact_force_scale=self.density,ordered_instantaneous_grip=self.density*np.mean(k),
             amplitude=args.amplitude,read_cut=middle,active_end=last,wall_samples=len(w))

    def rates(self,t,x,y,observe=False):
        coordinate=t/self.period-self.age[:,None]+self.phase[None,:]
        lower=np.floor(coordinate).astype(int)
        weight=coordinate-lower
        indices=lower+self.offset+self.labels
        if indices.min()<0 or indices.max()+1>=y.shape[1]:
            raise ArithmeticError('Material labels outside allocated stream')
        low=y[:,indices]; high=y[:,indices+1]
        blend=low+(high-low)*weight
        relaxed=np.sum((x+blend)*self.share,axis=2)
        stress=self.solve(self.force-self.scatter@relaxed.ravel())
        traction=(self.J@stress).reshape(2,self.points)
        slip=relaxed+traction/(self.density*self.mean_k)
        extension=slip[:,:,None]-x
        lower_force=self.k*(extension-low)*(1-weight)
        upper_force=self.k*(extension-high)*weight
        tip_force=lower_force+upper_force
        dx=tip_force-x
        dy=-self.grip*y
        for axis in range(2):
            for loc,force in ((indices,lower_force[axis]),(indices+1,upper_force[axis])):
                dy[axis]+=np.bincount(loc.ravel(),weights=(force*self.multiplicity[:,None]).ravel(),
                                      minlength=y.shape[1])
        if observe:
            # Actual wall support receives the tether and the wall-side drag.
            support=self.density*(x+dx).mean(axis=2)
            total=(self.xarea*support)[:,self.band].sum()/self.amplitude
            flux=self.cut@stress/self.amplitude
            share_error=np.max(abs(support-traction))/abs(self.amplitude)
            # Reconstruct action and reaction once for each tangent coefficient.
            sites=np.sum(tip_force*self.multiplicity[None,:,None],axis=1)
            packets=np.sum((dy+self.grip*y).reshape(2,self.families,self.parcels),axis=2)
            balance=np.max(abs(sites-packets))/abs(self.amplitude)
            at_entry=(lower-self.phase)>0
            fresh=(1-weight)*at_entry+weight*((lower+1-self.phase)>0)
            emit(test='reading',time=t,force_per_displacement=total,from_cut=flux,
                 fresh_fraction_min_in_band=float(np.min(fresh[self.band])),
                 force_assembly_error=share_error,force_exchange_error=balance)
            if not np.isfinite(total) or max(abs(total-flux),share_error,balance)>2e-7:
                raise ArithmeticError('Independent support, cut or contact force balance failed')
            return total
        return dx,dy


def run(args):
    s=Segment(args)
    x=np.zeros((2,s.points,s.families)); y=np.zeros((2,s.parcels*s.families))
    t=0.
    samples={0.:s.rates(0.,x,y,True)}
    for end in (.0625,1.,8.,16.,28.3125,s.transit):
        count=int(np.ceil((end-t)/args.dt-1e-12))
        dt=(end-t)/count
        start=t
        for i in range(count):
            t=start+i*dt
            dx,dy=s.rates(t,x,y)
            vx,vy=s.rates(t+dt/2,x+dx*(dt/2),y+dy*(dt/2))
            x+=dt*vx; y+=dt*vy
        t=end
        samples[end]=s.rates(end,x,y,True)
    emit(test='retention',reference_times=[.0625,28.3125],
         matched_ratio=samples[28.3125]/samples[.0625],
         transit_ratio=samples[s.transit]/samples[.0625],
         ratio_from_instant=samples[s.transit]/samples[0.])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=(26,34),default=34)
    p.add_argument('--sub',type=int,choices=(2,3,4),default=3)
    p.add_argument('--dt',type=float,default=1/128)
    p.add_argument('--brane-grip',type=float,default=1.)
    p.add_argument('--coverage',type=float,default=1.)
    p.add_argument('--phase',type=float,default=0.)
    p.add_argument('--amplitude',type=float,default=AMPLITUDE)
    args=p.parse_args()
    if args.coverage<=0 or args.dt<=0 or args.amplitude==0 or args.brane_grip<0:
        p.error('Invalid time step, coverage, anchoring or hold')
    run(args)
