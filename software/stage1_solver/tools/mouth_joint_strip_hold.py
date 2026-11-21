"""Weak hold of the first bend with an open, slipping strip coupled to the wall.

Twelve-site sections arrive from an open reservoir. Their longitudinal Morse
stretch and the transverse neighbor shear are separate terms in the material
energy. In this infinitesimal test the axial response is even in the bend
displacement, so zero and driven copies share the axial trajectory to first order.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.linalg import solve_banded

from mouth_connected_passage import WAIST, MIDDLE, arc_speed, bond_derivatives, integrals
from mouth_middle_tag_grip import cubic_trace
from mouth_order_wave_pulse import digest, emit
from mouth_segment_contact_hold import PITCH, Segment
from neutral_patch_sites import WITHIN, morse
from neutral_patch_turning_neighbor import AREA_PER_SITE, ENTRY, LAYERS, shear_coupon, sheet_modulus


class TravelingStrip:
    def __init__(self, trace, args):
        self.end = 2.
        knots = np.unique(np.r_[WAIST, MIDDLE, self.end,
                                trace.edges[(trace.edges > WAIST) & (trace.edges < self.end)]])
        w = np.r_[np.concatenate([np.linspace(a,b,48,endpoint=False)
                                  for a,b in zip(knots[:-1],knots[1:])]), self.end]
        self.depths = w
        self.arc = np.r_[0.,integrals(lambda v: np.hypot(1.,trace.geometry.radius(v)[1]),w)]/PITCH
        self.carrier_time = np.r_[0.,integrals(lambda v: 1/trace.speed(v),w)]
        self.gates = self.arc[np.searchsorted(w,[WAIST,MIDDLE,self.end])]
        self.speed = arc_speed(trace,w)/PITCH
        self.order = np.clip(trace.sample(trace.chi,w)[0],0.,1.)
        self.inlet, self.outlet = -args.lead, self.gates[-1]+args.exit_lead
        self.pos = np.arange(self.inlet,self.outlet,1.)
        self.id0 = 0
        self.hoop_at_gate = trace.geometry.radius(np.array([WAIST,MIDDLE,self.end]))[0]/trace.geometry.radius(WAIST)[0]

    def advance(self, dt):
        pos = self.pos
        tension, k = bond_derivatives(np.diff(pos))
        velocity = np.interp(pos,self.arc,self.speed)
        velocity[:-1] += tension
        velocity[1:] -= tension
        mat = np.zeros((3,len(pos)))
        mat[1] = 1
        mat[1,:-1] += dt*k
        mat[1,1:] += dt*k
        mat[0,1:] = mat[2,:-1] = -dt*k
        new = pos + dt*solve_banded((1,1),mat,velocity,check_finite=False)
        if not np.all(np.diff(new)>0) or not np.isfinite(new).all():
            raise ArithmeticError('Axial strip lost ordering')
        return new

    def exchange(self, new):
        n = np.searchsorted(new,self.outlet)
        add = new[0] >= self.inlet+1.
        if new[0] >= self.inlet+2. or n < 3:
            raise ArithmeticError('Unresolved material boundary event')
        self.pos = np.r_[new[0]-1.,new[:n]] if add else new[:n]
        self.id0 -= int(add)
        return n,add


class HeldWall:
    def __init__(self,args):
        self.wall = Segment(args)
        trace, intake = cubic_trace(args.mesh)
        self.strip = TravelingStrip(trace,args)
        s = self.strip
        np.testing.assert_allclose(s.carrier_time[-1],self.wall.transit,rtol=2e-9)
        # Recover the actual quadrature sites of the archived wall from their
        # emitted travel coordinates. Subsequent contact weights use live positions.
        self.stations = PchipInterpolator(s.carrier_time,s.arc)(self.wall.age*self.wall.period)
        if not np.all((self.stations<s.gates[1]) == self.wall.band):
            raise ArithmeticError('Wall support and material gate use different bands')
        with np.load(args.coupon) as data:
            self.rates = data['legacy_rates']
            self.port = data['legacy_B'][0]
            self.tags = data['labels'][:,1]
        prior = json.loads(Path('research/pde_ledger_v3/_measurements/mouth_mobile_neighbor_hold.json').read_text())
        if digest(args.coupon)!=prior['coupon_archive']['sha256']:
            raise ValueError('Wrong archived material response')
        if self.rates[0]!=0 or np.any(self.rates[1:]<=0):
            raise ValueError('Invalid isolated section modes')
        np.testing.assert_allclose(self.port.T@self.port,np.eye(self.wall.families),atol=2e-12)
        check = shear_coupon()
        self.shear = sheet_modulus(ENTRY)*AREA_PER_SITE/LAYERS if args.links else 0.
        emit(test='joint_setup',intake=intake,arc_gates=s.gates.tolist(),
             carrier_clocks=s.carrier_time[np.searchsorted(s.arc,s.gates)].tolist(),
             model='order-weighted reference shear per material bond, smooth axial sliding',
             shear_at_reference=self.shear,independent_shear_coupon=check['curvature_readings'][-1],
             input_sha256=digest(args.coupon),tags=self.tags.tolist(),lead=[args.lead,args.exit_lead],
             dt=args.dt,mesh=args.mesh,sub=args.sub,amplitude=args.amplitude)

    def rhs(self, pos, x, z, diagnostics=False):
        w = self.wall
        hi = np.searchsorted(pos,self.stations)
        lo = hi-1
        if lo.min()<0 or hi.max()>=len(pos):
            raise ArithmeticError('Material does not cover the active wall')
        hi_part = (self.stations-pos[lo])/(pos[hi]-pos[lo])
        tips = self.port.T @ z
        low = tips[:,:,lo].transpose(0,2,1)
        high = tips[:,:,hi].transpose(0,2,1)
        mean_tip = low+(high-low)*hi_part[None,:,None]
        relaxed = (x+mean_tip)@w.share
        stress = w.solve(w.force-w.scatter@relaxed.ravel())
        traction = (w.J@stress).reshape(2,w.points)
        slip = relaxed+traction/(w.density*w.mean_k)
        excess = slip[:,:,None]-x
        lower = w.k*(excess-low)*(1-hi_part[:,None])
        upper = w.k*(excess-high)*hi_part[:,None]
        applied = lower+upper
        dx = applied-x
        n = len(pos)
        offsets = n*np.arange(w.families)
        ext = np.empty((2,w.families,n))
        for axis in range(2):
            ext[axis] = (np.bincount((lo[:,None]+offsets).ravel(),weights=(lower[axis]*w.multiplicity[:,None]).ravel(),minlength=w.families*n)
                       + np.bincount((hi[:,None]+offsets).ravel(),weights=(upper[axis]*w.multiplicity[:,None]).ravel(),minlength=w.families*n)).reshape(w.families,n)
        order = np.interp((pos[1:]+pos[:-1])/2,self.strip.arc,self.strip.order)
        force = self.port @ ext
        difference = np.diff(z,axis=-1)
        tension = self.shear/ENTRY*order*difference
        force[:,:,:-1] += tension
        force[:,:,1:] -= tension
        net = float(w.cut@stress)
        support = float(np.sum((w.xarea*(w.density*(x+dx).mean(axis=2)))[:,w.band]))
        error = abs(support-net)
        if error>2e-8*max(abs(w.amplitude),1e-10) or not np.isfinite(net):
            raise ArithmeticError('Signed wall reaction disagrees with cut stress')
        if diagnostics:
            exchange = np.max(abs(ext.sum(axis=2)-(applied*w.multiplicity[None,:,None]).sum(axis=1)))
            # Neighbor work cancels the rate of stored shear at fixed material
            # geometry. Check with independently accumulated endpoint powers.
            dz = force-self.rates[None,:,None]*z
            store = np.sum(tension*np.diff(dz,axis=-1))+np.sum(self.rates[None,:,None]*z*dz)
            power = np.sum(ext*(self.port.T@dz))
            drag = np.sum(dz**2)
            residual = abs(power-store-drag)/max(abs(power),abs(store),drag,1e-28)
            if exchange>1e-7*max(abs(w.amplitude),1e-10) or residual>1e-7:
                raise ArithmeticError('Strip force transfer or frozen-geometry work balance failed')
            emit(test='balance',wall=net,integrated_support=support,exchange_error=exchange,
                 material_work_error=residual,material_drag=drag)
        return dx,force,net,error


def weights(rates, dt):
    integ=np.full_like(rates,dt)
    np.divide(-np.expm1(-dt*rates),rates,out=integ,where=rates!=0)
    return np.exp(-dt*rates)[None,:,None],integ[None,:,None]


def strained_cell(axial, hoop):
    # A separate native-lattice check at measured mean axial spacing. Shear
    # is a difference per reference section, as in the ribbon link energy.
    angle=np.arange(3)*np.pi/3
    ref_x=np.cos(angle)
    xx=axial*ref_x
    yy=hoop*np.sin(angle)
    r=np.hypot(xx,yy)
    tension,k=bond_derivatives(r)
    value=np.sum(ref_x**2*(tension/r+(k-tension/r)*(yy/r)**2))
    def energy(gamma):
        return np.sum(morse(np.hypot(xx,yy+gamma*ref_x),*WITHIN)[0])
    h=1e-4
    check=(energy(h)+energy(-h)-2*energy(0))/h**2
    np.testing.assert_allclose(value,check,rtol=3e-6,atol=2e-4)
    return float(value)


def run(args):
    system=HeldWall(args)
    row,w=system.strip,system.wall
    # Warm the open row while the exact unforced transverse solution is zero.
    warm=4*w.transit
    count=int(np.ceil(warm/args.dt))
    for _ in range(count):
        row.exchange(row.advance(warm/count))
    x=np.zeros((2,w.points,w.families))
    z=np.zeros((2,len(system.rates),len(row.pos)))
    events=[{}, {}, {}]
    tag_tips={}
    early=.0625
    time=0.
    history=[]
    wrst=0.
    tag_gap={}
    initial_inventory=[np.count_nonzero((row.pos>=0)&(row.pos<gate)) for gate in row.gates[1:]]
    # Choose the first twelve post-hold waist arrivals before inspecting forces.
    # Require the same tags to exit before the archived carrier-clock readout.
    cohort=[]
    read_idx=int(np.flatnonzero(system.tags==12)[0])
    for target in [early,28.3125,w.transit]:
        count=int(np.ceil((target-time)/args.dt-1e-12))
        step=(target-time)/count
        half,half_int=weights(system.rates,step/2)
        full,full_int=weights(system.rates,step)
        start=time
        for i in range(count):
            time=start+i*step
            old=row.pos
            new=row.advance(step)
            dx,f,force,error=system.rhs(old,x,z,diagnostics=(i==0))
            history.append([time,force])
            wrst=max(wrst,error)
            mid_x=x+step*dx/2
            mid_z=half*z+half_int*f
            vx,mid_f,_,_=system.rhs((old+new)/2,mid_x,mid_z)
            next_z=full*z+full_int*mid_f
            for gate_no,gate in enumerate(row.gates):
                if np.any((new<gate)&(old>=gate)):
                    raise ArithmeticError('Recrossing requires signed event treatment')
                crossed=np.flatnonzero((old<gate)&(new>=gate))
                for j in crossed:
                    tag=int(row.id0+j)
                    frac=(gate-old[j])/(new[j]-old[j])
                    if tag in events[gate_no]: raise ArithmeticError('Duplicate tag')
                    events[gate_no][tag]=time+frac*step
                    if gate_no==0 and len(cohort)<12: cohort.append(tag)
                    if tag in cohort:
                        tag_gap[tag,gate_no]=(1-frac)*(old[j+1]-old[j])+frac*(new[j+1]-new[j])
                        if j<8 or j+8>=len(new):
                            if gate_no==1: raise ArithmeticError('Missing comparison sites at middle')
                            continue
                        mix=(1-frac)*z+frac*next_z
                        tag_pos=mix[:,:,j]@system.port[:,read_idx]
                        ahead=mix[:,:,j+8]@system.port[:,read_idx]
                        behind=mix[:,:,j-8]@system.port[:,read_idx]
                        tag_tips[tag,gate_no]=(tag_pos-(ahead+behind)/2).tolist()
            x+=step*vx
            retained,added=row.exchange(new)
            z=next_z[:,:,:retained]
            if added:
                z=np.pad(z,((0,0),(0,0),(1,0)))
        time=target
    _,_,final,error=system.rhs(row.pos,x,z,diagnostics=True)
    history.append([time,final]);wrst=max(wrst,error)
    if len(cohort)!=12 or any(tag not in events[2] for tag in cohort):
        raise ArithmeticError('Prespecified tag cohort did not finish by the carrier clock')
    now,forces=np.array(history).T
    def reaction(t):return float(np.interp(t,now,forces))
    for j in range(3):
        gaps=[tag_gap[tag,j] for tag in cohort]
        avg=float(np.mean(gaps))
        emit(test='strain_check',gate=j,mean_forward_gap=avg,hoop=row.hoop_at_gate[j],
             neutral_curvature=strained_cell(1.,1.),local_affine_curvature=strained_cell(avg,row.hoop_at_gate[j]))

    initial=reaction(early)
    for j,end in enumerate(row.gates[1:],1):
        inv=int(np.count_nonzero((row.pos>=0)&(row.pos<end)))
        err=inv-initial_inventory[j-1]-len(events[0])+len(events[j])
        if err:raise ArithmeticError('Unbalanced material label inventory')
        trips=[events[j][tag]-events[0][tag] for tag in cohort]
        clock=float(np.mean(trips))
        sample_forces=[reaction(t) for t in trips]
        crossed_at=[events[j][tag] for tag in cohort]
        emit(entrance_times=[events[0][tag] for tag in cohort],gate_crossing_times=crossed_at,
             wall_forces_at_crossings=[reaction(t) for t in crossed_at],
             sample_force_range=[min(sample_forces),max(sample_forces)],test='tag_clock',gate=j,cohort=cohort,individual_transits=trips,clock=clock,
             entry_count=len(events[0]),exit_count=len(events[j]),initial_inventory=int(initial_inventory[j-1]),
             final_inventory=inv,force_at_clock=reaction(clock),early_force=initial,
             fraction=reaction(clock)/initial if initial else None)
    relative=np.array([tag_tips[tag,1] for tag in cohort])
    emit(test='middle_tag',section_site=12,neighbors_in_each_direction=8,cohort=cohort,
         individual_relative_displacements=relative.tolist(),
         mean_relative_by_component=np.mean(relative,axis=0).tolist(),amplitude=args.amplitude)
    emit(test='joint_force',early_time=early,early=initial,zero_time=reaction(0.),matched_time=28.3125,
         matched=reaction(28.3125),carrier_time=time,carrier_force=final,
         carrier_fraction=final/initial if initial else None,
         peak_cut_difference=wrst,peak_absolute_force=max(abs(forces)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--coupon',type=Path,default=Path('_scratch/s11c/mobile_neighbor_probe_02/response.npz'))
    p.add_argument('--lead',type=float,default=8.)
    p.add_argument('--exit-lead',type=float,default=8.)
    p.add_argument('--dt',type=float,default=1/256)
    p.add_argument('--sub',type=int,default=3)
    p.add_argument('--mesh',type=int,default=34)
    p.add_argument('--amplitude',type=float,default=1e-5)
    p.add_argument('--links',action=argparse.BooleanOptionalAction,default=True)
    p.set_defaults(coverage=1.,brane_grip=0.,phase=0.)
    run(p.parse_args())
