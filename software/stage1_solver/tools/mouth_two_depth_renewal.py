"""Free cubic mouth with quasi-static brane shear and renewing wall attachments.

Two distinct wall pushes use the established carrier state. The rest-slip
reference either travels with the measured wall flow or renews in place.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import linalg, sparse
from scipy.sparse.linalg import splu

from mouth_bend_renewal import ROOT, volume_cost, wall_data
from mouth_bend_static_shear import ANGULAR, archive_inputs
from mouth_carrier_lateral_inertia import taper
from mouth_carrier_lateral_response import prepare
from mouth_empty_bend_stiffness import INTERVALS
from mouth_free_bend_drive import Response
from mouth_free_bend_release import pencil
from mouth_order_wave_pulse import digest, emit
from mouth_two_depth_load import DEPTHS, PEAK_FORCE, polar


GAUGES = np.tri(3)
PORTS = GAUGES[[1, 2]].T


def side_traction(geo, edges, live, quad):
    """Integrate the x component of actual spring traction over the wall.

    For unit dual nodal stress m, the two orthonormal holding tractions
    are m/a² and m'/(2a sqrt(1+a'²)). Project onto motion of each station.
    """
    roots, weights = leggauss(quad+2)
    out = np.zeros((3,len(edges)))
    for j, (left,right) in enumerate(zip(edges[:-1],edges[1:])):
        cuts = [left, *(v for span in INTERVALS for v in span if left<v<right), right]
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            w = lo+(hi-lo)*(1+roots)/2
            a, slope = geo.radius(w)
            arc = np.hypot(1.,slope)
            t = (w-left)/(right-left)
            stress_m = np.column_stack((1-t,t))/a[:,None]**2
            stress_a = np.array([-1.,1.])/(2*a*arc*(right-left))[:,None]
            traction_x = slope[:,None]/arc[:,None]*stress_m+2*stress_a
            b = np.column_stack((np.ones_like(w),*(1-taper(w,*z)[0] for z in INTERVALS)))
            area = ANGULAR*a*a*arc*(hi-lo)*weights/2
            out[:,j:j+2] += b.T@(area[:,None]*traction_x)
    return -out[:,live]


class CoupledMouth:
    def __init__(self, args):
        source = ROOT/f'mouth_flared_cubic_n{args.mesh}'
        order_path = ROOT/f'mouth_free_order_cubic_n{args.mesh}_d04/order_seed1.55.npz'
        dom, field = prepare(SimpleNamespace(source=source))
        geo, chi, _ = archive_inputs(order_path,None)
        meta = json.loads((order_path.parent/'stdout').read_text().splitlines()[0])
        if digest(source/'field.npz')!=meta['source_sha256']:
            raise ValueError('Carrier and order belong to different preparations')
        np.testing.assert_array_equal(geo.points,dom.points)
        self.carrier = Response(dom,field,args.extension)
        self.car_energy,self.car_time = pencil(self.carrier,0.)
        self.nc = self.car_energy.shape[0]
        self.q = slice(self.nc-3,self.nc)
        edges,keep,bulk = volume_cost(geo,chi,args.sub,5,1e-8)
        contact,load,B,G,drift,self.transit,settings = wall_data(geo,chi,field,edges,keep,5)
        metric = bulk+contact
        norm = 1/np.sqrt(metric.diagonal())
        fac = linalg.cho_factor(norm[:,None]*metric*norm)
        def solve(rhs):
            return norm[:,None]*linalg.cho_solve(fac,norm[:,None]*rhs)
        self.dual_bend = solve(load)
        self.dual_ref = solve(B)
        self.side = side_traction(geo,edges,keep,5)
        if args.reference=='local':
            drift = np.zeros_like(drift)
        tau = self.transit
        # E y - T ydot = applied force. p is the two-component unloaded slip.
        # The incoming reference is zero for straight-prepared traveling contacts.
        self.history_bend = G@self.dual_bend/tau
        self.history_rate = drift+G@self.dual_ref/tau
        read = ANGULAR*load.T
        bend_stiff = read@self.dual_bend
        bend_memory = read@self.dual_ref
        nh = drift.shape[0]
        empty = sparse.csc_matrix((self.nc-2,nh))
        extra_k = sparse.block_diag((sparse.csc_matrix((self.nc-2,self.nc-2)),bend_stiff),format='csc')
        top = sparse.vstack((empty,-bend_memory),format='csc')
        lower = sparse.hstack((empty.T,-self.history_bend),format='csc')
        self.energy = sparse.bmat([[self.car_energy+extra_k,top],[lower,self.history_rate]],format='csc')
        self.time = sparse.block_diag((self.car_time,-sparse.eye(nh)),format='csc')
        self.forces = np.zeros((self.energy.shape[0],2))
        self.forces[self.q] = PORTS
        emit(test='material',reference=args.reference,sub=args.sub,dt=args.dt,clock=tau,
             field_hash=meta['source_sha256'],order_hash=digest(order_path),**settings,
             frozen_hold=bend_stiff.tolist(),load_quadrature_error=float(np.max(abs(self.side[1:]+read))),
             peak_force=PEAK_FORCE,station_depths=DEPTHS.tolist())

    def actual_traction(self,y):
        dual = self.dual_bend@y[self.nc-2:self.nc]-self.dual_ref@y[self.nc:]
        return self.side@dual


def pairs(label,value):
    return {label+'_real':value.real.tolist(),label+'_imag':value.imag.tolist()}


def phasor(model,dt):
    s = .30+.40j
    op = model.energy-s*model.time
    fix = splu(op.tocsc()).solve(model.forces)
    factor = splu((op/2-model.time/dt).tocsc())
    state = np.zeros_like(fix)
    for _ in range(round(96/dt)):
        state += factor.solve(model.forces-op@state)
    q = fix[model.q]
    z = model.carrier.force_matrix(s)
    # Material-on-wall force, distinct from additional holding effort.
    inferred = z@q-PORTS
    direct = model.actual_traction(fix)
    from_gauges = linalg.solve(GAUGES.T,inferred)
    through_wall = linalg.solve(GAUGES.T,direct)
    displacement = GAUGES@state[model.q]
    start_error = np.linalg.norm(displacement-GAUGES@q,axis=0)/np.linalg.norm(GAUGES@q,axis=0)
    err = np.max(abs(through_wall-from_gauges))
    response = linalg.solve(q[1:].T,-inferred[1:].T).T
    emit(test='driven',reference_rate=[s.real,s.imag],motion=polar(displacement),
         **pairs('motions',displacement),**pairs('particular_q',q),
         **pairs('station_force',through_wall),**pairs('carrier_subtraction',from_gauges),
         force_gain_per_peak=abs(through_wall).tolist(),
         finite_station_force_gain=abs(linalg.solve(GAUGES.T,model.actual_traction(state))).tolist(),
         switch_on_relative_error=start_error.tolist(),max_traction_error=float(err),
         bend_condition=float(np.linalg.cond(q[1:])),**pairs('bend_holding_response',response),
         **pairs('bare_carrier_matrix',z))
    if not np.isfinite(err) or err>2e-4 or max(start_error)>1e-4:
        raise ArithmeticError('Carrier subtraction, wall integration or finite switch-on failed')


def pulse(model,dt):
    n = model.nc-3
    a,t = model.car_energy,model.car_time
    field_a,field_t = a[:n,:n],t[:n,:n]
    replay = splu((field_a/2-field_t/dt).tocsc())
    # Replay receives only the observed gauge path and the bare carrier law.
    field_q,time_q = a[:n,n:],t[:n,n:]
    wall_f,wall_t,wall_q = a[n:,:n],t[n:,:n],a[n:,n:]
    advance = splu((model.energy/2-model.time/dt).tocsc())
    state = np.zeros_like(model.forces)
    bare = np.zeros((n,2))
    requested = {4.,8.,16.,round(model.transit/dt)*dt}
    checkpoints = {}
    worst = 0.
    for step in range(round(max(requested)/dt)+1):
        mid = (step+.5)*dt
        f = np.sin(np.pi*mid/8)**2 if mid<8 else 0.
        dy = advance.solve(f*model.forces-model.energy@state)
        ym = state+dy/2
        x = GAUGES@state[model.q]
        dx = GAUGES@dy[model.q]
        qm = linalg.solve(GAUGES,x+dx/2)
        dqdt = linalg.solve(GAUGES,dx/dt)
        dh = replay.solve(-field_a@bare-field_q@qm+time_q@dqdt)
        bare_required = wall_f@(bare+dh/2)+wall_q@qm-wall_t@(dh/dt)
        wall = linalg.solve(GAUGES.T,model.actual_traction(ym))
        subtracted = linalg.solve(GAUGES.T,bare_required-f*PORTS)
        err = np.max(abs(wall-subtracted))
        worst = max(worst,err)
        state += dy
        bare += dh
        for when in requested:
            if abs(mid-when)<dt*.51:
                checkpoints.setdefault(when,[]).append((wall,subtracted))
                if mid<when:
                    emit(test='pulse_motion',t=when,per_peak=(GAUGES@state[model.q]).tolist())
    for when,rows in sorted(checkpoints.items()):
        if len(rows)!=2:
            raise ValueError('Missing a midpoint bracketing a force reading')
        wall,subtracted=np.mean(rows,axis=0)
        emit(test='pulse_force',t=when,station_force_per_peak=wall.tolist(),
             from_bare_replay=subtracted.tolist(),error=float(np.max(abs(wall-subtracted))))
    emit(test='pulse_replay',end=max(requested),max_midpoint_force_error_per_peak=float(worst))
    if not np.isfinite(worst) or worst>2e-4:
        raise ArithmeticError('Prescribed-motion carrier replay failed to reproduce wall traction')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,required=True,choices=[26,34])
    p.add_argument('--extension',choices=['compact','wide'],required=True)
    p.add_argument('--dt',type=float,choices=[.25,.125],required=True)
    p.add_argument('--sub',type=int,choices=[1,2,3,4],required=True)
    p.add_argument('--reference',choices=['travelling','local'],required=True)
    args=p.parse_args()
    model=CoupledMouth(args)
    phasor(model,args.dt)
    pulse(model,args.dt)
