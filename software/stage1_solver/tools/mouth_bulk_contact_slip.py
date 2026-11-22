"""Differential tangential traction between two precontacted samples at fixed area."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_brace_grip import Material, clock
from mouth_bulk_contact import Tracking, converted, material_pair, observations
from mouth_middle_tag_grip import cubic_trace
from neutral_patch_structural_endpoints import PITCH
from neutral_patch_open_w import mix
from neutral_patch_sites import BETWEEN, TEMPERATURE, WITHIN, emit, morse


def traction(m, pos, sign, check):
    cell = m.cell
    cross = (cell.i < m.n//2) & (cell.j >= m.n//2)
    d = pos[cell.j[cross]] - pos[cell.i[cross]]
    d[:,:2] += cell.shifts[cross]-cell.box*np.rint(d[:,:2]/cell.box)
    r = np.linalg.norm(d,axis=1)
    d,r = d[r<cell.cutoff],r[r<cell.cutoff]
    f = mix(d[:,2],m.power)[0]
    radial = (1-f)*morse(r,*WITHIN)[1]+f*morse(r,*BETWEEN)[1]
    tension = np.sum(radial*d[:,0]/r)/cell.area
    if check:
        _,grad,_ = cell.evaluate(pos)
        np.testing.assert_allclose(sign@grad[:,0]/cell.area,tension,atol=2e-12)
        move = np.zeros_like(pos)
        move[:,0] = 2e-7*sign
        ep = cell.evaluate(pos+move)[0]
        em = cell.evaluate(pos-move)[0]
        np.testing.assert_allclose((ep-em)/(4e-7*cell.area),tension,atol=2e-7)
    return tension


def run(args):
    args.work.mkdir(parents=True,exist_ok=True)
    args.mesh,args.cutoff,args.side,args.brace,args.mobility,args.temperature = 34,4.5,4,'elastic',1.,TEMPERATURE
    depth,times,order,_ = clock(args.mesh)
    trace,_ = cubic_trace(args.mesh)
    transit = trace.transit(depth[-1],trace.edges[-1],8)
    base = Material(args)
    a = converted(base,319,times,order)
    b = converted(base,817,times,order)
    b[:,2] += max(a[:,2])-min(b[:,2])+PITCH
    b[:,:2] += PITCH*np.array([.5,np.sqrt(3)/6])
    pos = np.vstack([a,b]); pos[:,2] -= pos[:,2].mean()
    m = material_pair(base,True)
    rng = np.random.default_rng(319+23*817)
    count = int(np.ceil(transit/48/.001))
    dt = transit/(48*count)
    def blocks():
        for start in range(0,count,128):
            k = min(128,count-start)
            yield np.zeros(k),np.sqrt(2*TEMPERATURE*dt)*rng.normal(size=(k,m.n,3))
    for _ in range(48):
        for orders,noise in blocks():
            m.step(pos,dt,orders,noise)
    # Match the independent common-strain run's warm state before using it.
    # That trajectory supplies no tangential displacement or response here.
    folder = Path(f'_scratch/s11c/bulk_contact_p{args.power}_01')
    outcome = json.loads((folder/'outcome.json').read_text())
    assert outcome['limitsVerified'] and outcome['exitCode']==outcome['childOutcome']['exitCode']==0
    rows = map(json.loads,(folder/'stdout').read_text().splitlines())
    reference, = [r for r in rows if r['test']=='preparation' and r.get('contact')=='joined' and r['fraction']==1]
    whole,normal = observations(m.cell,pos)
    np.testing.assert_allclose([*whole,normal,np.ptp(pos[:,2])],
                               [*reference['rest_stress'],reference['bond_normal_traction'],reference['height_span']],atol=1e-10)
    sign = np.repeat([-.5,.5],m.n//2)
    steps = [0.,1e-4,-1e-4,5e-5,-5e-5]
    states = [pos.copy() for amp in steps]
    for state,amp in zip(states,steps):
        state[:,0] += sign*amp
    separation = (2/(m.n//2))*np.dot(sign,pos[:,0])
    isolated = material_pair(base,False)
    tracker = Tracking(pos,m.cell.box,True)
    emit(test='slip_setup',power=args.power,seeds=[319,817],transit=transit,area=m.cell.area,
         temperature=TEMPERATURE,restored_traction=normal,step=steps,axis='first tangential period')
    for frame in range(49):
        if frame:
            for orders,noise in blocks():
                for state in states:
                    m.step(state,dt,orders,noise,hold=True)
                tracker.update(states[0],m.cell.box)
        loads = np.array([traction(m,state,sign,frame in (0,48)) for state in states])
        error = np.array([2*np.dot(sign,x[:,0])/(m.n//2)-separation for x in states])-steps
        disconnected = traction(isolated,states[1],sign,frame in (0,48))
        emit(test='slip',fraction=frame/48,signed_traction=loads.tolist(),full=(loads[1]-loads[2])/2e-4,
             half=(loads[3]-loads[4])/1e-4,disconnected_same_positions=disconnected,
             centroid_error=max(abs(error)),**tracker.update(states[0],m.cell.box))
        if max(abs(error))>1.e-9 or disconnected!=0:
            raise ArithmeticError('Invalid held step or disconnected virtual traction')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--power',type=int,choices=(1,2),required=True)
    run(parser.parse_args())
