"""Measure an entrance-tag trajectory for the 0PN first-bend material test."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.interpolate import PchipInterpolator

from mouth_connected_passage import WAIST, MIDDLE
from mouth_empty_bend_stiffness import INTERVALS, shapes
from mouth_carrier_attachment_tension import area_stiffness
from mouth_joint_strip_hold import TravelingStrip
from mouth_middle_tag_grip import cubic_trace
from mouth_middepth_rearrangement import RECORD, wall_table
from mouth_order_wave_pulse import digest, emit
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import COUNT, WITHIN
from neutral_patch_turning_neighbor import ENTRY

PREVIOUS = Path('research/pde_ledger_v3/_measurements/mouth_middepth_bend_wall.json')


def prepare(directory):
    previous = json.loads(PREVIOUS.read_text())
    for file, expected in previous['hashes'].items():
        if digest(Path(file)) != expected:
            raise ValueError(f'Previous input changed: {file}')
    old = [json.loads(x) for x in json.loads(RECORD.read_text())['guarded_runs']['resolution']['literal_stdout']]
    setup = next(r for r in old if r['test']=='setup')
    gates = [r for r in old if r['test']=='tag_clock']
    mid_gap = next(r['mean_forward_gap'] for r in old if r['test']=='strain_check' and r['gate']==1)
    trace, intake = cubic_trace(34)
    strip = TravelingStrip(trace, SimpleNamespace(lead=8.,exit_lead=8.))
    dt = setup['dt']
    warm = 4*setup['transit']
    n = int(np.ceil(warm/dt))
    for _ in range(n):
        strip.exchange(strip.advance(warm/n))

    # Track an individual already just upstream of the fixed waist. All sites
    # in the transverse experiment use time translations of this measured path.
    i = int(np.searchsorted(strip.pos, -4.))
    tag = strip.id0+i
    t = 0.
    rows = []
    while True:
        i = tag-strip.id0
        s = strip.pos[i]
        rows.append((t,s,strip.pos[i+1]-s))
        if s>=strip.gates[-1]: break
        strip.exchange(strip.advance(dt))
        t+=dt
    tm, s, gaps = np.array(rows).T
    t_at = PchipInterpolator(s,tm)
    t0=float(t_at(0.)); clock=float(t_at(strip.gates[1])-t0)
    period=float(np.diff(gates[0]['gate_crossing_times']).mean())
    sampled_gap = float(np.interp(strip.gates[1],s,gaps))
    np.testing.assert_allclose(clock,gates[0]['clock'],rtol=.001)
    np.testing.assert_allclose(sampled_gap,mid_gap,rtol=.001)

    # A steady train injects one section every measured period at the waist.
    # Keep full measured entry-to-gauge passage and an open downstream buffer.
    reg=(tm-t0)/period
    sample=np.arange(0, int(np.floor(reg[-1]))*256+1)/256
    arc=PchipInterpolator(reg,s)(sample)
    gap=PchipInterpolator(reg,gaps)(sample)
    w=np.interp(arc,strip.arc,strip.depths)
    radius=trace.geometry.radius(w)[0]
    hoop=radius/trace.geometry.radius(np.array([WAIST]))[0][0]
    order=np.clip(trace.sample(trace.chi,w)[0],0,1)
    h=shapes(radius,w,'moving-box',trace.geometry.far,trace.geometry.roof)[0][:,0]
    wall=gap*hoop**2*order*sum(setup['neutral_tangents'])
    mid=strip.gates[1]
    bounds=[mid-6*mid_gap,mid+6*mid_gap]
    directory.mkdir(parents=True,exist_ok=True)
    path=directory/'path.dat'
    table=directory/'wall.dat'
    wall_table(table,setup['neutral_tangents'])
    with path.open('w') as f:
        header=[period,COUNT,ENTRY,*WITHIN,*bounds,mid,256,len(sample)]
        f.write(' '.join(map(str,header))+'\n')
        np.savetxt(f,np.column_stack((arc,gap,hoop,order,h,wall)),fmt='%.16g')
    wall_geometry = SimpleNamespace(blocks={2: (trace.edges,)}, radius=trace.geometry.radius)
    k,fd=area_stiffness(wall_geometry,INTERVALS[0])
    emit(test='tag_path',intake=intake,waist=WAIST,middle=MIDDLE,tag=tag,
         tag_clock=clock,archived_clock=gates[0]['clock'],arrival_period=period,
         gap_at_middle=sampled_gap,archived_gap=mid_gap,gauge_arc=bounds,middle_arc=mid,
         middle_registry=clock/period,gauge_registry=PchipInterpolator(arc,sample)(bounds).tolist(),
         max_sampled_registry=sample[-1],min_gap=min(gap),max_gap=max(gap),order_range=[min(order),max(order)],
         prior_sha256=digest(PREVIOUS),path_sha256=digest(path),wall_sha256=digest(table))
    emit(test='fixed_shape_tension',bend=INTERVALS[0],area_stiffness=k,finite_difference=fd,
         trial_tensions=[.02,.04],amplitude_in_radii=5e-6,incremental_holder_difference=(.04-.02)*k*5e-6)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory',type=Path)
    prepare(p.parse_args().directory)
