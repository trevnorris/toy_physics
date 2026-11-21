"""Run the local flowing patch at the archived mid-depth feed and wall coverage."""
import argparse
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from mouth_connected_passage import MIDDLE
from mouth_middle_tag_grip import cubic_trace
from mouth_order_wave_pulse import digest, emit
from mouth_segment_contact_hold import NEUTRAL, RECORD as CONTACT_RECORD
from neutral_patch_sites import BETWEEN, COUNT, GAP, WITHIN, morse
from neutral_patch_turning_neighbor import ENTRY

RECORD=Path('research/pde_ledger_v3/_measurements/mouth_joint_strip_hold.json')


def wall_table(path, tangents):
    # Spherical orientation average of the nine original offsets in three
    # wall tangents, with the same normal separation as the contact calibrator.
    ids=np.array(json.loads(CONTACT_RECORD.read_text())['independent_id_extraction']['laws']['central']['strict_surviving_ids'])
    xy=np.array(json.loads(NEUTRAL.read_text())['saved_measurements']['frames1321.json']['contents']['central']['rested_xy'])
    offsets=np.linalg.norm(xy[ids[:,1]]-xy[ids[:,0]],axis=1)
    cos,weights=np.polynomial.legendre.leggauss(16)
    step=.001
    disp=np.arange(5001)[:,None,None]*step
    projected=disp+offsets[None,:,None]*cos
    r=np.sqrt(GAP**2+offsets[None,:,None]**2+disp**2+2*disp*offsets[None,:,None]*cos)
    energy,tension=morse(r,*BETWEEN)
    norm=sum(tangents)
    e=((energy-energy[0])@(weights/2)).sum(axis=1)/norm
    slope=((tension*projected/r)@(weights/2)).sum(axis=1)/norm
    slope[0]=0.
    np.testing.assert_allclose(slope[1]/step,1.,rtol=1e-5)
    with path.open('w') as out:
        out.write(str(step)+'\n')
        np.savetxt(out,np.stack((e,slope),axis=1),fmt='%.17g')
    emit(test='wall_law',table_sha256=digest(path),slope_at_test_offset=slope[350],
         offset=.35,zero_curvature_check=slope[1]/step)


def run(args):
    records=json.loads(RECORD.read_text())
    rows=list(map(json.loads,records['guarded_runs']['resolution']['literal_stdout']))
    def get(test, **filters):
        return next(r for r in rows if r['test']==test and all(r[k]==v for k,v in filters.items()))
    mid=get('strain_check',gate=1)
    crossing=get('tag_clock',gate=1)['gate_crossing_times']
    period=float(np.mean(np.diff(crossing)))
    gap,hoop=mid['mean_forward_gap'],mid['hoop']
    trace,intake=cubic_trace(34)
    order=float(trace.sample(trace.chi,MIDDLE)[0])
    # Continuum wall quadrature per passing section: integral of its linear
    # overlap is one section spacing. Fixed tips, nine registered families.
    wall=gap*hoop**2*order*sum(get('setup')['neutral_tangents'])
    depth,alpha,rest=WITHIN
    vals=[args.length,args.lanes,args.sub,args.seed,args.graph,gap,hoop,period,
          order,wall,args.step,rest,alpha,depth,ENTRY,COUNT,args.cutoff,8,args.wave,args.sign]
    emit(test='local_setup',intake=intake,station=MIDDLE,reference_gap=gap,hoop=hoop,
         order=order,time_per_section=period,site_speed=gap/period,contact_per_section=wall,
         length=args.length,lanes=args.lanes,sub=args.sub,seed=args.seed,sign=args.sign,
         step=args.step,wave=args.wave,window_halfwidth=8,cutoff=args.cutoff,graph=args.graph,
         record=digest(RECORD),kernel=digest(args.binary),source=digest(Path(__file__).with_suffix('.cpp')))
    cmd=[str(args.binary.resolve()),*map(str,vals)]
    if args.wall=='morse':
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'wall.txt'
            wall_table(path,get('setup')['neutral_tangents'])
            subprocess.run([*cmd,str(path)],check=True)
    else:
        subprocess.run(cmd,check=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True)
    p.add_argument('--length',type=int,default=32)
    p.add_argument('--lanes',type=int,default=4)
    p.add_argument('--sub',type=int,default=128)
    p.add_argument('--seed',type=int,default=1)
    p.add_argument('--sign',type=int,choices=[-1,1],default=1)
    p.add_argument('--wall',choices=['harmonic','morse'],default='harmonic')
    p.add_argument('--graph',choices=['all','fixed'],default='all')
    p.add_argument('--step',type=float,default=1e-4)
    p.add_argument('--wave',type=float,default=8.)
    p.add_argument('--cutoff',type=float,default=3.)
    run(p.parse_args())
