"""Export authenticated 4D receiver endpoints, preserving the material positions and cell."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from mouth_bulk_relaxed_shape import PeriodicMaterial


SCRATCH = Path('_scratch/s11c')
LEDGER = Path('research/pde_ledger_v3/_measurements')


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def verify_guard(folder):
    outcome, child, limits, duration = [read(folder/name) for name in
        ('outcome.json','child-outcome.json','limit-validation.json','duration-validation.json')]
    if not (outcome['childOutcome']==child and outcome['exitCode']==child['exitCode']==0 and outcome['limitsVerified']
            and child['guardReason'] is None and limits['verified'] and int(limits['actual']['memory.swap.max'])==0
            and duration['verified'] and duration['actual']['RuntimeMaxUSec']=='infinity' and not (folder/'stderr').read_bytes()):
        raise ValueError(f'Unsuccessful source guard: {folder}')


def write_input(outdir,index,box,pos,temp,mobility,source,**details):
    n=128
    box=np.asarray(box,dtype=float);pos=np.asarray(pos,dtype=float)
    if box.shape!=(4,4) or pos.shape!=(n,4) or not np.isfinite(box).all() or not np.isfinite(pos).all() or not np.linalg.det(box)>0:
        raise ValueError('Invalid endpoint geometry')
    if not (0<=temp<np.inf and 0<mobility<np.inf):
        raise ValueError('Invalid source temperature or volume mobility')
    axis=np.zeros((4,4));axis[0,3]=axis[3,0]=.5
    outdir.mkdir(parents=True,exist_ok=True)
    dest=outdir/f'history{index}.txt'
    with dest.open('x') as stream:
        print(n,index,temp,mobility,2.5,digest(source),file=stream)
        for a in (box,np.linalg.inv(box),axis,pos):
            np.savetxt(stream,a,fmt='%.17g')
    # Seventeen significant digits are sufficient for an exact binary64 round trip.
    saved=np.loadtxt(dest,skiprows=1)
    np.testing.assert_array_equal(saved[:4],box)
    np.testing.assert_array_equal(saved[12:],pos)
    print(json.dumps(dict(history=index,source=str(source),source_sha256=digest(source),file=str(dest),
        export_sha256=digest(dest),source_temperature=temp,volume_mobility=mobility,axis='lab_xw',**details),allow_nan=False))


def export_endpoint(batch,task,index,outdir,guard=None,expected_temp=None):
    done=read(batch/'DONE.json')
    if guard is None:
        jobs=read(batch/'design.json')['jobs']
        group,=[k for k,group in jobs.items() if any(c['name']==task for c in group)]
        guard=batch/('guard_'+group)
    else:
        guard=batch/guard
    verify_guard(guard)
    entries=[json.loads(line) for line in (guard/'stdout').read_text().splitlines()]
    matched,=[line for line in entries if line.get('task')==task]
    command=matched['command']
    src=Path(command[2]);inp=Path(command[1]);exe=Path(command[0]);record=read(src)
    if (src.resolve()!=(batch/'runs'/(task+'.json')).resolve() or digest(src)!=done['results'][task] or
        digest(src)!=matched['result_sha256'] or digest(inp)!=matched['input_sha256'] or digest(exe)!=done['binary_sha256']):
        raise ValueError(f'Changed trajectory or incorrect source provenance for {task}')
    if record['n']!=128 or len(record['blocks'])!=16 or record['seed']!=int(command[3]) or record['mass']!=1 or record['gamma']!=1:
        raise ValueError('Incomplete native endpoint or changed dynamics')
    if expected_temp is not None and (record['temperature']!=expected_temp or record['pressure_external']!=0 or record['strain']!=0):
        raise ValueError('This branch requires an undriven zero-pressure precursor at its declared temperature')
    end=record['final'];last=record['blocks'][-1]['end_reading'];box=np.asarray(end['box']);pos=np.asarray(end['pos']);v=np.asarray(end['vel']);
    np.testing.assert_allclose(np.linalg.det(box),last['volume'],rtol=4.e-13)
    np.testing.assert_allclose(np.sum(v*v)/(4*128),last['temperature'],rtol=2.e-13)
    source_header=inp.read_text().splitlines()[0].split()
    if source_header[-1]!=record['source_sha256']:
        raise ValueError('Mismatched ancestral geometry')
    write_input(outdir,index,box,pos,record['temperature'],float(source_header[3]),src,
        prior_history=record['history'],prior_seed=record['seed'],source_guard=str(guard),
        prior_strain=record['strain'],prior_pressure=record['pressure_external'],
        prior_measured_kinetic_temperature=last['temperature'],coordinate_transfer='exact Cartesian endpoint; velocities resampled on launch')
    return float(source_header[3])


def export_single(src,guard,index,outdir,expected_sha256):
    verify_guard(guard)
    invocation=read(guard/'invocation.json');command=invocation['command'];cwd=Path(invocation['cwd'])
    if not 8<=len(command)<=14 or (cwd/command[2]).resolve()!=src.resolve():
        raise ValueError('The guard invocation did not directly run this native output')
    inp=cwd/command[1];exe=cwd/command[0];head=inp.read_text().splitlines()[0].split()
    if expected_sha256 is not None and digest(src)!=expected_sha256:
        raise ValueError('Raw output hash failed supervisor check')
    r=read(src);n=int(head[0]);source_temp=float(head[2]);base_mobility=float(head[3]);axis_idx=int(command[11]) if len(command)>11 else 0
    target=float(command[12]) if len(command)>12 else source_temp
    if not (r['source_sha256']==head[5] and r['history']==int(head[1]) and r['source_temperature']==source_temp
            and r['temperature']==target and r['n']==n==128 and r['mass']==r['gamma']==1 and r['axis_index']==axis_idx
            and r['seed']==int(command[3]) and r['requested_dt']==float(command[4]) and r['warm']==float(command[5])
            and r['measure']==float(command[6]) and r['strain']==float(command[7]) and r['volume_mobility']==base_mobility*(float(command[9]) if len(command)>9 else 1.)
            and r['pressure_external']==(float(command[10]) if len(command)>10 else 0.) and len(r['blocks'])==16):
        raise ValueError('Native input, settings or guard command disagree')
    np.testing.assert_allclose(sum(b['steps'] for b in r['blocks'])*r['dt'],r['measure'],atol=1.e-12)
    last=r['blocks'][-1]['end_reading'];f=r['final'];box=np.asarray(f['box']);x=np.asarray(f['pos']);v=np.asarray(f['vel']);volume=np.linalg.det(box)
    energy,_,de=PeriodicMaterial(n,2.5).read(x@np.linalg.inv(box),box)
    stress=(box.T@de-v.T@v)/volume
    np.testing.assert_allclose([volume, np.sum(v*v)/(4*n),energy/n,-np.trace(stress)/4,np.sum(stress*np.asarray(r['axis']))],
        [last['volume'],last['temperature'],last['potential_per_site'],last['pressure'],last['shear']],rtol=2.e-10,atol=2.e-9)
    np.testing.assert_allclose(stress,last['stress_tensile'],rtol=2.e-10,atol=2.e-9)
    write_input(outdir,index,box,x,target,base_mobility,src,prior_history=r['history'],prior_seed=r['seed'],source_guard=str(guard),
        checked_output_sha256=expected_sha256,source_input_sha256=digest(inp),native_sha256=digest(exe),
        prior_strain=r['strain'],prior_pressure=r['pressure_external'],recomputed_energy_per_site=energy/n,
        recomputed_pressure=-np.trace(stress)/4,prior_measured_kinetic_temperature=last['temperature'],
        coordinate_transfer='exact Cartesian endpoint; velocities resampled on launch')


def export_cold(outdir,mobility):
    g=SCRATCH/'iso4_pressure_glass41_02';p=g/'native/endpoint.npz'
    verify_guard(g)
    rec=read(LEDGER/'mouth_bulk_isotropic_4d.json')['runs'][g.name]
    literal=(g/'stdout').read_text()
    if literal!=rec['literal_stdout'] or not rec['guard_verified'] or rec['exit_code']!=0 or rec['child_exit']!=0:
        raise ValueError('Cold source guard or archived transcript differs')
    rows=list(map(json.loads,literal.splitlines())); origin,=[r for r in rows if r['test']=='input'];meas,=[r for r in rows if r['test']=='zero_load']
    if not (origin['temperature']==origin['applied_pressure']==0 and digest(Path(origin['source']))==origin['sha256']):
        raise ValueError('Cold source lineage differs')
    with np.load(p,allow_pickle=False) as state:
        assert int(state['seed'])==origin['seed']==meas['seed']==41 and float(state['rc'])==2.5
        y,box=state['y'],state['cell']
        n=len(y);x=y@box;v=np.linalg.det(box)
        np.testing.assert_array_equal(box,meas['cell'])
        energy,gx,de=PeriodicMaterial(n,2.5).read(y,box)
        pcheck=-np.trace(box.T@de)/(4*v)
        np.testing.assert_allclose(energy/n,meas['energy_per_site'],atol=2.e-12,rtol=0)
        np.testing.assert_allclose([pcheck,np.max(np.abs(gx)),n/v], [meas['pressure'],meas['max_site_force'],meas['density']],atol=3.e-11,rtol=0)
    write_input(outdir,141,box,x,0.,mobility,p,prior_history=41,prior_seed=41,source_guard=str(g),
        archived_energy_per_site=meas['energy_per_site'],recomputed_energy_per_site=energy/n,
        recomputed_pressure=pcheck,source_transcript_sha256=digest(g/'stdout'),cold_parent_sha256=origin['sha256'],
        coordinate_transfer='Cartesian positions from archived fractional sites and original lattice; velocities sampled on launch')


def branches(outdir):
    mobilities=[]
    for batch,offset,temp in [('receiver_temperature_return02',0,1.28),('receiver_temperature_bracket01',40,1.088)]:
        for h in range(3):
            mobilities.append(export_endpoint(SCRATCH/batch,f't0h{h}',offset+h,outdir,expected_temp=temp))
    if len(set(mobilities))!=1:
        raise ValueError('The original trials used different base volume-controller mobilities')
    export_cold(outdir,mobilities[0])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='mode',required=True)
    b=sub.add_parser('branches');b.add_argument('outdir',type=Path)
    e=sub.add_parser('endpoint');e.add_argument('batch',type=Path);e.add_argument('task');e.add_argument('index',type=int);e.add_argument('outdir',type=Path);e.add_argument('--guard')
    s=sub.add_parser('single');s.add_argument('raw',type=Path);s.add_argument('guard',type=Path);s.add_argument('index',type=int);s.add_argument('outdir',type=Path);s.add_argument('--expected-sha256')
    args=p.parse_args()
    if args.mode=='branches': branches(args.outdir)
    elif args.mode=='single': export_single(args.raw,args.guard,args.index,args.outdir,args.expected_sha256)
    else: export_endpoint(args.batch,args.task,args.index,args.outdir,args.guard)
