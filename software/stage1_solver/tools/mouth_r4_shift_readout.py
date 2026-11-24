"""Check and assemble the splice-free unfolded interior-grid displacement control."""
import argparse
from pathlib import Path

import numpy as np

from mouth_cp10_unfolded_depth import H,DEGREES,LEVELS,read,write_new
from mouth_exact_checkpoint_control import CHECKPOINT_SHA,digest


def reduce(root,area,out):
    prev=root/'r4_bucket_unfolded01'
    key=digest(prev/'control.npz')
    old=read(prev/'depth/result.json')
    source_hashes={}
    def checked(name,path):
        guard=area/f'guard_{name}'
        g=read(guard/'outcome.json')
        if (g['exitCode'] or not g['limitsVerified'] or g['childOutcome']['exitCode'] or
                int(read(guard/'effective-limits.json')['memory.swap.max'])!=0):
            raise RuntimeError('Guard failed '+str(guard))
        source_hashes[str(path)]=dict(sha256=digest(path),command=read(guard/'invocation.json')['command'])
        item=read(path)
        if item['parent_sha256']!=CHECKPOINT_SHA or item['control_sha256']!=key:
            raise ValueError('Frozen state provenance changed')
        return item
    labels=[checked(f'labels_{level}',area/f'labels_{level}/result.json') for level in LEVELS]
    for item in labels:
        for sign,term in item['cases']['unshifted']['tests'].items():
            # Signed depth banks use the same direct inverse. The static 1x
            # preparation includes a roundoff-sized registered harmonic trace snap.
            tolerance=np.finfo(float).eps if sign=='zero' else 0.
            if term['active_max_abs_difference']>tolerance:
                raise ValueError('Direct and spliced labels differ on ordered cell closure')
    slopes=checked('kinks',area/'kinks.json')
    cases=[checked(f'sector_{k}',area/f'sector_{k}/result.json') for k in range(6)]
    if [r['sector'] for r in cases]!=list(range(6)) or len({r['masks_sha256'] for r in cases})!=1:
        raise ValueError('Missing angular sector or inconsistent order front')
    table={};sensitivity={}
    for degree in DEGREES:
        fixed=checked(f'static_{degree}',area/f'static_{degree}/result.json')
        arr=np.array(fixed['raw_by_region'])
        for level in LEVELS:
            for kind in ('unshifted','shifted'):
                run=f'{kind}_{level}_{degree}';parts=[r['rows'][run] for r in cases]
                summed=lambda k:np.asarray(sum(np.asarray(x[k],dtype=np.longdouble) for x in parts),dtype=float)
                forces=summed('force'); local=summed('local');energies=summed('energies')
                if max(abs(forces-local.sum(axis=0)))>1e-13 or max(abs(forces-(energies[0]-energies[1])/(2*H)))>1e-8:
                    raise ArithmeticError('Stored split vs signed force')
                bank=fixed['axes_runs'].index(f'{kind}_{level}')
                native=arr[:,0,:];spatial=arr[:,bank,:];error=spatial-native
                energy,order=error.sum(axis=0)
                ftotal=forces[1]-forces[0]
                lf=local[:,1]-local[:,0]
                points=sum(r['points'] for r in parts)
                historical=old['rows'][f'{level}_{degree}']
                error_old=float(ftotal-historical['spatial_minus_native'])
                if kind=='unshifted' and (abs(error_old)>1e-12 or points!=historical['points']):
                    raise ArithmeticError(f'Published unshifted forces changed: {error_old}, {points}')
                table[run]=dict(energy=float(energy),order=float(order),remaining=float(ftotal-order),depth=float(ftotal),
                    raw={'force':(forces+old['common_bath_force']).tolist(),'E':arr.sum(axis=0)[[0,bank],0].tolist(),
                         'order':arr.sum(axis=0)[[0,bank],1].tolist()},
                    local=[dict(region=reg,energy=float(ee),order=float(oo),depth=float(ff),remaining=float(ff-oo))
                           for reg,(ee,oo),ff in zip(fixed['regions'],error,lf)],
                    archival_depth_delta=error_old,quadrature_points=points,
                    static_ordered_points=fixed['ordered_points'],min_sampled_Jinv=min(x['min_Jinv'] for x in parts))
            base,shift=[table[f'{k}_{level}_{degree}'] for k in ('unshifted','shifted')]
            sensitivity[f'{level}_{degree}']={metric:float(shift[metric]-base[metric]) for metric in ('energy','order','remaining','depth')}
    output=dict(parent_commit='344030b4',parent_sha256=CHECKPOINT_SHA,control_sha256=key,
        phase_drop=0.4510084282857874,fluid_potential='(n^5 - 5n + 4)/20',common_field_depth_force=old['common_bath_force'],
        convention='Eulerian minus original material, F = -dE/dlogL. Remaining = total error minus explicit order-shift error.',
        geometry='Same held CP10 geometry and order, half the native material displacements. One-quarter forward local cell width displacement of every interior normalized label-depth knot only. Order grid and domain endpoints identical for each signed motion.',
        h=H,labels=labels,native_face_kinks=slopes,
        source=source_hashes,ordered_regions=['outside seam neighborhood','within CP10 native χ-front cell neighborhood'],
        code_sha256={Path(__file__).with_name(n).name:digest(Path(__file__).with_name(n))
                     for n in ('mouth_r4_shift_inverse.py','mouth_r4_shift_work.cpp','mouth_r4_shift_work.py','mouth_r4_shift_depth.py','mouth_r4_shift_kinks.py')},
        rows=table,shift_minus_unshifted=sensitivity)
    write_new(out,output)
    for degree in DEGREES:
        for level in LEVELS:
            print(level,degree,sensitivity[f'{level}_{degree}'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('area',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();reduce(a.root.resolve(),a.area.resolve(),a.output.resolve())
