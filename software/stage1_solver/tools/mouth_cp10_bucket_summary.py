"""Combine fixed spatial masks, stored energy, order work and full matched depth work."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read


LEVELS = (1,2,4)
NAMES = ('A','B','C')


def grouped(values, tags, count):
    """Keep the sign and sum both sides of a cancellation in extended precision."""
    return np.stack([values[tags==i].sum(axis=0,dtype=np.longdouble) for i in range(count)])


def as_list(values):
    return np.asarray(values,dtype=float).tolist()


def collect(root, output):
    base = Path(__file__).resolve().parents[3]
    previous = base/'research/pde_ledger_v3/_measurements/mouth_r4_matched_knot.json'
    archive = read(previous)
    if archive['checkpoint_sha256']!=CHECKPOINT_SHA:
        raise ValueError('Unexpected checkpoint archive')
    sources={}
    def checked(folder, filename, guard_name):
        guard=folder/guard_name
        outcome=read(guard/'outcome.json')
        if (outcome['exitCode'] or outcome['childOutcome']['exitCode'] or not outcome['limitsVerified']
                or int(read(guard/'effective-limits.json')['memory.swap.max'])):
            raise ValueError(f'Unverified job {guard}')
        path=folder/filename
        data=read(path)
        if data['checkpoint_sha256']!=CHECKPOINT_SHA:
            raise ValueError(f'Incompatible checkpoint at {path}')
        sources[str(path.relative_to(root))]=dict(sha256=digest(path),command=read(guard/'invocation.json')['command'])
        return data
    geo_dir=root/'r4_bucket_geometry01'
    geo=checked(geo_dir,'geometry.json','guard')
    if digest(geo_dir/'masks.npz')!=geo['masks_sha256']:
        raise ValueError('Cell masks changed')
    with np.load(geo_dir/'masks.npz') as m:
        tags=m['bucket'];bits=m['bits']
        weld=m['weld'].astype(bool)
        fold=m['fold'].astype(bool)
        seam=m['seam'].astype(bool)
        direct=m['fold_hit'].astype(bool)
        if not np.array_equal(tags,np.where(fold|seam,2,np.where(weld,1,0))):
            raise ValueError('Bucket assignment does not implement the documented precedence')
    alt=np.where(weld,1,np.where(fold|seam,2,0))
    labels=[root/'r4_cp10_labels01/labels.npz']+[root/f'r4_nested_labels01/{i}/labels.npz' for i in (2,4)]
    spans=[]
    for f in labels:
        with np.load(f) as data:
            corners=np.lib.stride_tricks.sliding_window_view(data['labels'][...,0],(2,2,2))
            cut=geo['material_weld_chart_radius']
            spans.append((corners.min(axis=(3,4,5))<cut)&(corners.max(axis=(3,4,5))>cut))
    work_dir=root/'r4_bucket_work01'
    transport_dir=root/'r4_bucket_transport01'
    degrees={}
    for degree in (3,6,8):
        case_dir=work_dir/str(degree)
        work=checked(work_dir,f'{degree}/readout.json',f'guard_{degree}')
        if work['cells_sha256']!=digest(case_dir/'cells.npz') or work['label_sha256']!=list(map(digest,labels)):
            raise ValueError('Cell work input changed')
        with np.load(case_dir/'cells.npz') as cells:
            field=cells['spatial'];depth_side=cells['depth_side']
            columns=list(map(str,cells['columns']))
            sums=cells['total']
            primary=grouped(field,tags,3)
            other=grouped(field,alt,3)
            independent=grouped(field,bits,8)
            direct_fold=field[direct].sum(axis=0,dtype=np.longdouble)
            straddling=[as_list(cells[f'cells_{level}x'][mask].sum(axis=0,dtype=np.longdouble))
                        for mask,level in zip(spans,LEVELS)]
            depth_edges=as_list(cells['material_depth_edges'])
        closure=primary.sum(axis=0)-sums
        signed_idx=[columns.index('delta_E'),columns.index('delta_F')]
        max_closure=float(abs(closure[:,signed_idx]).max())
        if max_closure>1e-15 or primary.shape!=(3,3,len(columns)):
            raise ArithmeticError(f'Bucket sums disagree with whole integral: {max_closure}')
        top=[dict(p,bucket=NAMES[tags[tuple(p['native_cell'])]],bits=int(bits[tuple(p['native_cell'])]),
                  straddles_weld=bool(spans[-1][tuple(p['fine_cell'])])) for p in work['top10']]
        row=dict(columns=columns,totals=as_list(sums),order_levels=list(LEVELS),
                 primary=as_list(primary),weld_priority=as_list(other),by_mask_bits=as_list(independent),
                 direct_fold=as_list(direct_fold),actual_spanning_stencil_E_F=straddling,
                 work_arrays_sha256=work['cells_sha256'],max_signed_closure=max_closure,
                 work_checks=work['checks'],depth_edges=depth_edges,depth_side=as_list(depth_side),
                 ordered_points=work['ordered_points'],top10=top)
        degrees[str(degree)]=row
        if degree==3:
            continue
        full=[]
        for level in LEVELS:
            key=f'{level}_{degree}'
            candidates=[g for g in transport_dir.glob(f'*guard_{key}') if (g/'outcome.json').is_file()
                        and read(g/'outcome.json')['exitCode']==0]
            if len(candidates)!=1:
                raise ValueError(f'Expected one admitted transport result for {key}')
            transport=checked(transport_dir,f'{key}/readout.json',candidates[0].name)
            path=transport_dir/key/'cells.npz'
            if transport['fields_sha256']!=digest(path) or transport['label_sha256']!=work['label_sha256'][LEVELS.index(level)]:
                raise ValueError('Wrong transport or bond labels')
            with np.load(path) as cells:
                diffs=np.stack([cells[a]-cells[b] for a,b in (
                    ('force_spatial','force_native'),('cut_spatial','cut_native'),('straddling_spatial','straddling_native'))],axis=-1)
            fs=grouped(diffs,tags,3)
            ref=archive['summary'][key]
            mismatch=float(abs(fs[:,0].sum()-ref['depth_residual']))
            if mismatch>1e-13 or transport['points']!=archive['matched'][key]['steps'][-1]['shared_points']:
                raise ArithmeticError(f'Full depth force changed for {key}')
            full.append(dict(level=level,native_total=ref['native_total'],spatial_total=ref['spatial_total'],
                delta=ref['depth_residual'],primary=as_list(fs),by_mask_bits=as_list(grouped(diffs,bits,8)),
                weld_priority=as_list(grouped(diffs,alt,3)),cell_sum_error=mismatch,fields_sha256=transport['fields_sha256']))
        row['full_depth']=full
    witness=checked(work_dir,'witness_check.json','witness_guard')
    record=dict(input_commit='eda4ffed',checkpoint_sha256=CHECKPOINT_SHA,previous_measurement_sha256=digest(previous),
        phase_drop=archive['phase_drop'],fluid_potential='(n^5 - 5n + 4)/20',geometry=geo,sources=sources,
        bucket_order=list(NAMES),levels=list(LEVELS),depth_side_order=['core','collar'],
        full_depth_columns=['same-sample signed difference','fixed cut signed difference','stencil spanning weld, same-sample subset'],
        degrees=degrees,witness_check=witness)
    with output.open('x') as stream:
        json.dump(record,stream,indent=2,allow_nan=False)
        stream.write('\n')
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    record=collect(args.root.resolve(),args.output.resolve())
    for degree in ('6','8'):
        d=record['degrees'][degree]
        for i,level in enumerate(LEVELS):
            print(degree,level,'order',*[x[i][5] for x in d['primary']],
                  'full',*[x[0] for x in d['full_depth'][i]['primary']])
