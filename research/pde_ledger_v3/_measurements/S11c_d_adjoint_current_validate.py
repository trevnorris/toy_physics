#!/usr/bin/env python3
"""Validate an emitted adjoint-current packet and optionally publish it."""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import pickle
import shutil
import tempfile
import numpy as np
import sympy as sp
from S11c_d_modal_current_check import engine,digest
from S11c_d_joint_sheet_check import decoded_lines,_restore
from S11c_d_output_codec import restore_emission_index


def assoc(value):return {str(k):v for k,v in value}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--run-directory',type=Path,required=True)
    parser.add_argument('--modal-checkpoint',type=Path,required=True)
    parser.add_argument('--publish',action='store_true')
    args=parser.parse_args();base=args.run_directory;root=Path(engine.__file__).resolve().parents[1]
    summary=json.loads((base/'checks.json').read_text());provenance=summary['provenance']
    sources={'scripts/S11c_d_mixing_scattering_sympy_audit.py':'engineSha256',
             '_measurements/S11c_d_adjoint_current_check.py':'instrumentSha256'}
    for name,key in sources.items():
        if digest(root/name)!=provenance[key] or digest(base/'source'/name)!=provenance[key]:
            raise ValueError(('source drift',name))
    for name,sha in provenance['producerSources'].items():
        if name!='scripts/S11c_d_mixing_scattering_sympy_audit.py' and digest(root/name)!=sha:
            raise ValueError(('dependency drift',name))
    if digest(args.modal_checkpoint)!=provenance['modalCheckpointSha256']:
        raise ValueError('modal inventory pin')
    old_inventory=json.loads(args.modal_checkpoint.read_text());old_file=Path(old_inventory['runDirectory'])/'objects.pickle'
    if digest(old_file)!=provenance['modalCacheSha256']:raise ValueError('source modal payload pin')
    old,_=pickle.loads(old_file.read_bytes())
    if digest(base/'objects.pickle')!=summary['objectsSha256']:raise ValueError('new payload pin')
    result,known=pickle.loads((base/'objects.pickle').read_bytes())
    if (base/'stderr.txt').stat().st_size:raise ValueError('run stderr')
    progress=[json.loads(line) for line in (base/'progress.jsonl').read_text().splitlines()]
    if progress[-1]['stage']!='completed':raise ValueError('unfinished run')
    records={};tags=[]
    for line in decoded_lines(base/'full.out'):
        tag,sep,payload=line.rstrip('\n').partition(': ')
        if not sep or not tag.startswith('PY_S11CD_') or tag in records:raise ValueError(('tag',tag))
        records[tag]=_restore(payload);tags.append(tag)
    missing=[tag for tag in tags if not tag.startswith('PY_S11CD_METADATA_') and
             'PY_S11CD_METADATA_'+tag.removeprefix('PY_S11CD_') not in records]
    if missing:raise ValueError(('metadata missing',missing))
    if any('ZERO_MAP' in sp.srepr(v) for tag,v in records.items() if tag.startswith('PY_S11CD_METADATA_')):
        raise ValueError('unassigned zero dimension')
    index_tag='PY_S11CD_ADJOINT_CURRENT_MAP_EMISSION_LINES'
    indexed=restore_emission_index(assoc(records[index_tag]),tags[:tags.index(index_tag)])
    if records['PY_S11CD_ADJOINT_CURRENT_MAP_DIMENSION_CONSTRAINTS']:raise ValueError('dimension constraints')
    metadata_paths=0;checked=0;literal_scalars=0
    epsilon=next(s for s in known if str(s)=='epsilon_shape')
    def meta(tag,body,units=None):
        nonlocal metadata_paths
        paths=[]
        for descriptor in records['PY_S11CD_METADATA_'+tag.removeprefix('PY_S11CD_')]:
            item=assoc(descriptor)
            if len(item['DIMENSION_L_T_M'])!=3:raise ValueError('dimension arity')
            for path in item['PATHS']:
                restored=tuple(str(v) if isinstance(v,sp.core.symbol.Str) else int(v) for v in path)
                paths.append(restored)
                if units is not None and tuple(item['DIMENSION_L_T_M'])!=units[restored[0]]:
                    raise ValueError(('unit restoration',tag,restored))
        expected=[path for path,v in engine.leaves(engine.cas(body)) if not isinstance(v,sp.core.symbol.Str)]
        if Counter(paths)!=Counter(expected):raise ValueError(('metadata census',tag))
        metadata_paths+=len(paths)
    symbolic_zero=0
    for group in ('SYMBOLIC_OPERANDS','SYMBOLIC_RESIDUALS'):
        for name,body in result[group].items():
            tag='PY_S11CD_ADJOINT_CURRENT_MAP_'+group+'_'+name
            expected=engine.carrier_fingerprint(body) if group=='SYMBOLIC_OPERANDS' else body
            if records[tag]!=expected:raise ValueError(('symbolic payload',tag))
            meta(tag,body)
            if group=='SYMBOLIC_RESIDUALS':
                if any(v!=0 for v in body):raise ValueError(('symbolic residual',name))
                symbolic_zero+=len(body)
    if len(result['RECORDS'])!=len(old['RECORDS']):raise ValueError('missing source candidates')
    comparisons=[]
    for mode,source in zip(result['RECORDS'],old['RECORDS']):
        tag='PY_S11CD_ADJOINT_CURRENT_MAP_REFERENCE_LAB_HELD_RHO4_CONSTANT_'+str(mode['INDEX'])
        for key in ('INDEX','ROOT_DISK_INDEX','NORMAL_LIFT_SIGN','NULLITY','K','Q'):
            if mode[key]!=source[key]:raise ValueError(('source candidate join',key))
        if mode['SOURCE_CURRENT_NORMALIZATION_DEFINED']!=source['PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED']:
            raise ValueError('normalization gate changed')
        info={k:v for k,v in mode.items() if k!='ITEMS'}
        if records[tag+'_RECORD']!=engine.cas(info):raise ValueError('record payload')
        meta(tag+'_RECORD',info)
        lookup={(v['GROUP'],v['NAME']):v['VALUE'] for v in mode['ITEMS']}
        for item in mode['ITEMS']:
            array=item['VALUE'];name=tag+'_'+item['GROUP']+'_'+item['NAME']
            if not np.isfinite(array).all():raise ValueError(('nonfinite tensor',name))
            body=sp.ImmutableMatrix(*array.shape,[engine.FullPencilModes.number(v) for v in array.ravel()])*epsilon**item['EPSILON_POWER']
            if item['GROUP']=='RESIDUALS':
                if records[name]!=body:raise ValueError(('literal residual',name))
                norm=float(np.linalg.norm(array));literal_scalars+=array.size
                if norm!=mode['RESIDUAL_NORMS'][item['NAME']] or norm>1e-8:
                    raise ValueError(('numeric reconstruction diagnostic',name,norm))
            else:
                fingerprint=assoc(records[name]);expected_sha=hashlib.sha256(sp.srepr(engine.cas(body)).encode()).hexdigest()
                if str(fingerprint['OBJECT_SHA256'])!=expected_sha:raise ValueError(('tensor sha',name))
                samples=[]
                for index in range(3):
                    total=0j
                    for path,v in engine.leaves(body):
                        seed=hashlib.sha256((str(path)+':'+str(index)).encode()).digest()
                        total+=(int.from_bytes(seed[:2],'big')%97+1)/101*complex(v.subs(epsilon,sp.Rational(index+3,107)))
                    samples.append(engine.FullPencilModes.number(total))
                if tuple(fingerprint['NUMERIC_UNIT_FRAME_TENSOR_PIT'])!=tuple(samples):raise ValueError(('tensor PIT',name))
            meta(name,body,item['UNITS']);checked+=1
        if mode['INVERTIBLE_FIELD_MAP_DEFINED']:
            b=lookup['OPERANDS','POWER_MAP'];p=lookup['OPERANDS','PENCIL'];a=lookup['MAPS','ADJOINT_FIELD'];r=lookup['OPERANDS','RIGHT']
            n=np.linalg.matrix_rank(a,tol=1e-9)
            if n!=mode['NULLITY']:raise ValueError('incomplete adjoint field basis')
            norms={key:float(np.linalg.norm(lookup['RESIDUALS',key])) for key in
                ('ROW_TO_FIELD','PHYSICAL_FIELD_RECONSTRUCTION','CURRENT_FINITE_BRIDGE_RECONSTRUCTION',
                 'NORMAL_MIXED_CURRENT_ENERGY_RECONSTRUCTION','FREQUENCY_MIXED_CURRENT_ENERGY_RECONSTRUCTION')}
            mixed=lookup['FORMS','CURRENT_FINITE_MIXED']
            comparisons.append({'index':mode['INDEX'],'nullity':mode['NULLITY'],'adjointFieldRank':int(n),
                'diagnosticNorms':norms,'mixedCurrentReal':mixed.real.tolist(),'mixedCurrentImaginary':mixed.imag.tolist(),
                'bridgeReal':lookup['MAPS','PHYSICAL_BRIDGE'].real.tolist(),
                'bridgeImaginary':lookup['MAPS','PHYSICAL_BRIDGE'].imag.tolist(),
                'fieldDefectNorm':float(np.linalg.norm(lookup['MAPS','PHYSICAL_FIELD_DEFECT'])),
                'sourcePhysicalNormalization':mode['SOURCE_CURRENT_NORMALIZATION_DEFINED']})
    inventory={**summary,'runDirectory':str(base),'exitCode':0,'tagCount':len(tags),
        'metadataTagCount':sum(tag.startswith('PY_S11CD_METADATA_') for tag in tags),'sourceAssignments':len(indexed),
        'checkedTensorCount':checked,'checkedMetadataPaths':metadata_paths,'literalResidualScalars':literal_scalars,
        'checkedSymbolicZeroScalars':symbolic_zero,'sourceFiles':{name:digest(root/name) for name in sources},
        'validationSourceSha256':digest(Path(__file__)),'perModeMapSummary':comparisons,
        'artifacts':{name:{'bytes':(base/name).stat().st_size,'sha256':digest(base/name)} for name in
                     ('full.out','objects.pickle','checks.json','progress.jsonl','stderr.txt')},
        'scope':'regular reference root subspaces; power-weighted row representation and physical-current reconstruction',
        'pending':['BOTH_END_CURRENT_AND_MATCHING','FULL_CASE_INTEGRATION','GLOBAL_EXCEPTIONAL_COVERAGE',
                   'SCATTERING_PROFILE_FREQUENCY_POLES_BOOKKEEPING_EXPORT']}
    if args.publish:
        target=root/'scripts/out/S11c_d_adjoint_current_check.out'
        if target.exists() or target.is_symlink():raise ValueError('new publication path already exists')
        with tempfile.NamedTemporaryFile(dir=target.parent,prefix='.s11cd-adjoint-',delete=False) as stream:
            temporary=Path(stream.name)
            with (base/'full.out').open('rb') as src:shutil.copyfileobj(src,stream)
            stream.flush();os.fsync(stream.fileno())
        if digest(temporary)!=digest(base/'full.out'):raise ValueError('publication copy')
        os.replace(temporary,target)
        inventory['publication']={'path':str(target.relative_to(root)),'bytes':target.stat().st_size,
            'sha256':digest(target),'storage':'uncommitted transcript; annex at next requested checkpoint'}
        (root/'_measurements/S11c_d_adjoint_current_checkpoint.json').write_text(json.dumps(inventory,indent=2)+'\n')
    (base/'manifest.json').write_text(json.dumps(inventory,indent=2)+'\n')
    print(json.dumps({k:inventory[k] for k in ('recordCount','definedFieldMaps','fullSubspaceDirections',
        'tagCount','checkedTensorCount','literalResidualScalars','checkedSymbolicZeroScalars','checkedMetadataPaths','wallSeconds','peakRssKiB')},indent=2))


if __name__=='__main__':main()
