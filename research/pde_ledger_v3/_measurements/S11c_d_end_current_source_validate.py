#!/usr/bin/env python3
"""Validate saved end-source operands and publish completed diagnostic packets."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import pickle
import re
import shutil
import tempfile
from types import SimpleNamespace

import sympy as sp

from S11c_d_modal_current_check import digest,engine
from S11c_d_joint_sheet_check import decoded_lines,_restore
from S11c_d_output_codec import restore_emission_index
from S11c_d_current_runtime_metadata import cancellation_units, cancellation_packet


def assoc(value):return {str(k):v for k,v in value}


def validate(base):
    root=Path(engine.__file__).resolve().parents[1]
    checks=json.loads((base/'checks.json').read_text());provenance=checks['provenance']
    for name,key in (('scripts/S11c_d_mixing_scattering_sympy_audit.py','engineSha256'),
                     (provenance.get('instrumentPath','_measurements/S11c_d_end_current_source_check.py'),'instrumentSha256')):
        if digest(root/name)!=provenance[key] or digest(base/'source'/name)!=provenance[key]:
            raise ValueError(('end-source source pin',name))
    if digest(base/'objects.pickle')!=checks['objectsSha256']:raise ValueError('end-source payload pin')
    for name,pin in provenance.get('sourceFiles',{}).items():
        if digest(root/name)!=pin or digest(base/'source'/name)!=pin:raise ValueError(('end-source dependency pin',name))
    progress=[json.loads(line) for line in (base/'progress.jsonl').read_text().splitlines()]
    if progress[-1]['stage']!='completed':raise ValueError('end-source incomplete emission')
    data=pickle.loads((base/'objects.pickle').read_bytes())
    symbols={s.name:s for s in data['knownDimensions'] if isinstance(s,sp.Symbol)}
    r=SimpleNamespace(symbols=symbols,ell=symbols['L_W'])
    dimensions=engine.DimensionAnalysis.__new__(engine.DimensionAnalysis)
    dimensions.known,dimensions.unknown,dimensions.constraints,dimensions.solution=dict(data['knownDimensions']),{},set(),{}
    dimensions.zero=(sp.S.Zero,)*3
    engine.PHYSICAL_METADATA=engine.PhysicalMetadata(dimensions,r)
    records={};tags=[]
    for line in decoded_lines(base/'full.out'):
        tag,sep,payload=line.rstrip('\n').partition(': ')
        if not sep or not tag.startswith('PY_S11CD_') or tag in records:raise ValueError(('transcript tag',tag))
        records[tag]=_restore(payload);tags.append(tag)
    prefix='END_CURRENT_SOURCE_'+provenance['end']+'_LAB_HELD_RHO4_CONSTANT'
    index_tag='PY_S11CD_'+prefix+'_EMISSION_LINES'
    index=restore_emission_index(assoc(records[index_tag]),tags[:tags.index(index_tag)])
    for tag in tags:
        if not tag.startswith('PY_S11CD_METADATA_') and 'PY_S11CD_METADATA_'+tag.removeprefix('PY_S11CD_') not in records:
            raise ValueError(('missing metadata',tag))
    if any('ZERO_MAP' in sp.srepr(v) for tag,v in records.items() if tag.startswith('PY_S11CD_METADATA_')):
        raise ValueError('unassigned emitted dimension')
    if records['PY_S11CD_'+prefix+'_DIMENSION_CONSTRAINTS']:raise ValueError('dimension constraints')
    paths_checked=0;objects_checked=0;residuals=0
    for group,tag_prefix in (('conservative','CONSERVATIVE'),('slab','ENERGY_BALANCE'),('acoustic','CLOSED_ACOUSTIC')):
        for name,value in data[group].items():
            tag='PY_S11CD_'+tag_prefix+'_'+name+'_'+prefix
            body=engine.cas(value)
            expected=engine.carrier_fingerprint(body) if engine.dag_size(body)>1200 else body
            if records[tag]!=expected:raise ValueError(('source operand serialization',tag))
            metadata_body=value.lhs-value.rhs if name=='DEPTH_CONVERGENCE_DOMAIN' else body
            leaves=dict(engine.leaves(engine.cas(metadata_body)))
            entries=records['PY_S11CD_METADATA_'+tag.removeprefix('PY_S11CD_')]
            covered=[]
            for entry in entries:
                if len(entry)==2:
                    raw_path,raw_descriptor=entry;descriptor=assoc(raw_descriptor);paths=(raw_path,)
                else:
                    descriptor=assoc(entry);paths=descriptor['PATHS']
                if len(descriptor['DIMENSION_L_T_M'])!=3:raise ValueError('unit dimension arity')
                for raw_path in paths:
                    path=tuple(str(v) if isinstance(v,engine.Str) else int(v) for v in raw_path)
                    expression=leaves[path]
                    coefficients=engine.PHYSICAL_METADATA.coefficients(expression)
                    grades=tuple(tuple(int(v) for v in row) for row in descriptor['MULTIGRADE'])
                    if grades!=tuple(sorted(coefficients)):raise ValueError(('source multigrade',tag,path))
                    homotopy={}
                    for (e,a,b),c in coefficients.items():
                        homotopy[e,a+b]=homotopy.get((e,a+b),sp.S.Zero)+c*(symbols['W_0']/r.ell)**b
                    expected_support=tuple(sorted(g for g,c in homotopy.items() if c!=0))
                    support=tuple(tuple(int(v) for v in row) for row in descriptor['EPSILON_LAMBDA_SUPPORT'])
                    if support!=expected_support:raise ValueError(('source lambda support',tag,path))
                    covered.append(path);paths_checked+=1
            if Counter(covered)!=Counter(path for path,v in leaves.items() if not isinstance(v,engine.Str)):
                raise ValueError(('metadata path census',tag))
            objects_checked+=1
            if name.endswith('_RESIDUAL'):
                count=sum(not isinstance(v,engine.Str) for v in leaves.values())
                if count!=checks['residualInventory'][group+'_'+name]['scalars']:raise ValueError('residual census')
                residuals+=count
    proof_count=0
    if 'proofEmissionsSha256' in checks:
        proof_path=base/'proof-emissions.pickle'
        if digest(proof_path)!=checks['proofEmissionsSha256']:raise ValueError('cancellation proof payload pin')
        modes=engine.FullPencilModes.__new__(engine.FullPencilModes);modes.r=r
        restored_units=cancellation_units(engine,data['cancellationProof'],data['strong'],dimensions)
        for name,record in data['cancellationProof'].items():
            for key,value in record.items():
                if not key.endswith('_DOMAIN') and value!=0 and dimensions.measure(value)!=restored_units[name][key]:
                    raise ValueError(('cancellation operand dimension',name,key))
        proof_units={'CANCELLATION_'+name+'_'+key:unit for name,record in restored_units.items() for key,unit in record.items()}
        proof_keys={'CANCELLATION_'+name+'_'+key:(name,key) for name,record in data['cancellationProof'].items() for key in record}
        for item in pickle.loads(proof_path.read_bytes()):
            tag='PY_S11CD_'+prefix+'_'+item['name'];body=item['body']
            name,key=proof_keys[item['name']]
            expected_item=cancellation_packet(engine,modes,name,key,data['cancellationProof'][name],proof_units[item['name']])
            if item!=expected_item:raise ValueError(('cancellation restored packet',tag))
            expected=engine.carrier_fingerprint(body) if item['heavy'] else body
            if records[tag]!=expected:raise ValueError(('cancellation proof serialization',tag))
            metadata=expected_item['metadata']
            if records['PY_S11CD_METADATA_'+tag.removeprefix('PY_S11CD_')]!=metadata:
                raise ValueError(('cancellation proof metadata',tag))
            proof_count+=1
        if dimensions.constraints:raise ValueError('cancellation dimension constraints')
        if sum(v['CROSS_PRODUCT_RESIDUAL']!=0 for v in data['cancellationProof'].values())!=checks['nonzeroCancellationIdentities']:
            raise ValueError('cancellation residual census')
    return {**checks,'checkedCancellationObjects':proof_count,'runDirectory':str(base),'tagCount':len(tags),'sourceAssignments':len(index),
        'checkedSourceObjects':objects_checked,'checkedMetadataPaths':paths_checked,
        'literalSourceResidualScalars':residuals,'validationSourceSha256':digest(Path(__file__)),
        'artifacts':{name:{'bytes':(base/name).stat().st_size,'sha256':digest(base/name)}
                     for name in ('full.out','objects.pickle','checks.json','progress.jsonl','stderr.txt')}}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True)
    parser.add_argument('--publish',action='store_true')
    parser.add_argument('--publication-suffix',default='')
    parser.add_argument('--baseline-run',type=Path)
    args=parser.parse_args()
    inventory=validate(args.run_directory)
    if args.baseline_run:
        old_checks=json.loads((args.baseline_run/'checks.json').read_text())
        old_path=args.baseline_run/'objects.pickle'
        if digest(old_path)!=old_checks['objectsSha256']:raise ValueError('baseline payload pin')
        for key in ('end','inputSha256'):
            if inventory['provenance'].get(key)!=old_checks['provenance'].get(key):raise ValueError(('baseline input',key))
        old=pickle.loads(old_path.read_bytes());new=pickle.loads((args.run_directory/'objects.pickle').read_bytes())
        from S11c_d_end_current_resumable import polynomial_identity
        compared=0;exact_objects=0;algebraic_leaves=0;mapping_entries=0;nonzero=[]
        for group in ('conservative','slab','acoustic'):
            if old[group].keys()!=new[group].keys():raise ValueError('baseline object census')
            for name,before in old[group].items():
                after=new[group][name];compared+=1
                if engine.cas(before)==engine.cas(after):exact_objects+=1;continue
                # This native object is parameter_map.items(), whose order
                # follows an atom set. Join the complete map by symbolic key.
                if name=='SOURCE_PARAMETER_ALIGNMENT':
                    a,b=dict(before),dict(after)
                    if len(a)!=len(before) or len(b)!=len(after) or a.keys()!=b.keys():
                        raise ValueError('baseline parameter-map key census')
                    for key,value in a.items():
                        residual,_=polynomial_identity(value,b[key]);mapping_entries+=1
                        if residual!=0:nonzero.append((group,name,(sp.srepr(key),),sp.srepr(residual)))
                    continue
                a=dict(engine.leaves(engine.cas(before)));b=dict(engine.leaves(engine.cas(after)))
                if a.keys()!=b.keys():raise ValueError(('baseline leaf census',group,name))
                for path,value in a.items():
                    if value==b[path]:continue
                    residual,_=polynomial_identity(value,b[path]);algebraic_leaves+=1
                    if residual!=0:nonzero.append((group,name,path,sp.srepr(residual)))
        inventory['baselineRegression']={'runDirectory':str(args.baseline_run),'objectsSha256':digest(old_path),
            'comparedObjects':compared,'structurallyIdenticalObjects':exact_objects,
            'algebraicallyComparedLeaves':algebraic_leaves,'parameterBindingsComparedByKey':mapping_entries,
            'nonzeroResiduals':nonzero}
        (args.run_directory/'baseline-regression.json').write_text(json.dumps(inventory['baselineRegression'],indent=2)+'\n')
        if nonzero:raise ValueError('baseline reconstruction residual; inspect saved comparison')
    if args.publish:
        root=Path(engine.__file__).resolve().parents[1];end=inventory['provenance']['end'].lower()
        if args.publication_suffix and not re.fullmatch('[a-z0-9_]+',args.publication_suffix):raise ValueError('publication suffix')
        suffix='_'+args.publication_suffix if args.publication_suffix else ''
        stem=f'S11c_d_end_current_source_{end}{suffix}'
        target=root/'scripts/out'/(stem+'.out')
        if target.exists() or target.is_symlink():raise ValueError('publication already exists')
        with tempfile.NamedTemporaryFile(dir=target.parent,prefix='.s11cd-end-current-',delete=False) as stream:
            temporary=Path(stream.name)
            with (args.run_directory/'full.out').open('rb') as source:shutil.copyfileobj(source,stream)
            stream.flush();os.fsync(stream.fileno())
        if digest(temporary)!=inventory['artifacts']['full.out']['sha256']:raise ValueError('publication hash')
        os.replace(temporary,target)
        inventory['publication']={'path':str(target.relative_to(root)),'sha256':digest(target),'bytes':target.stat().st_size}
        (root/'_measurements'/(stem+'_checkpoint.json')).write_text(json.dumps(inventory,indent=2)+'\n')
    (args.run_directory/'manifest.json').write_text(json.dumps(inventory,indent=2)+'\n')
    print(json.dumps({k:inventory[k] for k in ('tagCount','checkedSourceObjects','checkedMetadataPaths',
        'literalSourceResidualScalars','retainedNonzeroScalars','wallSeconds','peakRssKiB')},indent=2))


if __name__=='__main__':main()
