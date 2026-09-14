#!/usr/bin/env python3
"""Validate a computed pairing packet, preserving nonzero diagnostic results."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import pickle
import re
import shutil
import tempfile

import sympy as sp
from S11c_d_joint_sheet_check import decoded_lines, _restore, engine
from S11c_d_modal_current_check import digest
from S11c_d_output_codec import restore_emission_index
from S11c_d_end_pairing_emit import diagnostic_record, worker_record, small_literal

ROOT=Path(__file__).resolve().parents[1]


def assoc(value):return {str(k):v for k,v in value}


def run():
    parser=argparse.ArgumentParser()
    parser.add_argument('--run-directory',type=Path,required=True)
    parser.add_argument('--publish',action='store_true')
    parser.add_argument('--publication-suffix',default='',
                        help='Optional fresh checkpoint suffix; existing publications are never replaced.')
    args=parser.parse_args();base=args.run_directory
    if args.publication_suffix and not re.fullmatch(r'[a-z][a-z0-9_]*',args.publication_suffix):
        raise ValueError('publication suffix must be a lowercase identifier')
    summary=json.loads((base/'checks.json').read_text());emission=json.loads((base/'emission.json').read_text())
    for name,sha in summary['sourceFiles'].items():
        if digest(ROOT/name)!=sha or digest(base/'source'/name)!=sha:raise ValueError(('source pin',name))
    if digest(ROOT/'_measurements/S11c_d_end_pairing_emit.py')!=emission['emitterSha256']:raise ValueError('emitter pin')
    if digest(base/'complete.pickle')!=summary['objectsSha256'] or summary['objectsSha256']!=emission['calculationSha256']:
        raise ValueError('calculation pin')
    if digest(base/'emissions.pickle')!=emission['emissionsSha256']:raise ValueError('emission pin')
    diagnostic=diagnostic_record(base)
    if diagnostic!=emission.get('diagnosticRuntime'):raise ValueError('diagnostic runtime provenance')
    workers=worker_record(base)
    if workers!=emission.get('workerRuntime'):raise ValueError('worker runtime provenance')
    packets=pickle.loads((base/'emissions.pickle').read_bytes())
    records={}
    for line in decoded_lines(base/'full.out'):
        tag,sep,body=line.rstrip('\n').partition(': ')
        if not sep or not tag.startswith('PY_S11CD_') or tag in records:raise ValueError(('transcript tag',tag))
        records[tag]=_restore(body)
    expected={prefix+p['tag'] for p in packets for prefix in ('PY_S11CD_','PY_S11CD_METADATA_')}
    if set(records)!=expected:raise ValueError('transcript packet census')
    index_packet=next(p for p in packets if p['tag'].endswith('_EMISSION_LINES'))
    index_tag='PY_S11CD_'+index_packet['tag']
    indexed=restore_emission_index(assoc(index_packet['body']),list(records)[:list(records).index(index_tag)])
    path_count=0;nonfinite=[];rational_paths=0
    for packet in packets:
        name=packet['tag'];body=packet['body']
        if packet['literal']!=(name.endswith('_RESIDUAL') and small_literal(body)):
            raise ValueError(('literal representation bound',name))
        actual=records['PY_S11CD_'+name]
        expected_body=body if packet['literal'] else engine.carrier_fingerprint(body)
        if actual!=expected_body or packet['fingerprint']!=engine.carrier_fingerprint(body):raise ValueError(('payload fingerprint',name))
        metadata=records['PY_S11CD_METADATA_'+name]
        if metadata!=packet['metadata']:raise ValueError(('metadata payload',name))
        if body.has(sp.nan,sp.zoo,sp.oo,-sp.oo):nonfinite.append(name)
        paths=[]
        for item in metadata:
            item=assoc(item);paths.append(tuple(str(v) if isinstance(v,engine.Str) else int(v) for v in item['OBJECT_PATH']))
            if len(item['VALUE_DIMENSION_L_T_M'])!=3 or any(v.free_symbols for v in item['VALUE_DIMENSION_L_T_M']):
                raise ValueError(('unresolved physical dimensions',name))
            rational_paths+=str(item['GRADE_REPRESENTATION'])=='EXACT_RATIONAL_NUMERATOR_DENOMINATOR'
            if not item['GRADE_DATA']:raise ValueError(('missing grade descriptor',name))
            for descriptor in item['GRADE_DATA']:
                descriptor=assoc(descriptor)
                if len(descriptor['DIMENSION_L_T_M'])!=3 or any(v.free_symbols for v in descriptor['DIMENSION_L_T_M']):
                    raise ValueError('grade dimensions unresolved')
                if not all(k in descriptor for k in ('MULTIGRADE','EPSILON_LAMBDA_SUPPORT','PATHS')):raise ValueError('grade census')
        body_paths=[path for path,v in engine.leaves(body) if not isinstance(v,engine.Str)]
        if Counter(paths)!=Counter(body_paths):raise ValueError(('metadata path census',name))
        path_count+=len(paths)
    if emission['dimensionConstraints'] or nonfinite:raise ValueError(('dimension/nonfinite',emission['dimensionConstraints'],nonfinite))
    computed,_=pickle.loads((base/'complete.pickle').read_bytes())
    scalar_count=nonzero_count=0
    nonzero_entries=[]
    for name,value in computed['retained'].items():
        leaves=list(engine.leaves(engine.cas(value)))
        values=[v for _,v in leaves]
        count=sum(v!=0 for v in values)
        if len(values)!=summary['records'][name]['scalars'] or count!=summary['records'][name]['retainedNonzeroScalars']:
            raise ValueError(('residual count',name))
        scalar_count+=len(values);nonzero_count+=count
        if count:
            source=next(p for p in packets if p['tag'].endswith('_RETAINED_'+name))
            for path,leaf in leaves:
                if leaf==0:continue
                metadata=next(item for item in source['metadata']
                              if tuple(assoc(item)['OBJECT_PATH'])==path)
                nonzero_entries.append({'family':name,'path':path,
                    'value':str(leaf),'factoredValue':str(sp.factor(leaf)),
                    'metadata':sp.srepr(metadata)})
    inventory={**summary,'runDirectory':str(base),'emitterSha256':emission['emitterSha256'],
        'validationSha256':digest(Path(__file__)),'objectCount':len(packets),'tagCount':len(records),
        'metadataPaths':path_count,'rationalMetadataPaths':rational_paths,'nonfiniteObjects':nonfinite,
        'sourceAssignments':len(indexed),
        'diagnosticRuntime':diagnostic,
        'workerRuntime':workers,
        'retainedResidualScalars':scalar_count,'retainedNonzeroScalars':nonzero_count,
        'retainedNonzeroEntries':nonzero_entries,
        'artifacts':{name:{'bytes':(base/name).stat().st_size,'sha256':digest(base/name)} for name in
                     ('full.out','complete.pickle','checks.json','emission.json','emissions.pickle')},
        'scope':'Development case and supplied endpoint limits; independent frequency legs; raw and retained balance records. No end-mode flux normalization.'}
    print(json.dumps({k:inventory[k] for k in ['end','objectCount','tagCount','metadataPaths','retainedResidualScalars','retainedNonzeroScalars']},indent=2))
    if args.publish:
        stem='S11c_d_end_pairing_'+summary['end'].lower()
        if args.publication_suffix:stem+='_'+args.publication_suffix
        target=ROOT/'scripts/out'/(stem+'.out')
        if target.exists() or target.is_symlink():raise ValueError('publication target exists')
        with tempfile.NamedTemporaryFile(dir=target.parent,prefix='.s11cd-end-pairing-',delete=False) as stream:
            temporary=Path(stream.name)
            with (base/'full.out').open('rb') as source:shutil.copyfileobj(source,stream)
            stream.flush();os.fsync(stream.fileno())
        if digest(temporary)!=inventory['artifacts']['full.out']['sha256']:raise ValueError('publication hash')
        os.replace(temporary,target)
        inventory['publication']={'path':str(target.relative_to(ROOT)),**inventory['artifacts']['full.out']}
        (ROOT/'_measurements'/(stem+'_checkpoint.json')).write_text(json.dumps(inventory,indent=2)+'\n')
    (base/'validation.json').write_text(json.dumps(inventory,indent=2)+'\n')


if __name__=='__main__':run()
