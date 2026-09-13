#!/usr/bin/env python3
"""Validate the saved end-frequency objects and atomically publish their transcript."""
import argparse
import json
import os
from pathlib import Path
import pickle
import tempfile
from types import SimpleNamespace
import numpy as np
import sympy as sp
from S11c_d_end_frequency_check import ROOT, load_packet, digest, engine, decoded_lines, _restore, check_certificate


def validate(base):
    summary=json.loads((base/'checks.json').read_text())
    for name,sha in summary['sourceFiles'].items():
        if digest(base/'source'/name)!=sha or digest(ROOT/name)!=sha:raise ValueError(('source pin',name))
    objects=base/'objects.pickle'
    if digest(objects)!=summary['objectsSha256']:raise ValueError('payload pin')
    args=SimpleNamespace(**{k:Path(v) if k!='end' else v for k,v in summary['arguments'].items()})
    modes,strong,units,inputs,provenance,native,coverage,packet=load_packet(args)
    result=pickle.loads(objects.read_bytes())
    if result['NATIVE_RECORDS']!=native or result['NATIVE_COVERAGE']!=coverage:raise ValueError('native packet cache join')
    engine.PHYSICAL_METADATA.dimensions.known.update(result['KNOWN_DIMENSIONS'])
    entries={}
    for line in decoded_lines(base/'full.out'):
        tag,sep,body=line.partition(': ')
        if not sep or tag in entries:raise ValueError('transcript syntax/duplicate tag')
        entries[tag]=_restore(body)
    expected=set();paths=0;residual_scalars=0
    for p in result['PACKETS']:
        tag='PY_S11CD_'+summary['prefix']+'_'+p['NAME'];meta='PY_S11CD_METADATA_'+summary['prefix']+'_'+p['NAME']
        expected.update((tag,meta));body=p['VALUE']
        fingerprint=engine.carrier_fingerprint(body) if p['HEAVY']=='carrier' else modes.compact_fingerprint(body) if p['HEAVY'] else body
        if entries.get(tag)!=fingerprint:raise ValueError(('object serialization',tag))
        metadata=modes.numeric_metadata(body,lambda path:p['UNITS'][path])
        if entries.get(meta)!=metadata:raise ValueError(('metadata serialization',tag))
        leaves={path for path,v in engine.leaves(body) if not isinstance(v,engine.Str)}
        emitted=[]
        for group in metadata:
            item={str(k):v for k,v in group};emitted.extend(tuple(str(x) if isinstance(x,engine.Str) else x for x in v) for v in item['PATHS'])
            if len(item['DIMENSION_L_T_M'])!=3:raise ValueError('dimension triple')
        if set(emitted)!=leaves or len(emitted)!=len(leaves) or set(p['UNITS'])!=leaves:
            raise ValueError(('metadata path census',tag))
        paths+=len(leaves)
        if p['NAME'].endswith('_RESIDUAL'):residual_scalars+=len(leaves)
    final='PY_S11CD_'+summary['prefix']+'_EMISSION_LINES';meta='PY_S11CD_METADATA_'+summary['prefix']+'_EMISSION_LINES'
    if set(entries)!=expected|{final,meta}:raise ValueError('final tag census')
    indexed=entries[final]
    from S11c_d_output_codec import restore_emission_index
    decoded_index={str(k):v for k,v in indexed}
    source_lines=restore_emission_index(decoded_index,[tag for tag in entries if tag in expected])
    if set(source_lines)!=expected:raise ValueError('emission line tag set')
    if entries[meta]!=modes.numeric_metadata(indexed,lambda p:(0,0,0)):raise ValueError('emission index metadata')
    reality=check_certificate(result['NORMAL_REALITY_COVERAGE'],result['RECORDS'])
    byname={p['NAME']:p['VALUE'] for p in result['PACKETS']}
    checks=[]
    for mode in result['RECORDS']:
        i=mode['INDEX'];n=mode.get('NULLITY',0);stem='MODE_'+str(i)+'_'
        if mode['NATIVE_NULLITY_RESIDUAL']!=0:raise ValueError('native subspace rank join')
        if not mode['FREQUENCY_NORMALIZATION_DEFINED']:continue
        array=lambda name:np.array(byname[stem+name],dtype=complex)
        r,l,dw,norm,p=array('RIGHT_BASIS'),array('LEFT_BASIS'),array('FREQUENCY_MATRIX'),array('FREQUENCY_NORMALIZED_LEFT'),array('FREQUENCY_FIELD_PROJECTOR')
        pairing=array('FREQUENCY_PAIRING')
        if np.linalg.matrix_rank(r,tol=1e-9)!=n or np.linalg.matrix_rank(l,tol=1e-9)!=n or np.linalg.matrix_rank(p,tol=1e-9)!=n:
            raise ValueError('incomplete subspace/projector rank')
        checks.append({'index':i,'nullity':n,'pairingReconstruction':float(np.linalg.norm(pairing-l.conj().T@dw@r)),
            'projectorReconstruction':float(np.linalg.norm(p-r@np.linalg.solve(pairing,l.conj().T@dw))),
            'normalization':float(np.linalg.norm(norm.conj().T@dw@r-np.eye(n)))})
    return {'runDirectory':str(base.resolve()),'summary':summary,'checks':checks,'normalReality':reality,
        'tagCount':len(entries),'objectCount':len(result['PACKETS']),'metadataPaths':paths,
        'residualScalars':residual_scalars,'sourceFiles':summary['sourceFiles'],
        'artifacts':{name:{'sha256':digest(base/name),'bytes':(base/name).stat().st_size}
                     for name in ('full.out','stderr.txt','objects.pickle','checks.json','progress.jsonl')}}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True)
    parser.add_argument('--publish',action='store_true');args=parser.parse_args()
    inventory=validate(args.run_directory)
    if args.publish:
        end=inventory['summary']['arguments']['end'].lower()
        destination=ROOT/'scripts/out'/('S11c_d_end_frequency_'+end+'.out')
        if destination.exists() or destination.is_symlink():raise ValueError('publication destination exists')
        handle,temp=tempfile.mkstemp(prefix='.s11cd-end-frequency-',dir=destination.parent)
        try:
            with os.fdopen(handle,'wb') as stream:stream.write((args.run_directory/'full.out').read_bytes())
            if digest(Path(temp))!=inventory['artifacts']['full.out']['sha256']:raise ValueError('publication copy digest')
            os.replace(temp,destination)
        finally:
            if Path(temp).exists():Path(temp).unlink()
        inventory['publishedTranscript']=str(destination.relative_to(ROOT))
        target=ROOT/'_measurements'/('S11c_d_end_frequency_'+end+'_checkpoint.json')
    else:target=args.run_directory/'validation.json'
    target.write_text(json.dumps(inventory,indent=2)+'\n')
    print(json.dumps({'inventory':str(target),'objects':inventory['objectCount'],'tags':inventory['tagCount'],
        'metadataPaths':inventory['metadataPaths'],'residualScalars':inventory['residualScalars']}))


if __name__=='__main__':main()
