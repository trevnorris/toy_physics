#!/usr/bin/env python3
"""Pinned native-root full-subspace current continuation."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import pickle
import resource
import time

import sympy as sp

from S11c_d_modal_current_check import build, digest, source_node, engine
from S11c_d_joint_sheet_check import _restore, decoded_lines


def load_current(args):
    pairing, inputs, provenance = build(args)
    checkpoint = json.loads(args.current_checkpoint.read_text())
    base = Path(checkpoint['runDirectory'])
    frozen = base/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py'
    if digest(frozen) != checkpoint['sourceFiles']['scripts/S11c_d_mixing_scattering_sympy_audit.py']:
        raise ValueError('current snapshot source mismatch')
    if source_node(frozen.read_text(), 'ClosedCurrentPairing') != source_node(Path(engine.__file__).read_text(), 'ClosedCurrentPairing'):
        raise ValueError('current constructor/proof source changed')
    for key in ('producerManifestSha256','cacheSha256','currentCacheSha256','currentManifestSha256','inputSha256'):
        if checkpoint['provenance'][key] != provenance[key]:
            raise ValueError(('current checkpoint input mismatch',key))
    objects = base/'objects.pickle'
    if digest(objects) != checkpoint['artifacts']['objects.pickle']['sha256']:
        raise ValueError('current checkpoint payload mismatch')
    result, known = pickle.loads(objects.read_bytes())
    engine.PHYSICAL_METADATA.dimensions.known.update(known)
    pairing.construct = lambda anchoring,end:result
    provenance.update({'currentCheckpointSha256':digest(args.current_checkpoint),
                       'pairingCacheSha256':digest(objects),
                       'instrumentSha256':digest(Path(__file__))})
    return pairing,result,inputs,provenance


def load_spectrum(args,pairing,inputs,provenance):
    manifest = json.loads(args.spectrum_manifest.read_text())
    base = Path(manifest['run_directory'])
    if manifest.get('exit_code') != 0 or manifest['source_hashes_before'] != manifest['source_hashes_after']:
        raise ValueError('completed stable native spectrum producer required')
    for name in ('full.out','symbols/REFERENCE_LAB_HELD_RHO4_CONSTANT.pickle'):
        if digest(base/name) != manifest['artifacts'][name]['sha256']:
            raise ValueError(('native spectrum artifact mismatch',name))
    for name,sha in manifest['source_hashes_after'].items():
        source = base/'source'/name if name.startswith('_measurements/') or name=='scripts/S11c_d_mixing_scattering_sympy_audit.py' else Path(engine.__file__).resolve().parents[1]/name
        if digest(source)!=sha:
            raise ValueError(('native spectrum source mismatch',name))
    if manifest['source_hashes_after']['_measurements/S11c_d_channel_preflight_input.json']!=digest(args.input):
        raise ValueError('native spectrum input changed')
    producer_strong = pickle.loads((base/'symbols/REFERENCE_LAB_HELD_RHO4_CONSTANT.pickle').read_bytes())[4]
    if producer_strong != pairing.acoustic.strong:
        raise ValueError('native spectrum/current physical pencil differs')
    prefix='PY_S11CD_END_SPECTRUM_INPUT_REFERENCE_LAB_HELD_RHO4_CONSTANT_0'
    records=[]; coverage=None; bound=None
    for line in decoded_lines(base/'full.out'):
        tag,_,payload=line.partition(': ')
        if tag.startswith(prefix+'_MODE_') and tag.endswith('_RECORD'):
            records.append({str(k):v for k,v in _restore(payload)})
        elif tag==prefix+'_ROOT_COVERAGE':
            coverage={str(k):v for k,v in _restore(payload)}
        elif tag==prefix+'_BOUND_CARRIERS':
            bound={str(k):v for k,v in _restore(payload)}
    if coverage is None or bound is None or not records:
        raise ValueError('missing native spectrum packet')
    if coverage['FINITE_POLYNOMIAL_ROOT_COVERAGE']!=sp.true:
        raise ValueError('native spectrum isolation unresolved')
    lifts={(int(r['ROOT_DISK_INDEX']),int(r['NORMAL_LIFT_SIGN'])) for r in records}
    expected={(i,sign) for i in range(int(coverage['DISTINCT_ROOT_COUNT'])) for sign in (-1,1)}
    if lifts!=expected or len(records)!=len(lifts):
        raise ValueError('incomplete or duplicate root/lift packet')
    for name,value in bound.items():
        if inputs.parameters.get(name)!=value:
            raise ValueError(('native packet material binding mismatch',name))
    provenance.update({'spectrumManifestSha256':digest(args.spectrum_manifest),
                       'spectrumTranscriptSha256':manifest['artifacts']['full.out']['sha256'],
                       'spectrumCacheSha256':manifest['artifacts']['symbols/REFERENCE_LAB_HELD_RHO4_CONSTANT.pickle']['sha256'],
                       'spectrumPacket':prefix,'rootDiskCount':int(coverage['DISTINCT_ROOT_COUNT']),
                       'rootLiftCount':len(records),'physicalPencilStructuralJoin':True})
    return records,coverage


def run():
    parser=argparse.ArgumentParser()
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--current-manifest',type=Path,required=True)
    parser.add_argument('--current-checkpoint',type=Path,required=True)
    parser.add_argument('--spectrum-manifest',type=Path,required=True)
    parser.add_argument('--run-directory',type=Path,required=True)
    parser.add_argument('--emit',action='store_true')
    args=parser.parse_args();args.end='REFERENCE'
    args.run_directory.mkdir(parents=True,exist_ok=True)
    started=time.monotonic()
    def progress(record):
        record['elapsedSeconds']=time.monotonic()-started
        if args.emit:
            with (args.run_directory/'progress.jsonl').open('a') as stream:stream.write(json.dumps(record)+'\n')
        else:print(json.dumps(record),flush=True)
    pairing,result,inputs,provenance=load_current(args)
    records,coverage=load_spectrum(args,pairing,inputs,provenance)
    modes=pairing.modes
    algebraic,relation,_=modes.analytic(pairing.acoustic.strong)
    mapping=inputs.mapping(algebraic,relation,(modes.k,modes.q,modes.eta,modes.sigma))
    mapping.update({modes.eta:0,modes.sigma:0})
    builder=engine.ModalCurrentSubspaces(pairing,result,mapping)
    progress({'stage':'loaded','provenance':provenance})
    computed=builder.construct(records,coverage,float(inputs.parameters['omega']),float(inputs.parameters['W_0']),progress)
    args.run_directory.mkdir(parents=True,exist_ok=True)
    objects=args.run_directory/'objects.pickle'
    objects.write_bytes(pickle.dumps((computed,engine.PHYSICAL_METADATA.dimensions.known),protocol=5))
    if args.emit:builder.emit(computed,provenance)
    summary={'provenance':provenance,'recordCount':len(computed['RECORDS']),
        'nullityCounts':{str(n):sum(r['NULLITY']==n for r in computed['RECORDS']) for n in sorted({r['NULLITY'] for r in computed['RECORDS']})},
        'physicalCurrentNormalizationCount':sum(r.get('PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED',False) for r in computed['RECORDS']),
        'residualNormMaxima':{name:max(r.get('RESIDUAL_NORMS',{}).get(name,0) for r in computed['RECORDS'])
            for name in sorted(set().union(*(r.get('RESIDUAL_NORMS',{}).keys() for r in computed['RECORDS'])))},
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'objectsSha256':digest(objects),'cutoffBinding':{'source':'W_0','value':str(inputs.parameters['W_0']),
            'role':'finite-depth balance diagnostic; infinite-depth current separately conditioned'}}
    (args.run_directory/'checks.json').write_text(json.dumps(summary,indent=2)+'\n')
    if args.emit:
        index=engine.emission_index(engine.EMISSION_LINES)
        engine.emit('MODAL_SUBSPACE_EMISSION_LINES',index)
        engine.emit('METADATA_MODAL_SUBSPACE_EMISSION_LINES',modes.numeric_metadata(engine.cas(index),lambda path:(0,0,0)))
    progress({'stage':'completed','summary':summary})


if __name__=='__main__':run()
