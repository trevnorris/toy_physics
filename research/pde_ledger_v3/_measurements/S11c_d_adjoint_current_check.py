#!/usr/bin/env python3
"""Pinned reference continuation of the physical-field / row-dual current map."""
import argparse
import ast
import json
from pathlib import Path
import pickle
import resource
import time
import numpy as np

from S11c_d_modal_subspace_check import load_current, engine
from S11c_d_modal_current_check import digest


def load_modal(args,provenance):
    checkpoint = json.loads(args.modal_checkpoint.read_text())
    base = Path(checkpoint['runDirectory'])
    frozen = base/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py'
    if digest(frozen)!=checkpoint['sourceFiles']['scripts/S11c_d_mixing_scattering_sympy_audit.py']:
        raise ValueError('modal producer snapshot pin')
    def existing_source(path):
        module = ast.parse(path.read_text())
        module.body = [node for node in module.body if not isinstance(node,ast.ClassDef) or node.name!='AdjointCurrentMap']
        return ast.dump(module,include_attributes=False)
    if existing_source(frozen)!=existing_source(Path(engine.__file__)):
        raise ValueError('modal producer changed outside the adjoint continuation')
    for key in ('producerManifestSha256','cacheSha256','currentCacheSha256','currentManifestSha256',
                'inputSha256','currentCheckpointSha256','pairingCacheSha256'):
        if checkpoint['provenance'][key]!=provenance[key]:raise ValueError(('modal input/proof mismatch',key))
    payload = base/'objects.pickle'
    if digest(payload)!=checkpoint['artifacts']['objects.pickle']['sha256']:
        raise ValueError('modal cache payload pin')
    result,known = pickle.loads(payload.read_bytes())
    engine.PHYSICAL_METADATA.dimensions.known.update(known)
    if len(result['RECORDS'])!=checkpoint['recordCount']:
        raise ValueError('modal candidate count mismatch')
    provenance.update({'modalCheckpointSha256':digest(args.modal_checkpoint),
        'modalCacheSha256':digest(payload),'modalProducerEngineSha256':digest(frozen),
        'rootLiftCount':len(result['RECORDS'])})
    return result


def main():
    parser = argparse.ArgumentParser()
    for name in ('manifest','input','current-manifest','current-checkpoint','modal-checkpoint','run-directory'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--emit',action='store_true')
    args = parser.parse_args();args.end='REFERENCE'
    args.run_directory.mkdir(parents=True,exist_ok=True)
    started = time.monotonic()
    def progress(record):
        record['elapsedSeconds'] = time.monotonic()-started
        with (args.run_directory/'progress.jsonl').open('a') as stream:stream.write(json.dumps(record)+'\n')
        if not args.emit:print(json.dumps(record),flush=True)
    pairing,current,inputs,provenance = load_current(args)
    source = load_modal(args,provenance)
    provenance['instrumentSha256'] = digest(Path(__file__))
    m = pairing.modes
    algebraic,relation,_ = m.analytic(pairing.acoustic.strong)
    binding = inputs.mapping(algebraic,relation,(m.k,m.q,m.eta,m.sigma))
    binding.update({m.eta:0,m.sigma:0})
    modal = engine.ModalCurrentSubspaces(pairing,current,binding)
    builder = engine.AdjointCurrentMap(modal)
    progress({'stage':'loaded','provenance':provenance})
    computed = builder.construct(source,progress)
    payload = args.run_directory/'objects.pickle'
    payload.write_bytes(pickle.dumps((computed,engine.PHYSICAL_METADATA.dimensions.known),protocol=5))
    progress({'stage':'objects_saved','sha256':digest(payload)})
    if args.emit:builder.emit(computed,provenance)
    residuals = {}
    for record in computed['RECORDS']:
        for item in record['ITEMS']:
            if item['GROUP']!='RESIDUALS':continue
            entry = residuals.setdefault(item['NAME'],{'matrices':0,'scalars':0,'maximumNorm':0.})
            entry['matrices']+=1;entry['scalars']+=item['VALUE'].size
            entry['maximumNorm']=max(entry['maximumNorm'],float(np.linalg.norm(item['VALUE'])))
    summary = {'provenance':provenance,'recordCount':len(computed['RECORDS']),
        'definedFieldMaps':sum(r['INVERTIBLE_FIELD_MAP_DEFINED'] for r in computed['RECORDS']),
        'fullSubspaceDirections':sum(r['NULLITY'] for r in computed['RECORDS']),
        'sourceNormalizedSubspaces':sum(r['SOURCE_CURRENT_NORMALIZATION_DEFINED'] for r in computed['RECORDS']),
        'powerMapRanks':dict(engine.Counter(str(r['POWER_MAP_RANK']) for r in computed['RECORDS'])),
        'symbolicResidualScalars':sum(len(v) for v in computed['SYMBOLIC_RESIDUALS'].values()),
        'symbolicNonzeroResidualScalars':sum(v!=0 for matrix in computed['SYMBOLIC_RESIDUALS'].values() for v in matrix),
        'numericalResiduals':residuals,'objectsSha256':digest(payload),
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (args.run_directory/'checks.json').write_text(json.dumps(summary,indent=2)+'\n')
    if args.emit:
        index = engine.emission_index(engine.EMISSION_LINES)
        engine.emit('ADJOINT_CURRENT_MAP_EMISSION_LINES',index)
        engine.emit('METADATA_ADJOINT_CURRENT_MAP_EMISSION_LINES',m.numeric_metadata(engine.cas(index),lambda path:(0,0,0)))
    progress({'stage':'completed','summary':summary})
    if summary['symbolicNonzeroResidualScalars'] or any(r['maximumNorm']>1e-8 for r in residuals.values()):
        raise ValueError('adjoint-map reconstruction residual; see emitted operands')


if __name__=='__main__':main()
