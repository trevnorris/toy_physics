#!/usr/bin/env python3
"""Source-pinned one-case energy-current construction from reduced operands."""
import argparse
import hashlib
import inspect
import json
from pathlib import Path
import pickle
import resource
import time

import sympy as sp

from S11c_d_joint_sheet_check import load, engine


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    started = time.monotonic()
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--end', choices=('REFERENCE', 'LEFT', 'RIGHT'), default='REFERENCE')
    parser.add_argument('--cache-result', type=Path)
    parser.add_argument('--acoustic', action='store_true')
    parser.add_argument('--slab-cache', type=Path)
    parser.add_argument('--require-zero-residuals', action='store_true',
                        help='guard after emitting every comparison; exit zero alone means emission completed')
    args = parser.parse_args()
    modes, strong, units, bindings, inputs, provenance = load(args)
    manifest = json.loads(args.manifest.read_text())
    cache = Path(manifest['run_directory'])/'symbols'/f'{args.end}_LAB_HELD_RHO4_CONSTANT.pickle'
    energy = pickle.loads(cache.read_bytes())[5]
    r = modes.r
    r.x = tuple(r.symbols['s11cc2X'+str(i)] for i in (1, 2, 3))
    r.t = r.symbols['s11cc2Time']
    r.z = r.symbols['s11cdNormalPosition']
    ends = engine.ConstantEndPencil.__new__(engine.ConstantEndPencil)
    ends.r, ends.kn = r, modes.k
    current = engine.UniformSlabCurrent(r, {'value':energy}, ends, strong[3, :])
    audit = engine.SlabEnergyBalance(current)
    end = {'REFERENCE':None, 'LEFT':-sp.oo, 'RIGHT':sp.oo}[args.end]
    cache_key = hashlib.sha256((digest(cache)+sp.__version__+args.end+
        inspect.getsource(engine.UniformSlabCurrent)+inspect.getsource(engine.SlabEnergyBalance)+
        inspect.getsource(engine.ConstantEndPencil.phase_terms)).encode()).hexdigest()
    if args.slab_cache and args.slab_cache.exists():
        saved = pickle.loads(args.slab_cache.read_bytes())
        if saved['key'] != cache_key:
            raise ValueError('development slab cache source/input mismatch')
        current.construct = lambda anchoring, endpoint:saved['conservative']
        audit.construct = lambda anchoring, endpoint:saved['balance']
    elif args.slab_cache:
        args.slab_cache.write_bytes(pickle.dumps({'key':cache_key,
            'balance':audit.construct('LAB_HELD', end),
            'conservative':current.construct('LAB_HELD', end)}, protocol=5))
    provenance['instrumentSha256'] = digest(Path(__file__))
    provenance['slabConstructionCacheKey'] = cache_key
    engine.physical('NONLOCAL_CURRENT_PREFLIGHT_PROVENANCE', provenance)
    result = audit.emit('LAB_HELD', end, args.end+'_LAB_HELD_RHO4_CONSTANT')
    acoustic = None
    if args.acoustic:
        acoustic_builder = engine.ClosedAcousticEnergy(audit, modes, strong)
        acoustic = acoustic_builder.construct('LAB_HELD', end)
    if args.cache_result:
        args.cache_result.write_bytes(pickle.dumps((result, current.construct('LAB_HELD', end), acoustic), protocol=5))
    if acoustic:
        acoustic_builder.emit('LAB_HELD', end, args.end+'_LAB_HELD_RHO4_CONSTANT')
    dims = engine.PHYSICAL_METADATA.dimensions
    engine.physical('NONLOCAL_CURRENT_PREFLIGHT_DIMENSION_CONSTRAINTS', tuple(dims.constraints))
    resources = {'wallSeconds':time.monotonic()-started,
                 'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    engine.emit('NONLOCAL_CURRENT_PREFLIGHT_RESOURCES', resources)
    engine.emit('METADATA_NONLOCAL_CURRENT_PREFLIGHT_RESOURCES', modes.numeric_metadata(
        engine.cas(resources), lambda p:(0, 1, 0) if p[0] == 'wallSeconds' else dims.zero))
    residuals = {key:value for key, value in result.items() if key.endswith('_RESIDUAL')}
    if acoustic:
        residuals.update({'ACOUSTIC_'+key:value for key, value in acoustic.items() if key.endswith('_RESIDUAL')})
        residuals.update({'FACE_'+str(i):face['CLOSURE_RESIDUAL'] for i, face in enumerate(acoustic['FACE_RECORDS'])})
    engine.physical('NONLOCAL_CURRENT_PREFLIGHT_RESIDUAL_CENSUS',
                    {key:sum(v != 0 for _, v in engine.leaves(engine.cas(value))) for key, value in residuals.items()},
                    zero_dimensions={(key,):dims.zero for key in residuals})
    if args.require_zero_residuals and any(v != 0 for value in residuals.values() for _, v in engine.leaves(engine.cas(value))):
        raise ValueError('energy-balance construction residual; see emitted operands')


if __name__ == '__main__':
    run()
