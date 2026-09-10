#!/usr/bin/env python3
"""Bounded diagnostic of a cached S11c-d mode's radical continuation.

Consumes a dev-symbol-cache and its explicit-input transcript. It prints
computed operands and residuals; it does not repair or classify the spectrum.
"""
import argparse
import cmath
import hashlib
import json
import pickle
from pathlib import Path
from types import SimpleNamespace

import sympy as sp
from sympy.core.symbol import Str
import S11c_d_mixing_scattering_sympy_audit as engine


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--symbol-cache', type=Path, required=True)
    parser.add_argument('--transcript', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--candidate-index', type=int, required=True)
    args = parser.parse_args()
    full, curl, units, known = pickle.loads(args.symbol_cache.read_bytes())[:4]
    symbols = {a.name:a for a in known if isinstance(a,sp.Symbol)}
    symbols.update({a.name:a for a in full.free_symbols})
    dimensions = engine.DimensionAnalysis.__new__(engine.DimensionAnalysis)
    dimensions.known = dict(known)
    dimensions.unknown, dimensions.constraints, dimensions.solution = {},set(),{}
    dimensions.zero = (sp.S.Zero,)*3
    r = SimpleNamespace(symbols=symbols,omega=symbols['omega'],ell=symbols['L_W'])
    engine.PHYSICAL_METADATA = engine.PhysicalMetadata(dimensions,r)
    ends = SimpleNamespace(r=r,kn=symbols['s11cdSpectralNormalMomentum'])
    modes = engine.FullPencilModes(ends,curl,units)
    _,relation,branch_join = modes.analytic(full)
    specification = json.loads(args.input.read_text())
    parameters = {k:sp.Rational(v) for k,v in specification['parameters'].items()}
    bindings = {a:parameters[a.name] for a in relation.free_symbols-{modes.k,modes.q}}
    q_squared = sp.solve(relation,modes.q**2)[0].xreplace(bindings)
    evaluate = sp.lambdify(modes.k,q_squared,'cmath')
    namespace = dict(vars(sp),Str=Str)
    for line_number,line in enumerate(args.transcript.open(),1):
        if line.startswith('PY_S11CD_FULL_PENCIL_MODE_INPUT_REFERENCE_'):
            candidates = eval(line.partition(': ')[2],{'__builtins__':{}},namespace)
            candidate = {str(k):v for k,v in candidates[args.candidate_index]}
            break
    else:
        raise ValueError('reference input mode packet is absent')
    k_end,q_engine = complex(candidate['K']),complex(candidate['Q'])
    k_start = complex(k_end.real)
    initial_square = complex(evaluate(k_start))
    engine.emit('SHEET_PROBE_REAL_AXIS_OPERAND',engine.cas(initial_square))
    engine.emit('METADATA_SHEET_PROBE_REAL_AXIS_OPERAND',modes.numeric_metadata(
        engine.cas(initial_square),lambda p:tuple(2*v for v in dimensions.measure(modes.q))))
    if parameters['omega']<=0 or initial_square.imag!=0 or initial_square.real>=0:
        raise ValueError('probe requires the positive-frequency evanescent real-momentum chart')
    initial_root = cmath.sqrt(initial_square)
    complex_momentum = sp.Dummy('complexNormalMomentum')
    branch_points = [complex(root.evalf()) for root in sp.Poly(
        q_squared.xreplace({modes.k:complex_momentum}),complex_momentum).all_roots()]
    tangent = k_end-k_start
    closest = [k_start+min(1.,max(0.,((point-k_start)*tangent.conjugate()).real/abs(tangent)**2))*tangent
               for point in branch_points]
    results=[]
    for steps in (16,64,256):
        continued = initial_root
        samples=[]
        for index in range(steps+1):
            k = k_start+tangent*index/steps
            square = complex(evaluate(k))
            root = cmath.sqrt(square)
            continued = min((root,-root),key=lambda value:abs(value-continued))
            samples.append((k,continued,continued**2-square))
        results.append({'steps':steps,'startK':k_start,'endK':k_end,
            'initialQ':initial_root,'engineQ':q_engine,'continuedQ':continued,
            'enginePhysicalSheetFlag':candidate['PHYSICAL_BULK_SHEET'],
            'rootDifference':q_engine-continued,'oppositeRootResidual':q_engine+continued,
            'maximumRadicalResidual':max(abs(v[2]) for v in samples),
            'minimumRadicalMagnitude':min(abs(v[1]) for v in samples),
            'minimumBranchPointDistance':min(abs(a-b) for a,b in zip(branch_points,closest))})
    provenance = {'candidateIndex':args.candidate_index,'sourceLine':line_number,
        **{name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in
           (('symbolCacheSha256',args.symbol_cache),('transcriptSha256',args.transcript),
            ('inputSha256',args.input),('engineSha256',Path(engine.__file__)))}}
    engine.emit('SHEET_PROBE_PROVENANCE',provenance)
    engine.emit('METADATA_SHEET_PROBE_PROVENANCE',modes.numeric_metadata(
        engine.cas(provenance),lambda p:dimensions.zero))
    engine.physical('SHEET_PROBE_SOURCE_RADICAL_RELATION',relation)
    engine.physical('SHEET_PROBE_SOURCE_BRANCH_JOIN_RESIDUAL',branch_join,
                    zero_dimensions={(i,):dimensions.zero for i in range(len(branch_join))})
    engine.emit('SHEET_PROBE_CONTINUATIONS',results)
    def unit(path):
        name=path[1]
        if name in ('startK','endK','minimumBranchPointDistance'): return dimensions.measure(modes.k)
        if name=='maximumRadicalResidual': return tuple(2*v for v in dimensions.measure(modes.q))
        if name in ('initialQ','engineQ','continuedQ','rootDifference','oppositeRootResidual','minimumRadicalMagnitude'):
            return dimensions.measure(modes.q)
        return dimensions.zero
    engine.emit('METADATA_SHEET_PROBE_CONTINUATIONS',modes.numeric_metadata(engine.cas(results),unit))
    engine.emit('SHEET_PROBE_DIMENSION_CONSTRAINTS',sorted(dimensions.constraints,key=str))


if __name__=='__main__':
    run()
