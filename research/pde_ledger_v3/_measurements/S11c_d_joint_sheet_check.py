#!/usr/bin/env python3
"""Source-pinned reduced-branch preflight and joint-path construction checks."""
import argparse
import hashlib
import json
from pathlib import Path
import pickle
import resource
import sys
import time
from types import SimpleNamespace

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import S11c_d_mixing_scattering_sympy_audit as engine
from ledger_fold import _restore
from S11c_d_output_codec import decoded_lines


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(args):
    manifest = json.loads(args.manifest.read_text())
    base = Path(manifest['run_directory'])
    cache = base/'symbols'/f'{args.end}_LAB_HELD_RHO4_CONSTANT.pickle'
    if manifest.get('exit_code') != 0 or manifest['source_hashes_before'] != manifest['source_hashes_after']:
        raise ValueError('completed stable producer required')
    for path in (cache, base/'full.out'):
        if digest(path) != manifest['artifacts'][str(path.relative_to(base))]['sha256']:
            raise ValueError(('producer artifact mismatch', str(path)))
    for name, expected in manifest['source_hashes_after'].items():
        # Producer instruments can acquire a lossless transcript reader.
        # Verify their frozen producer sources; physical inputs remain pinned
        # against the current files consumed by the new construction.
        source = base/'source'/name if name.startswith('_measurements/') else ROOT/name
        if name not in ('scripts/S11c_d_mixing_scattering_sympy_audit.py','scripts/S11c_d_output_codec.py') and digest(source) != expected:
            raise ValueError(('producer input mismatch', name))
    frozen = base/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py'
    if digest(frozen) != manifest['source_hashes_after']['scripts/S11c_d_mixing_scattering_sympy_audit.py']:
        raise ValueError('producer source snapshot mismatch')
    full,curl,units,known,strong = pickle.loads(cache.read_bytes())[:5]
    symbols = {s.name:s for s in known if isinstance(s,sp.Symbol)}
    symbols.update({s.name:s for s in full.free_symbols | strong.free_symbols})
    dims = engine.DimensionAnalysis.__new__(engine.DimensionAnalysis)
    dims.known,dims.unknown,dims.constraints,dims.solution = dict(known),{},set(),{}
    dims.zero = (sp.S.Zero,)*3
    r = SimpleNamespace(symbols=symbols,omega=symbols['omega'],ell=symbols['L_W'],
        tangents=tuple(symbols['s11cdTangentialMomentum'+str(i)] for i in (1,2)),
        xi=symbols['s11cdProfileCoordinate'])
    r.end_values = {(p,end):sp.Limit(sp.Function('s11cd'+p.upper()+'Profile')(r.xi),r.xi,end)
                    for p in ('w','m') for end in (-sp.oo,sp.oo)}
    engine.PHYSICAL_METADATA = engine.PhysicalMetadata(dims,r)
    modes = engine.FullPencilModes(SimpleNamespace(r=r,kn=symbols['s11cdSpectralNormalMomentum']),curl,units)
    specification = json.loads(args.input.read_text())
    inputs = engine.ChannelInput(r,specification)
    prefix='PY_S11CD_REDUCED_BINDING_OPERANDS_s11cc2ClosedSlabOperator_LAB_HELD_RHO4_CONSTANT_COMPUTED_BRANCH_BINDINGS_'
    bindings=[]
    for line in decoded_lines(base/'full.out'):
        if line.startswith(prefix):bindings.append(_restore(line.partition(': ')[2]))
    if not bindings:raise ValueError('missing reduced branch operands')
    field_units=[dims.known[sp.Function('s11cdReducedField'+name)] for name in ('u1','u2','u3','theta','eW')]
    row_units=[]
    for i in range(5):
        j=next(j for j in range(5) if strong[i,j]!=0)
        row_units.append(tuple(a+b for a,b in zip(dims.measure(strong[i,j]),field_units[j])))
    strong_units={(5*i+j,):tuple(a-b for a,b in zip(row_units[i],field_units[j])) for i in range(5) for j in range(5)}
    provenance={'producerManifestSha256':digest(args.manifest),'cacheSha256':digest(cache),
        'producerTranscriptSha256':digest(base/'full.out'),'producerSources':manifest['source_hashes_after'],
        'engineSha256':digest(Path(engine.__file__)),'instrumentSha256':digest(Path(__file__)),
        'inputSha256':digest(args.input),'end':args.end,'case':'LAB_HELD_RHO4_CONSTANT'}
    return modes,strong,strong_units,bindings,inputs,provenance


def run():
    started=time.monotonic()
    parser=argparse.ArgumentParser()
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--end',choices=('REFERENCE','LEFT','RIGHT'),default='REFERENCE')
    parser.add_argument('--paths',action='store_true')
    parser.add_argument('--pit',action='store_true')
    args=parser.parse_args()
    modes,strong,units,bindings,inputs,provenance=load(args)
    dims=engine.PHYSICAL_METADATA.dimensions
    def output(tag,value,unit=lambda p:dims.zero,heavy=False):
        engine.emit('JOINT_SHEET_PREFLIGHT_'+tag,modes.compact_fingerprint(value) if heavy else value)
        engine.emit('METADATA_JOINT_SHEET_PREFLIGHT_'+tag,modes.numeric_metadata(engine.cas(value),unit))
    output('PROVENANCE',provenance)
    if args.paths:
        manifest=json.loads(args.manifest.read_text())
        prefix='PY_S11CD_END_SPECTRUM_'+('PIT_' if args.pit else 'INPUT_')+args.end+'_LAB_HELD_RHO4_CONSTANT_0_MODE_'
        records=[]
        for line in decoded_lines(Path(manifest['run_directory'])/'full.out'):
            tag,_,payload=line.partition(': ')
            if tag.startswith(prefix) and tag.endswith('_RECORD'):
                records.append({str(k):v for k,v in _restore(payload)})
        if not records:raise ValueError('missing producer native spectrum')
        audit=engine.BulkContinuationAudit(modes,strong,units,bindings)
        audit.construct(args.end+'_LAB_HELD_RHO4_CONSTANT',channel_input=None if args.pit else inputs,
                        reference=args.end=='REFERENCE',spectrum={'RECORDS':records})
        output('DIMENSION_CONSTRAINTS',tuple(dims.constraints))
        output('RESOURCES',{'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
               lambda p:(0,1,0) if p[0]=='wallSeconds' else dims.zero)
        return
    algebraic,relation,joins=modes.analytic(strong)
    mapping=inputs.mapping(algebraic,relation,(modes.k,modes.q,modes.eta,modes.sigma,modes.r.omega))
    origin={modes.eta:0,modes.sigma:0} if args.end=='REFERENCE' else inputs.origin
    original=strong.xreplace(mapping).subs(origin)
    algebraic=algebraic.xreplace(mapping).subs(origin)
    relation=relation.xreplace(mapping)
    square=sp.solve(relation,modes.q**2)[0]
    scale=sp.sqrt(-sp.Poly(square,modes.k).nth(2))
    seeds=[scale*rhs.xreplace(dict(zip(lhs.args,(*modes.r.tangents,modes.k)))).xreplace(mapping)
           for lhs,rhs in bindings]
    output('SOURCE_BRANCH_JOINS',joins)
    output('RELATION',relation,lambda p:tuple(2*v for v in dims.measure(modes.q)))
    seed_differences=tuple(sp.simplify(seed-seeds[0]) for seed in seeds)
    output('REDUCED_SEED_DIFFERENCES',seed_differences,lambda p:dims.measure(modes.q))
    k=inputs.parameters[modes.r.tangents[0].name]
    cone=max(sp.solve(square.subs(modes.k,k),modes.r.omega))
    failures=[]
    for i,factor in enumerate((sp.Rational(1,2),sp.Integer(2),sp.Rational(-1,2),sp.Integer(-2))):
        w=factor*cone
        point={modes.r.omega:w,modes.k:k}
        q=sp.simplify(seeds[0].subs(point))
        source=original.subs(point).applyfunc(sp.simplify)
        continued=algebraic.subs({**point,modes.q:q}).applyfunc(sp.simplify)
        residual=(source-continued).applyfunc(sp.simplify)
        output('REAL_AXIS_'+str(i),{'OMEGA':w,'K':k,'Q':q},lambda p:
               dims.measure(modes.k) if p[0]=='K' else dims.measure(modes.q))
        output('SOURCE_MATRIX_'+str(i),source,lambda p:units[p],heavy=True)
        output('ALGEBRAIC_MATRIX_'+str(i),continued,lambda p:units[p],heavy=True)
        output('MATRIX_JOIN_RESIDUAL_'+str(i),residual,lambda p:units[p])
        output('RADICAL_RESIDUAL_'+str(i),sp.simplify(relation.subs({**point,modes.q:q})),lambda p:tuple(2*v for v in dims.measure(modes.q)))
        failures.extend(v for v in residual if v!=0)
    output('DIMENSION_CONSTRAINTS',tuple(dims.constraints))
    output('RESOURCES',{'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
           lambda p:(0,1,0) if p[0]=='wallSeconds' else dims.zero)
    if failures:raise ValueError('real-axis source/algebraic join differs')


if __name__=='__main__':run()
