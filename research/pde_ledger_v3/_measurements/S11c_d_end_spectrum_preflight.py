#!/usr/bin/env python3
"""Bound-input comparison of cached physical and potential-coordinate pencils."""
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


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    started = time.monotonic()
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--end', choices=('REFERENCE', 'LEFT', 'RIGHT'), default='RIGHT')
    parser.add_argument('--cache-result', type=Path)
    parser.add_argument('--coverage', action='store_true')
    parser.add_argument('--zero-tangent', action='store_true')
    parser.add_argument('--pit-sample', type=int)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    base = Path(manifest['run_directory'])
    cache = base/'symbols'/f'{args.end}_LAB_HELD_RHO4_CONSTANT.pickle'
    for path in (cache, base/'full.out'):
        if digest(path) != manifest['artifacts'][str(path.relative_to(base))]['sha256']:
            raise ValueError(('producer artifact changed', str(path)))
    for name, expected in manifest['source_hashes_after'].items():
        if name != 'scripts/S11c_d_mixing_scattering_sympy_audit.py' and digest(ROOT/name) != expected:
            raise ValueError(('producer input changed', name))
    full, curl, units, known, strong = pickle.loads(cache.read_bytes())[:5]
    symbols = {s.name:s for s in known if isinstance(s, sp.Symbol)}
    symbols.update({s.name:s for s in full.free_symbols | strong.free_symbols})
    dims = engine.DimensionAnalysis.__new__(engine.DimensionAnalysis)
    dims.known, dims.unknown, dims.constraints, dims.solution = dict(known), {}, set(), {}
    dims.zero = (sp.S.Zero,)*3
    r = SimpleNamespace(symbols=symbols, omega=symbols['omega'], ell=symbols['L_W'],
                        z=sp.Symbol('s11cdPreflightNormalCoordinate', real=True),
                        tangents=tuple(symbols['s11cdTangentialMomentum'+str(i)] for i in (1,2)))
    engine.PHYSICAL_METADATA = engine.PhysicalMetadata(dims, r)
    ends = SimpleNamespace(r=r, kn=symbols['s11cdSpectralNormalMomentum'])
    modes = engine.FullPencilModes(ends, curl, units)
    specification = json.loads(args.input.read_text())
    if args.zero_tangent:
        specification['parameters']['s11cdTangentialMomentum1'] = '0'
    parameters = {key:sp.Rational(value) for key,value in specification['parameters'].items()}
    xi = sp.Symbol('xi', real=True)
    profiles = {key:sp.sympify(value, locals={'xi':xi}) for key,value in specification['profiles'].items()}
    algebraic, relation, join = modes.analytic(full)
    physical, physical_relation, physical_join = modes.analytic(strong)
    live = {modes.k, modes.q, modes.eta, modes.sigma}
    mapping = {s:parameters[s.name] for s in (algebraic.free_symbols | physical.free_symbols | relation.free_symbols)-live}
    for limit in algebraic.atoms(sp.Limit) | physical.atoms(sp.Limit):
        # Bind only the actual profile-function names present in this cache.
        if limit.args[0].func.__name__ == 's11cdWProfile': profile = 'w'
        elif limit.args[0].func.__name__ == 's11cdMProfile': profile = 'm'
        else: raise ValueError(('unbound profile limit', limit))
        mapping[limit] = sp.limit(profiles[profile], xi, limit.args[2])
    origin = {modes.eta:sp.S.Zero, modes.sigma:sp.S.Zero} if args.end == 'REFERENCE' else {
        modes.eta:parameters['eta_bg'], modes.sigma:parameters['eta_bg']*parameters['W_0']/parameters['L_W']}
    quotient = algebraic.extract(modes.indices, modes.indices).xreplace(mapping).subs(origin).applyfunc(sp.cancel)
    physical = physical.xreplace(mapping).subs(origin).applyfunc(sp.cancel)
    bound_relation = relation.xreplace(mapping)

    def output(tag, value, unit=lambda path:dims.zero):
        engine.emit('END_SPECTRUM_PREFLIGHT_'+tag, value)
        engine.emit('METADATA_END_SPECTRUM_PREFLIGHT_'+tag, modes.numeric_metadata(engine.cas(value), unit))

    output('PROVENANCE', {'producerManifestSha256':digest(args.manifest), 'cacheSha256':digest(cache),
        'producerTranscriptSha256':digest(base/'full.out'), 'producerSources':manifest['source_hashes_after'],
        'engineSha256':digest(Path(engine.__file__)), 'instrumentSha256':digest(Path(__file__)),
        'inputSha256':digest(args.input), 'end':args.end, 'case':'LAB_HELD_RHO4_CONSTANT',
        'zeroTangentMutation':args.zero_tangent,'pitSample':args.pit_sample if args.pit_sample is not None else 'PHYSICAL_INPUT'})
    output('INPUT', specification)
    output('BOUND_CARRIERS', tuple((str(key),value) for key,value in mapping.items()),
           lambda path:dims.measure(next(key for key in mapping if str(key)==path[0])))
    output('GRADE_ORIGIN', tuple(origin.items()))
    output('BRANCH_JOIN_RESIDUALS', join+physical_join)
    output('RADICAL_RELATION_RESIDUAL', sp.expand(physical_relation-relation),
           lambda path:tuple(2*v for v in dims.measure(modes.q)))
    output('RADICAL_RELATION', bound_relation, lambda path:tuple(2*v for v in dims.measure(modes.q)))

    # Reuse the engine's field ansatz and differentiate it to construct the lift.
    # The tangential phase below is the harmonic ansatz, not a spectral answer.
    pencil = engine.ReducedPencil.__new__(engine.ReducedPencil)
    pencil.r = r
    pencil.fields = tuple(sp.Function('s11cdReducedField'+name) for name in ('u1','u2','u3','theta','eW'))
    tangent_coordinates = sp.symbols('s11cdPreflightX0:2', real=True)
    plane = sp.exp(sp.I*sum(k*x for k,x in zip(r.tangents,tangent_coordinates)))
    pencil.tangent_derivatives = tuple(sp.diff(plane,x)/plane for x in tangent_coordinates)
    lift_builder = engine.ConstantEndPencil.__new__(engine.ConstantEndPencil)
    lift_builder.r, lift_builder.kn = r, modes.k
    lift = lift_builder.field_lift(pencil).extract(range(5), modes.indices).xreplace(mapping)
    # The dual test is the harmonic character at reversed momentum, followed
    # by transpose, as in the engine's weak virtual-work projection.
    # Tangential substitution must precede binding, just as for the normal leg.
    unbound_lift = lift_builder.field_lift(pencil).extract(range(5), modes.indices)
    dual = unbound_lift.xreplace({s:-s for s in (*r.tangents,modes.k)}).xreplace(mapping).T
    if args.coverage:
        r.xi = symbols['s11cdProfileCoordinate']
        r.end_values = {(p,end):sp.Limit(sp.Function('s11cd'+p.upper()+'Profile')(r.xi),r.xi,end)
                        for p in ('w','m') for end in (-sp.oo,sp.oo)}
        channel_input = engine.ChannelInput(r,specification)
        field_units = [dims.known[f] for f in pencil.fields]
        row_units = []
        for i in range(5):
            j = next(j for j in range(5) if strong[i,j]!=0)
            row_units.append(tuple(a+b for a,b in zip(dims.measure(strong[i,j]),field_units[j])))
        strong_units = {(5*i+j,):tuple(a-b for a,b in zip(row_units[i],field_units[j]))
                        for i in range(5) for j in range(5)}
        coverage = engine.EndSpectrumCoverage(modes,strong,full,lift_builder.field_lift(pencil),strong_units)
        coverage.construct(args.end+'_LAB_HELD_RHO4_CONSTANT',-1 if args.end=='LEFT' else 1,
                           channel_input=channel_input if args.pit_sample is None else None,
                           reference=args.end=='REFERENCE',sample_index=args.pit_sample or 0)
        output('DIMENSION_CONSTRAINTS',tuple(dims.constraints))
        output('RESOURCES', {'wallSeconds':time.monotonic()-started,
                            'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
               lambda path:(0,1,0) if path[0]=='wallSeconds' else dims.zero)
        return
    route = (dual*physical*lift).applyfunc(sp.cancel)
    residual = (quotient-route).applyfunc(sp.cancel)
    output('LIFT_DETERMINANTS', (sp.factor(lift.det()), sp.factor(dual.det())),
           lambda path:dims.measure(unbound_lift.det()))
    output('PHYSICAL_PENCIL', engine.carrier_fingerprint(physical))
    output('QUOTIENT_PENCIL', engine.carrier_fingerprint(quotient))
    output('PULLBACK_PENCIL', engine.carrier_fingerprint(route))
    output('PULLBACK_RESIDUAL', residual,
           lambda path:units[(6*modes.indices[path[0]//5]+modes.indices[path[0]%5],)])
    k_square = sp.solve(bound_relation, modes.k**2)[0]
    divisor = sp.Poly(modes.k**2-k_square, modes.k)
    results = {}
    for label, matrix in (('PHYSICAL',physical), ('QUOTIENT',quotient)):
        (numerator,denominator), cleared, row_denominators = modes.rational_determinant(matrix)
        reduced = sp.rem(sp.Poly(numerator,modes.k),divisor).as_expr()
        output(label+'_ELIMINATION_REMAINDER', reduced.has(modes.k))
        polynomial = sp.Poly(reduced,modes.q)
        factorization = sp.factor_list(polynomial)
        output(label+'_DEGREE', polynomial.degree())
        output(label+'_FACTORS', tuple((factor.degree(), multiplicity) for factor,multiplicity in factorization[1]))
        output(label+'_DETERMINANT_DENOMINATOR', engine.carrier_fingerprint(denominator))
        output(label+'_ELIMINATION', engine.carrier_fingerprint(polynomial.as_expr()))
        roots = tuple((complex(root),multiplicity) for factor,multiplicity in factorization[1]
                      for root in sp.nroots(factor,n=40,maxsteps=300))
        output(label+'_RADICAL_ROOTS', roots, lambda path:dims.measure(modes.q) if path[1]==0 else dims.zero)
        results[label] = {'matrix':matrix,'numerator':numerator,'denominator':denominator,
                          'polynomial':polynomial,'roots':roots,'row_denominators':row_denominators}
    output('DETERMINANT_PULLBACK_RESIDUAL', sp.cancel(
        results['QUOTIENT']['numerator']/results['QUOTIENT']['denominator']
        -dual.det()*lift.det()*results['PHYSICAL']['numerator']/results['PHYSICAL']['denominator']))
    if args.cache_result:
        args.cache_result.write_bytes(pickle.dumps({'results':results,'relation':bound_relation,
            'k':modes.k,'q':modes.q,'mapping':mapping,'origin':origin,'lift':lift,'dual':dual}))
    output('DIMENSION_CONSTRAINTS', tuple(dims.constraints))
    output('RESOURCES', {'wallSeconds':time.monotonic()-started,
                        'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
           lambda path:(0,1,0) if path[0]=='wallSeconds' else dims.zero)


if __name__ == '__main__':
    run()
