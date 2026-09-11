#!/usr/bin/env python3
"""One-case source/image mutations of the computed Fourier operands.

Mathematical transform controls; these do not implement the physical profile
FORM or channel/flux controls in Section 5.
"""
import argparse
import hashlib
import json
from pathlib import Path
import pickle
import resource
import sys
import time

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import S11c_d_mixing_scattering_sympy_audit as engine
from ledger_fold import _restore


def run():
    started = time.monotonic()
    parser = argparse.ArgumentParser()
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--producer', type=Path, required=True)
    parser.add_argument('--production', action='store_true',
                        help='check every carrier and both action rows of the cached case')
    parser.add_argument('--carriers-only', action='store_true',
                        help='stop the production check before action-row contraction')
    args = parser.parse_args()
    prefix = 'PY_S11CD_BUILD_INPUT_DIGESTS: '
    pins = dict((str(k), str(v)) for k, v in next(
        _restore(line[len(prefix):]) for line in args.producer.open() if line.startswith(prefix)))
    for path, digest in pins.items():
        if path != 'scripts/S11c_d_mixing_scattering_sympy_audit.py':
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != digest:
                raise ValueError(('changed producer input', path))
    rows, known = pickle.loads(args.cache.read_bytes())
    r = engine.EdgeReduction(rows)
    dims = engine.DimensionAnalysis.__new__(engine.DimensionAnalysis)
    dims.known, dims.unknown, dims.constraints, dims.solution = known, {}, set(), {}
    dims.zero = (sp.S.Zero,)*3
    engine.PHYSICAL_METADATA = engine.PhysicalMetadata(dims, r)
    inverse = engine.FourierCarrierReconstruction(r)
    provenance = {'producer': str(args.producer), 'producer_pins': pins,
                  'cache_sha256': hashlib.sha256(args.cache.read_bytes()).hexdigest(),
                  'producer_sha256': hashlib.sha256(args.producer.read_bytes()).hexdigest(),
                  'engine_sha256': hashlib.sha256(Path(engine.__file__).read_bytes()).hexdigest(),
                  'instrument_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    engine.physical('INVERSE_CHECK_PROVENANCE', provenance)
    inverse.emit_kernel()
    if args.production or args.carriers_only:
        reconstruction = engine.EdgeReconstruction(r, inverse)
        for row_key in engine.CLOSED_KEYS:
            case, payload = next((case, payload) for case, payload in rows[row_key]['value']
                       if tuple(map(str, case)) == ('LAB_HELD', 'RHO4_CONSTANT'))
            suffix = '_'.join((row_key, *map(str, case)))
            definitions = tuple(engine.named(payload, 'FOURIER_PROFILE_BINDINGS'))
            for index, node in enumerate(sorted(payload.atoms(engine.AppliedUndef), key=sp.default_sort_key)):
                if node.func.__name__.startswith('s11cc2Fourier'):
                    inverse.carrier(node, definitions, suffix+'_'+str(index))
            if args.carriers_only:
                continue
            value = engine.named(payload, 'VALUE')
            units = engine.payload_units(value, engine.named(payload, 'DIMENSION_L_T_M'))
            reduced, records = r.value(value)
            reconstruction.row(value, reduced, records, suffix, units, definitions)
        engine.emit('INVERSE_CHECK_DIMENSION_CONSTRAINTS', tuple(dims.constraints))
        engine.emit('INVERSE_CHECK_DIMENSION_UNRESOLVED', tuple(dims.unknown))
        engine.emit('RESOURCE_MEASUREMENTS', (time.monotonic()-started,
                                             resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        engine.emit('PROCESS_COMPLETION', 'INVERSE_CHECK')
        return
    payload = next(payload for case, payload in rows[engine.CLOSED_KEYS[0]]['value']
                   if tuple(map(str, case)) == ('LAB_HELD', 'RHO4_CONSTANT'))
    definitions = tuple(engine.named(payload, 'FOURIER_PROFILE_BINDINGS'))
    coordinates = (inverse.p, -inverse.p, sp.S.Zero)
    for index, equation in enumerate(definitions):
        node = equation.lhs
        image = r.hat(node)[0]
        argument = node.args[2].xreplace(r.normal_map)
        reconstructed, _ = inverse.invert_image(image, argument)
        for label, multiplier, reverse in (('SOURCE_WEIGHT', 2, False), ('SOURCE_PHASE', 1, True)):
            source = equation.rhs
            phases = list(source.function.atoms(sp.exp))
            if len(phases) != 1:
                raise NotImplementedError(('mutation source phases', len(phases)))
            body = source.function*multiplier
            if reverse:
                body = body.xreplace({phases[0]: sp.exp(-phases[0].args[0])})
            changed = sp.Eq(node, sp.Integral(body, *source.limits), evaluate=False)
            mutated = tuple(changed if eq == equation else eq for eq in definitions)
            source_point, operands = inverse.source_inverse(node.func, mutated)
            source_values = tuple(source_point.subs(inverse.x[2], p) for p in coordinates)
            residual = tuple(sp.simplify(a-b) for a, b in zip(reconstructed, source_values))
            suffix = label+'_'+str(index)
            inverse.physical_zeros('INVERSE_CHECK_SOURCE_OPERANDS_'+suffix, (source, changed.rhs))
            inverse.physical_zeros('INVERSE_CHECK_VALUES_'+suffix, (reconstructed, source_values))
            inverse.physical_zeros('INVERSE_CHECK_RESIDUAL_'+suffix, residual)
        # Negate normal momenta in the computed image only. The source
        # definition and the inverse argument chart retain their original sign.
        reversal = {r.normal_map[g[2]]: -r.normal_map[g[2]] for g in r.momentum_groups}
        changed_image = image.xreplace(reversal)
        changed_inverse, _ = inverse.invert_image(changed_image, argument)
        residual = tuple(sp.simplify(a-b) for a, b in zip(changed_inverse, reconstructed))
        suffix = 'IMAGE_PHASE_'+str(index)
        inverse.physical_zeros('INVERSE_CHECK_IMAGE_OPERANDS_'+suffix, (image, changed_image),
                               {(0,): dims.measure(r.ell), (1,): dims.measure(r.ell)})
        inverse.physical_zeros('INVERSE_CHECK_VALUES_'+suffix, (reconstructed, changed_inverse))
        inverse.physical_zeros('INVERSE_CHECK_RESIDUAL_'+suffix, residual)
        if image.has(r.regulator):
            integrals = sorted(image.atoms(sp.Integral), key=sp.default_sort_key)
            localized_image = sp.Add(*(coefficient*sp.prod(v**p for v, p in zip(integrals, powers))
                for powers, coefficient in engine.polynomial_terms(image, integrals) if any(powers)))
            localized_inverse, _ = inverse.invert_image(localized_image, argument)
            residual = tuple(sp.simplify(a-b) for a, b in zip(localized_inverse, reconstructed))
            suffix = 'IMAGE_ABEL_REMOVAL_'+str(index)
            inverse.physical_zeros('INVERSE_CHECK_IMAGE_OPERANDS_'+suffix,
                                   (image, localized_image, image-localized_image))
            inverse.physical_zeros('INVERSE_CHECK_VALUES_'+suffix, (reconstructed, localized_inverse))
            inverse.physical_zeros('INVERSE_CHECK_RESIDUAL_'+suffix, residual)
    engine.emit('INVERSE_CHECK_DIMENSION_CONSTRAINTS', tuple(dims.constraints))
    engine.emit('INVERSE_CHECK_DIMENSION_UNRESOLVED', tuple(dims.unknown))
    engine.emit('RESOURCE_MEASUREMENTS', (time.monotonic()-started,
                                         resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    engine.emit('PROCESS_COMPLETION', 'INVERSE_CHECK')


if __name__ == '__main__':
    run()
