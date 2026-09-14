#!/usr/bin/env python3
"""Normalize the native four-case power residuals and fingerprint both controls."""
import argparse
import hashlib
import json
from pathlib import Path
import pickle
import resource
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
RUN = None
sys.path.insert(0, str(ROOT / 'scripts'))
import sympy as sp
from ledger_fold import _restore
import S11c_c2_selfenergy_fold_sympy_audit as c
import S11c_d_mixing_scattering_sympy_audit as d
from S11c_c2_trace_repair_check import TraceMetadata


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def object_orders(value, generators, homotopy_ratio, profiles):
    metadata = TraceMetadata.__new__(TraceMetadata)
    metadata.generators = generators
    grades=set();support=set()
    for _,leaf in d.leaves(value):
        if isinstance(leaf,c.Str):continue
        coefficients=metadata.coefficients(leaf.xreplace(profiles));grades.update(coefficients)
        homotopy={}
        for (e,a,b),coefficient in coefficients.items():
            term=coefficient*homotopy_ratio**b
            homotopy[e,a+b]=homotopy.get((e,a+b),sp.S.Zero)+term
        support.update(g for g,v in homotopy.items() if v!=0)
    return (sp.Tuple(*(sp.Tuple(*g) for g in sorted(grades))),
            sp.Tuple(*(sp.Tuple(*g) for g in sorted(support))))


def nonfinite_syntax(expression):
    """Count unbounded domains separately; this is not a convergence test."""
    domains = invalid = 0
    def visit(value):
        nonlocal domains, invalid
        if value in (sp.nan, sp.zoo, sp.oo, -sp.oo):
            invalid += 1
        elif isinstance(value, sp.Integral):
            visit(value.function)
            for limits in value.limits:
                visit(limits[0])
                for bound in limits[1:]:
                    if bound in (sp.oo, -sp.oo):
                        domains += 1
                    else:
                        visit(bound)
        elif isinstance(value, sp.Limit):
            visit(value.args[0])
            visit(value.args[1])
            if value.args[2] in (sp.oo, -sp.oo):
                domains += 1
            else:
                visit(value.args[2])
        elif isinstance(value, sp.Basic):
            for argument in value.args:
                visit(argument)
    visit(expression)
    return {'infiniteDomainEndpoints': domains, 'nonfiniteOutsideDomainEndpoints': invalid}


def canonical_residual(expression):
    carriers = expression.atoms(sp.Integral, sp.Derivative)
    mapping = {carrier: sp.Dummy('mechanicalRepairCarrier') for carrier in sorted(carriers, key=sp.default_sort_key)}
    shielded = expression.xreplace(mapping)
    active = tuple(atom for atom in mapping.values() if atom in shielded.free_symbols)
    if not active:
        return sp.cancel(shielded)
    reduced = sp.S.Zero
    terms = sp.Poly(shielded, *active).terms()
    for index, (powers, coefficient) in enumerate(terms):
        (RUN / 'progress.json').write_text(json.dumps({'carrierCoefficient': index+1, 'count': len(terms)})+'\n')
        reduced += sp.cancel(coefficient) * sp.Mul(*(atom**power for atom, power in zip(active, powers)))
    return reduced.xreplace({value: key for key, value in mapping.items()})


def run():
    global RUN
    parser=argparse.ArgumentParser()
    parser.add_argument('--run-directory',type=Path,required=True)
    parser.add_argument('--producer',type=Path,required=True)
    parser.add_argument('--inputs-cache',type=Path,required=True)
    args=parser.parse_args();RUN=args.run_directory
    started = time.monotonic()
    RUN.mkdir(parents=True,exist_ok=False)
    producer = args.producer
    manifest = json.loads((producer / 'manifest.json').read_text())
    transcript = producer / 'full.out'
    if manifest['exit_code'] != 0 or manifest['source_hashes_before'] != manifest['source_hashes_after']:
        raise ValueError('completed stable c2 producer required')
    if digest(transcript) != manifest['artifacts']['full.out']['sha256']:
        raise ValueError('c2 producer transcript pin')
    cache_meta=args.inputs_cache.with_suffix('.json')
    cache_pins=json.loads(cache_meta.read_text())
    if digest(args.inputs_cache)!=cache_pins['sha256']:raise ValueError('power input cache digest')
    for name,pin in cache_pins['sourcePins'].items():
        if manifest['source_hashes_after'][name]!=pin or digest(ROOT/name)!=pin:
            raise ValueError(('power input cache source',name))
    input_values,_=pickle.loads(args.inputs_cache.read_bytes())
    inputs=c.Inputs(input_values)
    generators=tuple(input_values[name] for name in ('epsilon_shape','eta_bg','sigma_W'))
    homotopy_ratio=input_values['W_0']/input_values['L_W']
    source_paths = (Path(__file__), ROOT / 'scripts/S11c_c2_selfenergy_fold_sympy_audit.py',
                    ROOT / 'scripts/S11c_d_mixing_scattering_sympy_audit.py',
                    ROOT / '_measurements/S11c_c2_trace_repair_check.py', producer / 'manifest.json',
                    args.inputs_cache,cache_meta)
    pins = {str(path): digest(path) for path in source_paths}
    keys = set()
    records = []
    def output(key, value, dimensions, heavy=False):
        if key in keys:
            raise ValueError(('duplicate write-key', key))
        keys.add(key)
        value = d.cas(value)
        if dimensions.has(sp.nan,sp.zoo):raise ValueError(('native power dimensions',key))
        orders,homotopy=object_orders(value,generators,homotopy_ratio,inputs.profiles)
        representation = d.carrier_fingerprint(value) if heavy else value
        print(json.dumps({'writeKey': key, 'value': sp.srepr(representation),
                          'dimensionLTM': sp.srepr(dimensions), 'epsEtaSigmaOrder': sp.srepr(orders),
                          'epsilonLambdaSupport':sp.srepr(homotopy)},
                         separators=(',', ':')), flush=True)
        records.append({'writeKey': key, 'scalarLeaves': sum(1 for _ in d.leaves(value)),
                        'nonzeroLeaves': sum(leaf != 0 for _, leaf in d.leaves(value)),
                        **nonfinite_syntax(value),
                        'nonfiniteRepresentationLeaves': sum(leaf.has(sp.nan, sp.zoo, sp.oo, -sp.oo)
                                                             for _, leaf in d.leaves(representation))})
        if heavy:
            records[-1]['pitSamples'] = [[sp.sstr(sample) for sample in row[2]] for row in representation[0][1]]
    wanted = ('TRACTION_SLAB_POWER_PAIRING_RESIDUAL_', 'TRACTION_SIGN_RESIDUAL_',
              'TRACTION_ROW_ROUTING_RESIDUAL_', 'TRACTION_CONSERVATIVE_VARIATION_')
    for line in transcript.open():
        tag, separator, payload = line.partition(': ')
        if not separator or not any(tag.startswith('PY_S11CC2_'+name) for name in wanted):
            continue
        obj = _restore(payload.rstrip())
        value = c.named(obj, 'VALUE')
        dimensions, orders = c.named(obj, 'DIMENSION_L_T_M'), c.named(obj, 'MULTIGRADE')
        key = 's11cc2TraceRepairFullPower' + ''.join(word.title() for word in tag.removeprefix('PY_S11CC2_').split('_'))
        if tag.startswith('PY_S11CC2_TRACTION_CONSERVATIVE_VARIATION_'):
            for field in ('ACTION_TO_ROW_MULTIPLIER', 'KINETIC_NORMALIZATION_RESIDUAL'):
                selected = c.named(value, field)
                output(key+''.join(word.title() for word in field.split('_')), selected,
                       c.named(dimensions, field))
        else:
            output(key+'Operand', value, dimensions, heavy=True)
            if tag.startswith('PY_S11CC2_TRACTION_SLAB_POWER_PAIRING_RESIDUAL_'):
                normalized = canonical_residual(value)
                output(key+'Canonical', normalized, dimensions, heavy=True)
    report = {'sourcePins': pins, 'sourcePinsAfter': {str(path): digest(path) for path in source_paths},
              'producerManifestSha256': digest(producer / 'manifest.json'), 'records': records,
              'wallSeconds': time.monotonic()-started, 'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (RUN / 'inventory.json').write_text(json.dumps(report, indent=2)+'\n')
    canonical=[row for row in records if row['writeKey'].endswith('Canonical')]
    if len(canonical)!=4 or any(row['nonzeroLeaves'] for row in canonical):
        raise ValueError('four-case canonical power residual census; inspect emitted operands')
    if report['sourcePins']!=report['sourcePinsAfter']:
        raise ValueError('power producer sources changed during execution')


if __name__ == '__main__':
    run()
