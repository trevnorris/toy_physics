#!/usr/bin/env python3
"""Normalize the native four-case power residuals and fingerprint both controls."""
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = Path('/tmp/s11c-mechanical-repair-20260912')
RUN = BASE / 'c2_checks'
sys.path.insert(0, str(ROOT / 'scripts'))
import sympy as sp
from ledger_fold import _restore
import S11c_c2_selfenergy_fold_sympy_audit as c
import S11c_d_mixing_scattering_sympy_audit as d


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def object_orders(value, source):
    symbols = {symbol.name: symbol for symbol in source.free_symbols}
    # An absent source knob cannot occur in its selected or normalized object.
    knobs = tuple(symbols.get(name) for name in ('epsilon_shape', 'eta_bg', 'sigma_W'))
    support = set().union(*(c.grades(leaf, *knobs) for _, leaf in d.leaves(value)))
    return sp.Tuple(*(sp.Tuple(*grade) for grade in sorted(support)))


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
    started = time.monotonic()
    RUN.mkdir(exist_ok=True)
    producer = BASE / 'c2_full'
    manifest = json.loads((producer / 'manifest.json').read_text())
    transcript = producer / 'full.out'
    if manifest['exit_code'] != 0 or manifest['source_hashes_before'] != manifest['source_hashes_after']:
        raise ValueError('completed stable c2 producer required')
    if digest(transcript) != manifest['artifacts']['full.out']['sha256']:
        raise ValueError('c2 producer transcript pin')
    source_paths = (Path(__file__), ROOT / 'scripts/S11c_c2_selfenergy_fold_sympy_audit.py',
                    ROOT / 'scripts/S11c_d_mixing_scattering_sympy_audit.py', producer / 'manifest.json')
    pins = {str(path): digest(path) for path in source_paths}
    keys = set()
    records = []
    def output(key, value, dimensions, orders, heavy=False):
        if key in keys:
            raise ValueError(('duplicate write-key', key))
        keys.add(key)
        value = d.cas(value)
        representation = d.carrier_fingerprint(value) if heavy else value
        print(json.dumps({'writeKey': key, 'value': sp.srepr(representation),
                          'dimensionLTM': sp.srepr(dimensions), 'epsEtaSigmaOrder': sp.srepr(orders)},
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
        key = 's11cMechanicalRepairFullC2' + ''.join(word.title() for word in tag.removeprefix('PY_S11CC2_').split('_'))
        if tag.startswith('PY_S11CC2_TRACTION_CONSERVATIVE_VARIATION_'):
            for field in ('ACTION_TO_ROW_MULTIPLIER', 'KINETIC_NORMALIZATION_RESIDUAL'):
                selected = c.named(value, field)
                output(key+''.join(word.title() for word in field.split('_')), selected,
                       c.named(dimensions, field), object_orders(selected, value))
        else:
            output(key+'Operand', value, dimensions, orders, heavy=True)
            if tag.startswith('PY_S11CC2_TRACTION_SLAB_POWER_PAIRING_RESIDUAL_'):
                normalized = canonical_residual(value)
                output(key+'Canonical', normalized, dimensions, object_orders(normalized, value), heavy=True)
    report = {'sourcePins': pins, 'sourcePinsAfter': {str(path): digest(path) for path in source_paths},
              'producerManifestSha256': digest(producer / 'manifest.json'), 'records': records,
              'wallSeconds': time.monotonic()-started, 'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (RUN / 'inventory.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    run()
