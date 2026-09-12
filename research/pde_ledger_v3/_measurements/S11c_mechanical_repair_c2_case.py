#!/usr/bin/env python3
"""One native repaired b-case operand, baseline c1 response, and c2 power checks."""
import hashlib
import json
from pathlib import Path
import pickle
import resource
import sys
import time

STARTED=time.monotonic()
ROOT=Path(__file__).resolve().parents[1]
BASE=Path('/tmp/s11c-mechanical-repair-20260912')
RUN=BASE/'c2_case'
sys.path.insert(0,str(ROOT/'scripts'))
import sympy as sp
from sympy.core.symbol import Str
from ledger_fold import load_model,check_consumer,assert_lookups_equal_manifest
import S11c_c2_selfenergy_fold_sympy_audit as c
import S11c_d_mixing_scattering_sympy_audit as fingerprints


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def step(name):
    (RUN/'progress.json').write_text(json.dumps({'stage':name,'wallSeconds':time.monotonic()-STARTED})+'\n')


def replace_case(root,case,value):
    # This is a declared native input-operand replacement in a preflight;
    # it is not a generated export or a claim that all four cases were rebuilt.
    return sp.Tuple(*(sp.Tuple(axes,sp.Tuple(*(sp.Tuple(key,c.cas(value) if str(key)=='VALUE' else member)
                                             for key,member in payload))) if tuple(map(str,axes))==case
                      else sp.Tuple(axes,payload) for axes,payload in root))


def output(name,value,inputs,*,units=None,heavy=False):
    value=c.cas(value)
    print(json.dumps({'writeKey':'s11cMechanicalRepairC2'+name,
        'value':sp.srepr(fingerprints.carrier_fingerprint(value) if heavy else value),
        'multigrade':sp.srepr(c.grade_object(value,inputs)),
        'dimensionLTM':sp.srepr(c.cas(c.dimensions(value) if units is None else units))},separators=(',',':')),flush=True)


def canonical_residual(expression):
    """Scalar carrier algebra; retain uncombined calculus carriers explicitly."""
    carriers=expression.atoms(sp.Integral,sp.Derivative)
    mapping={carrier:sp.Dummy('mechanicalRepairCarrier') for carrier in sorted(carriers,key=sp.default_sort_key)}
    shielded=expression.xreplace(mapping)
    active=tuple(atom for atom in mapping.values() if atom in shielded.free_symbols)
    if not active:
        return sp.cancel(shielded)
    terms=sp.Poly(shielded,*active).terms()
    reduced=sp.S.Zero
    for index,(powers,coefficient) in enumerate(terms):
        if index%10==0:
            step('canonical carrier coefficient '+str(index+1)+'/'+str(len(terms)))
        reduced+=sp.cancel(coefficient)*sp.Mul(*(atom**power for atom,power in zip(active,powers)))
    return reduced.xreplace({v:k for k,v in mapping.items()})


def run():
    RUN.mkdir(exist_ok=True)
    case=('LAB_HELD','RHO4_CONSTANT')
    paths=(ROOT/'scripts/S11c_c2_selfenergy_fold_sympy_audit.py',Path(__file__),
           ROOT/'scripts/S11c_b_brane_operator_sympy_audit.py',BASE/'b_case/objects.pickle',
           BASE/'b_case/manifest.json',BASE/'baseline/scripts/S11c_b_exports.py',
           BASE/'baseline/scripts/S11c_c1_exports.py',ROOT/'directives/S11c_c2_SHARED_PHYSICS.md')
    pins={str(path):digest(path) for path in paths}
    native=json.loads((BASE/'b_case/manifest.json').read_text())
    if native['sourcePins']!=native['sourcePinsAfter'] or native['cacheSha256']!=digest(BASE/'b_case/objects.pickle'):
        raise ValueError('native b case provenance')
    step('load baseline fold and native case')
    fold,fold_audit=load_model(str(BASE/'baseline/scripts/S11c_b_exports.py'),
                               str(BASE/'baseline/scripts/S11c_c1_exports.py'))
    closure=check_consumer(fold,c.IMPORT_KEYS)
    lookups=assert_lookups_equal_manifest(c.bind_inputs,fold,c.IMPORT_KEYS)
    original=lookups['result']
    values=dict(original.values)
    operator,origins,mu=pickle.loads((BASE/'b_case/objects.pickle').read_bytes())
    values['slab_operator']=replace_case(values['slab_operator'],case,operator)
    values['slab_operator_term_origins']=replace_case(values['slab_operator_term_origins'],case,origins)
    values['mu_theta_operator']=replace_case(values['mu_theta_operator'],case,mu)
    inputs=c.Inputs(values)
    output('Provenance',{'sourcePins':pins,'case':case,'fold':fold_audit,
                        'lookups':sorted(lookups['lookups']),'closure':sorted(closure['closure'])},inputs)
    step('close repaired native case')
    model=c.build_case(inputs,*case)
    (RUN/'model.pickle').write_bytes(pickle.dumps(model,protocol=5))
    step('independent conservative power and native traction')
    covectors,pairing,residual=c.traction_pairing(inputs,case,model)
    conservative=c.conservative_power_variation(inputs,case)
    power_units=c.dimensions(pairing['SLAB_POWER'])
    output('PowerPairing',pairing,inputs,heavy=True)
    output('KineticNormalizationResidual',conservative['KINETIC_NORMALIZATION_RESIDUAL'],inputs,
           units=c.dimensions(conservative['KINETIC_SOURCE_ROWS']))
    output('ActionToRowMultiplier',conservative['ACTION_TO_ROW_MULTIPLIER'],inputs)
    output('PowerResidual',residual,inputs,units=power_units,heavy=True)
    step('canonical carrier comparison')
    normalized=canonical_residual(residual)
    output('CanonicalPowerResidual',normalized,inputs,units=power_units,heavy=True)
    output('CanonicalResidualFreeSymbols',tuple(sorted(normalized.free_symbols,key=sp.default_sort_key)),inputs)
    step('traction sign control')
    _,flipped,flipped_residual=c.traction_pairing(inputs,case,model,flip=True)
    output('TractionSignOperand',flipped,inputs,heavy=True)
    output('TractionSignResidual',c.difference(flipped_residual,residual),inputs,units=power_units,heavy=True)
    step('incoming face-work routing control')
    routed=c.build_case(inputs,*case,face_work_routing=-sp.S.One)
    _,routing_pairing,routing_residual=c.traction_pairing(inputs,case,routed)
    output('FaceRoutingOperand',routing_pairing,inputs,heavy=True)
    output('FaceRoutingResidual',c.difference(routing_residual,residual),inputs,units=power_units,heavy=True)
    (RUN/'objects.pickle').write_bytes(pickle.dumps((conservative,pairing,residual,normalized,
                                                   flipped_residual,routing_residual),protocol=5))
    manifest={'case':case,'sourcePins':pins,'sourcePinsAfter':{str(p):digest(p) for p in paths},
              'wallSeconds':time.monotonic()-STARTED,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              'canonicalResidualZero':normalized==0,'modelSha256':digest(RUN/'model.pickle'),
              'objectsSha256':digest(RUN/'objects.pickle')}
    (RUN/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    step('emission complete')


if __name__=='__main__':
    run()
