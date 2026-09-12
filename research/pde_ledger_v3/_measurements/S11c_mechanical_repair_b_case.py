#!/usr/bin/env python3
"""Compute one native repaired b case and compare the pinned previous export."""
from pathlib import Path
import hashlib
import json
import pickle
import resource
import sys
import time

STARTED = time.monotonic()
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import sympy as sp
import S11c_b_brane_operator_sympy_audit as b
import S11c_d_mixing_scattering_sympy_audit as fingerprints
from ledger_fold import _restore
from S11c_inertia_artifact_audit import export_data

BASE = Path('/tmp/s11c-mechanical-repair-20260912')
RUN = BASE/'b_case'
RUN.mkdir(exist_ok=True)
CASE = ('LAB_HELD', 'RHO4_CONSTANT')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(operator):
    return sp.ImmutableMatrix([*b.named_tuple_row(b.named_tuple_row(operator, 'U_BODY_BALANCE'), 'EXPANDED'),
                               b.named_tuple_row(b.named_tuple_row(operator, 'E_W_BALANCE'), 'EXPANDED')])


def output(name, value, units=None, heavy=False):
    record = {'writeKey':'s11cMechanicalRepair'+name,
              'value':sp.srepr(fingerprints.carrier_fingerprint(fingerprints.cas(value)) if heavy else b.casify(value)),
              'epsEtaSigmaOrder':sp.srepr(b.multigrade(value)),
              'dimensionLTM':sp.srepr(b.casify(units if units is not None else b.dimension_object(value)))}
    print(json.dumps(record,separators=(',',':')),flush=True)


def run():
    pins = {name:digest(ROOT/name) for name in (
        'scripts/S11c_b_brane_operator_sympy_audit.py','scripts/S11c_a_exports.py',
        'directives/S11c_b_SHARED_PHYSICS.md','_measurements/S11c_mechanical_repair_b_case.py')}
    old_values,_,_ = export_data(BASE/'baseline/scripts/S11c_b_exports.py')
    before_cases = {tuple(map(str,axes)):b.named_tuple_row(payload,'VALUE')
                    for axes,payload in _restore(old_values['slab_operator'])}
    before = before_cases[CASE]
    output('Case',CASE)
    operator, origins, mu = b.build_operator(*CASE)
    operator, origins, mu = b.retained_grade(operator), b.retained_grade(origins), b.retained_grade(mu)
    (RUN/'objects.pickle').write_bytes(pickle.dumps((operator,origins,mu),protocol=5))
    after = b.casify(operator)
    old_rows,new_rows = rows(before),rows(after)
    source = b.named_tuple_row(after,'FACE_GENERALIZED_FORCE_ROWS')
    force = sp.ImmutableMatrix([*b.named_tuple_row(source,'U'),b.named_tuple_row(source,'E_W')])
    normalization=b.named_tuple_row(b.named_tuple_row(b.casify(origins),'FACE_VIRTUAL_WORK'),'ROW_NORMALIZATION')
    multiplier=b.named_tuple_row(normalization,'ACTION_TO_STORED_ROW_MULTIPLIER')
    control={b.INCOMING_LEDGER[name]['value']:sp.S.Zero for name in
             ('delta_p_plus','delta_p_minus','d_w_delta_p_plus','d_w_delta_p_minus','Lambda_X_0')}
    increment=(new_rows-new_rows.xreplace(control)).applyfunc(sp.expand)
    units=sp.Tuple(*[b.DIM_BODY_U]*3,b.DIM_ENERGY)
    output('RowNormalization',normalization)
    output('BeforeRows',old_rows,units,True)
    output('AfterRows',new_rows,units,True)
    output('ExternalForce',force,units,True)
    output('MechanicalFaceIncrement',increment,units,True)
    output('ActionLoad',multiplier*force,units,True)
    output('ActionLoadResidual',(increment-multiplier*force).applyfunc(sp.cancel),units)
    output('NonFacePreservationResidual',(new_rows.xreplace(control)-old_rows.xreplace(control)).applyfunc(sp.cancel),units)
    before_mass=b.named_tuple_row(b.named_tuple_row(before,'THETA_BALANCE'),'EXPANDED')
    after_mass=b.named_tuple_row(b.named_tuple_row(after,'THETA_BALANCE'),'EXPANDED')
    output('MassPreservationResidual',sp.cancel(after_mass-before_mass),b.DIM_RHOBR-b.DIM_T)
    old_force=b.named_tuple_row(before,'FACE_GENERALIZED_FORCE_ROWS')
    output('PhysicalForcePreservationResidual',sp.ImmutableMatrix([
        *[sp.cancel(x-y) for x,y in zip(b.named_tuple_row(source,'U'),b.named_tuple_row(old_force,'U'))],
        sp.cancel(b.named_tuple_row(source,'E_W')-b.named_tuple_row(old_force,'E_W'))]),units)
    old_mu=b.named_tuple_row(before,'MU_THETA_FACE_BINDING')[1]
    new_mu=b.named_tuple_row(after,'MU_THETA_FACE_BINDING')[1]
    output('ChemicalPreservationResidual',sp.cancel(new_mu-old_mu),b.DIM_ENERGY)
    output('SourcePins',pins)
    resources={'wallSeconds':time.monotonic()-STARTED,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    manifest={'case':CASE,'sourcePins':pins,'sourcePinsAfter':{name:digest(ROOT/name) for name in pins},
              'resources':resources,'cacheSha256':digest(RUN/'objects.pickle')}
    (RUN/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    output('WallSeconds',sp.Float(resources['wallSeconds']),b.DIM_T)
    output('PeakRssKiB',resources['peakRssKiB'],b.DIM_ZERO)


if __name__=='__main__':
    run()
