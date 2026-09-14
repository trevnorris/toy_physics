#!/usr/bin/env python3
"""Validate fresh full-subspace current/adjoint normalization and publication."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import pickle
import re
import shutil
from types import SimpleNamespace

import numpy as np
import sympy as sp

import S11c_d_end_normalization_check as runner
import S11c_d_end_normalization_remainder_check as remainder_checker
from S11c_d_end_normalization_check import ROOT, engine, digest, decoded_lines, _restore
from S11c_d_output_codec import restore_emission_index


def assoc(value):
    return {str(k):v for k,v in value}


def json_value(value):
    """Compare the in-memory value with its persisted JSON representation."""
    return json.loads(json.dumps(value))


def validate(base):
    summary = json.loads((base/'checks.json').read_text())
    validation_source = Path(__file__).resolve().relative_to(ROOT).as_posix()
    for name, sha in summary['sourceFiles'].items():
        if digest(base/'source'/name) != sha:
            raise ValueError(('frozen source pin', name))
        # A checker repair does not change the frozen calculation. Its new hash
        # is recorded separately; every current producer/helper pin stays exact.
        if name != validation_source and digest(ROOT/name) != sha:
            raise ValueError(('current producer source pin', name))
    for name, sha in summary['objects'].items():
        if digest(base/name) != sha:
            raise ValueError(('normalization payload pin', name))
    if (base/'stderr.txt').stat().st_size:
        raise ValueError('normalization run stderr')
    events = [json.loads(line) for line in (base/'progress.jsonl').read_text().splitlines()]
    if events[-1]['stage'] != 'completed':
        raise ValueError('normalization run incomplete')
    args = SimpleNamespace(**{k:Path(v) if v is not None and k!='end' else v
        for k,v in json.loads((base/'arguments.json').read_text()).items()})
    pairing, current, inputs, provenance = runner.load_pairing(args)
    native, coverage, binding, physical = runner.load_native(args,pairing,inputs,provenance)
    frequency = runner.load_frequency(args,provenance,physical)
    if json_value(provenance) != summary['provenance']:
        raise ValueError('normalization source/input joins')
    modal, known = pickle.loads((base/'modal.pickle').read_bytes())
    adjoint, adjoint_known = pickle.loads((base/'adjoint.pickle').read_bytes())
    engine.PHYSICAL_METADATA.dimensions.known.update(known)
    engine.PHYSICAL_METADATA.dimensions.known.update(adjoint_known)
    if modal['NATIVE_RECORDS'] != native or modal['NATIVE_COVERAGE'] != coverage:
        raise ValueError('native record and isolating disk join')
    entries = {}
    for line in decoded_lines(base/'full.out'):
        tag, sep, payload = line.rstrip('\n').partition(': ')
        if not sep or tag in entries or not tag.startswith('PY_S11CD_'):
            raise ValueError(('transcript tag', tag))
        entries[tag] = _restore(payload)
    prefix = summary['prefix']
    # Input mapping traverses a free-symbol set. Check every binding, then replay
    # its recorded association order so metadata grouping is reproducible.
    saved_bindings = assoc(entries['PY_S11CD_'+prefix+'_MATERIAL_AND_GRADE_BINDINGS'])
    names_to_symbols = {str(k):k for k in binding}
    if len(names_to_symbols) != len(binding) or saved_bindings != {str(k):v for k,v in binding.items()}:
        raise ValueError('material and grade binding census/value join')
    binding = {names_to_symbols[name]:binding[names_to_symbols[name]] for name in saved_bindings}
    final = 'PY_S11CD_'+prefix+'_EMISSION_LINES'
    if final not in entries:
        raise ValueError('missing emission index')
    indexed = restore_emission_index(assoc(entries[final]),list(entries)[:list(entries).index(final)])
    if set(indexed) != set(list(entries)[:list(entries).index(final)]):
        raise ValueError('emission index census')
    seen = set()
    original = engine.emit
    def check_emit(name, value):
        tag = 'PY_S11CD_'+name
        if tag in seen or entries.get(tag) != engine.cas(value):
            raise ValueError(('computed payload / metadata serialization', tag))
        seen.add(tag)
    engine.emit = check_emit
    try:
        builder = engine.ModalCurrentSubspaces(pairing,current,binding)
        builder.coefficient_residuals = modal['COEFFICIENT_RESIDUALS']
        runner.put(prefix,'EVALUATION_DOMAIN',{'END':args.end,'CASE':runner.CASE,'UNIT_FRAME':inputs.frame,
            'RETAINED_OPERATOR_FINITE_CONTRAST_EVALUATION':args.end!='REFERENCE',
            'CONTINUUM_REEXPANSION_PERFORMED':False,'GLOBAL_EXCEPTIONAL_COVERAGE_COMPUTED':False},pairing.modes)
        runner.put(prefix,'MATERIAL_AND_GRADE_BINDINGS',{str(k):v for k,v in binding.items()},pairing.modes,
            lambda p:engine.PHYSICAL_METADATA.dimensions.measure(next(k for k in binding if str(k)==p[0])))
        builder.emit(modal,provenance,prefix+'_MODAL',context=args.end+'_'+runner.CASE)
        joins = runner.frequency_joins(modal,frequency,builder,prefix)
        orientation = runner.orientations(modal,builder,args.end,prefix)
        engine.AdjointCurrentMap(builder).emit(adjoint,provenance,prefix+'_ADJOINT',context=args.end+'_'+runner.CASE)
        names = [tag.removeprefix('PY_S11CD_') for tag in entries if tag in seen and
                 not tag.startswith('PY_S11CD_METADATA_')]
        keys = {name:'s11cd'+''.join(v.title() for v in name.split('_')) for name in names}
        runner.put(prefix,'WRITE_KEYS',keys,pairing.modes)
        runner.put(prefix,'EMISSION_LINES',entries[final],pairing.modes)
    finally:
        engine.emit = original
    if seen != set(entries):
        raise ValueError(('unexpected/missing normalization tags', set(entries)-seen, seen-set(entries)))
    if len(keys) != len(set(keys.values())) or set(keys.values()) & set(engine.IMPORT_KEYS):
        raise ValueError('normalization key collision')
    paths = 0
    for name, body in entries.items():
        if name.startswith('PY_S11CD_METADATA_'):
            for value in body:
                descriptor = assoc(value)
                if len(descriptor['DIMENSION_L_T_M']) != 3 or any(v.free_symbols for v in descriptor['DIMENSION_L_T_M']):
                    raise ValueError(('unresolved dimension',name))
                if not all(k in descriptor for k in ('PATHS','MULTIGRADE','EPSILON_LAMBDA_SUPPORT')):
                    raise ValueError(('missing grades',name))
                paths += len(descriptor['PATHS'])
        elif name.endswith('_DIMENSION_CONSTRAINTS') and body:
            raise ValueError('nonempty dimension constraints')
    if joins != summary['frequencyJoins']:
        raise ValueError('independent frequency join summary')
    if json_value([{**v,'SHEET_MEMBERSHIP':str(v['SHEET_MEMBERSHIP'])} for v in orientation]) != summary['orientationRecords']:
        raise ValueError('outward orientation summary')
    if any(any(r['discreteResiduals'].values()) or max(r['residualNorms'].values())>1e-8 or
           max(r['liftDifferenceNorms'].values())>1e-8 for r in joins):
        raise ValueError('independent frequency/subspace discrepancy')
    if len(modal['RECORDS']) != len(native) or len(adjoint['RECORDS']) != len(native):
        raise ValueError('full native candidate census')
    residual_scalars = 0
    maxima = {}
    for mode, dual, old in zip(modal['RECORDS'],adjoint['RECORDS'],native):
        n = mode['NULLITY']
        for key in ('ROOT_DISK_INDEX','NORMAL_LIFT_SIGN','NULLITY'):
            if mode[key] != int(old[key]) or dual[key] != mode[key]:
                raise ValueError(('mode/adjoint/native candidate join',key))
        for side in ('RIGHT','LEFT'):
            if np.linalg.matrix_rank(mode['FORMS'][side],tol=1e-9) != n:
                raise ValueError('incomplete modal basis')
        for name, array in mode['RESIDUALS'].items():
            if not np.isfinite(array).all():
                raise ValueError(('nonfinite residual',name))
            norm = float(np.linalg.norm(array))
            if norm != mode['RESIDUAL_NORMS'][name]:
                raise ValueError('modal residual norm census')
            maxima['MODAL_'+name] = max(maxima.get('MODAL_'+name,0.),norm)
            residual_scalars += array.size
        for item in dual['ITEMS']:
            if not np.isfinite(item['VALUE']).all():
                raise ValueError(('nonfinite adjoint tensor',item['NAME']))
            if item['GROUP'] == 'RESIDUALS':
                name = item['NAME'];norm = float(np.linalg.norm(item['VALUE']))
                if norm != dual['RESIDUAL_NORMS'][name]:
                    raise ValueError('adjoint residual norm census')
                maxima['ADJOINT_'+name] = max(maxima.get('ADJOINT_'+name,0.),norm)
                residual_scalars += item['VALUE'].size
        if dual['INVERTIBLE_FIELD_MAP_DEFINED']:
            field = next(i['VALUE'] for i in dual['ITEMS'] if i['NAME']=='ADJOINT_FIELD')
            if np.linalg.matrix_rank(field,tol=1e-9) != n:
                raise ValueError('incomplete adjoint field basis')
    reality = runner.check_certificate(modal['NORMAL_REALITY_COVERAGE'],modal['RECORDS'])
    if reality != summary['normalReality']:
        raise ValueError('normal reality census')
    symbolic = [v for matrix in adjoint['SYMBOLIC_RESIDUALS'].values() for v in matrix]
    coefficients = [v for matrix in modal['COEFFICIENT_RESIDUALS'].values() for v in matrix]
    if any(v!=0 for v in symbolic+coefficients):
        raise ValueError('symbolic map/coefficient residual')
    diagnostic = {name:value for name,value in maxima.items() if value>1e-8}
    remainder = remainder_checker.compute(base,pairing,builder,modal,adjoint,binding)
    remainder_summary = remainder['summary']
    remainder_evidence = remainder_checker.emit_evidence(base,builder,remainder)
    (base/'remainders.pickle').write_bytes(pickle.dumps(remainder,protocol=5))
    if any(v['nonzeroScalars'] for v in remainder_summary['algebraicChecks'].values()):
        raise ValueError('pairing remainder decomposition/retained-grade discrepancy; inspect emitted evidence')
    if any(a<=1 and b<=1 for values in remainder_summary['coefficientGradesEtaSigma'].values() for a,b in values):
        raise ValueError('remainder contains retained background coefficient')
    if len(remainder_summary['records']) != len(native) or any(
            v['rightRank']!=v['nullity'] or v['adjointRank'] not in (None,v['nullity'])
            for v in remainder_summary['records']):
        raise ValueError('remainder full-subspace census')
    remainder_norms = remainder_summary['residualMinusRemainderNormMaxima']
    if any(value>1e-8 for value in remainder_norms.values()):
        raise ValueError('unexplained normalization balance residual; inspect emitted remainder operands')
    unaccounted = {name:value for name,value in diagnostic.items() if name not in remainder_norms}
    inventory = {**summary,'runDirectory':str(base),'tagCount':len(entries),
        'metadataPaths':paths,'numericResidualScalars':residual_scalars,
        'sourceAssignments':len(indexed),'validationSourceSha256':digest(Path(__file__)),
        'residualNormsAboveDiagnosticThreshold':diagnostic,
        'unaccountedResidualNormsAboveDiagnosticThreshold':unaccounted,
        'remainderAccounting':remainder_summary,'remainderEvidence':remainder_evidence,
        'diagnosticThreshold':1e-8,
        'artifacts':{name:{'bytes':(base/name).stat().st_size,'sha256':digest(base/name)} for name in
            ('full.out','stderr.txt','modal.pickle','adjoint.pickle','checks.json','progress.jsonl','arguments.json',
             'remainders.out','remainders.pickle')},
        'scope':'One supplied case, complete isolated root/lift subspaces, finite-contrast retained-operator evaluation; continuum re-expansion and global exceptional coverage remain open.'}
    (base/'validation.json').write_text(json.dumps(inventory,indent=2)+'\n')
    return inventory


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-directory',type=Path,required=True)
    parser.add_argument('--publish',action='store_true')
    parser.add_argument('--publication-suffix',default='thickness_repair')
    args = parser.parse_args();base = args.run_directory.resolve()
    if not re.fullmatch(r'[a-z][a-z0-9_]*',args.publication_suffix):
        raise ValueError('publication suffix')
    inventory = validate(base)
    print(json.dumps({key:inventory[key] for key in ('end','recordCount','basisDirections','definedFieldMaps',
        'physicalCurrentNormalizationCount','tagCount','metadataPaths','numericResidualScalars',
        'residualNormsAboveDiagnosticThreshold','unaccountedResidualNormsAboveDiagnosticThreshold')},indent=2))
    if inventory['unaccountedResidualNormsAboveDiagnosticThreshold']:
        raise ValueError('normalization diagnostics require investigation before publication')
    if args.publish:
        stem = 'S11c_d_end_normalization_'+inventory['end'].lower()+'_'+args.publication_suffix
        checkpoint = ROOT/'_measurements'/(stem+'_checkpoint.json')
        publications={'full.out':ROOT/'scripts/out'/(stem+'.out'),
                      'remainders.out':ROOT/'scripts/out'/(stem+'_remainders.out')}
        if checkpoint.exists() or any(target.exists() or target.is_symlink() for target in publications.values()):
            raise ValueError('publication target exists')
        inventory['publications']={}
        for name,target in publications.items():
            temporary = target.with_name('.'+target.name+'.new')
            with temporary.open('xb') as destination, (base/name).open('rb') as source:
                shutil.copyfileobj(source,destination);destination.flush();os.fsync(destination.fileno())
            if digest(temporary) != inventory['artifacts'][name]['sha256']:
                raise ValueError('publication copy hash')
            os.replace(temporary,target)
            inventory['publications'][name]={'path':str(target.relative_to(ROOT)),**inventory['artifacts'][name]}
        inventory['publication']=inventory['publications']['full.out']
        checkpoint.write_text(json.dumps(inventory,indent=2)+'\n')


if __name__ == '__main__':
    main()
