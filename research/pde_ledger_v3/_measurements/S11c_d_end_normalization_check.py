#!/usr/bin/env python3
"""Fresh endpoint/reference continuation of the native current and adjoint maps.

The input is a validated independent-frequency pairing, not a reference-mode
substitute. Full subspaces are joined to the independently computed frequency
packet. Finite-contrast evaluations retain their truncated-operator scope.
"""
import argparse
from collections import Counter
import faulthandler
import json
from pathlib import Path
import pickle
import resource
import time

import numpy as np
import sympy as sp

from S11c_d_end_pairing_check import build, ROOT, atomic
from S11c_d_modal_current_check import digest, source_node, engine
from S11c_d_joint_sheet_check import decoded_lines, _restore
from S11c_d_current_reality_check import check_certificate

CASE = 'LAB_HELD_RHO4_CONSTANT'
STACK_STREAM = None
SOURCES = ('scripts/S11c_d_mixing_scattering_sympy_audit.py',
    'scripts/S11c_d_output_codec.py', 'scripts/ledger_fold.py',
    '_measurements/S11c_d_end_normalization_check.py',
    '_measurements/S11c_d_end_normalization_validate.py',
    '_measurements/S11c_d_end_pairing_check.py',
    '_measurements/S11c_d_modal_current_check.py',
    '_measurements/S11c_d_joint_sheet_check.py',
    '_measurements/S11c_d_current_reality_check.py')


def checked_artifacts(checkpoint):
    base = Path(checkpoint['runDirectory']).resolve()
    for name, item in checkpoint['artifacts'].items():
        if digest(base/name) != item['sha256']:
            raise ValueError(('checkpoint artifact changed', str(base), name))
    for name, sha in checkpoint['sourceFiles'].items():
        if digest(base/'source'/name) != sha:
            raise ValueError(('checkpoint snapshot changed', name))
    return base


def load_pairing(args):
    checkpoint = json.loads(args.pairing_checkpoint.read_text())
    if checkpoint['end'] != args.end or checkpoint['retainedNonzeroScalars']:
        raise ValueError('matching endpoint pairing with no retained discrepancy required')
    base = checked_artifacts(checkpoint)
    pairing, inputs, provenance = build(args)
    old = checkpoint['provenance']
    # A new emitter or the downstream modal/adjoint classes may change; the
    # actual pairing and its source/input operands must remain identical.
    for key, value in provenance.items():
        if key in ('engineSha256', 'instrumentSha256'):
            continue
        if json.loads(json.dumps(value)) != old[key]:
            raise ValueError(('pairing producer join', key))
    frozen = (base/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py').read_text()
    for name in ('ClosedCurrentPairing', 'polynomial_terms'):
        if source_node(frozen, name) != source_node(Path(engine.__file__).read_text(), name):
            raise ValueError(('pairing construction changed', name))
    packet, known = pickle.loads((base/'complete.pickle').read_bytes())
    residual_leaves = {name:[v for _,v in engine.leaves(engine.cas(value))]
                       for name,value in packet['retained'].items()}
    counts = {name:sum(v != 0 for v in values) for name,values in residual_leaves.items()}
    if sum(counts.values()) or set(counts) != set(checkpoint['records']):
        raise ValueError('incomplete or nonzero retained pairing packet')
    if any(len(values)!=checkpoint['records'][name]['scalars'] for name,values in residual_leaves.items()):
        raise ValueError('retained pairing scalar census')
    engine.PHYSICAL_METADATA.dimensions.known.update(known)
    result = packet['result']
    def cached(anchoring, end):
        if (anchoring, end) != (pairing.anchoring, pairing.end):
            raise ValueError('pairing cache context mismatch')
        return result
    pairing.construct = cached
    provenance.update({'pairingCheckpointSha256':digest(args.pairing_checkpoint),
        'pairingCacheSha256':digest(base/'complete.pickle'),
        'pairingRetainedResidualScalars':checkpoint['retainedResidualScalars'],
        'pairingRetainedNonzeroScalars':sum(counts.values()),
        'instrumentSha256':digest(Path(__file__))})
    return pairing, result, inputs, provenance


def load_native(args, pairing, inputs, provenance):
    manifest = json.loads(args.manifest.read_text())
    base = Path(manifest['run_directory']).resolve()
    prefix = 'PY_S11CD_END_SPECTRUM_INPUT_'+args.end+'_'+CASE+'_0'
    packet, records = {}, []
    for line in decoded_lines(base/'full.out'):
        tag, _, payload = line.partition(': ')
        if not tag.startswith(prefix+'_'):
            continue
        name = tag[len(prefix)+1:]
        if name in ('ROOT_COVERAGE','BOUND_CARRIERS','GRADE_ORIGIN','PHYSICAL_PENCIL'):
            packet[name] = _restore(payload)
        elif name.startswith('MODE_') and name.endswith('_RECORD'):
            if int(name.split('_')[1]) != len(records):
                raise ValueError('native record sequence')
            records.append({str(k):v for k,v in _restore(payload)})
    coverage = {str(k):v for k,v in packet['ROOT_COVERAGE']}
    pairs = {(int(r['ROOT_DISK_INDEX']),int(r['NORMAL_LIFT_SIGN'])) for r in records}
    expected = {(i,s) for i in range(int(coverage['DISTINCT_ROOT_COUNT'])) for s in (-1,1)}
    if pairs != expected or len(records) != len(expected) or coverage['FINITE_POLYNOMIAL_ROOT_COVERAGE'] != sp.true:
        raise ValueError('native root/lift census unresolved')
    m = pairing.modes
    algebraic, relation, _ = m.analytic(pairing.acoustic.strong)
    mapping = inputs.mapping(algebraic, relation, (m.k,m.q,m.eta,m.sigma))
    origin = {m.eta:sp.S.Zero,m.sigma:sp.S.Zero} if args.end == 'REFERENCE' else inputs.origin
    if dict(packet['GRADE_ORIGIN']) != origin:
        raise ValueError('native grade origin differs')
    bound = {str(k):v for k,v in packet['BOUND_CARRIERS']}
    if any(bound.get(str(k)) != v for k,v in mapping.items()):
        raise ValueError('native material/profile bindings differ')
    mapping.update(origin)
    physical = algebraic.xreplace(mapping).applyfunc(sp.cancel)
    if engine.carrier_fingerprint(engine.cas(physical)) != packet['PHYSICAL_PENCIL']:
        raise ValueError('native and pairing physical-pencil fingerprint differs')
    provenance.update({'nativePacket':prefix,'rootDiskCount':int(coverage['DISTINCT_ROOT_COUNT']),
                       'rootLiftCount':len(records),'physicalPencilJoin':True})
    return records, coverage, mapping, physical


def load_frequency(args, provenance, physical):
    if args.end == 'REFERENCE':
        return None
    if args.frequency_checkpoint is None:
        raise ValueError('independent endpoint frequency checkpoint required')
    checkpoint = json.loads(args.frequency_checkpoint.read_text())
    base = checked_artifacts(checkpoint)
    old = checkpoint['summary']['provenance']
    for key in ('producerManifestSha256','cacheSha256','inputSha256','end','case'):
        if old[key] != provenance[key]:
            raise ValueError(('frequency checkpoint source/input join', key))
    frozen = (base/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py').read_text()
    if source_node(frozen,'EndModeFrequencyData') != source_node(Path(engine.__file__).read_text(),'EndModeFrequencyData'):
        raise ValueError('independent frequency constructor changed')
    result = pickle.loads((base/'objects.pickle').read_bytes())
    if result['PHYSICAL_PENCIL'] != physical:
        raise ValueError('independent frequency physical pencil differs')
    provenance.update({'frequencyCheckpointSha256':digest(args.frequency_checkpoint),
                       'frequencyCacheSha256':digest(base/'objects.pickle')})
    return result


def put(prefix, name, value, modes, unit=lambda path:(0,0,0), heavy=False):
    if isinstance(value,np.ndarray):
        value = sp.ImmutableMatrix(*value.shape,[modes.number(v) for v in value.ravel()])
    body = engine.cas(value)
    engine.emit(prefix+'_'+name,modes.compact_fingerprint(body) if heavy else body)
    engine.emit('METADATA_'+prefix+'_'+name,modes.numeric_metadata(body,unit))


def frequency_joins(modal, source, builder, prefix):
    if source is None:
        return []
    if len(modal['RECORDS']) != len(source['RECORDS']):
        raise ValueError('frequency/modal candidate census')
    packets = {p['NAME']:p['VALUE'] for p in source['PACKETS']}
    outputs = []
    for record, old in zip(modal['RECORDS'],source['RECORDS']):
        i = record['INDEX']; n = record['NULLITY']
        context = prefix+'_FREQUENCY_JOIN_'+str(i)
        discrete = {key:record[key]-old[key] for key in
                    ('INDEX','ROOT_DISK_INDEX','NORMAL_LIFT_SIGN','NULLITY','FREQUENCY_PAIRING_RANK')}
        numeric = {}
        for side in ('RIGHT','LEFT'):
            before = np.array(packets[f'MODE_{i}_{side}_BASIS'],dtype=complex)
            after = record['FORMS'][side]
            p0, p1 = before@before.conj().T, after@after.conj().T
            numeric[side+'_PROJECTOR'] = p1-p0
            # Coordinate projectors live in the numerical L/T/M unit frame;
            # restore each row/column's field or reciprocal-row unit.
            units = builder.field_units if side == 'RIGHT' else [tuple(-v for v in u) for u in builder.row_units]
            unit = lambda path:tuple(a-b for a,b in zip(units[path[0]//5],units[path[0]%5]))
            put(context,side+'_SOURCE_PROJECTOR',p0,builder.modes,unit,True)
            put(context,side+'_CURRENT_PROJECTOR',p1,builder.modes,unit,True)
            put(context,side+'_PROJECTOR_RESIDUAL',p1-p0,builder.modes,unit)
            discrete[side+'_BASIS_RANK'] = int(np.linalg.matrix_rank(after,tol=1e-9))-n
        before = np.array(packets[f'MODE_{i}_FREQUENCY_MATRIX'],dtype=complex)
        after = record['OPERANDS']['FREQUENCY_PENCIL_PLUS']
        numeric['FREQUENCY_MATRIX'] = after-before
        unit = lambda path:tuple(a-b-c for a,b,c in zip(builder.row_units[path[0]//5],builder.field_units[path[0]%5],builder.frequency_unit))
        put(context,'SOURCE_FREQUENCY_MATRIX',before,builder.modes,unit,True)
        put(context,'CURRENT_FREQUENCY_MATRIX',after,builder.modes,unit,True)
        put(context,'FREQUENCY_MATRIX_RESIDUAL',after-before,builder.modes,unit)
        put(context,'DISCRETE_RESIDUALS',discrete,builder.modes)
        # The low-degree exact-lift recovery and native numeric root must
        # identify the same isolated root; raw differences remain visible.
        coordinate = {'K':record['K']-old['K'],'Q':record['Q']-old['Q']}
        put(context,'LIFT_RESIDUALS',coordinate,builder.modes,
            lambda p:tuple(-v for v in builder.length_unit) if p[-1]=='K' else builder.frequency_unit)
        outputs.append({'index':i,'discreteResiduals':discrete,
            'residualNorms':{key:float(np.linalg.norm(v)) for key,v in numeric.items()},
            'liftDifferenceNorms':{key:float(abs(complex(v))) for key,v in coordinate.items()}})
    return outputs


def orientations(modal, builder, end, prefix):
    records = []
    orientation = {'LEFT':-1,'RIGHT':1}.get(end)
    put(prefix,'OUTWARD_END_ORIENTATION',{'END':end,'DEFINED':orientation is not None,
        'VALUE':() if orientation is None else orientation},builder.modes)
    for source in modal['RECORDS']:
        defined = orientation is not None and source.get('PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED',False)
        record = {'INDEX':source['INDEX'],'ROOT_DISK_INDEX':source['ROOT_DISK_INDEX'],
            'NORMAL_LIFT_SIGN':source['NORMAL_LIFT_SIGN'],'NULLITY':source['NULLITY'],
            'OUTWARD_FLUX_CLASSIFICATION_DEFINED':defined,'SHEET_MEMBERSHIP':source['SHEET_MEMBERSHIP'],
            'RETAINED_AS_MATCHING_CANDIDATE':True}
        if defined:
            signed = source['FORMS']['SIGNED_CURRENT']
            outward = orientation*signed
            incoming = np.flatnonzero(np.real(np.diag(outward))<0)
            outgoing = np.flatnonzero(np.real(np.diag(outward))>0)
            record.update({'INCOMING_BASIS_COLUMNS':tuple(map(int,incoming)),
                'OUTGOING_BASIS_COLUMNS':tuple(map(int,outgoing)),
                'UNCLASSIFIED_BASIS_DIRECTIONS':source['NULLITY']-len(incoming)-len(outgoing)})
            tag = prefix+'_ORIENTATION_'+str(source['INDEX'])
            put(tag,'NORMAL_SIGNED_CURRENT',builder.epsilon**2*sp.ImmutableMatrix(signed),builder.modes)
            put(tag,'OUTWARD_SIGNED_CURRENT',builder.epsilon**2*sp.ImmutableMatrix(outward),builder.modes)
        put(prefix,'ORIENTATION_'+str(source['INDEX'])+'_RECORD',record,builder.modes)
        records.append(record)
    return records


def main():
    global STACK_STREAM
    parser = argparse.ArgumentParser()
    for name in ('manifest','input','pairing-checkpoint','run-directory'):
        parser.add_argument('--'+name,type=Path,required=True)
    for name in ('source-checkpoint','current-manifest','frequency-checkpoint'):
        parser.add_argument('--'+name,type=Path)
    parser.add_argument('--end',choices=('LEFT','RIGHT','REFERENCE'),required=True)
    args = parser.parse_args()
    base = args.run_directory.resolve()
    base.relative_to(ROOT.parents[1]/'_scratch/s11c')
    base.mkdir(parents=True,exist_ok=False)
    started = time.monotonic()
    STACK_STREAM = (base/'stack-samples.txt').open('a')
    faulthandler.enable()
    faulthandler.dump_traceback_later(120,repeat=True,file=STACK_STREAM)
    pins = {name:digest(ROOT/name) for name in SOURCES}
    for name in SOURCES:
        target = base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes((ROOT/name).read_bytes())
    def progress(record):
        with (base/'progress.jsonl').open('a') as stream:
            stream.write(json.dumps({'elapsedSeconds':time.monotonic()-started,
                'cpuSeconds':time.process_time(),
                'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,**record})+'\n')
    (base/'arguments.json').write_text(json.dumps({k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},indent=2)+'\n')
    pairing, current, inputs, provenance = load_pairing(args)
    records, coverage, binding, physical = load_native(args,pairing,inputs,provenance)
    frequency = load_frequency(args,provenance,physical)
    prefix = 'END_NORMALIZATION_'+args.end+'_'+CASE
    builder = engine.ModalCurrentSubspaces(pairing,current,binding)
    put(prefix,'EVALUATION_DOMAIN',{'END':args.end,'CASE':CASE,'UNIT_FRAME':inputs.frame,
        'RETAINED_OPERATOR_FINITE_CONTRAST_EVALUATION':args.end!='REFERENCE',
        'CONTINUUM_REEXPANSION_PERFORMED':False,'GLOBAL_EXCEPTIONAL_COVERAGE_COMPUTED':False},pairing.modes)
    put(prefix,'MATERIAL_AND_GRADE_BINDINGS',{str(k):v for k,v in binding.items()},pairing.modes,
        lambda p:engine.PHYSICAL_METADATA.dimensions.measure(next(k for k in binding if str(k)==p[0])))
    progress({'stage':'loaded','provenance':provenance})
    modal = builder.construct(records,coverage,float(inputs.parameters['omega']),float(inputs.parameters['W_0']),progress)
    modal['COEFFICIENT_RESIDUALS'] = builder.coefficient_residuals
    modal['SYMBOLIC_OPERANDS'] = builder.symbolic_operands
    modal['SCALAR_OPERANDS'] = builder.scalar_operands
    modal['NATIVE_RECORDS'] = records
    modal['NATIVE_COVERAGE'] = coverage
    # Keep a construction packet before potentially lengthy serialization. The
    # established post-emission packet still includes the final unit registry.
    atomic(base/'modal-pre-emission.pickle',pickle.dumps(
        (modal,engine.PHYSICAL_METADATA.dimensions.known),protocol=5))
    progress({'stage':'modal_emission_started','constructionPacketSha256':digest(base/'modal-pre-emission.pickle')})
    builder.emit(modal,provenance,prefix+'_MODAL',context=args.end+'_'+CASE)
    (base/'modal.pickle').write_bytes(pickle.dumps((modal,engine.PHYSICAL_METADATA.dimensions.known),protocol=5))
    progress({'stage':'modal_saved','sha256':digest(base/'modal.pickle')})
    joins = frequency_joins(modal,frequency,builder,prefix)
    oriented = orientations(modal,builder,args.end,prefix)
    adjoint_builder = engine.AdjointCurrentMap(builder)
    adjoint = adjoint_builder.construct(modal,progress)
    atomic(base/'adjoint-pre-emission.pickle',pickle.dumps(
        (adjoint,engine.PHYSICAL_METADATA.dimensions.known),protocol=5))
    progress({'stage':'adjoint_emission_started','constructionPacketSha256':digest(base/'adjoint-pre-emission.pickle')})
    adjoint_builder.emit(adjoint,provenance,prefix+'_ADJOINT',context=args.end+'_'+CASE)
    (base/'adjoint.pickle').write_bytes(pickle.dumps((adjoint,engine.PHYSICAL_METADATA.dimensions.known),protocol=5))
    reality = check_certificate(modal['NORMAL_REALITY_COVERAGE'],modal['RECORDS'])
    residuals = lambda values:{key:max(v.get('RESIDUAL_NORMS',{}).get(key,0) for v in values)
        for key in sorted(set().union(*(v.get('RESIDUAL_NORMS',{}) for v in values)))}
    symbolic = [v for matrix in adjoint['SYMBOLIC_RESIDUALS'].values() for v in matrix]
    coefficients = [v for matrix in builder.coefficient_residuals.values() for v in matrix]
    summary = {'end':args.end,'case':CASE,'prefix':prefix,'provenance':provenance,'sourceFiles':pins,
        'recordCount':len(modal['RECORDS']),'basisDirections':sum(v['NULLITY'] for v in modal['RECORDS']),
        'nullityCounts':dict(Counter(v['NULLITY'] for v in modal['RECORDS'])),
        'definedFieldMaps':sum(v['INVERTIBLE_FIELD_MAP_DEFINED'] for v in adjoint['RECORDS']),
        'physicalCurrentNormalizationCount':sum(v.get('PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED',False) for v in modal['RECORDS']),
        'coefficientResidualScalars':len(coefficients),'coefficientNonzeroScalars':sum(v!=0 for v in coefficients),
        'symbolicResidualScalars':len(symbolic),'symbolicNonzeroScalars':sum(v!=0 for v in symbolic),
        'modalResidualNormMaxima':residuals(modal['RECORDS']),
        'adjointResidualNormMaxima':residuals(adjoint['RECORDS']),
        'frequencyJoins':joins,'orientationRecords':[{**v,'SHEET_MEMBERSHIP':str(v['SHEET_MEMBERSHIP'])} for v in oriented],
        'normalReality':reality,'wallSeconds':time.monotonic()-started,
        'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'objects':{name:digest(base/name) for name in
            ('modal.pickle','adjoint.pickle','modal-pre-emission.pickle','adjoint-pre-emission.pickle')}}
    # Map identities and finite-contrast physical balance diagnostics remain
    # separate. A nonzero raw balance cannot be silently accepted as closure;
    # retained/remainder accounting stays in the linked pairing checkpoint.
    (base/'checks.json').write_text(json.dumps(summary,indent=2)+'\n')
    names = [tag.removeprefix('PY_S11CD_') for tag in engine.EMISSION_LINES
             if not tag.startswith('PY_S11CD_METADATA_')]
    keys = {name:'s11cd'+''.join(v.title() for v in name.split('_')) for name in names}
    put(prefix,'WRITE_KEYS',keys,pairing.modes)
    index = engine.emission_index(engine.EMISSION_LINES)
    put(prefix,'EMISSION_LINES',index,pairing.modes)
    progress({'stage':'completed','summary':{k:summary[k] for k in ('recordCount','basisDirections','definedFieldMaps','wallSeconds')}})
    if pins != {name:digest(ROOT/name) for name in SOURCES}:
        raise ValueError('normalization source changed during run')
    if len(set(keys.values())) != len(keys) or set(keys.values()) & set(engine.IMPORT_KEYS):
        raise ValueError('normalization write-key collision')
    if summary['coefficientNonzeroScalars'] or summary['symbolicNonzeroScalars']:
        raise ValueError('symbolic reconstruction discrepancy; inspect emitted operands')
    if any(any(r['discreteResiduals'].values()) or max(r['residualNorms'].values())>1e-8 or
           max(r['liftDifferenceNorms'].values())>1e-8 for r in joins):
        raise ValueError('independent frequency/subspace join discrepancy; inspect emitted operands')


if __name__ == '__main__':
    try:
        main()
    finally:
        faulthandler.cancel_dump_traceback_later()
        if STACK_STREAM is not None:
            STACK_STREAM.close()
