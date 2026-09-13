#!/usr/bin/env python3
"""Source-pinned both-end frequency data without the physical-current dependency."""
import argparse
import ast
from collections import Counter
import json
from pathlib import Path
import pickle
import resource
import time

import sympy as sp
from S11c_d_joint_sheet_check import load, digest, engine, _restore, decoded_lines
from S11c_d_current_reality_check import check_certificate

ROOT = Path(engine.__file__).resolve().parents[1]
SOURCES = ('scripts/S11c_d_mixing_scattering_sympy_audit.py','scripts/S11c_d_output_codec.py',
    'scripts/ledger_fold.py','_measurements/S11c_d_end_frequency_check.py',
    '_measurements/S11c_d_end_frequency_validate.py','_measurements/S11c_d_joint_sheet_check.py',
    '_measurements/S11c_d_current_reality_check.py','_measurements/S11c_d_channel_preflight_input.json')


def method(source, cls, name):
    node=next(v for v in ast.parse(source).body if isinstance(v,ast.ClassDef) and v.name==cls)
    return ast.dump(next(v for v in node.body if isinstance(v,ast.FunctionDef) and v.name==name),include_attributes=False)


def load_packet(args):
    modes,strong,units,branches,inputs,provenance=load(args)
    manifest=json.loads(args.manifest.read_text());base=Path(manifest['run_directory'])
    frozen=(base/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py').read_text()
    current=Path(engine.__file__).read_text()
    joins={cls+'.'+name:method(frozen,cls,name)==method(current,cls,name)
        for cls,names in (('FullPencilModes',('analytic','rational_determinant')),
                          ('EndSpectrumCoverage',('construct','isolate'))) for name in names}
    if not all(joins.values()):raise ValueError(('native spectrum computation source changed',joins))
    prefix='PY_S11CD_END_SPECTRUM_INPUT_'+args.end+'_LAB_HELD_RHO4_CONSTANT_0'
    packet={};records=[]
    for line in decoded_lines(base/'full.out'):
        tag,_,payload=line.partition(': ')
        if not tag.startswith(prefix+'_'):continue
        name=tag[len(prefix)+1:]
        if name in ('ROOT_COVERAGE','BOUND_CARRIERS','GRADE_ORIGIN','PHYSICAL_PENCIL','ELIMINATION_OPERANDS'):
            packet[name]=_restore(payload)
        elif name.startswith('MODE_') and name.endswith('_RECORD'):
            if int(name.split('_')[1])!=len(records):raise ValueError('nonsequential native record')
            records.append({str(k):v for k,v in _restore(payload)})
    required={'ROOT_COVERAGE','BOUND_CARRIERS','GRADE_ORIGIN','PHYSICAL_PENCIL','ELIMINATION_OPERANDS'}
    if set(packet)!=required or not records:raise ValueError('incomplete native end packet')
    coverage={str(k):v for k,v in packet['ROOT_COVERAGE']}
    lifts={(int(v['ROOT_DISK_INDEX']),int(v['NORMAL_LIFT_SIGN'])) for v in records}
    expected={(i,s) for i in range(int(coverage['DISTINCT_ROOT_COUNT'])) for s in (-1,1)}
    if lifts!=expected or len(records)!=len(lifts) or coverage['FINITE_POLYNOMIAL_ROOT_COVERAGE']!=sp.true:
        raise ValueError('unresolved or duplicate native disk/lift coverage')
    algebraic,relation,_=modes.analytic(strong)
    mapping=inputs.mapping(algebraic,relation,(modes.k,modes.q,modes.eta,modes.sigma))
    bound={str(k):v for k,v in packet['BOUND_CARRIERS']}
    if any(bound.get(str(k))!=v for k,v in mapping.items()):raise ValueError('native end material/limit binding changed')
    if dict(packet['GRADE_ORIGIN'])!=inputs.origin:raise ValueError('native end grade origin changed')
    provenance.update({'nativePacket':prefix,'sourceMethodJoins':joins,
        'profileLimits':{str(k):str(v) for k,v in inputs.limits.items()},
        'instrumentSha256':digest(Path(__file__))})
    return modes,strong,units,inputs,provenance,records,coverage,packet


def main():
    parser=argparse.ArgumentParser()
    for name in ('manifest','input','run-directory'):parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--end',choices=('LEFT','RIGHT'),required=True)
    args=parser.parse_args();args.run_directory.mkdir(parents=True,exist_ok=True)
    started=time.monotonic()
    source_hashes={name:digest(ROOT/name) for name in SOURCES}
    for name in SOURCES:
        path=args.run_directory/'source'/name;path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes((ROOT/name).read_bytes())
    def progress(record):
        with (args.run_directory/'progress.jsonl').open('a') as stream:
            stream.write(json.dumps({'elapsedSeconds':time.monotonic()-started,**record})+'\n')
    modes,strong,units,inputs,provenance,records,coverage,native=load_packet(args)
    progress({'stage':'loaded','provenance':provenance})
    builder=engine.EndModeFrequencyData(modes,strong,units,inputs)
    result=builder.construct(records,coverage,progress)
    native_joins={
        'PHYSICAL_PENCIL_FINGERPRINT':engine.carrier_fingerprint(engine.cas(result['PHYSICAL_PENCIL']))==native['PHYSICAL_PENCIL'],
        'ELIMINATION_FINGERPRINT':engine.carrier_fingerprint(engine.cas(result['ELIMINATION']))==native['ELIMINATION_OPERANDS']}
    builder.put('NATIVE_PACKET_JOINS',native_joins)
    prefix='END_FREQUENCY_'+args.end+'_LAB_HELD_RHO4_CONSTANT'
    builder.emit(result,provenance,prefix)
    result['KNOWN_DIMENSIONS']=engine.PHYSICAL_METADATA.dimensions.known
    objects=args.run_directory/'objects.pickle';objects.write_bytes(pickle.dumps(result,protocol=5))
    index=engine.emission_index(engine.EMISSION_LINES)
    engine.emit(prefix+'_EMISSION_LINES',index)
    engine.emit('METADATA_'+prefix+'_EMISSION_LINES',modes.numeric_metadata(engine.cas(index),lambda p:(0,0,0)))
    reality=check_certificate(result['NORMAL_REALITY_COVERAGE'],result['RECORDS'])
    exact=[v for p in result['PACKETS'] if p['NAME'] in ('SOURCE_BRANCH_JOIN_RESIDUAL','FREQUENCY_WAVE_TANGENCY_RESIDUAL',
        'FREQUENCY_BINDING_RESIDUAL','NORMAL_REALITY_CHECKS') for _,v in engine.leaves(p['VALUE']) if not isinstance(v,engine.Str)]
    residuals={key:max(v['COEFFICIENT_FRAME_RESIDUAL_NORMS'].get(key,0) for v in result['RECORDS'])
        for key in sorted(set().union(*(v['COEFFICIENT_FRAME_RESIDUAL_NORMS'] for v in result['RECORDS'])))}
    summary={'provenance':provenance,'sourceFiles':source_hashes,'arguments':{k:str(v) for k,v in vars(args).items()},
        'prefix':prefix,'nativePacketJoins':native_joins,'recordCount':len(result['RECORDS']),
        'nullityCounts':dict(Counter(v.get('NULLITY',0) for v in result['RECORDS'])),
        'basisDirections':sum(v.get('NULLITY',0) for v in result['RECORDS']),
        'frequencyNormalizedCount':sum(v['FREQUENCY_NORMALIZATION_DEFINED'] for v in result['RECORDS']),
        'exactResidualScalars':len(exact),'nonzeroExactResidualScalars':sum(v!=0 for v in exact),
        'normalReality':reality,'residualNormMaxima':residuals,'objectsSha256':digest(objects),
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (args.run_directory/'checks.json').write_text(json.dumps(summary,indent=2)+'\n')
    progress({'stage':'completed','summary':summary})
    if source_hashes!={name:digest(ROOT/name) for name in SOURCES}:raise ValueError('source changed during run')
    if not all(native_joins.values()) or summary['nonzeroExactResidualScalars']:
        raise ValueError('source join residual; inspect emitted operands')


if __name__=='__main__':main()
