#!/usr/bin/env python3
"""Emit the computed endpoint pairing, raw/retained residuals and remainders."""
import argparse
import json
from pathlib import Path
import pickle
import resource
import time

import sympy as sp
from S11c_d_end_pairing_check import build, atomic, ROOT
from S11c_d_modal_current_check import digest, engine


def small_literal(value, maximum_nodes=1200):
    """Bound the printed expression tree, counting repeated shared nodes."""
    pending=[value]
    for _ in range(maximum_nodes):
        if not pending:return True
        pending.extend(pending.pop().args)
    return not pending


def diagnostic_record(base):
    """Bind optional scalar-cache instrumentation to the published provenance."""
    state=base/'scalar-state'
    if not state.exists():return None
    signature=json.loads((state/'signature.json').read_text())
    instrument=Path(signature['instrument'])
    if digest(instrument)!=signature['instrumentSha256'] or digest(state/instrument.name)!=signature['instrumentSha256']:
        raise ValueError('diagnostic instrument pin')
    if digest(base/'signature.json')!=signature['checkerSignatureSha256']:
        raise ValueError('diagnostic checker signature pin')
    operations=[json.loads(line) for line in (state/'operations.jsonl').read_text().splitlines()]
    if not operations or operations[-1]['stage']!='complete':
        raise ValueError('diagnostic replay is incomplete')
    packets={}
    for operation in operations:
        if operation['stage']!='saved':continue
        stem=operation['method']+'_'+operation['sourceSha256']
        packet=state/(stem+'.pickle');source=state/(stem+'.input.pickle')
        sha=digest(packet)
        if sha!=operation['sha256'] or sha!=(state/(stem+'.sha256')).read_text().strip():
            raise ValueError(('diagnostic result pin',stem))
        if digest(source)!=operation['sourceSha256']:
            raise ValueError(('diagnostic operand pin',stem))
        packets[stem]=sha
    return {'signature':signature,'signatureSha256':digest(state/'signature.json'),
            'operationLogSha256':digest(state/'operations.jsonl'),'operationCount':len(operations),
            'scalarPackets':packets}


def worker_record(base):
    """Read and verify the completed sequential-worker provenance."""
    state=base/'worker-state'
    if not state.exists():return None
    signature=json.loads((state/'signature.json').read_text())
    source=ROOT/'_measurements/S11c_d_end_pairing_workers.py'
    if digest(source)!=signature['environment']['workerSha256'] or digest(state/source.name)!=digest(source):
        raise ValueError('worker instrument source pin')
    if digest(Path(engine.__file__))!=signature['environment']['nativeSha256']:
        raise ValueError('worker native source pin')
    if digest(base/'signature.json')!=signature['checkerSignatureSha256']:
        raise ValueError('worker checker signature pin')
    seed=Path(signature['seedRun'])
    if digest(seed/'signature.json')!=signature['seedSignatureSha256']:
        raise ValueError('worker seed signature pin')
    operations=[json.loads(line) for line in (state/'operations.jsonl').read_text().splitlines()]
    if not operations or operations[-1]['stage']!='complete' or not operations[-1]['environmentStable']:
        raise ValueError('worker calculation incomplete or runtime changed')
    # This read-only method checks the complete request/input/result round trip.
    from S11c_d_end_pairing_workers import ScalarWorkers
    reader=ScalarWorkers.__new__(ScalarWorkers);reader.signature=signature
    workers={};seeds={}
    for operation in operations:
        if operation['stage']=='worker_saved':
            directory=Path(operation['directory'])
            _,_,record=reader.read_result(directory)
            if record!=operation['resource']:raise ValueError('worker operation/resource join')
            workers[directory.name]=record
        if operation['stage']=='reused':
            provenance=operation['provenance']
            if provenance['kind']=='diagnostic_seed':
                if digest(Path(provenance['packet']))!=provenance['sha256'] or digest(Path(provenance['input']))!=provenance['inputSha256']:
                    raise ValueError('reused scalar seed pin')
                before,_,_=pickle.loads(Path(provenance['packet']).read_bytes())
                if before!=pickle.loads(Path(provenance['input']).read_bytes()):raise ValueError('reused scalar operand join')
                seeds[provenance['packet']]=provenance
            elif provenance['kind']=='isolated_worker':
                _,_,record=reader.read_result(Path(provenance['directory']))
                if record['resultSha256']!=provenance['resultSha256']:raise ValueError('reused worker pin')
            else:raise ValueError('unknown scalar cache provenance')
    boundary=json.loads((state/'boundary.json').read_text())
    for item in boundary:
        path=Path(item['comparisonPacket'])
        if digest(path)!=item['comparisonSha256']:raise ValueError('worker boundary packet pin')
        _,before,after,residual=pickle.loads(path.read_bytes())
        if before!=after or residual!=after-before:raise ValueError('worker boundary comparison')
    return {'signature':signature,'signatureSha256':digest(state/'signature.json'),
            'operationLogSha256':digest(state/'operations.jsonl'),'operationCount':len(operations),
            'boundary':boundary,'workers':workers,'reusedDiagnosticPackets':seeds}


def run():
    parser=argparse.ArgumentParser()
    for name in ('manifest','input','run-directory'):parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--end',choices=('REFERENCE','LEFT','RIGHT'),required=True)
    parser.add_argument('--source-checkpoint',type=Path)
    parser.add_argument('--current-manifest',type=Path)
    args=parser.parse_args();base=args.run_directory;started=time.monotonic()
    summary=json.loads((base/'checks.json').read_text())
    if digest(base/'complete.pickle')!=summary['objectsSha256']:raise ValueError('pairing packet digest')
    for name,sha in summary['sourceFiles'].items():
        if digest(ROOT/name)!=sha:raise ValueError(('pairing source changed',name))
    builder,inputs,provenance=build(args)
    if json.loads(json.dumps(provenance))!=summary['provenance']:raise ValueError('pairing source/input provenance differs')
    packet,known=pickle.loads((base/'complete.pickle').read_bytes())
    d=engine.PHYSICAL_METADATA.dimensions;d.known.update(known)
    result,checks=packet['result'],packet['checks']
    diagnostic=diagnostic_record(base)
    workers=worker_record(base)
    prefix='END_CURRENT_PAIRING_'+args.end+'_LAB_HELD_RHO4_CONSTANT'
    records=[]
    def output(name,value,units=None,metadata_operand=None):
        body=engine.cas(value)
        source=body if metadata_operand is None else engine.cas(metadata_operand)
        units={} if units is None else units
        metadata=[]
        for path,leaf in engine.leaves(source):
            if isinstance(leaf,engine.Str):continue
            unit=units.get(path)
            if unit is None:unit=d.measure(leaf)
            if unit is None:raise ValueError(('untyped zero pairing leaf',name,path))
            try:
                engine.PHYSICAL_METADATA.coefficients(leaf)
                item=builder.modes.numeric_metadata(leaf,lambda p:unit)
                representation='POLYNOMIAL'
            except NotImplementedError:
                numerator,denominator=sp.together(leaf).as_numer_denom()
                denominator_unit=d.measure(denominator)
                numerator_unit=tuple(a+b for a,b in zip(unit,denominator_unit))
                parts={'EXACT_NUMERATOR':numerator,'EXACT_DENOMINATOR':denominator}
                item=builder.modes.numeric_metadata(engine.cas(parts),lambda p:
                    numerator_unit if p[0]=='EXACT_NUMERATOR' else denominator_unit)
                representation='EXACT_RATIONAL_NUMERATOR_DENOMINATOR'
            metadata.append({'OBJECT_PATH':path,'VALUE_DIMENSION_L_T_M':unit,
                             'GRADE_REPRESENTATION':representation,'GRADE_DATA':item})
        tag=prefix+'_'+name
        fingerprint=engine.carrier_fingerprint(body)
        literal=name.endswith('_RESIDUAL') and small_literal(body)
        engine.emit(tag,body if literal else fingerprint)
        engine.emit('METADATA_'+tag,metadata)
        records.append({'tag':tag,'body':body,'fingerprint':fingerprint,
                        'literal':literal,'metadata':engine.cas(metadata)})
    output('PROVENANCE',{**provenance,'emitterSha256':digest(Path(__file__)),
        'calculationSha256':summary['objectsSha256']})
    if diagnostic is not None:
        output('EXACT_ARITHMETIC_RUNTIME',diagnostic,
               {path:d.zero for path,_ in engine.leaves(engine.cas(diagnostic))})
    if workers is not None:
        output('ISOLATED_SCALAR_RUNTIME',workers,
               {path:(0,1,0) if path[-1] in ('wallSeconds','cpuSeconds') else d.zero
                for path,_ in engine.leaves(engine.cas(workers))})
    binding_values={str(key):value for key,value in packet['bindings'].items()}
    binding_units={(str(key),):d.measure(key) for key in packet['bindings']}
    output('MATERIAL_BINDINGS',binding_values,binding_units)
    for key,value in result.items():
        if key.endswith('_DOMAIN'):
            output('CONSTRUCTION_'+key,value,{():d.zero},value.lhs-value.rhs)
        else:output('CONSTRUCTION_'+key,value,builder.output_units(key,value))
    for key,value in checks.items():
        units=builder.check_output_units(key,value,checks)
        if key=='SLAB_CURRENT_DIAGONAL_JOIN_RESIDUAL':units=builder.output_units('SLAB_CURRENT_MATRIX',value)
        if key=='FINITE_COMPOSITION_DENOMINATOR_MATRIX':
            units={path:d.measure(leaf) for path,leaf in engine.leaves(engine.cas(value))}
        if key in ('EQUAL_DEPTH_WAVE_ROWS','EQUAL_DEPTH_WAVE_ELIMINANT'):
            unit=tuple(2*v for v in d.measure(builder.r.omega))
            units={(i,):unit for i in range(2)} if key.endswith('_ROWS') else {():unit}
        output('CHECK_'+key,value,units)
        if key in packet['retained']:
            output('RETAINED_'+key,packet['retained'][key],units)
            output('REMAINDER_'+key,packet['remainders'][key],units)
    output('RESIDUAL_CENSUS',summary['records'],
           {path:d.zero for path,_ in engine.leaves(engine.cas(summary['records']))})
    output('DIMENSION_CONSTRAINTS',tuple(d.constraints))
    output('RESOURCES',{'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
           {('wallSeconds',):(0,1,0),('peakRssKiB',):d.zero})
    # The source-line index is operational data, computed from actual emission.
    index=engine.emission_index(engine.EMISSION_LINES)
    output('EMISSION_LINES',index,{path:d.zero for path,_ in engine.leaves(engine.cas(index))})
    atomic(base/'emissions.pickle',pickle.dumps(records,protocol=5))
    atomic(base/'emission.json',(json.dumps({'end':args.end,'emitterSha256':digest(Path(__file__)),
        'calculationSha256':summary['objectsSha256'],'records':len(records),
        'emissionsSha256':digest(base/'emissions.pickle'),'dimensionConstraints':[sp.srepr(v) for v in d.constraints],
        'diagnosticRuntime':diagnostic,
        'workerRuntime':workers,
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},indent=2)+'\n').encode())
    if summary['retainedNonzeroScalars']:
        raise ValueError('retained endpoint balance residual; inspect emitted operands')


if __name__=='__main__':run()
