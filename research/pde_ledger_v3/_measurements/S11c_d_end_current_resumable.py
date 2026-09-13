#!/usr/bin/env python3
"""Source-pinned, resumable full symbolic end-current source construction."""
import argparse
from copy import copy
import faulthandler
import hashlib
import json
import os
from pathlib import Path
import pickle
import resource
import time
import threading

import sympy as sp
from sympy.polys.rings import ring
from S11c_d_joint_sheet_check import load, engine
from S11c_d_modal_current_check import digest, source_node
from S11c_d_current_runtime_metadata import cancellation_units, cancellation_packet


def atomic(path, body):
    temporary=path.with_name(path.name+'.partial')
    with temporary.open('wb') as stream:
        stream.write(body);stream.flush();os.fsync(stream.fileno())
    os.replace(temporary,path)


def expression_hash(value):return hashlib.sha256(sp.srepr(engine.cas(value)).encode()).hexdigest()


def polynomial_identity(left, right):
    """Combine exact fractions before sparse polynomial cross-multiplication.

    This avoids artificial denominator multiplicities in an uncombined sum.
    The caller also preserves the original, uncombined denominator domains.
    """
    a,b=sp.together(left).as_numer_denom();c,d=sp.together(right).as_numer_denom()
    symbols=tuple(sorted(a.free_symbols|b.free_symbols|c.free_symbols|d.free_symbols,key=sp.default_sort_key))
    polynomial_ring=ring(symbols,sp.QQ_I)[0]
    aa,bb,cc,dd=(polynomial_ring.from_expr(value) for value in (a,b,c,d))
    return (aa*dd-cc*bb).as_expr(),(a,b,c,d)


def main():
    parser=argparse.ArgumentParser()
    for name in ('manifest','input','run-directory'):parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--end',choices=('LEFT','RIGHT'),required=True)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--seed-run',type=Path)
    args=parser.parse_args();destination=args.run_directory
    destination.mkdir(parents=True,exist_ok=True);state=destination/'state';state.mkdir(exist_ok=True)
    started=time.monotonic();stage='loading';lock=threading.Lock();finished=threading.Event()
    def progress(name,**values):
        nonlocal stage
        if name!='resource':stage=name
        info={'stage':name,'activeStage':stage,'elapsedSeconds':time.monotonic()-started,
              'cpuSeconds':time.process_time(),'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,**values}
        with lock:
            with (destination/'progress.jsonl').open('a') as stream:stream.write(json.dumps(info)+'\n')
    def monitor():
        while not finished.wait(45):
            usage={}
            try:
                usage={line.split(':')[0]:line.split(':')[1].strip() for line in Path('/proc/self/status').read_text().splitlines()
                    if line.startswith(('VmRSS:','VmHWM:','VmSize:','Threads:'))}
            except OSError:pass
            progress('resource',processStatus=usage)
    faulthandler.enable()
    trace=(destination/'stack-samples.txt').open('a')
    faulthandler.dump_traceback_later(300,repeat=True,file=trace)
    threading.Thread(target=monitor,daemon=True).start()
    modes,source_strong,units,branches,inputs,provenance=load(args)
    manifest=json.loads(args.manifest.read_text());base=Path(manifest['run_directory'])
    root=Path(engine.__file__).resolve().parents[1]
    instrument=str(Path(__file__).resolve().relative_to(root))
    sources=('scripts/S11c_d_mixing_scattering_sympy_audit.py',instrument,
             '_measurements/S11c_d_joint_sheet_check.py','_measurements/S11c_d_modal_current_check.py',
             '_measurements/S11c_d_current_runtime_metadata.py',
             'scripts/S11c_d_output_codec.py','scripts/ledger_fold.py')
    pins={name:digest(root/name) for name in sources}
    signature={'sources':pins,'manifest':digest(args.manifest),'input':digest(args.input),'end':args.end}
    signature_path=state/'signature.json'
    if signature_path.exists():
        if not args.resume or json.loads(signature_path.read_text())!=signature:raise ValueError('resume source/input signature mismatch')
    elif args.resume:raise ValueError('no resumable state')
    else:atomic(signature_path,(json.dumps(signature,indent=2)+'\n').encode())
    for name in sources:
        target=destination/'source'/name;target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists() and digest(target)!=pins[name]:raise ValueError('source snapshot differs')
        if not target.exists():target.write_bytes((root/name).read_bytes())
    seed_records={}
    if args.seed_run:
        seed_state=args.seed_run/'state'
        seed_signature=json.loads((seed_state/'signature.json').read_text())
        if any(seed_signature[key]!=signature[key] for key in ('manifest','input','end')):
            raise ValueError('seed input mismatch')
        for name,pin in seed_signature['sources'].items():
            if digest(args.seed_run/'source'/name)!=pin:raise ValueError(('seed frozen source',name))
            # Only instrumentation and emission metadata may change. The complete
            # native construction, loader, codec and ledger must remain identical.
            if name not in (instrument,'_measurements/S11c_d_current_runtime_metadata.py') and pins.get(name)!=pin:
                raise ValueError(('seed construction source',name))
        seed_proofs=source_node((args.seed_run/'source'/instrument).read_text(),'polynomial_identity')==source_node(Path(__file__).read_text(),'polynomial_identity')
        for path in seed_state.glob('*.pickle'):
            name=path.stem
            if name.endswith('.input') or not (name in ('conservative','balance','acoustic') or name.startswith('rational_') or (seed_proofs and name.startswith('proof_'))):continue
            pin=path.with_suffix('.sha256')
            if not pin.exists() or digest(path)!=pin.read_text().strip():raise ValueError(('seed packet digest',name))
            packet=pickle.loads(path.read_bytes())
            if expression_hash(packet['source'])!=packet['sourceSha256']:raise ValueError(('seed packet operand',name))
            seed_records[name]=(packet,{'path':str(path),'sha256':digest(path)})
    frozen=base/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py'
    if source_node(frozen.read_text(),'UniformSlabCurrent')!=source_node(Path(engine.__file__).read_text(),'UniformSlabCurrent'):
        raise ValueError('end energy constructor differs from native producer')
    energy=pickle.loads((base/'symbols'/f'{args.end}_LAB_HELD_RHO4_CONSTANT.pickle').read_bytes())[5]
    r=copy(modes.r);r.x=tuple(r.symbols['s11cc2X'+str(i)] for i in (1,2,3));r.t=r.symbols['s11cc2Time'];r.z=r.symbols['s11cdNormalPosition']
    r.end_values={key:inputs.limits[value] for key,value in r.end_values.items()}
    endpoint=-sp.oo if args.end=='LEFT' else sp.oo
    strong=source_strong.xreplace(inputs.limits)
    ends=engine.ConstantEndPencil.__new__(engine.ConstantEndPencil);ends.r,ends.kn=r,modes.k
    modes=engine.FullPencilModes(ends,modes.curl,modes.units)
    current=engine.UniformSlabCurrent(r,{'value':energy},ends,strong[3,:])
    balance=engine.SlabEnergyBalance(current);acoustic=engine.ClosedAcousticEnergy(balance,modes,strong)
    provenance.update({'instrumentPath':instrument,'instrumentSha256':digest(Path(__file__)),
        'profileLimits':{str(k):str(v) for k,v in inputs.limits.items()},'endpoint':str(endpoint),
        'sourceEnergySha256':expression_hash(energy),'sourceFiles':pins,
        'runtime':{'sympy':sp.__version__,'materialParametersBound':False,'backgroundGradesTruncatedByCancellation':False}})
    if args.seed_run:
        provenance['seedRun']={'path':str(args.seed_run),'signatureSha256':digest(args.seed_run/'state/signature.json'),
                              'records':{name:pin for name,(_,pin) in seed_records.items()}}
    def saved(name,source,calculate):
        path=state/(name+'.pickle');pin=state/(name+'.sha256')
        source_hash=expression_hash(source)
        if path.exists():
            if not args.resume or not pin.exists() or digest(path)!=pin.read_text().strip():raise ValueError(('checkpoint digest',name))
            packet=pickle.loads(path.read_bytes())
            if packet['sourceSha256']!=source_hash:raise ValueError(('checkpoint operand differs',name))
            engine.PHYSICAL_METADATA.dimensions.known.update(packet['dimensions'])
            progress(name+'_resumed');return packet['value']
        atomic(state/(name+'.input.pickle'),pickle.dumps(source,protocol=5))
        if name in seed_records:
            seed,pin_record=seed_records[name]
            if seed['sourceSha256']!=source_hash:raise ValueError(('seed requested operand',name))
            engine.PHYSICAL_METADATA.dimensions.known.update(seed['dimensions'])
            value=seed['value'];progress(name+'_seeded',seed=pin_record)
        else:
            progress(name+'_started');value=calculate()
        packet={'sourceSha256':source_hash,'source':source,'value':value,'dimensions':engine.PHYSICAL_METADATA.dimensions.known}
        atomic(path,pickle.dumps(packet,protocol=5));atomic(pin,(digest(path)+'\n').encode())
        progress(name+'_saved',sha256=digest(path));return value
    coefficient_sources={}
    def rational_checkpoint(name,source,calculate):
        value=saved('rational_'+name,source,calculate);coefficient_sources[name]=(source,value);return value
    acoustic.rational_checkpoint=rational_checkpoint
    source=saved('conservative',strong,lambda:current.construct('LAB_HELD',endpoint))
    current.construct=lambda anchoring,end:source
    slab=saved('balance',strong,lambda:balance.construct('LAB_HELD',endpoint))
    balance.construct=lambda anchoring,end:slab
    bulk=saved('acoustic',strong,lambda:acoustic.construct('LAB_HELD',endpoint))
    acoustic.construct=lambda anchoring,end:bulk
    # An unchanged, seeded acoustic construction includes its complete arithmetic
    # census, even though its constructor did not need to execute again.
    if 'acoustic' in seed_records:
        for name,(packet,_) in seed_records.items():
            if name.startswith('rational_'):
                saved(name,packet['source'],lambda:None)
    # Every saved rational operation retains its actual source, including
    # the full reconstruction. Restore this proof census after a resumed stage.
    for path in sorted(state.glob('rational_*.pickle')):
        if path.name.endswith('.input.pickle'):continue
        pin=path.with_suffix('.sha256')
        if not pin.exists() or digest(path)!=pin.read_text().strip():raise ValueError('rational state digest')
        packet=pickle.loads(path.read_bytes())
        if expression_hash(packet['source'])!=packet['sourceSha256']:raise ValueError('rational source identity')
        coefficient_sources[path.stem.removeprefix('rational_')]=(packet['source'],packet['value'])
    progress('source_constructed')
    proof={}
    for name,(before,after) in coefficient_sources.items():
        residual,fractions=saved('proof_'+name,(before,after),lambda before=before,after=after:polynomial_identity(before,after))
        a,b,c,d=fractions
        source_denominator=before.as_numer_denom()[1]
        result_denominator=after.as_numer_denom()[1]
        proof[name]={'BEFORE':before,'AFTER':after,'SOURCE_NUMERATOR':a,'SOURCE_DENOMINATOR':b,
            'RESULT_NUMERATOR':c,'RESULT_DENOMINATOR':d,'CROSS_PRODUCT_RESIDUAL':residual,
            'SOURCE_DENOMINATOR_DOMAIN':sp.Ne(b,0),'RESULT_DENOMINATOR_DOMAIN':sp.Ne(d,0),
            'ORIGINAL_SOURCE_DENOMINATOR':source_denominator,'ORIGINAL_RESULT_DENOMINATOR':result_denominator,
            'ORIGINAL_SOURCE_DENOMINATOR_DOMAIN':sp.Ne(source_denominator,0),
            'ORIGINAL_RESULT_DENOMINATOR_DOMAIN':sp.Ne(result_denominator,0)}
    payload={'conservative':source,'slab':slab,'acoustic':bulk,'strong':strong,
        'profileBindings':inputs.limits,'knownDimensions':engine.PHYSICAL_METADATA.dimensions.known,
        'cancellationProof':proof}
    objects=destination/'objects.pickle';atomic(objects,pickle.dumps(payload,protocol=5))
    prefix='END_CURRENT_SOURCE_'+args.end+'_LAB_HELD_RHO4_CONSTANT'
    def output(name,value,unit=lambda path:(0,0,0),heavy=False,metadata_body=None,metadata=None):
        body=engine.cas(value);engine.emit(prefix+'_'+name,engine.carrier_fingerprint(body) if heavy else body)
        engine.emit('METADATA_'+prefix+'_'+name,metadata if metadata is not None else modes.numeric_metadata(body if metadata_body is None else metadata_body,unit))
    output('PROVENANCE',provenance);output('SOURCE_BOUND_PENCIL',strong,lambda path:units[path],True)
    current.emit('LAB_HELD',endpoint,prefix);balance.emit('LAB_HELD',endpoint,prefix);acoustic.emit('LAB_HELD',endpoint,prefix)
    proof_packets=[];dims=engine.PHYSICAL_METADATA.dimensions
    proof_units=cancellation_units(engine,proof,strong,dims)
    for name,record in proof.items():
        for key,value in record.items():
            item=cancellation_packet(engine,modes,name,key,record,proof_units[name][key])
            output(item['name'],item['body'],heavy=item['heavy'],metadata=item['metadata'])
            proof_packets.append(item)
    atomic(destination/'proof-emissions.pickle',pickle.dumps(proof_packets,protocol=5))
    inventory={}
    for group,values in (('conservative',source),('slab',slab),('acoustic',bulk)):
        for key,value in values.items():
            if not key.endswith('_RESIDUAL'):continue
            entries=[v for _,v in engine.leaves(engine.cas(value)) if not isinstance(v,engine.Str)]
            projected=[current.retained(v) for v in entries]
            inventory[group+'_'+key]={'scalars':len(entries),'rawNonzeroScalars':sum(v!=0 for v in entries),
                'retainedNonzeroScalars':sum(v!=0 for v in projected),
                'retainedGrades':sorted({g for v in projected for g in engine.PHYSICAL_METADATA.coefficients(v)})}
    output('RESIDUAL_GRADE_CENSUS',inventory);output('DIMENSION_CONSTRAINTS',tuple(dims.constraints))
    index=engine.emission_index(engine.EMISSION_LINES);engine.emit(prefix+'_EMISSION_LINES',index)
    engine.emit('METADATA_'+prefix+'_EMISSION_LINES',modes.numeric_metadata(engine.cas(index),lambda path:(0,0,0)))
    summary={'provenance':provenance,'residualInventory':inventory,'objectsSha256':digest(objects),
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'retainedNonzeroScalars':sum(v['retainedNonzeroScalars'] for v in inventory.values()),
        'cancellationIdentityCount':len(proof),'nonzeroCancellationIdentities':sum(v['CROSS_PRODUCT_RESIDUAL']!=0 for v in proof.values()),
        'proofEmissionsSha256':digest(destination/'proof-emissions.pickle')}
    atomic(destination/'checks.json',(json.dumps(summary,indent=2)+'\n').encode());progress('completed',summary=summary)
    finished.set();faulthandler.cancel_dump_traceback_later();trace.close()
    if pins!={name:digest(root/name) for name in sources}:raise ValueError('source changed during run')
    if summary['retainedNonzeroScalars'] or summary['nonzeroCancellationIdentities']:
        raise ValueError('computed source residual; inspect emitted operands')


if __name__=='__main__':main()
