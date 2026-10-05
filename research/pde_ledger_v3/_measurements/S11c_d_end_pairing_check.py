#!/usr/bin/env python3
"""Source-checked endpoint two-frequency pairing and retained-balance packets."""
import argparse
import ast
from copy import copy
import faulthandler
import hashlib
import json
import os
from pathlib import Path
import pickle
import resource
import threading
import time

import sympy as sp
from S11c_d_joint_sheet_check import load, engine
from S11c_d_modal_current_check import build as reference_build, digest, source_node

ROOT = Path(__file__).resolve().parents[1]


def atomic(path, value):
    temporary=path.with_name(path.name+'.partial')
    with temporary.open('wb') as stream:
        stream.write(value);stream.flush();os.fsync(stream.fileno())
    os.replace(temporary,path)


def build(args):
    if args.end=='REFERENCE':
        return reference_build(args)
    if args.source_checkpoint is None:
        raise ValueError('endpoint source checkpoint required')
    modes,strong,units,branches,inputs,provenance=load(args)
    manifest=json.loads(args.manifest.read_text());base=Path(manifest['run_directory'])
    energy=pickle.loads((base/'symbols'/f'{args.end}_LAB_HELD_RHO4_CONSTANT.pickle').read_bytes())[5]
    checkpoint=json.loads(args.source_checkpoint.read_text());source=Path(checkpoint['runDirectory'])
    if checkpoint['retainedNonzeroScalars'] or checkpoint['nonzeroCancellationIdentities']:
        raise ValueError('source checkpoint retains a discrepancy')
    old=checkpoint['provenance']
    if old['end']!=args.end or old['case']!=provenance['case'] or old['inputSha256']!=digest(args.input):
        raise ValueError('source checkpoint case/end/input mismatch')
    for name,record in checkpoint['artifacts'].items():
        if digest(source/name)!=record['sha256']:
            raise ValueError(('source artifact mismatch',name))
    for name,sha in old['sourceFiles'].items():
        if digest(source/'source'/name)!=sha:raise ValueError(('source snapshot mismatch',name))
        if name!='scripts/S11c_d_mixing_scattering_sympy_audit.py' and digest(ROOT/name)!=sha:
            raise ValueError(('source instrument changed',name))
    frozen=(source/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py').read_text()
    now=Path(engine.__file__).read_text()
    constructors=('UniformSlabCurrent','SlabEnergyBalance','ClosedAcousticEnergy','ConstantEndPencil','polynomial_terms')
    for name in constructors:
        if source_node(frozen,name)!=source_node(now,name):raise ValueError(('source constructor changed',name))
    for name,sha in old['producerSources'].items():
        if name.startswith('directives/') or name.endswith('_exports.py') or name=='scripts/ledger_fold.py':
            if digest(ROOT/name)!=sha:raise ValueError(('source physics input changed',name))
    packet=pickle.loads((source/'objects.pickle').read_bytes())
    strong=strong.xreplace(inputs.limits)
    if packet['strong']!=strong or packet['profileBindings']!=inputs.limits:
        raise ValueError('endpoint source/full pencil or profile join differs')
    energy_hash=hashlib.sha256(sp.srepr(engine.cas(energy)).encode()).hexdigest()
    if energy_hash!=old['sourceEnergySha256']:raise ValueError('endpoint energy source differs')
    r=copy(modes.r)
    r.x=tuple(r.symbols['s11cc2X'+str(i)] for i in (1,2,3));r.t=r.symbols['s11cc2Time'];r.z=r.symbols['s11cdNormalPosition']
    r.end_values={key:inputs.limits[value] for key,value in r.end_values.items()}
    endpoint=-sp.oo if args.end=='LEFT' else sp.oo
    ends=engine.ConstantEndPencil.__new__(engine.ConstantEndPencil);ends.r,ends.kn=r,modes.k
    modes=engine.FullPencilModes(ends,modes.curl,modes.units)
    current=engine.UniformSlabCurrent(r,{'value':energy},ends,strong[3,:])
    balance=engine.SlabEnergyBalance(current);acoustic=engine.ClosedAcousticEnergy(balance,modes,strong)
    engine.PHYSICAL_METADATA.dimensions.known.update(packet['knownDimensions'])
    def cached(value):
        def fetch(anchoring,end):
            if (anchoring,end)!=('LAB_HELD',endpoint):raise ValueError('cached endpoint source context mismatch')
            return value
        return fetch
    current.construct=cached(packet['conservative']);balance.construct=cached(packet['slab']);acoustic.construct=cached(packet['acoustic'])
    pairing=engine.ClosedCurrentPairing(acoustic,'LAB_HELD',endpoint)
    provenance.update({'instrumentSha256':digest(Path(__file__)),
        'sourceCheckpointSha256':digest(args.source_checkpoint),'sourcePacketSha256':digest(source/'objects.pickle'),
        'sourceConstructors':constructors,'sourceEnergySha256':energy_hash,
        'sourcePencilStructuralJoin':packet['strong']==strong,'sourceProfileStructuralJoin':packet['profileBindings']==inputs.limits})
    return pairing,inputs,provenance


def run():
    parser=argparse.ArgumentParser()
    for name in ('manifest','input','run-directory'):parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--end',choices=('REFERENCE','LEFT','RIGHT'),required=True)
    parser.add_argument('--source-checkpoint',type=Path)
    parser.add_argument('--current-manifest',type=Path)
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args();destination=args.run_directory
    destination.mkdir(parents=True,exist_ok=True)
    started=time.monotonic();finished=threading.Event();lock=threading.Lock();active='loading'
    def progress(record):
        nonlocal active
        if record['stage']!='resource':active=record['stage']
        record={**record,'activeStage':active,'elapsedSeconds':time.monotonic()-started,
                'cpuSeconds':time.process_time(),'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
        with lock:
            with (destination/'progress.jsonl').open('a') as stream:stream.write(json.dumps(record)+'\n')
    def monitor():
        while not finished.wait(45):progress({'stage':'resource'})
    threading.Thread(target=monitor,daemon=True).start()
    faulthandler.enable();trace=(destination/'stack-samples.txt').open('a')
    faulthandler.dump_traceback_later(300,repeat=True,file=trace)
    builder,inputs,provenance=build(args)
    files=('scripts/S11c_d_mixing_scattering_sympy_audit.py','_measurements/S11c_d_end_pairing_check.py',
           '_measurements/S11c_d_modal_current_check.py','_measurements/S11c_d_joint_sheet_check.py',
           'scripts/S11c_d_output_codec.py','scripts/ledger_fold.py')
    pins={name:digest(ROOT/name) for name in files}
    signature=json.loads(json.dumps({'sourceFiles':pins,'provenance':provenance,'end':args.end}))
    signature_path=destination/'signature.json'
    if signature_path.exists():
        if not args.resume or json.loads(signature_path.read_text())!=signature:raise ValueError('resume signature differs')
    else:
        if args.resume:raise ValueError('missing resume signature')
        atomic(signature_path,(json.dumps(signature,indent=2)+'\n').encode())
        for name in files:
            target=destination/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes((ROOT/name).read_bytes())
    def saved(name,calculate):
        path=destination/(name+'.pickle');pin=destination/(name+'.sha256')
        if path.exists():
            if not args.resume or digest(path)!=pin.read_text().strip():raise ValueError(('cache pin',name))
            value,known=pickle.loads(path.read_bytes());engine.PHYSICAL_METADATA.dimensions.known.update(known)
            progress({'stage':name+'_resumed'});return value
        progress({'stage':name+'_started'});value=calculate()
        atomic(path,pickle.dumps((value,engine.PHYSICAL_METADATA.dimensions.known),protocol=5))
        atomic(pin,(digest(path)+'\n').encode());progress({'stage':name+'_saved','sha256':digest(path)});return value
    progress({'stage':'loaded','provenance':provenance})
    result=saved('construction',lambda:builder.construct(builder.anchoring,builder.end))
    modes=builder.modes
    algebraic,relation,_=modes.analytic(builder.acoustic.strong)
    mapping=inputs.mapping(algebraic,relation,(modes.k,modes.q,modes.eta,modes.sigma,modes.r.omega))
    if args.end=='REFERENCE':mapping.update({modes.eta:0,modes.sigma:0})
    checks=saved('balance',lambda:builder.split_balance_checks(result,mapping,progress))
    derivatives=saved('derivatives',lambda:builder.split_derivative_checks(result,checks,mapping,progress))
    checks={**checks,**derivatives}
    same_frequency=dict.fromkeys(builder.frequencies,builder.r.omega)
    checks['SLAB_CURRENT_DIAGONAL_JOIN_RESIDUAL']=(result['SLAB_CURRENT_MATRIX'].xreplace(same_frequency)-
        builder.balance.construct(builder.anchoring,builder.end)['SLAB_CURRENT_MATRIX']).applyfunc(sp.expand)
    def project(value):
        if isinstance(value,sp.MatrixBase):return value.applyfunc(builder.c.retained)
        if isinstance(value,(tuple,list)):return tuple(project(v) for v in value)
        return builder.c.retained(value)
    projected=saved('retained',lambda:{key:project(value) for key,value in checks.items() if key.endswith('_RESIDUAL')})
    # Tuple containers require componentwise subtraction rather than tuple algebra.
    def subtract(a,b):
        if isinstance(a,(tuple,list)):return tuple(subtract(x,y) for x,y in zip(a,b))
        return a-b
    remainders={key:subtract(checks[key],value) for key,value in projected.items()}
    summary={'provenance':provenance,'sourceFiles':pins,'end':args.end,'records':{},
             'materialBindings':{str(k):str(v) for k,v in mapping.items()}}
    for key,value in projected.items():
        raw=[v for _,v in engine.leaves(engine.cas(checks[key]))]
        retained=[v for _,v in engine.leaves(engine.cas(value))]
        summary['records'][key]={'scalars':len(raw),'rawNonzeroScalars':sum(v!=0 for v in raw),
            'retainedNonzeroScalars':sum(v!=0 for v in retained)}
    summary.update({'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'retainedNonzeroScalars':sum(v['retainedNonzeroScalars'] for v in summary['records'].values()),
        'sourcePinsStable':pins=={name:digest(ROOT/name) for name in pins}})
    packet={'result':result,'checks':checks,'retained':projected,'remainders':remainders,'bindings':mapping,'summary':summary}
    atomic(destination/'complete.pickle',pickle.dumps((packet,engine.PHYSICAL_METADATA.dimensions.known),protocol=5))
    summary['objectsSha256']=digest(destination/'complete.pickle')
    atomic(destination/'checks.json',(json.dumps(summary,indent=2)+'\n').encode())
    progress({'stage':'complete','summary':summary})
    faulthandler.cancel_dump_traceback_later();finished.set()
    print(json.dumps({'end':args.end,'records':summary['records'],'wallSeconds':summary['wallSeconds'],
        'retainedNonzeroScalars':summary['retainedNonzeroScalars']},indent=2),flush=True)


if __name__=='__main__':run()
