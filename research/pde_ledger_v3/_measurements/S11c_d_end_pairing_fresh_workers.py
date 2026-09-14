#!/usr/bin/env python3
"""Fresh endpoint construction with compatible exact-scalar cache reuse.

Only pure native scalar operations cross the previous checkpoint boundary.
The endpoint construction, material binding and all physical checks run again.
"""
import argparse
import json
from pathlib import Path
import pickle
import resource
import shutil
import sys
import time

import S11c_d_end_pairing_workers as scalar
checker=scalar.checker
HERE=Path(__file__).resolve()


def checked_cache(checkpoint_path, environment):
    checkpoint=json.loads(checkpoint_path.read_text())
    runtime=checkpoint['workerRuntime']
    if runtime['signature']['environment']!=environment:raise ValueError('arithmetic cache runtime/source mismatch')
    base=Path(checkpoint['runDirectory']);state=base/'worker-state'
    if checker.digest(state/'signature.json')!=runtime['signatureSha256']:
        raise ValueError('arithmetic cache signature pin')
    if checker.digest(state/'operations.jsonl')!=runtime['operationLogSha256']:
        raise ValueError('arithmetic cache operation log pin')
    reader=scalar.ScalarWorkers.__new__(scalar.ScalarWorkers)
    reader.signature=runtime['signature'];values={}
    for name,expected in runtime['workers'].items():
        directory=state/name
        expression,value,record=reader.read_result(directory)
        if record!=expected:raise ValueError('arithmetic cache worker record differs')
        key=(record['method'],expression)
        if key in values and values[key][0]!=value:raise ValueError('conflicting arithmetic cache')
        values[key]=(value,{'kind':'external_worker','directory':str(directory),
            'checkpoint':str(checkpoint_path),'checkpointSha256':checker.digest(checkpoint_path),
            'resultSha256':record['resultSha256']},{})
    for name,expected in runtime['reusedDiagnosticPackets'].items():
        source=Path(expected['input']);packet=Path(name)
        if checker.digest(source)!=expected['inputSha256'] or checker.digest(packet)!=expected['sha256']:
            raise ValueError('arithmetic diagnostic cache pin')
        before,value,known=pickle.loads(packet.read_bytes())
        if before!=pickle.loads(source.read_bytes()):raise ValueError('arithmetic diagnostic operand join')
        method=next((method for method in scalar.METHODS if packet.name.startswith(method+'_')),None)
        if method is None:raise ValueError('arithmetic diagnostic method')
        key=(method,before)
        if key in values and values[key][0]!=value:raise ValueError('conflicting diagnostic arithmetic cache')
        values[key]=(value,{'kind':'external_diagnostic','packet':str(packet),'input':str(source),
            'checkpoint':str(checkpoint_path),'checkpointSha256':checker.digest(checkpoint_path),
            'sha256':expected['sha256'],'inputSha256':expected['inputSha256']},known)
    for item in runtime['boundary']:
        path=Path(item['comparisonPacket'])
        if checker.digest(path)!=item['comparisonSha256']:raise ValueError('boundary evidence pin')
        source,before,after,residual=pickle.loads(path.read_bytes())
        if before!=after or residual!=after-before:raise ValueError('boundary arithmetic comparison')
    if not runtime['boundary']:raise ValueError('no completed boundary comparisons')
    return values,runtime['boundary']


class FreshScalarWorkers(scalar.ScalarWorkers):
    def __init__(self,base,checkpoint):
        self.base=base;self.state=base/'worker-state';self.state.mkdir(parents=True,exist_ok=True)
        self.started=time.monotonic()
        self.signature={'mode':'FRESH_CONSTRUCTION_EXACT_SCALAR_CACHE','environment':scalar.environment(),
            'checkerSignatureSha256':checker.digest(base/'signature.json'),
            'orchestratorSha256':checker.digest(HERE),'arithmeticCheckpoint':str(checkpoint),
            'arithmeticCheckpointSha256':checker.digest(checkpoint)}
        pin=self.state/'signature.json'
        if pin.exists():
            if json.loads(pin.read_text())!=self.signature:raise ValueError('fresh-worker resume signature')
        else:
            scalar.write_json(pin,self.signature)
            shutil.copyfile(scalar.HERE,self.state/scalar.HERE.name)
            shutil.copyfile(HERE,self.state/HERE.name)
        self.cache,boundary=checked_cache(checkpoint,self.signature['environment'])
        scalar.write_json(self.state/'boundary.json',boundary)
        for path in sorted(self.state.glob('call_*/result.json')):
            expression,value,record=self.read_result(path.parent)
            key=(record['method'],expression)
            if key in self.cache and self.cache[key][0]!=value:raise ValueError('conflicting resumed scalar')
            self.cache[key]=(value,{'kind':'isolated_worker','directory':str(path.parent),
                                  'resultSha256':record['resultSha256']},{})
        self.progress({'stage':'arithmetic_cache_loaded','distinctExactOperands':len(self.cache),
                       'completedBoundaryComparisonsReused':len(boundary)})


def runtime_record(base):
    state=base/'worker-state';signature=json.loads((state/'signature.json').read_text())
    if checker.digest(HERE)!=signature['orchestratorSha256'] or checker.digest(state/HERE.name)!=checker.digest(HERE):
        raise ValueError('fresh worker orchestrator pin')
    if checker.digest(scalar.HERE)!=signature['environment']['workerSha256'] or checker.digest(state/scalar.HERE.name)!=checker.digest(scalar.HERE):
        raise ValueError('scalar worker source pin')
    if checker.digest(Path(checker.engine.__file__))!=signature['environment']['nativeSha256']:
        raise ValueError('fresh worker native pin')
    if checker.digest(base/'signature.json')!=signature['checkerSignatureSha256']:
        raise ValueError('fresh checker signature pin')
    checkpoint=Path(signature['arithmeticCheckpoint'])
    if checker.digest(checkpoint)!=signature['arithmeticCheckpointSha256']:raise ValueError('arithmetic checkpoint pin')
    cache,boundary=checked_cache(checkpoint,signature['environment'])
    if json.loads((state/'boundary.json').read_text())!=boundary:raise ValueError('boundary copy join')
    operations=[json.loads(v) for v in (state/'operations.jsonl').read_text().splitlines()]
    if not operations or operations[-1]['stage']!='complete' or not operations[-1]['environmentStable']:
        raise ValueError('fresh worker calculation incomplete')
    reader=scalar.ScalarWorkers.__new__(scalar.ScalarWorkers);reader.signature=signature
    workers={};reused={}
    cache_provenance={json.dumps(v[1],sort_keys=True) for v in cache.values()}
    for operation in operations:
        if operation['stage']=='worker_saved':
            directory=Path(operation['directory']);_,_,record=reader.read_result(directory)
            if record!=operation['resource']:raise ValueError('fresh worker resource join')
            workers[directory.name]=record
        elif operation['stage']=='reused':
            provenance=operation['provenance'];kind=provenance['kind']
            if kind=='isolated_worker':
                _,_,record=reader.read_result(Path(provenance['directory']))
                if record['resultSha256']!=provenance['resultSha256']:raise ValueError('resumed worker result pin')
            elif kind in ('external_worker','external_diagnostic'):
                key=json.dumps(provenance,sort_keys=True)
                if key not in cache_provenance:raise ValueError('external arithmetic provenance join')
                reused[provenance.get('directory',provenance.get('packet'))]=provenance
            else:raise ValueError('unknown fresh scalar provenance')
    return {'signature':signature,'signatureSha256':checker.digest(state/'signature.json'),
        'operationLogSha256':checker.digest(state/'operations.jsonl'),'operationCount':len(operations),
        'boundary':boundary,'workers':workers,'reusedExactScalarOperands':reused}


def main():
    parser=argparse.ArgumentParser(add_help=False)
    parser.add_argument('--arithmetic-checkpoint',type=Path,required=True)
    parser.add_argument('--inspect-arithmetic-cache',type=Path,
                        help='Write a read-only cache/provenance inventory and exit.')
    args,remaining=parser.parse_known_args()
    args.arithmetic_checkpoint=args.arithmetic_checkpoint.resolve()
    previous=json.loads(args.arithmetic_checkpoint.read_text())
    limits=previous['workerRuntime']['signature']['environment']['stackLimitBytes']
    resource.setrlimit(resource.RLIMIT_STACK,tuple(limits))
    if args.inspect_arithmetic_cache:
        started=time.monotonic()
        cache,boundary=checked_cache(args.arithmetic_checkpoint,scalar.environment())
        result={'sourceSha256':checker.digest(HERE),'checkpoint':str(args.arithmetic_checkpoint),
            'checkpointSha256':checker.digest(args.arithmetic_checkpoint),'distinctExactOperands':len(cache),
            'methods':sorted({key[0] for key in cache}),'completedBoundaryComparisons':len(boundary),
            'environment':scalar.environment(),'wallSeconds':time.monotonic()-started,
            'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'scope':'Read-only exact-operand cache and completed boundary validation; no fresh physical construction.'}
        scalar.write_json(args.inspect_arithmetic_cache,result)
        print(json.dumps(result,indent=2),flush=True)
        return
    base=Path(remaining[remaining.index('--run-directory')+1]);base.mkdir(parents=True,exist_ok=True)
    sys.excepthook=sys.__excepthook__
    coordinator=None
    def prepare():
        nonlocal coordinator
        coordinator=FreshScalarWorkers(base,args.arithmetic_checkpoint)
        for method in scalar.METHODS:
            def calculate(expression,method=method):return coordinator.calculate(method,expression)
            setattr(checker.engine.ClosedCurrentPairing,method,staticmethod(calculate))
    original_atomic=checker.atomic
    def atomic(path,body):
        original_atomic(path,body)
        if path==base/'signature.json' and coordinator is None:prepare()
    checker.atomic=atomic
    if (base/'signature.json').exists():prepare()
    sys.argv=[sys.argv[0],*remaining]
    checker.run()
    if coordinator is None:raise ValueError('worker setup was not reached')
    coordinator.progress({'stage':'complete','environmentStable':scalar.environment()==coordinator.signature['environment']})


if __name__=='__main__':main()
