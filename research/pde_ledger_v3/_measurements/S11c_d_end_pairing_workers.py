#!/usr/bin/env python3
"""Sequential isolated execution of unchanged native pairing scalar methods."""
import argparse
import faulthandler
import hashlib
import json
import os
from pathlib import Path
import pickle
import resource
import shutil
import subprocess
import sys
import time

import S11c_d_end_pairing_check as checker

METHODS=('rational_coefficient','carrier_expansion')
HERE=Path(__file__).resolve()


def write_json(path,value):
    checker.atomic(path,(json.dumps(value,indent=2)+'\n').encode())


def environment():
    return {'python':sys.version,'executable':str(Path(sys.executable).resolve()),
            'executableSha256':checker.digest(Path(sys.executable).resolve()),
            'sympy':checker.sp.__version__,'sympyInitSha256':checker.digest(Path(checker.sp.__file__)),
            'nativeSha256':checker.digest(Path(checker.engine.__file__)),
            'workerSha256':checker.digest(HERE),
            'stackLimitBytes':list(resource.getrlimit(resource.RLIMIT_STACK))}


def worker(request_path):
    sys.excepthook=sys.__excepthook__
    request=json.loads(request_path.read_text());directory=request_path.parent
    if request['environment']!=environment():raise ValueError('worker runtime/source differs')
    method=request['method']
    if method not in METHODS:raise ValueError('unknown scalar method')
    source=directory/'input.pickle'
    if checker.digest(source)!=request['sourceSha256']:raise ValueError('worker input digest')
    expression=pickle.loads(source.read_bytes())
    started=time.monotonic();faulthandler.enable()
    with (directory/'stack.txt').open('w') as trace:
        faulthandler.dump_traceback_later(300,repeat=True,file=trace)
        value=getattr(checker.engine.ClosedCurrentPairing,method)(expression)
        checker.atomic(directory/'result.pickle',pickle.dumps((expression,value),protocol=5))
        record={'requestSha256':checker.digest(request_path),'sourceSha256':request['sourceSha256'],
                'resultSha256':checker.digest(directory/'result.pickle'),'method':method,
                'wallSeconds':time.monotonic()-started,'cpuSeconds':time.process_time(),
                'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                'environmentStable':request['environment']==environment()}
        write_json(directory/'result.json',record)
        faulthandler.cancel_dump_traceback_later()
    print(json.dumps(record),flush=True)


class ScalarWorkers:
    def __init__(self,base,seed):
        self.base,self.seed=base,seed
        self.state=base/'worker-state';self.state.mkdir(parents=True,exist_ok=True)
        self.started=time.monotonic();self.cache={};self.original=[]
        self.signature={'environment':environment(),
                        'checkerSignatureSha256':checker.digest(base/'signature.json'),
                        'seedRun':str(seed),'seedSignatureSha256':checker.digest(seed/'signature.json')}
        pin=self.state/'signature.json'
        if pin.exists():
            if json.loads(pin.read_text())!=self.signature:raise ValueError('worker resume signature differs')
        else:
            write_json(pin,self.signature)
            checker.atomic(self.state/HERE.name,HERE.read_bytes())
        signature=json.loads((seed/'signature.json').read_text())
        if signature!=json.loads((base/'signature.json').read_text()):raise ValueError('seed/checker signature differs')
        for name,sha in signature['sourceFiles'].items():
            if checker.digest(checker.ROOT/name)!=sha or checker.digest(seed/'source'/name)!=sha:
                raise ValueError(('scalar seed source',name))
        diagnostic=json.loads((seed/'scalar-state/signature.json').read_text())
        instrument=Path(diagnostic['instrument'])
        if checker.digest(instrument)!=diagnostic['instrumentSha256'] or checker.digest(seed/'scalar-state'/instrument.name)!=diagnostic['instrumentSha256']:
            raise ValueError('scalar seed diagnostic instrument')
        if diagnostic['checkerSignatureSha256']!=checker.digest(seed/'signature.json'):
            raise ValueError('scalar seed diagnostic source signature')
        operations=[json.loads(line) for line in (seed/'scalar-state/operations.jsonl').read_text().splitlines()]
        for operation in operations:
            if operation['stage']!='saved':continue
            stem=operation['method']+'_'+operation['sourceSha256']
            packet=seed/'scalar-state'/(stem+'.pickle')
            source=seed/'scalar-state'/(stem+'.input.pickle')
            if checker.digest(packet)!=operation['sha256'] or checker.digest(packet)!=(packet.with_suffix('.sha256')).read_text().strip():
                raise ValueError(('scalar seed result',stem))
            if checker.digest(source)!=operation['sourceSha256']:raise ValueError(('scalar seed input',stem))
            before,value,known=pickle.loads(packet.read_bytes())
            if before!=pickle.loads(source.read_bytes()):raise ValueError(('scalar seed operand',stem))
            provenance={'kind':'diagnostic_seed','packet':str(packet),'sha256':checker.digest(packet),
                        'input':str(source),'inputSha256':checker.digest(source)}
            key=(operation['method'],before)
            if key in self.cache and self.cache[key][0]!=value:raise ValueError('conflicting exact scalar seeds')
            self.cache[key]=(value,provenance,known)
            self.original.append((operation['method'],before,value,provenance))
        for path in sorted(self.state.glob('call_*/result.json')):
            expression,value,record=self.read_result(path.parent)
            key=(record['method'],expression)
            if key in self.cache and self.cache[key][0]!=value:raise ValueError('conflicting worker scalar')
            self.cache[key]=(value,{'kind':'isolated_worker','directory':str(path.parent),
                                   'resultSha256':record['resultSha256']},{})
        self.progress({'stage':'cache_loaded','originalOperations':len(self.original),'distinctOperands':len(self.cache)})

    def progress(self,record):
        record={**record,'elapsedSeconds':time.monotonic()-self.started,
                'parentPeakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
        with (self.state/'operations.jsonl').open('a') as stream:stream.write(json.dumps(record)+'\n')

    def read_result(self,directory):
        request=json.loads((directory/'request.json').read_text())
        record=json.loads((directory/'result.json').read_text())
        if json.loads((directory/'exit.json').read_text())['exitCode']!=0:
            raise ValueError('worker did not exit successfully')
        if request['environment']!=self.signature['environment'] or not record['environmentStable']:
            raise ValueError('worker result environment')
        if checker.digest(directory/'request.json')!=record['requestSha256']:
            raise ValueError('worker result request pin')
        if checker.digest(directory/'input.pickle')!=record['sourceSha256'] or record['sourceSha256']!=request['sourceSha256']:
            raise ValueError('worker result input pin')
        if checker.digest(directory/'result.pickle')!=record['resultSha256']:
            raise ValueError('worker result payload pin')
        before,value=pickle.loads((directory/'result.pickle').read_bytes())
        if before!=pickle.loads((directory/'input.pickle').read_bytes()) or record['method']!=request['method']:
            raise ValueError('worker operand/method round trip')
        return before,value,record

    def calculate(self,method,expression,force=False):
        key=(method,expression)
        if not force and key in self.cache:
            value,provenance,known=self.cache[key]
            if checker.engine.PHYSICAL_METADATA is not None:
                checker.engine.PHYSICAL_METADATA.dimensions.known.update(known)
            self.progress({'stage':'reused','method':method,'provenance':provenance})
            return value
        body=pickle.dumps(expression,protocol=5);source_sha=hashlib.sha256(body).hexdigest()
        directory=self.state/('call_'+method+'_'+source_sha)
        if directory.exists():
            if not (directory/'result.json').exists():raise ValueError(('unfinished worker operation',str(directory)))
            before,value,record=self.read_result(directory)
            if before!=expression:raise ValueError('existing worker operand')
            self.cache[key]=(value,{'kind':'isolated_worker','directory':str(directory),
                                  'resultSha256':record['resultSha256']},{})
            return value
        directory.mkdir()
        checker.atomic(directory/'input.pickle',body)
        request={'method':method,'sourceSha256':source_sha,'environment':self.signature['environment']}
        write_json(directory/'request.json',request)
        self.progress({'stage':'worker_started','method':method,'directory':str(directory),'inputBytes':len(body)})
        with (directory/'stdout.txt').open('w') as stdout,(directory/'stderr.txt').open('w') as stderr:
            child=subprocess.run([sys.executable,'-X','faulthandler',str(HERE),'--worker',str(directory/'request.json')],
                                 stdout=stdout,stderr=stderr)
        write_json(directory/'exit.json',{'exitCode':child.returncode})
        self.progress({'stage':'worker_exited','method':method,'directory':str(directory),'exitCode':child.returncode})
        if child.returncode!=0:raise RuntimeError(('isolated scalar worker failed',child.returncode,str(directory)))
        before,value,record=self.read_result(directory)
        if before!=expression:raise ValueError('requested worker operand')
        self.cache[key]=(value,{'kind':'isolated_worker','directory':str(directory),
                               'resultSha256':record['resultSha256']},{})
        self.progress({'stage':'worker_saved','method':method,'directory':str(directory),'resource':record})
        return value

    def validate_boundary(self):
        validation=self.state/'boundary.json'
        if validation.exists():
            records=json.loads(validation.read_text())
            for item in records:
                packet=Path(item['comparisonPacket'])
                if checker.digest(packet)!=item['comparisonSha256']:raise ValueError('boundary comparison pin')
                source,before,after,residual=pickle.loads(packet.read_bytes())
                if before!=after or residual!=after-before:raise ValueError('boundary comparison changed')
            return records
        rational=next(item for item in self.original if item[0]=='rational_coefficient' and item[2]!=0)
        expansion=max((item for item in self.original if item[0]=='carrier_expansion'),key=lambda item:Path(item[3]['input']).stat().st_size)
        probe=self.seed.parent/'isolated-scalar.json'
        evidence=json.loads(probe.read_text());packet=self.seed.parent/'isolated-scalar.pickle'
        if checker.digest(packet)!=evidence['payloadSha256'] or evidence['nativeSha256']!=self.signature['environment']['nativeSha256']:
            raise ValueError('isolated probe provenance')
        source,expected=pickle.loads(packet.read_bytes())
        candidates=(rational,expansion,('carrier_expansion',source,expected,{'kind':'isolated_probe','path':str(packet),'sha256':checker.digest(packet)}))
        records=[]
        for i,(method,source,expected,provenance) in enumerate(candidates):
            value=self.calculate(method,source,force=True)
            residual=value-expected
            packet=self.state/f'boundary_{i}.pickle'
            checker.atomic(packet,pickle.dumps((source,expected,value,residual),protocol=5))
            record={'method':method,'baseline':provenance,'structuralEquality':value==expected,
                    'literalResidual':str(residual),'comparisonPacket':str(packet),'comparisonSha256':checker.digest(packet)}
            print(json.dumps(record),flush=True)
            records.append(record)
            if value!=expected:raise ValueError('isolated boundary result differs; comparison emitted')
        write_json(validation,records)
        self.progress({'stage':'boundary_validated','comparisons':len(records)})
        return records


def main():
    parser=argparse.ArgumentParser(add_help=False)
    parser.add_argument('--worker',type=Path)
    parser.add_argument('--seed-diagnostic',type=Path)
    parser.add_argument('--validate-boundary-only',action='store_true')
    args,remaining=parser.parse_known_args()
    if args.worker:
        worker(args.worker);return
    if args.seed_diagnostic is None:raise ValueError('diagnostic seed required')
    base=Path(remaining[remaining.index('--run-directory')+1]);seed=args.seed_diagnostic
    base.mkdir(parents=True,exist_ok=True)
    if not (base/'signature.json').exists():
        if checker.digest(seed/'construction.pickle')!=(seed/'construction.sha256').read_text().strip():
            raise ValueError('seed construction pin')
        for name in ('signature.json','construction.pickle','construction.sha256'):
            shutil.copy2(seed/name,base/name)
        shutil.copytree(seed/'source',base/'source')
    sys.excepthook=sys.__excepthook__
    coordinator=ScalarWorkers(base,seed)
    coordinator.validate_boundary()
    if args.validate_boundary_only:return
    for name in METHODS:
        def calculate(expression,name=name):return coordinator.calculate(name,expression)
        setattr(checker.engine.ClosedCurrentPairing,name,staticmethod(calculate))
    sys.argv=[sys.argv[0],*remaining]
    checker.run()
    coordinator.progress({'stage':'complete','environmentStable':environment()==coordinator.signature['environment']})


if __name__=='__main__':main()
