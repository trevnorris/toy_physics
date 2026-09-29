#!/usr/bin/env python3
"""Saved-return continuation: the same fixed transverse/current point and faces.

No producer import/replay, new root census, map, field or power construction.
This file is inert until the exact independent-review/authorization gate and
progress-dependent containment are verified. Interpretation belongs to a later report.
"""
import argparse
import builtins
from collections import OrderedDict
import importlib
import io
import ast
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import pickle
import resource
import signal
import time
import tempfile
import traceback

ROOT = Path('/var/projects/toy_physics')
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')
ENDS = ('REFERENCE', 'LEFT', 'RIGHT')
DRIVES = ('AMPLITUDE', 'PRESSURE', 'OUTWARD_VELOCITY', 'RELATIVE_MASS_FLUX',
          'AFFINITY', 'BULK_VELOCITY')


def require(value, reason):
    if not value:
        raise ValueError(reason)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()


def route(path):
    path = Path(path)
    return dict(path=str(path), canonicalPath=str(path.resolve(strict=True)),
                bytes=path.stat().st_size, sha256=sha(path))


class IntegrityError(RuntimeError):
    """A source or stored-object pin changed; never a local math fallback."""


def publish(path, payload):
    """Publish a complete fsynced byte string with no overwrite window."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.'+path.name+'.', suffix='.partial', dir=path.parent)
    with os.fdopen(fd, 'wb') as stream:
        stream.write(payload); stream.flush(); os.fsync(stream.fileno())
    os.link(temporary, path)  # Atomic exclusive publication; existing path fails.
    directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)
    os.unlink(temporary)


def save(path, value):
    publish(path, (json.dumps(value, indent=2, allow_nan=False)+'\n').encode())


class OperationBudget(BaseException):
    """Deliberately escapes ordinary algebra/domain exception handlers."""


class NativeDeadline(BaseException):
    """Fatal whole-job deadline; disarm before failure bookkeeping."""


class SavedCodec(pickle.Unpickler):
    """Only pinned data; no producer classes, arbitrary builtins, or eval.

    SymPy expression/matrix classes reconstruct their stored structure. The
    sole permitted SymPy helper reconstructs undefined symbolic functions.
    No solve, derivative, simplify, integral, or producer function is called.
    """
    def find_class(self, module, name):
        import sympy as sp
        import numpy as np
        safe = {
            ('collections', 'OrderedDict'): OrderedDict,
            ('builtins', 'complex'): builtins.complex,
            ('builtins', 'set'): builtins.set,
            ('builtins', 'frozenset'): builtins.frozenset,
            ('builtins', 'slice'): builtins.slice,
            ('numpy', 'ndarray'): np.ndarray,
            ('numpy', 'dtype'): np.dtype,
            ('numpy.core.multiarray', '_reconstruct'): np.core.multiarray._reconstruct,
            ('numpy.core.multiarray', 'scalar'): np.core.multiarray.scalar,
            ('numpy._core.multiarray', '_reconstruct'): np.core.multiarray._reconstruct,
            ('numpy._core.multiarray', 'scalar'): np.core.multiarray.scalar,
            ('numpy.core.numeric', '_frombuffer'): np.core.numeric._frombuffer,
            ('numpy._core.numeric', '_frombuffer'): np.core.numeric._frombuffer,
        }
        if (module, name) in safe:
            return safe[(module, name)]
        if module.startswith('sympy.') and not name.startswith('__'):
            item = getattr(importlib.import_module(module), name)
            if isinstance(item, type) and issubclass(item, (sp.Basic, sp.MatrixBase)):
                return item
            if (module, name) == ('sympy.core.function', '_rebuild_undef'):
                return item
        raise pickle.UnpicklingError('unsupported saved global: ' + module + '.' + name)

    def persistent_load(self, identifier):
        raise pickle.UnpicklingError('persistent IDs are unsupported')


class Journal:
    """Complete operands, atomic content-addressed values, local math failures."""
    def __init__(self, out, deadline, spec):
        self.out, self.deadline, self.spec = out, deadline, spec
        self.prior, self.prior_cache, self.reused = {}, {}, []
        self.prior_order, self.resumed, self.summary_reuse = [], [], []
        self.local_deadline, self.active_path = None, None
        self.armed_timer_kind = None
        self.same = None
        self.records, self.active, self.objects = [], None, {}
        self.progress_sequence = 0
        self.progress_paths = set()
        self.resume_suboperations = {}
        self.restored_suboperations = set()
        self.carrier_entry_number = 0
        def timeout(*_):
            signal.setitimer(signal.ITIMER_REAL, 0)
            raise NativeDeadline('no newly saved result for 3600 seconds')
        signal.signal(signal.SIGALRM, timeout)
        self.arm_native()

    def arm_native(self):
        signal.setitimer(signal.ITIMER_REAL, 0)
        remaining = self.deadline-time.monotonic()
        if remaining <= 0:
            raise NativeDeadline('no newly saved result for 3600 seconds')
        signal.setitimer(signal.ITIMER_REAL, remaining)

    def progress(self, receipt_path, kind):
        # Only call after exclusive, fsynced result publication. No heartbeat.
        receipt = route(receipt_path)
        require(receipt['path'] not in self.progress_paths,'result already counted as progress')
        self.progress_paths.add(receipt['path'])
        self.progress_sequence += 1
        event = dict(sequence=self.progress_sequence,kind=kind,receipt=receipt,
                     utc=datetime.now(timezone.utc).isoformat())
        save(self.out/'progress-events'/('%08d.json'%self.progress_sequence),event)
        payload=(json.dumps(event,indent=2)+'\n').encode()
        target=self.out.parent/'progress.json'
        fd,tmp=tempfile.mkstemp(prefix='.progress.',dir=target.parent)
        with os.fdopen(fd,'wb') as stream:
            stream.write(payload);stream.flush();os.fsync(stream.fileno())
        os.replace(tmp,target)
        self.deadline=time.monotonic()+3600
        self.arm_native()

    def configure_resume(self):
        root=Path(self.spec['resumeRoot'])
        expected={name:self.spec['inputs'][key] for name,key in self.spec['resumeFiles'].items()}
        require({str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()}==set(expected),
                'fixed-point partial file census changed')
        copied={}
        for relative,pin in expected.items():
            require(route(pin['path'])==pin,'fixed-point partial pin changed')
            target=self.out/'fixed-point-prior'/relative
            publish(target,Path(pin['path']).read_bytes())
            copied[relative]=route(target)
            require(copied[relative]['sha256']==pin['sha256'],'partial copy mismatch')
        save(self.out/'fixed-point-prior-copy-index.json',copied)
        index=json.loads((root/'operation-index.json').read_text())
        require(len(index)==14 and all(o['status']=='COMPLETE' for o in index[:11])
                and all(o['status']=='UNRESOLVED_NO_RETRY' for o in index[11:]),'partial journal census')
        self.resume_index={o['name']:o for o in index}
        for relative in expected:
            if '/suboperations/' in relative and relative.endswith('/complete.json'):
                record=json.loads((root/relative).read_text())
                key=(record['operation'],record['name'])
                require(key not in self.resume_suboperations,'duplicate saved suboperation')
                require(record==self.spec['resumeSuboperations'][record['operation']][record['name']],
                        'reviewed partial receipt differs')
                self.resume_suboperations[key]=(relative,record)
        require(len(self.resume_suboperations)==33,'expected 33 complete suboperations')
        self.progress(self.out/'fixed-point-prior-copy-index.json','PRESERVED_FIXED_POINT_RESULTS')

    def resume_has(self, name):
        return (self.active['name'],name) in self.resume_suboperations

    def restore_suboperation(self, name):
        key=(self.active['name'],name)
        require(key in self.resume_suboperations and key not in self.restored_suboperations,
                'missing or repeated saved suboperation')
        relative,record=self.resume_suboperations[key]
        pin=record['value']
        saved_relative=str(Path(pin['path']).relative_to(self.spec['resumeRoot']))
        source_pin=self.spec['inputs'][self.spec['resumeFiles'][saved_relative]]
        require(source_pin==pin,'suboperation object not in pinned input set')
        path=self.out/'fixed-point-prior'/saved_relative
        require(sha(path)==pin['sha256'],'copied suboperation object changed')
        stream=io.BytesIO(path.read_bytes()); value=SavedCodec(stream).load()
        require(not stream.read(1),'trailing saved object data')
        self.objects[id(value)]=(value,route(path))
        target=self.active_path/'suboperations'/name/'complete.json'
        saved=dict(name=name,value=route(path),operation=self.active['name'],
                   execution='RESTORED_PRIOR_COMPLETE_SUBOPERATION',priorReceipt=record)
        save(target,saved)
        self.restored_suboperations.add(key)
        self.progress(target,'RESTORED_PRIOR_COMPLETE_SUBOPERATION')
        return value

    def blob(self, value):
        existing = self.objects.get(id(value))
        if existing is not None and existing[0] is value:
            record = existing[1]
            if route(record['path']) != record:
                raise IntegrityError('previous object changed')
            return record
        payload = pickle.dumps(value, protocol=4)
        digest = hashlib.sha256(payload).hexdigest()
        path = self.out/'objects'/(digest+'.pickle')
        if path.exists():
            if path.stat().st_size != len(payload) or sha(path) != digest:
                raise IntegrityError('existing content-addressed object is incomplete or changed')
        else:
            publish(path, payload)
        record = route(path)
        if record['sha256'] != digest:
            raise IntegrityError('published object hash mismatch')
        self.objects[id(value)] = (value, record)
        return record

    def prior_json(self, relative):
        pin = self.spec['inputs'][self.spec['priorFiles'][relative]]
        require(route(pin['path']) == pin, 'prior JSON changed: '+relative)
        return json.loads(Path(pin['path']).read_text())

    def configure_prior(self):
        root = Path(self.spec['priorRoot'])
        observed = {name:self.spec['inputs'][key] for name,key in self.spec['priorFiles'].items()}
        index = self.prior_json('operation-index.json')
        checkpoint = json.loads(Path(self.spec['inputs']['priorCheckpoint']['path']).read_text())
        complete = [record for record in index if record['status']=='COMPLETE']
        unfinished = [record for record in index if record['status']=='UNRESOLVED_NO_RETRY']
        checks = dict(census={str(path.relative_to(root)) for path in root.rglob('*') if path.is_file()}==set(observed),
            checkpointArtifacts=observed==checkpoint['artifacts'],
            artifactIndex={key:value for key,value in observed.items() if key!='artifact-index.json'}==self.prior_json('artifact-index.json'),
            originalSpecification=self.prior_json('manifest.json')==self.spec['originalSpecification'],
            originalInputs=self.prior_json('posthashes.json')==self.spec['originalSpecification']['inputs'],
            counts=(len(complete)==11 and len(unfinished)==3 and len(index)==14),
            unfinishedNames={record['name'] for record in unfinished}=={label+'-source-bind' for label in ENDS})
        for ordinal,record in enumerate(index):
            folder = 'operations/%04d-%s/' % (ordinal,record['name'])
            key = 'complete.json' if record['status']=='COMPLETE' else 'unresolved.json'
            supplied = self.prior_json(folder+'input-references.json')
            checks['receipt-%02d'%ordinal] = (record==self.prior_json(folder+key)
                and supplied=={name:record[name] for name in ('name','startedUtc','operands')}
                and ('value' in record)==(record['status']=='COMPLETE'))
        save(self.out/'prior-metadata-joins.json',dict(checks=checks,index=index))
        require(all(checks.values()),'prior metadata mismatch; joins saved')
        self.prior = {record['name']:record for record in index}
        self.prior_order = [record['name'] for record in complete]
        require(len(self.prior)==14,'duplicate prior operations')
        # Every original file, including all unfinished operand bundles and
        # nonjournal evidence, is preserved byte-for-byte before restoration.
        copied = {}
        for relative,pin in observed.items():
            require(route(pin['path'])==pin,'prior artifact changed')
            payload = Path(pin['path']).read_bytes()
            require(hashlib.sha256(payload).hexdigest()==pin['sha256'],'prior read changed')
            target = self.out/'prior-complete'/relative
            publish(target,payload)
            copied[relative] = route(target)
            require(copied[relative]['sha256']==pin['sha256'] and copied[relative]['bytes']==pin['bytes'],
                    'prior exact copy mismatch')
        save(self.out/'prior-copy-index.json',copied)

    def prior_value(self, pin):
        allowed = {self.spec['inputs'][key]['sha256']:self.spec['inputs'][key]
                   for key in self.spec['priorRestoreKeys']}
        require(allowed.get(pin['sha256'])==pin,'unselected prior object')
        require(route(pin['path'])==pin,'prior object changed')
        digest = pin['sha256']
        if digest not in self.prior_cache:
            relative = str(Path(pin['path']).relative_to(self.spec['priorRoot']))
            copied = self.out/'prior-complete'/relative
            require(sha(copied)==digest,'copied object changed')
            stream = io.BytesIO(copied.read_bytes())
            value = SavedCodec(stream).load()
            require(not stream.read(1),'trailing bytes after saved object')
            self.prior_cache[digest] = value
            self.objects[id(value)] = (value,route(copied))
        return self.prior_cache[digest]

    def reuse_summary(self, relative, actual=None):
        require(relative in self.spec['priorSummaryFiles'],'unselected completed summary')
        pin = self.spec['inputs'][self.spec['priorFiles'][relative]]
        require(route(pin['path'])==pin,'summary changed')
        raw = Path(pin['path']).read_bytes()
        parsed = json.loads(raw)
        evidence = dict(source=pin,actualArgumentJoin=(parsed==actual) if actual is not None else None)
        save(self.out/'summary-joins'/(relative+'.join.json'),evidence)
        if actual is not None:
            require(evidence['actualArgumentJoin'],'completed nonjournal source/argument join changed')
        publish(self.out/relative,raw)
        self.summary_reuse.append(evidence)
        return parsed

    def checkpoint(self, name, value):
        # Local mathematical timers never cover serialization. Whole-job expiry
        # remains fatal; resuming uses the original end-binding deadline.
        self.arm_native()
        path = self.active_path/'suboperations'/name
        record = dict(name=name,value=self.blob(value),operation=self.active['name'])
        save(path/'complete.json',record)
        if not name.endswith('-input'):
            self.progress(path/'complete.json','NEW_COMPLETE_SUBOPERATION')
        return record

    def resume_local(self):
        self.arm_native()  # Preserve time since last saved result, never reset.

    def op(self, name, function, *args, seconds=20, required=False):
        path = self.out/'operations'/('%04d-%s' % (len(self.records),name))
        path.mkdir(parents=True,exist_ok=False)
        self.active_path = path
        self.active = dict(name=name,operands=[],startedUtc=datetime.now(timezone.utc).isoformat())
        save(path/'started.json',self.active)
        for arg in args: self.active['operands'].append(self.blob(arg))
        save(path/'input-references.json',self.active)
        previous = self.prior.get(name)
        if name in self.resume_index and name.endswith('-source-bind'):
            fixed=self.resume_index[name]
            require([p['sha256'] for p in fixed['operands']]==
                    [p['sha256'] for p in previous['operands']],
                    'fixed-point arguments differ from original joined operands')
            save(path/'fixed-point-incomplete-input-join.json',dict(
                priorFixedPoint=fixed,originalPrior=previous,sourceHashesEqual=True))
        if previous is not None:
            supplied = [self.prior_value(pin) for pin in previous['operands']]
            matches = [self.same(a,b) for a,b in zip(args,supplied)]
            evidence = dict(priorOperation=previous,suppliedOperands=self.active['operands'],
                count=len(args)==len(supplied),matches=matches,
                rule='RECURSIVE_SAVED_STRUCTURAL_IDENTITY_NO_ALGEBRA_OR_REPICKLE_CANONICALIZATION')
            save(path/'prior-argument-join.json',evidence)
            require(evidence['count'] and all(matches),'prior actual operand join changed')
            if previous['status']=='COMPLETE':
                if 'restore-' in name:
                    require(len(args)==1 and isinstance(args[0],str),'saved restore key schema')
                    source_key=args[0]
                    original_pin=self.spec['originalSpecification']['inputs'][source_key]
                    source_join=dict(key=source_key,originalSource=original_pin,
                        currentSource=self.spec['inputs'][source_key],observedSource=route(original_pin['path']))
                    save(path/'prior-restoration-source-join.json',source_join)
                    require(source_join['originalSource']==source_join['currentSource']==source_join['observedSource'],
                            'completed restore actual source route changed')
                require(name==self.prior_order[len(self.reused)],'complete return sequence changed')
                value = self.prior_value(previous['value'])
                self.reused.append(name)
                record = dict(self.active,status='COMPLETE',execution='RESTORED_PRIOR_COMPLETE_RETURN',
                              priorOperation=previous,argumentJoin=evidence,value=self.blob(value))
                save(path/'complete.json',record)
                self.records.append(record)
                self.progress(path/'complete.json','RESTORED_PRIOR_COMPLETE_RETURN')
                self.active=None;self.active_path=None
                return value  # The supplied completed function is never called.
            require(name not in self.resumed and len(self.resumed)<3 and seconds==180,
                    'one bounded attempt per unfinished end binding')
            self.resumed.append(name)
        require(len(self.reused)==11,'all original returns must be reused before new computation')
        remaining = self.deadline-time.monotonic()
        if remaining<=0: raise NativeDeadline('native whole-job deadline')
        self.local_deadline = None  # Old per-call caps superseded only for this approved continuation.
        self.resume_local()
        failed = None
        try:
            value = function(*args)
            signal.setitimer(signal.ITIMER_REAL,0)
        except NativeDeadline:
            signal.setitimer(signal.ITIMER_REAL,0)
            raise
        except OperationBudget:
            signal.setitimer(signal.ITIMER_REAL,0)
            if required: raise
            failed = dict(status='UNRESOLVED',reason='OPERATION_BUDGET',traceback=traceback.format_exc())
        except Exception as error:
            signal.setitimer(signal.ITIMER_REAL,0)
            if required or isinstance(error,(IntegrityError,OSError,pickle.UnpicklingError)): raise
            failed = dict(status='UNRESOLVED',reason='LOCAL_OPERATION_EXCEPTION',
                          exceptionType=type(error).__name__,traceback=traceback.format_exc())
        self.local_deadline=None
        if failed is not None: self.active['localFailureBeforePersistence']=failed
        self.arm_native()
        if failed is not None:
            failed['operationReceipt']=str(path/'unresolved.json')
            record = dict(self.active,status='UNRESOLVED_NO_RETRY',failure=failed)
            save(path/'unresolved.json',record)
        else:
            record = dict(self.active,status='COMPLETE',value=self.blob(value))
            save(path/'complete.json',record)
            self.progress(path/'complete.json','NEW_COMPLETE_OPERATION')
        self.records.append(record);self.active=None;self.active_path=None
        return failed if failed is not None else value


def containment():
    group = next(x[3:] for x in Path('/proc/self/cgroup').read_text().splitlines()
                 if x.startswith('0::'))
    cgroup = Path('/sys/fs/cgroup')/group.lstrip('/')
    limits = {key: (cgroup/key).read_text().strip()
              for key in ('memory.max', 'memory.swap.max', 'pids.max')}
    limits.update(nice=os.getpriority(os.PRIO_PROCESS, 0),
                  affinity=sorted(os.sched_getaffinity(0)),
                  threads={key: os.environ.get(key) for key in THREADS})
    require(limits['memory.max'] == '2147483648' and limits['memory.swap.max'] == '0'
            and limits['pids.max'] == '32' and limits['nice'] >= 15
            and len(limits['affinity']) == 1
            and all(x == '1' for x in limits['threads'].values()),
            'ordinary containment missing; science not imported')
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    def timeout(*_):
        raise NativeDeadline('native whole-job deadline')
    signal.signal(signal.SIGALRM, timeout)
    signal.setitimer(signal.ITIMER_REAL, 3600)
    return limits


def construct(spec, out, journal):
    import sympy as sp
    import numpy as np

    w = sp.Rational(spec['fixedPoint']['omega'])
    v = sp.Rational(spec['fixedPoint']['tangents'][1])
    k, kb = sp.symbols('premiseNormalMomentum premiseLeftNormalMomentum', real=True)
    q, qb = sp.symbols('premiseDepthMomentum premiseLeftDepthMomentum', complex=True)
    params = {name: sp.Rational(value) for name, value in spec['physicalInput']['parameters'].items()}

    def show(value):
        if isinstance(value, dict): return {str(a): show(b) for a, b in value.items()}
        if isinstance(value, (list, tuple)): return [show(x) for x in value]
        if isinstance(value, sp.MatrixBase): return show(value.tolist())
        if isinstance(value, np.ndarray): return show(value.tolist())
        if isinstance(value, complex):
            return dict(real=value.real if np.isfinite(value.real) else str(value.real),
                        imag=value.imag if np.isfinite(value.imag) else str(value.imag))
        if isinstance(value, (float, np.floating)):
            number = float(value)
            return number if np.isfinite(number) else str(value)
        if isinstance(value, (np.integer, np.bool_)): return value.item()
        if isinstance(value, sp.Basic): return str(value)
        return value

    def same(left,right):
        if left is right: return True
        if isinstance(left,np.ndarray) or isinstance(right,np.ndarray):
            return (isinstance(left,np.ndarray) and isinstance(right,np.ndarray)
                    and left.shape==right.shape and left.dtype==right.dtype
                    and bool(np.array_equal(left,right)))
        if isinstance(left,dict) or isinstance(right,dict):
            return (isinstance(left,dict) and isinstance(right,dict) and left.keys()==right.keys()
                    and all(same(left[key],right[key]) for key in left))
        if isinstance(left,(list,tuple)) or isinstance(right,(list,tuple)):
            return type(left) is type(right) and len(left)==len(right) and all(same(a,b) for a,b in zip(left,right))
        if isinstance(left,sp.MatrixBase) or isinstance(right,sp.MatrixBase):
            return (isinstance(left,sp.MatrixBase) and isinstance(right,sp.MatrixBase)
                    and left.shape==right.shape and all(a==b for a,b in zip(left,right)))
        return type(left) is type(right) and bool(left==right)
    journal.same=same

    def evidence(name,value):
        record=journal.checkpoint(name,value)
        journal.resume_local()
        return record

    def op(name, function, *args, seconds=20, required=False, render=True):
        value=journal.op(name,function,*args,seconds=seconds,required=required)
        if name in journal.reused:
            prior=journal.prior[name]
            old=journal.prior_json(name+'.json')
            if render:
                journal.reuse_summary(name+'.json',show(value))
            else:
                require(old.get('completeObjectReceipt')==prior,'original restore summary receipt join')
                journal.reuse_summary(name+'.json')
        else:
            save(out/(name+'.json'),show(value) if render else dict(operation=name,
                 completeObjectReceipt=journal.records[-1],wholeRestoredPayloadRendered=False))
        return value

    def restored(name, key):
        return op(name, restore, key, required=True, render=False)

    def restore(name):
        record = spec['inputs'][name]
        if route(record['path']) != record: raise IntegrityError('changed source '+name)
        with Path(record['path']).open('rb') as stream:
            return SavedCodec(stream).load()

    def flat(value):
        if isinstance(value, sp.MatrixBase): return list(value)
        if isinstance(value, dict): return [z for x in value.values() for z in flat(x)]
        if isinstance(value, (tuple, list)): return [z for x in value for z in flat(x)]
        return [value]

    def zero(value): return all(sp.cancel(x) == 0 for x in flat(value))
    def clean(matrix): return sp.ImmutableMatrix(matrix).applyfunc(sp.cancel)

    def raw_denominators(values):
        return tuple(dict.fromkeys(
            [sp.fraction(sp.together(value))[1] for value in values]+
            [power.base for value in values for power in value.atoms(sp.Pow)
             if power.exp.is_negative is True]))

    # This mapping is derived from ChannelInput's pinned source expression,
    # not a new physical closure. The expression reader admits scalar input
    # arithmetic only; no source function or producer is imported or called.
    engine_tree=ast.parse(Path(spec['inputs']['sourceEngine']['path']).read_text())
    channel=next(node for node in engine_tree.body if isinstance(node,ast.ClassDef) and node.name=='ChannelInput')
    initializer=next(node for node in channel.body if isinstance(node,ast.FunctionDef) and node.name=='__init__')
    origin_node=next(node for node in initializer.body if isinstance(node,ast.Assign)
        and any(isinstance(target,ast.Attribute) and target.attr=='origin' for target in node.targets))
    eta_node=next(node for node in initializer.body if isinstance(node,ast.Assign)
        and any(isinstance(target,ast.Name) and target.id=='eta' for target in node.targets))
    def input_arithmetic(node,env):
        if isinstance(node,ast.Constant): return node.value
        if isinstance(node,ast.Name): return env[node.id]
        if isinstance(node,ast.Attribute): return input_arithmetic(node.value,env)[node.attr]
        if isinstance(node,ast.Subscript): return input_arithmetic(node.value,env)[input_arithmetic(node.slice,env)]
        if isinstance(node,ast.BinOp):
            left,right=input_arithmetic(node.left,env),input_arithmetic(node.right,env)
            if isinstance(node.op,ast.Mult): return left*right
            if isinstance(node.op,ast.Div): return left/right
        raise ValueError('unsupported pinned grade-origin expression')
    origin_env={'self':{'parameters':params}}
    origin_env['eta']=input_arithmetic(eta_node.value,origin_env)
    physical_grades={ast.literal_eval(key.slice):input_arithmetic(value,origin_env)
                     for key,value in zip(origin_node.value.keys,origin_node.value.values)}
    require(set(physical_grades)=={'eta_bg','sigma_W'},'source grade origin schema')
    save(out/'physical-grade-origin-source.json',dict(source=spec['inputs']['sourceEngine'],
         assignment=ast.unparse(origin_node),values=show(physical_grades)))

    def bind(value,label,extra=None,carry=(),receipt=None):
        # No grade is created: carry contains existing source symbol names only.
        carried=set(carry)
        require(carried<={'epsilon_shape','eta_bg','sigma_W'},'explicit carried-symbol allowlist')
        endpoints=value.xreplace(uniform['profileEndpoints'])
        mapping=dict(extra or {})
        for symbol in endpoints.free_symbols-set(mapping):
            name=symbol.name
            if name in carried: mapping[symbol]=symbol
            elif name in ('omega','s11cdFrequency'): mapping[symbol]=w
            elif name in ('s11cdTangentialMomentum1','s11cdTangentialMomentum2'): mapping[symbol]=params[name]
            elif name=='s11cdSpectralNormalMomentum': mapping[symbol]=k
            elif name in physical_grades: mapping[symbol]=sp.S.Zero if label=='REFERENCE' else physical_grades[name]
            elif name in params: mapping[symbol]=params[name]
            elif symbol not in (k,kb,q,qb): raise ValueError('unsupported source symbol '+name)
        bound=endpoints.xreplace(mapping)
        if receipt is not None:
            receipt.update(original=value,afterEndpoints=endpoints,actualSubstitutions=mapping,
                carriedPresent=sorted(symbol.name for symbol in endpoints.free_symbols if symbol.name in carried),
                remainingSymbols=tuple(sorted(bound.free_symbols,key=str)),specialized=bound)
        require(not bound.atoms(sp.Integral,sp.Limit,sp.Derivative),'unresolved source atom')
        return bound

    def grade(matrix):
        epsilon = [s for s in matrix.free_symbols if s.name == 'epsilon_shape']
        if not epsilon and zero(matrix): return matrix, matrix
        require(len(epsilon) == 1, 'quadratic current amplitude carrier')
        e = epsilon[0]
        def extract(entry):
            index=journal.carrier_entry_number; journal.carrier_entry_number+=1
            name='carrier-entry-%06d'%index
            evidence(name+'-input',dict(entry=entry,epsilon=e))
            value=sp.Poly(entry,e).nth(2)
            evidence(name+'-return',dict(entry=entry,epsilon=e,coefficient=value))
            return value
        coefficient = matrix.applyfunc(extract)
        return coefficient, matrix-e**2*coefficient

    uniform = restored('restore-uniform','uniform')
    reduction = restored('restore-branch-context','reduction')
    units = restored('restore-unit-context','unitContext')
    common = restored('restore-uniform-common','uniformCommon')
    require(all(isinstance(x, dict) and x.get('reason') != 'OPERATION_BUDGET'
                for x in (uniform, reduction, units)), 'required source restoration incomplete')
    journal.reuse_summary('unit-context.json',show(dict(unitFrame=spec['physicalInput']['unit_frame'],
         fields=units['fieldUnits'],current=units['currentUnit'],strong=uniform['units']['strong'])))
    face_index = json.loads(Path(spec['inputs']['faceIndex']['path']).read_text())
    source_tree = ast.parse(Path(spec['inputs']['sourceEngine']['path']).read_text())
    acoustic_class = next(node for node in source_tree.body
                          if isinstance(node,ast.ClassDef) and node.name=='ClosedAcousticEnergy')
    acoustic_method = next(node for node in acoustic_class.body
                           if isinstance(node,ast.FunctionDef) and node.name=='construct')
    sign_loops = [node for node in ast.walk(acoustic_method)
                  if isinstance(node,ast.For) and isinstance(node.target,ast.Name)
                  and node.target.id=='sign' and isinstance(node.iter,ast.Tuple)]
    require(len(sign_loops)==1, 'source face-order schema')
    face_order = ast.literal_eval(sign_loops[0].iter)
    law_nodes = {name:next(node for node in acoustic_method.body if isinstance(node,ast.Assign)
                 and any(isinstance(target,ast.Name) and target.id==name for target in node.targets))
                 for name in ('pressure','velocity','local_flux','harmonic')}
    acoustic_law_provenance = {name:dict(line=node.lineno,assignment=ast.unparse(node)) for name,node in law_nodes.items()}
    journal.reuse_summary('acoustic-law-source.json',dict(source=spec['inputs']['sourceEngine'],laws=acoustic_law_provenance))
    journal.reuse_summary('native-face-order-source.json',dict(order=list(face_order),
        source=spec['inputs']['sourceEngine'],sourceLine=sign_loops[0].lineno,
        expression=ast.unparse(sign_loops[0].iter)))

    def branch_inventory(state):
        equations, mapping = state['branch_equations'], state['branch_map']
        residuals = [equation.rhs-mapping[equation.lhs] for equation in equations]
        return dict(equations=equations, mapping=mapping, residuals=residuals,
                    allEquationKeysPresent=all(eq.lhs in mapping for eq in equations),
                    joins=zero(residuals),scope='SUPPLIED_BRANCH_INVENTORY_INTEGRITY_NOT_INDEPENDENT_DERIVATION', tangents=state['tangents'],
                    groups=state['momentum_groups'], normalMap=state['normal_map'])

    branches = op('source-branch-map-joins', branch_inventory, reduction['reductionState'])
    require(branches.get('joins') is True, 'saved branch equation/map inconsistency')

    def native_acoustic_laws(phi,time_coordinate,depth_coordinate):
        env={'rho':params['rho_m'],'phi':phi,'positions':(depth_coordinate,),
             'r':{'t':time_coordinate}}
        def evaluate(node,values):
            if isinstance(node,ast.Name): return values[node.id]
            if isinstance(node,ast.Constant): return node.value
            if isinstance(node,ast.Attribute): return evaluate(node.value,values)[node.attr]
            if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub): return -evaluate(node.operand,values)
            if isinstance(node,ast.BinOp) and isinstance(node.op,ast.Mult):
                return evaluate(node.left,values)*evaluate(node.right,values)
            if isinstance(node,ast.ListComp):
                require(len(node.generators)==1 and not node.generators[0].ifs,'source velocity list schema')
                generator=node.generators[0]
                require(isinstance(generator.target,ast.Name),'source velocity coordinate schema')
                return [evaluate(node.elt,{**values,generator.target.id:item})
                        for item in evaluate(generator.iter,values)]
            if isinstance(node,ast.Call) and isinstance(node.func,ast.Attribute) and isinstance(node.func.value,ast.Name):
                require(node.func.value.id=='sp' and not node.keywords,'source acoustic call schema')
                operands=[evaluate(arg,values) for arg in node.args]
                if node.func.attr=='diff': return sp.diff(*operands)
                if node.func.attr=='ImmutableMatrix': return sp.ImmutableMatrix(*operands)
            raise ValueError('unsupported source acoustic law node')
        for name in ('pressure','velocity','local_flux'):
            env[name]=evaluate(law_nodes[name].value,env)
        return env['pressure'],env['velocity'][0],env['local_flux'][0]

    def carrier_law_check(label, result, mapping, depth):
        time_coordinate, depth_coordinate = sp.symbols('carrierTime carrierOutwardDepth',real=True)
        harmonic = sp.Symbol('carrierHarmonic',nonzero=True)
        epsilon = sp.Symbol('carrierRealFieldScale',real=True)
        minus_amplitude,plus_amplitude = sp.symbols('carrierConjugateAmplitude carrierAmplitude')
        plus = plus_amplitude*harmonic*sp.exp(sp.I*(q*depth_coordinate-w*time_coordinate))
        minus = minus_amplitude/harmonic*sp.exp(-sp.I*(qb*depth_coordinate-w*time_coordinate))
        real_ansatz = epsilon*(plus+minus)/2
        # Source-native pressure and outward velocity laws, applied to one
        # real harmonic carrier pair, not independent counterpropagating waves.
        pressure,velocity,flux=native_acoustic_laws(real_ansatz,time_coordinate,depth_coordinate)
        product=sp.expand(flux)
        average = product.coeff(harmonic,0)
        normalization = epsilon**2*minus_amplitude*plus_amplitude
        normalized_flux = sp.cancel(average/normalization).subs(depth_coordinate,0)
        plus_pressure,plus_velocity,_=native_acoustic_laws(plus,time_coordinate,depth_coordinate)
        pressure_scale=sp.cancel(plus_pressure/plus)
        velocity_scale=sp.cancel(plus_velocity/plus)
        saved_pressure = bind(result['OPEN_BULK_PRESSURE_SCALES'][0],label,mapping)
        saved_velocity = bind(result['OPEN_BULK_VELOCITY_SCALES'][0],label,mapping)
        minus_pressure,minus_velocity,_=native_acoustic_laws(minus,time_coordinate,depth_coordinate)
        minus_pressure_scale=sp.cancel(minus_pressure/minus)
        minus_velocity_scale=sp.cancel(minus_velocity/minus)
        saved_minus_pressure = bind(result['OPEN_BULK_PRESSURE_SCALES'][1],label,mapping)
        saved_minus_velocity = bind(result['OPEN_BULK_VELOCITY_SCALES'][1],label,mapping)
        residuals = dict(pressure=sp.cancel(saved_pressure-pressure_scale),
                         velocity=sp.cancel(saved_velocity-velocity_scale),
                         minusPressure=sp.cancel(saved_minus_pressure-minus_pressure_scale),
                         minusVelocity=sp.cancel(saved_minus_velocity-minus_velocity_scale),
                         depthFlux=sp.cancel(depth-normalized_flux))
        return dict(scope='SOURCE_LAW_AND_REAL_CARRIER_NORMALIZATION_CHECK_NOT_INDEPENDENT_CLOSURE',
            provenance=acoustic_law_provenance,plus=plus,conjugateLeg=minus,
            timeRates=(sp.cancel(sp.diff(plus,time_coordinate)/plus),sp.cancel(sp.diff(minus,time_coordinate)/minus)),
            realAnsatz=real_ansatz,pressure=pressure,velocity=velocity,product=product,
            zeroHarmonicAverage=average,amplitudeNormalization=normalization,
            normalizedFlux=normalized_flux,pressureScale=pressure_scale,velocityScale=velocity_scale,
            savedPressure=saved_pressure,savedVelocity=saved_velocity,savedFlux=depth,
            minusPressureScale=minus_pressure_scale,minusVelocityScale=minus_velocity_scale,
            savedMinusPressure=saved_minus_pressure,savedMinusVelocity=saved_minus_velocity,
            residuals=residuals,joins=zero(residuals))

    def bound_end(label, packet, pair_tuple):
        pair, dimensions = pair_tuple  # Authoritative end_pairing_check.py:165-166.
        result = pair['result']
        saved_faces = result['FACE_LEG_OBJECTS']
        source_addresses = [(index,leg) for index,face in enumerate(saved_faces)
                            for leg in range(len(face))]
        expected_addresses = [(0,0),(0,1),(1,0),(1,1)]
        source_face_census = dict(sourceFaceOrder=face_order,faceCount=len(saved_faces),
            legCounts=[len(face) for face in saved_faces],addresses=source_addresses,
            expectedAddresses=expected_addresses,expectedOrientations=(-1,1),
            complete=(len(face_order)==2 and set(face_order)=={-1,1} and len(saved_faces)==2
                      and len(source_addresses)==4 and set(source_addresses)==set(expected_addresses)))
        if not source_face_census['complete']:
            return dict(end=label,checks={'faceLegCensus':False},sourceFaceLegCensus=source_face_census,
                        savedFaceLegObjects=saved_faces,status='UNRESOLVED',reason='FACE_LEG_CENSUS_INCOMPLETE')
        wl, wr = result['FREQUENCY_LEGS']; kl, kr = result['NORMAL_LEGS']
        ql, qr = result['BULK_LEGS']
        mapping={wl:w,wr:w,kl:kb,kr:k,ql:qb,qr:q}
        matrix_specialization={}
        source_saved=journal.restore_suboperation('source-specialization')
        require(journal.same(source_saved['mapping'],mapping)
                and journal.same(source_saved['matrix']['original'],result['CLOSED_PENCIL_LEGS'][0])
                and journal.same(source_saved['sourceFaceCensus'],source_face_census),
                'saved source specialization argument join')
        matrix_specialization=source_saved['matrix']
        raw_matrix=matrix_specialization['specialized']
        wave_saved=journal.restore_suboperation('wave-polynomial')
        transverse=journal.restore_suboperation('transverse-restriction')
        require(journal.same(transverse['rawMatrix'],raw_matrix)
                and journal.same(transverse['wave'],wave_saved['wave'])
                and journal.same(transverse['branchResiduals'],packet['branchResiduals']),
                'saved transverse/source/branch bundle join')
        matrix=transverse['matrix']; wave=wave_saved['wave']
        original=transverse['original']; source_join=transverse['sourceJoin']
        wave_join=transverse['waveJoin']; lift=transverse['lift']; gram=transverse['gram']
        raw_weak=transverse['rawWeak']; raw_restricted=transverse['rawRestricted']
        invariant=transverse['invariantResidual']; common_curl_join=transverse['suppliedCurlJoin']
        require(wave_saved['degree']==2 and wave_saved['coefficientOne']==0,'saved depth polynomial shape')
        # These small context intermediates were not published in execution 4.
        # Reconstruct only this missing context; the three saved blocks and
        # every completed current-carrier result above/below are restored.
        q2=sp.cancel(-wave_saved['polynomial'].subs(q,0)/wave_saved['coefficientTwo'])
        oldq=packet['radical']; relation=bind(packet['originalRelation'],label,{oldq:oldq})
        scale=sp.sqrt(sp.cancel(sp.Poly(relation,oldq).nth(2)/wave_saved['coefficientTwo']))
        gauge=tuple(bind(value,label) for value in common['gauge'])
        longitudinal_probes=tuple(sp.ImmutableMatrix.vstack(value,sp.zeros(2,1)) for value in gauge)
        raw_gram=lift.H*lift
        weak,restricted=clean(raw_weak),clean(raw_restricted)
        raw_t_denominators=tuple(dict.fromkeys(
            [sp.fraction(sp.together(x))[1] for x in (*raw_weak,*raw_restricted)] +
            [power.base for x in (*raw_weak,*raw_restricted) for power in x.atoms(sp.Pow)
             if power.exp.is_negative is True]))
        evidence('recovered-unsaved-context',dict(q2=q2,relation=relation,scale=scale,gauge=gauge,
            longitudinalProbes=longitudinal_probes,rawGram=raw_gram,weak=weak,
            restricted=restricted,rawTDenominators=raw_t_denominators))
        currents,grade_residuals,current_specializations={},{},{}
        for name in ('SLAB_CURRENT_MATRIX','BULK_NORMAL_CURRENT_DENSITY_MATRIX',
                     'BULK_DEPTH_CURRENT_MATRIX','INTERFACE_POWER_MATRIX'):
            if journal.resume_has(name+'-specialized-input'):
                specialization=journal.restore_suboperation(name+'-specialized-input')
                require(journal.same(specialization['original'],result[name]),'saved current source join')
                specialized=specialization['specialized']
            else:
                specialization={}
                specialized=bind(result[name],label,mapping,
                    carry=('epsilon_shape','eta_bg','sigma_W'),receipt=specialization)
                evidence(name+'-specialized-input',specialization)
            if journal.resume_has(name+'-carrier-return'):
                carrier_saved=journal.restore_suboperation(name+'-carrier-return')
                coefficient=carrier_saved['coefficient']; residual=carrier_saved['reconstructionResidual']
            else:
                coefficient,residual=grade(specialized)
                evidence(name+'-carrier-return',dict(coefficient=coefficient,reconstructionResidual=residual,
                    remainingCoefficientSymbols=tuple(sorted(coefficient.free_symbols,key=str))))
            if journal.resume_has(name+'-end-specialization'):
                end_specialization=journal.restore_suboperation(name+'-end-specialization')
                require(journal.same(end_specialization['original'],coefficient),'saved end coefficient join')
                currents[name]=end_specialization['specialized']
            else:
                end_specialization={}
                currents[name]=bind(coefficient,label,receipt=end_specialization)
                evidence(name+'-end-specialization',end_specialization)
            grade_residuals[name]=residual
            current_specializations[name]=dict(carrier=specialization,end=end_specialization)
        amplitude_legs = tuple(tuple(sorted((symbol for symbol in dimensions
            if getattr(symbol,'name','').startswith('s11cdCurrent'+side+'Amplitude')),
            key=lambda symbol:symbol.name)) for side in ('Plus','Minus'))
        evidence('harmonic-amplitude-columns',dict(columns=amplitude_legs,counts=[len(x) for x in amplitude_legs]))
        require(all(len(amplitudes)==5 for amplitudes in amplitude_legs),
                'five physical amplitude columns in each harmonic leg')
        rows, reconstruction, raw_face_expressions = [], [], []
        flat_exterior = []
        for face_index_value,face in enumerate(result['FACE_LEG_OBJECTS']):
            require(len(face)==2,'both saved harmonic face legs required')
            face_rows,face_residuals,face_expressions = [],[],[]
            for leg,amplitudes in enumerate(amplitude_legs):
                row,residual,expressions = {},{},{}
                for name in DRIVES:
                    expression = face[leg][name]
                    coefficients=sp.ImmutableMatrix(1,5,lambda i,j:sp.diff(expression,amplitudes[j]))
                    row[name]=bind(coefficients,label,mapping)
                    residual[name] = expression-(coefficients*sp.ImmutableMatrix(amplitudes))[0]
                    expressions[name] = expression
                face_rows.append(row);face_residuals.append(residual);face_expressions.append(expressions)
                pressure_scale = bind(result['OPEN_BULK_PRESSURE_SCALES'][leg],label,mapping)
                velocity_scale = bind(result['OPEN_BULK_VELOCITY_SCALES'][leg],label,mapping)
                flat_exterior.append(dict(face=face_index_value,harmonicLeg=leg,
                    pressure=row['PRESSURE'],amplitude=row['AMPLITUDE'],bulkVelocity=row['BULK_VELOCITY'],
                    pressureScale=pressure_scale,velocityScale=velocity_scale,
                    pressureResidual=clean(row['PRESSURE']-pressure_scale*row['AMPLITUDE']),
                    velocityResidual=clean(row['BULK_VELOCITY']-velocity_scale*row['AMPLITUDE'])))
            rows.append(tuple(face_rows));reconstruction.append(tuple(face_residuals))
            raw_face_expressions.append(tuple(face_expressions))
        depth_source=result['OPEN_BULK_CURRENT_COEFFICIENTS'][3]
        free_amps={symbol:sp.S.One for symbol in depth_source.free_symbols
                   if symbol.name in ('s11cdAcousticLeftAmplitude','s11cdAcousticRightAmplitude')}
        depth_specialization={}
        depth_raw=bind(depth_source,label,{**mapping,**free_amps},
                       carry=('epsilon_shape','eta_bg','sigma_W'),receipt=depth_specialization)
        evidence('open-depth-specialized-input',depth_specialization)
        depth_matrix,depth_residual=grade(sp.ImmutableMatrix([depth_raw]))
        evidence('open-depth-carrier-return',dict(coefficient=depth_matrix,reconstructionResidual=depth_residual))
        depth_end_specialization={}
        depth=bind(depth_matrix[0],label,receipt=depth_end_specialization)
        evidence('open-depth-end-specialization',depth_end_specialization)
        carrier = carrier_law_check(label,result,mapping,depth)
        denominator_operands = [*raw_matrix,*raw_restricted,*raw_weak,*raw_gram,
            *[value for current in currents.values() for value in current],
            *[value for face in rows for row in face for coefficients in row.values() for value in coefficients]]
        raw_source_values = [*result['CLOSED_PENCIL_LEGS'][0],
            *[value for name in currents for value in result[name]],
            *[expression for face in raw_face_expressions for leg in face for expression in leg.values()]]
        carried_symbols = {symbol:symbol for columns in amplitude_legs for symbol in columns}
        carried_symbols.update({symbol:symbol for value in raw_source_values for symbol in value.free_symbols
                                if symbol.name=='epsilon_shape'})
        saved_source_factors = raw_denominators(raw_source_values)
        bound_saved_source_factors = tuple(bind(factor,label,{**mapping,**carried_symbols})
                                          for factor in saved_source_factors)
        denominators = tuple(dict.fromkeys((*raw_denominators(denominator_operands),*bound_saved_source_factors)))
        coupling=[uniform['records'][label]['coupling'][name] for name in ('TH','HT')]
        source_names={symbol.name for symbol in packet['originalAlgebraic'].free_symbols}
        specialization=dict(sourceFrequencyPresent='omega' in source_names,
            sourceTangentsPresent=all(name in source_names for name in
                ('s11cdTangentialMomentum1','s11cdTangentialMomentum2')),
            requestedFrequency=w,requestedTangents=tuple(params[name] for name in
                ('s11cdTangentialMomentum1','s11cdTangentialMomentum2')),
            originalSource=packet['originalAlgebraic'],actualBoundSource=original,
            actualMatrix=raw_matrix,sourceJoin=source_join,waveJoin=wave_join,
            remainingMatrixSymbols=tuple(sorted(matrix.free_symbols,key=str)),
            exactPoint=(w==params['omega'] and tuple(spec['fixedPoint']['tangents'])==
                tuple(str(params[name]) for name in ('s11cdTangentialMomentum1','s11cdTangentialMomentum2'))))
        evidence('source-current-face-residuals',dict(specialization=specialization,
            currentSpecializations=current_specializations,gradeResiduals=grade_residuals,
            depthEndSpecialization=depth_end_specialization,
            depthReconstructionResidual=depth_residual,carrierLaw=carrier,
            driveReconstruction=reconstruction,flatExteriorJoins=flat_exterior,
            sourceCoupling=coupling,pairingBranch=result['SOURCE_BRANCH_JOINS']))
        # CLOSED_PENCIL_LEGS[0] is the right single-leg pencil; only k and q
        # may remain there. Bilinear currents retain both normal/depth legs.
        checks=dict(faceLegCensus=source_face_census['complete'],sourceJoin=zero(source_join),
            waveJoin=zero(wave_join),suppliedCurlJoin=zero(common_curl_join),
            carrierLawAndNormalization=carrier['joins'],originalBranch=zero(packet['branchResiduals']),
            pairingBranch=zero(result['SOURCE_BRANCH_JOINS']),savedCoupling=zero(coupling),
            invariantT=zero(invariant),weakHermitian=zero(weak-weak.H),grade=zero(grade_residuals),
            depthGrade=zero(depth_residual),driveLinear=zero(reconstruction),
            flatExteriorJoin=zero([item[key] for item in flat_exterior for key in ('pressureResidual','velocityResidual')]),
            sourceFrequencyProvenance=specialization['sourceFrequencyPresent'],
            sourceTangentProvenance=specialization['sourceTangentsPresent'],
            exactPoint=specialization['exactPoint'],declaredSymbols=not matrix.free_symbols-{k,q})
        return dict(end=label,checks=checks,sourceFaceLegCensus=source_face_census,originalSource=packet['originalAlgebraic'],
            originalCurrentResult=result,sourceLegMapping=mapping,sourceSpecialization=specialization,
            currentSpecializations=current_specializations,
            depthEndSpecialization=depth_end_specialization,
            originalRelation=packet['originalRelation'], sourceJoin=source_join, waveJoin=wave_join,
            savedSourceBranch=packet['branchResiduals'], pairingSourceBranch=result['SOURCE_BRANCH_JOINS'],
            sourceCoupling=coupling, matrix=matrix, rawMatrix=raw_matrix, wave=wave, q2=q2, scale=scale,
            lift=lift, gram=gram, suppliedCurlJoin=common_curl_join, savedGauge=gauge,
            longitudinalDirectionProbes=longitudinal_probes, rawWeak=raw_weak, rawRestricted=raw_restricted,
            rawTDenominators=raw_t_denominators, chart=sp.factor(gram.det()), restricted=restricted, weak=weak,
            invariantResidual=invariant, weakHermitianResidual=weak-weak.H,
            currents=currents, gradeResiduals=grade_residuals, drives=rows,
            harmonicAmplitudeColumns=amplitude_legs,rawSavedFaceExpressions=raw_face_expressions,
            sourceDenominatorOperands=denominator_operands,
            rawSavedSourceDenominatorFactors=saved_source_factors,
            boundSavedSourceDenominatorFactors=bound_saved_source_factors,
            driveReconstruction=reconstruction, flatExteriorJoins=flat_exterior, sourceDenominators=denominators,
            outgoingDepthCoefficient=depth, carrierLawCheck=carrier,
            flatExteriorJoinScope='SOURCE_INTERNAL_CONSISTENCY_NOT_INDEPENDENT_CLOSURE',dimensions=dimensions)

    def face_address_census(records):
        addresses = [(item['faceOrdinal'],item['harmonicLeg']) for item in records]
        oriented = [(face_order[index] if 0<=index<len(face_order) else None,leg)
                    for index,leg in addresses]
        expected = [(0,0),(0,1),(1,0),(1,1)]
        expected_oriented = [(-1,0),(-1,1),(1,0),(1,1)]
        return dict(recordCount=len(records),distinctAddressCount=len(set(addresses)),
            addresses=addresses,orientedAddresses=oriented,sourceFaceOrder=face_order,
            expectedAddresses=expected,expectedOrientedAddresses=expected_oriented,
            complete=(len(face_order)==2 and set(face_order)=={-1,1} and len(records)==4
                      and set(addresses)==set(expected) and set(oriented)==set(expected_oriented)))

    def face_source_row(bound,face,leg,name,expression=None):
        source=bound['rawSavedFaceExpressions'][face][leg][name] if expression is None else expression
        amplitudes=bound['harmonicAmplitudeColumns'][leg]
        coefficients=sp.ImmutableMatrix(1,len(amplitudes),lambda i,j:sp.diff(source,amplitudes[j]))
        return bind(coefficients,bound['end'],bound['sourceLegMapping'])

    def face_omission(bound,face,leg,name,slot,term=None):
        source=bound['rawSavedFaceExpressions'][face][leg][name]
        amplitude=bound['harmonicAmplitudeColumns'][leg][slot]
        # A slot deletion acts on the native amplitude input; a term deletion
        # removes an actual additive term of that same source expression.
        changed=source.xreplace({amplitude:sp.S.Zero}) if term is None else source-term
        row=face_source_row(bound,face,leg,name,changed)
        return dict(sourceAddress=('FACE_LEG_OBJECTS',face,leg,name),source=source,
            amplitude=amplitude,omittedSourceTerm=term,changedSource=changed,recomputedRow=row)

    def faces(bound):
        plus_lift = bound['lift']; records = []
        literals = []
        if bound['end'] == 'REFERENCE':
            for item in face_index['records']:
                if item['tag'] != 'PY_S11CC2_FOLD_SYMBOL_MAP_LAB_HELD_RHO4_CONSTANT': continue
                for identity in item['velocityIdentifications']:
                    text = identity['savedValueSrepr']
                    if hashlib.sha256(text.encode()).hexdigest() != identity['literalSha256']:
                        raise IntegrityError('c2 index literal changed')
                    expression = sp.sympify(text)
                    field = next(symbol for symbol in expression.free_symbols if symbol.name=='e_W_t')
                    literals.append(dict(indexFaceLabel=item['face'],source=expression,
                        coefficientsByHarmonicLeg=tuple(bind(expression,'REFERENCE',{field:rate})
                            for rate in bound['carrierLawCheck']['timeRates']),literalSha256=identity['literalSha256']))
        for index,face_rows in enumerate(bound['drives']):
            for leg,row in enumerate(face_rows):
                lift=plus_lift if leg==0 else sp.conjugate(plus_lift).xreplace({k:kb})
                contractions = {name:clean(row[name]*lift) for name in DRIVES}
                loaded_controls = []
                for name in DRIVES:
                    urow = row[name][:,:3]
                    if zero(urow):
                        loaded_controls.append(dict(drive=name,status='NOT_APPLICABLE_NO_U_DEPENDENCE',
                            applicableCount=0,responsiveCount=0,nonapplicableCount=lift.cols,unresolvedCount=0,
                            originalRow=row[name],uRow=urow,baseline=contractions[name]))
                        continue
                    source=bound['rawSavedFaceExpressions'][index][leg][name]
                    amplitudes=bound['harmonicAmplitudeColumns'][leg]
                    source_terms=sp.Add.make_args(sp.expand(source))
                    terms=[]
                    for term in source_terms:
                        for slot in range(3):
                            if not term.has(amplitudes[slot]): continue
                            coefficient=bind(sp.diff(term,amplitudes[slot]),bound['end'],bound['sourceLegMapping'])
                            terms.append((slot,term,coefficient))
                    columns=[]
                    for column in range(lift.cols):
                        contributions=[(slot,term,sp.cancel(coefficient*lift[slot,column]))
                                       for slot,term,coefficient in terms]
                        rational=all(value.is_rational_function(k,kb,q,qb) is True for _,_,value in contributions)
                        loaded=[item for item in contributions if item[2]!=0] if rational else []
                        if not loaded:
                            structural=rational and all(value==0 for _,_,value in contributions)
                            columns.append(dict(column=column,contributions=contributions,
                                applicable=False if structural else None,
                                status='NOT_APPLICABLE_STRUCTURALLY_UNLOADED_COLUMN' if structural else 'UNRESOLVED'))
                            continue
                        slot,term,contribution=loaded[0]
                        route_control=face_omission(bound,index,leg,name,slot,term)
                        before=row[name]*lift[:,column]
                        after=route_control['recomputedRow']*lift[:,column]
                        raw_movement=after-before;movement=clean(raw_movement)
                        raw_control=dict(column=column,slot=slot,sourceRoute=route_control,
                            sourceContributions=contributions,originalRow=row[name],
                            fixedLoadedLift=lift[:,column],before=before,after=after,
                            rawMovement=raw_movement,reducedMovement=movement)
                        evidence('face-%d-leg-%d-%s-column-%d'%(index,leg,name,column),raw_control)
                        columns.append(dict(raw_control,applicable=True,
                            status='RATIONAL_NONIDENTITY_AT_FIXED_FREQUENCY_TANGENTS' if not zero(movement) else 'UNRESOLVED',
                            normalDepthPointwiseNonzeroClaim=False))
                    applicable_count=sum(item['applicable'] is True for item in columns)
                    responsive_count=sum(item['status']=='RATIONAL_NONIDENTITY_AT_FIXED_FREQUENCY_TANGENTS' for item in columns)
                    nonapplicable_count=sum(item['applicable'] is False for item in columns)
                    unresolved_count=sum(item['status']=='UNRESOLVED' for item in columns)
                    control_status=('UNRESOLVED' if unresolved_count or responsive_count!=applicable_count else
                        'RESPONSIVE_LOADED_COLUMNS' if applicable_count else 'NOT_APPLICABLE_STRUCTURALLY_UNLOADED_COLUMNS')
                    loaded_controls.append(dict(drive=name,status=control_status,
                        applicableCount=applicable_count,responsiveCount=responsive_count,
                        nonapplicableCount=nonapplicable_count,unresolvedCount=unresolved_count,
                        originalRow=row[name],baseline=contractions[name],columns=columns))
                slot_diagnostics = []
                for name in DRIVES:
                    for slot,field_name in ((3,'Theta'),(4,'E')):
                        original = row[name]; probe = sp.eye(5)[:,slot]
                        route_control=face_omission(bound,index,leg,name,slot)
                        changed=route_control['recomputedRow']
                        before,after=original*probe,changed*probe
                        raw_difference = after-before; reduced = clean(raw_difference)
                        coefficient = sp.cancel(original[0,slot])
                        rational = coefficient.is_rational_function(k,kb,q,qb) is True
                        status = ('NOT_APPLICABLE_STRUCTURAL_ZERO' if coefficient==0 else
                                  'GENERIC_COEFFICIENT_PRESENT' if rational else 'UNRESOLVED')
                        slot_diagnostics.append(dict(drive=name,field=field_name,slot=slot,status=status,
                            scope='NATIVE_SLOT_PRESENCE_ONLY_NOT_T_CANCELLATION_OR_POWER',
                            sourceRoute=route_control,originalRow=original,changedRow=changed,probe=probe,before=before,after=after,
                            rawDifference=raw_difference,reducedDifference=reduced,coefficient=coefficient))
                ew = next(item for item in slot_diagnostics
                          if item['drive']=='OUTWARD_VELOCITY' and item['field']=='E')
                evidence('face-%d-leg-%d-e-velocity-omission'%(index,leg),ew)
                ew_omission_responsive = not zero(ew['reducedDifference'])
                scaled_velocity_coefficient = sp.cancel(ew['coefficient']/w)
                ew_calibrated = (not scaled_velocity_coefficient.free_symbols
                                 and scaled_velocity_coefficient!=0 and w.is_positive is True)
                ew = dict(ew,frequencyNormalizedCoefficient=scaled_velocity_coefficient,
                          frequencyDomain=sp.Gt(w,0,evaluate=False),presenceCalibrationSupported=ew_calibrated,
                          omissionResponsive=ew_omission_responsive)
                probe_vectors = [probe if leg==0 else sp.conjugate(probe).xreplace({k:kb})
                                 for probe in bound['longitudinalDirectionProbes']]
                probes = [dict(kind='SAVED_CURL_GAUGE_DIRECTION_EMBEDDED_AS_U_PROBE',vector=probe,
                    normalization='Arbitrary coordinate probe; not a normalized physical longitudinal mode',
                    values={name:row[name]*probe for name in DRIVES}) for probe in probe_vectors]
                # Index values are orientation-blind. The two harmonic signs
                # come from their actual time characters, not face orientation.
                original_velocity = row['OUTWARD_VELOCITY']
                joins = [sp.cancel(original_velocity[0,4]-item['coefficientsByHarmonicLeg'][leg]) for item in literals]
                control_counts={key:sum(control[key] for control in loaded_controls)
                                for key in ('applicableCount','responsiveCount','nonapplicableCount','unresolvedCount')}
                evidence('face-%d-leg-%d-operands'%(index,leg),dict(nativeRows=row,lift=lift,
                    contractions=contractions,slotDiagnostics=slot_diagnostics,c2Residuals=joins,
                    nontransverseProbes=probes,loadedControls=loaded_controls,uControlCounts=control_counts))
                records.append(dict(faceOrdinal=index,nativeSourceOrientation=face_order[index],harmonicLeg=leg,
                    harmonicMeaning='POSITIVE_CHARACTER' if leg==0 else 'NEGATIVE_CONJUGATE_CHARACTER',
                    amplitudeColumns=bound['harmonicAmplitudeColumns'][leg],nativeRows=row,lift=lift,
                    contractions=contractions,driveZero=zero(contractions),loadedUTermControls=loaded_controls,
                    uControlCounts=control_counts,
                    nontransverseSourceProbes=probes,nativeThetaESlotDiagnostics=slot_diagnostics,
                    eWPresenceCalibration=ew,c2IndexLiteralJoins=joins,
                    c2IndexLiteralSupported=(len(literals)==2) if bound['end']=='REFERENCE' else None,
                    c2OrientationIndependentlyChecked=False))
        forms = {name:clean(plus_lift.xreplace({k:kb}).H*matrix*plus_lift) for name,matrix in bound['currents'].items()
                 if name!='SLAB_CURRENT_MATRIX'}
        face_census = face_address_census(records)
        complete_faces = face_census['complete']
        evidence('face-forms-and-census',dict(forms=forms,census=face_census))
        checks = dict(faceLegCensus=complete_faces,
            allHarmonicDrivesZero=complete_faces and all(item['driveZero'] for item in records),
            nativeEVelocityPresence=complete_faces and all(item['eWPresenceCalibration']['presenceCalibrationSupported'] for item in records),
            nativeEVelocityOmissionResponse=complete_faces and all(item['eWPresenceCalibration']['omissionResponsive'] for item in records),
            loadedUControls=complete_faces and all(control['unresolvedCount']==0
                and control['responsiveCount']==control['applicableCount']
                for item in records for control in item['loadedUTermControls']),
            c2IndexLiteral=complete_faces and (bound['end']!='REFERENCE' or all(item['c2IndexLiteralSupported'] and
                zero(item['c2IndexLiteralJoins']) for item in records)),lossSideFormsZero=zero(forms))
        return dict(end=bound['end'],faces=records,faceLegCensus=face_census,physicalForms=forms,checks=checks,
            harmonicLegScope='BOTH_SAVED_REAL_FIELD_LEGS_NOT_TWO_ACOUSTIC_INCIDENCE_WAVES',
            c2LiteralEvidence=literals,c2EvidenceScope='ORIENTATION_BLIND_INDEX_LITERAL_ONLY',
            flatExteriorSourceInternalJoins=bound['flatExteriorJoins'],
            c2PressureTraceFirstShapeJoin='NOT_SUPPORTED_BY_THIS_INPUT_INDEX')

    def numeric(matrix, mapping):
        return np.asarray(matrix.subs(mapping).evalf(40).tolist(), dtype=complex)

    def numerical_face_legs(bound,mapping,basis):
        records = []
        for face_index_value,face_rows in enumerate(bound['drives']):
            for leg,row in enumerate(face_rows):
                loaded_basis = basis if leg==0 else basis.conj()
                native_rows = {name:numeric(row[name],mapping) for name in DRIVES}
                contractions = {name:value@loaded_basis for name,value in native_rows.items()}
                e_coefficient = native_rows['OUTWARD_VELOCITY'][0,4]
                records.append(dict(faceOrdinal=face_index_value,harmonicLeg=leg,
                    nativeRows=native_rows,loadedBasis=loaded_basis,contractions=contractions,
                    nativeEVelocityCoefficient=e_coefficient,
                    nativeEVelocityPresent=abs(e_coefficient)>1e-12))
        return records

    def loaded_control(bound,source_kind,mapping,basis,scale):
        raw=bound['originalCurrentResult']['CLOSED_PENCIL_LEGS'][0]
        wl,wr=bound['originalCurrentResult']['FREQUENCY_LEGS']
        symbols={wl,wr} if source_kind=='frequency' else {
            symbol for symbol in raw.free_symbols if symbol.name in
            ('s11cdTangentialMomentum1','s11cdTangentialMomentum2')}
        baseline_matrix=numeric(bound['matrix'],mapping)
        original=(baseline_matrix/scale[:,None])@basis
        records=[]
        for column in range(basis.shape[1]):
            selected,best=None,0.
            for i in range(raw.rows):
                for j in range(raw.cols):
                    if abs(basis[j,column])<=1e-12 or not raw[i,j].has(*symbols): continue
                    for term in sp.Add.make_args(sp.expand(raw[i,j])):
                        if not term.has(*symbols): continue
                        fixed_term=bind(term,bound['end'],bound['sourceLegMapping'])
                        movement=abs(complex(fixed_term.subs(mapping).evalf(40))*basis[j,column])/scale[i]
                        if movement>best: best,selected=movement,(i,j,term,fixed_term)
            if selected is None:
                records.append(dict(column=column,status='UNRESOLVED',reason='NO_LOADED_NATIVE_SOURCE_TERM'))
                continue
            i,j,term,fixed_term=selected
            changed=sp.MutableDenseMatrix(raw);changed[i,j]-=term
            # Adapted live-parameter control: omit an addressed frequency- or
            # tangent-dependent native pencil term BEFORE fixed specialization.
            # Rebind the changed source; never replace the resulting matrix.
            altered_bound=bind(sp.ImmutableMatrix(changed),bound['end'],bound['sourceLegMapping'])
            altered_matrix=numeric(altered_bound,mapping)
            altered=(altered_matrix/scale[:,None])@basis
            difference=altered-original
            column_norm=float(np.linalg.norm(difference[:,column]))
            record=dict(column=column,sourceAddress=('CLOSED_PENCIL_LEGS',0,i,j),
                sourceEntry=raw[i,j],sourceSymbols=tuple(sorted(symbols,key=str)),
                omittedNativeTerm=term,specializedOmittedTerm=fixed_term,
                changedSourceEntry=changed[i,j],reboundMatrix=altered_bound,
                mutatedMatrix=altered_matrix,mutatedResidual=altered,difference=difference,
                loadedColumnMovement=column_norm)
            evidence('pencil-'+source_kind+'-column-'+str(column),record)
            records.append(dict(record,status='RESPONSIVE' if column_norm>1e-8 else 'UNRESOLVED'))
        return dict(status='RESPONSIVE' if records and all(item['status']=='RESPONSIVE' for item in records) else 'UNRESOLVED',
            scope='SOURCE_INPUT_OMISSION_AT_FIXED_POINT_NOT_LIVE_PARAMETER_COVERAGE',sourceKind=source_kind,
            baselineMatrix=baseline_matrix,fixedLoadedBasis=basis,baselineResidual=original,columns=records)

    def saved_seed(bound, old):
        info = old['info']
        k0, q0 = complex(info['K']), complex(info['PHYSICAL_Q'])
        w0 = sp.Rational(str(info['OMEGA'].real)) if isinstance(info['OMEGA'],complex) else sp.Rational(str(info['OMEGA']))
        v0 = params['s11cdTangentialMomentum2']
        mapping={k:k0,kb:k0.conjugate(),q:q0,qb:q0.conjugate()}
        matrix = numeric(bound['matrix'],mapping)
        lift = numeric(bound['lift'],mapping)
        # Projection through the 2-column geometric lift, not a new modal solve.
        gram = lift.conj().T@lift
        dual = np.linalg.inv(gram)@lift.conj().T
        right = np.asarray(old['right'],complex)
        t_residual = right-lift@dual@right
        scale = np.maximum(np.linalg.norm(matrix,axis=1),1.)
        pencil_join = (matrix-np.asarray(old['pencil'],complex))/scale[:,None]
        kernel = (matrix/scale[:,None])@right
        face = numerical_face_legs(bound,mapping,right)
        face_census = face_address_census(face)
        slab = numeric(bound['currents']['SLAB_CURRENT_MATRIX'],mapping)
        bulk = numeric(bound['currents']['BULK_NORMAL_CURRENT_DENSITY_MATRIX'],mapping)
        current = right.conj().T@slab@right
        projected_bulk = right.conj().T@bulk@right
        operand_joins = {'slab':slab-np.asarray(old['currentOperands']['CURRENT_SLAB'],complex),
                         'bulk':bulk-np.asarray(old['currentOperands']['CURRENT_BULK'],complex)}
        evidence('saved-seed-pencil-current-operands',dict(info=info,actualPointMapping=mapping,
            matrix=matrix,lift=lift,gram=gram,dual=dual,right=right,
            transverseResidual=t_residual,pencilResidual=pencil_join,kernelResidual=kernel,
            slab=slab,bulk=bulk,current=current,projectedBulk=projected_bulk,
            currentOperandResiduals=operand_joins,face=face,faceLegCensus=face_census))
        evidence('saved-seed-domain-operands',dict(sourceDenominators=bound['sourceDenominators'],
                 chart=bound['chart'],actualPointMapping=mapping))
        domain_values=[]
        for factor in (*bound['sourceDenominators'],bound['chart']):
            actual=factor.subs(mapping)
            remaining=tuple(sorted(actual.free_symbols,key=str))
            if remaining:
                domain_values.append(dict(sourceFactor=factor,pointValue=actual,remainingSymbols=remaining,
                     numeric=None,finite=False,magnitude=None))
                continue
            value=complex(actual.evalf(40))
            domain_values.append(dict(sourceFactor=factor,pointValue=actual,remainingSymbols=remaining,numeric=value,
                 finite=bool(np.isfinite(value)),magnitude=float(abs(value)) if np.isfinite(value) else None))
        evidence('saved-seed-source-domains',domain_values)
        domain_supported=all(item['finite'] and item['magnitude']>1e-12 for item in domain_values)
        controls=[loaded_control(bound,kind,mapping,right,scale) for kind in ('frequency','tangents')]
        wave = complex(bound['wave'].subs(mapping).evalf(40))
        normal_reality_residual = k0-k0.conjugate()
        normal_is_real = k0.imag == 0.0
        face_zero = face_census['complete'] and all(np.linalg.norm(value)<1e-8 for item in face
                        for value in item['contractions'].values())
        bulk_zero = np.linalg.norm(projected_bulk)<1e-8
        hermitian_residual = current-current.conj().T
        hermitian_part = (current+current.conj().T)/2
        current_eigenvalues = np.linalg.eigvalsh(hermitian_part)
        hermitian_supported = np.linalg.norm(hermitian_residual)<1e-8*max(1.,np.linalg.norm(current))
        current_rank_supported = bool(current_eigenvalues.size and np.min(np.abs(current_eigenvalues))>1e-9)
        transport = dict(scope='RECOMPUTED_SLAB_FORM_ONLY_ON_UNDRIVEN_T_ZERO_PROJECTED_BULK_DOMAIN',
            domainSupported=bool(face_zero and bulk_zero),right=right,slabMatrix=slab,bulkMatrix=bulk,
            currentForm=current,projectedBulk=projected_bulk,hermitianResidual=hermitian_residual,
            hermitianPart=hermitian_part,eigenvaluesOfHermitianPart=current_eigenvalues,
            hermitianTolerance=1e-8*max(1.,float(np.linalg.norm(current))),rankTolerance=1e-9,
            positiveCurrentRank=int(sum(current_eigenvalues>1e-9)),
            negativeCurrentRank=int(sum(current_eigenvalues < -1e-9)),
            interpretation='Tolerance-based current check on the saved basis; no new normalized basis')
        history = dict(currentDefined=old.get('currentDefined'),
            sourceSheetMembership=info.get('SHEET_MEMBERSHIP'),sourceExactRealNormal=info.get('EXACT_REAL_NORMAL'),
            sourceBulkDecayDiskCertified=info.get('BULK_DECAY_DISK_CERTIFIED'),
            sourcePhysicalCurrentNormalization=info.get('PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED'),
            scope='HISTORICAL_METADATA_NOT_RECOMPUTED_AVAILABILITY_GATE')
        missing = [name for name in ('depthIntegral','currentGram') if name not in old]
        if missing:
            historical_comparison = dict(status='UNAVAILABLE',missingFields=missing,
                availableSavedDepthIntegral=old.get('depthIntegral'),availableSavedGram=old.get('currentGram'),
                availabilityGated=False)
        else:
            depth = complex(old['depthIntegral'])
            weighted_matrix = slab+depth*bulk
            weighted_gram = right.conj().T@weighted_matrix@right
            saved_gram = np.asarray(old['currentGram'],complex)
            weighted_residual = weighted_gram-saved_gram
            comparison_supported = np.linalg.norm(weighted_residual)<1e-8
            historical_comparison = dict(status='COMPARED',right=right,slabMatrix=slab,bulkMatrix=bulk,
                savedDepthIntegral=depth,weightedMatrix=weighted_matrix,recomputedWeightedGram=weighted_gram,
                savedGram=saved_gram,residual=weighted_residual,residualNorm=float(np.linalg.norm(weighted_residual)),
                tolerance=1e-8,supported=bool(comparison_supported),availabilityGated=True,
                source='S11c_d_uniform_response.py:74-88; no depth integral is re-evaluated')
        depth_evidence = []
        for candidate in (q0,-q0):
            dm = dict(mapping); dm.update({q:candidate,qb:candidate.conjugate()})
            native_flux = complex(bound['outgoingDepthCoefficient'].subs(dm).evalf(40))
            carrier_flux = complex(bound['carrierLawCheck']['normalizedFlux'].subs(dm).evalf(40))
            normalization_residual = native_flux-carrier_flux
            decay = bool(abs(candidate.real)<1e-12 and candidate.imag>0)
            real_outward = bool(abs(candidate.imag)<1e-12 and abs(native_flux.imag)<1e-10 and native_flux.real>0)
            depth_evidence.append(dict(candidate=candidate,nativeFlux=native_flux,carrierFlux=carrier_flux,
                normalizationResidual=normalization_residual,decayingImaginaryDepth=decay,
                realOutwardDepth=real_outward,coordinateRealityTolerance=1e-12,
                historicalDecayDisk=info.get('BULK_DECAY_DISK_CERTIFIED'),
                supported=(abs(normalization_residual)<1e-10 and (decay or real_outward))))
        evidence('saved-seed-residuals',dict(sourceMatrix=matrix,savedPencil=old['pencil'],
            sourceLift=lift,savedRight=right,transverseResidual=t_residual,pencilResidual=pencil_join,
            kernelResidual=kernel,waveResidual=wave,currentOperandResiduals=operand_joins,
            current=current,projectedBulk=projected_bulk,hermitianResidual=hermitian_residual,
            normalRealityResidual=normal_reality_residual,face=face,faceCensus=face_census,
            physicalDepthEvidence=depth_evidence,historicalWeightedGram=historical_comparison,
            controls=controls,sourceDomain=domain_values))
        checks = dict(sourcePointDomain=domain_supported,faceLegCensus=face_census['complete'],savedFrequency=(w0==w==params['omega']),
            transverseMembership=np.linalg.norm(t_residual)<1e-8,
            physicalPencilJoin=np.linalg.norm(pencil_join)<1e-8,
            physicalKernel=np.linalg.norm(kernel)<1e-8, physicalWave=abs(wave)<1e-8,
            recomputedRealNormal=normal_is_real,
            actualPhysicalDepth=(depth_evidence[0]['supported'] and not depth_evidence[1]['supported']),
            slabCurrentOperands=np.linalg.norm(operand_joins['slab'])<1e-8,
            bulkCurrentOperands=np.linalg.norm(operand_joins['bulk'])<1e-8,
            nativeEVelocityPresence=face_census['complete'] and all(item['nativeEVelocityPresent'] for item in face),
            undriven=face_zero,projectedBulk=bulk_zero,
            transverseCurrentHermitian=hermitian_supported,transverseCurrentNonzeroRank=current_rank_supported,
            controls=all(item['status']=='RESPONSIVE' for item in controls))
        if historical_comparison['status']=='COMPARED':
            checks['historicalWeightedGramComparison'] = historical_comparison['supported']
        return dict(status='AVAILABLE' if all(checks.values()) else 'UNRESOLVED',
            scope='SELECTED_SAVED_SEED_SUBSPACE_JOIN_NOT_COMPLETE_CENSUS',info=info,
            historicalEligibility=history,historicalWeightedGramComparison=historical_comparison,
            checks=checks,sourceMatrix=matrix,savedPencil=old['pencil'],sourceLift=lift,
            savedRight=right,transverseResidual=t_residual,pencilResidual=pencil_join,
            kernelResidual=kernel,waveResidual=wave,face=face,faceLegCensus=face_census,recomputedTransverseTransport=transport,
            actualNormalMomentum=k0,normalRealityResidual=normal_reality_residual,
            normalRealityCriterion='Zero imaginary part of the saved numeric coefficient, checked directly',
            currentOperandResiduals=operand_joins,projectedBulk=projected_bulk,controls=controls,
            sourcePointDomainValues=domain_values,
            savedDepthSignEvidence=depth_evidence)

    packets = {}
    for label in ENDS:
        packets[label] = (restored(label+'-restore-original-symbol',label+'Frequency'),
                          restored(label+'-restore-current',label+'Pairing'))
    require(journal.reused==journal.prior_order and len(journal.reused)==11,'all original completed returns reused')
    save(out/'complete-return-reuse.json',dict(names=journal.reused,number=len(journal.reused),
         summaries=journal.summary_reuse))
    end_states, face_results, summaries = {}, {}, []
    seed_results={label:[dict(seedIndex=index,status='UNRESOLVED',reason='NOT_VISITED',attempted=False)
                         for index in (16,17)] for label in ENDS}
    for label in ENDS:
        if time.monotonic()>journal.deadline-60:
            end_states[label]='UNRESOLVED_BUDGET_NOT_ATTEMPTED'
            continue
        packet,pairing = packets[label]
        bound = op(label+'-source-bind',bound_end,label,packet,pairing,seconds=180,render=False)
        if 'checks' not in bound or not all(bound['checks'].values()):
            end_states[label]='UNRESOLVED_SOURCE_JOIN'
            continue
        face = op(label+'-face-premise',faces,bound,seconds=20)
        face_results[label]=face
        if 'checks' not in face or not all(face['checks'].values()):
            end_states[label]='UNRESOLVED_FACE_CHECK'
            continue
        for ordinal,index in enumerate((16,17)):
            if time.monotonic()>journal.deadline-30: break
            key=label+'Seed'+str(index)
            saved=restored(label+'-restore-seed-'+str(index),key)
            joined=op(label+'-seed-join-'+str(index),saved_seed,bound,saved,seconds=20)
            seed_results[label][ordinal]=dict(joined,seedIndex=index,attempted=True)
        statuses=[item.get('status','UNRESOLVED') for item in seed_results[label]]
        selected_premise_supported=len(statuses)==2 and all(
            item.get('status')=='AVAILABLE' and bool(item.get('checks')) and all(item['checks'].values())
            for item in seed_results[label])
        attempted_count=sum(item['attempted'] for item in seed_results[label])
        state='SELECTED_PREMISE_SUPPORTED' if selected_premise_supported else 'UNRESOLVED_SELECTED_PREMISE'
        end_states[label]=state
        summaries.append(dict(end=label,status=state,omega=str(w),
             tangents=[str(params['s11cdTangentialMomentum1']),str(params['s11cdTangentialMomentum2'])],
             scope='SELECTED_UNDRIVEN_TRANSVERSE_CURRENT_PREMISE_NOT_CHANNEL_CENSUS',
             coverage='SELECTED_SAVED_SUBSPACES_ONLY' if selected_premise_supported else 'PARTIAL_UNRESOLVED_SELECTED_SUBSPACES',
             completeInventory=False,selectedSeedIndices=[16,17],
             checkedSelectedSeeds=attempted_count,selectedSeedStatuses=statuses,
             selectedChecksComplete=selected_premise_supported))
    save(out/'source-end-states.json',end_states)
    save(out/'reference-face-drive.json',show(face_results.get('REFERENCE',{'status':'UNRESOLVED'})))
    seed_hashes={}
    for label in ENDS:
        for index in (16,17):
            key=label+'Seed'+str(index)
            seed_hashes.setdefault(spec['inputs'][key]['sha256'],[]).append(key)
    save(out/'saved-seed-comparison.json',show(dict(selectedRecords=seed_results,
         savedSeedProvenance=spec['savedSeedProvenance'],sharedAcceptedHashes=seed_hashes)))
    unvisited=[dict(end=label,status='UNRESOLVED',coverage='UNRESOLVED',reason=end_states[label])
               for label in ENDS if label not in {item['end'] for item in summaries}]
    save(out/'point-summary.json',summaries)
    save(out/'unvisited-point.json',unvisited)
    reference_face=face_results.get('REFERENCE',{})
    return dict(case=spec['case'],endStates=end_states,pointRecords=len(summaries),
                reusedCompleteCount=len(journal.reused),endBindingsAttempted=len(journal.resumed),
                selectedSeedChecks=sum(item['attempted'] for values in seed_results.values() for item in values),
                selectedSeedStatuses={label:[item['status'] for item in values] for label,values in seed_results.items()},
                referenceFaceAttempted='REFERENCE' in face_results,
                referenceFaceSupported=bool(reference_face.get('checks')) and all(reference_face['checks'].values()),
                pointStates={label:end_states[label] for label in ENDS})


def main():
    parser=argparse.ArgumentParser(__doc__)
    for name in ('input-manifest','gate-receipt','run-directory'):
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    spec=json.loads(args.input_manifest.read_text())
    gate=json.loads(args.gate_receipt.read_text())
    require(gate.get('status')=='READY_FOR_ONE_PROGRESS_GUARDED_FIXED_POINT_CONTINUATION','fresh exact execution-4 gate required')
    pins=dict(workerSha256=sha(__file__),inputManifestSha256=sha(args.input_manifest),
              scopeSha256=spec['inputs']['continuationScope']['sha256'],
              authorizationSha256=spec['inputs']['continuationAuthorization']['sha256'])
    require(all(gate.get(key)==value for key,value in pins.items()),'worker/manifest/scope/authorization gate pins')
    require(gate.get('independentBuildClearance') is True and gate.get('scienceJobOrdinal')==5
            and gate.get('maximumNewScientificExecutions')==1,'fresh review and one continuation only')
    require(gate.get('seconds')==0 and gate.get('nativeSeconds')==0
            and gate.get('progressStallSeconds')==3600 and gate.get('maximumEndBindings')==3
            and gate.get('automaticRetry') is False,'progress-dependent duration and no retry')
    require(gate.get('guardSha256')==sha(ROOT/'research/pde_ledger_v3/_measurements/S11c_d_transverse_face_fixed_point_continue_guard.py')
            and gate.get('sharedGuardSha256')==sha(ROOT/'scripts/s11c_guarded_run.py')
            and gate.get('supervisorSha256')==sha(ROOT/'research/pde_ledger_v3/_measurements/S11c_d_end_normalization_run.py'),
            'shared guard/supervisor pins')
    review=gate['reviewRecord']
    require(route(review['path'])==review,'review record pin')
    reviewed=json.loads(Path(review['path']).read_text())
    require(reviewed.get('independentBuildClearance') is True
            and all(reviewed.get(key)==value for key,value in pins.items()),'exact reviewed packet')
    legs=reviewed.get('reviews',[])
    require({leg.get('engine') for leg in legs}=={'claude','grok'} and len(legs)==2,'fresh independent Claude/Grok legs')
    for leg in legs:
        require(leg.get('verdict')=='CLEAR' and route(leg['output']['path'])==leg['output']
                and leg['output']['bytes']>0,'independent review output pin and recorded clearance')
    for name,record in spec['inputs'].items():
        require(route(record['path'])==record,'changed input '+name)
    approval=json.loads(Path(spec['inputs']['continuationAuthorization']['path']).read_text())
    require(approval.get('scientificExecutionOrdinal')==5 and approval.get('maximumNewScientificExecutions')==1
            and approval.get('progressStallSeconds')==3600 and approval.get('stopAfterResults') is True,
            'explicit longer continuation authorization scope')
    require(spec['case']=='LAB_HELD_RHO4_CONSTANT' and spec['fixedPoint']==dict(omega='1',tangents=['1/5','1/10']),
            'exact authorized one-point context')
    require(json.loads(Path(spec['inputs']['physicalInputFile']['path']).read_text())==spec['physicalInput'],
            'physical inputs unchanged')
    require({key:spec['inputs'][key] for key in spec['originalInputKeys']}==spec['originalSpecification']['inputs'],
            'original 33 source pins unchanged')
    for join in spec['checkpointJoins']:
        checkpoint=json.loads(Path(join['checkpoint']).read_text())
        require(checkpoint['artifacts'][join['artifact']]['sha256']==spec['inputs'][join['inputKey']]['sha256'],
                'accepted source checkpoint')
    out=args.run_directory.resolve()
    require(out==Path(spec['runRoot'])/'complete' and gate.get('runRoot')==spec['runRoot'],
            'one reserved run directory; no fallback')
    started=time.monotonic();limits=containment()
    out.mkdir(parents=True,exist_ok=False)
    journal=Journal(out,started+3600,spec)
    save(out/'native-limits.json',limits);save(out/'manifest.json',spec);save(out/'gate.json',gate)
    try:
        journal.configure_prior()
        journal.configure_resume()
        result=construct(spec,out,journal)
        post={name:route(record['path']) for name,record in spec['inputs'].items()}
        save(out/'posthashes.json',post);require(post==spec['inputs'],'posthash mismatch')
        result.update(wallSeconds=time.monotonic()-started,operations=len(journal.records),
                      allSourcesUnchanged=post==spec['inputs'])
        save(out/'operation-index.json',journal.records);save(out/'checks.json',result)
        # Preserve the established contract: one JSON document on stdout,
        # byte-identical to checks.json. Raw computed evidence has already
        # been emitted to immutable operation/suboperation artifacts.
        print(json.dumps(result,indent=2,allow_nan=False),flush=True)
    except BaseException:
        signal.setitimer(signal.ITIMER_REAL,0)
        save(out/'failure.json',dict(traceback=traceback.format_exc(),incompleteOperation=journal.active,
             wallSeconds=time.monotonic()-started,noRetry=True))
        save(out/'operation-index.json',journal.records)
        save(out/'failure-posthashes.json',{name:route(record['path']) for name,record in spec['inputs'].items()})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)
        save(out/'restored-suboperations.json',dict(restored=sorted(journal.restored_suboperations),
             expected=sorted(journal.resume_suboperations),count=len(journal.restored_suboperations)))
        save(out/'artifact-index.json',{str(path.relative_to(out)):route(path)
             for path in sorted(out.rglob('*')) if path.is_file()})


if __name__=='__main__':
    main()
