#!/usr/bin/env python3
"""Job-2 continuation: restore complete returns; finish the saved reader.

No producer imports or calls, binding, grading, derivatives, roots, modes,
solves, transforms, integrals, domain certification, or frequency substitution.
Scientific interpretation belongs to a separate report after this reader.
"""
import argparse
import builtins
from collections import OrderedDict
from datetime import datetime, timezone
import hashlib
import importlib
import io
import json
import os
from pathlib import Path
import pickle
import resource
import signal
import sys
import tempfile
import time
import traceback

ROOT = Path('/var/projects/toy_physics')
STORE = ROOT / '_scratch/s11c'
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')
ENDS = ('REFERENCE', 'LEFT', 'RIGHT')


def require(condition, reason):
    if not condition:
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


def publish(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.' + path.name + '.', dir=path.parent)
    with os.fdopen(fd, 'wb') as stream:
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())
    os.link(temporary, path)
    directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)
    os.unlink(temporary)


def encoded(value):
    return (json.dumps(value, indent=2, allow_nan=False) + '\n').encode()


def save(path, value):
    publish(path, encoded(value))


class NativeDeadline(BaseException):
    pass


def containment():
    group = next(line[3:] for line in Path('/proc/self/cgroup').read_text().splitlines()
                 if line.startswith('0::'))
    cgroup = Path('/sys/fs/cgroup') / group.lstrip('/')
    observed = {name: (cgroup / name).read_text().strip()
                for name in ('memory.max', 'memory.swap.max', 'pids.max')}
    observed.update(cgroup=str(cgroup), nice=os.getpriority(os.PRIO_PROCESS, 0),
                    affinity=sorted(os.sched_getaffinity(0)),
                    threads={name: os.environ.get(name) for name in THREADS})
    require(observed['memory.max'] == '2147483648'
            and observed['memory.swap.max'] == '0' and observed['pids.max'] == '32'
            and observed['nice'] >= 15 and len(observed['affinity']) == 1
            and all(value == '1' for value in observed['threads'].values()),
            'ordinary containment absent; no saved scientific objects restored')
    resource.setrlimit(resource.RLIMIT_AS, (2147483648, 2147483648))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    observed.update(rlimitAddressSpace=list(resource.getrlimit(resource.RLIMIT_AS)),
                    nativeSeconds=840)
    def timeout(*_):
        raise NativeDeadline('840-second native whole-reader deadline; no retry')
    signal.signal(signal.SIGALRM, timeout)
    signal.setitimer(signal.ITIMER_REAL, 840)
    return observed


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


def metadata_identity(left, right):
    """Object identity for restored science; exact equality for plain metadata.

    No re-pickle canonicalization or symbolic equivalence is used to decide
    whether a completed operation may be reused.
    """
    if left is right:
        return True
    scalar = (str, int, float, bool, type(None))
    if type(left) in scalar or type(right) in scalar:
        return type(left) is type(right) and left == right
    if type(left) is dict and type(right) is dict:
        if not all(type(key) in scalar for key in (*left, *right)):
            return False
        return (left.keys() == right.keys()
                and all(metadata_identity(left[key], right[key]) for key in left))
    if type(left) in (tuple, list) and type(left) is type(right):
        return len(left) == len(right) and all(metadata_identity(a, b) for a, b in zip(left, right))
    return False


class Journal:
    def __init__(self, out, spec):
        self.out, self.spec = out, spec
        self.records, self.active, self.objects = [], None, {}
        self.prior, self.prior_order, self.prior_cache = {}, [], {}
        self.reused, self.summary_reuse = [], []
        self.incomplete, self.resumed_incomplete = None, False

    def prior_json(self, relative):
        record = self.spec['inputs'][self.spec['priorFiles'][relative]]
        require(route(record['path']) == record, 'prior metadata changed: ' + relative)
        return json.loads(Path(record['path']).read_text())

    def configure_prior(self):
        root = Path(self.spec['priorRoot'])
        index = self.prior_json('operation-index.json')
        failure = self.prior_json('failure.json')
        checkpoint = json.loads(Path(self.spec['inputs']['priorFailureCheckpoint']['path']).read_text())
        observed = {name: self.spec['inputs'][key] for name, key in self.spec['priorFiles'].items()}
        names = {str(path.relative_to(root)) for path in root.rglob('*') if path.is_file()}
        checks = dict(artifactCensus=names == set(observed),
                      checkpointArtifacts=observed == checkpoint['artifacts'],
                      artifactIndex={name: item for name, item in observed.items()
                                     if name != 'artifact-index.json'} == self.prior_json('artifact-index.json'),
                      recordedFailure=failure == checkpoint['failure'],
                      completeCount=len(index) == self.spec['priorCompleteCount'],
                      originalManifest=self.prior_json('manifest.json') == self.spec['originalSpecification'])
        for ordinal, record in enumerate(index):
            folder = 'operations/%04d-%s/' % (ordinal, record['name'])
            prior_input = self.prior_json(folder + 'input-references.json')
            checks['receipt-%02d' % ordinal] = (
                record['status'] == 'COMPLETE' and record == self.prior_json(folder + 'complete.json')
                and prior_input == {key: record[key] for key in ('name', 'startedUtc', 'operands')})
        self.incomplete = failure['incompleteOperation']
        unfinished_folder = 'operations/%04d-%s/' % (len(index), self.incomplete['name'])
        checks['unfinishedReceipt'] = self.incomplete == self.prior_json(unfinished_folder + 'input-references.json')
        checks['unfinishedHasNoReturn'] = 'value' not in self.incomplete
        save(self.out / 'prior-continuation-metadata.json', dict(checks=checks, failure=failure,
             completeOperations=index, incompleteOperation=self.incomplete))
        require(all(checks.values()), 'prior continuation metadata mismatch; evidence saved')
        self.prior = {record['name']: record for record in index}
        self.prior_order = [record['name'] for record in index]
        require(len(self.prior) == len(index), 'duplicate prior operation names')
        # Preserve each exact prior journal operand/return byte string, including
        # the two incomplete operands, before any continuation restoration.
        for record in [*index, self.incomplete]:
            for item in [*record['operands'], *([record['value']] if 'value' in record else [])]:
                self.copy_prior_object(item)

    def copy_prior_object(self, record):
        allowed = [self.spec['inputs'][key] for key in self.spec['priorRestoreKeys']]
        require(record in allowed, 'unselected prior result object')
        require(route(record['path']) == record, 'prior result object changed')
        payload = Path(record['path']).read_bytes()
        require(hashlib.sha256(payload).hexdigest() == record['sha256'], 'prior bytes changed during read')
        path = self.out / 'prior-objects' / (record['sha256'] + '.pickle')
        if not path.exists():
            publish(path, payload)
        copied = route(path)
        require(copied['sha256'] == record['sha256'] and copied['bytes'] == record['bytes'],
                'prior object copy mismatch')
        return copied

    def prior_value(self, record):
        copied = self.copy_prior_object(record)
        digest = record['sha256']
        if digest not in self.prior_cache:
            stream = io.BytesIO(Path(copied['path']).read_bytes())
            value = SavedCodec(stream).load()
            require(not stream.read(1), 'trailing bytes after prior saved return')
            self.prior_cache[digest] = value
            self.objects[id(value)] = (value, copied)
        return self.prior_cache[digest]

    def has_prior(self, name):
        return name in self.prior

    def copy_summary(self, relative):
        require(relative in self.spec['priorSummaryFiles'], 'unselected complete summary')
        record = self.spec['inputs'][self.spec['priorFiles'][relative]]
        require(route(record['path']) == record, 'complete summary changed')
        payload = Path(record['path']).read_bytes()
        publish(self.out / relative, payload)
        copied = route(self.out / relative)
        require(copied['sha256'] == record['sha256'] and copied['bytes'] == record['bytes'],
                'complete summary copy mismatch')
        self.summary_reuse.append(dict(source=record, copy=copied,
                                       execution='COPIED_PRIOR_COMPLETE_SUMMARY'))
        return json.loads(payload)

    def blob(self, value):
        old = self.objects.get(id(value))
        if old is not None and old[0] is value:
            require(route(old[1]['path']) == old[1], 'journal object changed')
            return old[1]
        payload = pickle.dumps(value, protocol=4)
        digest = hashlib.sha256(payload).hexdigest()
        path = self.out / 'objects' / (digest + '.pickle')
        if not path.exists():
            publish(path, payload)
        record = route(path)
        require(record['sha256'] == digest and record['bytes'] == len(payload),
                'content-addressed publication mismatch')
        self.objects[id(value)] = (value, record)
        return record

    def op(self, name, function, *operands):
        directory = self.out / 'operations' / ('%04d-%s' % (len(self.records), name))
        self.active = dict(name=name, startedUtc=datetime.now(timezone.utc).isoformat(), operands=[])
        save(directory / 'started.json', self.active)
        for operand in operands:
            self.active['operands'].append(self.blob(operand))
        save(directory / 'input-references.json', self.active)
        if len(self.reused) < len(self.prior_order):
            require(name == self.prior_order[len(self.reused)], 'prior complete sequence changed')
            previous = self.prior[name]
            saved_operands = [self.prior_value(item) for item in previous['operands']]
            matches = [metadata_identity(a, b) for a, b in zip(operands, saved_operands)]
            evidence = dict(priorOperation=previous, suppliedOperands=self.active['operands'],
                            operandCount=len(operands) == len(saved_operands), matches=matches,
                            identityRule='RESTORED_OBJECT_IDENTITY_OR_PLAIN_METADATA_EQUALITY')
            save(directory / 'prior-argument-identity.json', evidence)
            require(evidence['operandCount'] and all(matches), 'prior completed arguments changed')
            value = self.prior_value(previous['value'])
            receipt = dict(self.active, status='COMPLETE', execution='RESTORED_PRIOR_COMPLETE_RETURN',
                           priorOperation=previous, argumentIdentity=evidence,
                           value=self.copy_prior_object(previous['value']))
            # Deliberately do not call function for any prior COMPLETE operation.
            self.reused.append(name)
        else:
            require(name not in self.prior, 'completed operation cannot be executed again')
            provenance = dict(execution='NEW_UNFINISHED_READER_OPERATION')
            if not self.resumed_incomplete:
                require(name == self.incomplete['name'], 'resume must begin at the recorded incomplete operation')
                saved_operands = [self.prior_value(item) for item in self.incomplete['operands']]
                matches = [metadata_identity(a, b) for a, b in zip(operands, saved_operands)]
                evidence = dict(priorIncomplete=self.incomplete, suppliedOperands=self.active['operands'],
                                operandCount=len(operands) == len(saved_operands), matches=matches,
                                identityRule='RESTORED_OBJECT_IDENTITY_OR_PLAIN_METADATA_EQUALITY')
                save(directory / 'prior-incomplete-identity.json', evidence)
                require(evidence['operandCount'] and all(matches), 'incomplete operands changed')
                provenance = dict(execution='RESUMED_PRIOR_INCOMPLETE_READER_OPERATION',
                                  priorIncomplete=self.incomplete, argumentIdentity=evidence)
            value = function(*operands)
            receipt = dict(self.active, status='COMPLETE', value=self.blob(value), **provenance)
            if name == self.incomplete['name']:
                self.resumed_incomplete = True
        save(directory / 'complete.json', receipt)
        self.records.append(receipt)
        self.active = None
        return value


def inspect_saved(spec, out, journal):
    import sympy as sp
    import numpy as np
    from sympy.core.function import UndefinedFunction

    inputs = spec['inputs']
    by_path = {record['path']: key for key, record in inputs.items()}
    cache = {}
    def read_json(key):
        record = inputs[key]
        require(route(record['path']) == record, 'JSON route changed: ' + key)
        return json.loads(Path(record['path']).read_text())

    def read_pickle(key):
        require(key in spec['restoreKeys'], 'unselected saved object: ' + key)
        if key not in cache:
            record = inputs[key]
            def restore(saved_route):
                require(route(saved_route['path']) == saved_route, 'saved route changed')
                payload = Path(saved_route['path']).read_bytes()
                require(hashlib.sha256(payload).hexdigest() == saved_route['sha256'],
                        'saved bytes changed during read')
                # Preserve the exact original bytes before deserialization.
                path = out / 'source-objects' / (saved_route['sha256'] + '.pickle')
                if not path.exists():
                    publish(path, payload)
                stream = io.BytesIO(payload)
                value = SavedCodec(stream).load()
                require(not stream.read(1), 'trailing bytes after saved object')
                return value
            cache[key] = journal.op('read-saved-%03d' % len(cache), restore, record)
        return cache[key]

    def same(left, right):
        if isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
            return (isinstance(left, np.ndarray) and isinstance(right, np.ndarray)
                    and left.shape == right.shape and left.dtype == right.dtype
                    and bool(np.array_equal(left, right)))
        if isinstance(left, dict) or isinstance(right, dict):
            return (isinstance(left, dict) and isinstance(right, dict)
                    and left.keys() == right.keys()
                    and all(same(left[key], right[key]) for key in left))
        if isinstance(left, (tuple, list)) or isinstance(right, (tuple, list)):
            return (type(left) is type(right) and len(left) == len(right)
                    and all(same(a, b) for a, b in zip(left, right)))
        return bool(left == right)

    def undefined_function_metadata(value):
        # Read the saved class namespace only; never apply/instantiate the class
        # or request its instance-only free_symbols/args properties.
        namespace = type.__getattribute__(value, '__dict__')
        return dict(type='SAVED_UNDEFINED_FUNCTION_CLASS', name=namespace['name'],
                    savedKeywordMetadata=display(namespace.get('_kwargs', {})))

    def display(value):
        if value is None or isinstance(value, (str, bool, int, float)):
            return value
        if isinstance(value, complex):
            return dict(real=value.real, imag=value.imag)
        if isinstance(value, np.generic):
            return display(value.item())
        if isinstance(value, dict):
            return {str(key): display(item) for key, item in value.items()}
        if isinstance(value, (tuple, list)):
            return [display(item) for item in value]
        if isinstance(value, (set, frozenset)):
            return dict(type=type(value).__name__,
                        entries=[display(item) for item in sorted(value, key=str)])
        if isinstance(value, slice):
            return dict(type='slice', start=display(value.start),
                        stop=display(value.stop), step=display(value.step))
        if isinstance(value, np.ndarray):
            return dict(shape=list(value.shape), dtype=str(value.dtype), entries=display(value.tolist()))
        if isinstance(value, sp.MatrixBase):
            return dict(shape=list(value.shape), entries=[[str(value[i, k])
                        for k in range(value.cols)] for i in range(value.rows)])
        if isinstance(value, UndefinedFunction):
            return undefined_function_metadata(value)
        if isinstance(value, sp.Basic):
            return str(value)
        raise TypeError('unsupported readable object ' + str(type(value)))

    def premise_show(value):
        # Exact source show() convention; used only for premise-rendered joins.
        if isinstance(value, dict): return {str(a): premise_show(b) for a, b in value.items()}
        if isinstance(value, (list, tuple)): return [premise_show(x) for x in value]
        if isinstance(value, sp.MatrixBase): return premise_show(value.tolist())
        if isinstance(value, np.ndarray): return premise_show(value.tolist())
        if isinstance(value, complex): return dict(real=value.real, imag=value.imag)
        if isinstance(value, (np.integer, np.floating, np.bool_)): return value.item()
        if isinstance(value, sp.Basic): return str(value)
        return value

    def inventory(value):
        if isinstance(value, UndefinedFunction):
            return undefined_function_metadata(value)
        if isinstance(value, dict):
            return {str(key): inventory(item) for key, item in value.items()}
        if isinstance(value, (tuple, list)):
            return [inventory(item) for item in value]
        if isinstance(value, (set, frozenset)):
            return dict(type=type(value).__name__,
                        entries=[inventory(item) for item in sorted(value, key=str)])
        if isinstance(value, slice):
            return dict(type='slice', start=inventory(value.start),
                        stop=inventory(value.stop), step=inventory(value.step))
        if isinstance(value, np.ndarray):
            return dict(type='ndarray', shape=list(value.shape), dtype=str(value.dtype))
        if isinstance(value, sp.MatrixBase):
            return dict(type=type(value).__name__, shape=list(value.shape),
                        symbols=sorted(str(symbol) for symbol in value.free_symbols))
        if isinstance(value, sp.Basic):
            return dict(type=type(value).__name__, symbols=sorted(str(s) for s in value.free_symbols))
        return display(value)

    def emit(name, value):
        relative = name + '.json'
        if relative in spec['priorSummaryFiles']:
            journal.copy_summary(relative)
        else:
            save(out / relative, display(value))

    # Plain JSON receipts are themselves journal operands, never success surrogates.
    parent = {name: read_json(key) for name, key in spec['job1Files'].items()
              if name.endswith('.json')}
    checkpoint = read_json('job1Checkpoint')
    parent_spec = parent['manifest.json']
    operations = parent['operation-index.json']
    def receipt_checks(metadata, cp):
        observed = {name: route(inputs[key]['path']) for name, key in spec['job1Files'].items()}
        root = Path(spec['job1Root'])
        actual_names = {str(path.relative_to(root)) for path in root.rglob('*') if path.is_file()}
        source_posthashes = metadata['posthashes.json']
        return dict(
            artifactCensus=actual_names == set(observed),
            checkpointArtifacts=observed == cp['artifacts'],
            checkpointChecks=metadata['checks.json'] == cp['checks'],
            workerArtifactIndex={key: record for key, record in observed.items()
                                 if key != 'artifact-index.json'} == metadata['artifact-index.json'],
            sourcePrePost=source_posthashes == metadata['manifest.json']['inputs'],
            actualSourceRoutes=all(route(record['path']) == record
                                   for record in source_posthashes.values()),
            stdoutChecksBytes=Path(inputs['job1Stdout']['path']).read_bytes()
                              == Path(inputs[spec['job1Files']['checks.json']]['path']).read_bytes(),
            strictStderrEmpty=all(Path(inputs[key]['path']).stat().st_size == 0
                                  for key in spec['job1StderrKeys']))
    receipt_result = journal.op('job1-receipts', receipt_checks, parent, checkpoint)
    emit('job1-receipts', receipt_result)
    require(all(receipt_result.values()), 'job1 receipt join failed; evidence saved')

    completed, failures, restored = {}, {}, {}
    for ordinal, previous in enumerate(operations):
        name = previous['name']
        folder = 'operations/%04d-%s/' % (ordinal, name)
        operand_values = [read_pickle(by_path[record['path']]) for record in previous['operands']]
        status = previous['status']
        receipt_name = 'complete.json' if status == 'COMPLETE' else 'unresolved.json'
        receipt = parent[folder + receipt_name]
        started = parent[folder + 'started.json']
        input_receipt = parent[folder + 'input-references.json']
        reused_detail = (journal.prior_json('prior-operation-%02d.json' % ordinal)
                         if journal.has_prior('inspect-prior-%02d' % ordinal) else None)
        joins = reused_detail['joins'].copy() if reused_detail is not None else dict(indexReceipt=previous == receipt,
                     inputReceipt=(input_receipt['name'] == name
                                   and input_receipt['operands'] == previous['operands']
                                   and input_receipt['startedUtc'] == previous['startedUtc']),
                     startedReceipt=(started['name'] == name
                                     and started['startedUtc'] == previous['startedUtc']))
        if status == 'COMPLETE':
            value = read_pickle(by_path[previous['value']['path']])
            completed[name] = value
            if 'restore-' in name:
                require(len(operand_values) == 1, 'saved restore operand schema')
                source_name = operand_values[0]
                source_key = spec['job1SourceInputs'][source_name]
                source_value = read_pickle(source_key)
                result = journal.op('source-join-' + name,
                                    lambda a, b: dict(structuralIdentity=same(a, b),
                                                      evidenceKind='SAVED_SOURCE_IDENTITY',
                                                      returnedInventory=inventory(a)),
                                    value, source_value)
                joins['restoredSourceStructuralIdentity'] = result['structuralIdentity']
                restored[name] = source_name
            if reused_detail is not None:
                detail = journal.op('inspect-prior-%02d' % ordinal, None,
                                    previous, operand_values, value, joins)
                emit('prior-operation-%02d' % ordinal, detail)
                require(all(detail['joins'].values()), 'saved prior operation joins; no recomputation')
                continue
            rendered = parent.get(name + '.json')
            if isinstance(rendered, dict) and 'completeObjectReceipt' in rendered:
                joins['renderedReceipt'] = rendered['completeObjectReceipt'] == previous
            elif rendered is not None:
                joins['renderedReturn'] = premise_show(value) == rendered
            detail = dict(name=name, status=status, sourceReceipt=previous,
                          returnedInventory=inventory(value), joins=joins)
        else:
            require(status == 'UNRESOLVED_NO_RETRY', 'unsupported prior operation status')
            failures[name] = dict(receipt=previous, operandInventory=inventory(operand_values),
                                  returnedValuePresent='value' in previous)
            joins['noInventedReturn'] = 'value' not in previous
            joins['failureReceiptRoute'] = previous['failure']['operationReceipt'] == inputs[
                spec['job1Files'][folder + 'unresolved.json']]['path']
            detail = dict(name=name, status=status, sourceReceipt=previous,
                          operandInventory=inventory(operand_values), joins=joins)
        # Save actual operands and returned evidence before checking any joins.
        journal.op('inspect-prior-%02d' % ordinal, lambda *args: detail,
                   previous, operand_values, completed.get(name), joins)
        emit('prior-operation-%02d' % ordinal, detail)
        require(all(joins.values()), 'prior operation join failed; evidence saved: ' + name)

    def context_joins(uniform, reduction, units, branch_input, branch_return, unit_json):
        state = reduction['reductionState']
        equations, mapping = state['branch_equations'], state['branch_map']
        equation_joins = [bool(eq.lhs in mapping and eq.rhs == mapping[eq.lhs]) for eq in equations]
        unit_operands = dict(unitFrame=parent_spec['physicalInput']['unit_frame'],
                             fields=units['fieldUnits'], current=units['currentUnit'],
                             strong=uniform['units']['strong'])
        branch_operands = dict(equations=equations, mapping=mapping,
                               tangents=state['tangents'], groups=state['momentum_groups'],
                               normalMap=state['normal_map'])
        joins = dict(branchOperationInput=same(branch_input, state),
                     recordedEquationCount=len(equations) == spec['expectedBranchEquationCount'],
                     branchReturnedOperands=all(same(branch_return[key], value)
                                               for key, value in branch_operands.items()),
                     branchEquationStructuralJoins=all(equation_joins),
                     unitRenderedJoin=premise_show(unit_operands) == unit_json,
                     returnedResidualCount=len(branch_return['residuals']) == len(equations))
        return dict(operands=branch_operands, savedBranchReturn=branch_return,
                    expectedEquationCount=spec['expectedBranchEquationCount'],
                    evidenceKind='SAVED_BRANCH_AND_UNIT_STRUCTURAL_IDENTITY',
                    unitOperands=unit_operands, equationStructuralJoins=equation_joins,
                    savedResidualZeroFlags=[bool(value == 0) for value in branch_return['residuals']],
                    joins=joins)
    branch_record = next(row for row in operations if row['name'] == 'source-branch-map-joins')
    context = journal.op('branch-unit-joins', context_joins,
                         completed['restore-uniform'], completed['restore-branch-context'],
                         completed['restore-unit-context'],
                         read_pickle(by_path[branch_record['operands'][0]['path']]),
                         completed['source-branch-map-joins'], parent['unit-context.json'])
    emit('branch-unit-joins', context)
    require(all(context['joins'].values()), 'branch/unit structural join failed; evidence saved')

    end_results = {}
    for end in ENDS:
        record = next(row for row in operations if row['name'] == end + '-source-bind')
        arguments = [read_pickle(by_path[item['path']]) for item in record['operands']]
        original, pair = completed[end + '-restore-original-symbol'], completed[end + '-restore-current']
        def inspect_end(label, args, source, current, receipt):
            joins = dict(label=args[0] == label, originalSymbolStructuralIdentity=same(args[1], source),
                         currentPairStructuralIdentity=same(args[2], current))
            # The saved pair is the source tuple, never a successfully bound map.
            return dict(joins=joins, sourceSymbolInventory=inventory(source),
                        evidenceKind='SAVED_OPERAND_IDENTITY_NOT_INDEPENDENT_PHYSICS',
                        currentPairInventory=inventory(current),
                        savedInputTupleLength=len(current),
                        sourceBindStatus=receipt['status'], sourceBindValuePresent='value' in receipt,
                        savedFailure=receipt.get('failure'),
                        sourceEndState=parent['source-end-states.json'].get(label))
        observed = journal.op('end-saved-input-' + end, inspect_end, end, arguments, original, pair, record)
        emit('end-saved-input-' + end, observed)
        require(all(observed['joins'].values()), 'end saved-operand join failed; evidence saved')
        end_results[end] = observed

    def coverage(metadata):
        planned = metadata['planned-grid.json']
        visited, unvisited = metadata['grid-summary.json'], metadata['unvisited-grid.json']
        expected = {(end, point['omega'], point['ray']) for end in planned['ends'] for point in planned['pairs']}
        seen = {(row['end'], row['omega'], row['ray']) for row in visited}
        missing = {(row['end'], row['omega'], row['ray']) for row in unvisited}
        selected = {}
        for omega in spec['comparisonFrequencies']:
            selected[omega] = {end: dict(
                visited=[row for row in visited if row['end'] == end and row['omega'] == omega
                         and row['ray'] == spec['comparisonRay']],
                unvisited=[row for row in unvisited if row['end'] == end and row['omega'] == omega
                           and row['ray'] == spec['comparisonRay']]) for end in ENDS}
        selected_counts = [dict(end=end, omega=omega, ray=spec['comparisonRay'],
                                planned=(end, omega, spec['comparisonRay']) in expected,
                                savedRowCount=len(rows['visited']) + len(rows['unvisited']))
                           for omega, ends in selected.items() for end, rows in ends.items()]
        return dict(joins=dict(disjoint=not bool(seen & missing), completeAccounting=seen | missing == expected,
                               uniqueVisited=len(seen) == len(visited), uniqueUnvisited=len(missing) == len(unvisited),
                               summaryCount=metadata['checks.json']['points'] == len(visited),
                               selectedTripleCount=len(selected_counts) == spec['expectedComparisonTripleCount'],
                               selectedTriplesPlanned=all(row['planned'] for row in selected_counts),
                               selectedTriplesAccountedOnce=all(row['savedRowCount'] == 1 for row in selected_counts)),
                    planned=len(expected), visited=len(visited), unvisited=len(unvisited),
                    selectedTripleCounts=selected_counts,
                    selectedSavedRows=selected, sourceEndStates=metadata['source-end-states.json'],
                    savedReferenceFace=metadata['reference-face-drive.json'],
                    savedSeeds=metadata['saved-seed-comparison.json'], savedLoci=metadata['threshold-loci.json'])
    grid = journal.op('saved-coverage', coverage, parent)
    emit('saved-coverage', grid)
    require(all(grid['joins'].values()), 'saved coverage accounting mismatch; evidence saved')

    # Selected historical numerical arrays, used only at their saved contrasts.
    response = read_pickle('oldResponse')
    remainder = read_pickle('oldRemainders')
    systems = read_pickle('oldSystems')
    def inspect_response(response, remainder, systems, unit_source):
        grades = list(response['solve']['coefficients'])
        vectors = []
        for contrast, entry in remainder.items():
            direct, retained, difference = (entry[key] for key in ('direct', 'retained', 'difference'))
            # Bookkeeping identity for a difference produced from these same operands.
            arithmetic = direct - retained - difference
            vectors.append(dict(contrast=contrast, direct=direct, retained=retained,
                                savedDifference=difference, subtractionBookkeepingIdentity=arithmetic,
                                savedEquationResidual=entry['directEquationResidual'],
                                savedMaximum=entry['maximumReferenceFrame']))
        joins = dict(remainderPacketStructuralIdentity=same(response['remainders'], remainder),
                     job1UnitSourceStructuralIdentity=same(response, unit_source),
                     coefficientGrades=set(grades) == set(systems['matrices']) == set(systems['rhs']),
                     coefficientSystemShapes=all(systems['matrices'][grade].shape == (5*response['size'], 5*response['size'])
                                                and systems['rhs'][grade].shape == response['solve']['coefficients'][grade].shape
                                                for grade in grades))
        return dict(joins=joins, vectors=vectors, grades=grades,
                    suppliedSourceConventions=spec['oldResponseSuppliedConventions'],
                    joinEvidenceKind='SAVED_PACKET_IDENTITY_AND_SHAPE_CHECKS',
                    coefficientIngredients={key: inventory(systems[key]) for key in ('matrices', 'rhs')},
                    savedCoefficients=response['solve']['coefficients'],
                    savedIndependentCoefficients=response['solve']['independentCoefficients'],
                    savedCoefficientResiduals=response['solve']['residual'],
                    savedScaledResiduals=response['solve']['scaledResidual'],
                    savedIndependentDifference=response['solve']['independentDifference'],
                    fieldUnits=response['fieldUnits'], rowUnits=response['rowUnits'],
                    currentUnit=response['currentUnit'], settings=response['settings'],
                    ratio=response['ratio'], savedScope=response['scope'],
                    labels=response['response']['labels'],
                    responseInventory=inventory(response['response']),
                    savedOpenFluxInventory=inventory(response['flux']),
                    savedOpenFluxScope=response['flux']['scope'],
                    savedInputPackets=response['inputPackets'],
                    selectedPacketTopLevelKeys=list(response),
                    topLevelNamedKeyProbe_notACensus={key: key in response
                                                     for key in spec['limitedComparatorKeyProbe']})
    anchor = journal.op('old-saved-contrast-vectors', inspect_response,
                        response, remainder, systems, completed['restore-unit-context'])
    emit('old-saved-contrast-vectors', anchor)
    require(all(anchor['joins'].values()), 'historical response packet join failed; evidence saved')

    # Already-readable historical operand exports are copied without parsing math.
    old_json = {role: read_json(key) for role, key in spec['oldReadable'].items()}
    def inspect_old_text(payloads, old_checkpoint, physical):
        context = payloads['context']
        return dict(physicalInput=context['physicalInput'],
                    physicalInputStructuralJoin=context['physicalInput'] == physical,
                    fieldUnits=context['fieldEntryUnits'], forcingUnits=context['forcingUnits'],
                    cases={name: dict(savedCase=row['case'], savedPlaneJetEstablished=row['planeJetEstablished'],
                                     savedForcingPairingEstablished=row['forcingPairingEstablished'],
                                     nativeTermsBySide=row['nativeTermsBySide'],
                                     regulatorPresentCount=row['regulatorPresentCount'],
                                     unreducedCoefficientZeroFlags=row['unreducedCoefficientZeroFlags'])
                           for name, row in payloads.items() if name.startswith('case-')},
                    ends={name: dict(savedEndField=row['endField'], fieldShape=row['field']['shape'],
                                    savedFieldNonzeroEntryCount=len(row['field']['entries']),
                                    savedPolynomialNonzeroEntryCount=len(row['polynomial']['entries']))
                          for name, row in payloads.items() if name.startswith('end-')},
                    originalReviewVerdicts=old_checkpoint['reviewVerdicts'],
                    originalAcceptanceScope=old_checkpoint['acceptanceScope'],
                    sourceReaderRoute=inputs['oldInspectionSource'])
    old_summary = journal.op('old-readable-inventory', inspect_old_text,
                             old_json, read_json('oldInspectionCheckpoint'), parent_spec['physicalInput'])
    for role, value in old_json.items():
        emit('old-readable-' + role, value)
    emit('old-readable-inventory', old_summary)
    require(old_summary['physicalInputStructuralJoin'], 'old saved binding mismatch; evidence saved')

    return dict(status='SAVED_READER_COMPLETED_REQUIRES_REPORT',
                completedParentOperations=len(completed), unresolvedParentOperations=len(failures),
                savedObjectsRead=len(cache), parentVisitedPoints=grid['visited'],
                parentUnvisitedPoints=grid['unvisited'], parentEndStates=grid['sourceEndStates'],
                sourceBindingReturnPresence={end: row['sourceBindValuePresent'] for end, row in end_results.items()},
                comparisonRows=grid['selectedSavedRows'], oldContrastCount=len(anchor['vectors']),
                oldReadableRecords=len(old_json), oldSavedReviewVerdicts=old_summary['originalReviewVerdicts'],
                physicalLightCalibration='OPEN_NOT_ESTABLISHED', scienceJobOrdinal=2,
                scientificExecutionOrdinal=3, restoredPriorCompleteOperations=len(journal.reused),
                resumedPriorIncompleteOperation=journal.resumed_incomplete,
                stopAfterJob2AndFrequencyComparison=True, automaticRetry=False)


def validate_gate(args, spec, gate):
    require(gate.get('status') == 'READY_FOR_ONE_GUARDED_TRANSVERSE_FACE_SAVED_READER_CONTINUATION',
            'new saved-reader gate required')
    require(gate.get('scienceJobOrdinal') == 2 and gate.get('independentBuildClearance') is True,
            'job-2 independent clearance required')
    continuation_route = gate['continuationAuthorizationRecord']
    require(continuation_route == spec['inputs']['continuationAuthorization']
            and route(continuation_route['path']) == continuation_route, 'continuation approval pin')
    continuation = json.loads(Path(continuation_route['path']).read_text())
    require(continuation['approvedSavedReaderContinuation'] is True
            and continuation['completedOperationReplayAuthorized'] is False
            and continuation['unfinishedPremiseBindingAuthorized'] is False
            and continuation['methodRouteDecisionAuthorized'] is False
            and continuation['fieldPowerConstructionAuthorized'] is False
            and continuation['stopAfterReaderAndFrequencyComparison'] is True
            and continuation['scientificExecutionOrdinal'] == gate.get('scientificExecutionOrdinal') == 3
            and continuation['maximumNewExecutions'] == 1
            and continuation['scientificExecutionsUsedBefore'] == 2
            and continuation['scientificExecutionCeiling'] == 4, 'bounded continuation scope')
    require(continuation['parentFailureCheckpoint'] == spec['inputs']['priorFailureCheckpoint'],
            'continuation failure target changed')
    approval = gate['userScienceApprovalRecord']
    require(approval == continuation['originalJobsApproval'], 'original approval join')
    require(route(approval['path']) == approval, 'approval receipt changed')
    authorization = json.loads(Path(approval['path']).read_text())
    require(authorization['authorizedScienceJobOrdinals'] == [1, 2]
            and authorization['stopAfterJob2AndFrequencyComparison'] is True
            and authorization['methodRouteDecisionAuthorized'] is False
            and authorization['fieldPowerConstructionAuthorized'] is False,
            'approved jobs-1-2 scope required')
    require(gate['workerSha256'] == sha(__file__)
            and gate['inputManifestSha256'] == sha(args.input_manifest), 'reader/source pins')
    review_route = gate['reviewRecord']
    require(route(review_route['path']) == review_route, 'new review record changed')
    review = json.loads(Path(review_route['path']).read_text())
    require(review.get('independentBuildClearance') is True, 'review clearance absent')
    packet_hashes = review['packetState']['fileHashes']
    for path in (Path(__file__), args.input_manifest):
        require(packet_hashes[str(path.resolve().relative_to(ROOT))] == sha(path),
                'review did not cover exact reader/manifest')
    limits = dict(seconds=900, nativeSeconds=840, memoryGiB=2, swapMax=0,
                  cpuCount=1, nice=15, tasksMax=32, nativeThreads=1, automaticRetry=False)
    require(all(gate.get(key) == value for key, value in limits.items()), 'ordinary reader limits')
    require(gate['guardSha256'] == sha(ROOT / 'scripts/s11c_guarded_run.py')
            and gate['supervisorSha256'] == sha(ROOT / 'research/pde_ledger_v3/_measurements/S11c_d_end_normalization_run.py'),
            'shared helpers changed')
    require(spec['scienceJobOrdinal'] == 2 and spec['scientificExecutionOrdinal'] == 3
            and spec['automaticRetry'] is False,
            'manifest scope')


def main():
    parser = argparse.ArgumentParser(__doc__)
    for name in ('input-manifest', 'gate-receipt', 'run-directory'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    args.input_manifest = args.input_manifest.resolve(strict=True)
    args.gate_receipt = args.gate_receipt.resolve(strict=True)
    spec = json.loads(args.input_manifest.read_text())
    gate = json.loads(args.gate_receipt.read_text())
    validate_gate(args, spec, gate)
    limits = containment()
    started = time.monotonic()
    out = args.run_directory.resolve()
    out.relative_to(STORE)
    require(out != STORE, 'dedicated reader output directory required')
    out.mkdir(parents=True, exist_ok=False)
    journal = Journal(out, spec)
    try:
        save(out / 'native-limits.json', limits)
        save(out / 'manifest.json', spec)
        save(out / 'gate.json', gate)
        launch_routes = {name: route(path) for name, path in {
            'worker': Path(__file__).resolve(), 'manifest': args.input_manifest,
            'gate': args.gate_receipt, 'review': Path(gate['reviewRecord']['path'])}.items()}
        save(out / 'launch-source-routes.json', launch_routes)
        for name, record in launch_routes.items():
            publish(out / 'launch-sources' / (name + Path(record['path']).suffix),
                    Path(record['path']).read_bytes())
        prehashes = {key: route(record['path']) for key, record in spec['inputs'].items()}
        save(out / 'prehashes.json', prehashes)
        require(prehashes == spec['inputs'], 'input changed before saved restoration')
        original = spec['originalSpecification']
        require(all(spec['inputs'][key] == value for key, value in original['inputs'].items()),
                'original 168 input routes changed')
        require(all(spec[key] == value for key, value in original.items() if key != 'inputs'),
                'original scientific selection or reader specification changed')
        journal.configure_prior()
        checkpoint_joins = journal.copy_summary('checkpoint-joins.json')
        require(all(row['hashMatches'] and row['sizeMatches'] for row in checkpoint_joins),
                'saved selected checkpoint joins; no recomputation')
        result = inspect_saved(spec, out, journal)
        require(journal.reused == journal.prior_order and journal.resumed_incomplete,
                'continuation did not consume its complete prefix and unfinished operation')
        save(out / 'complete-summary-reuse.json', journal.summary_reuse)
        post = {key: route(record['path']) for key, record in spec['inputs'].items()}
        save(out / 'posthashes.json', post)
        require(post == prehashes, 'source changed during saved inspection')
        launch_post = {name: route(record['path']) for name, record in launch_routes.items()}
        save(out / 'launch-source-posthashes.json', launch_post)
        require(launch_post == launch_routes, 'worker/manifest/gate/review changed during inspection')
        # Scientific inspection and pinning are complete; outer containment
        # still bounds final receipt/stdout publication.
        signal.setitimer(signal.ITIMER_REAL, 0)
        result.update(wallSeconds=time.monotonic() - started, operations=len(journal.records),
                      allSourceRoutesUnchanged=post == prehashes)
        save(out / 'operation-index.json', journal.records)
        payload = encoded(result)
        publish(out / 'checks.json', payload)
        # The exact bytes written to checks are the only successful stdout.
        sys.stdout.buffer.write(payload)
        sys.stdout.buffer.flush()
    except BaseException:
        signal.setitimer(signal.ITIMER_REAL, 0)
        save(out / 'failure.json', dict(traceback=traceback.format_exc(),
             incompleteOperation=journal.active, wallSeconds=time.monotonic() - started, noRetry=True))
        if not (out / 'operation-index.json').exists():
            save(out / 'operation-index.json', journal.records)
        post = {}
        for key, record in spec['inputs'].items():
            try:
                post[key] = route(record['path'])
            except Exception as error:
                post[key] = dict(errorType=type(error).__name__, error=str(error))
        save(out / 'failure-posthashes.json', post)
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        save(out / 'artifact-index.json', {str(path.relative_to(out)): route(path)
             for path in sorted(out.rglob('*')) if path.is_file()})


if __name__ == '__main__':
    main()
