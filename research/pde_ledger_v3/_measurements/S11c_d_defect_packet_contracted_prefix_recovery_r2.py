"""Proposed task-local, explicit single-parent recovery storage. NOT a worker.

Only stdlib JSON/SQLite byte transport. Scientific use requires an independently
assessed concrete worker, fresh exact gate, and the shared guarded supervisor.
No caller is launched here. Shared EvidenceStore and RequestIndex stay unchanged.
Tests may create manufactured databases; original scientific databases are only
read/copied as opaque bytes outside containment, never attached for writing.
"""
from collections import OrderedDict
import hashlib
import json
import os
from pathlib import Path
import shutil
import sqlite3

from S11c_d_defect_packet_evidence_store import EvidenceReader, EvidenceStore
from S11c_d_defect_packet_request_index import RequestIndex, canonical, require


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def panel_evidence(panels):
    """Only the metadata boundary changes. Numerical integer lookups stay intact."""
    require(type(panels) is dict and len(panels) == 2 and
            all(type(k) is int for k in panels) and set(panels) == {24, 48},
            'exact original panel-order keys; no generic stringification')
    result = {'24': panels[24], '48': panels[48]}
    require(len(result) == len(panels) and all(result[str(k)] is v for k, v in panels.items()),
            'lossless panel identity and collision refusal')
    return result


def readonly(path):
    db = sqlite3.connect(Path(path).resolve().as_uri() + '?mode=ro&immutable=1', uri=True)
    db.execute('PRAGMA query_only=ON')
    return db


def row_bytes(row):
    # Type-tagged stdlib fingerprint, no scientific restoration or arithmetic.
    return canonical([{'bytes': value.hex()} if type(value) is bytes else value for value in row])


def rows_digest(db, sql, args=()):
    h = hashlib.sha256(); count = 0
    for row in db.execute(sql, args):
        raw = row_bytes(row)
        h.update(str(len(raw)).encode() + b':' + raw); count += 1
    return {'rows': count, 'sha256': h.hexdigest()}


def exact(left, right):
    return canonical(left) == canonical(right)


def join_frame(expected, node, panels):
    """Strict encoded operand joins only; numerical geometry is not rederived."""
    parent = expected['descriptor']; ps = parent['settings']; pc = ps['context']
    require(exact(node['panel'], ps) and exact(node['originalA'], parent['arguments']['a']) and
            exact(node['originalB'], parent['arguments']['b']), 'actual parent settings and endpoints')
    require(type(node['side']) is int and node['side'] == 0 and
            type(node['ruleIndex']) is int and node['ruleIndex'] == 0,
            'first original node; no unpublished completed-loop prefix')
    require(ps['kind'] == 'outer' and ps['selectedRule'] == ps['ownRoute'] == 'A24' and
            ps['squareMapped'] is True and exact(pc['cell'], {'slabId': 0}) and
            type(pc['n']) is int and pc['n'] == 0 and pc['purpose'] == 'baseline' and
            type(pc['K']) is int and pc['K'] == 27 and type(pc['T']) is int and pc['T'] == 122,
            'exact first action and slab')
    require(len(panels) == 2, 'two complete paired panels')
    common = None; args = None; rules = []
    for witness in panels:
        d = witness['descriptor']; settings = d['settings']; context = settings['context']
        require(exact(d['context'], parent['context']) and d['route'] == parent['route'] == 'A24' and
                d['purpose'] == parent['purpose'] == 'baseline' and
                type(d['precision']) is int and d['precision'] == parent['precision'] == 30,
                'same entire mathematical context and own route')
        require(settings['kind'] == 'inner' and settings['ownRoute'] == 'A24' and
                settings['squareMapped'] is True and settings['selectedRule'] in ('A24', 'A48'),
                'actual paired inner rule')
        require(exact(settings['ruleReceipt'],
                      parent['context']['exactFamilyContexts']['ruleReceipts'][settings['selectedRule']]),
                'actual immutable selected-rule receipt')
        require(exact(context['outerParent'], node) and exact(context['m'], node['representedPoint']) and
                exact(context['cell'], {'cellId': 0, 'clipped': True}), 'same original point and inner cell')
        require(all(exact(context[k], pc[k]) for k in ('K', 'T', 'carrier', 'n')),
                'actual window/carrier/jet arguments')
        require(context['inputClips'] == 'J,Dh,Dq,Dr_wrong_root' and context['outputClips'] == 'Dr' and
                exact(context['primitiveKeys'], ['J', 'Dr', 'Dh', 'Dq']) and
                exact(context['semanticVariables'], {'X':'k', 'Y':'l'}) and
                context['sharedQuadratureCoordinateDoesNotIdentifyKAndL'] is True,
                'unaltered physical input/output contraction context')
        if common is not None:
            require(exact(common, context) and exact(args, d['arguments']), 'same full paired intervals/context')
        common = context; args = d['arguments']; rules.append(settings['selectedRule'])
    require(sorted(rules) == ['A24', 'A48'], 'one of each original paired rule')
    return {'parent': parent, 'node': node, 'panelArguments': args, 'context': common}


def audit_prefix(path, contract):
    """Complete byte/descriptor audit against an externally pinned prefix."""
    path = Path(path)
    require(path.is_file() and sha(path) == contract['databaseSha256'], 'fixed original database bytes')
    require(path.stat().st_size == contract['databaseBytes'], 'fixed original database length')
    require(all(not Path(str(path) + suffix).exists() for suffix in ('-wal', '-shm', '-journal')),
            'closed database, no sidecars or live writer')
    with EvidenceReader(path) as reader:
        evidence = reader.audit(contract['records'], contract['head'])
        db = reader.db
        require(db.execute('PRAGMA integrity_check').fetchone() == ('ok',), 'SQLite integrity')
        require(not list(db.execute('PRAGMA foreign_key_check')), 'SQLite foreign-key integrity')
        require(db.execute('SELECT count(*) FROM contracted_leaves').fetchone() == (0,),
                'no completed or active outer leaves in this recovery')
        # The inspection adapter has only RequestIndex.read_record/_request
        # methods and a read-only connection, not a scientific decoder.
        index = object.__new__(RequestIndex)
        index.db = db; index.maximum = contract['maxRecordBytes']; index.cache_limit = 0
        index.cache = OrderedDict(); index.cache_size = 0; index.cache_epoch = None
        expected = contract['pending']
        complete = 0; pending = []
        for ns, key, raw, state, inp, out in db.execute(
                'SELECT namespace,request_hash,descriptor,state,input_sequence,output_sequence FROM request_index'):
            value = json.loads(raw)
            got_key, got_raw = index._request(ns, value)
            require(got_key == key and got_raw == bytes(raw), 'full request and namespace identity')
            original, ir = index.read_record(inp)
            require(original == value and ir['namespace'] == ns and ir['record'] == 'request/' + key + '/input',
                    'exact request input receipt')
            if state == 'COMPLETE':
                _, rr = index.read_record(out)
                require(out > inp and rr['namespace'] == ns and rr['record'] == 'request/' + key + '/return',
                        'complete request return ancestry')
                complete += 1
            else:
                require(state == 'PENDING' and out is None, 'only one declared pending state')
                pending.append({'namespace': ns, 'requestHash': key,
                                'inputReceipt': ir, 'descriptor': value})
        require(complete == contract['completedRequests'] and pending == [expected],
                'entire completed census and single exact incomplete parent')
        require(expected['descriptor']['kind'] == 'outer-panel' and
                expected['descriptor']['route'] == 'A24' and
                expected['descriptor']['purpose'] == 'baseline' and
                expected['descriptor']['precision'] == 30,
                'declared first fixed-route outer parent only')
        mathematical, mr = index.read_record(contract['mathematicalContextReceipt']['sequence'])
        require(mr == contract['mathematicalContextReceipt'] and
                expected['descriptor']['context']['mathematicalContextReceipt'] == mr and
                exact(expected['descriptor']['context']['exactFamilyContexts'], mathematical['contexts']),
                'actual record-zero context and original receipt, not a new namespace substitute')
        actual_node, node_receipt = index.read_record(contract['firstOuterNode']['receipt']['sequence'])
        require(node_receipt == contract['firstOuterNode']['receipt'] and
                exact(actual_node, contract['firstOuterNode']['value']) and
                node_receipt['namespace'] == expected['namespace'] and
                node_receipt['record'].endswith('/outer-node-input'), 'full original node receipt and operands')
        join_frame(expected, actual_node, contract['panels'])
        for witness in contract['panels']:
            found = index.lookup(witness['namespace'], witness['descriptor'])
            require(found is not None and found['inputReceipt'] == witness['inputReceipt'] and
                    found['receipt'] == witness['returnReceipt'] and
                    found['descriptor']['kind'] == 'inner-panel', 'full completed panel witness')
            value = found['value']; d = found['descriptor']
            require(exact(value['a'], d['arguments']['a']) and exact(value['b'], d['arguments']['b']) and
                    exact(value['settings'], d['settings']) and
                    set(value['components']) == {'Y0','X1','X1T','X2T','X02T','C0T','C1T','C2T','Y0CT','Y1CT','Y0C','Y1C'},
                    'actual complete panel return belongs to original full input')
        actual_namespaces = [{'digest': digest, 'descriptor': json.loads(raw)} for digest, raw in
                             db.execute('SELECT digest,descriptor FROM namespaces ORDER BY digest')]
        require(exact(actual_namespaces, contract['originalNamespaces']), 'all actual original namespace definitions')
        serials = [int(name.split('/')[0]) for (name,) in db.execute('SELECT name FROM records')
                   if name.split('/')[0].isdigit()]
        require(bool(serials), 'published numeric serials exist')
        serial = max(serials)
        require(serial == contract['lastCommittedNumericSerial'], 'derive last committed serial from original names')
        return {'evidence': evidence, 'completeRequests': complete,
                'pending': pending[0], 'panels': contract['panels'], 'lastCommittedNumericSerial': serial,
                'recordsFingerprint': rows_digest(db, 'SELECT * FROM records ORDER BY sequence'),
                'chunksFingerprint': rows_digest(db, 'SELECT * FROM chunks ORDER BY sequence,ordinal'),
                'requestFingerprint': rows_digest(db, 'SELECT * FROM request_index ORDER BY namespace,request_hash'),
                'scientificAcceptance': False}


class ExplicitCloneStore(EvidenceStore):
    """New byte-identical branch of a closed pinned failure; never modify history."""
    def __init__(self, original, destination, contract, free_bytes=None):
        original = Path(original); destination = Path(destination)
        require(original.resolve() != destination.resolve(), 'new branch, not original database')
        require(not destination.exists(), 'exclusive new branch; no repeated recovery')
        self.original = original; self.contract = contract
        self.admission = audit_prefix(original, contract)
        free_bytes = free_bytes or (lambda: shutil.disk_usage(destination.parent).free)
        free = free_bytes(); needed = contract['reserveBytes'] + 2*contract['databaseBytes'] + 131072
        require(type(free) is int and free >= needed, 'clone disk reserve before any destination creation')
        self.copy_admission = {'freeBytes': free, 'requiredFreeBytes': needed}
        # Any copy/fsync/refusal failure preserves the partial exclusive branch.
        # It is never reopened here; no cleanup or automatic retry is attempted.
        with original.open('rb') as inp, destination.open('xb') as out:
            shutil.copyfileobj(inp, out, 1048576)
            out.flush()
            os.fsync(out.fileno())
        directory = os.open(destination.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
        require(sha(original) == sha(destination) == contract['databaseSha256'], 'byte-identical branch admission')
        self.path = destination; self.db = sqlite3.connect(str(destination)); self.closed = False
        self.db.execute('PRAGMA journal_mode=DELETE'); self.db.execute('PRAGMA synchronous=FULL')
        self.db.execute('PRAGMA foreign_keys=ON'); self.db.execute('PRAGMA cache_size=-4096')
        self.db.execute('PRAGMA temp_store=FILE')
        require(self.db.execute('PRAGMA synchronous').fetchone() == (2,) and
                self.db.execute('PRAGMA journal_mode').fetchone() == ('delete',) and
                self.db.execute('PRAGMA foreign_keys').fetchone() == (1,), 'same durable writer settings')
        self.sequence = contract['records']; self.previous = contract['head']

    def prefix_preserved(self):
        """All immutable old rows and complete request links remain exact."""
        require(sha(self.original) == self.contract['databaseSha256'], 'original still unchanged')
        old = readonly(self.original)
        try:
            count = self.contract['records']
            for sql in ('SELECT * FROM records WHERE sequence<? ORDER BY sequence',
                        'SELECT * FROM chunks WHERE sequence<? ORDER BY sequence,ordinal'):
                require(rows_digest(old, sql, (count,)) == rows_digest(self.db, sql, (count,)),
                        'complete immutable evidence/chunk prefix')
            for table, key in (('metadata', 'key'), ('namespaces', 'digest')):
                for row in old.execute('SELECT * FROM ' + table):
                    got = self.db.execute('SELECT * FROM ' + table + ' WHERE ' + key + '=?', (row[0],)).fetchone()
                    require(got == row, 'all original metadata/namespace bytes')
            expected = self.contract['pending']
            state = None
            for row in old.execute('SELECT * FROM request_index'):
                got = self.db.execute('SELECT * FROM request_index WHERE namespace=? AND request_hash=?', row[:2]).fetchone()
                require(got is not None, 'all original requests remain')
                if row[:2] != (expected['namespace'], expected['requestHash']):
                    require(got == row, 'every completed request is immutable')
                else:
                    require(got[:3] == row[:3] and got[4] == row[4], 'pending parent original operands immutable')
                    require((got[3] == 'PENDING' and got[5] is None) or
                            (got[3] == 'COMPLETE' and type(got[5]) is int and got[5] >= count),
                            'one-way explicitly resumed parent only')
                    state = got[3]
            return {'originalDatabaseSha256': self.contract['databaseSha256'],
                    'oldRecords': count, 'completedRequestsPreserved': self.contract['completedRequests'],
                    'parentState': state, 'currentLeafRows': self.db.execute('SELECT count(*) FROM contracted_leaves').fetchone()[0],
                    'originalLeafTableWasEmpty': True, 'scientificAcceptance': False}
        finally:
            old.close()


class ExplicitParentIndex(RequestIndex):
    """Only a caller with the exact reviewed parent can claim its unfinished work.

    Original lookup still refuses PENDING. The future worker must first call
    claim_parent with the actual full descriptor; this creates no new input and
    invokes no callback. It is deliberately not a generic automatic-resume index.
    """
    def __init__(self, store, contract, free_bytes=None):
        require(type(store) is ExplicitCloneStore and store.contract == contract and not store.closed,
                'fixed audited cloned store')
        self.store = store; self.db = store.db
        self.reserve = contract['reserveBytes']; self.maximum = contract['maxRecordBytes']
        self.cache_limit = contract['cacheBytes']
        require(type(self.reserve) is int and self.reserve >= 0 and
                type(self.maximum) is int and self.maximum > 0 and self.cache_limit == 0,
                'same explicit storage capacities, no payload cache')
        self.cache = OrderedDict(); self.cache_size = 0; self.cache_epoch = None
        self.free_bytes = free_bytes or (lambda: shutil.disk_usage(store.path.parent).free)
        self.contract = contract; self.claimed = False; self.completed_parent = False
        self._space(0)

    def claim_parent(self, namespace, descriptor):
        key, raw = self._request(namespace, descriptor)
        expected = self.contract['pending']
        require((namespace, key) == (expected['namespace'], expected['requestHash']),
                'initial parent mismatch; never fall through to new work')
        require(not self.claimed and not self.completed_parent, 'explicit parent can be claimed only once')
        require(canonical(expected['descriptor']) == raw, 'entire authorized parent descriptor')
        row = self.db.execute('SELECT descriptor,state,input_sequence,output_sequence FROM request_index WHERE namespace=? AND request_hash=?',
                              (namespace, key)).fetchone()
        require(row == (raw, 'PENDING', expected['inputReceipt']['sequence'], None),
                'same still-pending parent, not a general retry')
        inp, receipt = self.read_record(row[2])
        require(canonical(inp) == raw and receipt == expected['inputReceipt'], 'saved exact parent input')
        self.claimed = True
        return {'namespace': namespace, 'requestHash': key,
                'inputReceipt': receipt, 'descriptor': json.loads(raw)}

    def complete(self, ticket, value):
        parent = self.contract['pending']
        is_parent = (ticket['namespace'], ticket['requestHash']) == (parent['namespace'], parent['requestHash'])
        if is_parent:
            require(self.claimed and not self.completed_parent, 'explicit parent claim before completion')
        result = super().complete(ticket, value)
        if is_parent:
            self.completed_parent = True
        return result


class ExplicitRecoveryProtocol:
    """Source-level state machine used by a future fixed evaluator wrapper.

    Callbacks return encoded JSON; this class never decodes/evaluates science.
    The first outer callback is unfinished work. Both saved child callbacks are
    unreachable. before_put skips ONLY the exact saved first node record; pair
    publication must commit before normal new work is released.
    """
    def __init__(self, index, contract):
        require(type(index) is ExplicitParentIndex and index.contract == contract, 'fixed clone index')
        self.index = index; self.store = index.store; self.contract = contract
        self.state = 'PARENT'; self.events = []; self.attempts = []; self.pair_prepared = None
        self.saved_values = {}; self.parent_done = False
        self.witnesses = {x['descriptor']['settings']['selectedRule']:x for x in contract['panels']}
        require(set(self.witnesses) == {'A24','A48'}, 'two pre-registered child addresses')
        self.namespaces = {}
        for item in contract['originalNamespaces']:
            d = item['descriptor']
            if d['route'] == 'mathematical-inputs': continue
            key = (d['route'], d['settings']['purpose'])
            require(key not in self.namespaces, 'unique original route/purpose namespace')
            self.namespaces[key] = item
        require(set(self.namespaces) == {('A24','baseline'),('A24','formula-check:baseline')},
                'exact original request and formula-check namespaces')

    def _attempt(self, value):
        # Freeze encoded operands before any decision; caller mutation cannot
        # rewrite the failed attempt later reported from finally.
        self.attempts.append(json.loads(canonical(value)))

    def namespace(self, descriptor):
        self._attempt({'stage':'namespace','state':self.state,'descriptor':descriptor})
        key = (descriptor['route'], descriptor['settings']['purpose'])
        old = self.namespaces.get(key)
        if old is not None:
            require(exact(descriptor, old['descriptor']), 'existing route namespace cannot silently change')
            row = self.store.db.execute('SELECT descriptor FROM namespaces WHERE digest=?', (old['digest'],)).fetchone()
            require(row is not None and bytes(row[0]) == canonical(descriptor), 'actual existing namespace required')
            return old['digest']
        require(self.state == 'NEW', 'no new namespace before exact prefix recovery')
        return self.store.namespace(descriptor)

    def request(self, namespace, descriptor, callback):
        if self.state != 'NEW':
            self._attempt({'stage':'request','state':self.state,'namespace':namespace,'descriptor':descriptor})
        if self.state == 'PARENT':
            ticket = self.index.claim_parent(namespace, descriptor)
            self.events.append({'event':'explicit-parent-claimed','inputReceipt':ticket['inputReceipt']})
            self.state = 'NODE'
            value = callback()
            require(self.state == 'NEW' and not self.parent_done, 'all exact prefix transitions before parent completion')
            receipt = self.index.complete(ticket, value)
            self.parent_done = True
            self.events.append({'event':'original-parent-completed','returnReceipt':receipt,
                                'oldEvidenceRecords':self.contract['records'], 'mixedExecutionProvenance':True})
            return value
        if self.state in ('CHILD24', 'CHILD48'):
            rule = 'A24' if self.state == 'CHILD24' else 'A48'; old = self.witnesses[rule]
            require(namespace == old['namespace'] and exact(descriptor, old['descriptor']),
                    'required completed child descriptor mismatch; no begin or callback')
            found = self.index.lookup(namespace, descriptor)
            require(found is not None and found['receipt'] == old['returnReceipt'] and
                    found['inputReceipt'] == old['inputReceipt'],
                    'required completed child lookup miss; no begin or callback')
            self.saved_values[rule] = found['value']
            self.events.append({'event':'completed-child-reused','selectedRule':rule,
                                'inputReceipt':found['inputReceipt'],'returnReceipt':found['receipt'],
                                'callbackInvoked':False})
            self.state = 'CHILD48' if rule == 'A24' else 'PAIR'
            return found['value']
        require(self.state == 'NEW', 'no new request before saved node and both children plus committed pair')
        found = self.index.lookup(namespace, descriptor)
        if found is not None:
            self.events.append({'event':'ordinary-exact-reuse','inputReceipt':found['inputReceipt'],
                                'returnReceipt':found['receipt'],'callbackInvoked':False})
            return found['value']
        ticket = self.index.begin(namespace, descriptor)
        value = callback(); self.index.complete(ticket, value); return value

    def before_put(self, route, purpose, name, encoded):
        if self.state != 'NEW':
            self._attempt({'stage':'before-put','state':self.state,'route':route,
                           'purpose':purpose,'name':name,'encoded':encoded})
        if self.state == 'NODE':
            node = self.contract['firstOuterNode']
            require(route == 'A24' and purpose == 'baseline' and name == 'outer-node-input' and
                    exact(encoded, node['value']), 'exact original first node input required')
            actual, receipt = self.index.read_record(node['receipt']['sequence'])
            require(receipt == node['receipt'] and exact(actual, encoded), 'actual original first-node receipt')
            self.events.append({'event':'first-node-input-verified-and-skipped','receipt':receipt,
                                'duplicateRecordWritten':False})
            self.state = 'CHILD24'; return receipt
        if self.state == 'PAIR':
            require(self.pair_prepared is None and route == 'A24' and purpose == 'baseline' and
                    name == 'inner-positive-panel-pair', 'first new numeric put is the reconstructed pair')
            witness = self.witnesses['A24']['descriptor']
            require(exact(encoded['a'], witness['arguments']['a']) and
                    exact(encoded['b'], witness['arguments']['b']) and
                    exact(encoded['meta'], witness['settings']['context']['cell']) and
                    type(encoded['ownOrder']) is int and encoded['ownOrder'] == 24,
                    'same first pair coordinates/cell/own order')
            require(exact(encoded['pairedPanels'], {'24':self.saved_values['A24'],'48':self.saved_values['A48']}),
                    'pair contains the complete restored panel returns')
            require(set(encoded['positiveComponentContributions']) == set(self.saved_values['A24']['components']),
                    'all unpublished scalar accumulator components explicitly reconstructed')
            self.pair_prepared = canonical(encoded)
            return None
        require(self.state == 'NEW', 'no unclassified evidence before required prefix transitions')
        return None

    def after_put(self, route, purpose, name, encoded, receipt):
        if self.state == 'PAIR':
            require(self.pair_prepared is not None and self.pair_prepared == canonical(encoded) and
                    (route,purpose,name) == ('A24','baseline','inner-positive-panel-pair'), 'prepared pair only')
            actual, rr = self.index.read_record(receipt['sequence'])
            require(rr == receipt and exact(actual, encoded) and receipt['sequence'] >= self.contract['records'] and
                    receipt['namespace'] == self.contract['pending']['namespace'] and
                    receipt['record'] == str(self.contract['lastCommittedNumericSerial'] + 1) + '/inner-positive-panel-pair',
                    'actual new committed pair, original scalar formula under separate AST obligation')
            self.events.append({'event':'unpublished-parent-scalars-reconstructed','receipt':receipt,
                                'restoredPublishedReturn':False, 'completedCallbacksInvoked':False})
            self.pair_prepared = None; self.state = 'NEW'
        else:
            require(self.state == 'NEW', 'unexpected post-write transition')

    def final_record(self):
        # Future worker must persist this from finally, even when science fails.
        return {'state':self.state,'parentCompleted':self.parent_done,'events':self.events,'attempts':self.attempts,
                'prefix':self.store.prefix_preserved(),'scientificAcceptance':False,
                'oldSequenceRange':[0,self.contract['records']-1],
                'newSequenceRange':[self.contract['records'],self.store.sequence-1]}
