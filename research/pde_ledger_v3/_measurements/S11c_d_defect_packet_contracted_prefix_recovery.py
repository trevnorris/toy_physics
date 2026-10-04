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
        witnesses = []
        for witness in contract['panels']:
            found = index.lookup(witness['namespace'], witness['descriptor'])
            require(found is not None and found['inputReceipt'] == witness['inputReceipt'] and
                    found['receipt'] == witness['returnReceipt'] and
                    found['descriptor']['kind'] == 'inner-panel', 'full completed panel witness')
            require(found['descriptor']['context'] == expected['descriptor']['context'], 'same mathematical context')
            require(found['descriptor']['route'] == 'A24' and found['descriptor']['purpose'] == 'baseline' and
                    found['descriptor']['precision'] == 30, 'same actual own route and precision')
            # Selected rule A48 is a paired rule within OWN A24 route, never an
            # A48 outer-route result. Exact context/point joins are worker duties.
            witnesses.append(found['descriptor']['settings']['selectedRule'])
        require(sorted(witnesses) == ['A24', 'A48'], 'both original paired panels')
        return {'evidence': evidence, 'completeRequests': complete,
                'pending': pending[0], 'panels': contract['panels'],
                'recordsFingerprint': rows_digest(db, 'SELECT * FROM records ORDER BY sequence'),
                'chunksFingerprint': rows_digest(db, 'SELECT * FROM chunks ORDER BY sequence,ordinal'),
                'requestFingerprint': rows_digest(db, 'SELECT * FROM request_index ORDER BY namespace,request_hash'),
                'scientificAcceptance': False}


class ExplicitCloneStore(EvidenceStore):
    """New byte-identical branch of a closed pinned failure; never modify history."""
    def __init__(self, original, destination, contract):
        original = Path(original); destination = Path(destination)
        require(original.resolve() != destination.resolve(), 'new branch, not original database')
        require(not destination.exists(), 'exclusive new branch; no repeated recovery')
        self.original = original; self.contract = contract
        self.admission = audit_prefix(original, contract)
        with original.open('rb') as inp, destination.open('xb') as out:
            shutil.copyfileobj(inp, out, 1048576)
            out.flush()
            import os
            os.fsync(out.fileno())
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
                    'parentState': state, 'scientificAcceptance': False}
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
        if (namespace, key) != (expected['namespace'], expected['requestHash']):
            return None
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
