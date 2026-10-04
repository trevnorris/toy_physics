"""Unintegrated stdlib tooling for a future guarded numerical packet worker.

No scientific imports, old-bank reads, numerical evaluation or automatic resume.
The caller supplies a NEW EvidenceStore writer. Full scientific adequacy of the
request descriptor, node operands and returned values still needs build review.
An interrupted/pending request refuses reuse; it is never recomputed here.
"""
from collections import OrderedDict
import hashlib
import json
import shutil
import zlib


class RequestRefusal(ValueError):
    pass


def require(value, message):
    if value is not True:
        raise RequestRefusal(message)


def canonical(value):
    # Reuse the writer's JSON contract by rejecting implicit key/type conversion.
    stack = [(value, False)]; active = set()
    while stack:
        item, leaving = stack.pop()
        if type(item) in (dict, list):
            if leaving:
                active.remove(id(item)); continue
            require(id(item) not in active, 'cyclic JSON refused')
            active.add(id(item)); stack.append((item, True))
            if type(item) is dict:
                require(all(type(k) is str for k in item), 'string JSON keys required')
                stack.extend((v, False) for v in item.values())
            else:
                stack.extend((v, False) for v in item)
        else:
            require(type(item) in (str, int, float, bool, type(None)), 'encoded JSON only')
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def digest(data):
    return hashlib.sha256(data).hexdigest()


class RequestIndex:
    """One new writer, immutable full-operand requests and durable return links.

    States record a one-way transition; scientific evidence records are immutable.
    A crash between evidence commit and index update leaves a refused incomplete
    request. This intentionally provides no recovery/retry policy. Lookup checks
    the complete payload, receipt and immediate chain link; completion inspection
    must still audit the whole evidence chain. No scientific validity is inferred.
    """
    def __init__(self, store, *, reserve_bytes=20*1024**3,
                 max_record_bytes=8*1024**2, cache_bytes=0, free_bytes=None):
        require(type(reserve_bytes) is int and reserve_bytes >= 0, 'disk reserve')
        require(type(max_record_bytes) is int and max_record_bytes > 0, 'record capacity')
        require(type(cache_bytes) is int and 0 <= cache_bytes <= 64*1024**2, 'cache capacity')
        self.store = store
        self.db = store.db
        require(not store.closed, 'new writer must be open')
        require(self.db.execute("SELECT name FROM sqlite_master WHERE name='request_index'").fetchone() is None,
                'no existing request-index attachment or automatic resume')
        self.reserve = reserve_bytes
        self.maximum = max_record_bytes
        self.cache_limit = cache_bytes
        self.cache = OrderedDict()
        self.cache_size = 0
        self.cache_epoch = None
        self.free_bytes = free_bytes or (lambda: shutil.disk_usage(store.path.parent).free)
        with self.db:
            self.db.execute('''CREATE TABLE request_index (
 namespace TEXT NOT NULL REFERENCES namespaces(digest), request_hash TEXT NOT NULL,
 descriptor BLOB NOT NULL, state TEXT NOT NULL,
 input_sequence INTEGER REFERENCES records(sequence),
 output_sequence INTEGER REFERENCES records(sequence),
 PRIMARY KEY(namespace, request_hash))''')

    def _namespace(self, namespace):
        row = self.db.execute('SELECT descriptor FROM namespaces WHERE digest=?', (namespace,)).fetchone()
        require(row is not None, 'declared namespace required')
        raw = bytes(row[0])
        require(digest(raw) == namespace, 'namespace hash/operand integrity')
        value = json.loads(raw)
        require(canonical(value) == raw, 'canonical namespace operands')
        return value

    def _request(self, namespace, descriptor):
        ns = self._namespace(namespace)
        require(type(descriptor) is dict, 'full request object')
        require(set(('route', 'purpose', 'precision', 'kind', 'context', 'arguments', 'settings')) <= set(descriptor),
                'complete request-interface fields')
        require(descriptor['route'] in ('A24', 'A48', 'B50') and descriptor['route'] == ns['route'], 'route join')
        settings = ns['settings']
        require(type(descriptor['purpose']) is str and bool(descriptor['purpose']) and
                descriptor['purpose'] == settings.get('purpose'), 'purpose join')
        require(type(descriptor['precision']) is int and descriptor['precision'] > 0 and
                descriptor['precision'] == settings.get('precision'), 'actual precision join')
        require(type(descriptor['kind']) is str and bool(descriptor['kind']), 'request kind')
        require(all(type(descriptor[k]) is dict and bool(descriptor[k])
                    for k in ('context', 'arguments', 'settings')), 'explicit operands and settings')
        raw = canonical(descriptor)
        require(len(raw) <= self.maximum, 'request record capacity')
        return digest(raw), raw

    def _space(self, size):
        free = self.free_bytes()
        require(type(free) is int and free >= 0, 'valid disk-space observation')
        # Reserve space for this record plus rollback-journal/codec overhead.
        needed = self.reserve + 2*size + 131072
        require(free >= needed, 'disk reserve or bounded record space unavailable')
        return {'freeBytes': free, 'reserveBytes': self.reserve, 'nextRawBytes': size,
                'requiredFreeBytes': needed}

    def _put(self, namespace, name, value):
        raw = canonical(value)
        require(len(raw) <= self.maximum, 'evidence record capacity')
        self._space(len(raw))
        return self.store.put(namespace, name, value)

    def begin(self, namespace, descriptor):
        key, raw = self._request(namespace, descriptor)
        existing = self.db.execute('SELECT descriptor FROM request_index WHERE namespace=? AND request_hash=?',
                                   (namespace, key)).fetchone()
        if existing is not None:
            require(bytes(existing[0]) == raw, 'request hash collision')
            raise RequestRefusal('request already reserved; lookup complete returns, never rerun pending work')
        self._space(len(raw))
        # Reserve before the input write so a failed write cannot trigger replay.
        with self.db:
            self.db.execute('INSERT INTO request_index VALUES (?,?,?, ?,NULL,NULL)',
                            (namespace, key, raw, 'PREPARING'))
        receipt = self._put(namespace, 'request/'+key+'/input', descriptor)
        with self.db:
            self.db.execute('UPDATE request_index SET state=?,input_sequence=? WHERE namespace=? AND request_hash=? AND state=?',
                            ('PENDING', receipt['sequence'], namespace, key, 'PREPARING'))
        return {'namespace': namespace, 'requestHash': key, 'inputReceipt': receipt,
                'descriptor': json.loads(raw)}

    def complete(self, ticket, value):
        namespace = ticket['namespace']
        key, raw = self._request(namespace, ticket['descriptor'])
        require(key == ticket['requestHash'], 'ticket descriptor join')
        row = self.db.execute('SELECT descriptor,state,input_sequence FROM request_index WHERE namespace=? AND request_hash=?',
                              (namespace, key)).fetchone()
        require(row is not None and bytes(row[0]) == raw and row[1] == 'PENDING', 'one pending exact request')
        original, receipt = self.read_record(row[2])
        require(canonical(original) == raw and receipt == ticket['inputReceipt'], 'full input and receipt join')
        out = self._put(namespace, 'request/'+key+'/return', value)
        with self.db:
            changed = self.db.execute('UPDATE request_index SET state=?,output_sequence=? WHERE namespace=? AND request_hash=? AND state=?',
                                      ('COMPLETE', out['sequence'], namespace, key, 'PENDING')).rowcount
            require(changed == 1, 'single request completion')
        return out

    def lookup(self, namespace, descriptor):
        key, raw = self._request(namespace, descriptor)
        row = self.db.execute('SELECT descriptor,state,input_sequence,output_sequence FROM request_index WHERE namespace=? AND request_hash=?',
                              (namespace, key)).fetchone()
        if row is None:
            return None
        require(bytes(row[0]) == raw, 'request hash collision or changed full operands')
        require(row[1] == 'COMPLETE', 'incomplete request preserved; no automatic replay')
        inp, ir = self.read_record(row[2])
        require(canonical(inp) == raw and ir['namespace'] == namespace and ir['record'] == 'request/'+key+'/input',
                'full original input join')
        value, receipt = self.read_record(row[3])
        require(receipt['namespace'] == namespace and receipt['record'] == 'request/'+key+'/return', 'exact return link')
        return {'value': value, 'receipt': receipt, 'inputReceipt': ir, 'descriptor': json.loads(raw)}

    def read_record(self, sequence):
        require(type(sequence) is int and sequence >= 0, 'exact nonnegative sequence')
        # Any SQLite mutation invalidates payload cache, including writes through
        # another connection. Immutable indices do not authorize stale DB bytes.
        epoch = (self.db.total_changes, self.db.execute('PRAGMA data_version').fetchone()[0])
        if epoch != self.cache_epoch:
            self.cache.clear(); self.cache_size = 0; self.cache_epoch = epoch
        row = self.db.execute('SELECT namespace,name,digest,raw_bytes,compressed_bytes,chunk_count,previous,chain,codec FROM records WHERE sequence=?',
                              (sequence,)).fetchone()
        require(row is not None, 'existing full evidence record')
        ns, name, h, size, compressed, count, previous, chain, codec = row
        require(type(size) is int and 0 <= size <= self.maximum and codec == 'zlib-stream-v1', 'bounded known record codec')
        receipt = {'sequence': sequence, 'namespace': ns, 'record': name, 'sha256': h,
                   'bytes': size, 'compressedBytes': compressed, 'chunks': count,
                   'previous': previous, 'codec': codec}
        require(digest(canonical(receipt)) == chain, 'record receipt hash')
        old = self.db.execute('SELECT chain FROM records WHERE sequence=?', (sequence-1,)).fetchone() if sequence else None
        require(previous == (old[0] if old else '0'*64) and (sequence == 0 or old is not None), 'immediate chain link')
        # Cache contains canonical encoded bytes only. Every returned object is fresh.
        cache_key = (ns, sequence, h, chain)
        raw = self.cache.get(cache_key)
        if raw is None:
            decoder = zlib.decompressobj()
            parts = []; length = packed = number = 0
            for ordinal, blob in self.db.execute('SELECT ordinal,payload FROM chunks WHERE sequence=? ORDER BY ordinal', (sequence,)):
                require(ordinal == number, 'complete ordered chunks')
                piece = decoder.decompress(blob, self.maximum-length+1)
                length += len(piece); packed += len(blob); number += 1
                require(length <= self.maximum and not decoder.unconsumed_tail, 'bounded decompression')
                require(not decoder.unused_data, 'no trailing compressed data')
                parts.append(piece)
            raw = b''.join(parts)
            require(decoder.eof and length == size and packed == compressed and number == count, 'complete compressed record')
            require(digest(raw) == h, 'full evidence payload hash')
            if len(raw) <= self.cache_limit:
                while self.cache and self.cache_size+len(raw) > self.cache_limit:
                    _, evicted = self.cache.popitem(last=False); self.cache_size -= len(evicted)
                self.cache[cache_key] = raw; self.cache_size += len(raw)
        else:
            self.cache.move_to_end(cache_key)
        value = json.loads(raw)
        require(canonical(value) == raw, 'canonical full evidence value')
        return value, {**receipt, 'chainSha256': chain}

    def put_panel(self, namespace, name, panel):
        require(type(panel) is dict, 'full panel record')
        columns = ('points', 'weights', 'jacobians', 'values', 'operands')
        require(all(type(panel.get(k)) is list for k in columns), 'full node columns')
        count = len(panel['points'])
        require(0 < count <= 48 and all(len(panel[k]) == count for k in columns), 'complete bounded panel; no zip truncation')
        return self._put(namespace, name, panel)
