"""Proposed lossless JSON journal for a future guarded packet evaluator.

No numerical imports, computation, result cache, automatic resume, old-bank
migration or scientific acceptance. Tests use manufactured JSON only. Native
use requires the future evaluator's independently assessed build and guard.
"""
import hashlib
import json
import sqlite3
import zlib
from pathlib import Path

BLOCK=65536
ZERO='0'*64

def require(value,message):
    if value is not True:raise ValueError(message)


def json_value(value):
    """Refuse implicit key/type conversions; values already use the exact codec."""
    # The encoder detects container cycles. This walk uses enter/leave markers
    # so a repeated non-cyclic container is allowed without recursion.
    active=set();stack=[(value,False)]
    while stack:
        v,leaving=stack.pop()
        if type(v) in (dict,list):
            if leaving:active.remove(id(v));continue
            require(id(v) not in active,'cyclic JSON container');active.add(id(v));stack.append((v,True))
            if type(v) is dict:
                require(all(type(k) is str for k in v),'exact string JSON keys')
                stack.extend((x,False) for x in v.values())
            else:stack.extend((x,False) for x in v)
        else:require(type(v) in (type(None),str,int,float,bool),'pre-encoded JSON types only')


def canonical(value):
    json_value(value)
    return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode('utf-8')


def chain_digest(value):return hashlib.sha256(canonical(value)).hexdigest()


class EvidenceStore:
    """One writer, immutable route namespaces, one FULL transaction per record.

    Compression changes stored bytes only. The digest/length refer to the entire
    original canonical JSON stream. Every occurrence remains a separate record;
    no evaluated arrays or values are shared between numerical routes.
    """
    def __init__(self,path):
        self.path=Path(path)
        # Exclusive filesystem creation: an existing failed/completed journal
        # is never opened for writes, replaced, resumed or migrated.
        with self.path.open('xb'):pass
        self.db=sqlite3.connect(str(self.path));self.closed=False
        self.db.execute('PRAGMA journal_mode=DELETE')
        self.db.execute('PRAGMA synchronous=FULL')
        self.db.execute('PRAGMA foreign_keys=ON')
        self.db.execute('PRAGMA cache_size=-4096')
        self.db.execute('PRAGMA temp_store=FILE')
        require(self.db.execute('PRAGMA synchronous').fetchone()[0]==2,'FULL synchronous journal')
        require(self.db.execute('PRAGMA journal_mode').fetchone()[0]=='delete','DELETE journal')
        with self.db:
            self.db.executescript('''
CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE namespaces (digest TEXT PRIMARY KEY, descriptor BLOB NOT NULL UNIQUE);
CREATE TABLE records (
 sequence INTEGER PRIMARY KEY, namespace TEXT NOT NULL REFERENCES namespaces(digest),
 name TEXT NOT NULL, digest TEXT, raw_bytes INTEGER, compressed_bytes INTEGER,
 chunk_count INTEGER, previous TEXT NOT NULL, chain TEXT, codec TEXT NOT NULL,
 UNIQUE(namespace,name));
CREATE TABLE chunks (
 sequence INTEGER NOT NULL REFERENCES records(sequence), ordinal INTEGER NOT NULL,
 payload BLOB NOT NULL, PRIMARY KEY(sequence,ordinal));
''')
            self.db.execute('INSERT INTO metadata VALUES (?,?)',('format','packet-json-zlib-chunks-v1'))
        self.sequence=0;self.previous=ZERO

    def namespace(self,descriptor):
        require(not self.closed,'writer open')
        require(type(descriptor) is dict and descriptor.get('route') in ('A24','A48','B50','mathematical-inputs'),'explicit independent route namespace')
        require(type(descriptor.get('settings')) is dict,'explicit complete settings object; semantic adequacy needs build review')
        data=canonical(descriptor);digest=hashlib.sha256(data).hexdigest()
        existing=self.db.execute('SELECT descriptor FROM namespaces WHERE digest=?',(digest,)).fetchone()
        if existing is not None:require(existing[0]==data,'namespace hash collision; refuse')
        else:
            with self.db:self.db.execute('INSERT INTO namespaces VALUES (?,?)',(digest,data))
        return digest

    def put(self,namespace,name,value):
        require(not self.closed and type(name) is str and bool(name),'open writer and explicit record name')
        require(self.db.execute('SELECT 1 FROM namespaces WHERE digest=?',(namespace,)).fetchone() is not None,'declared exact namespace')
        require(self.db.execute('SELECT 1 FROM records WHERE namespace=? AND name=?',(namespace,name)).fetchone() is None,'immutable record name; no overwrite/retry')
        json_value(value)
        encoder=json.JSONEncoder(sort_keys=True,separators=(',',':'),allow_nan=False)
        compressor=zlib.compressobj(level=6);digest=hashlib.sha256();size=0;compressed=0;chunks=0;buffer=bytearray()
        sequence=self.sequence;previous=self.previous
        def append(data):
            nonlocal compressed,chunks
            buffer.extend(data)
            while len(buffer)>=BLOCK:
                block=bytes(buffer[:BLOCK]);del buffer[:BLOCK]
                self.db.execute('INSERT INTO chunks VALUES (?,?,?)',(sequence,chunks,block));chunks+=1;compressed+=len(block)
        with self.db:
            self.db.execute('INSERT INTO records(sequence,namespace,name,previous,codec) VALUES (?,?,?,?,?)',(sequence,namespace,name,previous,'zlib-stream-v1'))
            for piece in encoder.iterencode(value):
                # json encoder may materialize one large quoted string. This
                # is not a proof of bounded worker RSS; future record sizing
                # and containment must account for the original live object.
                for start in range(0,len(piece),BLOCK):
                    raw=piece[start:start+BLOCK].encode('utf-8');digest.update(raw);size+=len(raw);append(compressor.compress(raw))
            append(compressor.flush())
            if buffer:
                self.db.execute('INSERT INTO chunks VALUES (?,?,?)',(sequence,chunks,bytes(buffer)));chunks+=1;compressed+=len(buffer)
            receipt={'sequence':sequence,'namespace':namespace,'record':name,'sha256':digest.hexdigest(),'bytes':size,'compressedBytes':compressed,'chunks':chunks,'previous':previous,'codec':'zlib-stream-v1'}
            chain=chain_digest(receipt)
            self.db.execute('UPDATE records SET digest=?,raw_bytes=?,compressed_bytes=?,chunk_count=?,chain=? WHERE sequence=?',(receipt['sha256'],size,compressed,chunks,chain,sequence))
        # These advance only after transaction commit succeeds.
        self.sequence+=1;self.previous=chain
        return {**receipt,'chainSha256':chain}

    def close(self):
        if not self.closed:self.db.close();self.closed=True

    def __enter__(self):return self
    def __exit__(self,*args):self.close()


class EvidenceReader:
    """Read-only full-byte reconstruction. This does not restore scientific types."""
    def __init__(self,path):
        self.path=Path(path);require(self.path.is_file(),'existing journal required')
        self.db=sqlite3.connect(self.path.resolve().as_uri()+'?mode=ro',uri=True)
        self.db.execute('PRAGMA query_only=ON')
        require(self.db.execute("SELECT value FROM metadata WHERE key='format'").fetchone()==('packet-json-zlib-chunks-v1',),'exact journal format')

    def receipt(self,namespace,name):
        row=self.db.execute('SELECT sequence,namespace,name,digest,raw_bytes,compressed_bytes,chunk_count,previous,chain,codec FROM records WHERE namespace=? AND name=?',(namespace,name)).fetchone()
        require(row is not None and all(x is not None for x in row),'complete published record')
        seq,ns,n,digest,size,compressed,count,prev,chain,codec=row
        r={'sequence':seq,'namespace':ns,'record':n,'sha256':digest,'bytes':size,'compressedBytes':compressed,'chunks':count,'previous':prev,'codec':codec}
        require(codec=='zlib-stream-v1' and type(size) is int and size>=0 and type(count) is int and count>0,'valid record envelope')
        require(chain_digest(r)==chain,'record chain envelope matches')
        descriptor=self.db.execute('SELECT descriptor FROM namespaces WHERE digest=?',(namespace,)).fetchone()
        require(descriptor is not None and hashlib.sha256(descriptor[0]).hexdigest()==namespace,'original namespace bytes')
        return {**r,'chainSha256':chain}

    def iter_bytes(self,receipt):
        current=self.receipt(receipt['namespace'],receipt['record']);require(current==receipt,'complete actual receipt identity')
        inflater=zlib.decompressobj();h=hashlib.sha256();size=0;compressed=0;count=0
        for ordinal,payload in self.db.execute('SELECT ordinal,payload FROM chunks WHERE sequence=? ORDER BY ordinal',(receipt['sequence'],)):
            require(ordinal==count and type(payload) is bytes and 0<len(payload)<=BLOCK,'complete ordered chunk coverage');count+=1;compressed+=len(payload)
            pending=payload
            while pending:
                data=inflater.decompress(pending,BLOCK);pending=inflater.unconsumed_tail
                require(not inflater.unused_data,'no trailing compressed stream')
                size+=len(data);require(size<=receipt['bytes'],'decompressed length exceeds exact receipt')
                h.update(data)
                if data:yield data
        # A correct zlib stream emits all data above. flush must not conceal a
        # truncated stream or manufacture bytes outside the declared length.
        tail=inflater.flush();size+=len(tail);h.update(tail)
        require(inflater.eof and not inflater.unused_data and not inflater.unconsumed_tail,'complete single compressed stream')
        require((size,compressed,count,h.hexdigest())==(receipt['bytes'],receipt['compressedBytes'],receipt['chunks'],receipt['sha256']),'full payload byte identity')
        if tail:yield tail

    def read_bytes(self,receipt):return b''.join(self.iter_bytes(receipt))
    def read_json(self,receipt):return json.loads(self.read_bytes(receipt))

    def audit(self,expected_records=None,expected_head=None):
        require((expected_records is None)==(expected_head is None),'provide both external expected count and head')
        if expected_records is not None:require(type(expected_records) is int and expected_records>=0 and type(expected_head) is str,'exact external completion receipt')
        previous=ZERO;expected=0
        # Deliberately all records, not only complete-looking selected outputs.
        for ns,name in self.db.execute('SELECT namespace,name FROM records ORDER BY sequence'):
            r=self.receipt(ns,name);require(r['sequence']==expected and r['previous']==previous,'complete contiguous evidence chain')
            for _ in self.iter_bytes(r):pass
            previous=r['chainSha256'];expected+=1
        orphan=self.db.execute('SELECT COUNT(*) FROM chunks LEFT JOIN records USING(sequence) WHERE records.sequence IS NULL').fetchone()[0]
        require(orphan==0,'no orphan chunks')
        if expected_records is not None:require(expected==expected_records and previous==expected_head,'complete externally pinned chain, including final record')
        return {'records':expected,'lastChain':previous,'allBytesChecked':True,'scientificAcceptance':False}

    def close(self):self.db.close()
    def __enter__(self):return self
    def __exit__(self,*args):self.close()
