"""Synthetic JSON persistence/fault tests only; never open native science banks."""
import hashlib,importlib.util,json,math,os,random,sqlite3,tempfile,unittest
from pathlib import Path
P=Path(__file__).with_name('S11c_d_defect_packet_evidence_store.py')
s=importlib.util.spec_from_file_location('proposed_evidence_store',P);E=importlib.util.module_from_spec(s);s.loader.exec_module(E)
class StoreTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.path=Path(self.tmp.name)/'synthetic.sqlite';self.store=E.EvidenceStore(self.path)
        self.descriptor={'route':'A24','settings':{'precision':17,'manufactured':True,'source':'synthetic-only'}};self.ns=self.store.namespace(self.descriptor)
    def tearDown(self):self.store.close();self.tmp.cleanup()
    def read(self,r):
        with E.EvidenceReader(self.path) as reader:return reader.read_json(r)
    def test_json_types(self):
        x={'x':[None,True,False,0,-3,2**200,-0.0,1.25,'λ\x00\n'], 'empty':{}}
        r=self.store.put(self.ns,'types',x);self.assertEqual(self.read(r),x);self.assertEqual(r['sha256'],hashlib.sha256(E.canonical(x)).hexdigest());self.assertEqual(r['bytes'],len(E.canonical(x)))
    def test_exact_opaque_tuple_encoding(self):
        x={'mpf':[0,'7348327234732834728347283472834728347',-391,123],'decimal':'0.00000000001234567890123456789'}
        r=self.store.put(self.ns,'opaque',x);self.assertEqual(self.read(r),x)
    def test_multiple_routes_not_shared(self):
        x={'values':['2','3']};namespaces=[self.ns,self.store.namespace({'route':'A48','settings':self.descriptor['settings']}),self.store.namespace({'route':'B50','settings':self.descriptor['settings']})]
        receipts=[self.store.put(n,'same-name',x) for n in namespaces]
        self.assertEqual(len({r['namespace'] for r in receipts}),3);self.assertEqual(len({r['sequence'] for r in receipts}),3)
        self.assertEqual(self.store.db.execute('SELECT COUNT(DISTINCT sequence) FROM chunks').fetchone()[0],3)
    def test_namespace_settings_identity(self):
        other=self.store.namespace({'route':'A24','settings':{'precision':18,'manufactured':True,'source':'synthetic-only'}})
        self.assertNotEqual(self.ns,other);self.assertEqual(self.ns,self.store.namespace(self.descriptor))
    def test_unknown_namespace_refused(self):
        with self.assertRaises(ValueError):self.store.put('0'*64,'bad',{})
    def test_invalid_route_refused(self):
        with self.assertRaises(ValueError):self.store.namespace({'route':'mixed-A-B','settings':{}})
    def test_settings_required(self):
        with self.assertRaises(ValueError):self.store.namespace({'route':'A24'})
    def test_existing_file_not_modified(self):
        self.store.put(self.ns,'first',{'saved':True});before=self.path.read_bytes()
        with self.assertRaises(FileExistsError):E.EvidenceStore(self.path)
        self.assertEqual(before,self.path.read_bytes())
    def test_no_overwrite_even_equal(self):
        r=self.store.put(self.ns,'first',{'a':1})
        with self.assertRaises(ValueError):self.store.put(self.ns,'first',{'a':1})
        self.assertEqual(self.read(r),{'a':1});self.assertEqual(self.store.sequence,1)
    def test_no_overwrite_changed(self):
        self.store.put(self.ns,'first',{'a':1})
        with self.assertRaises(ValueError):self.store.put(self.ns,'first',{'a':2})
    def test_empty_payload_values(self):
        for i,x in enumerate([None,{},[],False,0,'']):self.assertEqual(self.read(self.store.put(self.ns,str(i),x)),x)
    def test_stream_spans_many_chunks(self):
        rng=random.Random(918);x={'opaque':[str(rng.getrandbits(128)) for _ in range(16000)]}
        r=self.store.put(self.ns,'many-chunks',x);self.assertGreater(r['chunks'],2)
        with E.EvidenceReader(self.path) as reader:
            pieces=list(reader.iter_bytes(r));self.assertTrue(all(len(v)<=E.BLOCK for v in pieces));self.assertEqual(b''.join(pieces),E.canonical(x))
    def test_repeated_large_content(self):
        x={'full':('abcλ\n'*500000)};r=self.store.put(self.ns,'compressible',x)
        self.assertEqual(self.read(r),x);self.assertLess(r['compressedBytes'],r['bytes'])
    def test_chain(self):
        records=[self.store.put(self.ns,'record/'+str(i),{'index':i,'previousResult':None}) for i in range(10)]
        self.assertEqual(records[0]['previous'],E.ZERO)
        for a,b in zip(records,records[1:]):self.assertEqual(b['previous'],a['chainSha256'])
        with E.EvidenceReader(self.path) as reader:self.assertEqual(reader.audit(),{'records':10,'lastChain':records[-1]['chainSha256'],'allBytesChecked':True,'scientificAcceptance':False})
    def test_external_head_detects_truncated_suffix(self):
        first=self.store.put(self.ns,'first',{});last=self.store.put(self.ns,'last',{})
        self.store.db.execute('DELETE FROM chunks WHERE sequence=?',(last['sequence'],));self.store.db.execute('DELETE FROM records WHERE sequence=?',(last['sequence'],));self.store.db.commit()
        with E.EvidenceReader(self.path) as reader:
            with self.assertRaises(ValueError):reader.audit(2,last['chainSha256'])
    def test_external_completion_receipt(self):
        r=self.store.put(self.ns,'first',{'synthetic':True})
        with E.EvidenceReader(self.path) as reader:self.assertEqual(reader.audit(1,r['chainSha256'])['records'],1)
    def test_incomplete_external_receipt_refused(self):
        with E.EvidenceReader(self.path) as reader:
            with self.assertRaises(ValueError):reader.audit(expected_records=0)
    def test_empty_audit(self):
        with E.EvidenceReader(self.path) as reader:self.assertEqual(reader.audit()['records'],0)
    def test_transaction_rollback_preserves_prefix(self):
        first=self.store.put(self.ns,'first',{'saved':'original'});chain=self.store.previous
        # Sorted key z is reached only after compressing the large earlier value.
        with self.assertRaises(ValueError):self.store.put(self.ns,'failing',{'a':'x'*500000,'z':float('nan')})
        self.assertEqual(self.store.previous,chain);self.assertEqual(self.store.sequence,1)
        self.assertEqual(self.store.db.execute('SELECT COUNT(*) FROM records').fetchone()[0],1)
        second=self.store.put(self.ns,'next-independent-tooling-record',{'afterFailure':'no scientific retry'})
        self.assertEqual(second['sequence'],1);self.assertEqual(second['previous'],first['chainSha256'])
        with E.EvidenceReader(self.path) as reader:self.assertEqual(reader.audit()['records'],2)
    def test_closed_writer_refuses(self):
        self.store.close()
        with self.assertRaises(ValueError):self.store.put(self.ns,'closed',{})
    def test_read_only_connection(self):
        with E.EvidenceReader(self.path) as reader:
            with self.assertRaises(sqlite3.OperationalError):reader.db.execute('DELETE FROM namespaces')
    def test_full_sync(self):self.assertEqual(self.store.db.execute('PRAGMA synchronous').fetchone()[0],2)
    def test_no_json_key_coercion(self):
        for x in [{3:'bad'},{'nested':{True:'bad'}}]:
            with self.assertRaises(ValueError):self.store.put(self.ns,'badkeys',x)
    def test_no_tuple_coercion(self):
        with self.assertRaises(ValueError):self.store.put(self.ns,'badtuple',{'a':(1,2)})
    def test_cycle_refused(self):
        x=[];x.append(x)
        with self.assertRaises(ValueError):self.store.put(self.ns,'cycle',x)
    def test_repeated_container_allowed(self):
        child=[1,2];x=[child,child];r=self.store.put(self.ns,'repeated-object',x);self.assertEqual(self.read(r),x)
    def test_nonfinite_rejected(self):
        for value in [float('nan'),float('inf'),-float('inf')]:
            with self.assertRaises(ValueError):self.store.put(self.ns,'nonfinite',value)
    def test_unknown_objects_refused(self):
        with self.assertRaises(ValueError):self.store.put(self.ns,'object',object())
    def test_receipt_change_refused(self):
        r=self.store.put(self.ns,'item',{'a':1});r['bytes']+=1
        with E.EvidenceReader(self.path) as reader:
            with self.assertRaises(ValueError):reader.read_bytes(r)
    def test_compressed_corruption_detected(self):
        r=self.store.put(self.ns,'item',{'a':'some value'});self.store.db.execute('UPDATE chunks SET payload=?',(b'broken',));self.store.db.commit()
        with E.EvidenceReader(self.path) as reader:
            with self.assertRaises((ValueError,zlib_error())):reader.read_bytes(r)
    def test_deleted_chunk_detected(self):
        r=self.store.put(self.ns,'item',{'a':1});self.store.db.execute('DELETE FROM chunks');self.store.db.commit()
        with E.EvidenceReader(self.path) as reader:
            with self.assertRaises(ValueError):reader.read_bytes(r)
    def test_trailing_stream_refused(self):
        import zlib
        r=self.store.put(self.ns,'item',{'a':1});old=self.store.db.execute('SELECT payload FROM chunks').fetchone()[0];self.store.db.execute('UPDATE chunks SET payload=?',(old+zlib.compress(b'extra'),));self.store.db.commit()
        with E.EvidenceReader(self.path) as reader:
            with self.assertRaises(ValueError):reader.read_bytes(r)
    def test_reordered_chunks_refused(self):
        rng=random.Random(1);r=self.store.put(self.ns,'item',[str(rng.getrandbits(128)) for _ in range(9000)])
        self.store.db.execute('UPDATE chunks SET ordinal=ordinal+10');self.store.db.commit()
        with E.EvidenceReader(self.path) as reader:
            with self.assertRaises(ValueError):reader.read_bytes(r)
    def test_record_deletion_chain_detected(self):
        first=self.store.put(self.ns,'first',{});self.store.put(self.ns,'second',{})
        self.store.db.execute('DELETE FROM chunks WHERE sequence=?',(first['sequence'],));self.store.db.execute('DELETE FROM records WHERE sequence=?',(first['sequence'],));self.store.db.commit()
        with E.EvidenceReader(self.path) as reader:
            with self.assertRaises(ValueError):reader.audit()
    def test_envelope_tamper_detected(self):
        r=self.store.put(self.ns,'item',{});self.store.db.execute('UPDATE records SET raw_bytes=raw_bytes+1');self.store.db.commit()
        with E.EvidenceReader(self.path) as reader:
            with self.assertRaises(ValueError):reader.receipt(self.ns,'item')
    def test_namespace_tamper_detected(self):
        r=self.store.put(self.ns,'item',{});self.store.db.execute('UPDATE namespaces SET descriptor=?',(b'{}',));self.store.db.commit()
        with E.EvidenceReader(self.path) as reader:
            with self.assertRaises(ValueError):reader.receipt(self.ns,'item')
    def test_nonexistent_reader_refuses(self):
        p=Path(self.tmp.name)/'absent'
        with self.assertRaises(ValueError):E.EvidenceReader(p)
        self.assertFalse(p.exists())

def zlib_error():
    import zlib
    return zlib.error

if __name__=='__main__':unittest.main()
