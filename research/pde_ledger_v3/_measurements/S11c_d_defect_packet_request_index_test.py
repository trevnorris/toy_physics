"""Manufactured JSON/SQLite only; no scientific payload, import or bank access."""
import copy
import importlib.util
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest.mock import patch
import zlib

HERE=Path(__file__).resolve().parent

def module(name):
    spec=importlib.util.spec_from_file_location(name,HERE/(name+'.py'))
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result

E=module('S11c_d_defect_packet_evidence_store')
Q=module('S11c_d_defect_packet_request_index')


class Requests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='packet-synthetic-requests-')
        self.store=E.EvidenceStore(Path(self.tmp.name)/'synthetic.sqlite')
        self.index=Q.RequestIndex(self.store,cache_bytes=4096,free_bytes=lambda:10**12)
        self.ns={r:self.store.namespace({'route':r,'settings':{'purpose':'baseline','precision':p}})
                 for r,p in [('A24',30),('A48',30),('B50',50)]}

    def tearDown(self):
        self.store.close();self.tmp.cleanup()

    def request(self,route='A48',**changes):
        d={'route':route,'purpose':'baseline','precision':50 if route=='B50' else 30,
           'kind':'synthetic-only','context':{'fixture':'no-native-science','units':'manufactured'},
           'arguments':{'node':{'mpf':[0,'3',-1,2],'decimal':'1.5'}},
           'settings':{'budget':'1/100','window':'manufactured','privateCheckPrecision':50}}
        d.update(changes);return d

    def finish(self,route='A48',value=None,descriptor=None):
        d=descriptor or self.request(route)
        ticket=self.index.begin(self.ns[route],d)
        value=value if value is not None else {'exact':{'mpf':[1,'5',-2,3],'decimal':'-1.25'},'synthetic':True}
        receipt=self.index.complete(ticket,value)
        return d,ticket,value,receipt

    def test_full_round_trip_and_receipt(self):
        d,t,v,r=self.finish();got=self.index.lookup(self.ns['A48'],d)
        self.assertEqual(got['value'],v);self.assertEqual(got['receipt'],r)
        self.assertEqual(got['inputReceipt'],t['inputReceipt']);self.assertEqual(got['descriptor'],d)

    def test_lookup_never_calls_evaluation_or_writer(self):
        d,_,v,_=self.finish()
        with patch.object(self.store,'put',side_effect=AssertionError('no new work')):
            self.assertEqual(self.index.lookup(self.ns['A48'],d)['value'],v)

    def test_pending_is_not_a_cache_miss(self):
        d=self.request();self.index.begin(self.ns['A48'],d)
        with self.assertRaisesRegex(Q.RequestRefusal,'incomplete'):self.index.lookup(self.ns['A48'],d)
        with self.assertRaisesRegex(Q.RequestRefusal,'already reserved'):self.index.begin(self.ns['A48'],d)

    def test_failed_input_write_leaves_reserved_refusal(self):
        d=self.request()
        with patch.object(self.store,'put',side_effect=OSError('manufactured write failure')):
            with self.assertRaises(OSError):self.index.begin(self.ns['A48'],d)
        self.assertEqual(self.store.db.execute('SELECT state FROM request_index').fetchone()[0],'PREPARING')
        with self.assertRaises(Q.RequestRefusal):self.index.lookup(self.ns['A48'],d)

    def test_failed_output_write_preserves_input(self):
        d=self.request();t=self.index.begin(self.ns['A48'],d)
        with patch.object(self.store,'put',side_effect=OSError('manufactured full disk')):
            with self.assertRaises(OSError):self.index.complete(t,{'result':1})
        self.assertEqual(self.index.read_record(t['inputReceipt']['sequence'])[0],d)
        with self.assertRaises(Q.RequestRefusal):self.index.lookup(self.ns['A48'],d)

    def test_failed_link_update_keeps_durable_orphan_and_refuses_replay(self):
        d=self.request();t=self.index.begin(self.ns['A48'],d)
        self.store.db.execute("CREATE TRIGGER refuse_complete BEFORE UPDATE OF state ON request_index WHEN NEW.state='COMPLETE' BEGIN SELECT RAISE(ABORT,'synthetic fault'); END")
        with self.assertRaises(sqlite3.IntegrityError):self.index.complete(t,{'result':'durable-before-link'})
        row=self.store.db.execute("SELECT sequence FROM records WHERE name LIKE '%/return'").fetchone()
        self.assertIsNotNone(row)
        self.assertEqual(self.index.read_record(row[0])[0],{'result':'durable-before-link'})
        with self.assertRaises(Q.RequestRefusal):self.index.lookup(self.ns['A48'],d)
        with self.assertRaises(Q.RequestRefusal):self.index.begin(self.ns['A48'],d)

    def test_duplicate_complete_is_refused(self):
        d,t,_,_=self.finish()
        with self.assertRaises(Q.RequestRefusal):self.index.complete(t,{'changed':1})
        with self.assertRaises(Q.RequestRefusal):self.index.begin(self.ns['A48'],d)

    def test_no_existing_index_reattachment(self):
        with self.assertRaisesRegex(Q.RequestRefusal,'automatic resume'):Q.RequestIndex(self.store)

    def test_different_full_argument_is_miss(self):
        d,_,_,_=self.finish();d['arguments']['node']['mpf'][1]='7'
        self.assertIsNone(self.index.lookup(self.ns['A48'],d))

    def test_error_budget_is_part_of_identity(self):
        d,_,_,_=self.finish();d['settings']['budget']='1/1000'
        self.assertIsNone(self.index.lookup(self.ns['A48'],d))

    def test_original_route_values_are_separate(self):
        receipts=[]
        for r in self.ns:receipts.append(self.finish(r)[3])
        self.assertEqual(len({r['namespace'] for r in receipts}),3)
        self.assertEqual(len({r['sequence'] for r in receipts}),3)

    def test_purpose_separates_control_from_baseline(self):
        d=self.request(purpose='mutant');ns=self.store.namespace({'route':'A48','settings':{'purpose':'mutant','precision':30}})
        t=self.index.begin(ns,d);self.index.complete(t,{'value':'mutant-only'})
        self.assertIsNone(self.index.lookup(self.ns['A48'],self.request()))
        with self.assertRaises(Q.RequestRefusal):self.index.lookup(self.ns['A48'],d)

    def test_route_precision_purpose_require_actual_namespace_join(self):
        for key,value in [('route','B50'),('precision',50),('purpose','check')]:
            with self.subTest(key=key),self.assertRaises(Q.RequestRefusal):
                self.index.begin(self.ns['A48'],self.request(**{key:value}))

    def test_incomplete_descriptor_refused(self):
        for key in self.request():
            d=self.request();del d[key]
            with self.subTest(key=key),self.assertRaises(Q.RequestRefusal):self.index.begin(self.ns['A48'],d)

    def test_empty_operands_refused(self):
        for key in ('context','arguments','settings'):
            with self.subTest(key=key),self.assertRaises(Q.RequestRefusal):self.index.begin(self.ns['A48'],self.request(**{key:{}}))

    def test_hash_alone_does_not_match_changed_descriptor(self):
        d,_,_,_=self.finish()
        with self.store.db:self.store.db.execute('UPDATE request_index SET descriptor=?',(b'{"changed":true}',))
        with self.assertRaisesRegex(Q.RequestRefusal,'full operands'):self.index.lookup(self.ns['A48'],d)

    def test_ticket_tampering_refused(self):
        d=self.request();t=self.index.begin(self.ns['A48'],d);t['descriptor']['settings']['budget']='changed'
        with self.assertRaises(Q.RequestRefusal):self.index.complete(t,{'x':1})

    def test_encoded_values_are_not_evaluated(self):
        v={'text':'__import__("os").system("not executed")','exact':[0,'9007199254740993',-4,54]}
        d,_,_,_=self.finish(value=v)
        self.assertEqual(self.index.lookup(self.ns['A48'],d)['value'],v)

    def test_implicit_tuple_or_key_conversion_refused(self):
        for v in [{'x':(1,2)},{1:'bad'}]:
            with self.subTest(v=v),self.assertRaises(Q.RequestRefusal):Q.canonical(v)

    def test_cycles_refused_repeated_acyclic_allowed(self):
        cycle=[];cycle.append(cycle)
        with self.assertRaisesRegex(Q.RequestRefusal,'cyclic'):Q.canonical(cycle)
        same={'a':[1,2]};self.assertEqual(json.loads(Q.canonical([same,same])),[same,same])

    def test_nonfinite_json_refused(self):
        for n in [float('nan'),float('inf'),-float('inf')]:
            with self.subTest(n=n),self.assertRaises(ValueError):Q.canonical({'n':n})

    def test_chunk_corruption_invalidates_cache_and_is_detected(self):
        d,_,_,r=self.finish();self.index.lookup(self.ns['A48'],d)
        with self.store.db:self.store.db.execute('UPDATE chunks SET payload=? WHERE sequence=?',(b'corrupt',r['sequence']))
        with self.assertRaises((Q.RequestRefusal,zlib.error)):self.index.lookup(self.ns['A48'],d)

    def test_external_sqlite_change_invalidates_cache(self):
        d,_,_,r=self.finish();self.index.lookup(self.ns['A48'],d)
        con=sqlite3.connect(self.store.path)
        with con:con.execute('UPDATE chunks SET payload=? WHERE sequence=?',(b'corrupt',r['sequence']))
        con.close()
        with self.assertRaises((Q.RequestRefusal,zlib.error)):self.index.lookup(self.ns['A48'],d)

    def test_receipt_corruption_refused(self):
        d,_,_,r=self.finish()
        with self.store.db:self.store.db.execute('UPDATE records SET digest=? WHERE sequence=?',('0'*64,r['sequence']))
        with self.assertRaisesRegex(Q.RequestRefusal,'receipt hash'):self.index.lookup(self.ns['A48'],d)

    def test_missing_chunk_refused(self):
        d,_,_,r=self.finish()
        with self.store.db:self.store.db.execute('DELETE FROM chunks WHERE sequence=?',(r['sequence'],))
        with self.assertRaises(Q.RequestRefusal):self.index.lookup(self.ns['A48'],d)

    def test_new_lookup_returns_fresh_mutable_object(self):
        d,_,v,_=self.finish();a=self.index.lookup(self.ns['A48'],d);a['value']['exact']['mpf'][1]='tampered-client'
        self.assertEqual(self.index.lookup(self.ns['A48'],d)['value'],v)

    def test_record_capacity_before_write(self):
        self.index.maximum=10
        with self.assertRaises(Q.RequestRefusal):self.index.begin(self.ns['A48'],self.request())
        self.assertEqual(self.store.db.execute('SELECT count(*) FROM request_index').fetchone()[0],0)

    def test_return_capacity_preserves_pending_input(self):
        d=self.request();t=self.index.begin(self.ns['A48'],d)
        with self.assertRaises(Q.RequestRefusal):self.index.complete(t,{'large':'x'*(self.index.maximum+1)})
        with self.assertRaises(Q.RequestRefusal):self.index.lookup(self.ns['A48'],d)
        self.assertEqual(self.index.read_record(t['inputReceipt']['sequence'])[0],d)

    def test_disk_reserve_before_reservation(self):
        self.index.free_bytes=lambda:0
        with self.assertRaisesRegex(Q.RequestRefusal,'disk reserve'):self.index.begin(self.ns['A48'],self.request())
        self.assertEqual(self.store.db.execute('SELECT count(*) FROM request_index').fetchone()[0],0)

    def test_disk_refusal_preserves_durable_input(self):
        d=self.request();t=self.index.begin(self.ns['A48'],d);self.index.free_bytes=lambda:0
        with self.assertRaises(Q.RequestRefusal):self.index.complete(t,{'x':1})
        self.assertEqual(self.index.read_record(t['inputReceipt']['sequence'])[0],d)

    def test_bounded_panel_completeness(self):
        p={key:[1,2] for key in ('points','weights','jacobians','values','operands')}
        r=self.index.put_panel(self.ns['A48'],'synthetic-panel',p)
        self.assertEqual(self.index.read_record(r['sequence'])[0],p)
        for k in p:
            bad=copy.deepcopy(p);bad[k].pop()
            with self.subTest(k=k),self.assertRaises(Q.RequestRefusal):self.index.put_panel(self.ns['A48'],'bad-'+k,bad)

    def test_no_more_than_48_nodes_in_one_record(self):
        for n in [0,49]:
            p={key:[1]*n for key in ('points','weights','jacobians','values','operands')}
            with self.subTest(n=n),self.assertRaises(Q.RequestRefusal):self.index.put_panel(self.ns['A48'],'bad-panel',p)

    def test_cache_limit_and_namespace_keys(self):
        self.index.cache_limit=700
        for i in range(4):
            d=self.request();d['arguments']['tag']=i
            self.finish(descriptor=d);self.index.lookup(self.ns['A48'],d)
            self.assertLessEqual(self.index.cache_size,700)
        self.assertEqual(self.index.cache_size,sum(len(b) for b in self.index.cache.values()))
        self.assertTrue(all(k[0]==self.ns['A48'] for k in self.index.cache))

    def test_unknown_record_sequence_refused(self):
        for n in [-1,False,1.2,999]:
            with self.subTest(n=n),self.assertRaises(Q.RequestRefusal):self.index.read_record(n)


if __name__=='__main__':unittest.main(verbosity=2)
