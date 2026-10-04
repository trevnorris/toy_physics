"""Manufactured JSON/SQLite and AST tests only. No scientific restoration."""
from copy import deepcopy
import ast
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest

import S11c_d_defect_packet_contracted_prefix_recovery as P
from S11c_d_defect_packet_evidence_store import EvidenceStore
from S11c_d_defect_packet_request_index import RequestIndex, canonical

M = Path(__file__).resolve().parent


class Fixture(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name); self.original = self.root/'original.sqlite'
        self.destination = self.root/'branch.sqlite'
        store = EvidenceStore(self.original); index = RequestIndex(store, reserve_bytes=0)
        self.ns = store.namespace({'route':'A24', 'settings':{'purpose':'baseline','precision':30}})
        context = {'meaning':'manufactured, no scientific values'}
        def desc(kind, order=24):
            return {'route':'A24','purpose':'baseline','precision':30,'kind':kind,
                    'context':context,'arguments':{'a':'-1','b':'0'},
                    'settings':{'selectedRule':'A'+str(order),'parent':'manufactured'}}
        self.pending_desc = desc('outer-panel')
        parent = index.begin(self.ns, self.pending_desc)
        self.panels=[]
        for order in (24,48):
            d=desc('inner-panel',order); ticket=index.begin(self.ns,d)
            out=index.complete(ticket,{'components':{'synthetic':{'value':'uninterpreted','error':'uninterpreted'}},'order':order})
            self.panels.append({'namespace':self.ns,'descriptor':d,'inputReceipt':ticket['inputReceipt'],'returnReceipt':out})
        with store.db:
            store.db.execute('CREATE TABLE contracted_leaves (job TEXT, ordinal INTEGER, active INTEGER, receipt BLOB NOT NULL, PRIMARY KEY(job,ordinal))')
        records,head=store.sequence,store.previous;store.close()
        self.contract={'databaseSha256':P.sha(self.original),'databaseBytes':self.original.stat().st_size,
                       'records':records,'head':head,'completedRequests':2,'pending':parent,'panels':self.panels,
                       'reserveBytes':0,'maxRecordBytes':8*1024**2,'cacheBytes':0}
        self.old_bytes=self.original.read_bytes()

    def branch(self):
        store=P.ExplicitCloneStore(self.original,self.destination,self.contract)
        self.addCleanup(store.close)
        index=P.ExplicitParentIndex(store,self.contract,free_bytes=lambda:10**12)
        return store,index

    def test_full_audit(self):
        got=P.audit_prefix(self.original,self.contract)
        self.assertEqual(got['completeRequests'],2); self.assertFalse(got['scientificAcceptance'])

    def test_branch_is_byte_identical_before_new_records(self):
        store,_=self.branch();self.assertEqual(self.destination.read_bytes(),self.old_bytes)
        self.assertEqual(store.sequence,self.contract['records'])
        self.assertEqual(store.previous,self.contract['head'])

    def test_original_unchanged_after_append(self):
        store,index=self.branch();index._put(self.ns,'new-evidence',{'notScience':True})
        self.assertEqual(self.original.read_bytes(),self.old_bytes)
        self.assertEqual(store.prefix_preserved()['parentState'],'PENDING')

    def test_ordinary_lookup_still_refuses_pending(self):
        _,index=self.branch()
        with self.assertRaisesRegex(ValueError,'incomplete request'):index.lookup(self.ns,self.pending_desc)

    def test_parent_exact_claim_and_one_way_completion(self):
        store,index=self.branch();ticket=index.claim_parent(self.ns,self.pending_desc)
        self.assertEqual(ticket,self.contract['pending'])
        index.complete(ticket,{'newUnfinishedReturn':'synthetic'})
        self.assertEqual(store.prefix_preserved()['parentState'],'COMPLETE')
        with self.assertRaises(ValueError):index.claim_parent(self.ns,self.pending_desc)
        with self.assertRaises(ValueError):index.complete(ticket,{'new':'repeat'})

    def test_claim_twice_refuses_without_reexecution(self):
        _,index=self.branch();index.claim_parent(self.ns,self.pending_desc)
        with self.assertRaises(ValueError):index.claim_parent(self.ns,self.pending_desc)

    def test_parent_completion_without_claim_refuses(self):
        _,index=self.branch()
        with self.assertRaisesRegex(ValueError,'claim'):index.complete(self.contract['pending'],{})

    def test_completed_lookup_returns_exact_operands_without_callback(self):
        _,index=self.branch()
        for witness in self.panels:
            got=index.lookup(witness['namespace'],witness['descriptor'])
            self.assertEqual(got['receipt'],witness['returnReceipt'])
            self.assertEqual(got['inputReceipt'],witness['inputReceipt'])
        self.assertFalse(index.claimed)

    def test_completed_request_cannot_be_started_again(self):
        _,index=self.branch()
        with self.assertRaises(ValueError):index.begin(self.ns,self.panels[0]['descriptor'])

    def test_new_request_works_normally(self):
        _,index=self.branch();d=deepcopy(self.pending_desc);d['arguments']['a']='new'
        self.assertIsNone(index.claim_parent(self.ns,d));self.assertIsNone(index.lookup(self.ns,d))
        ticket=index.begin(self.ns,d);index.complete(ticket,{'new':True})
        self.assertEqual(index.lookup(self.ns,d)['value'],{'new':True})

    def test_other_pending_request_has_no_claim_permission(self):
        _,index=self.branch();d=deepcopy(self.pending_desc);d['arguments']['a']='new'
        index.begin(self.ns,d);self.assertIsNone(index.claim_parent(self.ns,d))
        with self.assertRaises(ValueError):index.lookup(self.ns,d)

    def test_full_argument_change_cannot_claim_parent(self):
        _,index=self.branch();d=deepcopy(self.pending_desc);d['context']['meaning']='different'
        self.assertIsNone(index.claim_parent(self.ns,d));self.assertFalse(index.claimed)

    def test_precision_change_refuses(self):
        _,index=self.branch();d=deepcopy(self.pending_desc);d['precision']=50
        with self.assertRaises(ValueError):index.claim_parent(self.ns,d)

    def test_same_original_destination_refuses(self):
        with self.assertRaises(ValueError):P.ExplicitCloneStore(self.original,self.original,self.contract)
        self.assertEqual(self.original.read_bytes(),self.old_bytes)

    def test_existing_destination_refuses(self):
        self.destination.write_bytes(b'preserve')
        with self.assertRaises(ValueError):self.branch()
        self.assertEqual(self.destination.read_bytes(),b'preserve')

    def test_stale_hash_refuses(self):
        self.contract['databaseSha256']='0'*64
        with self.assertRaises(ValueError):self.branch()
        self.assertFalse(self.destination.exists())

    def test_sidecar_refuses(self):
        Path(str(self.original)+'-journal').write_bytes(b'pending')
        with self.assertRaises(ValueError):self.branch()

    def test_wrong_chain_head_refuses(self):
        self.contract['head']='0'*64
        with self.assertRaises(ValueError):self.branch()

    def test_wrong_complete_census_refuses(self):
        self.contract['completedRequests']=3
        with self.assertRaises(ValueError):self.branch()

    def test_wrong_panel_receipt_refuses(self):
        self.contract['panels'][0]['returnReceipt']['sha256']='0'*64
        with self.assertRaises(ValueError):self.branch()

    def test_missing_paired_order_refuses(self):
        self.contract['panels']=self.contract['panels'][:1]
        with self.assertRaises(ValueError):self.branch()

    def test_wrong_pending_descriptor_refuses(self):
        self.contract['pending']['descriptor']['arguments']['a']='different'
        with self.assertRaises(ValueError):self.branch()

    def test_prefix_chunk_mutation_detected(self):
        store,_=self.branch()
        with store.db:store.db.execute("UPDATE chunks SET payload=? WHERE sequence=0 AND ordinal=0",(b'changed',))
        with self.assertRaises(ValueError):store.prefix_preserved()

    def test_completed_return_link_mutation_detected(self):
        store,_=self.branch();h=P.hashlib.sha256(canonical(self.panels[0]['descriptor'])).hexdigest()
        with store.db:store.db.execute('UPDATE request_index SET output_sequence=0 WHERE request_hash=?',(h,))
        with self.assertRaises(ValueError):store.prefix_preserved()

    def test_original_mutation_detected(self):
        store,_=self.branch()
        with self.original.open('ab') as f:f.write(b'x')
        with self.assertRaises(ValueError):store.prefix_preserved()

    def test_namespace_mutation_detected(self):
        store,_=self.branch()
        with store.db:store.db.execute('UPDATE namespaces SET descriptor=? WHERE digest=?',(b'{}',self.ns))
        with self.assertRaises(ValueError):store.prefix_preserved()

    def test_new_request_storage_reserve_refusal(self):
        store,_=self.branch()
        with self.assertRaises(ValueError):P.ExplicitParentIndex(store,self.contract,free_bytes=lambda:0)


class EncodingTests(unittest.TestCase):
    def test_original_failure_and_exact_repair(self):
        left={'value':'opaque'};right={'value':'different'};old={24:left,48:right}
        with self.assertRaisesRegex(ValueError,'string JSON keys'):canonical({'pairedPanels':old})
        fixed=P.panel_evidence(old)
        self.assertIs(fixed['24'],left);self.assertIs(fixed['48'],right)
        self.assertEqual(json.loads(canonical(fixed)),{'24':left,'48':right})
        self.assertEqual(set(old),{24,48})

    def test_collision_and_unknown_keys_refuse(self):
        for value in ({24:1,'24':2,48:3},{'24':1,'48':2},{24:1},{24:1,48:2,96:3},{True:1,48:2},[24,48]):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):P.panel_evidence(value)

    def test_original_numeric_metadata_site_is_unique(self):
        source=(M/'S11c_d_defect_packet_contracted_numeric_lib.py').read_text();tree=ast.parse(source)
        sites=[]
        for n in ast.walk(tree):
            if isinstance(n,ast.Dict):
                for key,value in zip(n.keys,n.values):
                    if isinstance(key,ast.Constant) and key.value=='pairedPanels':sites.append(ast.dump(value))
        self.assertEqual(sites,["Name(id='panels', ctx=Load())"])
        self.assertIn('panels[24]',source);self.assertIn('panels[48]',source)

    def test_no_scientific_import_or_dynamic_execution_in_recovery(self):
        tree=ast.parse(Path(P.__file__).read_text())
        forbidden={'sympy','mpmath','numpy','scipy','pickle','dill'}
        for n in ast.walk(tree):
            if isinstance(n,ast.Import):self.assertFalse({a.name.split('.')[0] for a in n.names}&forbidden)
            if isinstance(n,ast.ImportFrom):self.assertNotIn((n.module or '').split('.')[0],forbidden)
            if isinstance(n,ast.Call) and isinstance(n.func,ast.Name):self.assertNotIn(n.func.id,{'eval','exec'})

    def test_only_evidence_boundary_changes_in_numeric_copy(self):
        original=ast.parse((M/'S11c_d_defect_packet_contracted_numeric_lib.py').read_text())
        candidate=ast.parse((M/'S11c_d_defect_packet_contracted_prefix_recovery_numeric.py').read_text())
        extra=[n for n in candidate.body if isinstance(n,ast.ImportFrom) and
               n.module=='S11c_d_defect_packet_contracted_prefix_recovery']
        self.assertEqual(len(extra),1);self.assertEqual([(a.name,a.asname) for a in extra[0].names],[('panel_evidence',None)])
        candidate.body.remove(extra[0]);sites=0
        for n in ast.walk(candidate):
            if isinstance(n,ast.Dict):
                for i,k in enumerate(n.keys):
                    if isinstance(k,ast.Constant) and k.value=='pairedPanels':
                        self.assertEqual(ast.dump(n.values[i]),"Call(func=Name(id='panel_evidence', ctx=Load()), args=[Name(id='panels', ctx=Load())], keywords=[])")
                        n.values[i]=ast.Name(id='panels',ctx=ast.Load());sites+=1
        self.assertEqual(sites,1);self.assertEqual(ast.dump(original),ast.dump(candidate))

    def test_actual_codec_transport_with_manufactured_values(self):
        tree=ast.parse((M/'S11c_d_defect_packet_contracted_numeric_lib.py').read_text())
        functions=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ('encode','encode_mpf','require')]
        class DummyEstimate:pass
        ns={'Estimate':DummyEstimate}
        exec(compile(ast.Module(body=functions,type_ignores=[]),'<source-extracted-stdlib-codec>','exec'),ns)
        raw={24:{'value':['synthetic',1]},48:{'value':['synthetic',2]}}
        with self.assertRaises(ValueError):canonical(ns['encode'](raw))
        self.assertEqual(json.loads(canonical(ns['encode'](P.panel_evidence(raw)))),
                         {'24':raw[24],'48':raw[48]})


if __name__=='__main__':unittest.main(verbosity=2)
