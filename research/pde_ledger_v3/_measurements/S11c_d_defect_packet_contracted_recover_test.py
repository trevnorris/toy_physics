"""Manufactured JSON/SQLite and AST tests only. No scientific restoration."""
from copy import deepcopy
import ast
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest import mock

import S11c_d_defect_packet_contracted_recover_protocol as P
from S11c_d_defect_packet_evidence_store import EvidenceStore
from S11c_d_defect_packet_request_index import RequestIndex, canonical

M = Path(__file__).resolve().parent


class Fixture(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name); self.original = self.root/'original.sqlite'
        self.destination = self.root/'branch.sqlite'
        store = EvidenceStore(self.original); index = RequestIndex(store, reserve_bytes=0)
        mathns=store.namespace({'route':'mathematical-inputs','settings':{'testOnly':True}})
        self.ns = store.namespace({'route':'A24', 'settings':{'purpose':'baseline','precision':30}})
        store.namespace({'route':'A24','settings':{'purpose':'formula-check:baseline','precision':30}})
        rules={'A24':{'path':'synthetic24','sha256':'a'*64,'bytes':24},
               'A48':{'path':'synthetic48','sha256':'b'*64,'bytes':48}}
        family={'meaning':'manufactured, no scientific values','ruleReceipts':rules}
        mathreceipt=index._put(mathns,'complete-numerical-context',{'contexts':family})
        context={'mathematicalContextReceipt':mathreceipt,'exactFamilyContexts':family}
        pc={'K':27,'T':122,'carrier':'synthetic-carrier','n':0,'purpose':'baseline',
            'cell':{'slabId':0},'epsilon':'synthetic-epsilon','geometryReceipt':'synthetic-geometry'}
        settings={'kind':'outer','ownRoute':'A24','selectedRule':'A24','squareMapped':True,
                  'ruleReceipt':rules['A24'],'context':pc}
        self.pending_desc={'route':'A24','purpose':'baseline','precision':30,'kind':'outer-panel',
                           'context':context,'arguments':{'a':'-149','b':'synthetic-b'},'settings':settings}
        parent=index.begin(self.ns,self.pending_desc)
        self.node={'panel':settings,'originalA':'-149','originalB':'synthetic-b','side':0,'ruleIndex':0,
                   'ruleNode':'synthetic-rule-node','representedPoint':'synthetic-outer-point','jacobian':'synthetic-jacobian'}
        nr=index._put(self.ns,'1/outer-node-input',self.node)
        self.panels=[];self.panel_values={}
        self.components={'Y0','X1','X1T','X2T','X02T','C0T','C1T','C2T','Y0CT','Y1CT','Y0C','Y1C'}
        for order in (24,48):
            innercontext={'n':0,'m':self.node['representedPoint'],'carrier':pc['carrier'],'K':27,'T':122,'M':149,
                          'outerParent':self.node,'primitiveKeys':['J','Dr','Dh','Dq'],'inputClips':'J,Dh,Dq,Dr_wrong_root',
                          'outputClips':'Dr','semanticVariables':{'X':'k','Y':'l'},'sharedQuadratureCoordinateDoesNotIdentifyKAndL':True,
                          'cell':{'cellId':0,'clipped':True}}
            inner={'kind':'inner','ownRoute':'A24','selectedRule':'A'+str(order),'squareMapped':True,
                   'ruleReceipt':rules['A'+str(order)],'context':innercontext}
            d={'route':'A24','purpose':'baseline','precision':30,'kind':'inner-panel','context':context,
               'arguments':{'a':'-27','b':'synthetic-inner-b'},'settings':inner}
            ticket=index.begin(self.ns,d)
            value={'components':{k:{'value':{'mpf':[0,'1',0,1]},'error':{'mpf':[0,'0',0,0]}} for k in self.components},
                   'a':d['arguments']['a'],'b':d['arguments']['b'],'settings':inner}
            self.panel_values['A'+str(order)]=value
            out=index.complete(ticket,value)
            self.panels.append({'namespace':self.ns,'descriptor':d,'inputReceipt':ticket['inputReceipt'],'returnReceipt':out})
        with store.db:
            store.db.execute('CREATE TABLE contracted_leaves (job TEXT, ordinal INTEGER, active INTEGER, receipt BLOB NOT NULL, PRIMARY KEY(job,ordinal))')
        index._put(self.ns,'580/manufactured-last-serial',{'synthetic':True})
        records,head=store.sequence,store.previous
        namespaces=[{'digest':d,'descriptor':json.loads(v)} for d,v in store.db.execute('SELECT digest,descriptor FROM namespaces ORDER BY digest')]
        store.close()
        self.contract={'databaseSha256':P.sha(self.original),'databaseBytes':self.original.stat().st_size,
                       'records':records,'head':head,'completedRequests':2,'pending':parent,'panels':self.panels,
                       'reserveBytes':0,'maxRecordBytes':8*1024**2,'cacheBytes':0,
                       'firstOuterNode':{'receipt':nr,'value':self.node},'lastCommittedNumericSerial':580,
                       'originalNamespaces':namespaces,'mathematicalContextReceipt':mathreceipt}
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
        with self.assertRaises(ValueError):index.claim_parent(self.ns,d)
        self.assertIsNone(index.lookup(self.ns,d))
        ticket=index.begin(self.ns,d);index.complete(ticket,{'new':True})
        self.assertEqual(index.lookup(self.ns,d)['value'],{'new':True})

    def test_other_pending_request_has_no_claim_permission(self):
        _,index=self.branch();d=deepcopy(self.pending_desc);d['arguments']['a']='new'
        index.begin(self.ns,d)
        with self.assertRaises(ValueError):index.claim_parent(self.ns,d)
        with self.assertRaises(ValueError):index.lookup(self.ns,d)

    def test_full_argument_change_cannot_claim_parent(self):
        _,index=self.branch();d=deepcopy(self.pending_desc);d['context']['meaning']='different'
        with self.assertRaises(ValueError):index.claim_parent(self.ns,d)
        self.assertFalse(index.claimed)

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
        with self.assertRaises(Exception):store.prefix_preserved()

    def test_completed_return_link_mutation_detected(self):
        store,_=self.branch();h=P.hashlib.sha256(canonical(self.panels[0]['descriptor'])).hexdigest()
        with store.db:store.db.execute('UPDATE request_index SET output_sequence=0 WHERE request_hash=?',(h,))
        with self.assertRaises(Exception):store.prefix_preserved()

    def test_original_mutation_detected(self):
        store,_=self.branch()
        with self.original.open('ab') as f:f.write(b'x')
        with self.assertRaises(Exception):store.prefix_preserved()

    def test_namespace_mutation_detected(self):
        store,_=self.branch()
        with store.db:store.db.execute('UPDATE namespaces SET descriptor=? WHERE digest=?',(b'{}',self.ns))
        with self.assertRaises(Exception):store.prefix_preserved()

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


class ProtocolTests(unittest.TestCase):
    # Share only fixture setup, not duplicate test-method executions.
    setUp = Fixture.setUp
    branch = Fixture.branch
    def protocol(self):
        store,index=self.branch()
        ns=store.namespace({'route':'mathematical-inputs','settings':{'execution':'synthetic-recovery'}})
        receipt=index._put(ns,'recovery-execution-provenance',{'originalDatabaseSha256':self.contract['databaseSha256'],'mathematicalNamespaceUnchanged':True})
        return store,index,P.ExplicitRecoveryProtocol(index,self.contract,receipt)

    def dead_callback(self):
        raise AssertionError('completed child callback MUST NOT execute')

    def enter_children(self,p):
        got=p.before_put('A24','baseline','outer-node-input',self.node)
        self.assertEqual(got,self.contract['firstOuterNode']['receipt'])
        values={}
        for witness in self.panels:
            rule=witness['descriptor']['settings']['selectedRule']
            values[rule]=p.request(self.ns,witness['descriptor'],self.dead_callback)
        return values

    def publish_pair(self,p,index,values):
        d=self.panels[0]['descriptor']
        encoded={'a':d['arguments']['a'],'b':d['arguments']['b'],'meta':d['settings']['context']['cell'],
                 'pairedPanels':{'24':values['A24'],'48':values['A48']},'ownOrder':24,
                 'positiveComponentContributions':{k:dict(values['A24']['components'][k]) for k in self.components}}
        for k,v in encoded['positiveComponentContributions'].items():
            p.observe_pair(k,values['A24']['components'][k],values['A48']['components'][k],24,v)
        self.assertIsNone(p.before_put('A24','baseline','inner-positive-panel-pair',encoded))
        rr=index._put(self.ns,'581/inner-positive-panel-pair',encoded)
        p.after_put('A24','baseline','inner-positive-panel-pair',encoded,rr)
        return rr

    def test_complete_protocol_has_zero_child_callback_calls(self):
        store,index,p=self.protocol();before=store.sequence
        def unfinished():
            values=self.enter_children(p);self.publish_pair(p,index,values)
            return {'newUnfinishedParentReturn':'synthetic'}
        result=p.request(self.ns,self.pending_desc,unfinished)
        self.assertEqual(result,{'newUnfinishedParentReturn':'synthetic'})
        self.assertEqual(store.sequence,before+3) # Pair, durable recovery ledger, parent return; no duplicate node.
        final=p.final_record();self.assertTrue(final['parentCompleted'])
        self.assertEqual(final['prefix']['parentState'],'COMPLETE')
        self.assertEqual([e['event'] for e in final['events']],['explicit-parent-claimed',
            'first-node-input-verified-and-skipped','completed-child-reused','completed-child-reused',
            'unpublished-parent-scalars-reconstructed','original-parent-completed'])
        self.assertEqual(self.original.read_bytes(),self.old_bytes)

    def test_parent_mismatch_refuses_before_begin_and_callback(self):
        store,index,p=self.protocol();d=deepcopy(self.pending_desc);d['arguments']['a']='nearby'
        with mock.patch.object(index,'begin',side_effect=AssertionError('must not begin')):
            with self.assertRaises(ValueError):p.request(self.ns,d,self.dead_callback)
        self.assertFalse(index.claimed);self.assertEqual(store.sequence,self.contract['records']+1)

    def test_child_mismatch_refuses_before_begin_and_callback(self):
        store,index,p=self.protocol()
        def unfinished():
            p.before_put('A24','baseline','outer-node-input',self.node)
            d=deepcopy(self.panels[0]['descriptor']);d['arguments']['b']='nearby'
            return p.request(self.ns,d,self.dead_callback)
        with mock.patch.object(index,'begin',side_effect=AssertionError('must not begin')):
            with self.assertRaises(ValueError):p.request(self.ns,self.pending_desc,unfinished)
        self.assertEqual(p.state,'FAILED');self.assertEqual(store.sequence,self.contract['records']+1)

    def test_child_lookup_miss_refuses_before_begin_and_callback(self):
        store,index,p=self.protocol()
        def unfinished():
            p.before_put('A24','baseline','outer-node-input',self.node)
            with mock.patch.object(index,'lookup',return_value=None):
                return p.request(self.ns,self.panels[0]['descriptor'],self.dead_callback)
        with mock.patch.object(index,'begin',side_effect=AssertionError('must not begin')):
            with self.assertRaisesRegex(ValueError,'lookup miss'):p.request(self.ns,self.pending_desc,unfinished)
        self.assertEqual(store.sequence,self.contract['records']+1)

    def test_wrong_child_order_refuses(self):
        _,_,p=self.protocol()
        def unfinished():
            p.before_put('A24','baseline','outer-node-input',self.node)
            return p.request(self.ns,self.panels[1]['descriptor'],self.dead_callback)
        with self.assertRaises(ValueError):p.request(self.ns,self.pending_desc,unfinished)

    def test_no_request_before_saved_node_join(self):
        _,_,p=self.protocol()
        with self.assertRaises(ValueError):
            p.request(self.ns,self.pending_desc,lambda:p.request(self.ns,self.panels[0]['descriptor'],self.dead_callback))

    def test_changed_node_input_refuses_no_duplicate_put(self):
        store,_,p=self.protocol();node=deepcopy(self.node);node['representedPoint']='changed'
        with self.assertRaises(ValueError):
            p.request(self.ns,self.pending_desc,lambda:p.before_put('A24','baseline','outer-node-input',node))
        self.assertEqual(store.sequence,self.contract['records']+1)

    def test_no_new_work_before_pair_commit(self):
        _,_,p=self.protocol()
        def unfinished():
            self.enter_children(p)
            d=deepcopy(self.pending_desc);d['arguments']['a']='new'
            return p.request(self.ns,d,self.dead_callback)
        with self.assertRaises(ValueError):p.request(self.ns,self.pending_desc,unfinished)
        self.assertEqual(p.state,'FAILED')

    def test_incomplete_parent_callback_cannot_publish_return(self):
        store,_,p=self.protocol()
        with self.assertRaises(ValueError):p.request(self.ns,self.pending_desc,lambda:{'tooSoon':True})
        self.assertEqual(store.sequence,self.contract['records']+1);self.assertFalse(p.parent_done)

    def test_existing_namespace_mismatch_never_calls_store_insertion(self):
        store,_,p=self.protocol();d=next(x['descriptor'] for x in self.contract['originalNamespaces']
                                      if x['digest']==self.ns);d=deepcopy(d);d['settings']['precision']=31
        with mock.patch.object(store,'namespace',side_effect=AssertionError('must not insert')):
            with self.assertRaises(ValueError):p.namespace(d)

    def test_existing_namespace_uses_original_digest_without_insertion(self):
        store,_,p=self.protocol();d=next(x['descriptor'] for x in self.contract['originalNamespaces'] if x['digest']==self.ns)
        with mock.patch.object(store,'namespace',side_effect=AssertionError('must not insert')):
            self.assertEqual(p.namespace(d),self.ns)

    def test_new_namespace_refuses_before_prefix_done(self):
        store,_,p=self.protocol()
        with mock.patch.object(store,'namespace',side_effect=AssertionError('must not insert')):
            with self.assertRaises(ValueError):p.namespace({'route':'A48','settings':{'purpose':'baseline','precision':30}})

    def test_frame_mutants_refuse(self):
        for location in ('node','point','cell','carrier','window','rule','own-route','input-endpoint','side'):
            parent=deepcopy(self.contract['pending']);node=deepcopy(self.node);panels=deepcopy(self.panels)
            c=panels[0]['descriptor']['settings']['context']
            if location=='node':c['outerParent']['ruleIndex']=1
            if location=='point':c['m']='wrong'
            if location=='cell':c['cell']['cellId']=1
            if location=='carrier':c['carrier']='wrong'
            if location=='window':c['T']=124
            if location=='rule':panels[0]['descriptor']['settings']['ruleReceipt']['path']='wrong'
            if location=='own-route':panels[0]['descriptor']['route']='A48'
            if location=='input-endpoint':panels[0]['descriptor']['arguments']['a']='wrong'
            if location=='side':node['side']=1
            with self.subTest(location=location):
                with self.assertRaises(ValueError):P.join_frame(parent,node,panels)

    def test_pair_return_must_be_actual_committed_record(self):
        _,index,p=self.protocol()
        def unfinished():
            values=self.enter_children(p);d=self.panels[0]['descriptor']
            encoded={'a':d['arguments']['a'],'b':d['arguments']['b'],'meta':d['settings']['context']['cell'],
                     'pairedPanels':{'24':values['A24'],'48':values['A48']},'ownOrder':24,
                     'positiveComponentContributions':{k:'synthetic' for k in self.components}}
            p.before_put('A24','baseline','inner-positive-panel-pair',encoded)
            p.after_put('A24','baseline','inner-positive-panel-pair',encoded,self.panels[0]['returnReceipt'])
        with self.assertRaises(ValueError):p.request(self.ns,self.pending_desc,unfinished)
        self.assertEqual(p.state,'FAILED')

    def test_final_record_preserves_failed_attempt(self):
        _,_,p=self.protocol();d=deepcopy(self.pending_desc);d['arguments']['a']='wrong'
        with self.assertRaises(ValueError):p.request(self.ns,d,self.dead_callback)
        final=p.final_record();self.assertEqual(final['attempts'][0]['descriptor'],d)
        self.assertEqual(final['prefix']['parentState'],'PENDING');self.assertFalse(final['parentCompleted'])

    def test_failed_attempt_is_frozen_before_caller_mutation(self):
        _,_,p=self.protocol();d=deepcopy(self.pending_desc);d['arguments']['a']='wrong'
        with self.assertRaises(ValueError):p.request(self.ns,d,self.dead_callback)
        d['arguments']['a']='mutated-after-refusal'
        self.assertEqual(p.final_record()['attempts'][0]['descriptor']['arguments']['a'],'wrong')

    def test_namespace_attempt_is_frozen_before_caller_mutation(self):
        _,_,p=self.protocol();d={'route':'A48','settings':{'purpose':'baseline','precision':30}}
        with self.assertRaises(ValueError):p.namespace(d)
        d['settings']['purpose']='changed'
        self.assertEqual(p.final_record()['attempts'][0]['descriptor']['settings']['purpose'],'baseline')

    def test_copy_disk_refusal_before_destination_creation(self):
        with self.assertRaisesRegex(ValueError,'clone disk reserve'):
            P.ExplicitCloneStore(self.original,self.destination,self.contract,free_bytes=lambda:0)
        self.assertFalse(self.destination.exists())

    def test_partial_copy_is_preserved_exclusively(self):
        def broken(inp,out,length):
            out.write(b'partial');raise OSError('manufactured write failure')
        with mock.patch.object(P.shutil,'copyfileobj',side_effect=broken):
            with self.assertRaises(OSError):P.ExplicitCloneStore(self.original,self.destination,self.contract)
        self.assertEqual(self.destination.read_bytes(),b'partial')
        with self.assertRaises(ValueError):self.branch()
        self.assertEqual(self.original.read_bytes(),self.old_bytes)

    def test_file_and_directory_fsync(self):
        with mock.patch.object(P.os,'fsync',wraps=P.os.fsync) as fs:
            self.branch();self.assertEqual(fs.call_count,2)

    def test_serial_is_derived_not_accepted_as_arbitrary_pin(self):
        self.contract['lastCommittedNumericSerial']=579
        with self.assertRaises(ValueError):self.branch()



class BuildProtocolTests(unittest.TestCase):
    setUp=Fixture.setUp
    branch=Fixture.branch
    protocol=ProtocolTests.protocol
    enter_children=ProtocolTests.enter_children
    publish_pair=ProtocolTests.publish_pair
    dead_callback=ProtocolTests.dead_callback

    def pair_data(self,p):
        values=self.enter_children(p);d=self.panels[0]['descriptor']
        out={'a':d['arguments']['a'],'b':d['arguments']['b'],'meta':d['settings']['context']['cell'],
             'pairedPanels':{'24':values['A24'],'48':values['A48']},'ownOrder':24,
             'positiveComponentContributions':{k:dict(values['A24']['components'][k]) for k in self.components}}
        for k,v in out['positiveComponentContributions'].items():
            p.observe_pair(k,values['A24']['components'][k],values['A48']['components'][k],24,v)
        return out

    def test_terminal_refusal_and_callback_cannot_retry(self):
        _,_,p=self.protocol()
        with self.assertRaises(RuntimeError):p.request(self.ns,self.pending_desc,lambda:(_ for _ in ()).throw(RuntimeError('synthetic stop')))
        self.assertEqual(p.state,'FAILED')
        with self.assertRaisesRegex(ValueError,'terminal'):p.request(self.ns,self.pending_desc,self.dead_callback)

    def test_caught_child_refusal_latches(self):
        _,_,p=self.protocol()
        def callback():
            p.before_put('A24','baseline','outer-node-input',self.node)
            d=deepcopy(self.panels[0]['descriptor']);d['arguments']['a']='wrong'
            with self.assertRaises(ValueError):p.request(self.ns,d,self.dead_callback)
            return p.request(self.ns,self.panels[0]['descriptor'],self.dead_callback)
        with self.assertRaisesRegex(ValueError,'terminal'):p.request(self.ns,self.pending_desc,callback)

    def assert_pair_mutant(self,mutate):
        _,_,p=self.protocol()
        def callback():
            out=self.pair_data(p);mutate(out)
            p.before_put('A24','baseline','inner-positive-panel-pair',out)
        with self.assertRaises((ValueError,KeyError)):p.request(self.ns,self.pending_desc,callback)
        self.assertEqual(p.state,'FAILED')

    def test_pair_mutated_panels(self):self.assert_pair_mutant(lambda x:x['pairedPanels']['24'].update(a='wrong'))
    def test_pair_mutated_order(self):self.assert_pair_mutant(lambda x:x.update(ownOrder=48))
    def test_pair_missing_component(self):self.assert_pair_mutant(lambda x:x['positiveComponentContributions'].pop('Y0'))
    def test_pair_mutated_own_value(self):self.assert_pair_mutant(lambda x:x['positiveComponentContributions']['Y0'].update(value={'mpf':[0,'3',0,2]}))
    def test_pair_mutated_error(self):self.assert_pair_mutant(lambda x:x['positiveComponentContributions']['Y0'].update(error={'mpf':[0,'1',0,1]}))
    def test_pair_bad_shape(self):self.assert_pair_mutant(lambda x:x['positiveComponentContributions']['Y0'].update(extra=True))
    def test_pair_bad_error_atom(self):self.assert_pair_mutant(lambda x:x['positiveComponentContributions']['Y0'].update(error={'mpf':[0,'NaN',0,-1]}))

    def bad_pair_receipt(self,kind):
        store,index,p=self.protocol()
        def callback():
            out=self.pair_data(p);p.before_put('A24','baseline','inner-positive-panel-pair',out)
            if kind=='sequence':index._put(self.ns,'unwanted-intervening-record',{})
            name={'name':'581/wrong-name','serial':'582/inner-positive-panel-pair'}.get(kind,'581/inner-positive-panel-pair')
            receipt=index._put(self.ns,name,out)
            p.after_put('A24','baseline','inner-positive-panel-pair',out,receipt)
        with self.assertRaises(ValueError):p.request(self.ns,self.pending_desc,callback)
        self.assertEqual(p.state,'FAILED')
    def test_pair_wrong_sequence(self):self.bad_pair_receipt('sequence')
    def test_pair_wrong_serial(self):self.bad_pair_receipt('serial')
    def test_pair_wrong_name(self):self.bad_pair_receipt('name')

    def test_interruption_after_pair_commit_preserves_pending(self):
        store,index,p=self.protocol()
        def callback():
            out=self.pair_data(p);p.before_put('A24','baseline','inner-positive-panel-pair',out)
            index._put(self.ns,'581/inner-positive-panel-pair',out)
            raise KeyboardInterrupt('manufactured interruption before after_put')
        with self.assertRaises(KeyboardInterrupt):p.request(self.ns,self.pending_desc,callback)
        final=p.final_record()
        self.assertEqual(final['state'],'FAILED');self.assertIsNone(final['recoveryLedgerReceipt'])
        self.assertEqual(final['prefix']['parentState'],'PENDING')
        self.assertEqual(store.sequence,self.contract['records']+2)
        with self.assertRaises(ValueError):p.request(self.ns,self.pending_desc,self.dead_callback)

    def test_exact_sequence_durable_ledger_then_new_reuse(self):
        store,index,p=self.protocol();reused=[]
        def callback():
            values=self.enter_children(p);rr=self.publish_pair(p,index,values)
            self.assertEqual(rr['sequence'],self.contract['records']+1)
            self.assertEqual(p.recovery_ledger_receipt['sequence'],self.contract['records']+2)
            self.assertEqual(p.state,'NEW')
            value=p.request(self.ns,self.panels[0]['descriptor'],self.dead_callback,on_reuse=lambda r:reused.append(r['receipt']))
            self.assertEqual(value,values['A24'])
            return {'unfinishedParent':'now-complete'}
        p.request(self.ns,self.pending_desc,callback)
        self.assertEqual(reused,[self.panels[0]['returnReceipt']]);self.assertTrue(p.parent_done)
        self.assertEqual(p.final_record()['prefix']['originalLeafRows'],0)

    def test_schema_mutation_detected(self):
        store,_=self.branch();store.db.execute('CREATE TABLE extra(x)');store.db.commit()
        with self.assertRaisesRegex(ValueError,'schema'):store.prefix_preserved()

    def test_new_suffix_chunk_corruption_detected(self):
        store,index=self.branch();rr=index._put(self.ns,'new',{'synthetic':'data'})
        store.db.execute('UPDATE chunks SET payload=? WHERE sequence=?',(b'corrupt',rr['sequence']));store.db.commit()
        with self.assertRaises(Exception):store.prefix_preserved()

    def test_new_pending_request_reported(self):
        store,index=self.branch();d=deepcopy(self.pending_desc);d['arguments']['a']='new'
        ticket=index.begin(self.ns,d);out=store.prefix_preserved()
        self.assertEqual(out['allRequestStates']['PENDING'],2)
        self.assertIn(ticket['requestHash'],[x['requestHash'] for x in out['pendingRequests']])

if __name__=='__main__':unittest.main(verbosity=2)
