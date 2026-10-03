#!/usr/bin/env python3
"""Stdlib metadata/storage tests only. No scientific imports or evaluation."""
import ast
import hashlib
import json
from pathlib import Path
import runpy
import sqlite3
import tempfile
import unittest

M=Path(__file__).resolve().parent
PREFIX='S11c_d_defect_packet_fourier'
WORKER=M/(PREFIX+'.py');LIB=M/(PREFIX+'_lib.py');MANIFEST=M/(PREFIX+'_inputs.json')
W=runpy.run_path(str(WORKER),run_name='inert_worker_tests')
L=runpy.run_path(str(LIB),run_name='inert_library_tests')
m=json.loads(MANIFEST.read_text())
read=lambda p:json.loads(Path(p).read_text())

def function(tree,name):return next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==name)

class StorageTests(unittest.TestCase):
    def test_full_transaction_and_opaque_receipt(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'evidence.sqlite';s=L['DurableStore'](p)
            r=s.put('synthetic/input',{'label':'tooling only','operand':[1,'x',None]});s.close()
            with sqlite3.connect('file:'+str(p)+'?mode=ro',uri=True) as db:
                name,h,raw=db.execute('SELECT name,sha256,payload FROM records').fetchone()
            self.assertEqual(name,r['record']);self.assertEqual(hashlib.sha256(raw).hexdigest(),h);self.assertEqual(h,r['sha256']);self.assertEqual(len(raw),r['bytes'])
    def test_immutable_duplicate(self):
        with tempfile.TemporaryDirectory() as d:
            s=L['DurableStore'](Path(d)/'db');s.put('x',{'a':1})
            with self.assertRaises(sqlite3.IntegrityError):s.put('x',{'a':2})
            self.assertEqual(json.loads(s.db.execute('SELECT payload FROM records').fetchone()[0]),{'a':1});s.close()
    def test_existing_store_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'db';p.touch()
            with self.assertRaises(ValueError):L['DurableStore'](p)
    def test_nan_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            s=L['DurableStore'](Path(d)/'db')
            with self.assertRaises(ValueError):s.put('x',{'value':float('nan')})
            self.assertEqual(s.db.execute('SELECT count(*) FROM records').fetchone()[0],0);s.close()
    def test_complete_request_key(self):
        baseline={'fieldId':'native','coefficientsAscending':[['1','0']],'center':'-5/2','carrier':'0','role':'X','argumentDerivative':0,'spatialOrders':[1,0,0]}
        h=L['polynomial_request_key'](baseline)
        for key,value in [('fieldId','different'),('coefficientsAscending',[['2','0']]),('center','5/2'),('carrier','sqrt(595)/10'),('role','Y'),('argumentDerivative',1),('spatialOrders',[0,0,0])]:
            with self.subTest(key=key):self.assertNotEqual(h,L['polynomial_request_key']({**baseline,key:value}))
    def test_key_order_stable(self):self.assertEqual(L['digest']({'a':1,'b':2}),L['digest']({'b':2,'a':1}))
    def test_exact_true_only(self):
        L['require'](True,'ok')
        for x in (1,False,None,'true'):
            with self.assertRaises(ValueError):L['require'](x,'refuse')

class SourceTests(unittest.TestCase):
    def test_top_level_imports_stdlib(self):
        allowed={'argparse','ast','hashlib','itertools','json','math','os','pathlib','resource','shutil','sys','time','traceback','sqlite3','fractions'}
        for path in (WORKER,LIB):
            tree=ast.parse(path.read_text())
            for n in tree.body:
                if isinstance(n,ast.Import):self.assertTrue(all(a.name in allowed for a in n.names))
                if isinstance(n,ast.ImportFrom):self.assertIn(n.module,allowed)
    def test_no_timer_or_old_producer_calls(self):
        for path in (WORKER,LIB):
            tree=ast.parse(path.read_text());calls={ast.unparse(n.func) for n in ast.walk(tree) if isinstance(n,ast.Call)}
            self.assertFalse(calls&{'signal.alarm','signal.setitimer','build_face','run_science_old','sp.Integral','np.linalg.solve','time.sleep'})
    def test_science_import_follows_containment(self):
        main=function(ast.parse(WORKER.read_text()),'main');s=ast.unparse(main)
        self.assertLess(s.index("ns['containment']()"),s.index('import sympy'))
    def test_unchanged_journal_and_containment_helper_selection(self):
        self.assertEqual(W['HELPERS'],('require','sha','save','replace_json','containment','Journal','decode','one_symbol'))
    def test_failed_numerical_journal_finalized(self):
        run=function(ast.parse(WORKER.read_text()),'run_science')
        final=next(n.finalbody for n in ast.walk(run) if isinstance(n,ast.Try) and any('store.close' in ast.unparse(a) for a in n.finalbody))
        self.assertIn('numerical-journal-receipt',ast.unparse(ast.Module(body=final,type_ignores=[])))
    def test_gate_requires_actual_record_and_stage(self):
        s=ast.unparse(function(ast.parse(WORKER.read_text()),'verify_gate'))
        for text in ['READY_FOR_ONE_PACKET_FOURIER_BANK','CLEAR FOR THIS PACKET-ACTION FOURIER BUILD','independentBuildClearance','scientificRunsAuthorized','durationLimits']:
            self.assertIn(text,s)
    def test_negative_test_carrier_and_Y_derivative(self):
        s=ast.unparse(function(ast.parse(LIB.read_text()),'product'));self.assertIn("role == 'X' else -1",s);self.assertIn('argumentDerivative',s)
        worker=ast.unparse(function(ast.parse(WORKER.read_text()),'run_science'));self.assertIn('-sqrt(595)/10',worker)
    def test_independent_precision_and_rules(self):
        s=LIB.read_text()
        self.assertIn('self.A=mpmath.mp.clone();self.A.dps=30',s);self.assertIn('self.B=mpmath.mp.clone();self.B.dps=50',s)
        self.assertIn('for order in (24,48)',s);self.assertIn('kronrod15(self.B,store)',s)
    def test_no_relative_accuracy_promotion(self):
        self.assertIn("'relativeAccuracyClaim':False",LIB.read_text());self.assertIn("'fullActionAccuracyClaim':False",LIB.read_text())
    def test_source_hashes(self):
        for p,h in m['sourcePins'].items():
            with self.subTest(path=p):self.assertEqual(W['sha'](p),h)
    def test_original_complete_selection(self):
        saved=m['savedInputs'];a=read(saved['selected/pressure-addresses.json']['path'])['selected'];orig=read(saved['inventory/THETA_BALANCE-ordered-addresses.json']['path'])
        self.assertEqual(a,[x for x in orig if x['jet']['channel']=='e_W'])
        fam=W['family_plan'](a);self.assertEqual(len(fam),19)
        self.assertEqual(sum(f['spec']['role']=='X' for f in fam),13)
        self.assertEqual(sum(f['spec']['role']=='Y' and f['spec']['argumentDerivative']==1 for f in fam),3)
    def test_mutated_source_order_new_family(self):
        a=read(m['savedInputs']['selected/pressure-addresses.json']['path'])['selected'];live=next(x for x in a if x['status']=='FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED')
        other=json.loads(json.dumps(live));other['jet']['timeOrder']+=1;other['addressId']=-1
        self.assertEqual(len(W['family_plan']([live,other])),4)
    def test_authority_scope(self):
        a=read(m['executionAuthority']);self.assertEqual(a['scope'],m['scope']);self.assertEqual(a['scienceExecutionsAuthorized'],1);self.assertTrue(a['noDeadline']);self.assertFalse(a['automaticScientificRetry'])
    def test_native_input_and_helper_unchanged(self):
        old=read(M/'S11c_d_defect_packet_preflight_inputs.json')
        self.assertEqual(m['helperSource'],old['helperSource'])
        for n,r in old['savedInputs'].items():self.assertEqual(m['savedInputs'][n],r)
    def test_all_saved_files_full_receipts(self):
        for n,r in m['savedInputs'].items():
            with self.subTest(alias=n):self.assertEqual(Path(r['path']).stat().st_size,r['bytes']);self.assertEqual(W['sha'](r['path']),r['sha256'])
    def test_no_runtime_ready_gate(self):self.assertFalse((M/(PREFIX+'_gate.json')).exists())
    def test_launcher_external_runtime_snapshot(self):
        source=(M/(PREFIX+'_launch.py')).read_text();self.assertIn('external-runtime',source);self.assertIn("'--pool','s11c-near-unity'",source);self.assertIn("'--memory-gib','4'",source)

if __name__=='__main__':unittest.main()
