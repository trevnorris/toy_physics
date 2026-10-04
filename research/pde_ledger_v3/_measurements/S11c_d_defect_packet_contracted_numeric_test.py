"""Stdlib-only synthetic/metadata tests. No native plan or scientific input restored."""
import ast
import cmath
from copy import deepcopy
from fractions import Fraction as F
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import sqlite3

import S11c_d_defect_packet_contracted_geometry as G
import S11c_d_defect_packet_contracted_numeric_lib as N
import S11c_d_defect_packet_contracted_validation as V
from S11c_d_defect_packet_evidence_store import EvidenceStore
from S11c_d_defect_packet_request_index import RequestIndex

HERE=Path(__file__).parent


class Arithmetic(unittest.TestCase):
    def test_outer_leaf_ranks_binding_nested_guard(self):
        total={'a':N.Estimate(0,F(9,10))};nested={'a':F(4,5)}
        done,key,metric,details=N.outer_decision(total,nested,1)
        self.assertFalse(done);self.assertEqual(details['bindingGuard'],'nested')
        left=metric({'a':N.Estimate(0,F(6,10))},{'a':F(1,5)},key)
        right=metric({'a':N.Estimate(0,F(3,10))},{'a':F(3,5)},key)
        self.assertGreater(right,left)
    def test_outer_leaf_ranks_total_when_binding(self):
        done,key,metric,details=N.outer_decision({'a':N.Estimate(0,4)},{'a':F(1,2)},1)
        self.assertEqual(details['bindingGuard'],'total');self.assertEqual(metric({'a':N.Estimate(0,3)},{'a':F(1,9)},key),3)
    def test_outer_checks_all_components(self):
        done,key,_,_=N.outer_decision({'a':N.Estimate(0,0),'b':N.Estimate(0,F(1,3))},{'a':0,'b':F(1,3)},1)
        self.assertFalse(done);self.assertEqual(key,'b')
    def test_early_leaf_receipt(self):
        with sqlite3.connect(':memory:') as db:self.assertEqual(N.leaf_count(db),{'initialized':False,'active':0})
    def test_initialized_leaf_receipt(self):
        with sqlite3.connect(':memory:') as db:
            db.execute('CREATE TABLE contracted_leaves(active INTEGER)');db.executemany('INSERT INTO contracted_leaves VALUES (?)',[(1,),(0,),(1,)])
            self.assertEqual(N.leaf_count(db),{'initialized':True,'active':2})
    def test_wrong_root_only_has_its_dr_label(self):
        e=N.Evaluator.__new__(N.Evaluator);e.contexts={'entries':[{'addressId':8347,'n':0,'alpha':['1','0'],'primitives':['Dr','Dh','Dq']}]}
        class C:
            j=1j
            mpf=staticmethod(F)
        self.assertEqual(e.addressed(C(),0,'wrong-root-mutant',{'Dr_wrong_root':N.Estimate(5,F(1,4))}),{'8347/Dr':N.Estimate(5+0j,F(1,4))})
    def test_positive_product_errors(self):
        a,b=N.Estimate(F(-3),F(1,10)),N.Estimate(F(7),F(1,5))
        v=a*b;self.assertEqual(v.value,-21);self.assertEqual(v.error,F(66,50))
    def test_addition_does_not_cancel_errors(self):
        v=N.Estimate(2,F(1,10))+N.Estimate(-2,F(1,5));self.assertEqual((v.value,v.error),(0,F(3,10)))
    def test_negative_error_refuses(self):
        with self.assertRaises(ValueError):N.Estimate(1,-1)
    def test_fixed_divisor(self):
        v=N.Estimate(4,F(1,5))/-2;self.assertEqual((v.value,v.error),(-2,F(1,10)))
    def test_uncertain_divisor_refuses(self):
        with self.assertRaises(ValueError):N.Estimate(1,0)/N.Estimate(2,0)
    def test_zero_divisor_refuses(self):
        with self.assertRaises(ValueError):N.Estimate(1,0)/0
    def test_combine_uses_distinct_prefactors(self):
        z={n:N.Estimate(F(i+1),F(1,100)) for i,n in enumerate(['Y0','X1','X1T','X2T','X02T','C0T','C1T','C2T','Y0CT','Y1CT','Y0C','Y1C'])}
        a=N.combine(F(2),F(3),F(5),F(7),F(11),z);b=N.combine(F(2),F(3),F(5),F(14),F(11),z)
        self.assertEqual(b['J'].value,2*a['J'].value)
        for k in ('Dr','Dh','Dq','Dr_wrong_root'):self.assertEqual(a[k],b[k])
    def test_reflected_control_changes_assembly(self):
        z={n:N.Estimate(F(i+1),0) for i,n in enumerate(['Y0','X1','X1T','X2T','X02T','C0T','C1T','C2T','Y0CT','Y1CT','Y0C','Y1C'])}
        v=N.combine(F(2),F(3),F(5),F(7),F(11),z);self.assertNotEqual(v['Dr'].value,v['Dr_wrong_root'].value)

    def test_panel_discrepancies_cannot_cancel(self):
        first=N.paired_indicator({'value':7,'error':0},{'value':3,'error':0},48)
        second=N.paired_indicator({'value':-7,'error':0},{'value':-3,'error':0},48)
        self.assertEqual((first+second).value,0)
        self.assertEqual((first+second).error,8)
    def test_pair_keeps_own_route_value(self):
        left={'value':3,'error':F(1,10)};right={'value':4,'error':F(1,5)}
        self.assertEqual(N.paired_indicator(left,right,24).value,3)
        self.assertEqual(N.paired_indicator(left,right,48).value,4)
        self.assertEqual(N.paired_indicator(left,right,48).error,F(13,10))


class ToyContext:
    one=1;zero=0;j=1j;pi=cmath.pi
    exp=staticmethod(cmath.exp);sqrt=staticmethod(cmath.sqrt)


class ExactCodec(unittest.TestCase):
    class C:
        @staticmethod
        def make_mpf(t):return type('Number',(),{'_mpf_':t})()
        @staticmethod
        def isfinite(v):return True
    def test_exact_tuple_roundtrip(self):
        value={'mpf':[1,'123456789',-19,27]}
        self.assertEqual(N.encode(N.decode(self.C(),value)),value)
    def test_tuple_normalization_refuses(self):
        class Bad(self.C):
            @staticmethod
            def make_mpf(t):return ExactCodec.C.make_mpf((t[0],t[1]+1,t[2],t[3]))
        with self.assertRaises(ValueError):N.decode(Bad(),{'mpf':[0,'5',-2,3]})
    def test_bad_mantissa_refuses(self):
        with self.assertRaises(ValueError):N.decode(self.C(),{'mpf':[0,'not-an-integer',-2,3]})


class Selector(unittest.TestCase):
    def setUp(self):
        self.address={'addressId':17,'sourceTransform':{'coefficientId':'fixture-X'},'consumerTransform':{'coefficientId':'fixture-Y'},'jet':{'spatialOrders':[2,0,0],'timeOrder':0}}
        self.x={'role':'X','argumentDerivative':0,'fieldId':'fixture-X','timeOrder':0,'spatialOrders':[2,0,0],'center':'-5/2','width':'8','profileLength':'10','coefficientsAscending':[['7','2']]}
        self.y={**self.x,'role':'Y','fieldId':'fixture-Y','spatialOrders':[0,0,0],'center':'5/2'}
        self.units={'matchingXFamilyInterfaces':[self.x],'matchingYFamilyInterfaces':[{'spec':self.y,'addresses':[17]}]}
    def test_actual_selector(self):self.assertIs(V.original_selector(self.units,self.address,'source',[['7','2']]),self.x)
    def test_nonzero_selector_refuses(self):
        self.x['argumentDerivative']=1
        with self.assertRaises(ValueError):V.original_selector(self.units,self.address,'source',[['7','2']])
    def test_boolean_selector_refuses(self):
        self.x['argumentDerivative']=False
        with self.assertRaises(ValueError):V.original_selector(self.units,self.address,'source',[['7','2']])
    def test_missing_interface_refuses(self):
        self.units['matchingXFamilyInterfaces']=[]
        with self.assertRaises(ValueError):V.original_selector(self.units,self.address,'source',[['7','2']])
    def test_wrong_y_address_refuses(self):
        self.units['matchingYFamilyInterfaces'][0]['addresses']=[18]
        with self.assertRaises(ValueError):V.original_selector(self.units,self.address,'consumer',[['7','2']])
    def test_wrong_coefficient_refuses(self):
        with self.assertRaises(ValueError):V.original_selector(self.units,self.address,'consumer',[['8','2']])
    def test_wrong_order_refuses(self):
        self.x['spatialOrders']=[0,0,0]
        with self.assertRaises(ValueError):V.original_selector(self.units,self.address,'source',[['7','2']])


class ManufacturedGaussian(unittest.TestCase):
    def test_moment_vs_derivative_unrelated_data(self):
        c=ToyContext()
        for n in (0,1,2,3):
            for p in (-.7,.3,1.1):
                a=N.gaussian(c,p,.2,-.4,2,n,'X','A');b=N.gaussian(c,p,.2,-.4,2,n,'X','B')
                self.assertLess(abs(a[0]-b[0]),1e-12)
    def test_mutant_actual_b_route(self):
        c=ToyContext();args=(c,.6,.2,-.4,2,2,'X')
        a=N.gaussian(*args,'A',True);b=N.gaussian(*args,'B',True);base=N.gaussian(*args,'B',False)
        self.assertLess(abs(a[0]-b[0]),1e-14);self.assertGreater(abs(base[0]-b[0]),1e-4)
    def test_zero_carrier_mutant(self):
        c=ToyContext();self.assertEqual(N.gaussian(c,.6,0,-.4,2,2,'X','B',True)[0],0)
    def test_y_negative_argument_phase(self):
        c=ToyContext();p=.6;carrier=.2;center=.4;width=2
        v=N.gaussian(c,p,carrier,center,width,0,'Y','B')[0]
        expected=width*cmath.sqrt(2*cmath.pi)*cmath.exp(-width**2*(carrier-p)**2/2+1j*(p-carrier)*center)
        self.assertLess(abs(v-expected),1e-14)
    def test_no_extra_conjugation(self):
        c=ToyContext();v=N.gaussian(c,.6,.2,.4,2,0,'Y','A')[0]
        self.assertGreater(abs(v-v.conjugate()),.1)
    def test_stripped_amplitude_reconstruction(self):
        for kind in ('X','Y'):
            f,a,g=N.gaussian(ToyContext(),.7,.2,-.4,2,2,kind,'B');self.assertEqual(f,g*a)


class SyntheticGeometry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.k=G.Quad(0,1,F(2));cls.plan=G.mesh(cls.k,3,7,cls.k*0,5,3)
    def test_six_slabs(self):self.assertEqual(len(self.plan['affinity']),6)
    def test_full_wings(self):self.assertEqual((self.plan['slabs'][0]['left'],self.plan['slabs'][-1]['right']),(-10,10))
    def test_true_window(self):self.assertEqual(G.true_window(3,7,F(-9)),(-3,F(-2)))
    def test_right_window(self):self.assertEqual(G.true_window(3,7,F(9)),(F(2),3))
    def test_audit(self):self.assertTrue(G.audit(self.plan))
    def test_left_corruption(self):self.assertEqual(G.wing_control(self.plan,-1)['refusal'],'actual max/min clipping window')
    def test_right_corruption(self):self.assertEqual(G.wing_control(self.plan,1)['refusal'],'actual max/min clipping window')
    def test_missing_cell_refuses(self):
        p=deepcopy(self.plan);p['slabs'][5]['cells'].pop()
        with self.assertRaises(ValueError):G.audit(p)
    def test_reversed_interval_refuses(self):
        p=deepcopy(self.plan);s=p['slabs'][2];s['left'],s['right']=s['right'],s['left']
        with self.assertRaises(ValueError):G.audit(p)
    def test_missing_resolution_cut_refuses(self):
        p=deepcopy(self.plan);p['cuts'].pop(2)
        with self.assertRaises(ValueError):G.audit(p)
    def test_wrong_clip_membership_refuses(self):
        p=deepcopy(self.plan);p['slabs'][0]['cells'][0]['clipped']=not p['slabs'][0]['cells'][0]['clipped']
        with self.assertRaises(ValueError):G.audit(p)
    def test_no_float_inputs(self):
        with self.assertRaises(ValueError):G.mesh(self.k,3.0,7,self.k*0,5,3)
    def test_carrier_plan(self):self.assertTrue(G.audit(G.mesh(self.k,3,7,self.k,5,3)))
    def test_labels_preserved(self):
        for rec in self.plan['affinity']:
            labels=[l for line in rec['graphs'] for l in line.labels]
            self.assertEqual(len(labels),len(set(labels)));self.assertIn('clip-low',labels);self.assertIn('clip-high',labels)
    def test_serialization_complete(self):
        encoded=json.loads(json.dumps(G.packed(self.plan)));self.assertEqual(len(encoded['slabs']),len(self.plan['slabs']))


class RequestOrchestration(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.store=EvidenceStore(Path(self.tmp.name)/'test.sqlite')
        self.index=RequestIndex(self.store,reserve_bytes=0,free_bytes=lambda:10**12)
        self.e=N.Evaluator.__new__(N.Evaluator);self.e.store=self.store;self.e.index=self.index;self.e.contexts={'manufactured':'no native values'};self.e.math_receipt={'synthetic':'only'};self.e.namespaces={};self.e.serial=0
        self.e.rules={r:(type('Context',(),{'dps':p})(),[],[],{}) for r,p in [('A24',30),('A48',30),('B50',50)]}
    def tearDown(self):self.store.close();self.tmp.cleanup()
    def test_exact_completed_lookup_no_callback(self):
        calls=[]
        fn=lambda:calls.append(1) or {'v':3}
        for _ in range(2):self.assertEqual(self.e.request('A24','baseline','synthetic',{'x':1},{'setting':2},fn),{'v':3})
        self.assertEqual(calls,[1])
    def test_route_purpose_separation(self):
        calls=[]
        for route,purpose in [('A24','baseline'),('A48','baseline'),('A24','mutant')]:self.e.request(route,purpose,'synthetic',{'x':1},{'s':2},lambda:calls.append(1) or 3)
        self.assertEqual(len(calls),3)
    def test_pending_failure_no_retry(self):
        def fail():raise RuntimeError('manufactured')
        with self.assertRaises(RuntimeError):self.e.request('A24','baseline','synthetic',{'x':1},{'s':2},fail)
        with self.assertRaises(ValueError):self.e.request('A24','baseline','synthetic',{'x':1},{'s':2},lambda:3)
        self.assertEqual(self.store.db.execute('SELECT state FROM request_index').fetchone()[0],'PENDING')
    def test_full_changed_arguments(self):
        calls=[]
        for x in (1,2):self.e.request('A24','baseline','synthetic',{'x':x},{'s':2},lambda:calls.append(1) or 3)
        self.assertEqual(len(calls),2)


class SourceContracts(unittest.TestCase):
    def test_actual_tail_and_profile_ast_only(self):
        contracts=json.loads((HERE/'S11c_d_defect_packet_contracted_numeric_build_source_contracts.json').read_text())
        class C:
            @staticmethod
            def fragment(record):return ast.parse(record['text']).body[0]
        class J:
            @staticmethod
            def emit(name,value):pass
        V.source_ast_joins(contracts,ast.parse((HERE/'S11c_d_defect_packet_contracted_numeric_lib.py').read_text()),ast.parse((HERE/'S11c_d_defect_packet_contracted_prepare.py').read_text()),C(),J())
    def test_no_science_import_at_library_top(self):
        for name in ('contracted_numeric_lib','contracted_geometry','contracted_prepare'):
            tree=ast.parse((HERE/('S11c_d_defect_packet_'+name+'.py')).read_text())
            imports=[n for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom))]
            self.assertFalse(any('sympy' in ast.unparse(n) or 'mpmath' in ast.unparse(n) for n in imports))
    def test_no_timer_or_rule_constructor(self):
        text=(HERE/'S11c_d_defect_packet_contracted_numeric_lib.py').read_text()
        for forbidden in ('gauss_quadrature','quadgl(','quadts(','signal.alarm','setitimer','perf_counter'):
            self.assertNotIn(forbidden,text)
    def test_strict_control_selection(self):
        text=(HERE/'S11c_d_defect_packet_contracted_prepare.py').read_text();self.assertIn("e['component']==COMPONENTS[0] and e['n']==2",text)
    def test_fixed_rule_no_automatic_order_growth(self):
        text=(HERE/'S11c_d_defect_packet_contracted_numeric_lib.py').read_text();self.assertIn('for order in (24,48)',text);self.assertNotIn('order+=',text)
    def test_disk_backed_leaf_schema(self):
        text=(HERE/'S11c_d_defect_packet_contracted_numeric_lib.py').read_text();self.assertIn('CREATE TABLE contracted_leaves',text);self.assertIn('SELECT ordinal,receipt',text)
    def test_record_before_comparison(self):
        tree=ast.parse((HERE/'S11c_d_defect_packet_contracted_numeric.py').read_text());self.assertTrue(any(isinstance(n,ast.FunctionDef) and n.name=='compare' for n in ast.walk(tree)))
    def test_exact_tail_relations_converted_without_relaxing_predicate(self):
        tree=ast.parse((HERE/'S11c_d_defect_packet_contracted_prepare.py').read_text())
        messages={'same positive original tail constants','positive actual source/test bound','positive new tail domain','positive all-address/primitive tail allocation'}
        calls=[n for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='require' and len(n.args)>1 and isinstance(n.args[1],ast.Constant) and n.args[1].value in messages]
        self.assertEqual(len(calls),4)
        for n in calls:self.assertIsInstance(n.args[0],ast.Call);self.assertEqual(n.args[0].func.id,'bool')
    def test_symbolic_unknown_bool_refuses(self):
        class Unknown:
            def __bool__(self):raise TypeError('unknown exact relation')
        with self.assertRaises(TypeError):N.require(bool(Unknown()),'new symbolic relation')
    def test_all_sources_parse(self):
        for p in HERE.glob('S11c_d_defect_packet_contracted*.py'):ast.parse(p.read_text())


class BuildMetadata(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.manifest=json.loads((HERE/'S11c_d_defect_packet_contracted_numeric_inputs.json').read_text())
    def test_scope_does_not_claim_full_action(self):
        self.assertFalse(self.manifest['scope']['completePacketAction']);self.assertTrue(self.manifest['scope']['noLeakageOrCurrent'])
    def test_original_windows(self):self.assertEqual(self.manifest['scope']['originalWindowPairs'],[[27,122],[29,124]])
    def test_guard_memory_no_deadline(self):
        r=self.manifest['resources'];self.assertEqual((r['memoryBytes'],r['swapBytes'],r['tasksMax']),(4*1024**3,0,32));self.assertIsNone(r['durationLimits'])
    def test_every_input_complete_path_exists(self):
        for rec in self.manifest['savedInputs'].values():self.assertEqual(Path(rec['path']).stat().st_size,rec['bytes'])
    def test_native_imports_follow_containment(self):
        text=(HERE/'S11c_d_defect_packet_contracted_numeric.py').read_text()
        self.assertLess(text.index("enforced=ns['containment']()"),text.index('import sympy as sp'))
    def test_actual_command_scope(self):
        text=(HERE/'S11c_d_defect_packet_contracted_numeric_launch.py').read_text()
        self.assertIn("'defect_packet_contracted_numeric'",text);self.assertIn("'--pool','s11c-near-unity'",text)
    def test_hook_thread(self):
        self.assertIn('01a0e01b-ef84-7192-817f-584cda5d339b',(HERE/'S11c_d_defect_packet_contracted_numeric_launch.py').read_text())
    def test_original_branch_ast(self):
        old=ast.parse((HERE/'S11c_d_defect_packet_inner_lib.py').read_text());new=ast.parse((HERE/'S11c_d_defect_packet_contracted_numeric_lib.py').read_text())
        oldq=next(n for n in ast.walk(old) if isinstance(n,ast.FunctionDef) and n.name=='q')
        newq=next(n for n in ast.walk(new) if isinstance(n,ast.FunctionDef) and n.name=='q')
        a=next(n.value for n in oldq.body if isinstance(n,ast.Assign));b=next(n.value for n in newq.body if isinstance(n,ast.Assign));self.assertEqual(ast.dump(a),ast.dump(b))


if __name__=='__main__':unittest.main()
