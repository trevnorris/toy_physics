#!/usr/bin/env python3
"""Stdlib metadata and synthetic routing tests. No scientific restoration."""
import ast
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import runpy
import unittest

M=Path(__file__).resolve().parent;P='S11c_d_defect_packet_inner'
W=runpy.run_path(str(M/(P+'.py')),run_name='inert_worker_tests')
L=runpy.run_path(str(M/(P+'_lib.py')),run_name='inert_library_tests')
m=json.loads((M/(P+'_inputs.json')).read_text())
read=lambda p:json.loads(Path(p).read_text())
raw=lambda n:read(m['savedInputs'][n]['path'])

class FakeContext:
    @staticmethod
    def fsum(values):return sum(values,F(0))

class MetadataTests(unittest.TestCase):
    def test_declared_point_census(self):
        p=L['point_plan']();self.assertEqual(len(p),38);self.assertEqual(len({(x['k'],x['l']) for x in p}),38)
        self.assertTrue(all(F(x['k']) not in [-1,1] and F(x['l']) not in [-1,1] for x in p))
    def test_collisions_exact_and_both_sides(self):
        p=L['point_plan']()
        for label,target in [('opposite',0),('sum-plus',2),('sum-minus',-2),('profile',0)]:
            rows=[r for r in p if r['label']==label];self.assertEqual(len(rows),5)
            values=[F(r['l'])-F(r['k']) if label=='profile' else F(r['k'])+F(r['l']) for r in rows]
            self.assertIn(F(target),values);self.assertTrue(min(values)<target<max(values))
    def test_exact_root_label_merging(self):
        for point in L['point_plan']():
            if point['label'] in ['opposite','sum-plus','sum-minus'] and point['shift']=='0':
                cuts=L['exact_cuts'](point['k'],point['l'],122)
                self.assertTrue(any(any(v.startswith('height-root') for v in r['labels']) and any(v.startswith('reflected-root') for v in r['labels']) for r in cuts))
    def test_near_collisions_are_not_merged(self):
        cuts=L['exact_cuts']('1/3',str(F(-1,3)+F(1,128)),122)
        self.assertFalse(any(any(v.startswith('height-root') for v in r['labels']) and any(v.startswith('reflected-root') for v in r['labels']) for r in cuts))
    def test_profile_labels_and_radius_preserved(self):
        cuts=L['exact_cuts']('0','0',122);center=next(r for r in cuts if r['a']=='0' and r['b']=='0')
        self.assertIn('t=0',center['labels']);self.assertIn('t=Q',center['labels'])
        self.assertIn({'a':'0','b':'-122','labels':['left-radius']},cuts)
    def test_squared_orientation_constant_synthetic(self):
        cls=L['InnerEvaluator'];obj=object.__new__(cls);obj.rules={24:([F(-1,2),F(1,2)],[F(1),F(1)])};records=[]
        obj.emit=lambda c,n,v:records.append((n,v));c=FakeContext()
        result=obj.gauss('toy',c,[(F(-2),F(7))],lambda x:([F(1)],{'toyX':x}),24)
        self.assertEqual(result,[9]);self.assertEqual(len(records),3)
        for _,r in records[:-1]:
            self.assertTrue(all(-2<x<7 for x in r['points']));self.assertTrue(all(j>0 for j in r['jacobians']))
    def test_failed_panel_prefix_preserved_synthetic(self):
        cls=L['InnerEvaluator'];obj=object.__new__(cls);obj.rules={24:([F(-1,2),F(1,2)],[F(1),F(1)])};records=[];calls=[]
        obj.emit=lambda c,n,v:records.append((n,v))
        def failing(x):
            calls.append(x)
            if len(calls)==2:raise ValueError('synthetic')
            return [F(1)],{'toy':x}
        with self.assertRaises(ValueError):obj.gauss('toy',FakeContext(),[(F(0),F(1))],failing,24)
        self.assertTrue(records[-1][0].endswith('/failed-panel'));self.assertEqual(len(records[-1][1]['valuesPrefix']),1)
    def test_global_singular_error_oracle_terminates(self):
        cls=L['InnerEvaluator'];obj=object.__new__(cls);records=[]
        obj.emit=lambda c,n,v:records.append((n,v))
        # Manufactured bookkeeping oracle, not a quadrature rule or native kernel.
        def panel(key,c,a,b,f,label):
            width=b-a;depth=width.denominator.bit_length()-1
            error=F(1,2**(depth//2)) if a==0 else F(0)
            return {'a':a,'b':b,'K':[width,2*width],'error':[error,error/2]}
        obj.adaptive_panel=panel;budget=F(1,256)
        value=obj.adaptive('toy-root',FakeContext(),[(F(0),F(1))],None,budget)
        self.assertEqual(value,[1,2]);last=records[-1][1]
        self.assertTrue(all(e<=budget for e in last['summedEmpiricalErrors']))
        self.assertEqual(last['refinements'],16);self.assertEqual(len(last['activeLeaves']),17)
        self.assertTrue(any(n.endswith('/global-sum/16') for n,_ in records))
    def test_old_local_halving_cannot_meet_manufactured_root_budget(self):
        for depth in range(33):self.assertGreater(F(1,2**(depth//2)),F(1,256*2**depth))
    def test_global_vector_budget_and_final_partition(self):
        cls=L['InnerEvaluator'];obj=object.__new__(cls);records=[];evaluated={}
        obj.emit=lambda c,n,v:records.append((n,v))
        def panel(key,c,a,b,f,label):
            width=b-a;depth=width.denominator.bit_length()-1
            record={'a':a,'b':b,'K':[width],'error':[F(1,2**(depth//2)) if a==0 else F(0),width**2/100]}
            # Error component count and K component count must agree in real vector panels.
            record['K']=[width,width];evaluated[label]=record;return record
        obj.adaptive_panel=panel;budget=F(1,256)
        self.assertEqual(obj.adaptive('toy-vector',FakeContext(),[(F(0),F(1,2)),(F(1,2),F(1))],None,budget),[1,1])
        final=records[-1][1];active=[evaluated[n] for n in final['activeLeaves']];ordered=sorted(active,key=lambda v:v['a'])
        self.assertEqual(ordered[0]['a'],0);self.assertEqual(ordered[-1]['b'],1)
        self.assertTrue(all(a['b']==b['a'] for a,b in zip(ordered,ordered[1:])))
        sums=[sum((v['error'][j] for v in active),F(0)) for j in range(2)]
        self.assertEqual(sums,final['summedEmpiricalErrors']);self.assertTrue(all(v<=budget for v in sums))
    def test_global_zero_budget_refused(self):
        obj=object.__new__(L['InnerEvaluator'])
        with self.assertRaises(ValueError):obj.adaptive('toy',FakeContext(),[(F(0),F(1))],None,F(0))
    def test_subdivision_preserves_original_endpoints_under_rounding(self):
        with localcontext() as ctx:
            ctx.prec=10;a=Decimal('-1e20');b=Decimal('0.1');n=3
            self.assertNotEqual(a+(b-a)*n/n,b)
            points=L['partition_interval'](a,b,n)
            self.assertIs(points[0],a);self.assertIs(points[-1],b)
            self.assertTrue(all(x<y for x,y in zip(points,points[1:])))
    def test_subdivision_bad_inputs_refuse(self):
        for a,b,n in [(F(1),F(1),2),(F(2),F(1),2),(F(0),F(1),0),(F(0),F(1),True)]:
            with self.assertRaises(ValueError):L['partition_interval'](a,b,n)
    def test_native_profile_scale_actual_operand(self):
        r=raw('inventory/native-profile-scale-join.json');self.assertEqual(W['native_profile_scale'](r,'10'),10)
    def test_native_profile_scale_mutations_refuse(self):
        baseline=raw('inventory/native-profile-scale-join.json')
        for key,value in [('physicalLength','1'),('declaredLength',1),('savedLength',{'text':'1','srepr':'Integer(1)'})]:
            changed=json.loads(json.dumps(baseline));changed[key]=value
            with self.assertRaises(ValueError):W['native_profile_scale'](changed,'10')
        for key,value in [('rule','L_W**0'),('source',"value*=1"),('functionExecuted',True)]:
            changed=json.loads(json.dumps(baseline));changed['nativeRule'][key]=value
            with self.assertRaises(ValueError):W['native_profile_scale'](changed,'10')
    def test_require_exact_true(self):
        L['require'](True,'ok')
        for value in [1,None,False,'true']:
            with self.assertRaises(ValueError):L['require'](value,'refuse')
    def test_all_complete_saved_receipts(self):
        for alias,r in m['savedInputs'].items():
            with self.subTest(alias=alias):self.assertEqual(Path(r['path']).stat().st_size,r['bytes']);self.assertEqual(W['sha'](r['path']),r['sha256'])
    def test_all_source_pins(self):
        for path,h in m['sourcePins'].items():
            with self.subTest(path=path):self.assertEqual(W['sha'](path),h)
    def test_complete_original_selection(self):
        self.assertEqual(raw('selected/pressure-addresses.json')['selected'],[a for a in raw('inventory/THETA_BALANCE-ordered-addresses.json') if a['jet']['channel']=='e_W'])
    def test_every_original_factor_receipt(self):
        defs=raw('preflight/numeric-factor-adapters.json')['definitions']
        for a in raw('selected/pressure-addresses.json')['selected']:
            entry=defs['-'.join([a['face'],a['slot'],a['component']])];op=raw('factors/'+a['fullFactorProof']['proof']+'-operands.json')
            self.assertEqual(entry['original'],op['mappedAddressFactor'])
    def test_every_template_proof_operand(self):
        defs=raw('preflight/numeric-factor-adapters.json')['definitions'];self.assertEqual(len(defs),20)
        for label,e in defs.items():
            op=raw('preflight/numeric-factor-'+label+'-input.json');ret=raw('preflight/numeric-factor-'+label+'-return.json')
            self.assertEqual(op['left'],e['mapped']);self.assertEqual(op['right'],e['template']);self.assertEqual(ret['cancelled'],W['ZERO'])
    def test_native_profile_same_physical_input(self):self.assertEqual(raw('native/binding.json')['physicalInput'],raw('physical-input.json'))
    def test_all_whole_definitions_original_receipts(self):
        tags=raw('pressure/whole-tags.json')
        for name,alias in m['wholeDefinitionInputs'].items():self.assertEqual(tags[name]['savedDefinition'],raw(alias));self.assertEqual(tags[name]['sha256'],m['savedInputs'][alias]['sha256'])
    def test_rule_row_full_byte_receipts(self):
        ex=raw('rules/extraction.json');self.assertEqual(ex['database'],raw('fourier/journal-receipt.json'))
        for name,r in ex['records'].items():self.assertEqual(r['sha256'],m['savedInputs']['rules/'+name+'.json']['sha256']);self.assertEqual(r['bytes'],m['savedInputs']['rules/'+name+'.json']['bytes'])
    def test_rule_census_metadata(self):
        for name,n in [('A-GL24',24),('A-GL48',48)]:
            r=raw('rules/'+name+'.json');self.assertEqual(len(r['nodes']),n);self.assertEqual(len(r['weights']),n);self.assertEqual(len(r['momentResiduals']),2*n)
        r=raw('rules/B-G7-K15.json');self.assertEqual(len(r['kronrodNodes']),15);self.assertEqual(len(r['gaussNodes']),7);self.assertEqual(len(r['momentResidualsThrough23']),24)
    def test_authority_one_no_deadline(self):
        a=read(m['executionAuthority']);self.assertEqual(a['scope'],m['scope']);self.assertEqual(a['scienceExecutionsAuthorized'],1);self.assertIs(a['noDeadline'],True);self.assertIs(a['automaticScientificRetry'],False)

class SourceTests(unittest.TestCase):
    def test_inert_top_imports(self):
        allowed={'argparse','ast','hashlib','itertools','json','math','os','pathlib','resource','shutil','sys','time','traceback','fractions','heapq'}
        for name in [P+'.py',P+'_lib.py']:
            for n in ast.parse((M/name).read_text()).body:
                if isinstance(n,ast.Import):self.assertTrue(all(a.name in allowed for a in n.names))
                if isinstance(n,ast.ImportFrom):self.assertIn(n.module,allowed)
    def test_no_completed_numerical_constructor_or_producer(self):
        for name in [P+'.py',P+'_lib.py']:
            calls={ast.unparse(n.func) for n in ast.walk(ast.parse((M/name).read_text())) if isinstance(n,ast.Call)}
            self.assertFalse(calls&{'FourierEvaluator','oldlib.FourierEvaluator','kronrod15','c.gauss_quadrature','sp.Integral','build_face','signal.alarm','signal.setitimer','time.sleep'})
    def test_science_import_after_containment(self):
        source=(M/(P+'.py')).read_text();self.assertLess(source.index("ns['containment']()"),source.index('import sympy'))
    def test_journal_finalized_on_failure(self):
        source=(M/(P+'.py')).read_text();self.assertIn("finally:\n        store.close();J.emit('numerical-journal-receipt'",source)
    def test_fixed_rule_reuse_not_reconstruction(self):
        tree=ast.parse((M/(P+'_lib.py')).read_text());s=ast.unparse(next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='__init__'))
        self.assertIn('self.restore',s);self.assertNotIn('gauss_quadrature',s);self.assertNotIn('kronrod15(',s)
    def test_ready_gate_not_present(self):self.assertFalse((M/(P+'_gate.json')).exists())
    def test_pool_and_hook_no_computation_timer(self):
        source=(M/(P+'_launch.py')).read_text();self.assertIn("'--pool','s11c-near-unity'",source);self.assertIn("'--memory-gib','4'",source);self.assertIn('completion-watcher/state.json',source)
    def test_no_scientific_stage_asserted_by_metadata(self):
        self.assertFalse(m['completedFunctionsReplayed']);self.assertEqual(m['status'],'PREPARED_FOR_BUILD_REVIEW_NOT_EXECUTED')

if __name__=='__main__':unittest.main()
