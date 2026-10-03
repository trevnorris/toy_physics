#!/usr/bin/env python3
"""Standard-library source/metadata tests. Never restore scientific payloads."""
import ast
import copy
import hashlib
import json
from pathlib import Path
import runpy
import tempfile
from types import SimpleNamespace
import unittest

HERE=Path(__file__).resolve().parent
PREFIX='S11c_d_defect_packet_preflight'
SOURCE=HERE/(PREFIX+'.py')
NS=runpy.run_path(str(SOURCE),run_name='source_tests_only')
TREE=ast.parse(SOURCE.read_text())
MANIFEST=json.loads((HERE/(PREFIX+'_inputs.json')).read_text())
ZERO={'text':'0','srepr':'Integer(0)'}

def read(alias):
    return json.loads(Path(MANIFEST['savedInputs'][alias]['path']).read_text())

def top(name):
    return next(n for n in TREE.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name==name)

class Tooling(unittest.TestCase):
    def test_no_scientific_import_on_load(self):
        for n in TREE.body:
            if isinstance(n,ast.Import):self.assertFalse(any(v.name.split('.')[0] in ('sympy','numpy','scipy','mpmath') for v in n.names))
            if isinstance(n,ast.ImportFrom):self.assertFalse(n.module.split('.')[0] in ('sympy','numpy','scipy','mpmath'))
        self.assertNotIn('sp',NS)
    def test_containment_precedes_scientific_import(self):
        text=ast.get_source_segment(SOURCE.read_text(),top('main'))
        self.assertLess(text.index("ns['containment']()"),text.index('import sympy'))
        self.assertLess(text.index('verify_invocation('),text.index('args.out.mkdir'))
    def test_predicate_rejects_truthiness(self):
        class Poison:
            def __bool__(self):raise AssertionError('Must not coerce unknown')
        for x in (False,None,0,1,'True',Poison()):
            with self.assertRaises(ValueError):NS['require'](x,'refused')
        NS['require'](True,'accepted')
    def test_exact_true_atom_interface(self):
        true_atom=object();g=NS['require'].__globals__;before=g.get('sp')
        try:
            g['sp']=SimpleNamespace(S=SimpleNamespace(true=true_atom))
            NS['require'](true_atom,'known exact atom')
            with self.assertRaises(ValueError):NS['require'](object(),'unknown object')
        finally:
            if before is None:g.pop('sp',None)
            else:g['sp']=before
    def test_selection_actual_originals(self):
        a=read('inventory/THETA_BALANCE-ordered-addresses.json');c=read('local/all-local-cells.json')
        selected,local,status=NS['validate_selection'](a,c)
        self.assertEqual(selected,read('selected/pressure-addresses.json')['selected'])
        self.assertEqual(local,read('selected/local-cells.json')['selected'])
        self.assertEqual(sum(status.values()),544)
    def test_selection_rejects_duplicate(self):
        a=read('inventory/THETA_BALANCE-ordered-addresses.json');c=read('local/all-local-cells.json')
        ids=[i for i,v in enumerate(a) if v['jet']['channel']=='e_W'];a[ids[1]]=copy.deepcopy(a[ids[0]])
        with self.assertRaises(ValueError):NS['validate_selection'](a,c)
    def test_selection_rejects_grade_corruption(self):
        a=read('inventory/THETA_BALANCE-ordered-addresses.json');c=read('local/all-local-cells.json')
        next(v for v in a if v['jet']['channel']=='e_W')['targetGrade']=[2,0]
        with self.assertRaises(ValueError):NS['validate_selection'](a,c)
    def test_selection_rejects_missing_local(self):
        a=read('inventory/THETA_BALANCE-ordered-addresses.json');c=read('local/all-local-cells.json')
        c.remove(next(v for v in c if v['row']=='THETA_BALANCE' and v['field']=='e_W'))
        with self.assertRaises(ValueError):NS['validate_selection'](a,c)
    def test_exact_all_input_byte_receipts(self):
        for alias,r in MANIFEST['savedInputs'].items():
            p=Path(r['path']);self.assertEqual(NS['sha'](p),r['sha256'],alias);self.assertEqual(p.stat().st_size,r['bytes'])
    def test_current_source_pins(self):
        for p,h in MANIFEST['sourcePins'].items():self.assertEqual(NS['sha'](p),h,p)
    def test_whole_definitions_are_actual_full_files(self):
        tags=read('pressure/whole-tags.json')
        for tag,alias in MANIFEST['wholeDefinitionInputs'].items():
            self.assertEqual(tags[tag]['savedDefinition'],read(alias));self.assertEqual(tags[tag]['sha256'],MANIFEST['savedInputs'][alias]['sha256'])
    def test_inherited_proof_operand_schemas(self):
        a=read('selected/pressure-addresses.json')['selected']
        for x in a:
            name=x['fullFactorProof']['proof'];op=read('factors/'+name+'-operands.json')
            for k in ('sha256','original','mapped','map','symbolAssumptions','flatSupport','frequency','positiveRegulatorContinuation'):
                self.assertEqual(op['actualResponseMap'][k],x['responseMap'][k])
            self.assertEqual(op['addressNormalOriginal'],x['normalOriginal']);self.assertEqual(op['requiredMap'],x['fullFactorProof']['completeNormalMap'])
            inp=read('factors/'+name+'-full-mapped-residual-input.json');self.assertEqual(inp['right'],op['mappedAddressFactor'])
            inp=read('factors/'+name+'-normal-source-join-input.json');self.assertEqual(inp['left'],op['addressNormalOriginal']);self.assertEqual(inp['right'],op['savedNormal'])
            for t in ('-full-mapped-residual-return.json','-normal-source-join-return.json'):self.assertEqual(read('factors/'+name+t)['cancelled'],ZERO)
    def test_field_operand_and_quotient_schema(self):
        fields=read('pressure/fields.json');certs=read('pressure/coefficient-certificates.json')
        self.assertEqual(set(fields),set(certs));self.assertEqual(len(fields),34)
        for fid,v in fields.items():
            p=read('field/'+fid+'-polynomial.json');i=read('field/'+fid+'-reconstruction-input.json');r=read('field/'+fid+'-reconstruction-return.json')
            self.assertEqual(i['left'],v['field']);self.assertEqual(r['cancelled'],ZERO)
            self.assertEqual(len(p['coefficients']),p['degree']+1);self.assertIn('denominator',p)
    def test_local_proof_and_ancestry_not_replayed(self):
        for c in read('selected/local-cells.json')['selected']:
            self.assertEqual([i['name'] for i in c['identities']],['cell-source-sum','cell-physical-polynomial','cell-derivative'])
            self.assertTrue(all(i['cancelled']==ZERO for i in c['identities']))
            self.assertEqual(c['identities'][0]['left'],c['coefficient']);self.assertEqual(c['polynomial']['original'],c['coefficient'])
            self.assertEqual(len(c['sourceChildren']),len(c['summands']))
        self.assertNotIn("sum(v['summands'])",SOURCE.read_text())
    def test_physical_schema(self):
        c=read('local/context.json');p=read('physical-input.json')
        self.assertEqual(c['physical'],p);self.assertEqual(c['saved']['physicalInput'],p)
        self.assertEqual(c['frequencyOverride'],{'old':'1','actual':3});self.assertEqual(c['numeric']['omega'],{'text':'3','srepr':'Integer(3)'})
        self.assertTrue(c['effectiveSpeedOnlyInPressure'])
        names=[p[0]['text'] for p in read('ends/left-match.json')['point']['mappingPairs']]
        self.assertEqual(set(names),{'weak_end_p','weak_end_cs','weak_end_q'})
    def test_control_candidates_and_normal_height_zeros(self):
        a=read('selected/pressure-addresses.json')['selected'];f=read('pressure/fields.json')
        height=[v for v in a if v['slot']=='normal' and v['component']=='NATIVE_HEIGHT']
        self.assertEqual(len(height),48);self.assertTrue(all(v['status']=='EXACT_ZERO_CONSUMER' for v in height))
        live=[v for v in a if v['status']=='FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED']
        self.assertTrue(any(v['component']=='NATIVE_MIXED_ITERATION' for v in live))
        self.assertTrue(any(v['component']=='INHERITED_DIRECT_WHOLE_OFF_DIAGONAL' for v in live))
        self.assertTrue(any(v['slot']=='normal' and v['component']=='NATIVE_SLOPE' and not f[v['consumerTransform']['coefficientId']]['constant'] for v in live))
        self.assertTrue(any(v['jet']['spatialOrders'][0]>0 and not f[v['sourceTransform']['coefficientId']]['constant'] for v in live))
    def test_readiness_fails_before_output_or_import(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'gate.json';p.write_text('{"status":"NOT_READY"}')
            with self.assertRaises(ValueError):NS['verify_gate'](p,HERE/(PREFIX+'_inputs.json'),MANIFEST)
            self.assertEqual(sorted(x.name for x in Path(t).iterdir()),['gate.json'])
    def test_output_and_argv_route_refusal(self):
        a=SimpleNamespace(out=Path('/tmp/x'),inputs=Path('/tmp/i'),gate=Path('/tmp/g'))
        argv=[str(SOURCE),'--out','/tmp/x','--inputs','/tmp/i','--gate','/tmp/g']
        NS['verify_invocation'](a,{'command':argv,'outputDirectory':'/tmp/x'},argv)
        with self.assertRaises(ValueError):NS['verify_invocation'](a,{'command':argv,'outputDirectory':'/tmp/y'},argv)
        with self.assertRaises(ValueError):NS['verify_invocation'](a,{'command':argv,'outputDirectory':'/tmp/x'},argv+['extra'])
    def test_evidence_exclusive_and_durable_interface(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'evidence.json';NS['save'](p,{'old':1});before=p.read_bytes()
            with self.assertRaises(FileExistsError):NS['save'](p,{'new':2})
            self.assertEqual(p.read_bytes(),before)
    def test_partial_tail_record_precedes_capacity_refusal(self):
        text=ast.get_source_segment(SOURCE.read_text(),top('run_science'))
        self.assertLess(text.index("J.emit('tail-selection-step-'"),text.index("require(K<=m['preflightCapacity']"))
        self.assertIn("'quadratureReady':False",text);self.assertIn("'numericalAction':None",text)
    def test_launcher_guard_and_hook_routes(self):
        t=(HERE/(PREFIX+'_launch.py')).read_text();ast.parse(t)
        self.assertIn("READY_FOR_ONE_PACKET_PREFLIGHT",t);self.assertIn("'defect_packet_preflight'",t)
        self.assertIn("--memory-gib','4'",t);self.assertIn('s11c_guarded_run.py',t);self.assertIn('S11c_d_end_normalization_run.py',t)
        self.assertLess(t.index("path = RUN / 'completion-watcher/state.json'"),t.index("os.write(writer, b'1')"))
    def test_no_computational_deadlines_or_integral_calls(self):
        forbidden={'quad','quadgl','quadts','integrate','Integral','solve','alarm','setitimer','sleep'}
        for n in ast.walk(TREE):
            if isinstance(n,ast.Call):
                name=n.func.id if isinstance(n.func,ast.Name) else n.func.attr if isinstance(n.func,ast.Attribute) else None
                self.assertNotIn(name,forbidden)
        self.assertIsNone(MANIFEST['resources']['durationLimits']);self.assertFalse((HERE/(PREFIX+'_gate.json')).exists())

if __name__=='__main__':unittest.main()
