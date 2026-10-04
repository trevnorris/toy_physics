"""Finite stdlib tooling and saved-operand metadata tests. No scientific restoration."""
import ast,importlib.util,json,sys,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
HERE=Path(__file__).resolve().parent;PREFIX='S11c_d_first_order_receiving_regular'
spec=importlib.util.spec_from_file_location('regular_build',HERE/(PREFIX+'.py'));w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
def manifest():return json.loads((HERE/(PREFIX+'_inputs.json')).read_text())
def raw(alias):return json.loads(Path(manifest()['savedFiles'][alias]['path']).read_text())
class Tooling(unittest.TestCase):
    def test_import_and_compile_without_science(self):
        self.assertNotIn('sympy',sys.modules)
        for suffix in ['.py','_launch.py']:compile((HERE/(PREFIX+suffix)).read_text(),PREFIX+suffix,'exec')
    def test_strict_predicate(self):
        w.require(True,'ok')
        for v in [1,0,None,False,[],[True],'true']:
            with self.assertRaises(ValueError):w.require(v,'refuse')
    def test_variations_with_internal_zeros(self):
        for signs,count in [([],0),([1,0,1],0),([1,0,-1,0,1],2),([-1,-1,0,1],1),([0,0],0)]:self.assertEqual(w.variations(signs),count)
    def test_variations_refuses_nonexact_signs(self):
        for signs in [[True],[1.0],[None],[2],[-2]]:
            with self.assertRaises(ValueError):w.variations(signs)
    def test_unknown_sign_stops_without_numeric_fallback(self):
        for z,p,n,expected in [(True,False,False,0),(False,True,False,1),(False,False,True,-1)]:
            self.assertEqual(w.exact_sign(SimpleNamespace(is_zero=z,is_positive=p,is_negative=n)),expected)
        with self.assertRaises(ValueError):w.exact_sign(SimpleNamespace(is_zero=None,is_positive=None,is_negative=None))
    def test_exact_invocation(self):
        a=SimpleNamespace(out=Path('/tmp/regular-out'),inputs=Path('/tmp/regular-input'),gate=Path('/tmp/regular-gate'))
        argv=[str(Path(w.__file__).resolve()),'--out',str(a.out),'--inputs',str(a.inputs),'--gate',str(a.gate)];g={'command':['guard']+argv,'outputDirectory':str(a.out)}
        w.verify_invocation(a,g,argv)
        with self.assertRaises(ValueError):w.verify_invocation(a,g,argv+['retry'])
        with self.assertRaises(ValueError):w.verify_invocation(a,{**g,'outputDirectory':'/tmp/else'},argv)
    def test_append_only(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'return';w.save(p,{'value':1})
            with self.assertRaises(FileExistsError):w.save(p,{'value':2})
            self.assertEqual(json.loads(p.read_text()),{'value':1})
    def test_guard_before_import_and_creation(self):
        t=ast.parse((HERE/(PREFIX+'.py')).read_text());main=ast.unparse(next(n for n in t.body if isinstance(n,ast.FunctionDef) and n.name=='main'))
        self.assertLess(main.index('verify_gate('),main.index('args.out.mkdir'))
        self.assertLess(main.index("ns['containment']()"),main.index('import sympy'))
    def test_no_deadlines_or_old_scientific_function_calls(self):
        t=ast.parse((HERE/(PREFIX+'.py')).read_text());calls={ast.unparse(n.func) for n in ast.walk(t) if isinstance(n,ast.Call)}
        self.assertFalse(calls&{'signal.alarm','signal.setitimer','sp.integrate','sp.limit','sp.solve','sp.nroots','sp.N','pickle.loads','subprocess.run'})
        imports=[ast.unparse(n) for n in ast.walk(t) if isinstance(n,(ast.Import,ast.ImportFrom))]
        self.assertFalse(any('S11c_d_' in n for n in imports))
    def test_original_saved_bytes(self):
        for alias,v in manifest()['savedFiles'].items():
            self.assertEqual(Path(v['path']).stat().st_size,v['bytes'],alias);self.assertEqual(w.sha(v['path']),v['sha256'],alias)
    def test_pressure_joins_have_complete_native_rows_and_faces(self):
        assembly=raw('source/first-order-pressure-assembly.json');receiving=raw('receiving/pressure-receiving-assembly.json')
        self.assertEqual(len(assembly),10);self.assertEqual(len(receiving['entries']),20)
        for row in w.ROWS:
            for col in [0,1]:
                recs=[v for v in assembly if v['row']==row and v['column']==col];self.assertEqual(len(recs),1)
                self.assertEqual({(v['face'],v['slot']) for v in recs[0]['pieces']},{(f,s) for f in ['plus','minus'] for s in ['pressure','normal']})
                for piece in recs[0]['pieces']:
                    factor=raw('source/'+piece['face']+'-'+piece['slot']+'-flat-arguments.json')
                    self.assertEqual(factor['factor'],piece['flatNormalFactor'])
                    self.assertEqual(piece['source01Envelope'],raw('source/'+piece['face']+'-source-01-column-'+str(col)+'-return.json')['value'])
    def test_all_local_force_ancestry_is_supplied(self):
        cells=raw('source-input/local/all-local-cells.json');self.assertEqual(len(cells),400)
        for row in w.ROWS:
            for grade in ['10','01']:
                for col in [0,1]:
                    op=raw('source/'+row+'-'+grade+'-column-'+str(col)+'-input.json')
                    self.assertEqual(len(op['cells']),20)
                    for v in op['cells']:self.assertIn(v['cell'],cells)
    def test_inherited_matrix_proof_census_and_arguments(self):
        for name,n,m,folder in [('geometric-dual',5,5,'receiving'),('transverse-to-scalars-coupling',3,2,'receiving'),('scalars-to-transverse-coupling',2,3,'receiving'),('transverse-block',2,2,'receiving'),('incident-chart-match',5,2,'receiving'),('FULL-offwave-source-correspondence',5,5,'receiving'),('LEFT-raw-vs-receiving',5,5,'receiving'),('native-baseline-current-join',5,5,'flux'),('LEFT-slab-epsilon-once',5,5,'flux')]:
            op=raw(folder+'/'+name+'-matrix-operands.json');self.assertEqual(set(op),{'left','right'})
            for i in range(n):
                for j in range(m):
                    v=raw(folder+'/'+name+'-'+str(i)+'-'+str(j)+'-return.json');self.assertEqual(v['cancelled']['srepr'],'Integer(0)')
                    self.assertEqual(set(raw(folder+'/'+name+'-'+str(i)+'-'+str(j)+'-input.json')),{'left','right'})
    def test_original_profile_returns_and_length(self):
        profiles=raw('source/restored-profile-rules.json');self.assertEqual(len(profiles['certificates']),40)
        for v in profiles['certificates']:
            self.assertEqual(v['L']['srepr'],'Integer(10)')
            for proof in v['identities']:self.assertEqual(proof['cancelled']['srepr'],'Integer(0)')
    def test_original_current_source_selection(self):
        native=raw('original/native/LEFT-original.json')['actual']['slab']['SLAB_CURRENT_MATRIX']
        self.assertEqual(native,raw('flux/LEFT-native-current-source.json')['slabCurrent'])
        self.assertEqual(native,raw('flux/LEFT-slab-binding.json')['original'])
        self.assertEqual(raw('flux/native-baseline-current-join-matrix-operands.json')['right'],raw('flux/physical-linear-survival-operands.json')['leftUnboundCurrent'])
    def test_original_both_threshold_arguments_and_block(self):
        a=raw('receiving/grazing-minus-determinant-and-domains.json');b=raw('receiving/grazing-plus-determinant-and-domains.json')
        self.assertEqual(a['originalBlock'],b['originalBlock']);self.assertEqual(a['determinant'],b['determinant'])
        self.assertNotEqual(a['point'][0],b['point'][0]);self.assertEqual(a['point'][1],b['point'][1])
    def test_full_native_affine_map_origins(self):
        native=raw('original/native/LEFT-original.json')['actual']['acoustic']
        self.assertEqual(native,raw('receiving/native-face-original-operands.json'))
        for face in ['plus','minus']:
            for col in [0,1]:
                a=raw('receiving/'+face+'-affine-face-column-'+str(col)+'.json')
                self.assertIn(a['nativeFace'],native['FACE_RECORDS'])
                self.assertEqual(a['inheritedSourceNormalization']['operands']['left'],raw('source/'+face+'-source-01-column-'+str(col)+'-return.json')['value'])
                self.assertEqual(a['inheritedSourceNormalization']['return']['cancelled']['srepr'],'Integer(0)')
    def test_actual_source_not_relabelled_as_loss(self):
        for col in [0,1]:
            f=raw('source/full-local-force-column-'+str(col)+'.json');self.assertFalse(f['diagnosticIsPhysicalLossProjection'])
            self.assertTrue(f['pressureAddition'].startswith('negative C00 R00 S01'))
        self.assertEqual(raw('source/incident-columns.json')['fieldOrder'],list(w.FIELDS))
    def test_staged_adjugate_and_no_field_solver(self):
        s=(HERE/(PREFIX+'.py')).read_text()
        self.assertLess(s.index("certificates=[exclude_ray"),s.index('adj=Cq.adjugate'))
        self.assertIn('REAL_AXIS_RECEIVING_POLE_OR_ENDPOINT_ROOT',s)
        self.assertIn('UNRESOLVED_EXACT_SIGN',s)
        self.assertIn('chart_extra=1',s)
        self.assertNotIn('.inv(',s)
    def test_method_build_execution_separation(self):
        m=manifest();r=json.loads(Path(m['methodRecord']).read_text());self.assertTrue(r['methodAssessed'])
        self.assertEqual(r['literalVerdict'],'CLEAR FOR THIS REAL-AXIS RECEIVING AND CROSS-FLUX METHOD')
        self.assertTrue(m['scope']['physicalPowerBalancePending']);self.assertIsNone(m['scope']['leakageFactor'])
        self.assertFalse((HERE/(PREFIX+'_gate.json')).exists())
if __name__=='__main__':unittest.main()
