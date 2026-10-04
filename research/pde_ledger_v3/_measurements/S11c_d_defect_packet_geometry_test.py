"""Standard-library tooling tests. All computed geometry is SYNTHETIC Q(sqrt(2)).

No actual physical plan, scientific symbol restoration or quadrature is run.
"""
import ast
import copy
from fractions import Fraction as F
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

M=Path(__file__).resolve().parent

def load(name):
    spec=importlib.util.spec_from_file_location(name,M/(name+'.py'));m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m

G=load('S11c_d_defect_packet_geometry_lib');W=load('S11c_d_defect_packet_geometry')

def q(a=0,b=0):return G.Quad(a,b,2)

def fixture():
    lines=[G.Line(0,q(-3),('bottom',)),G.Line(0,q(3),('top',)),
           G.Line(1,q(),('diagonal',)),G.Line(-1,q(),('anti',)),
           G.Line(0,q(0,1),('root',)),G.Line(1,q(),('duplicate-geometry',))]
    cuts=[(q(0,-1),'vertical-root'),(q(),'zero')]
    return lines,cuts,[q(-3),q(3),q(-3),q(3)]


class ExactTests(unittest.TestCase):
    def test_arithmetic(self):
        self.assertEqual(q(1,1)*q(1,-1),-1);self.assertEqual(q(1,1)/q(1,-1),q(-3,-2))
        self.assertEqual(q(4,3)+q(-4,-3),0);self.assertEqual(q(2)/2,1)
    def test_order_mixed_signs(self):
        self.assertTrue(q(1,-1)<0<q(-1,1));self.assertTrue(q(3,-2)>0);self.assertTrue(q(-3,2)<0)
        self.assertTrue(q(0,-1)<q()<q(0,1));self.assertEqual(sorted([q(0,1),q(1),q()]),[q(),q(1),q(0,1)])
    def test_hash_contract(self):self.assertEqual(hash(q(2)),hash(2));self.assertEqual({q(2):'a'}[2],'a')
    def test_no_floats(self):
        for v in [1.0,True,float('nan'),float('inf')]:
            with self.assertRaises(ValueError):G.Quad(v,0,2)
    def test_no_bad_field(self):
        for d in [0,-2,4,'9/16']:
            with self.assertRaises(ValueError):G.Quad(0,1,d)
    def test_field_mismatch(self):
        with self.assertRaises(ValueError):q(1)+G.Quad(1,0,3)
    def test_zero_divisor(self):
        with self.assertRaises(ValueError):q(1)/q()
    def test_line_type(self):
        with self.assertRaises(ValueError):G.Line(.5,q(),('x',))
    def test_duplicate_label(self):
        with self.assertRaises(ValueError):G.canonical_lines([G.Line(0,q(),('x',)),G.Line(1,q(),('x',))])
    def test_coalesced_labels(self):
        lines,_,_=fixture();actual=G.canonical_lines(lines)
        self.assertEqual(len(actual),5);self.assertIn(('diagonal','duplicate-geometry'),[v.labels for v in actual])
    def test_exact_coverage(self):
        lines,cuts,box=fixture();p=G.arrangement(lines,*box,cuts);a=G.audit(p,lines,cuts)
        self.assertEqual(a['area'],36);self.assertTrue(a['disjointInteriors']);self.assertEqual(a['numericalNodesEvaluated'],0)
    def test_missing_line(self):
        lines,cuts,box=fixture();bad=[l for l in lines if l.labels!=('anti',)];p=G.arrangement(bad,*box,cuts)
        self.assertEqual(G.audit(p,bad,cuts)['area'],36)
        with self.assertRaisesRegex(ValueError,'required collision'):G.audit(p,lines,cuts)
    def test_reversed_cell(self):
        lines,cuts,box=fixture();p=G.arrangement(lines,*box,cuts);c=p['slabs'][0]['cells'][0]
        c['lowerLine'],c['upperLine']=c['upperLine'],c['lowerLine']
        with self.assertRaisesRegex(ValueError,'oriented adjacent'):G.audit(p,lines,cuts)
    def test_missing_cell(self):
        lines,cuts,box=fixture();p=G.arrangement(lines,*box,cuts);p['slabs'][0]['cells'].pop()
        with self.assertRaisesRegex(ValueError,'cell census'):G.audit(p,lines,cuts)
    def test_changed_area(self):
        lines,cuts,box=fixture();p=G.arrangement(lines,*box,cuts);p['slabs'][0]['cells'][0]['area']+=1
        with self.assertRaisesRegex(ValueError,'cell area'):G.audit(p,lines,cuts)
    def test_missing_crossing_record(self):
        lines,cuts,box=fixture();p=G.arrangement(lines,*box,cuts);p['intersections'].pop()
        with self.assertRaisesRegex(ValueError,'crossing provenance'):G.audit(p,lines,cuts)
    def test_missing_crossing_label(self):
        lines,cuts,box=fixture();p=G.arrangement(lines,*box,cuts)
        r=next(v for v in p['cutLabels'] if any(t.startswith('intersection:') for t in v['labels']))
        r['labels']=[t for t in r['labels'] if not t.startswith('intersection:')]
        with self.assertRaisesRegex(ValueError,'crossing label'):G.audit(p,lines,cuts)
    def test_synthetic_full_specs(self):
        specs=G.specifications(q(0,1),4,7,q(0,1),F(2),F(3))
        for name,(lines,cuts,box) in specs.items():
            audit=G.audit(G.arrangement(lines,*box,cuts),lines,cuts)
            self.assertEqual(audit['area'],64 if name=='square' else 56)
    def test_height_both_signs(self):
        lines,_,_=G.specifications(q(0,1),4,7,q(),F(2),F(3))['height']
        for v in ['branch:l-','branch:l+']:
            pair=[l for l in lines if any(x.endswith(v) for x in l.labels)]
            self.assertEqual(sorted(l.slope for l in pair),[-1,1]);self.assertEqual(len(pair),2)
    def test_quadrature_occurrence_count(self):
        c=G.static_counts({'cells':7,'slabs':3})
        self.assertEqual(c['A24']['outerNodeOccurrences'],7*48**2)
        self.assertEqual(c['A48']['outerNodeOccurrences'],7*96**2)
        self.assertEqual(c['B_initial']['outerNodeOccurrencesLowerBound'],7*225)
        self.assertEqual(c['oldBankRequestsReused'],0)
    def test_exact_json(self):
        self.assertEqual(W.packed({'x':q(1,2),'f':F(2,3)}),{'x':{'a':'1','b':'2','basisSquare':'2'},'f':'2/3'})
        self.assertEqual(W.packed(2.0),2.0)  # old elapsed-time metadata only
        with self.assertRaises(ValueError):W.packed(float('nan'))
    def test_journal_before_return(self):
        with tempfile.TemporaryDirectory() as d:
            j=W.Journal(Path(d));j.start('synthetic',{'x':q(1,1)});self.assertEqual(j.active,'synthetic')
            j.finish({'x':q(2,2)});rows=[json.loads(s) for s in (Path(d)/'evidence-chain.jsonl').read_text().splitlines()]
            self.assertEqual(rows[1]['previous'],rows[0]['chainSha256']);self.assertEqual(j.completed,['synthetic'])
            with self.assertRaises(FileExistsError):j.emit('synthetic-input',{})
    def test_no_evaluators_imported(self):
        tree=ast.parse((M/'S11c_d_defect_packet_geometry.py').read_text())
        imports=[a.name for n in ast.walk(tree) if isinstance(n,(ast.Import,ast.ImportFrom)) for a in n.names]
        self.assertFalse(set(imports)&{'sympy','mpmath','numpy','scipy'})
    def test_no_unprotected_entrypoint(self):
        source=(M/'S11c_d_defect_packet_geometry.py').read_text()
        self.assertLess(source.index("gate=verify_gate(args.gate"),source.index("G=load_geometry(m)"))
        self.assertLess(source.index("enforced=ns['containment']()"),source.index("G=load_geometry(m)"))
    def test_standing_scope(self):
        source=(M/'S11c_d_defect_packet_geometry.py').read_text()
        self.assertIn("'pressureSummandUnits':'PENDING_SEPARATE_REQUIRED_CERTIFICATE'",source)
        self.assertIn("'numericalEvaluatorReady':False",source)
    def test_launcher_no_deadline(self):
        source=(M/'S11c_d_defect_packet_geometry_launch.py').read_text()
        self.assertIn("'--stage','defect_packet_geometry'",source);self.assertNotIn('timeout=',source)
        self.assertIn("'--pool','s11c-near-unity','--memory-gib','4'",source)
    def test_manifest_saved_receipts(self):
        p=M/'S11c_d_defect_packet_geometry_inputs.json'
        if not p.exists():self.skipTest('manifest prepared after initial synthetic tests')
        m=json.loads(p.read_text());self.assertFalse(m['completedFunctionsReplayed'])
        for receipt in m['savedInputs'].values():
            self.assertEqual(W.sha(receipt['path']),receipt['sha256']);self.assertEqual(Path(receipt['path']).stat().st_size,receipt['bytes'])
    def test_saved_rule_counts_only(self):
        p=M/'S11c_d_defect_packet_geometry_inputs.json'
        if not p.exists():self.skipTest('manifest prepared after initial synthetic tests')
        m=json.loads(p.read_text())
        for name,n in [('A-GL24',24),('A-GL48',48),('B-G7-K15',15)]:
            v=json.loads(Path(m['savedInputs']['rules/'+name+'.json']['path']).read_text())
            self.assertEqual(len(v['kronrodNodes' if n==15 else 'nodes']),n)


    def test_paired_branch_work_counts(self):
        entries=[{'addressId':i,'component':'NATIVE_HEIGHT','status':W.LIVE,
                  'targetGrade':[1,0],'families':{'X':'x','Y':'y','response':'r'}} for i in range(2)]
        result=W.count_work(entries,{}, {'slabs':3,'cells':7},'height',G)
        record=result['primitives'][0]
        for route in record['routes'].values():
            nodes=route['nodeOccurrencesPerFamily']
            self.assertEqual(route['addressProductOccurrences'],nodes*2)
            self.assertEqual(route['addressBranchProductOccurrences'],nodes*4)
            self.assertEqual(route['candidateResponseBranchOccurrences'],nodes*2)
            self.assertEqual(route['candidateYBranchOccurrencesBeforeCoordinateReuse'],nodes*2)
            self.assertEqual(route['candidateXFamilyOccurrencesBeforeCoordinateReuse'],nodes)
        self.assertEqual(record['addressedBranchScalarReturnOccurrencesA24A48WithoutSharing'],
                         2*record['addressedScalarReturnOccurrencesA24A48WithoutSharing'])
    def test_flat_counts_single_outer_integral(self):
        e={'addressId':0,'component':'NATIVE_FLAT','status':W.LIVE,'targetGrade':[0,0],
           'families':{'X':'x','Y':'y','response':'r'}}
        r=W.count_work([e],{}, {'slabs':3,'cells':7},'square',G)['primitives'][0]
        self.assertEqual(r['routes']['A24']['nodeOccurrencesPerFamily'],3*48)
        self.assertEqual(r['routes']['A24']['addressBranchProductOccurrences'],3*48)
    def test_empty_primitive_counts(self):
        for r in W.count_work([],{}, {'slabs':3,'cells':7},'height',G)['primitives']:
            self.assertEqual(r['addressedBranchScalarReturnOccurrencesA24A48WithoutSharing'],0)
            self.assertEqual(r['addressedScalarReturnOccurrencesA24A48WithoutSharing'],0)
    def test_wrong_box_explicit_join_refuses(self):
        lines,cuts,box=fixture();actual_box=[q(-4),q(4),box[2],box[3]]
        plan=G.arrangement(lines,*actual_box,cuts)
        self.assertTrue(G.audit(plan,lines,cuts)['coverage'])
        with self.assertRaisesRegex(ValueError,'actual arrangement box'):
            W.require(plan['box']==box,'actual arrangement box matches specification')
    def test_missing_saved_directory_scan_is_empty(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(list((Path(d)/'does-not-exist').rglob('*')),[])
    def test_evidence_precedes_new_guards(self):
        source=(M/'S11c_d_defect_packet_geometry.py').read_text()
        self.assertLess(source.index("J.emit('geometry-basis'"),source.index("require(kappa*kappa"))
        self.assertLess(source.index("J.emit(name+'-full-arrangement'"),source.index("require(plan['box']==box"))
        self.assertLess(source.index("require(plan['box']==box"),source.index('result=G.audit(plan,lines,cuts)'))


if __name__=='__main__':unittest.main()
