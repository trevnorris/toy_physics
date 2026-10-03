"""Stdlib regression tests only. Never import/restore scientific payloads."""
import ast
import json
from pathlib import Path
import runpy
import tempfile
import types
import unittest
from unittest.mock import patch

M=Path(__file__).resolve().parent
W=M/'S11c_d_defect_end_uniform_continue.py'
OLD=M/'S11c_d_defect_end_uniform.py'
ns=runpy.run_path(str(W),run_name='inert_tooling_tests')

class Sentinel:pass
class Test(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.out=Path(self.tmp.name)
    def refs(self,a,b):
        blobs={'a':{'value':a},'b':{'value':b}}
        receipts={k:{'path':'prior/'+k,'sha256':k*64,'bytes':10} for k in blobs}
        ref=ns['OperandReferences'](blobs,receipts)
        for k in blobs:ref.register(k,('value',))
        return ref
    def evidence(self,refs,predicate=None):
        class Base:
            def __init__(self,out):self.out=out;self.active=None
            def emit(self,name,value):
                p=self.out/(name+'.json');ns['save'](p,value)
                return {'path':p.name,'sha256':ns['sha'](p),'bytes':p.stat().st_size}
        if predicate is None:
            base={'np':types.SimpleNamespace(ndarray=Sentinel),'sp':types.SimpleNamespace(S=types.SimpleNamespace(true=object()))}
            exec(compile(ns['definitions'](OLD.read_text(),('exact_structure',)),'original-predicate','exec'),base)
            predicate=base['exact_structure']
        return ns['make_evidence'](Base,predicate,refs)(self.out)
    def test_old_failure_location_is_materialized_json(self):
        base={'Path':Path,'json':json,'os':ns['os']}
        exec(compile(ns['definitions'](OLD.read_text(),('save',)),'old-writer','exec'),base)
        with patch.object(json,'dumps',side_effect=MemoryError('synthetic allocation refusal')):
            with self.assertRaises(MemoryError):base['save'](self.out/'old.json',{'value':[1,2]})
            ns['save'](self.out/'new.json',{'value':[1,2]})
        self.assertEqual(json.loads((self.out/'new.json').read_text()),{'value':[1,2]})
    def test_writer_exclusive_and_rejects_nonfinite(self):
        ns['save'](self.out/'a.json',{'x':1})
        with self.assertRaises(FileExistsError):ns['save'](self.out/'a.json',{'x':2})
        with self.assertRaises(ValueError):ns['save'](self.out/'bad.json',{'x':float('nan')})
        self.assertTrue((self.out/'bad.json').exists())
        self.assertEqual(json.loads((self.out/'a.json').read_text()),{'x':1})
    def test_large_operand_is_not_rendered(self):
        class Unprintable:
            def __str__(self):raise AssertionError('must not render native source')
            def __repr__(self):raise AssertionError('must not render native source')
        a,b=Unprintable(),Unprintable();j=self.evidence(self.refs(a,b),lambda x,y:True)
        j.join('join',a,b)
        saved=json.loads((self.out/'join-operands.json').read_text())
        self.assertEqual(saved['actual']['selector'],['value']);self.assertEqual(saved['expected']['member'],'b')
    def test_unknown_operand_refused(self):
        a,b={'n':1},{'n':1};refs=self.refs(a,b)
        with self.assertRaises(ValueError):refs.reference({'n':1})
    def test_replaced_selector_refused(self):
        a,b={'n':1},{'n':1};refs=self.refs(a,b);refs.blobs['a']['value']={'n':1}
        with self.assertRaises(ValueError):refs.reference(a)
    def test_same_value_different_type_still_refused(self):
        a,b=[1],(1,);j=self.evidence(self.refs(a,b))
        with self.assertRaises(ValueError):j.join('bad',a,b)
        self.assertFalse(json.loads((self.out/'bad.json').read_text())['passed'])
        self.assertTrue((self.out/'bad-operands.json').exists())
    def test_exact_reuse_does_not_call_comparison_twice(self):
        calls=[];a,b={'x':[1,2]},{'x':[1,2]}
        j=self.evidence(self.refs(a,b),lambda x,y:calls.append((x,y)) or True)
        j.join('first',a,b);j.join('second',a,b)
        self.assertEqual(len(calls),1)
        second=json.loads((self.out/'second.json').read_text())
        self.assertFalse(second['predicateCalled']);self.assertEqual(second['reusedReturn']['path'],'first.json')
    def test_failure_keeps_operands_before_predicate(self):
        a,b=[1],[1]
        def fail(x,y):raise MemoryError('synthetic comparison allocation')
        j=self.evidence(self.refs(a,b),fail)
        with self.assertRaises(MemoryError):j.join('pending',a,b)
        self.assertTrue((self.out/'pending-operands.json').exists());self.assertFalse((self.out/'pending.json').exists())
        self.assertEqual(j.active,'pending')
    def test_actual_input_uses_full_blob_selector(self):
        a,b={'nested':[1,2]},{'nested':[1,2]};j=self.evidence(self.refs(a,b))
        j.emit('RIGHT-limit-1-actual-input',a)
        obj=json.loads((self.out/'RIGHT-limit-1-actual-input.json').read_text())
        self.assertEqual(obj['operand']['selector'],['value']);self.assertFalse(obj['printedScientificSummary'])
    def test_tail_starts_exactly_at_unsaved_join(self):
        text=OLD.read_text();parts=ns['source_parts'](text)
        self.assertIn("native_args['pairing']",ast.unparse(parts['grazingPending'][0]))
        self.assertEqual(ast.dump(ast.Module(body=parts['grazingFull'][5:],type_ignores=[])),ast.dump(ast.Module(body=parts['grazingPending'],type_ignores=[])))
        self.assertNotIn('wave_test(',ast.unparse(ast.Module(body=parts['controls'],type_ignores=[])))
        self.assertIn("record_control('opposite-sheet'",ast.unparse(ast.Module(body=parts['controls'],type_ignores=[])))
    def test_changed_tail_detection(self):
        old=ns['source_parts'](OLD.read_text());changed=ns['source_parts'](OLD.read_text().replace("p:sp.S.One,q:sp.Integer(2)","p:sp.Integer(7),q:sp.Integer(2)"))
        self.assertNotEqual(ns['ast_sha'](old['controls']),ns['ast_sha'](changed['controls']))
    def test_no_scientific_import_before_containment(self):
        tree=ast.parse(W.read_text());topimports=[n for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom))]
        self.assertNotIn('sympy',' '.join(ast.unparse(n) for n in topimports));self.assertNotIn('numpy',' '.join(ast.unparse(n) for n in topimports))
        body=W.read_text();self.assertLess(body.index("ns['containment']()"),body.index('import sympy as sp'))
        calls=[n for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name)]
        self.assertNotIn('run_science',[n.func.id for n in calls]);self.assertNotIn('wave_test',[n.func.id for n in calls])
    def test_exact_prior_snapshot_census(self):
        data=json.loads((M/'S11c_d_defect_end_uniform_files.json').read_text())
        self.assertEqual(len(data),6885);self.assertEqual(sum(r['bytes'] for r in data.values()),540133852)
        self.assertIn('LEFT-comparison.json',data);self.assertIn('RIGHT-limit--1-native.json',data)
        self.assertNotIn('RIGHT-native--1-pairing.json',data);self.assertNotIn('RIGHT-comparison.json',data)
    def test_original_worker_and_fragment_hashes_match_saved_failure(self):
        rec=json.loads((M/'S11c_d_defect_end_uniform_completion.json').read_text())
        self.assertEqual(rec['workerSha256'],ns['sha'](OLD))
        self.assertEqual(rec['literalBuildVerdicts']['claude'],'NEEDS REVISION')
        self.assertFalse(rec['scientificOverallAcceptance'])

if __name__=='__main__':unittest.main()
