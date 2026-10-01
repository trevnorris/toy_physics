#!/usr/bin/env python3
"""Standard-library source/serialization tests only; no scientific execution."""
import ast
import json
from pathlib import Path
import tempfile
import unittest

HERE = Path(__file__).parent
WORKER = HERE/'S11c_upstream_mixed_trace.py'
SOURCE = WORKER.read_text()
TREE = ast.parse(SOURCE)
ns = {'ast':ast,'json':json,'Path':Path}
for name in ('require','extract_function','literal_record','assignment','selected_constructor'):
    node = next(x for x in TREE.body if isinstance(x,ast.FunctionDef) and x.name == name)
    exec(compile(ast.Module(body=[node],type_ignores=[]),str(WORKER),'exec'),ns)

class Tooling(unittest.TestCase):
    def test_native_text_roundtrip_all_selected(self):
        for record in json.loads((HERE/'S11c_upstream_mixed_trace_native.json').read_text()).values():
            actual=ns['selected_constructor'](Path(record['source']).read_text(),record['recordKey'],
                                               record['case'],record['outer'])
            self.assertEqual(actual,record['constructor'])
    def test_negative_case_label(self):
        raw="Tuple(Tuple(Tuple(Str('B'), Integer(-1)), Tuple(Tuple(Str('VALUE'), Integer(7)))))"
        src="x={'row':{'value':_restore("+repr(raw)+")}}"
        self.assertEqual(ns['selected_constructor'](src,'row',['B',-1]),'Integer(7)')
    def test_duplicate_case_refused(self):
        row="Tuple(Tuple(Str('B'),Integer(1)),Tuple(Tuple(Str('VALUE'),Integer(0))))"
        src="x={'row':{'value':_restore("+repr('Tuple('+row+','+row+')')+")}}"
        with self.assertRaises(ValueError):ns['selected_constructor'](src,'row',['B',1])
    def test_missing_case_refused(self):
        src="x={'row':{'value':_restore(\"Tuple(Tuple(Tuple(Str('B'),Integer(1)),Tuple(Tuple(Str('VALUE'),Integer(0)))))\")}}"
        with self.assertRaises(ValueError):ns['selected_constructor'](src,'row',['C',1])
    def test_native_assignment_extraction(self):
        src=(HERE.parent/'scripts/S11c_c2_selfenergy_fold_sympy_audit.py').read_text()
        for fn,name in [('kernel_bridge','z_three'),('kernel_bridge','three_inverse'),
                        ('reference_pressure_kernels','trace_three')]:
            fragment=ns['assignment'](src,fn,name)
            self.assertIsInstance(ast.parse(fragment).body[0],ast.Assign)
        fragment=ns['assignment'](src,'kernel_bridge','z_three')
        slot=ast.parse(fragment).body[0].value.args[0].elts[0].elts[2]
        self.assertEqual(ast.literal_eval(slot),0)
    def test_scientific_import_after_containment(self):
        top_imports=[x for x in TREE.body if isinstance(x,(ast.Import,ast.ImportFrom))]
        self.assertFalse(any('sympy' in ast.unparse(x) for x in top_imports))
        main=next(x for x in TREE.body if isinstance(x,ast.FunctionDef) and x.name=='main')
        text=ast.get_source_segment(SOURCE,main)
        self.assertLess(text.index('containment()'),text.index('import sympy'))
    def test_no_replay_entrypoints(self):
        # Worker embeds neither a producer import nor a call to prior science.
        calls=[x.func.id for x in ast.walk(TREE) if isinstance(x,ast.Call) and isinstance(x.func,ast.Name)]
        self.assertNotIn('native',calls)
        self.assertNotIn('boundary',calls)
        self.assertNotIn('action',calls)
        self.assertNotIn('load_model',calls)
    def test_no_deadline_calls(self):
        calls=[ast.unparse(x.func) for x in ast.walk(TREE) if isinstance(x,ast.Call)]
        self.assertFalse(any(x.endswith(('.alarm','.setitimer')) for x in calls))
        self.assertNotIn('timeout=',SOURCE)

if __name__=='__main__':unittest.main()
