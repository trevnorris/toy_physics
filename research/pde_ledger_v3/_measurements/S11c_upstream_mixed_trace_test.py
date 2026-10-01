#!/usr/bin/env python3
"""Standard-library source/serialization tests only; no scientific execution."""
import ast
import json
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace

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

    def test_native_namespace_old_failure_and_corrected_execution(self):
        src=(HERE.parent/'scripts/S11c_c2_selfenergy_fold_sympy_audit.py').read_text()
        fragment=ns['assignment'](src,'kernel_bridge','z_three')
        left,right=object(),object()
        calls=[]
        class SourceMatrix:
            def __getitem__(self,index):return ('diagonal',index)
        class Transfer:
            def xreplace(self,mapping):
                calls.append(mapping)
                return ('transfer',mapping)
        common=dict(sp=SimpleNamespace(Matrix=lambda rows:rows),z_matrix=SourceMatrix(),
                    transfer=Transfer(),z_middle='middle')
        old=dict(common,leftmap=left,rightmap=right)
        with self.assertRaisesRegex(NameError,'left_map'):exec(fragment,old)
        corrected=dict(common,left_map=left,right_map=right)
        exec(fragment,corrected)
        self.assertEqual(calls,[left,right])
        self.assertEqual(corrected['z_three'][0][2],0)
        self.assertIs(corrected['z_three'][0][1][1],left)
        self.assertIs(corrected['z_three'][1][2][1],right)
        science=next(x for x in TREE.body if isinstance(x,ast.FunctionDef) and x.name=='science')
        update=next(x for x in ast.walk(science) if isinstance(x,ast.Call)
                    and isinstance(x.func,ast.Attribute) and x.func.attr=='update'
                    and any(k.arg=='z_matrix' for k in x.keywords))
        keys={k.arg for k in update.keywords}
        self.assertTrue({'left_map','right_map'}<=keys)
        self.assertFalse({'leftmap','rightmap'}&keys)
    def test_all_native_fragment_name_loads(self):
        src=(HERE.parent/'scripts/S11c_c2_selfenergy_fold_sympy_audit.py').read_text()
        cases=[('kernel_bridge','z_three',{'sp','z_matrix','transfer','left_map','right_map','z_middle'}),
               ('kernel_bridge','three_inverse',{'sp','coefficient','z_three'}),
               ('reference_pressure_kernels','trace_three',{'sp','value_coefficient',
                 'height_constant','normal_output','height_kernel','left','normal_middle','right','normal_input'})]
        for fn,target,keys in cases:
            tree=ast.parse(ns['assignment'](src,fn,target))
            loads={n.id for n in ast.walk(tree) if isinstance(n,ast.Name) and isinstance(n.ctx,ast.Load)}
            self.assertEqual(loads,keys)
    def test_final_evidence_precedes_related_guards(self):
        science=next(x for x in TREE.body if isinstance(x,ast.FunctionDef) and x.name=='science')
        emits={}
        guards={}
        for n in ast.walk(science):
            if not isinstance(n,ast.Call):continue
            if isinstance(n.func,ast.Attribute) and n.func.attr=='emit' and isinstance(n.args[0],ast.Constant):
                emits[n.args[0].value]=n.lineno
            if isinstance(n.func,ast.Name) and n.func.id=='require' and isinstance(n.args[-1],ast.Constant):
                guards[n.args[-1].value]=n.lineno
        for key in ['fixed external factors independent of middle momentum and shape grades',
                    'external closure domain','finite response factor','native outgoing middle closure domain']:
            self.assertLess(max(emits[n] for n in ['physical-factors','response-status','controls','middle-domain']),guards[key])

if __name__=='__main__':unittest.main()
