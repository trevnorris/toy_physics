"""Stdlib-only label/refusal and exact-input transport tests; no science."""
import ast, hashlib, json, runpy, unittest
from pathlib import Path
M=Path(__file__).resolve().parent
P=M/'S11c_d_background_force_pairing.py'
NS=runpy.run_path(str(P),run_name='metadata_test_only')
MAN=json.loads((M/'S11c_d_background_force_pairing_inputs.json').read_text())
class ToolingTests(unittest.TestCase):
    def test_exact_label_selection(self):
        self.assertEqual(NS['named']([('U',(0,0,0)),('E_W','unchanged')],'E_W'),'unchanged')
    def test_duplicate_label_refused(self):
        with self.assertRaises(ValueError):NS['named']([('E_W',1),('E_W',2)],'E_W')
    def test_missing_label_refused(self):
        with self.assertRaises(ValueError):NS['named']([('U',0)],'E_W')
    def test_strict_predicate(self):
        NS['require'](True,'ok')
        for v in [False,1,0,None,[],[1],'true']:
            with self.assertRaises(ValueError):NS['require'](v,'refuse')
    def test_every_input_receipt(self):
        for r in MAN['savedFiles'].values():
            p=Path(r['path']);self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),r['sha256']);self.assertEqual(p.stat().st_size,r['bytes'])
    def test_literal_original_join(self):
        for rec in MAN['registryInputs'].values():
            v=json.loads(Path(MAN['savedFiles'][rec['savedAlias']]['path']).read_text());p=Path(rec['originalSource']);lines=p.read_text().splitlines()
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),v['sourceSha256'])
            self.assertEqual(ast.literal_eval(lines[v['displaySourceLine']-1].strip().removeprefix("'display': ").removesuffix(',')),v['completePublishedDisplay'])
            self.assertEqual(ast.literal_eval(lines[v['constructorSourceLine']-1].strip().removeprefix("'value': _restore(").removesuffix('),')),v['completePublishedConstructorString'])
    def test_helper_selection_only(self):
        tree=NS['definitions'](Path(MAN['helperSource']).read_text());self.assertEqual({n.name for n in tree.body},set(NS['HELPERS']))
    def test_no_ready_gate(self):
        self.assertFalse((M/'S11c_d_background_force_pairing_gate.json').exists())
if __name__=='__main__':unittest.main()
