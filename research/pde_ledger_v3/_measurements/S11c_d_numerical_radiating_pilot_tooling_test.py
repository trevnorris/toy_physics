"""Local stdlib-only launch/storage tests. No scientific imports or payloads."""
import ast
import json
from pathlib import Path
import tempfile
import unittest
import S11c_d_numerical_radiating_end_maps_v2 as storage
import S11c_d_numerical_radiating_integration as integration

M=Path(__file__).resolve().parent

class ToolingTests(unittest.TestCase):
    def test_failure_saves_exact_inputs_and_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            journal=storage.Journal(Path(directory))
            def fail():integration.require(False,'synthetic failed check',{'computed':[1,2]})
            with self.assertRaises(ValueError):journal.call('integration/check',{'operand':[3,4]},fail)
            active=json.loads((Path(directory)/'active-operation.json').read_text())
            evidence=json.loads((Path(directory)/'failed-check-evidence.json').read_text())
            self.assertEqual(storage.decode(journal.store.get(active['input'])),{'operand':[3,4]})
            self.assertEqual(storage.decode(journal.store.get(evidence['evidence'])),{'computed':[1,2]})
            self.assertEqual(journal.count,0);journal.store.close()
    def test_immutable_saved_result_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            journal=storage.Journal(Path(directory));receipt=journal.blob('saved',{'a':3})
            with self.assertRaises(Exception):journal.blob('saved',{'a':4})
            self.assertEqual(storage.decode(journal.store.get(receipt)),{'a':3});journal.store.integrity_check();journal.store.close()
    def test_source_only_imports_and_no_deadlines(self):
        names=('pilot','integration','finite')
        for name in names:
            tree=ast.parse((M/f'S11c_d_numerical_radiating_{name}.py').read_text())
            for node in tree.body:
                if isinstance(node,ast.Import):self.assertFalse(any(x.name.startswith(('numpy','sympy','scipy')) for x in node.names))
                if isinstance(node,ast.ImportFrom):self.assertFalse((node.module or '').startswith(('numpy','sympy','scipy')))
            calls=[ast.unparse(n.func) for n in ast.walk(tree) if isinstance(n,ast.Call)]
            self.assertFalse(any(c.endswith(('.alarm','.setitimer','.nroots','.rsolve')) for c in calls))
    def test_selected_helper_definitions_and_dependencies_exist(self):
        engine=M.parent/'scripts/S11c_d_mixing_scattering_sympy_audit.py'
        names={getattr(n,'name',None) for n in ast.parse(engine.read_text()).body}
        self.assertTrue({'memo_xreplace','dag_free_symbols','BoundedActionQuadrature','BoundedSourceFourierQuadrature'}<=names)
        finite={getattr(n,'name',None) for n in ast.parse((M/'S11c_d_finite_scattering.py').read_text()).body}
        self.assertTrue({'polynomial_basis','BasisMomentum'}<=finite)
    def test_fixed_contrasts_and_signed_deficit_not_clipped(self):
        text=(M/'S11c_d_numerical_radiating_finite.py').read_text();tree=ast.parse(text)
        functions={getattr(n,'name',None):n for n in tree.body if isinstance(n,ast.FunctionDef)}
        calls=[ast.unparse(n.func) for n in ast.walk(functions['solve']) if isinstance(n,ast.Call)]
        self.assertNotIn('np.clip',calls)
        self.assertIn('1-output.sum(axis=0)/incident',text)
        self.assertIn('NEGATIVE_DEFICIT_UNRESOLVED',text)
        self.assertIn("('base','refinement'),('refinement','domain'),('domain','regulator')",text)

if __name__=='__main__':unittest.main()
