"""Stdlib-only bookkeeping and expression-walker tests; no scientific payloads."""
import ast
from dataclasses import dataclass
import json
from pathlib import Path
import tempfile
import unittest
import S11c_d_numerical_radiating_end_maps as worker

@dataclass(frozen=True)
class Node:
    kind: str
    args: tuple=()
    @property
    def is_Add(self):return self.kind=='add'
    @property
    def is_Mul(self):return self.kind=='mul'
    @property
    def is_Pow(self):return self.kind=='pow'
    @property
    def base(self):return self.args[0]
    @property
    def exp(self):return self.args[1]
    def has(self,item):return self==item or any(v.has(item) for v in self.args if isinstance(v,Node))
class Integer(int):
    is_Integer=True

class ToolingTests(unittest.TestCase):
    def test_carrier_walk(self):
        e=Node('epsilon');c=Node('opaque');square=Node('pow',(e,Integer(2)))
        cases=[(Node('mul',(e,e,c)),{2}),(square,{2}),(Node('add',(square,Node('mul',(square,c)))),{2}),(Node('add',(square,c)),{0,2}),(c,{0})]
        for value,expected in cases:self.assertEqual(worker.carrier_degrees(value,e,{}),expected)
        for value in (Node('function',(e,)),Node('pow',(e,Integer(-1)))):
            with self.assertRaises(ValueError):worker.carrier_degrees(value,e,{})
    def test_journal_failure_preserves_operands(self):
        with tempfile.TemporaryDirectory() as name:
            j=worker.Journal(Path(name));self.assertEqual(j.call('good',{'source':1},lambda:7),7)
            def fail():worker.require(False,'test failure',{'residual':3})
            with self.assertRaises(ValueError):j.call('bad',{'source':2},fail)
            self.assertEqual(j.count,1);self.assertEqual(j.active,'bad')
            evidence=json.loads((Path(name)/'failed-check-evidence.json').read_text())
            self.assertEqual(evidence['operation'],'bad');self.assertEqual(j.store.count(),4)
            j.store.get(evidence['evidence']);j.store.integrity_check();j.store.close()
    def test_no_producer_or_timer(self):
        tree=ast.parse(Path(worker.__file__).read_text());imports=[];calls=[]
        for node in ast.walk(tree):
            if isinstance(node,ast.Import):imports.extend(x.name for x in node.names)
            elif isinstance(node,ast.ImportFrom):imports.append(node.module)
            elif isinstance(node,ast.Call):calls.append(ast.unparse(node.func))
        self.assertEqual([n for n in imports if n.startswith('S11')],['S11c_d_numerical_radiating_blob_store'])
        self.assertFalse(any(n.endswith(('.alarm','.setitimer','.Poly','.diff','.integrate','.nroots')) for n in calls))
        self.assertNotIn('signal',imports)

if __name__=='__main__':unittest.main()
