#!/usr/bin/env python3
"""Standard-library source/persistence tests only; no scientific import or CAS."""
import ast
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

HERE=Path(__file__).resolve().parent
WORKER=HERE/'S11c_upstream_mixed_tanh.py'
spec=importlib.util.spec_from_file_location('mixed_tanh_tooling',WORKER)
worker=importlib.util.module_from_spec(spec);spec.loader.exec_module(worker)


class Tooling(unittest.TestCase):
    def test_import_is_inert(self):
        self.assertFalse(hasattr(worker,'sp'))
        self.assertFalse(hasattr(worker,'Str'))

    def test_literal_only(self):
        source="x={'dtn_kernel':{'value':_restore(\"Tuple(Integer(1))\")}}"
        self.assertEqual(worker.literal_record(source,'dtn_kernel'),'Tuple(Integer(1))')
        with self.assertRaises(ValueError):
            worker.literal_record(source+'\ny='+source[2:],'dtn_kernel')
        with self.assertRaises(ValueError):
            worker.literal_record("x={'dtn_kernel':{'value':other('x')}}",'dtn_kernel')

    def test_function_extraction_does_not_import_producer(self):
        text="raise RuntimeError('producer must not execute')\ndef small(x):\n    return x+1\n"
        ns={};exec(worker.extract_function(text,'small'),ns)
        self.assertEqual(ns['small'](2),3)

    def test_immutable_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'x.json';worker.save(p,{'a':[1,2]})
            self.assertEqual(json.loads(p.read_text()),{'a':[1,2]})
            with self.assertRaises(FileExistsError):worker.save(p,{'a':0})

    def test_scientific_imports_are_after_containment(self):
        tree=ast.parse(WORKER.read_text())
        main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
        contained=min(n.lineno for n in ast.walk(main) if isinstance(n,ast.Call)
            and isinstance(n.func,ast.Name) and n.func.id=='containment')
        imported=min(n.lineno for n in ast.walk(main) if isinstance(n,ast.Import)
            and any(a.name=='sympy' for a in n.names))
        self.assertLess(contained,imported)
        self.assertFalse(any(isinstance(n,(ast.Import,ast.ImportFrom)) and
            ('sympy' in ast.unparse(n) or 'pickle' in ast.unparse(n)) for n in tree.body))


if __name__=='__main__':unittest.main()
