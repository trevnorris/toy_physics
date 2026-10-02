#!/usr/bin/env python3
"""Stdlib-only implementation tests. No science import or payload restoration."""
import ast
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import re
import unittest

HERE=Path(__file__).resolve().parent
WORKER=HERE/'S11c_d_defect_source_composition.py'
NAMES={'require','triples','jet_spec','broad_row_census','quotient_recurrence','definitions'}
NS={'ast':ast,'hashlib':hashlib,'itertools':itertools,'re':re,
    'G':((0,0),(1,0),(0,1),(1,1))}
TREE=ast.parse(WORKER.read_text())
exec(compile(ast.Module(body=[n for n in TREE.body if isinstance(n,ast.FunctionDef) and n.name in NAMES],type_ignores=[]),'stdlib-worker-functions','exec'),NS)

class Tests(unittest.TestCase):
    def test_grade_counts(self):
        self.assertEqual([len(NS['triples'](g)) for g in NS['G']],[1,3,3,9])
    def test_direct_support(self):
        self.assertEqual([(a,c) for a,b,c in NS['triples']((1,1)) if b==(1,1)],[((0,0),(0,0))])
    def test_grade_order_noncommutes(self):
        self.assertIn(((1,0),(0,1),(0,0)),NS['triples']((1,1)))
        self.assertIn(((0,0),(0,1),(1,0)),NS['triples']((1,1)))
    def test_jet_parser(self):
        self.assertEqual(NS['jet_spec']('u_2_ttd1d3')['spatialOrders'],[1,0,1])
        self.assertEqual(NS['jet_spec']('u_2_ttd1d3')['timeOrder'],2)
        self.assertEqual(NS['jet_spec']('grad_theta_2')['channel'],'theta')
        self.assertIsNone(NS['jet_spec']('e_W_d4'))
    def test_broad_names_not_whitelist(self):
        rows=NS['broad_row_census']("Add(Symbol('other'), Mul(Symbol('d_w_unknown'), Symbol('delta_p_plus')))")
        self.assertEqual([h['name'] for h in rows[1]['hits']],['d_w_unknown','delta_p_plus'])
        self.assertEqual(rows[0]['hits'],[])
    def test_full_child_partition(self):
        rows=NS['broad_row_census']("Add(Symbol('d_w_delta_p_minus'), Symbol('u_1'), Symbol('delta_p_plus'))")
        self.assertEqual([r['childIndex'] for r in rows],[0,1,2])
        self.assertTrue(all(hashlib.sha256(r['constructorText'].encode()).hexdigest()==r['sha256'] for r in rows))
    def test_rational_quotient_recurrence(self):
        # (2 + 3 eta + 5 sigma)/(1 + eta): exact discarded eta^2 must survive.
        out=NS['quotient_recurrence']({(0,0):Fraction(2),(1,0):Fraction(3),(0,1):Fraction(5)},
            {(0,0):Fraction(1),(1,0):Fraction(1)},tuple(itertools.product(range(3),repeat=2)),lambda a:a)
        self.assertEqual([out[g] for g in NS['G']],[2,1,5,-5])
        self.assertEqual(out[2,0],-1)
        self.assertEqual(out[2,1],5)
    def test_singular_denominator_refuses(self):
        with self.assertRaises(ZeroDivisionError):
            NS['quotient_recurrence']({(0,0):Fraction(1)},{(0,0):Fraction(0)},[(0,0)],lambda a:a)
    def test_no_top_level_scientific_import(self):
        for node in TREE.body:
            if isinstance(node,(ast.Import,ast.ImportFrom)):
                names=[x.name for x in node.names] if isinstance(node,ast.Import) else [node.module]
                self.assertFalse(any(n.startswith(('sympy','numpy','scipy')) for n in names))
    def test_no_old_producer_calls(self):
        calls={n.func.id for n in ast.walk(TREE) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name)}
        self.assertFalse(calls & {'Inputs','build_face','build_case','kernel_bridge','reference_pressure_kernels','shape_coefficients'})
        self.assertFalse(any(isinstance(n,ast.Attribute) and n.attr in {'integrate','solve','upper_triangular_solve','doit','Poly'} for n in ast.walk(TREE)))
    def test_gate_needs_concrete_assessment(self):
        gate=next(n for n in TREE.body if isinstance(n,ast.FunctionDef) and n.name=='verify_gate')
        text=ast.get_source_segment(WORKER.read_text(),gate)
        for fragment in ('correctedMethodAssessed','independentBuildClearance','workerSha256','launcherSha256','methodSha256'):
            self.assertIn(fragment,text)
    def test_native_assignment_targets_exist(self):
        native=(HERE.parent/'scripts/S11c_c2_selfenergy_fold_sympy_audit.py').read_text()
        fn=next(n for n in ast.parse(native).body if isinstance(n,ast.FunctionDef) and n.name=='build_face')
        targets=[t.id for n in fn.body if isinstance(n,ast.Assign) for t in n.targets if isinstance(t,ast.Name)]
        for target in ('extension','jet_transfer','normal_jet','reference_pressure'):
            self.assertEqual(targets.count(target),1)
    def test_metadata_route_class_supported(self):
        native=(HERE.parent/'scripts/S11c_c2_selfenergy_fold_sympy_audit.py').read_text()
        entries=[n for n in ast.parse(native).body if isinstance(n,(ast.ClassDef,ast.FunctionDef)) and n.name=='Inputs']
        self.assertEqual(len(entries),1)
        self.assertIsInstance(entries[0],ast.ClassDef)

if __name__=='__main__':unittest.main()
