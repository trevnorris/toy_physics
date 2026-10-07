#!/usr/bin/env python3
"""Synthetic serialized-stream controls. Never open production inputs/import engines."""
from dataclasses import replace
import io
import json
from pathlib import Path
import tempfile
import unittest

import O2_cross_engine_comparator as c


def sy(name):
    return f'Symbol({name!r})'


def fn(name,arg):
    return f'Function({name!r})({arg})'


class Controls(unittest.TestCase):
    def fixture(self,py,wl,*,names=c.NAME_TABLE,py_tag='LEFT',wl_tag='RIGHT',key='item'):
        # Native serialization is part of every control's path, including repoints.
        with tempfile.TemporaryDirectory() as directory:
            p,w = Path(directory)/'py.out',Path(directory)/'wl.out'
            p.write_text(f'PY_O2_{py_tag}: Tuple(Tuple(Str("value"), {py}))\n')
            w.write_text(f'WL_O2_{wl_tag}: <|\n"{key}" -> {wl}\n|>\n')
            row = c.JoinRow('synthetic',(py_tag,'value'),(wl_tag,key),1,1,'synthetic')
            output,accounting = io.StringIO(),io.StringIO()
            c.run(p,w,output,accounting,joins=(row,),names=names)
        parsed = [json.loads(line) for line in output.getvalue().splitlines()]
        residual = [line['comparison'] for line in parsed if line['kind']=='residual']
        return residual,parsed,[json.loads(line) for line in accounting.getvalue().splitlines()]

    def changes(self,py,wl,mutation):
        before = self.fixture(py,wl)[0]
        after = self.fixture(py,mutation)[0]
        self.assertNotEqual(before,after,'mutation must move the printed residual, not just operand text')

    def test_one_sided_operand(self):
        self.changes('Add(Integer(2), Symbol("x1", real=True))','2+x1','3+x1')

    def test_form_change(self):
        self.changes('Mul(Symbol("x1"), Symbol("x2"))','x1*x2','x1+x2^2')

    def test_each_name_binding_repoint(self):
        for row in c.NAME_TABLE:
            with self.subTest(name=row.py):
                if row.kind == 'profile':
                    py,wl = fn(row.py,'Integer(7)'),f'{row.wl}[7]'
                else:
                    py,wl = f'Add(Integer(7), {sy(row.py)})',f'7+{row.wl}'
                before = self.fixture(py,wl)[0]
                # Exchange two declared object bindings. The destination is an
                # existing, independently named object, not an absent sentinel.
                other = next((r for r in c.NAME_TABLE if r != row and r.kind == row.kind),
                             next(r for r in c.NAME_TABLE if r != row))
                rows = tuple(replace(r,wl=other.wl) if r == row else
                             replace(r,wl=row.wl) if r == other else r for r in c.NAME_TABLE)
                after = self.fixture(py,wl,names=rows)[0]
                self.assertNotEqual(before,after,'repoint must move residual content')
                if row.kind in ('profile','coordinate','parameter'):
                    self.assertEqual(before[0]['outcome'],'exact')
                    self.assertEqual(after[0]['outcome'],'exact')
                    self.assertNotEqual(before[0]['residual'],after[0]['residual'])

    def test_applied_argument_stripped(self):
        self.changes(fn('V_r','Symbol("x1")'),'VR[x1]','VR')

    def test_live_profile_frozen(self):
        self.changes(fn('o2_rho_br_live','Symbol("x1")'),'RhoBr[x1]','11')

    def test_open_head(self):
        self.changes(fn('OPEN_Synthetic','Symbol("I_br_live"), Symbol("x1")'),
                     'OpenAction[Role,{IBrLive},x1]','OpenAction[OtherRole,{IBrLive},x1]')

    def test_open_named_operand(self):
        self.changes(fn('OPEN_Synthetic','Symbol("I_br_live"), Symbol("x1")'),
                     'OpenAction[Role,{IBrLive},x1]','OpenAction[Role,{NBrLive},x1]')

    def test_open_argument(self):
        self.changes(fn('OPEN_Synthetic','Symbol("I_br_live"), Symbol("x1")'),
                     'OpenAction[Role,{IBrLive},VR[x1]]','OpenAction[Role,{IBrLive},VR[x2]]')

    def test_open_orientation(self):
        self.changes(fn('OPEN_Synthetic','Symbol("I_br_live"), Symbol("x1")'),
                     'OpenAction[Role,{IBrLive},x1]','-OpenAction[Role,{IBrLive},x1]')

    def test_derivative_order(self):
        py = 'Subs(Derivative(Function("V_r")(Dummy("z")), Tuple(Dummy("z"), Integer(1))), Tuple(Dummy("z")), Tuple(Symbol("x1")))'
        self.changes(py,'Derivative[1][VR][x1]','Derivative[2][VR][x1]')

    def test_binder_structure(self):
        self.changes('Lambda(Tuple(Symbol("o2_s_face")), Function("OPEN_Synthetic")(Symbol("o2_s_face")))',
                     'Function[{s},OpenAction[Role,s]]','Function[{t},OpenAction[Role,s]]')

    def test_nested_sibling_removed(self):
        self.changes('Tuple(Tuple(Str("a"), Integer(2)), Tuple(Str("b"), Integer(3)))',
                     '<|"a" -> 2,"b" -> 3|>','<|"a" -> 2|>')

    def test_moved_tag_and_key_stays_joined(self):
        before = self.fixture('Integer(13)','13')
        after = self.fixture('Integer(13)','13',wl_tag='RENAMED_CONTAINER',key='relocated')
        self.assertNotEqual(before[1],after[1])
        self.assertEqual(after[2][0]['accounting'],'joined')
        self.assertFalse(any(r.get('accounting')=='unjoined' for r in after[2]))
        self.assertEqual(before[0],after[0])

    def test_duplicate_join_rejected(self):
        r = c.J('a','L','R/k',1,1,'fixture')
        with self.assertRaises(c.InputError):
            c.validate_tables((r,r),())

    def test_parent_child_join_rejected(self):
        with self.assertRaises(c.InputError):
            c.validate_tables((c.J('a','L','R/a',1,1,'fixture'),
                               c.J('b','L/x','R/b',1,1,'fixture')),())

    def test_duplicate_name_rejected_both_sides(self):
        r = c.NAME_TABLE[0]
        for extra in (replace(r,wl='Other'),replace(r,py='Other')):
            with self.assertRaises(c.InputError):
                c.validate_tables((),(r,extra))

    def test_boolean_does_not_hide_algebraic_sibling(self):
        residual,_,accounting = self.fixture('Tuple(true, Add(Integer(3), Symbol("x1")))',
                                            '{True,3+x1}')
        parts = residual[0]['children']
        self.assertEqual(parts['0']['outcome'],'not_formed')
        self.assertEqual(parts['1']['outcome'],'exact')
        self.assertIn('residual',parts['1'])
        self.assertLess(accounting[0]['compared_leaves'][0],accounting[0]['parsed_leaves'][0])

    def test_lossless_function_metadata_survives(self):
        a = c.py_parse("Function('native', **{'real': True, 'finite': True})(Symbol('x1'))")
        self.assertIn(('real','True'),a.args[0].options)
        before = self.fixture("Function('native', **{'real': True})(Symbol('x1'))",'native[x1]')[0]
        after = self.fixture("Function('native', **{'real': False})(Symbol('x1'))",'native[x1]')[0]
        self.assertNotEqual(before,after)

    def test_empty_containers_are_accounted(self):
        _,_,rows = self.fixture('Tuple()','{}')
        self.assertGreater(rows[0]['parsed_leaves'][0],0)
        self.assertEqual(rows[0]['accounting'],'joined')

    def test_matching_text_and_names_not_zero_evidence(self):
        for py,wl in [('Str("same")','"same"'),('Symbol("x1")','x1')]:
            r = self.fixture(py,wl)[0][0]
            self.assertEqual(r['outcome'],'not_formed')
            self.assertNotIn('residual',r)

    def test_grammatical_unknown_head_is_coverage(self):
        r = self.fixture('Mystery(Integer(2))','Mystery[2]')[0][0]
        self.assertEqual(r['outcome'],'not_formed')

    def test_malformed_input_rejected(self):
        for source,parser in [('Function(',c.py_parse),('<|"a" -> x',c.wl_parse),
                              ('<|"a"->x,"a"->y|>',c.wl_parse),('__import__("os").system("true")',c.py_parse)]:
            with self.assertRaises(c.InputError):
                parser(source)

    def test_declared_tables_are_injective(self):
        c.validate_tables(c.JOIN_TABLE,c.NAME_TABLE)

    def test_independent_census_accounts_for_nested_sibling(self):
        _,_,rows = self.fixture('Tuple(Tuple(Str("a"), Integer(2)), Tuple(Str("b"), true))',
                               '<|"a"->2,"b"->True,"extra"->{7,9}|>')
        summaries = [r for r in rows if r.get('kind') == 'stream_object_count']
        self.assertEqual(len(summaries),2)
        for row in summaries:
            self.assertEqual(row['parsed_leaves'],row['accounted_leaves'])
            self.assertLess(row['compared_leaves'],row['parsed_leaves'])

    def test_profile_derivative_reaches_scalar_subtraction(self):
        py = 'Subs(Derivative(Function("V_r")(Dummy("z")), Tuple(Dummy("z"), Integer(1))), Tuple(Dummy("z")), Tuple(Symbol("x1")))'
        result = self.fixture(py,'Derivative[1][VR][x1]')[0][0]
        self.assertEqual(result['outcome'],'exact')
        self.assertIn('residual',result)

    def test_unpaired_open_subtree_mutation_moves_residual(self):
        self.changes(fn('OPEN_Synthetic','Symbol("I_br_live")'),
                     'OpenAction[Role,{IBrLive},Nested[VR[x1]]]',
                     'OpenAction[Role,{IBrLive},Nested[VR[x2]]]')


if __name__ == '__main__':
    unittest.main(verbosity=2)
