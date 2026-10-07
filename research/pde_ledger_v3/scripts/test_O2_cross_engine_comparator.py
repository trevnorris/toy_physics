#!/usr/bin/env python3
"""Synthetic serialized-stream controls. Never open production inputs/import engines."""
from dataclasses import replace
import io
import json
from pathlib import Path
import tempfile
import unittest
import sympy as sp

import O2_cross_engine_comparator as c


def sy(name):
    return f'Symbol({name!r})'


def fn(name,arg):
    return f'Function({name!r})({arg})'


def products(value):
    """Only computed outcomes/residuals/deltas; never echoed scalar operands."""
    if isinstance(value,list):
        return [products(v) for v in value]
    if not isinstance(value,dict):
        return value
    if 'outcome' in value:
        return {key:products(item) for key,item in value.items()
                if key in ('outcome','residual','structural_delta','children','operands',
                           'head','py_present','wl_present')}
    # Child keys and structural-delta fields are comparison products themselves.
    return {key:products(item) for key,item in value.items()
            if key not in ('py_value','wl_value','points','point_policy')}


class Controls(unittest.TestCase):
    def streams(self,py_text,wl_text,joins,*,names=c.NAME_TABLE):
        # Native serialization is part of every control's path, including repoints.
        with tempfile.TemporaryDirectory() as directory:
            p,w = Path(directory)/'py.out',Path(directory)/'wl.out'
            p.write_text(py_text)
            w.write_text(wl_text)
            output,accounting = io.StringIO(),io.StringIO()
            c.run(p,w,output,accounting,joins=joins,names=names)
        parsed = [json.loads(line) for line in output.getvalue().splitlines()]
        residual = [line['comparison'] for line in parsed if line['kind']=='residual']
        return residual,parsed,[json.loads(line) for line in accounting.getvalue().splitlines()]

    def fixture(self,py,wl,*,names=c.NAME_TABLE,py_tag='LEFT',wl_tag='RIGHT',key='item',layout=''):
        row = c.JoinRow('synthetic',(py_tag,'value'),(wl_tag,key),1,1,'synthetic',layout)
        return self.streams(f'PY_O2_{py_tag}: Tuple(Tuple(Str("value"), {py}))\n',
                            f'WL_O2_{wl_tag}: <|\n"{key}" -> {wl}\n|>\n',(row,),names=names)

    def changes(self,py,wl,mutation):
        before = self.fixture(py,wl)
        after = self.fixture(py,mutation)
        self.assertNotEqual(products(before[0]),products(after[0]),
                            'mutation must move outcome/residual/delta, not operand text')
        return before,after

    def census(self,result,side,field):
        return next(r[side][field] for r in result[1] if r['kind']=='structure')

    def changes_open(self,py,wl,mutation,field):
        before,after = self.changes(py,wl,mutation)
        self.assertNotEqual(self.census(before,'wl',field),self.census(after,'wl',field),
                            'OPEN control must also move its computed census')

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
                self.assertNotEqual(products(before),products(after),'repoint must move residual content')
                if row.kind in ('profile','coordinate','parameter'):
                    self.assertEqual(before[0]['outcome'],'exact')
                    self.assertEqual(after[0]['outcome'],'exact')
                    self.assertNotEqual(before[0]['residual'],after[0]['residual'])

    def test_applied_argument_stripped(self):
        self.changes(fn('V_r','Symbol("x1")'),'VR[x1]','VR')

    def test_live_profile_frozen(self):
        self.changes(fn('o2_rho_br_live','Symbol("x1")'),'RhoBr[x1]','11')

    def test_open_head(self):
        self.changes_open(fn('OPEN_Synthetic','Symbol("I_br_live"), Symbol("x1")'),
                     'OpenAction[Role,{IBrLive},x1]','OpenAction[OtherRole,{IBrLive},x1]',
                     'open_action_orientation_counts')

    def test_open_named_operand(self):
        self.changes_open(fn('OPEN_Synthetic','Symbol("I_br_live"), Symbol("x1")'),
                     'OpenAction[Role,{IBrLive},x1]','OpenAction[Role,{NBrLive},x1]',
                     'named_operands_and_heads')

    def test_open_argument(self):
        self.changes_open(fn('OPEN_Synthetic','Symbol("I_br_live"), Symbol("x1")'),
                     'OpenAction[Role,{IBrLive},VR[x1]]','OpenAction[Role,{IBrLive},VR[x2]]',
                     'named_operands_and_heads')

    def test_open_orientation(self):
        self.changes_open(fn('OPEN_Synthetic','Symbol("I_br_live"), Symbol("x1")'),
                     'OpenAction[Role,{IBrLive},x1]','-OpenAction[Role,{IBrLive},x1]',
                     'open_action_orientation_counts')

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
        self.assertNotEqual(before[2],after[2])
        self.assertEqual(after[2][0]['accounting'],'joined')
        self.assertFalse(any(r.get('accounting')=='unjoined' for r in after[2]))
        self.assertEqual(products(before[0]),products(after[0]))

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
        self.assertNotEqual(products(before),products(after))

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
        self.changes_open(fn('OPEN_Synthetic','Symbol("I_br_live")'),
                     'OpenAction[Role,{IBrLive},Nested[VR[x1]]]',
                     'OpenAction[Role,{IBrLive},Nested[VR[x2]]]', 'named_operands_and_heads')

    def test_exact_argument_content(self):
        before,after = self.changes(fn('V_r','Symbol("x1")'),'VR[x1]','VR[x2]')
        for result in (before,after):
            self.assertEqual(result[0][0]['outcome'],'exact')
        self.assertNotEqual(before[0][0]['residual'],after[0][0]['residual'])

    def test_common_translation_preserves_residual(self):
        base = self.fixture('Add(Pow(Symbol("x1"),Integer(2)),Integer(3))','x1*x2+5')
        translated = self.fixture('Add(Pow(Symbol("x1"),Integer(2)),Integer(3),Symbol("x3"))',
                                  'x1*x2+5+x3')
        self.assertEqual(products(base[0]),products(translated[0]))

    def test_residual_reconstructs_left_operand(self):
        # Independent synthetic construction: the residual must reconstruct the
        # left operand when added to the right. No residual value is prescribed.
        x1,x2 = sp.symbols('x1 x2',real=True)
        left,right = x1**2+2*x2+3,x1*x2-4
        result = self.fixture('Add(Pow(Symbol("x1"),Integer(2)),Mul(Integer(2),Symbol("x2")),Integer(3))',
                              'x1*x2-4')[0][0]
        self.assertEqual(result['outcome'],'exact')
        residual = sp.sympify(result['residual'],locals={'x1':x1,'x2':x2})
        self.assertEqual(sp.expand(right+residual),sp.expand(left))

    def test_transpose_layout_matches_explicit_transpose(self):
        py = 'ImmutableDenseMatrix([[Integer(2),Integer(3),Integer(5)],[Integer(7),Integer(11),Integer(13)]])'
        transpose = 'ImmutableDenseMatrix([[Integer(2),Integer(7)],[Integer(3),Integer(11)],[Integer(5),Integer(13)]])'
        wl = '{{2,7},{3,11},{5,13}}'
        mapped = self.fixture(py,wl,layout='transpose')
        explicit = self.fixture(transpose,wl)
        self.assertEqual(products(mapped[0]),products(explicit[0]))
        changed = self.fixture(py,'{{2,7},{3,11},{5,17}}',layout='transpose')
        self.assertNotEqual(products(mapped[0]),products(changed[0]))

    def test_held_aggregate_orientation(self):
        action = 'OpenAction[NativeToCoordinateDensity,{JMap},x1]'
        held = f'Inactive[Total][Inactive[Map][Function[{{s}},{action}],OpenAction[FaceSet,{{}},x2]]]'
        direct = self.fixture('Integer(2)',f'-({action})')
        aggregate = self.fixture('Integer(2)',f'-({held})')
        opposite = self.fixture('Integer(2)',held)
        field = 'open_action_orientation_counts'
        self.assertEqual(self.census(direct,'wl',field),self.census(aggregate,'wl',field))
        self.assertNotEqual(self.census(aggregate,'wl',field),self.census(opposite,'wl',field))
        # Moving the polarity inside the bound body must preserve its census.
        inside = f'Inactive[Total][Inactive[Map][Function[{{s}},-({action})],OpenAction[FaceSet,{{}},x2]]]'
        self.assertEqual(self.census(aggregate,'wl',field),
                         self.census(self.fixture('Integer(2)',inside),'wl',field))

    def test_sympy_aggregate_body_orientation(self):
        action = fn('OPEN_ApplyNativeFaceReduction','Symbol("J_map"), Symbol("x1")')
        held = fn('OPEN_SumOverAllNativeFaces',
                  'Symbol("J_map"), Lambda(Tuple(Symbol("o2_s_face")), '+action+')')
        negative = self.fixture('Mul(Integer(-1),'+held+')','2')
        positive = self.fixture(held,'2')
        field = 'open_action_orientation_counts'
        def body(result):
            return [r for r in self.census(result,'py',field) if r[0]=='py::OPEN_ApplyNativeFaceReduction']
        self.assertTrue(body(negative))
        self.assertNotEqual(body(negative),body(positive))
        direct = self.fixture('Mul(Integer(-1),'+action+')','2')
        self.assertEqual(body(negative),body(direct))

    def test_open_arguments_are_not_signed_balance_terms(self):
        outer = 'OpenAction[Outer,{OpenAction[Dependency,{},x1]},x2]'
        positive = self.fixture('Integer(2)',outer)
        negative = self.fixture('Integer(2)','-('+outer+')')
        term_field,arg_field = 'open_action_orientation_counts','open_argument_occurrence_counts'
        self.assertNotEqual(self.census(positive,'wl',term_field),self.census(negative,'wl',term_field))
        self.assertEqual(self.census(positive,'wl',arg_field),self.census(negative,'wl',arg_field))
        term_roles = {r[0] for r in self.census(negative,'wl',term_field)}
        argument_roles = {r[0] for r in self.census(negative,'wl',arg_field)}
        self.assertTrue(argument_roles)
        self.assertFalse(term_roles & argument_roles)

    def test_mass_loss_join_keeps_rhs_separate(self):
        row = next(r for r in c.JOIN_TABLE if r.py==('MASS_INPUT','outward_loss'))
        # This is an object/role assertion, not an expected residual. It rejects
        # a different-role occurrence even if its current payload happens to be
        # identical. The source citations identify this outward-loss occurrence.
        self.assertEqual(row.wl,('EXCHANGE_MOMENTUM','MaterialIdentification','@1'))
        other_rows = tuple(r for r in c.JOIN_TABLE
                           if r.label in ('mass_equation','mass_divergence'))
        self.assertEqual({r.label for r in other_rows},{'mass_equation','mass_divergence'})
        joins = (row,*other_rows)
        py = ('PY_O2_MASS_INPUT: Tuple('
              'Tuple(Str("outward_loss"),Function("j_n")(Symbol("x1"))),'
              'Tuple(Str("divergence"),Add(Symbol("x1"),Integer(3))),'
              'Tuple(Str("supplied_equation"),Equality(Add(Symbol("x1"),Integer(3)),'
              'Mul(Integer(-1),Function("j_n")(Symbol("x1"))))))\n')
        def stream(loss,rhs):
            return ('WL_O2_EXCHANGE_MOMENTUM: <|"MaterialIdentification" -> SameExchangedMaterial['+loss+',PiN,JMap]|>\n'
                    'WL_O2_MASS_INPUT: <|"Divergence" -> x1+3,'
                    '"Equation" -> x1+3 == -Jn[x1],"RHS" -> '+rhs+'|>\n')
        base = self.streams(py,stream('Jn[x1]','-Jn[x1]'),joins)
        rhs_changed = self.streams(py,stream('Jn[x1]','-Jn[x2]'),joins)
        loss_changed = self.streams(py,stream('Jn[x2]','-Jn[x1]'),joins)
        for result in (base,rhs_changed,loss_changed):
            resolved = {r['row'] for r in result[2] if r.get('accounting')=='joined'}
            self.assertEqual(resolved,{r.label for r in joins})
        self.assertEqual(products(base[0]),products(rhs_changed[0]))
        self.assertNotEqual(products(base[0]),products(loss_changed[0]))
        rhs = next(r for r in base[2] if r.get('path')==['MASS_INPUT','RHS'])
        self.assertEqual(rhs['accounting'],'unjoined')
        self.assertIn('same-role',rhs['reason'])

    def test_open_operand_binding_consistency(self):
        for py,wl in (
            ('material_action_compatibility','UnresolvedStressInertiaNormalIdentifications'),
            ('energy_accounting_overlap','UnresolvedEnergyOccurrenceIdentifications'),
            ('S12_reaction_system','OPENReactionSystem'),
            ('face_support_partition','UnresolvedSupportPartition')):
            with self.subTest(operand=py):
                result = self.fixture(fn('OPEN_Synthetic',sy(py)),f'OpenAction[Role,{{{wl}}},x1]')
                for side in ('py','wl'):
                    names = self.census(result,side,'named_operands_and_heads')
                    self.assertIn(py,names)
                    self.assertNotIn(side+'::'+(py if side=='py' else wl),names)


if __name__ == '__main__':
    unittest.main(verbosity=2)
