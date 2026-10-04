"""Stdlib constructor transport/refusal tests on manufactured expressions only."""
import ast
import unittest
import S11c_d_defect_packet_contracted_restore as R
import S11c_d_defect_packet_contraction_lib as OLD


class Factory:
    I=('constant','I');pi=('constant','pi')
    def __getattr__(self,name):
        if name=='Function':return lambda label:lambda *args:('function',label,args)
        return lambda *args,**kw:(name,args,kw)


class ConstructorTransport(unittest.TestCase):
    def restore(self,text):return R.restore_scalar(Factory(),{'srepr':text})
    def test_old_integer(self):self.assertEqual(self.restore('Integer(-17)'),('Integer',(-17,),{}))
    def test_old_rational(self):self.assertEqual(self.restore('Rational(-5,13)'),('Rational',(-5,13),{}))
    def test_old_assumptions(self):self.assertEqual(self.restore("Symbol('toy', real=True, finite=False)"),('Symbol',('toy',),{'real':True,'finite':False}))
    def test_original_subset_identical_constructor_calls(self):
        for text in ["Mul(Integer(7), Pow(Symbol('toy', positive=True), Rational(-1,3)))","Add(I, Rational(2,7), Integer(-9))"]:
            self.assertEqual(self.restore(text),OLD.restore_scalar(Factory(),{'srepr':text}))
    def test_sinh_and_pi(self):
        result=self.restore("sinh(Mul(pi,Symbol('toy',real=True)))")
        self.assertEqual(result,('sinh',(('Mul',(('constant','pi'),('Symbol',('toy',),{'real':True})),{}),),{}))
    def test_named_function_preserves_arguments(self):
        result=self.restore("Function('Hwhole')(Symbol('toy'),Integer(7),Integer(13))")
        self.assertEqual(result,('function','Hwhole',(('Symbol',('toy',),{}),('Integer',(7,),{}),('Integer',(13,),{}))))
    def test_all_named_arities(self):
        for name,n in R.FUNCTION_ARITIES.items():
            text="Function(%r)(%s)"%(name,','.join('Integer(%s)'%(j+2) for j in range(n)))
            self.assertEqual(len(self.restore(text)[2]),n)
    def test_old_reader_refuses_new_sinh(self):
        with self.assertRaises(ValueError):OLD.restore_scalar(Factory(),{'srepr':"sinh(Symbol('toy'))"})
    def test_old_reader_refuses_applied_function(self):
        with self.assertRaises(ValueError):OLD.restore_scalar(Factory(),{'srepr':"Function('common_outgoing_q')(Symbol('toy'))"})
    def test_plan_has_no_scientific_objects(self):
        import json
        plan=R.constructor_plan({'srepr':"Mul(pi,Function('common_outgoing_q')(Symbol('toy',real=True)))"})
        self.assertEqual(json.loads(json.dumps(plan)),plan)
    def test_unknown_function(self):
        with self.assertRaises(ValueError):self.restore("Function('unapproved')(Integer(1))")
    def test_wrong_function_arity(self):
        with self.assertRaises(ValueError):self.restore("Function('Hwhole')(Integer(1))")
    def test_function_keyword(self):
        with self.assertRaises(ValueError):self.restore("Function('Hwhole',real=True)(Integer(1),Integer(2),Integer(3))")
    def test_applied_keyword(self):
        with self.assertRaises(ValueError):self.restore("Function('common_outgoing_q')(Integer(1),foo=True)")
    def test_arbitrary_call_chain(self):
        with self.assertRaises(ValueError):self.restore("Symbol('x')(Integer(1))")
    def test_attribute_call(self):
        with self.assertRaises(ValueError):self.restore("module.sinh(Integer(1))")
    def test_attribute_constant(self):
        with self.assertRaises(ValueError):self.restore('module.pi')
    def test_unknown_constant(self):
        with self.assertRaises(ValueError):self.restore('E')
    def test_unknown_constructor(self):
        with self.assertRaises(ValueError):self.restore('sin(Integer(1))')
    def test_sinh_wrong_arity(self):
        with self.assertRaises(ValueError):self.restore('sinh(Integer(1),Integer(2))')
    def test_float(self):
        with self.assertRaises(ValueError):self.restore('Integer(1.2)')
    def test_float_constructor(self):
        with self.assertRaises(ValueError):self.restore("Float('1.2')")
    def test_boolean_integer(self):
        with self.assertRaises(ValueError):self.restore('Integer(True)')
    def test_zero_denominator(self):
        with self.assertRaises(ValueError):self.restore('Rational(1,0)')
    def test_wrong_rational_arity(self):
        with self.assertRaises(ValueError):self.restore('Rational(1)')
    def test_string_inside_add(self):
        with self.assertRaises(ValueError):self.restore("Add('pi',Integer(1))")
    def test_empty_add(self):
        with self.assertRaises(ValueError):self.restore('Add()')
    def test_wrong_power_arity(self):
        with self.assertRaises(ValueError):self.restore('Pow(Integer(2))')
    def test_unexpected_keyword(self):
        with self.assertRaises(ValueError):self.restore('Add(Integer(1),evaluate=False)')
    def test_assumption_unknown(self):
        with self.assertRaises(ValueError):self.restore("Symbol('x',extended_real=True)")
    def test_assumption_nonboolean(self):
        with self.assertRaises(ValueError):self.restore("Symbol('x',real=1)")
    def test_duplicate_assumption(self):
        with self.assertRaises(ValueError):self.restore("Symbol('x',real=True,real=False)")
    def test_symbol_arity(self):
        with self.assertRaises(ValueError):self.restore("Symbol('x','y')")
    def test_unpacking(self):
        with self.assertRaises(ValueError):self.restore("Add(*values)")
    def test_comprehension(self):
        with self.assertRaises(ValueError):self.restore('[Integer(x) for x in values]')
    def test_no_dynamic_evaluator(self):
        from pathlib import Path
        source=Path(R.__file__).read_text();tree=ast.parse(source)
        self.assertFalse(any(isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id in ('eval','exec','sympify','compile') for n in ast.walk(tree)))
    def test_no_scientific_import(self):
        from pathlib import Path
        tree=ast.parse(Path(R.__file__).read_text())
        self.assertEqual([n.names[0].name for n in tree.body if isinstance(n,ast.Import)],['ast'])


if __name__=='__main__':unittest.main()
