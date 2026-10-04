"""Stdlib manufactured data and AST checks only; no native restoration/science."""
import ast,copy,json,unittest
from pathlib import Path
from collections import Counter
from S11c_d_defect_packet_contracted_continue_resume import (
    require,UNITS,ZERO,canonical_slots,check_census,group_indices,guard_inventory,prefix_return,exact_constant)
M=Path(__file__).resolve().parent

def entries():
    return [{'addressId':i,'face':'plus' if i<10 else 'minus','grade':[1,1],'unit':UNITS.copy(),
             'component':'NATIVE_MIXED_ITERATION' if i%2==0 else 'INHERITED_DIRECT_WHOLE_OFF_DIAGONAL',
             'primitives':['J'] if i%2==0 else ['Dr','Dh','Dq']} for i in range(20)]

class Census(unittest.TestCase):
    def setUp(self):self.entries=entries();self.expected=canonical_slots(self.entries,27,122,'synthetic')
    def test_full(self):self.assertTrue(check_census(self.expected,self.expected));self.assertEqual(len(self.expected),40)
    def test_twenty_unique(self):
        x=entries();x[1]['addressId']=0
        with self.assertRaises(ValueError):canonical_slots(x,27,122,'x')
    def test_wrong_primitive_selection(self):
        x=entries();x[1]['primitives']=['Dr']
        with self.assertRaises(ValueError):canonical_slots(x,27,122,'x')
    def test_unit_refusal(self):
        x=entries();x[0]['unit']=['0','0','0']
        with self.assertRaises(ValueError):canonical_slots(x,27,122,'x')
    def test_grade_refusal(self):
        x=entries();x[0]['grade']=[0,0]
        with self.assertRaises(ValueError):canonical_slots(x,27,122,'x')
    def test_missing(self):
        with self.assertRaises(ValueError):check_census(self.expected[:-1],self.expected)
    def test_duplicate(self):
        with self.assertRaises(ValueError):check_census(self.expected+[self.expected[0]],self.expected)
    def test_duplicate_replaces_missing(self):
        x=copy.deepcopy(self.expected);x[-1]=x[0]
        with self.assertRaises(ValueError):check_census(x,self.expected)
    def test_wrong_face(self):
        x=copy.deepcopy(self.expected);x[0]['face']='minus'
        with self.assertRaises(ValueError):check_census(x,self.expected)
    def test_wrong_K(self):
        x=copy.deepcopy(self.expected);x[0]['K']=29
        with self.assertRaises(ValueError):check_census(x,self.expected)
    def test_wrong_T(self):
        x=copy.deepcopy(self.expected);x[0]['T']=124
        with self.assertRaises(ValueError):check_census(x,self.expected)
    def test_wrong_carrier(self):
        x=copy.deepcopy(self.expected);x[0]['carrier']='other'
        with self.assertRaises(ValueError):check_census(x,self.expected)
    def test_wrong_purpose(self):
        x=copy.deepcopy(self.expected);x[0]['purpose']='mutant'
        with self.assertRaises(ValueError):check_census(x,self.expected)
    def test_order_not_identity(self):self.assertTrue(check_census(self.expected[::-1],self.expected))
    def test_direct_three_times(self):
        g=group_indices(self.expected);self.assertEqual(len(g['component/D']),30);self.assertEqual(len(g['component/J']),10)
        self.assertEqual(Counter(self.expected[i]['addressId'] for i in g['component/D']),{i:3 for i in range(1,20,2)})
    def test_positive_groups_all_retained(self):
        g=group_indices(self.expected);self.assertEqual(g['full'],list(range(40)));self.assertEqual(len(g['face/plus']),20);self.assertEqual(len(g['grade/1,1']),40)
    def test_same_census_mutant_path(self):
        # Zero numerical values cannot conceal a metadata mismatch.
        x=copy.deepcopy(self.expected)
        for row in x:row['total']=0
        self.assertTrue(check_census(x,self.expected));x[-1]['addressId']=0
        with self.assertRaises(ValueError):check_census(x,self.expected)

class Prefix(unittest.TestCase):
    def zero(self):return {'input':{'left':{'srepr':'Symbol("x")'},'right':{'srepr':'Symbol("x")'},'context':{},'newDerivation':True},'raw':{'residual':ZERO},'decision':{'cancelled':ZERO},'return':{'residual':ZERO}}
    def test_zero(self):self.assertEqual(prefix_return('new-fixture',self.zero()),'inherited-new-zero')
    def test_nonzero(self):
        x=self.zero();x['return']['residual']={'text':'1','srepr':'Integer(1)'}
        with self.assertRaises(ValueError):prefix_return('new-fixture',x)
    def test_missing_input(self):
        x=self.zero();del x['input']
        with self.assertRaises(ValueError):prefix_return('new-fixture',x)
    def test_changed_decision(self):
        x=self.zero();x['decision']['cancelled']={'text':'unknown','srepr':'Symbol("unknown")'}
        with self.assertRaises(ValueError):prefix_return('new-fixture',x)
    def test_changed_arguments(self):
        x=self.zero();del x['input']['left']
        with self.assertRaises(ValueError):prefix_return('new-fixture',x)
    def test_enlarged(self):
        x={'input':{'K':29,'T':124},'decision':{'outer':ZERO,'middle':ZERO,'originalDomain':True,'sourceFunctionCalled':False},'return':{'outer':ZERO,'middle':ZERO}}
        self.assertEqual(prefix_return('new-enlarged-tail-synthetic',x),'enlarged-tail-observation');x['input']['K']=27
        with self.assertRaises(ValueError):prefix_return('new-enlarged-tail-synthetic',x)
    def test_source_call_refused(self):
        x={'input':{'K':29,'T':124},'decision':{'outer':ZERO,'middle':ZERO,'originalDomain':True,'sourceFunctionCalled':True},'return':{'outer':ZERO,'middle':ZERO}}
        with self.assertRaises(ValueError):prefix_return('new-enlarged-tail-synthetic',x)
    def test_strict_boolean(self):
        for v in [1,0,None,[],{},'True']:
            with self.assertRaises(ValueError):require(v,'synthetic')
    def test_guard_inventory(self):
        text='def prepare():\n    require(True,"ok")\n    if enabled:\n        require(True,"branch")\n    require(False,"failed")\n    require(True,"unreached")\n'
        rows=guard_inventory(text,5);self.assertEqual([r['line'] for r in rows],[2,4,5,6]);self.assertEqual(rows[2]['disposition'],'failed');self.assertEqual(rows[3]['disposition'],'unreached');self.assertFalse(any(r['originalPredicateFunctionReexecuted'] for r in rows))

class ASTAndContracts(unittest.TestCase):
    def tree(self,n):return ast.parse((M/(n+'.py')).read_text())
    def function(self,n,fn):return next(x for x in self.tree(n).body if isinstance(x,ast.FunctionDef) and x.name==fn)
    def test_unchanged_numerical_worker_tail(self):
        old=self.function('S11c_d_defect_packet_contracted_numeric','run');new=self.function('S11c_d_defect_packet_contracted_continue','run');self.assertEqual(ast.dump(old),ast.dump(new))
    def test_unchanged_geometry_tail(self):
        old=self.function('S11c_d_defect_packet_contracted_prepare','prepare');start=next(i for i,n in enumerate(old.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='plans' for t in n.targets));new=self.function('S11c_d_defect_packet_contracted_continue_tail','continue_geometry')
        self.assertEqual(ast.dump(ast.Module(body=old.body[start:],type_ignores=[])),ast.dump(ast.Module(body=new.body[2:],type_ignores=[])))
    def test_no_completed_function_calls(self):
        tree=self.tree('S11c_d_defect_packet_contracted_continue_prepare')
        forbidden={'weighted_tail','exponential_moment','whole_densities','response_selection','source_ast_joins','combine','params','original_selector'}
        calls=[n.func.id if isinstance(n.func,ast.Name) else n.func.attr if isinstance(n.func,ast.Attribute) else '' for n in ast.walk(tree) if isinstance(n,ast.Call)]
        self.assertFalse(forbidden&set(calls))
    def test_no_science_import(self):
        for f in ['resume','prepare','tail']:
            tree=self.tree('S11c_d_defect_packet_contracted_continue_'+f)
            for n in ast.walk(tree):
                if isinstance(n,ast.Import):self.assertFalse(any(x.name.split('.')[0] in ('sympy','mpmath','numpy') for x in n.names))
                if isinstance(n,ast.ImportFrom):self.assertNotIn(n.module,('sympy','mpmath','numpy'))
    def test_H_formula_source_AST(self):
        orig=self.function('S11c_d_defect_packet_preflight','run') if any(isinstance(n,ast.FunctionDef) and n.name=='run' for n in self.tree('S11c_d_defect_packet_preflight').body) else self.tree('S11c_d_defect_packet_preflight')
        terms=[n.value for n in ast.walk(orig) if isinstance(n,ast.AugAssign) and isinstance(n.target,ast.Name) and n.target.id=='middle']
        target=ast.parse('(4/b)*cx*cy*E30*F[2]**2*sp.Rational(55,3)*sp.Rational(1,2)**Tlim',mode='eval').body
        self.assertEqual(len(terms),1);self.assertEqual(ast.dump(terms[0]),ast.dump(target))
    def test_exact_comparison_wrapped_for_strict_predicate(self):
        tree=self.tree('S11c_d_defect_packet_contracted_continue_prepare')
        for n in ast.walk(tree):
            if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='require':
                operand=n.args[0]
                scientific=any(isinstance(v,ast.Call) and isinstance(v.func,ast.Name) and v.func.id=='const' for v in ast.walk(operand))
                if scientific:self.assertTrue(isinstance(operand,ast.Call) and isinstance(operand.func,ast.Name) and operand.func.id=='bool')
    def test_saved_no_unit_conversion(self):
        tree=self.tree('S11c_d_defect_packet_contracted_continue_prepare')
        # Exact unit strings join originals; no new unit-inference callable.
        self.assertFalse(any(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr in ('dimension','unit','scale') for n in ast.walk(tree)))
    def test_no_gate_created(self):self.assertFalse((M/'S11c_d_defect_packet_contracted_continue_gate.json').exists())
    def test_all_modules_parse(self):
        for suffix in ['', '_prepare','_resume','_tail','_launch']:self.tree('S11c_d_defect_packet_contracted_continue'+suffix)

if __name__=='__main__':unittest.main()
