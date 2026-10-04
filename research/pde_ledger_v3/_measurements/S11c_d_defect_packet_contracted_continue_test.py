"""Stdlib manufactured data and AST checks only; no native restoration/science."""
import ast,copy,json,tempfile,unittest
from pathlib import Path
from collections import Counter
from S11c_d_defect_packet_contracted_continue_resume import (
    require,UNITS,ZERO,canonical_slots,check_census,group_indices,guard_inventory,prefix_return,exact_constant,tree_inventory,check_tree,join_base_census,bound_source,independent_slots,tail_formula_join)
M=Path(__file__).resolve().parent

def entries():
    return [{'addressId':i,'face':'plus' if i<10 else 'minus','grade':[1,1],'unit':UNITS.copy(),
             'component':'NATIVE_MIXED_ITERATION' if i%2==0 else 'INHERITED_DIRECT_WHOLE_OFF_DIAGONAL',
             'primitives':['J'] if i%2==0 else ['Dr','Dh','Dq']} for i in range(20)]

class Census(unittest.TestCase):
    def setUp(self):
        self.entries=entries();self.expected=canonical_slots(self.entries,27,122,'synthetic')
        for row in self.expected:row.update(boundSource={'fixture':'original','K':27,'T':122},outerOperand={'text':'1','srepr':'Integer(1)'},middleOperand=ZERO)
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
        self.assertFalse(any(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and isinstance(n.func.value,ast.Name) and n.func.value.id in ('N','V','G') for n in ast.walk(tree)))
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

class IndependentAncestry(unittest.TestCase):
    def fixture(self):
        es=entries();addresses=[{'addressId':e['addressId'],'face':e['face'],'targetGrade':e['grade'],'component':e['component'],'status':'FORMAL'} for e in es]
        addresses += [{'addressId':i,'face':'plus','targetGrade':[0,0],'component':'OTHER','status':'EXACT_ZERO_SOURCE'} for i in range(20,544)]
        one={'text':'1','srepr':'Integer(1)'}
        plan={'K':27,'T':122,'allAddresses':[{'addressId':e['addressId'],'component':e['component'],'grade':e['grade'],'outer':one,'middle':ZERO} for e in es]}
        enlarged={e['addressId']:{'input':{'K':29,'T':124,'address':e},'return':{'outer':ZERO,'middle':ZERO}} for e in es}
        return es,addresses,plan,enlarged
    def test_native_base_census(self):
        es,a,p,x=self.fixture();self.assertEqual(len(join_base_census(p,a,es)),20)
    def test_unlisted_extra_live(self):
        es,a,p,x=self.fixture();a[-1]['status']='FORMAL';a[-1]['component']=es[0]['component']
        with self.assertRaises(ValueError):join_base_census(p,a,es)
    def test_extra_base_entry(self):
        es,a,p,x=self.fixture();p['allAddresses'].append(copy.deepcopy(p['allAddresses'][0]))
        with self.assertRaises(ValueError):join_base_census(p,a,es)
    def test_unbalanced_faces(self):
        es,a,p,x=self.fixture();a[0]['face']='minus'
        with self.assertRaises(ValueError):join_base_census(p,a,es)
    def test_wrong_native_grade(self):
        es,a,p,x=self.fixture();a[0]['targetGrade']=[0,0]
        with self.assertRaises(ValueError):join_base_census(p,a,es)
    def test_independent_windows(self):
        es,a,p,x=self.fixture();small=independent_slots(p,x,a,27,122,'x');large=independent_slots(p,x,a,29,124,'x')
        self.assertNotEqual(small[0]['outerOperand'],large[0]['outerOperand']);self.assertNotEqual(small[0]['boundSource'],large[0]['boundSource'])
    def test_wrong_saved_base_window(self):
        es,a,p,x=self.fixture();p['T']=124
        with self.assertRaises(ValueError):independent_slots(p,x,a,27,122,'x')
    def test_wrong_saved_enlarged_window(self):
        es,a,p,x=self.fixture();x[0]['input']['K']=27
        with self.assertRaises(ValueError):independent_slots(p,x,a,29,124,'x')
    def test_wrong_saved_enlarged_id(self):
        es,a,p,x=self.fixture();x[0]['input']['address']['addressId']=100
        with self.assertRaises(ValueError):independent_slots(p,x,a,29,124,'x')
    def test_swapped_value_refused(self):
        es,a,p,x=self.fixture();expected=independent_slots(p,x,a,27,122,'x');mut=copy.deepcopy(expected);mut[0]['outerOperand']=ZERO
        with self.assertRaises(ValueError):check_census(mut,expected)
    def test_wrong_bound_source_refused(self):
        es,a,p,x=self.fixture();expected=independent_slots(p,x,a,27,122,'x');mut=copy.deepcopy(expected);mut[0]['boundSource']['addressId']=19
        with self.assertRaises(ValueError):check_census(mut,expected)
    def test_prior_tree(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'a').write_bytes(b'abc');expected={'a':{'bytes':3}}
            self.assertTrue(check_tree(tree_inventory(p),expected,1,3));(p/'extra').write_bytes(b'x')
            with self.assertRaises(ValueError):check_tree(tree_inventory(p),expected,1,3)
    def test_prior_missing(self):
        with self.assertRaises(ValueError):check_tree({}, {'a':{'bytes':3}},1,3)
    def test_prior_wrong_bytes(self):
        with self.assertRaises(ValueError):check_tree({'a':{'bytes':2,'symlink':False}},{'a':{'bytes':3}},1,3)
    def test_prior_symlink(self):
        with self.assertRaises(ValueError):check_tree({'a':{'bytes':3,'symlink':True}},{'a':{'bytes':3}},1,3)
    def source_pair(self):
        old=(M/'S11c_d_defect_packet_contracted_prepare.py').read_text();base=(M/'S11c_d_defect_packet_preflight.py').read_text()
        node=next(v for v in ast.walk(ast.parse(base)) if isinstance(v,ast.FunctionDef) and v.name=='contributions')
        return old,ast.get_source_segment(base,node)
    def test_exact_actual_source_AST_pairs(self):
        old,base=self.source_pair();pairs=tail_formula_join(old,base);self.assertEqual(set(pairs),{'outer','J','H-overcount','D'});self.assertTrue(all(v['same'] for v in pairs.values()))
    def test_changed_H_source_refused(self):
        old,base=self.source_pair();pairs=tail_formula_join(old.replace('middlebound+=(4/b)','middlebound+=(5/b)'),base);self.assertFalse(pairs['H-overcount']['same'])
    def test_changed_D_source_refused(self):
        old,base=self.source_pair();pairs=tail_formula_join(old.replace('(36*121/b**2)','(37*121/b**2)'),base);self.assertFalse(pairs['D']['same'])
    def test_changed_outer_source_refused(self):
        old,base=self.source_pair();pairs=tail_formula_join(old.replace('outerbound=2*Cordinary','outerbound=3*Cordinary'),base);self.assertFalse(pairs['outer']['same'])
    def test_all_groups_before_budget_guard(self):
        text=(M/'S11c_d_defect_packet_contracted_continue_prepare.py').read_text();self.assertLess(text.index("J.emit('new-all-analytic-tail-group-decisions'"),text.index("'actual separate analytic subset aggregate budget'"))
    def test_gate_requires_literal_pair(self):
        text=(M/'S11c_d_defect_packet_contracted_continue.py').read_text();self.assertIn("g['literalBuildVerdicts']",text);self.assertIn("r['reports'][e]['literalVerdict']",text)

if __name__=='__main__':unittest.main()
