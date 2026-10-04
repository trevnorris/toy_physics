"""Synthetic helper tests only; no native unit/source calculation."""
import ast,copy,importlib.util,json,unittest
from pathlib import Path
P=Path(__file__).with_name('S11c_d_defect_packet_source_units_lib.py')
s=importlib.util.spec_from_file_location('source_units_test_library',P);U=importlib.util.module_from_spec(s);s.loader.exec_module(U)
def e(s):return {'srepr':s,'text':'synthetic'}
class Sink:
    def __init__(self):self.values=[]
    def emit(self,n,v):self.values.append((n,v))
class Tests(unittest.TestCase):
    def test_same_ast_not_printed_text(self):self.assertTrue(U.equal(e('Integer(1)'),{'srepr':'Integer( 1 )','text':'different'}))
    def test_false_identity(self):self.assertFalse(U.equal(e('Integer(1)'),e('Integer(2)')))
    def test_no_constant_restoration(self):self.assertEqual(U.symbols("Add(Integer(3),I)"),set())
    def test_symbol_census(self):self.assertEqual(U.symbols("Add(Symbol('x',real=True), Symbol('y'))"),{'x','y'})
    def test_tuple_leaf(self):self.assertTrue(U.equal(U.encoded_leaf(e('Tuple(Integer(1),Integer(2))'),1),e('Integer(2)')))
    def test_matrix_unit(self):self.assertEqual(U.matrix_unit(e('MutableDenseMatrix([[Integer(-1)],[Integer(2)],[Integer(0)]])')),(-1,2,0))
    def test_wrong_matrix_shape(self):
        with self.assertRaises(ValueError):U.matrix_unit(e('MutableDenseMatrix([[Integer(1)]])'))
    def test_unknown_unit_refuses(self):
        with self.assertRaises(ValueError):U.Units({}).dimension(U.parse("Symbol('unknown')"))
    def test_bad_sum_refuses(self):
        with self.assertRaises(ValueError):U.Units({'x':(1,0,0)}).dimension(U.parse("Add(Symbol('x'),Integer(1))"))
    def test_dimensional_substitution(self):self.assertEqual(U.Units({'L':(1,0,0),'g':(0,0,0)}).dimension(U.parse("Mul(Symbol('L'),Add(Integer(1),Symbol('g')))")),(1,0,0))
    def test_zero_no_unit(self):self.assertIsNone(U.Units({}).dimension(U.parse('Integer(0)')))
    def test_reciprocal_actual_base(self):
        text="Pow(Symbol('x'),Integer(-1))";walk=U.Units({'x':(2,0,0)});walk.dimension(U.parse(text));r=U.reciprocal_bases(text,walk.events);self.assertEqual(r[0]['baseUnit'],['2','0','0'])
    def test_reciprocal_missing_path(self):
        text="Pow(Symbol('x'),Integer(-1))";walk=U.Units({'x':(2,0,0)});walk.dimension(U.parse(text))
        with self.assertRaises(ValueError):U.reciprocal_bases(text,walk.events[1:])
    def test_reciprocal_wrong_root(self):
        walk=U.Units({'x':(2,0,0)});walk.dimension(U.parse("Symbol('x')"))
        with self.assertRaises(ValueError):U.reciprocal_bases('Integer(1)',walk.events)
    def test_strip_epsilon(self):self.assertTrue(U.equal_text(U.strip_epsilon(e("Mul(Integer(2),Symbol('epsilon_shape'),Symbol('x'))"),e("Symbol('epsilon_shape')")),"Mul(Integer(2),Symbol('x'))"))
    def test_strip_epsilon_zero(self):self.assertEqual(U.strip_epsilon(e('Integer(0)'),e("Symbol('epsilon_shape')")),'Integer(0)')
    def test_strip_epsilon_missing(self):
        with self.assertRaises(ValueError):U.strip_epsilon(e("Symbol('x')"),e("Symbol('epsilon_shape')"))
    def test_strip_epsilon_twice(self):
        with self.assertRaises(ValueError):U.strip_epsilon(e("Mul(Symbol('epsilon_shape'),Symbol('epsilon_shape'))"),e("Symbol('epsilon_shape')"))
    def test_jet_unit(self):self.assertEqual(U.jet_unit({'name':'e_W_d1_t','channel':'e_W','timeOrder':1,'spatialOrders':[1,0,0]},{'e_W_d1_t':(-1,-1,0)},e("Symbol('e_W_d1_t')")),(-1,-1,0))
    def test_jet_schema_disagreement(self):
        with self.assertRaises(ValueError):U.jet_unit({'name':'e_W','channel':'e_W','timeOrder':0,'spatialOrders':[0,0,0]},{'e_W':(1,0,0)},e("Symbol('e_W')"))
    def profile(self):return {'original':e("Symbol('w1_profile_d1')"),'map':[[e("Symbol('w1_profile_d1')"),e('Integer(2)')]],'oneDimensional':e('Integer(2)'),'L':e('Integer(10)'),'independentSigma':True,'profileTransformEvaluated':False}
    def context(self):return {'numeric':{'L_W':e('Integer(10)')}}
    def registry(self):return {'L_W':(1,0,0),'w1_profile_d1':(0,0,0),'w1_profile_d2':(0,0,0)}
    def test_profile_scale(self):self.assertEqual(U.profile_transport(self.profile(),self.context(),self.registry())['entries'][0]['formalMappedUnit'],(0,0,0))
    def test_profile_missing_map(self):
        p=self.profile();p['map']=[]
        with self.assertRaises(ValueError):U.profile_transport(p,self.context(),self.registry())
    def test_profile_wrong_L(self):
        p=self.profile();p['L']=e('Integer(1)')
        with self.assertRaises(ValueError):U.profile_transport(p,self.context(),self.registry())
    def test_transverse_zero(self):
        p=self.profile();p['original']=e("Symbol('w1_profile_d2')");p['map']=[[p['original'],e('Integer(0)')]]
        self.assertTrue(U.profile_transport(p,self.context(),self.registry())['entries'][0]['transverseZero'])
    def test_transverse_nonzero_refuses(self):
        p=self.profile();p['original']=e("Symbol('w1_profile_d2')");p['map'][0][0]=p['original']
        with self.assertRaises(ValueError):U.profile_transport(p,self.context(),self.registry())
    def test_regular_flags(self):
        d={'p-input.json':{'value':e('Integer(2)')},'p-state.json':{'finite':True,'nonzero':True}};s=Sink();U.regularity(d,'p',e('Integer(2)'),s);self.assertEqual(len(s.values),1)
    def test_regular_unknown_refuses(self):
        d={'p-input.json':{'value':e('Integer(2)')},'p-state.json':{'finite':None,'nonzero':True}}
        with self.assertRaises(ValueError):U.regularity(d,'p',e('Integer(2)'),Sink())
    def test_regular_other_operand_refuses(self):
        d={'p-input.json':{'value':e('Integer(2)')},'p-state.json':{'finite':True,'nonzero':True}}
        with self.assertRaises(ValueError):U.regularity(d,'p',e('Integer(3)'),Sink())
    def test_field_proof(self):self.assertTrue(U.field_proof({'degree':0,'coefficients':[e('Integer(2)')],'denominator':e('Integer(1)')},{'left':e('Integer(2)'),'right':e('Integer(2)')},{'cancelled':e('Integer(0)')},e('Integer(2)'))['fieldReconstructionInherited'])
    def test_field_proof_wrong_input(self):
        with self.assertRaises(ValueError):U.field_proof({'degree':0,'coefficients':[e('Integer(2)')],'denominator':e('Integer(1)')},{'left':e('Integer(3)')},{'cancelled':e('Integer(0)')},e('Integer(2)'))
    def test_no_sympy_import(self):
        for name in ('S11c_d_defect_packet_source_units.py','S11c_d_defect_packet_source_units_lib.py'):
            t=ast.parse(P.with_name(name).read_text());self.assertFalse(any(isinstance(n,ast.Import) and any(a.name.startswith('sympy') for a in n.names) or isinstance(n,ast.ImportFrom) and str(n.module).startswith('sympy') for n in ast.walk(t)))
if __name__=='__main__':unittest.main()
