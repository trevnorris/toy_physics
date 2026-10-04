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
    def profile(self):return {'original':e("Symbol('w1_profile_d1')"),'map':[[e("Symbol('w1_profile_d1')"),e("Add(Rational(1,2),Mul(Rational(-1,2),Pow(tanh(Mul(Rational(1,10),Symbol('composition_x',real=True))),Integer(2))))")]],'oneDimensional':e('Integer(2)'),'L':e('Integer(10)'),'independentSigma':True,'profileTransformEvaluated':False}
    def context(self):return {'numeric':{'L_W':e('Integer(10)')}}
    def registry(self):return {'L_W':(1,0,0),'w1_profile_d1':(0,0,0),'w1_profile_d2':(0,0,0)}
    def test_profile_scale(self):self.assertEqual(U.profile_transport(self.profile(),self.context(),self.registry(),U.ScaleVerifier({},self.context(),Sink()))['entries'][0]['formalMappedUnit'],(0,0,0))
    def test_profile_missing_map(self):
        p=self.profile();p['map']=[]
        with self.assertRaises(ValueError):U.profile_transport(p,self.context(),self.registry(),U.ScaleVerifier({},self.context(),Sink()))
    def test_profile_wrong_L(self):
        p=self.profile();p['L']=e('Integer(1)')
        with self.assertRaises(ValueError):U.profile_transport(p,self.context(),self.registry(),U.ScaleVerifier({},self.context(),Sink()))
    def test_transverse_zero(self):
        p=self.profile();p['original']=e("Symbol('w1_profile_d2')");p['map']=[[p['original'],e('Integer(0)')]]
        self.assertTrue(U.profile_transport(p,self.context(),self.registry(),U.ScaleVerifier({},self.context(),Sink()))['entries'][0]['transverseZero'])
    def test_transverse_nonzero_refuses(self):
        p=self.profile();p['original']=e("Symbol('w1_profile_d2')");p['map'][0][0]=p['original']
        with self.assertRaises(ValueError):U.profile_transport(p,self.context(),self.registry(),U.ScaleVerifier({},self.context(),Sink()))
    def test_regular_flags_alone_refuse(self):
        d={'p-input.json':{'value':e('Integer(2)')},'p-state.json':{'finite':True,'nonzero':True}}
        with self.assertRaises(ValueError):U.regularity(d,'p',e('Integer(2)'),Sink())
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
    def test_ring_distributive_assembly(self):
        self.assertEqual(U.ring(U.parse("Mul(Integer(2),Add(Symbol('x'),Integer(1)))")),U.ring(U.parse("Add(Mul(Integer(2),Symbol('x')),Integer(2))")))
    def test_ring_unequal_operands(self):self.assertNotEqual(U.ring(U.parse("Symbol('x')")),U.ring(U.parse("Mul(Integer(2),Symbol('x'))")))
    def test_ring_inverse_scalar_normalization(self):self.assertEqual(U.ring(U.parse("Pow(Add(Integer(2),Mul(Integer(2),Symbol('x'))),Integer(-1))")),U.ring(U.parse("Mul(Rational(1,2),Pow(Add(Integer(1),Symbol('x')),Integer(-1)))")))
    def test_ring_imaginary_constant_inverse(self):self.assertEqual(U.ring(U.parse("Pow(I,Integer(-1))")),U.ring(U.parse('Mul(Integer(-1),I)')))
    def test_ring_zero_denominator(self):
        with self.assertRaises(ValueError):U.ring(U.parse('Pow(Integer(0),Integer(-1))'))
    def test_ring_no_float(self):
        with self.assertRaises(ValueError):U.ring(U.parse('Float(0.5)'))
    def test_ring_fractional_symbolic_power_refuses(self):
        with self.assertRaises(ValueError):U.ring(U.parse("Pow(Symbol('x'),Rational(1,2))"))
    def test_new_binding_operand_join(self):
        ctx={'densityMap':{'live':e("Mul(Symbol('rho'),Symbol('profile'))")},'profileEqualities':{'profile':e("Add(Integer(1),Symbol('eta'))")},'numeric':{'rho':e('Integer(2)')}}
        node,steps=U.binding_operand("Symbol('live')",ctx);self.assertEqual(U.ring(node),U.ring(U.parse("Add(Integer(2),Mul(Integer(2),Symbol('eta')))")));self.assertEqual(len(steps),3)
    def test_binding_cycle_refuses(self):
        ctx={'densityMap':{},'profileEqualities':{'a':e("Symbol('b')"),'b':e("Symbol('a')")},'numeric':{}}
        with self.assertRaises(ValueError):U.binding_operand("Symbol('a')",ctx)
    def test_assembly_remainder_uses_full_table(self):
        gs=[e("Symbol('eta')"),e("Symbol('sig')")];split={'retained':{'(0, 0)':e('Integer(1)'),'(1, 0)':e('Integer(2)'),'(0, 1)':e('Integer(0)'),'(1, 1)':e('Integer(0)')},'numerator':e("Add(Integer(1),Mul(Integer(2),Symbol('eta')))"),'denominator':e('Integer(1)')}
        out=U.full_remainder_operands(split['numerator'],split,gs);self.assertEqual(U.ring(out['fullHigherRemainder']),{});self.assertEqual(U.ring(out['quotientRingNumeratorRemainder']),{})
        split['retained']['(1, 0)']=e('Integer(3)');self.assertNotEqual(U.ring(U.full_remainder_operands(split['numerator'],split,gs)['fullHigherRemainder']),{})
    def test_profile_recurrence_actual_operand(self):
        p=self.profile();cert=U.scale_certificate('w',1,p['map'][0][1],{},self.context());self.assertTrue(cert['matches'])
    def test_profile_missing_scale_actual_path(self):
        p=self.profile();original=p['map'][0][1];p['map'][0][1]=e('Mul(Rational(1,10),'+original['srepr']+')');sink=Sink()
        with self.assertRaisesRegex(ValueError,'actual mapped profile'):U.profile_transport(p,self.context(),self.registry(),U.ScaleVerifier({},self.context(),sink))
        self.assertTrue(any(n.endswith('-decision-operands') and not v['matches'] for n,v in sink.values))
    def test_profile_argument_mismatch(self):
        p=self.profile();p['map'][0][1]=e(p['map'][0][1]['srepr'].replace('Rational(1,10)','Rational(1,11)'))
        with self.assertRaisesRegex(ValueError,'x/L'):U.profile_transport(p,self.context(),self.registry(),U.ScaleVerifier({},self.context(),Sink()))
    def test_scale_certificate_reuse_requires_actual_operands(self):
        p=self.profile();sink=Sink();v=U.ScaleVerifier({},self.context(),sink);v.verify('w',1,p['map'][0][1]);count=len(sink.values);v.verify('w',1,p['map'][0][1]);self.assertEqual(len(sink.values),count)
    def test_scale_catalog_collision_refuses(self):
        p=self.profile();q=copy.deepcopy(p);q['map'][0][1]=e('Integer(3)')
        with self.assertRaises(ValueError):U.scale_catalog({'p':p,'q':q},self.context())
    def test_higher_profile_requires_saved_predecessor(self):
        with self.assertRaises(ValueError):U.scale_certificate('w',2,e('Integer(1)'),{},self.context())
    def fraction_fixture(self):
        f={'numerator':e('Integer(1)'),'denominator':e('Integer(1)'),'finite':[[True,True],[True,True]],'signedNonzero':[[True,False],[True,False]]}
        d={'p-input.json':{'value':e('Integer(1)')},'p-fraction.json':f}
        for n in ('fraction-reconstruction','numerator-components','denominator-components'):
            d['p-'+n+'-input.json']={'left':e('Integer(1)'),'right':e('Integer(1)')};d['p-'+n+'-return.json']={'cancelled':e('Integer(0)')}
        return d
    def test_regular_one_nonzero_component_suffices(self):self.assertIn('fraction',U.regularity(self.fraction_fixture(),'p',e('Integer(1)'),Sink()))
    def test_regular_both_zero_components_refuse(self):
        d=self.fraction_fixture();d['p-fraction.json']['signedNonzero'][1]=[False,False]
        with self.assertRaises(ValueError):U.regularity(d,'p',e('Integer(1)'),Sink())
    def test_worker_collection_guards_are_explicit_boolean(self):
        t=ast.parse(P.with_name('S11c_d_defect_packet_source_units.py').read_text());hits=[]
        for n in ast.walk(t):
            if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='require' and n.args and isinstance(n.args[0],ast.Call) and isinstance(n.args[0].func,ast.Name) and n.args[0].func.id=='bool':
                hits.append(n.args[0].args[0].id)
        self.assertTrue({'maps','matching','candidates'}<=set(hits))
if __name__=='__main__':unittest.main()
