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
NAMES={'require','triples','jet_spec','broad_row_census','quotient_recurrence','definitions','native_jet_dimensions','select_certified_candidate','epsilon_power','address_sum_for','complete_factor_map','join_source_input','verify_helper_paths'}
NS={'ast':ast,'hashlib':hashlib,'itertools':itertools,'re':re,
    'ROOT':Path('/var/projects/toy_physics'),
    'G':((0,0),(1,0),(0,1),(1,1))}
TREE=ast.parse(WORKER.read_text())
exec(compile(ast.Module(body=[n for n in TREE.body if isinstance(n,ast.FunctionDef) and n.name in NAMES],type_ignores=[]),'stdlib-worker-functions','exec'),NS)


class SourceFixture:
    """Synthetic scalar substitution model, not a restored scientific payload."""
    class Raw:
        def __init__(self,face,density=2):self.face,self.density=face,density
        def subs(self,mapping,simultaneous=False):
            assert simultaneous
            return ('identified',mapping['s11cc1_V_lab_held_'+self.face],
                    mapping['s11cc1_mu_theta_lab_held_'+self.face],self.density)
    def __init__(self,face):
        self.events=[];self.face=face
        self.inp={'raw':self.Raw(face),'velocityAmplitude':Fraction(5),'chemicalAmplitude':Fraction(7),
                  'velocityCoefficient':Fraction(2),'combined':Fraction(41,2)}
        self.chemical={'raw':('mu',Fraction(14)),'amplitude':Fraction(7)}
        self.normalization={'velocityCoefficient':Fraction(2)}
    def emit(self,name,value):self.events.append(('input',name,value))
    def zero(self,name,left,right):
        self.events.append(('zero',name,(left,right)))
        if left!=right:raise ValueError(name)
    def bind(self,value):
        if isinstance(value,tuple):
            _,v,m,d=value
            return 2*v+3*m/d
        return value
    def one(self,objects,name):
        if name not in ('s11cc1_V_lab_held_'+self.face,'s11cc1_mu_theta_lab_held_'+self.face,'rho_br_bg_rho4_constant'):
            raise ValueError('wrong native face atom')
        return name
    def run(self):
        NS['join_source_input'](self,self.face,self.inp,self.chemical,self.normalization,self.bind,self.one,Fraction(2))

class Tests(unittest.TestCase):
    def test_source_join_routes_both_faces_and_saves_before_checks(self):
        for face in ('plus','minus'):
            fixture=SourceFixture(face);fixture.run()
            self.assertEqual([v[0] for v in fixture.events],['input','input','zero','zero','zero','zero'])
            self.assertIn('stage2Map',fixture.events[0][2])
            self.assertEqual(fixture.events[1][2]['bound'],Fraction(41,2))
    def test_source_join_rejects_native_chemical_corruption(self):
        fixture=SourceFixture('plus');fixture.chemical['raw']=('mu',Fraction(16))
        with self.assertRaisesRegex(ValueError,'chemical-epsilon'):fixture.run()
        self.assertEqual(fixture.events[0][0],'input')
    def test_source_join_rejects_disconnected_chemical_amplitude(self):
        fixture=SourceFixture('minus');fixture.inp['chemicalAmplitude']=Fraction(9)
        with self.assertRaisesRegex(ValueError,'chemical-amplitude'):fixture.run()
    def test_source_join_rejects_missing_live_density_binding(self):
        fixture=SourceFixture('minus');fixture.inp['raw'].density=1
        with self.assertRaisesRegex(ValueError,'live-density'):fixture.run()
    def test_source_join_rejects_wrong_normalization_and_combined(self):
        for key in ('velocityCoefficient','combined'):
            fixture=SourceFixture('plus');fixture.inp[key]+=1
            with self.assertRaises(ValueError):fixture.run()
    def test_helper_paths_must_match_actual_command(self):
        paths={'sharedGuard':str(NS['ROOT']/'scripts/s11c_guarded_run.py'),
               'supervisor':str(HERE/'S11c_d_end_normalization_run.py')}
        NS['verify_helper_paths'](paths)
        for key in paths:
            with self.assertRaises(ValueError):NS['verify_helper_paths']({**paths,key:'/tmp/hash-matching-copy.py'})
            missing=dict(paths);missing.pop(key)
            with self.assertRaises(KeyError):NS['verify_helper_paths'](missing)
        gate=next(n for n in TREE.body if isinstance(n,ast.FunctionDef) and n.name=='verify_gate')
        self.assertTrue(any(isinstance(n,ast.Call) and getattr(n.func,'id',None)=='verify_helper_paths' for n in ast.walk(gate)))
    def test_saved_source_dependencies_are_consumed(self):
        calls={n.args[0].value for n in ast.walk(TREE) if isinstance(n,ast.Call)
               and getattr(n.func,'id',None)=='load' and n.args and isinstance(n.args[0],ast.Constant)}
        for name in ('native-chemical-amplitude','chemical-amplitude-domain','inherited-source-normalization','selected-fourier-contraction'):
            self.assertIn('consumer/'+name+'.json',calls)
        self.assertNotIn("'densityRestored':True",WORKER.read_text())
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
    def test_native_dimension_rule_is_read_from_source(self):
        native=(HERE.parent/'scripts/S11c_c2_selfenergy_fold_sympy_audit.py').read_text()
        fn=next(n for n in ast.parse(native).body if isinstance(n,ast.FunctionDef) and n.name=='wave_jet')
        text=ast.get_source_segment(native,fn)
        result=NS['native_jet_dimensions'](text)
        self.assertEqual(result['velocity'],[1,0,0])
        self.assertEqual(result['scalar'],[0,0,0])
        self.assertEqual(result['sourceSha256'],hashlib.sha256(text.encode()).hexdigest())
        with self.assertRaises(ValueError):NS['native_jet_dimensions'](text.replace('(1, 0, 0) if','(2, 0, 0) if'))
    def test_candidate_search_skips_zero_and_unknown(self):
        entries=[{'eligible':True,'zero':True,'finite':True,'jet':'zero'},
                 {'eligible':True,'zero':None,'finite':True,'jet':'unknown'},
                 {'eligible':True,'zero':False,'finite':None,'jet':'unproved-finite'},
                 {'eligible':False,'zero':False,'finite':True,'jet':'ineligible'},
                 {'eligible':True,'zero':False,'finite':True,'jet':'applicable'}]
        self.assertEqual(NS['select_certified_candidate'](entries)['jet'],'applicable')
        self.assertIsNone(NS['select_certified_candidate'](entries[:-1]))
    def test_actual_address_sum_and_corruption(self):
        # Independent affine row C(eta)*R(eta)*S(eta), source jet=7.
        # Target eta coefficient: C1 R0 S0 + C0 R1 S0 + C0 R0 S1.
        row_terms=[(3,5,11),(2,13,11),(2,5,17)]
        expected=7*(3*5*11+2*13*11+2*5*17)
        rows=[{'row':'A','face':'minus','targetGrade':(1,0),'consumerOriginal':Fraction(c),
               'responsePlaceholder':Fraction(r),'sourceOriginal':Fraction(v),'sourceAtom':7}
              for c,r,v in row_terms]
        unrelated={**rows[0],'row':'B','consumerOriginal':999}
        selected,actual=NS['address_sum_for'](rows+[unrelated],'A','minus',(1,0))
        self.assertEqual(actual,expected);self.assertEqual(len(selected),3)
        self.assertNotEqual(NS['address_sum_for'](rows[:-1],'A','minus',(1,0))[1],expected)
        self.assertNotEqual(NS['address_sum_for'](rows+[rows[0]],'A','minus',(1,0))[1],expected)
        self.assertEqual(NS['address_sum_for'](rows,'A','plus',(1,0)),([],0))
    def test_zero_addresses_have_no_epsilon_power(self):
        self.assertEqual(NS['epsilon_power'](True),0)
        self.assertEqual(NS['epsilon_power'](False),1)
    def test_correction_obligations_have_persistent_inputs(self):
        source=WORKER.read_text()
        for name in ('actual-row-placeholder-input','actual-address-reconstruction',
                     'assembled-mixed-components','whole-native-trace-direct-injection',
                     'inherited-isolated-factor','inherited-full-row-density-proof',
                     'historical-pressure-ablation-scope','physical-depth-branch-conditions',
                     'source-control-candidates','control-address-selection','constant-end-and-zero-profile-scope'):
            self.assertIn(name,source)
        self.assertIn("'reference':sign*context['numeric']['W_0']/2",source)
        self.assertNotIn("'tag':'Dwhole_'+face+'(2,3/2)'",source)
    def test_response_adapter_native_reference_assignment(self):
        native=(HERE.parent/'scripts/S11c_c2_selfenergy_fold_sympy_audit.py').read_text()
        fn=next(n for n in ast.parse(native).body if isinstance(n,ast.FunctionDef) and n.name=='build_face')
        assignments=[n for n in fn.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='reference' for t in n.targets)]
        self.assertEqual(len(assignments),1)
        class Inputs:values={'W_0':Fraction(6,5)}
        namespace={'face':-1,'inputs':Inputs()}
        exec(compile(ast.Module(body=assignments,type_ignores=[]),'native-reference-standin','exec'),namespace)
        self.assertEqual(namespace['reference'],Fraction(-3,5))

    def test_broad_scan_refuses_attribute_and_dynamic_routes(self):
        for text in ("sp.Symbol('delta_p_plus')", "Function(name)", "symbols('delta_p_plus d_w_delta_p_plus')",
                     "Symbol(f'delta_p_{face}')", "sp.Add(Symbol('delta_p_plus'),Integer(0))"):
            rows=NS['broad_row_census'](text)
            self.assertTrue(any(not v['completeNameCoverage'] for v in rows),text)
            self.assertTrue(any(v['unsupportedConstructors'] for v in rows),text)
    def test_broad_scan_accounts_overlapping_name_substrings(self):
        rows=NS['broad_row_census']("Add(Symbol('d_w_delta_p_minus'),Symbol('delta_p_plus'))")
        self.assertEqual([v['rawPressureSubstringCount'] for v in rows],[2,1])
        self.assertTrue(all(v['completeNameCoverage'] for v in rows))
    def test_broad_scan_refuses_hidden_substring(self):
        rows=NS['broad_row_census']("Add(Symbol('u_1'),Str('delta_p_hidden'))")
        self.assertFalse(rows[1]['completeNameCoverage'])
        self.assertEqual(rows[1]['hits'],[])
    def test_full_factor_map_includes_normal_output(self):
        roles={'qi':'q(k)','qo':'q(l)','omega':3}
        mapped=NS['complete_factor_map']({'qi'},{'qo'},roles,[['qi','q(k)']])
        self.assertEqual(mapped,{'qi':'q(k)','qo':'q(l)'})
        # The direct whole response has no raw depth; its normal prefactor still does.
        self.assertEqual(NS['complete_factor_map'](set(),{'qo'},roles,[]),{'qo':'q(l)'})
    def test_full_factor_map_rejects_wrong_or_missing_role(self):
        roles={'qi':'q(k)','qo':'q(l)','omega':3}
        for pairs in ([],[['qi','q(l)']],[['qi','q(k)'],['qo','q(l)']],[['qi','q(k)'],['qi','q(k)']]):
            with self.assertRaises(ValueError):NS['complete_factor_map']({'qi'},{'qo'},roles,pairs)
    def test_off_delta_controls_do_not_use_profile_at_zero(self):
        text=WORKER.read_text();start=text.index('    def choose_source(')
        controls=text[start:text.index("    J.emit('responsive-formal-controls'",start)]
        self.assertNotIn('subs(x,0)',controls)
        self.assertIn("x in item['field'].free_symbols",controls)
        self.assertIn("cncert['status']=='NONCONSTANT_POLYNOMIAL_IN_TANH_CERTIFIED'",controls)
        self.assertIn("'transformValue':'NOT_EVALUATED_OR_CERTIFIED_NONZERO'",controls)
    def test_source_unit_expectation_not_claimed_as_new_proof(self):
        text=WORKER.read_text()
        self.assertNotIn('coefficientDimensionFromSource',text)
        self.assertIn("'coefficientDimensionIndependentlyVerified':False",text)

    def test_metadata_route_class_supported(self):
        native=(HERE.parent/'scripts/S11c_c2_selfenergy_fold_sympy_audit.py').read_text()
        entries=[n for n in ast.parse(native).body if isinstance(n,(ast.ClassDef,ast.FunctionDef)) and n.name=='Inputs']
        self.assertEqual(len(entries),1)
        self.assertIsInstance(entries[0],ast.ClassDef)

if __name__=='__main__':unittest.main()
