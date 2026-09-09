#!/usr/bin/env python3
"""Synthetic-only N6 comparator fixtures derived from the two emitter schemas.

No measured stream is read; no measured builder is imported or executed.
"""
import contextlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import sympy as sp
import S11c_c2_N6_cross_engine_comparator as c

CASE=(('ANCHORING','LAB_HELD'),('DENSITY','RHO4_CONSTANT'))
DIM={'computed':[[1,-1,0]],'consistent':True,'epsilon_support':[0],'target':[1,-1,0]}
WDIM='<|"COEFFICIENT_SUPPORT" -> {{1,-1,0}}, "MEASURE_SUPPORT" -> {{0,0,0}}, "TRIAL_TEST_SUPPORT" -> {0,0,0}, "PAIRED_SUPPORT" -> {{1,-1,0}}, "ADDITION_CONSISTENCY" -> Inactive[Equal][1,1]|>'


def sparse(values):
    nonzero=[(i+1,v) for i,v in enumerate(values) if v]
    ptr=[0]; n=0
    for v in values:
        n+=bool(v); ptr.append(n)
    cols=','.join('{1}' for _ in nonzero)
    vals=','.join(str(v) for _,v in nonzero)
    return 'SparseArray[Automatic,{%d,1},0,{1,{{%s},{%s}},{%s}}]'%(len(values),','.join(map(str,ptr)),cols,vals)


def slot(column,expr='x',values=(7,0)):
    return '<|'+column+' -> <|"ARITHMETIC" -> {'+expr+'}, "DIMENSIONS" -> {'+WDIM+'}, "COMPONENT_AXES" -> {1}, "PROBE_NUMERATORS" -> '+sparse(values)+', "PROBE_DENOMINATORS" -> SparseArray[Automatic,{2,1},1,{1,{{0,0,0},{}},{}}], "SAMPLE_INDEX" -> {"own points"}|>|>'


def schema(family):
    if family in c.CARRIERS:
        return ['THETA','delta_p_plus',1,[0,0]], '{"THETA",1,deltaPPlus,{1,0,0}}'
    if family in c.RC_SOURCES: return [1,[0,0]], '{1,"theta",{1,0,0}}'
    if family in c.COV_SOURCES: return [1,'theta',[0,0]], '{1,"theta",{1,0,0}}'
    if family in c.GUARDS:
        if '_SLOT_' in family: return ['THETA','linear',[0,0]], '{"EULERIAN",4,1,{1,0,0}}'
        return ['THETA',1,[0,0]], '{"EULERIAN",4,1,"FOURIER_KOUT_Y",{1,0,0}}'
    return [['THICKNESS_TO_TRANSVERSE','THETA'],[0,0],6,1], '{"THETA_FROM_TRANSVERSE","FOURIER_KOUT_Y",{1,0,0},1}'


def records(family,column,expr='x',observation=True,probe=None):
    node={'op':'variable','args':[{'literal':[expr,'global']}],'degree':[1,0]}
    row={'object':'S11CC2_'+family,'anchoring':'LAB_HELD','density':'RHO4_CONSTANT',
         'probe':probe,'dimension':[DIM],'data':{'columns':[column],
         'nonzero_modular_numerator':[observation],'numerator_denominator':[[[7 if observation else 0,1]]]}}
    nd=dict(row,object='S11CC2_'+family+'_NODES',data={'columns':[column],'root_ids':[0],'root_nodes':[node]})
    dag=dict(row,object='S11CC2_'+family.split('_')[0]+'_ARITHMETIC_DAG',data={'nodes':[node]})
    return row,nd,dag


class Mechanics(unittest.TestCase):
    def setUp(self):
        self.saved=c.ACTIVE_NAMES
        self.saved_points=c.POINT_NAMES
        c.POINT_NAMES={}
        c.ACTIVE_NAMES=c.base.checked_mechanical_symbol_map({'delta_p_plus','e_W','x'})
    def tearDown(self):
        c.ACTIVE_NAMES=self.saved
        c.POINT_NAMES=self.saved_points

    def test_axes(self):
        self.assertEqual(dict(c.make_key(CASE)),dict(CASE))
        for items in ([('ANCHORING','LAB_HELD')]*2,[('UNKNOWN','x')],[('DENSITY','LAB_HELD')]):
            with self.assertRaises(c.AxisError): c.make_key(items)

    def test_wl_equals_and_lazy_address(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'wl'
            p.write_text('WL_S11CC2_N6COV_R_COV = <|{"LAB_HELD","RHO4_CONSTANT"} -> '+slot(schema('N6COV_R_COV')[1])+'|>\n')
            idx=c.load_wl(p); key,ref=next(c.wl_case_addresses(idx['N6COV_R_COV'][0]))
            self.assertEqual(key,c.make_key(CASE)); self.assertTrue(ref.read().startswith('<|'))
            p.write_text('WL_S11CC2_N6COV_R_COV: <||>\n')
            with self.assertRaises(c.InputError): c.load_wl(p)

    def test_guard_probe_duplicate_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'py'
            rows=[records('N6_SLOT_GUARD_NATIVE',['THETA','linear',[0,0]],probe=probe)[0]
                  for probe in ('EULERIAN','MATERIAL')]
            p.write_text(''.join(json.dumps(r)+'\n' for r in rows))
            idx=c.load_py_jsonl([p])
            self.assertEqual(len(idx.records),2)
            self.assertTrue(all(len(x)==1 for x in idx.records.values()))

    def test_dag_refs_and_formal_jet(self):
        x=sp.Symbol('x',real=True); f=sp.Function('f'); expr=sp.Derivative(f(x),x,evaluate=False)
        nodes=[{'op':'variable','args':[{'literal':[sp.srepr(expr),'formal_jet']}]},
               {'op':'number','args':[{'literal':2}]},
               {'op':'pow','args':[{'ref':1},{'literal':-1}]},
               {'op':'mul','args':[{'ref':0},{'ref':2}]},
               {'op':'add','args':[{'ref':3},{'ref':1}]}]
        self.assertEqual(c.DAG({'nodes':nodes}).get(4),expr/2+2)
        with self.assertRaisesRegex(c.InputError,'ARITHMETIC_DAG'): c.DAG(None)
        for literal in (['algebraic_i'],['w1_profile','X'],['formal_integral_guard',6]):
            self.assertIsInstance(c.dag_literal(literal),sp.Basic)

    def test_every_numeric_family_schema_extracts(self):
        for family in sorted(c.SHARED-c.META):
            with self.subTest(family=family),contextlib.redirect_stdout(io.StringIO()):
                col,wcol=schema(family); row,nd,dag=records(family,col,probe='EULERIAN' if family in c.GUARDS else None)
                a=c.extract_py_numeric(family,row,nd,c.DAG(dag['data']),CASE)
                b=c.extract_wl_numeric(family,slot(wcol),CASE)
                self.assertEqual(len(a),1); self.assertEqual(len(b),1)
                self.assertTrue(c.materialize(a[0])); self.assertTrue(c.materialize(b[0]))

    def test_every_metadata_schema_extracts(self):
        fixtures={
          'N6COV_FROZEN_PHI':({'substitution_map':[['theta','theta+x']], 'a_rho':0},'<|"MAP" -> {theta -> theta+x}, "A_RHO" -> 0|>'),
          'N6COV_PHI_DOMAIN_CENSUS':({'coverage':[['theta',True,['theta',[],None,None]]],'uncovered':[],'imported_mu_sha256':'a'*64},'<|"COVERAGE" -> <|theta -> Inactive[Equal][1,1]|>, "UNCOVERED" -> {}|>'),
          'N6COV_ACTUAL_CONTROL_PARAMETERS':({'kappa_a':1,'kappa_j':0,'junk_dimension':[-1,-2,1]},'<|"ADVECTION" -> 1, "JUNK" -> 0, "JUNK_DIMENSIONS" -> {{-1,-2,1}}|>'),
          'N6RC_ADVECTION_ABSENCE':({'a_rho':0,'density_gradient':[0,0,0],'material_mu_tag_derivatives':[0,0]},'<|"ADVECTION" -> 0, "DENSITY_GRADIENT" -> {0,0,0}, "MATERIAL_MU_TAG_DERIVATIVE" -> 0|>'),
          'N6RC_FROZEN_RELATIONS':({'field_and_derivative_maps':[['theta','theta+x']],'fixed_axes':['LAB_HELD','RHO4_CONSTANT']},'<|"FIELD_MAP" -> {theta -> theta+x}, "MATERIAL_COVECTOR_MAP" -> {{1,0},{0,1}}|>'),
          'N6RC_DIMENSIONS':({'eulerian':{"('THETA', 'delta_p_plus')":{'computed':[1,-1,0],'consistent':True,'unknown':[]}}},'<|"OBJECTS" -> <|"RC:CARRIER_EULERIAN" -> <|{"THETA",1,deltaPPlus,{1,0,0}} -> {'+WDIM+'}|>|>, "BASE_OPERAND" -> deltaPPlus|>')}
        self.assertEqual(set(fixtures),c.META)
        for family,(a,b) in fixtures.items():
            with self.subTest(family=family),contextlib.redirect_stdout(io.StringIO()):
                aa=c.extract_meta('SymPy',family,a,CASE); bb=c.extract_meta('WL',family,b,CASE)
                self.assertGreater(len(aa),0); self.assertGreater(len(bb),0)
                for leaf in aa+bb: self.assertTrue(c.materialize(leaf),leaf.error)

    def test_spelling_only_bridge(self):
        with self.assertRaises(c.InputError): c.base.checked_mechanical_symbol_map({'alpha_beta','alphaBeta'})
        col,wcol=schema('N6RC_R_N6')
        self.assertNotEqual(c.decode_py_key('N6RC_R_N6',col),c.decode_wl_key('N6RC_R_N6',wcol))
        p=c.decode_py_key('N6COV_R_COV',[0,'theta',[0,1]])
        w=c.decode_wl_key('N6COV_R_COV','{"SUM","theta",{1,0,1}}')
        self.assertEqual(p,w)

    def test_no_applied_head_collapse_or_predicate_evaluation(self):
        self.assertEqual(c.parse_wl_value('f[x]'),sp.Function('f')(sp.Symbol('x')))
        self.assertNotEqual(c.parse_wl_value('muThetaOperand[x]'),sp.Symbol('mu_theta_drive'))
        self.assertNotIsInstance(c.parse_wl_value('Inactive[Greater][1,0]'),bool)
        self.assertNotIsInstance(c.parse_wl_value('Inactive[Equal][1,1]'),bool)
        self.assertEqual(c.parse_wl_value('Inactive[Plus][x,Inactive[Times][-1,y]]'),sp.Symbol('x')-sp.Symbol('y'))

    def test_pit_seal_and_one_sided_support(self):
        family='N6COV_R_COV'; col,wcol=schema(family); row,nd,dag=records(family,col,observation=False)
        out=io.StringIO()
        with contextlib.redirect_stdout(out):
            c.seal('SymPy',family,CASE,'primes',[11,17,29]); c.seal('WL',family,CASE,'primes',[11,17,41])
            a=c.extract_py_numeric(family,row,nd,c.DAG(dag['data']),CASE)
            b=c.extract_wl_numeric(family,slot(wcol,values=(3,4,5)),CASE)
            acc=c.compare_family(family,a,b)
        self.assertEqual(acc.join,1)
        self.assertIsInstance(c.support_residual(c.support(False),c.support(True)),c.UndecidedResidual)
        self.assertEqual(c.wl_support(sparse((0,0)),1),[c.NO_NONZERO_FOUND])
        seals=[json.loads(line.split(' ',1)[1]) for line in out.getvalue().splitlines() if line.startswith('SURFACED ')]
        self.assertTrue(seals); self.assertTrue(all('A_minus_B' not in s for s in seals))
        self.assertIn('pit_sealed',out.getvalue()); self.assertIn('NO_NONZERO_FOUND',out.getvalue())

    def test_three_valued_and_carrier_repoint_form(self):
        self.assertIsInstance(c.measured_residual(False,sp.S.Zero,1),c.BooleanNotResidualable)
        self.assertIsInstance(c.measured_residual(c.UndecidedResidual(None,None,'budget'),sp.S.Zero,1),c.UndecidedResidual)
        f='N6RC_CARRIER_EULERIAN'; col,wcol=schema(f)
        outputs=[]
        for expr in ('x','Inactive[Plus][Inactive[Times][-1,x],y*z]'):
            row,nd,dag=records(f,col); out=io.StringIO()
            with contextlib.redirect_stdout(out):
                aa=c.extract_py_numeric(f,row,nd,c.DAG(dag['data']),CASE)
                bb=c.extract_wl_numeric(f,slot(wcol,expr=expr),CASE)
                c.compare_family(f,aa,bb)
            outputs.append([json.loads(x.split(' ',1)[1])['A_minus_B'] for x in out.getvalue().splitlines() if x.startswith('CASE ')][0])
        self.assertEqual(outputs[0],'Integer(0)'); self.assertNotEqual(outputs[0],outputs[1])

    def test_zero_extraction_operational_guard(self):
        a=c.Accounting()
        with contextlib.redirect_stdout(io.StringIO()): c.zero_extract_guard([],[],True,True,a)
        self.assertEqual(a.zero_extract_failures,1)

    def test_sympy_only_and_missing_dag(self):
        self.assertEqual(c.excluded_reason('N6_MU_RECONSTRUCTION_RESIDUAL'),'no WL sibling')
        self.assertNotIn('N6_MU_RECONSTRUCTION_RESIDUAL',c.SHARED)
        f='N6COV_R_COV'; row,nd,dag=records(f,schema(f)[0])
        with contextlib.redirect_stdout(io.StringIO()),self.assertRaisesRegex(c.InputError,'ARITHMETIC_DAG'):
            c.extract_py_numeric(f,row,nd,None,CASE)

    def test_cli_disagreement_exits_zero_and_missing_dag_nonzero(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); p=root/'py'; w=root/'wl'; f='N6COV_R_COV'
            row,nd,dag=records(f,schema(f)[0]); p.write_text(''.join(json.dumps(r)+'\n' for r in (row,nd,dag)))
            w.write_text('WL_S11CC2_'+f+' = <|{"LAB_HELD","RHO4_CONSTANT"} -> '+slot(schema(f)[1],expr='y*z-x')+'|>\n')
            args=[sys.executable,c.__file__,'--py',str(p),'--wl',str(w),'--family',f,'--deferred-ledger',str(root/'deferred')]
            result=subprocess.run(args,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout[-2000:]+result.stderr)
            for word in ('operand_A','operand_B','A_minus_B','ACCOUNTING','RUN_ACCOUNTING'):
                self.assertIn(word,result.stdout)
            p.write_text(''.join(json.dumps(r)+'\n' for r in (row,nd)))
            result=subprocess.run(args,capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0); self.assertIn('parse_failed',result.stdout)

    def test_row_index_origin_and_form_repoint(self):
        f='N6RC_CARRIER_EULERIAN'
        for i in range(3):
            py=c.decode_py_key(f,[f'U{i}','delta_p_plus',1,[0,0]],CASE)
            wl=c.decode_wl_key(f,[f'U{i+1}',1,'deltaPPlus',[1,0,0]],CASE)
            self.assertEqual(py,wl)
            self.assertEqual(dict(py)['ROW_COMPONENT'],str(i+1))
        self.assertNotEqual(c.decode_py_key(f,['U1','delta_p_plus',1,[0,0]]),
                            c.decode_wl_key(f,['U1',1,'deltaPPlus',[1,0,0]]))
        target={'object':'N6RC_CARRIER_EULERIAN','arithmetic':'x'}
        other={'object':'N6RC_CARRIER_MATERIAL','arithmetic':'Inactive[Plus][Inactive[Times][-1,x],y*z]'}
        measured=[]
        for replacement in (target,other):
            out=io.StringIO(); row,nd,dag=records(f,['U0','delta_p_plus',1,[0,0]])
            with contextlib.redirect_stdout(out):
                a=c.extract_py_numeric(f,row,nd,c.DAG(dag['data']),CASE)
                b=c.extract_wl_numeric(f,slot('{"U1",1,deltaPPlus,{1,0,0}}',replacement['arithmetic']),CASE)
                acc=c.compare_family(target['object'],a,b)
            self.assertEqual(acc.join,1)
            measured.append([json.loads(x.split(' ',1)[1])['A_minus_B'] for x in out.getvalue().splitlines() if x.startswith('CASE ')][0])
        self.assertEqual(measured[0],'Integer(0)'); self.assertNotEqual(measured[0],measured[1])

    def test_disjoint_axis_sets_are_counted(self):
        f='N6RC_R_N6'; col,wcol=schema(f);row,nd,dag=records(f,col)
        with contextlib.redirect_stdout(io.StringIO()):
            a=c.extract_py_numeric(f,row,nd,c.DAG(dag['data']),CASE)
            b=c.extract_wl_numeric(f,slot(wcol),CASE)
            acc=c.compare_family(f,a,b)
        self.assertEqual(acc.axis_set_mismatch,2); self.assertEqual(acc.unmatched_key,2)
        self.assertEqual(acc.join,0)

    def test_metadata_and_bare_jet_spelling_injectivity(self):
        index=c.PyIndex()
        index.symbols={'W_bg','W_bg_d1','theta_d1d1','w1_profile_d1','u_1'}
        index.local_literals={('w1_profile_d1','X'),('w1_profile_d1','Y')}
        c.ACTIVE_NAMES,c.POINT_NAMES=c.verified_spelling_maps(index)
        self.assertEqual(c.parse_wl_value('WBgJet1'),sp.Symbol('W_bg_d1'))
        self.assertEqual(c.parse_wl_value('thetaJet11'),sp.Symbol('theta_d1d1'))
        self.assertEqual(c.parse_wl_value('w1ProfileJet1AtY'),c.dag_literal(['w1_profile_d1','Y']))
        self.assertNotEqual(c.parse_wl_value('w1ProfileJet1AtY'),c.parse_wl_value('w1ProfileJet1AtX'))
        a=c.parse_py_expr('W_bg_d1*u_1/W_bg'); b=c.parse_wl_value('WBgJet1*u1/WBg')
        self.assertEqual(c.measured_residual(a,b,1),sp.S.Zero)
        self.assertEqual(c.parse_wl_value('thetaJet11[x]'),sp.Function('thetaJet11')(sp.Symbol('x')))
        index.symbols.add('WBgJet1')
        with self.assertRaisesRegex(c.InputError,'non-injective'):c.verified_spelling_maps(index)
        self.assertIn('W_bg',c.metadata_symbols({'h_alpha':'W_bg_d1*u_1/W_bg'}))

    def test_uncommon_held_heads_are_not_placeholder_equalities(self):
        a=c.parse_wl_value('Inactive[SameQ][x]'); b=c.parse_wl_value('Inactive[Exp][x]')
        self.assertEqual(a.func.__name__,'HeldInactiveSameQ')
        self.assertEqual(b.func.__name__,'HeldInactiveExp')
        self.assertNotEqual(c.measured_residual(a,b,1),sp.S.Zero)
        self.assertNotIn('N6Shield',str(a)+str(b))

    def test_domain_coverage_is_keyed_on_atom(self):
        f='N6COV_PHI_DOMAIN_CENSUS'
        with contextlib.redirect_stdout(io.StringIO()):
            a=c.extract_meta('SymPy',f,{'coverage':[['theta',True,['theta',[],None,None]]]},CASE)
            b=c.extract_meta('WL',f,'<|"COVERAGE" -> <|theta -> Inactive[Equal][1,1]|>|>',CASE)
        aa=next(x for x in a if dict(x.key)['FIELD_PATH']=='["coverage"]')
        self.assertEqual(aa.key,b[0].key)
        self.assertEqual(dict(aa.key)['DOMAIN_ATOM'],'theta')
        c.materialize(aa);c.materialize(b[0])
        self.assertIsInstance(c.measured_residual(aa.value,b[0].value,1),c.BooleanNotResidualable)

    def test_support_uses_components_not_nonzero_count(self):
        a='SparseArray[Automatic,{2,2},0,{1,{{0,1,2},{{1},{2}}},{7,8}}]'
        b='SparseArray[Automatic,{2,2},0,{1,{{0,1,2},{{1},{1}}},{7,8}}]'
        self.assertEqual(c.wl_support(a,2),[c.NONZERO_WITNESSED,c.NONZERO_WITNESSED])
        self.assertEqual(c.wl_support(b,2),[c.NONZERO_WITNESSED,c.NO_NONZERO_FOUND])

    def test_dag_rejects_boolean_number_and_extra_arguments(self):
        for args in ([{'literal':True}],[{'literal':1},{'literal':2}]):
            with self.assertRaises(c.InputError):c.DAG({'nodes':[{'op':'number','args':args}]})

    def test_budget_emission_is_distinct_from_parse_failure(self):
        from types import SimpleNamespace
        import time
        def sleeper(conn,*args): time.sleep(.5)
        with tempfile.TemporaryDirectory() as tmp:
            args=SimpleNamespace(object_seconds=.005,object_rss_mib=8192.,residual_leaf_seconds=1.,
                                 deferred_ledger=Path(tmp)/'ledger.md')
            out=io.StringIO()
            with patch.object(c,'_worker',sleeper),contextlib.redirect_stdout(out):
                acc,peak=c.bounded_object('N6COV_R_COV',('LAB_HELD','RHO4_CONSTANT'),c.PyIndex(),[],args)
            self.assertEqual(acc.deferred_oversize,1); self.assertEqual(acc.parse_failed,0)
            self.assertIn('UndecidedResidual',out.getvalue())
            for word in ('operand_A','operand_B','A_minus_B'):self.assertIn(word,out.getvalue())
            self.assertIn('wallclock ceiling',args.deferred_ledger.read_text())

    def test_nonempty_broken_extractor_and_all_deferred_exit_nonzero(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=root/'py';w=root/'wl';f='N6COV_R_COV'
            row,nd,dag=records(f,schema(f)[0]);p.write_text(''.join(json.dumps(r)+'\n' for r in (row,nd,dag)))
            w.write_text('WL_S11CC2_'+f+' = <|{"LAB_HELD","RHO4_CONSTANT"} -> '+slot(schema(f)[1])+'|>\n')
            args=['--py',str(p),'--wl',str(w),'--family',f,'--deferred-ledger',str(root/'ledger')]
            with patch.object(c,'extract_py_numeric',return_value=[]),patch.object(c,'extract_wl_numeric',return_value=[]),contextlib.redirect_stdout(io.StringIO()):
                self.assertNotEqual(c.run(args),0)
            result=subprocess.run([sys.executable,c.__file__,*args,'--object-seconds','0.000001'],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            summary=json.loads(next(x.split(' ',1)[1] for x in result.stdout.splitlines() if x.startswith('RUN_ACCOUNTING ')))
            self.assertGreater(summary['deferred_oversize'],0); self.assertEqual(summary['parse_failed'],0)

    def test_provenance_is_surfaced_and_diagnostic_is_only_accounted(self):
        self.assertNotEqual(c.excluded_reason('N6COV_PROVENANCE'),'no WL sibling')
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=root/'py';w=root/'wl';f='N6COV_R_COV'
            row,nd,dag=records(f,schema(f)[0]);internal=dict(row,object='S11CC2_N6_MU_RECONSTRUCTION_RESIDUAL')
            p.write_text(''.join(json.dumps(r)+'\n' for r in (row,nd,dag,internal)))
            w.write_text('WL_S11CC2_'+f+' = <|{"LAB_HELD","RHO4_CONSTANT"} -> '+slot(schema(f)[1])+'|>\n')
            result=subprocess.run([sys.executable,c.__file__,'--py',str(p),'--wl',str(w),'--family',f],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            rows=[json.loads(x.split(' ',1)[1]) for x in result.stdout.splitlines() if x.startswith('ACCOUNTING_ROW ')]
            finding=next(x for x in rows if x['family']=='N6_MU_RECONSTRUCTION_RESIDUAL')
            self.assertEqual(finding['reason'],'no WL sibling');self.assertEqual(finding['finding'],'sympy_only')
            self.assertFalse(any('N6_MU_RECONSTRUCTION_RESIDUAL' in x for x in result.stdout.splitlines() if x.startswith('CASE ')))


if __name__=='__main__': unittest.main()
