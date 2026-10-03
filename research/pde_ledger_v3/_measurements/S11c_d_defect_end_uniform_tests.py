#!/usr/bin/env python3
"""Stdlib metadata/stand-in tests only, never scientific payload decoding."""
import ast,copy,hashlib,json,re,tempfile,unittest
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace
M=Path(__file__).resolve().parent;P=M/'S11c_d_defect_end_uniform.py';source=P.read_text();tree=ast.parse(source)
names=('require','convolution','rotated_name','definitions','exact_structure','save','sha','Evidence','source_fragment','SourceStage','SourceMapUnresolved','attribution_usable','local_key','normalization_routes','merge_native_dimensions')
class Basic:pass
class Matrix:pass
class FunctionClass:pass
class Array:
 def __init__(self,raw,shape=(2,),dtype='int64'):self.raw,self.shape,self.dtype=raw,shape,dtype
 def tobytes(self):return self.raw
 def __eq__(self,other):raise ValueError('ambiguous array truth')
class NpInt(int):
 def item(self):return int(self)
class NpFloat(float):
 def item(self):return float(self)
class NpBool:pass
np=SimpleNamespace(ndarray=Array,complexfloating=complex,integer=NpInt,bool_=NpBool,floating=NpFloat)
sp=SimpleNamespace(Basic=Basic,MatrixBase=Matrix,S=SimpleNamespace(true=object()),srepr=lambda v:'Integer(0)')
ns={'re':re,'np':np,'sp':sp,'FunctionClass':FunctionClass,'Path':Path,'json':json,'hashlib':hashlib,'os':__import__('os'),'base64':__import__('base64'),'ast':ast,'traceback':__import__('traceback'),'ROWS':('U0','U1','U2','THETA_BALANCE','E_W_BALANCE'),'FIELDS':('u_1','u_2','u_3','theta','e_W'),'G':((0,0),(1,0),(0,1),(1,1))}
exec(compile(ast.Module(body=[n for n in tree.body if getattr(n,'name',None) in names],type_ignores=[]),'stdlib-standins','exec'),ns)
MAN=json.loads((M/'S11c_d_defect_end_uniform_inputs.json').read_text())
class Checks(unittest.TestCase):
 def test_convolution_retains_cross_grade(self):
  a={(0,0):Fraction(2),(1,0):Fraction(3),(0,1):Fraction(5)}
  b={(0,0):Fraction(7),(1,0):Fraction(11),(0,1):Fraction(13)}
  self.assertEqual(ns['convolution'](a,b),{(0,0):14,(1,0):43,(0,1):61,(1,1):94})
 def test_full_remainder_terms_not_erased(self):
  full=ns['convolution']({(1,0):2,(0,1):3},{(1,0):5,(0,1):7},False)
  self.assertEqual(full,{(2,0):10,(1,1):29,(0,2):21})
 def test_cyclic_coordinate_and_vector_map(self):
  for a,b in [('u_1','u_3'),('u_2_d1d3','u_1_d2d3'),('theta_d1d2','theta_d1d3'),('grad_theta_3','grad_theta_2'),('u_1_tt','u_3_tt'),('u_2_t_d1','u_1_t_d3'),('w1_profile_d1','w1_profile_d3')]:self.assertEqual(ns['rotated_name'](a,set()),b)
 def test_scalar_not_guessed_from_suffix(self):
  self.assertEqual(ns['rotated_name']('B_rho_3',{'B_rho_3'}),'B_rho_3')
  with self.assertRaises(ValueError):ns['rotated_name']('unclassified_tensor_1',set())
 def test_bad_derivative_refuses(self):
  with self.assertRaises(ValueError):ns['rotated_name']('u_1_d4',set())
 def test_array_identity_never_boolean_conversion(self):
  same=ns['exact_structure'];self.assertTrue(same({'x':Array(b'ab')},{'x':Array(b'ab')}));self.assertFalse(same(Array(b'ab'),Array(b'ac')));self.assertFalse(same(Array(b'ab'),Array(b'ab',dtype='uint64')));self.assertFalse(same([1],(1,)))
 def test_ambiguous_predicate_is_not_true(self):
  with self.assertRaises(ValueError):ns['require'](1,'not literal bool')
 def test_encoded_zero_final_json(self):
  with tempfile.TemporaryDirectory() as d:
   ev=ns['Evidence'](Path(d));value=ev.encode({'zero':Basic(),'dict':{Basic():Basic()}});ns['save'](Path(d)/'x.json',value);self.assertEqual(json.loads((Path(d)/'x.json').read_text()),value)
 def test_append_only_evidence_preserves_before_failure(self):
  with tempfile.TemporaryDirectory() as d:
   ev=ns['Evidence'](Path(d));ev.emit('operand',{'x':7})
   with self.assertRaises(FileExistsError):ev.emit('operand',{'x':8})
   self.assertEqual(json.loads((Path(d)/'operand.json').read_text()),{'x':7});self.assertEqual(ev.count,1)
 def test_actual_native_symbol_classifier_is_bijective(self):
  names=set()
  for row in ('U0','U1','U2','THETA_BALANCE','E_W_BALANCE'):
   raw=json.loads(Path(MAN['savedInputs']['native/'+row+'.json']['path']).read_text())
   nodes=ast.walk(ast.parse(raw['fullConstructor'],mode='eval'))
   names.update(n.args[0].value for n in nodes if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='Symbol')
  physical=json.loads(Path(MAN['physicalInput']).read_text())
  scalars=set(physical['parameters'])|{'epsilon_shape','sigma_W','delta_p_plus','delta_p_minus','d_w_delta_p_plus','d_w_delta_p_minus'}
  mapping={n:ns['rotated_name'](n,scalars) for n in names}
  self.assertEqual(len(names),209);self.assertEqual(set(mapping.values()),names)
  for n in names:self.assertEqual(mapping[mapping[mapping[n]]],n)
 def test_actual_native_unit_and_wave_source_contract(self):
  contract=json.loads(Path(MAN['savedInputs']['native/wave-profile.json']['path']).read_text());c2=Path(MAN['c2Source']).read_text();t=ast.parse(c2)
  dims=next(n.value for n in t.body if isinstance(n,ast.Assign) and any(isinstance(v,ast.Name) and v.id=='DIMENSION_SCHEMA' for v in n.targets))
  self.assertEqual(contract['dimensions'],ast.literal_eval(dims))
  node=next(n for n in t.body if isinstance(n,ast.FunctionDef) and n.name=='wave_jet');self.assertEqual(contract['waveSource'],ast.get_source_segment(c2,node))
 def test_all_current_pins(self):
  for path,h in MAN['sourcePins'].items():self.assertEqual(ns['sha'](path),h,path)
 def test_source_and_returns_are_pinned_before_decode(self):
  self.assertLess(source.index("'opaque old blob '"),source.index('value=codec(payload)'))
  self.assertEqual(len(MAN['uniformReceipts']['completedOperations']),18)
  self.assertEqual(len([r for r in MAN['uniformReceipts']['completedOperations'] if '/exact-limits/' in r['name']]),4)
 def test_only_inert_codec_extracted(self):
  self.assertIn("('SavedCodec','decode')",source);self.assertNotIn('EndBinding(',source)
  executed_calls={ast.unparse(n.func) for n in ast.walk(tree) if isinstance(n,ast.Call)}
  for forbidden in ['sp.Poly','sp.solve','sp.limit','signal.alarm','setitimer','lambdify']:
   self.assertNotIn(forbidden,executed_calls)
 def test_imports_follow_containment(self):
  main=ast.get_source_segment(source,next(n for n in tree.body if getattr(n,'name',None)=='main'))
  self.assertLess(main.index("ns['containment']()"),main.index('import sympy'))
  for n in tree.body:
   if isinstance(n,ast.Import):self.assertFalse(any(a.name in ('sympy','numpy','scipy') for a in n.names))
 def test_every_cell_route_once_and_closed_values_saved(self):
  cells=json.loads(Path(MAN['savedInputs']['ends/symbols.json']['path']).read_text());keys={(c['side'],c['row'],c['field'],tuple(c['grade'])) for c in cells};self.assertEqual(len(keys),200)
  for c in cells:self.assertIn('closedGrazingValue',c);self.assertEqual(c['sumIdentity']['cancelled'],{'text':'0','srepr':'Integer(0)'})
 def test_review_and_scope_authority(self):
  method=json.loads(Path(MAN['methodRecord']).read_text());self.assertTrue(method['jointIndependentMethodClearance'])
  authority=json.loads(Path(MAN['executionAuthority']).read_text());self.assertEqual(authority['scope'],MAN['scope']);self.assertTrue(authority['noDeadline']);self.assertFalse(authority['automaticScientificRetry'])
  self.assertEqual(MAN['resources']['memoryBytes'],4*1024**3);self.assertIsNone(MAN['resources']['durationLimits'])
 def test_all_executable_native_contract_fragments(self):
  sources={'engine':Path(MAN['engineSource']).read_text(),'uniform':Path(MAN['uniformWorker']).read_text(),'ends':Path(MAN['endsWorker']).read_text(),'full':Path(MAN['fullWeakWorker']).read_text()}
  function=next(n for n in tree.body if getattr(n,'name',None)=='source_chart_and_scale');calls=[n for n in ast.walk(function) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='contract']
  self.assertEqual(len(calls),20)
  for n in calls:ns['source_fragment'](sources[n.args[0].id],ast.literal_eval(n.args[1]),ast.literal_eval(n.args[2]))
 def test_source_contract_refuses_changed_time_or_factor(self):
  src='class C:\n def f(self):\n  return x / 2\n'
  ns['source_fragment'](src,('C','f'),'return x / 2')
  with self.assertRaises(ValueError):ns['source_fragment'](src,('C','f'),'return x')
  src='def f():\n phase = -3j * t\n'
  with self.assertRaises(ValueError):ns['source_fragment'](src,('f',),'phase = 3j * t')
 def test_actual_phase_metadata_joins(self):
  phase=json.loads(Path(MAN['savedInputs']['ends/phase-arguments.json']['path']).read_text());duality=json.loads(Path(MAN['savedInputs']['native/weak-duality.json']['path']).read_text())
  self.assertEqual(phase['sourceExpression'],duality['inheritedConvention']['sourcePhase']);self.assertEqual(phase['profileExpression'],duality['inheritedConvention']['profilePhase'])
  for k in ('kout','kin','ko','ki'):self.assertEqual([x['srepr'] for x in phase['bindings'][k][1:]],['Rational(1, 5)','Rational(1, 10)'])
 def test_coordinate_cycle_discriminates_tangent_swap(self):
  old=(Fraction(1,5),Fraction(1,10),Fraction(7));expected=(Fraction(7),Fraction(1,5),Fraction(1,10))
  self.assertEqual((old[2],old[0],old[1]),expected);self.assertNotEqual((old[2],old[1],old[0]),expected);self.assertNotEqual((-old[2],old[0],old[1]),expected)
 def test_half_height_side_bijection(self):
  heights={'minus':Fraction(0),'plus':Fraction(1,2)};old={'LEFT':Fraction(0),'RIGHT':Fraction(1)}
  self.assertEqual({e:next(s for s,h in heights.items() if h==v/2) for e,v in old.items()},{'LEFT':'minus','RIGHT':'plus'})
  self.assertNotEqual({e:next(s for s,h in heights.items() if h==(1-v)/2) for e,v in old.items()},{'LEFT':'minus','RIGHT':'plus'})
 def test_attribution_uses_wave_not_offwave_zero(self):
  f=next(n for n in tree.body if getattr(n,'name',None)=='run_science');text=ast.unparse(f)
  self.assertIn('attribution = wave_test(',text);self.assertNotIn("J.zero(nm + '-attribution'",text)
  self.assertIn("'newLawRemainder': newc['remainder']",text)
 def test_unavailable_results_and_fallback_are_persisted(self):
  self.assertIn("UNRESOLVED_UNSUPPORTED_DEPTH_DEPENDENCE",source);self.assertIn("UNRESOLVED_OPERAND_FINITE_DOMAIN",source);self.assertIn("INAPPLICABLE_NO_CERTIFIED_MOVEMENT",source)
  self.assertIn("J.emit(end + '/control-curl-entry-ablation'",ast.unparse(tree))
  self.assertNotIn("require(regular,",source)
 def test_memory_evidence_uses_actual_constructor_sizes(self):
  sizes=[len(json.loads(Path(MAN['savedInputs']['native/'+r+'.json']['path']).read_text())['fullConstructor'].encode()) for r in ('U0','U1','U2','THETA_BALANCE','E_W_BALANCE')]
  self.assertEqual(sizes,[118565,118565,118565,199688,340239]);self.assertEqual(sum(sizes),895622)
  self.assertNotIn("'sourceRows': native",ast.unparse(tree));self.assertIn('peakRssKiBBeforeCovariance',source)
 def test_raw_invariant_and_grazing_scope(self):
  self.assertIn("Iold=map(old",source);self.assertIn("inv['rawResidual']",source);self.assertIn("R['onWaveInvariant']",source)
  self.assertIn("g=C*LL-LL*DD",source);self.assertIn("UNAVAILABLE_NO_NEW_REMAINDER_LIMIT_COMPUTED",source)
 def test_actual_saved_normalization_routes(self):
  raw={a:json.loads(Path(v['path']).read_text()) for a,v in MAN['savedInputs'].items() if a.startswith('normalization/') or a=='native/U0.json'}
  local,ends,returns=ns['normalization_routes'](raw)
  self.assertEqual((len(local),len(ends),len(returns)),(400,400,16))
  self.assertEqual([returns[i]['mappedGrades'][0]['value']['srepr'] for i in (122,113,114)],['Integer(-9)','Rational(1, 25)','Rational(1, 100)'])
 def test_changed_saved_coefficient_is_refused(self):
  raw={a:json.loads(Path(v['path']).read_text()) for a,v in MAN['savedInputs'].items() if a.startswith('normalization/') or a=='native/U0.json'}
  raw['normalization/local-cells.json'][0]['coefficient']={'text':'-357/20','srepr':'Rational(-357, 20)'}
  with self.assertRaisesRegex(ValueError,'actual saved cell'):ns['normalization_routes'](raw)
 def test_changed_native_witness_is_refused(self):
  raw={a:json.loads(Path(v['path']).read_text()) for a,v in MAN['savedInputs'].items() if a.startswith('normalization/') or a=='native/U0.json'}
  raw['normalization/batch-return.json'][1]['sourceConstructor']='Integer(1)'
  with self.assertRaisesRegex(ValueError,'actual child'):ns['normalization_routes'](raw)
 def test_saved_batch_receipts_are_exact(self):
  raw=json.loads(Path(MAN['savedInputs']['normalization/operation-index.json']['path']).read_text());r=next(x for x in raw if x['name']=='U0-local-batch-007')
  for slot,alias in [('input','batch-input'),('result','batch-return')]:
   pin=MAN['savedInputs']['normalization/'+alias+'.json'];self.assertEqual(r[slot],{'path':Path(pin['path']).name,'sha256':pin['sha256'],'bytes':pin['bytes']})
 def test_all_actual_end_local_terms_have_saved_ancestry(self):
  local=json.loads(Path(MAN['savedInputs']['normalization/local-end-terms.json']['path']).read_text());ends=json.loads(Path(MAN['savedInputs']['ends/symbols.json']['path']).read_text())
  lm={ns['local_key'](x['savedCell']):x for x in local}
  for c in ends:
   keys=[(c['row'],c['field'],n,tuple(c['grade'])) for n in range(4)]
   self.assertEqual(c['local'],[lm[k]['newTerms'][c['side']] for k in keys])
   self.assertEqual(c['localAncestry'],[{'row':k[0],'field':k[1],'xOrder':k[2],'grade':list(k[3])} for k in keys])
 def test_witness_uses_actual_original_epsilon_derivative(self):
  t=ast.unparse(next(n for n in tree.body if getattr(n,'name',None)=='saved_normalization_joins'))
  self.assertIn("original = v['original']",t);self.assertIn('bound_child = original.xreplace(',t);self.assertIn('extracted = sp.diff(bound_child, eps)',t)
  self.assertIn("cell['summands'][position]",t);self.assertNotIn('normalization_trial',source)
 def test_nonzero_attribution_is_fatal(self):
  for status in ('NONZERO_CERTIFIED','UNRESOLVED_NONZERO_SYMBOLIC_REMAINDER','UNKNOWN'):
   for remainder in (True,False):
    with self.assertRaisesRegex(ValueError,'reconstruction failed'):ns['attribution_usable'](status,remainder)
 def test_only_declared_attribution_unavailability_survives(self):
  self.assertFalse(ns['attribution_usable']('UNRESOLVED_UNSUPPORTED_DEPTH_DEPENDENCE',True))
  self.assertFalse(ns['attribution_usable']('ZERO',False));self.assertTrue(ns['attribution_usable']('ZERO',True))
 def test_source_failure_classification_preserves_exception(self):
  class Journal:
   active='source';records=[]
   def emit(self,name,record):self.records.append((name,record))
  j=Journal()
  with self.assertRaises(ns['SourceMapUnresolved']):
   with ns['SourceStage'](j,'units'):raise ValueError('unit mismatch')
  self.assertEqual(j.records[-1][1]['status'],'SOURCE_MAP_UNRESOLVED');self.assertEqual(j.records[-1][1]['exceptionType'],'ValueError')
 def test_source_stage_does_not_swallow_resource_failure(self):
  with self.assertRaises(MemoryError):
   with ns['SourceStage'](SimpleNamespace(),'units'):raise MemoryError('cap')
 def test_source_native_checks_are_inside_classifying_stage(self):
  f=next(n for n in tree.body if getattr(n,'name',None)=='run_science')
  blocks=[ast.unparse(n) for n in ast.walk(f) if isinstance(n,ast.With)]
  main=next(t for t in blocks if "'native-source-chart-units'" in t)
  self.assertIn('rotated_name(',main);self.assertIn('rowunits =',main);self.assertIn("'actual native unit schema'",main)
  raw=next(t for t in blocks if "'-native-source-binding-units'" in t);self.assertIn('originjoin =',raw);self.assertIn("'strong native row/field units'",raw)

 def unit_fixture(self):
  class Atom:
   is_Number=False;is_Symbol=True
   def __init__(self,name):self.name=name
  class Journal:
   def __init__(self):self.records=[]
   def emit(self,n,v):self.records.append((n,v))
  return Atom,Journal
 def test_old_dimension_walker_refuses_actual_gamma_atom(self):
  Atom,Journal=self.unit_fixture();fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='dimension')
  scope={'sp':SimpleNamespace(I=complex(0,1),Rational=Fraction),'dimensions':{'epsilon_shape':(0,0,0)},'zero_dim':(0,0,0),'require':ns['require']}
  exec(compile(ast.Module(body=[fn],type_ignores=[]),'unchanged-unit-predicate','exec'),scope)
  with self.assertRaisesRegex(ValueError,'native dimensions gamma_'):scope['dimension'](Atom('gamma_s11cb_mu_r_bg_05'))
 def test_unit_merge_covers_actual_names_without_guessing_values(self):
  Atom,Journal=self.unit_fixture();base=json.loads(Path(MAN['savedInputs']['native/wave-profile.json']['path']).read_text())['dimensions'];names=set()
  for row in ('U0','U1','U2','THETA_BALANCE','E_W_BALANCE'):
   raw=json.loads(Path(MAN['savedInputs']['native/'+row+'.json']['path']).read_text());names.update(n.args[0].value for n in ast.walk(ast.parse(raw['fullConstructor'],mode='eval')) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='Symbol')
  missing=names-set(base);self.assertTrue(missing);self.assertTrue(all(n.startswith('gamma_') for n in missing))
  atoms=[Atom(n) for n in sorted(names)];units={a:(i+1,-2,3) for i,a in enumerate(atoms) if a.name in missing};j=Journal()
  merged=ns['merge_native_dimensions'](base,atoms,{'saved1':units,'saved2':dict(units)},j)
  fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='dimension');scope={'sp':SimpleNamespace(I=complex(0,1),Rational=Fraction),'dimensions':merged,'zero_dim':(0,0,0),'require':ns['require']};exec(compile(ast.Module(body=[fn],type_ignores=[]),'unchanged-unit-predicate','exec'),scope)
  for a in atoms:self.assertEqual(scope['dimension'](a),tuple(units[a] if a.name in missing else base[a.name]))
  self.assertEqual(len(j.records),len(missing)+1)
 def test_unit_merge_missing_or_assumption_mismatch_refuses(self):
  Atom,Journal=self.unit_fixture();a=Atom('gamma_x');b=Atom('gamma_x')
  for registry in ({},{b:(1,2,3)}):
   with self.assertRaisesRegex(ValueError,'missing actual saved'):ns['merge_native_dimensions']({},[a],{'saved':registry},Journal())
 def test_unit_merge_conflicting_or_inexact_refuses(self):
  Atom,Journal=self.unit_fixture();a=Atom('gamma_x')
  with self.assertRaisesRegex(ValueError,'conflicting'):ns['merge_native_dimensions']({},[a],{'left':{a:(1,2,3)},'right':{a:(2,2,3)}},Journal())
  for v in ((1.0,2,3),(True,2,3),(1,2),None):
   with self.assertRaises(ValueError):ns['merge_native_dimensions']({},[a],{'saved':{a:v}},Journal())
 def test_no_unknown_non_gamma_exemption(self):
  Atom,Journal=self.unit_fixture();a=Atom('unknown_material')
  with self.assertRaisesRegex(ValueError,'non-generated'):ns['merge_native_dimensions']({},[a],{'saved':{a:(1,2,3)}},Journal())
 def test_actual_saved_unit_registry_receipts_and_source_routes(self):
  for name,o in MAN['unitRegistryOrigins'].items():
   d=json.loads(Path(MAN['savedInputs'][o['checksAlias']]['path']).read_text());self.assertEqual(d['objectsSha256'],MAN['originalPackets'][name]['sha256'])
   source=Path(o['producerSource']).read_text();ns['source_fragment'](source,o['scope'],o['statement'])
   self.assertEqual(d['provenance']['producerSources']['scripts/S11c_b_exports.py'],ns['sha'](MAN['bSource']))
 def test_dimension_predicate_and_comparison_functions_unchanged(self):
  old=ast.parse((Path(MAN['worker']).parents[3]/'_scratch/s11c/s11c-defect-end-uniform-20261003/build-review-r3/packet/worker.py').read_text())
  for name in ('dimension','polynomial','coefficients','wave_test','finite_constant','inspect_matrix','saved_normalization_joins'):
   a=next(n for n in ast.walk(old) if isinstance(n,ast.FunctionDef) and n.name==name);b=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==name)
   self.assertEqual(ast.dump(a,include_attributes=False),ast.dump(b,include_attributes=False))

if __name__=='__main__':unittest.main()
