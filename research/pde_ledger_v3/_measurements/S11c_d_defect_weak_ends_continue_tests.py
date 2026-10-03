"""Stdlib metadata/source and exact-true interface tests; no scientific restoration."""
from pathlib import Path
import ast,json,runpy,symtable,unittest
M=Path(__file__).resolve().parent;P='S11c_d_defect_weak_ends_continue';B='S11c_d_defect_weak_ends'
ns=runpy.run_path(str(M/(P+'.py')),run_name='stdlib_test_only');source=(M/(B+'.py')).read_text()
prior=Path('/var/projects/toy_physics/_scratch/s11c/s11c-defect-weak-ends-20261002/diagnostic-01/complete')
read=lambda n:json.loads((prior/n).read_text())
class TrueAtom:
 def __bool__(self):return True
class Unknown:
 def __bool__(self):raise AssertionError('truth coercion forbidden')
 def __eq__(self,other):raise AssertionError('equality coercion forbidden')
class LyingEqual:
 def __eq__(self,other):return True
 def __bool__(self):return True
class Tests(unittest.TestCase):
 def test_python_true(self):self.assertTrue(ns['known_true'](True,TrueAtom()))
 def test_exact_true_atom(self):
  atom=TrueAtom();self.assertTrue(ns['known_true'](atom,atom))
 def test_false_numeric_and_text_refused(self):
  atom=TrueAtom()
  for x in (False,None,0,1,'True',[],{},TrueAtom(),LyingEqual()):self.assertFalse(ns['known_true'](x,atom))
 def test_unknown_not_coerced(self):self.assertFalse(ns['known_true'](Unknown(),TrueAtom()))
 def test_original_failure_reproduced(self):
  atom=TrueAtom();self.assertTrue(bool(atom))
  with self.assertRaises(ValueError):ns['require'](atom and atom,'control speed in scoped interval')
  ns['require'](ns['known_true'](atom,atom),'same exact-true decision')
 def test_original_tail_identical(self):
  setup,tail,nested=ns['source_parts'](source);run=next(n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef) and n.name=='run_science')
  start=next(i for i,n in enumerate(run.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and len(n.value.args)==2 and isinstance(n.value.args[1],ast.Constant) and n.value.args[1].value=='control speed in scoped interval')
  self.assertEqual(ast.dump(ast.Module(body=tail,type_ignores=[])),ast.dump(ast.Module(body=run.body[start+1:],type_ignores=[])))
  self.assertEqual({n.name for n in nested},{'zero','constant'})
 def test_unpublished_setup_only(self):
  setup,_,_=ns['source_parts'](source)
  self.assertEqual([ast.unparse(n) for n in setup],['controls = []','point = {p: sp.Integer(1), q: sp.Integer(2)}','cs2 = sp.Rational(180, 101)'])
 def test_no_dispersion_or_prefix_replay_in_tail(self):
  _,tail,_=ns['source_parts'](source);text=ast.unparse(ast.Module(body=tail,type_ignores=[]))
  for word in ('control-outgoing-physical-dispersion','new-field-endpoints','run_science','native_phase_exponent'):self.assertNotIn(word,text)
 def test_tail_namespace(self):
  _,tail,_=ns['source_parts'](source)
  text='def tail():\n'+''.join('    '+l+'\n' for l in ast.unparse(ast.Module(body=tail,type_ignores=[])).splitlines())
  st=symtable.symtable(text,'tail','exec').get_children()[0]
  names={s.get_name() for s in st.get_symbols() if s.is_global() and s.is_referenced()}
  allowed={'possible_controls','sp','W','Bh','reversed_limits','point','audit','constant','tuple','next','symbols','zero','J','len','list','require','controls','cs2','join_records','sorted','used','copies','seen','H','manifest'}
  self.assertLessEqual(names,allowed)
 def test_actual_completed_operation_joins(self):
  idx=read('artifact-index.json');ops=read('operation-index.json');self.assertEqual(len(ops),34)
  for o in ops:
   self.assertTrue(o['name'].startswith('new-field-endpoints-'))
   for k in ('input','result'):self.assertEqual(o[k],idx[o[k]['path']])
   self.assertEqual(read(o['input']['path']),read(o['result']['path'])['inherited'])
 def test_actual_saved_candidate_counts(self):
  tree=ast.parse(source);nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='control_eligible'];local={'require':ns['require']};exec(compile(ast.Module(body=nodes,type_ignores=[]),'metadata-only','exec'),local)
  counts={k:0 for k in ('omit-height-contact','reverse-translation-phase','omit-lower-normal-sign')}
  for row in ns['ROWS']:
   for r in read(row+'-new-end-addresses.json'):
    for k in counts:counts[k]+=r['endContributions']['plus']!=ns['ZERO'] and local['control_eligible'](r['inheritedAddress'],k)
  self.assertEqual(counts,{'omit-height-contact':48,'reverse-translation-phase':48,'omit-lower-normal-sign':24})
 def test_all_saved_cells_remain(self):
  cells=read('new-complete-retained-weak-end-symbols.json');self.assertEqual(len(cells),200)
  self.assertTrue(all(c['sumIdentity']['cancelled']==ns['ZERO'] for c in cells))
 def test_no_science_at_import(self):
  text=(M/(P+'.py')).read_text();tree=ast.parse(text)
  for n in tree.body:
   if isinstance(n,ast.Import):self.assertFalse(any(x.name in ('sympy','numpy','scipy') for x in n.names))
  main=ast.get_source_segment(text,next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main'))
  self.assertLess(main.index("ns['containment']()"),main.index('import sympy'))
 def test_interval_persisted_before_guard(self):
  text=(M/(P+'.py')).read_text()
  self.assertLess(text.index("J.emit('interval-predicate-operands'"),text.index("require(known_true(lower"))
  self.assertNotIn('bool(lower)',text);self.assertNotIn('bool(upper)',text)
 def test_fresh_hook_and_guard(self):
  text=(M/(P+'_launch.py')).read_text();self.assertIn('20261002/continuation-01',text);self.assertIn('s11c_guarded_run.py',text);self.assertIn('defect_weak_ends_continue',text)
  self.assertLess(text.index("if path.exists() and read(path)['status'] == 'waiting':"),text.index("os.write(writer, b'1')"))
 def test_no_deadline_or_retry(self):
  text=(M/(P+'.py')).read_text();self.assertNotIn('signal.alarm',text);self.assertNotIn('setitimer',text)
  self.assertIn("'automaticRetry':False",text);self.assertIn("continuationIndependentBuildClearance",text)
if __name__=='__main__':unittest.main()
