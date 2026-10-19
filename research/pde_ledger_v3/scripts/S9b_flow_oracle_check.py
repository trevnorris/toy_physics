#!/usr/bin/env python3
"""Mechanical transcript reader and launcher. Run only through s11c_guarded_run.py.
No CAS is performed in this driver; the single math child inherits the guard.
The AST conversion below is syntax-only, never an independently typed operand.
"""
import ast
import argparse
import hashlib
import json
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path('/var/projects/toy_physics')
HERE = ROOT / '_scratch/s9b_oracle_check'
SCRIPT = pathlib.Path(__file__).resolve()
ASSETS = SCRIPT.parent / 'S9b_flow_oracle_check_data'
SOURCES = {
 'PY': ROOT / 'research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out',
 'WL': ROOT / 'research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out',
}
PY_TAGS = ['PY_LOCAL_S9B_GENERAL_PROFILES', 'PY_LOCAL_S9B_COUNTING',
 'PY_LOCAL_S9B_VECTOR_DISPERSION', 'PY_LOCAL_S9B_VECTOR_ADVECTION',
 'PY_LOCAL_S9B_SPEED_IDENTIFICATION', 'PY_LOCAL_S9B_MASS_BALANCE']
WL_TAGS = ['WL_LOCAL_S9B_SUPPLIED_INPUTS', 'WL_LOCAL_S9B_PROFILE_ANSATZ',
 'WL_LOCAL_S9B_RETAINED_GRADES', 'WL_LOCAL_S9B_RESTRICTION_FLOW_ONLY']
for route in ['DEFLECTION', 'RADAR']:
 PY_TAGS.append('PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_' + route)
 WL_TAGS.append('WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_' + route)


def sha(path):
 h = hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda: f.read(1024*1024), b''): h.update(b)
 return h.hexdigest()


def read_tags(path, wanted):
 result = {}
 with path.open() as f:
  for n, line in enumerate(f, 1):
   tag, sep, raw = line.rstrip('\n').partition(': ')
   if sep and tag in wanted:
    if tag in result: raise ValueError('duplicate tag ' + tag)
    result[tag] = {'line': n, 'raw': raw, 'source': str(path)}
 if set(result) != set(wanted): raise ValueError('missing tags ' + str(set(wanted)-set(result)))
 return result


# Wolfram symbol names cannot contain the SymPy underscores. This reversible
# lexical encoding does not identify either engine's physical objects.
def name(s):
 return 'py`' + s.replace('_', '$')


def wl(n):
 if isinstance(n, ast.Constant):
  if isinstance(n.value, str): return json.dumps(n.value)
  return str(n.value)
 if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub):
  return 'Minus['+wl(n.operand)+']'
 if isinstance(n, (ast.List, ast.Tuple)):
  return '{' + ','.join(wl(a) for a in n.elts) + '}'
 if isinstance(n, ast.Name):
  return {'oo':'Infinity', 'true':'True', 'false':'False'}.get(n.id, name(n.id))
 if not isinstance(n, ast.Call): raise ValueError(ast.dump(n))
 if isinstance(n.func, ast.Call): return wl(n.func) + '[' + ','.join(wl(a) for a in n.args) + ']'
 h=n.func.id
 a=n.args
 if h in ['Symbol', 'Function']: return name(ast.literal_eval(a[0]))
 if h=='Integer': return str(ast.literal_eval(a[0]))
 if h=='Str': return json.dumps(ast.literal_eval(a[0]))
 if h=='Limit':
  direction=ast.literal_eval(a[3].args[0]) if len(a)>3 else '+'
  return 'Inactive[Limit]['+wl(a[0])+','+wl(a[1])+'->'+wl(a[2])+',Direction->'+('1' if direction=='-' else '-1')+']'
 if h=='Subs':
  return '('+wl(a[0])+'/.Thread['+wl(a[1])+'->'+wl(a[2])+'])'
 if h=='Derivative': return 'D['+','.join(wl(x) for x in a)+']'
 heads={'Tuple':'List','Add':'Plus','Mul':'Times','Pow':'Power','Rational':'Rational',
 'Equality':'Equal','Unequality':'Unequal','StrictGreaterThan':'Greater',
 'StrictLessThan':'Less','GreaterThan':'GreaterEqual','LessThan':'LessEqual',
 'ImmutableDenseMatrix':'Identity','Integral':'Inactive[Integrate]', 'acos':'ArcCos',
 'ConditionSet':'sourceConditionSet','Interval':'sourceInterval','ProductSet':'sourceProductSet',
 'Intersection':'sourceIntersection','Union':'sourceUnion','Contains':'sourceContains'}
 return heads.get(h,h)+'['+','.join(wl(x) for x in a)+']'


def node_record(n):
 return {'raw': ast.unparse(n), 'wl': wl(n)}


def declared_assumptions(tree):
 # Interpret the actual srepr keyword declarations, never infer from names.
 declarations={}
 for n in ast.walk(tree):
  constructor=n.func if isinstance(n,ast.Call) else None
  if isinstance(constructor,ast.Call) and isinstance(constructor.func,ast.Name) and constructor.func.id=='Function':
   attrs={}
   for kw in constructor.keywords:
    if kw.arg is None: attrs.update(ast.literal_eval(kw.value))
    else: attrs[kw.arg]=ast.literal_eval(kw.value)
  elif isinstance(constructor,ast.Name) and constructor.id=='Symbol':
   attrs={kw.arg:ast.literal_eval(kw.value) for kw in n.keywords}
  else: continue
  expression=wl(n)
  predicates=[]
  for prop,head in [('positive','Greater'),('negative','Less'),('nonnegative','GreaterEqual'),('nonpositive','LessEqual'),('zero','Equal'),('nonzero','Unequal')]:
   if attrs.get(prop) is True: predicates.append(head+'['+expression+',0]')
  if attrs.get('real') is True or attrs.get('positive') is True or attrs.get('negative') is True:
   predicates.append('Element['+expression+',Reals]')
  for predicate in predicates: declarations[predicate]={'wl':predicate,'raw':ast.unparse(n)}
 return list(declarations.values())


def attached_domain(tree):
 # First tuple contains the exterior interval and the general ray predicate.
 radial=tree.args[0].args[1]
 interval=tree.args[0].args[2].args[1]
 lower,upper,left_open,right_open=interval.args
 relation='Greater' if left_open.id=='true' else 'GreaterEqual'
 result=[{'wl':relation+'['+wl(radial)+','+wl(lower)+']','raw':ast.unparse(interval)}]
 if not isinstance(upper,ast.Name) or upper.id!='oo':
  relation='Less' if right_open.id=='true' else 'LessEqual'
  result.append({'wl':relation+'['+wl(radial)+','+wl(upper)+']','raw':ast.unparse(interval)})
 result.append(node_record(tree.args[0].args[3].args[1]))
 return result


def extract_py(record):
 a=ast.parse(record['raw'], mode='eval').body
 meta={ast.literal_eval(x.args[0].args[0]):x.args[1] for x in a.args[2].args}
 cases=[]
 for case in a.args[1].args:
  fields=case.args
  eqs=[x for x in fields if isinstance(x,ast.Call) and isinstance(x.func,ast.Name) and x.func.id=='Equality']
  cases.append({'label':ast.literal_eval(fields[0].args[0]),
                'domain_raw':ast.unparse(fields[1]),
                'raw':ast.unparse(case),
                'velocity':node_record(eqs[0]), 'exchange':node_record(eqs[-1]),
                'domain':node_record(fields[1].args[1]),
                'velocity_derivative':node_record(eqs[1]) if len(eqs)>2 else None})
 return {'raw':record['raw'], 'line':record['line'], 'source':record['source'],
         'metadata_raw':ast.unparse(a.args[2]), 'domain_raw':ast.unparse(a.args[0]),
         'declared':declared_assumptions(a),'attached':attached_domain(a),
         'constants':node_record(meta['constants']),
         'definition':node_record(meta['definition']),
         'squared':node_record(meta['squared']),
         'general':node_record(meta['general_solution']), 'cases':cases}


def main():
 parser=argparse.ArgumentParser()
 parser.add_argument('--control',choices=['baseline','visser_form','hl_spatial_form'],default='baseline')
 args=parser.parse_args()
 data={'provenance':{}, 'py':{}, 'wl':{}, 'quotes':json.loads((ASSETS/'oracle_quotes.json').read_text()),
       'map':json.loads((ASSETS/'identification_map.json').read_text()),'control':args.control}
 for engine, path in SOURCES.items():
  annex_target=subprocess.check_output(['git','show','503736ff:'+str(path.relative_to(ROOT))],cwd=ROOT,text=True).strip()
  if path.readlink().as_posix()!=annex_target:
   raise ValueError('source annex target differs from accepted revision 503736ff: '+str(path))
  data['provenance'][engine]={'path':str(path),'sha256':sha(path),'accepted_revision':'503736ff','annex_target':annex_target}
  records=read_tags(path,PY_TAGS if engine=='PY' else WL_TAGS)
  for tag, record in records.items():
   if engine=='PY':
    if 'IMPLIED_JN' in tag: data['py'][tag]=extract_py(record)
    else:
     tree=ast.parse(record['raw'],mode='eval').body
     record['wl']=wl(tree)
     record['declared']=declared_assumptions(tree)
     data['py'][tag]=record
   else: data['wl'][tag]=record
 for f in ['gr-qc_0411060.pdf','gr-qc_9712010.pdf']:
  data['provenance'][f]={'path':str(HERE/'sources'/f),'sha256':sha(HERE/'sources'/f)}
 for f in ['S9b_flow_oracle_check.py','S9b_flow_oracle_check.wl','oracle_quotes.json','identification_map.json']:
  artifact=SCRIPT if f==SCRIPT.name else ASSETS/f
  data['provenance'][f]={'path':str(artifact),'sha256':sha(artifact)}
 (HERE/'inputs.json').write_text(json.dumps(data,indent=2)+'\n')
 print('WORKING_ON: native Wolfram computation from both transcript representations',flush=True)
 messages=[]
 completed=False
 with (HERE/('latest_kernel_'+args.control+'.out')).open('w') as output:
  proc=subprocess.Popen(['math','-script',str(ASSETS/'S9b_flow_oracle_check.wl')],stdin=subprocess.DEVNULL,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
  for line in proc.stdout:
   print(line,end='',flush=True)
   output.write(line)
   output.flush()
   if re.match(r'^[A-Za-z$][A-Za-z0-9$]*::',line): messages.append(line.rstrip())
   if line.startswith('INSTRUMENT_COMPLETION: '): completed=True
  code=proc.wait()
 print('DRIVER_COMPLETION_GUARD:',json.dumps({'kernel_exit':code,'completion_record':completed,'messages':messages}),flush=True)
 return code or (0 if completed and not messages else 2)

if __name__=='__main__':
 try: sys.exit(main())
 except Exception as e:
  print('SOURCE_OR_DRIVER_ERROR:',repr(e),flush=True)
  print('OPERAND: UNMATCHED; RESIDUAL: UNMATCHED; GUARD: source/driver failure',flush=True)
  raise
