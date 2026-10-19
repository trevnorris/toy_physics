#!/usr/bin/env python3
"""S9b record lookups only. No engine imports, evaluation, or symbolic/numeric work.

Native field retrieval uses balanced delimiters and returns original substrings.
JSON walks retrieve stored fields and count literal encodings, never test mathematical
predicates. All temporary files belong in _scratch/s9b_record. --write regenerates
only the measurements file; it does not author the physics record or register.
"""
from pathlib import Path
import argparse
import collections
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'research/pde_ledger_v3'
SCRIPT = 'research/pde_ledger_v3/scripts/S9b_record_measurements.py'
NATIVE = {'py': BASE/'scripts/out/S9b_light_bending_sympy_audit.out',
          'wl': BASE/'mathematica/out/S9b_light_bending_mathematica_audit.out'}
CMP = BASE/'scripts/out/S9b_cross_engine_comparator.out'
ACCOUNT = BASE/'scripts/out/S9b_cross_engine_comparator_accounting.jsonl'


def split(s):
    """Split commas outside balanced syntax; do not interpret expressions."""
    out=[]; stack=[]; quote=None; escaped=False; start=0; i=0
    while i<len(s):
        c=s[i]
        if quote:
            if escaped: escaped=False
            elif c=='\\': escaped=True
            elif c==quote: quote=None
        elif c in '\"\'': quote=c
        elif s[i:i+2]=='<|': stack.append('|>'); i+=1
        elif s[i:i+2]=='|>' and stack and stack[-1]=='|>': stack.pop(); i+=1
        elif c in '([{': stack.append({'(':')','[':']','{':'}'}[c])
        elif stack and c==stack[-1]: stack.pop()
        elif c==',' and not stack: out.append(s[start:i].strip()); start=i+1
        i+=1
    if stack or quote: raise ValueError('unbalanced stored syntax')
    if s[start:].strip(): out.append(s[start:].strip())
    return out


def children(s,engine):
    if engine=='wl' and s.startswith('<|') and s.endswith('|>'):
        return {part.split(' -> ',1)[0].strip('"'):part.split(' -> ',1)[1] for part in split(s[2:-2])}
    if engine=='wl' and s.startswith('{') and s.endswith('}'):
        return dict(enumerate(split(s[1:-1])))
    for head,op,cl in [('Tuple','(',')'),('ConditionalExpression','[',']')]:
        if s.startswith(head+op) and s.endswith(cl):
            return dict(enumerate(split(s[len(head)+1:-1])))
    return {}


def native(engine,tag):
    prefix=engine.upper()+('_LOCAL_S9B_' if tag.startswith('LOCAL_') else '_S9B_')
    if tag.startswith('LOCAL_'):tag=tag[6:]
    prefix+=tag+': '
    with NATIVE[engine].open() as f:
        for line in f:
            if line.startswith(prefix):return line[len(prefix):].rstrip('\n')
    raise KeyError(prefix)


def native_fields(engine):
    """Compact original fields. Full gates are separately retrieved, never replaced."""
    with NATIVE[engine].open() as f:
        for lineno,line in enumerate(f,1):
            if ': ' not in line:continue
            tag,s=line.rstrip('\n').split(': ',1)
            if '_LOCAL_' in tag and not any(k in tag for k in ['COUNTING','RETAINED_GRADES','MASS_BALANCE','C_RESPONSES','C_CONSTANT_RESPONSE','C_FIXED_RATIO_RESPONSE','C_POWER_RESPONSE','RADAR_ANGULAR_OPERATOR_RATIO','ABEL_REFERENCE_INVERSES','PROFILE_ANSATZ','GENERAL_PROFILES','SPEED_IDENTIFICATION']):continue
            def walk(v,path):
                ch=children(v,engine)
                if engine=='py' and (v.startswith("Tuple(Str('ATTACHED_DOMAIN')") or v.startswith("Tuple(Str('branch_complement')")):
                    return
                if ch:
                    for k,val in ch.items():
                        if engine=='wl' and k in ['Domain','BranchDomain','AlongRayDomain','OutsideDomain','AdmissibleRadialDomain','FunctionSpace']:continue
                        if engine=='wl' and v.startswith('ConditionalExpression[') and k==1:continue
                        walk(val,path+[str(k)])
                elif len(v)<20000:
                    print(f'{tag} line={lineno} path={"/".join(path)}\n{v}')
            walk(s,[])


def comparator():
    counts=collections.Counter(); samples=collections.Counter(); reps=collections.Counter(); exceptions=collections.Counter(); symbolic=collections.Counter(); boolean=collections.Counter(); residual30_rows=collections.Counter()
    with CMP.open() as f:
        for lineno,line in enumerate(f,1):
            d=json.loads(line)
            if d.get('kind')=='configuration':
                print('CONFIGURATION',line.rstrip());continue
            if d.get('kind')!='comparison':
                if d.get('kind')=='completion':print('COMPLETION',line.rstrip())
                continue
            name=d['name']; local=collections.Counter(); differences=[]; representations=[]
            def walk(x,path=''):
                if isinstance(x,dict):
                    if 'same_bound_serialization' in x:
                        reps[str(x['same_bound_serialization'])]+=1
                        if x['same_bound_serialization'] is False:representations.append(path)
                    if x.get('outcome')=='symbolic':
                        symbolic[(x.get('role',''),x.get('zero_class',''),x.get('sample_policy','sampled'))]+=1
                        print('SYMBOLIC_VALUE',name,path,json.dumps({k:v for k,v in x.items() if k in ['outcome','role','py_value','wl_value','residual','zero_class','sample_policy','comparison_formed']},ensure_ascii=False,sort_keys=True))
                    if x.get('outcome')=='predicate':
                        boolean[(x.get('role',''),x.get('boolean_residual','<absent>') if x.get('boolean_residual') in ('False','True') else '<stored expression>')]+=1
                        print('PREDICATE_LOCATION',name,path,json.dumps({k:v for k,v in x.items() if k in ['role','comparison_formed']},sort_keys=True))
                    if 'samples' in x:
                        for i,s in enumerate(x['samples']):
                            c=s.get('comparison',{}); r=c.get('residual',{}); enc=json.dumps({k:v for k,v in c.items() if k!='residual'},sort_keys=True)
                            samples[(x.get('role',''),enc)]+=1; local[(x.get('role',''),enc)]+=1
                            # Literal exceptions, not a numeric/tolerance comparison.
                            if (isinstance(r,dict) and any(r.get(k,'(0.0 + 0.0j)') != '(0.0 + 0.0j)' for k in ('computed_18','computed_30'))) or c.get('boolean_difference') is True or c.get('same_exact_representation') is False:
                                differences.append({'path':path+'/samples/'+str(i),'sample':s}); exceptions[name]+=1
                    if 'samples' in x:
                        for sample in x['samples']:
                            residual=sample.get('comparison',{}).get('residual',{})
                            if isinstance(residual,dict) and 'computed_30' in residual:
                                counts[(x.get('role',''),residual['computed_30'])]+=1
                                if residual['computed_30']!='(0.0 + 0.0j)':residual30_rows[name]+=1
                    if x.get('outcome') in ['label','marker']:
                        print('LABEL_OR_MARKER',name,path,json.dumps({k:v for k,v in x.items() if k in ['outcome','comparison','py_text','wl_text']},ensure_ascii=False,sort_keys=True))
                    for k,v in x.items():
                        if k not in ['samples','py_operand','wl_operand']:
                            walk(v,path+'/'+k)
                elif isinstance(x,list):
                    for i,v in enumerate(x):walk(v,path+'/'+str(i))
            walk(d['comparison'])
            print('ROW',name,'line',lineno,'SAMPLE_LITERAL_COUNTS',json.dumps([[*k,v] for k,v in local.items()],ensure_ascii=False))
            print('REPRESENTATION_FALSE_PATHS',name,json.dumps(representations))
            for diff in differences:print('SAMPLE_LITERAL_EXCEPTION',name,json.dumps(diff,ensure_ascii=False,sort_keys=True))
    print('SAMPLE_LITERAL_TOTALS',json.dumps([[*k,v] for k,v in samples.items()],ensure_ascii=False))
    print('RESIDUAL_LITERAL_TOTALS',json.dumps([[*k,v] for k,v in counts.items()],ensure_ascii=False))
    print('REPRESENTATION_LITERAL_TOTALS',json.dumps(reps,sort_keys=True))
    print('SAMPLE_LITERAL_EXCEPTION_COUNTS',json.dumps(exceptions,sort_keys=True))
    print('SYMBOLIC_CLASS_LITERAL_COUNTS',json.dumps([[*k,v] for k,v in symbolic.items()]))
    print('BOOLEAN_RESIDUAL_LITERAL_COUNTS',json.dumps([[*k,v] for k,v in boolean.items()]))
    print('COMPUTED30_NONZERO_LITERAL_ROWS',json.dumps(residual30_rows,sort_keys=True))
    print('COMPUTED30_NONZERO_LITERAL_COUNTS',sum(residual30_rows.values()),len(residual30_rows))


def section(name):
    if name=='M1':
        paths=['directives/S9b_record_directive.md','directives/S9b_repair_decision_list.md',
               'directives/_measurements/S9b_production_runs.md','directives/_measurements/S9b_comparator_production_run.md']
        for p in paths:print('FILE',p);print((BASE/p).read_text())
        print('SPEC A-C');print((BASE/'directives/S9b_SHARED_PHYSICS.md').read_text().split('## Part D inputs:')[0])
    elif name=='M2':
        for engine in ('py','wl'):
            for tag in ('BRANCH_EXISTENCE','PATH_TRAVERSAL','BRANCH_TYPE'):
                print(engine,tag);print(native(engine,tag))
    elif name=='M3':native_fields('wl')
    elif name=='M4':native_fields('py')
    elif name=='M5':comparator()
    elif name=='M6':print(ACCOUNT.read_text(),end='')
    elif name=='M7':
        for p in ['S9b_comparator_build_report.md','S9b_comparator_build_report_coverage_reasons.md','S9b_comparator_build_report_function_name_table.md']:
            print('FILE',p);print((BASE/'directives/_measurements'/p).read_text())
    elif name=='M8':
        for stem,n in [('S9b_repair_build_r',12),('S9b_comparator_build_r',5)]:
            for i in range(n):
                p=BASE/'directives/_measurements'/f'{stem}{i}_review_disposition.md'
                print('FILE',p.relative_to(ROOT));print(p.read_text())
    elif name=='M9':
        p=BASE/'steps/O2_steady_brane_balance.md';print('FILE',p.relative_to(ROOT));print(p.read_text())
    elif name=='M12':
        pins=[('34a994d5','directives/S9b_SHARED_PHYSICS.md'),('79055918','directives/S9b_repair_decision_list.md'),('72866fcf','steps/O2_steady_brane_balance.md'),('503736ff','scripts/S9b_light_bending_sympy_audit.py'),('503736ff','mathematica/S9b_light_bending_mathematica_audit.wl'),('88bb1f41','scripts/S9b_cross_engine_comparator.py')]
        for rev,path in pins:
            old=subprocess.run(['git','show',rev+':research/pde_ledger_v3/'+path],cwd=ROOT,check=True,capture_output=True).stdout
            print(rev,path,'literal_bytes_match',old==(BASE/path).read_bytes())
        print('PY B_GAMMA_DIFFERENCE on-domain value (original substring)')
        print(children(children(native('py','B_GAMMA_DIFFERENCE'),'py')[0],'py')[1])
    elif name=='M11':
        print(NATIVE['wl'].read_text(),end='')
    elif name=='M10':
        p='research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md'
        old=subprocess.run(['git','show','HEAD:'+p],cwd=ROOT,check=True,capture_output=True,text=True).stdout
        new=(ROOT/p).read_text()
        for label,s in [('HEAD',old),('WORKTREE',new)]:
            rows=[l for l in s.splitlines() if l.startswith('### R-')]
            print(label,'ENTRY_COUNT',len(rows));print('\n'.join(rows))
        print('WORKTREE REGISTER');print(new)
    else:raise ValueError(name)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--section');ap.add_argument('--native',nargs='+');ap.add_argument('--write',action='store_true');a=ap.parse_args()
    if a.native:
        engine,tag,*path=a.native;v=native(engine,tag)
        for key in path:
            ch=children(v,engine);v=ch[int(key) if key.isdecimal() else key]
        print(v)
    elif a.section:section(a.section)
    elif a.write:
        dst=BASE/'steps/_measurements/S9b_record_measurements.md'
        scratch=ROOT/'_scratch/s9b_record';scratch.mkdir(parents=True,exist_ok=True)
        temp=scratch/'S9b_record_measurements.md'
        with temp.open('w') as out:
            out.write('# S9b record measurements\n\nGenerated by `scripts/S9b_record_measurements.py`; retrieval only.\n'
                      'Native compact listings preserve retrieved substrings; Named domain fields and scalar fields of 20000 characters or more are omitted from the compact views, not evaluated or asserted; full native retrieval remains available. '
                      'M2 retains the general gates in full; every complete tagged object remains retrievable with `--native ENGINE TAG`. '
                      'M5 counts literal JSON fields, without numeric interpretation, tolerance or algebra. '
                      'M11 preserves the full Wolfram native stream, including every restricted/forward/response domain; compact fields do not replace those gates. M6/M7 preserve the complete row accounting and remaining-coverage inventory.\n\n')
            for i in range(1,13):
                key='M'+str(i);cmd=[sys.executable,SCRIPT,'--section',key]
                out.write(f'## {key}\n\n~~~~~~text\n$ python3 {SCRIPT} --section {key}\n');out.flush()
                subprocess.run(cmd,cwd=ROOT,stdout=out,check=True)
                out.write('~~~~~~\n\n')
        temp.replace(dst)
    else:ap.error('choose --section, --native, or --write')
if __name__=='__main__':main()
