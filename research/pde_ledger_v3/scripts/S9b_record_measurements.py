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
import functools
import gzip
import hashlib
import json
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'research/pde_ledger_v3'
SCRIPT = 'research/pde_ledger_v3/scripts/S9b_record_measurements.py'
BUILD = 'fb2d6d25'
COMPARATOR = '88bb1f41'
BASELINE = '36f59ff1'
REGISTER_BEFORE = 'caa0710d'
ROUND_BEFORE = 'd1511b20'


@functools.lru_cache(None)
def revision_bytes(rev, path):
    """Read the named revision, resolving annex content only in this checkout."""
    data = subprocess.run(['git', 'show', rev+':'+path], cwd=ROOT,
                          check=True, capture_output=True).stdout
    mode = subprocess.run(['git', 'ls-tree', rev, '--', path], cwd=ROOT,
                          check=True, capture_output=True, text=True).stdout.split()[0]
    if mode == '120000':
        target = (ROOT/path).parent/data.decode().strip()
        target = target.resolve()
        if not target.is_relative_to(ROOT/'.git/annex'):
            raise ValueError('annex target outside this checkout: '+path)
        if not target.is_file():
            subprocess.run(['datalad', 'get', path], cwd=ROOT, check=True)
        return target
    return data


def revision_text(rev, path):
    data = revision_bytes(rev, path)
    return data.read_text() if isinstance(data, Path) else data.decode()


def revision_file(rev, path):
    data = revision_bytes(rev, path)
    if isinstance(data, Path):
        return data
    scratch = ROOT/'_scratch/s9b_record/pinned'
    scratch.mkdir(parents=True, exist_ok=True)
    target = scratch/(rev+'_'+Path(path).name)
    target.write_bytes(data)
    return target


def native_file(engine):
    path = {'py': 'scripts/out/S9b_light_bending_sympy_audit.out',
            'wl': 'mathematica/out/S9b_light_bending_mathematica_audit.out'}[engine]
    return revision_file(BUILD, 'research/pde_ledger_v3/'+path)


def source(rev, path):
    print('SOURCE git show '+rev+':'+path)
    print(revision_text(rev, path))


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
    with native_file(engine).open() as f:
        for line in f:
            if line.startswith(prefix):return line[len(prefix):].rstrip('\n')
    raise KeyError(prefix)


def native_fields(engine):
    """Compact original fields. Full gates are separately retrieved, never replaced."""
    with native_file(engine).open() as f:
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
    with revision_file(COMPARATOR, 'research/pde_ledger_v3/scripts/out/S9b_cross_engine_comparator.out').open() as f:
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


def condition_pairs():
    """Compare only original strings at stored syntax paths, never conditions."""
    for engine in ('py', 'wl'):
        objects = {}
        with native_file(engine).open() as stream:
            for line in stream:
                tag, sep, payload = line.rstrip('\n').partition(': ')
                if sep and ('_CONDITION_DEFLECTION' in tag or '_CONDITION_RADAR' in tag):
                    objects[tag] = payload
        for tag, left in objects.items():
            if not tag.endswith('_CONDITION_DEFLECTION'):continue
            other = tag[:-len('DEFLECTION')]+'RADAR'
            right = objects[other]
            print('PAIR', tag, other, 'whole_payload_literal_equal', left == right)
            equal=[]; different=[]
            def walk(a,b,path):
                if a==b:
                    equal.append(path)
                    print('IDENTICAL_FIELD',path,'bytes',len(a.encode()),'sha256',hashlib.sha256(a.encode()).hexdigest())
                    if len(a)<20000:print(a)
                    return
                ac,bc=children(a,engine),children(b,engine)
                if ac and bc and ac.keys()==bc.keys():
                    for k in ac:walk(ac[k],bc[k],path+'/'+str(k))
                else:
                    different.append(path)
                    print('DIFFERING_FIELD',path,'DEFLECTION',a,'RADAR',b,sep='\n')
            # Even identical whole payloads get their immediate stored fields listed.
            lc,rc=children(left,engine),children(right,engine)
            if lc and rc and lc.keys()==rc.keys():
                for k in lc:walk(lc[k],rc[k],str(k))
            else:walk(left,right,'<payload>')
            print('PAIR_SUMMARY',tag,json.dumps({'identical_paths':equal,'differing_paths':different},sort_keys=True))


def harness_retrieval():
    prefix = 'research/pde_ledger_v3/directives/_measurements/'
    path = prefix+'S9b_repair_sympy_tag_inventory.json'
    print('SOURCE git show '+BUILD+':'+path)
    inventory = json.loads(revision_text(BUILD, path))
    print('PY flag_qualification', inventory['flag_qualification'])
    print('PY source_pins', json.dumps(inventory['source_pins'], sort_keys=True))
    print('PY K11', json.dumps(inventory['K11'], sort_keys=True))
    for run in inventory['runs']:
        if run['id']=='full_harness':
            print('PY harness provenance',json.dumps({k:v for k,v in run.items() if k!='tags'},sort_keys=True))
            for row in run['tags']:
                if 'A_NONRECIPROCAL_' in row['name']:print('PY HARNESS_TAG',json.dumps(row,sort_keys=True))
    for knife, data in inventory['comparison_counts_by_knife'].items():
        grouped=collections.defaultdict(list)
        for row in data['tags']:
            if 'A_NONRECIPROCAL_' in row['tag']:
                print('PY STORED_MARK',knife,json.dumps(row,sort_keys=True))
                grouped[row.get('mark', 'NOT_EVALUATED' if row.get('evaluated') is False else '<absent>')].append(row['tag'])
        print('PY SHARED_NONRECIPROCAL_MARKS',knife,json.dumps(grouped,sort_keys=True))
    for knife in inventory['knives']:
        print('PY STORED_WITNESS_LOCATION',json.dumps(knife,sort_keys=True))
    path=prefix+'S9b_repair_wolfram_round11_harness.stdout.gz'
    print('SOURCE gzip text of git show '+BUILD+':'+path+' (annex content)')
    marks=collections.defaultdict(lambda:collections.defaultdict(list))
    with gzip.open(revision_file(BUILD,path),'rt') as stream:
        for lineno,line in enumerate(stream,1):
            tag,sep,payload=line.rstrip('\n').partition(': ')
            match=re.fullmatch(r'WL_LOCAL_S9B_ABLATION_(K[0-9]+[abc]?)_(WL_S9B_(?:F_FLOW_|F_LIVE_)?A_NONRECIPROCAL_.+)',tag)
            if match:
                knife,obj=match.groups()
                # Last printed structural mark; no mathematical comparison is performed.
                found=re.findall(r'"Mark" -> "([^"]+)"',payload)
                mark=found[-1] if found else '<absent>'
                marks[knife][mark].append(obj)
                suffix=payload[payload.rfind('"StructuralComparison"'):]
                print('WL HARNESS_TAG',lineno,tag,'Baseline_field', '"Baseline" ->' in payload,'Corrupted_field','"Corrupted" ->' in payload,'NOT_EVALUATED_literal','NOT_EVALUATED' in payload)
                print(suffix if found else payload)
            elif 'K11' in tag and ('COVERAGE' in tag or 'MANIFEST' in tag):
                print('WL K11 COVERAGE',lineno,line.rstrip())
            elif ('WITNESS' in tag or 'SEARCH' in tag) and 'A_NONRECIPROCAL_' in payload:
                print('WL STORED_WITNESS',lineno,line.rstrip())
    for knife,grouped in marks.items():
        print('WL SHARED_NONRECIPROCAL_MARKS',knife,json.dumps(grouped,sort_keys=True))



def source_lines(rev, path, first, last):
    print('SOURCE git show '+rev+':'+path, 'lines', first, last)
    for n,line in enumerate(revision_text(rev,path).splitlines(),1):
        if first<=n<=last:print(str(n)+(':' if not line else ': '+line))


def gamma_strata():
    source_lines(BASELINE, 'research/pde_ledger_v3/directives/S9b_record_directive.md',109,113)
    for engine in ('py','wl'):
        for prefix in ('','R_FLOW_ONLY_','R_SPEED_ONLY_','R_TILT_ONLY_','F_FLOW_','F_LIVE_'):
            for suffix in ('DEFLECTION','RADAR','DIFFERENCE'):
                tag=prefix+'B_GAMMA_'+suffix
                raw=native(engine,tag)
                print(engine,tag,'COMPLETE_OBJECT');print(raw)
                obj=children(raw,engine)
                if engine=='py':value=children(obj[0],engine)[1]
                else:
                    if 'Object' in obj:obj=children(obj['Object'],engine)
                    value=children(obj['OnDomain'],engine)[0]
                print(engine,tag,'ON_DOMAIN_VALUE');print(value)
    source_lines(BUILD,'research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py',610,663)
    source_lines(BUILD,'research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl',430,460)


def restriction_objects():
    for engine in ('py','wl'):
        for prefix in ('R_FLOW_ONLY_B_','R_SPEED_ONLY_B_','R_TILT_ONLY_B_',
                       'R_BULK_ONLY_C_FIXED_RATIO_','R_BULK_ONLY_C_POWER_'):
            for object_name in ('CONDITION','IMPLIED_JN'):
                for observable in ('DEFLECTION','RADAR'):
                    tag=prefix+object_name+'_'+observable
                    print(engine,tag,'COMPLETE_OBJECT');print(native(engine,tag))
    for path,first,last in [('scripts/S9b_light_bending_sympy_audit.py',435,490),
                            ('scripts/S9b_light_bending_sympy_audit.py',910,982),
                            ('mathematica/S9b_light_bending_mathematica_audit.wl',488,578)]:
        source_lines(BUILD,'research/pde_ledger_v3/'+path,first,last)


def branch_objects():
    for path,first,last in [('scripts/S9b_light_bending_sympy_audit.py',78,87),
                            ('scripts/S9b_light_bending_sympy_audit.py',527,575),
                            ('scripts/S9b_light_bending_sympy_audit.py',666,724),
                            ('scripts/S9b_light_bending_sympy_audit.py',770,773),
                            ('scripts/S9b_light_bending_sympy_audit.py',910,923),
                            ('scripts/S9b_light_bending_sympy_audit.py',1041,1062),
                            ('scripts/S9b_light_bending_sympy_audit.py',1070,1077),
                            ('mathematica/S9b_light_bending_mathematica_audit.wl',120,185),
                            ('mathematica/S9b_light_bending_mathematica_audit.wl',240,282),
                            ('mathematica/S9b_light_bending_mathematica_audit.wl',650,670)]:
        source_lines(BUILD,'research/pde_ledger_v3/'+path,first,last)
    names=[]
    for prefix in ('','F_FLOW_','F_LIVE_'):
        for suffix in ('BRANCH_EXISTENCE','PATH_TRAVERSAL','BRANCH_TYPE'):
            tag=prefix+suffix;names.append(tag)
            pair=[]
            for engine in ('py','wl'):
                raw=native(engine,tag);pair.append(raw)
                print(engine,tag,'COMPLETE_OBJECT');print(raw)
            print('SAME_NAMED_NATIVE_PAYLOADS',tag,'literal_equal',pair[0]==pair[1])
    with revision_file(COMPARATOR,'research/pde_ledger_v3/scripts/out/S9b_cross_engine_comparator_accounting.jsonl').open() as stream:
        for line in stream:
            row=json.loads(line)
            if row.get('name') in names:print('ACCOUNTING',line.rstrip())


def section(name):
    print('PINNED_READS native/harness='+BUILD+' comparator='+COMPARATOR+' register_before='+REGISTER_BEFORE+' repository_context='+BASELINE)
    if name=='M1':
        paths=['directives/S9b_record_directive.md','directives/S9b_repair_decision_list.md',
               'directives/_measurements/S9b_production_runs.md','directives/_measurements/S9b_comparator_production_run.md']
        for p in paths:
            rev = '79055918' if p.endswith('S9b_repair_decision_list.md') else (BUILD if p.endswith('S9b_production_runs.md') else (COMPARATOR if p.endswith('S9b_comparator_production_run.md') else BASELINE))
            source(rev, 'research/pde_ledger_v3/'+p)
        print('SPEC A-C git show 34a994d5:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md');print(revision_text('34a994d5', 'research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md').split('## Part D inputs:')[0])
    elif name=='M2':
        for engine in ('py','wl'):
            for tag in ('BRANCH_EXISTENCE','PATH_TRAVERSAL','BRANCH_TYPE'):
                print(engine,tag);print(native(engine,tag))
    elif name=='M3':native_fields('wl')
    elif name=='M4':native_fields('py')
    elif name=='M5':comparator()
    elif name=='M6':source(COMPARATOR, 'research/pde_ledger_v3/scripts/out/S9b_cross_engine_comparator_accounting.jsonl')
    elif name=='M7':
        for p in ['S9b_comparator_build_report.md','S9b_comparator_build_report_coverage_reasons.md','S9b_comparator_build_report_function_name_table.md']:
            source(COMPARATOR, 'research/pde_ledger_v3/directives/_measurements/'+p)
    elif name=='M8':
        for stem,n in [('S9b_repair_build_r',12),('S9b_comparator_build_r',5)]:
            for i in range(n):
                p=BASE/'directives/_measurements'/f'{stem}{i}_review_disposition.md'
                source(BUILD if stem=='S9b_repair_build_r' else COMPARATOR, str(p.relative_to(ROOT)))
    elif name=='M9':
        source('72866fcf', 'research/pde_ledger_v3/steps/O2_steady_brane_balance.md')
    elif name=='M12':
        pins=[('34a994d5','directives/S9b_SHARED_PHYSICS.md'),('79055918','directives/S9b_repair_decision_list.md'),('72866fcf','steps/O2_steady_brane_balance.md'),('503736ff','scripts/S9b_light_bending_sympy_audit.py'),('503736ff','mathematica/S9b_light_bending_mathematica_audit.wl'),('88bb1f41','scripts/S9b_cross_engine_comparator.py')]
        for rev,path in pins:
            other = BUILD if path.endswith(('_audit.py', '_audit.wl')) else BASELINE
            old=revision_text(rev, 'research/pde_ledger_v3/'+path)
            new=revision_text(other, 'research/pde_ledger_v3/'+path)
            print('git show',rev+':research/pde_ledger_v3/'+path,'versus git show',other+':research/pde_ledger_v3/'+path,'literal_bytes_match',old==new)
        print('PY B_GAMMA_DIFFERENCE on-domain value (original substring)')
        print(children(children(native('py','B_GAMMA_DIFFERENCE'),'py')[0],'py')[1])
    elif name=='M11':
        print(native_file('wl').read_text(),end='')
    elif name=='M10':
        p='research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md'
        old=revision_text(REGISTER_BEFORE, p)
        new=(ROOT/p).read_text()
        prior=revision_text(ROUND_BEFORE,p)
        for label,s in [(REGISTER_BEFORE,old),(ROUND_BEFORE,prior),('WORKTREE',new)]:
            rows=[l for l in s.splitlines() if l.startswith('### R-')]
            print(label,'ENTRY_COUNT',len(rows));print('\n'.join(rows))
        print('BEFORE REGISTER git show '+REGISTER_BEFORE+':'+p);print(old)
        print('WORKTREE REGISTER');print(new)
        def sourced_entries(text):
            found=[]
            for block in re.split(r'(?m)^### ',text)[1:]:
                if not block.startswith('R-'):continue
                entry=block.split(' ',1)[0]
                source_block=block.split('- **source**',1)[1].split('- **requirement**',1)[0]
                if 'S9b' in source_block:found.append(entry)
            return found
        for label,text in [(REGISTER_BEFORE,old),(ROUND_BEFORE,prior),('WORKTREE',new)]:
            print(label,'S9B_SOURCE_ENTRIES',json.dumps(sourced_entries(text)))
        print('S9B_SOURCE_GAINS_SINCE_'+REGISTER_BEFORE,json.dumps([x for x in sourced_entries(new) if x not in sourced_entries(old)]))
        print('S9B_SOURCE_GAINS_SINCE_'+ROUND_BEFORE,json.dumps([x for x in sourced_entries(new) if x not in sourced_entries(prior)]))
    elif name=='M13':condition_pairs()
    elif name=='M14':harness_retrieval()
    elif name=='M15':
        source(BASELINE, 'docs/medium_requirements_and_prior_art.md')
        source(BASELINE, 'research/pde_ledger_v3/directives/O2_steady_brane_balance_scoping.md')
        source(BASELINE, 'research/pde_ledger_v3/directives/prior_art_4d_bulk_search_prompt.md')
    elif name=='M16':
        for engine in ('py','wl'):
            for tag in ('F_MASS_SOLUTION','C_CONSTANT_N_DEPENDENCE','C_FIXED_RATIO_N_DEPENDENCE','C_POWER_N_DEPENDENCE'):
                print(engine,tag);print(native(engine,tag))
        for path,first,last in [('scripts/S9b_light_bending_sympy_audit.py',1031,1050),
                                ('mathematica/S9b_light_bending_mathematica_audit.wl',633,636),
                                ('mathematica/S9b_light_bending_mathematica_audit.wl',651,664)]:
            full='research/pde_ledger_v3/'+path
            print('SOURCE git show '+BUILD+':'+full, 'lines', first, last)
            for n,line in enumerate(revision_text(BUILD,full).splitlines(),1):
                if first<=n<=last:print(str(n)+': '+line)
    elif name=='M17':gamma_strata()
    elif name=='M18':restriction_objects()
    elif name=='M19':branch_objects()
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
                      'M11 preserves the full Wolfram native stream, including every restricted/forward/response domain; compact fields do not replace those gates. M6/M7 preserve the complete row accounting and remaining-coverage inventory. '
                      'M13 performs literal within-engine condition-payload comparisons only; M14 retrieves harness coverage/marks; M15 preserves repository prior-art sources; M16 retrieves mass and N-dependence objects by engine. '
                      'M17 retrieves gamma objects and their strata; M18 retrieves complete single-mechanism conditions/exchange; M19 retrieves branch constructors and same-named native rows. '
                      'All revision-bearing reads use explicit git revisions, including caa0710d for the pre-S9b register; only the after-pass register is read from the worktree. HEAD is not an input.\n\n')
            for i in range(1,20):
                key='M'+str(i);cmd=[sys.executable,SCRIPT,'--section',key]
                out.write(f'## {key}\n\n~~~~~~text\n$ python3 {SCRIPT} --section {key}\n');out.flush()
                subprocess.run(cmd,cwd=ROOT,stdout=out,check=True)
                out.write('~~~~~~\n\n')
        temp.replace(dst)
    else:ap.error('choose --section, --native, or --write')
if __name__=='__main__':main()
