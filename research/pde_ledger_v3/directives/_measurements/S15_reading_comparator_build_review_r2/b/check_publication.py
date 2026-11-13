"""Replay the reporting stage and perturb one native rule, plane, point or source argument at a time."""
from collections import Counter
import copy
import gzip
import hashlib
import io
import json
from pathlib import Path
import shlex
import subprocess
import sys
import tarfile
import tempfile

import reference as ref

sys.path.insert(0, str(ref.ROOT/'scripts'))
import S15_reading_coverage as candidate

def copy_archive(original, tmp, updates):
    # Keep unmodified members byte-for-byte. An update replaces one decoded blob.
    with tarfile.open(original) as old, tarfile.open(tmp,'w') as new:
        for info in old:
            data=old.extractfile(info).read()
            if info.name in updates:
                data=updates[info.name]
            info=copy.copy(info)
            info.size=len(data)
            new.addfile(info,io.BytesIO(data))

def gzjson(value):
    return gzip.compress(json.dumps(value).encode())

def summary(records, name):
    row,=[r for r in records if r['object']==name]
    return row, len(row['sampled']),len(row['untested'])

def statuses(items):
    return dict(Counter(r['status'] for r in items if 'velocity_change' in r['argument'][0]))

def change_rule(original, output):
    with tarfile.open(original) as t:
        family,=(row for row in ref.member(t,'context.json.gz') if row['field']=='family/centre_velocity_changes')
        records=ref.member(t,'coverage.json.gz')
        native,=(r for r in records if r['id']==family['right'][0])
        oldhash=native['payload']
        with gzip.open(t.extractfile('operands.jsonl.gz'),'rt') as source:
            buffer=io.BytesIO()
            hits=0
            with gzip.open(buffer,'wt') as dest:
                for line in source:
                    row=json.loads(line)
                    if row['sha256']==oldhash:
                        operand=row['native_operand']
                        changed=operand.replace('Derivative[1][xx[1, 2]][t]','Derivative[1][xx[2, 2]][t]',1)
                        assert changed!=operand
                        newhash=hashlib.sha256(changed.encode()).hexdigest()
                        row.update(native_operand=changed,sha256=newhash)
                        line=json.dumps(row)+'\n'
                        hits+=1
                    dest.write(line)
        assert hits==1
        native['payload']=newhash
    copy_archive(original,output,{'operands.jsonl.gz':buffer.getvalue(),'coverage.json.gz':gzjson(records)})

def report_replay(tmp):
    primary=ref.ROOT/'scripts/out/S15_reading_comparator.out'
    man_path=ref.ROOT/'directives/_measurements/S15_reading_comparator_run_manifest.json'
    archived=json.loads(man_path.read_text())
    args=shlex.split(archived['report_command'])
    assert args[-2]=='--publication'
    args[0:2]=[sys.executable,'-B',str(ref.ROOT/'scripts/S15_reading_report.py')]
    args[-1]=str(tmp/'replay')
    r=subprocess.run(args,capture_output=True,text=True,check=True)
    print('normal report stdout\n'+r.stdout.strip())
    rebuilt=tmp/'replay'
    audited=json.loads((rebuilt/'directives/_measurements'/man_path.name).read_text())
    covered=json.loads(ref.COV.read_text())
    lines=primary.read_text().splitlines()
    records=list(map(json.loads,lines))
    for report,expected in [('item3_exponent_coverage',covered['counts']),('observed_argument_checks',covered['argument_counts'])]:
        row,=(r for r in records if r['record']==report)
        values=row.get('counts',{k:v for k,v in row.items() if k!='record'})
        assert values==expected
        assert [line for line in lines if json.loads(line)['record']==report][0] in archived['transcript_summary_lines_verbatim']
        assert [line for line in lines if json.loads(line)['record']==report][0] in audited['transcript_summary_lines_verbatim']
    for identity,key in [(primary,'transcript'), (ref.ARCH,'detail_archive'), (ref.COV,'final_coverage')]:
        rebuilt_file=rebuilt/identity.relative_to(ref.ROOT)
        sha=ref.digest(identity)
        print('published and normal replay',key,'bytes',identity.stat().st_size,'sha256',sha)
        assert sha==ref.digest(rebuilt_file)==archived[key]['sha256']==audited[key]['sha256']
        if key!='final_coverage': assert archived[key]['bytes']==audited[key]['bytes']==identity.stat().st_size
    assert audited['final_coverage']['source_sha256']==archived['final_coverage']['source_sha256']==ref.digest(ref.ROOT/'scripts/S15_reading_coverage.py')
    assert archived['transcript_summary_lines_verbatim']==[x for x in lines if json.loads(x)['record']!='exception']
    assert audited['transcript_summary_lines_verbatim']==archived['transcript_summary_lines_verbatim']
    assert covered['counts']==archived['final_coverage']['exponent_coverage']==audited['final_coverage']['exponent_coverage']
    assert covered['argument_counts']==archived['final_coverage']['argument_counts']==audited['final_coverage']['argument_counts']
    print('normal report JSON transcript rows',len(records),'verbatim manifest summary rows',len(archived['transcript_summary_lines_verbatim']),'both revised rows present once')

def main():
    mode='remaining' if sys.argv[1:]==['--remaining'] else 'all'
    assert mode in ('all','remaining')
    with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as tmpname:
        tmp=Path(tmpname)
        if mode=='all': report_replay(tmp)
        with tarfile.open(ref.ARCH) as t:
            cfg=ref.member(t,'configuration.json')
            src=(ref.ROOT/'mathematica/S15_separation_departure.wl').read_text()
            source_info=dict(path='research/pde_ledger_v3/mathematica/S15_separation_departure.wl',revision='cc62dc9f',sha256=hashlib.sha256(src.encode()).hexdigest())
            if mode=='all': velocity=candidate.argument_slots(t,src,source_info)
            name='inner/carried/component_1/labels_drain1+drain2'
            k=cfg['required_objects'].index(name)
            file=f'objects/{k:03d}.json.gz'
            obj=ref.member(t,file)
            before,n,m=summary(candidate.exponent_coverage(t,cfg),name)
            if mode=='all':
                print('positive controls',name,(n,m),'velocity',statuses(velocity))
                changed=src.replace('Do[register[centreVelocityChange[j],lenDim-timeDim,{timeDim}]','Do[register[centreVelocityChange[j],lenDim-timeDim,{lenDim}]',1)
                assert changed!=src
                path=tmp/'Wolfram-length-argument.wl'
                path.write_text(changed)
                bad=candidate.argument_slots(t,path.read_text(),source_info|{'sha256':ref.digest(path)})
                print('one-sided native registration timeDim to lenDim',[(r['wl_head'],r['observed'],r['wl_expected'],r['status']) for r in bad if 'velocity_change' in r['argument'][0]])
                assert statuses(bad)=={'slot_dimension_difference':3}
        py,pr=ref.native('py');wl,wr=ref.native('wl')
        native=tuple(sorted(py[name][0]))
        feasible=ref.arrangement(native)
        assert all(row in wl[name][0] for row in native)
        with tarfile.open(ref.ARCH) as t:
            ids,points=ref.observation_set(t,file,obj,Counter())
        if mode=='all':
            index=next(i for i,p in enumerate(obj['exponent_loci']) if p[-1]!='0')
            mutant=copy.deepcopy(obj)
            original=mutant['exponent_loci'][index][:]
            mutant['exponent_loci'][index][-1]=str(-ref.F(original[-1]))
            af=tmp/'plane-sign.tar'
            copy_archive(ref.ARCH,af,{file:gzjson(mutant)})
            with tarfile.open(af) as t:
                after,_,_=summary(candidate.exponent_coverage(t,cfg),name)
            try: ref.compare_sets(after,native,feasible,points)
            except AssertionError as exc:
                assert exc.args==('native plane mismatch',)
                print('one-sided checkpoint plane sign',original,'->',mutant['exponent_loci'][index],'native object comparison:',str(exc))
            else: raise AssertionError('plane sign control passed')
        changed=copy.deepcopy(obj)
        row=next(r for r in changed['rows'] if r['field']==ref.FIELD)
        s=row['samples'][0]
        assert s['case'][3]=='generic'
        point=s['exponents'][:]
        s['exponents'][0]='-3';s['exponents'][2]='-2'
        af=tmp/'sample.tar'
        copy_archive(ref.ARCH,af,{file:gzjson(changed)})
        with tarfile.open(af) as t:
            updated,mc,nc=summary(candidate.exponent_coverage(t,cfg),name)
            ids,newpoints=ref.observation_set(t,file,changed,Counter())
        # JSON converts integer codimension keys to strings during publication.
        ref.compare_sets(json.loads(json.dumps(updated)),native,feasible,newpoints)
        newsel=set(map(tuple,updated['sampled']))-set(map(tuple,before['sampled']))
        assert len(newsel)==1 and len(next(iter(newsel)))==3 and (mc,nc)==(n+1,m-1)
        print('one archived generic case only',s['case'],point,'->',s['exponents'],'added simultaneous plane indices',sorted(newsel),'sampled/untested',(n,m),'->',(mc,nc))
        af=tmp/'rule-slot.tar'
        change_rule(ref.ARCH,af)
        with tarfile.open(af) as t:
            try: candidate.argument_slots(t,src,source_info)
            except ValueError as exc:
                assert str(exc)=='source centre or velocity slot differs from declared binding'
                print('one-sided Wolfram derivative centre [1,2] -> [2,2], all right-hand responses retained:',str(exc))
            else: raise AssertionError('velocity slot mutation passed')
        # Independent non-affine mutation on the native condition syntax alone.
        try: ref.affine('a_1**2 + b_1','-1')
        except AssertionError:
            print('native condition form a_1 + b_1 == -1 -> a_1**2 + b_1 == -1: rejected as non-affine')
        else: raise AssertionError('non-affine control passed')

if __name__=='__main__': main()
