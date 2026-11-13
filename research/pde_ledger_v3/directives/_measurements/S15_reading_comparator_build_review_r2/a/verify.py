"""Compare candidate publication to the independent reference on the original packet."""
import gzip
import hashlib
import io
import json
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from collections import Counter
from pathlib import Path
from reference import ARCHIVE, COVERAGE, HERE, MANIFEST, P, PACKET, SUMMARY, member

sys.path.insert(0,str(PACKET/P/'scripts'))
import S15_reading_coverage as candidate


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()


def normalize(row):
    """Translate ordered index lists to native plane coefficients as identities."""
    planes = [tuple(s) for s in row['planes']]
    assert len(planes)==len(set(planes))
    result = {}
    for name in ('sampled','untested'):
        groups = []
        for indices in row[name]:
            assert all(type(i) is int and 0<=i<len(planes) for i in indices)
            assert list(indices)==sorted(set(indices))
            groups.append(frozenset(planes[i] for i in indices))
        assert len(groups)==len(set(groups))
        result[name] = set(groups)
    assert not result['sampled'] & result['untested']
    return set(planes), result


def write_copy(original, out, replacements):
    with tarfile.open(original) as src, tarfile.open(out,'w') as dest:
        for info in src:
            if info.name in replacements:
                data = replacements[info.name]
                info.size = len(data)
                dest.addfile(info,io.BytesIO(data))
            else:
                dest.addfile(info,src.extractfile(info))
    return out


def enc(obj):
    return gzip.compress(json.dumps(obj,separators=(',',':')).encode(),compresslevel=1)


def entry(archive, index):
    with tarfile.open(archive) as ar:
        cfg=member(ar,'configuration.json')
        rows=candidate.exponent_coverage(ar,cfg)
    return next(r for r in rows if r['archive_index']==index)


def main():
    ref=json.loads((HERE/'reference.json').read_text())
    baseline=json.loads(COVERAGE.read_text())
    manifest=json.loads(MANIFEST.read_text())
    entries=baseline['objects']
    assert [r['object'] for r in entries]==[r['object'] for r in ref['objects']]
    for r,gold in zip(entries,ref['objects']):
        assert r['archive_index']==gold['index'] and r['operand_ids']==gold['operand_ids']
        assert normalize(r)==normalize(gold) and r['untested_by_codimension']==gold['ranks'], r['object']
    print('per-object plane sets, IDs, every sampled/untested selection and gap codimension:',len(entries),'equal')
    with tarfile.open(ARCHIVE) as ar:
        archived=member(ar,'context_checks.json.gz')['argument_slots']
        cov={r['id']:r for r in member(ar,'coverage.json.gz')}
        def operand_keys(row):
            return [(i.split(':')[0],cov[i]['payload']) for i in row['left']+row['right']]
        alias_count=0
        for ix in range(336):
            obj=member(ar,f'objects/{ix:03}.json.gz')
            r,=[r for r in obj['rows'] if r['field']=='item3/object']
            while 'same_comparison_as' in r:
                old=r
                r=obj['rows'][r['same_comparison_as']]
                assert operand_keys(old)==operand_keys(r)
                alias_count+=1
        print('item3/object alias edges with identical per-engine payload digests:',alias_count)
        mixed=member(ar,'objects/288.json.gz')
        base_obj=member(ar,'objects/010.json.gz')
    originals={tuple(r['argument']):r for r in archived}
    published={tuple(r['argument']):r for r in baseline['argument_slots']}
    assert len(originals)==len(archived)==len(published)==len(baseline['argument_slots']) and originals.keys()==published.keys()
    slot_refs=ref['slots']
    corrected=[]
    for key,row in published.items():
        if key[0] in slot_refs:
            assert row==slot_refs[key[0]]
            corrected.append(key)
        else:
            assert row==originals[key]
    assert len(corrected)==3
    counts=dict(Counter(s['status'] for s in published.values()))
    assert baseline['argument_counts']==counts
    print('source-derived velocity entries equal:',len(corrected),'unchanged other slot entries:',len(published)-len(corrected),counts)
    lines=SUMMARY.read_text().splitlines()
    by_kind={}
    for line in lines:
        item=json.loads(line)
        by_kind.setdefault(item['record'],[]).append(item)
    tc,=by_kind['item3_exponent_coverage']
    sc,=by_kind['observed_argument_checks']
    assert {k:v for k,v in tc.items() if k!='record'}==baseline['counts']==manifest['final_coverage']['exponent_coverage']
    assert sc['counts']==counts==manifest['final_coverage']['argument_counts']
    assert sc['python_empty']==[s['argument'] for s in published.values() if not s['observed']]
    assert manifest['transcript_summary_lines_verbatim']==[l for l in lines if json.loads(l)['record']!='exception']
    exceptions={l for l in lines if json.loads(l)['record']=='exception'}
    quoted=manifest['representative_exceptions_verbatim']
    assert len(quoted)==4 and len(exceptions)==8 and set(quoted)<=exceptions
    for k,path in [('final_coverage',COVERAGE),('detail_archive',ARCHIVE),('transcript',SUMMARY)]:
        meta=manifest[k]
        assert sha(path)==meta['sha256']==sha(Path(meta['path']))
        if k!='final_coverage':
            assert path.stat().st_size == meta['bytes']
    assert baseline['archive_sha256']==sha(ARCHIVE)
    assert manifest['source_sha256']['S15_reading_report.py']==sha(PACKET/P/'scripts/S15_reading_report.py')
    assert manifest['final_coverage']['source_sha256']==sha(PACKET/P/'scripts/S15_reading_coverage.py')
    for loc in [P/'mathematica/S15_separation_departure.wl',P/'scripts/S15_separation_departure.py']:
        committed=subprocess.check_output(['git','show','cc62dc9f:'+str(loc)])
        assert committed==(PACKET/loc).read_bytes()
    print('normal transcript and manifest:',len(lines),'lines with all summary lines and four representative exceptions verbatim; 3 publication file hashes and 2 publisher code hashes equal')
    source=(PACKET/P/'mathematica/S15_separation_departure.wl').read_text()
    sinfo=ref['slots']['X_1_1_velocity_change']['wl_source']
    sinfo={k:sinfo[k] for k in ('path','revision','sha256')}
    with tempfile.TemporaryDirectory() as temp:
        folder=Path(temp)
        produced=folder/'coverage.json'
        candidate.publish(manifest,ARCHIVE,PACKET,produced)
        assert COVERAGE.read_bytes()==produced.read_bytes()
        print('fresh normal publisher coverage sha256:',sha(produced),'byte identical')
        # Change one archived plane constant for one object. All samples and other
        # operands, including the velocity slot records, stay live in the archive.
        key='objects/010.json.gz'
        normal=next(q for q in entries if q['archive_index']==10)
        old=base_obj['exponent_loci'][1].copy()
        base_obj['exponent_loci'][1][-1]='-1'
        assert old!=base_obj['exponent_loci'][1]
        changed=entry(write_copy(ARCHIVE,folder/'sign.tar',{key:enc(base_obj)}),10)
        pn,sn=normalize(normal)
        pm,sm=normalize(changed)
        assert pn-pm=={tuple(old)} and len(pm-pn)==1
        assert len(sn['sampled'])==len(sm['sampled'])+1 and len(sm['untested'])==len(sn['untested'])+1
        print('one-sided input plane constant',old,'->',base_obj['exponent_loci'][1], 'sampled/untested',len(sn['sampled']),len(sn['untested']),'->', len(sm['sampled']),len(sm['untested']))
        # Change only one actual sample tuple, leaving its case label unmodified.
        row,=[r for r in mixed['rows'] if r['field']=='item3/object']
        sample=row['samples'][0]
        old=sample['exponents']
        sample['exponents']=['-1','-1','-2','-2']
        normal=next(q for q in entries if q['archive_index']==288)
        changed=entry(write_copy(ARCHIVE,folder/'sample.tar',{'objects/288.json.gz':enc(mixed)}),288)
        added=set(map(tuple,changed['sampled']))-set(map(tuple,normal['sampled']))
        assert added=={(0,2,4,6)}
        assert set(map(tuple,normal['untested']))-set(map(tuple,changed['untested']))==added
        print('one-sided tuple',old,'->',sample['exponents'],'unchanged',sample['case'],'newly witnessed',sorted(added),'sampled/untested',len(changed['sampled']),len(changed['untested']))
        # Change the form of the source declaration: the argument takes length.
        modified=source.replace('centreVelocityChange[j],lenDim-timeDim,{timeDim}','centreVelocityChange[j],lenDim-timeDim,{lenDim}')
        assert source!=modified
        sourcepath=folder/'S15_separation_departure.wl'
        sourcepath.write_text(modified)
        with tarfile.open(ARCHIVE) as ar:
            slots=candidate.argument_slots(ar,sourcepath.read_text(),{**sinfo,'sha256':sha(sourcepath)})
        diffs=[(s['wl_head'],s['wl_expected'],s['status']) for s in slots if s['argument'][0] in slot_refs]
        assert all(d==[0,1,0] and st=='slot_dimension_difference' for h,d,st in diffs) and len(diffs)==3
        print('one-sided native input dimension timeDim -> lenDim:',diffs)
        # Change only the SymPy observation of component 2, leaving source and
        # native rule positions intact. Compare the complete 406-row results.
        row,=[r for r in archived if r['argument']==['X_1_2_velocity_change',0,'t']]
        row['observed']=[[0,-1,1]]
        with tarfile.open(ARCHIVE) as ar:
            ctx=member(ar,'context_checks.json.gz')
        ctx['argument_slots']=archived
        m=write_copy(ARCHIVE,folder/'slot.tar',{'context_checks.json.gz':enc(ctx)})
        with tarfile.open(m) as ar:
            slots=candidate.argument_slots(ar,source,sinfo)
        changed=[s for s in slots if s!=published[tuple(s['argument'])]]
        assert len(changed)==1 and changed[0]['wl_head']=='centreVelocityChange[2]' and changed[0]['status']=='slot_dimension_difference'
        print('one-sided observed component 2 [0,0,1] -> [0,-1,1]: only',changed[0]['wl_head'],changed[0]['status'])
        # Test an actual Wolfram rule rebind in a content-addressed temporary
        # catalogue with the new payload hash referenced by the same occurrence.
        with tarfile.open(ARCHIVE) as ar:
            records=member(ar,'coverage.json.gz')
            rule,=[r for r in records if r['id']==slot_refs['X_1_1_velocity_change']['wl_velocity_rules']]
            old=rule['payload']
            compressed=folder/'operands.jsonl.gz'
            with gzip.open(ar.extractfile('operands.jsonl.gz'),'rb') as fr, gzip.open(compressed,'wb',compresslevel=1) as fw:
                for line in fr:
                    if old.encode() in line:
                        value=json.loads(line)
                        assert value['sha256']==old
                        native=value['native_operand']
                        replacement=native.replace('xx[1, 2]','xx[2, 2]')
                        assert replacement!=native
                        value['native_operand']=replacement
                        value['sha256']=hashlib.sha256(replacement.encode()).hexdigest()
                        rule['payload']=value['sha256']
                        line=(json.dumps(value,separators=(',',':'))+'\n').encode()
                    fw.write(line)
            assert old!=rule['payload']
        m=write_copy(ARCHIVE,folder/'rule.tar',{'coverage.json.gz':enc(records),'operands.jsonl.gz':compressed.read_bytes()})
        with tarfile.open(m) as ar:
            try:
                candidate.argument_slots(ar,source,sinfo)
            except ValueError as error:
                assert str(error)=='source centre or velocity slot differs from declared binding'
                print('one-sided native rule xx[1, 2] -> xx[2, 2]:',str(error))
            else:
                raise AssertionError('Wrong source centre accepted')
    print('all scoped publication and mutation comparisons passed')

if __name__=='__main__':
    main()
