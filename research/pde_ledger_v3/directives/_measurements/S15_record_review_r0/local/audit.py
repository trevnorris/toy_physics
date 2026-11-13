from collections import Counter, defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tarfile

root = Path.cwd()
local = Path(__file__).resolve().parent
base = Path("research/pde_ledger_v3")
frozen = local.parent / "candidate"
evidence = (frozen / base / "steps/_measurements/S15_record_measurements.md").read_text()
m3, m4 = evidence.split("## M3.")[1].split("## M4.")
m4 = m4.split("## M5.")[0]
manifest = json.loads(subprocess.check_output(["git", "show", f"41adb784:{base}/directives/_measurements/S15_separation_departure_wolfram_run_manifest.json"]))
eligible = {n for n,s in manifest['item4_statuses'].items() if 'Not applicable' not in s}
fd = hashlib.sha256()
latest = defaultdict(list)
completed = {}
ap = b"/item4/"
with open(manifest['transcript'],'rb') as f:
    for line in f:
        fd.update(line)
        if line.startswith(b"END_OF_ITEM4_PASS:"):
            obj = json.loads(line.split(b": ", 1)[1])
            completed[obj] = latest.pop(obj)
        elif line.startswith(b'WORKING_ON: "reading/') and line.endswith(b' / object"\n'):
            obj = line[len(b'WORKING_ON: "reading/'):-len(b' / object"\n')].decode()
            latest.pop(obj,None)
        elif line.startswith(b"END_OF_BASE_PASS:"):
            completed.update({k:v for k,v in latest.items() if v[-1][1] == manifest['item4_statuses'][k]})
            latest.clear()
        elif line.startswith(b"reading/") and ap in line[:170]:
            name, val = line.rstrip(b'\n').split(b': ',1)
            obj, tag = (x.decode() for x in name[len(b'reading/'):].split(ap))
            if obj not in eligible or not re.fullmatch(r'(?:condition|Q_under_condition)(?:_\d+)?|status|unresolved_parameter_branches',tag):
                continue
            latest[obj].append((tag,val.decode()))
assert fd.hexdigest() == manifest['transcript_sha256']
assert set(completed) == eligible
sections = re.split(r"^### (.+)$",m4,flags=re.M)[1:]
published = {n:[tuple(s.split(': ',1)) for s in section.splitlines() if s and not s.startswith('```')] for n,section in zip(sections[::2],sections[1::2])}
assert set(published) == eligible
count = Counter()
for obj, rows in published.items():
    assert Counter(rows) == Counter(completed[obj]), obj
    for k,v in rows:
        count[re.sub(r'_\d+$','',k)]+=1
    for k in {k for k, _ in rows}:
        assert [v for key,v in rows if key==k] == [v for key,v in completed[obj] if key==k], (obj,k)
print("M4",len(published),"eligible native objects", dict(count))
print("M4 Wolfram full sha256",fd.hexdigest())
sections = re.split(r"^### (\d+)\. (.+)$",m3,flags=re.M)[1:]
m3_rows = {int(i):(name,rows) for i,name,rows in zip(sections[::3],sections[1::3],sections[2::3])}
counts=Counter()
corrections=Counter()
raw_counts=Counter()
baseline='research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md'
with tarfile.open(base/'scripts/out/S15_reading_comparator_details.out') as ar, (base/'scripts/out/S15_separation_departure.out').open('rb') as py:
    def read(n):return json.loads(gzip.decompress(ar.extractfile(n).read()))
    locations={r['id']:r for r in read('coverage.json.gz')}
    for i,(name,report) in m3_rows.items():
        src = read(f'objects/{i:03d}.json.gz')
        assert src['object']==name
        data={(r['role'],r['field']):r for r in src['rows']}
        status = {}
        refpath=f'refinement/{i:03d}.json.gz'
        if refpath in ar.getnames():
            ref=read(refpath)
            assert ref['object']==name
            if ref['branch']:
                assert ref['branch']['same_elementary_boundaries']
                status['item3/exponent_branch_inventory']='same_elementary_loci'
            logpath=f'log_branches/{i:03d}.json.gz'
            log=read(logpath) if logpath in ar.getnames() else {'object':name,'rechecks':[]}
            assert log['object']==name
            replacements={(r['field'],tuple(r['case']),tuple(r['exponents'])):r for r in log['rechecks']}
            for fld in ref['fields']:
                k=fld['field']
                values=Counter()
                for p in fld['additional_samples']:
                    p=replacements.pop((k, tuple(p['case']),tuple(p['exponents'])),p)
                    assert p['status']=='within_estimate'
                    values[p['status']]+=1
                    corrections['branch_checks']+=1
                for p in fld['rechecks']:
                    assert p['status']=='within_precision_estimate'
                    values[p['status']]+=1
                    corrections['precision_checks']+=1
                assert f"{k} {dict(values)}" in report
                initial=data['cross',k]['status']
                status[k]=('sample_agreement' if initial not in ('exact_identity','sample_coverage_incomplete') and values else initial)
            assert not replacements
        matched=re.findall(r'^- `(\w+) / ([^`]+)` (\w+) → (\w+)\. (.+) ↔ (.+)$',report,re.M)
        for role,k,initial,final,left,right in matched:
            source=data[role,k]
            assert source['status']==initial
            assert (','.join(left.split(', ')),','.join(right.split(', ')))==(','.join(source['left']),','.join(source['right']))
            assert (status.get(k,initial) if role=='cross' else initial)==final, (i,k)
            counts[(role,k)]+=1
        expected_fields=['Q','E_geo','E_resp','separation_derivative','departure_partition/residual','item3/family_integral','item3/object','item3/Q_closed_operand','item3/exponent_branch_inventory','velocity_part_1/residual','velocity_part_2/residual']
        assert {k for role,k,*_ in matched if role=='cross'}=={k for k in expected_fields if ('cross',k) in data}
        assert {role for role,*_ in matched}=={'cross','py_same_engine','wl_same_engine'}
        wolfram=locations[data['cross','Q']['right'][0]]['tag'].removeprefix('reading/').removesuffix('/Q')
        assert f"Single-engine item-4 status for `{wolfram}` {manifest['item4_statuses'][wolfram]}" in report
        for k in (1,2):
            loc=locations[data['cross',f'velocity_part_{k}/residual']['left'][0]]
            span=loc['span'];py.seek(span['offset']);raw=py.read(span['length']).removesuffix(b'\n')
            assert hashlib.sha256(raw).hexdigest()==loc['payload']
            if i < 3:
                print(name,f"drain{k}",loc['id'], '0' if raw==b'0' else f"{len(raw)}-byte expression")
            if name.startswith(('both_surfaces/','inner/','outer/')) and raw==b'0':
                raw_counts[f'literal-zero-surface-drain{k}']+=1
print("M3",len(m3_rows),"objects",sum(counts.values()),"role/field/operand rows",dict(corrections),dict(raw_counts))
def census(s):
    parts=re.split(r'^### (R-\S+).*$',s.split('## Population passes')[0],flags=re.M)[1:]
    return {k:re.search(r'\*\*status\*\* (\w+)',v).group(1) for k,v in zip(parts[::2],parts[1::2])}
prior=census(subprocess.check_output(['git','show',f'83237f4a:{baseline}'],text=True))
now=census((frozen/baseline).read_text())
assert now==prior
print("register",dict(Counter(prior.values())),"->",dict(Counter(now.values())))
