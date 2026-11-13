"""Reconcile all initial, precision and extra branch sample counts and locators."""
from collections import Counter,defaultdict
from packet_io import *
d=Records()
m=json.loads((S.parent/'directives/_measurements/S15_reading_comparator_run_manifest.json').read_text())
counts=Counter();record_types=Counter();parents=set();strata=Counter();absent=[];sampled=defaultdict(set)
for oidx,name in enumerate(d.cfg['required_objects']):
    o=d.obj(oidx)
    for j,r in enumerate(o['rows']):
        if r['role']!='cross':continue
        rr=o['rows'][r.get('same_comparison_as',j)]
        if 'samples' not in rr:continue
        cases=[tuple(s['case']) for s in rr['samples']]
        expected={(f,t,p,b) for f in d.cfg['forms'] for t in d.cfg['strata'] for p in range(3)
                  for b in (['generic']+['locus_'+str(i) for i in range(len(o['exponent_loci']))] if r['field'] in ('item3/family_integral','item3/object','item3/Q_closed_operand') else ['generic'])}
        assert set(cases)==expected and len(cases)==len(expected)
        if any(s['status']=='separated' for s in rr['samples']):sampled[(name,r['field'])]={tuple(s['case']) for s in rr['samples'] if s['status']=='separated'}
        parents.add((oidx,r.get('same_comparison_as',j)))
    refpath=f'refinement/{oidx:03d}.json.gz'
    if refpath not in d.tar.getnames():continue
    ref=member(d.tar,refpath)
    logpath=f'log_branches/{oidx:03d}.json.gz'
    logs=member(d.tar,logpath)['rechecks'] if logpath in d.tar.getnames() else []
    replace={(r['field'],tuple(r['case'])):r for r in logs}
    for field in ref['fields']:
        nm=field['field']
        for row in field['rechecks']:
            if (name,nm) in sampled:sampled[(name,nm)].remove(tuple(row['case']))
            counts['high precision',row['status']]+=1
        if not field['additional_samples']:continue
        b=ref['branch']
        labels=['single_'+str(i) for i in range(len(b['single_planes']))]+['joint_'+str(i) for i in range(len(b['joint_planes']))]
        expected={(f,t,p,l) for f in d.cfg['forms'] for t in d.cfg['strata'] for p in range(3) for l in labels}
        actual=set()
        for s0 in field['additional_samples']:
            case=tuple(s0['case']);actual.add(case)
            rr=replace.pop((nm,case),None)
            if rr is not None:
                assert rr['exponents']==s0['exponents'] and s0['status']=='unresolved'
                record_types['log replacements']+=1
                s0=rr
            counts['additional',s0['status']]+=1
            if 'values' in s0:assert s0['values'][-1][0]-s0['values'][-1][1]==s0['residual']
        assert actual==expected and len(actual)==len(field['additional_samples'])
    assert not replace
print('unique primary arithmetic sample rows',len(parents),'precision tallies',dict(counts),'source of logged rechecks',dict(record_types),'unresolved early separated',sum(map(len,sampled.values())))
import hashlib,subprocess
for lang,prefix in [('scripts','py'),('mathematica','wl')]:
    ext='py' if prefix=='py' else 'wl'
    for part,revision in [('two_drain_region','3192123b'),('separation_departure','cc62dc9f')]:
        rel=f'research/pde_ledger_v3/{lang}/S15_{part}.{ext}'
        dig=sha((P/rel).read_bytes());git=sha(subprocess.check_output(['git','show',revision+':'+rel]))
        print('native source',rel,revision,dig,'match',git==dig);assert dig==git
print('supplemental source declaration')
for ln in (P/'research/pde_ledger_v3/mathematica/S15_separation_departure.wl').read_text().splitlines()[224:229]:print(ln)
print('all_sources_match_manifest',all(sha((S/name).read_bytes())==h for name,h in m['source_sha256'].items()))
types=Counter();missing=[];used=set()
for i in range(len(d.cfg['required_objects'])):
    o=d.obj(i)
    for r in o['rows']:
        ids=r['left']+r['right'];used.update(ids)
        if r['status']=='unpaired':types[r['reason']]+=1
for r in d.idx.values():
    if r['role']=='scoped' and r['id'] not in used:missing.append(r['id'])
print('unpaired reasons',dict(types),'scoped missing ids',missing)
