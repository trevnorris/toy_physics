"""Validate every saved comparison, then publish the transcript and guarded run manifest."""
from argparse import ArgumentParser
from collections import Counter, defaultdict
import gzip
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import sys
import tarfile

REPO=next(p for p in Path(__file__).resolve().parents if (p/'research/pde_ledger_v3').is_dir() and (p/'.git').exists())
LEDGER=REPO/'research/pde_ledger_v3'
MEASURE=LEDGER/'directives/_measurements'


def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()


def read(path):
    open_file=gzip.open if path.suffix=='.gz' else open
    with open_file(path,'rt') as f:
        return json.load(f)


def relative(path):
    return str(path.resolve().relative_to(REPO))


def guarded(path):
    inv=read(path/'invocation.json')
    validation=read(path/'limit-validation.json')
    actual=validation['actual']
    out=read(path/'outcome.json')
    if not (out['limitsVerified'] and validation['verified'] and out['exitCode']==out['childOutcome']['exitCode']==0
            and out['childOutcome']['guardReason'] is None and actual['memory.swap.max']=='0'
            and inv['wallDeadlineSeconds'] is None and inv['inactivityDeadlineSeconds'] is None):
        raise ValueError(f'Guard failed: {path}')
    with (path/'resource-samples.jsonl').open() as f:
        samples=list(map(json.loads,f))
    return dict(directory=relative(path),unit=inv['unit'],actual_limits=actual,
                guard_child_command=inv['command'],pool=inv['pool'],started_utc=inv['startedUtc'],
                last_resource_sample_utc=samples[-1]['utc'],guard_runtime_seconds=out['wallSeconds'],
                exit_code=out['exitCode'],child_exit_code=out['childOutcome']['exitCode'],
                sampled_memory_peak_bytes=max(int(r['memory.peak']) for r in samples),
                guard_outcome=out,guard_receipt_sha256=digest(path/'outcome.json'))


def main():
    cli=ArgumentParser(description=__doc__)
    cli.add_argument('root',type=Path)
    cli.add_argument('--review-controls',type=Path,help='guard for the one-sided native mutations')
    cli.add_argument('--unit-controls',type=Path,help='guard for the unit controls')
    cli.add_argument('--publication',type=Path,help='optional alternate root for candidate outputs')
    cli.add_argument('--guard-parent',type=Path,required=True)
    cli.add_argument('--context-run',required=True)
    opt=cli.parse_args()
    root=opt.root
    config=read(root/'configuration.json')
    summary=read(root/'summary.json')
    coverage=read(root/'coverage.json.gz')
    context=read(root/'context.json.gz')
    aux=read(root/'context_checks.json.gz')
    required=config['required_objects']
    if len(required)!=336 or summary['objects']!=len(required) or any(any(e.values()) for e in config['roster'].values()):
        raise ValueError('roster or checkpoint count mismatch')
    by_id={x['id']:x for x in coverage}
    if len(by_id)!=len(coverage):raise ValueError('duplicate occurrence ID')
    # Verify content-addressed operands, even when a later record references an earlier comparison.
    needed={x['payload'] for x in coverage if 'payload' in x}
    with gzip.open(root/'operands.jsonl.gz','rt') as f:
        for line in f:
            entry=json.loads(line)
            if hashlib.sha256(entry['native_operand'].encode()).hexdigest()!=entry['sha256']:
                raise ValueError('catalogue blob mismatch')
            needed.remove(entry['sha256'])
    if needed:raise ValueError('missing native operands')
    lines=[]
    def report(kind,**data):
        lines.append(json.dumps({'record':kind,**data},ensure_ascii=False,sort_keys=True,separators=(',',':')))
    counts=defaultdict(Counter)
    by_field=defaultdict(Counter)
    links=defaultdict(Counter)
    family_forms=defaultdict(Counter)
    zero_spellings=defaultdict(Counter)
    literal_zero=hashlib.sha256(b'0').hexdigest()
    widths=defaultdict(Counter)
    examples=defaultdict(list)
    numsamples=Counter()
    sampled_residuals=defaultdict(Counter)
    exceptions=[]
    paired=set()
    initial_status={}
    initial_counts={}
    original_flags={}
    native_forms={}
    for k,name in enumerate(required):
        obj=read(root/'objects'/f'{k:03d}.json.gz')
        if (obj['index'],obj['object'],obj['build_id'])!=(k,name,config['build_id']):
            raise ValueError('wrong saved object')
        byrolefield={}
        for index,row in enumerate(obj['rows']):
            role=row['role'];field=row['field'];st=row['status']
            key=role+' / '+field
            counts[str(row['item'])][st]+=1
            by_field[key][st]+=1
            refs=row['left']+row['right']
            if any(r not in by_id or 'payload' not in by_id[r] for r in refs):raise ValueError('bad native reference')
            paired.update(refs)
            origin=row
            if 'same_comparison_as' in row:
                ancestor=row['same_comparison_as']
                if not 0<=ancestor<index:raise ValueError('forward comparison reference')
                origin=obj['rows'][ancestor]
                if st!=origin['status']:raise ValueError('reference status mismatch')
                lhs=[(r[:2],by_id[r]['payload']) for r in refs]
                rhs=[(r[:2],by_id[r]['payload']) for r in origin['left']+origin['right']]
                if lhs!=rhs:raise ValueError('operand changed in reference')
            byrolefield[key]=st
            if role=='cross':
                initial_status[name,field]=st
                if 'sample_counts' in origin:
                    initial_counts[name,field]=origin['sample_counts']
                    if st=='sample_difference':
                        original_flags[name,field]={tuple(v['case']) for v in origin['samples'] if v['status']=='separated'}
            if field=='item3/object' and role=='cross':
                native_form=('closed' if any(r['field']=='item3/Q_closed_operand' for r in obj['rows']) else 'integral')
                family_forms[native_form][st]+=1
                native_forms[name]=native_form
            if field in ('E_geo','E_resp') and len(refs)==2:
                az,bz=[by_id[r]['payload']==literal_zero for r in refs]
                if az!=bz:zero_spellings[field][st]+=1
            if role.endswith('_same_engine'):
                links[role][('direct_' if len(row['right'])==1 else 'sum_')+st]+=1
            if 'samples' in origin:
                sc=Counter(s['status'] for s in origin['samples'])
                if sc!=origin['sample_counts']:raise ValueError('sample status miscount')
                if origin is row:numsamples.update(sc)
                sampled_residuals[key].update(sc)
                present={tuple(s['case']) for s in origin['samples']}
                expected={(f,t,p,b) for f in config['forms'] for t in config['strata']
                          for p in range(len(config['points'])) for b in
                          (['generic']+[f'locus_{l}' for l in range(len(obj['exponent_loci']))]
                           if field in ('item3/object','item3/Q_closed_operand','item3/family_integral') else ['generic'])}
                # Identical native operands can refer to an earlier row; branch sampling is part of its identity.
                if present!=expected or len(present)!=len(origin['samples']):raise ValueError('case grid mismatch '+key)
                for sample in origin['samples']:
                    if sample['status']!='unresolved':
                        a,b=sample['values'][-1]
                        if a-b!=sample['residual']:raise ValueError('stored sample residual mismatch')
            if st in ('sample_difference','exact_difference','dimension_difference','unresolved','sample_coverage_incomplete',
                      'different_exponent_loci'):
                info=dict(object=name,role=role,field=field,status=st)
                if 'residual' in origin and 'samples' not in origin:info['residual']=origin['residual']
                if 'reason' in origin:info['reason']=origin['reason']
                if 'samples' in origin:
                    info['sample_counts']=origin['sample_counts']
                    bad=[q for q in origin['samples'] if q['status']!='within_estimate']
                    if bad:info['first_case']=bad[0]
                exceptions.append(info)
                if len(examples[key])<2:examples[key].append(info)
            if field=='item3/exponent_branch_inventory':widths['exponent_loci'][len(obj['exponent_loci'])]+=1
        for field in ('Q','E_geo','E_resp','item3/object','velocity_part_1/residual','velocity_part_2/residual'):
            if 'cross / '+field not in byrolefield:raise ValueError('missing main outcome')
    missing_refs=[r['id'] for r in coverage if r['role']=='scoped' and r['id'] not in paired]
    if missing_refs:raise ValueError('unaccounted per-object native fields: '+str(missing_refs[:12]))
    if counts!=summary['by_item'] or by_field!=summary['by_field']:raise ValueError('independent tally mismatch')
    report('roster',objects=len(required),gaps=config['roster'],counts=config['occurrences'])
    for stream in config['occurrences']:
        report('native_roles',stream=stream,counts=Counter(x['role'] for x in coverage if x['id'].startswith(stream+':')))
    for role,v in sorted(links.items()):report('reading_to_original',role=role,counts=v)
    for k,v in sorted(by_field.items()):report('field',field=k,counts=v)
    for k,v in sorted(family_forms.items()):report('family_result',form=k,counts=v)
    for k,v in sorted(zero_spellings.items()):report('different_literal_zero_spelling',field=k,counts=v)
    report('separate_numerical_cases',counts=numsamples)
    report('exponent_boundaries_per_object',counts=widths['exponent_loci'])
    for key,v in sorted(sampled_residuals.items()):report('sample_cases_per_field',field=key,counts=v)
    for row in context:
        if 'normalized_py' in row or row['field'] in ('family/choice','profile_family/chart_note'):
            report('family_context',field=row['field'],status=row['status'])
    report('kernel_checks',counts=Counter(q['status'] for q in aux['kernels']))
    ds=aux['wolfram_divisor_loci']
    report('wolfram_native_divisor_audit',records=len(ds),supports=sum(v['branch_supports'] for v in ds),
           nonfalse=sum(v['other_returns'] for v in ds))
    refined=defaultdict(dict)
    log_replacements={}
    for path in sorted((root/'log_branches').glob('*.gz')):
        obj=read(path)
        for row in obj['rechecks']:
            key=(obj['object'],row['field'],tuple(row['case']))
            if key in log_replacements:raise ValueError('duplicate log recheck')
            log_replacements[key]=row
    targeted=Counter()
    precision_examples={}
    added=Counter()
    joint_counts=Counter()
    for path in sorted((root/'refinement').glob('*.gz')):
        obj=read(path)
        if obj['build_id']!=config['build_id']:raise ValueError('refinement input mismatch')
        name=obj['object']
        branch=obj['branch']
        if branch:
            joint_counts[(len(branch['single_planes']),len(branch['joint_planes']))]+=1
            refined[name]['item3/exponent_branch_inventory']=('same_elementary_loci' if branch['same_elementary_boundaries']
                                                              else 'different_exponent_loci')
        for field in obj['fields']:
            tag=field['field']
            initial=initial_status[name,tag]
            evidence=[]
            for sample in field['additional_samples']:
                key=(name,tag,tuple(sample['case']))
                if key in log_replacements:
                    if sample['status']!='unresolved':raise ValueError('unexpected log replacement')
                    new=log_replacements.pop(key)
                    if new['exponents']!=sample['exponents']:raise ValueError('log branch exponents changed')
                    sample=new
                added[sample['status']]+=1
                evidence.append(sample['status'])
            past=initial_counts.get((name,tag),{})
            rechecks=field['rechecks']
            old_flags=original_flags.get((name,tag),set())
            if old_flags!={tuple(q['case']) for q in rechecks if tuple(q['case']) in old_flags}:
                raise ValueError('unchecked initial disagreement')
            for r in rechecks:
                targeted[r['status']]+=1
                precision_examples.setdefault(r['status'],dict(object=name,field=tag,**r))
            evidence.extend(r['status'] for r in rechecks)
            if any(v in ('separated','persistent_difference') for v in evidence):
                status='sample_difference'
            elif any(v in ('unresolved','poorly_resolved') for v in evidence) or past.get('poorly_resolved') or past.get('unresolved'):
                status='sample_coverage_incomplete'
            elif evidence and initial not in ('exact_identity','sample_coverage_incomplete'):
                status='sample_agreement'
            else:status=initial
            refined[name][tag]=status
            if tag=='velocity_part_2/residual':
                refined[name]['velocity_part_2/object']=status
    if log_replacements:raise ValueError('orphan log recheck')
    if sum(targeted.values())!=sum(len(flags) for (name,field),flags in original_flags.items() if field!='velocity_part_2/object'):
        raise ValueError('initial precision flag tally differs')
    final_counts=defaultdict(Counter)
    for (name,field),old in initial_status.items():
        final_counts[field][refined[name].get(field,old)]+=1
    for n in required:
        if n.startswith('both_surfaces/sum_crossing/'):
            report('signed_crossing',object=n,counts={field:refined[n].get(field,initial_status[n,field])
                   for field in ('Q','E_geo','E_resp','velocity_part_1/residual','velocity_part_2/residual')})
    report('precision_rechecks',counts=targeted)
    for r in precision_examples.values():report('precision_example',**r)
    report('unpaired_internal_kernel_occurrences',occurrences=sum(len(v) for v in aux['unpaired_internal_kernel_occurrences'].values()))

    report('additional_branch_samples',counts=added)
    for (n,j),objects in joint_counts.items():
        report('compound_branch_inventory',objects=objects,elementary_boundaries=n,explicit_wolfram_intersections=j)
    for field in ['Q','E_geo','E_resp','separation_derivative','departure_partition/residual',
                  'item3/object','item3/Q_closed_operand','item3/family_integral','item3/exponent_branch_inventory',
                  'velocity_part_1/residual','velocity_part_2/residual']:
        report('final_field',field=field,counts=final_counts[field])
    for k,v in sorted(family_forms.items()):
        chosen=[n for n in required if native_forms[n]==k]
        report('final_family_result',form=k,counts=Counter(refined[n].get('item3/object',initial_status[n,'item3/object']) for n in chosen))
    for field,pairs in sorted(final_counts.items()):
        if field not in ('Q','E_geo','E_resp','item3/family_integral','velocity_part_1/residual','velocity_part_2/residual'):
            continue
        if pairs.get('sample_coverage_incomplete'):
            report('still_incomplete',field=field,objects=[n for n in required if refined[n].get(field,initial_status[n,field])=='sample_coverage_incomplete'])
    dimension_fields=[k for k in final_counts if '/dimension_' in k]
    report('all_dimension_records',total=sum(sum(final_counts[k].values()) for k in dimension_fields),
           exact=sum(final_counts[k]['exact_identity'] for k in dimension_fields),fields=len(dimension_fields))
    for q in exceptions:report('exception',**q)
    publish = opt.publication or LEDGER
    (publish/'scripts/out').mkdir(parents=True,exist_ok=True)
    (publish/'directives/_measurements').mkdir(parents=True,exist_ok=True)
    if opt.review_controls:
        run = guarded(opt.review_controls)
        report('native_mutations',verbatim=(opt.review_controls/'stdout').read_text().splitlines())
    outpath=publish/'scripts/out/S15_reading_comparator.out'
    # Operand archive includes the comparison referring to it and each source occurrence. Local regenerable output.
    archive=outpath.with_name('S15_reading_comparator_details.out')
    tmp=archive.with_suffix('.out.tmp')
    with tarfile.open(tmp,'w') as tf:
        for p in sorted(root.rglob('*')):
            if p.is_file() and p.name not in ('manifest.json',):
                tf.add(p,arcname=str(p.relative_to(root)),recursive=False)
    tmp.replace(archive)
    scripts={**config['sources']}
    source_dir=Path(__file__).resolve().parent
    def source_path(name):
        candidate=source_dir/name
        return candidate if candidate.exists() else LEDGER/'scripts'/name
    if opt.review_controls:
        name='S15_reading_review_controls.py'
        scripts[name]=digest(source_path(name))
    for fn in ('S15_reading_refinement.py','S15_reading_log_branches.py','S15_reading_context_checks.py','test_S15_reading_comparator.py',Path(__file__).name):
        scripts[fn]=digest(source_path(fn))
    for fn,h in scripts.items():
        if digest(source_path(fn))!=h:raise ValueError('comparator source changed: '+fn)
    runs={name:guarded(opt.guard_parent/name) for name in (root.name,opt.context_run,'precision01','jointlogs01')}
    unit_guard=opt.unit_controls or opt.guard_parent/'controls01'
    runs['controls01']=guarded(unit_guard)
    if opt.review_controls:
        runs['native_controls']=guarded(opt.review_controls)
    if (opt.guard_parent/'divisors01').exists():
        runs['divisors01']=guarded(opt.guard_parent/'divisors01')
    controls=(unit_guard/'stderr').read_text().splitlines()
    if not controls or controls[-1]!='OK':raise ValueError('builder controls not successful')
    manifest=dict(status='Builder run; independent build review pending',implementation=config['build_id'],
                  decision_commit='2603dc7c',engine_commits={'reading':'cc62dc9f','original':'3192123b'},
                  previous_comparator_commit='5a019394',repository_inputs=config['inputs'],source_sha256=scripts,
                  guard_source_sha256=digest(REPO/'scripts/s11c_guarded_run.py'),runs=runs,
                  engine_sources={},reproduction_commands=[],controls_verbatim=controls,
                  versions=dict(python=sys.version,sympy=version('sympy'),scipy=version('scipy'),numpy=version('numpy')),
                  detail_archive=dict(path=relative(archive),bytes=archive.stat().st_size,sha256=digest(archive),
                                      format='POSIX tar of JSON and gzip JSON records'),
                  representative_exceptions_verbatim=[json.dumps({'record':'exception',**ex},ensure_ascii=False,
                      sort_keys=True,separators=(',',':')) for issues in sorted(examples.items()) for ex in issues[1]])
    import subprocess
    for label,p,rev in [('decision',LEDGER/'directives/S15_reading_comparator_decision_list.md','2603dc7c'),
                        ('original manifest',MEASURE/'S15_build_run_record.md','3192123b'),
                        ('SymPy manifest',MEASURE/'S15_separation_departure_sympy_run_manifest.json','cc62dc9f'),
                        ('Wolfram manifest',MEASURE/'S15_separation_departure_wolfram_run_manifest.json','cc62dc9f')]:
        locator=rev+':'+relative(p)
        content=subprocess.check_output(['git','show',locator],cwd=REPO)
        manifest.setdefault('accepted_artifacts',{})[label]=dict(git_object=locator,sha256=hashlib.sha256(content).hexdigest())
    for name,run in runs.items():
        part=' --all-cpus-exclusive' if len(run['actual_limits']['affinity'])>1 else ''
        limit=int(run['actual_limits']['memory.max'])//2**30
        import shlex
        manifest['reproduction_commands'].append('python3 scripts/s11c_guarded_run.py --pool s15-reading-comparator'+
            f' --memory-gib {limit}{part} --log-directory '+run['directory']+' -- '+shlex.join(run['guard_child_command']))
    for folder,ext in [('scripts','py'),('mathematica','wl')]:
        for fn in ('S15_separation_departure','S15_two_drain_region'):
            p=LEDGER/folder/f'{fn}.{ext}'
            manifest['engine_sources'][relative(p)]=digest(p)
    manifest['guarded_stage_source_notes']={
        'context01':'Initial kernel and argument checks. Its initial Boolean-tag tally is superseded by divisors01.',
        'divisors01':'Reads native Boolean False nodes; rechecks only divisor records.',
        'precision01':'45- and 75-decimal closed forms, exact stationary substitution and compound exponent branches.',
        'jointlogs01':'Li0(z)=z/(1-z) and Li1(z)=-log(1-z) for real |z|<1, previously unexercised log cases only.'}
    baseline=opt.root.parent/'pre_precision/S15_reading_context_checks.py'
    if baseline.exists():
        manifest['runs'][opt.context_run]['source_sha256_at_launch']=digest(baseline)
    manifest['report_command']=shlex.join([sys.executable,*sys.argv])
    if opt.review_controls:
        manifest['guarded_stage_source_notes']={
            opt.context_run:'Signed endpoints, full azimuth bounds, argument slots and native divisor flags.',
            'precision01':'Closed-form precision refinements and compound exponent branches on the repaired profiles.',
            'jointlogs01':'Remaining evaluable log-branch cases.',
            'native_controls':'One-sided centre, domain, endpoint, azimuth and three native total response edits.',
            'controls01':'Unit controls including 78 exact stationary velocity checks in two backends.'}
        manifest['reused_exact']=config['reused_exact']
    manifest['runs'][root.name]['end_receipt_utc']=read(root.parent/(root.name+'.DONE.json'))['ended']
    manifest['versions']['mpmath']=version('mpmath')
    from S15_reading_coverage import publish as publish_coverage
    coverage_path=publish/'directives/_measurements/S15_reading_comparator_coverage.json'
    coverage_counts, slot_counts=publish_coverage(manifest, archive, REPO, coverage_path)
    manifest['final_coverage']=dict(path=relative(coverage_path),sha256=digest(coverage_path),
                                    exponent_coverage=coverage_counts,argument_counts=slot_counts,
                                    source_sha256=digest(source_path('S15_reading_coverage.py')))
    report('item3_exponent_coverage',**coverage_counts)
    report('observed_argument_checks',counts=slot_counts,
           python_empty=[q['argument'] for q in aux['argument_slots'] if not q['observed']])
    tmp=outpath.with_suffix('.out.tmp')
    with tmp.open('w') as f:
        for line in lines:f.write(line+'\n')
    tmp.replace(outpath)
    manifest['transcript']=dict(path=relative(outpath),bytes=outpath.stat().st_size,sha256=digest(outpath))
    manifest['transcript_summary_lines_verbatim']=[line for line in lines if json.loads(line)['record']!='exception']
    path=publish/'directives/_measurements/S15_reading_comparator_run_manifest.json'
    path.write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
    for kind in ('roster','reading_to_original','kernel_checks','observed_argument_checks'):
        for line in lines:
            if json.loads(line)['record']==kind:print(line)
    print('validated',len(required),'object checkpoints;',len(lines),'transcript lines')


if __name__=='__main__':
    main()
