#!/usr/bin/env python3
"""Collect literal guard/output evidence; performs no symbolic computation.

The input is a JSON list of attempts with id, phase, guard_directory, optional
work_directory and supervisor_file. An attempt record may additionally carry
an explicitly sourced historical stop account. The output is an evidence
packet, not an inference that the build is complete.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()

def read(path):
    path=Path(path)
    return json.loads(path.read_text()) if path.exists() else None

def collect(spec):
    guard=Path(spec['guard_directory'])
    record={k:v for k,v in spec.items() if k not in
            ('guard_directory','work_directory','supervisor_file')}
    record['guard_invocation']=read(guard/'invocation.json')
    record['guard_outcome']=read(guard/'outcome.json')
    record['child_outcome']=read(guard/'child-outcome.json')
    for name in ('limit-validation','duration-validation','effective-limits'):
        record[name.replace('-','_')]=read(guard/(name+'.json'))
    record['supervisor']=read(spec['supervisor_file']) if spec.get('supervisor_file') else None
    samples=guard/'resource-samples.jsonl'
    first=last=peak_sample=None
    peak=-1;count=0
    if samples.exists():
        with samples.open() as f:
            for line in f:
                sample=json.loads(line);count+=1
                if first is None: first=line.rstrip('\n')
                last=line.rstrip('\n')
                value=int(sample['memory.peak'])
                if value>peak: peak=value;peak_sample=last
    record['resources']={'sample_count':count,'peak_bytes':None if peak<0 else peak,
                         'first_sample_literal':first,'last_sample_literal':last,
                         'peak_sample_literal':peak_sample,
                         'samples_sha256':sha(samples) if samples.exists() else None}
    record['transcripts']={}
    record['quoted_stdout_lines']=[]
    record['object_end_records']=[]
    for name in ('stdout','stderr'):
        path=guard/name
        if not path.exists():
            record['transcripts'][name]=None;continue
        record['transcripts'][name]={'size_bytes':path.stat().st_size,'sha256':sha(path)}
        if name=='stdout':
            with path.open() as f:
                for line in f:
                    line=line.rstrip('\n')
                    if line.startswith(('OBJECT_END:','ENGINE_END:','ATTEMPT_ERROR:',
                                        'WORKING_ON:','BASIS_SIZE:','ADDITION_COUNT:',
                                        'SELECTED_INDICES:','CONDITION_COUNT:',
                                        'SINGULAR_COMPLETED:','SINGULAR_VERSION:')):
                        record['quoted_stdout_lines'].append(line)
                    if line.startswith('OBJECT_END: '):
                        record['object_end_records'].append(json.loads(line.partition(': ')[2]))
    work=Path(spec['work_directory']) if spec.get('work_directory') else None
    record['engine_run_data']=read(work/'run_data.json') if work else None
    record['generated_inputs']=[]
    if work and work.exists():
        record['generated_inputs']=[{'name':p.name,'size_bytes':p.stat().st_size,'sha256':sha(p)}
                                    for p in sorted(work.glob('object_*.sing'))]
    record['missing_records']=[name for name,key in
            [('invocation.json','guard_invocation'),('outcome.json','guard_outcome'),
             ('child-outcome.json','child_outcome'),('run_data.json','engine_run_data')]
             if record[key] is None]
    return record

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--attempt-list',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    args=p.parse_args()
    specs=json.loads(args.attempt_list.read_text())
    packet={'schema':'s15-singular-attempt-evidence-v1',
            'collected_utc':datetime.now(timezone.utc).isoformat(),
            'collector_sha256':sha(__file__),
            'attempts':[collect(spec) for spec in specs]}
    args.output.write_text(json.dumps(packet,indent=2)+'\n')

if __name__=='__main__':
    main()
