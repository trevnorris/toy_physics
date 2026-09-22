#!/usr/bin/env python3
"""Stdlib-only final metadata handoff; completed scientific inputs stay saved."""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import time

REPO = Path('/var/projects/toy_physics')
M = REPO/'research/pde_ledger_v3/_measurements'
F = REPO/'_scratch/s11c/s11c-remaining-case-frequency-20260921'
OLD = F/'lab-operator-inputs'
NAME = 'S11c_d_remaining_case_frequency_lab_operator_inputs_finish'
CHECKS_SHA = '589fa276597ff5a74f07e8efde3d91e78ba052e7782fcfb3e9e04dca52f13e15'
SOURCE_SHA = '9ee5165cce2654ecbe7112420221a2c39c3b3ad8b5a33405d0da631f11d4e39d'
INVENTORY_SHA = '44f63d5a299d8d4565c2308d1d18c5adc56f72ada4879015d571d8f45d033b9d'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''): h.update(block)
        os.posix_fadvise(stream.fileno(), 0, 0, os.POSIX_FADV_DONTNEED)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--run-directory', type=Path, required=True)
    base = parser.parse_args().run_directory.resolve(); base.relative_to(F); base.mkdir(parents=True, exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);start=time.monotonic()
    paths = {}; artifacts = {}
    def retain(path, expected=None):
        p = Path(path); h = digest(p)
        if expected is not None: assert h == expected, str(p)
        record = {'path': str(p), 'canonical': str(p.resolve()), 'sha256': h, 'bytes': p.stat().st_size,
                  'directLink': str(p.readlink()) if p.is_symlink() else None}
        paths[str(p)] = record; return record
    def read(path, expected=None): retain(path, expected); return json.loads(Path(path).read_text())
    def save(name, value):
        payload = json.dumps(value, indent=2)+'\n'; p = base/name
        assert not p.exists() and not p.is_symlink();p.parent.mkdir(parents=True,exist_ok=True)
        with p.open('x') as out:out.write(payload);out.flush();os.fsync(out.fileno())
        record = retain(p);artifacts[name] = record;return record
    inventory_path = F/'lab-operator-inputs-finish/original-file-inventory.json'
    inventory = read(inventory_path, INVENTORY_SHA)
    for path, record in inventory.items(): assert retain(path, record['sha256']) == record
    source_path = M/'S11c_d_remaining_case_frequency_lab_operator_inputs.py'
    source = retain(source_path, SOURCE_SHA); tree = ast.parse(source_path.read_text())
    function = next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
    # Preserve exact completed native metadata order and identify the alias.
    tail = [ast.unparse(n) for n in function.body[-5:]]
    assert "'artifacts': journal.artifacts" in ast.unparse(function)
    assert tail[-3:] == ["journal.json('checks.json', checks)", 'signal.alarm(0)', 'print(json.dumps(checks, indent=2))']
    checks = read(OLD/'complete/checks.json', CHECKS_SHA)
    printed = read(OLD/'frequency_lab_operator_inputs.stdout')
    self_record = printed['artifacts'].pop('checks.json')
    assert printed == checks and self_record['sha256'] == CHECKS_SHA and self_record['path'] == str(OLD/'complete/checks.json')
    outcomes = [read(OLD/n) for n in ('resource-guard/outcome.json','resource-guard/child-outcome.json','frequency_lab_operator_inputs.invocation.json')]
    assert all(v['exitCode']==0 for v in outcomes) and outcomes[1]['guardReason'] is None
    limits = read(OLD/'resource-guard/effective-limits.json')
    assert limits['memory.max']=='2147483648' and limits['memory.swap.max']=='0' and limits['pids.max']=='32' and limits['nice']==15 and len(limits['affinity'])==1 and set(limits['threads'].values())=={'1'}
    for name in ('frequency_lab_operator_inputs.stderr','guard.stderr','resource-guard/stderr'): assert retain(OLD/name)['bytes']==0
    for line in (OLD/'resource-guard/resource-samples.jsonl').read_text().splitlines():
        sample=json.loads(line);events=dict(v.split() for v in sample['memory.events'].splitlines())
        assert int(sample['memory.swap.current'])==0 and int(sample['memory.peak'])<=2*1024**3 and all(int(events[k])==0 for k in ('max','oom','oom_kill'))
    for path, record in read(OLD/'complete/inputs.json')['consumedRoutes'].items():
        current = retain(path,record['sha256'])
        assert all(current[k]==record[k] for k in ('canonical','directLink','bytes'))
    for name, record in checks['artifacts'].items():
        original = retain(record['path'],record['sha256']); target = base/name
        target.parent.mkdir(parents=True,exist_ok=True);assert not target.exists() and not target.is_symlink()
        target.symlink_to(Path(original['canonical']));artifacts[name] = retain(target,record['sha256'])
    save('completed-input-emission-join.json',{'originalChecks':retain(OLD/'complete/checks.json',CHECKS_SHA),
        'originalStdout':retain(OLD/'frequency_lab_operator_inputs.stdout'),'source':source,'sourceTail':tail,
        'printedSelfEntry':self_record,'soleOriginalDifference':'Journal appended its checks record to the aliased artifacts dictionary after writing checks bytes.',
        'originalIdentityPassed':False,'originalOutcomes':outcomes,'originalInventory':retain(inventory_path,INVENTORY_SHA),
        'newScientificImports':0,'newScientificReads':0,'newScientificCalls':0,'originalCompleteSummariesDirectlyReferenced':True})
    for path in (Path(__file__).resolve(),M/(NAME+'_plan.md')):
        rec=retain(path);target=base/'source'/path.name;target.parent.mkdir(exist_ok=True)
        with target.open('xb') as out:out.write(Path(rec['canonical']).read_bytes())
        retain(target,rec['sha256'])
    for path, record in list(paths.items()): assert retain(path,record['sha256'])==record
    save('final-consumed-paths.json',dict(paths))
    result={'status':'COMPLETED_SAVED_LAB_OPERATOR_INPUT_HANDOFF','cases':checks['cases'],
        'originalRunDirectory':str(OLD/'complete'),'originalChecksSha256':CHECKS_SHA,'originalPrintedSelfEntryPreserved':True,
        'newScientificCalls':0,'newScientificReads':0,'allConsumedHashesUnchanged':True,'consumedLogicalPaths':len(paths),
        'artifacts':dict(artifacts),'wallSeconds':time.monotonic()-start,'scope':checks['scope']}
    save('checks.json',result);signal.alarm(0)
    # Print the exact durable bytes; save() cannot alter result.artifacts.
    print((base/'checks.json').read_text(),end='')


if __name__=='__main__': main()
