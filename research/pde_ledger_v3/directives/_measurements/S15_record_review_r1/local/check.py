"""Read-only checks of candidate M2 and M3 against accepted native payloads and archive."""
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tarfile


base = Path('research/pde_ledger_v3')
candidate = Path('_scratch/s15/record_build/r1/candidate') / base
file = candidate / 'steps/_measurements/S15_record_measurements.md'
text = file.read_text()
m2 = text.split('## M2.')[1].split('## M3.')[0]
printed = dict((addr, (field, val)) for field, addr, val in re.findall(
    r'^(\S+) ((?:py|wl)_reading:\d+) (.+)$', m2, re.M))
by_tag = re.findall(r'^### (\d{3})\. (.+)$', text, re.M)
manifest = json.loads(subprocess.check_output(['git','show',
    f'41adb784:{base}/directives/_measurements/S15_reading_comparator_run_manifest.json']))
archive = base / 'scripts/out/S15_reading_comparator_details.out'
assert hashlib.sha256(archive.read_bytes()).hexdigest() == manifest['detail_archive']['sha256']
native = {k:open(manifest['repository_inputs'][k]['path'],'rb') for k in ('wl_reading','py_reading')}
counts, exact = Counter(), Counter()
checked = set()
with tarfile.open(archive) as tar:
    spans = {entry['id']: entry for entry in json.loads(gzip.decompress(tar.extractfile('coverage.json.gz').read()))}
    for idx, name in by_tag:
        obj = json.loads(gzip.decompress(tar.extractfile(f'objects/{idx}.json.gz').read()))
        assert (int(idx), name) == (obj['index'], obj['object'])
        for row in obj['rows']:
            field = row['field']
            if row['role'] != 'cross' or field not in ('item3/Q_closed_operand','item3/object','velocity_part_1/residual','velocity_part_2/residual'):
                continue
            vals=[]
            for addr in (*row['left'], *row['right']):
                engine, offset = addr.split(':')
                entry = spans[addr]
                native[engine].seek(int(offset))
                val = native[engine].read(entry['span']['length']).removesuffix(b'\n')
                assert hashlib.sha256(val).hexdigest() == entry['payload']
                if addr in printed:
                    assert (field, val.decode()) == printed[addr]
                    checked.add(addr)
                vals.append(val)
            if field.startswith('velocity'):
                literal = [val == b'0' for val in vals]
                assert literal[0] == literal[1]
                counts[(field,literal[0])] += 1
                if name.startswith('both_surfaces/sum_crossing'):
                    exact[(field,row['status'])] += 1
assert len(by_tag) == 336
assert len(checked) == len(printed) == 52
print('M2 native values agreeing with original accepted archive spans and stream bytes:',len(checked))
print('M3 indexed object names and indices matching original archive:',len(by_tag))
print('Paired literal zeros vs other native spellings:', sorted(counts.items()))
print('Net velocity cross-engine statuses:', sorted(exact.items()))
