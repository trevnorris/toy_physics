"""Locate the accepted reading records for the S15 comparator decision list."""
import argparse
import ast
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[3]
LEDGER = ROOT / 'research/pde_ledger_v3'
MEASUREMENTS = LEDGER / 'directives/_measurements'
TAG = re.compile(rb'^([A-Za-z][A-Za-z0-9_{},+ /-]*): ')
BASE = '3192123b'
READING = 'cc62dc9f'


def accepted(path, revision):
    return subprocess.check_output(['git', 'show', f'{revision}:{path.relative_to(ROOT)}'], cwd=ROOT)


def record_index(path, stop_tag):
    digest = sha256()
    records = {}
    manifest = []
    end = None
    previous = None
    with path.open('rb') as stream:
        while True:
            start = stream.tell()
            line = stream.readline()
            if not line:
                break
            digest.update(line)
            match = TAG.match(line)
            if not match:
                if previous is None:
                    raise ValueError(f'unattached continuation in {path} at {start}')
                previous['length'] += len(line)
                continue
            tag = match[1].decode()
            body_start = start + match.end()
            previous = {'offset': body_start, 'length': len(line) - match.end()}
            records.setdefault(tag, []).append(previous)
            if tag == 'Q_manifest':
                if manifest:
                    raise ValueError(f'duplicate Q_manifest in {path}')
                body = line[match.end():].decode().strip()
                manifest = (ast.literal_eval(body) if body.startswith('[')
                            else json.loads('[' + body[1:-1] + ']'))
            if tag == stop_tag:
                end = stream.tell()
                break
        if end is None or (stop_tag == 'END_OF_TRANSCRIPT' and stream.read(1)):
            raise ValueError(f'bad end marker in {path}')
    return records, manifest, digest.hexdigest(), end


def payload(stream, entry):
    stream.seek(entry['offset'])
    data = stream.read(entry['length'])
    if not data.endswith(b'\n'):
        raise ValueError('incomplete payload')
    return data[:-1]


def unique(records, key):
    entries = records.get(key, [])
    if len(entries) != 1:
        raise ValueError(f'{key!r} has {len(entries)} occurrences')
    return entries[0]


def sum_operands(name):
    match = re.fullmatch(r'(inner|outer|both_surfaces)/sum_(rest|transport|crossing)/(component_[1-3])', name)
    if not match:
        return []
    location, kind, component = match.groups()
    surfaces = ('inner', 'outer') if location == 'both_surfaces' else (location,)
    rest = ('carried', 'pressure', 'viscous', 'other')
    transport = ('transport_centre', 'transport_radius')
    pieces = {'rest': rest, 'transport': transport, 'crossing': rest + transport}[kind]
    return [f'{surface}/{piece}/{component}' for surface in surfaces for piece in pieces]


def build_index(side):
    wl = side == 'wolfram'
    manifest_file = MEASUREMENTS / f'S15_separation_departure_{side}_run_manifest.json'
    manifest = json.loads(accepted(manifest_file, READING))
    reading = ROOT / manifest['transcript']
    original = reading.with_name('S15_two_drain_region.out')
    record_file = MEASUREMENTS / 'S15_build_run_record.md'
    section = accepted(record_file, BASE).decode().split(f'## {manifest["engine"]}\n', 1)[1]
    expected_hash, = re.findall(r'- Transcript SHA256: `([0-9a-f]{64})`', section.split('\n## ', 1)[0])
    old, _, old_hash, old_bytes = record_index(original, 'END_OF_TRANSCRIPT')
    if old_hash != expected_hash:
        raise ValueError(f'{original} differs from the accepted transcript')
    stop = 'END_OF_BASE_PASS' if wl else 'END_OF_TRANSCRIPT'
    tags, names, digest, end = record_index(reading, stop)
    if stop == 'END_OF_BASE_PASS':
        verified = manifest['final_transcript_production']['verification']
        if (digest, end) != (verified['prefix_sha256'], verified['prefix_bytes']):
            raise ValueError('Wolfram base-pass prefix differs from the accepted manifest')
    elif digest != manifest['transcript_sha256']:
        raise ValueError('SymPy reading transcript differs from the accepted manifest')
    if not names or len(names) != len(set(names)):
        raise ValueError('invalid Q_manifest')
    printed_names = {tag.removeprefix('reading/').removesuffix('/Q')
                     for tag in tags if tag.startswith('reading/') and tag.endswith('/Q')}
    if set(names) != printed_names:
        raise ValueError('Q records differ from Q_manifest')
    fields = {
        '1': ['Q'],
        '2': ['E_geo', 'E_resp', 'separation_derivative', 'departure_existence_condition'],
        '3': ['item3/object', 'item3/family_integral', 'item3/existence_condition',
              'item3/exponent_branch_inventory', 'item3/status'],
        '5': [f'velocity_part_{k}/{field}' for k in (1, 2)
              for field in ('operand_a', 'operand_b', 'object' if wl else 'residual')],
        '6': [f'{field}/dimension_reference_residuals'
              for field in ('E_geo', 'E_resp', 'velocity_part_1', 'velocity_part_2')]
    }
    rows = {}
    groups = Counter()
    statuses = Counter()
    with reading.open('rb') as new_stream, original.open('rb') as old_stream:
        for name in names:
            key = re.sub(r'labels_\{(drain1|drain2|remainder)(?:, (?:drain1|drain2|remainder))*\}',
                         lambda m: 'labels_' + '+'.join(m[0].removeprefix('labels_{').removesuffix('}').split(', ')), name)
            if key in rows:
                raise ValueError(f'object-name collision for {key}')
            prefix = f'reading/{name}/'
            spans = {item: {f: unique(tags, prefix + f) for f in required}
                     for item, required in fields.items()}
            item2_literals = {field: payload(new_stream, spans['2'][field]) == b'0'
                              for field in ('E_geo', 'E_resp')}
            status = payload(new_stream, spans['3']['item3/status']).decode()
            statuses[status] += 1
            row = {'native_name': name, 'selected_records': spans, 'item3_native_status': status,
                   'item2_prints_literal_zero': item2_literals}
            if name in old:
                original_span = unique(old, name)
                identical = payload(new_stream, spans['1']['Q']) == payload(old_stream, original_span)
                row['original_span'] = original_span
                row['original_lookup'] = 'identical payload' if identical else 'different payload'
            elif operands := sum_operands(name):
                row['original_operand_spans'] = {tag: unique(old, tag) for tag in operands}
                row['original_lookup'] = 'new sum with all original operand tags located'
            else:
                row['original_lookup'] = 'no identical tag'
            groups[row['original_lookup']] += 1
            rows[key] = row
    return {
        'engine': side, 'original': str(original.relative_to(ROOT)), 'original_bytes': old_bytes,
        'original_sha256': old_hash, 'reading': str(reading.relative_to(ROOT)), 'read_through': stop,
        'reading_bytes_checked': end, 'reading_bytes_sha256': digest, 'original_lookup': dict(groups),
        'item3_native_status_counts': dict(statuses), 'objects': rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index', type=Path, help='optional local file for native tags and byte offsets')
    args = parser.parse_args()
    sympy, wolfram = (build_index(side) for side in ('sympy', 'wolfram'))
    report = {
        'accepted_reading_commit': READING, 'accepted_original_commit': BASE,
        'index_scope': 'Selected native records, byte lookups and print-form differences only.',
        'sympy': sympy, 'wolfram': wolfram,
    }
    if args.index is not None:
        with args.index.open('x') as output:
            json.dump(report, output, indent=2)
            output.write('\n')
    for item in (sympy, wolfram):
        print(item['engine'], 'original bytes and sha256:', item['original_bytes'], item['original_sha256'])
        print(item['engine'], item['read_through'], 'bytes and sha256:',
              item['reading_bytes_checked'], item['reading_bytes_sha256'])
        print('manifest objects with selected tags:', len(item['objects']))
        print('item 1 same-engine lookup:', item['original_lookup'])
        print('item 3 status tag text:', item['item3_native_status_counts'])
    py_names, wl_names = set(sympy['objects']), set(wolfram['objects'])
    print('common object names after label spelling normalization:', len(py_names & wl_names))
    print('SymPy-only names:', sorted(py_names - wl_names))
    print('Wolfram-only names:', sorted(wl_names - py_names))
    for field in ('E_geo', 'E_resp'):
        differences = [
            name for name in sorted(py_names & wl_names)
            if sympy['objects'][name]['item2_prints_literal_zero'][field]
            != wolfram['objects'][name]['item2_prints_literal_zero'][field]
        ]
        print(field, 'literal-zero spelling differs across engines:', len(differences))
    for item in (sympy, wolfram):
        differences = [v['native_name'] for v in item['objects'].values() if v['original_lookup'] == 'different payload']
        print(item['engine'], 'changed native payloads with same tag:', differences)
        sums = [row for row in item['objects'].values() if 'original_operand_spans' in row]
        print(item['engine'], 'sum rows:', len(sums), 'original operand references:',
              sum(len(row['original_operand_spans']) for row in sums))

if __name__ == '__main__':
    main()
