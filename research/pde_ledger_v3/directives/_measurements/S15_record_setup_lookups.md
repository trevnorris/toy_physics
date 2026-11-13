# S15 record decision list evidence lookups

Run each command from the repository root. The first checks the local products against the manifest filed at `41adb784`. The Wolfram reading input here ends at the accepted base-pass byte boundary. Later Wolfram passes concern item 4. The second retrieves one component of the signed shell crossing and the two native expressions for a single drain-2 centre-transport piece. The item-3 rows refer to the original quantity on the profile family and its geometric departure, respectively. [Accepted comparator report](S15_reading_comparator_run_record.md) · [Definition of the family and its departure](../S15_reading_decision_list.md)

## Accepted comparator products and reading inputs

```sh
python3 - <<'PY'
from hashlib import sha256
import json
from pathlib import Path
import subprocess

base = Path('research/pde_ledger_v3')
manifest_path = base / 'directives/_measurements/S15_reading_comparator_run_manifest.json'
manifest = json.loads(subprocess.check_output(['git', 'show', f'41adb784:{manifest_path}']))
for key in ('transcript', 'detail_archive', 'final_coverage'):
    entry = manifest[key]
    directory = 'directives/_measurements' if key == 'final_coverage' else 'scripts/out'
    path = base / directory / Path(entry['path']).name
    digest = sha256(path.read_bytes()).hexdigest()
    assert digest == entry['sha256']
    print(path.name, digest)
for key in ('py_reading', 'wl_reading'):
    entry = manifest['repository_inputs'][key]
    remaining = entry['bytes']
    digest = sha256()
    with Path(entry['path']).open('rb') as stream:
        while remaining:
            data = stream.read(min(remaining, 1024 * 1024))
            assert data
            digest.update(data)
            remaining -= len(data)
    assert digest.hexdigest() == entry['sha256']
    print(key, entry['bytes'], digest.hexdigest())
PY
```

```text
S15_reading_comparator.out 473997968d83faee9843fff93a07bf45c6ea9c77fbf1251567d41f5c19f33b54
S15_reading_comparator_details.out 206109b070a649e8bbcd9c6a77766acd1e8a98290a17c7ec9408af9609741ae4
S15_reading_comparator_coverage.json 141d8a7fc93cfc28d2b6fab6a99e3cc61fd99965a3e9c016355363d7838c76e6
py_reading 192065031 10160f89704d910139b0401fbbcc1419078e45e915a2f1e663f24807647182d6
wl_reading 333550018 75622128d04da2d934c2a52468505b4438fff9f17a902d7a6cd35993d769b6c4
```

## Signed shell crossing and one native transport term

```sh
python3 - <<'PY'
import gzip
import json
from pathlib import Path
import tarfile

base = Path('research/pde_ledger_v3')
streams = {
    'py_reading': base / 'scripts/out/S15_separation_departure.out',
    'wl_reading': base / 'mathematica/out/S15_separation_departure.out',
}
cases = {
    '000': ('Q', 'E_geo', 'E_resp', 'velocity_part_1/residual', 'velocity_part_2/residual'),
    '080': ('item3/Q_closed_operand', 'item3/object'),
}
with tarfile.open(base / 'scripts/out/S15_reading_comparator_details.out') as archive:
    for index, fields in cases.items():
        with gzip.open(archive.extractfile(f'objects/{index}.json.gz')) as records:
            obj = json.load(records)
        print(obj['object'])
        for field in fields:
            row, = (r for r in obj['rows'] if r['role'] == 'cross' and r['field'] == field)
            left, = row['left']
            right, = row['right']
            print(field, row['status'], left, right)
            if index == '080':
                for locator in (left, right):
                    engine, offset = locator.split(':')
                    with streams[engine].open('rb') as native:
                        native.seek(int(offset))
                        print(native.readline().decode().rstrip())
PY
```

```text
both_surfaces/sum_crossing/component_1
Q exact_identity py_reading:175986249 wl_reading:189014051
E_geo exact_identity py_reading:176242291 wl_reading:191886167
E_resp sample_agreement py_reading:176425206 wl_reading:194089093
velocity_part_1/residual exact_identity py_reading:177374202 wl_reading:201177044
velocity_part_2/residual exact_identity py_reading:177466790 wl_reading:201850880
inner/transport_centre/component_1/labels_drain2
item3/Q_closed_operand exact_identity py_reading:83732737 wl_reading:60330854
-4*pi*A_2(t)*B_2(t)*r_in(t)**(a_2 + b_2 + 2)*Derivative(X_2_1(t), t)/3
-((2*(Pi - Pi*directionSeparation[1]^2) + (2*(-Pi + 3*Pi*directionSeparation[1]^2))/3)*rin[t]^(3 + 2*(-1/2 + aa[2]/2 + bb[2]/2))*AA[2][t]*BB[2][t]*Derivative[1][xx[2, 1]][t])
item3/object exact_identity py_reading:83733127 wl_reading:60331358
-8*pi*A_2(t)*B_2(t)*r_in(t)**(a_2 + b_2 + 2)*Derivative(X_2_1(t), t)/3
-2*(2*(Pi - Pi*directionSeparation[1]^2) + (2*(-Pi + 3*Pi*directionSeparation[1]^2))/3)*rin[t]^(3 + 2*(-1/2 + aa[2]/2 + bb[2]/2))*AA[2][t]*BB[2][t]*Derivative[1][xx[2, 1]][t]
```
