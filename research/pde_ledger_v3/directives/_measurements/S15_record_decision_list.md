# S15 record decision list measurements

Reviewed baseline. The [setup lookups](S15_record_setup_lookups.md) give the input fingerprints, native expressions and measured comparison statuses. Run from the repository root.

```sh
python3 - <<'PY'
from hashlib import sha256
from pathlib import Path
p = Path('research/pde_ledger_v3/directives')
for name in ('S15_record_decision_list.md', '_measurements/S15_record_setup_lookups.md'):
    print(name, sha256((p / name).read_bytes()).hexdigest())
PY
```

```text
S15_record_decision_list.md 4d16e5346522b84c13187282c983d9d920601d046ded862a6a36f2e3c554c944
_measurements/S15_record_setup_lookups.md 3c39ab92b79d1c3160edae7f98c58003798ba704654dbb55e02ca7133f09dbdf
```
