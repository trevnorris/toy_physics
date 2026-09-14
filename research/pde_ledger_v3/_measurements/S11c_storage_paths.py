#!/usr/bin/env python3
"""Check or restore historical path aliases for repository-resident S11c runs.

No calculation payload is written. New runs should use the durable paths
directly; aliases preserve old absolute paths in frozen manifests and packets.
"""
import argparse
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
DEFAULT = Path(__file__).with_name('S11c_storage_migration_checkpoint.json')


def inspect(entries, repository):
    durable = repository/'_scratch/s11c'
    rows = []
    for item in entries:
        old, new = Path(item['oldPath']), Path(item['durablePath'])
        if old.parent != Path('/tmp') or 's11c' not in old.name.lower():
            raise ValueError(('unexpected compatibility path', str(old)))
        if new.parent != durable or new.name != old.name:
            raise ValueError(('unexpected durable path', str(new)))
        if not new.exists() or new.is_symlink():
            raise ValueError(('durable payload missing or redirected', str(new)))
        if old.is_symlink():
            if old.resolve() != new:
                raise ValueError(('conflicting compatibility link', str(old)))
            status = 'present'
        elif old.exists():
            raise ValueError(('existing file would be overwritten', str(old)))
        else:
            status = 'missing'
        rows.append((old, new, status))
    if len({old for old, _, _ in rows}) != len(rows):
        raise ValueError('duplicate compatibility paths')
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', type=Path, default=DEFAULT)
    parser.add_argument('--restore', action='store_true',
                        help='Recreate missing aliases; refuse every conflicting path.')
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    rows = inspect(manifest['entries'], REPO)
    created = []
    if args.restore:
        for old, new, status in rows:
            if status == 'missing':
                old.symlink_to(new, target_is_directory=new.is_dir())
                created.append(str(old))
        rows = inspect(manifest['entries'], REPO)
    missing = [str(old) for old, _, status in rows if status == 'missing']
    print(json.dumps({'durableRoot': str(REPO/'_scratch/s11c'),
                      'checkedPaths': len(rows), 'createdLinks': created,
                      'missingLinks': missing}, indent=2))
    if missing:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
