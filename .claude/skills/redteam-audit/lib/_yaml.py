#!/usr/bin/env python3
"""
YAML read/write helper for the redteam-audit skill's bash CLI.

Subcommands:
  get FILE KEY              Print value at dotted KEY (e.g. project.name)
  get-list FILE KEY         Print list items at KEY, one per line
  set FILE KEY VALUE        Set scalar VALUE at KEY (atomic write).
                            VALUE is parsed as YAML; pass --raw to force string.
  set-yaml FILE KEY         Read YAML from stdin, set at KEY (atomic write)
  has FILE KEY              Exit 0 if KEY exists, 1 otherwise
  keys FILE KEY             Print top-level keys of mapping at KEY
  dump FILE                 Pretty-print full document
  init FILE                 Write {} to FILE if it doesn't exist
  batches FILE              Print "id<TAB>start<TAB>end<TAB>label" per batch
  stages-in-batch FILE ID   Print stage numbers (zero-padded) in batch ID

Dotted keys descend into mappings and lists (integer index for lists).
List indices are written as `[i]` or accessed by integer key.
"""
import os
import sys
import yaml
from pathlib import Path

def load(path):
    if not Path(path).exists():
        return {}
    with open(path) as f:
        return yaml.safe_load(f) or {}

def atomic_write(path, data):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + ".tmp")
    with open(tmp, 'w') as f:
        yaml.safe_dump(data, f, default_flow_style=False, sort_keys=False, width=120)
    tmp.replace(p)

def descend(obj, key_path):
    """Walk dotted key path through dicts/lists. Return value or raise KeyError."""
    parts = key_path.split('.')
    cur = obj
    for p in parts:
        if isinstance(cur, list):
            cur = cur[int(p)]
        elif isinstance(cur, dict):
            if p not in cur:
                raise KeyError(key_path)
            cur = cur[p]
        else:
            raise KeyError(key_path)
    return cur

def set_path(obj, key_path, value):
    """Set value at dotted key path, creating intermediate dicts as needed."""
    parts = key_path.split('.')
    cur = obj
    for p in parts[:-1]:
        if isinstance(cur, list):
            cur = cur[int(p)]
        else:
            if p not in cur or not isinstance(cur[p], (dict, list)):
                cur[p] = {}
            cur = cur[p]
    last = parts[-1]
    if isinstance(cur, list):
        cur[int(last)] = value
    else:
        cur[last] = value

def cmd_get(args):
    data = load(args[0])
    try:
        v = descend(data, args[1])
    except (KeyError, IndexError, TypeError):
        sys.exit(2)
    if v is None:
        print("")
    elif isinstance(v, (dict, list)):
        print(yaml.safe_dump(v, default_flow_style=False, sort_keys=False).rstrip())
    else:
        print(v)

def cmd_get_list(args):
    data = load(args[0])
    try:
        v = descend(data, args[1])
    except (KeyError, IndexError, TypeError):
        sys.exit(2)
    if not isinstance(v, list):
        sys.exit(2)
    for item in v:
        if isinstance(item, (dict, list)):
            print(yaml.safe_dump(item, default_flow_style=True, sort_keys=False).strip())
        else:
            print(item)

def cmd_has(args):
    data = load(args[0])
    try:
        descend(data, args[1])
        sys.exit(0)
    except (KeyError, IndexError, TypeError):
        sys.exit(1)

def cmd_keys(args):
    data = load(args[0])
    try:
        v = data if args[1] == "" else descend(data, args[1])
    except (KeyError, IndexError, TypeError):
        sys.exit(2)
    if isinstance(v, dict):
        for k in v.keys():
            print(k)
    elif isinstance(v, list):
        for i in range(len(v)):
            print(i)

def cmd_set(args):
    fp, key, val_str = args[0], args[1], args[2]
    raw = "--raw" in args[3:]
    data = load(fp)
    if raw:
        val = val_str
    else:
        try:
            val = yaml.safe_load(val_str)
        except yaml.YAMLError:
            val = val_str
    set_path(data, key, val)
    atomic_write(fp, data)

def cmd_set_yaml(args):
    fp, key = args[0], args[1]
    payload = yaml.safe_load(sys.stdin.read())
    data = load(fp)
    set_path(data, key, payload)
    atomic_write(fp, data)

def cmd_dump(args):
    data = load(args[0])
    print(yaml.safe_dump(data, default_flow_style=False, sort_keys=False, width=120).rstrip())

def cmd_init(args):
    fp = args[0]
    if not Path(fp).exists():
        atomic_write(fp, {})

def cmd_batches(args):
    cfg = load(args[0])
    batches = cfg.get("batches", []) or []
    for b in batches:
        bid = b.get("id", "?")
        rng = b.get("range", [0, 0])
        label = b.get("label", "")
        print(f"{bid}\t{rng[0]}\t{rng[1]}\t{label}")

def cmd_stages_in_batch(args):
    cfg = load(args[0])
    pad = int(cfg.get("stage_pad", 3))
    target_id = args[1]
    excludes = set(cfg.get("stages", {}).get("exclude", []) or [])
    for b in cfg.get("batches", []) or []:
        if str(b.get("id")) == target_id:
            rng = b.get("range", [0, 0])
            for n in range(int(rng[0]), int(rng[1]) + 1):
                if n in excludes:
                    continue
                print(f"{n:0{pad}d}")
            return
    sys.exit(2)

CMDS = {
    "get": cmd_get,
    "get-list": cmd_get_list,
    "has": cmd_has,
    "keys": cmd_keys,
    "set": cmd_set,
    "set-yaml": cmd_set_yaml,
    "dump": cmd_dump,
    "init": cmd_init,
    "batches": cmd_batches,
    "stages-in-batch": cmd_stages_in_batch,
}

def main():
    if len(sys.argv) < 2 or sys.argv[1] not in CMDS:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    CMDS[sys.argv[1]](sys.argv[2:])

if __name__ == "__main__":
    main()
