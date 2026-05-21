#!/usr/bin/env python3
"""
Fast manifest queries — read MANIFEST.yaml once per invocation.

Used by the bash CLI to avoid spawning hundreds of _yaml.py processes for
operations that touch all stages.

Subcommands:
  summary               Per-batch status table (stdout).
  batch-status <ID>     Stage table for one batch (stdout).
  state-list <STATE>    List "stage<TAB>batch<TAB>status" lines, one per stage in STATE.
  blocked               List blocked_* stages.
  next-batch            Print the earliest batch with any audit-eligible unit
                        (status in {pending, upstream_stale}).
  render-batches        Regenerate <redteam-dir>/BATCHES.md from manifest.
  mark-stale-downstream <FROM>  Set upstream_stale=true on every stage > FROM.
  set-stage-fields <NNN> [--file -]  Read YAML mapping from stdin/file and merge
                          into stages.<NNN>. Atomic write.
  append-stage-list <NNN> <FIELD>  Append the stdin string as a new element
                          of stages.<NNN>.<FIELD> (creates the list if absent).

Usage:
  manifest.py <config-path> <subcommand> [args...]
"""
import os
import sys
import yaml
from collections import defaultdict
from datetime import datetime
from pathlib import Path

def load_yaml(path):
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

def manifest_path(cfg, cfg_path):
    project_root = Path(cfg_path).parent
    rd = cfg.get("project", {}).get("redteam_dir", "redteam")
    return project_root / rd / "MANIFEST.yaml"

def batches_md_path(cfg, cfg_path):
    project_root = Path(cfg_path).parent
    rd = cfg.get("project", {}).get("redteam_dir", "redteam")
    return project_root / rd / "BATCHES.md"

def iter_batches(cfg):
    """Yield (id, start, end, label) tuples preserving config order."""
    for b in cfg.get("batches", []) or []:
        rng = b.get("range", [0, 0])
        yield b.get("id"), int(rng[0]), int(rng[1]), b.get("label", "")

def stage_pad(cfg):
    return int(cfg.get("stage_pad", 3))

def stages_in_batch(cfg, bid):
    pad = stage_pad(cfg)
    for b in cfg.get("batches", []) or []:
        if str(b.get("id")) == str(bid):
            rng = b.get("range", [0, 0])
            return [f"{n:0{pad}d}" for n in range(int(rng[0]), int(rng[1]) + 1)]
    return []

def batch_counts(manifest, cfg):
    """Return [(bid, total, {state: count}, label), ...]"""
    stages = manifest.get("stages", {}) or {}
    out = []
    for bid, start, end, label in iter_batches(cfg):
        pad = stage_pad(cfg)
        counts = defaultdict(int)
        total = 0
        for n in range(start, end + 1):
            key = f"{n:0{pad}d}"
            entry = stages.get(key)
            if entry is None:
                continue
            counts[entry.get("status", "(none)")] += 1
            total += 1
        out.append((bid, total, dict(counts), label))
    return out

# --- subcommands ------------------------------------------------------------

def cmd_summary(cfg, cfg_path, args):
    m = load_yaml(manifest_path(cfg, cfg_path))
    rows = batch_counts(m, cfg)
    print(f"Project: {m.get('project_name', '?')}")
    print(f"Manifest: {manifest_path(cfg, cfg_path)}")
    print()
    print(f"{'BATCH':<8}  {'TOTAL':>5}  STATES")
    for bid, total, counts, _ in rows:
        cells = " ".join(f"{k}={v}" for k, v in counts.items()) or "(empty)"
        print(f"{bid:<8}  {total:>5}  {cells}")

def cmd_batch_status(cfg, cfg_path, args):
    if not args:
        print("usage: batch-status <ID>", file=sys.stderr)
        sys.exit(2)
    bid = args[0]
    m = load_yaml(manifest_path(cfg, cfg_path))
    stages = m.get("stages", {}) or {}
    keys = stages_in_batch(cfg, bid)
    if not keys:
        print(f"unknown batch: {bid}", file=sys.stderr)
        sys.exit(2)
    print(f"Batch {bid}")
    print()
    print(f"{'STAGE':<6}  {'STATUS':<30}  {'SYMPY':<8}  {'MATH':<8}  {'CHECKPT':<8}")
    for key in keys:
        e = stages.get(key, {})
        status = e.get("status", "(none)")
        files = e.get("files", {}) or {}
        sympy_exists = (files.get("sympy") or {}).get("exists", "?")
        math_exists = (files.get("mathematica") or {}).get("exists", "?")
        cp = "yes" if e.get("is_checkpoint") else ""
        print(f"{key:<6}  {status:<30}  {str(sympy_exists):<8}  {str(math_exists):<8}  {cp:<8}")

def cmd_state_list(cfg, cfg_path, args):
    if not args:
        print("usage: state-list <STATE>", file=sys.stderr)
        sys.exit(2)
    target = args[0]
    m = load_yaml(manifest_path(cfg, cfg_path))
    for key, e in (m.get("stages", {}) or {}).items():
        if e.get("status") == target:
            print(f"{key}\t{e.get('batch_id', '')}\t{e.get('status')}")

def cmd_blocked(cfg, cfg_path, args):
    m = load_yaml(manifest_path(cfg, cfg_path))
    rd = cfg.get("project", {}).get("redteam_dir", "redteam")
    for key, e in (m.get("stages", {}) or {}).items():
        s = e.get("status", "")
        if s.startswith("blocked_"):
            print(f"{key}\t{e.get('batch_id', '')}\t{s}\t{rd}/reports/stage_{key}.md")

def cmd_next_batch(cfg, cfg_path, args):
    """Print the next batch with at least one audit-eligible unit.

    Audit-eligible = status in {pending, upstream_stale}. All other states
    are either already done (verified) or awaiting human input (blocked_*)
    or mid-loop (auditing/fixing/codex_applied/verifying/etc.). We pick the
    EARLIEST batch by config order so cascade-affected upstream batches
    surface before later ones — a fresh material change in batch I.1 forces
    the orchestrator to re-audit the stale parts of I.1 / II.1 before
    advancing.
    """
    m = load_yaml(manifest_path(cfg, cfg_path))
    stages = m.get("stages", {}) or {}
    pad = stage_pad(cfg)
    ELIGIBLE = {"pending", "upstream_stale"}
    for bid, start, end, _ in iter_batches(cfg):
        for n in range(start, end + 1):
            key = f"{n:0{pad}d}"
            e = stages.get(key)
            if e is None:
                continue
            if e.get("status") in ELIGIBLE:
                print(bid)
                return
    print("(no batch with audit-eligible units)")

def cmd_render_batches(cfg, cfg_path, args):
    m = load_yaml(manifest_path(cfg, cfg_path))
    rows = batch_counts(m, cfg)
    md_path = batches_md_path(cfg, cfg_path)
    lines = [
        "# Red-Team Batch Status",
        "",
        f"Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"Project: {m.get('project_name', '?')}",
        "",
        "| Batch | Range | Stages | States | Label |",
        "|---|---|---:|---|---|",
    ]
    for bid, start, end, label in iter_batches(cfg):
        # find matching row
        match = next((r for r in rows if r[0] == bid), (bid, 0, {}, label))
        total = match[1]
        counts = match[2]
        cells = " ".join(f"{k}={v}" for k, v in counts.items()) or "(empty)"
        lines.append(f"| {bid} | {start}–{end} | {total} | {cells} | {label} |")
    md_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.write_text("\n".join(lines) + "\n")
    print(f"rendered: {md_path}")

def cmd_mark_stale_downstream(cfg, cfg_path, args):
    """Mark every unit > FROM as upstream_stale.

    A unit's status is the source of truth for what audit/verify see. Setting
    only the upstream_stale flag would leave a previously `verified` unit
    still counting as verified in status queries and next-batch — the very
    thing the cascade is supposed to prevent. So we also DEMOTE the status:

      pending                       → stays pending (will be audited anyway)
      blocked_unfixable             → stays blocked (human decision pending)
      blocked_critical_downstream   → stays blocked
      upstream_stale                → stays upstream_stale (already there)
      anything else (incl. verified, audited, directive_ready, codex_applied,
                     verifying, needs_rework, fixing, auditing)
                                    → upstream_stale (must be re-audited)

    For demoted units, also clear codex_session — the old session's context
    is about a now-stale premise.
    """
    if not args:
        print("usage: mark-stale-downstream <FROM>", file=sys.stderr)
        sys.exit(2)
    from_n = int(args[0])
    mp = manifest_path(cfg, cfg_path)
    m = load_yaml(mp)
    stages = m.get("stages", {}) or {}
    KEEP = {"pending", "blocked_unfixable", "blocked_critical_downstream", "upstream_stale"}
    demoted = 0
    untouched = 0
    for key, e in stages.items():
        if e.get("stage", 0) > from_n:
            cur = e.get("status", "")
            if cur not in KEEP:
                e["status"] = "upstream_stale"
                e["upstream_stale"] = True
                e["codex_session"] = ""
                demoted += 1
            else:
                # Pending units will get a fresh audit anyway; blocked units
                # need a human decision regardless. We leave both status and
                # the flag alone so flag stays synced with status everywhere.
                untouched += 1
    atomic_write(mp, m)
    print(f"mark-stale-downstream: demoted {demoted} unit(s) to upstream_stale "
          f"(codex sessions cleared); {untouched} unit(s) left as-is")

def cmd_set_stage_fields(cfg, cfg_path, args):
    if not args:
        print("usage: set-stage-fields <NNN> [--file PATH]", file=sys.stderr)
        sys.exit(2)
    key = args[0]
    src = sys.stdin
    if len(args) >= 3 and args[1] == "--file":
        src = open(args[2])
    payload = yaml.safe_load(src.read())
    if not isinstance(payload, dict):
        print("payload must be a YAML mapping", file=sys.stderr)
        sys.exit(2)
    mp = manifest_path(cfg, cfg_path)
    m = load_yaml(mp)
    stages = m.setdefault("stages", {})
    entry = stages.setdefault(key, {})
    entry.update(payload)
    atomic_write(mp, m)

def cmd_append_stage_list(cfg, cfg_path, args):
    if len(args) < 2:
        print("usage: append-stage-list <NNN> <FIELD>", file=sys.stderr)
        sys.exit(2)
    key, field = args[0], args[1]
    value = sys.stdin.read().strip()
    if not value:
        return
    mp = manifest_path(cfg, cfg_path)
    m = load_yaml(mp)
    stages = m.setdefault("stages", {})
    entry = stages.setdefault(key, {})
    lst = entry.get(field) or []
    if not isinstance(lst, list):
        lst = []
    lst.append(value)
    entry[field] = lst
    atomic_write(mp, m)


CMDS = {
    "summary": cmd_summary,
    "batch-status": cmd_batch_status,
    "state-list": cmd_state_list,
    "blocked": cmd_blocked,
    "next-batch": cmd_next_batch,
    "render-batches": cmd_render_batches,
    "mark-stale-downstream": cmd_mark_stale_downstream,
    "set-stage-fields": cmd_set_stage_fields,
    "append-stage-list": cmd_append_stage_list,
}

def main():
    if len(sys.argv) < 3 or sys.argv[2] not in CMDS:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    cfg_path = sys.argv[1]
    cfg = load_yaml(cfg_path)
    CMDS[sys.argv[2]](cfg, cfg_path, sys.argv[3:])

if __name__ == "__main__":
    main()
