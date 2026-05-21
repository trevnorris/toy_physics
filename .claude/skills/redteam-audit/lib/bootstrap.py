#!/usr/bin/env python3
"""
Generate a starter .redteam-config.yaml from a detect.py scan.

Usage:
  bootstrap.py <project-root> [--batch-size N] [--project-name NAME]
                              [--prefer-role role=glob ...] [--out PATH]
                              [--force]

Internally calls detect.py to scan the project, picks the strongest cluster
for each role (most files wins; ties broken by directory depth), generates a
batch map (uniform-N-per-batch by default), and writes a config to
<project-root>/.redteam-config.yaml unless --out overrides.

Refuses to overwrite an existing config unless --force is passed.

Always prints a "review and customize" checklist to stdout — the generated
config is a STARTING POINT, not finished.
"""

import argparse
import subprocess
import sys
from pathlib import Path

import yaml


def detect(project_root, detect_script):
    out = subprocess.check_output(
        ['python3', str(detect_script), str(project_root)],
        text=True,
    )
    return yaml.safe_load(out)


def parse_prefer(items):
    """Parse --prefer-role role=glob entries."""
    out = {}
    for kv in items or []:
        if '=' not in kv:
            print(f"warning: ignoring malformed --prefer-role: {kv}", file=sys.stderr)
            continue
        role, glob = kv.split('=', 1)
        out[role.strip()] = glob.strip()
    return out


def pick_role_globs(clusters, preferred):
    """For each role, pick the cluster with the most files. Honor --prefer-role overrides."""
    by_role = {}
    # Sort clusters: highest count first, shallower dirs first as tiebreaker.
    sorted_clusters = sorted(
        clusters,
        key=lambda c: (-c.get('count', 0), len(Path(c.get('rel_dir', '')).parts)),
    )
    for c in sorted_clusters:
        role = c.get('role')
        if not role:
            continue
        if role not in by_role:
            by_role[role] = c['glob']
    by_role.update(preferred)
    return by_role


def build_batch_map(stage_range, batch_size):
    """Uniform batches of `batch_size` stages each, ID'd B1, B2, ..."""
    start, end = stage_range
    batches = []
    n = start
    idx = 1
    while n <= end:
        batch_end = min(n + batch_size - 1, end)
        batches.append({
            'id': f"B{idx}",
            'label': f"Batch {idx} — stages {n}-{batch_end} (rename to reflect content)",
            'range': [n, batch_end],
        })
        n = batch_end + 1
        idx += 1
    return batches


def build_config(report, args):
    role_globs = pick_role_globs(report.get('clusters', []), parse_prefer(args.prefer_role))

    stage_range = report.get('stage_range') or [1, 1]

    # Pad: most common pad observed.
    pad_counts = report.get('stage_pads') or {}
    if pad_counts:
        pad = int(max(pad_counts.items(), key=lambda kv: kv[1])[0])
    else:
        pad = 3

    batches = build_batch_map(stage_range, max(1, args.batch_size))

    # The red-team only tracks script-related roles. tex/notes are detected
    # (and shown in the detection summary) but not emitted into the manifest —
    # doc alignment is handled manually, out-of-band.
    paths = {}
    for role in ['sympy', 'mathematica', 'sympy_output', 'mathematica_output']:
        if role in role_globs:
            paths[role] = [role_globs[role]]

    runners = {}
    if 'sympy' in role_globs:
        runners['sympy'] = 'python3 {script}'
    if 'mathematica' in role_globs:
        runners['mathematica'] = 'math -script {script}'

    return {
        'project': {
            'name': args.project_name or report.get('project_name') or 'unknown',
            'redteam_dir': 'redteam',
        },
        'stage_pad': pad,
        'paths': paths,
        'runners': runners,
        'stages': {
            'start': stage_range[0],
            'end': stage_range[1],
            'exclude': [],
        },
        'batches': batches,
        'checkpoints': [],
        'limits': {
            'max_iterations': 2,
            'parallel_audit_max': 10,
            'parallel_verify_max': 10,
            'fix_parallelism': 1,
        },
        'codex': {
            'chat_wrapper': '~/.claude/hooks/codex-chat/codex-chat',
            'sandbox': 'workspace-write',
            'dry_run': False,
        },
        'status_only_candidates': [],
    }


def print_review_summary(out_path, report, config):
    print(f"wrote: {out_path}")
    print()
    print("Detection summary:")
    print(f"  files found:    {report.get('file_count', 0)}")
    sr = report.get('stage_range')
    if sr:
        print(f"  stage range:    {sr[0]}-{sr[1]}")
    print(f"  tracked script roles: {sorted(config['paths'].keys())}")
    # Surface tex/notes detection separately — they are NOT tracked by the
    # red-team (scripts-only audit), but knowing they exist helps the user
    # confirm the layout is what they expected.
    untracked = sorted(
        c.get('role') for c in report.get('clusters', [])
        if c.get('role') in ('tex', 'notes') and c.get('role')
    )
    if untracked:
        print(f"  detected (not tracked): {sorted(set(untracked))}")
    print(f"  batches:        {len(config['batches'])} (size {config['batches'][0]['range'][1] - config['batches'][0]['range'][0] + 1 if config['batches'] else 0})")
    parts_dirs = report.get('parts_dirs') or []
    if parts_dirs:
        print()
        print("  Candidate parts/chapters dirs found (consider using these for batch boundaries):")
        for d in parts_dirs:
            print(f"    - {d['path']}: {d['tex_count']} tex file(s)")
    print()
    print("Review and customize the generated config:")
    print(f"  1. Open {out_path}")
    print("  2. Verify the path globs match the actual file layout")
    print("  3. Add checkpoint stage numbers to `checkpoints:`")
    print("  4. Add stages with intentionally no executable audit to `status_only_candidates:`")
    print("  5. Rename batch IDs / labels to reflect content (currently auto-numbered)")
    print("  6. Adjust `runners:` if you use non-standard interpreters")
    print()
    print("Then initialize the manifest:")
    print(f"  cd {Path(out_path).parent} && <skill-path>/lib/redteam.sh init")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('project_root')
    p.add_argument('--batch-size', type=int, default=12,
                   help='stages per batch (default 12)')
    p.add_argument('--project-name', help='override project name (default: dir name)')
    p.add_argument('--prefer-role', action='append', default=[],
                   help='override role->glob mapping, e.g. --prefer-role sympy=scripts/audits/*_stage{N}*.py')
    p.add_argument('--out', help='output path (default: <project-root>/.redteam-config.yaml)')
    p.add_argument('--force', action='store_true', help='overwrite existing config')
    args = p.parse_args()

    skill_lib = Path(__file__).parent
    project_root = Path(args.project_root).resolve()

    if not project_root.is_dir():
        print(f"error: not a directory: {project_root}", file=sys.stderr)
        sys.exit(2)

    detect_script = skill_lib / 'detect.py'
    if not detect_script.exists():
        print(f"error: detect.py not found at {detect_script}", file=sys.stderr)
        sys.exit(2)

    report = detect(project_root, detect_script)

    out_path = Path(args.out) if args.out else project_root / '.redteam-config.yaml'
    if out_path.exists() and not args.force:
        print(f"error: config already exists at {out_path}", file=sys.stderr)
        print("use --force to overwrite, or --out to write elsewhere", file=sys.stderr)
        sys.exit(1)

    config = build_config(report, args)

    if not config['paths']:
        print("error: no stage files detected. Nothing to bootstrap.", file=sys.stderr)
        print("Project tree should contain files with `stage<NUM>` in their names.", file=sys.stderr)
        sys.exit(1)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(yaml.safe_dump(config, default_flow_style=False, sort_keys=False, width=120))

    print_review_summary(out_path, report, config)


if __name__ == '__main__':
    main()
