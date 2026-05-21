#!/usr/bin/env python3
"""
Project layout detector for redteam-audit bootstrap.

Walks a project tree, finds files with `stage<NUM>` in their name, clusters
them by (parent_dir, extension), derives glob templates with {N} replacing the
stage number and `*` replacing variable descriptive middles, and emits a YAML
detection report on stdout.

Usage:
  detect.py <project-root>
"""

import re
import sys
import yaml
from collections import defaultdict
from pathlib import Path

STAGE_RE = re.compile(r'stage[_]?(\d+)', re.IGNORECASE)

# Dirs we never descend into.
# Note: 'output'/'outputs' are NOT skipped — they often hold audit transcripts
# that we want to detect as sympy_output / mathematica_output roles.
SKIP_DIRS = {
    '.git', '.svn', '.hg',
    'node_modules', '__pycache__', '.pytest_cache', '.mypy_cache',
    '.venv', 'venv', 'env',
    '_archive', '_backup', 'archive', 'backup',
    'redteam',  # don't recursively re-detect our own scaffold
}

# Max directory depth from project root.
MAX_DEPTH = 6


def common_prefix(strings):
    if not strings:
        return ""
    s = sorted(strings)
    s1, s2 = s[0], s[-1]
    for i, c in enumerate(s1):
        if i >= len(s2) or c != s2[i]:
            return s1[:i]
    return s1


def common_suffix(strings):
    return common_prefix([s[::-1] for s in strings])[::-1]


def dominant_subset(templates):
    """Return the largest subset of templates that share the same prefix-before-{N}.
    A single outlier with a different prefix-before-{N} can otherwise collapse the
    common prefix to empty and produce a useless glob."""
    groups = defaultdict(list)
    for t in templates:
        idx = t.find('{N}')
        before = t[:idx] if idx >= 0 else t
        groups[before].append(t)
    if not groups:
        return list(templates)
    # Pick the group with the most templates; tie broken by longest prefix.
    best = max(groups.values(), key=lambda g: (len(g), len(g[0][:g[0].find('{N}')])))
    return best


def derive_glob(templates):
    """Given a set of filename templates (each with {N} for stage number),
    derive the smallest glob that matches all of them, using `*` for any
    variable middle text. Drops outliers whose before-{N} prefix differs
    from the dominant pattern."""
    templates = list(set(templates))
    if not templates:
        return None
    if len(templates) == 1:
        return templates[0]
    # Filter to the dominant before-{N} subset to avoid one outlier
    # collapsing the common prefix to empty.
    templates = dominant_subset(templates)
    if len(templates) == 1:
        return templates[0]
    prefix = common_prefix(templates)
    suffix = common_suffix(templates)
    sample = templates[0]
    if len(prefix) + len(suffix) >= len(sample):
        return prefix
    return prefix + "*" + suffix


def classify_role(rel_dir, ext, glob, sample_files=None):
    """Heuristic role assignment based on extension, path, glob, and sample
    filenames. Sample names are inspected because the derived glob may collapse
    the signature (one outlier file can force the common suffix down to just
    '.py', erasing keywords that were in every other name)."""
    g = glob.lower()
    d = rel_dir.lower()
    samples = ' '.join((sample_files or [])).lower()
    blob = g + ' ' + samples

    if ext == '.tex':
        return 'tex'
    if ext == '.wl':
        return 'mathematica'
    if ext == '.md':
        if 'notes' in d or 'stages' in d or 'docs' in d:
            return 'notes'
        return None
    if ext == '.py':
        if 'audit' in blob or 'sympy' in blob or 'verify' in blob:
            return 'sympy'
        return None
    if ext == '.txt':
        if 'sympy' in blob:
            return 'sympy_output'
        if 'math' in blob:
            return 'mathematica_output'
        return None
    return None


def walk_project(root):
    """Yield records for every file under `root` whose name contains a stage number."""
    root = Path(root).resolve()
    for fpath in root.rglob('*'):
        if not fpath.is_file():
            continue
        rel_parts = fpath.relative_to(root).parts
        if any(p in SKIP_DIRS for p in rel_parts):
            continue
        if len(rel_parts) > MAX_DEPTH:
            continue
        m = STAGE_RE.search(fpath.name)
        if not m:
            continue
        digits = m.group(1)
        start, end = m.start(1), m.end(1)
        template = fpath.name[:start] + '{N}' + fpath.name[end:]
        yield {
            'rel_dir': str(fpath.parent.relative_to(root)),
            'name': fpath.name,
            'template': template,
            'pad': len(digits),
            'stage_num': int(digits),
            'ext': fpath.suffix,
        }


def find_parts_dirs(root):
    """Find candidate 'parts' / 'chapters' directories that hold tex files."""
    root = Path(root).resolve()
    candidates = ['parts', 'chapters', 'sections', 'chapter', 'part']
    found = []
    for c in candidates:
        for p in root.rglob(c):
            if not p.is_dir():
                continue
            rel = p.relative_to(root)
            if any(part in SKIP_DIRS for part in rel.parts):
                continue
            if len(rel.parts) > MAX_DEPTH:
                continue
            tex_count = sum(1 for _ in p.glob('*.tex'))
            if tex_count > 0:
                found.append({
                    'path': str(rel),
                    'tex_count': tex_count,
                    'files': sorted(f.name for f in p.glob('*.tex'))[:6],
                })
    return found


def main():
    if len(sys.argv) != 2:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    root = Path(sys.argv[1]).resolve()
    if not root.is_dir():
        print(f"error: not a directory: {root}", file=sys.stderr)
        sys.exit(2)

    files = list(walk_project(root))

    # Cluster by (rel_dir, ext).
    clusters_by_key = defaultdict(list)
    for f in files:
        clusters_by_key[(f['rel_dir'], f['ext'])].append(f)

    report = {
        'project_root': str(root),
        'project_name': root.name,
        'file_count': len(files),
        'stage_range': None,
        'stage_pads': {},
        'parts_dirs': find_parts_dirs(root),
        'clusters': [],
    }

    if files:
        nums = [f['stage_num'] for f in files]
        report['stage_range'] = [min(nums), max(nums)]
        pad_counter = defaultdict(int)
        for f in files:
            pad_counter[f['pad']] += 1
        report['stage_pads'] = dict(pad_counter)

    for (rel_dir, ext), members in sorted(clusters_by_key.items()):
        templates = [m['template'] for m in members]
        glob = derive_glob(templates)
        full_glob = glob if rel_dir == '.' else f"{rel_dir}/{glob}"
        sample_names = sorted({m['name'] for m in members})[:6]
        role = classify_role(rel_dir, ext, full_glob, sample_names)
        report['clusters'].append({
            'rel_dir': rel_dir,
            'ext': ext,
            'count': len(members),
            'glob': full_glob,
            'role': role,
            'sample_files': sample_names[:3],
        })

    yaml.safe_dump(report, sys.stdout, default_flow_style=False, sort_keys=False, width=120)


if __name__ == '__main__':
    main()
