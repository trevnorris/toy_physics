#!/usr/bin/env python3
"""Inventory the frozen S11c-era tree; metadata only, no science or annex fetch.

Run from any directory: python3 path/to/inventory.py [--check]
Writes only INVENTORY.tsv beside this script; prints the workstream/kind table.
--check compares existing TSV bytes without writing. Rules classify file roles,
not scientific validity. First matching rule wins. No rename inference is used.
"""
import argparse
import collections
import csv
import io
from pathlib import Path
import re
import subprocess

BASE = "dada3b7d"
END = "archive/pre-cleanup-2026-10-04"
ROOT = Path(__file__).resolve().parents[3]
DEST = Path(__file__).with_name("INVENTORY.tsv")
STREAMS = {
    1: "01 Physics specification and contracts",
    2: "02 Symbolic build and omega1 premises",
    3: "03 Upstream repairs and diagnostics",
    4: "04 Numerical scattering and omega3 benchmark",
    5: "05 Near-unity and defect prerequisites",
    6: "06 Clean-condition packet",
    7: "07 Lean and S9-S11 formalization audits",
    8: "08 Exploratory throat and EM",
    9: "09 Muonium gravity corrections",
    10: "10 Execution infrastructure",
    11: "11 Ledger front matter and other notes",
}
KINDS = (
    "spec/amendment/contract", "decision list or build directive", "result report",
    "review report or disposition", "script that produces a reported result",
    "worker/continuation/tooling", "process record", "output (.out/.json result)",
    "Lean proof/contract", "Lean CAS bridge", "doc/note", "ledger front matter",
)


def git(*args):
    return subprocess.check_output(["git", "-C", str(ROOT), *args])


def workstream(path):
    p = path.lower()
    n = Path(p).name
    if "muonium" in p:
        return 9
    if p.startswith("docs/"):
        return 11 if n == "s11_maccullagh_differentiation.md" else 8
    if any(s in n for s in ("core_response", "conversion_face_work", "elastic_reference_comparison")):
        return 8
    if "clean_condition" in p or "zinvariant_operator_blocks" in p:
        return 6
    if "/lean/" in p or re.match(r"s(?:9|10|11)(?:_|\.)", n):
        return 7
    if (p.startswith("scripts/") or n == "agents.md" or n == ".gitattributes"
            or n == ".gitignore" or any(s in n for s in
                ("guarded_run", "review_cas", "desktop_freeze", "resource_guard",
                 "freeze_incident", "host_freeze", "review_policy", "no_deadline",
                 "s11c_parallel_", "s11c_storage_"))):
        return 10
    if any(s in n for s in ("near_unity", "s11c_d_defect_", "s11c_d_first_order_",
                            "physics_refocus", "nonuniform_work", "background_force",
                            "nonuniform_background", "centre_drive")):
        return 5
    if (re.match(r"s11c_[abc][12]?_", n) or any(s in n for s in
            ("upstream_", "inertia_", "mechanical_", "thickness_coordinate",
             "trace_repair", "wolfram_pressure", "wolfram_repair", "sheet_continuation",
             "sheet_path", "sheet_repair", "fourier_sheet", "transverse_sign"))):
        return 3
    if any(s in n for s in ("numerical_", "finite_balance", "saved_balance")):
        return 4
    if any(s in n for s in ("shared_physics", "spec_v", "nonlinear_pole",
                            "scattering_form", "total_transverse_loss",
                            "channel_map_born", "channel_map_review", "exploratory_acceptance")):
        return 1
    if n.startswith("s11c_d_"):
        return 2
    return 11


def kind(path):
    p = path.lower()
    n = Path(p).name
    ext = Path(p).suffix
    # Bridge role wins for its bindings/generators as well as Lean source.
    if "/s10audit/cas/" in p or "s10_lean_cas_" in n or "cas_bridge" in n:
        return KINDS[9]
    if ext == ".lean" or ("/lean/" in p and any(s in n for s in
            ("coverage", "fidelity.md", "lakefile", "lean-toolchain"))):
        return KINDS[8]
    if ext in (".py", ".wl", ".sh"):
        # Primary engines/audits and their directly named companion reports.
        # Launchers/tests/readers/continuations remain tooling even if cited.
        tooling = ("test", "launch", "runner", "supervisor", "continue", "recover",
                   "inspect", "read", "watch", "guard", "review", "codec", "exports",
                   "prepare", "manifest", "dispatch", "protocol", "resume", "_lib")
        if not any(s in n for s in tooling) and (
                "_audit." in n or "_probe." in n or "_diagnostic." in n or
                "_comparator." in n or n[:-len(ext)] + "_report.md" in NAMES or
                n[:-len(ext)] + "_result.md" in NAMES):
            return KINDS[4]
        return KINDS[5]
    if ext == ".out" or n.endswith((".stdout", ".stdout.txt")):
        return KINDS[7]
    if ext == ".json":
        if "_run_" in n or "_progress" in n:
            return KINDS[6]
        # Scientific inputs/outputs/checks are not decoded. Process-role names
        # take precedence; JSON result_record is a receipt, not the return.
        if any(s in n for s in ("record", "readiness", "authority", "authorization",
                "launch", "receipt", "checkpoint", "manifest", "packet", "gate",
                "preparation", "completion", "execution", "resource", "snapshot",
                "hash", "_run.", "_runs.", "_state.", "review", "registration",
                "source_inspection", "static_", "tooling", "disposition")):
            # 'packet' in S11c_d_defect_packet is a scientific workstream name.
            stripped = n.replace("s11c_d_defect_packet_", "")
            if any(s in stripped for s in ("record", "readiness", "authority", "authorization",
                    "launch", "receipt", "checkpoint", "manifest", "packet", "gate",
                    "preparation", "completion", "execution", "resource", "snapshot",
                    "hash", "_run.", "_runs.", "_state.", "review", "registration",
                    "source_inspection", "static_", "tooling", "disposition")):
                return KINDS[6]
        return KINDS[7]
    if ext in (".stderr", ".log", ".gz", ".kdl") or any(s in n for s in
            ("_prompt", "evidence_guide", "review_guide", "_index.txt", "source_excerpt",
             ".stderr.txt", ".time.txt", "_time.txt", "completion_message", "rebuild_stack")):
        return KINDS[6]
    if any(s in n for s in ("disposition", "adjudication", "assessment", "review",
                            "_claude", "_grok", "_opus", "_codex")):
        return KINDS[3]
    if any(s in n for s in ("report", "result", "verification", "regeneration_status",
                            "endpoint_status", "applicability_status")):
        return KINDS[2]
    if any(s in n for s in ("directive", "decision", "program_brief", "build.md",
                            "implementation", "_plan.md", "_plan.txt", "build_notes")) and n != "v3_step_plan.md":
        return KINDS[1]
    if any(s in n for s in ("shared_physics", "amendment", "contract", "method",
                            "clean_condition.md", "formalization_policy")):
        return KINDS[0]
    if (n in ("status.md", "readme.md", "claude.md", "agents.md", ".gitignore",
              ".gitattributes", "v3_step_plan.md", "deferred_heavy_runs.md")
            or "/steps/" in p):
        return KINDS[11]
    return KINDS[10]


def inventory():
    parts = git("diff", "--name-status", "--no-renames", "-z", BASE, END).split(b"\0")[:-1]
    changes = list(zip(parts[::2], parts[1::2]))
    assert len(parts) == 2 * len(changes)
    assert collections.Counter(s for s, _ in changes) == {b"A": 4434, b"M": 24}
    tree = {}
    for record in git("ls-tree", "-rlz", END).split(b"\0"):
        if record:
            metadata, path = record.split(b"\t", 1)
            mode, typ, oid, size = metadata.split()
            tree[path] = (mode, typ, oid, int(size))
    rows = []
    annex = 0
    for status, rawpath in changes:
        path = rawpath.decode()
        mode, typ, oid, size = tree[rawpath]
        assert typ == b"blob"
        if mode == b"120000":
            target = git("cat-file", "blob", oid.decode()).decode()
            key_size = re.search(r"/[^/]*-s(\d+)--", target)
            if "/annex/objects/" in target and key_size:
                size = int(key_size.group(1))
                annex += 1
            else:
                raise ValueError(f"Unclassified symlink byte size: {path}")
        rows.append((path, status.decode(), size, STREAMS[workstream(path)], kind(path)))
    assert len({r[0] for r in rows}) == len(rows)
    return sorted(rows), annex


def check_selection(live_exports=False):
    """Phase 3 metadata only; H3 keeps result outputs, not a cache replay package."""
    import ast
    import hashlib
    import json
    import posixpath
    folder = DEST.parent
    def table(name):
        with (folder / name).open() as stream:
            return list(csv.DictReader(stream, delimiter='\t'))
    kept_rows, pruned_rows = table('KEEP.tsv'), table('PRUNE.tsv')
    kept = {r['path']: r for r in kept_rows}
    assert len(kept) == len(kept_rows), 'Duplicate KEEP paths'
    pruned = {r['path'] for r in pruned_rows}
    original = {r['path'] for r in table('INVENTORY.tsv')}
    archived = set(git('ls-tree', '-r', '--name-only', END).decode().splitlines())
    tracked = set(git('ls-files').decode().splitlines())
    assert not set(kept) & pruned
    assert original == (set(kept) | pruned) & original
    assert pruned <= archived
    def sha(path):
        with path.open('rb') as stream:
            h = hashlib.sha256()
            for block in iter(lambda: stream.read(1024*1024), b''):
                h.update(block)
            return h.hexdigest()
    output_prefix = 'research/pde_ledger_v3/_measurements/S11c_d_near_unity_uniform_output/'
    outputs = 0
    for name, row in kept.items():
        assert (ROOT/name).is_file(), name
        assert not name.startswith('research/pde_ledger_v3/_measurements/retained_uniform/'), name
        assert 'S11c_uniform_retained_' not in name, name
        if row['sha256']:
            assert sha(ROOT/name) == row['sha256'], name
        if name.startswith(output_prefix):
            assert name.endswith('.json') and row['sha256'], name
            json.loads((ROOT/name).read_text())  # JSON syntax only, no expression evaluation.
            outputs += 1
    assert outputs == 31
    def candidates(name, owner):
        path = Path(name)
        choices = [path] if path.is_absolute() else [owner.parent/path, ROOT/'research/pde_ledger_v3'/path,
            ROOT/'research/pde_ledger_v3/scripts'/path, ROOT/'research/pde_ledger_v3/directives'/path, ROOT/path]
        found = []
        for choice in choices:
            choice = Path(posixpath.normpath(str(choice)))
            try: rel = str(choice.relative_to(ROOT))
            except ValueError: continue
            if choice.is_file() and rel not in pruned and (rel in kept or rel in tracked): found.append(choice)
        return found
    pins = []
    def walk(value, owner):
        if isinstance(value, dict):
            if isinstance(value.get('path'), str) and re.fullmatch('[0-9a-f]{64}', str(value.get('sha256', ''))):
                pins.append((owner, value['path'], value['sha256'], value))
            for key, item in value.items():
                if '/' in str(key) and isinstance(item, str) and re.fullmatch('[0-9a-f]{64}', item):
                    pins.append((owner, key, item, {}))
                elif '/' in str(key) and isinstance(item, dict) and re.fullmatch('[0-9a-f]{64}', str(item.get('sha256',''))):
                    pins.append((owner, key, item['sha256'], item))
                walk(item, owner)
        elif isinstance(value, list):
            for item in value: walk(item, owner)
    for name, row in kept.items():
        # Cited outputs retain historical metadata, not a promise to preserve runtime caches (H3).
        if row['tier']=='1' and name.endswith('.json') and not name.startswith(output_prefix):
            walk(json.loads((ROOT/name).read_text()), ROOT/name)
    for owner, name, expected, detail in pins:
        def matches(path):
            if 'byteOffset' in detail:
                with path.open('rb') as stream:
                    stream.seek(detail['byteOffset'])
                    return hashlib.sha256(stream.read(detail['bytes'])).hexdigest()==expected
            return sha(path)==expected
        assert any(matches(p) for p in candidates(name, owner)), (owner,name,expected)
    exports = 0
    live_bad = []
    known_edits = {'S10_brane_mode_spectrum_sympy_audit.py':'56595cf7',
                   'S11_stray_longitudinal_sympy_audit.py':'035bb654'}
    for path in sorted((ROOT/'research/pde_ledger_v3/scripts').glob('*_exports.py')):
        lines = path.read_text().splitlines(True)
        start = next((i for i,line in enumerate(lines) if line.startswith('BUILD_INPUT_DIGESTS = ')), None)
        if start is None: continue
        for end in range(start+1, len(lines)+1):
            try: value = ast.parse(''.join(lines[start:end])).body[0].value; break
            except SyntaxError: continue
        mapping = ast.literal_eval(value.args[0] if isinstance(value, ast.Call) else value)
        for name, expected in mapping.items():
            exports += 1
            if not any(sha(p)==expected for p in candidates(name,path)):
                assert name in known_edits, (path,name,expected)
                producer = 'research/pde_ledger_v3/scripts/'+name
                assert hashlib.sha256(git('show', BASE+':'+producer)).hexdigest()==expected
                assert sha(ROOT/producer)==hashlib.sha256(git('show',known_edits[name]+':'+producer)).hexdigest()
                live_bad.append((path.name,name))
    assert len(live_bad)==2, live_bad
    for engine, directory in [('PY','scripts'),('WL','mathematica')]:
        source = ROOT/f'research/pde_ledger_v3/lean/s10/S10Audit/CAS/{engine}.lean'
        expected = source.read_text().splitlines()[1].split()[-1]
        producer = 'sympy' if engine=='PY' else 'mathematica'
        assert sha(ROOT/f'research/pde_ledger_v3/{directory}/out/S10_anisotropic_strata_{producer}_audit.out')==expected
    # H5: exact frozen bytes, with historical links interpreted at the archive revision.
    frozen = {'research/pde_ledger_v3/directives/S11c_d_SCATTERING_FORM_AMENDMENT.md'}
    for name in frozen: assert (ROOT/name).read_bytes()==git('show',END+':'+name)
    local, archive_links = 0, 0
    for name,row in kept.items():
        if name in frozen or not name.endswith('.md') or re.search(r'(claude|grok|opus)(_|\.|/)|/(_legs|_reviews)/',name,re.I): continue
        text = (ROOT/name).read_text()
        links = re.findall(r'\]\(([^)]+)\)',text) + re.findall(r'^\[[^\]]+\]:\s*(\S+)',text,re.M)
        for target in links:
            target=re.sub(r':\d+(?:[–-]\d+)?$','',target.strip('<>').split('#')[0])
            if not target or target.startswith(('http','app:','codex:','mailto:')): continue
            if target.startswith(END+':'):
                assert target[len(END)+1:] in archived, (name,target);archive_links+=1;continue
            if '/' not in target and '.' not in target: continue
            assert (ROOT/Path(posixpath.normpath(posixpath.join(str(Path(name).parent),target)))).exists(), (name,target)
            local+=1
        for target in re.findall(re.escape(END)+r':([^`\s)]+)',text):
            target=re.sub(r':\d+(?:[–-]\d+)?$','',target)
            if target.startswith('<'): continue
            assert target in archived or any(p.startswith(target.rstrip('/')+'/') for p in archived), (name,target)
            archive_links+=1
    print(f'PASS: inventory partition {len(original & set(kept))} kept + {len(original & pruned)} pruned = {len(original)}.')
    print(f'PASS: {len(kept)} kept paths; {outputs} byte-identical uniform JSON outputs; no relocation-map rows.')
    print(f'PASS: {len(pins)} other Tier1 JSON file/range pins; {exports-len(live_bad)}/{exports} current export pins; 2 compact Lean input hashes.')
    print(f'PASS: {local} local record links and {archive_links} archive references; {len(frozen)} exact frozen method exempt from live-link checks.')
    print(f'KNOWN OPEN (H4): 2 export/producer mismatches introduced by 56595cf7 and 035bb654: {live_bad}')
    if live_exports and live_bad: raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--check-selection", action="store_true", help="Check Phase 3 path/hash metadata only")
    parser.add_argument("--check-live-exports", action="store_true", help="Also fail on live export/producer drift")
    args = parser.parse_args()
    if args.check_selection or args.check_live_exports:
        check_selection(args.check_live_exports)
        raise SystemExit(0)
    NAMES = {Path(p).name.lower() for p in git("ls-tree", "-r", "--name-only", END).decode().splitlines()}
    rows, annex = inventory()
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(("path", "A/M", "bytes", "workstream", "kind"))
    writer.writerows(rows)
    data = buffer.getvalue().encode()
    if args.check:
        assert DEST.read_bytes() == data, "INVENTORY.tsv differs from regenerated metadata"
    else:
        DEST.write_bytes(data)
    print(f"{len(rows)} files: 4434 A, 24 M; {annex} annex sizes read from keys; no content fetched.")
    print(f"Frozen endpoint: {git('rev-parse', END + '^{commit}').decode().strip()}")
    counts = collections.Counter((r[3], r[4]) for r in rows)
    print("| Workstream | " + " | ".join(f"K{i+1}" for i in range(len(KINDS))) + " | Total |")
    print("| --- | " + " | ".join("---:" for _ in range(len(KINDS)+1)) + " |")
    for stream in STREAMS.values():
        values = [counts[stream, k] for k in KINDS]
        print("| " + stream + " | " + " | ".join(map(str, values)) + f" | {sum(values)} |")
    print("| Total | " + " | ".join(str(sum(counts[s, k] for s in STREAMS.values())) for k in KINDS) + f" | {len(rows)} |")
