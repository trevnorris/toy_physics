#!/usr/bin/env python3
"""
Cosmetic 3-digit padding of GENUINE 2-digit `Stage NN` neighbor refs in the
070-089 explicit-branch band notes (user-approved 2026-06-04).

Two read-only classifiers (+ the old-form-vs-new-form / +17 trap) confirmed EVERY
2-digit Stage reference in stage070..089 is a genuine current-canonical neighbor ref
(zero stale survivors, zero non-refs). So padding every 2-digit number that sits inside
a Stage-citation token is safe — nothing to exclude.

Method: match whole Stage-citation tokens (lead + optional compound tail covering
`/`, en-dash, hyphen, ` and ` separators), then zero-pad every standalone 2-digit run
inside the matched token to 3 digits. Idempotent (3-digit runs untouched). Label-only:
only digit characters are inserted; no other byte changes. Prints every change for eyeball.
"""
import os, re, glob

HERE = os.path.dirname(os.path.abspath(__file__))
NOTES = os.path.join(HERE, "..", "notes", "stages")

# A Stage-citation token: lead "Stage"/"Stages" + space-or-hyphen + number,
# then zero or more compound continuations joined by / en-dash hyphen or " and ".
TOKEN = re.compile(
    r"(?<![A-Za-z])Stages?[ \-]\d{1,3}(?:(?:[–\-/]|\s+and\s+)\d{1,3})*"
)
SHORT_NUM = re.compile(r"(?<!\d)\d{1,2}(?!\d)")  # 1- or 2-digit runs; 3-digit runs untouched

def pad_token(tok):
    return SHORT_NUM.sub(lambda m: m.group(0).zfill(3), tok)

total_changes = 0
files_touched = 0
for n in range(70, 90):
    hits = glob.glob(os.path.join(NOTES, f"moving_throat_pde_stage0{n}_*.md"))
    if not hits:
        continue
    p = hits[0]
    content = open(p, encoding="utf-8").read()
    changes = []  # (old_tok, new_tok)

    def repl(m):
        tok = m.group(0)
        new = pad_token(tok)
        if new != tok:
            changes.append((tok, new))
        return new

    new_content = TOKEN.sub(repl, content)
    if new_content != content:
        open(p, "w", encoding="utf-8").write(new_content)
        files_touched += 1
        total_changes += len(changes)
        print(f"\n=== stage0{n} ({os.path.basename(p)}) : {len(changes)} token(s) padded ===")
        for old, new in changes:
            print(f"   {old!r}  ->  {new!r}")

print(f"\n--- files touched: {files_touched}   token-instances padded: {total_changes} ---")
