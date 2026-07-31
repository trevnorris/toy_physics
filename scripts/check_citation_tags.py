#!/usr/bin/env python3
"""Verify every `FILE.md#tag` citation resolves to a `{#tag}` that exists.

Line-number citations in this repo go stale silently whenever a cited file
changes length -- six measured instances in one session, every gate green
throughout.  Tags do not move, but they can still be renamed or deleted, and
that is what this checks.

Exit 0 = every citation resolves and every tag is unique within its file.
Exit 1 = at least one dangling citation, duplicate tag, or unresolvable path.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# `{#tag}` as it appears where a tag is DEFINED.
TAG_DEF = re.compile(r"\{#([A-Za-z0-9_-]+)\}")
# `some/path/FILE.md#tag` as it appears where a tag is CITED.
TAG_CITE = re.compile(r"([A-Za-z0-9_/.-]+\.md)#([A-Za-z0-9_-]+)")

SKIP_DIRS = {".git", "archive", "__pycache__", "node_modules"}


def markdown_files() -> list[Path]:
    return [
        p
        for p in REPO.rglob("*.md")
        if not any(part in SKIP_DIRS for part in p.relative_to(REPO).parts)
    ]


def resolve(cited_path: str, citing_file: Path) -> Path | None:
    """Resolve a cited path repo-relative, then citing-file-relative, then by
    unique basename.  Ambiguous basenames resolve to None -- an ambiguous
    citation is a real defect, not something to guess at."""
    candidate = REPO / cited_path
    if candidate.is_file():
        return candidate
    candidate = citing_file.parent / cited_path
    if candidate.is_file():
        return candidate
    matches = [p for p in markdown_files() if p.name == Path(cited_path).name]
    return matches[0] if len(matches) == 1 else None


def main() -> int:
    files = markdown_files()

    # Where each tag is defined, and how many times.
    definitions: dict[Path, dict[str, int]] = {}
    for path in files:
        counts: dict[str, int] = {}
        for tag in TAG_DEF.findall(path.read_text(encoding="utf-8", errors="replace")):
            counts[tag] = counts.get(tag, 0) + 1
        if counts:
            definitions[path] = counts

    failures: list[str] = []

    for path, counts in definitions.items():
        for tag, n in counts.items():
            if n > 1:
                failures.append(
                    f"DUPLICATE TAG: {{#{tag}}} defined {n} times in "
                    f"{path.relative_to(REPO)} -- a citation to it is ambiguous"
                )

    citations = 0
    for path in files:
        text = path.read_text(encoding="utf-8", errors="replace")
        for cited_path, tag in TAG_CITE.findall(text):
            citations += 1
            target = resolve(cited_path, path)
            if target is None:
                failures.append(
                    f"UNRESOLVED PATH: {path.relative_to(REPO)} cites "
                    f"{cited_path}#{tag}, but that path is missing or ambiguous"
                )
                continue
            if tag not in definitions.get(target, {}):
                failures.append(
                    f"DANGLING CITATION: {path.relative_to(REPO)} cites "
                    f"{cited_path}#{tag}, but {{#{tag}}} is not defined in "
                    f"{target.relative_to(REPO)}"
                )

    defined = sum(len(c) for c in definitions.values())
    print(f"TAGS_DEFINED: {defined}   CITATIONS_CHECKED: {citations}")

    if failures:
        for line in failures:
            print(f"  {line}")
        print(f"CITATION_TAG_CHECK: FAIL ({len(failures)} problem(s))")
        return 1

    print("CITATION_TAG_CHECK: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
