#!/usr/bin/env python3
"""Report classified SymPy declaration sites from explicitly named engines."""

import argparse
import ast
import re
from pathlib import Path


CLASS_TAGS = ("KNOB", "STRUCTURAL", "COORDINATE", "CONTROL", "PREMISE", "DERIVED")
ANNOTATION = re.compile(
    rf"# ({'|'.join(CLASS_TAGS)}) · ([^#\n]+)$"
)
CONSTRUCTORS = {"Symbol", "symbols", "Function"}


class DeclarationVisitor(ast.NodeVisitor):
    def __init__(self):
        self.lines = set()

    def visit_Call(self, node):
        function = node.func
        if (
            isinstance(function, ast.Attribute)
            and isinstance(function.value, ast.Name)
            and function.value.id == "sp"
            and function.attr in CONSTRUCTORS
        ):
            self.lines.add(node.lineno)
        self.generic_visit(node)


def scan(path):
    source = path.read_text(encoding="utf-8")
    source_lines = source.splitlines()
    tree = ast.parse(source, filename=str(path))
    visitor = DeclarationVisitor()
    visitor.visit(tree)
    records = []
    findings = []
    for line_number in sorted(visitor.lines):
        line = source_lines[line_number - 1]
        match = ANNOTATION.search(line)
        declaration = line.split("#", 1)[0].strip()
        if match is None:
            findings.append(f"{path}:{line_number}: declaration has no valid class tag: {declaration}")
            continue
        records.append((match.group(1), line_number, declaration, match.group(2)))
    return records, findings


def main():
    parser = argparse.ArgumentParser(
        description="Extract annotated declaration sites from an explicit engine-file list."
    )
    parser.add_argument("engines", nargs="+", type=Path)
    arguments = parser.parse_args()

    print("INVENTORY_UNIT: annotated declaration site; a class-uniform grouped declaration is one site")
    print(
        "SOURCE_LIMIT: runtime symbols built from dynamic expressions cannot be enumerated from source; "
        "the annotated sp.Symbol(f\"dim_{coefficient}_{axis}\") declaration site is reported"
    )
    all_findings = []
    for path in arguments.engines:
        records, findings = scan(path)
        all_findings.extend(findings)
        print(f"FILE: {path}")
        for class_tag in CLASS_TAGS:
            entries = [
                f"line {line_number}: {declaration} -- {description}"
                for record_class, line_number, declaration, description in records
                if record_class == class_tag
            ]
            print(f"{class_tag}: {'; '.join(entries) if entries else '(none)'}")
    if all_findings:
        print("FINDINGS:")
        for finding in all_findings:
            print(f"  {finding}")
        raise SystemExit(1)
    print("FINDINGS: none")


if __name__ == "__main__":
    main()
