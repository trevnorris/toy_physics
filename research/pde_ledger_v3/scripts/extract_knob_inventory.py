#!/usr/bin/env python3
"""Report classified SymPy declaration sites from explicitly named engines."""

import argparse
import ast
import hashlib
import re
from pathlib import Path

import sympy as sp
from sympy.assumptions.assume import AppliedPredicate
from sympy.core.symbol import Str


CLASS_TAGS = ("KNOB", "STRUCTURAL", "COORDINATE", "CONTROL", "PREMISE", "DERIVED")
ANNOTATION = re.compile(
    rf"# ({'|'.join(CLASS_TAGS)}) · (\S(?:[^#\n]*\S)?)$"
)
CONSTRUCTORS = {"Symbol", "symbols", "Function"}


def file_sha256(path):
    """Return the digest of one generator input exactly as read from disk."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def assumption_channel_operands_and_residual(ledger, assumption_objects):
    """Compare exported Symbol assumptions with unary Q predicates used by an engine."""
    predicates = {}
    for assumption_object in assumption_objects:
        for predicate in assumption_object.atoms(AppliedPredicate):
            if len(predicate.arguments) != 1:
                continue
            symbol = predicate.arguments[0]
            if not isinstance(symbol, sp.Symbol):
                continue
            predicates[(symbol.name, sp.srepr(predicate))] = predicate

    rows = []
    residual = sp.Integer(0)
    for (symbol_name, _), predicate in sorted(predicates.items()):
        binding = ledger.get(symbol_name, {}).get("value")
        if isinstance(binding, sp.Symbol):
            inherited_predicate = predicate.function(binding)
            inherited_result = sp.ask(inherited_predicate)
        else:
            inherited_predicate = Str("missing_symbol_binding")
            inherited_result = None
        result_object = (
            sp.true if inherited_result is True
            else sp.false if inherited_result is False
            else Str("unknown")
        )
        row_residual = sp.Integer(inherited_result is not True)
        rows.append(
            sp.Tuple(
                Str(symbol_name),
                Str(sp.srepr(binding)),
                predicate,
                inherited_predicate,
                result_object,
                row_residual,
            )
        )
        residual += row_residual
    return sp.Tuple(*rows), residual


class DeclarationVisitor(ast.NodeVisitor):
    def __init__(self):
        self.constructor_lines = set()
        self.assignment_names = {}
        self.all_assignment_names = {}
        self.module_scope = True

    @staticmethod
    def _target_names(target):
        if isinstance(target, ast.Name):
            return [target.id]
        if isinstance(target, (ast.Tuple, ast.List)):
            return [
                name
                for element in target.elts
                for name in DeclarationVisitor._target_names(element)
            ]
        return []

    def visit_Assign(self, node):
        names = [
            name
            for target in node.targets
            for name in self._target_names(target)
        ]
        self.all_assignment_names[node.lineno] = names
        if self.module_scope:
            self.assignment_names[node.lineno] = names
        self.generic_visit(node)

    def visit_AnnAssign(self, node):
        names = self._target_names(node.target)
        self.all_assignment_names[node.lineno] = names
        if self.module_scope:
            self.assignment_names[node.lineno] = names
        self.generic_visit(node)

    def visit_AugAssign(self, node):
        names = self._target_names(node.target)
        self.all_assignment_names[node.lineno] = names
        if self.module_scope:
            self.assignment_names[node.lineno] = names
        self.generic_visit(node)

    def _visit_nested_scope(self, node):
        was_module_scope = self.module_scope
        self.module_scope = False
        self.generic_visit(node)
        self.module_scope = was_module_scope

    visit_FunctionDef = _visit_nested_scope
    visit_AsyncFunctionDef = _visit_nested_scope
    visit_ClassDef = _visit_nested_scope
    visit_Lambda = _visit_nested_scope

    def visit_Call(self, node):
        function = node.func
        if (
            isinstance(function, ast.Attribute)
            and isinstance(function.value, ast.Name)
            and function.value.id == "sp"
            and function.attr in CONSTRUCTORS
        ):
            self.constructor_lines.add(node.lineno)
        self.generic_visit(node)


def scan(path):
    source = path.read_text(encoding="utf-8")
    source_lines = source.splitlines()
    tree = ast.parse(source, filename=str(path))
    visitor = DeclarationVisitor()
    visitor.visit(tree)
    records = []
    findings = []
    declaration_lines = visitor.constructor_lines | set(visitor.assignment_names)
    for line_number in sorted(declaration_lines):
        line = source_lines[line_number - 1]
        match = ANNOTATION.search(line)
        declaration = line.split("#", 1)[0].strip()
        if match is None:
            findings.append(f"{path}:{line_number}: declaration has no valid class tag: {declaration}")
            continue
        records.append(
            (
                match.group(1),
                line_number,
                declaration,
                match.group(2),
                tuple(visitor.assignment_names.get(line_number, ())),
            )
        )
    return records, findings


def declaration_classes(path):
    _, findings = scan(path)
    if findings:
        raise ValueError("\n".join(findings))
    source = path.read_text(encoding="utf-8")
    source_lines = source.splitlines()
    tree = ast.parse(source, filename=str(path))
    visitor = DeclarationVisitor()
    visitor.visit(tree)
    classes = {}
    annotated_assignments = [
        (match.group(1), line_number, names)
        for line_number, names in visitor.all_assignment_names.items()
        if (match := ANNOTATION.search(source_lines[line_number - 1])) is not None
    ]
    for class_tag, line_number, names in annotated_assignments:
        for name in names:
            if name in classes and classes[name] != class_tag:
                raise ValueError(f"{path}:{line_number}: conflicting class tags for {name}")
            classes[name] = class_tag
    return classes


def main():
    parser = argparse.ArgumentParser(
        description="Extract annotated declaration sites from an explicit engine-file list."
    )
    parser.add_argument("engines", nargs="+", type=Path)
    arguments = parser.parse_args()

    print("INVENTORY_UNIT: annotated declaration site; a class-uniform grouped declaration is one site")
    print(
        "SOURCE_LIMIT: runtime symbols built from dynamic expressions cannot be enumerated from source; "
        "the annotated coefficient-dimension Symbol declaration site is reported"
    )
    all_findings = []
    for path in arguments.engines:
        records, findings = scan(path)
        all_findings.extend(findings)
        print(f"FILE: {path}")
        for class_tag in CLASS_TAGS:
            entries = [
                f"line {line_number}: {declaration} -- {description}"
                for record_class, line_number, declaration, description, _ in records
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
