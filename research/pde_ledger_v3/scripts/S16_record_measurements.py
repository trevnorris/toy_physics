#!/usr/bin/env python3
"""Retrieve S16 record evidence; never import or run an engine/comparator.

Only text retrieval, JSON field selection, literal-match counts and stored-object
shape/cardinality are used. No expression is parsed as mathematics or evaluated.
Run from any directory. Working output stays in _scratch/s16/rec/work/.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import contextlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[3]
BASE = "4d4d22cb"
V3 = "research/pde_ledger_v3/"
SCRIPT = V3 + "scripts/S16_record_measurements.py"
CMP = V3 + "scripts/out/S16_cross_engine_comparator.out"
MEASUREMENTS = V3 + "steps/_measurements/S16_record_measurements.md"
ANCILLARY = [V3 + "SUBSTRATE_REQUIREMENTS.md", V3 + "steps/S14a_2_slow_limit.md", "STATUS.md"]


def emit(value):
    print(json.dumps(value, ensure_ascii=False, sort_keys=True))


def command(argv):
    print("$ " + " ".join(argv))
    result = subprocess.run(argv, cwd=ROOT, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, check=True)
    print(result.stdout, end="" if result.stdout.endswith("\n") else "\n")


def lines(path, start=1, end=None):
    print(f"FILE {path} lines {start}:{end or 'EOF'}")
    for n, line in enumerate((ROOT / path).read_text().splitlines(), 1):
        if n >= start and (end is None or n <= end):
            print(f"{n}: {line}")


def records():
    with (ROOT / CMP).open() as handle:
        for n, line in enumerate(handle, 1):
            yield n, json.loads(line)


def source_set(names):
    for name in names:
        lines(V3 + "directives/_measurements/" + name + ".md")


def m01():
    command(["git", "log", "-1", "--format=%H %s", BASE])
    for commit in ("0464b80a", "1d990ee6", "f08db34a", "9d069620", "79511272"):
        command(["git", "log", "-1", "--format=%H %s", commit])
    lines(V3 + "directives/S16_record_decision_list.md")
    lines(V3 + "directives/S16_decision_list.md")
    for path, commit in (("directives/S16_decision_list.md", "0464b80a"),
                         ("scripts/S16_moving_region.py", "1d990ee6"),
                         ("mathematica/S16_moving_region.wl", "1d990ee6"),
                         ("scripts/S16_cross_engine_comparator.py", "f08db34a")):
        command(["git", "diff", "--exit-code", commit, "--", V3 + path])


def m02():
    source_set([f"S16_decision_list_review_p{p}_disposition" for p in (1, 2, 3, 4)])
    source_set(["S16_build_review_r0_disposition", "S16_build_review_r1_disposition",
                "S16_build_run_record"])


def m03():
    lines(V3 + "directives/S16_comparator_decision_list.md")
    source_set([f"S16_comparator_review_r{r}_disposition" for r in (0, 1, 2)])
    source_set(["S16_comparator_run_record", "S16_comparator_builder_report",
                "S16_comparator_production_run"])


def m04():
    for path, spans in [
        (V3 + "V3_STEP_PLAN.md", [(670, 746), (766, 788), (1296, 1312)]),
        ("research/4d_1pn_full/paper/4d_1pn_full.tex", [(695, 720), (761, 765), (885, 892)]),
        ("research/4d_2_5pn/paper/4d_2_5pn.tex", [(605, 614)]),
        (V3 + "CHARTER.md", [(74, 80), (109, 119)]),
        (V3 + "DEFECT_REGISTER.md", [(166, 172)]),
        (V3 + "directives/O2_SHARED_PHYSICS.md", [(125, 140), (607, 607)]),
        (V3 + "directives/O2_steady_brane_balance_scoping.md", [(164, 164)]),
        ("notes/muonium_gravity_corrections.md", [(21, 25)]),
        (V3 + "steps/S14a_1_projection.md", [(35, 103)]),
        (V3 + "steps/S12_order_source.md", [(1, 24), (280, 303)]),
    ]:
        for start, end in spans:
            lines(path, start, end)


def m05():
    lines("docs/medium_requirements_and_prior_art.md", 370, 474)
    lines("docs/_measurements/gravity_prior_art_sources.md")
    for r in (0, 1, 2):
        lines(f"docs/_measurements/gravity_prior_art_review_r{r}_disposition.md")


def m06():
    lines(V3 + "scripts/S16_moving_region.py")
    lines(V3 + "mathematica/S16_moving_region.wl")


def m07():
    # Whole native payloads, except the separately accounted native sum audits.
    for path in (V3 + "scripts/out/S16_moving_region.out",
                 V3 + "mathematica/out/S16_moving_region.out"):
        print("FILE", path, "literal lines excluding tags containing .dimensionSum.")
        for n, line in enumerate((ROOT / path).read_text().splitlines(), 1):
            if ".dimensionSum." not in line.split(" = ", 1)[0]:
                print(f"{n}: {line}")
        print("LITERAL Missing[\"zeroTerm\"] LINE COUNT",
              sum('Missing["zeroTerm"]' in s for s in (ROOT / path).read_text().splitlines()))


def m08():
    for n, row in records():
        if row["kind"] in ("configuration", "summary", "resources"):
            print("JSONL_LINE", n)
            emit(row)


def m09():
    # Every comparison's stored residual fields, with its own binding text.
    bindings = {}
    for n, row in records():
        if row["kind"] == "configuration":
            bindings = {b["key"]: b for b in row["binding_table"]}
        if row["kind"] == "comparison":
            fields = {k: v for k, v in row.items() if k not in (
                "samples", "py_value", "wl_value", "py_alias_resolved", "wl_alias_resolved")}
            emit({"jsonl_line": n, "binding_object": bindings[row["binding"]]["object"],
                  "stored_fields": fields})


def m10():
    # Exact literal frequencies only: strings are never interpreted as values.
    kinds = Counter()
    fields = Counter()
    by_rule = Counter()
    by_rule_residual = Counter()
    exact = defaultdict(list)
    conditions = defaultdict(list)
    texts = defaultdict(list)
    dimensions = defaultdict(list)
    samples = Counter()
    for n, row in records():
        kinds[row["kind"]] += 1
        if row["kind"] != "comparison":
            continue
        fields.update(row.keys())
        by_rule[row["binding"].split(":", 1)[0]] += 1
        if "exact_residual" in row:
            pair = ("literal 0" if row["exact_residual"] == "0" else "other literal",
                    "literal 0" if row["resolved_exact_residual"] == "0" else "other literal")
            exact[pair].append(row["binding"])
            by_rule_residual[(row["binding"].split(":", 1)[0], *pair)] += 1
            conditions[row["condition_residual"]].append(row["binding"])
        if "text_difference" in row:
            texts[json.dumps(row["text_difference"], ensure_ascii=False, sort_keys=True)].append(row["binding"])
        if "dimension_difference" in row:
            dimensions[row["dimension_difference"]].append(row["binding"])
        for sample in row.get("samples", []):
            samples[(sample.get("coverage"), sample.get("residual"))] += 1
    emit({"row_class_literal_counts": dict(kinds), "comparison_field_counts": dict(fields),
          "comparison_rule_literal_counts": dict(by_rule)})
    for (rule, exact_literal, resolved_literal), count in by_rule_residual.items():
        emit({"binding_rule": rule, "exact_residual_literal_category": exact_literal,
              "resolved_residual_literal_category": resolved_literal, "count": count})
    for title, groups in (("exact/resolved literal category", exact),
                          ("condition_residual literal", conditions),
                          ("text_difference literal", texts),
                          ("dimension_difference literal", dimensions)):
        for literal, members in groups.items():
            emit({"group": title, "literal": literal, "count": len(members), "bindings": members})
    for (coverage, residual), count in samples.items():
        emit({"sample_coverage_literal": coverage, "sample_residual_literal": residual, "count": count})


def m11():
    # Native sample rows and their exact operands for domain comparisons, and
    # every exclusion/unresolved sample. Other samples get literal frequencies
    # PER BINDING, so no binding's zero samples are borrowed from another.
    for n, row in records():
        if row["kind"] != "comparison" or "samples" not in row:
            continue
        if row["binding"] in ("domains:0", "domains:1"):
            print("JSONL_LINE", n)
            emit(row)
        groups = Counter((s.get("coverage"), s.get("residual"), s.get("condition_residual"))
                         for s in row["samples"])
        emit({"binding": row["binding"], "literal_sample_frequencies": [
            {"coverage": k[0], "residual": k[1], "condition_residual": k[2], "count": v}
            for k, v in groups.items()]})
        for sample in row["samples"]:
            if sample.get("coverage") != "evaluated":
                emit({"binding": row["binding"], "stored_sample": sample})
    # The complete test functions are already stored; retrieval makes no new
    # profiles, quadratures, substitutions or derivative checks.
    for n, row in records():
        if row["kind"] == "test_functions":
            print("JSONL_LINE", n)
            emit(row)


def m12():
    # Exhaustive membership by exact stored reason; no mathematical grouping.
    unjoined = defaultdict(list)
    unpaired = defaultdict(list)
    counts = Counter()
    for n, row in records():
        if row["kind"] != "tag_accounting":
            continue
        counts[(row["engine"], row["accounting"])] += 1
        if row["accounting"] == "unjoined":
            unjoined[(row["engine"], row["reason"])].append(
                {k: row[k] for k in ("tag", "line", "occurrence", "parsed_leaves", "actually_compared_leaves")})
        else:
            for component in row["unpaired_components"]:
                unpaired[(row["engine"], component["reason"])].append(
                    {"tag": row["tag"], "line": row["line"], "component": component})
    for key, count in counts.items():
        emit({"engine": key[0], "accounting": key[1], "tag_count": count})
    for title, groups in (("unjoined tag", unjoined), ("unpaired bound component", unpaired)):
        for (engine, reason), members in groups.items():
            emit({"group": title, "engine": engine, "literal_reason": reason, "count": len(members)})
            for member in members:
                emit(member)


def m13():
    # All binding counts, including every below-count binding named by summary.
    for n, row in records():
        if row["kind"] == "binding_accounting":
            print("JSONL_LINE", n)
            emit(row)


def m14():
    for version in (BASE, "4ae7e6cf", "fa1fa413", None):
        path = V3 + "SUBSTRATE_REQUIREMENTS.md"
        text = (subprocess.run(["git", "show", f"{version}:{path}"], cwd=ROOT,
                               text=True, capture_output=True, check=True).stdout
                if version else (ROOT / path).read_text())
        print("REGISTER", version or "working tree")
        print("ENTRY_HEADING_LITERAL_COUNT", sum(s.startswith("### R-") for s in text.splitlines()))
        for n, line in enumerate(text.splitlines(), 1):
            if line.startswith("### R-") or "**target**" in line or "**status**" in line:
                print(f"{n}: {line}")
        if version:
            print(text)
    lines(path)


def m15():
    source_set(["S14a_2_record_r2_review_disposition"])
    for path in (V3 + "steps/S14a_2_slow_limit.md", "STATUS.md"):
        lines(path)


def m16():
    command(["git", "diff", "--no-ext-diff", "--no-color", BASE, "--", *ANCILLARY])
    command(["git", "diff", "--name-only", BASE])
    command(["git", "status", "--short", "--", V3 + "steps/S16_moving_region.md",
             MEASUREMENTS, SCRIPT, *ANCILLARY, ".claude/skills/prep-for-compact/SKILL.md"])


def m17():
    command(["git", "log", "-1", "--format=%H %s", "4ae7e6cf"])
    source_set(["S16_record_r0_review_disposition", "S16_record_r0_review_lookups"])
    command(["git", "show", "--format=", "--name-only", "4ae7e6cf", "--",
             V3 + "directives/_measurements/S16_record_r0_review_disposition.md",
             V3 + "directives/_measurements/S16_record_r0_review_lookups.md"])


def m18():
    # Y4: select only the seven literal rule names; retrieve their stored
    # operands, domains, residuals and samples without evaluating an expression.
    rules = {"supplied_material", "supplied_momentum", "integral_definitions",
             "geometry_surfaces", "geometry_normal", "geometry_areaMeasure",
             "geometry_volumeMeasure"}
    counts = Counter()
    for n, row in records():
        if row["kind"] == "configuration":
            for binding in row["binding_table"]:
                if binding["rule"] in rules:
                    emit({"configuration_line": n, "stored_binding": binding})
        if row["kind"] == "comparison":
            rule = row["binding"].split(":", 1)[0]
            if rule in rules:
                counts[rule] += 1
                print("JSONL_LINE", n)
                emit(row)
    emit({"literal_rule_counts": dict(counts)})


def m19():
    command(["git", "log", "-1", "--format=%H %s", "fa1fa413"])
    command(["git", "log", "-1", "--format=%H %s", "4e5a94bc"])
    source_set(["S16_record_r1_review_disposition", "S16_record_r1_review_lookups"])
    paths = [V3 + "directives/_measurements/S16_record_r1_review_disposition.md",
             V3 + "directives/_measurements/S16_record_r1_review_lookups.md"]
    command(["git", "show", "--format=", "--name-only", "fa1fa413", "--", *paths])
    command(["git", "show", "--format=", "--name-only", "4e5a94bc", "--", *paths])


SECTIONS = {
    "M01": ("Contract, decisions and accepted versions", m01),
    "M02": ("Step-list and build dispositions, with build run record", m02),
    "M03": ("Comparator contract, acceptance and provenance", m03),
    "M04": ("Separate source-row inventories and inherited interfaces", m04),
    "M05": ("Filed gravity prior art, conditions, access limits and reviews", m05),
    "M06": ("Accepted construction source lines", m06),
    "M07": ("Native printed objects and per-object dimensions", m07),
    "M08": ("Comparator configuration, summary and resources verbatim as JSON", m08),
    "M09": ("Every binding's own object and stored comparison fields", m09),
    "M10": ("Literal counts and exhaustive binding membership", m10),
    "M11": ("Per-binding sample literals, domain probes, exclusions and live test functions", m11),
    "M12": ("Every unjoined tag and unpaired bound component, by stored reason", m12),
    "M13": ("Every binding's stored parsed/compared counts", m13),
    "M14": ("Register before/after, schema, rest-on method and both routes", m14),
    "M15": ("S14a part 2 acceptance and current provenance/status text", m15),
    "M16": ("Exact ancillary diff against the starting commit and edit inventory", m16),
    "M17": ("Preserved round-0 record review and repair scope Y1–Y4", m17),
    "M18": ("Supplied-equation, integral-definition and geometry comparisons", m18),
    "M19": ("Round-1 W1 dispositions and the user's dated admission decision", m19),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lookup", choices=SECTIONS)
    args = parser.parse_args()
    if args.lookup:
        SECTIONS[args.lookup][1]()
        return
    work = ROOT / "_scratch/s16/rec/work"
    work.mkdir(parents=True, exist_ok=True)
    temporary = work / "S16_record_measurements.md"
    with temporary.open("w") as output:
        output.write("# S16 record measurements\n\n")
        output.write(f"Generator: `{SCRIPT}`. Baseline: `{BASE}`. Run `python3 {SCRIPT}` to regenerate.\n\n")
        output.write("Only retrieval, literal-match counts and named stored-object shape/cardinality. "
                     "No engine, comparator, CAS, numeric evaluation or new comparison is run. "
                     "JSON objects are selected/re-serialized without interpreting their mathematical strings. "
                     "Every block below is the named command and its complete literal stdout. "
                     "Line numbers in JSONL refer to the filed comparator output.\n\n")
        for key, (title, action) in SECTIONS.items():
            output.write(f"## {key}. {title}\n\n````text\n$ python3 {SCRIPT} --lookup {key}\n")
            with contextlib.redirect_stdout(output):
                action()
            output.write("````\n\n")
    (ROOT / MEASUREMENTS).write_bytes(temporary.read_bytes())


if __name__ == "__main__":
    main()
