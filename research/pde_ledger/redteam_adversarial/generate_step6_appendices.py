#!/usr/bin/env python3
"""Emit Step 6 Part A audit appendices from redteam_adversarial YAML.

The script is deterministic: it reads only repository artifacts and never reads
the wall clock.  Re-running it should reproduce byte-identical appendix files.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
import re
import unicodedata

import yaml


REDTEAM_ROOT = Path(__file__).resolve().parent
PDE_ROOT = REDTEAM_ROOT.parent
PAPER_APPENDICES = PDE_ROOT / "paper" / "appendices"

FIT_INDEX_OUT = PAPER_APPENDICES / "fit_insertion_index.tex"
PROVENANCE_OUT = PAPER_APPENDICES / "parameter_provenance_audit.tex"

CONSTRAINT_ORDER = [
    "internal_consistency",
    "free_choice",
    "published_target",
    "no_constraint_record",
]

ASCII_REPLACEMENTS = {
    "\u2018": "'",
    "\u2019": "'",
    "\u201c": '"',
    "\u201d": '"',
    "\u2013": "-",
    "\u2014": "-",
    "\u2212": "-",
    "\u00d7": "x",
    "\u2026": "...",
    "\u03b1": "alpha",
    "\u0391": "Alpha",
    "\u03b2": "beta",
    "\u0392": "Beta",
    "\u03b3": "gamma",
    "\u0393": "Gamma",
    "\u03b4": "delta",
    "\u0394": "Delta",
    "\u03b5": "epsilon",
    "\u03b7": "eta",
    "\u03b8": "theta",
    "\u03ba": "kappa",
    "\u03bb": "lambda",
    "\u039b": "Lambda",
    "\u03bc": "mu",
    "\u03c0": "pi",
    "\u03a0": "Pi",
    "\u03c1": "rho",
    "\u03c3": "sigma",
    "\u03a3": "Sigma",
    "\u03c4": "tau",
    "\u03c7": "chi",
    "\u03a7": "Chi",
    "\u03be": "xi",
    "\u039e": "Xi",
    "\u03c9": "omega",
    "\u03a9": "Omega",
    "\u2113": "ell",
    "\U0001d52f": "r",
    "\u221a": "sqrt",
    "\u2264": "<=",
    "\u2265": ">=",
    "\u2248": "approx",
    "\u2260": "!=",
}


def load_yaml(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def to_ascii(value) -> str:
    text = "" if value is None else str(value)
    for src, dst in ASCII_REPLACEMENTS.items():
        text = text.replace(src, dst)
    text = unicodedata.normalize("NFKD", text)
    text = text.encode("ascii", "ignore").decode("ascii")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def tex_escape(value) -> str:
    text = to_ascii(value)
    text = text.replace("\\", r"\textbackslash{}")
    replacements = {
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    for src, dst in replacements.items():
        text = text.replace(src, dst)
    return text


def tex_path(value) -> str:
    text = to_ascii(value)
    text = text.replace("{", "").replace("}", "")
    return r"\AuditPath{" + text + "}"


def tex_id(value) -> str:
    text = to_ascii(value)
    text = text.replace("{", "").replace("}", "")
    return r"\AuditID{" + text + "}"


def tex_id_list(values) -> str:
    values = [value for value in values if value is not None]
    if not values:
        return tex_escape("none recorded")
    return ", ".join(tex_id(value) for value in values)


def citation_ref(citation: dict | None) -> str:
    if not citation:
        return tex_escape("not recorded")
    path = citation.get("path", "not recorded")
    line = citation.get("line")
    if line is None:
        return tex_path(path)
    return tex_path(path) + ":" + tex_escape(line)


def first_citation(citations: list[dict] | None) -> tuple[dict | None, int]:
    citations = citations or []
    return (citations[0] if citations else None, len(citations))


def latest_status_note(candidate: dict) -> str | None:
    notes = candidate.get("status_notes") or []
    if not notes:
        return None
    return notes[-1].get("note")


def manifest_provenance_path(path_text: str) -> Path:
    rel = path_text
    prefix = "redteam_adversarial/"
    if rel.startswith(prefix):
        rel = rel[len(prefix) :]
    return REDTEAM_ROOT / rel


def stage_sort_key(stage_values) -> tuple[int, str]:
    if not stage_values:
        return (9999, "")
    stage = str(stage_values[0])
    try:
        return (int(stage), stage)
    except ValueError:
        return (9999, stage)


def split_label(kinds: list[str]) -> str:
    if not kinds:
        return "no_constraint_record"
    return "+".join(kinds)


def sorted_constraint_labels(counter: Counter) -> list[str]:
    labels = list(counter)
    labels.sort(key=lambda label: ([CONSTRAINT_ORDER.index(part) if part in CONSTRAINT_ORDER else 99 for part in label.split("+")], label))
    return labels


def make_table_row(cells: list[str]) -> str:
    return " & ".join(cells) + r" \\"


def load_sources():
    fit = load_yaml(REDTEAM_ROOT / "fit_insertion_points.yaml")
    manifest = load_yaml(REDTEAM_ROOT / "MANIFEST.yaml")
    benchmarks = load_yaml(REDTEAM_ROOT / "benchmarks.yaml")
    family_map = load_yaml(REDTEAM_ROOT / "provenance" / "_family_map.yaml")
    dedup = load_yaml(REDTEAM_ROOT / "provenance" / "_dedup_proposal.yaml")

    alias_to_canonical: dict[str, str] = {}
    for group in dedup.get("alias_groups", []) or []:
        canonical_id = group.get("canonical_id")
        for alias in group.get("aliases", []) or []:
            alias_to_canonical[alias] = canonical_id

    benchmark_by_family: dict[str, list[dict]] = defaultdict(list)
    for entry in benchmarks.get("entries", []) or []:
        benchmark_by_family[entry.get("family_id")].append(entry)

    provenance_records: list[dict] = []
    provenance_by_candidate: dict[str, list[dict]] = defaultdict(list)
    for path in sorted((REDTEAM_ROOT / "provenance").glob("*.yaml")):
        if path.name.startswith("_"):
            continue
        data = load_yaml(path)
        record = {"path": path, "data": data}
        provenance_records.append(record)
        provenance_by_candidate[data.get("candidate_id")].append(record)

    manifest_candidates = manifest.get("candidates", {}) or {}
    fit_candidates = fit.get("candidates", []) or []
    fit_ids = {candidate["id"] for candidate in fit_candidates}
    manifest_ids = set(manifest_candidates)
    if fit_ids != manifest_ids:
        missing_manifest = sorted(fit_ids - manifest_ids)
        missing_fit = sorted(manifest_ids - fit_ids)
        raise SystemExit(
            "fit/MANIFEST candidate mismatch: "
            f"missing_manifest={missing_manifest[:10]} missing_fit={missing_fit[:10]}"
        )

    referenced_paths: set[Path] = set()
    for candidate_id, candidate in manifest_candidates.items():
        for rel in ((candidate.get("paths") or {}).get("provenance") or []):
            path = manifest_provenance_path(rel)
            if not path.exists():
                raise SystemExit(f"missing provenance path for {candidate_id}: {rel}")
            referenced_paths.add(path.resolve())

    unreferenced = [
        record["path"].name
        for record in provenance_records
        if record["path"].resolve() not in referenced_paths
    ]
    if unreferenced:
        raise SystemExit(f"top-level provenance records not referenced by MANIFEST: {unreferenced[:10]}")

    return {
        "fit": fit,
        "manifest": manifest,
        "benchmarks": benchmarks,
        "family_map": family_map,
        "dedup": dedup,
        "alias_to_canonical": alias_to_canonical,
        "benchmark_by_family": benchmark_by_family,
        "provenance_records": provenance_records,
        "provenance_by_candidate": provenance_by_candidate,
    }


def record_constraint_kinds(record: dict) -> list[str]:
    constraints = record["data"].get("constraints") or []
    return sorted({item.get("constraint_kind") for item in constraints if item.get("constraint_kind")})


def candidate_constraint_kinds(candidate_id: str, provenance_by_candidate: dict[str, list[dict]], alias_to_canonical: dict[str, str]) -> list[str]:
    records = provenance_by_candidate.get(candidate_id, [])
    kinds = sorted({kind for record in records for kind in record_constraint_kinds(record)})
    if kinds:
        return kinds
    canonical = alias_to_canonical.get(candidate_id)
    if canonical:
        return sorted({kind for record in provenance_by_candidate.get(canonical, []) for kind in record_constraint_kinds(record)})
    return []


def candidate_disposition(candidate_id: str, manifest_candidate: dict, alias_to_canonical: dict[str, str]) -> str:
    status = manifest_candidate.get("status", "unknown")
    note = latest_status_note(manifest_candidate)
    if status == "verdict_logged" and note:
        return tex_escape(f"{status}: {note}")
    if candidate_id in alias_to_canonical:
        return tex_escape(f"{status}; alias folded to {alias_to_canonical[candidate_id]}")
    return tex_escape(status)


def provenance_family_ids(candidate_id: str, family_map: dict) -> tuple[list[str], list[str]]:
    primary = list((family_map.get("primary_candidate_family_map") or {}).get(candidate_id, []) or [])
    all_families = list((family_map.get("candidate_family_map") or {}).get(candidate_id, []) or [])
    overlays = [family for family in all_families if family not in primary]
    return primary, overlays


def basis_for_record(record: dict, family_ids: list[str], benchmark_by_family: dict[str, list[dict]]) -> str:
    data = record["data"]
    constraints = data.get("constraints") or []
    kinds = record_constraint_kinds(record)
    kind = kinds[0] if kinds else "no_constraint_record"
    evidence = constraints[0].get("evidence_citation") if constraints else None

    if kind == "internal_consistency":
        return tex_escape("located derivation step: ") + citation_ref(evidence)
    if kind == "free_choice":
        return tex_escape("posited choice: ") + citation_ref(evidence)
    if kind == "published_target":
        entries = []
        for family_id in family_ids:
            entries.extend(benchmark_by_family.get(family_id, []))
        if entries:
            pieces = []
            seen = set()
            for entry in entries:
                key = entry.get("id")
                if key in seen:
                    continue
                seen.add(key)
                citation = entry.get("source_citation", "benchmark citation not recorded")
                pieces.append(f"{entry.get('id')}: {entry.get('source_type')} - {citation}")
            return tex_escape("benchmark: " + " | ".join(pieces))
        return tex_escape("published_target; no benchmarks.yaml entry joined for listed family; evidence: ") + citation_ref(evidence)
    return tex_escape("no constraints block in provenance record; surfaced as scanner/non-parameter accounting")


def downstream_cell(data: dict) -> str:
    dependents = data.get("downstream_dependents") or []
    if not dependents:
        return tex_escape("none recorded")
    return tex_escape(", ".join(str(item) for item in dependents))


def generated_preamble(generator_name: str) -> list[str]:
    return [
        "% GENERATED FILE -- do not edit by hand.",
        f"% Generated by redteam_adversarial/{generator_name}.",
        "% Deterministic source: redteam_adversarial YAML artifacts; no wall-clock read.",
        r"\providecommand{\AuditPath}[1]{\begingroup\footnotesize\nolinkurl{#1}\endgroup}",
        r"\providecommand{\AuditID}[1]{\begingroup\footnotesize\nolinkurl{#1}\endgroup}",
        "",
    ]


def emit_fit_index(sources: dict) -> tuple[str, Counter, Counter]:
    fit_candidates = list(sources["fit"].get("candidates", []) or [])
    manifest_candidates = sources["manifest"].get("candidates", {}) or {}
    provenance_by_candidate = sources["provenance_by_candidate"]
    alias_to_canonical = sources["alias_to_canonical"]

    status_split = Counter(manifest_candidates[candidate["id"]].get("status", "unknown") for candidate in fit_candidates)
    kind_split = Counter()
    rows = []

    for candidate in fit_candidates:
        candidate_id = candidate["id"]
        manifest_candidate = manifest_candidates[candidate_id]
        kinds = candidate_constraint_kinds(candidate_id, provenance_by_candidate, alias_to_canonical)
        label = split_label(kinds)
        kind_split[label] += 1
        first_anchor, anchor_count = first_citation(manifest_candidate.get("file_line_citations") or candidate.get("file_line_citations"))
        anchor = citation_ref(first_anchor)
        if anchor_count > 1:
            anchor += tex_escape(f" ({anchor_count} manifest anchors)")
        params = candidate.get("parameter_names") or manifest_candidate.get("parameter_names") or []
        parameter_cell = tex_id(candidate_id) + r"\newline " + tex_id_list(params)
        if candidate_id in alias_to_canonical:
            kind_cell = tex_id(label) + r"\newline " + tex_escape("alias of ") + tex_id(alias_to_canonical[candidate_id])
        else:
            kind_cell = tex_id(label)
        rows.append(
            (
                stage_sort_key(candidate.get("anchor_stages")),
                candidate_id,
                make_table_row(
                    [
                        tex_escape(", ".join(str(stage) for stage in candidate.get("anchor_stages", []))),
                        parameter_cell,
                        anchor,
                        kind_cell,
                        candidate_disposition(candidate_id, manifest_candidate, alias_to_canonical),
                    ]
                ),
            )
        )

    rows.sort(key=lambda item: (item[0], item[1]))

    lines = generated_preamble("generate_step6_appendices.py")
    lines.extend(
        [
            r"\chapter{Layer-2 Fit-Insertion Index}",
            r"\label{app:fit-insertion-index}",
            "",
            r"\section{Generated Coverage Summary}",
            "",
            tex_escape(
                "This generated appendix indexes every Phase-A fit-insertion candidate from the layer-2 adversarial audit. "
                f"It covers {len(fit_candidates)} candidates. The audit result is operational: no surviving fatal flaw was found "
                "in this toy-analog ledger; this is not a truth claim about the model."
            ),
            "",
            r"\begin{longtable}{@{}L{0.45\textwidth}L{0.45\textwidth}@{}}",
            r"\caption{Fit-insertion index generated coverage totals.}\label{tab:fit-insertion-generated-totals}\\",
            r"\toprule",
            r"Item & Count \\",
            r"\midrule",
            r"\endfirsthead",
            r"\toprule",
            r"Item & Count \\",
            r"\midrule",
            r"\endhead",
            make_table_row([tex_escape("Phase-A candidates"), tex_escape(len(fit_candidates))]),
        ]
    )
    for label in sorted_constraint_labels(kind_split):
        lines.append(make_table_row([tex_escape(f"constraint_kind: {label}"), tex_escape(kind_split[label])]))
    for status in sorted(status_split):
        lines.append(make_table_row([tex_escape(f"MANIFEST status: {status}"), tex_escape(status_split[status])]))
    lines.extend(
        [
            make_table_row([tex_escape("Alias-folded candidates surfaced in rows"), tex_escape(len(alias_to_canonical))]),
            r"\bottomrule",
            r"\end{longtable}",
            "",
            r"\section{Phase-A Candidate Rows}",
            "",
            r"\begingroup\scriptsize",
            r"\begin{longtable}{@{}L{0.06\textwidth}L{0.27\textwidth}L{0.24\textwidth}L{0.18\textwidth}L{0.17\textwidth}@{}}",
            r"\caption{Every Phase-A candidate with source anchor, constraint-kind accounting, and audit disposition.}\label{tab:fit-insertion-candidates}\\",
            r"\toprule",
            r"Stage & Candidate / parameter & Source anchor & constraint\_kind & Audit disposition \\",
            r"\midrule",
            r"\endfirsthead",
            r"\toprule",
            r"Stage & Candidate / parameter & Source anchor & constraint\_kind & Audit disposition \\",
            r"\midrule",
            r"\endhead",
        ]
    )
    lines.extend(row for _, _, row in rows)
    lines.extend(
        [
            r"\bottomrule",
            r"\end{longtable}",
            r"\endgroup",
            "",
            r"\section{Visible Grouping and Exclusion Accounting}",
            "",
            tex_escape(
                "No Phase-A candidate is silently dropped. Alias-folded scanned candidates remain listed in the main table and are also enumerated here. "
                "Records whose provenance YAML has no constraints block are listed with constraint_kind no_constraint_record rather than being forced into a fit-vs-derive class."
            ),
            "",
            r"\begin{longtable}{@{}L{0.38\textwidth}L{0.52\textwidth}@{}}",
            r"\caption{Alias accounting from the dedup proposal.}\label{tab:fit-insertion-alias-accounting}\\",
            r"\toprule",
            r"Alias candidate & Canonical candidate \\",
            r"\midrule",
            r"\endfirsthead",
            r"\toprule",
            r"Alias candidate & Canonical candidate \\",
            r"\midrule",
            r"\endhead",
        ]
    )
    for alias, canonical in sorted(alias_to_canonical.items()):
        lines.append(make_table_row([tex_id(alias), tex_id(canonical)]))
    lines.extend([r"\bottomrule", r"\end{longtable}", ""])
    return "\n".join(lines) + "\n", kind_split, status_split


def emit_parameter_provenance(sources: dict) -> tuple[str, Counter]:
    provenance_records = list(sources["provenance_records"])
    family_map = sources["family_map"]
    benchmark_by_family = sources["benchmark_by_family"]
    alias_to_canonical = sources["alias_to_canonical"]
    benchmarks = sources["benchmarks"].get("entries", []) or []

    kind_split = Counter(split_label(record_constraint_kinds(record)) for record in provenance_records)

    def row_key(record: dict):
        data = record["data"]
        return (
            stage_sort_key(data.get("anchor_stages")),
            data.get("candidate_id", ""),
            data.get("parameter_name", ""),
            record["path"].name,
        )

    rows = []
    for record in sorted(provenance_records, key=row_key):
        data = record["data"]
        candidate_id = data.get("candidate_id")
        primary_families, overlay_families = provenance_family_ids(candidate_id, family_map)
        family_text = tex_id_list(primary_families) if primary_families else tex_escape("no primary family")
        if overlay_families:
            family_text += tex_escape("; overlays: ") + tex_id_list(overlay_families)
        kinds = split_label(record_constraint_kinds(record))
        origin_claims = data.get("origin_claims") or []
        origin = origin_claims[0] if origin_claims else {}
        origin_stage = origin.get("introduced_at_stage") or ", ".join(str(stage) for stage in data.get("anchor_stages", []))
        origin_citation = citation_ref(origin.get("citation"))
        parameter_cell = (
            tex_id(data.get("parameter_name"))
            + r"\newline "
            + tex_id(candidate_id)
            + r"\newline "
            + family_text
        )
        rows.append(
            make_table_row(
                [
                    tex_escape(", ".join(str(stage) for stage in data.get("anchor_stages", []))),
                    parameter_cell,
                    tex_escape(f"origin stage {origin_stage}") + r"\newline " + origin_citation,
                    tex_id(kinds),
                    basis_for_record(record, primary_families, benchmark_by_family),
                    downstream_cell(data),
                ]
            )
        )

    lines = generated_preamble("generate_step6_appendices.py")
    lines.extend(
        [
            r"\chapter{Layer-2 Parameter-Provenance Audit}",
            r"\label{app:parameter-provenance-audit}",
            "",
            r"\section{Generated Coverage Summary}",
            "",
            tex_escape(
                "This generated appendix exposes the parameter-level genealogy used by the layer-2 fit-vs-derive audit. "
                f"It covers {len(provenance_records)} top-level provenance records. The audit result is operational: no surviving fatal flaw was found "
                "in this toy-analog ledger; this is not a truth claim about the model."
            ),
            "",
            r"\begin{longtable}{@{}L{0.45\textwidth}L{0.45\textwidth}@{}}",
            r"\caption{Parameter-provenance generated coverage totals.}\label{tab:parameter-provenance-generated-totals}\\",
            r"\toprule",
            r"Item & Count \\",
            r"\midrule",
            r"\endfirsthead",
            r"\toprule",
            r"Item & Count \\",
            r"\midrule",
            r"\endhead",
            make_table_row([tex_escape("Top-level provenance parameter records"), tex_escape(len(provenance_records))]),
            make_table_row([tex_escape("Primary canonical families"), tex_escape((family_map.get("summary") or {}).get("family_count", "not recorded"))]),
            make_table_row([tex_escape("Canonical candidates"), tex_escape((family_map.get("summary") or {}).get("canonical_candidate_count", "not recorded"))]),
            make_table_row([tex_escape("Alias candidates excluded from canonical count but surfaced below"), tex_escape((family_map.get("summary") or {}).get("alias_count_excluded", len(alias_to_canonical)))]),
        ]
    )
    for label in sorted_constraint_labels(kind_split):
        lines.append(make_table_row([tex_escape(f"constraint_kind: {label}"), tex_escape(kind_split[label])]))
    lines.extend([r"\bottomrule", r"\end{longtable}", ""])

    lines.extend(
        [
            r"\section{Published-Target Benchmark Sources}",
            "",
            tex_escape(
                "The rows below are the benchmark citations loaded from benchmarks.yaml. Published-target provenance rows cite these entries where their family map joins to one of these family identifiers; if no join exists, the row states that absence explicitly."
            ),
            "",
            r"\begingroup\scriptsize",
            r"\begin{longtable}{@{}L{0.16\textwidth}L{0.14\textwidth}L{0.62\textwidth}@{}}",
            r"\caption{Benchmark citations available to published-target provenance rows.}\label{tab:parameter-provenance-benchmarks}\\",
            r"\toprule",
            r"Family & Type & Citation \\",
            r"\midrule",
            r"\endfirsthead",
            r"\toprule",
            r"Family & Type & Citation \\",
            r"\midrule",
            r"\endhead",
        ]
    )
    for entry in sorted(benchmarks, key=lambda item: (item.get("family_id", ""), item.get("id", ""))):
        lines.append(
            make_table_row(
                [
                    tex_id(entry.get("family_id")),
                    tex_escape(entry.get("source_type")),
                    tex_escape(f"{entry.get('id')}: {entry.get('source_citation')}"),
                ]
            )
        )
    lines.extend([r"\bottomrule", r"\end{longtable}", r"\endgroup", ""])

    lines.extend(
        [
            r"\section{Parameter Genealogy Rows}",
            "",
            r"\begingroup\scriptsize",
            r"\begin{longtable}{@{}L{0.05\textwidth}L{0.24\textwidth}L{0.16\textwidth}L{0.13\textwidth}L{0.28\textwidth}L{0.08\textwidth}@{}}",
            r"\caption{Parameter-level provenance records with basis and downstream dependents.}\label{tab:parameter-provenance-records}\\",
            r"\toprule",
            r"Stage & Parameter / candidate / family & Origin & constraint\_kind & Basis & Downstream \\",
            r"\midrule",
            r"\endfirsthead",
            r"\toprule",
            r"Stage & Parameter / candidate / family & Origin & constraint\_kind & Basis & Downstream \\",
            r"\midrule",
            r"\endhead",
        ]
    )
    lines.extend(rows)
    lines.extend([r"\bottomrule", r"\end{longtable}", r"\endgroup", ""])

    lines.extend(
        [
            r"\section{Visible Alias Accounting}",
            "",
            tex_escape(
                "The family map canonicalizes 915 candidates and excludes 7 alias candidates from the canonical count. Those aliases are visible here so the appendix does not silently drop them."
            ),
            "",
            r"\begin{longtable}{@{}L{0.38\textwidth}L{0.52\textwidth}@{}}",
            r"\caption{Alias accounting used by the parameter-provenance family map.}\label{tab:parameter-provenance-alias-accounting}\\",
            r"\toprule",
            r"Alias candidate & Canonical candidate \\",
            r"\midrule",
            r"\endfirsthead",
            r"\toprule",
            r"Alias candidate & Canonical candidate \\",
            r"\midrule",
            r"\endhead",
        ]
    )
    for alias, canonical in sorted(alias_to_canonical.items()):
        lines.append(make_table_row([tex_id(alias), tex_id(canonical)]))
    lines.extend([r"\bottomrule", r"\end{longtable}", ""])
    return "\n".join(lines) + "\n", kind_split


def main() -> None:
    sources = load_sources()
    fit_tex, fit_kind_split, fit_status_split = emit_fit_index(sources)
    provenance_tex, provenance_kind_split = emit_parameter_provenance(sources)

    FIT_INDEX_OUT.write_text(fit_tex, encoding="utf-8")
    PROVENANCE_OUT.write_text(provenance_tex, encoding="utf-8")

    print(f"wrote {FIT_INDEX_OUT.relative_to(PDE_ROOT)}")
    print(f"wrote {PROVENANCE_OUT.relative_to(PDE_ROOT)}")
    print("fit rows:", len(sources["fit"].get("candidates", []) or []))
    print("fit constraint_kind split:", dict(sorted(fit_kind_split.items())))
    print("fit status split:", dict(sorted(fit_status_split.items())))
    print("parameter rows:", len(sources["provenance_records"]))
    print("parameter constraint_kind split:", dict(sorted(provenance_kind_split.items())))
    print("aliases surfaced:", len(sources["alias_to_canonical"]))


if __name__ == "__main__":
    main()
