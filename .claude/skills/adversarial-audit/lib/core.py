#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml


STATE_MACHINE = [
    "pending",
    "scanned",
    "provenance_built",
    "audit_pending",
    "audited",
    "defense_pending",
    "defended",
    "adjudicating",
    "verdict_logged",
    "blocked",
]


class NoAliasDumper(yaml.SafeDumper):
    def ignore_aliases(self, data: Any) -> bool:
        return True

WRITE_COMMANDS = {
    "init",
    "phase-a-scan",
    "phase-b-build",
    "phase-c-render",
    "set-status",
    "record-codex-defense",
    "dry-run",
    "purge-dry-run",
}

REQUIRED_DIRS = [
    "phase_a_fragments",
    "provenance",
    "reports",
    "defenses",
    "verdicts",
    "tmp_prompts",
    "codex_logs",
]

MODALITIES = [
    "numeric_literal",
    "claim_label",
    "graph",
    "existing_provenance",
]

CLAIM_KEYWORDS = [
    "fit",
    "fixed to",
    "matched",
    "matches",
    "matching",
    "benchmark",
    "calibrated",
    "calibration",
    "canonical condition",
    "determine",
    "determined",
    "external",
    "fingerprint",
    "normalization mismatch",
    "not free",
    "obstruction",
    "pin",
    "sourced benchmark",
]

PARAMETER_DENYLIST = {
    "a",
    "b",
    "c",
    "d",
    "g",
    "i",
    "j",
    "k",
    "l",
    "m",
    "n",
    "q",
    "r",
    "t",
    "w",
    "x",
    "y",
    "z",
    "A",
    "C",
    "I",
    "K",
    "L",
    "M",
    "N",
    "O",
    "P",
    "Q",
    "X",
    "Y",
    "c_s",
    "h_2",
    "j_2",
    "Lam_out_series",
    "Lambda",
    "Lambda_2",
    "m_0",
    "N_Q",
    "Omega_Q",
    "Omega_Q_3c_s_2a",
    "sigma_can",
    "Y_2",
    "Y_Q",
    "Yret_series",
    "y_2",
    "alpha",
    "beta",
    "gamma",
    "delta",
    "eta",
    "omega",
    "varpi",
    "atom_work",
    "claimstatus",
    "emph",
    "eqref",
    "expect_zero",
    "insufficient_verification",
    "item",
    "material_change",
    "mathematica_transliteration",
    "paper_alignment",
    "paper_misalignment",
    "script_missing_paper_claim",
    "section",
    "stagefield",
    "symbol_assumption_error",
    "widehat",
}

TEX_COMMAND_DENYLIST = {
    "begin",
    "bigl",
    "bigr",
    "boxed",
    "cdot",
    "claimstatus",
    "cos",
    "dfrac",
    "dot",
    "emph",
    "end",
    "eqref",
    "frac",
    "in",
    "item",
    "label",
    "left",
    "ln",
    "operatorname",
    "qquad",
    "quad",
    "ref",
    "resultanchor",
    "right",
    "rm",
    "section",
    "sin",
    "sqrt",
    "stagefield",
    "StageFile",
    "StatusExactClosure",
    "StatusReduced",
    "sum",
    "text",
    "widehat",
}


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def slug(text: str) -> str:
    out = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_").lower()
    return out or "item"


def load_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with path.open(encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def write_yaml(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        yaml.dump(data, f, Dumper=NoAliasDumper, default_flow_style=False, sort_keys=False, width=120, allow_unicode=True)
    tmp.replace(path)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def require_manifest_lock(command: str) -> None:
    if command in WRITE_COMMANDS and os.environ.get("ADVERSARIAL_MANIFEST_LOCKED") != "1":
        print(
            f"error: {command} writes MANIFEST.yaml and must be invoked through adversarial.sh flock locking",
            file=sys.stderr,
        )
        raise SystemExit(2)


def nested_get(data: dict[str, Any], dotted: str, default: Any = None) -> Any:
    cur: Any = data
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur


class Env:
    def __init__(self, config_path: Path):
        self.config_path = config_path.resolve()
        self.project_root = self.config_path.parent
        self.cfg = load_yaml(self.config_path)
        self.adv = self.cfg.get("adversarial") or {}
        if not self.adv:
            raise SystemExit("error: .redteam-config.yaml has no adversarial section")
        self.repo_root = self.resolve_project_path(self.adv.get("repo_root", "../.."))
        self.skill_dir = Path(__file__).resolve().parent.parent

    def adv_value(self, key: str, default: Any = None) -> Any:
        if key.startswith("codex."):
            return nested_get(self.adv, key, nested_get(self.cfg, key, default))
        return nested_get(self.adv, key, default)

    def resolve_project_path(self, value: str | Path) -> Path:
        text = str(value)
        text = os.path.expanduser(text)
        p = Path(text)
        if p.is_absolute():
            return p.resolve()
        return (self.project_root / p).resolve()

    @property
    def artifact_root(self) -> Path:
        return self.resolve_project_path(self.adv.get("artifact_root", "redteam_adversarial"))

    @property
    def manifest_path(self) -> Path:
        return self.artifact_root / "MANIFEST.yaml"

    @property
    def fit_path(self) -> Path:
        return self.artifact_root / "fit_insertion_points.yaml"

    @property
    def benchmarks_path(self) -> Path:
        return self.artifact_root / "benchmarks.yaml"

    @property
    def batches_path(self) -> Path:
        return self.artifact_root / "BATCHES.md"

    def config_path_value(self, key: str) -> Path:
        value = self.adv_value(key)
        if value is None:
            raise SystemExit(f"error: missing adversarial config key: {key}")
        return self.resolve_project_path(str(value))

    def rel(self, path: Path) -> str:
        path = path.resolve()
        try:
            return str(path.relative_to(self.project_root))
        except ValueError:
            try:
                return str(path.relative_to(self.repo_root))
            except ValueError:
                return str(path)

    def abs_from_rel(self, path_text: str) -> Path:
        p = Path(path_text)
        if p.is_absolute():
            return p
        candidate = self.project_root / p
        if candidate.exists() or not str(path_text).startswith(("docs/", "graph/")):
            return candidate
        return self.repo_root / p


def ensure_tree(env: Env) -> None:
    env.artifact_root.mkdir(parents=True, exist_ok=True)
    for name in REQUIRED_DIRS:
        (env.artifact_root / name).mkdir(parents=True, exist_ok=True)


def initial_manifest(env: Env) -> dict[str, Any]:
    ts = now_iso()
    return {
        "schema_version": 1,
        "project_name": env.adv.get("project_name", env.cfg.get("project", {}).get("name", "pde_ledger")),
        "artifact_root": env.adv.get("artifact_root", "redteam_adversarial"),
        "created_at": ts,
        "updated_at": ts,
        "state_machine": STATE_MACHINE,
        "candidates": {},
        "dry_runs": {},
    }


def load_manifest(env: Env) -> dict[str, Any]:
    m = load_yaml(env.manifest_path)
    if not m:
        return initial_manifest(env)
    m.setdefault("state_machine", STATE_MACHINE)
    m.setdefault("candidates", {})
    m.setdefault("dry_runs", {})
    return m


def save_manifest(env: Env, manifest: dict[str, Any]) -> None:
    manifest["updated_at"] = now_iso()
    write_yaml(env.manifest_path, manifest)


def transition(entry: dict[str, Any], status: str) -> None:
    if status not in STATE_MACHINE:
        raise SystemExit(f"error: invalid candidate status: {status}")
    entry["status"] = status
    entry.setdefault("status_timestamps", {})[status] = now_iso()


def validate_manual_transition(current: str | None, target: str) -> None:
    if target not in STATE_MACHINE:
        raise SystemExit(f"error: invalid candidate status: {target}")
    current = current or "pending"
    if current not in STATE_MACHINE:
        raise SystemExit(f"error: candidate has invalid current status: {current}")
    if target == current:
        return
    if target == "blocked":
        return
    current_i = STATE_MACHINE.index(current)
    target_i = STATE_MACHINE.index(target)
    if target_i == current_i + 1:
        return
    raise SystemExit(f"error: invalid status transition: {current} -> {target}")


def source_files_for_stage(env: Env, stage: str) -> list[tuple[str, Path]]:
    stage = f"{int(stage):03d}"
    paths: list[tuple[str, Path]] = []
    candidates = [
        ("paper_stage_tex", env.project_root / f"paper/stages/stage_{stage}.tex"),
    ]
    for role, path in candidates:
        if path.exists():
            paths.append((role, path))
    for path in sorted((env.project_root / "notes/stages").glob(f"moving_throat_pde_stage{stage}_*.md")):
        paths.append(("notes_stage", path))
    for path in sorted((env.project_root / "scripts").glob(f"moving_throat_pde_stage{stage}_*_sympy_audit.py")):
        paths.append(("sympy_script", path))
    for path in sorted((env.project_root / "mathematica").glob(f"moving_throat_pde_stage{stage}_*_mathematica_audit.wl")):
        paths.append(("mathematica_script", path))
    return paths


def pass2_report_path(env: Env, stage: str) -> Path:
    pattern = nested_get(env.adv, "phase_a_seeds.pass2_stage_reports_glob")
    if not pattern:
        return env.project_root / f"redteam/pass2/reports/stage_{int(stage):03d}.md"
    return env.resolve_project_path(pattern.format(stage=f"{int(stage):03d}"))


def checkpoint_provenance_path(env: Env) -> Path:
    return env.config_path_value("phase_a_seeds.checkpoint_provenance")


def has_fit_signal(text: str) -> bool:
    lower = text.lower()
    if any(keyword in lower for keyword in CLAIM_KEYWORDS):
        return True
    return bool(
        re.search(
            r"(?i)\b(?:value|constant|scalar|parameter|coefficient|normalization|branch|benchmark)\b.*\b(?:fixed|matched|matching|determined|pinned|not free)\b",
            text,
        )
    )


def is_scaffolding_line(text: str) -> bool:
    lower = text.lower()
    markers = (
        "\\stagefield{purpose}",
        "\\stagefield{downstream use}",
        "\\stagefield{verification}",
        "\\claimstatus",
        "\\section",
        "documentation index-garble",
        "paper_misalignment",
        "deliverable forms checked",
        "not a fully solved",
        "resolve before fix_loop",
    )
    return any(marker in lower for marker in markers)


def mentioned_stage(text: str, stages: list[str]) -> str | None:
    for stage in stages:
        n = int(stage)
        patterns = [
            rf"(?i)stage[\s_~`:-]*0*{n}(?!\d)",
        ]
        if any(re.search(pattern, text) for pattern in patterns):
            return f"{n:03d}"
    return None


def clean_param_piece(text: str) -> str:
    text = text.replace("\\rm", "")
    text = text.replace("\\mathrm", "")
    text = re.sub(r"[^A-Za-z0-9]+", "_", text)
    return re.sub(r"_+", "_", text).strip("_")


def normalized_param_name(text: str) -> str | None:
    text = text.strip().strip("`")
    text = text.replace("\\", "")
    text = text.replace("{", "").replace("}", "")
    text = text.replace("^", "_")
    text = clean_param_piece(text)
    if not text or text in PARAMETER_DENYLIST:
        return None
    if re.search(r"_Q[A-Za-z0-9_]", text):
        return None
    if len(text) == 1 and text.isalpha():
        return None
    return text


def add_param(params: list[str], value: str | None) -> None:
    if value and value not in PARAMETER_DENYLIST:
        params.append(value)


def detect_parameters(text: str) -> list[str]:
    params: list[str] = []
    for match in re.finditer(r"`([A-Za-z][A-Za-z0-9]*(?:_[A-Za-z0-9]+)+)`", text):
        add_param(params, normalized_param_name(match.group(1)))
    for match in re.finditer(r"\b([A-Za-z][A-Za-z0-9]*(?:_[A-Za-z0-9]+)+)\b", text):
        add_param(params, normalized_param_name(match.group(1)))
    for match in re.finditer(r"\\([A-Za-z]+)(?:_\{?([^}\s^]+)\}?)?", text):
        command, subscript = match.groups()
        if command in TEX_COMMAND_DENYLIST:
            continue
        if not subscript and command in PARAMETER_DENYLIST:
            continue
        if subscript:
            add_param(params, normalized_param_name(f"{command}_{subscript}"))
        elif command not in PARAMETER_DENYLIST:
            add_param(params, normalized_param_name(command))
    return sorted(set(params))


def candidate_parameters(text: str) -> list[str]:
    params = detect_parameters(text)
    lower = text.lower()
    has_dtn = "dtn" in lower or "fingerprint" in lower
    has_value = bool(re.search(r"(\\frac\{|\d+\s*/\s*\d+|\b\d+(?:\.\d+)?\b|=)", text))
    is_value_claim = has_fit_signal(text)
    if has_dtn and has_value and is_value_claim:
        params.append("matched_fingerprint_value")
    if params and is_value_claim:
        return sorted(set(params))
    return []


def candidate_key_for(stage: str, text: str) -> str:
    params = candidate_parameters(text)
    if params:
        return f"stage_{int(stage):03d}_{slug(params[0])}"
    raise ValueError("cannot create candidate without parameter/value class")


def dry_run_fixture_group_key(frag: dict[str, Any]) -> str:
    """Fixture-only grouping for the mandated 104/105 chi_Q dry-run episode."""
    params = set(frag.get("parameter_names") or [])
    citation = frag.get("citation") or {}
    text = str(citation.get("excerpt") or "").lower()
    has_chi_episode_param = bool(params & {"chi_Q", "Gamma_5", "sigma_Q", "matched_fingerprint_value", "outgoing_l2_DtN_fingerprint"})
    if has_chi_episode_param and ("dtn" in text or "fingerprint" in text or "chi" in text or "sigma" in text):
        return "chi_q_outgoing_dtn"
    return str(frag["candidate_key"])


def make_fragment(
    env: Env,
    *,
    modality: str,
    stage: str,
    role: str,
    path: Path,
    line: int,
    text: str,
    reason: str,
) -> dict[str, Any]:
    params = candidate_parameters(text)
    if not params:
        raise ValueError("cannot create fragment without parameter/value class")
    return {
        "candidate_key": candidate_key_for(stage, text),
        "modality": modality,
        "anchor_stage": f"{int(stage):03d}",
        "parameter_names": params,
        "reason": reason,
        "citation": {
            "path": env.rel(path),
            "line": line,
            "role": role,
            "stage": f"{int(stage):03d}",
            "excerpt": text.strip(),
        },
    }


def iter_file_lines(path: Path) -> list[tuple[int, str]]:
    try:
        lines = read_text(path).splitlines()
    except OSError:
        return []
    return list(enumerate(lines, 1))


def numeric_literal_scan(env: Env, stages: list[str]) -> list[dict[str, Any]]:
    literal_re = re.compile(r"(\\frac\{|\d+\s*/\s*\d+|\b\d+(?:\.\d+)?\b|=)")
    out: list[dict[str, Any]] = []
    for stage in stages:
        for role, path in source_files_for_stage(env, stage):
            if role not in {"notes_stage", "paper_stage_tex"}:
                continue
            for line_no, line in iter_file_lines(path):
                if is_scaffolding_line(line):
                    continue
                if not literal_re.search(line):
                    continue
                if not has_fit_signal(line):
                    continue
                if not candidate_parameters(line):
                    continue
                out.append(
                    make_fragment(
                        env,
                        modality="numeric_literal",
                        stage=stage,
                        role=role,
                        path=path,
                        line=line_no,
                        text=line,
                        reason="target-related numerical literal or closed-form coefficient",
                    )
                )
    return out


def claim_label_scan(env: Env, stages: list[str]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for stage in stages:
        for role, path in source_files_for_stage(env, stage):
            if role not in {"notes_stage", "paper_stage_tex"}:
                continue
            for line_no, line in iter_file_lines(path):
                if is_scaffolding_line(line):
                    continue
                lower = line.lower()
                if not any(keyword in lower for keyword in CLAIM_KEYWORDS):
                    continue
                if not candidate_parameters(line):
                    continue
                out.append(
                    make_fragment(
                        env,
                        modality="claim_label",
                        stage=stage,
                        role=role,
                        path=path,
                        line=line_no,
                        text=line,
                        reason="claim label or status wording near target-related parameter",
                    )
                )
    return out


def existing_provenance_scan(env: Env, stages: list[str]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for stage in stages:
        report = pass2_report_path(env, stage)
        if report.exists():
            for line_no, line in iter_file_lines(report):
                if is_scaffolding_line(line):
                    continue
                if "moving_throat_pde_stage" in line and any(ext in line for ext in (".py", ".wl", ".md", ".tex", ".txt")):
                    continue
                if not any(k in line.lower() for k in ("f1", "pin", "fingerprint", "provenance", "fixed", "constant", "normalization", "benchmark")):
                    continue
                if not candidate_parameters(line):
                    continue
                out.append(
                    make_fragment(
                        env,
                        modality="existing_provenance",
                        stage=stage,
                        role="pass2_stage_report",
                        path=report,
                        line=line_no,
                        text=line,
                        reason="pass-2 reconciliation or red-team provenance seed mentions a candidate value",
                    )
                )
    cp = checkpoint_provenance_path(env)
    if cp.exists():
        for line_no, line in iter_file_lines(cp):
            if not candidate_parameters(line):
                continue
            anchor = mentioned_stage(line, stages)
            if anchor is None:
                continue
            if any(k in line.lower() for k in ("provenance", "stage", "derived", "fingerprint", "dtn", "parameter", "constant")):
                out.append(
                    make_fragment(
                        env,
                        modality="existing_provenance",
                        stage=anchor,
                        role="checkpoint_constant_provenance",
                        path=cp,
                        line=line_no,
                        text=line,
                        reason="checkpoint provenance seed mentions a candidate parameter",
                    )
                )
    return out


def source_query_variants(env: Env, path: Path, role: str, stage: str) -> list[str]:
    variants: list[str] = []

    def add(value: str) -> None:
        if value and value not in variants:
            variants.append(value)

    add(env.rel(path))
    try:
        add(str(path.resolve().relative_to(env.repo_root)))
    except ValueError:
        pass
    if role == "notes_stage":
        paper = env.project_root / f"paper/stages/stage_{int(stage):03d}.tex"
        if paper.exists():
            add(env.rel(paper))
            try:
                add(str(paper.resolve().relative_to(env.repo_root)))
            except ValueError:
                pass
    return variants


def query_graph_source(env: Env, source_queries: list[str]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    query_graph = env.config_path_value("query_graph")
    graph_path = env.config_path_value("atlas_graph")
    errors: list[dict[str, Any]] = []
    for query in source_queries:
        cmd = ["timeout", "600", "python3", str(query_graph), "--graph", str(graph_path), "source", query, "--format", "json"]
        proc = subprocess.run(cmd, cwd=env.repo_root, text=True, capture_output=True, check=False)
        if proc.returncode != 0:
            errors.append({"source": query, "reason": f"query_graph exit {proc.returncode}", "stderr": proc.stderr.strip()[:500]})
            continue
        payload = json.loads(proc.stdout or "{}")
        nodes = payload.get("nodes") or []
        if nodes:
            return nodes, []
        errors.append({"source": query, "reason": "no atlas node tied to this source path"})
    return [], errors


def graph_scan(env: Env, stages: list[str]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    fragments: list[dict[str, Any]] = []
    gaps: list[dict[str, Any]] = []
    for stage in stages:
        for role, path in source_files_for_stage(env, stage):
            if role not in {"paper_stage_tex", "notes_stage"}:
                continue
            source_queries = source_query_variants(env, path, role, stage)
            nodes, errors = query_graph_source(env, source_queries)
            query = source_queries[0]
            if not nodes:
                gaps.append(
                    {
                        "stage": f"{int(stage):03d}",
                        "source": query,
                        "graph_gap": True,
                        "attempted_sources": source_queries,
                        "reason": "; ".join(sorted({e["reason"] for e in errors})) or "no atlas node tied to this source path",
                        "stderr": next((e.get("stderr") for e in errors if e.get("stderr")), ""),
                    }
                )
                continue
            for node in nodes:
                text = " ".join(str(v) for v in node.values() if not isinstance(v, (dict, list)))
                if candidate_parameters(text):
                    fragments.append(
                        {
                            "candidate_key": candidate_key_for(stage, text),
                            "modality": "graph",
                            "anchor_stage": f"{int(stage):03d}",
                            "parameter_names": candidate_parameters(text),
                            "reason": "atlas graph node carries target-related text",
                            "citation": {
                                "path": query,
                                "line": None,
                                "role": role,
                                "excerpt": node.get("id", ""),
                            },
                            "graph_node": node.get("id"),
                        }
                    )
    return fragments, gaps


def merge_fragments(fragments: list[dict[str, Any]], *, dry_run: bool, dry_run_id: str | None) -> list[dict[str, Any]]:
    grouped: dict[str, dict[str, Any]] = {}
    for frag in fragments:
        key = dry_run_fixture_group_key(frag) if dry_run else frag["candidate_key"]
        entry = grouped.setdefault(
            key,
            {
                "id": "",
                "candidate_key": key,
                "dry_run": dry_run,
                "dry_run_id": dry_run_id,
                "anchor_stages": [],
                "parameter_names": [],
                "citations": [],
                "modality_attribution": [],
                "modality_fragments": [],
            },
        )
        stage = frag["anchor_stage"]
        if stage not in entry["anchor_stages"]:
            entry["anchor_stages"].append(stage)
        for param in frag.get("parameter_names") or []:
            if param not in entry["parameter_names"]:
                entry["parameter_names"].append(param)
        if frag["modality"] not in entry["modality_attribution"]:
            entry["modality_attribution"].append(frag["modality"])
        cit = frag.get("citation")
        if cit and cit not in entry["citations"]:
            entry["citations"].append(cit)
        entry["modality_fragments"].append(frag)
    out = []
    for key, entry in grouped.items():
        entry["anchor_stages"] = sorted(entry["anchor_stages"])
        entry["parameter_names"] = sorted(entry["parameter_names"]) or ["unclassified_target_parameter"]
        entry["modality_attribution"] = sorted(entry["modality_attribution"])
        prefix = f"dryrun_{dry_run_id}_" if dry_run else "fit_"
        entry["id"] = prefix + key
        out.append(entry)
    out.sort(key=lambda item: item["id"])
    return out


def initial_candidate_entry(candidate: dict[str, Any]) -> dict[str, Any]:
    ts = now_iso()
    params = candidate.get("parameter_names") or []
    return {
        "id": candidate["id"],
        "candidate_key": candidate["candidate_key"],
        "dry_run": bool(candidate.get("dry_run")),
        "dry_run_id": candidate.get("dry_run_id"),
        "anchor_stages": candidate.get("anchor_stages") or [],
        "file_line_citations": candidate.get("citations") or [],
        "parameter_names": params,
        "batch_id": None,
        "status": "pending",
        "status_timestamps": {"pending": ts},
        "codex_session": {
            "by_parameter": {
                param: {
                    "session_id": None,
                    "log_paths": [],
                    "last_iter": None,
                    "last_exit": None,
                    "last_run": None,
                }
                for param in params
            }
        },
        "paths": {
            "report": None,
            "defense": None,
            "verdict": None,
            "provenance": [],
            "phase_c_prompt": None,
        },
        "verdict": {
            "adversarial": None,
            "adjudication": None,
        },
        "modality_attribution": candidate.get("modality_attribution") or [],
        "modality_fragments": candidate.get("modality_fragments") or [],
    }


def render_fit_file(env: Env, manifest: dict[str, Any], stage_results: dict[str, Any] | None = None) -> None:
    candidates = []
    for entry in (manifest.get("candidates") or {}).values():
        candidates.append(
            {
                "id": entry["id"],
                "dry_run": entry.get("dry_run", False),
                "dry_run_id": entry.get("dry_run_id"),
                "anchor_stages": entry.get("anchor_stages", []),
                "parameter_names": entry.get("parameter_names", []),
                "status": entry.get("status"),
                "file_line_citations": entry.get("file_line_citations", []),
                "modality_attribution": entry.get("modality_attribution", []),
                "batch_id": entry.get("batch_id"),
            }
        )
    candidates.sort(key=lambda item: item["id"])
    if stage_results is None:
        stage_results = {}
        for dry_run in (manifest.get("dry_runs") or {}).values():
            stage_results.update(dry_run.get("stage_results") or {})
    payload = {
        "schema_version": 1,
        "generated_at": now_iso(),
        "dry_run_artifacts_present": any(c.get("dry_run") for c in candidates),
        "stage_results": stage_results or {},
        "candidates": candidates,
    }
    write_yaml(env.fit_path, payload)


def render_batches(env: Env, manifest: dict[str, Any]) -> None:
    counts = defaultdict(int)
    dry_counts = defaultdict(int)
    verdict_count = 0
    for entry in (manifest.get("candidates") or {}).values():
        counts[entry.get("status", "unknown")] += 1
        if entry.get("dry_run"):
            dry_counts[entry.get("status", "unknown")] += 1
        verdict = entry.get("verdict") or {}
        if verdict.get("adversarial") or verdict.get("adjudication"):
            verdict_count += 1
    lines = [
        "# Adversarial Audit Status",
        "",
        f"Generated: {now_iso()}",
        f"Project: {manifest.get('project_name', '?')}",
        "",
        "| Scope | Counts |",
        "|---|---|",
        f"| all candidates | {' '.join(f'{k}={v}' for k, v in sorted(counts.items())) or 'none'} |",
        f"| dry-run candidates | {' '.join(f'{k}={v}' for k, v in sorted(dry_counts.items())) or 'none'} |",
        f"| binding verdict fields populated | {verdict_count} |",
        "",
        "Phase C batches are intentionally unset until the Step 3 batching consult.",
    ]
    write_text(env.batches_path, "\n".join(lines) + "\n")


def ensure_benchmarks_file(env: Env) -> dict[str, Any]:
    data = load_yaml(env.benchmarks_path)
    if not data:
        data = {"schema_version": 1, "entries": []}
    data.setdefault("entries", [])
    return data


def update_benchmarks_for_candidate(env: Env, candidate: dict[str, Any]) -> None:
    data = ensure_benchmarks_file(env)
    cid = candidate["id"]
    entries = [e for e in data.get("entries", []) if e.get("candidate_id") != cid]
    citations = [
        c
        for c in candidate.get("file_line_citations", [])
        if "104" in c.get("path", "") or "105" in c.get("path", "")
    ][:8]
    params = set(candidate.get("parameter_names") or [])
    if candidate.get("dry_run") and "chi_Q" in params and citations:
        entries.append(
            {
                "id": f"{cid}__dry_run_benchmark_placeholder",
                "candidate_id": cid,
                "dry_run": True,
                "dry_run_id": candidate.get("dry_run_id"),
                "non_binding": True,
                "fixture_type": "dry-run source-anchor placeholder",
                "claim": "DRY-RUN PLACEHOLDER: Phase C must use sourced benchmark entries, not model memory.",
                "value": None,
                "source_type": "dry-run project-source anchor list",
                "source_citations": citations,
                "requires_agent_sourcing": True,
                "external_match_policy": "For real campaign use, replace this placeholder with an agent-built sourced benchmark entry.",
            }
        )
    data["entries"] = entries
    write_yaml(env.benchmarks_path, data)


def render_template(env: Env, template_name: str, replacements: dict[str, str]) -> str:
    template_path = env.skill_dir / "prompts" / template_name
    text = read_text(template_path)
    for key, value in replacements.items():
        text = text.replace("{" + key + "}", value)
    return text


def yaml_block(data: Any) -> str:
    return yaml.safe_dump(data, default_flow_style=False, sort_keys=False, width=120, allow_unicode=True).rstrip()


def stage_source_context(env: Env, stages: list[str]) -> list[dict[str, Any]]:
    sources: list[dict[str, Any]] = []
    for stage in stages:
        for role, path in source_files_for_stage(env, stage):
            sources.append({"stage": f"{int(stage):03d}", "role": role, "path": env.rel(path)})
    return sources


def phase_a_stage_reports(env: Env, stages: list[str]) -> list[str]:
    return [env.rel(pass2_report_path(env, stage)) for stage in stages if pass2_report_path(env, stage).exists()]


def render_phase_a_modality_prompts(env: Env, stages: list[str], prefix: str, dry_run: bool = False) -> dict[str, str]:
    ensure_tree(env)
    stages = [f"{int(s):03d}" for s in stages]
    source_files = stage_source_context(env, stages)
    graph_cmd = " ".join(
        [
            "timeout",
            "600",
            "python3",
            env.rel(env.config_path_value("query_graph")),
            "--graph",
            env.rel(env.config_path_value("atlas_graph")),
        ]
    )
    common = {
        "STAGES": ", ".join(stages),
        "SOURCE_FILES": yaml_block(source_files),
        "PASS2_RECONCILIATION": env.rel(env.config_path_value("phase_a_seeds.pass2_reconciliation")),
        "PASS2_STAGE_REPORTS": yaml_block(phase_a_stage_reports(env, stages)),
        "CHECKPOINT_PROVENANCE": env.rel(checkpoint_provenance_path(env)),
        "GRAPH_WRAPPER": graph_cmd,
    }
    templates = {
        "numeric_literal": "phase_a_numeric_literal.md",
        "claim_label": "phase_a_claim_label.md",
        "graph": "phase_a_graph.md",
        "existing_provenance": "phase_a_existing_provenance.md",
    }
    rendered: dict[str, str] = {}
    for modality, template in templates.items():
        prompt = render_template(env, template, common)
        if dry_run:
            prompt = "dry_run: true\nnon_binding: true\n\n" + prompt
        path = env.artifact_root / "tmp_prompts" / f"{prefix}_phase_a_{modality}.md"
        write_text(path, prompt)
        rendered[modality] = env.rel(path)
    return rendered


def run_phase_a(env: Env, stages: list[str], dry_run: bool, dry_run_id: str | None) -> dict[str, Any]:
    ensure_tree(env)
    stages = [f"{int(s):03d}" for s in stages]
    prefix = dry_run_id or "phase_a"
    modality_prompts = render_phase_a_modality_prompts(env, stages, prefix=prefix, dry_run=dry_run)
    fragments_by_modality: dict[str, list[dict[str, Any]]] = {
        "numeric_literal": numeric_literal_scan(env, stages),
        "claim_label": claim_label_scan(env, stages),
        "existing_provenance": existing_provenance_scan(env, stages),
    }
    graph_fragments, graph_gaps = graph_scan(env, stages)
    fragments_by_modality["graph"] = graph_fragments

    for modality in MODALITIES:
        payload = {
            "schema_version": 1,
            "dry_run": dry_run,
            "dry_run_id": dry_run_id,
            "modality": modality,
            "blind_to_modalities": [m for m in MODALITIES if m != modality],
            "stages": stages,
            "candidates": fragments_by_modality.get(modality, []),
        }
        if modality == "graph":
            payload["graph_gaps"] = graph_gaps
        write_yaml(env.artifact_root / "phase_a_fragments" / f"{prefix}_{modality}.yaml", payload)

    all_fragments = [frag for frags in fragments_by_modality.values() for frag in frags]
    candidates = merge_fragments(all_fragments, dry_run=dry_run, dry_run_id=dry_run_id)

    manifest = load_manifest(env)
    manifest.setdefault("dry_runs", {})
    if dry_run and dry_run_id:
        manifest["dry_runs"][dry_run_id] = {
            "dry_run": True,
            "stages": stages,
            "started_at": manifest["dry_runs"].get(dry_run_id, {}).get("started_at", now_iso()),
            "last_phase": "phase_a_scan",
        }

    for candidate in candidates:
        entry = manifest["candidates"].get(candidate["id"]) or initial_candidate_entry(candidate)
        entry["file_line_citations"] = candidate.get("citations", [])
        entry["parameter_names"] = candidate.get("parameter_names", [])
        entry["modality_attribution"] = candidate.get("modality_attribution", [])
        entry["modality_fragments"] = candidate.get("modality_fragments", [])
        transition(entry, "scanned")
        manifest["candidates"][candidate["id"]] = entry

    stage_results = {}
    for stage in stages:
        count = sum(1 for c in candidates if stage in c.get("anchor_stages", []))
        stage_results[stage] = {
            "candidate_count": count,
            "structurally_vacuous": count == 0,
            "dry_run": dry_run,
        }
    if dry_run and dry_run_id:
        manifest["dry_runs"][dry_run_id]["stage_results"] = stage_results
    render_fit_file(env, manifest, stage_results=stage_results)
    render_batches(env, manifest)
    save_manifest(env, manifest)

    critic_prompt = render_template(
        env,
        "phase_a_completeness_critic.md",
        {
            "STAGES": ", ".join(stages),
            "FIT_INSERTION_POINTS_YAML": yaml_block({"candidates": candidates, "stage_results": stage_results}),
            "MODALITY_FRAGMENT_PATHS": "\n".join(
                f"- {env.rel(env.artifact_root / 'phase_a_fragments' / f'{prefix}_{m}.yaml')}" for m in MODALITIES
            ),
        },
    )
    critic_path = env.artifact_root / "tmp_prompts" / f"{prefix}_phase_a_completeness_critic.md"
    write_text(critic_path, critic_prompt)

    return {
        "stages": stages,
        "candidate_ids": [c["id"] for c in candidates],
        "stage_results": stage_results,
        "graph_gaps": graph_gaps,
        "modality_prompts": modality_prompts,
        "critic_prompt": env.rel(critic_path),
    }


def source_evidence_for_candidate(env: Env, entry: dict[str, Any]) -> list[dict[str, Any]]:
    evidence: list[dict[str, Any]] = []
    params = [str(p) for p in entry.get("parameter_names") or []]
    param_needles = {p for p in params}
    param_needles.update(p.replace("_", "\\_") for p in params)
    param_needles.update(p.replace("_", "") for p in params)
    token_needles: set[str] = set()
    for param in params:
        for token in re.split(r"[^A-Za-z0-9]+", param):
            if len(token) >= 3:
                token_needles.add(token.lower())
    cited_paths = {str(c.get("path")) for c in entry.get("file_line_citations", []) if c.get("path")}
    for stage in entry.get("anchor_stages", []):
        for role, path in source_files_for_stage(env, stage):
            if role not in {"notes_stage", "sympy_script", "paper_stage_tex"}:
                continue
            for line_no, line in iter_file_lines(path):
                mentions_param = any(needle and needle in line for needle in param_needles)
                lower = line.lower()
                mentions_token = any(token in lower for token in token_needles)
                same_source_fit_line = env.rel(path) in cited_paths and has_fit_signal(line)
                if not mentions_param and not mentions_token and not same_source_fit_line:
                    continue
                evidence.append(
                    {
                        "path": env.rel(path),
                        "line": line_no,
                        "role": role,
                        "excerpt": line.strip(),
                    }
                )
    return evidence[:40]


def graph_context_for_candidate(env: Env, entry: dict[str, Any]) -> dict[str, Any]:
    contexts = []
    gaps = []
    seen_sources: set[tuple[str, ...]] = set()
    for citation in entry.get("file_line_citations", []):
        if citation.get("role") not in {"notes_stage", "paper_stage_tex"}:
            continue
        query = citation.get("path")
        if not query:
            continue
        role = str(citation.get("role"))
        stage = str(citation.get("stage") or (entry.get("anchor_stages") or ["0"])[0])
        source_path = env.abs_from_rel(str(query))
        source_queries = source_query_variants(env, source_path, role, stage)
        key = tuple(source_queries)
        if key in seen_sources:
            continue
        seen_sources.add(key)
        nodes, errors = query_graph_source(env, source_queries)
        if not nodes:
            gaps.append(
                {
                    "source": query,
                    "attempted_sources": source_queries,
                    "graph_gap": True,
                    "reason": "; ".join(sorted({e["reason"] for e in errors})) or "no atlas node tied to this source path",
                }
            )
        else:
            contexts.append({"source": query, "source_queries": source_queries, "nodes": nodes})
    return {"contexts": contexts, "graph_gaps": gaps}


def build_provenance(env: Env, candidate_id: str) -> dict[str, Any]:
    manifest = load_manifest(env)
    entry = (manifest.get("candidates") or {}).get(candidate_id)
    if not entry:
        raise SystemExit(f"error: unknown candidate: {candidate_id}")

    evidence = source_evidence_for_candidate(env, entry)
    graph_context = graph_context_for_candidate(env, entry)
    provenance_paths: list[str] = []
    for param in entry.get("parameter_names") or ["unclassified_target_parameter"]:
        prompt_payload = dict(entry)
        prompt_payload["phase_b_parameter_name"] = param
        phase_b_prompt = render_template(
            env,
            "phase_b_provenance_builder.md",
            {
                "CANDIDATE_ID": candidate_id,
                "CANDIDATE_YAML": yaml_block(prompt_payload),
            },
        )
        prompt_path = env.artifact_root / "tmp_prompts" / f"{candidate_id}__phase_b__{slug(param)}.md"
        write_text(prompt_path, phase_b_prompt)
        payload = {
            "schema_version": 1,
            "candidate_id": candidate_id,
            "parameter_name": param,
            "dry_run": bool(entry.get("dry_run")),
            "dry_run_id": entry.get("dry_run_id"),
            "non_binding": bool(entry.get("dry_run")),
            "generated_content_kind": "mechanical_evidence_bundle",
            "agent_prompt_path": env.rel(prompt_path),
            "agent_synthesis_required": True,
            "dry_run_fixture": bool(entry.get("dry_run")),
            "fixture_note": (
                "Dry-run provenance stores only opened source evidence and graph context. "
                "Origin claims, constraints, and contradiction findings must be produced by the Phase B agent path."
            )
            if entry.get("dry_run")
            else None,
            "taxonomy": {
                "ledger_scope": "parameter-value provenance only",
                "excluded_layers": ["file/source provenance", "result/stage provenance without value genealogy"],
            },
            "anchor_stages": entry.get("anchor_stages", []),
            "notes_source_opened": True,
            "source_evidence": evidence,
            "origin_claims": [],
            "constraints": [],
            "graph_context": graph_context,
            "provenance_findings": [],
        }
        if graph_context.get("graph_gaps"):
            payload["provenance_findings"].append(
                {
                    "type": "graph_gap",
                    "severity": "dry_run_non_binding" if entry.get("dry_run") else "needs_triage",
                    "summary": "The atlas graph wrapper returned no source nodes for one or more primary sources.",
                    "gaps": graph_context["graph_gaps"],
                }
            )
        path = env.artifact_root / "provenance" / f"{candidate_id}__{slug(param)}.yaml"
        write_yaml(path, payload)
        provenance_paths.append(env.rel(path))

    entry.setdefault("paths", {})["provenance"] = provenance_paths
    transition(entry, "provenance_built")
    manifest["candidates"][candidate_id] = entry
    update_benchmarks_for_candidate(env, entry)
    render_fit_file(env, manifest)
    render_batches(env, manifest)
    save_manifest(env, manifest)
    return {"candidate_id": candidate_id, "provenance_paths": provenance_paths, "status": entry["status"]}


def load_benchmarks_for_candidate(env: Env, candidate_id: str) -> list[dict[str, Any]]:
    data = ensure_benchmarks_file(env)
    return [e for e in data.get("entries", []) if e.get("candidate_id") == candidate_id]


def render_phase_c(env: Env, candidate_id: str) -> dict[str, Any]:
    manifest = load_manifest(env)
    entry = (manifest.get("candidates") or {}).get(candidate_id)
    if not entry:
        raise SystemExit(f"error: unknown candidate: {candidate_id}")
    frozen = read_text(env.config_path_value("frozen_directive"))
    provenance_docs = []
    for path_text in entry.get("paths", {}).get("provenance") or []:
        path = env.abs_from_rel(path_text)
        if path.exists():
            provenance_docs.append(load_yaml(path))
    graph_context = graph_context_for_candidate(env, entry)
    primary_sources = sorted({c.get("path") for c in entry.get("file_line_citations", []) if c.get("path")})
    prompt = render_template(
        env,
        "phase_c_adversarial.md",
        {
            "CANDIDATE_ID": candidate_id,
            "DRY_RUN_BLOCK": yaml_block({"dry_run": bool(entry.get("dry_run")), "non_binding": bool(entry.get("dry_run"))}),
            "FROZEN_DIRECTIVE": frozen.rstrip(),
            "CANDIDATE_YAML": yaml_block(entry),
            "PRIMARY_SOURCES": "\n".join(f"- {p}" for p in primary_sources),
            "PROVENANCE_SLICE": yaml_block(provenance_docs),
            "BENCHMARKS": yaml_block(load_benchmarks_for_candidate(env, candidate_id)),
            "GRAPH_CONTEXT": yaml_block(graph_context),
        },
    )
    prompt_path = env.artifact_root / "tmp_prompts" / f"{candidate_id}__phase_c_adversarial.md"
    write_text(prompt_path, prompt)
    report_path = env.artifact_root / "reports" / f"{candidate_id}__dry_run_prompt_render.md"
    if entry.get("dry_run"):
        write_text(
            report_path,
            "\n".join(
                [
                    "---",
                    "dry_run: true",
                    "non_binding: true",
                    f"candidate_id: {candidate_id}",
                    "binding_verdict: null",
                    "adversarial_verdict: null",
                    "---",
                    "",
                    "# Dry-Run Phase C Prompt Render",
                    "",
                    "This file records that Phase C prompt assembly completed. It is not an adversarial report and records no verdict.",
                    "",
                    f"- prompt: {env.rel(prompt_path)}",
                ]
            )
            + "\n",
        )
        entry.setdefault("paths", {})["report"] = env.rel(report_path)
    entry.setdefault("paths", {})["phase_c_prompt"] = env.rel(prompt_path)
    transition(entry, "audit_pending")
    manifest["candidates"][candidate_id] = entry
    render_fit_file(env, manifest)
    render_batches(env, manifest)
    save_manifest(env, manifest)
    return {"candidate_id": candidate_id, "phase_c_prompt": env.rel(prompt_path), "status": entry["status"]}


def render_defense_prompt(env: Env, candidate_id: str, parameter: str) -> str:
    manifest = load_manifest(env)
    entry = (manifest.get("candidates") or {}).get(candidate_id)
    if not entry:
        raise SystemExit(f"error: unknown candidate: {candidate_id}")
    report_path = entry.get("paths", {}).get("report")
    report_text = "(no adversarial report path is recorded)"
    if report_path:
        p = env.abs_from_rel(report_path)
        if p.exists():
            report_text = read_text(p)
    provenance_docs = []
    for path_text in entry.get("paths", {}).get("provenance") or []:
        if slug(parameter) in path_text or parameter in path_text:
            p = env.abs_from_rel(path_text)
            if p.exists():
                provenance_docs.append(load_yaml(p))
    prompt = render_template(
        env,
        "codex_defense.md",
        {
            "CANDIDATE_ID": candidate_id,
            "PARAMETER": parameter,
            "CANDIDATE_YAML": yaml_block(entry),
            "ADVERSARIAL_REPORT": report_text.rstrip(),
            "PROVENANCE_SLICE": yaml_block(provenance_docs),
            "BENCHMARKS": yaml_block(load_benchmarks_for_candidate(env, candidate_id)),
        },
    )
    prompt_path = env.artifact_root / "tmp_prompts" / f"{candidate_id}__defense__{slug(parameter)}.md"
    write_text(prompt_path, prompt)
    return str(prompt_path)


def codex_log_path(env: Env, candidate_id: str, parameter: str, iteration: str) -> str:
    return str(env.artifact_root / "codex_logs" / f"{candidate_id}__{slug(parameter)}__iter{iteration}.txt")


def parameter_session(env: Env, candidate_id: str, parameter: str) -> str:
    manifest = load_manifest(env)
    entry = (manifest.get("candidates") or {}).get(candidate_id) or {}
    by_param = (((entry.get("codex_session") or {}).get("by_parameter")) or {})
    value = (by_param.get(parameter) or by_param.get(slug(parameter)) or {}).get("session_id")
    return "" if value is None else str(value)


def parse_defense_outcomes(log_path: Path) -> list[dict[str, Any]]:
    text = read_text(log_path) if log_path.exists() else ""
    outcomes: list[dict[str, Any]] = []
    for line_no, line in enumerate(text.splitlines(), 1):
        upper = line.upper()
        if "|" in line or all(token in upper for token in ("DEFEND", "CONCEDE", "PARTIAL")):
            continue
        field_match = re.search(
            r"(?i)^\s*(?:defense_verdict|verdict|outcome|finding[_ -]?\d*[_ -]?outcome)\s*:\s*(DEFEND|CONCEDE|PARTIAL)\b",
            line,
        )
        bullet_match = re.search(r"(?i)^\s*[-*]\s*(?:finding[^:]*:\s*)?(DEFEND|CONCEDE|PARTIAL)\b", line)
        match = field_match or bullet_match
        if match:
            outcomes.append({"line": line_no, "outcome": match.group(1).upper(), "excerpt": line.strip()[:300]})
    return outcomes


def record_codex_defense(
    env: Env,
    candidate_id: str,
    parameter: str,
    iteration: str,
    rc: str,
    log_path: str,
    new_session: str,
) -> dict[str, Any]:
    manifest = load_manifest(env)
    entry = (manifest.get("candidates") or {}).get(candidate_id)
    if not entry:
        raise SystemExit(f"error: unknown candidate: {candidate_id}")
    by_param = entry.setdefault("codex_session", {}).setdefault("by_parameter", {})
    record = by_param.setdefault(parameter, {"session_id": None, "log_paths": []})
    if new_session:
        record["session_id"] = new_session
    record.setdefault("log_paths", []).append(env.rel(Path(log_path)))
    record["last_iter"] = int(iteration)
    record["last_exit"] = int(rc)
    record["last_run"] = now_iso()
    outcomes = parse_defense_outcomes(Path(log_path)) if int(rc) == 0 else []
    if int(rc) == 0 and not outcomes:
        transition(entry, "defense_pending")
        save_manifest(env, manifest)
        raise SystemExit(
            f"error: defense log contains no per-finding DEFEND/CONCEDE/PARTIAL outcome: {env.rel(Path(log_path))}"
        )
    if outcomes:
        record["last_outcomes"] = outcomes
        entry.setdefault("defense_outcomes", {})[parameter] = outcomes
    transition(entry, "defended" if int(rc) == 0 else "defense_pending")
    save_manifest(env, manifest)
    return {"candidate_id": candidate_id, "parameter": parameter, "session_id": record.get("session_id"), "outcomes": outcomes}


def cmd_init(env: Env, _args: argparse.Namespace) -> None:
    ensure_tree(env)
    manifest = load_manifest(env)
    save_manifest(env, manifest)
    if not env.fit_path.exists():
        render_fit_file(env, manifest)
    if not env.benchmarks_path.exists():
        write_yaml(env.benchmarks_path, {"schema_version": 1, "entries": []})
    render_batches(env, manifest)
    print(f"init: {env.rel(env.artifact_root)}")


def cmd_status(env: Env, _args: argparse.Namespace) -> None:
    manifest = load_manifest(env)
    counts = defaultdict(int)
    dry_counts = defaultdict(int)
    verdict_count = 0
    for entry in (manifest.get("candidates") or {}).values():
        counts[entry.get("status", "unknown")] += 1
        if entry.get("dry_run"):
            dry_counts[entry.get("status", "unknown")] += 1
        verdict = entry.get("verdict") or {}
        if verdict.get("adversarial") or verdict.get("adjudication"):
            verdict_count += 1
    payload = {
        "project": manifest.get("project_name"),
        "manifest": env.rel(env.manifest_path),
        "state_counts": dict(sorted(counts.items())),
        "dry_run_state_counts": dict(sorted(dry_counts.items())),
        "binding_verdict_fields_populated": verdict_count,
        "dry_runs": manifest.get("dry_runs", {}),
    }
    print(yaml_block(payload))


def cmd_candidate_info(env: Env, args: argparse.Namespace) -> None:
    entry = (load_manifest(env).get("candidates") or {}).get(args.candidate_id)
    if not entry:
        raise SystemExit(f"error: unknown candidate: {args.candidate_id}")
    print(yaml_block(entry))


def cmd_render_phase_a_prompts(env: Env, args: argparse.Namespace) -> None:
    if not args.stages:
        raise SystemExit("error: render-phase-a-prompts requires --stages")
    prefix = args.prefix or f"phase_a_{args.stages[0]}_{args.stages[-1]}"
    paths = render_phase_a_modality_prompts(env, args.stages, prefix=prefix, dry_run=args.dry_run)
    print(yaml_block({"stages": [f"{int(s):03d}" for s in args.stages], "modality_prompts": paths}))


def cmd_phase_a_scan(env: Env, args: argparse.Namespace) -> None:
    if not args.stages:
        raise SystemExit("error: phase-a-scan requires --stages")
    result = run_phase_a(env, args.stages, dry_run=args.dry_run, dry_run_id=args.dry_run_id)
    print(yaml_block(result))


def cmd_phase_b_build(env: Env, args: argparse.Namespace) -> None:
    print(yaml_block(build_provenance(env, args.candidate_id)))


def cmd_phase_c_render(env: Env, args: argparse.Namespace) -> None:
    print(yaml_block(render_phase_c(env, args.candidate_id)))


def cmd_set_status(env: Env, args: argparse.Namespace) -> None:
    manifest = load_manifest(env)
    entry = (manifest.get("candidates") or {}).get(args.candidate_id)
    if not entry:
        raise SystemExit(f"error: unknown candidate: {args.candidate_id}")
    current = entry.get("status")
    validate_manual_transition(current, args.status)
    transition(entry, args.status)
    if args.note:
        entry.setdefault("status_notes", []).append({"status": args.status, "note": args.note, "at": now_iso()})
    manifest["candidates"][args.candidate_id] = entry
    render_fit_file(env, manifest)
    render_batches(env, manifest)
    save_manifest(env, manifest)
    print(yaml_block({"candidate_id": args.candidate_id, "from": current, "to": args.status}))


def cmd_render_defense_prompt(env: Env, args: argparse.Namespace) -> None:
    print(render_defense_prompt(env, args.candidate_id, args.parameter))


def cmd_dry_run(env: Env, args: argparse.Namespace) -> None:
    stages = args.stages or ["003", "104", "105"]
    dry_run_id = args.dry_run_id or "stages_003_104_105"
    result = run_phase_a(env, stages, dry_run=True, dry_run_id=dry_run_id)
    manifest = load_manifest(env)
    candidate_ids = result["candidate_ids"]
    for cid in candidate_ids:
        build_provenance(env, cid)
        render_phase_c(env, cid)
    manifest = load_manifest(env)
    if dry_run_id in manifest.get("dry_runs", {}):
        manifest["dry_runs"][dry_run_id]["completed_at"] = now_iso()
        manifest["dry_runs"][dry_run_id]["last_phase"] = "phase_c_prompt_rendered"
        manifest["dry_runs"][dry_run_id]["candidate_ids"] = candidate_ids
        save_manifest(env, manifest)
    status = {
        "dry_run": True,
        "dry_run_id": dry_run_id,
        "stages": [f"{int(s):03d}" for s in stages],
        "candidate_ids": candidate_ids,
        "stage_results": result["stage_results"],
        "graph_gaps": result["graph_gaps"],
        "critic_prompt": result["critic_prompt"],
        "binding_verdicts_recorded": sum(
            1
            for entry in (manifest.get("candidates") or {}).values()
            if (entry.get("verdict") or {}).get("adversarial") or (entry.get("verdict") or {}).get("adjudication")
        ),
    }
    print(yaml_block(status))


def remove_path(path: Path) -> None:
    if not path.exists():
        return
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()


def collect_dry_run_ids(env: Env, manifest: dict[str, Any]) -> set[str]:
    ids = set(str(k) for k in (manifest.get("dry_runs") or {}).keys())
    for entry in (manifest.get("candidates") or {}).values():
        if entry.get("dry_run") and entry.get("dry_run_id"):
            ids.add(str(entry["dry_run_id"]))
    for path in (env.artifact_root / "phase_a_fragments").glob("*.yaml"):
        data = load_yaml(path)
        if data.get("dry_run") and data.get("dry_run_id"):
            ids.add(str(data["dry_run_id"]))
    for path in (env.artifact_root / "tmp_prompts").glob("*_phase_a_*.md"):
        if not artifact_has_dry_run_marker(path):
            continue
        name = path.name
        for suffix in ("_phase_a_numeric_literal.md", "_phase_a_claim_label.md", "_phase_a_graph.md", "_phase_a_existing_provenance.md"):
            if name.endswith(suffix):
                ids.add(name[: -len(suffix)])
    for directory in ("provenance", "reports", "defenses", "verdicts", "codex_logs", "tmp_prompts"):
        for path in (env.artifact_root / directory).glob("dryrun_*"):
            try:
                text = read_text(path)[:20000]
            except OSError:
                continue
            match = re.search(r"(?m)^\s*dry_run_id:\s*['\"]?([^'\"\s]+)", text)
            if match:
                ids.add(match.group(1))
    return ids


def artifact_has_dry_run_marker(path: Path) -> bool:
    if "codex_logs" in path.parts:
        return False
    if not path.is_file() or path.suffix not in {".yaml", ".yml", ".md", ".txt"}:
        return False
    try:
        text = read_text(path)[:20000]
    except OSError:
        return False
    lines = text.lstrip().splitlines()
    if lines and lines[0].strip() == "---":
        marker_lines = []
        for line in lines[1:80]:
            if line.strip() == "---":
                break
            marker_lines.append(line)
    else:
        marker_lines = lines[:40]
    marker_text = "\n".join(marker_lines)
    return bool(re.search(r"(?m)^\s*(dry_run:\s*true|dry_run_id:|non_binding:\s*true)\b", marker_text))


def cmd_purge_dry_run(env: Env, args: argparse.Namespace) -> None:
    manifest = load_manifest(env)
    target = args.dry_run_id
    targets = sorted(collect_dry_run_ids(env, manifest)) if target == "all" else [target]
    removed_candidates = []
    for cid, entry in list((manifest.get("candidates") or {}).items()):
        if not entry.get("dry_run"):
            continue
        if entry.get("dry_run_id") not in targets:
            continue
        for value in (entry.get("paths") or {}).values():
            if isinstance(value, list):
                for item in value:
                    remove_path(env.abs_from_rel(item))
            elif value:
                remove_path(env.abs_from_rel(value))
        removed_candidates.append(cid)
        del manifest["candidates"][cid]
    for dry_run_id in targets:
        manifest.get("dry_runs", {}).pop(dry_run_id, None)
    if target == "all":
        manifest["dry_runs"] = {}
    for directory in ("phase_a_fragments", "tmp_prompts", "provenance", "reports", "defenses", "verdicts", "codex_logs"):
        root = env.artifact_root / directory
        for dry_run_id in targets:
            for path in root.glob(f"{dry_run_id}_*"):
                remove_path(path)
            for path in root.glob(f"*{dry_run_id}*"):
                if path.name.startswith("dryrun_") or artifact_has_dry_run_marker(path):
                    remove_path(path)
        for cid in removed_candidates:
            for path in root.glob(f"{cid}*"):
                remove_path(path)
    if target == "all":
        for directory in ("phase_a_fragments", "tmp_prompts", "provenance", "reports", "defenses", "verdicts", "codex_logs"):
            root = env.artifact_root / directory
            if not root.exists():
                continue
            for path in root.iterdir():
                if artifact_has_dry_run_marker(path):
                    remove_path(path)
    data = ensure_benchmarks_file(env)
    data["entries"] = [
        e
        for e in data.get("entries", [])
        if not e.get("dry_run") or (target != "all" and e.get("dry_run_id") not in targets)
    ]
    write_yaml(env.benchmarks_path, data)
    render_fit_file(env, manifest)
    render_batches(env, manifest)
    save_manifest(env, manifest)
    print(yaml_block({"removed_dry_run": target, "enumerated_dry_run_ids": targets, "removed_candidates": removed_candidates}))


def cmd_artifact_root(env: Env, _args: argparse.Namespace) -> None:
    print(env.artifact_root)


def cmd_config_path(env: Env, args: argparse.Namespace) -> None:
    print(env.config_path_value(args.key))


def cmd_config_value(env: Env, args: argparse.Namespace) -> None:
    value = env.adv_value(args.key)
    print("" if value is None else value)


def cmd_codex_log_path(env: Env, args: argparse.Namespace) -> None:
    print(codex_log_path(env, args.candidate_id, args.parameter, args.iteration))


def cmd_parameter_session(env: Env, args: argparse.Namespace) -> None:
    print(parameter_session(env, args.candidate_id, args.parameter))


def cmd_record_codex_defense(env: Env, args: argparse.Namespace) -> None:
    print(
        yaml_block(
            record_codex_defense(
                env,
                args.candidate_id,
                args.parameter,
                args.iteration,
                args.return_code,
                args.log_path,
                args.new_session,
            )
        )
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Adversarial audit skill core helper")
    parser.add_argument("config_path", type=Path)
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("artifact-root").set_defaults(func=cmd_artifact_root)
    sub.add_parser("init").set_defaults(func=cmd_init)
    sub.add_parser("status").set_defaults(func=cmd_status)

    ci = sub.add_parser("candidate-info")
    ci.add_argument("candidate_id")
    ci.set_defaults(func=cmd_candidate_info)

    rp = sub.add_parser("render-phase-a-prompts")
    rp.add_argument("--stages", nargs="+", required=True)
    rp.add_argument("--prefix")
    rp.add_argument("--dry-run", action="store_true")
    rp.set_defaults(func=cmd_render_phase_a_prompts)

    pa = sub.add_parser("phase-a-scan")
    pa.add_argument("--stages", nargs="+", required=True)
    pa.add_argument("--dry-run", action="store_true")
    pa.add_argument("--dry-run-id")
    pa.set_defaults(func=cmd_phase_a_scan)

    pb = sub.add_parser("phase-b-build")
    pb.add_argument("candidate_id")
    pb.set_defaults(func=cmd_phase_b_build)

    pc = sub.add_parser("phase-c-render")
    pc.add_argument("candidate_id")
    pc.set_defaults(func=cmd_phase_c_render)

    ss = sub.add_parser("set-status")
    ss.add_argument("candidate_id")
    ss.add_argument("status", choices=STATE_MACHINE)
    ss.add_argument("--note")
    ss.set_defaults(func=cmd_set_status)

    rd = sub.add_parser("render-defense-prompt")
    rd.add_argument("candidate_id")
    rd.add_argument("parameter")
    rd.set_defaults(func=cmd_render_defense_prompt)

    dr = sub.add_parser("dry-run")
    dr.add_argument("--stages", nargs="+")
    dr.add_argument("--dry-run-id")
    dr.set_defaults(func=cmd_dry_run)

    purge = sub.add_parser("purge-dry-run")
    purge.add_argument("dry_run_id")
    purge.set_defaults(func=cmd_purge_dry_run)

    cp = sub.add_parser("config-path")
    cp.add_argument("key")
    cp.set_defaults(func=cmd_config_path)

    cv = sub.add_parser("config-value")
    cv.add_argument("key")
    cv.set_defaults(func=cmd_config_value)

    cl = sub.add_parser("codex-log-path")
    cl.add_argument("candidate_id")
    cl.add_argument("parameter")
    cl.add_argument("iteration")
    cl.set_defaults(func=cmd_codex_log_path)

    ps = sub.add_parser("parameter-session")
    ps.add_argument("candidate_id")
    ps.add_argument("parameter")
    ps.set_defaults(func=cmd_parameter_session)

    rc = sub.add_parser("record-codex-defense")
    rc.add_argument("candidate_id")
    rc.add_argument("parameter")
    rc.add_argument("iteration")
    rc.add_argument("return_code")
    rc.add_argument("log_path")
    rc.add_argument("new_session", nargs="?", default="")
    rc.set_defaults(func=cmd_record_codex_defense)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    require_manifest_lock(args.command)
    env = Env(args.config_path)
    args.func(env, args)


if __name__ == "__main__":
    main()
