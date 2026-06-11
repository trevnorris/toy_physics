#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable

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


_LOADER = yaml.CSafeLoader if hasattr(yaml, "CSafeLoader") else yaml.SafeLoader


class NoAliasDumper(yaml.SafeDumper):
    def ignore_aliases(self, data: Any) -> bool:
        return True


if hasattr(yaml, "CSafeDumper"):
    class NoAliasCDumper(yaml.CSafeDumper):
        def ignore_aliases(self, data: Any) -> bool:
            return True

    _C_DUMPER = NoAliasCDumper
else:
    _C_DUMPER = None

# CSafeDumper is not content-neutral for this campaign manifest: it double-quotes
# and escapes some astral Unicode scalars where SafeDumper emits single-quoted
# Unicode text. Keep artifact writes on the proven no-alias SafeDumper.
_DUMPER = NoAliasDumper


WRITE_COMMANDS = {
    "init",
    "phase-a-scan",
    "phase-a-ingest",
    "apply-alias-map",
    "phase-b-build",
    "phase-b-build-all",
    "phase-b-ingest",
    "benchmark-ingest",
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

INGEST_MODALITIES = MODALITIES + ["completeness_critic"]

ALIAS_ADJACENT_LINE_WINDOW = 3

TARGET_LAYER_DEFAULT = "provenance/_target_layer.yaml"
CONCEPT_ALIASES_DEFAULT = "provenance/_concept_aliases.yaml"

PHASE_B_CONSTRAINT_KINDS = {
    "internal_consistency",
    "published_target",
    "free_choice",
}

BENCHMARK_SOURCE_TYPES = {
    "web_lookup",
    "textbook",
    "CODATA",
}

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


PHASE_B_FILENAME_STEM_MAX = 220


def bounded_phase_b_stem(candidate_id: str, parameter_name: str, marker: str) -> str:
    candidate_slug = slug(candidate_id)
    parameter_slug = slug(parameter_name)
    stem = f"{candidate_slug}{marker}{parameter_slug}"
    if len(stem) <= PHASE_B_FILENAME_STEM_MAX:
        return stem
    digest = hashlib.sha1(f"{candidate_id}\0{parameter_name}".encode("utf-8")).hexdigest()[:12]
    prefix = f"{candidate_slug}{marker}"
    available = PHASE_B_FILENAME_STEM_MAX - len(prefix) - len(digest) - 1
    if available >= 16:
        shortened_parameter = parameter_slug[:available].rstrip("_") or parameter_slug[:available]
        return f"{prefix}{shortened_parameter}_{digest}"

    available = PHASE_B_FILENAME_STEM_MAX - len(marker) - len(digest) - 2
    candidate_available = max(16, min(len(candidate_slug), available // 2))
    parameter_available = max(16, available - candidate_available)
    shortened_candidate = candidate_slug[:candidate_available].rstrip("_") or candidate_slug[:candidate_available]
    shortened_parameter = parameter_slug[:parameter_available].rstrip("_") or parameter_slug[:parameter_available]
    return f"{shortened_candidate}{marker}{shortened_parameter}_{digest}"


def phase_b_prompt_path(env: Env, candidate_id: str, parameter_name: str) -> Path:
    return env.artifact_root / "tmp_prompts" / f"{bounded_phase_b_stem(candidate_id, parameter_name, '__phase_b__')}.md"


def phase_b_bundle_path(env: Env, candidate_id: str, parameter_name: str) -> Path:
    return env.artifact_root / "provenance" / f"{bounded_phase_b_stem(candidate_id, parameter_name, '__')}.yaml"


def load_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with path.open(encoding="utf-8") as f:
        return yaml.load(f, Loader=_LOADER) or {}


def write_yaml(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        yaml.dump(data, f, Dumper=_DUMPER, default_flow_style=False, sort_keys=False, width=120, allow_unicode=True)
    tmp.replace(path)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def require_manifest_lock(command: str, *, writes_manifest: bool | None = None) -> None:
    writes = command in WRITE_COMMANDS if writes_manifest is None else writes_manifest
    if writes and os.environ.get("ADVERSARIAL_MANIFEST_LOCKED") != "1":
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


def normalize_family_param(text: Any) -> str:
    value = str(text or "").strip()
    value = value.replace("\\", "")
    value = value.replace("{", "").replace("}", "")
    value = value.replace("^", "_")
    value = clean_param_piece(value)
    return value.replace("_", "").lower()


def normalized_parameter_family(parameter_names: Iterable[Any]) -> list[str]:
    family = sorted({p for p in (normalize_family_param(param) for param in parameter_names) if p})
    return family or ["unclassified_target_parameter"]


def normalized_parameter_family_key(parameter_names: Iterable[Any]) -> str:
    return "+".join(normalized_parameter_family(parameter_names))


DEFAULT_CONCEPT_ALIAS_DATA: dict[str, Any] = {
    "schema_version": 1,
    "description": (
        "Reviewable concept-variant aliases for adversarial target resolution. "
        "Variants map to one canonical concept id; seeded chains are non-primary family overlays."
    ),
    "aliases": {
        "Sigma_0_can": "Sigma0_can",
        "Sigma0_can": "Sigma0_can",
        "Sigma0_num": "Sigma0_can",
        "Sigma0_can_expected": "Sigma0_can",
        "Sigma_0_can_expected": "Sigma0_can",
        "Sigma_0_star": "Sigma0_star",
        "Sigma0_hat": "Sigma0_hat",
        "chi_q": "chi_Q",
        "chiQ": "chi_Q",
        "chi_Q": "chi_Q",
        "chi_Q_R": "chi_Q",
        "chi_Q_hyb": "chi_Q",
        "m0_hat": "mhat_0",
        "mhat0": "mhat_0",
        "mhat_0": "mhat_0",
        "m_hat_0": "mhat_0",
        "m_hat0": "mhat_0",
        "mhat_minus": "mhat_0",
        "P0target": "P0_target",
        "P_0_target": "P0_target",
        "P0_target": "P0_target",
        "P_0": "P0_target",
        "P0_compact": "P0_target",
        "NQ_target": "N_Q",
        "N_Q_target": "N_Q",
        "NQ_from_def": "N_Q",
        "N_Q": "N_Q",
        "Gamma5_target": "Gamma_5",
        "Gamma5_port": "Gamma_5",
        "Gamma_5": "Gamma_5",
        "gamma_GR": "Gamma_5",
        "normalization_54_5": "Gamma_5",
        "coefficient_54_over_5": "Gamma_5",
        "target_coefficient_54_5": "Gamma_5",
        "L_over_a": "L_over_a",
        "Lambda_ell": "L_over_a",
        "ell_over_a": "L_over_a",
        "aspect_ratio": "L_over_a",
        "aspect_ratio_37_20": "L_over_a",
        "aspect_ratio_37_over_20": "L_over_a",
        "rF1": "r_F1",
        "r_F1": "r_F1",
        "r_f1": "r_F1",
        "r_star": "r_F1",
        "mathfrak_r_F1": "r_F1",
        "mathfrak_r_geom": "r_F1",
        "g_F1": "g_star",
        "g_minus_F1": "g_star",
        "mathfrak_g_minus_F1": "g_star",
        "T_can": "T_can",
        "T_m_can": "T_can",
        "T_hat_can": "T_can",
        "T_hat_can_expected": "T_can",
        "T_hat_m": "T_can",
        "S_can": "S_can",
        "S_can_expected": "S_can",
        "Pi_can": "Pi_can",
        "Pi_can_expected": "Pi_can",
        "V_known": "V_known",
        "V_known_x1": "V_known",
        "DeltaV_req": "DeltaV_req",
        "lambda_L": "lambda_L",
        "lambda_L_soft": "lambda_L",
        "lambda_L_paper": "lambda_L",
        "Lvar_soft": "lambda_L",
        "K_eta": "K_eta",
        "T_Omega": "T_Omega",
        "T_w": "T_w",
        "mu_eta": "mu_eta",
        "K_Sigma": "K_Sigma",
        "Theta_w": "Theta_w",
        "Theta_w_chi": "Theta_w",
        "Theta_w_J": "Theta_w",
    },
    "key_rules": [
        {"pattern": r"(?:^|_)sigma0_num(?:_|$)", "target": "Sigma0_num"},
        {"pattern": r"(?:^|_)sigma0_can(?:_|$)", "target": "Sigma0_can"},
        {"pattern": r"(?:^|_)sigma_0_can(?:_|$)", "target": "Sigma_0_can"},
        {"pattern": r"fixedpoint_expected_targets_json", "target": "Sigma0_can", "multi_target": True},
        {"pattern": r"canonical_quartet_literals", "target": "Sigma0_can", "multi_target": True},
        {"pattern": r"canonical_mouth_block", "target": "Sigma0_can", "multi_target": True},
        {"pattern": r"(?:^|_)chi_?q(?:_|$)|chiq|chiQ", "target": "chi_Q"},
        {"pattern": r"unit_product_m0hat2_chiQ_NQ", "target": "chi_Q", "multi_target": True},
        {"pattern": r"(?:^|_)m0hat|mhat0|mhat_0|m_hat_0", "target": "mhat_0"},
        {"pattern": r"(?:^|_)nq(?:_|$)|nq_target|NQ_target", "target": "N_Q"},
        {"pattern": r"(?:^|_)p0target|p0_target|p_0_target|universal_p0", "target": "P0_target"},
        {"pattern": r"54_5|54_over_5|quadrupole_normalization|normalization_target", "target": "P0_target"},
        {"pattern": r"aspect_ratio|37_20|37_over_20", "target": "L_over_a"},
        {"pattern": r"rexact_closed_form_4107|radical_constant_forced_by_rf1", "target": "r_F1"},
        {"pattern": r"(?:^|_)lambda_l(?:_|$)|lambda_L", "target": "lambda_L"},
        {"pattern": r"vknown|V_known|barrier_benchmark", "target": "V_known"},
        {"pattern": r"wall_action_constitutive_coefficients", "target": "K_eta", "multi_target": True},
        {"pattern": r"calibration_slice|benchmark_masses|stage252_slice_inputs", "target": "mu_eta", "multi_target": True},
    ],
    "multi_target_key_patterns": [
        r"block",
        r"bundle",
        r"calibration_slice",
        r"canonical_quartet",
        r"constants_block",
        r"expected_targets",
        r"inputs",
        r"packet",
        r"quintuple",
        r"slice",
    ],
    "chains": {
        "chain_quad_54_5": {
            "description": "54/5 plus gamma_GR quadrupole normalization target.",
            "concepts": ["P0_target", "N_Q", "Gamma_5", "mhat_0"],
            "expected_stages": ["019", "022", "023", "025", "030", "189", "193", "195", "197", "223"],
            "expected_candidate_ids": [
                "fit_stage019_P0_target_upstream_pin",
                "fit_stage019_p0target_54_5_raw_retype",
                "fit_stage022_p0_quadrupole_normalization_target",
                "fit_stage023_nq_target_gr_quadrupole_match",
                "fit_stage025_target_coefficient_54_5",
                "fit_stage030_nq_target_retype_54_5",
                "fit_stage189_mhat0_quadrupole_normalization_target",
                "fit_stage189_normalization_target_2g_5c5",
                "fit_stage195_p0_target_canonical_value",
                "fit_stage197_p0_target_quadrupole_normalization",
                "fit_stage223_p0_target_calibration_coefficient",
            ],
        },
        "chain_aspect_37_20": {
            "description": "37/20 aspect-ratio cascade into Family-1, including the 4107 radical.",
            "concepts": ["L_over_a", "r_F1", "g_star", "epsilon_r", "eta"],
            "expected_stages": ["073", "121", "131", "139", "146", "148", "168"],
            "expected_candidate_ids": [
                "fit_stage073_Lambda_ell_fixed_37",
                "fit_stage073_aspect_ratio_L_over_a",
                "fit_stage073_wall_fraction_epsilon_r",
                "fit_stage121_aspect_ratio_37_20_reference_freeze",
                "fit_stage121_aspect_ratio_37_over_20",
                "fit_stage146_aspect_ratio_37_20_primitive",
                "fit_stage148_radical_constant_forced_by_rf1",
                "fit_stage168_rexact_closed_form_4107",
            ],
        },
        "chain_chi_Q_norm": {
            "description": "chi_Q / mhat_0 / N_Q outgoing normalization chain.",
            "concepts": ["chi_Q", "mhat_0", "N_Q", "sigma_Q_can", "xi_Q"],
            "expected_stages": ["100", "101", "104", "105", "107", "108", "192", "194", "195", "196", "197", "232"],
            "expected_candidate_ids": [
                "fit_stage100_chiQ_fixed_downstream_card",
                "fit_stage104_unit_product_m0hat2_chiQ_NQ",
                "fit_stage105_chiQ_eq_1_by_matching_card",
                "fit_stage105_chi_q_canonical_fix",
                "fit_stage105_chi_q_canonical_unity",
                "fit_stage_105_chi_q",
                "fit_stage194_chi_Q_canonical_fixing",
                "fit_stage195_mhat0_natural_source_map",
                "fit_stage196_sigma_Q_canonical_outgoing_scale",
                "fit_stage197_chi_q_unity_packet_a_target",
                "fit_stage232_nq_canonical_normalization",
            ],
        },
        "chain_sigma0_transport": {
            "description": "Sigma0_can / T_can transport chain, including 867/876 digit drift.",
            "concepts": ["Sigma0_can", "T_can", "S_can", "Pi_can"],
            "expected_stages": ["155", "156", "157", "158", "163", "168"],
            "expected_candidate_ids": [
                "fit_stage155_fixedpoint_expected_targets_json",
                "fit_stage156_unique_traction_renormalization",
                "fit_stage157_expected_values_sidecar_json",
                "fit_stage158_canonical_quartet_literals",
                "fit_stage163_sigma0_can_digit_variant",
                "fit_stage168_canonical_mouth_block",
                "fit_stage168_sigma0_num_digit_variant",
            ],
        },
        "chain_barrier_222_224": {
            "description": "Barrier 222-224 chain, including V_known.",
            "concepts": ["V_known", "DeltaV_req", "epsilon_barrier", "lambda_20", "lambda_21", "lambda_22", "barP0_compat", "T_quad"],
            "expected_stages": ["222", "223", "224"],
            "expected_candidate_ids": [
                "fit_stage222_illustrative_barrier_benchmark",
                "fit_stage222_vknown_barrier_benchmark",
                "fit_stage223_barrier_benchmark_and_eta_window",
                "fit_stage224_grouped_signature_exact",
                "fit_stage224_kill_test_budgets_slice_anchored",
                "fit_stage224_t_quad_target_scale",
            ],
        },
        "chain_calibration_245_253": {
            "description": "Calibration 245-253 chain, including lambda_L and CODATA mass-ratio carriers.",
            "concepts": ["lambda_L", "epsilon_eta", "K_turn", "f_lat", "gamma_lattice_red", "gamma_lattice_legacy", "mu_eta", "Upsilon_lat", "m_s"],
            "expected_stages": ["245", "247", "248", "250", "253"],
            "expected_candidate_ids": [
                "fit_stage245_eps_eta_session1_match",
                "fit_stage247_lambda_L_backsolve",
                "fit_stage247_lambda_L_closure",
                "fit_stage247_lambda_l_fixed_by_recorded_point",
                "fit_stage248_xi_turn_lambda_th_carried_hardcodes",
                "fit_stage250_benchmark_masses_and_window",
                "fit_stage253_calibration_slice_nominal_constants",
                "fit_stage253_K_turn_force_match",
                "fit_stage253_upsilon_lat_calibration",
            ],
        },
        "chain_wall_action": {
            "description": "Wall-action coefficients plus mhat/P0 normalization carriers.",
            "concepts": ["K_eta", "T_Omega", "T_w", "mu_eta", "mhat_0", "P0_target", "K_Sigma"],
            "expected_stages": ["001", "019", "022", "038", "193", "223"],
            "expected_candidate_ids": [
                "fit_stage001_wall_action_constitutive_coefficients",
                "fit_stage019_K_Sigma_fixed_by_one_pole",
                "fit_stage019_K_Sigma_fixed_by_outgoing_normalization",
                "fit_stage019_P0_target_upstream_pin",
                "fit_stage022_mhat0_unity_branch",
                "fit_stage022_p0_quadrupole_normalization_target",
                "fit_stage038_mhat_p0_fitted_scales",
                "fit_stage193_p0_target_surface_frozen_without_derivation_edge",
                "fit_stage223_p0_target_import_unprovenanced",
                "fit_stage223_universal_p0_target",
            ],
        },
    },
}


def concept_alias_path(env: Env) -> Path:
    return env.artifact_root / CONCEPT_ALIASES_DEFAULT


def target_layer_path(env: Env) -> Path:
    return env.artifact_root / TARGET_LAYER_DEFAULT


def ensure_concept_alias_file(env: Env) -> Path:
    path = concept_alias_path(env)
    if not path.exists():
        write_yaml(path, DEFAULT_CONCEPT_ALIAS_DATA)
    return path


def load_concept_aliases(env: Env) -> dict[str, Any]:
    path = ensure_concept_alias_file(env)
    data = load_yaml(path)
    if not isinstance(data, dict):
        raise SystemExit(f"error: concept alias table must be a YAML mapping: {env.rel(path)}")
    data.setdefault("aliases", {})
    data.setdefault("key_rules", [])
    data.setdefault("chains", {})
    return data


GREEK_NAME_MAP = {
    "alpha": "alpha",
    "beta": "beta",
    "chi": "chi",
    "delta": "delta",
    "Delta": "Delta",
    "epsilon": "epsilon",
    "eta": "eta",
    "Gamma": "Gamma",
    "gamma": "gamma",
    "kappa": "kappa",
    "lambda": "lambda",
    "Lambda": "Lambda",
    "mu": "mu",
    "Omega": "Omega",
    "omega": "omega",
    "Pi": "Pi",
    "rho": "rho",
    "Sigma": "Sigma",
    "sigma": "sigma",
    "Theta": "Theta",
    "theta": "theta",
    "xi": "xi",
    "zeta": "zeta",
}


def normalize_symbol_match_key(text: Any) -> str:
    value = str(text or "")
    value = value.replace("\\widehat", "widehat")
    value = value.replace("\\hat", "hat")
    for tex, name in GREEK_NAME_MAP.items():
        value = value.replace("\\" + tex, name)
    value = value.replace("{", "").replace("}", "")
    value = value.replace("^", "_")
    value = value.replace("-", "_")
    value = re.sub(r"[^A-Za-z0-9]+", "_", value)
    return re.sub(r"_+", "_", value).strip("_").lower()


def concept_alias_lookup(alias_data: dict[str, Any]) -> dict[str, str]:
    cached = alias_data.get("_compiled_alias_lookup")
    if isinstance(cached, dict):
        return cached
    lookup: dict[str, str] = {}
    for variant, canonical in (alias_data.get("aliases") or {}).items():
        canonical_text = str(canonical)
        lookup[str(variant)] = canonical_text
        lookup[normalize_symbol_match_key(variant)] = canonical_text
        lookup[normalize_family_param(variant)] = canonical_text
        lookup[canonical_text] = canonical_text
        lookup[normalize_symbol_match_key(canonical_text)] = canonical_text
        lookup[normalize_family_param(canonical_text)] = canonical_text
    alias_data["_compiled_alias_lookup"] = lookup
    return lookup


def canonicalize_concept(value: Any, alias_data: dict[str, Any]) -> str:
    text = str(value or "").strip()
    if not text:
        return "unclassified_target_parameter"
    lookup = concept_alias_lookup(alias_data)
    return (
        lookup.get(text)
        or lookup.get(normalize_symbol_match_key(text))
        or lookup.get(normalize_family_param(text))
        or text
    )


def aliases_for_canonical(canonical: str, alias_data: dict[str, Any]) -> list[str]:
    cache = alias_data.setdefault("_compiled_aliases_by_canonical", {})
    if isinstance(cache, dict) and canonical in cache:
        return list(cache[canonical])
    variants = {canonical}
    lookup = concept_alias_lookup(alias_data)
    for variant, target in (alias_data.get("aliases") or {}).items():
        if (lookup.get(str(target)) or lookup.get(normalize_symbol_match_key(target)) or str(target)) == canonical:
            variants.add(str(variant))
    result = sorted(variants)
    if isinstance(cache, dict):
        cache[canonical] = result
    return result


def target_search_variants(canonical: str, raw: str, alias_data: dict[str, Any]) -> list[str]:
    variants = set(aliases_for_canonical(canonical, alias_data))
    variants.add(raw)
    variants.add(canonical)
    expanded = set()
    for variant in variants:
        value = str(variant)
        expanded.add(value)
        expanded.add(value.replace("_", ""))
        expanded.add(value.replace("_", "\\_"))
        expanded.add(value.replace("_", " "))
        if value == "L_over_a":
            expanded.update({"L/a", "L}{a", "L\\over a", "\\frac{L}{a}"})
        if value in {"mhat_0", "m0_hat"}:
            expanded.update({"mhat0", "m_hat_0", "m_{\\hat 0}", "\\widehat m_0", "\\widehat m_0^{\\,2}"})
        if value == "chi_Q":
            expanded.update({"chiQ", "\\chi_Q", "\\chi_{Q}"})
        if value == "Sigma0_can":
            expanded.update({"Sigma0", "Sigma0_can", "Sigma_0_can", "\\Sigma_0", "\\Sigma_0^{\\rm can}"})
        if value == "P0_target":
            expanded.update({"P0target", "P_0", "P0", "P_0_target"})
        if value == "N_Q":
            expanded.update({"NQ", "NQ_target", "N_Q_target", "N_Q"})
    return sorted({item for item in expanded if item}, key=lambda item: (-len(item), item))


def strip_stage_key_prefix(candidate_key: Any) -> str:
    key = str(candidate_key or "")
    key = re.sub(r"^fit_", "", key)
    key = re.sub(r"^stage_?0*\d+_?", "", key)
    return key


def plausible_parameter_names(entry: dict[str, Any], alias_data: dict[str, Any]) -> list[str]:
    values = []
    for param in entry.get("parameter_names") or []:
        text = str(param or "").strip()
        if not text:
            continue
        if text in {"matched_fingerprint_value", "stale_output", "assert_nonzero", "assert_zero", "paragraph"}:
            continue
        canonical = canonicalize_concept(text, alias_data)
        if canonical not in values:
            values.append(canonical)
    return values


def key_rule_target(candidate_key: str, alias_data: dict[str, Any]) -> tuple[str | None, str | None, bool]:
    key = strip_stage_key_prefix(candidate_key)
    for rule in alias_data.get("key_rules") or []:
        pattern = str(rule.get("pattern") or "")
        if pattern and re.search(pattern, key, re.IGNORECASE):
            target = str(rule.get("target") or "")
            if target:
                return target, f"candidate_key rule /{pattern}/", bool(rule.get("multi_target"))
    return None, None, False


def parse_primary_target(entry: dict[str, Any], alias_data: dict[str, Any]) -> dict[str, Any]:
    candidate_key = str(entry.get("candidate_key") or entry.get("id") or "")
    key_suffix = strip_stage_key_prefix(candidate_key)
    raw_target, basis, rule_multi = key_rule_target(candidate_key, alias_data)
    resolution_confidence = "high"
    low_reasons: list[str] = []
    if raw_target:
        raw = raw_target
    else:
        normalized_key = normalize_symbol_match_key(key_suffix)
        matches = []
        for param in entry.get("parameter_names") or []:
            param_text = str(param or "")
            candidate_forms = {normalize_symbol_match_key(param_text), normalize_family_param(param_text)}
            canonical = canonicalize_concept(param_text, alias_data)
            candidate_forms.add(normalize_symbol_match_key(canonical))
            candidate_forms.add(normalize_family_param(canonical))
            if any(form and re.search(rf"(?:^|_){re.escape(form)}(?:_|$)", normalized_key) for form in candidate_forms):
                matches.append((param_text, canonical))
        canonical_matches = []
        for _raw, canonical in matches:
            if canonical not in canonical_matches:
                canonical_matches.append(canonical)
        if len(canonical_matches) == 1:
            raw = matches[0][0]
            basis = "candidate_key matched manifest parameter"
        else:
            params = plausible_parameter_names(entry, alias_data)
            if params:
                raw = params[0]
                basis = "parameter_names fallback because candidate_key was uninformative"
                resolution_confidence = "low"
                low_reasons.append("key_uninformative")
                if len(params) > 1:
                    low_reasons.append("multiple_plausible_targets")
            else:
                raw = "unclassified_target_parameter"
                basis = "no resolvable candidate_key or parameter_names target"
                resolution_confidence = "low"
                low_reasons.append("no_plausible_target")
    primary = canonicalize_concept(raw, alias_data)
    additional: list[str] = []
    key_patterns = alias_data.get("multi_target_key_patterns") or []
    key_multi = rule_multi or any(re.search(str(pattern), key_suffix, re.IGNORECASE) for pattern in key_patterns)
    if key_multi:
        for param in plausible_parameter_names(entry, alias_data):
            canonical = canonicalize_concept(param, alias_data)
            if canonical != primary and canonical not in additional:
                additional.append(canonical)
    return {
        "raw_primary_target_parameter": raw,
        "primary_target_parameter": primary,
        "resolution_basis": basis or "candidate_key parse",
        "target_resolution_confidence": resolution_confidence,
        "low_reasons": low_reasons,
        "multi_target": bool(additional),
        "additional_target_parameters": additional,
    }


def normalize_literal_value(text: str) -> str:
    value = str(text).strip()
    frac = re.fullmatch(r"\\frac\{([^{}]+)\}\{([^{}]+)\}", value)
    if frac:
        value = f"{frac.group(1)}/{frac.group(2)}"
    value = value.replace(" ", "")
    value = value.replace("\\,", "")
    value = value.replace("{", "").replace("}", "")
    value = value.replace("−", "-")
    return value


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def normalize_decimal_text(text: str) -> str:
    value = text.strip().replace("_", "").replace("−", "-")
    try:
        decimal_value = Decimal(value)
    except InvalidOperation:
        return value
    if decimal_value == decimal_value.to_integral_value():
        return str(decimal_value.quantize(Decimal(1)))
    frac = Fraction(decimal_value)
    if frac.denominator <= 1000:
        return fraction_text(frac)
    return format(decimal_value.normalize(), "f").rstrip("0").rstrip(".")


def normalize_extracted_value(raw: str) -> str:
    value = normalize_literal_value(raw).strip().strip("`")
    value = value.replace("'", '"')
    sp_float = re.fullmatch(r"""sp\.Float\(\s*"([^"]+)"(?:\s*,\s*\d+)?\s*\)""", value)
    if sp_float:
        return normalize_decimal_text(sp_float.group(1))
    sp_int = re.fullmatch(r"sp\.Integer\(\s*([-+]?\d+)\s*\)", value)
    if sp_int:
        return str(int(sp_int.group(1)))
    sp_rat = re.fullmatch(r"sp\.Rational\(\s*([-+]?\d+)\s*,\s*([-+]?\d+)\s*\)", value)
    if sp_rat:
        return fraction_text(Fraction(int(sp_rat.group(1)), int(sp_rat.group(2))))
    m_rat = re.fullmatch(r"Rational\[\s*([-+]?\d+)\s*,\s*([-+]?\d+)\s*\]", value)
    if m_rat:
        return fraction_text(Fraction(int(m_rat.group(1)), int(m_rat.group(2))))
    frac = re.fullmatch(r"([-+]?\d+)\s*/\s*([-+]?\d+)", value)
    if frac:
        return fraction_text(Fraction(int(frac.group(1)), int(frac.group(2))))
    if re.fullmatch(r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?", value):
        return normalize_decimal_text(value)
    if re.fullmatch(r"[-+]?\d+(?:[eE][-+]?\d+)?", value):
        return normalize_decimal_text(value)
    value = re.sub(r"\s+", "", value)
    return value


SP_FLOAT_RE = re.compile(r"""sp\.Float\(\s*["']([^"']+)["'](?:\s*,\s*\d+)?\s*\)""")
SP_RATIONAL_RE = re.compile(r"sp\.Rational\(\s*([-+]?\d+)\s*,\s*([-+]?\d+)\s*\)")
SP_INTEGER_RE = re.compile(r"sp\.Integer\(\s*([-+]?\d+)\s*\)")
M_RATIONAL_RE = re.compile(r"Rational\[\s*([-+]?\d+)\s*,\s*([-+]?\d+)\s*\]")
CLOSED_FORM_RE = re.compile(
    r"(?:Sqrt\[[^\]\n]+\]|sqrt\([^\)\n]+\)|sp\.sqrt\([^\)\n]+\))(?:\s*/\s*(?:\([^)\n]+\)|[A-Za-z0-9_.*^+-]+))?"
)
VALUE_EXPR_RE = re.compile(
    r"(sp\.Float\(\s*['\"][^'\"]+['\"](?:\s*,\s*\d+)?\s*\)|sp\.Rational\(\s*[-+]?\d+\s*,\s*[-+]?\d+\s*\)|sp\.Integer\(\s*[-+]?\d+\s*\)|Rational\[\s*[-+]?\d+\s*,\s*[-+]?\d+\s*\]|\\frac\{[-+]?\d+\}\{[-+]?\d+\}|[-+]?\d+\s*/\s*[-+]?\d+|[-+]?\d+\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+(?:[eE][-+]?\d+)?)"
)


def add_value_candidate(values: list[dict[str, str]], raw: str, rule: str) -> None:
    normalized = normalize_extracted_value(raw)
    if not normalized:
        return
    item = {"raw": str(raw).strip(), "normalized": normalized, "comparison_key": normalized, "extraction_rule": rule}
    if item not in values:
        values.append(item)


def extract_values_from_rhs(rhs: str) -> list[dict[str, str]]:
    values: list[dict[str, str]] = []
    masked = list(rhs)

    def mask_span(start: int, end: int) -> None:
        for i in range(start, end):
            masked[i] = " "

    for regex, rule, formatter in (
        (SP_FLOAT_RE, "sp.Float", lambda m: f'sp.Float("{m.group(1)}")'),
        (SP_RATIONAL_RE, "sp.Rational", lambda m: f"sp.Rational({m.group(1)},{m.group(2)})"),
        (SP_INTEGER_RE, "sp.Integer", lambda m: f"sp.Integer({m.group(1)})"),
        (M_RATIONAL_RE, "Rational[]", lambda m: f"Rational[{m.group(1)},{m.group(2)}]"),
    ):
        for match in regex.finditer(rhs):
            add_value_candidate(values, formatter(match), rule)
            mask_span(match.start(), match.end())

    masked_text = "".join(masked)
    for match in CLOSED_FORM_RE.finditer(masked_text):
        add_value_candidate(values, match.group(0), "closed_form")
        mask_span(match.start(), match.end())
    masked_text = "".join(masked)

    frac_54_5 = re.search(r"(?<!\d)54(?:\s*\*[^/\n;]+)?/\s*\(?\s*5(?!\d)", masked_text)
    if frac_54_5:
        add_value_candidate(values, "54/5", "coefficient_fraction_54_5")
    for match in re.finditer(r"(?:^|=|:|->|\\to|\\mapsto)\s*" + VALUE_EXPR_RE.pattern, masked_text):
        add_value_candidate(values, match.group(1), "relation_value")
    for match in re.finditer(r"\\frac\{[-+]?\d+\}\{[-+]?\d+\}", masked_text):
        add_value_candidate(values, match.group(0), "tex_fraction")
    return values


def extract_leading_values_from_rhs(rhs: str) -> list[dict[str, str]]:
    values: list[dict[str, str]] = []
    text = str(rhs or "").strip()
    if not text:
        return values
    frac_54_5 = re.search(r"(?<!\d)54(?:\s*\*[^/\n;]+)?/\s*\(?\s*5(?!\d)", text[:160])
    if frac_54_5:
        add_value_candidate(values, "54/5", "coefficient_fraction_54_5")
        return values
    closed = CLOSED_FORM_RE.match(text)
    if closed:
        add_value_candidate(values, closed.group(0), "closed_form")
        text = text[closed.end() :]
    else:
        first = VALUE_EXPR_RE.match(text)
        if first:
            add_value_candidate(values, first.group(1), "relation_value")
            text = text[first.end() :]
    while True:
        chained = re.match(r"^\s*=\s*" + VALUE_EXPR_RE.pattern, text)
        if not chained:
            break
        add_value_candidate(values, chained.group(1), "chained_relation_value")
        text = text[chained.end() :]
    return values


def line_mentions_target(text: str, canonical: str, raw: str, alias_data: dict[str, Any]) -> bool:
    normalized_line = normalize_symbol_match_key(text)
    for variant in target_search_variants(canonical, raw, alias_data):
        key = normalize_symbol_match_key(variant)
        if key and re.search(rf"(?:^|_){re.escape(key)}(?:_|$)", normalized_line):
            return True
    return False


def lhs_matches_target(lhs: str, canonical: str, alias_data: dict[str, Any]) -> bool:
    return canonicalize_concept(lhs.strip().strip('"').strip("'"), alias_data) == canonical


def extract_target_values_from_text(text: str, canonical: str, raw: str, alias_data: dict[str, Any]) -> list[dict[str, str]]:
    source = str(text or "")
    values: list[dict[str, str]] = []
    for match in re.finditer(r"""["']([^"']+)["']\s*:\s*([^,\n#]+)""", source):
        if lhs_matches_target(match.group(1), canonical, alias_data):
            append_unique_items(values, extract_leading_values_from_rhs(match.group(2)))
    assignment_re = re.compile(r"([\\A-Za-z][\\A-Za-z0-9_{}^]*)\s*(?::=|=|:)\s*([^#;\n]+)")
    for match in assignment_re.finditer(source):
        lhs = match.group(1)
        if lhs_matches_target(lhs, canonical, alias_data) or line_mentions_target(lhs, canonical, raw, alias_data):
            append_unique_items(values, extract_leading_values_from_rhs(match.group(2)))
    if line_mentions_target(source, canonical, raw, alias_data):
        for match in re.finditer(r"(?P<prefix>.{0,120}?)(?:=|:=|:|->|\\to|\\mapsto)\s*(?P<rhs>[^#;\n]+)", source):
            if line_mentions_target(match.group("prefix"), canonical, raw, alias_data):
                append_unique_items(values, extract_leading_values_from_rhs(match.group("rhs")))
        if any(marker in strip_stage_key_prefix(raw).lower() for marker in ("37_20", "54_5")):
            append_unique_items(values, extract_values_from_rhs(source))
    return values


VALUE_TOKEN_RE = r"(\\frac\{[^{}]+\}\{[^{}]+\}|[-+]?\d+\s*/\s*[-+]?\d+|[-+]?\d+\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+(?:[eE][-+]?\d+)?)"


def extract_literal_values(text: Any) -> list[str]:
    source = str(text or "")
    values: list[str] = []
    for match in re.finditer(rf"(?:=|:=|\\to|->|\\mapsto)\s*{VALUE_TOKEN_RE}", source):
        value = normalize_literal_value(match.group(1))
        if value not in values:
            values.append(value)
    for match in re.finditer(r"\\frac\{[^{}]+\}\{[^{}]+\}", source):
        value = normalize_literal_value(match.group(0))
        if value not in values:
            values.append(value)
    for match in re.finditer(r"(?<![A-Za-z0-9_])[-+]?\d+\s*/\s*[-+]?\d+(?![A-Za-z0-9_])", source):
        value = normalize_literal_value(match.group(0))
        if value not in values:
            values.append(value)
    return values


def candidate_literal_values(entry: dict[str, Any]) -> list[str]:
    values: list[str] = []
    for citation in entry.get("file_line_citations") or []:
        append_unique_items(values, extract_literal_values(citation.get("excerpt")))
    for fragment in entry.get("modality_fragments") or []:
        citation = fragment.get("citation") or {}
        append_unique_items(values, extract_literal_values(citation.get("excerpt")))
    return sorted(values)


def resolve_citation_source_path(env: Env, path_text: Any) -> Path:
    text = str(path_text or "").strip()
    p = Path(text)
    if p.is_absolute():
        return p
    candidates = [
        env.project_root / text,
        env.repo_root / text,
    ]
    if text.startswith("research/pde_ledger/"):
        candidates.insert(0, env.repo_root / text)
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


SOURCE_LINE_CACHE: dict[Path, list[str] | None] = {}


def read_citation_source_line(env: Env, citation: dict[str, Any]) -> str | None:
    line_no = citation_line(citation)
    if line_no is None:
        return None
    source_path = resolve_citation_source_path(env, citation.get("path"))
    if source_path not in SOURCE_LINE_CACHE:
        try:
            SOURCE_LINE_CACHE[source_path] = read_text(source_path).splitlines()
        except OSError:
            SOURCE_LINE_CACHE[source_path] = None
    lines = SOURCE_LINE_CACHE.get(source_path)
    if lines is None:
        return None
    if 1 <= line_no <= len(lines):
        return lines[line_no - 1].strip()
    return None


def target_value_sources(env: Env, entry: dict[str, Any]) -> list[dict[str, Any]]:
    sources: list[dict[str, Any]] = []
    for citation in entry.get("file_line_citations") or []:
        if citation.get("excerpt"):
            sources.append(
                {
                    "source_kind": "citation_excerpt",
                    "path": citation.get("path"),
                    "line": citation.get("line"),
                    "text": str(citation.get("excerpt")),
                }
            )
        source_line = read_citation_source_line(env, citation)
        if source_line and source_line != citation.get("excerpt"):
            sources.append(
                {
                    "source_kind": "cited_source_line",
                    "path": citation.get("path"),
                    "line": citation.get("line"),
                    "text": source_line,
                }
            )
    for fragment in entry.get("modality_fragments") or []:
        citation = fragment.get("citation") or {}
        if citation.get("excerpt"):
            sources.append(
                {
                    "source_kind": "modality_fragment_excerpt",
                    "path": citation.get("path"),
                    "line": citation.get("line"),
                    "text": str(citation.get("excerpt")),
                }
            )
    deduped: list[dict[str, Any]] = []
    seen = set()
    for source in sources:
        key = (source.get("source_kind"), source.get("path"), source.get("line"), source.get("text"))
        if key not in seen:
            seen.add(key)
            deduped.append(source)
    return deduped


def target_value_evidence(env: Env, entry: dict[str, Any], target_info: dict[str, Any], alias_data: dict[str, Any]) -> list[dict[str, Any]]:
    evidence: list[dict[str, Any]] = []
    canonical = target_info["primary_target_parameter"]
    raw = target_info["raw_primary_target_parameter"]
    candidate_key = str(entry.get("candidate_key") or "")
    for source in target_value_sources(env, entry):
        values = extract_target_values_from_text(source.get("text") or "", canonical, raw, alias_data)
        if not values and re.search(r"(?:37_20|37_over_20|54_5|54_over_5|4107)", candidate_key, re.IGNORECASE):
            values = extract_values_from_rhs(str(source.get("text") or ""))
        for value in values:
            item = dict(value)
            item["source_kind"] = source.get("source_kind")
            item["path"] = source.get("path")
            item["line"] = source.get("line")
            if item not in evidence:
                evidence.append(item)
    evidence.sort(key=lambda item: (str(item.get("path")), str(item.get("line")), item.get("normalized", ""), item.get("raw", "")))
    return evidence


def target_layer_record(env: Env, entry: dict[str, Any], alias_data: dict[str, Any]) -> dict[str, Any]:
    target_info = parse_primary_target(entry, alias_data)
    value_evidence = target_value_evidence(env, entry, target_info, alias_data)
    target_values = sorted({item["normalized"] for item in value_evidence if item.get("normalized")})
    target_value_keys = sorted({item["comparison_key"] for item in value_evidence if item.get("comparison_key")})
    low_reasons = list(target_info.get("low_reasons") or [])
    if not target_values:
        low_reasons.append("target_value_unresolved")
    confidence = "low" if low_reasons or target_info.get("target_resolution_confidence") == "low" else "high"
    basis_parts = [str(target_info.get("resolution_basis") or "candidate_key parse")]
    if value_evidence:
        value_sources = sorted(
            {
                f"{item.get('source_kind')}:{normalize_citation_path(item.get('path'))}:{item.get('line')}"
                for item in value_evidence
            }
        )
        basis_parts.append("values from " + ", ".join(value_sources[:4]))
    else:
        basis_parts.append("no confident target value extracted from citation/source line")
    if low_reasons:
        basis_parts.append("low: " + ", ".join(sorted(set(low_reasons))))
    anchor_stages = [str(stage) for stage in entry.get("anchor_stages") or []]
    return {
        "candidate_id": str(entry.get("id")),
        "candidate_key": str(entry.get("candidate_key") or ""),
        "anchor_stages": anchor_stages,
        "anchor_stage": anchor_stages[0] if len(anchor_stages) == 1 else "+".join(anchor_stages),
        "raw_primary_target_parameter": target_info["raw_primary_target_parameter"],
        "primary_target_parameter": target_info["primary_target_parameter"],
        "additional_target_parameters": target_info.get("additional_target_parameters") or [],
        "multi_target": bool(target_info.get("multi_target")),
        "target_values": target_values,
        "target_value_keys": target_value_keys,
        "target_value_evidence": value_evidence,
        "target_confidence": confidence,
        "resolution_basis": "; ".join(basis_parts),
    }


def build_target_layer(env: Env, out_path: Path) -> dict[str, Any]:
    manifest = load_manifest(env)
    alias_data = load_concept_aliases(env)
    scanned_entries = [
        entry
        for entry in (manifest.get("candidates") or {}).values()
        if entry.get("status") == "scanned" and not entry.get("dry_run")
    ]
    records = [target_layer_record(env, entry, alias_data) for entry in sorted(scanned_entries, key=lambda item: str(item.get("id") or ""))]
    confidence_counts = Counter(record.get("target_confidence") for record in records)
    payload = {
        "schema_version": 1,
        "generated_at": now_iso(),
        "source_manifest": env.rel(env.manifest_path),
        "concept_aliases": env.rel(concept_alias_path(env)),
        "read_only_manifest": True,
        "rule": {
            "primary_target_resolution": [
                "candidate_key key_rules from _concept_aliases.yaml",
                "candidate_key token match against manifest parameter_names when no key_rule fires",
                "parameter_names fallback only when the key is uninformative",
                "concept aliases canonicalize variants to one primary_target_parameter",
            ],
            "target_confidence": "low when key fallback/multiple plausible targets/no confident target value; high otherwise",
            "value_extraction": [
                "citation excerpt and cited source line are both scanned",
                "sp.Float/sp.Rational/sp.Integer and Mathematica Rational[]",
                "JSON colon values",
                "chained equality/relation values",
                "TeX/plain fractions, decimals, and integers",
                "closed-form Sqrt/sqrt expressions",
                "exact decimal/fraction normalization such as 0.05 == 1/20",
            ],
        },
        "summary": {
            "scanned_candidate_count": len(records),
            "target_confidence_counts": dict(sorted(confidence_counts.items())),
            "value_extracted_count": sum(1 for record in records if record.get("target_values")),
            "multi_target_count": sum(1 for record in records if record.get("multi_target")),
        },
        "records": records,
    }
    write_yaml(out_path, payload)
    return payload


def load_or_build_target_layer(env: Env, *, build_missing: bool = True) -> dict[str, Any]:
    path = target_layer_path(env)
    if not path.exists():
        if not build_missing:
            raise SystemExit(f"error: target layer missing for read-only path: {env.rel(path)}")
        return build_target_layer(env, path)
    data = load_yaml(path)
    if not isinstance(data, dict) or not isinstance(data.get("records"), list):
        if not build_missing:
            raise SystemExit(f"error: target layer malformed for read-only path: {env.rel(path)}")
        return build_target_layer(env, path)
    return data


def target_layer_by_candidate(env: Env, *, build_missing: bool = True) -> dict[str, dict[str, Any]]:
    data = load_or_build_target_layer(env, build_missing=build_missing)
    return {str(record.get("candidate_id")): record for record in data.get("records") or []}


def citation_line(citation: dict[str, Any]) -> int | None:
    value = citation.get("line")
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def normalize_citation_path(path_text: Any) -> str:
    path = str(path_text or "").strip().replace("\\", "/")
    match = re.search(r"(?:^|/)research/pde_ledger/(.+)$", path)
    if match:
        return match.group(1)
    return path


def citations_adjacent(left: dict[str, Any], right: dict[str, Any]) -> bool:
    left_path = normalize_citation_path(left.get("path"))
    right_path = normalize_citation_path(right.get("path"))
    if not left_path or left_path != right_path:
        return False
    left_line = citation_line(left)
    right_line = citation_line(right)
    if left_line is not None and right_line is not None:
        return abs(left_line - right_line) <= ALIAS_ADJACENT_LINE_WINDOW
    if left_line is None and right_line is None:
        return str(left.get("excerpt") or "") == str(right.get("excerpt") or "")
    return False


def entries_have_adjacent_citation(left: dict[str, Any], right: dict[str, Any]) -> bool:
    for left_cit in left.get("file_line_citations") or []:
        for right_cit in right.get("file_line_citations") or []:
            if citations_adjacent(left_cit, right_cit):
                return True
    return False


def citation_summary(citations: list[dict[str, Any]]) -> dict[str, Any]:
    paths = sorted({normalize_citation_path(c.get("path")) for c in citations if c.get("path")})
    lines = sorted({citation_line(c) for c in citations if citation_line(c) is not None})
    summary: dict[str, Any] = {
        "paths": paths,
        "adjacency_rule": f"same path and line distance <= {ALIAS_ADJACENT_LINE_WINDOW}; missing lines require identical excerpt",
    }
    if len(paths) == 1:
        summary["path"] = paths[0]
    if lines:
        summary["line_min"] = min(lines)
        summary["line_max"] = max(lines)
    return summary


def canonical_candidate_entries(manifest: dict[str, Any], *, scanned_only: bool = True) -> list[dict[str, Any]]:
    entries = []
    for entry in (manifest.get("candidates") or {}).values():
        if entry.get("duplicate_of"):
            continue
        if scanned_only and entry.get("status") != "scanned":
            continue
        entries.append(entry)
    return sorted(entries, key=lambda item: str(item.get("id") or ""))


def stable_entry_id(prefix: str, parts: Iterable[Any]) -> str:
    raw = "|".join(str(part) for part in parts)
    digest = hashlib.sha1(raw.encode("utf-8")).hexdigest()[:12]
    return f"{prefix}_{digest}"


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


def graph_flatten_text(value: Any) -> Iterable[str]:
    if value is None:
        return
    if isinstance(value, dict):
        for key, child in value.items():
            yield str(key)
            yield from graph_flatten_text(child)
        return
    if isinstance(value, list):
        for child in value:
            yield from graph_flatten_text(child)
        return
    yield str(value)


def graph_normalize_terms(query: str) -> list[str]:
    return [term.lower() for term in re.split(r"\s+", query.strip()) if term]


class GraphSourceIndex:
    def __init__(self, graph: dict[str, Any]) -> None:
        self.graph = graph
        self.nodes = graph.get("nodes") or []

    @classmethod
    def load(cls, path: Path) -> "GraphSourceIndex":
        return cls(load_yaml(path) or {})

    def source_nodes(self, query: str) -> list[dict[str, Any]]:
        terms = graph_normalize_terms(query)
        results = []
        for node in self.nodes:
            text = " ".join(
                graph_flatten_text(
                    {
                        "file": node.get("file"),
                        "files": node.get("files"),
                        "sources": node.get("sources"),
                        "legacy_sources": node.get("legacy_sources"),
                        "legacy_file": node.get("legacy_file"),
                        "source_note_files": node.get("source_note_files"),
                        "canonical_target_files": node.get("canonical_target_files"),
                    }
                )
            ).lower()
            if text and all(term in text for term in terms):
                results.append(node)
        results.sort(key=lambda node: (node.get("layer", ""), node["id"]))
        return results


GRAPH_SOURCE_INDEX_CACHE: dict[Path, GraphSourceIndex] = {}


def graph_source_index(env: Env) -> GraphSourceIndex:
    graph_path = env.config_path_value("atlas_graph").resolve()
    cached = GRAPH_SOURCE_INDEX_CACHE.get(graph_path)
    if cached is None:
        try:
            cached = GraphSourceIndex.load(graph_path)
        except OSError as exc:
            raise SystemExit(f"error: cannot read atlas graph for source lookup: {graph_path}: {exc}") from exc
        GRAPH_SOURCE_INDEX_CACHE[graph_path] = cached
    return cached


def query_graph_source(
    env: Env,
    source_queries: list[str],
    graph_index: GraphSourceIndex | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    index = graph_index or graph_source_index(env)
    errors: list[dict[str, Any]] = []
    for query in source_queries:
        nodes = index.source_nodes(query)
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


def new_codex_session_record() -> dict[str, Any]:
    return {
        "session_id": None,
        "log_paths": [],
        "last_iter": None,
        "last_exit": None,
        "last_run": None,
    }


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
                param: new_codex_session_record()
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


def append_unique_items(target: list[Any], values: list[Any]) -> None:
    for value in values:
        if value not in target:
            target.append(value)


def merge_candidate_into_manifest(manifest: dict[str, Any], candidate: dict[str, Any]) -> tuple[dict[str, Any], bool, str | None]:
    candidates = manifest.setdefault("candidates", {})
    entry = candidates.get(candidate["id"])
    created = entry is None
    old_status = None if created else entry.get("status")
    if created:
        entry = initial_candidate_entry(candidate)
    else:
        append_unique_items(entry.setdefault("anchor_stages", []), candidate.get("anchor_stages") or [])
        append_unique_items(entry.setdefault("file_line_citations", []), candidate.get("citations") or [])
        append_unique_items(entry.setdefault("parameter_names", []), candidate.get("parameter_names") or [])
        append_unique_items(entry.setdefault("modality_attribution", []), candidate.get("modality_attribution") or [])
        append_unique_items(entry.setdefault("modality_fragments", []), candidate.get("modality_fragments") or [])
        entry["anchor_stages"] = sorted(entry.get("anchor_stages") or [])
        entry["parameter_names"] = sorted(entry.get("parameter_names") or [])
        entry["modality_attribution"] = sorted(entry.get("modality_attribution") or [])
        by_param = entry.setdefault("codex_session", {}).setdefault("by_parameter", {})
        for param in entry.get("parameter_names") or []:
            by_param.setdefault(param, new_codex_session_record())
    current = entry.get("status") or "pending"
    if current not in STATE_MACHINE or STATE_MACHINE.index(current) < STATE_MACHINE.index("scanned"):
        transition(entry, "scanned")
    candidates[candidate["id"]] = entry
    return entry, created, old_status


def stage_results_from_manifest(manifest: dict[str, Any]) -> dict[str, Any]:
    stage_to_entries: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for entry in (manifest.get("candidates") or {}).values():
        for stage in entry.get("anchor_stages") or []:
            stage_to_entries[str(stage)].append(entry)
    return {
        stage: {
            "candidate_count": len(entries),
            "structurally_vacuous": len(entries) == 0,
            "dry_run": all(bool(entry.get("dry_run")) for entry in entries),
        }
        for stage, entries in sorted(stage_to_entries.items())
    }


def render_fit_file(env: Env, manifest: dict[str, Any], stage_results: dict[str, Any] | None = None) -> None:
    candidates = []
    for entry in (manifest.get("candidates") or {}).values():
        item = {
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
        if entry.get("phase_b_status"):
            item["phase_b_status"] = entry.get("phase_b_status")
        if entry.get("duplicate_of"):
            item["duplicate_of"] = entry.get("duplicate_of")
        candidates.append(item)
    candidates.sort(key=lambda item: item["id"])
    if stage_results is None:
        stage_results = stage_results_from_manifest(manifest)
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
    alias_count = 0
    canonical_count = 0
    alias_group_count = 0
    batch_counts = defaultdict(int)
    for entry in (manifest.get("candidates") or {}).values():
        counts[entry.get("status", "unknown")] += 1
        if entry.get("dry_run"):
            dry_counts[entry.get("status", "unknown")] += 1
        if entry.get("duplicate_of"):
            alias_count += 1
        else:
            canonical_count += 1
        if entry.get("alias_ids"):
            alias_group_count += 1
        batch = entry.get("batch") or entry.get("batch_id")
        if batch:
            batch_counts[str(batch)] += 1
        verdict = entry.get("verdict") or {}
        if verdict.get("adversarial") or verdict.get("adjudication"):
            verdict_count += 1
    family_map_path = env.artifact_root / "provenance" / "_family_map.yaml"
    family_map = load_yaml(family_map_path) if family_map_path.exists() else {}
    family_summary = family_map.get("summary") or {}
    family_rows = []
    for family in family_map.get("families") or []:
        member_count = len(family.get("member_candidate_ids") or [])
        if member_count <= 1:
            continue
        family_rows.append(
            [
                str(family.get("family_id") or ""),
                str(member_count),
                ", ".join(str(p) for p in family.get("parameter_names") or []) or "-",
                ", ".join(str(v) for v in family.get("representative_values") or []) or "-",
                ", ".join(str(s) for s in family.get("stages_touched") or []) or "-",
            ]
        )
    lines = [
        "# Adversarial Audit Status",
        "",
        f"Generated: {now_iso()}",
        f"Project: {manifest.get('project_name', '?')}",
        "",
        "Authoritative consult record: `BATCHING_DECISIONS.md`.",
        "",
        "| Scope | Counts |",
        "|---|---|",
        f"| all candidates | {' '.join(f'{k}={v}' for k, v in sorted(counts.items())) or 'none'} |",
        f"| dry-run candidates | {' '.join(f'{k}={v}' for k, v in sorted(dry_counts.items())) or 'none'} |",
        f"| binding verdict fields populated | {verdict_count} |",
        f"| dedup canonicals / aliases | canonical={canonical_count} aliases={alias_count} alias_groups={alias_group_count} |",
        "",
        "## Phase B Dedup And Families",
        "",
        f"- Dedup state: {canonical_count} canonical entries, {alias_count} aliased entries.",
        f"- Family map: `{env.rel(family_map_path)}`" if family_map else "- Family map: not rendered yet.",
    ]
    if family_map:
        lines.extend(
            [
                f"- Families: {family_summary.get('family_count', len(family_rows))}; singletons: {family_summary.get('singleton_count', 'unknown')}; unmapped canonicals: {family_summary.get('unmapped_canonical_count', 'unknown')}.",
                "",
                "Non-singleton family groupings:",
                "",
                "| Family | Members | Parameters | Values | Stages |",
                "|---|---:|---|---|---|",
            ]
        )
        if family_rows:
            for family_id, member_count, params, values, stages in family_rows:
                lines.append(f"| {family_id} | {member_count} | {params} | {values} | {stages} |")
        else:
            lines.append("| none | 0 | - | - | - |")
    lines.extend(["", "## Phase C Batch Assignment", ""])
    if batch_counts:
        lines.extend(["| Batch | Candidate Count |", "|---|---:|"])
        for batch, count in sorted(batch_counts.items()):
            lines.append(f"| {batch} | {count} |")
    else:
        lines.append("batch assignment pending (Step 5)")
    write_text(env.batches_path, "\n".join(lines) + "\n")


def ensure_benchmarks_file(env: Env) -> dict[str, Any]:
    data = load_yaml(env.benchmarks_path)
    if not data:
        data = {"schema_version": 1, "entries": []}
    data.setdefault("entries", [])
    return data


def update_benchmarks_for_candidate(env: Env, candidate: dict[str, Any]) -> None:
    """Maintain only the dry-run placeholder owned by phase-b-build.

    Real sourced benchmark entries are created by benchmark-ingest and must be preserved
    across later provenance rebuilds for the same candidate.
    """
    if not candidate.get("dry_run"):
        return
    data = ensure_benchmarks_file(env)
    cid = candidate["id"]
    placeholder_id = f"{cid}__dry_run_benchmark_placeholder"
    entries = [e for e in data.get("entries", []) if e.get("id") != placeholder_id]
    citations = [
        c
        for c in candidate.get("file_line_citations", [])
        if "104" in c.get("path", "") or "105" in c.get("path", "")
    ][:8]
    params = set(candidate.get("parameter_names") or [])
    if candidate.get("dry_run") and "chi_Q" in params and citations:
        entries.append(
            {
                "id": placeholder_id,
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
    return yaml.dump(data, Dumper=_DUMPER, default_flow_style=False, sort_keys=False, width=120, allow_unicode=True).rstrip()


def stage_source_context(env: Env, stages: list[str]) -> list[dict[str, Any]]:
    sources: list[dict[str, Any]] = []
    for stage in stages:
        for role, path in source_files_for_stage(env, stage):
            sources.append({"stage": f"{int(stage):03d}", "role": role, "path": env.rel(path)})
    return sources


def phase_a_stage_reports(env: Env, stages: list[str]) -> list[str]:
    return [env.rel(pass2_report_path(env, stage)) for stage in stages if pass2_report_path(env, stage).exists()]


def critic_candidate_payload(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    candidates = []
    for entry in (manifest.get("candidates") or {}).values():
        candidates.append(
            {
                "id": entry.get("id"),
                "candidate_key": entry.get("candidate_key"),
                "dry_run": bool(entry.get("dry_run")),
                "dry_run_id": entry.get("dry_run_id"),
                "anchor_stages": entry.get("anchor_stages", []),
                "parameter_names": entry.get("parameter_names", []),
                "status": entry.get("status"),
                "file_line_citations": entry.get("file_line_citations", []),
                "modality_attribution": entry.get("modality_attribution", []),
            }
        )
    candidates.sort(key=lambda item: str(item.get("id") or ""))
    return candidates


def fragment_paths_for_critic(env: Env) -> list[str]:
    root = env.artifact_root / "phase_a_fragments"
    return [env.rel(path) for path in sorted(root.glob("*.yaml"))]


def stages_for_critic(fit_payload: dict[str, Any]) -> list[str]:
    stages = set()
    for candidate in fit_payload.get("candidates") or []:
        stages.update(str(stage) for stage in candidate.get("anchor_stages") or [])
    stages.update(str(stage) for stage in (fit_payload.get("stage_results") or {}).keys())
    return sorted(stages)


def render_critic_prompt(
    env: Env,
    *,
    prefix: str,
    manifest: dict[str, Any] | None = None,
    fit_payload: dict[str, Any] | None = None,
) -> str:
    ensure_tree(env)
    if manifest is None:
        manifest = load_manifest(env)
    if fit_payload is None:
        fit_payload = {
            "candidates": critic_candidate_payload(manifest),
            "stage_results": stage_results_from_manifest(manifest),
        }
    stages = stages_for_critic(fit_payload)
    critic_prompt = render_template(
        env,
        "phase_a_completeness_critic.md",
        {
            "STAGES": ", ".join(stages) if stages else "(none)",
            "FIT_INSERTION_POINTS_YAML": yaml_block(fit_payload),
            "MODALITY_FRAGMENT_PATHS": "\n".join(f"- {path}" for path in fragment_paths_for_critic(env)) or "- (none)",
        },
    )
    critic_path = env.artifact_root / "tmp_prompts" / f"{prefix}_phase_a_completeness_critic.md"
    write_text(critic_path, critic_prompt)
    return env.rel(critic_path)


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


def run_phase_a(env: Env, stages: list[str], dry_run: bool, dry_run_id: str | None, prefix: str | None = None) -> dict[str, Any]:
    ensure_tree(env)
    stages = [f"{int(s):03d}" for s in stages]
    prefix = prefix or "phase_a"
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
        merge_candidate_into_manifest(manifest, candidate)

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

    critic_prompt = render_critic_prompt(
        env,
        prefix=prefix,
        manifest=manifest,
        fit_payload={"candidates": candidates, "stage_results": stage_results},
    )

    return {
        "stages": stages,
        "candidate_ids": [c["id"] for c in candidates],
        "stage_results": stage_results,
        "graph_gaps": graph_gaps,
        "modality_prompts": modality_prompts,
        "critic_prompt": critic_prompt,
    }


def resolve_input_path(env: Env, path_text: str) -> Path:
    path = Path(path_text).expanduser()
    if path.is_absolute():
        return path
    cwd_candidate = (Path.cwd() / path).resolve()
    if cwd_candidate.exists():
        return cwd_candidate
    project_candidate = (env.project_root / path).resolve()
    if project_candidate.exists():
        return project_candidate
    repo_candidate = (env.repo_root / path).resolve()
    if repo_candidate.exists():
        return repo_candidate
    return cwd_candidate


def candidate_error_label(index: int, candidate: Any) -> str:
    key = None
    if isinstance(candidate, dict):
        key = candidate.get("candidate_key")
    suffix = f" ({key})" if key else ""
    return f"candidate[{index}]{suffix}"


def require_candidate_field(
    errors: list[str],
    file_label: str,
    candidate: dict[str, Any],
    index: int,
    field: str,
) -> Any:
    value = candidate.get(field)
    if value is None or value == "":
        errors.append(f"{file_label}: {candidate_error_label(index, candidate)} missing required field: {field}")
    return value


def normalize_agent_fragment_file(env: Env, source: Path) -> tuple[dict[str, Any] | None, list[str]]:
    file_label = env.rel(source)
    errors: list[str] = []
    if not source.exists():
        return None, [f"{file_label}: fragment file does not exist"]
    try:
        data = load_yaml(source)
    except yaml.YAMLError as exc:
        return None, [f"{file_label}: invalid YAML: {exc}"]
    if not isinstance(data, dict):
        return None, [f"{file_label}: top-level YAML must be a mapping"]
    modality = data.get("modality")
    if modality not in INGEST_MODALITIES:
        errors.append(
            f"{file_label}: modality must be one of {', '.join(INGEST_MODALITIES)}; got {modality!r}"
        )
    candidates_raw = data.get("candidates")
    if candidates_raw is None:
        errors.append(f"{file_label}: missing required top-level field: candidates")
        candidates_raw = []
    if not isinstance(candidates_raw, list):
        errors.append(f"{file_label}: candidates must be a list")
        candidates_raw = []
    dry_run = bool(data.get("dry_run", False))
    dry_run_id = data.get("dry_run_id")
    if dry_run and not dry_run_id:
        errors.append(f"{file_label}: dry_run: true fragments must include dry_run_id for purgeability")
    dry_run_id_text = str(dry_run_id) if dry_run_id else None

    normalized_candidates: list[dict[str, Any]] = []
    for index, candidate in enumerate(candidates_raw, 1):
        error_count_before_candidate = len(errors)
        if not isinstance(candidate, dict):
            errors.append(f"{file_label}: candidate[{index}] must be a mapping")
            continue
        candidate_key = require_candidate_field(errors, file_label, candidate, index, "candidate_key")
        anchor_stage = require_candidate_field(errors, file_label, candidate, index, "anchor_stage")
        parameter_names = candidate.get("parameter_names")
        if not isinstance(parameter_names, list):
            errors.append(f"{file_label}: {candidate_error_label(index, candidate)} parameter_names must be a non-empty list")
            parameter_names = []
        elif not parameter_names:
            errors.append(f"{file_label}: {candidate_error_label(index, candidate)} parameter_names must be non-empty")
        elif any(str(param).strip() == "" for param in parameter_names):
            errors.append(f"{file_label}: {candidate_error_label(index, candidate)} parameter_names contains an empty value")
        citation = candidate.get("citation")
        if not isinstance(citation, dict):
            errors.append(f"{file_label}: {candidate_error_label(index, candidate)} citation must be a mapping")
            citation = {}
        citation_path = citation.get("path")
        if not citation_path:
            errors.append(f"{file_label}: {candidate_error_label(index, candidate)} citation missing required field: path")
        if "line" not in citation:
            errors.append(f"{file_label}: {candidate_error_label(index, candidate)} citation missing required field: line")
        citation_excerpt = citation.get("excerpt")
        if not citation_excerpt:
            errors.append(f"{file_label}: {candidate_error_label(index, candidate)} citation missing required field: excerpt")
        reason = require_candidate_field(errors, file_label, candidate, index, "reason")
        try:
            stage = f"{int(str(anchor_stage)):03d}"
        except (TypeError, ValueError):
            errors.append(f"{file_label}: {candidate_error_label(index, candidate)} anchor_stage is not a stage number: {anchor_stage!r}")
            stage = str(anchor_stage or "")
        if len(errors) != error_count_before_candidate:
            continue
        normalized_citation = dict(citation)
        normalized_citation["path"] = str(citation_path)
        normalized_citation["stage"] = str(normalized_citation.get("stage") or stage)
        normalized = {
            "candidate_key": str(candidate_key),
            "modality": str(modality),
            "anchor_stage": stage,
            "parameter_names": [str(param) for param in parameter_names],
            "reason": str(reason),
            "citation": normalized_citation,
        }
        for optional_key in ("graph_node",):
            if optional_key in candidate:
                normalized[optional_key] = candidate[optional_key]
        normalized_candidates.append(normalized)

    if errors:
        return None, errors

    stages = sorted({candidate["anchor_stage"] for candidate in normalized_candidates})
    payload: dict[str, Any] = {
        "schema_version": int(data.get("schema_version") or 1),
        "dry_run": dry_run,
        "dry_run_id": dry_run_id_text,
        "modality": str(modality),
        "ingested_from": file_label,
        "ingested_at": now_iso(),
        "stages": stages,
        "candidates": normalized_candidates,
    }
    if modality in MODALITIES:
        payload["blind_to_modalities"] = [m for m in MODALITIES if m != modality]
    if "graph_gaps" in data:
        payload["graph_gaps"] = data.get("graph_gaps") or []
    if "pure_identities" in data:
        payload["pure_identities"] = data.get("pure_identities") or []
    return payload, []


def persisted_agent_fragment_path(env: Env, source: Path, payload: dict[str, Any]) -> Path:
    fragment_root = env.artifact_root / "phase_a_fragments"
    try:
        if source.resolve().parent == fragment_root.resolve():
            return source.resolve()
    except OSError:
        pass
    dry_run_id = payload.get("dry_run_id")
    prefix = str(dry_run_id) if payload.get("dry_run") and dry_run_id else "agent"
    return fragment_root / f"{prefix}_{slug(source.stem)}_{payload['modality']}.yaml"


def ingest_phase_a_fragments(env: Env, fragment_paths: list[str]) -> dict[str, Any]:
    ensure_tree(env)
    normalized_payloads: list[tuple[Path, Path, dict[str, Any]]] = []
    errors: list[str] = []
    for path_text in fragment_paths:
        source = resolve_input_path(env, path_text)
        payload, file_errors = normalize_agent_fragment_file(env, source)
        errors.extend(file_errors)
        if payload is not None:
            normalized_payloads.append((source, persisted_agent_fragment_path(env, source, payload), payload))
    if errors:
        raise SystemExit("error: malformed Phase A fragment input:\n- " + "\n- ".join(errors))

    for _source, dest, payload in normalized_payloads:
        write_yaml(dest, payload)

    manifest = load_manifest(env)
    grouped_fragments: dict[tuple[bool, str | None], list[dict[str, Any]]] = defaultdict(list)
    dry_run_fragments: dict[str, dict[str, Any]] = {}
    persisted_paths = []
    for _source, dest, payload in normalized_payloads:
        persisted_paths.append(env.rel(dest))
        dry_run = bool(payload.get("dry_run"))
        dry_run_id = payload.get("dry_run_id")
        grouped_fragments[(dry_run, dry_run_id)].extend(payload.get("candidates") or [])
        if dry_run and dry_run_id:
            record = dry_run_fragments.setdefault(
                str(dry_run_id),
                {
                    "dry_run": True,
                    "stages": set(),
                    "ingested_fragments": [],
                },
            )
            record["stages"].update(payload.get("stages") or [])
            record["ingested_fragments"].append(env.rel(dest))

    created_candidate_ids: list[str] = []
    updated_candidate_ids: list[str] = []
    status_preserved: dict[str, str] = {}
    merged_candidate_ids: list[str] = []
    dry_run_candidate_ids: dict[str, list[str]] = defaultdict(list)
    for (dry_run, dry_run_id), fragments in grouped_fragments.items():
        for candidate in merge_fragments(fragments, dry_run=dry_run, dry_run_id=dry_run_id):
            entry, created, old_status = merge_candidate_into_manifest(manifest, candidate)
            merged_candidate_ids.append(entry["id"])
            if dry_run and dry_run_id:
                dry_run_candidate_ids[str(dry_run_id)].append(entry["id"])
            if created:
                created_candidate_ids.append(entry["id"])
            else:
                updated_candidate_ids.append(entry["id"])
                if old_status and old_status == entry.get("status"):
                    status_preserved[entry["id"]] = old_status

    for dry_run_id, record in dry_run_fragments.items():
        manifest_record = manifest.setdefault("dry_runs", {}).setdefault(
            dry_run_id,
            {
                "dry_run": True,
                "started_at": now_iso(),
            },
        )
        existing_stages = set(manifest_record.get("stages") or [])
        existing_stages.update(record["stages"])
        manifest_record["stages"] = sorted(existing_stages)
        manifest_record["last_phase"] = "phase_a_ingest"
        append_unique_items(manifest_record.setdefault("ingested_fragments", []), record["ingested_fragments"])
        append_unique_items(manifest_record.setdefault("candidate_ids", []), dry_run_candidate_ids.get(dry_run_id, []))
        manifest_record["stage_results"] = {
            stage: {
                "candidate_count": sum(
                    1
                    for cid in manifest_record.get("candidate_ids", [])
                    if stage in ((manifest.get("candidates") or {}).get(cid, {}).get("anchor_stages") or [])
                ),
                "structurally_vacuous": False,
                "dry_run": True,
            }
            for stage in manifest_record["stages"]
        }

    render_fit_file(env, manifest, stage_results=stage_results_from_manifest(manifest))
    render_batches(env, manifest)
    save_manifest(env, manifest)
    return {
        "ingested_fragment_paths": persisted_paths,
        "candidate_ids": sorted(set(merged_candidate_ids)),
        "created_candidate_ids": sorted(created_candidate_ids),
        "updated_candidate_ids": sorted(set(updated_candidate_ids)),
        "status_preserved": dict(sorted(status_preserved.items())),
    }


def alias_member_evidence(entry: dict[str, Any], target_record: dict[str, Any] | None = None) -> dict[str, Any]:
    target_record = target_record or {}
    return {
        "id": entry.get("id"),
        "candidate_key": entry.get("candidate_key"),
        "primary_target_parameter": target_record.get("primary_target_parameter"),
        "raw_primary_target_parameter": target_record.get("raw_primary_target_parameter"),
        "target_values": target_record.get("target_values") or [],
        "target_value_keys": target_record.get("target_value_keys") or [],
        "target_confidence": target_record.get("target_confidence"),
        "citation": citation_summary(entry.get("file_line_citations") or []),
        "citations": entry.get("file_line_citations") or [],
        "modalities": entry.get("modality_attribution") or [],
    }


def dedup_record(entry: dict[str, Any], target_record: dict[str, Any]) -> dict[str, Any]:
    anchor_stages = [str(stage) for stage in entry.get("anchor_stages") or []]
    return {
        "id": str(entry.get("id")),
        "entry": entry,
        "anchor_stage": anchor_stages[0] if len(anchor_stages) == 1 else "+".join(anchor_stages),
        "primary_target_parameter": target_record.get("primary_target_parameter"),
        "target_values": target_record.get("target_values") or [],
        "target_value_keys": target_record.get("target_value_keys") or [],
        "target_confidence": target_record.get("target_confidence"),
        "target_record": target_record,
        "citations": entry.get("file_line_citations") or [],
        "modality_count": len(entry.get("modality_attribution") or []),
    }


def connected_components(
    records: list[dict[str, Any]],
    adjacent: Any,
) -> list[list[dict[str, Any]]]:
    remaining = {record["id"]: record for record in records}
    components: list[list[dict[str, Any]]] = []
    while remaining:
        first_id = sorted(remaining)[0]
        stack = [remaining.pop(first_id)]
        component = []
        while stack:
            current = stack.pop()
            component.append(current)
            adjacent_ids = [
                other_id
                for other_id, other in remaining.items()
                if adjacent(current, other)
            ]
            for other_id in sorted(adjacent_ids):
                stack.append(remaining.pop(other_id))
        components.append(sorted(component, key=lambda record: record["id"]))
    return components


def component_parameter_family(component: list[dict[str, Any]]) -> list[str]:
    param_sets = [set(record.get("parameter_family") or []) for record in component]
    common = set.intersection(*param_sets) if param_sets else set()
    if common:
        return sorted(common)
    return sorted(set().union(*param_sets)) if param_sets else []


def component_value_conflict(component: list[dict[str, Any]]) -> bool:
    value_sets = [set(record.get("literal_values") or []) for record in component if record.get("literal_values")]
    if len(value_sets) < 2:
        return False
    return not set.intersection(*value_sets) and len(set().union(*value_sets)) > 1


def split_component_by_literal_value(component: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    if not component_value_conflict(component):
        return [component]
    buckets: dict[str, list[dict[str, Any]]] = defaultdict(list)
    unknown: list[dict[str, Any]] = []
    value_counts = Counter(value for record in component for value in (record.get("literal_values") or []))
    for record in component:
        values = record.get("literal_values") or []
        if not values:
            unknown.append(record)
            continue
        key = sorted(values, key=lambda value: (-value_counts[value], value))[0]
        buckets[key].append(record)
    if buckets and unknown:
        largest_key = sorted(buckets, key=lambda key: (-len(buckets[key]), key))[0]
        buckets[largest_key].extend(unknown)
    elif unknown:
        buckets["__no_literal_value__"].extend(unknown)
    return [sorted(records, key=lambda record: record["id"]) for _key, records in sorted(buckets.items())]


def choose_canonical(component: list[dict[str, Any]]) -> dict[str, Any]:
    return sorted(component, key=lambda record: (-record["modality_count"], record["id"]))[0]


def pure_path_prefix_variant(left_path: str, right_path: str) -> bool:
    if not left_path or not right_path:
        return False
    left = left_path.strip("/")
    right = right_path.strip("/")
    return left == right or left.endswith("/" + right) or right.endswith("/" + left)


def citation_pair_alias_eligible(left: dict[str, Any], right: dict[str, Any]) -> bool:
    left_raw = str(left.get("path") or "")
    right_raw = str(right.get("path") or "")
    left_path = normalize_citation_path(left_raw)
    right_path = normalize_citation_path(right_raw)
    left_line = citation_line(left)
    right_line = citation_line(right)
    if left_raw == right_raw and left_line == right_line:
        return True
    if left_path and left_path == right_path and left_line is not None and right_line is not None:
        return abs(left_line - right_line) <= ALIAS_ADJACENT_LINE_WINDOW
    if left_line is not None and right_line is not None and left_line == right_line:
        return pure_path_prefix_variant(left_raw, right_raw) or pure_path_prefix_variant(left_path, right_path)
    return False


def records_have_alias_eligible_citation(left: dict[str, Any], right: dict[str, Any]) -> bool:
    for left_cit in left.get("citations") or []:
        for right_cit in right.get("citations") or []:
            if citation_pair_alias_eligible(left_cit, right_cit):
                return True
    return False


def strict_target_alias_adjacent(left: dict[str, Any], right: dict[str, Any]) -> bool:
    return (
        left.get("anchor_stage") == right.get("anchor_stage")
        and left.get("primary_target_parameter") == right.get("primary_target_parameter")
        and left.get("target_confidence") == "high"
        and right.get("target_confidence") == "high"
        and len(set(left.get("target_value_keys") or [])) == 1
        and set(left.get("target_value_keys") or []) == set(right.get("target_value_keys") or [])
        and records_have_alias_eligible_citation(left, right)
    )


def ambiguous_group_record(
    stage: str,
    primary_target_parameter: str,
    records: list[dict[str, Any]],
    reason: str,
    components: list[list[dict[str, Any]]] | None = None,
) -> dict[str, Any]:
    return {
        "anchor_stage": stage,
        "primary_target_parameter": primary_target_parameter,
        "reason": reason,
        "component_count": len(components or []),
        "members": [
            alias_member_evidence(record["entry"], record.get("target_record"))
            for record in sorted(records, key=lambda item: item["id"])
        ],
    }


def build_dedup_proposal(env: Env, out_path: Path) -> dict[str, Any]:
    manifest = load_manifest(env)
    targets_by_id = target_layer_by_candidate(env)
    scanned_entries = [
        entry
        for entry in (manifest.get("candidates") or {}).values()
        if entry.get("status") == "scanned" and not entry.get("duplicate_of") and not entry.get("dry_run")
    ]
    records: list[dict[str, Any]] = []
    for entry in scanned_entries:
        target_record = targets_by_id.get(str(entry.get("id")))
        if not target_record:
            target_record = target_layer_record(env, entry, load_concept_aliases(env))
        records.append(dedup_record(entry, target_record))

    alias_groups: list[dict[str, Any]] = []
    alias_member_ids: set[str] = set()
    value_buckets: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        value_keys = sorted(set(record.get("target_value_keys") or []))
        if record.get("target_confidence") != "high" or len(value_keys) != 1:
            continue
        value_buckets[(record["anchor_stage"], str(record.get("primary_target_parameter")), value_keys[0])].append(record)

    for (stage, target, value_key), bucket in sorted(value_buckets.items()):
        if len(bucket) < 2:
            continue
        for component in connected_components(bucket, strict_target_alias_adjacent):
            if len(component) < 2:
                continue
            canonical = choose_canonical(component)
            aliases = [record["id"] for record in component if record["id"] != canonical["id"]]
            if not aliases:
                continue
            all_citations = [citation for record in component for citation in record.get("citations") or []]
            alias_member_ids.update(record["id"] for record in component)
            alias_groups.append(
                {
                    "canonical_id": canonical["id"],
                    "aliases": sorted(aliases),
                    "shared": {
                        "anchor_stage": stage,
                        "primary_target_parameter": target,
                        "target_value": value_key,
                        "citation": citation_summary(all_citations),
                    },
                    "canonical_selection_rule": "most modality attributions; ties broken by candidate id sort",
                    "members": [alias_member_evidence(record["entry"], record.get("target_record")) for record in component],
                }
            )

    ambiguous: list[dict[str, Any]] = []
    ambiguous_member_ids: set[str] = set()
    grouped_by_stage_target: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        if record["id"] in alias_member_ids:
            continue
        grouped_by_stage_target[(record["anchor_stage"], str(record.get("primary_target_parameter")))].append(record)

    for (stage, target), bucket in sorted(grouped_by_stage_target.items()):
        if len(bucket) < 2:
            continue
        value_sets = {tuple(sorted(set(record.get("target_value_keys") or []))) for record in bucket}
        confident_single_value = [
            record
            for record in bucket
            if record.get("target_confidence") == "high" and len(set(record.get("target_value_keys") or [])) == 1
        ]
        all_same_site = False
        if len(confident_single_value) >= 2 and len(value_sets) == 1:
            all_same_site = any(
                records_have_alias_eligible_citation(left, right)
                for i, left in enumerate(confident_single_value)
                for right in confident_single_value[i + 1 :]
            )
        reason_bits = []
        if len(value_sets) > 1:
            reason_bits.append("differing target values")
        if any(record.get("target_confidence") != "high" for record in bucket):
            reason_bits.append("one or more low-confidence target resolutions")
        if not all_same_site:
            reason_bits.append("citations are outside strict alias adjacency")
        ambiguous_member_ids.update(record["id"] for record in bucket)
        ambiguous.append(
            ambiguous_group_record(
                stage,
                target,
                bucket,
                "; ".join(reason_bits) or "same stage and target but not strict same-value same-site aliases",
                connected_components(bucket, lambda left, right: records_have_alias_eligible_citation(left, right)),
            )
        )

    alias_groups.sort(key=lambda group: str(group["canonical_id"]))
    ambiguous.sort(key=lambda group: (str(group.get("anchor_stage")), str(group.get("primary_target_parameter"))))
    overlap = sorted(alias_member_ids & ambiguous_member_ids)
    if overlap:
        raise SystemExit("error: dedup proposal invariant failed; alias/ambiguous overlap: " + ", ".join(overlap[:20]))
    total_aliases = sum(len(group.get("aliases") or []) for group in alias_groups)
    standalone_ids = sorted({record["id"] for record in records} - alias_member_ids - ambiguous_member_ids)
    payload = {
        "schema_version": 1,
        "generated_at": now_iso(),
        "source_manifest": env.rel(env.manifest_path),
        "target_layer": env.rel(target_layer_path(env)),
        "read_only_manifest": True,
        "rule": {
            "alias": (
                "same anchor_stage, same concept-normalized primary_target_parameter, exactly one non-empty normalized target_value "
                "with high confidence on both sides, and strict citation identity/adjacency"
            ),
            "adjacent_line_window": ALIAS_ADJACENT_LINE_WINDOW,
            "removed_fallbacks": (
                "the value-compatibility / no-conflicting-value fallback and parameter-family union merging are disabled"
            ),
            "canonical_selection": "highest modality-attribution count; ties by lexicographic candidate id",
            "deletion_policy": "no deletions; aliases keep manifest entries and point duplicate_of to the canonical when applied",
            "partition_invariant": "candidate ids may appear in exactly one of alias_groups, ambiguous, or standalone_canonical_ids",
        },
        "summary": {
            "scanned_candidate_count": len(scanned_entries),
            "alias_group_count": len(alias_groups),
            "total_aliases": total_aliases,
            "resulting_canonical_count": len(scanned_entries) - total_aliases,
            "ambiguous_group_count": len(ambiguous),
            "alias_ambiguous_id_overlap": len(overlap),
            "standalone_canonical_count": len(standalone_ids),
        },
        "alias_groups": alias_groups,
        "ambiguous": ambiguous,
        "standalone_canonical_ids": standalone_ids,
    }
    write_yaml(out_path, payload)
    return payload["summary"] | {"out": env.rel(out_path)}


def resolve_output_path(env: Env, path_text: str | None, default_rel: str) -> Path:
    if not path_text:
        return env.artifact_root / default_rel
    path = Path(path_text).expanduser()
    if path.is_absolute():
        return path
    return (Path.cwd() / path).resolve()


def merge_alias_into_canonical(canonical: dict[str, Any], alias: dict[str, Any]) -> None:
    append_unique_items(canonical.setdefault("anchor_stages", []), alias.get("anchor_stages") or [])
    append_unique_items(canonical.setdefault("file_line_citations", []), alias.get("file_line_citations") or [])
    append_unique_items(canonical.setdefault("parameter_names", []), alias.get("parameter_names") or [])
    append_unique_items(canonical.setdefault("modality_attribution", []), alias.get("modality_attribution") or [])
    append_unique_items(canonical.setdefault("modality_fragments", []), alias.get("modality_fragments") or [])
    append_unique_items(canonical.setdefault("alias_ids", []), [alias["id"]])
    canonical["anchor_stages"] = sorted(canonical.get("anchor_stages") or [])
    canonical["parameter_names"] = sorted(canonical.get("parameter_names") or [])
    canonical["modality_attribution"] = sorted(canonical.get("modality_attribution") or [])
    by_param = canonical.setdefault("codex_session", {}).setdefault("by_parameter", {})
    for param in canonical.get("parameter_names") or []:
        by_param.setdefault(param, new_codex_session_record())


def apply_alias_map(env: Env, map_path_text: str) -> dict[str, Any]:
    source = resolve_input_path(env, map_path_text)
    data = load_yaml(source)
    if not isinstance(data, dict):
        raise SystemExit(f"error: alias map must be a YAML mapping: {env.rel(source)}")
    groups = data.get("alias_groups")
    if groups is None:
        groups = data.get("groups")
    if not isinstance(groups, list):
        raise SystemExit(f"error: alias map missing list field: alias_groups")

    manifest = load_manifest(env)
    candidates = manifest.setdefault("candidates", {})
    errors: list[str] = []
    normalized_groups: list[tuple[str, list[str]]] = []
    seen_aliases: dict[str, str] = {}
    for index, group in enumerate(groups, 1):
        if not isinstance(group, dict):
            errors.append(f"group[{index}] must be a mapping")
            continue
        canonical_id = str(group.get("canonical_id") or "")
        aliases = group.get("aliases")
        if not canonical_id:
            errors.append(f"group[{index}] missing canonical_id")
            continue
        if canonical_id not in candidates:
            errors.append(f"group[{index}] unknown canonical_id: {canonical_id}")
        canonical = candidates.get(canonical_id) or {}
        if canonical.get("duplicate_of"):
            errors.append(f"group[{index}] canonical_id is already an alias: {canonical_id} -> {canonical.get('duplicate_of')}")
        if not isinstance(aliases, list) or not aliases:
            errors.append(f"group[{index}] aliases must be a non-empty list")
            continue
        normalized_aliases = []
        for alias in aliases:
            alias_id = str(alias)
            if alias_id == canonical_id:
                errors.append(f"group[{index}] alias repeats canonical_id: {alias_id}")
                continue
            if alias_id not in candidates:
                errors.append(f"group[{index}] unknown alias id: {alias_id}")
                continue
            existing_target = candidates[alias_id].get("duplicate_of")
            if existing_target and existing_target != canonical_id:
                errors.append(f"group[{index}] alias already points elsewhere: {alias_id} -> {existing_target}")
                continue
            if seen_aliases.get(alias_id) and seen_aliases[alias_id] != canonical_id:
                errors.append(f"group[{index}] alias appears under multiple canonicals: {alias_id}")
                continue
            seen_aliases[alias_id] = canonical_id
            normalized_aliases.append(alias_id)
        normalized_groups.append((canonical_id, sorted(set(normalized_aliases))))
    if errors:
        raise SystemExit("error: malformed alias map:\n- " + "\n- ".join(errors))

    applied = 0
    already_applied = 0
    for canonical_id, aliases in normalized_groups:
        canonical = candidates[canonical_id]
        for alias_id in aliases:
            alias = candidates[alias_id]
            if alias.get("duplicate_of") == canonical_id:
                already_applied += 1
                continue
            merge_alias_into_canonical(canonical, alias)
            alias["duplicate_of"] = canonical_id
            alias.setdefault("status_notes", []).append(
                {
                    "status": alias.get("status"),
                    "note": f"Marked duplicate_of {canonical_id} by alias map {env.rel(source)}",
                    "at": now_iso(),
                }
            )
            applied += 1
        candidates[canonical_id] = canonical
    render_fit_file(env, manifest)
    render_batches(env, manifest)
    save_manifest(env, manifest)
    return {
        "alias_map": env.rel(source),
        "groups": len(normalized_groups),
        "aliases_applied": applied,
        "aliases_already_applied": already_applied,
    }


def family_member_record(entry: dict[str, Any], target_record: dict[str, Any], family_targets: list[str] | None = None) -> dict[str, Any]:
    return {
        "candidate_id": entry.get("id"),
        "candidate_key": entry.get("candidate_key"),
        "anchor_stages": entry.get("anchor_stages") or [],
        "primary_target_parameter": target_record.get("primary_target_parameter"),
        "family_targets": family_targets or [target_record.get("primary_target_parameter")],
        "target_values": target_record.get("target_values") or [],
        "target_value_keys": target_record.get("target_value_keys") or [],
        "target_confidence": target_record.get("target_confidence"),
        "resolution_basis": target_record.get("resolution_basis"),
    }


def family_value_divergence_finding(members: list[dict[str, Any]]) -> dict[str, Any] | None:
    by_value: dict[str, list[dict[str, Any]]] = defaultdict(list)
    display_for_key: dict[str, str] = {}
    for member in members:
        keys = member.get("target_value_keys") or []
        values = member.get("target_values") or []
        for index, key in enumerate(keys):
            if not key:
                continue
            display_for_key.setdefault(str(key), str(values[index] if index < len(values) else key))
            by_value[str(key)].append(
                {
                    "candidate_id": member.get("candidate_id"),
                    "stages": member.get("anchor_stages") or [],
                    "primary_target_parameter": member.get("primary_target_parameter"),
                }
            )
    if len(by_value) <= 1:
        return None
    return {
        "type": "value_divergence",
        "severity": "needs_agent_review",
        "summary": "Same conceptual family carries divergent normalized target values.",
        "distinct_values": [
            {
                "value_key": key,
                "display_value": display_for_key.get(key, key),
                "occurrences": by_value[key],
            }
            for key in sorted(by_value)
        ],
    }


def family_record_for_members(
    family_id: str,
    family_kind: str,
    members: list[tuple[dict[str, Any], dict[str, Any]]],
    *,
    primary_target_parameter: str | None = None,
    concepts: list[str] | None = None,
    primary_family: bool = True,
    findings: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    member_records = [
        family_member_record(entry, target_record, concepts or ([primary_target_parameter] if primary_target_parameter else None))
        for entry, target_record in sorted(members, key=lambda item: str(item[0].get("id")))
    ]
    values = sorted({value for member in member_records for value in member.get("target_values") or []})
    all_findings = list(findings or [])
    divergence = family_value_divergence_finding(member_records)
    if divergence:
        all_findings.append(divergence)
    return {
        "family_id": family_id,
        "family_kind": family_kind,
        "primary_family": primary_family,
        "primary_target_parameter": primary_target_parameter,
        "concepts": concepts or ([primary_target_parameter] if primary_target_parameter else []),
        "parameter_names": concepts or ([primary_target_parameter] if primary_target_parameter else []),
        "representative_values": values,
        "member_candidate_ids": [member["candidate_id"] for member in member_records],
        "stages_touched": sorted({str(stage) for member in member_records for stage in member.get("anchor_stages") or []}),
        "members": member_records,
        "findings": all_findings,
    }


def build_seeded_chain_family(
    chain_id: str,
    chain: dict[str, Any],
    canonical_entries: list[dict[str, Any]],
    target_records: dict[str, dict[str, Any]],
    alias_data: dict[str, Any],
) -> dict[str, Any]:
    concepts = sorted({canonicalize_concept(concept, alias_data) for concept in chain.get("concepts") or []})
    members: list[tuple[dict[str, Any], dict[str, Any]]] = []
    concept_set = set(concepts)
    for entry in canonical_entries:
        record = target_records.get(str(entry.get("id")))
        if not record:
            continue
        record_concepts = {record.get("primary_target_parameter")}
        record_concepts.update(record.get("additional_target_parameters") or [])
        if record_concepts & concept_set:
            members.append((entry, record))
    member_ids = {str(entry.get("id")) for entry, _record in members}
    member_stages = {str(stage) for entry, _record in members for stage in entry.get("anchor_stages") or []}
    warnings = []
    missing_ids = [cid for cid in chain.get("expected_candidate_ids") or [] if cid not in member_ids]
    missing_stages = [stage for stage in chain.get("expected_stages") or [] if str(stage) not in member_stages]
    if missing_ids or missing_stages:
        warnings.append(
            {
                "type": "coverage_warning",
                "severity": "hard",
                "summary": "Seeded chain mechanical concept matching missed expected coverage.",
                "missing_expected_candidate_ids": missing_ids,
                "missing_expected_stages": missing_stages,
            }
        )
    return family_record_for_members(
        chain_id,
        "seeded_d3_chain",
        members,
        concepts=concepts,
        primary_family=False,
        findings=warnings,
    ) | {"description": chain.get("description")}


def build_family_map(env: Env, out_path: Path) -> dict[str, Any]:
    manifest = load_manifest(env)
    alias_data = load_concept_aliases(env)
    targets_by_id = target_layer_by_candidate(env)
    all_candidates = manifest.get("candidates") or {}
    alias_count = sum(1 for entry in all_candidates.values() if entry.get("duplicate_of"))
    canonical_entries = [
        entry
        for entry in all_candidates.values()
        if not entry.get("duplicate_of") and not entry.get("dry_run")
    ]
    grouped: dict[str, list[tuple[dict[str, Any], dict[str, Any]]]] = defaultdict(list)
    primary_candidate_family_map: dict[str, list[str]] = defaultdict(list)
    candidate_family_map: dict[str, list[str]] = defaultdict(list)
    target_records: dict[str, dict[str, Any]] = {}
    for entry in canonical_entries:
        record = targets_by_id.get(str(entry.get("id")))
        if not record:
            record = target_layer_record(env, entry, alias_data)
        target_records[str(entry.get("id"))] = record
        targets = [record.get("primary_target_parameter") or "unclassified_target_parameter"]
        if record.get("multi_target"):
            append_unique_items(targets, record.get("additional_target_parameters") or [])
        for target in targets:
            grouped[str(target)].append((entry, record))

    families: list[dict[str, Any]] = []
    for index, (target, members) in enumerate(sorted(grouped.items()), 1):
        family_id = f"fam_{index:04d}_{slug(target)}"
        family_kind = "conceptual_parameter" if len(members) > 1 else "explicit_singleton_conceptual_parameter"
        family = family_record_for_members(
            family_id,
            family_kind,
            members,
            primary_target_parameter=target,
            concepts=[target],
            primary_family=True,
        )
        families.append(family)
        for cid in family["member_candidate_ids"]:
            primary_candidate_family_map[cid].append(family_id)
            candidate_family_map[cid].append(family_id)

    seeded_chain_families = [
        build_seeded_chain_family(chain_id, chain, canonical_entries, target_records, alias_data)
        for chain_id, chain in sorted((alias_data.get("chains") or {}).items())
    ]
    for chain_family in seeded_chain_families:
        families.append(chain_family)
        for cid in chain_family.get("member_candidate_ids") or []:
            candidate_family_map[cid].append(chain_family["family_id"])

    canonical_ids = sorted(str(entry["id"]) for entry in canonical_entries)
    unmapped = [cid for cid in canonical_ids if not primary_candidate_family_map.get(cid)]
    chain_summaries = {
        family["family_id"]: {
            "member_count": len(family.get("member_candidate_ids") or []),
            "value_divergence": any(f.get("type") == "value_divergence" for f in family.get("findings") or []),
            "coverage_warning": any(f.get("type") == "coverage_warning" for f in family.get("findings") or []),
            "coverage_findings": [f for f in family.get("findings") or [] if f.get("type") == "coverage_warning"],
        }
        for family in seeded_chain_families
    }
    payload = {
        "schema_version": 1,
        "generated_at": now_iso(),
        "source_manifest": env.rel(env.manifest_path),
        "target_layer": env.rel(target_layer_path(env)),
        "concept_aliases": env.rel(concept_alias_path(env)),
        "read_only_manifest": True,
        "rule": {
            "family": (
                "primary families are keyed by concept-normalized primary_target_parameter alone; target values are per-member "
                "attributes and divergent values produce value_divergence findings"
            ),
            "alias_handling": "entries with duplicate_of are excluded; if no aliases are present this is a pre-apply map and should be re-run after apply-alias-map",
            "seeded_chains": "D3 headline chains are emitted as non-primary named overlay families populated by concept matching with hard coverage warnings for expected misses",
        },
        "summary": {
            "canonical_candidate_count": len(canonical_entries),
            "alias_count_excluded": alias_count,
            "pre_alias_map": alias_count == 0,
            "rerun_after_alias_apply": alias_count == 0,
            "family_count": len(families),
            "concept_family_count": len(grouped),
            "seeded_chain_count": len(seeded_chain_families),
            "singleton_count": sum(
                1
                for family in families
                if family.get("primary_family") and len(family.get("member_candidate_ids") or []) == 1
            ),
            "unmapped_canonical_count": len(unmapped),
            "unmapped_canonical_ids": unmapped,
            "seeded_chains": chain_summaries,
        },
        "primary_candidate_family_map": dict(sorted((cid, sorted(ids)) for cid, ids in primary_candidate_family_map.items())),
        "candidate_family_map": dict(sorted((cid, sorted(ids)) for cid, ids in candidate_family_map.items())),
        "seeded_chain_families": seeded_chain_families,
        "families": families,
    }
    write_yaml(out_path, payload)
    return payload["summary"] | {"out": env.rel(out_path)}


def source_evidence_for_candidate(env: Env, entry: dict[str, Any], target_parameters: list[str] | None = None) -> list[dict[str, Any]]:
    evidence: list[dict[str, Any]] = []
    params = [str(p) for p in (target_parameters if target_parameters is not None else entry.get("parameter_names") or [])]
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


def graph_context_for_candidate(
    env: Env,
    entry: dict[str, Any],
    graph_index: GraphSourceIndex | None = None,
) -> dict[str, Any]:
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
        nodes, errors = query_graph_source(env, source_queries, graph_index)
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


def phase_b_target_parameters(
    env: Env,
    entry: dict[str, Any],
    target_layer: dict[str, dict[str, Any]] | None = None,
) -> list[str]:
    if entry.get("dry_run"):
        return [str(p) for p in entry.get("parameter_names") or ["unclassified_target_parameter"]]
    targets_by_id = target_layer if target_layer is not None else target_layer_by_candidate(env)
    record = targets_by_id.get(str(entry.get("id")))
    if not record:
        return [str(p) for p in entry.get("parameter_names") or ["unclassified_target_parameter"]]
    targets = [str(record.get("primary_target_parameter") or "unclassified_target_parameter")]
    if record.get("multi_target"):
        append_unique_items(targets, [str(p) for p in record.get("additional_target_parameters") or []])
    return [target for target in targets if target]


def phase_b_bundle_paths_written(env: Env, entry: dict[str, Any], params: list[str]) -> bool:
    if not params:
        return False
    candidate_id = str(entry.get("id"))
    return all(phase_b_bundle_path(env, candidate_id, param).exists() for param in params)


def build_provenance_for_entry(
    env: Env,
    manifest: dict[str, Any],
    candidate_id: str,
    *,
    graph_index: GraphSourceIndex | None = None,
    target_layer: dict[str, dict[str, Any]] | None = None,
    save: bool = False,
    render: bool = False,
) -> dict[str, Any]:
    entry = (manifest.get("candidates") or {}).get(candidate_id)
    if not entry:
        raise SystemExit(f"error: unknown candidate: {candidate_id}")
    if entry.get("duplicate_of"):
        raise SystemExit(f"error: {candidate_id} is an alias of {entry.get('duplicate_of')}; build provenance on the canonical")

    params = phase_b_target_parameters(env, entry, target_layer)
    evidence = source_evidence_for_candidate(env, entry, params)
    graph_context = graph_context_for_candidate(env, entry, graph_index)
    provenance_paths: list[str] = []
    for param in params:
        prompt_payload = dict(entry)
        prompt_payload["phase_b_parameter_name"] = param
        prompt_payload["phase_b_target_parameters"] = params
        phase_b_prompt = render_template(
            env,
            "phase_b_provenance_builder.md",
            {
                "CANDIDATE_ID": candidate_id,
                "CANDIDATE_YAML": yaml_block(prompt_payload),
            },
        )
        prompt_path = phase_b_prompt_path(env, candidate_id, param)
        write_text(prompt_path, phase_b_prompt)
        payload = {
            "schema_version": 1,
            "candidate_id": candidate_id,
            "parameter_name": param,
            "dry_run": bool(entry.get("dry_run")),
            "dry_run_id": entry.get("dry_run_id"),
            "non_binding": bool(entry.get("dry_run")),
            "generated_content_kind": "mechanical_evidence_bundle",
            "synthesis_status": "pending",
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
        path = phase_b_bundle_path(env, candidate_id, param)
        write_yaml(path, payload)
        provenance_paths.append(env.rel(path))

    entry.setdefault("paths", {})["provenance"] = provenance_paths
    if entry.get("dry_run"):
        transition(entry, "provenance_built")
        entry["phase_b_status"] = "dry_run_mechanical_complete"
    else:
        entry["phase_b_status"] = "synthesis_pending"
    manifest["candidates"][candidate_id] = entry
    update_benchmarks_for_candidate(env, entry)
    if render:
        render_fit_file(env, manifest)
        render_batches(env, manifest)
    if save:
        save_manifest(env, manifest)
    return {
        "candidate_id": candidate_id,
        "provenance_paths": provenance_paths,
        "status": entry["status"],
        "bundle_count": len(provenance_paths),
    }


def build_provenance(env: Env, candidate_id: str) -> dict[str, Any]:
    manifest = load_manifest(env)
    return build_provenance_for_entry(
        env,
        manifest,
        candidate_id,
        graph_index=graph_source_index(env),
        target_layer=target_layer_by_candidate(env),
        save=True,
        render=True,
    )


def read_ids_file(path_text: str) -> tuple[list[str], list[str]]:
    path = Path(path_text).expanduser()
    if not path.is_absolute():
        path = (Path.cwd() / path).resolve()
    ids: list[str] = []
    duplicate_ids: list[str] = []
    seen: set[str] = set()
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        candidate_id = raw_line.split("#", 1)[0].strip()
        if not candidate_id:
            continue
        if candidate_id in seen:
            duplicate_ids.append(candidate_id)
            continue
        seen.add(candidate_id)
        ids.append(candidate_id)
    return ids, duplicate_ids


def phase_b_not_eligible_reason(entry: dict[str, Any] | None) -> str | None:
    if entry is None:
        return "unknown_candidate"
    if entry.get("duplicate_of"):
        return f"duplicate_of:{entry.get('duplicate_of')}"
    if entry.get("status") != "scanned":
        return f"status:{entry.get('status', 'unknown')}"
    return None


def manifest_count_summary(manifest: dict[str, Any]) -> dict[str, dict[str, int]]:
    candidates = (manifest.get("candidates") or {}).values()
    return {
        "status_counts": dict(sorted(Counter(str(entry.get("status", "unknown")) for entry in candidates).items())),
        "phase_b_status_counts": dict(
            sorted(Counter(str(entry.get("phase_b_status", "unset")) for entry in candidates).items())
        ),
    }


def phase_b_build_all(env: Env, args: argparse.Namespace) -> dict[str, Any]:
    manifest = load_manifest(env)
    candidates = manifest.get("candidates") or {}
    if args.ids_file:
        candidate_ids, duplicate_ids = read_ids_file(args.ids_file)
        ids_source = str(Path(args.ids_file).expanduser())
    else:
        candidate_ids = sorted(str(candidate_id) for candidate_id in candidates)
        duplicate_ids = []
        ids_source = "auto"

    target_layer = target_layer_by_candidate(env, build_missing=not args.dry)
    graph_index = None if args.dry else graph_source_index(env)
    target_items: list[tuple[str, list[str]]] = []
    skipped_already_built: list[str] = []
    skipped_not_eligible: list[dict[str, str]] = []
    target_errors: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []

    for candidate_id in candidate_ids:
        entry = candidates.get(candidate_id)
        reason = phase_b_not_eligible_reason(entry)
        if reason:
            skipped_not_eligible.append({"candidate_id": candidate_id, "reason": reason})
            continue
        try:
            params = phase_b_target_parameters(env, entry, target_layer)
            already_built = phase_b_bundle_paths_written(env, entry, params)
        except (Exception, SystemExit) as exc:
            error = {"candidate_id": candidate_id, "error": str(exc)}
            target_errors.append(error)
            errors.append(error)
            continue
        if not args.force and already_built:
            skipped_already_built.append(candidate_id)
            continue
        target_items.append((candidate_id, params))

    limit = args.limit
    if limit is not None and limit < 0:
        raise SystemExit("error: --limit must be non-negative")
    selected_items = target_items[:limit] if limit is not None else target_items
    limited_out_count = max(0, len(target_items) - len(selected_items))
    target_distribution = Counter(len(params) for _candidate_id, params in target_items)
    selected_distribution = Counter(len(params) for _candidate_id, params in selected_items)

    built = 0
    total_bundles_written = 0
    built_distribution: Counter[int] = Counter()
    built_ids: list[str] = []
    if not args.dry:
        for candidate_id, _params in selected_items:
            try:
                result = build_provenance_for_entry(
                    env,
                    manifest,
                    candidate_id,
                    graph_index=graph_index,
                    target_layer=target_layer,
                    save=False,
                    render=False,
                )
            except SystemExit as exc:
                errors.append({"candidate_id": candidate_id, "error": str(exc)})
                continue
            except Exception as exc:
                errors.append({"candidate_id": candidate_id, "error": f"{type(exc).__name__}: {exc}"})
                continue
            bundle_count = int(result.get("bundle_count") or len(result.get("provenance_paths") or []))
            built += 1
            total_bundles_written += bundle_count
            built_distribution[bundle_count] += 1
            built_ids.append(candidate_id)

        if built:
            render_fit_file(env, manifest)
            render_batches(env, manifest)
            save_manifest(env, manifest)

    target_parameter_counts = {candidate_id: len(params) for candidate_id, params in selected_items}
    report: dict[str, Any] = {
        "dry": bool(args.dry),
        "force": bool(args.force),
        "ids_source": ids_source,
        "ids_considered": len(candidate_ids),
        "duplicate_ids_ignored": duplicate_ids,
        "limit": limit,
        "target_set_size": len(target_items) + len(target_errors),
        "buildable_target_count": len(target_items),
        "target_error_count": len(target_errors),
        "selected_target_count": len(selected_items),
        "limited_out_count": limited_out_count,
        "built": built,
        "would_build": len(selected_items) if args.dry else 0,
        "skipped_already_built": len(skipped_already_built),
        "skipped_not_eligible": len(skipped_not_eligible),
        "total_bundles_written": total_bundles_written,
        "target_bundles_per_candidate_distribution": dict(sorted(target_distribution.items())),
        "selected_bundles_per_candidate_distribution": dict(sorted(selected_distribution.items())),
        "built_bundles_per_candidate_distribution": dict(sorted(built_distribution.items())),
        "target_parameter_counts": target_parameter_counts,
        "errors": errors,
    }
    if skipped_already_built:
        report["skipped_already_built_sample"] = skipped_already_built[:20]
    if skipped_not_eligible:
        report["skipped_not_eligible_sample"] = skipped_not_eligible[:20]
    if built_ids:
        report["built_candidate_ids"] = built_ids
    if not args.dry:
        report["saved_manifest"] = bool(built)
        report["final_manifest_counts"] = manifest_count_summary(manifest)
    return report


def citation_path_under_notes(path_text: Any) -> bool:
    path = str(path_text or "").strip().replace("\\", "/")
    if path.startswith("notes/"):
        return True
    return "/notes/" in path


def validate_notes_citation(
    errors: list[str],
    file_label: str,
    where: str,
    citation: Any,
    *,
    citation_field: str,
) -> dict[str, Any]:
    if not isinstance(citation, dict):
        errors.append(f"{file_label}: {where} {citation_field} must be a mapping")
        return {}
    path = citation.get("path")
    line = citation.get("line")
    excerpt = citation.get("excerpt")
    if not path:
        errors.append(f"{file_label}: {where} {citation_field} missing path")
    elif not citation_path_under_notes(path):
        errors.append(f"{file_label}: {where} {citation_field}.path must be under notes/: {path}")
    if line is None:
        errors.append(f"{file_label}: {where} {citation_field} missing line")
    if not excerpt:
        errors.append(f"{file_label}: {where} {citation_field} missing excerpt")
    return dict(citation)


def validate_optional_citation(errors: list[str], file_label: str, where: str, citation: Any) -> dict[str, Any]:
    if not isinstance(citation, dict):
        errors.append(f"{file_label}: {where} citation must be a mapping")
        return {}
    path = citation.get("path")
    if not path:
        errors.append(f"{file_label}: {where} citation missing path")
    if "line" not in citation:
        errors.append(f"{file_label}: {where} citation missing line")
    if not citation.get("excerpt"):
        errors.append(f"{file_label}: {where} citation missing excerpt")
    return dict(citation)


def normalize_phase_b_synthesis_file(env: Env, source: Path) -> tuple[dict[str, Any] | None, list[str]]:
    file_label = env.rel(source)
    errors: list[str] = []
    if not source.exists():
        return None, [f"{file_label}: synthesis file does not exist"]
    try:
        data = load_yaml(source)
    except yaml.YAMLError as exc:
        return None, [f"{file_label}: invalid YAML: {exc}"]
    if not isinstance(data, dict):
        return None, [f"{file_label}: top-level YAML must be a mapping"]
    candidate_id = str(data.get("candidate_id") or "")
    parameter_name = str(data.get("parameter_name") or data.get("parameter") or "")
    if not candidate_id:
        errors.append(f"{file_label}: missing required field: candidate_id")
    if not parameter_name:
        errors.append(f"{file_label}: missing required field: parameter_name")

    origin_claims_raw = data.get("origin_claims")
    constraints_raw = data.get("constraints")
    downstream_raw = data.get("downstream_dependents", [])
    findings_raw = data.get("provenance_findings", [])
    if not isinstance(origin_claims_raw, list):
        errors.append(f"{file_label}: origin_claims must be a list")
        origin_claims_raw = []
    if not isinstance(constraints_raw, list):
        errors.append(f"{file_label}: constraints must be a list")
        constraints_raw = []
    if not isinstance(downstream_raw, list):
        errors.append(f"{file_label}: downstream_dependents must be a list")
        downstream_raw = []
    if not isinstance(findings_raw, list):
        errors.append(f"{file_label}: provenance_findings must be a list")
        findings_raw = []

    origin_claims: list[dict[str, Any]] = []
    for index, claim in enumerate(origin_claims_raw, 1):
        where = f"origin_claims[{index}]"
        if not isinstance(claim, dict):
            errors.append(f"{file_label}: {where} must be a mapping")
            continue
        for field in ("parameter", "introduced_at_stage", "introduced_at_line"):
            if claim.get(field) in (None, ""):
                errors.append(f"{file_label}: {where} missing {field}")
        normalized = dict(claim)
        normalized["citation"] = validate_notes_citation(
            errors,
            file_label,
            where,
            claim.get("citation"),
            citation_field="citation",
        )
        origin_claims.append(normalized)

    constraints: list[dict[str, Any]] = []
    for index, constraint in enumerate(constraints_raw, 1):
        where = f"constraints[{index}]"
        if not isinstance(constraint, dict):
            errors.append(f"{file_label}: {where} must be a mapping")
            continue
        if not constraint.get("parameter"):
            errors.append(f"{file_label}: {where} missing parameter")
        kind = constraint.get("constraint_kind")
        if kind not in PHASE_B_CONSTRAINT_KINDS:
            errors.append(
                f"{file_label}: {where} constraint_kind must be one of {', '.join(sorted(PHASE_B_CONSTRAINT_KINDS))}; got {kind!r}"
            )
        normalized = dict(constraint)
        normalized["evidence_citation"] = validate_notes_citation(
            errors,
            file_label,
            where,
            constraint.get("evidence_citation"),
            citation_field="evidence_citation",
        )
        constraints.append(normalized)

    downstream_dependents: list[str] = []
    for index, stage in enumerate(downstream_raw, 1):
        try:
            downstream_dependents.append(f"{int(str(stage)):03d}")
        except (TypeError, ValueError):
            errors.append(f"{file_label}: downstream_dependents[{index}] is not a stage number: {stage!r}")

    findings: list[dict[str, Any]] = []
    for index, finding in enumerate(findings_raw, 1):
        where = f"provenance_findings[{index}]"
        if not isinstance(finding, dict):
            errors.append(f"{file_label}: {where} must be a mapping")
            continue
        for field in ("type", "severity", "summary"):
            if not finding.get(field):
                errors.append(f"{file_label}: {where} missing {field}")
        citations = finding.get("citations", [])
        if not isinstance(citations, list):
            errors.append(f"{file_label}: {where} citations must be a list")
            citations = []
        normalized = dict(finding)
        normalized["citations"] = [
            validate_optional_citation(errors, file_label, f"{where}.citations[{c_index}]", citation)
            for c_index, citation in enumerate(citations, 1)
        ]
        findings.append(normalized)

    if errors:
        return None, errors
    return (
        {
            "candidate_id": candidate_id,
            "parameter_name": parameter_name,
            "origin_claims": origin_claims,
            "constraints": constraints,
            "downstream_dependents": sorted(set(downstream_dependents)),
            "provenance_findings": findings,
            "source_path": env.rel(source),
        },
        [],
    )


def provenance_path_for_parameter(env: Env, entry: dict[str, Any], parameter_name: str) -> Path:
    expected = phase_b_bundle_path(env, str(entry["id"]), parameter_name)
    expected_rel = env.rel(expected)
    for path_text in entry.get("paths", {}).get("provenance") or []:
        if path_text == expected_rel or Path(path_text).name == expected.name:
            return env.abs_from_rel(path_text)
    return expected


def provenance_complete_for_candidate(env: Env, entry: dict[str, Any]) -> bool:
    paths = entry.get("paths", {}).get("provenance") or []
    if not paths:
        return False
    for path_text in paths:
        path = env.abs_from_rel(path_text)
        if not path.exists():
            return False
        data = load_yaml(path)
        if data.get("synthesis_status") != "complete":
            return False
    return True


def ingest_phase_b_syntheses(env: Env, synthesis_paths: list[str]) -> dict[str, Any]:
    normalized_payloads: list[dict[str, Any]] = []
    errors: list[str] = []
    for path_text in synthesis_paths:
        source = resolve_input_path(env, path_text)
        payload, file_errors = normalize_phase_b_synthesis_file(env, source)
        errors.extend(file_errors)
        if payload is not None:
            normalized_payloads.append(payload)
    if errors:
        raise SystemExit("error: malformed Phase B synthesis input:\n- " + "\n- ".join(errors))

    manifest = load_manifest(env)
    candidates = manifest.get("candidates") or {}
    plans: list[tuple[dict[str, Any], dict[str, Any], Path, dict[str, Any]]] = []
    updated_paths: list[str] = []
    completed_candidates: list[str] = []
    for payload in normalized_payloads:
        candidate_id = payload["candidate_id"]
        entry = candidates.get(candidate_id)
        if not entry:
            raise SystemExit(f"error: unknown candidate in synthesis: {candidate_id}")
        if entry.get("duplicate_of"):
            raise SystemExit(f"error: {candidate_id} is an alias of {entry.get('duplicate_of')}; ingest synthesis for the canonical")
        path = provenance_path_for_parameter(env, entry, payload["parameter_name"])
        if not path.exists():
            raise SystemExit(
                f"error: provenance bundle does not exist for {candidate_id} parameter {payload['parameter_name']}: {env.rel(path)}"
            )
        existing = load_yaml(path)
        if existing.get("candidate_id") != candidate_id:
            raise SystemExit(f"error: provenance bundle candidate_id mismatch: {env.rel(path)}")
        if str(existing.get("parameter_name")) != payload["parameter_name"]:
            raise SystemExit(
                f"error: provenance bundle parameter_name mismatch for {env.rel(path)}: {existing.get('parameter_name')!r}"
            )
        plans.append((payload, entry, path, existing))

    for payload, entry, path, existing in plans:
        append_unique_items(existing.setdefault("origin_claims", []), payload["origin_claims"])
        append_unique_items(existing.setdefault("constraints", []), payload["constraints"])
        append_unique_items(existing.setdefault("downstream_dependents", []), payload["downstream_dependents"])
        append_unique_items(existing.setdefault("provenance_findings", []), payload["provenance_findings"])
        existing["synthesis_status"] = "complete"
        existing["synthesis_ingested_at"] = now_iso()
        existing["agent_synthesis_path"] = payload["source_path"]
        write_yaml(path, existing)
        updated_paths.append(env.rel(path))
        entry["phase_b_status"] = "synthesis_partial"
        if provenance_complete_for_candidate(env, entry):
            transition(entry, "provenance_built")
            entry["phase_b_status"] = "synthesis_complete"
            completed_candidates.append(candidate_id)
        candidates[candidate_id] = entry

    render_fit_file(env, manifest)
    render_batches(env, manifest)
    save_manifest(env, manifest)
    return {
        "ingested_syntheses": [payload["source_path"] for payload in normalized_payloads],
        "updated_provenance_paths": updated_paths,
        "completed_candidates": sorted(set(completed_candidates)),
    }


def citation_is_real(value: Any) -> bool:
    placeholder = {"", "none", "n/a", "na", "tbd", "todo", "model memory", "memory", "from memory"}
    if isinstance(value, str):
        return value.strip().lower() not in placeholder
    if isinstance(value, dict):
        flattened = " ".join(str(v).strip() for v in value.values() if v not in (None, ""))
        return flattened.strip().lower() not in placeholder
    return False


def normalize_benchmark_entries(data: Any, file_label: str) -> tuple[list[dict[str, Any]], list[str]]:
    errors: list[str] = []
    if isinstance(data, dict):
        entries_raw = data.get("entries")
        if entries_raw is None and any(key in data for key in ("claim", "value", "source_type", "source_citation")):
            entries_raw = [data]
    elif isinstance(data, list):
        entries_raw = data
    else:
        return [], [f"{file_label}: top-level YAML must be a mapping or list"]
    if not isinstance(entries_raw, list):
        return [], [f"{file_label}: missing list field: entries"]
    entries: list[dict[str, Any]] = []
    for index, entry in enumerate(entries_raw, 1):
        where = f"entries[{index}]"
        if not isinstance(entry, dict):
            errors.append(f"{file_label}: {where} must be a mapping")
            continue
        candidate_id = entry.get("candidate_id")
        family_id = entry.get("family_id")
        if not candidate_id and not family_id:
            errors.append(f"{file_label}: {where} must include candidate_id or family_id")
        if candidate_id and family_id:
            errors.append(f"{file_label}: {where} must key by either candidate_id or family_id, not both")
        for field in ("claim", "value", "source_type", "source_citation", "obtained_note"):
            if entry.get(field) in (None, ""):
                errors.append(f"{file_label}: {where} missing {field}")
        if entry.get("source_type") not in BENCHMARK_SOURCE_TYPES:
            errors.append(
                f"{file_label}: {where} source_type must be one of {', '.join(sorted(BENCHMARK_SOURCE_TYPES))}; got {entry.get('source_type')!r}"
            )
        if not citation_is_real(entry.get("source_citation")):
            errors.append(f"{file_label}: {where} source_citation must be a real URL, DOI, CODATA identifier, or textbook ref+page")
        normalized = dict(entry)
        if not normalized.get("id"):
            normalized["id"] = stable_entry_id(
                "bench",
                [
                    normalized.get("family_id") or normalized.get("candidate_id"),
                    normalized.get("claim"),
                    normalized.get("value"),
                    normalized.get("source_type"),
                    normalized.get("source_citation"),
                ],
            )
        normalized["ingested_at"] = now_iso()
        entries.append(normalized)
    return ([], errors) if errors else (entries, [])


def ingest_benchmarks(env: Env, benchmark_paths: list[str]) -> dict[str, Any]:
    normalized_entries: list[dict[str, Any]] = []
    errors: list[str] = []
    source_paths: list[str] = []
    for path_text in benchmark_paths:
        source = resolve_input_path(env, path_text)
        file_label = env.rel(source)
        if not source.exists():
            errors.append(f"{file_label}: benchmark file does not exist")
            continue
        try:
            data = load_yaml(source)
        except yaml.YAMLError as exc:
            errors.append(f"{file_label}: invalid YAML: {exc}")
            continue
        entries, file_errors = normalize_benchmark_entries(data, file_label)
        errors.extend(file_errors)
        for entry in entries:
            entry["ingested_from"] = file_label
        normalized_entries.extend(entries)
        source_paths.append(file_label)
    if errors:
        raise SystemExit("error: malformed benchmark input:\n- " + "\n- ".join(errors))

    data = ensure_benchmarks_file(env)
    existing_by_id = {str(entry.get("id")): entry for entry in data.get("entries", []) if entry.get("id")}
    inserted = 0
    updated = 0
    unchanged = 0
    for entry in normalized_entries:
        entry_id = str(entry["id"])
        if entry_id not in existing_by_id:
            existing_by_id[entry_id] = entry
            inserted += 1
        elif existing_by_id[entry_id] == entry:
            unchanged += 1
        else:
            existing_by_id[entry_id] = entry
            updated += 1
    data["entries"] = sorted(existing_by_id.values(), key=lambda item: str(item.get("id")))
    write_yaml(env.benchmarks_path, data)
    return {
        "ingested_benchmark_files": source_paths,
        "entries_inserted": inserted,
        "entries_updated": updated,
        "entries_unchanged": unchanged,
        "benchmark_path": env.rel(env.benchmarks_path),
    }


def load_benchmarks_for_candidate(env: Env, candidate_id: str) -> list[dict[str, Any]]:
    data = ensure_benchmarks_file(env)
    return [e for e in data.get("entries", []) if e.get("candidate_id") == candidate_id]


def render_phase_c(env: Env, candidate_id: str) -> dict[str, Any]:
    manifest = load_manifest(env)
    entry = (manifest.get("candidates") or {}).get(candidate_id)
    if not entry:
        raise SystemExit(f"error: unknown candidate: {candidate_id}")
    if entry.get("duplicate_of"):
        raise SystemExit(f"error: {candidate_id} is an alias of {entry.get('duplicate_of')}; render Phase C on the canonical")
    if not entry.get("dry_run") and entry.get("status") != "provenance_built":
        raise SystemExit(
            f"error: {candidate_id} status is {entry.get('status')}; run phase-b-ingest before phase-c-render"
        )
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
    result = run_phase_a(env, args.stages, dry_run=args.dry_run, dry_run_id=args.dry_run_id, prefix=args.prefix)
    print(yaml_block(result))


def cmd_phase_a_ingest(env: Env, args: argparse.Namespace) -> None:
    print(yaml_block(ingest_phase_a_fragments(env, args.fragments)))


def cmd_target_resolve(env: Env, args: argparse.Namespace) -> None:
    out_path = resolve_output_path(env, args.out, TARGET_LAYER_DEFAULT)
    payload = build_target_layer(env, out_path)
    print(yaml_block(payload["summary"] | {"out": env.rel(out_path), "concept_aliases": env.rel(concept_alias_path(env))}))


def cmd_dedup_propose(env: Env, args: argparse.Namespace) -> None:
    out_path = resolve_output_path(env, args.out, "provenance/_dedup_proposal.yaml")
    print(yaml_block(build_dedup_proposal(env, out_path)))


def cmd_apply_alias_map(env: Env, args: argparse.Namespace) -> None:
    print(yaml_block(apply_alias_map(env, args.map_path)))


def cmd_family_build(env: Env, args: argparse.Namespace) -> None:
    out_path = resolve_output_path(env, args.out, "provenance/_family_map.yaml")
    print(yaml_block(build_family_map(env, out_path)))


def cmd_render_critic(env: Env, args: argparse.Namespace) -> None:
    manifest = load_manifest(env)
    path = render_critic_prompt(env, prefix=args.prefix or "phase_a", manifest=manifest)
    print(
        yaml_block(
            {
                "critic_prompt": path,
                "fragment_paths": fragment_paths_for_critic(env),
                "candidate_count": len((manifest.get("candidates") or {})),
            }
        )
    )


def cmd_phase_b_build(env: Env, args: argparse.Namespace) -> None:
    print(yaml_block(build_provenance(env, args.candidate_id)))


def cmd_phase_b_build_all(env: Env, args: argparse.Namespace) -> None:
    print(yaml_block(phase_b_build_all(env, args)))


def cmd_phase_b_ingest(env: Env, args: argparse.Namespace) -> None:
    print(yaml_block(ingest_phase_b_syntheses(env, args.syntheses)))


def cmd_benchmark_ingest(env: Env, args: argparse.Namespace) -> None:
    print(yaml_block(ingest_benchmarks(env, args.benchmarks)))


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
    result = run_phase_a(env, stages, dry_run=True, dry_run_id=dry_run_id, prefix=dry_run_id)
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
        try:
            data = load_yaml(path)
        except (OSError, yaml.YAMLError):
            continue
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
    pa.add_argument("--prefix", default="phase_a")
    pa.add_argument("--dry-run", action="store_true")
    pa.add_argument("--dry-run-id")
    pa.set_defaults(func=cmd_phase_a_scan)

    ingest = sub.add_parser("phase-a-ingest")
    ingest.add_argument("fragments", nargs="+")
    ingest.set_defaults(func=cmd_phase_a_ingest)

    target = sub.add_parser("target-resolve")
    target.add_argument("--out")
    target.set_defaults(func=cmd_target_resolve)

    dedup = sub.add_parser("dedup-propose")
    dedup.add_argument("--out")
    dedup.set_defaults(func=cmd_dedup_propose)

    alias = sub.add_parser("apply-alias-map")
    alias.add_argument("map_path")
    alias.set_defaults(func=cmd_apply_alias_map)

    family = sub.add_parser("family-build")
    family.add_argument("--out")
    family.set_defaults(func=cmd_family_build)

    critic = sub.add_parser("render-critic")
    critic.add_argument("--prefix", default="phase_a")
    critic.set_defaults(func=cmd_render_critic)

    pb = sub.add_parser("phase-b-build")
    pb.add_argument("candidate_id")
    pb.set_defaults(func=cmd_phase_b_build)

    pba = sub.add_parser("phase-b-build-all")
    pba.add_argument("--ids-file")
    pba.add_argument("--limit", type=int)
    pba.add_argument("--force", action="store_true")
    pba.add_argument("--dry", action="store_true")
    pba.set_defaults(func=cmd_phase_b_build_all)

    pbi = sub.add_parser("phase-b-ingest")
    pbi.add_argument("syntheses", nargs="+")
    pbi.set_defaults(func=cmd_phase_b_ingest)

    bi = sub.add_parser("benchmark-ingest")
    bi.add_argument("benchmarks", nargs="+")
    bi.set_defaults(func=cmd_benchmark_ingest)

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
    writes_manifest = args.command in WRITE_COMMANDS
    if args.command == "phase-b-build-all" and getattr(args, "dry", False):
        writes_manifest = False
    require_manifest_lock(args.command, writes_manifest=writes_manifest)
    env = Env(args.config_path)
    args.func(env, args)


if __name__ == "__main__":
    main()
