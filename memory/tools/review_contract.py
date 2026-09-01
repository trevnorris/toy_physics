"""Versioned constants shared by the Grok reviewer and deterministic verifier."""

from __future__ import annotations

import hashlib
import json


GROK_REVIEW_MODEL = "grok-4.6"

REVIEW_PROMPT = (
    "Act only as an independent reviewer; do not rewrite the candidate. The complete immutable task, "
    "schema, applicable frozen prompts, sealed semantic inputs, and candidate are appended below. "
    "Check factual fidelity, citation support and anchor scope, lifecycle and "
    "evidence labels, stable IDs, omissions, conflicts, and any legacy-migration requirements. The first "
    "result must use the required JSON schema. Give concise findings ordered by severity, citing candidate "
    "headings/lines and sealed source paths. A PASS may include only non-blocking suggestions."
)

REVIEW_SCHEMA = json.dumps({
    "type": "object",
    "properties": {
        "verdict": {"type": "string", "enum": ["PASS", "FAIL"]},
        "findings": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "severity": {"type": "string", "enum": ["blocking", "major", "minor", "note"]},
                    "summary": {"type": "string"},
                    "candidate_location": {"type": "string"},
                    "source_evidence": {"type": "string"},
                },
                "required": ["severity", "summary", "candidate_location", "source_evidence"],
                "additionalProperties": False,
            },
        },
    },
    "required": ["verdict", "findings"],
    "additionalProperties": False,
}, separators=(",", ":"), sort_keys=True)


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


REVIEW_PROMPT_SHA256 = sha256_text(REVIEW_PROMPT)
REVIEW_SCHEMA_SHA256 = sha256_text(REVIEW_SCHEMA)
REVIEW_CONTRACT_SHA256 = hashlib.sha256(json.dumps({
    "model": GROK_REVIEW_MODEL,
    "prompt_sha256": REVIEW_PROMPT_SHA256,
    "schema_sha256": REVIEW_SCHEMA_SHA256,
}, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
