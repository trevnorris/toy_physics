#!/usr/bin/env python3
"""Check or explicitly accept the whole-file ledger-dimensions pin."""

from __future__ import annotations

import argparse
import hashlib
import re
import sys
from pathlib import Path
from typing import Iterable


SCRIPT_DIR = Path(__file__).resolve().parent
LEDGER_DIR = SCRIPT_DIR.parent
LEDGER_DIMENSIONS_MODULE = SCRIPT_DIR / "ledger_dimensions.py"
ACCEPTED_SHA256_PATH = SCRIPT_DIR / "ledger_dimensions.accepted.sha256"
AUTHORITY_LINE_RE = re.compile(
    r"(?P<digest>[0-9a-f]{64})  ledger_dimensions\.py\n?\Z"
)


class ModulePinError(ValueError):
    """The module pin could not be evaluated."""


class ModulePinMismatch(ModulePinError):
    """The current module bytes are not the deliberately accepted bytes."""


class ModulePinArgumentParser(argparse.ArgumentParser):
    """Report command-line failures with the control's structured marker."""

    def error(self, message: str) -> None:
        self.print_usage(sys.stderr)
        self.exit(2, f"MODULE_PIN_ERROR|argument_error={message}\n")


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(LEDGER_DIR))
    except ValueError:
        return str(path)


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of one file, computed in-process."""

    try:
        content = path.read_bytes()
    except OSError as exc:
        raise ModulePinError(
            f"MODULE_PIN_ERROR|module_read={display_path(path)}|detail={exc}"
        ) from exc
    return hashlib.sha256(content).hexdigest()


def accepted_sha256() -> str:
    try:
        authority_text = ACCEPTED_SHA256_PATH.read_text(encoding="ascii")
    except (OSError, UnicodeError) as exc:
        raise ModulePinError(
            f"MODULE_PIN_ERROR|authority={display_path(ACCEPTED_SHA256_PATH)}"
            f"|detail={exc}"
        ) from exc
    match = AUTHORITY_LINE_RE.fullmatch(authority_text)
    if match is None:
        raise ModulePinError(
            f"MODULE_PIN_ERROR|malformed_authority="
            f"{display_path(ACCEPTED_SHA256_PATH)}"
        )
    return match.group("digest")


def require_accepted_ledger_dimensions(
    consumer_artifacts: Iterable[Path] = (),
) -> str:
    """Reject any byte difference from the separately accepted module digest."""

    artifacts = tuple(Path(path) for path in consumer_artifacts)
    missing = tuple(path for path in artifacts if not path.is_file())
    if missing:
        raise ModulePinError(
            "MODULE_PIN_ERROR|missing_consumer_artifact="
            + ",".join(display_path(path) for path in missing)
        )

    expected = accepted_sha256()
    actual = sha256_file(LEDGER_DIMENSIONS_MODULE)
    if actual != expected:
        artifact_text = ",".join(display_path(path) for path in artifacts)
        raise ModulePinMismatch(
            f"MODULE_PIN_MISMATCH|module={display_path(LEDGER_DIMENSIONS_MODULE)}"
            f"|authority={display_path(ACCEPTED_SHA256_PATH)}"
            f"|expected_sha256={expected}|actual_sha256={actual}"
            f"|consumer_artifacts={artifact_text or '(none)'}"
        )
    return actual


def accept_current_module() -> str:
    """Explicitly replace the accepted digest with the current module digest."""

    digest = sha256_file(LEDGER_DIMENSIONS_MODULE)
    try:
        ACCEPTED_SHA256_PATH.write_text(
            f"{digest}  {LEDGER_DIMENSIONS_MODULE.name}\n",
            encoding="ascii",
        )
    except OSError as exc:
        raise ModulePinError(
            f"MODULE_PIN_ERROR|authority_write="
            f"{display_path(ACCEPTED_SHA256_PATH)}|detail={exc}"
        ) from exc
    return digest


def main() -> int:
    parser = ModulePinArgumentParser(allow_abbrev=False)
    parser.add_argument(
        "--accept",
        action="store_true",
        help="explicitly accept the current ledger_dimensions.py bytes",
    )
    parser.add_argument(
        "--consumer-artifact",
        action="append",
        type=Path,
        default=[],
        help="artifact whose producer consumed ledger_dimensions.py",
    )
    args = parser.parse_args()
    if args.accept and args.consumer_artifact:
        parser.error("--accept cannot be combined with --consumer-artifact")

    try:
        if args.accept:
            digest = accept_current_module()
            print(
                f"MODULE_PIN_ACCEPTED|module="
                f"{display_path(LEDGER_DIMENSIONS_MODULE)}"
                f"|authority={display_path(ACCEPTED_SHA256_PATH)}"
                f"|accepted_sha256={digest}"
            )
        else:
            digest = require_accepted_ledger_dimensions(args.consumer_artifact)
            artifacts = ",".join(
                display_path(path) for path in args.consumer_artifact
            )
            print(
                f"MODULE_PIN_OK|module={display_path(LEDGER_DIMENSIONS_MODULE)}"
                f"|authority={display_path(ACCEPTED_SHA256_PATH)}"
                f"|accepted_sha256={digest}"
                f"|consumer_artifacts={artifacts or '(none)'}"
            )
        return 0
    except ModulePinError as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
