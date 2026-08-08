"""Regression tests for the consumer-derived engine fence."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from registry_import_fence import discover_registry_importing_engines


HERE = Path(__file__).resolve().parent


def test_registry_import_fence_is_computed_from_import_syntax(tmp_path: Path) -> None:
    scripts = tmp_path / "scripts"
    mathematica = tmp_path / "mathematica"
    scripts.mkdir()
    mathematica.mkdir()
    imported = scripts / "imported_audit.py"
    imported.write_text("import registry_read\n", encoding="utf-8")
    from_imported = scripts / "from_imported_audit.py"
    from_imported.write_text("from registry_read import load_registry\n", encoding="utf-8")
    comment_only = scripts / "comment_only_audit.py"
    comment_only.write_text("# import registry_read\n", encoding="utf-8")
    unrelated = mathematica / "unrelated_audit.wl"
    unrelated.write_text("Print[registry_read]\n", encoding="utf-8")

    observed = discover_registry_importing_engines(tmp_path)
    assert set(observed) == {imported.resolve(), from_imported.resolve()}


def test_registry_import_fence_command_answers_the_current_set() -> None:
    completed = subprocess.run(
        [sys.executable, str(HERE / "registry_import_fence.py"), "--list"],
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert completed.returncode == 0
    assert not completed.stderr
    assert tuple(Path(line) for line in completed.stdout.splitlines()) == (
        discover_registry_importing_engines()
    )
