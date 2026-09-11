#!/usr/bin/env python3
"""Regression and adversarial checks for the focused S10 stratum comparison."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import tempfile

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "scripts"))
from S10_anisotropic_strata_comparator import compare
from S10_cross_engine_comparator import read_transcript


def main():
    py = BASE / "scripts/out/S10_anisotropic_strata_sympy_audit.out"
    wl = BASE / "mathematica/out/S10_anisotropic_strata_mathematica_audit.out"
    checks = []
    result = compare([py, wl])
    assert result["status"] == "PASS" and len(result["joined"]) == 10
    checks.append({"name": "canonical_focused_pair", "outcome": "PASS", "root_cases": 10})

    with tempfile.TemporaryDirectory(prefix="s10-strata-") as temp:
        changed = Path(temp) / "changed.out"
        text = wl.read_text()
        changed.write_text(text.replace("STRATUM1", "STRATUM_TEMP")
                           .replace("STRATUM2", "STRATUM1").replace("STRATUM_TEMP", "STRATUM2"))
        assert compare([py, changed])["status"] == "PASS"
        checks.append({"name": "stratum_number_permutation", "outcome": "PASS"})
        original = py.read_text()
        target = "PY_S10_XFORM_ANISO_D3_Q8_STRATUM2_ROOT3_N3_TRANSVERSE_NULLITY: 1"
        assert target in original
        variants = {
            "corrupt_transverse_count": original.replace(target, target[:-1] + "0"),
            "omit_perpendicular_rerun": "\n".join(line for line in original.splitlines()
                if not line.startswith("PY_S10_XFORM_ANISO_D3_Q8_STRATUM2_")) + "\n",
            "corrupt_root_list": original.replace(
                "PY_S10_XFORM_ANISO_D3_Q8_STRATUM2_ROOT_ORDERING: (0, mu_R/rho_br, mu_R/(rho_br*s_rho))",
                "PY_S10_XFORM_ANISO_D3_Q8_STRATUM2_ROOT_ORDERING: (0, mu_R/rho_br, -mu_R/(rho_br*s_rho))"),
        }
        for name, mutated in variants.items():
            assert mutated != original, name
            changed.write_text(mutated)
            try:
                compare([changed, wl])
            except ValueError as error:
                checks.append({"name": name, "outcome": "REJECTED", "reason": str(error)})
            else:
                raise AssertionError(f"mutation accepted: {name}")

    for engine, current, old in [
        ("PY", py, BASE / "scripts/out/S10_brane_mode_spectrum_sympy_audit.out"),
        ("WL", wl, BASE / "mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out"),
    ]:
        now, before = read_transcript(current), read_transcript(old)
        checked = 0
        for d in (3, 4):
            for root in (1, 2, 3):
                for field in ("N1_MATRIX", "N2_RANK", "N2_NULLITY", "N3_STACKED_MATRIX",
                              "N3_STACKED_RANK", "N3_TRANSVERSE_NULLITY", "N4_NULLITY_DIFFERENCE", "N7_BASIS_COUNT"):
                    name = f"S10_XFORM_ANISO_D{d}_ROOT{root}_{field}"
                    assert now.records[name].raw == before.records[name].raw, name
                    checked += 1
        checks.append({"name": engine + "_generic_outputs_unchanged", "outcome": "PASS", "rows": checked})

    export = BASE / "scripts/S10_exports.py"
    before = hashlib.sha256(export.read_bytes()).hexdigest()
    guard = subprocess.run([sys.executable, str(BASE / "scripts/S10_brane_mode_spectrum_sympy_audit.py"),
                            "--package", "XFORM_ANISO"], capture_output=True, text=True)
    assert guard.returncode == 2 and "requires --no-export" in guard.stderr
    assert hashlib.sha256(export.read_bytes()).hexdigest() == before
    checks.append({"name": "subset_export_guard", "outcome": "REJECTED", "export_sha256": before})
    print(json.dumps({"status": "PASS", "checks": checks}, indent=2))


if __name__ == "__main__":
    main()
