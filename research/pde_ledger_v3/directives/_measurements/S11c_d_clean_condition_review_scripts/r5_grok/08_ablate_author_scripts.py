#!/usr/bin/env python3
"""FORM-ablate the v5 author scripts on copies. A typed script stays byte-identical."""
import shutil
import subprocess
import sys
from pathlib import Path

src = Path(
    "/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/"
    "S11c_d_clean_condition_v5_author_scripts"
)
here = Path("/tmp/s11cd_clean_review_r5_grok/author_ablation")
here.mkdir(parents=True, exist_ok=True)

# Copy
for name in [
    "bulk_pullback_pairing_audit.py",
    "round_spin_carrier_audit.py",
]:
    shutil.copy(src / name, here / name)

# Ablation 1: flip the wave-equation sign (FORM of the acoustics)
bulk = (here / "bulk_pullback_pairing_audit.py").read_text()
bulk_abl = bulk.replace(
    "wave_equation = sp.diff(phi, t, 2) - c_s0**2 * sum(",
    "wave_equation = sp.diff(phi, t, 2) + c_s0**2 * sum(",
)
(here / "bulk_pullback_pairing_audit_ablated.py").write_text(bulk_abl)

# Ablation 2: drop the polar-vector minus in inversion_parity (FORM of the
# transformation law)
rnd = (here / "round_spin_carrier_audit.py").read_text()
rnd_abl = rnd.replace(
    "sp.expand(-component.subs({x: -x, y: -y, z: -z}, simultaneous=True))",
    "sp.expand(component.subs({x: -x, y: -y, z: -z}, simultaneous=True))",
)
(here / "round_spin_carrier_audit_ablated.py").write_text(rnd_abl)


def run(path):
    r = subprocess.run(
        [sys.executable, str(path)],
        capture_output=True,
        text=True,
        timeout=120,
    )
    return r.returncode, r.stdout, r.stderr


for label, path in [
    ("BULK_BASE", here / "bulk_pullback_pairing_audit.py"),
    ("BULK_ABL", here / "bulk_pullback_pairing_audit_ablated.py"),
    ("ROUND_BASE", here / "round_spin_carrier_audit.py"),
    ("ROUND_ABL", here / "round_spin_carrier_audit_ablated.py"),
]:
    code, out, err = run(path)
    out_path = here / f"{label}.stdout.txt"
    out_path.write_text(out)
    print(f"=== {label} exit={code} ===")
    print(out)
    if err.strip():
        print("STDERR", err[:500])

base_b = (here / "BULK_BASE.stdout.txt").read_text()
abl_b = (here / "BULK_ABL.stdout.txt").read_text()
base_r = (here / "ROUND_BASE.stdout.txt").read_text()
abl_r = (here / "ROUND_ABL.stdout.txt").read_text()
print("BULK_BYTE_IDENTICAL", base_b == abl_b)
print("ROUND_BYTE_IDENTICAL", base_r == abl_r)
