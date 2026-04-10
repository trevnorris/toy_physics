#!/usr/bin/env python3
"""
3pn_referee_master_sympy_audit.py

Standalone end-to-end SymPy referee audit for the full conservative 3PN derivation chain.

What this file is for
---------------------
This is the top-level SymPy referee runner for the 4d_3pn paper bundle. It replays the
full stage-audit chain used to derive the paper's formulas, executes the stages in isolated
namespaces, and stops immediately on the first failed symbolic identity.

Coverage map
------------
This master audit covers the derivation chain corresponding to the current 3PN summary/writeup:

  • one-body 3PN Schwarzschild gate and repair ledger
  • grouped-P2 kickoff / isotropy bookkeeping
  • comparable-mass scaffold and cubic Legendre machinery
  • COM ordinary/Hamiltonian compiler and exact GR COM target
  • generic-frame COM projection, seed repair, COM-null ideal, and contact/gauge orbit
  • generic-frame target import, Hamiltonian-first lift, and fixed-chart uniqueness
  • grouped real P2 target pack, minimal-lift no-go, and richer exact middle-block closure
  • static P0/geometry completion, sigma-collapse, and pure-kinetic collapse
  • consolidated conservative 3PN theorem ledger

Important scope note
--------------------
This verifies the full encoded symbolic derivation chain inside the declared closure framework
used by the 3PN paper. It is therefore the correct referee script for the paper as written,
but it does not claim to prove statements outside that framework.

Reproducibility / integrity
---------------------------
Each extracted stage carries a pinned SHA256 digest inside the master runner. The script uses
only the Python standard library plus SymPy and does not import any local modules.
"""

from __future__ import annotations

import hashlib
import sys
import time
import traceback
from pathlib import Path
from typing import List, Tuple

ROOT = Path(__file__).resolve().parent
STAGES: List[Tuple[str, str]] = [('3pn_onebody_audit.py', 'fa3887f9994dfe4a30a81000aa6551d097ea27ed0039670cded88f43a9f221c6'), ('3pn_grouped_p2_audit.py', '9eda9aff4a120717cdab2b3b4d32348e02e25290b0fe2be7362a01359185db18'), ('3pn_comparable_mass_audit.py', 'a2c24af46549eb1fb17dc563aac96cd2faf2009f2702c4ebf92d7bf677bbe218'), ('3pn_com_linear_map_audit.py', 'ffe6f85da870d09d59178358b062b5accd4792e7a17a1dcb98f98ece9e0ea331'), ('3pn_com_gr_target_audit.py', '69efe2d65f00da96f902409135ea84a11ec4b4c88014694108b13d42bda62913'), ('3pn_generic_frame_com_projection_audit.py', '578ed9231f3d225aa678a30a3d0164547e77e817c37f5243fb32dd94ad7309d5'), ('3pn_seed_repair_and_com_null_ideal_audit.py', '558a1741ad0f361a1e132d04fdba33147febb94838a2d3cc042b14a98259b06c'), ('3pn_contact_generator_and_gauge_orbit_audit.py', 'fbe8d1d6215c0a86207cdf4a57c1fad9f32865e377db105558505e209aa15325'), ('3pn_generic_frame_target_import_audit.py', 'a05ee0283e107068398ad8165a634fbd835cf2270b4ac13e4c5a78cf478c6d28'), ('3pn_hamiltonian_level_lift_audit.py', 'b5f68dc6e4f9b6e7ab5d012e02abb071cc786a4721bb4f8a05d5f1bfa2ba60a9'), ('3pn_generic_frame_hamiltonian_compiler_audit.py', 'e4776bd914c2d7b5c2c9dc3705cee8a2017e3e830efd5153202232f93f71a415'), ('3pn_grouped_p2_target_pack_audit.py', '796309e827ac3c5d76dc3634dd6891ca2b12f38281a4b4f3516b75149ff47cdf'), ('3pn_grouped_p2_middle_block_test_audit.py', '3311e142bc60496f780515df157664b8b1315794a8db6530fd01b5371444f3bc'), ('3pn_grouped_p2_richer_lift_audit.py', '7bc1e7a2e365369cf8804f6dd012812c41b757d0c00146780717ec6e97d247fe'), ('3pn_static_p0_geometry_counterterm_audit.py', '2324a01d6eaf9d3c005676bda5d79aedb57d16de90bdfa8e30c3f414a1e99011'), ('3pn_sigma_collapse_and_unique_geometry_completion_audit.py', 'f7bb035c045f8feb2a073a8c2306a2fb7972ae0a33da43dc281ba1ffd5d3bee0'), ('3pn_pure_kinetic_collapse_audit.py', '02edbdea3df9a2e703cbd3c6ce434883d6587b9553a12bb6b314cefdc0eec8a9'), ('3pn_conservative_theorem_audit.py', 'e2249faeb8ef62e355fdb95a6ae665989d6a1ba3812b70e047811dc43835e807')]


def banner(title: str) -> None:
    line = "=" * 100
    print("\n" + line)
    print(title)
    print(line)


def run_stage(index: int, total: int, filename: str, expected_sha: str) -> float:
    path = ROOT / filename
    source = path.read_text()
    actual_sha = hashlib.sha256(source.encode()).hexdigest()
    if actual_sha != expected_sha:
        raise RuntimeError(f"Stage hash mismatch for {filename}: {actual_sha} != {expected_sha}")
    banner(f"STAGE {index}/{total}: {filename}\nSHA256: {expected_sha}")
    ns = {"__name__": "__main__", "__file__": str(path), "__package__": None}
    t0 = time.perf_counter()
    code = compile(source, str(path), "exec")
    exec(code, ns, ns)
    dt = time.perf_counter() - t0
    print(f"\n[stage ok] {filename} finished in {dt:.3f} s")
    return dt


def main() -> None:
    total = len(STAGES)
    banner("3PN REFEREE MASTER SYMPY AUDIT — START")
    print(f"Stage files: {total}")
    print("This run will stop immediately on the first failed symbolic identity.")
    timings = []
    t0 = time.perf_counter()
    for i, (filename, sha) in enumerate(STAGES, start=1):
        try:
            dt = run_stage(i, total, filename, sha)
            timings.append((filename, dt))
        except Exception as exc:  # pragma: no cover
            print(f"\n[stage failed] {filename}: {exc}")
            traceback.print_exc()
            raise
    total_dt = time.perf_counter() - t0
    banner("3PN REFEREE MASTER SYMPY AUDIT — ALL STAGES PASSED")
    print(f"Total wall time: {total_dt:.3f} s")
    print("Stage timings:")
    for filename, dt in timings:
        print(f"  - {filename}: {dt:.3f} s")
    print("\nFinal status: every stage file completed without assertion failure.")
    print("This is the stand-alone referee audit for the full conservative 3PN derivation chain.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("\n[FAIL] Master referee audit stopped on an error:", exc, file=sys.stderr)
        raise
