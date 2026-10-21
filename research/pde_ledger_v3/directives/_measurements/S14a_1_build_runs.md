# S14a part 1 build runs — 2026-10-21

Deliverable: two independently constructed CAS engines, their final-code guarded transcripts, and this run record. Contract: `directives/S14a_1_decision_list.md`, revision 3, at HEAD `e9eed6145b4b170719d0c852e2199d44e9a86da8`. No commits were made.

The Wolfram source was written before the SymPy source. Neither engine reads the other engine or its output. No cross-engine comparison or scientific interpretation was performed. Formal integrals retain the live fields; velocity objects carry their nonzero-density domain. All contract OPEN items remain OPEN.

The requested `_scratch/s14a/guard` scratch path is a symlink to `_scratch/s11c/s14a`, satisfying the unchanged guard's log-root restriction. Runs were sequential, without time limits; each used pool `s14a`, a 50 GiB cap, zero swap and the default 32-task cap.

**SymPy**

- Repair scope: C1 only, against reviewed baseline `eaccd417`. The comparison description now states that the residual tests the time and in-plane endpoint pieces, while `bulk_parts` and `atoms` are shared by its operands and are not tested. No physics object changes.
- Script: `research/pde_ledger_v3/scripts/S14a_1_projection.py`; SHA-256 `3d70f411dc894f22dc416b5cae7f7a09cbf8a3603b92a9a8e70ac7639da3e2f9`.
- Transcript: `research/pde_ledger_v3/scripts/out/S14a_1_projection.out`; SHA-256 `3bfb17aa5205e7f0ddec1eea44c920d4f0bb4e3dffe4ccf6fcc74abe88fbe0db`.
- Guard outcome: exit `0`; child exit `0`; `limitsVerified = true`; `guardReason = null`; stderr `0` bytes.
- Runtime: child `2.0343507020152174` s; guard launcher `2.0996876590070315` s. Peak cgroup memory: `62517248` bytes (whole guarded cgroup, not process RSS).
- Requested objects not produced: none. The transcript is byte-identical to its captured stdout. Compared with the reviewed SymPy transcript, only the four comparison-description lines changed; every other byte is unchanged.
- No commit, review round or cross-engine comparison was performed in this repair; contract OPEN items remain OPEN.

```bash
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s14a --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s14a/guard/sympy_C1_01 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S14a_1_projection.py
```

**Wolfram**

- Repair scope: G1 only, against reviewed baseline `eaccd417`. The order entry's `"source"` now holds `suppliedSource`, the projected `n Γ_B` source plus the full face parts that appears on the right of the projected order balance. Order lost and the face channel minus order lost use that entry. Divergence objects keep their expressions and are printed under divergence names; the projected source and its face parts are named explicitly. The material source entry and every other object retain their previous expressions.
- Script: `research/pde_ledger_v3/mathematica/S14a_1_projection.wl`; SHA-256 `365be5d41d3c1f34ddb727bef8a7b49b02bc39e2122a85644aa3dcd2a5870905`.
- Transcript: `research/pde_ledger_v3/mathematica/out/S14a_1_projection.out`; SHA-256 `7911ff19f667a388a2a54a3eb0d076ed06143d791ddb2b5bccfe397588652d5b`. It is byte-identical to the final-code guarded stdout.
- Guard outcome: exit `0`; child exit `0`; `limitsVerified = true`; `guardReason = null`; stderr `0` bytes. Pool `s14a`; memory cap `53687091200` bytes; zero swap; 32 tasks; both Wolfram thread controls retained; no time limit.
- Runtime: child `5.219873964990256` s; guard launcher `5.2969460879976396` s. Peak cgroup memory: `229904384` bytes (whole guarded cgroup, not process RSS).
- Requested objects not produced: none. Live fields, moving faces, the spherical specialization and the nonzero-density conditions are retained.
- Validation: after accounting for the source/divergence label changes, only the order-lost and face-channel-minus-order-lost payloads differ from the reviewed Wolfram transcript; every other line is byte-identical. The order-lost payload is byte-identical to the retained supplied-source form. No Wolfram message lines were found. No commit, review round or cross-engine comparison was performed; all contract OPEN items remain OPEN.

```bash
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s14a --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s14a/guard/wolfram_G1_01 -- math -script /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S14a_1_projection.wl
```

Validation: Python source parsed; Wolfram delimiters checked; both final sources ran through the guard. Each transcript is nonempty, includes the requested object tags and dimensional objects, and is byte-identical to its own captured stdout. No error-token matches were found. The repository-referenced `derived_or_declared.py` and `engine_output_checks.py` tools are absent from this checkout; no substitute review was launched.
