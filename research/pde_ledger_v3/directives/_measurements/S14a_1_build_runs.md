# S14a part 1 build runs — 2026-10-21

Deliverable: two independently constructed CAS engines, their final-code guarded transcripts, and this run record. Contract: `directives/S14a_1_decision_list.md`, revision 3, at HEAD `e9eed6145b4b170719d0c852e2199d44e9a86da8`. No commits were made.

The Wolfram source was written before the SymPy source. Neither engine reads the other engine or its output. No cross-engine comparison or scientific interpretation was performed. Formal integrals retain the live fields; velocity objects carry their nonzero-density domain. All contract OPEN items remain OPEN.

The requested `_scratch/s14a/guard` scratch path is a symlink to `_scratch/s11c/s14a`, satisfying the unchanged guard's log-root restriction. Runs were sequential, without time limits; each used pool `s14a`, a 50 GiB cap, zero swap and the default 32-task cap.

**SymPy**

- Script: `research/pde_ledger_v3/scripts/S14a_1_projection.py`; SHA-256 `de2cca647e1eb42e14cce81e0cad295f5a5637d7b130dfdfae5e27a50876db76`.
- Transcript: `research/pde_ledger_v3/scripts/out/S14a_1_projection.out`; SHA-256 `1eefeb1d1a6eeed8f2fb697d8d7c659909b46075ae0ba1de832558d7e9997f24`.
- Guard outcome: exit `0`; child exit `0`; `limitsVerified = true`; `guardReason = null`; stderr `0` bytes.
- Runtime: child `1.9185772690107115` s; guard launcher `1.9864753239962738` s. Peak cgroup memory: `61476864` bytes (whole guarded cgroup, not process RSS).
- Requested objects not produced: none.

```bash
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s14a --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s14a/guard/sympy_01 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S14a_1_projection.py
```

**Wolfram**

- Script: `research/pde_ledger_v3/mathematica/S14a_1_projection.wl`; SHA-256 `1022c6d0a8d73f6901f6767df9c48fb83d5971b93f65864fb2aef5b5a2f768c5`.
- Transcript: `research/pde_ledger_v3/mathematica/out/S14a_1_projection.out`; SHA-256 `2d7bedef5cbd90d8105b3729f37492116bf5ec2c9d624d75184e30794e91940f`.
- Guard outcome: exit `0`; child exit `0`; `limitsVerified = true`; `guardReason = null`; stderr `0` bytes.
- Runtime: child `4.825878250005189` s; guard launcher `4.896841159003088` s. Peak cgroup memory: `227004416` bytes (whole guarded cgroup, not process RSS).
- Requested objects not produced: none.

```bash
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s14a --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s14a/guard/wolfram_01 -- math -script /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S14a_1_projection.wl
```

Validation: Python source parsed; Wolfram delimiters checked; both final sources ran through the guard. Each transcript is nonempty, includes the requested object tags and dimensional objects, and is byte-identical to its own captured stdout. No error-token matches were found. The repository-referenced `derived_or_declared.py` and `engine_output_checks.py` tools are absent from this checkout; no substitute review was launched.
