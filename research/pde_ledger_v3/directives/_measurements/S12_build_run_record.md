# S12 build run record

2026-10-20. Contract: `../S12_decision_list.md`, authorized after review at `a48ce830`; continuation with full access authorized by the user. Deliverable: two independently constructed engines, final-code guarded transcripts, and this run record. No comparison, review, or commit was performed.

Paths below are relative to `research/pde_ledger_v3/` unless absolute. Both final runs used the requested pooled guard, 50 GiB cap, zero swap, 32 tasks, verified limits, and no execution deadline. Jobs ran serially, Wolfram first. Each transcript is byte-identical to its own guard stdout; coverage checks found no missing requested output tags. These checks concern production of artifacts, not their scientific agreement or review status.

## Wolfram

Script: `mathematica/S12_order_source.wl`.
Source SHA256: `12e24d6888fa31088c6ce6ceb433e269050abef7a6c7832181c80c192398e1bf`.

```bash
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s12 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s12/guard/wolfram_03 -- math -script /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S12_order_source.wl
```

Guard outcome: exit **0**, child exit **0**, `limitsVerified: true`, `guardReason: null`, stderr **0 bytes**. Guard wall time: **5.040406816988252 s**; child workload wall time: **4.925288975005969 s**. Peak guarded-cgroup memory: **222760960 bytes**, from the maximum recorded `memory.peak`.

Transcript: `mathematica/out/S12_order_source.out`, **113412 bytes**.
Transcript SHA256: `b85dc6565c77bf0df279f458db34aba30b765a73dfe46a3de8830ac7fdf3753a`.

## SymPy

Script: `scripts/S12_order_source.py`.
Source SHA256: `60618f4bfb4b663d6e973afcb20562d4a1a62098dd3448ca721e5ca5f23a2300`.

```bash
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s12 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s12/guard/sympy_02 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S12_order_source.py
```

Guard outcome: exit **0**, child exit **0**, `limitsVerified: true`, `guardReason: null`, stderr **0 bytes**. Guard wall time: **0.8368255689856596 s**; child workload wall time: **0.7669781669974327 s**. Peak guarded-cgroup memory: **53604352 bytes**, from the maximum recorded `memory.peak`.

Transcript: `scripts/out/S12_order_source.out`, **117204 bytes**.
Transcript SHA256: `c3ca625fc035168ab520ebd0eb0fa20c5830e71d0c139302990a88a5fc556646`.

## Coverage and retained scope

Both engines emitted the general required source in the two exterior regions, the interior and both faces, including face-concentrated contributions; the labelled conditional two-offset `G` specialization on the same regions and faces; and dimensional operands and residuals for both balances and source objects. **Requested objects not produced: none.** No forms were chosen for the live fields. Objects without a second route are printed with an explicit comparison-availability note. Contract OPEN items remain untouched; inventories were not undertaken.

Working files and mechanical coverage/hash metadata: repository-relative `_scratch/s12/build_work/`. Guard receipts and resource samples: the two log directories in the commands. `_scratch/s12/guard` resolves to `_scratch/s11c/s12_guard` to satisfy the existing guard's log-path restriction; the guard was not changed. The earlier sandbox-blocked attempt and the admitted `wolfram_02` code-error attempt are superseded, not final-code transcripts. The Wolfram bug was use of `UnitStep` instead of the distributional `HeavisideTheta`; its partial transcript remains in `_scratch/s12/build_work/wolfram_02.partial.out`. No completed run was replayed.

Commits: **none**.
