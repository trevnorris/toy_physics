# S12 build run record — decision-list revision 4

2026-10-20. Contract: `../S12_decision_list.md`, revision 4 cleared at `b08673a2`; reviewed build baseline preserved at `04a24c36`. This is the authorized build update, with no new review, engine-output comparison, or commit. Paths below are relative to `research/pde_ledger_v3/` unless absolute.

Both engines give n, every velocity component, and every order-flux component independent live functions and one-sided traces in all three open regions. Order remains a live interior field and zero outside under U1. Conserved products are formed within each region before distributional differentiation. The computed source includes its face-concentrated terms, and total-material face requirements are printed separately. No input field has a concentrated face part: those terms from the previous build were removed, as revision 4 requires. Their possible inclusion remains OPEN outside this object.

Wolfram eliminates the three regional density time derivatives using their material balances. SymPy extracts and removes regional bulk material-balance multiples. Neither multiplies a jumping order field by the global material-balance distribution. Both retain the general interior order, the labelled conditional two-offset G specialization, and dimensional operands/residuals, without choosing jumps or field laws.

Final-code runs were serial, Wolfram first, through the requested pooled guard: 50 GiB cap, zero swap, 32 tasks, verified limits and no execution deadline.

## Wolfram

Script: `mathematica/S12_order_source.wl`. Source SHA256: `f116251b1f47fb31d2b67bac6867f653683fa6a76216b09feab72c4a73bcdda9`.

```bash
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s12 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s12/guard/wolfram_revision4_01 -- math -script /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S12_order_source.wl
```

Guard/child exit **0**, `limitsVerified: true`, `guardReason: null`, stderr **0 bytes**. Workload wall time: **4.618431990005774 s**; guard wall time: **4.701620332984021 s**. Peak guarded-cgroup memory: **223019008 bytes**, from the maximum recorded `memory.peak`.

Transcript: `mathematica/out/S12_order_source.out`, **303204 bytes**. SHA256: `ec0d1acdd0666456cb50ce785b0c5188208b55e88f19b339535a97bcaf7c9055`.

## SymPy

Script: `scripts/S12_order_source.py`. Source SHA256: `871495c15867462eb766da203bd338786d11fdb45474c4abc86b2df79e1cf926`.

```bash
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s12 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s12/guard/sympy_revision4_01 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S12_order_source.py
```

Guard/child exit **0**, `limitsVerified: true`, `guardReason: null`, stderr **0 bytes**. Workload wall time: **1.5678965460101608 s**; guard wall time: **1.6370442509942222 s**. Peak guarded-cgroup memory: **55709696 bytes**, from the maximum recorded `memory.peak`.

Transcript: `scripts/out/S12_order_source.out`, **366976 bytes**. SHA256: `988149e6eb2162c78329eb6cbfa90762bf505fea448c716f7de2f33686036a97`.

## Production checks and scope

Each transcript equals its own guard stdout. Separate structural checks found all requested source, regional, face, material-face-requirement and dimensional output tags; all regional density, velocity and flux identifiers occur in the face operands. Excluded concentrated-input identifiers and output tags are absent. These checks concern output coverage, not a scientific verdict. Requested objects not produced: **none**. No guard refusal, kill, or code-error attempt occurred in this update.

Working scripts and hash/coverage metadata: repository-relative `_scratch/s12/build_work/revision4/`. Guard receipts and resource samples are in the two command log directories; `_scratch/s12/guard` remains the existing alias to `_scratch/s11c/s12_guard`. The guard and committed annex objects were unchanged. Inventories and the remaining OPEN questions were not undertaken. Commits: **none**.
