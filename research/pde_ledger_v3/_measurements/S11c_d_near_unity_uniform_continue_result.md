# Selected uniform near-unity check completed — 2026-10-01

**All scheduled selected-uniform checks passed.** Both transverse polarizations
remain finite, carry nonzero current and have zero selected face drives at the
tested effective sound speeds on both sides of, and at, their modal/acoustic
matches. This removes the observed denominator-decision blocker. It does not
measure defect leakage or establish the flowing, primitively calibrated model.

The completion record (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11c_d_near_unity_uniform_continue_completion.json`)
contains the inspected receipts, preservation checks, resources, certificate
summary and per-point coverage. The immutable continuation is
`_scratch/s11c/s11c-parallel-near-unity-20261001/uniform-continuation-01`.
The earlier partial result (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11c_d_near_unity_uniform_result.md`) and its 44
domain refusals remain unchanged as history.

## Saved outputs and replay limit

The completed outputs below are kept as their original JSON bytes. They support this report; they are not a
self-contained replay package. Intermediate caches remain only on local disk under `_scratch/`. A new run
needs those caches regenerated from the kept b/c2 engines and the S10/S11 export chain repaired; see
[STATUS, open questions](../../../STATUS.md#open-questions-and-owners). No rerun occurred during cleanup.

- [Aggregate checks and results](S11c_d_near_unity_uniform_output/checks.json).
- [Exact denominator certificates](S11c_d_near_unity_uniform_output/nonzero-certificate-summary.json) and
  [certificate self-checks](S11c_d_near_unity_uniform_output/certificate-self-checks.json).
- Saved grazing limits: [LEFT minus](S11c_d_near_unity_uniform_output/LEFT-limits-sign--1.json),
  [LEFT plus](S11c_d_near_unity_uniform_output/LEFT-limits-sign-1.json),
  [RIGHT minus](S11c_d_near_unity_uniform_output/RIGHT-limits-sign--1.json),
  [RIGHT plus](S11c_d_near_unity_uniform_output/RIGHT-limits-sign-1.json).

| Schedule index | LEFT point output | RIGHT point output |
| --- | --- | --- |
| 0 | [LEFT](S11c_d_near_unity_uniform_output/point-00-LEFT.json) | [RIGHT](S11c_d_near_unity_uniform_output/point-00-RIGHT.json) |
| 1 | [LEFT](S11c_d_near_unity_uniform_output/point-01-LEFT.json) | [RIGHT](S11c_d_near_unity_uniform_output/point-01-RIGHT.json) |
| 2 | [LEFT](S11c_d_near_unity_uniform_output/point-02-LEFT.json) | [RIGHT](S11c_d_near_unity_uniform_output/point-02-RIGHT.json) |
| 3 | [LEFT](S11c_d_near_unity_uniform_output/point-03-LEFT.json) | [RIGHT](S11c_d_near_unity_uniform_output/point-03-RIGHT.json) |
| 4 | [LEFT](S11c_d_near_unity_uniform_output/point-04-LEFT.json) | [RIGHT](S11c_d_near_unity_uniform_output/point-04-RIGHT.json) |
| 5 | [LEFT](S11c_d_near_unity_uniform_output/point-05-LEFT.json) | [RIGHT](S11c_d_near_unity_uniform_output/point-05-RIGHT.json) |
| 6 | [LEFT](S11c_d_near_unity_uniform_output/point-06-LEFT.json) | [RIGHT](S11c_d_near_unity_uniform_output/point-06-RIGHT.json) |
| 7 | [LEFT](S11c_d_near_unity_uniform_output/point-07-LEFT.json) | [RIGHT](S11c_d_near_unity_uniform_output/point-07-RIGHT.json) |
| 8 | [LEFT](S11c_d_near_unity_uniform_output/point-08-LEFT.json) | [RIGHT](S11c_d_near_unity_uniform_output/point-08-RIGHT.json) |
| 9 | [LEFT](S11c_d_near_unity_uniform_output/point-09-LEFT.json) | [RIGHT](S11c_d_near_unity_uniform_output/point-09-RIGHT.json) |
| 10 | [LEFT](S11c_d_near_unity_uniform_output/point-10-LEFT.json) | [RIGHT](S11c_d_near_unity_uniform_output/point-10-RIGHT.json) |
| 11 | [LEFT](S11c_d_near_unity_uniform_output/point-11-LEFT.json) | [RIGHT](S11c_d_near_unity_uniform_output/point-11-RIGHT.json) |

## Actual coverage and result

This is LAB_HELD/RHO4_CONSTANT, omega=3, tangents (1/5,1/10), strict rest bulk,
with the original physical inputs and a selected effective-bulk-speed family.
The completed schedule has 12 exact speeds, each evaluated at both ends and in
both normal directions: 24 end evaluations and 48 sign checks. Each sign check
uses the two source-selected transverse polarizations. Four matching sign
checks were restored from the original run; 44 unfinished checks are new.

| End | Exact modal/acoustic matching speed | Selected normal momentum | Current eigenvalue magnitudes |
| --- | --- | --- | --- |
| LEFT | sqrt(3/2) | +/-sqrt(595)/10 | 0.2195335965177084; 32.930039477656265 |
| RIGHT | sqrt(150/101) | +/-sqrt(601)/10 | 0.22063771209836272; 33.42661338290195 |

The selected lift has rank two throughout. The current eigenvalue signs follow
the independently differentiated dispersion direction, with normalized Grams
equal to +/-identity to floating roundoff (about 2.2e-16). Current condition
numbers are 150 and 151.5; lift condition numbers are approximately 12.2474 and
12.3085. Those current values and conditions are unchanged across the sampled
sound speeds within each end. Full-pencil rank three/nullity two is a numerical
estimate with its saved tolerance, not an exact complete-mode classification.

All 1,536 named face-channel records have literal exact-zero selected rows and
zero numerical relative residual. Their full native rows are finite. The three
selected normalized bulk-normal/depth/interface forms have zero numerical
entries; recorded exact form differences from the corresponding grazing limits
are zero. Selected harmonic-leg conjugacy residuals are also zero.

The four original end/sign grazing calculations are reused. Their radiating
and evanescent selected limits agree exactly in the saved residuals. Raw zero
denominators and `zoo` from direct substitution at grazing remain recorded;
the source-local limits supply those points. The continuation does not turn
raw singular evaluations into finite values by ignoring them.

Across sign checks, LEFT has 8 evanescent-depth, 2 exact-grazing and 14
radiating-depth points; RIGHT has 14, 2 and 8 respectively. These labels describe
the bulk depth root at the selected transverse mode, not a defect-scattering
channel census. Finite samples and the selected limits are not an interval-wide
smoothness or absence-of-leakage theorem.

## Controls actually reached

- Native-pencil omission changes its residual by 0.04136–0.04159; omission of a
  contributing current cross term changes the normalized Gram by
  0.004107–0.004148.
- All 48 orientation reversals disagree with the independently derived group
  direction. All 192 addressed eW outward-velocity omissions move by 1.
- Both previously unvisited native sheet tests ran at schedule index 0,
  positive normal sign, an evanescent-depth point. The two faces, two harmonic
  legs and eW/theta unit columns give eight probes per end. Reversing the depth
  sheet moves the actual native row actions by 1.5319–1.8167; the reassembled
  closure residuals are literal zero. These unit columns are controls, not
  modes; this is not sheet coverage at every point or a retarded-equivalence
  proof.

The controls demonstrate sensitivity of the addressed source/current/face
operations. They do not convert zero uniform drives into zero defect loss.

## Exact domain certificates and preservation

All 1,354 distinct saved denominator values were inspected across 6,314 uses:
705 have an explicit original nonzero property; 649 previously undecided values
have saved exact real/imaginary components, literal zero reconstruction, real
components and a strictly signed nonzero component. For example, the original
refusal

`2250 + sqrt(266)*(-1365 + 900*I)/4 + 7500*I`

has saved imaginary component `225*sqrt(266)+7500`, positive, and reconstruction
residual `0`. No certificate uses numerical smallness. All 649 signed witnesses
have both component signs resolved; the self-checks reject zero, infinity and
an unbound symbol. Saved certificate metadata and the fail-closed implementation
were inspected without independently re-evaluating the scientific expressions.

All 27 old complete returns were restored without calling their functions.
All 32 prior files, totalling 46,689,296 bytes, were copied byte-identically,
including the 44 original refusals. The four completed matching sign results
are unchanged in the new point JSON. The original source reconstructions,
profile/native joins and exact limits were not replayed.

There are 45 new complete journal operations: 44 unfinished sign checks plus
the aggregate continuation. The 1,916 new and 5,383 preserved prior opaque blobs
pass their stored length/hash checks; every journal and certificate reference
joins. All 70 actual argument records report exact structural identity:
12 also have identical serialized bytes, while 58 use structural equality.
Inspection read only the literal Boolean serialization metadata for these
records, without unpickling their scientific operands. This is saved provenance
evidence, not independent recomputation of the restored expressions.

All 705 metadata/byte checks passed, including 81 posthash routes and 40 source
snapshots. Strict scientific/guard/coordinator/watcher stderr is empty; stdout
is byte-identical to checks.json. No incomplete or failed continuation operation
remains. The result tree has 64 files / 93,592,964 bytes including the old copy.

## Cost, limits and next dependency

Worker time was 1,845.891 seconds (30.8 minutes); peak cgroup memory was
458,420,224 bytes (about 437 MiB). All 924 resource samples have zero swap and
memory-event counters. Host available memory stayed above 21,497,450,496 bytes.
Actual enforcement verified 8 GiB cgroup/native caps in the 16 GiB aggregate
pool, CPU 15, one thread, 32 tasks, a 4 GiB host reserve and desktop priority.
`RuntimeMaxUSec=infinity` and `Restart=no`; no deadline or resource stop occurred.

Both original method verdicts remain **CLEAR FOR THIS SELECTED UNIFORM METHOD**,
packet `bf55b921706cc1d2c3df225fe87ec8a20405df782de771f6e50f8d60506f3752`.
There is no new independent worker/result review. Completion inspection used
JSON, source, readable output, serialization metadata and opaque-byte hashing;
no scientific module was imported and no payload was restored outside the job.

The speed match here uses the actual selected transverse modal speed. It is
distinct from setting the register's bare speed ratio to one or carrying out a
primitive one-medium calibration. No draining background is present, and the
small-flow approximation is nonuniform at grazing for fixed nonzero flow.

The delivered [upstream joint assessment](S11c_upstream_repair_review_joint_disposition.md)
still leaves direct height–slope/nonuniform composition unresolved. The uniform
result retains that conditional applicability; it cannot clear the defect
operator. The parallel preparation (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11c_upstream_parallel_next_work.md`)
identifies a small same-profile direct-term diagnostic and affected downstream
objects. That diagnostic is not launched by this completion. Stop here before
any defect sweep, producer regeneration, new validator or loss calculation.
Full Green/FORM/A11/A12, physical calibration and draining-model leakage remain
open. Scratch stays ignored; old failures, review literals, scientific inputs,
Lean and the protected builder remain intact.
