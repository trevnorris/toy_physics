# S15 comparator build review, round 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S15_comparator_build_review_r1_lookups.md`, beside this file
(generator `_scratch/s15/gen/s15_cmp_r1_lookups.sh`). "M-n" cites a section of that file. Every computed result quoted
here is a leg's saved stdout; the legs' scripts were not re-run.

## The round

- **Artifact.** Repair 1 of the round-0 findings, uncommitted on `0ad8773e`:
  - the tests `scripts/test_S15_cross_engine_comparator.py`;
  - the run record `directives/_measurements/S15_comparator_run_record.md`;
  - the builder's final report `directives/_measurements/S15_comparator_builder_report.md`, copied verbatim.

  The comparator and the guard are byte-identical to `0ad8773e`. The repair changed only the three files above (M2).
  The production output `_scratch/s15/cmp_work/full_07.jsonl` and its audits are unchanged. Everything matched the
  baseline after the round. The working tree was unchanged apart from this file and the lookups file, and no guard
  unit was still running (M0).
- **Gate.** G1: a Codex-written script, so its legs are fresh Claude + Grok, reviewed until clear.
- **Prompt and legs.** `_scratch/s15/cmp_review/s15_cmp_build_review_prompt_r1.md`, sha256 `546ba3e7…` (M0). Both legs
  received it identically, in parallel:
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s15-cmp-r1-grok/REVIEW.md`;
  - fresh Claude (Opus agent): `/var/projects/_review_scratch/s15-cmp-r1-claude/REVIEW.md`.
- **Verdicts.** Grok: CLEAR. Claude: CLEAR. Neither leg filed a finding (M1).

## The round-0 findings

**C1, the global disable switch: repaired.**
- The tests no longer set `disable_comparison`. It is named only in a comment at test line 96. The run record states
  "No global `disable_comparison` option is used." (M3).
- *Claude:* in a copy of the comparator, the leg altered one step at a time and left the rest of the comparison on.
  Each of eight steps then failed its own group with that step's assertion. The steps were the exact residual, the
  samples, the text difference, the domain, the dimensions, formation, leaf credit and the record ordinal (M3).
  Restoring round 0's global switch fails with `residual fault disabled evaluation`: the harness refuses a fault that
  is not targeted (M3).
- *Grok:* restoring the global switch gives `{'not_compared': 'comparison disabled'}` and an inactive-baseline
  failure, not a targeted one (M3).

**C2, the sparse-moment oracle: repaired.**
- The oracle takes its moments from its own `fourier_moment` (test lines 915, 928–929). The three remaining uses of
  `c.angular_moment` are the halving ablation, which replaces it and then restores it (lines 946–952) (M4).
- *Claude:* halving the full-range moments, or the half-range ones, fails `sparse_integral_oracle` with `sparse moments
  differ from independent Fourier expansion`. Under the same halving with round 0's oracle restored, the group passes
  (`ERRORS []`) (M4).
- *Grok:* on 98 exponent pairs the Fourier moment equals `angular_moment` (`mismatches 0`). Under halving, the
  production value is `2*pi/3` and the Fourier reference `4*pi/3`. An oracle built from the halved `angular_moment`
  matches production (`form_ablate_02.stdout`, lines 3–5) (M4).
- One line in Grok's `c1c2.stdout`, `round0_shared_moment_would_match False`, compares a product at exponents
  (3,0)(0,2) with a reference at (1,0)(0,2) (probe lines 132, 143, 148; M4). The exponents differ, so that line says
  nothing about round 0's oracle. Grok's report rests on `form_ablate_02.stdout`, which agrees with Claude's restored-
  oracle run. This is not a finding.

## The zero-operand question (round 0, binding `region_momentum_groups:1373`)

**Settled by computation in both legs: the printed `0` is the integral's value.**
- The comparator prints SymPy `0`, Wolfram `0` and residual `0` at all 45 samples (M5).
- *Claude:* the leg read the SymPy operand with its own reader and the Wolfram operand with Wolfram's parser. It
  evaluated them by quadrature at the comparator's 45 samples and printed each value over ∫|integrand|. For 1373 that ratio is
  at most `3.25e-16` (SymPy) and `3.76e-17` (Wolfram). Across the six integrals that print zero it is at most
  `7.36e-16`. The nonzero drain-1 siblings give `1.00e+00` and `8.07e-01`, a check that the quadrature sees a value
  when one is there (M5).
- *Grok:* on both sides the integrand reduces to a factor independent of φ times `sin(φ)`. The φ integral is `0`, the
  absolute integral is nonzero, and `value_over_scale 0` (M5).
- Round 0's unscaled quadrature values (`1.9e-07`, `1.0e-04`) are therefore not evidence against the printed zero.

## The output is the final code's output

- The comparator's bytes are the same at the round-0 and round-1 baselines (M2), so `full_07` is the output of the
  code accepted here.
- *Claude:* a rerun of 1470 bindings matches `full_07` (`identical=1470 differing=0`). Workers 1, 2 (reversed) and
  3 (forward) give identical records on 18 bindings with group, nested and linearity integrals (M6).
- *Grok:* six rerun bindings match `full_07`, and a worker split gives equal records (M6). Binding ends are
  injective at both engines (1616/1616), and name rows too (45/45) (M6).

## Examined and not filed

- **Table-content ablations move production and no test** (Claude f01, f02, f03 and f06; M7). Grok's ablations of the
  same kinds moved production (its report).
  - This is the round-0 disposition, unchanged. The repoint controls are built from the loaded tables, so they test
    detection, not the tables' content.
  - The content is checked by both legs' compare-first. Claude's pairing agrees with the comparator object by
    object. In Grok's, twelve selected pairs have zero residual and distinct objects stay distinct (their reports).
  - The strata are printed: 45 samples on each of the 712 sampled bindings.
- **Grok's density collision.**
  - Grok replaced `profile_cached` so that `rho_2` receives `rho_1`'s test function. The pressure-piece sample moved
    and its exact residual stayed `0` (M8). Grok did not run the tests under that replacement.
  - The control `test_functions` asserts distinctness for every declared profile, `rho_1` and `rho_2` included,
    through `profile`, which calls `profile_cached` (M8). That is a reading of the source, not a run.
  - In the unaltered code, the two profiles differ (`rho_profiles_differ True`, M8). Not a finding.
- **The switch itself.** `disable_comparison` remains in the comparator at line 1180 (M3). Nothing in the controls
  relies on it any longer.
- **The legs' failed jobs** were their own scripts (M9):
  - Claude's `item10_py_01` failed on a ZeroDivisionError in its own script; `item10_py_03` is the run it reports.
  - Grok's `form_ablate_02` exited 1 after printing its lines, in its own `clear()`. Its `own_compare_05` was
    superseded by `_06`.
  - Every one of these jobs ran with verified limits.
- **Coverage of the exclusive-mode control.** Neither leg ran `exclusive_guard_controls`, which needs the exclusive
  guard mode. Grok checked worker splitting under an affinity override instead. The builder's guarded test run
  `tests_20261025_15` completed all 27 groups (run record).

## Decision

**Accepted at this version.** Nothing outstanding changes what the comparator computes or what may be claimed from
`full_07` (G4).
- Both legs are CLEAR. C1 and C2 are repaired, as both legs' runs show, and the zero-operand question is settled by
  computation.
- The commit includes the repaired tests, the run record, the builder report, this file and its lookups.
- The production output is `full_07.jsonl`, the run of the unchanged comparator that rounds 0 and 1 reviewed. It is
  filed as `scripts/out/S15_cross_engine_comparator.out` (git-annex, GIN), with `S15_comparator_production_run.md`
  (identity with the reviewed run, inputs, guard outcomes), as S16 did at `f08db34a`.
- What the output may be taken to show is the S15 record's question. The record's decision list carries three
  notes:
  - `c_comp` appears only as `p_br′`;
  - `balance_equation` is an unreduced `Eq`, not a statement that the balance holds;
  - prose conditions compared as text are textual evidence only.
