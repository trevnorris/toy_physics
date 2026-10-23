# S16 comparator build review, round 0: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S16_comparator_review_r0_lookups.md`, beside this file
(generator `_scratch/s16/cmp_review/gen_r0_lookups.sh`). "P-n" cites a section of that file.

## The round

- **Artifact.** The comparator `research/pde_ledger_v3/scripts/S16_cross_engine_comparator.py`, its tests, the run
  record, the builder's report, and the production output `_scratch/s16/cmp_work/full_05.jsonl` with its summary and
  the last test output. Codex (`gpt-6-astra`, high) wrote them, to the comparator decision list at `5994e9ac`.
  - The build was interrupted twice: Codex's model API returned "at capacity", and the machine rebooted during the
    `full_04` run. Codex resumed its own session both times. The production output is the `full_05` run.
  - The builder's report is its final message, copied verbatim from its transcript.
- **Prompt and legs.** `_scratch/s16/cmp_review/s16_cmp_build_review_prompt_r0.md` (P1). Codex wrote the artifact, so
  the legs are a fresh Claude agent (Opus) and Grok (`grok-4.7`, high) (G1). They had the identical prompt and ran in
  parallel. Grok exited 0.
- **Baseline.** All seven baseline files matched after the round. The working tree was unchanged apart from the
  lookups file (P1).
- **Verdicts.** Claude: NEEDS REVISION, three findings. Grok: NEEDS REVISION, one finding, the same as Claude's first
  (P2).

## Findings

**K1. A volume integrand's partial time derivative is read as the derivative along the moving chart.** (Claude 1,
Grok 1.) **CONFIRMED.**
- The comparator substitutes the chart into a volume integrand with SymPy's `subs(…, simultaneous=True)`
  (`Names.integral`, l.763; P3).
- Both engines define the volume head with `x` a bound Cartesian argument (P3):
  - the Wolfram schema: `"V" -> "volume integral; x remains a bound Cartesian argument"`;
  - the SymPy engine: "Volume arguments retain bound x".
  
  So a derivative in `t` inside a volume integrand is the partial derivative at fixed `x`, evaluated at the chart
  point.
- That substitution turns `Derivative(F(x1,x2,x3,t), t)` into `Derivative(F(rho + Xc1(t), x2, x3, t), t)`, a
  derivative along the chart. A non-simultaneous substitution keeps `Subs(Derivative(F(x1,x2,x3,t), t), x1, …)`
  (Claude's probe; Grok's `REPRO` lines; P3).
- Claude counts 44 value bindings whose volume integrand carries such a derivative: on 28 both operands are read
  this way, on 16 only the Wolfram operand (P3).
- In `full_05`, those 16 bindings, 8 `localEvolution` pieces and 8 `partsSum` partitions of the bulk face terms of
  `Nrate` and `Drate0`–`Drate2`, are the only ones with numeric sample residuals other than 0. On binding 182 these
  occur in the generic, radius-velocity-zero and zero-count strata, and not where the centre velocity vanishes (P3).
  The two other bindings the scan lists, `domains:0` and `domains:1`, are condition comparisons that print `False`
  (P3).
- Claude's probe repair of l.763 alone moved all 44 bindings. All 13 test groups passed both before and after it
  (P3), so no control detects the misreading.
- The list: item 1, "A different object is still read as different. Examples: another derivative slot, order or
  point". Item 4: "Each engine's `integral.*` lines and source define its heads." Item 5: a live function's derivative
  taken at another point inside a volume integral is a different object, and evaluation sees it.
- **What must be true.**
  - In a volume integral, both operands evaluate the integrand at the chart point. A derivative in `t` there stays
    the partial derivative at fixed `x`, on both sides.
  - The 44 bindings are compared as those objects.
  - A control fails when a volume integrand's partial time derivative is read as the derivative along the chart.

**K2. The zero-count stratum sets `N` to zero but leaves the printed count nonzero.** (Claude 2.) **CONFIRMED.**
- The stratum changes only the test value of `N` (`profile()`, l.845; P4). `Σ`'s test function does not depend on
  the stratum.
- In `full_05`, at form `mixed_quadratic`, point 0 (P4):
  - `definitions:37`, `N(t)`, is `19/13` in the generic stratum and `0` in the zero-count stratum;
  - `definitions:38`, the printed count `∫Σ`, is the same nonzero value in both.
- Claude counts 7785 zero-count samples identical to their generic sample, and 18 that differ. The 18 are in
  `definitions:37` and `domains:1` (P4).
- The test control for this coverage (test l.376) uses a witness in the symbol `N(t)`, so it passes (P4).
- Grok found the stratum present, with normalized objects marked outside their domain (P4). That is consistent: the
  stratum exists, but the unnormalized objects in it are evaluated at a nonzero count.
- The list (P4): item 5, Strata: "a stratum where `N` vanishes at the point, with the fields and the geometry live.
  There the unnormalized objects are compared". Domains: "The unnormalized balances and moments keep the full domain,
  including where `N` is zero".
- **What may not be claimed now.** That any unnormalized object was compared where the count vanishes.
- **What must be true.**
  - In that stratum the live fields make the printed count vanish at the point, with `Σ` still live. `N` is zero
    there too.
  - The unnormalized objects are evaluated there.
  - The coverage control for this stratum fails when the stratum's fields leave the printed count nonzero.

**K3. Exact residuals that contain the comparator's own unevaluated substitution wrappers are printed unresolved.**
(Claude 3.) **CONFIRMED.**
- `exact_difference` (l.956–970) expands and groups the difference. It does not evaluate the comparator's canonical
  `Subs(…, (_slot0, …, _slot5), (Xc1(t), …, t))` wrappers (P5).
- In `full_05`, 655 subtraction residuals print as `0` and 230 do not. Of the 230, 114 contain the literal
  `Subs(0,`. An example is `derivative_piece_centerMotion:183`, whose residual is a sum of integrals of
  `Derivative(Xc1(t), t)*Subs(0, …)` (P5).
- Claude's own canonicalization resolves 228 of these (P5).
- The list, item 5: "Exact first. Use exact symbolic subtraction where feasible, and print the exact residual as
  computed. Where it is feasible, evaluation at points is added beside it, never in its place."
- **What may not be claimed now.** An exact result for these bindings. Only the sampled comparison stands for them.
- **What must be true.**
  - Where the comparator's own wrappers can be evaluated, the exact residual is also printed in its resolved form,
    beside the as-computed one, never in its place.
  - The resolution keeps every derivative's slot, order and point (item 1), so K1's distinction survives it.

**K4. The name table cites the wrong SymPy lines for `N`, `D` and `P`.** (Claude, under "Name table".) **CONFIRMED;
minor.**
- The table cites SymPy lines 102, 103 and 104. The engine introduces `N`, `D` and `P` at lines 101, 102 and 103
  (P6). The Wolfram citations (120–122) are the definition lines.
- The list, item 3: "Each entry cites the object both names denote and the line in each engine that introduces it."
- It changes nothing computed. It is corrected with the repair, because the printed table is the grounding a record
  cites. On its own it would not hold acceptance.
- **What must be true.** Each cited line introduces the object.

## Not findings

- **The tests stay green when the production tables are ablated** (P7). Claude: pairing surfaces by index, transposing
  `Π`, merging pieces into totals, dropping the zero-count stratum, pairing components by position. Grok: transposing
  `Π`, retargeting piece bindings, removing the zero-count stratum.
  - Under each of these, the production output moves (both legs, P7).
  - Item 7's controls test that a mutation made in a fixture changes the comparison. They do not pin the production
    tables, which both legs checked against their own mappings.
  - K2's control is the exception, and K2 covers it.
- **One-sided corruption** (P8). Claude's corrupted SymPy line moved only its own binding's comparison and leaf
  counts. Grok's moved only its own tag.
- **Reproducibility** (P8). Each leg's guarded rerun of the final code reproduced `full_05`, apart from runtime and
  memory fields.
- Both legs found the reading, the binding table, the name table's injectivity and index order, the integral heads,
  the conditions, the geometry, the dimensions and labels, and the absence of a target sound. The run record matched
  the guard records.

## Outcome

Four findings confirmed: K1–K3 change what is computed or what may be claimed, and K4 is a citation. The reviewed
version is preserved with these dispositions. Codex repairs it to what must be true. The repair gets a round 1 review
by a fresh Claude agent and Grok.
