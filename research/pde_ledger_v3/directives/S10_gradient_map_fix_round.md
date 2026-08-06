# Close two gate holes, correct one comment, and put the cited evidence where it survives

**Do not commit.** Four items. Both legs cleared the nine declarations themselves — they derived the map
independently from the two engines and agreed — so ⛔ **do not change any declaration.**

## The measurement all four items rest on

Over all `9! = 362 880` bijections of the nine gradient symbols, **288 leave every emitted `§Q7` density
invariant**. Both legs measured this independently and agree on the number and on its structure:

```
S3 on the three diagonal symbols  x  S3 on the three off-diagonal PAIRS  x  2^3 within-pair swaps
       6                          x          6                          x    8      =  288
```

⇒ the payloads pin only **which symbols are diagonal** and **how the off-diagonal ones pair**. Nothing
else. A declaration wrong by any of those 288 produces the identical verdict on every declared row —
under the transpose, the harness's **entire counter line** is byte-identical.

**Measured consequence (a leg, eight mutations):** two wrong declarations pass **both** existing gates.
`SWAP_g22_g33` — a wrong map — passes the unit suite and `ACCEPTANCE 19`. `DROP_g33_ONLY` — a silently
removed declaration — passes both, while four `§Q7` rows quietly fall out of verdict and item 19 prints
`AS_DECLARED=[AGREE=8,NAMING_MISMATCH=4]` and still reports `PASS`.

---

## D1 — `reduction/test_engine_output_checks.py`: pin the diagonal entries

`test_s10_declares_only_verified_dimension_scale_and_gradient_spellings` pins the applied-rename set on
`main_d3_q7_stiffness`. That payload is the curl density, in which `g11`, `g22`, `g33` **cancel** — so the
three diagonal declarations are pinned by nothing.

Add assertions for the one declared row whose payload carries **all nine** symbols:
`xform_fullgrad_d3_q7_stiffness`. Pin its status, its **complete** `naming_applied` set, and that its
`undeclared_spellings` is empty — the rename set written out **by hand** from

```
WL  g{r}x{c}   ==   PY  g{r}{c}        for r, c in 1..3        both standing for  d(u_c) / d(x_r)
```

⛔ **not** pasted from a harness run. ⭐ Writing it by hand is what lets this assertion catch a mutation
**no payload can detect** — a leg confirmed it catches the full transpose, which `ACCEPTANCE 19` cannot.

⛔ Do not remove or weaken any existing assertion in that test.

## D2 — `reduction/harness_ablation.py`: assert that the declaration is COMPLETE

`ACCEPTANCE 19` compares `AS_DECLARED` against `DERANGED` and `ABSENT` and never pins `AS_DECLARED`
itself, which is why a dropped entry passes.

Add one further requirement: **under `AS_DECLARED`, no declared `§Q7` row is `NAMING_MISMATCH`.**

⭐ That is a statement about the declaration being **applied to every symbol that occurs**, ⛔ not about
which map is right — so it does not reintroduce an assertion the payloads cannot support, and it leaves
`TRANSPOSE` unasserted as before.

⚠ **And print which rows moved**, per variant, not only the status tally. A reader currently cannot see
that under `DERANGED` the `xform_fullgrad` stiffness row stays `AGREE` — its `Σ g²` is invariant under a
within-row column cycle. Other items in the module already print names (`ACCEPTANCE 4 partial_rows`,
`5 new_named_gaps`, `6 main_missing`); follow that convention.

## D3 — `reduction/checks_S10.yaml`: correct the comment, ⛔ not the declarations

Three corrections to the comment above the nine entries. **Change nothing else in the file.**

1. It says the payloads cannot separate the declared map **from its transpose**. True, and it understates
   the result by a factor of 144. Replace with the measurement above: the size of the group, its
   structure, and that the payloads fix only the diagonal/off-diagonal-pair partition.
2. It cites *"measured on all 18 emitted D=3 pairs"*. **Four of those 18 are `0` on both sides** and carry
   no information about any map; the informative base is **14**. Say 14.
3. Its closing sentence says a per-row search could absorb a disagreement. A leg tested that claim by
   applying a deranged relabelling to one row: the fixed map does ⛔ **not** turn that row into
   `DISAGREE` — it stays `NAMING_MISMATCH`, i.e. a **coverage gap**, and becomes visible because it is
   then the **only** non-`AGREE` row among the twelve. Correct the sentence to say that.

⭐ Keep the provenance block with its engine line numbers exactly as it is.

## D4 — the cited evidence must survive the commit

`.gitignore:96` ignores `research/**/_scratch/`, so **every script the config comment and the step record
cite will never be committed.** A committed artifact that cites evidence which evaporates is worse than
one that cites nothing.

Move these nine measurement scripts from `research/pde_ledger_v3/_scratch/` to a **tracked** directory,
`research/pde_ledger_v3/reduction/measurements/`:

```
q7_index_map_discrimination.py     q7_payload_invariance_group.py    q7_map_diagonal_gap.py
q7_map_declaration_ablation.py     declaration_load_ablation.py      stiffness_coefficient_recovery.py
q8_stratum_manual_comparison.py    bare_word_payload_scan.py         q3_sign_adjudication.py
route_row_information_content.py
```

(That is ten; move all ten.) Each locates the tree by `Path(__file__).resolve().parents[N]` — **fix `N`
for the new depth** and change nothing else about what any of them computes. Add a short `README.md`
listing what each one measures, in one line each.

Then update the citation in `checks_S10.yaml`'s comment to the new path.

⛔ **Do not rewrite these scripts.** Both legs ran them and reported them clean; a rewrite would discard
that review.

---

## Out of scope — report, do not fix

Either engine · either committed `.out` · `steps/` · `paper/` · `directives/` · `REBUILD_HANDOFF.md`
(known stale in two places; the orchestrator owns it) · any other test or acceptance item.

## Acceptance — run these and paste literal stdout

1. The full unit suite from `reduction/`. A one-line cause for any failure.
2. **Show D1's new assertion can fail:** on a copy under `/tmp/`, swap the `g22` and `g33` declarations,
   run that one test, paste the failure. Confirm the committed configuration passes.
3. **Show D2's new requirement can fail:** on a copy under `/tmp/`, delete the `g33` declaration, run the
   battery, paste the failure. Confirm the committed configuration passes.
4. `harness_ablation.py` in full, literal.
5. **For every moved script: its literal stdout after the move, and confirmation that it is byte-identical
   to its stdout before the move.** Show the comparison command and its result.
6. The harness on the S10 config — the literal `CROSS_ENGINE:` line only.
7. `git status --short` and `git diff --stat`.

## Constraints

- No script over 600 seconds; a timeout means reformulate, never raise the limit.
  ⚠ `declaration_load_ablation.py` on the S10 config takes about 13 minutes — run it in the background
  and wait for it, do not wrap it in a shorter timeout and do not skip it.
- **Do not launch Mathematica or `wolframscript`.**
- No `git add`, no `git commit`, no other git write.

## Report back — under 20 lines

1. `D1`–`D4`: done / not, with line numbers and the new paths.
2. The acceptance output.
3. Anything you measured that contradicts this directive.
