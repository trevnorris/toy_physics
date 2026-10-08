# O2 comparator build r4: review dispositions (orchestrator)

**Artifact:** the r4 comparator, its tests and its ablation harness. Hashes are in the baseline
`_scratch/s9b_build/o2_comparator_build_review_baseline_r4.sha256`:

| File | sha (first 8) |
|---|---|
| comparator | `70c7e7c9…` |
| tests | `3c2ca9d0…` |
| harness | `f0961070…` |

The fresh Codex author (session `01a118a9…`, its second revision) produced it to amended directive `44f78617` from
brief `_scratch/s9b_build/o2_comparator_repair4_prompt.md`. Builder run:

- 70 tests OK;
- 33 harness ablations, each failing its designated test;
- the full guarded run exited 0 in 171 s, with sources last edited before every run.

The artifact is Codex-written, so it gets fresh Claude + Grok (G1), reviewed until clear (G4).

**Legs (identical prompt `_scratch/s9b_build/o2_comparator_build_review_prompt_r4.md`; review round 5).** Both
reported before adjudication.

- **Grok:** CLEAR, no findings (`_scratch/s9b_build/o2_comparator_build_review_r4_grok.txt`).
- **Fresh Claude (opus):** NEEDS REVISION, 3 findings (`_scratch/s9b_build/o2_comparator_build_review_r4_claude.md`).
  Its evidence files are copied from `/tmp/o2r4_fc/` to `_scratch/s9b_build/o2_comparator_build_review_r4_claude_evidence/`.

**Agreed sound by both:**
- coverage is complete, with 160 joined, 204 unjoined and 0 unaccounted;
- the tables are injective and the joins checked are correct;
- the readers are lossless;
- every algebraic residual recomputes to 0 from the raw transcripts by independent scripts;
- unformable residuals print as `not_formed` with distinct reasons;
- the one-sided corruptions are caught;
- there is no verdict token and exit is 0.

Each verification below is a mechanical lookup: sed, grep and sha256sum on the artifact, plus verbatim retrieval of
the leg's filed stdout. Commands and literal output are in `O2_comparator_build_r4_review_disposition_lookups.md`.

| # | Finding (leg) | Disposition | What must be true after repair |
|---|---|---|---|
| W1 | The named-OPEN inventory counts WL provenance strings and held-operator symbols as operands, and the computed limit miscounts those occurrences. (Claude 1) | **ACCEPT.** Lookups W1. Effect on output: the production output prints `"wl::6"` and `"wl::3,8"` on 6 lines each, and `"wl::D"` and `"wl::Map"` on 20 lines each. These come from WL `operand[UnspecifiedRotationalNormalRates,"6","3,8"]` (WL L266) and from `Inactive[D]`/`Inactive[Map]` heads. The comparator's own `leaf_policy` (L995) puts provenance and syntax outside the inventory. The printed named difference therefore reports WL-only "operands" that are not operands. This is a false difference in the output, not a false agreement. It sits in material changed in r4 (the labels repair and the computed limit). | Only operands and labels that an OPEN occurrence names enter the named inventory. Provenance and held-operator spelling stay outside it, wherever the action sits in the term. Its own role head is consumed, not counted as outside, at any depth. Every printed limit count follows the same leaf policy. |
| W2 | Three balance-entry orientation branches have no control. (Claude 2) | **ACCEPT.** Lookups W2. Balance entries take their orientation from `term_orientation` (L1142–1161, used at L1187/1189). The harness's orientation knife edits `signed_actions` (harness L85–86), and no harness line names `term_orientation`. The leg's copies removed the Lambda-body, inactive-aggregate and held-body branches, and all 70 tests still passed (F1–F3 `OK`). On the production streams the same copies print byte-identical output (`differing lines 0` ×3), so today's output is faithful. Even so, the amended item-6 *Balances* bullets ("held or inactive aggregates"; "SymPy `Derivative` and `Lambda` bodies") have no control that fails. That meets the review filter. | Each item-6 *Balances* orientation bullet has a control through the balance-entry path, and a harness ablation of the branch it covers fails a test. |
| W3 | The printed limit does not name the coordinate of a held derivative, or the variable of a variation, that wraps an OPEN action, nor a held aggregate's set or binder. Changing it moves no comparison product. (Claude 3) | **ACCEPT, as a statement of the limit.** Lookups W3. `not_compared` (L992–993) names argument placement, multiplicity and "argument content outside these inventories". A wrapper's variable is not an argument of the action. The leg's one-sided corruption D (PY `HOLD_INPLANE`, the derivative coordinate of `OPEN_MomentumFlux_0_0` changed from `x1` to `x2`) moves only operand text and diagnostic hashes (`differing lines 5`, all in operands, hashes and entry operand text). A reader could take the transport term's empty entry difference to cover its divergence direction. The leg proposes no widening, in line with the user's three-inventory scope. Whether to compare wrapper variables is a scope question for the user. | The printed limit names every kind of content in a paired OPEN occurrence or balance entry that is not compared, including wrapper derivative or variation variables and held-aggregate sets or binders. An empty difference cannot be read as covering them. |

**Note N1 (not a comparator defect; carried to the record).** PY `OPEN_UnfixedNormalGeneralizedRates` (audit
L277, normal only) is paired with WL `UnspecifiedRotationalNormalRates` (WL L266, rotational and normal), citing §6
generalized rates (comparator L452, L557). PY's separate `RotationalGeneralizedWork` (audit L279) is unjoined with
the reason that WL keeps rotational work inside its joint material-work action (accounting). That accounting is
true, and the difference in grain is printed. The record (sub-step 7) states this pairing's scope difference, and
does not read the empty role difference on that row as identity of scope.

**Grok's clear.** Grok missed W1–W3. Its checks of coverage, residuals, the stated limit for `NativeFaceMeasure`,
and controls are consistent with the Claude leg, which found these three by targeted form ablation, a one-sided
corruption, and inventory recomputation.

**Pattern (G4).** All three findings sit in material r4 changed:
- W1 in the labels and computed-limit repair;
- W2 in the item-6 control repair;
- W3 in the printed-limit repair.

No residual moved. Across five rounds and two authors, every leg has recomputed the closed results as exact
zeros; this round the Claude leg recomputed 137 algebraic rows and found 0 nonzero. The open findings concern what
the OPEN-content difference prints and claims, plus one contract-control gap. The memory rule says that more than two
rounds means stop and put the premise to the user. The user owns that decision.

**Next.** Preserve r4 as the reviewed baseline (not accepted). The repair route and the W3 scope question go to
the user.
