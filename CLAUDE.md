# How we work

Seventeen rules. Every one exists because ignoring it cost a session.
Full process: `docs/development_pipeline.md`. What we're building: `docs/development_plan.md`.
Where we are: `STATUS.md`.

## The method — the CAS answers, not me

1. **Two engines exist so they can disagree.** Independent construction, not hidden answers. The
   disagreement is the measurement.
2. **A script prints computed objects. It never states conclusions.** Emit both operands and the residual,
   then guard — a residual asserted zero always prints `0` and carries no information. Interpretation
   belongs to the step record.
   **This binds me too, and that is the half I keep exempting myself from.** Every review-leg prompt I
   write says *"a prose derivation is worth nothing; show the script and its literal stdout, or the claim
   is discarded."* The same standard applies to anything I write. **A claim about an artifact carries the
   command that produced it**, and the commands and their literal output go in a `_measurements/` file
   beside the document — the same thing I make a leg do when I demand its script paths.
   Measured 2026-08-12: four export-chain designs and eight legs died on a question one `len(LEDGER)`
   answered; *"the imported action is in the LEDGER"* (there is no action row), *"no change to the publish
   guard"* (false under execution) and a class citation pointing at two sections that declare no classes
   each cost a full round. **All four were typed conclusions with no computation behind them — the exact
   defect this rebuild exists to remove, relocated into the orchestrator.**
3. **Name the object. Do not specify the recipe.** If a review is arguing about *how* to derive something —
   is this quotient well-defined, is this weight unique — the question was manufactured by specifying a
   derivation path. Ask for the object; let the engine hand over what it built.
4. **If I am deciding in prose what the engines should compute, I have inverted the method.** The tell:
   many turns reasoning toward an answer a script would settle in one. The cause, measured: the instrument
   was broken, so no measurement was available. Fix the instrument; don't reason around it.
5. **The spec says what to compute — never what anything equals**, is expected, or was measured. Withhold
   exactly one thing: an acceptance criterion referencing an expected value. A builder iterating to exit 0
   converges on any target it can see. Supply everything else, as equations.
6. **A disagreement is a finding.** Don't try to make divergence impossible with more careful prose; that
   defeats the reason there are two engines.

## The gates

7. **Whatever writes does not review.** Two legs on anything physics-bearing — and a spec both engines read
   is physics-bearing, because an error there makes both engines agree on the same wrong thing.
   Orchestrator-written → Codex + Grok. Codex-written → fresh Claude agent + Grok.
   **TRIGGER — no builder launches until its decision list has had two legs.** The list is
   orchestrator-written and is the one artifact the *builder* trusts: everything downstream is checked
   twice, the list is checked zero times. One pass, then fold and go — never iterated to green.
   Measured 2026-08-09: six spec defects, each costing a build round *plus* two legs, when two legs before
   the build would have caught them — three "level-above" misses, one exception named instead of the
   property (which bred a regression), measured counts stated four lines above "do not state the counts",
   and **an acceptance test that would have passed with the defect still in place.**
8. **Launch legs on sight, before I look at the result.** A self-check discharges the felt need for an
   independent one, and it is most convincing when it finds things.
9. **No commit before both legs report.** The commit is the last step. Reviewing the directive does not pay
   the tax for the build. My own verification is not a leg.
10. **Stop when nothing outstanding changes what is computed or what may be claimed** — not when both legs
    are green. "A leg that finds nothing is weak evidence" is my prior; put it in a leg's prompt and it
    becomes a quota.

## The traps

11. **Correctness is king; cost is never a reason** to drop a control, narrow a check, or skip a leg.
    Scaling work down is the user's call.
12. **A prohibition is not a control.** Blindness is enforced by *absence* — by bounding what the builder is
    handed — never by a sentence forbidding a read. A do-not-read list is a denylist, and a denylist means
    the architecture is wrong. The measured failure is absence of computation, not anchoring; quarantine is
    cut and rule 2 replaced it. **This applies to these rules as well: every one is prose I drift from
    under load, so the ones that hold are the ones that leave an artifact whose absence you can see.**
13. **A finding is not a mandate — verify it myself.** Legs have been wrong in both directions.
14. **Ablate to test; don't read.** A form control tests physics; a coefficient control tests arithmetic.
    Demand a script and its literal stdout — a prose re-derivation is the same defect relocated into the
    review. **Controls written against a script that does not exist cannot be ablated**, so they get
    reviewed by reading, which is the weakest instrument we have: put them in the build review, not in a
    document about a build.
15. **If successive revisions keep breeding defects in the material just changed, change the author.** Don't
    fold a fourth time.
16. **Prior art is an oracle, never a premise.** Check our computed result against it; never assume its
    result for our object. Its conditions may not be ours.
17. **A freeze is a red flag; a required freeze is the finding.** The one move that has cost this program
    a round, nearly every time, is treating a quantity that VARIES — a background field, a field's argument
    dependence, a jet/derivative order, a rate, a shape — as if it were constant. It manufactures a
    wrong-but-consistent answer both engines can share, so the comparator reads agreement and the defect
    hides (26=26 was a coincidence of two frozen mechanisms). Keep every varying quantity LIVE and
    differentiate it; when a method seems to require holding one fixed to proceed, that requirement is the
    measurement, not a step. Caught only by a ground-truth anchor or a variable-coefficient/form ablation.

## Repository infrastructure — the `.out` transcripts live in git-annex + GIN (set up 2026-09-01)

The v3 CAS audit transcripts (`research/pde_ledger_v3/scripts/out/*.out` and
`research/pde_ledger_v3/mathematica/out/*.out`, ~370 MB) are **git-annex content backed by GIN**
(`gin.g-node.org/trevnorris/toy_physics`, public), NOT plain git blobs — one exceeded GitHub's 100 MB/file cap.
The policy is the root `.gitattributes` (last-match-wins): everything is `annex.largefiles=nothing` (plain git)
**except** those two `out/*.out` globs, which are `anything` (annex). GitHub is `annex-ignore` by design (git +
tiny pointers only); GIN holds the bytes.

- **Generating/updating a v3 `.out`**: just `datalad save -m "…"` — the policy annexes `out/*.out` automatically
  and keeps everything else in git. Then publish with **both** `datalad push --to gin` (content → GIN) **and**
  `git push origin <branch>` (git + pointers → GitHub). ⛔ Never `git add -f` a big `.out`. ⛔ Never annex an
  `*_exports.py` — they are hash-chained plain-git inputs the next step imports.
- **After any fresh checkout/clone the `.out` are annex symlinks with NO content** until `datalad get <path>`.
  `grep`-by-line still resolves once content is present (the symlink is followed), so a script/directive/leg
  that reads a committed `.out` must `datalad get` it first — otherwise it reads a dangling link.
- The GIN token is stored in the datalad keyring under credential name **`gin`**; use `--credential gin` with any
  `create-sibling-gin`/publish. Full record, exact commands, and open follow-ups: the
  `project-datalad-gin-out-storage` auto-memory.
