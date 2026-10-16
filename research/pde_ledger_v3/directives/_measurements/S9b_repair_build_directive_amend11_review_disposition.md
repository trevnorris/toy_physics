# S9b repair build directive, amendment 11: review dispositions (orchestrator)

Amendment 11 answers the build review r7 (`S9b_repair_build_r7_review_disposition.md`, rows R1–R4) and carries the
user's choice of 2026-10-16 about what the harness must show. It changes what two ablation harnesses print, so it is
reviewed as physics until clear. The orchestrator wrote it, so the legs are Codex and Grok (`CLAUDE.md` G1). Lookups:
`S9b_repair_build_directive_amend11_review_disposition_lookups.md` (mechanical; generator
`_scratch/s9b_build/gen/s9b_amend11_r0_lookups.sh`).

**The user's choice (2026-10-16).** I put the build review's eighth round to the user: what must round 8 clear?
The options were all of item 14 on every output, or a narrowed bar. The user chose the narrowed bar:
- no false zero anywhere;
- at each knife's own site, a shown change in both engines;
- away from the sites, undecided is allowed.

The amendment states this as item 14's "What the harness must show".

## Round 0

**Artifact.** The working-tree directive at sha256 `0329a03d…`. I folded round 1 into the working tree before saving
round 0's bytes. They were rebuilt exactly: the Codex leg's captured `git diff 4b4db191` (`amendment.stdout`),
applied to the directive at `4b4db191`. The rebuilt copy's hash equals the round-0 baseline entry (lookups). It is
preserved in `_scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r0.md`. Prompt:
`_scratch/s9b_build/s9b_repair_build_directive_amend11_review_prompt_r0.md` (`161a6ea8…`).

**Legs.**
- **Codex (gpt-6.1-sol, xhigh): CLEAR, no findings.** Final message:
  `_scratch/s9b_build/s9b_amend11_review_r0_codex_final.txt`. Full report and scripts:
  `_scratch/s9b_build/s9b_amend11_review_r0_codex_evidence/` (`report.md`).
- **Grok: NEEDS REVISION, three findings.** Report: `_scratch/s9b_build/s9b_amend11_review_r0_grok.txt`; the report
  starts at "**Verdict**". Script and stdout: `_scratch/s9b_build/s9b_amend11_review_r0_grok_evidence/`.

Both legs found:
- no leakage;
- the engines unchanged;
- the text consistent with "It prints no verdict".

| # | Finding (leg) | Disposition | Round-1 text |
|---|---|---|---|
| G1 | A part the engine has shown changed may print as a change or as undecided. Rule 2 prints a change "where the engine shows the part changed". The new bullet lets "any other part that differs" print undecided after a witness search. The header dropped the search condition ("Any other part may print undecided"). Two builders can read this two ways, which is R3's case. (Grok 1) | **ACCEPT IN PART.** The header's dropped condition is a defect: it reads as a second rule that needs no search. The text must also say which rule governs a part away from the sites. Grok's fix makes rule 2 bind everywhere. I reject that part. The user's choice allows undecided away from the sites, and making rule 2 bind there would put R3's off-site leaves back into round 8's bar. | The header points at item 14 and restates nothing. The bullet says that a part that is not a site part of the knife, and that differs, may print undecided once the harness has looked for a witness, even where the engine has shown it changed. At site parts, the bullet requires rule 2. |
| G2 | "Where both values are defined", with the point free to ignore the attached domains, can print a change for a value equal on its domain. Grok's script: `sqrt(x**2)` against `x` agree at `x=1` and differ at `x=-1`, where both are real; the attached domain `x>0` was what made them equal. Item 9 says a response "has a value only where both squares are real and non-negative" and "A principal root of a negative square is not a value". The Wolfram engine carries that condition on the response's `ConditionalExpression` (L97–102), which the harness collects as an attached domain. (Grok 2) | **ACCEPT.** The defect is mine. Round 0's "A value is compared apart from its domain" was my reading of R2, and it lets a witness fall where an object has no value. That would make a change at a site print where none exists. R1 and R2 do not need it. R1's defect is the zero, and R2's is a gate that accepts only literal `True`. Both builders read "defined" to include the attached domains, and that reading is safe: a value equal on its common domain cannot then print as a change. So the witness stays inside both domains. What changes is how a condition is shown to hold. | For a value, a witness is a real point of its variables at which both values are defined and both attached domains hold. A condition holds at a point only where the harness shows that it holds there, whether it is stated as a formula or in words (such as that the required derivatives exist). The spec's order counting is not such a condition. If the two attached domains share no point, that shows nothing about the values. |
| G3 | "Parts computed from the mutated construction" is not fixed for K10 or K11. K10: `BRANCH_EXISTENCE` is emitted from `real_domain` before the gate is assigned (SymPy L719–720, L770; Wolfram L166). K11: the method's residual does not involve the mutation, tags the method does not reach never meet the condition, and nothing makes the method reach the mutated construction. Grok also cites `docs/development_pipeline.md` §4: zero difference under a FORM knife is a fail, and on a two-route object the other route must stay. (Grok 3) | **ACCEPT IN PART.** Name the set, including K10's and K11's. Two parts are rejected. On §4: the harness prints and the review judges whether a knife bites. A conditional requirement on the harness does not make zero difference a pass, and stating the expected bite in the builder's directive would hand the builder the outcome (M2). On two-route objects: the directive names no S9b knife as acting on one (lookups: no "route" in item 14's knives). | A knife's site parts are the parts of the emitted objects whose computation uses the mutated construction, including the object emitted at the construction site itself. They do not include the build-input digest, or a part computed before the mutation enters. For K10 they are the parts gated on the negated condition. For K11 they are the compact differences of the reached tags whose computation the mutation enters, not the method's residual, and the method reaches at least one such tag. |

**R2's stated outcome is revised.** The build review r7 disposition's R2 row says "The attached domains need not hold
there." G2 shows that outcome is wrong, so round 1 supersedes it. R2 is now met this way:
- the witness stays inside both attached domains;
- a condition holds where the harness shows it holds, not only where it evaluates to a literal `True`.

**Author.** G2 is the first defect in material this amendment changed, and it is mine. Under G4, a further defect in
material a fold just changed moves the text to a new author.

**Process.** Grok's guarded runs are under `_scratch/s11c/s9b-amend11-review-r0/defs-*`, and Codex's under `k11-*`,
`logic-wl-*` and `semantic-*`. All completed.

## Round 1

**Artifact.** The working-tree directive at sha256 `75f4d562…`, saved before launch in
`_scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r1.md`. Prompt:
`_scratch/s9b_build/s9b_repair_build_directive_amend11_review_prompt_r1.md` (`22bf4bfe…`). Lookups:
`S9b_repair_build_directive_amend11_review_r1_lookups.md` (generator `_scratch/s9b_build/gen/s9b_amend11_r1_lookups.sh`).

**Legs.**
- **Codex (gpt-6.1-sol, xhigh): CLEAR, no findings.** Final message:
  `_scratch/s9b_build/s9b_amend11_review_r1_codex_final.txt`. Report and evidence:
  `_scratch/s9b_build/s9b_amend11_review_r1_codex_evidence/` (`report.md`, `evidence.md`).
- **Grok: NEEDS REVISION, two findings.** Report: `_scratch/s9b_build/s9b_amend11_review_r1_grok.txt`; the report
  starts at "**Verdict**". Scripts and stdout: `_scratch/s9b_build/s9b_amend11_review_r1_grok_evidence/`.

Grok confirms round 0's G2 fix on its own cases. `sqrt(x**2)` against `x` on `x>0` no longer prints a change, and
`x` against `2x` on disjoint domains no longer prints zero.

| # | Finding (leg) | Disposition | What must be true |
|---|---|---|---|
| H1 | "A condition holds at a point only where the harness shows that it holds there" does not fix what showing is, so two builders print differently. Wolfram attaches structured domains, not formulas: the reduction's `regularityDomain` (L400–404) has a formula radius condition, three conditions stated in words, and the spec's order counting. The propagation gate is an `Inactive[ForAll]` (L262). `TrueQ` is false on a bare quantifier that `Resolve` reduces to `True` (`resolve_wl.stdout`). Grok's classifier gives `x` against `2x` under such a domain three readings: undecided under two of them, a change under the third (summaries: 3 CHANGED, 9 UNDECIDED). (Grok 1) | **ACCEPT.** Round 1's sentence, which I wrote, names an obligation without its test. A literal-`True` gate, which is R2, stays a legal reading. | Two builders apply the same test of whether an attached domain holds at a witness point. The test covers a formula, including a quantified or inactive one; a condition stated in words; and a structured domain that mixes them. The spec's order counting is never tested. A condition the harness cannot decide at a point leaves that point unusable. It never makes a part print zero or a change. |
| H2 | K10's site parts are two different sets. The general sentence takes the parts whose computation uses the mutated construction. The K10 sentence takes "the parts gated on the negated condition". In Wolfram, `branchGateLocal` (L166) feeds `LOCAL_BRANCH_TYPES` (L182, L184), the emitted propagation gate `BRANCH_EXISTENCE` (L275) and the local nonreciprocal construction's domain, none of which is passed through `gated`. Grok's script marks those three `DIFFER` between the two sentences, and agrees on the rest. SymPy `BRANCH_EXISTENCE` is emitted from `real_domain` (L720) before the gate is assigned (L770), so neither sentence includes it. (Grok 2) | **ACCEPT.** My round-1 K10 sentence narrowed the general rule. | One definition of a knife's site parts that gives the same set in both engines. For K10, that set follows from the engines' constructions after the negation, not from a narrower K10-only phrase. SymPy `BRANCH_EXISTENCE`, computed before the mutation, is outside it. |

**Author change (G4).** H1 and H2 are both in sentences I wrote in round 1, the fold of round 0. Round 0's G2 was in a
sentence I wrote too. Under the rule this disposition set in round 0, the text moves to a new author. Codex
(gpt-6.1-sol, xhigh) writes item 14's amended text for round 2. It works from:
- the outcomes required above and in round 0 (the "Round-1 text" column, read as what must be true);
- the build review r7's R1, R2 and R4, with R2's outcome as revised in round 0;
- the user's choice;
- both legs' reports from both rounds;
- the engines.

The author is given what must be true, not my wording. The orchestrator then copies the text in unchanged after
checking that the diff stays in item 14's amendment-11 text, the header and the Part 2 and 3 lines.

**Legs for round 2.** The round-2 text is Codex-written, so its legs are a fresh Claude agent and Grok (G1). My
round-0 and round-1 sentences that survive have already had Codex and Grok. Authorship is recorded as mixed:
orchestrator-drafted in rounds 0 and 1, Codex-written in round 2. Codex, the author, does not review.

**Process.** Both legs' guarded runs are under `_scratch/s11c/s9b-amend11-review-r1/` (Codex `codex-8S0Pna-*`, Grok `defs-r1-*`, `resolve-wl-r1-*`, `shows-wl-r1-*`). All completed.

## Round 2

**Artifact.** The working-tree directive at sha256 `d8a54eda…`, Codex-written: author session note
`_scratch/s9b_build/amend11_author/author_note.md`, brief `brief.md` beside it. The orchestrator copied the author's
text in unchanged. The diff against round 1 touches only the amendment-11 header and item 14's amendment-11
paragraphs. Saved before launch in `_scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r2.md`. Prompt:
`_scratch/s9b_build/s9b_repair_build_directive_amend11_review_prompt_r2.md` (`6ac2be30…`).

**Legs** (the text is Codex-written, so a fresh Claude agent and Grok, G1):
- **Fresh Claude (opus): NEEDS REVISION, five findings.** Report: `_scratch/s9b_build/s9b_amend11_review_r2_claude.md`,
  saved verbatim from the hand-back. Scripts and stdout: `_scratch/s9b_build/s9b_amend11_review_r2_claude_evidence/`.
- **Grok: NEEDS REVISION, four findings.** Report: `_scratch/s9b_build/s9b_amend11_review_r2_grok.txt`. Script and
  stdout: `_scratch/s9b_build/s9b_amend11_review_r2_grok_evidence/`.

Both legs confirm what rounds 0 and 1 required:
- disjoint domains no longer print zero;
- a value equal on its domain no longer prints as a change;
- K10's trace starts after the negation;
- K11's residual and unreached tags cannot meet a knife's obligation;
- no leakage.

**Findings.**
- **Both legs (Claude 1, Grok 3): identical worded conditions print zero when their translated predicates differ.**
  Rule 1 prints zero for "the same condition on the same variables" (L289). Round 2 makes a worded condition a
  predicate computed from each copy's own constructions (L325–326). Claude's K4a case: `WHOLE_DOMAIN_SAMEQ True`,
  while "required derivatives exist" for `δ` is `False` on the baseline and `True` on the knifed copy. Grok's
  classifier: carrier equality prints `ZERO`, and the translated predicates print `CHANGED`. This is a defect in
  material round 2 added.
- **Claude only:**
  - 2: "unchanged" for a value that carries a domain has two readings, everywhere or on the common domain; the site
    obligation's "differs" is undefined.
  - 3: being "defined" has no test, so a regularized value could stand where an ordinary one diverges.
  - 4: "its containing case" misses domains Wolfram attaches one level above the case.
  - 5: the print obligation is unbounded across witness candidates.
- **Grok only:**
  - 1: the witness point leaves domain parameters such as `bFar` unassigned.
  - 2: "the exterior half-line" is not one set (`(0,∞)` against `r > bFar`).
  - 4: worded entries outside the listed tests have no shared predicate.

**Disposition: held for the user.** This is the amendment's third round. The findings have gone from 3 to 2 to 9 as
the text grew, and round 2's longest addition bred the defect both legs found. Each round has found new cases in a
general rule for deciding whether any engine object changed. The review is past the point where the user asked to
be consulted (2026-10-11), and the pattern points at the premise, not the wording. So I am putting the premise to the
user before any fold.

**Process.** The legs' guarded runs are under `_scratch/s11c/s9b-amend11-review-r2/` (Claude `claudeleg-*`). All
completed. The Claude leg's guarded copy of the Wolfram engine reproduces the production transcript byte for byte.
