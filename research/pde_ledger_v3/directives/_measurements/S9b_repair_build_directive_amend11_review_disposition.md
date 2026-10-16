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

**The user's decision (2026-10-16).** Two choices, both made through the question tool. The option texts are mine.
Both are quoted from the session transcript.
- 09:27 UTC. Question: "What should round 8 have to clear before the harness is accepted?" Chosen: "Narrowed bar
  (Recommended)", whose description reads: "Round 8 must clear R1 everywhere (no false zeros) and R4 at every
  knife's target in both engines. 'Undecided' is allowed away from the targets, and the record claims nothing from
  those labels." Rounds 0–2 were written to this choice.
- 11:13 UTC, after the round-2 findings. Question: "How should we get the harness finished?" Chosen: "Smaller
  harness (Recommended)", whose description reads: "Replace the general 'did it change' rule with: 'identical' only
  for literally identical outputs, otherwise 'not compared', plus one concrete, printed witness per knife per engine
  that the reviewers check. A new author writes a short amendment, reviewed until clear". The other option was to
  keep refining amendment 11's general rule.

So round 2's findings are not folded. Round 3 replaces the general rule. Item 14's per-part classification (zero,
change, undecided) is removed. Every tag is marked identical, when its two payloads are structurally the same
expression, or not compared. Each knife has one concrete printed witness in each engine. R1 is answered because
nothing prints a zero unless the payloads are identical. R4 is answered by the witness, which must lie at a site
part. R2 is answered because the witness's attached conditions must be shown to hold at its point. Round 2's
findings bear on round 3 only where they concern the witness test: free symbols such as `bFar`, definedness,
conditions at every payload level, and worded conditions.

## Round 3

**Author.** The author changes again (G4): a fresh Codex author (gpt-6.1-sol, xhigh; session `01a1446c…`). It
starts from the last accepted directive (`4b4db191`), so none of rounds 0–2's text carries over. Brief and note:
`_scratch/s9b_build/amend11_author2/brief.md` and `author_note.md`. The brief states what must be true, from the
user's choice.

**Artifact.** The working-tree directive at sha256 `7d3283a4…`, Codex-written. The orchestrator copied it in
unchanged after checking that the diff against `4b4db191` touches only these:
- the amendment-11 header;
- the earlier amendment headers that named the removed classification (amendments 3, 5, 6 and 7);
- item 14;
- K11's printing lines;
- the Part 2 and Part 3 round-8 lines.

Saved before launch in `_scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r3.md`. Prompt:
`_scratch/s9b_build/s9b_repair_build_directive_amend11_review_prompt_r3.md` (`dd3e3146…`). Baseline: 17 files, all
unchanged after both legs. Lookups: `S9b_repair_build_directive_amend11_review_r3_lookups.md` beside this file
(gen `_scratch/s9b_build/gen/s9b_amend11_r3_lookups.sh`).

**Legs** (Codex-written, so a fresh Claude agent and Grok, G1):
- **Fresh Claude (opus): NEEDS REVISION, six findings.** Report: `_scratch/s9b_build/s9b_amend11_review_r3_claude.md`,
  saved verbatim from the hand-back. Scripts and stdout: `_scratch/s9b_build/s9b_amend11_review_r3_claude_evidence/`.
- **Grok: NEEDS REVISION, three findings.** Report: `_scratch/s9b_build/s9b_amend11_review_r3_grok.txt`. Scripts and
  stdout: `_scratch/s9b_build/s9b_amend11_review_r3_grok_evidence/`.

Both legs confirm these:
- the marks claim only structural identity, and nothing prints zero, change or undecided at a tag;
- one site rule covers K6, K10 and K11, and excludes the digest and K11's residual;
- the point assigns domain parameters such as `bFar`;
- finiteness is ordinary;
- a witness's point and profile stay out of every engine object and out of item 13;
- "The engine is unchanged" and "Print no verdict" are consistent with item 14;
- there is no leakage.

Claude confirms that no text outside K11 relies on the removed classification.

**Findings and dispositions.** All are accepted. Each changes what a harness prints, or whether a printed witness
shows a change in the emitted object. The column says what must be true; the wording is the author's.

| # | Finding (evidence in the lookups file) | Must be true in round 4 |
|---|---|---|
| C1 | K11 still prints "their difference" and a residual, and nothing governs them. The text that held them to item 14's rule was removed. Both harnesses compute them with the r7 classifier, which prints `ValueDifference -> 0` for `x` on `x>0` against `2x` on `x<0` (WL C1–C3; SymPy L1020). R1 survives at K11. | No K11 difference or residual prints or implies a zero or a change. K11 pairs carry the same marks as tags. Any K11 claim of change is a witness. |
| C2 | "Structurally the same under a structural comparison the harness names" admits a normalizing comparison. `SameQ[Together[e1],e2]` is True while `e1` is `Indeterminate` at `x=1` and `e2` is 2. `cancel` gives the same in SymPy. Machine `SameQ` is True across a 2⁻⁵² difference. | Identical means the payloads are the same expression as emitted. Nothing is evaluated, simplified or canonicalized first. Approximate numbers are identical only when they are exactly equal. |
| C3 | A domain witness needs one copy's condition true and the other false, but nothing about the conditions that contain it. Container `x>−1`, `x>0` against `x>0 ∨ x<−5`, at `x=−6`: the part conditions differ and the effective sets are equal. | In a domain witness, every condition at the containing levels, other than the compared one, is computed true in both copies. |
| C4 | Which entries are attached conditions is not stated. SymPy's `RHO_POSITIVE`/`RHO_NONPOSITIVE` are branch labels at the level containing every C-row object; read as conditions, no SymPy K9 value witness can exist (C2 `EmptySet`). SymPy attaches a set, `Interval.open(b_min, ∞)`, to `RADIAL_DIFFERENTIAL_IDENTITY` over `r` (L966). Read pointwise, it splits from the every-`r` reading at K7 and K1–K4 (D1/D2). Labels such as `"Measure"` and `"Profiles"` are left to the harness. | The text says which entries are attached conditions. They are the entries stating where the part's value holds: domains, gates and premises. A branch label qualifies only its own branch. Descriptive strings and entries that are themselves computed values (e.g. `Closedness`) are not conditions. The harness prints which entries it evaluated as conditions. A set attached to a relation over a variable is that variable's range. Such a relation is decided over the whole range, as the every-`b` reduction (item 6) states it. |
| C5 | The admissible-profile list omits `j_n`. The spec makes `j_n` live (L80). Both engines attach `j_n = 0` to every forward object (SymPy L1050, L1070; WL `"NormalExchange"`, L659–661). | `j_n` is on the list where it enters. |
| C6 | "A printed bound covering both errors" and "computed true" are met by a quadrature error estimate (E2: the values differ by 1.8×10⁻⁴; the estimates sum to 10⁻³⁶) or by sampling a quantifier (F1 all true; F2 counterexample `π`). | A numeric bound is a proved enclosure: exact, or outward-rounded interval or ball arithmetic. A condition is computed true or false only by an exact decision over its whole stated range. |
| G1 | "Computed true" is not one test. `TrueQ` is False on a true `Inactive[ForAll]`, which `Resolve` decides as True. SymPy's `.doit()` leaves a ConditionSet emptiness undecided. `regularityDomain` mixes a formula, worded clauses and order-counting strings. `Closedness` is a value under one reading and a condition under another. | One decision procedure for every condition: an exact decision of the condition as stated, with its quantifiers active and its bound variables bound. A test that leaves a form inactive or unevaluated decides nothing. A worded condition is decided as its mathematical content for that object, as the directive and spec define it. The exterior half-line is fixed as one set. Order-counting strings (`OpticalCounting`, the spec's order counting) are not evaluated. |
| G2 | "Every free symbol" by a symbol scan includes bound variables: a quantifier variable, an integral variable, a ConditionSet binder. Substituting a binder makes a false emptiness statement `True` (SymPy) or breaks the integral (WL). | The point assigns only the symbols that remain free once the bound variables are set aside. Bound variables are never substituted. |
| G3 | "It need not meet the condition being compared" names no condition now that zero, change and undecided are gone. It reads as letting a value witness's profile leave its domain: `sqrt(x²)` against `x` at `x=−1` gives finite, differing values while `x>0` is false. | The phrase goes. In a value witness the profile makes every attached condition true. In a domain witness the compared condition is true in one copy and false in the other. |

Round 2's witness-related findings are met or carried here: definedness and free symbols by both legs' checks; worded
entries and the exterior half-line by G1; containing levels by C3.

**Author for round 4.** Round 3's findings are gaps in new text, not defects bred by a fold, so the author stays
Codex (gpt-6.1-sol, xhigh). The fold goes to a fresh session with this table, both reports and both evidence
directories. The text is to stay short, saying what must be true and not how to compute it.

**Process.** The legs' guarded runs are under `_scratch/s11c/s9b-amend11-review-r3/` (Claude `claudeleg-*`, Grok
`bind-wl-1`, `cond-*`, `trueq-wl-1`). All completed.

## Round 4

**Artifact.** The working-tree directive at sha256 `166c108d…`, written by a fresh Codex session (gpt-6.1-sol,
xhigh; session `01a1448c…`). Brief and note: `_scratch/s9b_build/amend11_author3/brief.md`, `author_note.md`. The
orchestrator copied it in unchanged after checking that the diff against round 3 (`a65f752c`) touches only item 14
and K11's printing lines. Saved before launch in
`_scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r4.md`. Prompt:
`_scratch/s9b_build/s9b_repair_build_directive_amend11_review_prompt_r4.md` (`d85f21f7…`). Baseline: 18 files, all
unchanged after both legs. Lookups: `S9b_repair_build_directive_amend11_review_r4_lookups.md` beside this file
(gen `_scratch/s9b_build/gen/s9b_amend11_r4_lookups.sh`).

**Legs** (Codex-written, so a fresh Claude agent and Grok, G1):
- **Fresh Claude (opus): NEEDS REVISION, three findings.** Report: `_scratch/s9b_build/s9b_amend11_review_r4_claude.md`,
  saved verbatim from the hand-back. Scripts and stdout: `_scratch/s9b_build/s9b_amend11_review_r4_claude_evidence/`.
- **Grok: NEEDS REVISION, three findings.** Report: `_scratch/s9b_build/s9b_amend11_review_r4_grok.txt`. Scripts and
  stdout: `_scratch/s9b_build/s9b_amend11_review_r4_grok_evidence/`.

Both legs confirm these:
- round 3's C1, C3, C5, G2 and G3 are met; Claude also finds C2, C6 and G1 met, and C4 met as far as it reaches;
- K11 prints ordered operand pairs with the tag marks and no subtraction;
- what round 3's legs confirmed still holds;
- there is no leakage.

**Findings and dispositions.** Five are accepted and one is rejected. Two of the accepted ones (C4′ and G1′ below)
trace to wording in round 3's table that the author used as written: "a set attached to a relation" and "a branch
label". So this table states outcomes only. The wording is the author's.

| # | Finding (evidence in the lookups file) | Disposition: what must be true in round 5 |
|---|---|---|
| C1′ | Ordered cases hide a domain witness. Wolfram's `branchTypes` is a `Piecewise`: a later case holds only where every earlier case fails. K10 makes the third case's condition true in one copy and false in the other, with no containing condition, while the emitted classification at that point is the same in both (A1–A4). | **Accepted.** An attached condition includes whatever the engine's own structure makes a case depend on, including the failure of earlier cases in an ordered case object. A domain witness shows that the emitted object differs at its point. |
| C2′ | A symbol each copy's own construction fixes gets one shared value. Item 5 fixes each turning radius by that ray's own construction (L157–158). Wolfram's gates leave `rTurnFlyby`, `rTurnRadar` and `impactRadar` free, each fixed by `turningCondition` or `radarEndpointEquation`. Under K3 the roots differ: (3+√21)/2 against 3 (B2), so a shared value meets only one copy's construction (B3). SymPy's free `s9b_velocity_zero_radius` (76 occurrences in `B_IMPLIED_JN_DEFLECTION`) is the same kind of symbol. | **Accepted.** Here a quantity that varies between copies would be held fixed (M3). A symbol that a copy's own construction fixes takes, in each copy, the value that copy's construction gives, printed with its computation. The point assigns only the independent inputs. |
| C4′ | The range rule covers only SymPy's attached set. Wolfram emits the reduced identity as `"ProfileEquation"` with `"ProfileEquationRadius" -> r`, and its range is the domain formula `r > bFar`. Read pointwise against every-`r`, the two readings split (P1). P2 shows that the every-`r` reading can furnish a witness. | **Accepted.** In both engines, the radius of a reduced every-`b` condition is decided over its stated range, whatever carries that range. |
| G1′ | The density-sign labels are not tied to one part. Item 9 makes each fixed-ratio and power response one object on both parts. SymPy's tuple carries `RHO_POSITIVE` and `RHO_NONPOSITIVE` beside `OBJECT_ON_BOTH_PARTS`. Their conjunction is `EmptySet`, so a reader who counts both as conditions finds no SymPy K9 witness. Wolfram nests one `LabelCondition` in each part. "Branch label" elsewhere means branch existence. | **Accepted.** Item 9's density-sign labels are not attached conditions of the shared object. A witness may lie on either part, and it prints which part its point lies on. Parts at `ρ ≤ 0` are kept, not excluded. |
| G2′ | Wolfram's worded domain still has more than one test. The engine's `abelInverse` uses `GenerateConditions -> False`, which returns `0` for a kernel whose ordinary limit is `Infinity` (ABEL). The half-line's coordinate is open: `r > bFar` (`Radius`) against `b > bFar` (`FarZone`, K7's domain) (HALFLINE). "Not evaluated" for `OpticalCounting`, together with "Unresolved evaluations furnish no witness", can block every witness under that domain. | **Accepted.** A worded "exists" or "converges" means ordinary existence or convergence. Each worded condition names the variable it ranges over in each engine. Order-counting entries are not attached conditions, so they cannot leave a witness unresolved. |
| G3′ | Unequal approximate payloads can still be marked identical: machine `SameQ[1., 1.+2^-52]` is `True`. | **Rejected: already met.** Round 4 says that approximate numbers "must be exactly equal" and that "their printing preserves every distinction" (L282–283). A comparison that marks that pair identical breaks the text as written. Whether the named comparison is exact on approximate numbers is a check for the build review r8, and its prompt will ask for it. |

**Author for round 5.** The same model, in a fresh Codex session, with this table, both reports and both evidence
directories. Its brief says to keep item 14 short and to state each outcome in its own words.

**Loop record.** Amendment 11 has had five rounds; the smaller harness, two. Findings in the smaller harness went
from 9 to 5 accepted, and each round has taken about 35 minutes. If round 5 does not clear, or its findings are again
new cases of the same witness test, the loop goes to the user before another fold.

**Process.** The legs' guarded runs are under `_scratch/s11c/s9b-amend11-review-r4/` (Claude `claudeleg-*`; Grok `py1`, `py2`,
`wl1`). All completed; two Claude runs (`claudeleg-py1`, `claudeleg-wl1`) failed on the leg's own script and were superseded.

## Round 5

**Artifact.** The working-tree directive at sha256 `efa14a51…`, written by a fresh Codex session (gpt-6.1-sol,
xhigh; session `01a144aa…`). Brief and note: `_scratch/s9b_build/amend11_author4/brief.md`, `author_note.md`. The
orchestrator copied it in unchanged after checking that the diff against round 4 (`46d32255`) touches only item 14.
Saved before launch in `_scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r5.md`. Prompt:
`_scratch/s9b_build/s9b_repair_build_directive_amend11_review_prompt_r5.md` (`d2673c87…`). Baseline: 19 files, all
unchanged after both legs. Lookups: `S9b_repair_build_directive_amend11_review_r5_lookups.md` beside this file
(gen `_scratch/s9b_build/gen/s9b_amend11_r5_lookups.sh`).

**Legs** (Codex-written, so a fresh Claude agent and Grok, G1):
- **Fresh Claude (opus): NEEDS REVISION, two findings.** Report: `_scratch/s9b_build/s9b_amend11_review_r5_claude.md`,
  saved verbatim from the hand-back. Scripts and stdout: `_scratch/s9b_build/s9b_amend11_review_r5_claude_evidence/`.
- **Grok: NEEDS REVISION, one finding.** Report: `_scratch/s9b_build/s9b_amend11_review_r5_grok.txt`. Scripts and
  stdout: `_scratch/s9b_build/s9b_amend11_review_r5_grok_evidence/`.

Both legs confirm these:
- the marks;
- the site rule;
- C1′ and G1′ in form;
- G2′'s half-line, ordinary existence and order counting;
- the K11 pairs;
- the rejected row's ground;
- engines unchanged, no verdict, no leakage.

Grok also finds C2′ and C4′ met.

**Findings.** All three are supported by the legs' scripts and accepted.
- **Claude 1 (C2′ only partly met).** "Takes that copy's own computed value" assumes exactly one value. Under K3, the
  baseline turning condition has no solution where K3's has one (F1a). The flyby turning condition has two roots,
  of which only one meets the rest of the construction (F1b). Wolfram fixes `impactRadar` by two equations
  (L253–260, L348). SymPy's free `s9b_velocity_zero_radius` has the zero set `{2, 4}` on one branch, so a builder
  who takes a different zero in each copy compares `j_n` at different radii, `−4` against `−8`, on a branch the knife
  does not change.
- **Claude 2 (C4′ generalized too far).** Nothing says which inequalities in a domain formula state a range. Wolfram's
  baseline `r > bFar && rAnchor > bFar` and K7's `bFixed > bFar` have the same form (L411–413, L421). Read as a range,
  K7's `bFixed` gives back the baseline's every-`b` quantifier, which undoes K7 in Wolfram alone (F2). Read as a range,
  SymPy's `Gt(z, s9b_b_min)` takes the join branch out of every witness. Where the radius is ranged, a density label
  holds on part of the range only, so "print which part contains its point" has no answer.
- **Grok 1 (G2′).** Wolfram's worded `"required derivatives exist"` has two ordinary tests: at the witness point, or
  at every radius on `r > bFar`. On `Abs(r−2)` at `r = 3` the K4a witness passes the first test and fails the second.

**Disposition: held for the user.** The loop record in round 4 set this condition: if round 5 does not clear, or its
findings are again new cases of the same witness test, the loop goes to the user before another fold. All three
findings are such cases.

**Process.** The legs' guarded runs are under `_scratch/s11c/s9b-amend11-review-r5/` (Claude `claudeleg-*`, of which
`py4`, `wl4` and `wl5` carry the verdict; Grok `py-witness-attempt-*`, `wl-witness-attempt-1`). All completed.
