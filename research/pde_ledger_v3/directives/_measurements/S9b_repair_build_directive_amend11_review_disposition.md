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
