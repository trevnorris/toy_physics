# Independent physics review — S11b-B BUILD DIRECTIVES, REVISION 4 (before any build runs)

## Artifact

- `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11bB_wl_directive.md`
- `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11bB_py_directive.md`

Engine-specific header + a shared physics specification verified byte-identical in both
(`sha256 54e53adee5460c1f803dec7ef1f433c5eb7bbe11112bc85e53f272f132064827`).

## Why this review exists

Two independent engines will each implement these directives, and their agreement is what certifies the
result. **The directive is the one artifact both engines share.** An error in it lands in both, they agree,
and dual-engine certifies wrong physics.

⚠ **That is not hypothetical here.** On the immediately preceding sub-step a task asked for *"the order in
`v₀/c_s0`"*; both engines dutifully returned a power of `v₀/c_s0`, and the correct answer was a different
quantity altogether which exceeds first order in the very regime that mattered. The specification had
**presupposed the form of the answer**, and three review rounds missed it because each asked whether the
script matched the spec.

## ⛔⛔ WHAT MAKES THIS REVISION DIFFERENT — READ THIS BEFORE ANYTHING ELSE

**Three previous revisions were rejected by independent review before any build ran.** Rev 1 and rev 2
failed *structurally*. ⚠⚠ **Rev 3's findings were then all located in the parts that had been CHANGED to
fix rev 2** — each rewrite bred the next round's defects, of the same class.

⇒ ⭐⭐ **THE PRIMARY QUESTION FOR YOU IS WHETHER REV 4'S CHANGES INTRODUCED NEW DEFECTS.** Weight your
attention toward the changed material, listed below. A fix that creates a fresh way for the physics to be
wrong is the outcome this review exists to prevent.

⭐ **The sharpest instance from last round:** rev 2's branch-of-`q` rule was *vague but safe*; rev 3
replaced it with a *precise and wrong* one, which would have flipped the sign of the deliverable.
⇒ **A wrong precise answer is worse than an acknowledged ambiguity.** ⚠ Rev 4 responded by removing the
rule and asking the engines to derive it — **judge whether that trade was made correctly** (see check 2).

### What changed in rev 4, by location — ⛔ the *findings* behind these are deliberately withheld

`§0` scope · `§1` complex-frequency block · `§2` permeable-face paragraph, plus a new limitation note and
the validity note · `§3` energy-basis instruction · `§3b` derivation prescription and a new accounting
rule · `B1` · `B4`/`B5`/`B6` reporting requirements · `B7` list · `B8` control A and the transverse note.

## What to check

**1. ⭐⭐ DID ANY REV-4 CHANGE INTRODUCE A NEW DEFECT?** Highest priority. For each changed location, judge
it *on its own merits* as if seeing it fresh: does it presuppose an answer, assert something that should be
computed, under-determine something two engines must agree on, or close off an outcome?

**2. ⭐⭐ IS THE COMPLEX-FREQUENCY / BRANCH SECTION NOW UNDER-DETERMINED?** §1 supplies **no** branch rule
and instead gives three physical requirements and asks each engine to derive, state and justify the
prescription. ⚠ **Judge whether requirements 1–3 actually determine `q_out(ω,k)` uniquely.** If they do
not, two engines can derive different prescriptions and the deliverable's imaginary part is unresolved —
**say so, and say what minimal additional statement would fix it without supplying the answer.** ⛔ Do not
supply the prescription itself in your report.

**3. ⭐⭐ DOES THE SPECIFICATION PRESUPPOSE THE FORM OF ANY ANSWER?** Go task by task, including the
unchanged ones. Flag any place where phrasing constrains the shape of a result — a named expansion
parameter, an assumed proportionality, a fixed number of regimes, a dimension only one candidate could
have.

**4. ⭐⭐ IS THE NEW "SINGLE ACCOUNTING RULE" IN §3b SOUND, AND IS IT GENUINELY ABLE TO FAIL?** It removes
the Rayleigh dissipation function, declares `J_±` determined by the closure plus the bulk solution, and
frames "`Z_perm` is the only place the interface enters the mechanical equations" as a **hypothesis to
test** with a loss-channel enumeration and a power-balance cross-check. ⚠ **Judge: (a) is the physics
right; (b) is dissipation now counted exactly once; (c) can the stated checks actually refute the
hypothesis, or would an engine that assumed it pass them trivially; (d) is there a loss channel the
framing makes invisible?**

**5. ⭐ IS THE PHYSICS CLOSED, GIVEN THE CHANGED DERIVATION ROUTE?** The route changed, so **recount**:
equations against unknowns, explicitly. Does varying `u`, `δW`, `θ`, `λ` plus the external face force
determine the system? Is anything required and absent?

**6. ⭐ IS `B1` CORRECT AS WRITTEN?** The exact balance is now supplied rather than left as a task.
⚠ **Verify the balance itself is right**, and that linearising it **necessarily** produces every term it
should — in particular that no term can be dropped as higher-order that is not. ⛔ If the supplied balance
is wrong, both engines inherit it.

**7. ⭐ IS THE `§3` ENERGY-BASIS TASK WELL-POSED?** Engines must construct the closed basis of quadratic
invariants and classify omissions as independent or redundant *modulo `B1`'s constraint*. ⚠ Judge whether
that has a determinate answer, whether the symmetry statement given is sufficient to pin the basis, and
whether "redundant modulo the constraint" is a test an engine can actually apply.

**8. ⭐ IS THE INSTABILITY CHANNEL NOW GENUINELY OPEN?** §0 and B5 admit a growing root as a first-class
outcome. ⚠ Check that nothing **elsewhere** in the specification still forecloses it — a positivity
assumption, a branch instruction, a passivity requirement, or a task phrased so a growing root reads as an
error.

**9. DO THE TWO DIRECTIVES SPECIFY THE SAME PHYSICS?** The shared block is byte-identical, so the risk is
the headers. Does either add, remove, reweight or reinterpret a task, or give one engine information the
other lacks?

**10. ⛔ Is any task ill-posed or fabrication-forcing?** A task whose only honest outcome is an invented
value is worse than no task. ⚠ Several tasks were changed to *require a normalization statement before a
value* — check that none of those can only be answered by inventing a convention.

**11. ⛔ Does either directive LEAK an expected answer?** The author pre-registered predictions and must not
have supplied them.

## ⛔ CLEARED LAST ROUND BY BOTH INDEPENDENT LEGS — do not re-open

Scope boundary explicit and not leaked · header physics symmetric and `reduction/` barred from both · no
pre-registration leak · `B8` controls B/C/D are genuine form cuts · the closed count is square *given a
correct `B1`* · tractability and orphan check clean.

⚠ Re-opening these costs a round. ⭐ **Unless a rev-4 change has newly broken one** — in that case report it
and say which change did.

## Do not read

- Any file whose name contains `PREREGISTERED` or `PREREG`, and blobs
  `e6d40a292968d89adcc9b17962186ca4f509da4a`, `0acb69cfa703d075415e4874569270c6c06cf413`,
  `d5732c46fc6229f3924863547483ff690addea49`, `a68245b4fc5eb402dc29aa245a335d2caa453460`.
- Any other reviewer's output, and `research/pde_ledger_v3/_scratch/`.
- `research/pde_audit/`.
- ⛔ The **git history** of the directive files. Judge rev 4 as it stands; do not diff it against rev 3.

⭐ You **may and should** read `research/pde_ledger_v3/steps/`, `V3_STEP_PLAN.md`, and `reduction/` to judge
completeness. Those are barred from the *builders*, not from you.

## Required method

DOCUMENT review. Form your own view of what this step must establish **before** judging whether the
directive establishes it. Derive independently where you can; quote both the directive and your own
reasoning for every finding.
⚠ Where you believe a task will produce a specific answer, ⛔ do not report that answer — report only
whether the task as written can produce it **and its negation**.

## Physics filter

Report a finding only if it catches a way the physics could be wrong, or a way the dual-engine check could
be defeated. ⛔ No style, no formatting, no check-strength auditing — hardening is a separate phase.

## Output

Findings, most severe first: file and line, one sentence stating the defect, and a concrete scenario
showing how it lets wrong physics through or defeats dual-engine. ⭐ **Mark each finding NEW-IN-REV-4 or
PRE-EXISTING.** If nothing survives the filter, say so plainly. Under 60 lines.
