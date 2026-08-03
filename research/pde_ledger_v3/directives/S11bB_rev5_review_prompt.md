# Independent physics review — S11b-B BUILD DIRECTIVES, REVISION 5 (before any build runs)

## Artifact

- `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11bB_wl_directive.md`
- `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11bB_py_directive.md`

Engine-specific header + a shared physics specification verified byte-identical in both
(`sha256 b80b55680b5a7feac664d170047188bd176f1de8f1c4ec9b22759b7a4a2266ce`).

## Why this review exists

Two independent engines will each implement these directives, and their agreement is what certifies the
result. **The directive is the one artifact both engines share.** An error in it lands in both, they agree,
and dual-engine certifies wrong physics.

## ⛔⛔ THE PATTERN YOU ARE HERE TO BREAK

**Four revisions have been rejected by independent review before any build ran.** ⚠⚠ **In each of the last
three rounds, EVERY finding was located in the material that had just been CHANGED to fix the previous
round.** The fixes keep breeding the next round's defects.

⇒ ⭐⭐ **Weight your attention toward what rev 5 changed.** A fix that creates a fresh way for the physics to
be wrong is the outcome this review exists to prevent.

### What changed in rev 5 — ⛔ the findings behind these are deliberately withheld

- **§1b** — the complex-frequency section now **supplies** a branch prescription. Previous revisions either
  asserted a different one or supplied none.
- **§2b** — **new**. The interface closure gains a chemical-potential term; the face response is no longer
  supplied but must be derived; an acceptance check against an externally derived `Z_perm` is added;
  passivity is split into two computed questions.
- **§3b** — the derivation route is **replaced**. An action principle with a Lagrange multiplier is gone;
  balance laws replace it, three routes are forbidden, and the previous energy-accounting check is deleted
  and replaced.
- **§0 / B5** — a growing root now requires reported diagnostics.
- **§3** — the symmetry group is stated in full; the independence test is restated.
- **B8** — a new control **E**; control **C** widened.

## What to check

**1. ⭐⭐ IS §1b's SUPPLIED BRANCH PRESCRIPTION CORRECT, AND DOES IT UNIQUELY DETERMINE `q_out`?**
Highest priority. ⚠ **Do not take it on trust — derive it and TEST IT NUMERICALLY.** `q_out` has branch
points on the real `ω` axis and the modes of interest are complex. Specifically judge:
(a) does the prescription actually single-value `q_out` where the calculation needs it;
(b) is it consistent with the two real-axis requirements;
(c) is the instruction *not to re-impose those requirements at complex `ω`* correct, or does it discard
    something physical;
(d) could two engines still implement it differently?
⛔ **If the supplied prescription is wrong, that is the most important finding you can return** — a
previous revision supplied a precise rule that was wrong and it would have flipped the deliverable's sign.

**2. ⭐⭐ IS §3b's SUPPLIED DERIVATION ROUTE CORRECT?** ⚠ Judge independently, and **check the claims it
makes**: that a retarded kernel cannot be varied; that each forbidden route is genuinely wrong; that the
stated signatures of the forbidden routes are real; and that the mandated balance-law route is
**unambiguous** — two engines following steps 1–6 must not be able to write different equations.
⭐ Also judge whether the **causality diagnostic** and the **energy-sink attribution** can actually fail, or
whether an engine following the prescription passes them by construction. ⚠ The check they replaced was
un-failable, which is why it was removed; ⛔ do not let the replacement have the same defect.

**3. ⭐⭐ IS §2b RIGHT?** New material, so judge it cold:
(a) is the **interfacial mass balance** `v_bulk,± = V_± + J_±/ρ_m` correct, including signs?
(b) is the modified closure `J_± = Λ_p δp − Λ_μ μ_θ + Λ_V V_±` thermodynamically sensible, and is `μ_θ` the
    right conjugate variable?
(c) ⭐ is the `Λ_μ⁰ = 0` **acceptance check** against the supplied `Z_perm` genuinely **able to fail**, or
    could an engine satisfy it trivially?
(d) are the two passivity questions actually computable from what is supplied?

**4. ⭐⭐ DOES THE SPECIFICATION PRESUPPOSE THE FORM OF ANY ANSWER?** Go task by task, including unchanged
ones. ⚠ Rev 5 supplies more than its predecessors did — **judge whether anything now supplied is really a
result being smuggled in as setup.**

**5. ⭐ IS THE PHYSICS CLOSED UNDER THE NEW ROUTE?** The route changed, so **recount** equations against
unknowns explicitly. With `θ` never eliminated and the mass balance an independent evolution equation, is
the system determined? Is anything required and absent?

**6. ⭐ IS THE INSTABILITY CHANNEL BOTH OPEN AND PROTECTED?** §0 and B5 admit a growing root as a
first-class outcome, and now require diagnostics to separate a real one from an artifact. ⚠ Check that
nothing elsewhere forecloses growth, **and** that the diagnostics do not amount to a licence to discard it.

**7. ⭐ IS `B1` STILL CORRECT** given §2b's changes, and does linearising the supplied balance necessarily
produce every term it should?

**8. DO THE TWO DIRECTIVES SPECIFY THE SAME PHYSICS?** The shared block is byte-identical, so the risk is
the headers.

**9. ⛔ Is any task ill-posed or fabrication-forcing?** ⛔ Does either directive LEAK an expected answer?

## ⛔ CLEARED BY BOTH INDEPENDENT LEGS IN EARLIER ROUNDS — do not re-open

Scope boundary explicit and not leaked · header physics symmetric and `reduction/` barred from both · no
pre-registration leak · `B8` controls B/C/D are genuine form cuts · tractability and orphan check clean ·
`B1`'s supplied balance is the correct integrated mass balance.

⚠ Re-opening these costs a round. ⭐ **Unless a rev-5 change has newly broken one** — then report it and say
which change did.

## Do not read

- Any file whose name contains `PREREGISTERED` or `PREREG`, and blob
  `a68245b4fc5eb402dc29aa245a335d2caa453460`.
- Any other reviewer's output, and `research/pde_ledger_v3/_scratch/`.
- `research/pde_audit/`.
- ⛔ The **git history** of the directive files, and ⛔ `research/pde_ledger_v3/steps/S11b_HANDOFF.md` —
  both record the reasoning behind rev 5's changes, and this review is worth nothing if it reads that
  reasoning instead of forming its own. **Judge rev 5 as it stands.**

⭐ You **may and should** read `research/pde_ledger_v3/steps/S11bA_interface_response.md`, `V3_STEP_PLAN.md`,
and `reduction/` to judge completeness. Those are barred from the *builders*, not from you.

## Required method

DOCUMENT review, but ⚠ **not document-only**: checks 1 and 2 require you to **compute**. Use SymPy/numerics
freely. Form your own view of what this step must establish **before** judging whether the directive
establishes it. Quote both the directive and your own reasoning for every finding.
⚠ Where you believe a task will produce a specific answer, ⛔ do not report that answer — report only
whether the task as written can produce it **and its negation**.

## Physics filter

Report a finding only if it catches a way the physics could be wrong, or a way the dual-engine check could
be defeated. ⛔ No style, no formatting, no check-strength auditing — hardening is a separate phase.

## Output

Findings, most severe first: file and line, one sentence stating the defect, and a concrete scenario
showing how it lets wrong physics through or defeats dual-engine. ⭐ **Mark each finding NEW-IN-REV-5 or
PRE-EXISTING.** ⭐ For checks 1, 2 and 3 state explicitly whether you **verified** the supplied material or
only read it. If nothing survives the filter, say so plainly. Under 60 lines.
