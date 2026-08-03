# Independent physics review — S11b-B BUILD DIRECTIVES, REVISION 6 (before any build runs)

## Artifact

- `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11bB_wl_directive.md`
- `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11bB_py_directive.md`

Engine-specific header + a shared physics specification verified byte-identical in both
(`sha256 18c3e9faafcc0c1b6036e8cc658743a0dd15914ba2108459290177e28b3b0f33`).

## Why this review exists

Two independent engines will each implement these directives, and their agreement is what certifies the
result. **The directive is the one artifact both engines share.** An error in it lands in both, they agree,
and dual-engine certifies wrong physics.

## ⛔⛔ THE PATTERN YOU ARE HERE TO BREAK

**Five revisions have been rejected by independent review before any build ran.** ⚠⚠ **In each of the last
three rounds, EVERY finding was located in the material that had just been CHANGED to fix the previous
round.** The fixes keep breeding the next round's defects.

⇒ ⭐⭐ **Weight your attention toward what rev 6 changed.** A fix that creates a fresh way for the physics to
be wrong is the outcome this review exists to prevent.

### What changed in rev 6 — ⛔ the findings behind these are deliberately withheld

- **§1b** — two sentences removed; a supplied branch prescription otherwise unchanged from the previous
  revision, which two independent legs verified by computation.
- **§2b** — the closure is **re-parameterized**: the flux is now driven by an **affinity `𝒜` that the
  engines must CONSTRUCT** rather than a supplied combination. The passivity requirement is rewritten as a
  **two-port** question and both parts may now return `NOT_ESTABLISHED`.
- **§3b** — **substantially rewritten.** A Lagrange multiplier is now permitted where the previous revision
  banned it; a binding **virtual-displacement rule** is added; constitutive quantities become **variational**
  derivatives; explicit geometry/sign conventions are supplied for the bulk load; an overclaim about
  stationarity is narrowed and a doubled-variable cross-check is now permitted; **two convention
  cross-checks** are added with an explicitly bounded scope.
- **B7** — a units-restored homogeneity check on every final equation, plus a mandatory ablation of it.
- **B6 / B8** — symbol updates following §2b's re-parameterization; control **E** re-scoped.

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

**3. ⭐⭐ IS §2b RIGHT?** Re-parameterized, so judge it cold:
(a) is the **interfacial mass balance** `v_bulk,± = V_± + J_±/ρ_m` correct, including signs?
(b) ⭐⭐ **is the affinity `𝒜` well-posed as a task** — can an engine construct it uniquely from entropy
    production given what is supplied, or could two engines build different affinities and both comply?
    ⚠ Check the **normalization** point specifically: `U` is per in-plane 3-volume with `θ` fractional,
    while the bulk pressure is a 4-volume quantity.
(c) ⭐ is the reduction **acceptance check** genuinely **able to fail**, given the engines now also report
    the parameter mapping used to reach it — or can that mapping be chosen to force a pass?
(d) is the **two-port** passivity treatment correct, and is allowing `NOT_ESTABLISHED` for both parts right,
    or does it let engines dodge a question that is actually answerable?

**3b. ⭐⭐ IS §3b's VIRTUAL-DISPLACEMENT RULE CORRECT AND SUFFICIENT?** ⚠ **Verify it by computation.**
It asserts that a virtual variation transfers no mass, hence `δΣ = 0` binds the variations, hence `U` must
**not** be varied at fixed `θ`. Judge: is that right; does it leave the derivation **exactly one** outcome;
and do the two supplied **convention cross-checks** actually discriminate the readings, ⛔ or can a wrong
derivation pass them? ⚠ Also judge whether their **stated scope limitation** is correctly drawn — they must
not be usable to reject a growing root of the full problem.

**4. ⭐⭐ DOES THE SPECIFICATION PRESUPPOSE THE FORM OF ANY ANSWER?** Go task by task, including unchanged
ones. ⚠ Rev 6 supplies more than its predecessors did — **judge whether anything now supplied is really a
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
  both record the reasoning behind rev 6's changes, and this review is worth nothing if it reads that
  reasoning instead of forming its own. **Judge rev 6 as it stands.**

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
showing how it lets wrong physics through or defeats dual-engine. ⭐ **Mark each finding NEW-IN-REV-6 or
PRE-EXISTING.** ⭐ For checks 1, 2 and 3 state explicitly whether you **verified** the supplied material or
only read it. If nothing survives the filter, say so plainly. Under 60 lines.
