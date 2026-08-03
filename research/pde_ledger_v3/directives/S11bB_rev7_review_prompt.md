# Independent physics review — S11b-B BUILD DIRECTIVES, REVISION 7 (before any build runs)

## Artifact

- `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11bB_wl_directive.md`
- `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11bB_py_directive.md`

Engine-specific header + a shared physics specification, byte-identical in both
(`sha256 dd4887d5facdade1ece25c0cbdb4449d428ee962db1962c6cec0c07b4823f1da`).

## Why this review exists

Two independent engines will each implement these directives, and their agreement is what certifies the
result. **The directive is the one artifact both engines share.** An error in it lands in both, they agree,
and dual-engine certifies wrong physics.

⚠ **Six revisions have been rejected by independent review before any build ran**, and in each of the last
four rounds every finding sat in material that had just been changed to fix the previous round.
⇒ ⭐⭐ **Weight your attention toward what rev 7 changed**, listed below.

### What rev 7 changed

- **§3b · virtual-displacement rule** — rewritten. New `δ_v` notation, an explicit material-area mass
  `Σ_mat`, and a displayed VIRTUAL CONSTRAINT with a guard.
- **§3b · role of `J_±`** — a direct virtual-work contribution from `J_±` is now explicitly set to zero.
- **§3b · energy diagnostic** — signed source/sink accounting, plus a **new pressure-work sign check**.
- **§2b · two-port power** — the power identity was rewritten; port dissipativity restated.
- **§2b · reduction check** — the coefficient mapping is now fixed *before* the face solve, refitting is
  prohibited, and a new **affinity-power check** was added.
- **B8 control A** — now states that zero face velocity does not imply zero flux.

## ⭐⭐ FOUR SPECIFIC QUESTIONS — settle each by COMPUTATION, not by reading

These are contested. ⛔ **No position is supplied and you should not infer one.** ⭐ Decide each yourself and
show the work.

**Q1. Is the VIRTUAL-DISPLACEMENT RULE in §3b correct?** It distinguishes an Eulerian slab density from a
material-area mass, imposes `δ_v Σ_mat = 0`, and derives a constraint relating `δ_vθ`, `δ_ve_W` and
`∇_x·δ_v u`. ⚠ **Derive it independently.** Is the constraint right? Is the claim that the *Eulerian*
version must not be used correct? Could two engines still comply and get different equations?

**Q2. Is the TWO-PORT POWER IDENTITY in §2b complete?** With the bulk's normal velocity related to the face
velocity by `v_bulk,± = V_± + J_±/ρ_m`, ⭐ **construct the total energy flux into the bulk yourself** and
compare it against what the directive specifies. Is any term missing? Is `NOT_ESTABLISHED` appropriately
available, or does it license dodging a question that is actually answerable?

**Q3. Can the `Z_perm` REDUCTION CHECK be forced to pass?** Engines construct an affinity `𝒜` and report a
coefficient mapping. ⚠ Judge whether an affinity with a **wrong bulk-side sign or normalization** could
still reproduce the supplied `Z_perm` by choosing the mapping — and whether the new ordering constraint and
the affinity-power check actually prevent that.

**Q4. Does the new PRESSURE-WORK SIGN CHECK catch a reversed traction sign?** ⚠ **Test it:** reverse the
sign of the bulk traction in a minimal model with positive radiation resistance, and determine whether the
directive's checks — causality, energy accounting, and the new pressure-work check — would flag it, or
whether all of them pass while the roots move from decay to growth.

## Also check

**5. ⭐⭐ DOES THE SPECIFICATION PRESUPPOSE THE FORM OF ANY ANSWER?** Go task by task. ⛔ An earlier revision
was caught asserting where the deliverable's roots lie. Flag any statement that constrains what the
engines may find.

**6. ⭐ IS THE INSTABILITY CHANNEL STILL OPEN?** A **growing** root must remain an admissible, first-class
outcome. ⚠ Several diagnostics were added this round — check that none of them, alone or combined, can be
used to reject a genuine growing root, and that each one that could is **explicitly scope-bounded** to a
sub-case where growth would necessarily be a derivation error.

**7. ⭐ IS EVERY CHECK ABLE TO FAIL?** The directive now states, for each check, what wrong derivation it
catches. ⚠ **Verify those claims** — a check that passes by construction is worse than no check.

**8. ⭐ IS THE PHYSICS CLOSED?** Count equations against unknowns explicitly.

**9. Do the two directives specify the same physics?** The shared block is byte-identical; the risk is the
headers. **10. ⛔ Does either directive leak an expected answer, or contain an ill-posed task?**

## ⛔ CLEARED BY INDEPENDENT LEGS IN EARLIER ROUNDS — do not re-open

The scope boundary · header symmetry and the `reduction/` bar · `B1`'s supplied mass balance · the `A/B/C`
step split · `B8` controls B/C/D as form cuts · **§1b's branch prescription as such** (verified
computationally by three independent parties — ⭐ but *do* report any defect in how it is now stated).

⚠ Re-opening these costs a round. ⭐ Unless a rev-7 change has newly broken one — then say which change did.

## Do not read

- Anything named `PREREGISTERED`/`PREREG`, and blob `a68245b4fc5eb402dc29aa245a335d2caa453460`.
- ⛔⛔ `research/pde_ledger_v3/_scratch/` — it contains **prior reviewers' reports and the authoring
  transcript for this very revision.** Reading it would replace your independent judgement with theirs.
- ⛔ `research/pde_ledger_v3/steps/S11b_HANDOFF.md`, and the **git history** of the directive files — both
  record the reasoning behind these changes. **Judge rev 7 as it stands.**
- `research/pde_audit/`.

⭐ You **may and should** read `research/pde_ledger_v3/steps/S11bA_interface_response.md`, `V3_STEP_PLAN.md`
and `reduction/` to judge completeness. Those are barred from the *builders*, not from you.

## Required method

⚠ **This is not a document-only review.** Q1–Q4 require computation; use SymPy/numerics freely. Form your
own view of what this step must establish **before** judging whether the directive establishes it.
⚠ Where you believe a task will produce a specific answer, ⛔ do not report that answer — report only
whether the task as written can produce it **and its negation**.

## Physics filter

Report a finding only if it catches a way the physics could be wrong, or a way the dual-engine check could
be defeated. ⛔ No style, no formatting, no check-strength auditing beyond what checks 6 and 7 ask.

## Output

⭐ **Answer Q1–Q4 explicitly and first**, each with a one-word verdict and the computation behind it.
Then any further findings, most severe first: file and line, one sentence stating the defect, and a
concrete scenario. ⭐ Mark each NEW-IN-REV-7 or PRE-EXISTING. ⭐ State for each of Q1–Q4 whether you
**verified by computation** or only read. If nothing survives the filter, say so plainly. Under 60 lines.
