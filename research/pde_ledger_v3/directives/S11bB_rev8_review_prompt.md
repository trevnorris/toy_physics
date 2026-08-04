# Independent physics review — S11b-B BUILD DIRECTIVES, REVISION 8 (before any build runs)

## Artifact

- `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11bB_wl_directive.md`
- `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11bB_py_directive.md`

Engine-specific header + a shared physics specification, byte-identical in both
(`sha256 92c139f39fb058bd6387c5ac983a47838d066bfbd5e6f46f577375b0aa60b760`).

## Why this review exists

Two independent engines will each implement these directives, and their agreement is what certifies the
result. **The directive is the one artifact both engines share.** An error in it lands in both, they agree,
and dual-engine certifies wrong physics.

⚠ **Seven revisions have been rejected by independent review before any build ran**, and in each of the last
four rounds every finding sat in material that had just been changed to fix the previous round.
⇒ ⭐⭐ **Weight your attention toward what rev 8 changed**, listed below.

### What rev 8 changed

- **§2b · the affinity is now SUPPLIED** as `𝒜 = μ_s − δp/ρ_m`, its construction task and derivation tag
  removed, and the affinity-power check **deleted**. ⚠ The author records the cost: this "sacrifices
  independent falsification of its sign and normalization."
- **§2b/§3b · a reciprocal traction `Λ_X(ω)𝒜_±`** is added with a free coefficient, propagating into the
  traction, the virtual work, the two-port power identity, the Onsager matrix, the dimension list, and a
  new FORM control **F**.
- **§1b** — the continuation is extended to the **upper** half-plane.
- **TASKS rule (4)** — scoped.
- **Diagnostics** — the conservative and pressure-work checks are now explicitly restricted to `Λ_X⁰ = 0`.

## ⭐⭐ FOUR SPECIFIC QUESTIONS — settle each by COMPUTATION, not by reading

These are contested. ⛔ **No position is supplied and you should not infer one.** ⭐ Decide each yourself and
show the work.

**Q1. ⭐⭐ IS THE SUPPLIED AFFINITY `𝒜 = μ_s − δp/ρ_m` CORRECT?** ⚠⚠ **This is now the ONLY check on it.**
It was previously an engine output; it is now supplied, so no dual-engine disagreement can surface an error
in it. ⭐ **Derive the affinity independently from interfacial entropy production** and compare — sign,
normalization, and which density each term is divided by. ⛔ Do not verify it against the directive's own
power identity; that is circular. ⚠ Its bulk-side **sign** determines whether the evanescent channel
dissipates or gains, hence whether a growing root appears — so an error here is not recoverable downstream.

**Q2. IS THE RECIPROCAL TRACTION `Λ_X` CORRECTLY AND COMPLETELY IMPLEMENTED?** ⭐ Verify by computation:
(a) does `Λ_X⁰ = 0` recover the previous specification **exactly** — traction, mechanics, power identity?
(b) is the traction's **sign and placement** right, given the outward-normal conventions?
(c) is the Onsager matrix as displayed the correct flux–force structure, and is the reciprocity relation it
    would force between `Λ_X` and `Λ_V` left as an **output** rather than asserted?
(d) does `Λ_X` appear everywhere it must — power, entropy production, thickness equation, dimensions,
    controls — or is there a place it was missed?

**Q3. ⭐⭐ IS THERE NOW A DIAGNOSTIC GAP AT `Λ_X⁰ ≠ 0`?** The conservative and pressure-work checks are
restricted to `Λ_X⁰ = 0`. ⚠ **Determine whether any check still constrains the derivation when
`Λ_X⁰ ≠ 0`** — or whether the very regime the new term was introduced to explore is now unpoliced. ⛔ If it
is unpoliced, say what minimal check would cover it **without** foreclosing a growing root.

**Q4. IS THE UPPER-HALF-PLANE CONTINUATION CORRECT AND CONSISTENT?** ⭐ Verify that the extended definition
agrees with the lower-half-plane prescription on the real axis and on the imaginary axis, and that a root
with `Im ω > 0` now has an unambiguous `q_out`.

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

The scope boundary · header symmetry and the `reduction/` bar · `B1`'s supplied mass balance · §3b's virtual-displacement rule (verified correct by computation) · the `A/B/C`
step split · `B8` controls B/C/D as form cuts · **§1b's branch prescription as such** (verified
computationally by three independent parties — ⭐ but *do* report any defect in how it is now stated).

⚠ Re-opening these costs a round. ⭐ Unless a rev-7 change has newly broken one — then say which change did.

## Do not read

- Anything named `PREREGISTERED`/`PREREG`, and blob `a68245b4fc5eb402dc29aa245a335d2caa453460`.
- ⛔⛔ `research/pde_ledger_v3/_scratch/` — it contains **prior reviewers' reports and the authoring
  transcript for this very revision.** Reading it would replace your independent judgement with theirs.
- ⛔ `research/pde_ledger_v3/steps/S11b_HANDOFF.md`, and the **git history** of the directive files — both
  record the reasoning behind these changes. **Judge rev 8 as it stands.**
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
concrete scenario. ⭐ Mark each NEW-IN-REV-8 or PRE-EXISTING. ⭐ State for each of Q1–Q4 whether you
**verified by computation** or only read. If nothing survives the filter, say so plainly. Under 60 lines.
