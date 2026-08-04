# Independent physics review — S11b-B BUILD DIRECTIVES, REVISION 10 (before any build runs)

## Artifact

- `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11bB_wl_directive.md`
- `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11bB_py_directive.md`

Engine-specific header + a shared physics block, byte-identical in both
(`sha256 6e5516a792cdb0fe472bd3215beb689ac4984eb2d4abcb8f566bbd68e1a55d92`).

## Why this review exists

Two independent engines will each implement these directives, and their agreement is what certifies the
result. **The directive is the one artifact both engines share.** An error in it lands in both, they agree,
and dual-engine certifies wrong physics.

⚠ **Nine revisions have been rejected before any build ran**, and in most rounds every finding sat in
material that had just been changed. ⇒ ⭐⭐ **weight your attention toward what rev 10 changed.**

### What rev 10 changed

- **§1b** — the framing around the branch prescription was made **conditional**; the prescription itself is
  untouched and is in the cleared set.
- **§2b** — the admissibility task now asks for a condition "in whatever form it takes"; a sentence
  privileging one relaxation-time case was removed; the revision-history paragraph was deleted entirely.
- **§3b** — the energy-balance check was declared an **off-shell coefficient check** with two degenerate
  readings explicitly forbidden, and its blind spot recorded; the **causality diagnostic** was restated in
  terms of **inherited pole location after removable cancellations** rather than symbol appearance.
- **§0b** — the **git history, object store, refs, reflogs, index, commit messages and patches** are now
  barred, as are `paper/` and named process documents. The **py header's** registry hint was removed.
- **§3/§3b** — divergence-invariant instruction made plural. **`Z_perm`** check's advertised catch narrowed.

## ⭐⭐ SPECIFIC QUESTIONS — settle each by COMPUTATION where computation applies

**Q1. ⭐⭐ CONTAMINATION SWEEP — the highest-priority check this round.**
§2b asks the engines to compute an **admissibility / reciprocity** result about the closure coefficients
`Λ_A⁰`, `Λ_V⁰`, `Λ_X⁰` and the relaxation times. ⚠⚠ **A leak was found and removed from this exact area in
the previous draft**: wording that named a particular candidate outcome and framed it as the status quo.
⇒ ⭐ **Sweep the whole directive for any remaining statement, phrasing, ordering, example, tag name, or
emphasis that names, hints at, presupposes, or makes easier to reach any particular outcome** of:
(a) which coefficients admissibility forces to vanish, if any;
(b) what relation reciprocity forces between `Λ_X` and `Λ_V`, if any;
(c) what relation among `τ_A`, `τ_V`, `τ_X` is forced, if any;
(d) the fate of the longitudinal mode, or the sign of any imaginary part.
⛔ **Report the leak; ⛔ do NOT report the answer you believe is correct.** ⚠ A phrase like "still requires"
or a tag naming one coefficient is exactly the failure mode — look for that class, not just literal values.

**Q2. ⭐⭐ IS THE OFF-SHELL ENERGY-BALANCE CHECK NOW UNAMBIGUOUS, AND DOES IT STILL DISCRIMINATE?**
⚠ It previously admitted three readings, two of which were degenerate — one false-FAILing a correct
derivation, one passing vacuously in both engines. ⭐ **Verify by computation:** (a) does the prescribed
computation now have exactly one meaning; (b) does it still flag a traction sign error, a factor error, an
omitted face, and a wrong slab-side conjugate; (c) is the recorded blind spot correctly and completely
stated, or does it hide further insensitivity; (d) does it classify any root or comment on any imaginary
part? ⛔ It must not.

**Q3. IS THE RESTATED CAUSALITY DIAGNOSTIC SOUND?** It now keys on **inherited pole location after
removable cancellations** instead of the appearance of a conjugated symbol. ⭐ Verify: (a) it no longer
false-FAILs when a correct derivation rationalizes a complex denominator; (b) it still catches a transposed
kernel in each of the three coupling channels — ⚠ including the one whose pole is **displaced by the
interface coupling**; (c) it does not fire on the mandated Hermitian power forms.

**Q4. ARE THE THREE RELAXATION TIMES GENUINELY FREE?** ⭐ Confirm no wording now privileges any relation
among `τ_A`, `τ_V`, `τ_X`, and that no check or acceptance value silently requires one. ⚠ The single
supplied acceptance value is known to live in the equal-time slice — verify the directive **says so
plainly** rather than leaving the general case looking covered.

## Also check

**6. ⭐ IS THE INSTABILITY CHANNEL STILL OPEN?** A growing root must remain first-class. ⚠ Confirm no
diagnostic, alone or combined, can reject one, and that each that could is explicitly scope-bounded.
**7. ⭐ IS EVERY CHECK ABLE TO FAIL?** The directive states per-check what wrong derivation it catches —
⚠ verify those claims. **8. ⭐ IS THE PHYSICS CLOSED?** Count equations against unknowns explicitly.
**9. Do the two directives specify the same physics?** **10. ⛔ Any ill-posed or fabrication-forcing task?**

## ⛔ CLEARED BY INDEPENDENT LEGS — do not re-open

Scope boundary · `B1`'s mass balance · the `A/B/C` split · §1b's branch prescription **including its
upper-half-plane extension** · §3b's virtual-displacement rule · the **supplied affinity `𝒜`**
(independently derived and confirmed by two separate parties) · `B8` controls B/C/D · closure count.

⚠ Unless a rev-9 change newly broke one — then say which change did.

## Do not read

- Anything named `PREREGISTERED`/`PREREG`.
- ⛔⛔ `research/pde_ledger_v3/_scratch/` — it holds prior reviewers' reports **and the authoring transcript
  for this revision**, including a derived candidate answer to the question Q1 polices. Reading it would
  both replace your judgement with another's and destroy your ability to run Q1 honestly.
- ⛔ `research/pde_ledger_v3/steps/S11b_HANDOFF.md`, the other `S11bB_rev*_fix_directive.md` files, and the
  **git history** of the directive files — all record the reasoning behind these changes.
- `research/pde_audit/`.

⭐ You may read `steps/S11bA_interface_response.md`, `V3_STEP_PLAN.md` and `reduction/` to judge
completeness. Those are barred from the *builders*, not from you.

## Required method

⚠ Not document-only: Q2–Q4 require computation; use SymPy/numerics freely. Q1 and Q5 require reading the
whole directive and, for Q5, the filesystem. Form your own view of what this step must establish **before**
judging whether the directive establishes it.
⚠ Where you believe a task will produce a specific answer, ⛔ do not report that answer.

## Physics filter

Report a finding only if it catches a way the physics could be wrong, or a way the dual-engine check could
be defeated. ⭐ Q1 and Q5 findings count under the second clause.

## Output

⭐ Answer Q1–Q5 explicitly and first, each with a one-word verdict, and state whether you **verified by
computation** or only read. Then further findings, most severe first, each marked NEW-IN-REV-10 or
PRE-EXISTING. If nothing survives the filter, say so plainly. Under 60 lines.
