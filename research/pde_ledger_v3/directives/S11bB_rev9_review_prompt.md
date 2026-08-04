# Independent physics review — S11b-B BUILD DIRECTIVES, REVISION 9 (before any build runs)

## Artifact

- `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11bB_wl_directive.md`
- `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11bB_py_directive.md`

Engine-specific header + a shared physics block, byte-identical in both
(`sha256 9bace451fc436eb4626b09d188a69a1f852fec6ccd437a2f4725c6e3d19e45c3`).

## Why this review exists

Two independent engines will each implement these directives, and their agreement is what certifies the
result. **The directive is the one artifact both engines share.** An error in it lands in both, they agree,
and dual-engine certifies wrong physics.

⚠ **Eight revisions have been rejected before any build ran**, and in most rounds every finding sat in
material that had just been changed. ⇒ ⭐⭐ **weight your attention toward what rev 9 changed.**

### What rev 9 changed

- **§0b** — the read-bar was rewritten from prefix globs to absolute paths; both **headers** were made to
  state the identical external-read restriction.
- **§2b** — separate relaxation times `τ_A`, `τ_V`, `τ_X` replace a single `τ`, propagated throughout.
- **§3b** — the **causality diagnostic was scoped** to a named set of objects; a **new two-port energy
  balance check resolved order-by-order** in the reciprocal coefficient `Λ_X⁰` replaces a check that was
  restricted to `Λ_X⁰ = 0`.
- **§3b/§3** — the divergence-invariant instruction was made plural.
- **B8** — Control F's acceptance criterion was made self-contained.

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

**Q2. Does the new ORDER-BY-ORDER two-port balance check actually cover the `Λ_X⁰ ≠ 0` regime, and can it
fail?** ⭐ Verify by computation. The previous check was restricted to `Λ_X⁰ = 0`, which meant the engine's
energy balance was compared against an independent standard **only where the interface does nothing**.
⚠ Judge: (a) does the new check genuinely constrain the derivation at `Λ_X⁰ ≠ 0`; (b) construct a wrong
derivation — a traction sign or factor error, a one-face omission, or eliminating the wrong variable
against `B1` — and confirm the check flags it; (c) does it classify any root or comment on the sign of any
imaginary part? ⛔ It must not.

**Q3. Are `τ_A`, `τ_V`, `τ_X` correctly and completely propagated?** ⭐ Verify: does `τ_A = τ_V = τ_X`
recover the single-`τ` specification **exactly**? Is any place still assuming a shared relaxation time
implicitly — the face response, the power identity, the Onsager array, the validity conditions, the
controls? ⚠ An implicit shared `τ` surviving anywhere would silently re-decide the question §2b asks.

**Q4. Is the SCOPED causality diagnostic still able to catch what it must?** It was narrowed so it no longer
fires on the Hermitian power forms the directive itself requires. ⭐ Confirm both halves: (a) it no longer
produces a false FAIL on the mandated passivity computation; (b) it **still** catches an advanced or
conjugated kernel appearing in the equations of motion, closure, face response or determinant.

**Q5. Is the read-bar now actually complete?** ⭐ Check it against the real filesystem: is there any file
reachable from the engines' working directory that contains this step's reasoning, prior reviews, or
another revision of the directive, and is **not** covered? ⚠ Answer for **both** engines and confirm the
bar is symmetric.

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
computation** or only read. Then further findings, most severe first, each marked NEW-IN-REV-9 or
PRE-EXISTING. If nothing survives the filter, say so plainly. Under 60 lines.
