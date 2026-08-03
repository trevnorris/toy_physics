# S11b-B — WRITE REVISION 7 OF THE BUILD DIRECTIVE

⛔⛔ **THIS IS A DOCUMENT TASK. ⛔ Do NOT write, run, or modify any `.wl` or `.py` script.**

## Your role has changed, and the reason matters

**You are the AUTHOR of rev 7.** The previous six revisions were authored by someone else and each was
rejected by independent review before any build ran. ⚠⚠ **In the last four rounds, every finding was located
in material that had just been changed to fix the previous round.**

⭐⭐ **The author's characteristic failure has been diagnosed and you must not repeat it:**
**every time a binding kinematic relation was DESCRIBED IN PROSE instead of WRITTEN AS AN EQUATION, the
in-plane divergence term `∇·u` fell out of it.** That has now happened **four times**:
- the original whole-interface spec (both engines would have found a free longitudinal and agreed),
- twice more in `B1`, which only stopped failing once the exact balance law was **written out as an
  equation**,
- and again in rev 6's virtual-displacement rule, where the prose *"conserves the slab's material content
  per unit in-plane area, `δΣ = 0`"* read literally gives `δθ = −δe_W` with **no `∇·δu`**.

⇒ ⭐⭐⭐ **RULE FOR REV 7: every binding kinematic relation MUST appear as an explicit equation.** ⛔ If you
find yourself writing a sentence that a reader must translate into an equation, write the equation instead.
⭐ Where a term has been dropped before, add an explicit guard naming what must not be absent.

## Artifacts

**Edit only:** `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11bB_SHARED_PHYSICS.md`

Then **reassemble both directives** so the shared block is byte-identical in each:
```
cd /var/projects/toy_physics/research/pde_ledger_v3/directives
cat S11bB_wl_header.md S11bB_SHARED_PHYSICS.md > S11bB_wl_directive.md
cat S11bB_py_header.md S11bB_SHARED_PHYSICS.md > S11bB_py_directive.md
```
⛔ **Do not edit the headers.** Both legs have cleared them repeatedly.
⭐ Verify byte-identity of the shared block in both outputs and **report the `sha256` of
`S11bB_SHARED_PHYSICS.md`**.

## The findings you are fixing — ⭐ READ THE REVIEW IN FULL

`/var/projects/toy_physics/research/pde_ledger_v3/_scratch/S11bB_dir_rev6_codex.txt`

⚠ **Read the FINDINGS section at the end of that file** — it is a long transcript and the findings are at
the bottom. Fix **every** finding that survives your own judgement. ⭐ **If you judge a finding to be wrong,
say so and say why — do not silently skip it.**

⚠⚠ **A SECOND, INDEPENDENT REVIEW OF REV 6 IS STILL RUNNING** and is not available to you. ⇒ ⛔ **Do not
assume the findings above are the complete set**, and ⛔ **do not write anything that presents rev 7 as
final or fully reviewed.** Its findings will be applied in the following round.

## ⭐ Derive the fixes; ⛔ do not transcribe them

The reviews describe defects and sometimes gesture at remedies. ⭐ **Work out the correct physics yourself**
and verify it — **use SymPy/numerics freely.** ⛔ Adopting a suggested remedy without checking it is how a
previous revision introduced a *precise and wrong* branch rule that would have flipped the deliverable's
sign.

⭐ In particular, **the virtual-displacement relation must be derived and written as an equation**, with the
`Σ` it refers to unambiguously defined (⚠ Eulerian versus material is exactly what went wrong).

## ⛔⛔ CONSTRAINTS — violating any of these is worse than leaving a finding unfixed

1. ⛔ **MINIMAL AND SURGICAL.** Fix the findings; change nothing else. ⚠ Wholesale rewrites are what bred
   each round's new defects.
2. ⛔ **NEVER SUPPLY AN EXPECTED RESULT OR ITS REASON.** The engines that will implement this must be able to
   return **any** answer, including the negation of what anyone expects. ⚠ Rev 5 was caught asserting where
   B5's roots lie; that prejudged the deliverable. ⛔ Do not state or hint at what the longitudinal mode
   will do, whether the transverse coupling vanishes, or what sign any imaginary part takes.
3. ⛔⛔ **KEEP THE FALSIFICATION CHANNEL OPEN.** A **growing** root is an admissible, first-class outcome.
   ⛔ Do not add any stability assumption, positivity requirement, or acceptance gate that would suppress
   one. ⚠ Any diagnostic that could reject a growing root **must** have its scope explicitly bounded to a
   sub-case where growth would necessarily be a derivation error.
4. ⛔ **EVERY CHECK MUST BE ABLE TO FAIL.** ⚠ A previous revision shipped a check that enumerated loss
   channels *from the equations its own prescribed route produced* — an identity. ⭐ For each check you add
   or modify, state in one line **what wrong derivation it would catch.** ⛔ If you cannot, delete it rather
   than ship it.
5. ⛔ **DO NOT RE-OPEN what independent legs have cleared:** the scope boundary; header symmetry and the
   `reduction/` bar; `B1`'s supplied mass balance; the `A/B/C` step split; `B8` controls B/C/D as form cuts;
   §1b's branch prescription **as such** (⚠ verified computationally by three independent parties —
   ⭐ but *do* fix any defect a rev-6 finding identifies in how it is stated).
6. ⛔ **TWO ENGINES MUST NOT BE ABLE TO DIVERGE.** For every instruction you write, ask whether two competent
   independent implementers could comply and still produce different equations. ⚠ That is the failure mode
   this whole document exists to prevent.
7. ⛔ **No new scope.** Do not add tasks beyond what a finding requires. ⚠ The step is already large.

## Output

1. The edited `S11bB_SHARED_PHYSICS.md` and both reassembled directives.
2. A report, **under 60 lines**, listing each finding and **what you changed for it** — or why you judged it
   wrong. ⭐ Include the `sha256`, and for every check you added or modified, the one-line statement of what
   wrong derivation it catches.

⛔ Then stop and exit. ⛔ Do not write a script, do not run a build, do not write a summary document.
