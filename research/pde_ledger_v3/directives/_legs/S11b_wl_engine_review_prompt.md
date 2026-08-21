# Independent physics review — S11b blind Wolfram engine

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11b_interface_coupling_law_mathematica_audit.wl`

This is a **Mathematica script** (`wolframscript`). It is the *blind* engine of a dual-engine step: it
imports nothing and re-derives the S11b interface coupling law from the shared physics spec alone, emitting
one tag per named object as `WL_S11B_<QUANTITY>: <value>` in `InputForm`. There is deliberately **no
VERDICT/PASS/FAIL line** — a boolean-valued test is emitted as the symbolic object the test returned.

## What to check
The `.wl` re-derives the interface coupling law from the spec ALONE and emits **computed** objects for:
bulk face response / acoustic impedance and its frequency regimes, added (entrained) mass, the
mass-balance constraint, the equations of motion, the thickness and compressional responses, passivity /
reciprocity / admissibility, the longitudinal mode fate + grazing threshold, the breathing-mode stability
slice, the transverse mode, parity, kernel orientation (the memory-kernel/causal sign), and
dimensions / homogeneity. Your job: **does each emitted object follow from the physics, or is any of them a
hand-typed answer, a tautological residual, or an object whose value does not depend on the derivation?**

## What you are handed
- The artifact above (read it, run it, ablate copies of it).
- The **sole physics authority** to derive against:
  `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md`.
  Every equation, ansatz, supplied law, derivation route (§§0–13) and tag obligation (§10) is there.
- `/var/projects/toy_physics/CLAUDE.md` (how we work).
- ⛔ You are **NOT** handed, and must **NOT** seek out, the sibling SymPy engine
  (`scripts/S11b_interface_coupling_law_sympy_audit.py`), its `scripts/S11b_exports.py`, or its `.out`.
  The two engines are independent by construction; reading the sibling's output would turn your
  independent derivation into a confirmation. If you find the `.wl` reads any of them (or any repository
  file) at run time, that is itself a finding — the blind engine must import nothing.

## Required method — THIS IS A SCRIPT
Derive independently. **Write your own derivation script BEFORE opening the artifact in detail, and save
BOTH the script and its literal stdout to named absolute paths** — report those paths. A prose
re-derivation is worth nothing; without a runnable script and its stdout, your derivation claims will be
discarded. Two hand-derivations agreeing is not evidence.

Ablate every load-bearing check and report its **literal** output; code-reading alone has repeatedly missed
real defects here. Probe for: a value verified using the predicate/definition that produced it
(`c := sqrt(x)` then emitting `c^2 - x`); a conclusion emitted as an unconditional literal; a check whose
expected value lives inside the artifact it checks; and any `assert`/`Abort` that PRECEDES the value it
guards (report it — a perturbation strong enough to flip such a guard kills the run, so you would see only
pass-or-crash).

⛔⛔ **A FORM ABLATION IS MANDATORY, not optional — it is the only thing that has ever caught the worst
defect here.** Change the **structure** of a load-bearing object — e.g. flip a sign *and* an off-diagonal
term together, or collapse two independent symbols into one — then re-run and report the **literal** diff.
A COEFFICIENT rescale tests only arithmetic; only a FORM change tests physics, because a rescale never
leaves the family. If an emitted tag is byte-identical before and after a structural change to the object
it supposedly reports, that tag is not computed from the derivation — that is the finding.

⇒ Ask of every emitted object: **WHICH LINE COMPUTED THIS?** Give the line number, or report it as
uncomputed.

### The withheld acceptance criteria — DERIVE THESE YOURSELF, do not look for them stated
The spec deliberately withholds the numeric/relational answers to several checks so the engine cannot fit
to them. Independently derive, from the spec's setup, and compare to what the `.wl` emits:
- the `μ_s = 0` reduction of the compressional / face-closure response (B2c's acceptance obligation): the
  relation it forces between the pressure-channel and affinity-channel interface coefficients;
- the transverse mode on a **uniform** background (B6) — is its coupling computed or asserted, and what is
  its dispersion / growth rate;
- the breathing-mode stability quantity on the `k = 0`, impermeable slice (B4);
- the count of independent quadratic invariants in the closed energy basis (B5/§5) — report the number the
  engine reaches and how it reaches it. (⚠ This count is a known open question across engines; report the
  `.wl`'s number and its construction, do not assume a value.)
Report each as "engine emits X; my independent derivation gives Y" with your script path.

### ONE-SIDED CORRUPTION — the only test that two routes are INDEPENDENT
Where the `.wl` claims two independent routes to one object and emits their residual, a zero residual proves
nothing on its own. **Corrupt ONE route at a time and report which objects moved.** If breaking route A
also moves route B's object, they were never independent and the residual is decoration. ⛔ Do NOT verify a
supplied object against the artifact's own identity — that is circular.

## Physics filter
Report a finding only if it catches a way the **physics** could be wrong. Do not report "the script would be
wrong on a different input" or style. A leg that returns "nothing survives the filter" is weak evidence —
if you find nothing, show the FORM ablations you ran and their literal diffs so the null result is measured.

## Ablation sandbox — MATHEMATICA, TWO-SEAT LICENCE (binding, in this prompt so both legs run identically)
⛔ Copy the artifact to `/tmp` and ablate the **COPY**. ⛔ Never modify the working tree.
⛔ Wrap EVERY kernel run in `timeout 600`. A 600s hit is a FAILED ablation — report it and move on.
⛔ NEVER raise the timeout. ⛔ NEVER run more than one kernel at a time (the licence has TWO seats and the
orchestrator may be holding one).
⛔ Kill any kernel whose resident memory passes **6 GB**. If a kernel is killed, run `free -h` FIRST before
anything else — an orphaned kernel leaks memory and will OOM unrelated jobs; then
`ps -eo pid,rss,pcpu,etime,comm --sort=-rss | head` to find and kill the orphan.
⭐ Save every ablation script AND its literal stdout to named absolute paths, and report those paths.

## Report
Under ~30 lines. For every finding: the emitted tag, the line that (should) compute it, your independent
result, the literal ablation diff, and why it is a physics error (not an input-sensitivity nit). Keep
findings separately numbered.
