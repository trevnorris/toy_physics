# Independent review — the S9 export-chain build directive

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S9_export_chain_rebuild_directive.md`

It was written by the orchestrator and **nothing has been built from it yet.** This review decides whether
it gets built. A Codex builder will execute it as written.

## Why this review exists

The directive is the one artifact the builder reads, so an error in it is not caught by anything
downstream. Measured in this project over the last two days: the orchestrator's W0 directive carried six
verified defects; two orchestrator-written decision lists stated things about the engines that were false
and a builder correctly refused one of them; the S11 pair of engine directives went two review rounds and
seven blocking defects. **The base rate of defects in orchestrator-written directives here is high.**

## Method — read the repository FIRST, the directive SECOND

Form your own view of what S9's two engines currently compute and emit **before** you open the directive.
Reading the directive first anchors you to its framing, and its framing is the thing under test.

- `/var/projects/toy_physics/research/pde_ledger_v3/S9_REWRITE_PLAN.md` — the governing plan the directive
  claims to apply. Decisions **D1–D12** are settled there and a builder does not re-choose them.
- `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.py`
- `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S9_light_requires_shear_mathematica_audit.wl`
- the two committed outputs under `scripts/out/` and `mathematica/out/`
- `/var/projects/toy_physics/CLAUDE.md` — how this project works, sixteen rules

The registry directory `research/pde_ledger_v3/reduction/` was deleted at commit `fb29bba2`. Retrieve the
old comparison config with `git show 67dd3ce2:research/pde_ledger_v3/reduction/checks_S9.yaml` if you need
it. **Do not restore it.**

**Do not modify the repository. Do not build anything. Do not commit.**

## What to attack

**① The twelve-row rename table — is every pair real?** The directive asserts twelve
Wolfram-tag/SymPy-tag correspondences and instructs both engines to rename them to one shared name. Check
each against the committed outputs and the engine sources. A wrong pair renames two *different* objects to
one name, which silently converts a real disagreement into an apparent agreement — the worst failure this
architecture can have. Report per row: correct, wrong, or unverifiable.

**② The container-versus-component claim.** The directive asserts that for several pairs one engine emits
the object and the other emits a container holding it, and that this falls "almost entirely on the SymPy
side". **Determine, per row, which side is the container** and give the evidence. The directive forbids
fixing a Wolfram-side container (D12 restricts the `.wl` to name strings) and requires reporting it
instead. If Wolfram-side containers are more common than the directive claims, the instruction to rename
those rows produces a shared name over unlike payloads and must change.

**③ Is D1's export boundary well-defined against the actual script?** D1 says export "every object S9's
`MAIN` package emits", on the grounds that package membership is mechanical. Does the SymPy engine have a
`MAIN` package that a script can enumerate without judgement calls, or is that a property of the tag
strings only? What exactly would a builder export, and is that set determinate?

**④ Does the directive leak an expected result anywhere?** Rule 5 of `CLAUDE.md`: the spec says what to
compute, never what anything equals, is expected, or was measured. A builder that can read a target
iterates until it matches. Check every table, example and verification item — including whether naming an
object (`TRANSVERSE_SPEED_SQUARED`, `COEFFICIENT_DIMENSION_DIFFERENCE`) tells the builder something about
the answer that it should have had to compute.

**⑤ Is verification item 1 able to fail?** It requires building an old→new name map, applying it to the
committed baseline, and diffing. Construct the case where a real physics change survives that check. Is
"byte-identical payload under a renamed tag" strong enough to establish the restructure moved no physics —
and what class of change does it miss? Say plainly if the answer is that it misses nothing.

**⑥ Is anything here unbuildable, or already built?** The directive asks for an in-run exact round-trip
check on `S9_exports.py` (D3), a dimension emitted "in whatever shape the engine's own dimension solve
produces" (D4), and a knob extractor over annotated declarations (D6). Are those buildable against this
engine as written? Give the concrete blocker if not.

**⑦ Does the directive contradict the plan, or quietly extend it?** The plan is closed and says explicitly
that nothing may be added to it by the orchestrator or a builder. Name anything the directive adds,
narrows, or reinterprets — including the rename scheme in §4, which the plan does not spell out at this
level of detail. **A finding that a piece of this is unnecessary is worth as much as one that it is
insufficient.**

## Physics filter

> Report a finding only if it changes what the project computes, what it may claim, or what method it
> adopts.

Worth most: *"pair N is wrong and would fuse two different objects"* · *"the Wolfram engine is the
container for pair N, so that rename is unsound"* · *"D1's boundary is not determinate in this script"* ·
*"verification item 1 passes on a real physics change, here is the case"* · *"the directive supplies a
result the builder should have had to compute."*

**This project's failure mode is adding machinery, not omitting it.** Do not propose new checks, new
files, or new layers unless their absence lets a specific wrong answer through — and then say which one.

A leg returning "nothing survives the filter" is weak evidence on its own: state what you checked, what you
could not, and what would have had to be true for you to find something.

## Output

Rank most-severe first. For each finding: **the claim · the evidence (a quotation with `file:line`, or a
literal command and its output) · what must change in the directive.**

Give literal commands and literal output for every claim about what an engine does or does not emit. A
prose assertion that you checked something is not evidence in this project.
