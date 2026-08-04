# S11b-B — WRITE REVISION 11 OF THE BUILD DIRECTIVE

⛔⛔ **DOCUMENT TASK. ⛔ Do NOT write, run, or modify any `.wl` or `.py` script.**

You authored revisions 7–10. ⭐ **You are the author again for rev 11.**

## Artifacts

**Edit:** `directives/S11bB_SHARED_PHYSICS.md` and `S11bB_py_header.md` (title only — see below).
Reassemble both directives, verify byte-identity of the shared block, report the new `sha256`.

## The review to fold

`/var/projects/toy_physics/research/pde_ledger_v3/_scratch/S11bB_dir_rev10_agent.txt`

⚠ **F2 and F3 are NOT yours** — the read-bar is being fixed operationally at build time, by moving
artifacts out of the tree rather than by enumerating them in a prohibition. ⛔ **Do not edit §0b.**
⭐ Fix **F1, F4, F5, F6, F7** and the header title. Say so if you judge a finding wrong.
⚠ The other reviewer found nothing; ⛔ it has now cleared four revisions containing real defects — **do not
read its silence as corroboration.**

## ⭐⭐ F1 IS A REGRESSION, AND IT IS THE ONE THAT MATTERS

The causality diagnostic was restated last revision in terms of **pole location after removable
cancellations**. That restatement **created a hole**: the A-channel pole is displaced by the interface
coupling, so the test **passes a transposed kernel** in part of the declared parameter space and
**false-FAILs a correct derivation** in another part — and symbolically returns no verdict at all.

⚠⚠ **Both failure modes must be handled by whatever you write:**
- the **symbol-appearance** test it replaced caught all three channels, but false-FAILed a correct
  derivation that rationalizes a complex denominator by its conjugate;
- the **pole-location** test fixed that, but is indeterminate exactly where the pole is displaced.

⭐ **Derive a criterion that catches a transposed kernel in ALL THREE channels without firing on a correct
derivation.** ⛔ Do not adopt the reviewer's suggested direction without verifying it — ⚠ **the previous
restatement was adopted from a reviewer's suggestion unverified, and that is what produced this
regression.** ⭐ Verify with SymPy over the declared parameter ranges, ⛔ not on one convenient sign.
⚠ If no single criterion covers all three, ⭐ **say so plainly and state what each variant does and does not
catch** — an honestly bounded check beats a check credited with a catch it cannot make.

## ⭐⭐ F4 — THE ASYMMETRY THAT QUIETLY LOSES A FALSIFIER

Four places say a mistake manufactures **growth**; **none** names a mistake that manufactures spurious
**decay**, though each has a mirror that does. And only growth carries a procedural burden — diagnostics to
run first, a three-part report in B5.

⚠ The gates themselves are correctly scope-bounded and the channel is not closed. ⛔ **But the attention
budget is one-sided, and it biases an engine toward reporting decay** — which is the quiet way to lose a
falsification channel that survives every explicit gate.

⭐ **Balance it.** For each named growth-manufacturing mistake, name the mirror that manufactures spurious
decay, and make the reporting burden **symmetric in outcome** — whatever an engine must establish before
reporting growth, it must equally establish before reporting decay. ⛔ Do not do this by weakening the
growth requirements; ⭐ raise the decay side to match.

## The remaining findings

- **F5 — two structural presuppositions in §2b.** (i) "unconditional admissibility" versus any "**smaller**
  region that also imposes reciprocity" presupposes **strict nesting**; ⚠ coincidence of the two regions is
  one of the candidate outcomes, and an engine primed this way will report nesting rather than test for it.
  ⭐ Reword so coincidence, nesting and incomparability are all reachable. (ii) "identifying relaxation
  times before computing can give a false reciprocity relation" leans away from one candidate answer.
  ⭐ Reword to forbid the *procedure* without characterising the *outcome*.
- **F6 — "its corresponding velocity" is ambiguous** in the off-shell balance check: contracting the
  thickness equation with the face velocity instead of the thickness rate breaks a cancellation and yields a
  false FAIL. ⭐ Name each equation's conjugate rate explicitly. ⚠ Also correct the blind-spot summary: the
  check tests the traction's **sign, factor and face count**, ⛔ not its analytic structure — it is blind to
  a transposed kernel used identically on both sides.
- **F7 — the Onsager–Casimir step is under-specified for a mixed representation.** Parities are supplied
  but no rule for handling a mixed flux/force pair, so two engines can legitimately reach opposite signs.
  ⚠ `NOT_ESTABLISHED` is forbidden for this object, so an under-specified procedure is not survivable.
  ⭐ Either supply the conversion rule, or lift the `NOT_ESTABLISHED` prohibition for this part and say why.
  ⛔ Do not state what the resulting relation is.
- **Header title** — the PY header's **title** still says "registry insertion" while its body forbids
  registry work, and the WL header says nothing of the kind. ⭐ Fix the title only.

## ⛔⛔ CONSTRAINTS

1. ⛔ **MINIMAL AND SURGICAL.** ⛔ **Do not edit §0b.**
2. ⛔ **NEVER SUPPLY AN EXPECTED RESULT, ITS REASON, OR ITS SHAPE**, and ⛔ **never justify anything by
   reference to a previous revision or by naming a baseline case** — the engines cannot read either.
3. ⛔⛔ **KEEP THE FALSIFICATION CHANNEL OPEN**, and after F4, keep it **symmetric**.
4. ⛔ **EVERY CHECK MUST BE ABLE TO FAIL**, with a **TRUE** one-line statement of what it catches. ⚠ If a
   check's coverage is partial, ⭐ state the boundary in the directive.
5. ⛔ **DO NOT RE-OPEN** what independent legs cleared: scope boundary · `B1`'s mass balance · the `A/B/C`
   split · §1b's branch prescription and its upper-half-plane extension · §3b's virtual-displacement rule ·
   the supplied affinity `𝒜` · `B8` controls B/C/D · the closure count · the off-shell balance check's
   **structure** (fix only F6's ambiguity and its blind-spot wording).
6. ⛔ **TWO ENGINES MUST NOT BE ABLE TO DIVERGE.** 7. ⛔ **No new scope.**

## Output

The edited files, plus a report **under 60 lines**: each finding, what you changed, the `sha256`, and for
every check added/kept/modified a **true** one-line statement of what wrong derivation it catches — ⭐ with
its coverage boundary where partial. ⭐ For F1, state explicitly what you verified and over what parameter
range.

⛔ Then stop and exit.
