# Independent physics review

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/paper/steps/S11_stray_longitudinal.tex`
(working-tree, uncommitted — inspect with `git status` / `git diff`)

## What to check

**FIDELITY-ONLY review of a TeX card against its step record.** The record is the **source of truth**:
`/var/projects/toy_physics/research/pde_ledger_v3/steps/S11_stray_longitudinal.md`.

Report: **overclaims · dropped caveats · false artifacts · number mismatches.**

⛔⛔ **THE SPECIFIC RISK THIS CARD HAS.** Its honest content is unusually **negative** — it establishes a
spectrum and then states three things it **cannot** settle. ⇒ **The failure mode is a card that reads as
a clean result with the limits in small print.** That is drift produced by making a card read well, and
it would misrepresent the step. Check these specifically:

1. ⭐⭐ **Does the departure section LEAD** with the unified statement — *the stray mode's entire physical
   status (radiative, observable, Lorentz-breaking, or none) reduces to ONE unbuilt object, the
   brane–bulk interface coupling law* — **before** any table of what was computed? ⛔ Or is it appended
   as a caveat after the spectrum?
2. ⭐ Is **"what S11 does NOT deliver" a VISIBLE section**, not a footnote, naming: bound-vs-radiating ·
   that light's confinement is unconditional · whether the longitudinal is observable · any leakage rate?
3. ⭐ Are **BOTH corrections to the walk rendered AS corrections** — the invariant-count claim
   (**including that the METHOD was the error**, not just the arithmetic), and that the bulk-matching
   reading was overstated in **BOTH** directions?
4. Is `B_comp` labelled **postulated** everywhere it appears, ⛔ never derived, with its retirement
   condition named?
5. Are the **known limits** reported — including that two of six dimensional checks in the `.wl` and
   assertion 12 in the `.py` **cannot fail** and are **not counted as coverage**?
6. Is the engines' `VERDICT: PASS` kept as an **internal-consistency** statement and ⛔ **not** presented
   as a verdict on the physics?
7. ⚠ **MACRO HAZARD.** Read `/var/projects/toy_physics/research/pde_ledger_v3/paper/macros.tex`. Some
   `\stagefield`-style macros are **SUPPRESSED in the default build**. Anything a reader must see that
   sits in a suppressed field is **invisible in the PDF** — that is a real finding, and it has happened
   before on a gap disclosure.

Also: confirm the paper still builds clean (`pdflatex` twice) and report the page count.

## Do not read

- `/var/projects/toy_physics/research/pde_ledger_v3/_scratch/S11_tex_directive.md` — the build
  instructions. ⛔ Judge the card against the **record**, not against what the builder was told. A card
  can satisfy its directive and still misrepresent the record, and that is exactly what this leg is for.
- any prior review report or reviewer output for this step.

## Required method

Read the record **first**, form your own view of what the step establishes and what it does not, and
**then** read the card. ⛔ Do not read them in the other order — reading the card first anchors you to
its framing, which is the thing under test.

Quote the card and the record side by side for every finding.

## Physics filter

**"Report a finding only if it catches a way the physics could be MISREPRESENTED to a reader.**
⛔ Do not report LaTeX style, wording preference, formatting, or typography." ⚠ A short report with one
real finding beats a long one. If the card is faithful, say so plainly and stop.

## Quarantine rule

Not applicable — nothing is quarantined for this review.

## Output format

```
VERDICT: CLEAN | FINDINGS
OVERCLAIMS        (list, or "none")
DROPPED CAVEATS   (list, or "none")
FALSE ARTIFACTS   (list, or "none")
NUMBER MISMATCHES (list, or "none")
FRAMING           (items 1-3: does the negative content lead, or is it buried?)
SUPPRESSED-FIELD  (item 7: anything a reader must see that will not render)
BUILD             (clean? page count?)
WHAT IS FAITHFUL  (brief)
```
