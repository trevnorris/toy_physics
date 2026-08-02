# DIRECTIVE — the S11 TeX card

Write **one** file: `research/pde_ledger_v3/paper/steps/S11_stray_longitudinal.tex`, and add its
`\input` to the light part alongside the S9 and S10 cards. ⛔ Do not commit. ⛔ Do not edit any other file.

**Source of truth:** `research/pde_ledger_v3/steps/S11_stray_longitudinal.md`. The card is a faithful
rendering of that record — ⛔ do not add physics, do not soften it, do not "improve" a claim.

**Match the existing cards' structure and macros:** read `paper/steps/S10_two_transverse_photons.tex`
and `paper/steps/S9_light_requires_shear.tex` and follow their conventions exactly.
⚠ Check `paper/macros.tex` before using any `\stagefield`-style macro — **some are suppressed in the
default build**, and a disclosure placed in a suppressed field is invisible to a reader. Anything a
reader must see goes in visible body text.

---

## ⛔⛔ THE ONE WAY THIS CARD GOES WRONG

This step's honest content is unusually **negative**: it establishes a spectrum and then states three
things it **cannot** settle. ⛔ **The failure mode is a card that reads as a clean result with the limits
in small print at the bottom.** That is drift produced by trying to make a card read well, and it would
misrepresent the step.

⇒ **Requirements, not preferences:**

1. ⭐⭐ **The departure section LEADS with the unified statement** — *the second mode's entire physical
   status (radiative, observable, Lorentz-breaking, or none of the above) reduces to one unbuilt object,
   the brane–bulk interface coupling law* — **before** any table of what was computed. ⛔ It must not be
   a caveat appended after the spectrum.
2. ⭐ **"What S11 does NOT deliver" is a VISIBLE section**, not a footnote: bound-vs-radiating · that
   light's confinement is unconditional · whether the longitudinal is observable · any leakage rate.
3. ⭐ **Both corrections to the walk appear as corrections**, with what was wrong stated plainly: the
   invariant-count claim (and that the *method* was the error), and that the bulk-matching reading was
   overstated in **both** directions.
4. ⭐ `B_comp` is labelled **postulated** everywhere it appears, never derived, with its retirement
   condition named.
5. ⚠ **The known limits are reported**, including that two of six dimensional checks in the `.wl` and
   assertion 12 in the `.py` cannot fail and are not counted as coverage.
6. ⛔ The engines' `VERDICT: PASS` must **not** be presented as a verdict on the physics.

## Acceptance

- `cd research/pde_ledger_v3/paper && pdflatex` (twice) builds clean, **no new warnings**.
- Report the page count before and after.
- Every number, dimension vector and equation in the card traces to the step record. ⛔ If the record
  and your rendering disagree, **stop and report it** — do not reconcile.

## Report

Under 20 lines: where the departure statement sits relative to the spectrum, how you rendered the two
corrections, the page count delta, and anything in the record you could not render faithfully.
⛔ Do not summarise this directive back to me.
