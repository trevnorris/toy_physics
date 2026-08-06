# Bring the S10 ledger card up to the rebuilt record

**File:** `research/pde_ledger_v3/paper/steps/S10_two_transverse_photons.tex` — it already exists.
**Update it in place. Do not commit. Do not modify any other file except as D6 requires.**

## The single source of truth

`research/pde_ledger_v3/steps/S10_two_transverse_photons.md`, as it now stands. **Read it in full before
writing a line.** ⛔ Do not take anything from the card's current text as still true, and ⛔ do not consult
the engines, the harness output, or any other document for a number — every figure in the card must come
from that record.

⚠ **The card is currently wrong in ways a reader cannot detect.** Its Verification paragraph names two
engine files that **do not exist** (`S10_two_transverse_photons_*_audit`; the engines are
`S10_brane_mode_spectrum_*_audit`), describes a **quarantine procedure that has been retired**, and cites
gate results from an instrument that has since been rebuilt. Its Equations paragraph states that the
`D=3` reduction to `|∇×u|²` is *"computed, with residual zero, rather than assumed"* — which the record
now shows was true of **one** engine only until the `§Q7` repair.

## What the card must become

⭐ **Match the existing card format exactly** — `\stagefield`, `\claimstatus`, `\resultanchor`,
`\StatusText`, `\StageFile`, the table style. Compare against `paper/steps/S9_light_requires_shear.tex`
and `S11_stray_longitudinal.tex`. ⛔ Do not invent a new format and ⛔ do not restructure the sections that
are still correct.

**D1 — keep, unchanged in substance:** the three-move argument, the mode-count table, the
not-a-codimension paragraph, the general-`D` dimensions and the identity that explains the S9 loose end,
the what-is-new table, inputs, regime, and the departure paragraph. These survive the rebuild intact.

**D2 — the mode-count table gains a column or a companion paragraph** for the **transverse nullity**,
which the engines compute separately from the plain nullity, and a statement of the `§Q8` stratum result:
which packages have an allowed stratum, which do not, and what happens on the one that does.
⚠ The record is explicit that the stratum comparison is a **hand** pairing and not a harness result. The
card must not present it as machine-verified.

**D3 — a new prior-art paragraph.** The record's `PRIOR ART` section, condensed. ⛔ It must carry the
over-claim refusal: this step linearises about `v₀ = 0`, which is MacCullagh's own regime, so agreement
with him is agreement inside his domain of validity and is **not** evidence the medium differs from a
nineteenth-century elastic aether. ⭐ Then what is genuinely ours, as the record lists it.

**D4 — rewrite the Verification paragraph from scratch.** Correct engine paths, the harness counter block,
the action comparison, and the `§Q7` history — *including* that the SymPy `§Q7` was a check that could
not fail until it was repaired, and that a zero `§Q7` residual is a claim about the **form** and never
about the **term**.

**D5 — a new limitations paragraph**, from the record's `WHAT THIS STEP STILL DOES NOT ESTABLISH`. ⛔ Do
not summarise it into optimism. ⭐ In particular the card must state, in the paper's own voice: that the
`D=3` count of two rests on an in-plane/out-of-plane decoupling **nobody has computed**; that
`D_brane = 3` is postulated; and that cross-engine agreement did **not** improve when `§Q7` was repaired.

**D6 —** confirm the card's `\input` is already present in `paper/parts/part01_light.tex`; add it only if
it is missing. Then **build the paper** and report that it builds clean, with the page count.
⚠ `latexmk` is **not installed**; the working recipe is
`pdflatex -interaction=nonstopmode pde_ledger_v3.tex` run **twice**, from `paper/` — see
`paper/README.md`.

## ⛔ The rule that governs every number

**Every figure in the card must appear in the record, with the same value.** ⛔ Do not round, do not
recompute, do not "improve" a number, and ⛔ do not carry forward any number from the card's current text
without finding it in the record first. ⭐ If the record does not give a number the card seems to need,
**leave it out and report that** — do not go looking for it elsewhere.

## Acceptance — run these and paste literal stdout

1. **Every numeric figure in your new card, in a list, each with the record line number it came from.**
   Any figure you cannot source that way is a failure — report it rather than sourcing it elsewhere.
2. The paper build: the command, that it exits clean, and the page count.
3. `grep -n "two_transverse_photons_.*audit\|quarantin" paper/steps/S10_two_transverse_photons.tex` —
   paste the output. It should be empty.
4. Old and new line counts.
5. `git status --short`.

## Constraints

- **Do not launch Mathematica or `wolframscript`.**
- Do not run either engine and do not run the harness. Everything you need is in the record.
- No `git add`, no `git commit`, no other git write.

## Report back — under 20 lines

1. `D1`–`D6`: done / not.
2. The acceptance output.
3. Anything in the record you found unclear, self-contradictory, or impossible to render — ⭐ **this is
   wanted**, it is a review of the record as much as a build of the card.
