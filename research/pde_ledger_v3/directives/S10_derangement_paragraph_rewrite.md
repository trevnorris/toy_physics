# Rewrite one paragraph that has been wrong three times, and finish three propagations

**Do not commit.** Five items. ⭐ **`F1` is the reason a different author is writing this.**

---

## F1 — the single-row-derangement paragraph. ⛔ THE ORCHESTRATOR MUST NOT WRITE IT A FOURTH TIME.

**Where it lives, in two places that must agree:**
- `steps/S10_two_transverse_photons.md`, the table and the sentences around it (~`:495-508`)
- `reduction/checks_S10.yaml`, the parallel comment above the nine gradient declarations

**Its history, and this is why the instruction below is shaped as it is.** The paragraph describes what
happens when the nine gradient symbols inside **one row's payload** are deranged, under the committed
declaration. It has been written three times and been wrong three times, each time the same way: a leg
measured **one row under one permutation**, and the orchestrator wrote that branch as the general rule.

1. First draft: *"the corrupted row stays `NAMING_MISMATCH`, a coverage gap that stands out."*
2. Second: *"the corrupted row becomes `DISAGREE`."*
3. Third, the current text: a three-row table pinning one status per row.

**All three are branch generalisations.** Two legs have now measured it. The second measured **every**
declared `§Q7` row under **27 permutations** (11 named + 16 seeded random derangements) and computed each
payload's **stabiliser** — the number of the `9!` relabellings that leave it unchanged:

| declared row | stabiliser | statuses reached over the 27 |
|---|---|---|
| the four curl-form stiffness rows (`main`, `signflip`, `aniso`, `xcoef`) | 288 | AGREE 7 · **DISAGREE 16** · NAMING_MISMATCH 4 |
| `xform_divonly_d3_q7_stiffness` | 4320 | AGREE 7 · NAMING_MISMATCH 20 · **never DISAGREE** |
| `xform_divonly_d3_q7_difference`, `xform_fullgrad_d3_q7_difference` | 288 | AGREE 7 · DISAGREE 20 |
| **`xform_fullgrad_d3_q7_stiffness`** | **362880** | **AGREE 27/27** |
| the four `_q7_difference` rows whose payload is `0` | 362880 (vacuously) | AGREE 27/27 |

⇒ the current table's `main_d3_q7_stiffness → NAMING_MISMATCH` cell names the **minority** branch, 4
against 16.

⛔⛔ **And one sentence in the current text is simply false.** It says the six emitted
`ORDINARY_CURL_NORM` payloads are invariant under all `9!` relabellings, as `Σ g_ij²` is. **They are
not** — their stabiliser is **288**, and the relabelling `g3x2↔g3x3` changes every one of them.
⭐ **Exactly one emitted payload-bearing object is `9!`-invariant: `XFORM_FULLGRAD`'s stiffness density.**

⭐⭐ **A leg supplied the permutation-free statement the paragraph has been groping for, and it is short:**

> A derangement drawn from the **288-element joint stabiliser** produces no signal on **any** row.
> Outside it, a payload-bearing row moves — to **`DISAGREE`** when the corrupted payload cannot be
> re-reconciled by a per-row bijection, and to **`NAMING_MISMATCH`** when it can.

### What to write

⭐ **Build it on the STABILISER, not on a sample of permutations.** That is what makes a statement
permutation-free, and it is why the three previous drafts failed: each described a sample.

- ⭐ Use the leg's permutation-free statement above, or your own equivalent of it.
- ⭐ State that exactly one payload-bearing row, `xform_fullgrad_d3_q7_stiffness`, has stabiliser `9!` —
  ⇒ **no payload check can police it**, which is why the unit test pins its nine renames by hand.
- ⛔ **Delete the per-row status table.** A status reached under some permutations and not others is an
  instrument detail; it belongs in the stdout of
  `reduction/measurements/q7_map_declaration_ablation.py`, ⛔ not in a step record or a config comment.
- ⛔ **Write no sentence of the form "the corrupted row becomes X"** unless X holds for **every**
  permutation outside the stabiliser.
- ⛔ **Do not repeat the `ORDINARY_CURL_NORM` invariance claim in any form.**

⭐ **Keep** the surrounding claims, which are not in question: that one map is applied to every row; that
the 288-element group is what the payloads cannot see; that the declaration is sourced from the engine
definitions.

⚠ **Record the history in one sentence** — that this paragraph was wrong three times by generalising from
a single measured branch. It is the most useful thing on the page for whoever writes the next one.

---

## F2 — the paper card still claims what the record retracted

`paper/steps/S10_two_transverse_photons.tex` (~`:193-196`) still lists the **in-plane/out-of-plane split**
and the **`D`-sweep** under *"what is genuinely ours"*. The record no longer does: branons are filed as
**KNOWN** prior art, and the `D`-sweep is an evidence exercise over a standard identity. Bring the card
into line with `steps/S10_two_transverse_photons.md`.

⚠ The record's citation for the branon filing is a loose anchor. The KNOWN filing is at
`docs/medium_requirements_and_prior_art.md:142`. Point both at that.

## F3 — the paper card still carries the `D=3`-only figure

`paper/steps/S10_two_transverse_photons.tex` (~`:98-100`) states the generic anisotropic transverse
nullity as *"one at the propagating root"*. The engines emit **1 at `D=3` and 2 at `D=4`**, on the two
nonzero roots. The record's table was corrected; the card was not.

## F4 — "`§8` has no stratum scope token" is stated absolutely in two more places

Corrected once in the record's stratum section, and **not** in the record's limitations list
(~`:558`) nor in the card (~`:103-104`). It is true of the **as-built** grammar (`s10-as-built`) and
**false** of the repaired `directives/S10_SHARED_PHYSICS.md`, which adds `STRATUM<s>`. Scope both.

## F5 — six further corrections the same leg measured

Each is a one-line factual fix. ⭐ Verify each yourself before applying it.

1. `steps/S10…md` — *"of the 18 emitted `D=3` pairs, **14** carry information about the map"*. **13** do.
   The 14th is `XFORM_FULLGRAD`'s stiffness pair, which the same page calls invisible.
   ⚠ `checks_S10.yaml`'s own wording — *"on the 14 informative pairs, both the declared map and its
   transpose give residual 0"* — is **correct** and must not be changed; only the record's paraphrase
   drifted.
2. `steps/S10…md` — *"one row is INVISIBLE"*. **Five** of the twelve are: the `9!`-invariant one **plus
   the four whose payload is literally `0`**. Say which kind is meant.
3. `steps/S10…md` — *"the battery now requires the declaration to leave **no** `§Q7` row
   `NAMING_MISMATCH`"*. The code requires **all twelve to be `AGREE`**, which is strictly stronger, and
   its own comment explains why the weaker predicate is insufficient. ⇒ the record describes the
   predicate the code **rejects**.
4. `steps/S10…md` — the *"genuinely ours"* list still contains *"the longitudinal zero kept deliberately
   as the charge anchor"*. The prior-art register records that only as **NOT-FOUND**, and the same page
   says twelve lines later that NOT-FOUND is not originality. Apply that to this item too.
5. `SUBSTRATE_REQUIREMENTS.md` — the *"Honest sizing"* paragraph still says pass 1 produced **six**
   entries, *"three times what the file held before"*. It produced **seven** (nine total, two
   pre-existing), which the paragraph above it already states correctly.
6. `SUBSTRATE_REQUIREMENTS.md` `R-S8-03` — *"`B2` shuts only the route from a polar `P` to `μ_R`'s
   **magnitude**"*. `DEFECT_REGISTER.md` scopes `B2` as a route-failure for **deriving `μ_R` from a polar
   substructure**, which closes that route for the sign too. ⭐ The entry's **conclusion** — that the
   requirement is live via other routes — is unaffected; fix only the characterisation.

⚠ Items 5 and 6 are in `SUBSTRATE_REQUIREMENTS.md`, which the out-of-scope list below otherwise excludes.
**These two are in scope; nothing else in that file is.**

## F6 — sweep for the same three defects anywhere else

Each of `F2`–`F4` is a fix applied in one place and missed in another. **Grep both the record and the card
for every other occurrence** of: a novelty claim, the generic anisotropic nullity, and the stratum-token
claim. Report what you find, fixed or not.

---

## Out of scope

Either engine · either committed `.out` · `reduction/engine_output_checks.py` ·
`reduction/harness_ablation.py` · `DEFECT_REGISTER.md` · `steps/S9_*` · all of
`SUBSTRATE_REQUIREMENTS.md` **except** the two lines named in `F5.5` and `F5.6`.

## Acceptance — run these and paste literal stdout

1. Your rewritten `F1` paragraph, both copies, quoted in full. **Then read your own text back and list
   every status-per-row claim it makes.** ⛔ Any that does not name a permutation is a failure.
2. `grep -n` over the record and the card for the three swept claims (`F6`), with output.
3. The full unit suite from `reduction/`, and `harness_ablation.py`'s `ACCEPTANCE 19` line.
4. The paper build: two `pdflatex -interaction=nonstopmode pde_ledger_v3.tex` passes from `paper/`,
   the error count and the page count.
5. `git status --short`.

## Constraints

- No script over 600 seconds. **Do not launch Mathematica or `wolframscript`.**
- No `git add`, no `git commit`, no other git write.

## Report back — under 20 lines

1. `F1`–`F6`: done / not.
2. The acceptance output.
3. Anything you measured that contradicts this directive — ⭐ including if the leg's table above is
   itself wrong. **Check it before you rely on it.**
