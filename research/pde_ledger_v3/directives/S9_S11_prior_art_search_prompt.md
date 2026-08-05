# Prior-art search — is any of this known, or is it new?

You are doing a **literature and prior-art search**, ⛔ not a review. Nobody is asking you to check whether
this work is correct. The question is narrower and purely factual:

> ⭐ **Which parts of this mathematics already exist in the published literature, and which parts appear to
> be new — or at least, which you cannot find anywhere?**

## Read these first, in full

**The three steps in question (light sector, S9 → S11):**
- `/var/projects/toy_physics/research/pde_ledger_v3/steps/S9_light_requires_shear.md`
- `/var/projects/toy_physics/research/pde_ledger_v3/steps/S10_two_transverse_photons.md`
- `/var/projects/toy_physics/research/pde_ledger_v3/steps/S11_stray_longitudinal.md`
- `/var/projects/toy_physics/research/pde_ledger_v3/steps/S11bA_interface_response.md`
- `/var/projects/toy_physics/research/pde_ledger_v3/steps/S11bB_interface_assembly.md`

**Context for what the whole thing is trying to be:**
- `/var/projects/toy_physics/research/pde_ledger_v3/CHARTER.md`
- `/var/projects/toy_physics/docs/model_map.md`

⚠ **These documents are mid-revision and some contain known errors.** ⛔ Do not report on their quality,
their internal consistency, or their defects — that is being handled separately and is not your task.
⭐ Read them to extract **the mathematical content**, then go looking for that content in the literature.

---

## What to do

### 1 · State the mathematics in neutral, standard language

⭐ **First, and this determines whether the rest of the search is any good:** strip the project's own
vocabulary and restate what is actually being done in the language a physicist or mathematician outside
this project would use. The project says "brane", "defect", "throat", "branon", "drum head" — ⭐ translate
those into whatever the standard terms are, if standard terms exist.

⛔ Do not accept the project's framing as the right description. If what is really going on is, say,
ordinary linear elasticity in `D` dimensions with a particular stiffness tensor, ⭐ **say that plainly.**

### 2 · Search for each distinguishable piece, separately

⭐ **Search piece by piece, ⛔ not for the whole thing at once.** A combination can be novel while every
component is textbook, and that distinction is the entire point of this exercise. At minimum, search for:

- an **elastic Lagrangian whose stiffness is built only from the antisymmetric derivative**
  (`(∂_i u_j − ∂_j u_i)²`), and what mode content that gives in `D` dimensions;
- **elastic / solid-state analogues of electromagnetism** — media whose transverse modes are identified
  with photons, and how far those analogies are normally pushed;
- **brane-world models with an elastic or fluid brane**, and **branon** literature (the brane's own
  transverse fluctuation as a field);
- **elastic and hydrodynamic analogue-gravity** models — acoustic metrics, analogue horizons, and any
  programme that tries to get **both** gravity-like and electromagnetism-like behaviour from **one** medium;
- **topological defects in ordered media** carrying a conserved, sign-valued charge — and whether anyone
  derives a **Coulomb-like interaction from the energetics of the deformation** rather than by introducing
  a gauge field;
- **superfluid / Bose-condensate models of emergent electrodynamics** (the Volovik programme and anything
  adjacent);
- the **interface / matching-condition** mathematics in the later steps — an impedance-style or
  Onsager-style response matrix at a boundary between two such media.

### 3 · For each piece, report where it sits

Use exactly these labels:

| label | meaning |
|---|---|
| **STANDARD** | textbook or long-established; give the standard name and a canonical reference |
| **KNOWN** | appears in the research literature; give specific papers, authors, years |
| **ADJACENT** | something close exists but differs in a way you can state precisely |
| **NOT FOUND** | you searched and could not find it — ⭐ say **what** you searched for |

⛔ **`NOT FOUND` means "I could not find it", ⛔ NOT "it does not exist."** ⭐ Say so plainly, and say how
hard you looked. ⚠ A negative result from a search is weak evidence, and overstating it would be worse than
useless here.

### 4 · The honest bottom line

⭐ Answer directly, and ⛔ do not soften it in either direction:

- Which parts are **completely standard** and would be recognised immediately by someone in the field?
- Which parts are **known but obscure** — real literature exists, but it is not widely known?
- Is the **combination** — one medium, one action, aiming at both a gravity-like and an EM-like sector,
  with charge as a topological defect — something that has been **tried before**? By whom? ⭐ And if it has,
  **how did it go**, and where did it run into trouble?
- Is there anything here you genuinely **could not find any precedent for**?

⚠ ⭐ **If the honest answer is "essentially all of this is known and the novel part is small", say that.**
That is a useful answer and it is the one this project would rather have early than late. ⛔ Do not inflate
novelty to be encouraging.

⚠ Equally, ⛔ do not dismiss it by pattern-matching to "another aether model" without checking — ⭐ say
specifically **which** prior work it resembles and **where it differs**.

---

## Report format — ⛔ under 60 lines

1. **The mathematics, restated in standard language** — a short paragraph, no project vocabulary.
2. **Piece-by-piece table**: piece → `STANDARD`/`KNOWN`/`ADJACENT`/`NOT FOUND` → references or what you
   searched.
3. **The combination** — tried before or not, by whom, how it went.
4. **What you could not find precedent for**, with how hard you looked.
5. ⭐ **Anything you think this project should read**, most useful first — papers, books, review articles.
   ⭐ This is the most valuable part of your answer.

⛔ Write nothing to the repository. ⛔ Make no edits. Report only.
