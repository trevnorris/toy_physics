# Author v5 of the S11c-d clean-condition packet

You are the **author** of version 5 of a two-file packet. v1–v3 were written by one author. v4 was written by a
previous Codex author session, from which you inherit authorship. Four rounds of independent review have produced
verified, accepted findings. Round 4 found three new defects in material that v4 changed (R4-2, R4-4, R4-6), so the
standard for v5 is that **every change is grounded in the code and the governing sources, not in a paraphrase**.

Your job is: **revise → ground → write the change log → stop.**
- ⛔ Do not launch, call or spawn any other AI, agent, reviewer or review process. Review is not your job.
- ⛔ Do not commit.
- Write only the files named under *Deliverables*.
- ⛔ Do not run the S11c engines or any Wolfram kernel.
- ⛔ Do not signal, resume or otherwise touch PID 4097233.

## The packet (v4, committed at `5693e861`)

All paths are under `/var/projects/toy_physics/research/pde_ledger_v3/directives/`.

- **A.** `S11c_d_clean_condition.md`
  - **Role:** a physics-bearing proposal/record. It states a symmetry-based "clean condition" (hypothesis H) under
    which light does not leak from the brane slab into the bulk at linear order. It proposes R-LEAK-1 and
    R-LEAK-2 and re-scopes S11c-d.
  - A is read by the user and by reviewers, **not** by the builder.
- **B.** `S11c_d_zinvariant_operator_blocks_directive.md`
  - **Role:** a builder directive. A Codex builder will receive **only B** and execute it.

## What you are handed

- **The review record:** `_measurements/S11c_d_clean_condition_review_disposition.md`. It has four rounds. Each
  finding is listed with its verification and the resolution owed. The **round-4 section (R4-1 … R4-10)** is what
  v5 must resolve. The round-1, round-2 and round-3 resolutions must not regress.
- **The round-4 reports**, with the legs' scripts and literal stdout:
  - `_measurements/S11c_d_clean_condition_review_r4_claude.txt`
  - `_measurements/S11c_d_clean_condition_review_r4_grok.txt`
  - `_measurements/S11c_d_clean_condition_review_scripts/r4_claude/` and `…/r4_grok/`
- **The v4 change log and author scripts:** `_measurements/S11c_d_clean_condition_v4_changelog.md` and
  `_measurements/S11c_d_clean_condition_v4_author_scripts/`. Round 4 found that those scripts are typed (R4-6).
- **The current grounding files:**
  - `_measurements/S11c_d_clean_condition.md` (for A)
  - `_measurements/S11c_d_zinvariant_operator_blocks_directive.md` (for B)
- **The governing sources**, all under `/var/projects/toy_physics/`:
  - `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md`
  - `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`
  - `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md`
  - `research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md`
  - `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py`
  - `research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py`
  - `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`
  - `research/pde_ledger_v3/scripts/S11c_b_exports.py`
  - `research/pde_ledger_v3/CHARTER.md`
  - `research/pde_ledger_v3/V3_STEP_PLAN.md`
  - `docs/toy_model_ontology_summary.md`
  - `docs/native_light_em_and_vortex_throat_interpretation.md`
  - `AGENTS.md`
  - `scripts/s11c_guarded_run.py`
- The rest of the repository is readable.

## Standards v5 must meet

1. **B names the object, not a recipe.**
   - Where B must list fields, inputs, outputs, rows or symbols, **enumerate them from the engines' code**, with
     file:line, not from a paraphrase. The previous author's lists were paraphrased, and each list missed or
     mis-typed something.
   - Where a requirement can be stated as a **symmetry** or a **domain** rather than as a list of operations,
     state it that way, and have the builder report how it was implemented, with file:line.
2. **B supplies verified physics as equations, labelled as supplied and unfalsifiable within that build.**
   - B never states what any entry equals, is expected to be, or was measured to be.
   - B contains no acceptance criterion that references an expected value.
   - B **must not encode H's prediction**: neither which blocks are expected to vanish, nor which may be nonzero.
     Printing every entry, including zeros, is the object.
   - Tag names name objects (field, row, block, control), never a value, sign or shape of a result.
3. **B keeps everything the review record lists as resolved.** That includes:
   - §0 execution safeguards (guarded runner, the PID 4097233 stop, one Wolfram kernel);
   - the declared freezes, listed first;
   - the import whitelist;
   - the blind Wolfram engine;
   - the raw-first comparator;
   - FORM controls K1 and K2 that enter at the stored energy;
   - the three clauses and the structural rule;
   - the Part 2 scope note that surfaces, and does not invent, the missing drain-flow specification;
   - build → run → report → stop.
4. **A stays a proposal.** Its claims stay conditional, scoped to linear order and to the stated premises. It
   claims nothing that Part 1's computation could not support.
5. **Every repo citation is grounded verbatim.**
   - Regenerate both grounding files from **commands**: mechanical lookups such as `sed -n` and `grep -n`, each
     shown with its literal output, run from the repo root.
   - Keep the `§P` prior-art table in A's grounding file verbatim.
   - Use fences longer than any fence inside the quoted text.
6. **Any physics claim you add or change carries a script and its literal stdout.** Save both under
   `_measurements/S11c_d_clean_condition_v5_author_scripts/`. A prose derivation is not evidence.
   - A script must **compute** what it prints from its inputs. A value that would survive deleting its inputs is
     typed, not computed (R4-6).
   - If a statement cannot be computed, present it as supplied or as bookkeeping, ⛔ not as evidence.
7. **The change log claims only what v5 demonstrably does.** For each finding, "resolved" requires the change and
   its location. Where resolution is partial, say so.

## Deliverables

- v5 of `S11c_d_clean_condition.md` and of `S11c_d_zinvariant_operator_blocks_directive.md`, edited in place.
  - Status line: **v5**, Codex-revised (v4–v5) from the orchestrator-written v1–v3, ⛔ not governing until
    review-cleared.
  - Name the v4 baseline commit `5693e861`.
- Both grounding files, regenerated for v5.
- `_measurements/S11c_d_clean_condition_v5_changelog.md`. For each of R4-1 … R4-10 it gives:
  - what changed;
  - where (file:line in v5);
  - the evidence path, if a script was used.

  It also states any finding you could not resolve, and why.
- Any scripts and stdout, under `_measurements/S11c_d_clean_condition_v5_author_scripts/`.

Stop after writing these.
