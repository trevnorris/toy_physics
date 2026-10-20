# S12 comparator build, review round 0: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S12_comparator_build_review_r0_lookups.md`, beside this file
(generator `_scratch/s12/gen/cmp_build_r0_lookups.sh`). "D-n" cites a section of that file.

## The round

- **Artifact.** The S12 cross-engine comparator as first built, uncommitted at review on top of `3822b64e`:
  - `scripts/S12_cross_engine_comparator.py` and `scripts/test_S12_cross_engine_comparator.py`;
  - `directives/_measurements/S12_comparator_run_record.md` and `S12_comparator_builder_report.md`;
  - the output in `_scratch/s12/cmp_work/`, which is scratch and is not committed.

  The four repository files are committed unchanged with this file. All 22 baseline hashes matched after both legs,
  and the working tree was unchanged apart from this file and the lookups (D1).
- **Prompt and legs.** `_scratch/s12/cmp/review/s12_cmp_build_review_prompt_r0.md`, sha256 `6cbd54cc…` (D1). Both
  legs received it identically and ran in parallel, each under a 25 GiB pooled cap.
  - **Claude:** a fresh opus agent. Report `/var/projects/_review_scratch/s12-cmp-r0-claude/REVIEW.md`.
  - **Grok:** `grok-4.7`, effort high. Report `/var/projects/_review_scratch/s12-cmp-r0-grok/grok_stdout.txt`.
- **Verdicts.** Claude: NEEDS REVISION, 4 findings. Grok: CLEAR, no findings (D2).

## Findings

**F1. The same derivative, printed in the two engines' notations, is read as two different objects.** (Claude F1.)
**CONFIRMED.**
- **The list.** Item 4: "Use exact symbolic subtraction where feasible, and print the exact residual as computed.
  Where it is feasible, evaluation at points is added beside it, never in its place."
- **The gap.**
  - The reader builds a Wolfram derivative at composite arguments as a substitution over every argument slot
    (l.418). It reads SymPy's printed derivative at a composite argument in the form SymPy prints. The exact
    difference only expands (l.726–729).
  - Both engines print such derivatives (D3). Where they are the same derivative of the same function at the same
    arguments, the two readings stay different trees. On those rows the exact residual measures the reader, not the
    engines, and the row rests on evaluation alone, which F2 weakens further.
  - Each leg's own comparison read those derivatives as one object (Claude's probe, D3). Grok's report notes that its
    exact comparison rewrote the derivative heads, and that the comparator leaves the difference in the residual (D3).
  - A second instance is an inequality printed with its sides exchanged, whose difference stays an unreduced `Xor`
    (Claude, D3).
- **What must be true.**
  - Where the two engines print the same object in different notations, the reader builds the same object. Examples
    are a derivative at composite arguments and an inequality with its sides exchanged. The exact residual is then the
    residual of the printed objects.
  - A different derivative (another slot, order or point) or a different relation is still read as different. Item 6's
    derivative-structure controls cover this.
  - This translates each engine's notation; it does not reconcile two representations. Content that one engine prints
    and the other does not stays in the residual.

**F2. Distinct live functions can share a test function, so evaluation cannot tell them apart.** (Claude F2.)
**CONFIRMED.**
- **The list.** Item 4: "every live function is replaced by test functions of several different forms". Item 2: "A
  mapped name is a substitution, never by itself a zero."
- **The gap.**
  - A live function's test function depends on its name only through a seed taken mod 13 (l.640, l.657). Otherwise it
    depends only on the form and the number of arguments (l.612–622).
  - Distinct live functions of the same arity can therefore share a test function in one form, or in every form.
    Evaluation cannot tell them apart, and a name pointed at the wrong function gives zero samples.
  - Claude's one-sided corruption C2 exchanged two such functions in one engine's rows, and no sample moved (D4).
  - The name-table repoint control checks only exact residuals (test l.134–136), so it does not see this.
- **What must be true.**
  - Distinct live functions get distinct test functions, in every form and stratum. The exception is the deliberately
    equal one-sided values of a jump stratum.
  - A control fails when two coincide. It checks the evaluation as well as the exact residual.

**F3. No control goes through the current rows' component pairing.** (Claude F3.) **CONFIRMED.**
- **The list.** Item 6 requires "a repoint of a component correspondence (item 1)", and says "Each fixture … goes
  through the same reader, extraction, join and comparison path as the production run".
- **The gap.**
  - The current rows are paired in `currents()` (l.970), not through `tree()`'s explicit paths (l.957). The only
    correspondence-repoint control goes through `tree()` (test l.260–261).
  - Changing only l.970 to pair by position left every test passing (Claude, D5).
  - Grok's position ablation changed both l.957 and l.970. It failed that control through the l.957 change (D5), so
    the two legs' results agree.
  - Both legs found the production pairing correct.
- **What must be true.** A control through the current rows' pairing, the unpaired time-density included, fails when
  that pairing is changed.

**F4. Dimension-ledger terms are bound by being identical, then compared with each other and counted as compared.**
(Claude F4, minor.) **CONFIRMED.**
- **The list.** Item 6 asks for two leaf counts per joined row, one of them "the count actually compared", and says "A
  row compared below its parsed count stays visible."
- **The gap.**
  - Terms are bound across engines by their parsed object (l.997). The bound terms are then compared with each other
    (l.1014), a comparison that cannot differ.
  - Those leaves are credited as compared (l.935–937). The "actually compared" counts on dimension rows therefore
    include comparisons that carry no information, which can hide a row compared below its parsed count.
- **What must be true.** A binding key is neither printed nor counted as a compared value.

## Notes for the repair (no separate finding)

- **N1.** The component-correspondence citations point at the lines that emit each object, not the lines that build
  it (Claude, D7). Item 1 asks for the lines that build them. Both legs found the pairing itself correct.
- **N2.** Under Claude's ablation A2c, which keeps only forms depending on the sum of the arguments, the suite fails.
  But the dedicated unseparated-arguments coverage control, run alone, still passes (D7). Item 6 asks that each
  coverage control fail under its own reduction.
- **N3.** In the Claude leg, a guarded job of a few seconds overlapped another guarded job. The pooled budget admitted
  both, neither was killed, and no kernel ran (D7). It has no effect on the evidence.

## On the CLEAR leg

Grok recorded F1's underlying fact without reporting it (D3). It did not test whether distinct live functions get
distinct test functions (F2). Its position ablation changed two lines at once, so its control failure came from the
other change (F3). Legs have been wrong in both directions. These four findings rest on the comparator's source lines
and on the legs' scripts and stdout.

## Outcome

Four findings, confirmed. The reviewed build is preserved unchanged in the commit that adds this file. It is not
accepted. The same builder repairs it against F1–F4 and N1–N2 (its first repair). Review round 1 follows, with a fresh
Claude agent and Grok.
