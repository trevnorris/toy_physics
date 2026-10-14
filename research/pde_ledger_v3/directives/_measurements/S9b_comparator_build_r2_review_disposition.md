# S9b cross-engine comparator, build review round 2: dispositions (orchestrator)

The comparator is Codex-written, so the legs are a fresh Claude agent and Grok (`CLAUDE.md` G1). It is an instrument,
reviewed until clear (G4). This round reviews the second repair (round-1 dispositions:
`S9b_comparator_build_r1_review_disposition.md`). Lookups: `S9b_comparator_build_r2_review_disposition_lookups.md`
(mechanical; generator `_scratch/s9b_build/gen/s9b_comparator_r2_lookups.sh`). Operand values the legs computed stay
in their evidence and are not copied here.

**Artifacts.** The repaired comparator (sha256 `a9071a94…`), its tests (`7c2b2159…`), the builder's report
(`278d646d…`) with its two appendices, and the guarded runs `repair-02/full-01` (`comparison.jsonl` sha256
`8ab9e5ad…`; `accounting.jsonl` `0e6f002f…`) and `repair-02/tests-04`, staged under `_scratch/s9b_build/comparator_r2/`
and `_scratch/s9b_build/comparator_output/repair-02/`. Baseline: `_scratch/s9b_build/s9b_comparator_build_review_baseline_r2.sha256`
(all 22 entries verify). Prompt: `_scratch/s9b_build/s9b_comparator_build_review_prompt_r2.md`.

**Legs.**
- **Fresh Claude (opus): NEEDS REVISION, three findings.** Report: `_scratch/s9b_build/s9b_comparator_build_review_r2_claude.md`.
  Evidence: `_scratch/s9b_build/s9b_comparator_build_review_r2_claude_evidence/`.
- **Grok: CLEAR.** Report: `_scratch/s9b_build/s9b_comparator_build_review_r2_grok.txt`. Evidence:
  `_scratch/s9b_build/s9b_comparator_build_review_r2_grok_evidence/`.

Both legs found sound: their own enumerations of the two transcripts (330 joined, 49 and 44 local, none unjoined)
match the accounting; both name tables are injective and their cited lines denote the same spec objects; the
integration-constant row is gone and both names stay engine-qualified; whether a residual forms no longer depends on
its value; the sample points are generic; a one-sided corruption of a formed value row moves the printed output; each
leg's recomputation of formed samples from the raw transcripts matches the printed lines; no target token appears,
the run exits 0, and no ablation transcript or export delta is joined. The three findings below concern what the
output can be cited as having compared, not a wrong printed value.

| # | Finding (leg) | Disposition |
|---|---|---|
| F1 | The row-level `comparison_formed` flag (line 1468, `any(used)`) counts a formed `NOT_ESTABLISHED` marker, density label or label condition as a formed comparison, so the report's "Rows with no formed comparison: 37" overstates what was compared. 76 rows are flagged formed with no value comparison formed: 36 mixed-grade one-way and round-trip rows, 18 `B_GAMMA`-family rows and 22 Part C rows. A one-sided corruption of a Part C object changed no computed product and left the row flagged formed. (Claude 1) | **ACCEPT.** Lookups: lines 1214 and 1468–1469; 293 rows flagged `true` and 37 `false` in `accounting.jsonl`; the `A_ONE_WAY_ER_G021` row (one compared leaf per engine, flagged formed); the report line; Claude's per-kind census and its Part C ablation line. Grok quoted the same 37 without classifying it. Item 6 asks for "the count actually compared". Required: the accounting states, per row, which kinds of comparison formed (value, domain, marker, label), and the rows with no formed value comparison are counted and reported as such. No report line counts a marker or label comparison as a compared value. A control fails if a row whose only formed comparison is a marker or a label is counted as having a formed value comparison. The flag is new in the second repair (the round-1 comparator has none; lookup) and is its answer to round 1's R3, so this is a defect in the material that repair changed. |
| F2 | The exponential test profile never forms on the 34 improper-integral rows (every `A_DEFLECTION`, `A_RADAR_LOG` and `B_RESIDUAL` row and their `F_LIVE`/`R_*` copies). `integral_screen` (lines 1015–1073) calls `as_leading_term` (line 1064), which raises `NotImplementedError` when `exp` is present; the fallback needs a finite plain limit, which an integrable inverse-square-root endpoint does not have, so it prints "endpoint convergence not established" (line 1073). Under the inverse and shifted profiles the same rows form. On these rows the only profile families compared are the two rational ones. (Claude 2) | **ACCEPT.** Lookups: the cited lines; 34 rows of `comparison.jsonl` carry the reason; Claude's per-row census (exponential formed 0 of 3 on all 34 rows, inverse and shifted 3 of 3) and its screen check (the screen passes inverse and shifted, rejects exponential at each point). Grok found the unformed outcomes printed faithfully, which this does not contradict: the printing is faithful, but the screen gives up for a limit of the comparator's own construction, not because convergence is undetermined. Item 4 asks for test profiles "of different functional forms" for an object that is a functional of a live profile. Required: an improper integral forms under every test profile for which its convergence at the endpoint can be established; unresolved coverage is printed only where convergence is genuinely not established, and its reason says which. A control fails unless a convergent improper integral with an integrable endpoint singularity forms under a non-rational test profile. |
| F3 | No sample forms on 36 mixed-grade one-way and round-trip value rows (`A_ONE_WAY_ER`, `A_ONE_WAY_RE`, `A_ROUND_TRIP` at G021, G101, G120, G121, with their `F_` copies). The SymPy operand holds a `Subs` of a `b`-derivative of an integral whose limit is the derivative variable, printed as "held derivative, substitution or limit" (line 954); the Wolfram operand integrates trigonometric functions over an angle range whose pole set the screen does not establish (line 1026). Both operands are finite functionals of the same profiles; the printed outcome is faithful, but nothing is measured for these objects. (Claude 3) | **ACCEPT.** Lookups: the cited lines; 36 rows carry the held-derivative reason (1,260 occurrences) and 33 the pole-set reason (813); Claude's census shows all 36 rows' values with only unresolved coverage. Round 1 routed these rows to the record as a caveat. This round shows that the cause is how the comparator evaluates the two engines' representations, not an object without a value at the sample points, and the user chose to repair all three findings (2026-10-14). Required: wherever both operands are defined at a sample point, the comparison forms; unresolved coverage is printed only where an operand's value there is genuinely not established, and its reason names what is not established about the object. Neither representation is reconciled into the other (item 4). A control fails unless a held derivative of an integral with a variable limit, and an integral of a trigonometric integrand over a range with no pole, each form. |
| O1 | Claude, checked and found sound: narrowing the leaf cache key to `(lv, role)` passes all 52 tests and, on corrupted streams, masks real residuals; the production key `(lv, rv, role)` is correct, so it is a test-suite gap, not a current defect. | **NO REPAIR.** Lookups: Claude's cache-key ablation lines. The production key is correct and the output is unchanged. The third repair changes the evaluation path the cache serves, and its review ablates it again. |
| O2 | Claude, minor: a stratum whose condition appears under two labels can switch between "relation true" and "union piece 0" between runs on identical input; coverage is unchanged. | **NO REPAIR.** No computed comparison or count changes. |

**Record caveats, not repairs.** These stand as round 1 routed them, apart from F3: no along-ray domain predicate
forms, because the Wolfram domains carry free turning radii. Claude's one-sided domain corruption of
`A_DEFLECTION_G010` moved no formed comparison and showed only in the printed condition text. Whatever remains
unformed after the third repair limits what the S9b record may claim.

**Route.** The author changes (G4). Round 1's R1 was a defect in the material the first repair added, and F1 is a
defect in the material the second repair changed. That is two successive repairs breeding defects in what they
changed. The user chose a full third repair of all three findings by a new author, followed by review round 3
(2026-10-14). The new author is a fresh Codex session, working in a new export, and is reviewed by fresh legs.
