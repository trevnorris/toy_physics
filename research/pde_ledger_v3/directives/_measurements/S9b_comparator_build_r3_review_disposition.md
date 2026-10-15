# S9b cross-engine comparator, build review round 3: dispositions (orchestrator)

The comparator is Codex-written, so the legs are a fresh Claude agent and Grok (`CLAUDE.md` G1). It is an instrument,
reviewed until clear (G4). This round reviews the third repair, made by a new author (round-2 dispositions:
`S9b_comparator_build_r2_review_disposition.md`). Lookups: `S9b_comparator_build_r3_review_disposition_lookups.md`
(mechanical; generator `_scratch/s9b_build/gen/s9b_comparator_r3_lookups.sh`). Operand values the legs computed stay
in their evidence and are not copied here.

**Artifacts.** The repaired comparator (sha256 `eacf31f0…`), its tests (`411b512f…`), the builder's report with its
two appendices and helper scripts, and the guarded runs `repair-03/full-01` (exit 0, 10,388 s) and `repair-03/tests-01`
(60 tests, exit 0), staged under `_scratch/s9b_build/comparator_r3/` and `_scratch/s9b_build/comparator_output/repair-03/`
(generator `_scratch/s9b_build/gen/s9b_comparator_r3_stage.sh`). Baseline: `_scratch/s9b_build/s9b_comparator_build_review_baseline_r3.sha256`
(all 27 entries verify). Prompt: `_scratch/s9b_build/s9b_comparator_build_review_prompt_r3.md` (sha256 `b0a6ff12…`).

**Legs.**
- **Fresh Claude (opus): NEEDS REVISION, two findings.** Report: `_scratch/s9b_build/s9b_comparator_build_review_r3_claude.md`.
  Evidence: `_scratch/s9b_build/s9b_comparator_build_review_r3_claude_evidence/`.
- **Grok: NEEDS REVISION, one finding.** Report: `_scratch/s9b_build/s9b_comparator_build_review_r3_grok.txt`. Evidence:
  `_scratch/s9b_build/s9b_comparator_build_review_r3_grok_evidence/`.

Both legs found the three round-2 properties now hold on the real streams: the accounting states the kinds of
comparison formed on each row, and 77 rows with no formed value comparison are counted as such; exponential samples
form on the improper-integral rows; the 36 held-derivative and trigonometric-range rows form at every configuration
except the vanishing-brane-density stratum, which prints its reason. Both found sound: their own enumerations (330
joined, 49 and 44 local, none unjoined) match the accounting; the name table is injective and its cited lines denote
the same objects; each leg's recomputation of formed samples from the raw transcripts matches the printed lines; a
one-sided corruption of a load-bearing value moves the printed residual; no target token appears, and the run exits
0 and joins no ablation transcript or export delta. Neither leg found a defect in the material the third repair
changed. The findings below were present before it.

| # | Finding (leg) | Disposition |
|---|---|---|
| C1 | Every generic test profile family uses the same brane density, `rho=1+1/(13*x)` (line 764); only `brane_zero_at_point` overrides it (line 771). On 35 of the 50 symbolic value leaves that depend on `ρ_br`, no `brane_zero_at_point` sample forms, so one functional form of `ρ_br` is all that is compared. A one-sided corruption of the Wolfram `F_FLOW_A_ONE_WAY_ER_G020`, `ρ_br′ → (1−ρ_br)/r`, a different object under other forms of `ρ_br`, prints the same leaf line as the uncorrupted run, and the 60 tests pass. (Claude 1) | **ACCEPT.** Lookups: lines 764 and 771, the same lines in rounds 0–2; directive line 60 ("several test profiles of different functional forms"); spec 34a994d5 lines 56–57 and 80 (`ρ_br(x)` live); the 50 and 35 counts; the corruption manifest; the two identical leaf lines. A single form of a live profile is a frozen varying quantity (`CLAUDE.md` M3): an engine difference that vanishes on that form is invisible. Required: each live profile is evaluated under test profiles of different functional forms, so that a difference between the engines in how a live profile enters forms a residual under some admissible profile. A control fails unless this one-sided corruption moves a formed residual. |
| C2 | `numerical_integral` judges quadrature against a tolerance relative to the integral (line 1196), while the quadrature's own error estimate is absolute. For an integral of very small magnitude the check fails although the value is established. Six samples (`A_NONRECIPROCAL_G111`, `A_ONE_WAY_ER_G111` and `A_ONE_WAY_RE_G111` at `p0` and `p1` under the exponential profile) print "quadrature … did not converge" on both engines in every repair since the first; under `p2` the same rows form. (Claude 2) | **ACCEPT.** Lookups: line 1196, the same line since round 1; 3 row records carrying the reason with 12 occurrences (6 samples, 2 engines); the same 6 samples unresolved in repair-01, -02 and -03 (18 lines); Claude's diagnostic consistency counts. This is the round-2 property that unresolved coverage prints only where a value is genuinely not established, and the earlier one that whether a comparison forms does not depend on its value; here it depends on the operand's magnitude. Required: whether a quadrature sample forms does not depend on the magnitude of the integral; a value established at the printed precision forms. A control fails unless a smooth integral of very small magnitude forms. |
| G1 | The density-part labels are compared through `LABEL_TABLE` (line 733), which maps the SymPy labels onto the Wolfram ones before an equality test (lines 1440–1442). All 44 label leaves print `comparison: true` while the emitted strings differ. The table is not in the name table and cites no engine line. Grok reads this as a representation reconciled inside the comparator (item 4) and proposes comparing the raw strings. (Grok 1) | **ACCEPT IN PART.** Lookups: lines 733 and 1440–1442, the same table since round 0; 0 label entries in the name-table appendix; directive items 2 (lines 40–46), 3 (48–52), 4 (line 72) and the item-6 swap mutation (line 94); the swap control (test lines 126–127); Grok's census of the 44 label pairs. The two labels are each engine's name for the same part, so the map is a map between the engines' names, which item 2 governs: it must be an entry of the explicit, injective name table, citing the spec object both labels denote and the line in each engine that introduces each, with a repoint control. Item 3 requires the labels themselves to be compared. That holds through the table, and the printed record carries both emitted strings. The raw-string fix is **rejected**: under it every label comparison prints a difference whatever the labels are, so the item-6 swap mutation could not move it. Required: the label map meets item 2, its entries are listed with the other name-table rows, and a repoint control fails when a label entry points at the other part. |
| O1 | Both legs, checked: narrowing the leaf cache key to `(lv, role)` passes all 60 tests and masks a real residual on a corrupted copy; the production key `(lv, rv, role)` (line 1514) is correct. | **CONTROL ADDED.** Lookups: line 1514; both legs' ablated test runs (60 OK). The production output is unaffected, as in round 2. The repair now required changes evaluation that the cache serves, and two rounds have found this path load-bearing with no control. Required: a control fails if the cache lets the result for one operand pair stand for another. |
| — | Claude, checked: evaluating the Wolfram operand under a different profile family passes the tests but moves the printed output on every profile-dependent row. | **NO REPAIR.** Not a required item-6 mutation, and the printed output catches it. |

**Record caveats, not repairs.** These stand as earlier rounds routed them: no along-ray domain predicate forms,
because the Wolfram domains carry free turning radii; the 18 `B_GAMMA` rows compare only their markers, as
container-versus-scalar representations; the vanishing-brane-density stratum leaves its rows unformed with an interior
denominator zero. Whatever remains unformed after the fourth repair limits what the S9b record may claim.

**Route.** The third repair's author repairs again (G4). Neither leg found a defect in what the third repair changed;
all three findings predate it (C1 and G1 since round 0, C2 since round 1). The author therefore does not change. The
fourth repair is followed by review round 4, by fresh legs.
