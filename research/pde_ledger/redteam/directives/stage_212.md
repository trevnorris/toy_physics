---
unit_id: 212
batch: VI.1
created_at: 2026-06-01T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-06-02T11:39:51-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 212

Apply EVERY finding below in order (F1, F2) AND the now-RESOLVED notes edits (R1, R2 — see the `## RESOLVED` block). After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

F2 is a SCRIPT-side banner fix (`.py` line 35) despite its "paper_misalignment" title — apply it. The two NOTES discrepancies (R1, R2) have been RESOLVED by the user (direction (a) for both); you are now EXPLICITLY AUTHORIZED to edit the notes file as specified in the `## RESOLVED` block. The SymPy script itself is correct and must NOT change (other than the F2 banner).

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — missing_verification_script (missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_mathematica_audit.wl`

**Issue:** Stage 212 is non-status-only and non-checkpoint, yet has only a SymPy audit; the stage card states "Mathematica audit: none yet." A second engine is required because every claim is independently verifiable in Mathematica. Create the `.wl` so that it independently verifies the claim manifest below and exits 0.

**Required change:**
Write a new Mathematica script at the Target path that verifies M1–M4. Each verifiable result must be checked with a guard that exits nonzero on failure (e.g. `If[result =!= expected, Print["FAIL: M_n"]; Exit[1]]`, or a `Throw`/`Assert` pattern that aborts with nonzero exit). Print each result before the guard. The script must derive each claim INDEPENDENTLY using native Mathematica primitives — `Subsets`, `Binomial`, `Count`, `Min`, `Reduce`, `Resolve`, `ForAll`, `Element[..., Reals]` — via a DIFFERENT decomposition than the `.py`. In particular, the four interval/order theorems (M3) must be proved by symbolic quantifier elimination over the reals (`Resolve[ForAll[...], Reals]` returning `True`), NOT by the `.py`'s finite integer enumeration. A line-by-line port of the SymPy `itertools` loops or the integer-range `for`-loops is rejected as transliteration.

**Claim manifest** (the new `.wl` must independently verify each):

- **M1 (combinatorial ledger).** Over the five primitive axes \(\mathfrak I_5=\{\lambda,c,\gamma,U,W\}\):
  - \(|\{\{i,j\}\subset\mathfrak I_5\}| = \binom{5}{2} = 10\) and \(|\{\{i,j,k\}\subset\mathfrak I_5\}| = \binom{5}{3} = 10\).
  - Every primitive pair is contained in exactly 3 primitive triples.
  - (consistency, as the `.py` also checks) every primitive axis is contained in exactly 6 primitive triples.
  Route: build the pair/triple sets with `Subsets[axes,{2}]` and `Subsets[axes,{3}]`; compare `Length` against `Binomial[5,2]`/`Binomial[5,3]`; for incidence use `Count[triples, t_ /; SubsetQ[t, pair]]` and `Count[triples, t_ /; MemberQ[t, axis]]`.

- **M2 (boundary-splice nested-min identity).** For symbolic reals \(\iota, a, b, c\):
  \[ \min(\iota,\ \min(a,b,c)) = \min(\iota, a, b, c). \]
  Verify both the lo-form (with \(\beta^{\rm lo}=\min(\tau_{ij}^{\rm lo},\tau_{ik}^{\rm lo},\tau_{jk}^{\rm lo})\)) and hi-form. Route: assert `Min[iota, Min[a,b,c]] === Min[iota,a,b,c]` (Mathematica auto-flattens `Min`); this confirms the splice definition \(\tau_{T,\min}^{\rm lo,\triangle}=\min(\beta_T^{\rm lo},\tau_{ijk,\min}^{\rm lo,int})\) collapses to a single flat minimum over the interior bound and the three edge cones.

- **M3 (interval / order theorems, symbolic).** For real ordered bounds, prove each as a universally-quantified inequality via `Resolve`/`Reduce` over `Reals` (this is the independent-decomposition requirement — do NOT enumerate integers):
  - **M3a local full-simplex sandwich:** for all \(\beta^{\rm lo}\le\beta^{\rm hi}\), \(\iota^{\rm lo}\le\iota^{\rm hi}\), \(b\in[\beta^{\rm lo},\beta^{\rm hi}]\), \(i\in[\iota^{\rm lo},\iota^{\rm hi}]\):
    \[ \min(\beta^{\rm lo},\iota^{\rm lo})\ \le\ \min(b,i)\ \le\ \min(\beta^{\rm hi},\iota^{\rm hi}). \]
  - **M3b classification orders:** if \(\iota^{\rm hi}<\beta^{\rm lo}\) then for all admissible \(i,b\): \(i<b\) (interior-certified); if \(\iota^{\rm lo}>\beta^{\rm hi}\) then for all admissible \(i,b\): \(b<i\) (boundary-certified).
  - **M3c primitive-triple ranking:** if \(U_1<L_2\) then for all \(x\in[L_1,U_1]\), \(y\in[L_2,U_2]\) (with \(L_1\le U_1\), \(L_2\le U_2\)): \(x<y\).
  - **M3d global splice sandwich + improvement/no-improvement:** the support-\(\le3\) sandwich \(\min(p^{\rm lo},t^{\rm lo})\le\min(p_*,t_*)\le\min(p^{\rm hi},t^{\rm hi})\); and if \(t^{\rm hi}<p^{\rm lo}\) then \(t_*<p_*\) (improvement), if \(t^{\rm lo}>p^{\rm hi}\) then \(p_*<t_*\) (no-improvement).
  Route: `Resolve[ForAll[{vars}, hypotheses \[Implies] conclusion], Reals]` must return `True` for each; guard `If[res =!= True, Print["FAIL M3x"]; Exit[1]]`.

- **M4 (finite evaluation budget).** With 10 primitive pairs at 12 evaluations each (6 per envelope × 2 envelopes) and 10 primitive triples at 48 interior evaluations each (24 per envelope × 2 envelopes):
  \[ 10\times 12 = 120, \qquad 10\times 48 = 480, \qquad 120+480 = 600. \]
  Verify each literal equality. (Note: use 120 for the pairwise total, matching the script line 202, the stage card "budget 600", and the appendix `10\times12+10\times48=600`. Do NOT use the notes literal `188` — see the Resolve block; the correct value is 120 and is the side the script already uses.)

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 212` and confirm the new checks (M1–M4) appear AND the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_mathematica_audit.wl`
- summary: Added the independent Mathematica audit covering M1-M4 with native combinatorics, Min flattening, quantified real-order theorems, and budget checks.
- deviation: none

## F2 — paper_misalignment (notes_contradicts_script) — script-side label fix

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_sympy_audit.py:35`

**Issue:** The banner prints `STAGE 195` although this is the stage-212 audit and the closing line (206) already says "All Stage 212 identities... verified." This is a copy-from-template residue in the script (safely script-side correctable). It is unrelated to the notes naming question in the Resolve block below.

**Required change:**
Edit line 35 only. Before:
```
banner("STAGE 195 — FULL PRIMITIVE-TRIPLE RANKING AND THE UP-TO-THREE-COORDINATE SIEVE")
```
After:
```
banner("STAGE 212 — FULL PRIMITIVE-TRIPLE RANKING AND THE UP-TO-THREE-COORDINATE SIEVE")
```
Do not change any other line. Re-run `python3 <path>`; it must still exit 0.

**Verification command:**
After Codex applies, `redteam exec-sympy 212` re-runs; the saved output header banner reads `STAGE 212 — FULL PRIMITIVE-TRIPLE RANKING...` and exit code is 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_sympy_audit.py`
- summary: Corrected the Stage 212 SymPy audit banner from `STAGE 195` to `STAGE 212`.
- deviation: none

## RESOLVED (user direction: (a) for BOTH R1 and R2 — notes edits AUTHORIZED)

The user confirmed the script is correct and AUTHORIZED the following notes-side prose edits. Codex applies the notes edit; Claude reviews the diff (file-ownership rule). The SymPy script is NOT changed for R1/R2 (it already uses the correct values). Target notes file: `notes/stages/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem.md` (the stage-212 notes; confirm the exact filename by glob if needed).

**R1 — pairwise-budget typo `188` → `120` (CONFIRMED, direction (a)).** In the notes, correct §7.1 `10\times 12 = 188` → `10\times 12 = 120`, and §7.3 `188+480=600` → `120+480=600`. (10×12=120 is arithmetically forced and matches the script line 202, the card budget 600, and the appendix `10×12+10×48=600`.) No script change.

**R2 — stale body stage-labels (direction (a) — renumber to current scheme).** The notes body calls this stage "Stage 246" and its imports "Stage 245 / Stage 243"; the notes title and the stage card both name it Stage 212. Renumber:
- All self-references "Stage 246" → "Stage 212" (unambiguous — confirmed by the notes title + card).
- The import references "Stage 245" and "Stage 243" → their CURRENT canonical upstream stage numbers. Determine these from what Stage 212 actually imports (read the import context in the notes/script and map to the canonical stage that owns each imported result). The offset 246→212 is −34, which would map 245→211 and 243→209; apply that mapping ONLY if 211/209 are in fact the stages that own the imported results — verify before writing. If an import target is ambiguous, renumber the unambiguous self-references (246→212) and append a `## Blocked: R2-import` note naming the ambiguity rather than guessing. This is notes-prose only; no verified identity changes.

After the notes edits, re-run the scripts (`python3` / `math -script`); they must still exit 0 (the notes edits do not touch any script, so the script results are unchanged — this is just the F1/F2 re-run confirmation).

## Applied: R1/R2

- files_changed:
  - `notes/stages/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem.md`
- summary: Corrected the pairwise budget typo to `120` and renumbered the Stage 212 notes self/import references to canonical Stage 212, Stage 211, and Stage 209 labels.
- deviation: none

---

## Iteration 2 (orchestrator) — R2 follow-up: one missed forward-reference

The R2 renumber (iter 1) correctly fixed the self-references (246→212) and the imported-upstream references (243→209, 245→211), but ONE stale forward-reference was missed because the iter-1 RESOLVED block only enumerated 246/245/243. Fix it now (notes-prose only, user-authorized direction (a), same −34 scheme):

- In `notes/stages/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem.md`, line ~498 (content-anchor on the sentence), the forward pointer reads: "**Stage 247** should build the exact four-coordinate positive simplex, reduce its full boundary to the Stage 212 triple packets, ...". Renumber **`Stage 247` → `Stage 213`** (247 − 34 = 213; Stage 213 is the four-coordinate mixed-simplex stage, which is the actual next stage in the current scheme). Change ONLY that one token; leave the rest of the sentence (including the already-correct "Stage 212 triple packets") intact.

No script change; do not touch any other line. After the edit, grep the notes for `Stage 24` / `Stage 23` to confirm no other stale ≥240 stage token remains (there should be none).
