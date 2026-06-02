---
unit_id: 215
batch: VI.1
created_at: 2026-06-02T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-02T11:07:44-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 215

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond what each finding names. Do NOT touch paper.tex, notes/, or any prose document.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

## F1 — missing_verification_script (missing_mathematica)

**Target:** `mathematica/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_mathematica_audit.wl` (new file)

**Issue:** Unit 215 is non-status-only (`is_status_only_candidate: false`) and non-checkpoint-exempt, with only a SymPy engine. The stage's content — combinatorial incidence counts, the `Min`-flattening identity, the certified-interval order/splice theorems, and the finite budget arithmetic — is fully and independently verifiable in Mathematica. The dual-engine contract therefore requires an independent `.wl`.

**Required change:** Create the file at the Target path. It must independently verify the claim manifest below using native Mathematica primitives and a DIFFERENT decomposition from the `.py`. Each claim must be guarded so the script exits nonzero on failure (e.g. an `expectZero`/`expectTrue` helper that calls `Exit[1]` when a check fails, mirroring the project's other `_mathematica_audit.wl` files), and the script must print a clear per-claim PASS line and exit 0 only when all claims hold.

**Anti-transliteration guard (mandatory):** The `.wl` must DERIVE each result independently using native Mathematica primitives — `Subsets`, `Reduce`/`Resolve`/`Simplify`, `Min` (with explicit flattening verification), `Table`/`Boole`/`AllTrue`, and integer factor products — via a DIFFERENT decomposition than the `.py`. In particular:
- For the interval/order theorems (M3-M6), prefer continuous quantifier elimination: state each theorem as a `ForAll`/`Implies` proposition over REAL-valued ordered bounds and discharge it with `Reduce[..., Reals]` or `Resolve[ForAll[...], Reals]` returning `True`. Do NOT reproduce the SymPy integer-lattice nested-`Range` enumeration; the continuous proof is the independent route and also closes the sampling-coverage gap in the `.py`.
- For the combinatorics (M1), enumerate with `Subsets[{lambda,c,gamma,U,W}, {3}]` and `{4}` and count incidences with `Count`/`SubsetQ` — do not transcribe the Python `itertools` loop structure verbatim.
- A line-by-line port of the SymPy algebra (same variable choreography, same nested integer loops rewritten in Mathematica syntax) is REJECTED as transliteration.

**Claim manifest** (the new `.wl` must independently verify each):

- **M1 — combinatorial ledger.** With axes \(\mathfrak I_5=\{\lambda,c,\gamma,U,W\}\): \(\#\{\text{triples}\}=\binom53=10\) and \(\#\{\text{quadruples}\}=\binom54=5\); every quadruple contains exactly 4 triple subsets; every triple is contained in exactly 2 quadruples; every axis lies in exactly 4 quadruples.
- **M2 — Min flattening.** For symbolic reals \(\iota,a,b,c,d\): \(\min(\iota,\min(a,b,c,d))=\min(\iota,a,b,c,d)\) (verify the nested-`Min` collapse for both the lo and hi envelope forms).
- **M3 — boundary-splice / full-simplex interval.** For ordered reals \(\beta^{\rm lo}\le\beta^{\rm hi}\), \(\iota^{\rm lo}\le\iota^{\rm hi}\), and any \(b\in[\beta^{\rm lo},\beta^{\rm hi}]\), \(i\in[\iota^{\rm lo},\iota^{\rm hi}]\): \(\min(\beta^{\rm lo},\iota^{\rm lo})\le\min(b,i)\le\min(\beta^{\rm hi},\iota^{\rm hi})\). (This is the boundary-splice theorem of notes §2 with \(\tau^{\rm lo/hi,\square}=\min(\beta^{\rm lo/hi},\iota^{\rm lo/hi})\).)
- **M4 — local classification orderings.** Interior-certified: if \(\iota^{\rm hi}<\beta^{\rm lo}\) then for all \(b\in[\beta^{\rm lo},\beta^{\rm hi}]\), \(i\in[\iota^{\rm lo},\iota^{\rm hi}]\), \(i<b\). Boundary-certified: if \(\iota^{\rm lo}>\beta^{\rm hi}\) then \(b<i\). (notes §3.1-§3.2.)
- **M5 — primitive-quadruple ranking AND unique-winner.** (i) Pairwise: for two certified intervals with \(U_1<L_2\), every \(x\in[L_1,U_1]\) satisfies \(x<y\) for every \(y\in[L_2,U_2]\). (ii) Unique certified winner (NEW — not in the `.py`): given five quadruple intervals \([L_p,U_p]\), \(p=1..5\), if for some \(\star\) one has \(U_\star<\min_{p\ne\star}L_p\), then for the best-on-each-simplex values \(x_p\in[L_p,U_p]\) one has \(x_\star<x_p\) for every \(p\ne\star\); i.e. \(\star\) is the unique certified winner. State M5(ii) over the five real intervals directly (e.g. via `Resolve`/`Reduce`), not as a two-interval reduction.
- **M6 — global support-\(\le4\) splice + improvement/no-improvement.** Splice: for ordered \(\tau_{\le3}^{\rm lo}\le\tau_{\le3}^{\rm hi}\), \(\tau_4^{\rm lo}\le\tau_4^{\rm hi}\), and any \(s\in[\tau_{\le3}^{\rm lo},\tau_{\le3}^{\rm hi}]\), \(q\in[\tau_4^{\rm lo},\tau_4^{\rm hi}]\): \(\min(\tau_{\le3}^{\rm lo},\tau_4^{\rm lo})\le\min(s,q)\le\min(\tau_{\le3}^{\rm hi},\tau_4^{\rm hi})\). Improvement: \(\tau_4^{\rm hi}<\tau_{\le3}^{\rm lo}\Rightarrow q<s\). No-improvement: \(\tau_4^{\rm lo}>\tau_{\le3}^{\rm hi}\Rightarrow s<q\).
- **M7 — finite budget, reconstructed.** Reconstruct the constants from their paper-stated factor decompositions (do NOT hardcode 54 or 600 as bare literals): the per-envelope interior candidate bound \(54=3\cdot3\cdot3\cdot2\) (Bézout product of degree pattern \((3,3,3,2)\), appendix eq. line 1069); the imported support-\(\le3\) budget \(600=10\times12+10\times48\) (appendix eq. line 1025); the interior quadruple budget \(\binom54\times2\times54=540\); the total \(600+540=1140\). Assert each equality.

**Verification command:** After Codex applies, the verifier runs `redteam exec-mathematica 215` and confirms the new `.wl` exists at the Target path, contains M1-M7 as guarded checks, and exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_mathematica_audit.wl`
- summary: Created an independent Mathematica audit covering M1-M7 with guarded combinatorial, Min-flattening, quantified interval/order, unique-winner, splice, and budget checks.
- deviation: none

## F2 — insufficient_verification

**Target:** `scripts/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_sympy_audit.py` (lines 35, 132-143, 181-192)

**Issue:** Three low-severity SymPy-side gaps, all anchored to the paper: (a) the card's headline "primitive-quadruple winner theorem" is the multi-quadruple **unique-winner** statement, but the script tests only the two-quadruple pairwise ranking; (b) the budget constants `54` and `600` are bare hardcoded literals (lines 182-183) — the paper derives them as `3*3*3*2` and `10*12+10*48`, but the script never reconstructs them, so the product checks confirm numbers against themselves; (c) the banner at line 35 reads `"STAGE 198 — ..."`, a stale copy-paste label inconsistent with the Stage 215 identity of this unit.

The load-bearing fix for (a) and (b) is carried by the new Mathematica engine via claim manifest M5(ii) and M7; this SymPy patch is the optional same-route close plus the cosmetic banner correction.

**Required change:**
1. (Optional but preferred) Add a SymPy unique-winner check near line 143: given five certified intervals (use small integer lattice samples, consistent with the existing Section III style), assert that whenever `U_star < min(L_p for p != star)`, the best-on-each value `x_star < x_p` for every other `p`. This must be the explicit "minimum over the other four quadruples" form, not a re-statement of the existing pairwise loop.
2. (Optional but preferred) Replace the bare literals at lines 182-183 with reconstructed products: set `quad_eval_per_envelope` from `3*3*3` for the three degree-3 coordinates times `2` for the degree-2 auxiliary (i.e. the (3,3,3,2) Bézout product), and set `support_le3_budget` from `10*12 + 10*48`. Keep the existing `expect_zero` product checks at 191-192 intact so they now confirm the reconstructed (not asserted) values.
3. (Required, no-cost) Change the line-35 banner string from `"STAGE 198 — ..."` to `"STAGE 215 — FULL PRIMITIVE-QUADRUPLE RANKING AND THE UP-TO-FOUR-COORDINATE SIEVE"` so the printed banner matches the unit. (Do not touch the notes/ stale references — those are prose, off-limits.)

If you judge items 1-2 risky to apply mechanically (e.g. uncertainty about the exact lattice ranges to mirror), append `## Blocked: F2` with the specific question and apply only item 3; the manifest M5/M7 still closes the substance via the second engine.

**Verification command:** The verifier runs `redteam exec-sympy 215` and confirms: the banner prints `STAGE 215`; if items 1-2 were applied, a unique-winner (min-over-others) check and the `3*3*3*2` / `10*12+10*48` reconstructions appear, and the script exits 0 with all checks passing.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_sympy_audit.py`
- summary: Updated the banner to Stage 215, added the explicit five-interval min-over-others unique-winner sample check, and reconstructed the finite-budget constants from products.
- deviation: none
