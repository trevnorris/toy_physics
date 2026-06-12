# Phase C — Batch 2 verdicts (targeted free_choice sweep, stages 000-079)

Batch: `batch_c2_freechoice_000_079` (PC6 targeted sweep). Roster: `_phase_c_batch2_roster.yaml`.
31 free_choice candidates (provenance_built) -> 14 disjoint stage-cluster/family units, per-member sub-verdicts.
Adversarial framing = the **posit-dressed-as-derived** test (the stage_192 failure mode): is a *chosen* value labeled
derived/forced/exact/non-tunable/a-theorem in any paper card or per-stage notes line? `Benchmarks: []` is correct/expected
here (free_choice makes no external-match claim; the GR/CODATA test does not bite).

## HEADLINE: 14/14 units NO FATAL FLAW. 31/31 members NO. No YES verdicts -> no defense, no adjudication, no user gate.

Disclosure is consistently honest across the band: posited values are labeled "(new ansatz)" / "not yet frozen" /
"not derived here" / `\StatusOpen` / "choose" / "imported" / "natural ... not a universal theorem" / "frozen reference value" /
"selected earlier" / "carried ... reference-branch numerical freeze, not a new theorem". Every `\StatusExactClosure` /
"exact" / "lock" badge attaches to an **algebraic identity or closure-relative consequence**, never to the posited number
itself — confirmed against the paper's own legend (`paper/appendices/stage_ledger.tex:17`: ExactClosure = "exact AFTER the
stated reduced closure has been declared ... does NOT mean exact for every possible nonlinear moving-throat PDE branch").
This is the stage_192 vocabulary used CORRECTLY.

## Per-unit verdict table

| unit | stages | members | verdict | one-line basis |
|------|--------|---------|---------|----------------|
| u01_stage001_wall_action_coeffs | 001 | K_eta,T_w,mu_eta,T_Omega (1 cand) | NO | introduced under "(new ansatz)/not yet frozen" (notes:194/196); "fixed effective constitutive functions of the chosen reference throat ... not to be refit to rescue the chain" (notes:214-217) |
| u02_stage002_gamma_ab | 002 | Gamma_AB | NO | value-less deferred open-system completion; card says "not derived here" in 3 places (stage_002.tex:5/71/115) |
| u03_b0_normalization_chain | 011/020/023 | B_0 (3 cands) | NO | B_0 takes NO numeric value at any anchor -> no re-fit drift possible; honestly carried symbol |
| u04_stage015_s_em | 015 | S_EM | NO | imported parent-theory localized-Maxwell action term (stage_001.tex:9 "imports"); StatusExact scoped to declared parent data |
| u05_stage021_adapter_couplings | 021 | C_a,Omega_A_l (2 cands) | NO | structural gauge-invariant definition + "smallest useful reduced model"/"only the one-mode reduced prototype" (notes:82/314) |
| u06_stage022_023_coeffs | 022/023 | D_0,D_A,mhat_0,N_n (4 cands) | NO | all values `\StatusOpen`; "exact closure" scoped to bridge algebra; "branch test, not a proof of branch realization" |
| u07_mid_loading_source | 028/049 | alpha,sigma (2 cands) | NO | alpha an explicit Input (sign only posited; 029 derives magnitude); sigma "choose the uniform source density"; exact attaches to downstream conditional results |
| u08_n5_eos_exponent | 062/070 | K_eos,n (2 cands) | NO | n=5 consistently "frozen n=5 EOS" Input at both anchors; "deriving" applies to shell T_X/K_X coefficients, never to the exponent |
| u09_stage065_canonical_normalization | 065 | V_0,V_0f,f_prime_0 (3 cands) | NO | V_0 an open wall-amplitude solved-for; **V_0f is a scanner-fused token** (`V_0 f(...)`, stage_065.tex:7), not a real param; f'(0)=1 labeled "canonical normalization" (script:118) |
| u10_km_mouth_closure | 071/073 | K_m (3 cands) | NO | self-flagged "not a universal theorem ... first natural local Robin closure" (notes071:95); identical K_m=T_X/ell at every anchor |
| u11_aspect_37_20_origin | 073 | L_over_a (4 cands) | NO | classification already adjudicated free_choice; "reference-branch numerical freeze, not a new theorem" (notes073:56); ExactClosure badge closure-scoped |
| u12_epsilon_r | 073 | epsilon_r | NO | "frozen inputs ... epsilon_r=0.05=1/20" (notes073:18/30); grep for {deriv,forced,exact,non-tunable,theorem} near epsilon_r = 0 hits |
| u13_stage074_healing_ell | 074 | ell (2 cands) | NO | "a controlled local closure, not yet a theorem ... instead of leaving it free" (notes074:55); ExactClosure legend-scoped |
| u14_stage075_076_wall_norm | 075/076 | alpha_r,lambda_mu (2 cands) | NO | alpha_r=10 "frozen reference value"/"already carries"; lambda_mu an open convention ("whether one chooses lambda_mu=1 ..."), carried symbolically |

## Orchestrator spot-checks (distrust-all-clean backstop)
- **u09 V_0f**: VERIFIED against `paper/stages/stage_065.tex:7` = `V_{\rm conf}(r;a)=V_0f((r-a)/\ell)` -> "V_0f" is `V_0 * f(...)` fused by the scanner; agent's tokenization-artifact call CONFIRMED (no distinct parameter exists).
- **u11 stage_073**: VERIFIED `stage_073.tex:5` claimstatus closure-scoped, `:7` Inputs "selected earlier", `:11` `\subsection*{Derivation}`, `:17` `L/a=37/20` -> 37/20 is the INPUT to the derivation (37/20 / (1/20) = 37), not itself claimed-derived; agent's "presentational tension, not derive-vs-fit mismatch" CONFIRMED.
Both ratify the agents; the all-NO is earned (agents also surfaced genuine non-fatal nuances, not rubber-stamps).

## OPTIONAL close-out items (NON-fatal, NOT verdicts — cosmetic/metadata polish; defer to Step-6 at user discretion)
None of these is an overclaim that makes the program assert anything false; all are metadata/presentation nits.
1. **u09 V_0f dedup**: `fit_stage_065_v_0f` is a tokenization artifact (`V_0 f`) — mark as alias/dedup of V_0 rather than a standalone free_choice candidate (parallels batch-1's scanner-mislabel cleanups).
2. **u09 f'(0)=1**: documented only in the sympy script (line 118), absent from notes/ — optional notes mention (disclosure completeness, not an overclaim).
3. **u01**: card stage_001.tex:149 omits the word "ansatz"/"not yet frozen" the notes carry — under-disclosure asymmetry (bundles self-flag severity:low); optional card hedge.
4. **u08**: stage 070 establishes no own notes-level provenance for n=5 (inherits 062/Phase-0) — stale-anchor metadata; optional.
5. **u05**: Omega_A_l bundle `introduced_at_stage:17` — benign EM-extension renumber artifact.
6. **u11**: stage_073.tex 37/20 boxed under `\subsection*{Derivation}` — mild presentational tension, already self-flagged in the `fit_stage073_paper_fixes_geometry_ratios` bundle; optional reword.
7. **u12**: card boxes `ell/a=1/20` without an inline "posited" tag (surrounding "selected/declared/reference" framing carries disclosure); optional.
8. **u14**: provenance `downstream_dependents` (075->'062','063'; 076->'063','078') appear stale/pre-renumber (cards route values to 076/077/078) — metadata artifact.

## Overall
Batch 2 extends the layer-2 result into the foundations/early-geometry/EOS/geometry-freeze band: the free_choice ansatz
surface there is honestly disclosed as posited, with no posit-dressed-as-derived overclaim. NO program revision required;
the 8 optional items are cosmetic/metadata. The three ANSATZ_LEDGER branch-determinable values that originate in this band
(37/20 origin, epsilon_r=1/20, K_m) are all confirmed honestly-posited free_choice at their origin — consistent with the
batch-1 37/20 adjudication.
