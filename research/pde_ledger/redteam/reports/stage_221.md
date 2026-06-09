---
unit_id: 221
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [notes/stages/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 221 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_221.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows/blocks at lines 54, 369, 461-513, 1431)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_mathematica_audit.txt`

## What the paper claims

The card's `\stagefield{Output}` states the dispersive/no-free-lunch theorem: "$|\Re\chi_*|/|\Im\chi_*|=|\delta|/\gamma_*$, so resonant conservative gain is inseparable from absorptive loading, with maximal balanced leverage at equal conservative and absorptive magnitudes." The notes (section 9) and the Part VII appendix (lines 461-513) enumerate nine concrete deliverables: (1) the generic simple-pole Breit-Wigner normal form `chi_s ≈ A_*/(δ−iγ_*)` from the linearized denominator `F_Π = F_0'(ω_*)δ − Π(ω_*)Z_*`; (2) the Stage-220 wall derivative identity `∂_Π D_Π = −N` with `N = (A·G_W + R·G_U)²/Δ_Π²`; (3) the wall-like pole specialization `chi_qq ≈ (1/D_0')·1/(δ−iγ_wall)`; (4) the exact Re/Im line-shape forms and ratio `|Re|/|Im| = r`; (5) the resonance optimum `max|Re χ_*| = |A_*|/(2γ_*)` at `r=1`; (6) the low-loss factorization identity for `r/(1+r²) − η/(1+η²)`; (7) the barrier/absorbed-power ratio `P̄/(ω|U_disp|) = γ_*/|δ| = 1/r`; (8) the quality-factor detuning bound `|ω−ω_*|/ω_* ≥ 1/(2Q_*η)`; and (9) the linear survival window `(|A_j|/γ_*)·η/(1+η²)·S_j² ≥ 2ΔV_req`. The card status is `\StatusExactClosure{}`; the appendix flags only "whether a completed branch supplies a useful pole" as `\StatusOpen{}`, not any of the line-shape algebra.

Note on numbering: the notes file internally calls this "Stage 238" and attributes the derivative identity to "Stage 237," whereas the card, appendix, and both scripts call it Stage 221 and attribute the derivative identity to "Stage 220." This is a renumbering-era drift in the notes prose; the math is identical. I treat it as a low-severity notes-vs-paper label mismatch, not a math defect (see F4).

## What the script claims to verify

The SymPy script verifies, via 12 `assert simplify(...) == 0` checks: the simple-pole normal-form factoring (I), the passive substitution `Π→iΓ` form of `γ_*` (I), the Stage-220 derivative identity (II), the wall-like specialization factoring (III), the Re/Im decomposition and `|Re|/|Im|=r` on the positive-detuning slice `δ=rγ` (IV), the maximum-identity and low-loss factorization identities (IV.b), the barrier/absorbed-power ratio `=1/r` (V), and the two quality-factor detuning relations (VI). It then defines but only *prints* (no assertion) the low-loss `|U_disp|_max` and the necessary residue/linewidth ratio (the survival-window quantities), plus a probe-only numeric slice. The `.wl` mirrors the same 12 checks (as 12 `expectZero` calls) plus 2 additional `expectZero` calls (lines 160, 161-164) that purport to verify the survival window.

## Paper ↔ script cross-check

| Paper deliverable | SymPy | Mathematica | Status |
|---|---|---|---|
| (1) simple-pole Breit-Wigner normal form | L18-28 | L64-76 | match (algebraic factoring) |
| (2) `∂_Π D_Π = −N` derivative identity | L52 | L93 | match (substantive) |
| (3) wall-like specialization | L62 | L102-105 | match (algebraic factoring) |
| (4) Re/Im forms + `|Re|/|Im|=r` | L77-79 | L119-121 | match (substantive) |
| (5) optimum at r=1 (`max=A_*/2γ_*`) | L91 | L124-127 | partial (factorization identity only; sup/optimum value not asserted) |
| (6) low-loss factorization identity | L96 | L128-131 | match (substantive) |
| (7) barrier/absorbed-power ratio `=1/r` | L109 | L139 | match (substantive) |
| (8) Q-factor detuning bound | L124, L127 | L151, L152 | partial (equality forms asserted; the ≥ inequality direction not asserted) |
| (9) linear survival window inequality | L137-146 (print-only) | L160, L161-164 (tautological) | missing/mismatch |

Dominant pattern: most deliverables are faithfully exercised; deliverable (9) is unverified in both engines (print-only in SymPy, tautological in `.wl`), and the `.wl` is a transliteration. Hence `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 24 | `simplify(chi_lin - chi_bw)==0` | (1) | partial (denominator factoring) |
| A2 | sympy | 28 | `simplify(chi_passive - chi_passive_expected)==0` | (1) | partial (Π→iΓ factoring) |
| A3 | sympy | 52 | `simplify(diff(D_Pi,Pi)+Nfun)==0` | (2) | yes |
| A4 | sympy | 62 | wall-like factoring `==0` | (3) | partial |
| A5 | sympy | 77 | `simplify(re_r - re_expected)==0` | (4) | yes |
| A6 | sympy | 78 | `simplify(im_r - im_expected)==0` | (4) | yes |
| A7 | sympy | 79 | `simplify(re_r/im_r - r)==0` | (4) | yes |
| A8 | sympy | 91 | max-identity factorization `==0` | (5) | yes (identity), partial (sup) |
| A9 | sympy | 96 | low-loss factorization `==0` | (6) | yes |
| A10 | sympy | 109 | `simplify(ratio_barrier - 1/r)==0` | (7) | yes |
| A11 | sympy | 124 | detuning Q-form `==0` | (8) | yes (equality) |
| A12 | sympy | 127 | low-loss detuning boundary `==0` | (8) | yes (equality) |
| — | sympy | 137-146 | (no assertion; print-only) | (9) | NO |
| B1-B12 | mathematica | 71,76,93,102,119,120,121,124,128,139,151,152 | `expectZero` mirroring A1-A12 | (1)-(8) | same as A1-A12 |
| B13 | mathematica | 160 | `expectZero[survivalLeft - 2 UdispLowLossMax]` | (9) | NO — tautological |
| B14 | mathematica | 161-164 | `expectZero[residueRequirement·η/(1+η²)·S² - 2 DeltaVreq]` | (9) | NO — tautological |

Checkpoint PASS count: SymPy emits 12 asserts (all pass, exit 0). `.wl` emits 14 `expectZero` calls (line 29 is the definition, not a call) → 14 PASS lines in the output, all present. No missing/short PASS lines, but B13/B14 are tautological (see F1) so the effective substantive `.wl` count is 12, of which deliverable (9) has zero genuine coverage.

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_mathematica_audit.wl:154-164`

**What's wrong:**
The two `.wl` checks that purport to verify the linear survival window (paper deliverable #9) are algebraic round-trips of the immediately-preceding definitions and cannot fail for any physics. At L154-155:
```
UdispLowLossMax = normalizeExpr[(1/2) Aabs/gamma * eta/(1 + eta^2) * Sfam^2];
survivalLeft    = normalizeExpr[Aabs/gamma * eta/(1 + eta^2) * Sfam^2];
```
`survivalLeft` is defined as exactly `2*UdispLowLossMax`. The check at L160 `survivalLeft - 2 UdispLowLossMax` is therefore identically `0` by construction. Likewise at L156-158 `residueRequirement = 2 DeltaVreq/Sfam^2 * (1 + eta^2)/eta`, and the check at L161-164 `residueRequirement * eta/(1 + eta^2) * Sfam^2 - 2 DeltaVreq` just inverts that definition: `(2 DeltaVreq/Sfam²·(1+η²)/η)·η/(1+η²)·Sfam² − 2 DeltaVreq ≡ 0`. Neither check exercises any relation between the survival-window left side, the low-loss bound, and `ΔV_req`; they verify only that `x = 2·(x/2)` and `y·k = y·k`.

**Why this matters:**
Deliverable #9 — the central low-loss survival inequality `(|A_j|/γ_*)·η/(1+η²)·S_j² ≥ 2ΔV_req` — is the paper card's bottom-line gate ("the linear resonant corridor lives or dies on the residue-to-linewidth ratio"). These two PASS lines create the false impression the `.wl` covers it. As a checkpoint, this stage anchors downstream trust (MTDC-T11.1); a tautological "PASS" on the headline gate is a trust hazard.

**Required change:**
Replace the two tautological checks with substantive ones that tie the survival-window saturation to the low-loss line-shape bound derived earlier in the script (not to a re-statement of the same definition). Concretely: at the boundary `r = 1/η`, the bound `|U_disp| ≤ (1/2)(|A|/γ)·η/(1+η²)·S²` must follow from the actual line-shape `U_disp = -(1/2)·reExpected·Sfam²` evaluated with `reExpected = Aabs·r/(gamma(1+r²))` at `r → 1/η`. Verify `(-Udisp /. r -> 1/eta) - UdispLowLossMax == 0` (this is non-trivial: it requires `reExpected` at `r=1/η` to equal `(Aabs/gamma)·η/(1+η²)`). Then verify the survival inequality is *saturated* at equality by checking that substituting `Aabs/gamma -> residueRequirement` into `survivalLeft` yields exactly `2 DeltaVreq`, i.e. `(survivalLeft /. Aabs -> residueRequirement gamma) - 2 DeltaVreq == 0` — this connects `residueRequirement` to `survivalLeft` through the line-shape, not through its own definition.

**Verification:** New `expectZero` checks at the survival-window section should reference `Udisp`/`reExpected` (the earlier-derived line shape) rather than re-stating `UdispLowLossMax`/`residueRequirement`; the verifier confirms the checks involve the line-shape expression and that the `.wl` still exits 0.

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.py:136-146`

**What's wrong:**
The SymPy script computes the two survival-window quantities `Udisp_lowloss_max` (L137-139) and `residue_requirement` (L140-142) and only *prints* them (L144-146). There is no `assert` for paper deliverable #9 in the SymPy engine at all. So the headline survival-window gate has zero assertion in SymPy and only tautological assertions in the `.wl` (F1). On a checkpoint with `is_checkpoint: True`, the bottom-line `\stagefield{Output}` deliverable must be exercised by a real assertion in both engines.

**Why this matters:**
The card's claim status is `\StatusExactClosure{}` for the line-shape and "low-loss survival inequality" (appendix L513). A closure claim with no assertion behind it in either engine for that specific deliverable is unverified scaffolding masquerading as proof.

**Required change:**
Add a substantive SymPy assertion connecting the low-loss line-shape bound to the survival quantities. Mirror the non-tautological route prescribed in F1: assert `simplify((-U_disp).subs(r, 1/eta) - Udisp_lowloss_max) == 0` (this exercises that the general line shape `-U_disp = (1/2)·Aabs·r/(gamma(1+r²))·Sfam²` evaluated at the low-loss boundary `r=1/η` reproduces `Udisp_lowloss_max`), and assert that the residue requirement saturates equality through the line shape: `simplify((Aabs/gamma * eta/(1+eta**2) * Sfam**2).subs(Aabs, residue_requirement*gamma) - 2*DeltaV_req) == 0`.

**Verification:** New asserts appear after L142; the verifier confirms both reference the earlier line-shape symbol `re_expected`/`U_disp` (not a re-statement of `Udisp_lowloss_max`), and the SymPy script still exits 0.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_mathematica_audit.wl` (whole file; representative L64-76, L84-93, L112-131)

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent re-derivation. The variable choreography, intermediate definitions, ordering, and hand-supplied target forms are identical; only the syntax differs. Three corresponding sections:

Section I — both build the SAME intermediate `chi_bw`/`chiBW` rather than deriving the Breit-Wigner form independently:
- py L18-23: `F_lin = Fprime*delta - Pi*Z_star; chi_lin = together(Num_star/F_lin); A_star = Num_star/Fprime; chi_bw = together(A_star/(delta - Pi*Z_star/Fprime))`
- wl L64-69: `Flin = Fprime delta - portPi Zstar; chiLin = Together[Numstar/Flin]; Astar = ...; chiBW = Together[Astar/(delta - portPi Zstar/Fprime)]`

Section II — both write the SAME closed-form target `N = (A G_W + R G_U)²/Δ²` and check `∂D + N`, rather than the `.wl` independently computing `D[Q_Pi/DeltaPi, portPi]` and letting Mathematica discover/`FullSimplify` the residue form:
- py L48-52: `Q_Pi = ...; D_Pi = together(K_B - Q_Pi/Delta_Pi); Nfun = together((A*G_W + R*G_U)**2/Delta_Pi**2); assert diff(D_Pi,Pi)+Nfun==0`
- wl L89-93: `QPi = ...; DPi = Together[KB - QPi/DeltaPi]; Nfun = Together[(Afun GW + R GU)^2/DeltaPi^2]; expectZero[D[DPi,portPi] + Nfun]`

Section IV — both hand-write the same `re_expected`/`im_expected` and the same factorizations `(f-1/2)+(r-1)²/(2(1+r²))`, `(f-η/(1+η²))+(r-η)(ηr-1)/((1+r²)(1+η²))`:
- py L75-96 vs wl L116-131 — identical expressions, identical assertion order.

The only structural difference is the 2 extra (tautological, see F1) survival checks. There is no independent route in the `.wl`: no `Residue`/`Series`/`Limit` to obtain the pole reduction, no alternative decomposition of `N`, no independent computation of the line shape from `1/(δ−iγ)`. This violates the dual-engine policy for a checkpoint, which requires both engines to derive results independently.

**Why this matters:**
A transliterated second engine does not provide an independent confirmation; a transcription error in the `.py` target form would be copied verbatim into the `.wl` and both would "agree" on the wrong thing. For a trust-anchoring checkpoint this is the dual-engine policy's exact failure mode.

**Required change:**
Re-derive at least the load-bearing results in the `.wl` via native Mathematica primitives and a DIFFERENT decomposition (see directive F3 claim manifest). Keep the existing engine cross-check, but the `.wl` must reach the same results independently. Do NOT merely rename variables.

**Verification:** The verifier confirms the `.wl` derives (M1) the `N` form via `Together[D[QPi/DeltaPi, portPi]]` simplified to a perfect square WITHOUT pre-writing `(Afun GW + R GU)^2` as the target, (M2) the Breit-Wigner reduction via `Series`/`Normal` about `delta=0` (or `Limit`), and (M3) the line shape via `ComplexExpand` of `Aabs/(delta - I gamma)` then substitution `delta -> r gamma`, rather than the pre-collapsed `Aabs/(r gamma - I gamma)`; and that the `.wl` still exits 0.

### F4 — paper_misalignment (notes_contradicts_script)

**Severity:** low
**Files:**
- paper/notes side: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.md` (title + §3 + §9 repeatedly say "Stage 238" and "Stage 237 derivative identity")
- script side: `/var/projects/toy_physics/.../scripts/...sympy_audit.py:36-37` comment "Exact Stage 220 wall derivative identity"; `.wl:82` "II. Exact Stage 220 derivative identity"; card `paper/stages/stage_221.tex:9` "Imports Stage~220".

**What's wrong:**
The notes prose names this stage "Stage 238" and attributes the wall derivative identity to "Stage 237," while the card, appendix, and both scripts call it Stage 221 and attribute the identity to Stage 220. This is a renumbering-era inconsistency in the notes prose, not a math discrepancy (the identity `∂_Π D_Π = −N` and all formulas are the same on both sides). Per the auditor contract this is a `paper_misalignment` (notes prose vs. card/scripts) that I do not resolve — but it is low severity and purely a label.

**Why this matters:**
A reader cross-referencing the notes against the scripts will see "Stage 237/238" vs "Stage 220/221" and may suspect a deeper mismatch. For a checkpoint it is worth correcting the notes numbering, but this is a prose edit (notes/), routed to the user — Codex must not touch it.

**Required change (routed to user, NOT Codex):** see directive `## Resolve before fix_loop`. The math needs no change.

**Verification:** N/A for scripts; user decides whether to renumber the notes prose.

## Independent-derivation check (Mathematica)

The `.wl` is a transliteration of the `.py` (see F3). Every load-bearing definition (`Flin/chiBW`, `DPi/Nfun` with the pre-written square, `reExpected/imExpected`, the two factorization identities, `U_disp/P_abs`, the detuning substitutions) appears in the same order with the same hand-supplied target forms, differing only in Mathematica vs Python syntax. No native independent route (`Residue`, `Series`, `Limit`, alternative `N` decomposition) is used. Verdict: transliteration.

## Engine cross-check

Both engines produce identical symbolic results and identical probe numerics:
- `chi(omega)`: py `Num_star/(Fprime*delta - I*Gamma_out*Z_star)` ↔ wl `Numstar/(delta*Fprime - I*GammaOut*Zstar)`
- `A_*`, `gamma_*`, `gamma_wall`: identical.
- Re/Im, `|Re|/|Im|=r`, `f(r)`, max/low-loss identities, barrier ratio `1/r`, detuning `r/(2Qfac)` and `1/(2 Qfac eta)`, survival formulas: identical.
- Probe slice: Re=1.05, Im=0.35, ratio=3, |U|=13.125, P_abs=48.125, P/(ω|U|)=0.3333…, low-loss|U|=10.294…, required|A|/γ=0.204 — match to all printed digits.
No `engine_disagreement`. (Note this agreement is weakened as independent evidence by F3 — agreement of a transliteration is expected.)

## Verdict justification

The substantive line-shape algebra holds up under attack: the Stage-220 derivative identity (A3/B3) is a genuine non-trivial computation whose claimed closed form `N = (A G_W + R G_U)²/Δ²` I re-derived by hand and confirmed; the Re/Im decomposition and `|Re|/|Im|=r` (A5-A7), the max/low-loss factorizations (A8/A9), the barrier ratio (A10), and the detuning relations (A11/A12) are all real and correctly anchored. The Breit-Wigner/wall-factoring checks (A1/A2/A4) are low-content algebraic round-trips but acceptable as normal-form verifications. The defects are: (F1) the `.wl`'s two survival-window checks are tautological round-trips of their own definitions; (F2) the SymPy survival window has no assertion at all (print-only) — so the card's headline deliverable #9 is unverified in both engines; (F3) the `.wl` is a transliteration of the `.py`, not an independent derivation, which the checkpoint dual-engine bar forbids; (F4) the notes prose carries stale "Stage 237/238" numbering vs the card/scripts' "Stage 220/221" (low, prose-only, user-routed). Verdict: `findings`. Not `stop_cold`: the math is internally consistent and correct, and no downstream-propagating constant changes — F1/F2/F3 are added-coverage and independence fixes, F4 is a label.

## Self-test notes

I checked: (1) variable-independence — the only derivative is `∂_Π D_Π` (py L52 / wl L93); `D_Pi` genuinely depends on `Pi` through `Wfun = Omega_W²-omega²-Pi`, so the derivative is non-trivially nonzero and the identity is a real test (not an identically-zero trap). My prescribed new check `(-Udisp /. r->1/eta) - UdispLowLossMax` substitutes a value, no spurious derivative. (2) Trivial-case: for F1/F2 I confirmed `reExpected` at `r=1/η` equals `(Aabs/gamma)·(1/η)/(1+1/η²) = (Aabs/gamma)·η/(1+η²)`, so the proposed `(-Udisp at r=1/η) - UdispLowLossMax` reduces to 0 non-trivially (it would FAIL if the line shape were wrong, unlike the current tautology). (3) Sign/parity: positive-detuning slice `δ=rγ>0` with `A_*>0` gives `Im χ>0`, consistent with the paper's `A_*/(δ−iγ_*)`; no symmetric-domain integral here. (4) No fabricated load-bearing literals — only the explicitly probe-only numeric slice, matching across engines. (5) Paper round-trip: my prescribed fixes introduce no constant not already in the paper (`η`, `ΔV_req`, `A`, `γ`, `S_j` all appear in card/appendix deliverable #9).
