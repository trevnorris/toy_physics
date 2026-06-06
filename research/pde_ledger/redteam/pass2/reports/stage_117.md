---
unit_id: 117
batch: IV.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage117_outlet_core_status.md]
  paper_appendix: present
---

# Audit unit 117 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_117.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage117_outlet_core_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (row at L1268 `\input{stages/stage_117}`; load-bearing narrative at L375-565, L1073-1077)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage117_outlet_core_status_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage117_outlet_core_status_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.txt`

## What the paper claims

Stage 117 ("Concrete Outlet-Core Status") is a status-consolidation ledger card. It has no explicit `\stagefield{Output}` equation; the card carries `\claimstatus{\StatusExactClosure / \StatusOpen}`, names the step a "core outlet realization ledger step," and lists three diagnostic Checks (Schur-complement signs; D/N mixed-tube length normalization vs `L/a`; parent overlap ratios before interpreting them as branch values). The authoritative source is the notes file plus the appendix block `subsec:app-part04-compensated-robin-mixed` / `subsec:app-part04-core-balance`. The notes state the concrete shell/mixed core lands on the canonical compensated outgoing `l=2` branch iff (i) the coupling-balance law `g_s²(K_sK_q+λ²)=4(K_sg_q−λg_s)²` holds, (ii) the auxiliary D/N tube has length `L_W=(πa/2)√((1+r_c)/3)` with `r_c=λ²/(K_sK_q)`, and that on that surface the outlet reduces to `Λ_eff=Λ₂^out+4σ_*−σ_*/(1−z²/3−iz⁵/9)` with `σ_*=g_s²/(4K_s)`, preserving `Ŷ₂(z)=1+z²/9+4z⁴/81+iz⁵/27+O(z⁶)`. The appendix additionally records the full deformation classification: scale class (β=±1 only), pure Robin (`χ_Q^R=3/(3−ρ_R)`, canonical only at ρ_R=0), standalone mixed pole (`κ_W=−1/9`, then `σ_W=0`), and the two hybrid even-solutions `(ρ_R,κ_W)=(σ_W,0)` and `(4σ_W,1/3)` with `χ_Q^hyb=(1−9σ_Wγ_W)/(1−σ_W)` canonical iff `γ_W=1/9`. The surviving open question is microscopic: does the actual moving-throat core realize this surface? (The notes/card do not claim it does — `\StatusOpen`.)

Note: the card's `\stagefield{Purpose}` text reads "Stage~134 is a core outlet realization ledger step" — a stale self-label (117+17=134, the known +17 pre-renumber drift; neighbors 118→135, 119→136, 120→137, 121→138, 122→139, 123→140 are all exactly +17). This is a cosmetic numbering-band artifact already scoped to the dedicated numbering pass, not a math-side finding. Recorded as a side observation only.

## What the script claims to verify

The script (SymPy + a mirrored Mathematica `.wl`) runs a six-section low-frequency classification of outlet deformations: (1) the scale/argument class admits only β=±1, and the β=1 branch preserves χ_Q=1; (2) pure Robin preserves the canonical even fingerprint only at ρ=0; (3) a standalone mixed pole forces a formal κ=−1/9 and then σ=0 (the pole disappears); (4) the hybrid outlet has exactly two even-matching branches — the trivial cancellation `(ρ,κ)=(σ,0)` and the nontrivial compensated `(4σ,1/3)` — with the compensated branch collapsing to `(1−σ)Λ_out` at γ=1/9; (5) the concrete two-channel shell/mixed core obeys `ρ_c=4σ_c` on the balance surface with `σ_c|_balance=σ_*=g_s²/(4K_s)`, and routes κ₀ through the stage-116 FORWARD tube-length law `L_W=πa√((1+r_c)/3)/2 ⟹ κ₀=4L_W²/(π²a²)`, with γ₀=(1+r_c)/9 carried as a postulated ansatz, so that `δΛ_core` collapses to the compensated form; (6) a capstone wires booleans from sections 1-5 residuals and asserts the unique nontrivial survivor is the compensated Robin-mixed core realization. The script does NOT assert the microscopic realization (it prints the open question).

## Paper ↔ script cross-check

| Paper / notes / appendix deliverable | Script-side check | Status |
|---|---|---|
| Scale class β=±1, β=1 preserves χ_Q (app L375-397) | §1 `beta_solutions=={±1}`; `chi_arg(β=1)−1=0` | match |
| Pure Robin canonical only at ρ_R=0 (app L402-412) | §2 `robin_even_solutions==[{ρ:0}]`; `chi_R(ρ=0)−1=0` | match |
| Standalone mixed pole κ_W=−1/9 then σ_W=0 (app L414-423) | §3 `kappa_match+1/9=0`; `sigma_match=0`; `chi_mix(σ=0)−1=0` | match |
| Hybrid even-solutions (σ,0) & (4σ,1/3) (app L437-443) | §4 `hybrid_solutions=[{κ:0,ρ:σ},{κ:1/3,ρ:4σ}]` | match |
| `χ_Q^hyb=(1−9σγ)/(1−σ)`, canonical iff γ=1/9 (app L444-456) | §4 `chi_cancel−(1−9σγ)=0`; `chi_comp(γ=1/9)−1=0`; collapse to `(1−σ)Λ_out` | match |
| Coupling-balance `g_s²(K_sK_q+λ²)=4(K_sg_q−λg_s)²` (notes L33-37; app L524-528) | §5 `rho_c−4σ_c=0` solved for g_q (two roots) | match |
| `σ_*=g_s²/(4K_s)` (notes L51; app—via σ_c on balance) | §5 `σ_c|_balance−σ_*=0` (both roots) | match |
| `L_W=(πa/2)√((1+r_c)/3)`, κ_c=1/3 (notes L36-39; app L514-521, L544-548) | §5 `kappa0_from_tube=4L_W_forward²/(π²a²)`; load-bearing residual `δ_core−δ_core_exp=0` | match (de-tautologized; see Verdict) |
| `Λ_eff=4σ_*−σ_*/(1−z²/3−iz⁵/9)` (notes L42-51; app L497-509) | §5 `delta_core−delta_core_expected=0` at O(z⁶) | match |
| γ_c=1/9 via γ₀=(1+r_c)/9 (notes; app L508, L1077) | §5 γ₀ carried as ANSATZ; feeds the §5 residual (no separate tautological assertion) | match (ansatz-conditioned) |
| Capstone: unique nontrivial survivor (status tag) | §6 booleans wired from §1-5 residuals; `survivors==["compensated…"]` | match |
| Schur-complement sign Check (card Checks item 1) | §5 implicit via squared `(K_sg_q−λg_s)²`; no explicit positivity test | partial |
| Parent overlap ratios Check (card Checks item 3) | not tested in script | missing (diagnostic-only) |

Paper alignment: **aligned** on every load-bearing deliverable (the notes/appendix equations are all exercised by non-tautological checks). The two unexercised items (Schur sign positivity, parent overlaps) are diagnostic Checks-list reminders, not `\stagefield{Output}` deliverables, and were explicitly out of scope in the pass-1 directive; they do not pull alignment below `aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 64-65 | `if {beta sols} != {±1}: raise` | scale class structure | yes |
| A2 | sympy | 66 | `expect_zero(chi_arg(β=1)−1)` | scale class odd norm | yes |
| A3 | sympy | 80-81 | `if robin_even != [{ρ:0}]: raise` | pure Robin uniqueness | yes |
| A4 | sympy | 82 | `expect_zero(chi_R(ρ=0)−1)` | pure Robin odd norm | yes |
| A5 | sympy | 98 | `expect_zero(kappa_match+1/9)` | standalone pole κ=−1/9 | yes (solved, non-trivial) |
| A6 | sympy | 99 | `expect_zero(sigma_match)` | standalone pole vanishes | yes |
| A7 | sympy | 100 | `expect_zero(chi_mix(σ=0)−1)` | trivial at σ=0 | yes |
| A8 | sympy | 122 | `expect_zero(chi_cancel−(1−9σγ))` | cancellation branch odd-norm | yes |
| A9 | sympy | 123-126 | `expect_zero(... gamma=0 trivial)` | cancellation trivial at γ=0 | yes |
| A10 | sympy | 127 | `expect_zero(chi_comp(γ=1/9)−1)` | compensated branch odd norm | yes |
| A11 | sympy | 128-138 | `expect_zero(... collapses to (1−σ)Λ_out)` | compensated == pure scale | yes |
| A12 | sympy | 164-167 | `expect_zero(σ_c(root0)−σ_c(root1))` | both g_q roots → same σ | yes |
| A13 | sympy | 168 | `expect_zero(σ_c(root0)−σ_*)` | σ_c on balance = g_s²/(4K_s) | yes |
| A14 | sympy | 190 | `expect_zero(delta_core−delta_core_expected)` | **load-bearing: full Λ_eff collapse + κ_c=1/3 via forward L_W** | yes (falsifiable; de-tautologized) |
| A15 | sympy | 229-230 | `if survivors != ["compensated…"]: raise` | capstone uniqueness | yes (flags computed from §1-5) |
| B1-B15 | mathematica | 41-178 | mirrors A1-A15 1:1 | same set | mirror (see F1) |

All SymPy rows are "yes." There is no longer a standalone tautological κ_c=1/3 / γ_c=1/9 assertion (the pass-1 F2/F3 `Lw_required` back-solve and the bare-γ identity have been removed/folded into A14). The capstone flags (A15) are computed, not hardcoded (pass-1 F4 fixed).

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl:34-178`
- (parallel) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage117_outlet_core_status_sympy_audit.py:50-230`

**What's wrong:**
The `.wl` is a line-by-line transliteration of the `.py`, not an independent second-engine derivation. Corresponding sections, in identical order, with identical choreography and verbatim check-label strings:

- Generating function — sympy L50 `Lambda_out = -3 + z**2/3 + z**4/9 + I*z**5/9`; wl L34 `lambdaOut = -3 + z^2/3 + z^4/9 + I z^5/9` — identical.
- §1 — sympy L54 `Y_arg = sp.series((-3*S)/(S*Lambda_out.subs(z, beta*z)), z, 0, 6)`; wl L37 `yArg = Series[(-3 s)/(s (lambdaOut /. z -> beta z)), {z,0,5}]`; both then `sp.solve([m2==1/9, m4==4/81], beta)` ↔ `Solve[{m2Arg==1/9, m4Arg==4/81}, beta]` — same solve, same RHS literals.
- §4 — sympy L103 `Lambda_hyb = sp.series(Lambda_out + rho - sigma/(1 - kappa*z**2 - I*gamma*z**5), …)`; wl L75 `lambdaHyb = Series[lambdaOut + rho - sigma/(1 - kappa z^2 - I gamma z^5), …]` — identical.
- §5 (the heart) — sympy L156-160 define `r_c, rho_c, sigma_c, kappa_c, gamma_c`; wl L109-113 define `rC, rhoC, sigmaC, kappaC, gammaC` with byte-identical formulas. The `.wl` does **not** derive these from the 2×2 shell/mixed Schur system (no `Solve` of the core matrix, no explicit `adj/det` elimination of the side-channel `q`) — it re-types the SymPy script's already-reduced closed forms `(K_s g_q − λ g_s)²/(K_s² K_q (1+r_c))` etc.
- §5 F2 path — sympy L177-179 `L_W_forward`, `kappa0_from_tube`; wl L128-130 `lWForward`, `kappa0FromTube` — identical construction.
- Every `expect_zero(...)` label string equals its `expectZero(...)` counterpart verbatim, and §6's classification rows/`nontrivial*` wiring mirror 1:1.

The only "independence" is that two CAS each Taylor-expand the same scalar rational functions and solve the same matching equations. A shared algebraic/transcription error in a reduced closed form (e.g. in `sigma_c` or in the hand-written `delta_core_expected`) would be replicated identically in both files and pass both engines.

**Why this matters:**
This is exactly the defect the IV.2 transliteration-watch band (105-175) is being re-scrutinized for. The IV.2 batch precedent is direct and recent: stage 114's `.wl` re-typed the same family of reduced `ρ_c/σ_c/κ_c/γ_c` closed forms via a mirrored matrix-inverse, and the user AUTHORIZED re-authoring the `.wl` to a genuinely independent route (Schur elimination by `Solve`-ing the 2×2 core system) rather than accepting the mirror. Stage 117 §5 carries the identical re-typed reduced forms and is feasible to fix the same way — derive `δΛ_core` in Mathematica by independently eliminating the side-channel `q` from the core matrix `{{K_s,λ},{λ,−K_q D_W}}` with `D_W=1−κ₀z²−iγ₀z⁵`, then assert the same `delta_core_expected`. This converts the second engine from "Mathematica's `Series`/`Solve` agree with SymPy's" into an independent corroboration of the Schur reduction.

(Pass-1 raised this same finding F1 and it was dispositioned "accepted policy-mirror" under the older Mathematica-mirror policy. The IV.2 calibration superseded that for this band: line-by-line ports where an independent primitive is feasible are now flagged for re-authoring, not accepted. Direction confirmed by the stage-114 IV.2 directive. Severity medium: the emitted values are all correct, so this is a second-engine-strength defect, not a correctness defect.)

**Required change:**
Re-author the `.wl` §5 so the core correction `δΛ_core(z)` is derived by an independent route — explicitly eliminate the side-channel coordinate `q` from the core system `{{K_s,λ},{λ,−K_q D_W(z)}}·{s,q}=u·{g_s,g_q}` (e.g. `Solve` the second row for `q`, back-substitute into the first for `s`, then `δΛ_core = g_s s + g_q q` at `u=1`), with `D_W(z)=1−κ₀z²−iγ₀z⁵`. Keep `kappa0FromTube`/`lWForward` (the §5 F2 forward law) and the same final `deltaCoreExpected`; keep all `expectZero` labels and the §1-4/§6 structure. The acceptance criterion is that the `.wl` reaches `delta_core_expected` through an algebraically distinct path (a `Solve`/elimination on the core matrix, not a re-typed reduced `sigmaC`), so a wrong reduced closed form cannot pass both engines. (Codex designs the exact route; this states requirement + acceptance only.)

**Verification:**
After the rewrite, the `.wl` §5 must contain a `Solve[...]` (or explicit `adj/det`) elimination of `q` from the 2×2 core matrix in place of the re-typed `sigmaC = (ks gq - lam gs)^2/(ks^2 kq (1+rC))` line (wl L111). The variable names / intermediate steps must differ from the `.py`'s; all `expectZero` checks must still print `= 0` / `PASS`; the printed core-balance branches, `σ_*`, and the final `concrete core collapses…` residual must be byte-identical to the current transcript. Re-run `math -script` to confirm `Exit[0]`.

## Independent-derivation check (Mathematica)

Not independent — see F1. The `.wl` mirrors the `.py` step for step across all six sections: identical `lambdaOut` definition, identical `Series[…,{z,0,5}]` expansions of the same scalar generating functions, identical `Solve` matching equations with the same RHS literals (`1/9`, `4/81`), identical re-typed §5 reduced forms (`rC, rhoC, sigmaC, kappaC, gammaC`), identical `lWForward`/`kappa0FromTube` forward-law construction, identical `deltaCoreExpected`, and identical §6 boolean wiring. Every check-label string matches verbatim. The Mathematica run therefore corroborates only "Mathematica's `Series`/`Solve`/`FullSimplify` agree with SymPy's on the same algebra," not the physics by an independent path.

## Engine cross-check

Both engines exit 0 and agree at the level they claim (they are mirrors, so engine_disagreement is n/a). Every asserted residual prints `= 0` in both transcripts:
- §1 `scale/argument solutions = {{β→−1},{β→1}}` (both); `chi_Q=1` residual `= 0` (both).
- §3 `kappa match = -1/9`, `sigma match = 0` (both); all three residuals `= 0`.
- §4 `hybrid branches = {{κ:0,ρ:σ},{κ:1/3,ρ:4σ}}` (both); all four residuals `= 0`.
- §5 core-balance branches `g_q = (g_s/2K_s)(2λ ± √(K_qK_s+λ²))` (both); `σ_*` residuals `= 0`; load-bearing `concrete core collapses… = 0` (both).
- §6 capstone: both print exactly one nontrivial survivor `compensated Robin-mixed core realization (True,True,True)`.
The Mathematica log ends `Stage 117 Mathematica audit passed.` (12 `PASS:` lines, reconciling with the 12 SymPy `expect_zero` assertions per the pass-1 verification count).

## Verdict justification

**Verdict: findings (1, medium).** I re-confirmed — independently, not by assuming the pass-1/verification record — the three historically flagged items, all of which hold up in the current scripts:

1. **κ₀ is NOT target-inverted (pass-1 F2 genuinely fixed).** The `Lw_required = solve(4Lw²/(π²a²)=(1+r_c)/3, Lw)` back-solve is GONE from both engines (the X−X pitfall variable is absent — no `Lw_required`/`lWRequired` anywhere). κ₀ is built FORWARD: `L_W_forward = πa√((1+r_c)/3)/2 ⟹ kappa0_from_tube = 4L_W_forward²/(π²a²)`, which simplifies to `(1+r_c)/3`, so `κ_c=1/3`. Adversarial falsifiability check (mine): if the law were `√((1+r_c)/2)`, `kappa0_from_tube=(1+r_c)/2 ⟹ κ_c=1/2`, and the O(z²) pole of `δ_core` would no longer match `δ_core_expected`'s `z²/3`, so the load-bearing A14 residual would be nonzero and FAIL. The check genuinely exercises (and could falsify) the in-script tube-length coefficient — acceptable de-tautologization for a consolidation card (it is a consistency check against the written forward law, not an in-stage re-solution of the D/N half-wave BVP, which lives upstream at stage 116).

2. **γ₀=(1+r_c)/9 is a POSTULATED pure-scale ANSATZ, correctly labeled (pass-1 F3 fixed).** The §5 comment (sympy L148-153 / wl L98-108) and the status print (sympy L180 / wl L131) explicitly say γ₀ is "a pure-scale ANSATZ of the canonical compact outgoing l=2 branch … not derived," and that `γ_c=1/9` is "a consistency-of-assumption check, not an independent derivation." The stale "Stage 119" / "derived upstream" attribution is gone (no `119` in either script). There is no standalone tautological `γ_c=1/9` assertion; γ₀ feeds only the A14 residual. This ansatz is flagged below for the catalog.

3. **Capstone flags are computed, not hardcoded (pass-1 F4 fixed).** §6 `even_ok_*`/`odd_ok_*`/`nontrivial_*` reference genuinely-upstream §1-5 variables (`beta_solutions`, `robin_even_solutions`, `kappa_match`, `sigma_match`, `chi_mix`, `hybrid_solutions`, `chi_cancel`, `chi_comp`, and the §5 residual for `nontrivial_compensated`); the capstone yields exactly one survivor.

Other attacks that failed to produce findings: (a) §4 `χ_comp` algebra — `(1−9σγ)/(1−σ)` at γ=1/9, σ-free in the canonical limit, evaluates to 1 ✓; (b) §5 `σ_c|_balance` — both g_q roots give `(K_sg_q−λg_s)²=g_s²(K_sK_q+λ²)/4`, so `σ_c=g_s²/(4K_s)=σ_*` ✓; (c) coupling-balance form `g_s²(K_sK_q+λ²)=4(K_sg_q−λg_s)²` matches the notes/appendix verbatim ✓; (d) symbol domains — `Ks,Kq,lam,kappa0,gamma0,a,Lw` positive, `rho,sigma,kappa,gamma,gs,gq` real-only; positivity is justified by the physical setup and does not enable a false simplification (the residuals are identically-zero rational functions); (e) value reconciliation — every emitted deliverable maps to the notes and/or appendix (see below), zero misaligned.

The lone defect is F1: the `.wl` is a verbatim transliteration of the `.py`, and §5 re-types the reduced Schur closed forms rather than deriving them independently — precisely the IV.2-band defect that was re-flagged and re-authored at stage 114. No `stop_cold`: F1 changes no emitted value (the re-authored route reaches the same residuals), so no downstream propagation; `material_change` would be false.

## Self-test notes

Checked: (1) **variable independence** — the F1 required change introduces a `Solve`-based elimination of `q` from the 2×2 core matrix, which legitimately depends on `s,q,K_s,K_q,λ,g_s,g_q,D_W(z)`; no `diff`/`D[expr,var]`-against-absent-variable trap is introduced (no derivatives in the prescribed route). (2) **symmetry/parity** — no integrals over unbounded domains anywhere in this unit; n/a. (3) **trivial-case** — the prescribed independent `δΛ_core = g_s s + g_q q` must reduce to the same `ρ_c − σ_c/(…)` rational form (verified mentally via the 2×2 Cramer solution: `s,q` are rational in `D_W`, and the mouth feedback reproduces the Schur complement), so the existing `expectZero` residuals stay exact zeros. (4) **path specification** — F1 targets the existing `.wl` under `mathematica/`; path verified, no new file. (5) **paper round-trip** — the F1 rewrite keeps `lWForward`/`kappa0FromTube` and `deltaCoreExpected` unchanged and changes only the derivation route of `σ_c/δΛ_core`, so it introduces no new constant and cannot create a fresh `paper_misalignment`; every emitted value stays byte-identical.

## Value Reconciliation (pass-2 augmentation)

Outputs are FRESH (sympy `.txt` mtime 2026-05-29 13:25:44 > `.py` 13:23:26; mathematica `.txt` 13:25:44 > `.wl` 13:20:34), so reconciliation rests on script source + committed transcripts. All §-emitted deliverables are symbolic / pure-rational; the only literal constants are the canonical fingerprint rationals (1/3, 1/9, 4/81, etc.), which are physics targets the notes/appendix state, not free hardcodes.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| scale class roots `β = ±1` | py L58-63 / wl L41 (out L9) | app L375-397 (β=1 preserves) | MATCH |
| pure-Robin canonical only at `ρ_R=0` | py L74-78 / wl L52 (out L15) | app L402-412 (`χ_Q^R=3/(3−ρ_R)`) | MATCH |
| standalone mixed pole `κ_W=−1/9` | py L90,98 / wl L63,70 (out L21,23) | app L423 ("force κ_W=−1/9") | MATCH |
| standalone mixed pole `σ_W=0` | py L91,99 / wl L64,71 (out L22,24) | app L423 ("then σ_W=0") | MATCH |
| hybrid even-solutions `(σ,0)` & `(4σ,1/3)` | py L108-118 / wl L80-83 (out L30/L35) | app L437-443 (eq hybrid-even-solutions) | MATCH |
| `χ_Q^hyb=(1−9σγ)/(1−σ)`; canonical iff `γ_W=1/9` | py L120-127 / wl L84-91 | app L444-456 (eq hybrid-chiQ, gammaW-canonical) | MATCH |
| compensated branch `→ (1−σ)Λ_out` | py L128-138 / wl L92-95 (out L34) | notes L43-55; app L456 ("collapses to … pure scale") | MATCH |
| coupling-balance `g_s²(K_sK_q+λ²)=4(K_sg_q−λg_s)²` / `g_q=(g_s/2K_s)(2λ±√(K_qK_s+λ²))` | py L161 (out L39) / wl L114 (out L48) | notes L33-37; app L524-528 (eq core-coupling-balance) | MATCH |
| `σ_* = g_s²/(4K_s)` | py L163,168 / wl L116,121 (out L41) | notes L51; app L502/L517 (via `ρ_c=4σ_c`) | MATCH |
| `r_c = λ²/(K_sK_q)` | py L156 / wl L109 | notes L37-39; app L490 (eq rc-def) | MATCH |
| `ρ_c=g_s²/K_s`, `σ_c=(K_sg_q−λg_s)²/(K_s²K_q(1+r_c))`, `κ_c=κ₀/(1+r_c)`, `γ_c=γ₀/(1+r_c)` | py L157-160 / wl L110-113 | app L500-509 (eq core-identifications) | MATCH |
| `L_W = πa√((1+r_c)/3)/2` (forward law) → `κ₀=(1+r_c)/3` → `κ_c=1/3` | py L177-179 (out L42) / wl L128-130 (out L53) | notes L36-39; app L544-548 (eq LW-compensation), L1075 (`κ₀=(1+r_c)/3`) | MATCH |
| `γ₀=(1+r_c)/9` (ANSATZ) → `γ_c=1/9` | py L180,188 (out L43) / wl L131,135 (out L54) | notes (compensated form, L43-55); app L508, L1077 (`γ₀=(1+r_c)/9`) | MATCH |
| `Λ_eff/δΛ_core = 4σ_*−σ_*/(1−z²/3−iz⁵/9)` (load-bearing collapse) | py L189-190 (out L44) / wl L136-139 (out L55) | notes L42-51; app L497-509 | MATCH |
| capstone: unique nontrivial survivor = compensated core realization | py L218-230 (out L49-53) / wl L167-178 (out L61-65) | app claimstatus L37; narrative L78-79 | MATCH |

INTERNAL scaffolding (accounted for, no finding): `Lambda_out/lambdaOut`, the per-section series coefficients (`m2_arg, m4_arg, c2_R, c4_R, L0/L2/L4/L5_mix/_hyb`), `chi_arg/chi_R/chi_mix/chi_cancel/chi_comp`, `kappa_match/sigma_match` intermediates, `gq_solutions` list, the `kappa0_from_tube` intermediate, and all `even_ok_*/odd_ok_*/nontrivial_*` boolean flags.

Postulated ansatz value: **γ₀ = (1+r_c)/9** — pure-scale ANSATZ of the canonical compact outgoing l=2 branch (NOT derived; postulated in the stage-116 note "Bare outgoing normalization," carried as a hardcoded input at stages 115/116). Correctly labeled as such in both scripts (sympy L150-153, L180 / wl L104-107, L131). Flag for the ansatz catalog. (κ₀=(1+r_c)/3, by contrast, IS forward-derived upstream at stage 116 from the D/N half-wave eigenvalue and is exercised here via the forward tube-length law — not an ansatz.)

reconciliation: complete; 16 deliverable values checked, 0 misaligned.
