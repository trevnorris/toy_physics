---
batch: III.3
created: 2026-05-26
status: awaiting_codex_recommendations
items:
  - Q1: 062 F1 — second equality and Cauchy-Schwarz bound for boxed G_micro
  - Q2: 062 F2 — sign convention on σφ coupling in F_red / S_parent
---

# Batch III.3 paper-alignment resolutions

Two substantive `paper_misalignment` items, both on stage 062. The two other auditor-flagged `paper_misalignment` findings (067 F1, 072 F1) are pure script-banner relabels with unambiguous direction (paper card and filename are authoritative); they are routed directly to Codex as script-side label fixes and do not appear here.

The Codex recommendation pass writes `## Recommendation` blocks below each question with `direction:` and rationale. The user reviews and approves; then a separate Codex apply pass executes the approved direction.

---

## Q1 — Stage 062 F1: boxed `G_micro` second equality + Cauchy-Schwarz bound

**Auditor finding (full):** `redteam/reports/stage_062.md` F1 (subtype `script_missing_paper_claim`).

### Paper claim

`paper/stages/stage_062.tex:33-46` boxes eq. `app-stage062-Gmicro` with **two** equalities:

```
G_micro = ρ_* g_φ² O_{σφ}² / (m c_{s,*}² K_X N_{σσ})                          (i)
        = [ρ_* g_φ² N_{φφ} / (m c_{s,*}² K_X)] · C_{σφ}²                       (ii)
```

with `C_{σφ}² := O_{σφ}² / (N_{σσ} N_{φφ})` and `0 ≤ C_{σφ}² ≤ 1` (eq. `app-stage062-Csigma`).

### Script coverage

- Both engines verify the **first equality (i)** via on-shell elimination + closed-form check (sympy:86, wl:80).
- The **second equality (ii)** is only printed as a definition (sympy:88, wl:82). No assertion.
- The **Cauchy-Schwarz bound** `0 ≤ C_{σφ}² ≤ 1` is not exercised anywhere.

### Options

- **(a) Extend scripts** — add `C_sp_sq = Osp²/(Nss·Npp)` and a `G_micro_factored = (ρ_* g_φ² N_{φφ}/(m c_{s,*}² K_X)) · C_sp_sq` closed form; `expect_zero(G_micro_closed - G_micro_factored)`. For the bound, substitute `O_{σφ} = cos(θ)·√(N_σσ·N_φφ)` and confirm `C_sp² = cos²(θ) ∈ [0,1]`. Symmetric edit in `.wl`.
- **(b) Trim paper card** — drop the second equality and the Cauchy-Schwarz bound from `paper/stages/stage_062.tex`, leaving only `G_micro = ρ_* g_φ² O_{σφ}²/(m c_{s,*}² K_X N_σσ)`.
- **(c) Acknowledge** — document in script comments that the second equality is a paper-side rewrite using the definition of `C_{σφ}²` and the Cauchy-Schwarz bound is a standard inequality on overlap integrals; neither requires script-side verification. No assertion change.

### Destination-verification questions for Codex

- Do any downstream stages (063-072 in this batch, or any later stage) consume the `N_{φφ}·C_{σφ}²` factored form of `G_micro` rather than the `O_{σφ}²/N_{σσ}` form?
- Does any downstream stage import the Cauchy-Schwarz bound `0 ≤ C_{σφ}² ≤ 1` (e.g., the Cauchy no-go inequality in stage 063 — but verify whether 063's check actually relies on this bound being verified at 062, or whether 063 verifies the bound internally)?
- If downstream stages consume only the first equality form, (b) is justified.
- If downstream stages cite `C_{σφ}²` and the bound is load-bearing, (a) is justified.

## Recommendation

- direction: a
- rationale: Stage 063 is a direct downstream consumer of the factored coherence form and the Cauchy bound: its paper card lists as inputs "the parent gain formula" and "the Cauchy bound on \(C_{\sigma\phi}^2\)" at `paper/stages/stage_063.tex:6-7`, then boxes fixed-\(g_\phi\) coherence thresholds and a Cauchy no-go at `paper/stages/stage_063.tex:23-36`. The Stage 063 SymPy audit defines `C_sp_sq` at `scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py:35-36`, substitutes `Osp**2 -> C2 * Nss * Npp` in the threshold checks at `scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py:50-60`, asserts the factored relation `G_micro - G_max*C^2` at `scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py:69-73`, and uses Cauchy saturation at `scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py:117-121`; the Mathematica audit mirrors this with `c2Def = oSP^2/(nSS*nPP)` and `G_micro - G_max C^2` at `mathematica/moving_throat_pde_stage063_parent_thresholds_mathematica_audit.wl:41-58` plus Cauchy saturation at `mathematica/moving_throat_pde_stage063_parent_thresholds_mathematica_audit.wl:101-104`. Because the disputed Stage 062 second equality and bound are load-bearing in the next stage, Stage 062 should verify them rather than merely print the definition at `scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:88` and `mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:82`.
- downstream_impact: Stage 063 directly consumes the factored `N_{\phi\phi} C_{\sigma\phi}^2` form and Cauchy no-go; Stage 064 independently rederives an equilibrium \(C_{\sigma\phi}^2\le1\) profile bound; Stages 067-068 consume later benchmark/generic coherence penalties rather than the raw Stage 062 overlap identity.
- notes: Stage 063 does check the Cauchy saturation downstream, but that does not close the missing verification on the Stage 062 paper-card claim itself.

## Apply log

- applied_at: 2026-05-26T18:40:40Z
- direction: a (Q1)
- files_changed:
  - scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py
  - mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl
- sympy_exit: 0
- mathematica_exit: 0
- summary: Added script-side verification for the boxed G_micro factored equality and the Cauchy-Schwarz parameterization in both engines.

---

## Q2 — Stage 062 F2: σφ coupling sign in F_red / S_parent

**Auditor finding (full):** `redteam/reports/stage_062.md` F2 (subtype `notes_contradicts_script`).

### Notes claim

`notes/stages/moving_throat_pde_stage062_parent_action_gain.md:44`:
```
F_red[σ,φ] = ∫₀^L ds [ (Θ_σ/2)σ² − Λ_φ σ φ + (T_X/2)φ_s² + (K_X/2)φ² ] + (K_m/2)φ(0)²
```

Note the **minus** sign on the `σφ` coupling — derived from the physical perturbation `δV_conf = −g_φ χ_φ(y) φ(s)`.

### Script implementation

Both engines build `S_parent` with **plus** sign on the coupling:
- sympy:71-75: `S_parent = ... + Λ_φ σ φ + ...`
- wl:65: `sParent = ... + λ_φ σ φ + ...`

### Effect

- The gain `G_micro = Λ²/(Θ K_X)` is invariant under this sign flip (quadratic in `Λ`), so all assertions pass.
- On-shell `σ_*(φ)` prints with opposite sign:
  ```
  script:  σ_* = −O_sp g_φ φ ρ_*/(N_ss c_{s,*}² m)
  notes:   σ_* = +O_sp g_φ φ ρ_*/(N_ss c_{s,*}² m)
  ```
- If a downstream stage imports `σ_*(φ)` itself (not just the gain), it picks up the wrong sign.

### Options

- **(a) Script-side fix** — change `S_parent` / `sParent` to use `−Λ_φ σ φ` to match notes. Gain assertion still passes; printed `σ_*` flips sign to match notes.
- **(b) Paper/notes-side fix** — update notes to use the script's `+Λ_φ σ φ` convention, with a remark on the physical derivation route. Requires editing `notes/stages/moving_throat_pde_stage062_parent_action_gain.md`.

### Destination-verification questions for Codex

- Does any downstream stage (063-072 + later) import the **expression** `σ_*(φ)` from stage 062 (not just the gain `G_micro`)?
- Stage 064 ("equilibrium source profile") in particular — does it construct `σ_*(φ)` from notes-derived premises and would it conflict with the script's sign convention?
- If any downstream consumer reads `σ_*(φ)` with the notes' sign convention, (a) is required.
- If `σ_*(φ)` is purely a stage-062-internal display quantity (gain alone propagates), either (a) or (b) is fine — but (a) is cheaper.

## Recommendation

- direction: a
- rationale: A search for the exact Stage 062 on-shell expression found no downstream import of `sigma_star = ±O_sp*g_phi*phi*rho_star/(N_ss*cs_star_sq*m)` outside the Stage 062 audit report, so the printed \(\sigma_*(\phi)\) is not itself propagated. However, Stage 064 is sign-sensitive and follows the notes' minus-coupling convention: the paper boxes the positive equilibrium source profile \(\chi_\sigma(y)=g_\phi\chi_\phi(y)/H(y)\) at `paper/stages/stage_064.tex:13-19`; the SymPy audit derives it from `F_loc = (1/2) H sigma^2 - g_phi chi_phi sigma` and checks the closure at `scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py:56-68`, then uses the same reduced static sign `F = ... - Lambda * sigma * phi + ...` at `scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py:143-148`; the Mathematica audit mirrors those premises at `mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl:33-42` and `mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl:110-121`. Since Stage 062 currently builds `S_parent` with the opposite plus sign at `scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:71-76` and `mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:63-66`, the cheaper and more aligned repair is to flip the scripts to the notes/Stage 064 convention; the gain checks remain invariant.
- downstream_impact: No downstream stage appears to import the literal Stage 062 \(\sigma_*(\phi)\) expression, but Stage 064 constructs the equilibrium source profile with the notes' sign convention and would remain conceptually inconsistent with a plus-sign Stage 062 action.
- notes: Apply pass should update the Stage 062 scripts and transcripts only; the paper/notes sign convention is already the one used by Stage 064.

## Apply log

- applied_at: 2026-05-26T18:40:40Z
- direction: a (Q2)
- files_changed:
  - scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py
  - mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl
- sympy_exit: 0
- mathematica_exit: 0
- summary: Flipped the Stage 062 parent-action sigma-phi coupling to the notes' minus-sign convention in both engines while preserving the invariant gain checks.
