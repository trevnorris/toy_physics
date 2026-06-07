---
unit_id: 129
batch: IV.4
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage129_mouth_boundary_layer.md]
  paper_appendix: present
---

# Audit unit 129 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_129.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage129_mouth_boundary_layer.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (line 1292 is just `\input{stages/stage_129}`; the appendix carries no separate summary row — the card IS the appendix content)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage129_mouth_boundary_layer_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage129_mouth_boundary_layer_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage129_mouth_boundary_layer_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage129_mouth_boundary_layer_mathematica_audit.txt`

## What the paper claims

The `.tex` card (terse) states the audit target via a quote: "Electrochemical zero-flux law yields the exponential source profile \(\sigma_\Pi\)." The notes carry the full derivation: starting from a minimal positive source free-energy functional `F_mouth[σ]` (notes l.26-35), a linearized mouth potential `V_m(z) ≈ V₁z` with `V₁ := ∂_z δV_conf|_m − q_* ∂_z A₀|_m` (l.47-55), and an electrochemical potential `μ_σ^chem(z) = Θ_σ ln(σ/σ_*) + V₁z` (l.64-70), the standard positive Onsager current `J_σ = −M_σ σ ∂_z μ_σ^chem = −M_σ(Θ_σ σ' + V₁σ)` (l.77-89) is set to zero on the stationary recirculating branch, yielding the ODE `Θ_σ σ' + V₁σ = 0` and hence `σ = C e^{−V₁z/Θ_σ}`. Normalizing `∫₀^L σ dz = 1` gives the boxed exact positive family `σ_Π(z) = Π e^{−Πz/L} / (L(1−e^{−Π}))` with `Π := V₁L/Θ_σ > 0` (l.113-119), and the single remaining shape parameter `Π_m = V₁L/Θ_σ` (boxed, l.132-141). Distinct deliverables: (D1) the boxed normalized profile `σ_Π`; (D2) its unit normalization on `[0,L]`; (D3) zero-flux/stationary-ODE satisfaction; (D4) the chemical-potential → Onsager-current law that PRODUCES the ODE; (D5) the bias parameter `Π_m = V₁L/Θ_σ`. The card status is `Reduced / ExactClosure` (a derivation-ledger entry, not an unconditional theorem).

## What the script claims to verify

Both scripts hard-define the closed-form profile `σ = Π e^{−Πz/L}/(L(1−e^{−Π}))` (sympy l.13; wl l.33) and the reduced current `J = −M(Θσ' + V₁σ)` (sympy l.14; wl l.34). They then verify three identities: (1) the profile integrates to 1 on `[0,L]` (sympy l.17,26; wl l.38); (2) after substituting `V₁ → ΠΘ/L`, the current `J` vanishes (sympy l.19,28; wl l.40,42); (3) under the same substitution the ODE residual `Θσ' + V₁σ` vanishes (sympy l.23-24,30; wl l.44,46). They print `Π_m = V₁L/Θ`. The SymPy script also computes `mu = Θ log(σ/σ_*) + V₁z` (l.22) but never uses it in any assertion — it is dead code; the `.wl` does not define `mu` at all.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| D1 boxed profile `σ_Π` | profile hard-defined identical to boxed form (sympy l.13 / wl l.33) | match (definition, not derivation) |
| D2 unit normalization | `∫₀^L σ = 1` asserted (sympy l.26 / wl l.38) | match |
| D3 zero-flux / stationary ODE | `J=0` and residual `=0` after `V₁→ΠΘ/L` (sympy l.28,30 / wl l.42,46) | match |
| D4 chemical-potential → current law `J_σ = −M_σ σ ∂_z μ_σ^chem` | `J` defined directly in reduced form; `μ` computed (sympy l.22) but never linked to `J` by any assertion; never present in `.wl` | missing |
| D5 `Π_m = V₁L/Θ_σ` | printed only (sympy l.34 / wl l.50) | match (print, no assertion needed — it is a definition) |

Dominant pattern: the headline boxed deliverables (D1/D2/D3/D5) are faithfully exercised; one stated derivation link (D4, μ→J) is unexercised. No identity mismatch and no value mismatch exist — `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 26 | `integrate(sigma,(z,0,L)) - 1 == 0` | D2 normalization | yes (denominator `L(1−e^{−Π})` must be exact) |
| A2 | sympy | 28 | `simplify(J.subs(V1→ΠΘ/L)) == 0` | D3 zero-flux | yes (ties profile decay rate to `Π`) |
| A3 | sympy | 30 | `(Θσ'+V₁σ).subs(V1→ΠΘ/L) == 0` | D3 stationary ODE | yes (same content as A2; redundant but non-tautological) |
| A4 | math | 38 | `expectZero[Integrate[sigma,{z,0,lM}]-1]` | D2 normalization | yes |
| A5 | math | 42 | `expectZero[jSub]` (`v1→piM thetaSigma/lM`) | D3 zero-flux | yes |
| A6 | math | 46 | `expectZero[residual]` (same subst) | D3 stationary ODE | yes |
| — | sympy | 22 | `mu = Θ log(σ/σ_*) + V₁z` | D4 μ→J link | NO — computed, never asserted (dead code) |

A2/A3 (and A5/A6) are not tautologies-by-construction: `σ` is defined with decay rate `Π/L`, and the checks confirm that this rate is consistent with the independently-imposed relation `V₁ = ΠΘ/L`. A wrong exponent in `σ` would make the residual nonzero. They are, however, mutually redundant (A3 is A2 with the `−M` factor stripped). The normalization (A1/A4) is the most substantive check (the `L(1−e^{−Π})` denominator must be exactly right). No assertion covers D4.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage129_mouth_boundary_layer_mathematica_audit.wl:33-46`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage129_mouth_boundary_layer_sympy_audit.py:13-30`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent re-derivation. Corresponding sections:
- profile definition — sympy l.13 `sigma = Pi*sp.exp(-Pi*z/L)/(L*(1-sp.exp(-Pi)))`; wl l.33 `sigma = piM*Exp[-piM*z/lM]/(lM*(1 - Exp[-piM]));` — identical closed form, just renamed symbols.
- current — sympy l.14 `J = -M*(Theta*sp.diff(sigma, z) + V1*sigma)`; wl l.34 `jSigma = -mobility*(thetaSigma*D[sigma, z] + v1*sigma);` — identical structure.
- zero-flux subst — sympy l.19 `J.subs({V1: Pi*Theta/L})`; wl l.40 `jSigma /. v1 -> piM*thetaSigma/lM` — identical substitution choreography.
- ODE residual subst — sympy l.23-24 / wl l.44 — identical.
Both engines take the SAME pre-written answer `σ_Π` as a given and verify the SAME three identities with the SAME substitution. Neither engine re-derives `σ` from the zero-flux ODE (`Θσ' + V₁σ = 0`) plus normalization — which is exactly the derivation the notes perform (l.98-119). The second engine echoes the first engine's algebra rather than confirming the result by an independent route.

**Why this matters:**
The dual-engine policy exists so a second engine catches errors the first cannot — but a transliteration shares the first engine's assumptions and choreography, so it can only catch typos, not derivation errors. The whole boxed result (D1) rests on a definition both scripts simply copy; nothing independently confirms the exponential is the actual solution of the stated ODE.

**Required change:**
In the `.wl`, derive the profile independently rather than copying it. Replace the hard-coded `sigma` (wl l.33) with a `DSolve`-based recovery: solve `thetaSigma*sigmaFn'[z] + v1*sigmaFn[z] == 0` for `sigmaFn[z]`, impose the normalization `Integrate[sigmaFn[z], {z,0,lM}] == 1` to fix the integration constant, substitute `v1 -> piM*thetaSigma/lM`, and `expectZero` that the resulting closed form minus `piM*Exp[-piM*z/lM]/(lM*(1 - Exp[-piM]))` is identically 0. Keep the existing normalization / zero-flux / residual checks as confirmations. This makes the `.wl` an independent derivation of D1 rather than a copy.

**Verification:**
The verifier confirms a `DSolve[...]` (or equivalent independent solve) appears in the `.wl`, that a new `expectZero` comparing the derived profile to the boxed form is present and passes, and that `math -script` exits 0. The SymPy script is unchanged.

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage129_mouth_boundary_layer_sympy_audit.py:22`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage129_mouth_boundary_layer_mathematica_audit.wl` (μ→J link absent entirely)

**What's wrong:**
The notes' central derivation is that the reduced current `J_σ = −M_σ(Θ_σ σ' + V₁σ)` is the consequence of the Onsager law `J_σ = −M_σ σ ∂_z μ_σ^chem` applied to the electrochemical potential `μ_σ^chem = Θ_σ ln(σ/σ_*) + V₁z` (notes l.64-89). The scripts assert the *reduced* form by fiat: sympy l.14 and wl l.34 define `J = −M(Θσ' + V₁σ)` directly. The SymPy script computes `mu = Theta*sp.log(sigma/sigma_star) + V1*z` on l.22 but never uses it in any assertion — `mu` is dead code; the `.wl` omits `μ` entirely. The identity that should be checked — `−M σ ∂_z μ = −M(Θσ' + V₁σ)` (true because `∂_z μ = Θ σ'/σ + V₁`) — connects the chemical-potential definition to the reduced current and is the actual content of deliverable D4. It is never exercised in either engine.

**Why this matters:**
A reader (and the audit) is told the boundary-layer law derives from a GNLS/localized-Maxwell electrochemical balance, but the scripts only verify properties of the final profile, not the chemical-potential → current step that motivates the whole stage. The unused `mu` variable signals an intended but abandoned check. Without it, the `Θ_σ ln(σ/σ_*)` (entropic) and `V₁z` (potential) structure of `μ_σ^chem` is never tested at all.

**Required change:**
Add one assertion in each engine linking `μ_σ^chem` to the reduced current. In SymPy, after l.22, add:
`assert sp.simplify((-M*sigma*sp.diff(mu, z)) - (-M*(Theta*sp.diff(sigma, z) + V1*sigma))) == 0` (and a corresponding `print`). This uses the already-computed `mu`, so no new symbols. In the `.wl`, define `mu = thetaSigma*Log[sigma/sigmaStar] + v1*z;` and add `expectZero["chemical-potential current identity", (-mobility*sigma*D[mu, z]) - (-mobility*(thetaSigma*D[sigma, z] + v1*sigma))];`.

**Verification:**
The verifier confirms a new μ→J assertion appears in both scripts, that it passes (residual exactly 0), and that both scripts exit 0. The new check is non-tautological: it would fail if `μ_σ^chem` were defined with a wrong coefficient on either term.

## Independent-derivation check (Mathematica)

The `.wl` is a faithful transliteration of the `.py` — see F1. Same symbol set (`z, lM, piM, thetaSigma, v1, mobility, sigmaStar` ↔ `z, L, Pi, Theta, V1, M, sigma_star`), same hard-coded `sigma` (wl l.33 ↔ py l.13), same `jSigma` (wl l.34 ↔ py l.14), same three checks in the same order with the same `v1 -> piM*thetaSigma/lM` substitution (wl l.40,44 ↔ py l.19,24), and same trailing prints. No independent solve (`DSolve`) is performed; the profile is copied, not re-derived.

## Engine cross-check

Both engines produce identical results (outputs fresh — see below):
- SymPy output: `sigma_Pi(z) = Pi*exp(Pi*(L - z)/L)/(L*(exp(Pi) - 1))` (= boxed form ×e^Π/e^Π), `Normalization = 1`, `Zero-flux current J_sigma = 0`, `Stationary zero-flux ODE residual = 0`.
- Mathematica output: `sigma_Pi(z) = piM/(E^((piM*z)/lM)*(1 - E^(-piM))*lM)` (= boxed form), `Normalization = 1`, all three `PASS` with residuals `= 0`.
The two profile printings are algebraically identical (multiply num+denom of one by `e^Π`). The engines agree fully. (This agreement is unsurprising precisely because the `.wl` is a transliteration — F1.)

Freshness: sympy script mtime `2026-05-11 11:56:51` < output `2026-05-11 12:45:34`; wl script mtime `2026-05-27 17:20:55` < output `2026-05-27 17:46:31`. Both outputs are NEWER than their scripts and their content matches the current scripts → fresh; no `stale_output` finding.

## Verdict justification

The headline boxed deliverables hold up: the normalization check (A1/A4) is substantive (the exact `L(1−e^{−Π})` denominator is required), and the zero-flux/ODE checks (A2/A3/A5/A6) are non-tautological because they tie the profile's decay rate to the imposed `V₁ = ΠΘ/L`. I attacked the normalization integral (re-derived `∫=1` by substitution — holds), the symbol domains (all positivity matches the notes' `Π,L,Θ_σ > 0`; no invalid `simplify`), and the profile printings across engines (algebraically identical). Two real defects remain: (F1) the `.wl` is a line-by-line transliteration that never independently derives the profile from the ODE, weakening the second-engine guarantee; and (F2) the chemical-potential → current link (deliverable D4) — the whole physical motivation of the stage — is never asserted, with the SymPy `mu` variable left as dead code and `μ` absent from the `.wl`. Both are script-side and mechanically fixable; neither changes the paper's claim, so `verdict: findings`, `stop_cold: null`, `paper_alignment: aligned`. I confirm I read the paper card, the notes, and the appendix; the script's verified profile matches the paper's boxed claim exactly.

## Value Reconciliation (pass-2 augmentation)

All emitted deliverable values are symbolic (no numeric constants are computed). Enumerating every RESULT/labeled value the scripts emit:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `σ_Π(z) = Π e^{−Πz/L}/(L(1−e^{−Π}))` | py l.13/out l.9; wl l.33/out l.5 | notes l.113-119 (boxed); `.tex` quote l.16 ("exponential source profile σ_Π") | MATCH |
| Normalization `∫₀^L σ dz = 1` | py l.17/out l.10; wl l.37-38/out l.6 | notes l.108-110 | MATCH |
| Zero-flux current `J_σ = 0` | py l.20/out l.11; wl l.41/out l.9 | notes l.92-100 | MATCH |
| Stationary zero-flux ODE residual `= 0` (`Θσ'+V₁σ=0`) | py l.24/out l.12; wl l.45/out l.12 | notes l.98-100 | MATCH |
| `Π_m = V₁L/Θ_σ` | py l.34/out l.15; wl l.50/out l.17 | notes l.132-141 (boxed); l.118 (`Π := V₁L/Θ_σ`) | MATCH |

INTERNAL (scaffolding, not expected in prose, no finding): `J = −M(Θσ'+V₁σ)` reduced current definition; `mu = Θ log(σ/σ_*) + V₁z` (computed, unused — see F2); pass/fail flags; the `v1 → piM*thetaSigma/lM` substitution constant.

Note: the `.tex` card is intentionally terse and carries only the prose-level claim; every deliverable value lives correctly in the `.md` notes, which is the natural carrier — so all are MATCH. No MISMATCH and no MISSING-DELIVERABLE. The F1/F2 findings are about verification adequacy and independence, NOT value discrepancies, so they do not create reconciliation `paper_misalignment`s.

`reconciliation: complete; 5 values checked, 0 misaligned`

## Self-test notes

I checked: (1) Variable independence — `D[sigma, z]`/`sp.diff(sigma, z)` is correct since `σ` genuinely depends on `z`; the proposed F2 `D[mu, z]` likewise depends on `z` (via `log σ` and the `V₁z` term) so the derivative is not identically zero. (2) Trivial-case pre-check — the proposed F2 residual `−Mσ∂_zμ − (−M(Θσ'+V₁σ))` reduces to exactly 0 analytically (`∂_zμ = Θσ'/σ + V₁` ⇒ `−Mσ∂_zμ = −M(Θσ'+V₁σ)`), confirming it is a passing-but-non-tautological check. (3) Path specifications — F1/F2 target the existing `.wl` in `mathematica/` and `.py` in `scripts/`; no new files. (4) Paper round-trip — neither F1 nor F2 introduces a new value or constant, so no new `paper_misalignment` is created; both keep the boxed `σ_Π` and `Π_m = V₁L/Θ_σ` exactly as the notes state.
