---
unit_id: 165
batch: V.1
created_at: 2026-06-08T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-08T16:08:34-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 165

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage165_exact_branch_drifts_mathematica_audit.wl:35-37`

**Issue:** The `.wl` solver block (eqR/eqG + Solve) is a verbatim coefficient-for-coefficient port of the SymPy block (`scripts/...stage165...py:47-49`); the second engine adds no independent route to the result. The achievable independence is limited (single-path linear algebra), so this is a low-severity policy-fidelity flag, not a correctness flag — engines already agree exactly.

**Required change:**
Give the Mathematica script its own derivation of the two linear conditions `eqR`/`eqG` instead of writing the coefficient vectors as literals at wl:35-36. Build them from the Stage-164 imbalance channels expressed as log-derivative sums, then set each channel to zero to OBTAIN the linear equations:

1. Above the current `eqR`/`eqG` definitions (wl:35), insert symbolic log-drift definitions matching notes §3 (`K_s K_q/λ²`) and §4 (`g_q K_s/(g_s λ)`). Using the notes scalings — `λ ∝ a²·ℓ·L_W^{1/2}·v_w0` with `ℓ ∝ 1/c_{s,w}`; `K_q ∝ Z_q·c_s²/L_W²`; `K_s ∝ a²/(ρ_w·ℓ) ∝ a²·c_{s,w}/ρ_w`; `J_s ∝ a²·ℓ ∝ a²/c_{s,w}` — assemble the two channel log-drifts as linear combinations of `{dZ, drho, dcsw, dcs, dv, da, dLW}` (and `dT` for the fixed-g channel via `T_m`'s law in notes §4). Define them so that:
   - `chanR := δln(K_s K_q/λ²)` reduces to `dZ + 2*dcs + 3*dcsw - drho - 2*dv - 2*da - 3*dLW`
   - `chanG := δln(g_q K_s/(g_s λ))` reduces to `dZ + 3*dcsw - drho - dT - dv - 2*da - 2*dLW`
2. Set `eqR = (chanR == 0); eqG = (chanG == 0);` — i.e. derive, do not transcribe.
3. Add a one-line assertion that the assembled coefficient vectors match the SymPy ones (e.g. `expectZero["eqR matches Stage-164 fixed-r channel", chanR - (dZ + 2*dcs + 3*dcsw - drho - 2*dv - 2*da - 3*dLW)]` and the analogous one for `chanG`) so the verifier can confirm the derived form is correct.
4. Leave the rest of the script (Solve, DN substitution, ratio/product/n=5 checks) unchanged.

If reconstructing the channel scalings from the notes is judged too underspecified to transcribe safely (the notes give proportionalities, not normalized constants — but only the LINEAR log-derivative coefficients matter here, and those ARE fully determined), append `## Blocked: F1` with the specific ambiguity rather than guessing a wrong coefficient.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 165` and confirms (a) wl:35-36 no longer hold the literal coefficient sums verbatim, (b) the two new "eqR/eqG matches channel" expectZero checks PASS, (c) all four original expectZero checks still PASS, and (d) the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage165_exact_branch_drifts_mathematica_audit.wl`
- summary: Replaced the literal fixed-r/fixed-g equations with channel log-drifts assembled from the Stage-164 scaling definitions and checked those channels against the expected coefficient forms.
- deviation: none

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.py:97` (append a new numeric block before the final carry-forward prints, or after them, but before EOF)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage165_exact_branch_drifts_mathematica_audit.wl:83` (append a numeric block before the final `Exit[0]` at wl:88)

**Issue:** Notes §5-6 box closed-form amplitude laws and report five numeric branch constants as stage results, but neither engine evaluates them; a transcription slip in any quoted decimal would go undetected.

**Required change (both engines, mirrored):**
Add a numeric verification block that pins the lower-branch constants and checks each closed-form prefactor against the notes decimal to tolerance `1e-12`.

1. Pin: `rstar` from the canonical Family-1 closed form `Sqrt[4107 - 100*Pi^2]/(10*Pi)` (= `1.77799353547498…`; do NOT paste the bare decimal — derive it, mirroring stages 168/169), and `gstar = 0.758035078944663` (notes §5 lines 198-200; `g_*` has no simpler closed form in this card's inputs, so it legitimately remains a decimal). Using the exact symbolic `rstar` makes the four prefactors exact, so their `N[...]` match the notes decimals to full precision.
2. Compute the four prefactors (use high-precision arithmetic — SymPy `sp.Float(..., 30)` / `sp.nsimplify`-free direct float on exact symbolic forms, or Mathematica `N[..., 30]`):
   - `Tm_pref   = 3*Sqrt[10]*3^(3/4) / (5*Pi*gstar*(1+rstar^2)^(1/4))`  → expect `1.2715890393387603` (notes line 232)
   - `v_pref    = 9*Sqrt[10]*3^(1/4)*rstar / (20*(1+rstar^2)^(3/4))`    → expect `1.1428896163056477` (notes line 241; magnitude — the law carries a leading minus)
   - `ratio_pref = Sqrt[3]*Pi*gstar*rstar / (4*Sqrt[1+rstar^2])`        → expect `0.8987885086678338` (notes line 310; magnitude)
   - `prod_pref  = 81*rstar / (10*Pi*gstar*(1+rstar^2))`                → expect `1.4532859092683434` (notes line 318; magnitude)
3. Assert each `Abs[pref - expected] < 1e-12` (SymPy: `expect_zero`-style with a tolerance variant, or `assert abs(float(pref) - expected) < 1e-12`; Mathematica: `If[Abs[pref-expected] < 10^-12, pass[...], fail[..., pref-expected]]`). Print each computed prefactor as a labeled result.

Use the constants EXACTLY as quoted in the notes — do not introduce any new constant or re-derive `rstar`/`gstar` from upstream (that is out of scope for this stage; they are carried-forward inputs).

**Self-test confirmation (already done by auditor):** the prefactor forms are copied verbatim from notes §5-6 boxes; they depend on `rstar`/`gstar` (live, not identically constant); substituting the quoted `rstar=1.77799…`, `gstar=0.758035…` into `Sqrt[3]*Pi*gstar*rstar/(4*Sqrt[1+rstar^2])` yields ≈0.8988, matching the notes — so the assertions are non-trivially true, not silent-pass.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 165` and `redteam exec-mathematica 165` and confirms four new labeled prefactor prints appear in each output, all four numeric comparisons PASS in both engines, and both scripts exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage165_exact_branch_drifts_mathematica_audit.wl`
- summary: Added mirrored lower-branch numeric prefactor checks for Tm_pref, v_pref, ratio_pref, and prod_pref against the notes decimals at 1e-12 tolerance.
- deviation: none
