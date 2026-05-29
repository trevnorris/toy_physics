---
unit_id: 160
batch: IV.6
created_at: 2026-05-27T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 160

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line range named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

Do NOT modify the SymPy script. Only the Mathematica `.wl` is in scope for this directive.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage160_bare_mixed_port_slippage_mathematica_audit.wl:28-46`

**Issue:**
The Mathematica audit is a line-by-line port of the SymPy audit: identical symbol names (`eps, rc, drc, dk0, dg0`), identical intermediate names (`k0Star, g0Star, kW, gW, dkW, dgW`), identical recipe (`Series[..., {eps, 0, 1}]` followed by `Coefficient[..., eps, 1]`), identical identity assembly, identical gate substitution, and identical hardcoded numerics in the print-only tangential block (`0.832409471081635`, `1.16275838754222`). The second-engine policy requires both engines to derive the claim independently from the physical premises. The current `.wl` only confirms that Mathematica can execute SymPy's algorithm.

**Required change:**

Rewrite lines 28-46 of `mathematica/moving_throat_pde_stage160_bare_mixed_port_slippage_mathematica_audit.wl` so the derivation reaches the same residual via a structurally different path. Specifics:

1. Rename internal symbols away from SymPy's spellings. Suggested renames (use these exactly so the verifier can confirm the diff):
   - `eps` → remove entirely; no perturbation parameter
   - `rc` → `rStar`
   - `drc` → `deltaR`
   - `dk0` → `dKappa0`
   - `dg0` → `dGamma0`
   - `k0Star` → `kappa0Canon`
   - `g0Star` → `gamma0Canon`
   - `kW` → `kappaWperturbed`
   - `gW` → `gammaWperturbed`
   - `dkW` → `deltaKappaW`
   - `dgW` → `deltaGammaW`

2. Replace the `Series` + `Coefficient` recipe with a direct first-derivative extraction at the canonical point. Concretely:

   Before (lines 31-36):
   ```
   k0Star = (1 + rc)/3;
   g0Star = (1 + rc)/9;
   kW = Normal[Series[(k0Star + eps*dk0)/(1 + rc + eps*drc), {eps, 0, 1}]];
   gW = Normal[Series[(g0Star + eps*dg0)/(1 + rc + eps*drc), {eps, 0, 1}]];
   dkW = Expand[Coefficient[kW, eps, 1]];
   dgW = Expand[Coefficient[gW, eps, 1]];
   ```

   After:
   ```
   kappa0Canon = (1 + rStar)/3;
   gamma0Canon = (1 + rStar)/9;
   (* Total differential of kappa_W = kappa_0 / (1 + r_c) at the canonical point. *)
   (* d(kappa_W) = (1/(1+r_c)) d(kappa_0) - (kappa_0/(1+r_c)^2) d(r_c) *)
   deltaKappaW = Together[
     dKappa0/(1 + rStar) - (kappa0Canon/(1 + rStar)^2) * deltaR
   ];
   deltaGammaW = Together[
     dGamma0/(1 + rStar) - (gamma0Canon/(1 + rStar)^2) * deltaR
   ];
   ```

3. The `Clear[...]` and `$Assumptions` lines must be updated to reference the new symbol names. The two `expectZero` calls and the `dgWGate` line must be retargeted to the new symbol names:

   Before (lines 28-29):
   ```
   Clear[eps, rc, drc, dk0, dg0];
   $Assumptions = Element[{eps, rc, drc, dk0, dg0}, Reals];
   ```
   After:
   ```
   Clear[rStar, deltaR, dKappa0, dGamma0];
   $Assumptions = Element[{rStar, deltaR, dKappa0, dGamma0}, Reals];
   ```

   Before (line 41):
   ```
   identity = dgW - 1/3*dkW - (dg0 - 1/3*dk0)/(1 + rc);
   ```
   After:
   ```
   identity = deltaGammaW - (1/3)*deltaKappaW - (dGamma0 - (1/3)*dKappa0)/(1 + rStar);
   ```

   Before (lines 44-46):
   ```
   dgWGate = FullSimplify[(dg0 - 1/3*dk0)/(1 + rc)];
   Print["dγ_W under dκ_W = 0 = ", fmt[dgWGate]];
   expectZero["pure-scale harmlessness", dgWGate /. dg0 -> dk0/3];
   ```
   After:
   ```
   deltaGammaWGate = FullSimplify[(dGamma0 - (1/3)*dKappa0)/(1 + rStar)];
   Print["dγ_W under dκ_W = 0 = ", fmt[deltaGammaWGate]];
   expectZero["pure-scale harmlessness", deltaGammaWGate /. dGamma0 -> dKappa0/3];
   ```

   The two `Print["dκ_W = ", ...]` and `Print["dγ_W = ", ...]` lines (lines 38-39) must also be retargeted to `deltaKappaW` / `deltaGammaW`.

4. Do not modify the tangential-block (lines 48-69) or the final `Print` of carry-forward formulas. Those are downstream print-only; the transliteration concern was only about the load-bearing derivation block.

5. Leave the SymPy file untouched.

**Self-test reasoning (already performed by auditor):**
- The derivative form `d(kappa_W) = d(kappa_0)/(1+r_c) - kappa_0/(1+r_c)^2 d(r_c)` is the chain rule applied directly at the canonical point — algebraically equivalent to extracting the eps^1 coefficient of `Series[(k0Star + eps dk0)/(1+rc+eps drc), {eps,0,1}]`, but expressed without an explicit perturbation parameter or `Coefficient[]` extraction.
- Substituting `kappa0Canon = (1+rStar)/3`, the residual `deltaGammaW - (1/3) deltaKappaW - (dGamma0 - (1/3) dKappa0)/(1+rStar)` reduces to:
  `(dGamma0/(1+rStar) - (1+rStar)/9 / (1+rStar)^2 deltaR) - (1/3)(dKappa0/(1+rStar) - (1+rStar)/3/(1+rStar)^2 deltaR) - (dGamma0 - dKappa0/3)/(1+rStar)`
  = `dGamma0/(1+rStar) - deltaR/(9(1+rStar)) - dKappa0/(3(1+rStar)) + deltaR/(9(1+rStar)) - dGamma0/(1+rStar) + dKappa0/(3(1+rStar))`
  = `0`. Confirmed by hand.
- Pure-scale residual after `dGamma0 -> dKappa0/3`: `(dKappa0/3 - dKappa0/3)/(1+rStar) = 0`. Confirmed.
- Symbol-name diff: none of the new names (`rStar`, `deltaR`, `dKappa0`, `dGamma0`, `kappa0Canon`, `gamma0Canon`, `kappaWperturbed`, `gammaWperturbed`, `deltaKappaW`, `deltaGammaW`) appear in the SymPy script. Verifier can grep for `Series` and `Coefficient` in the new `.wl` and find neither.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 160` and confirm:
- the script exits 0
- both `expectZero` calls still pass
- the literal substrings `Series[`, `Coefficient[`, `k0Star`, `g0Star`, `dk0`, `dg0`, and `, eps,` no longer appear in the `.wl` file (i.e., the derivation block was actually restructured, not just superficially renamed)
- the residual still simplifies to `0`
