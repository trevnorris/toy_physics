---
unit_id: 193
batch: V.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface.md]
  paper_appendix: present
---

# Audit unit 193 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_193.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 117, 265, 1216–1261, 1456)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

Stage 193 freezes the exact conservative front-end target surface for the grouped-real `P2` (20/21/22) bundle and proves a scalar/geometry firewall. `\stagefield{Output}`: "Freezes the isotropic one-pole conservative target surface and proves the scalar/geometry lane enters the grouped carrier only at quadratic anisotropy order." The notes and appendix (eqs. `app-part05-conservative-target-surface`, `-one-pole-surface-D`, `-conservative-one-pole-module`, `-geometry-firewall-schur`, `-geometry-firewall-result`) enumerate four deliverables: (1) the isotropy surface `a_2=b_2=a_4=b_4=0`, equivalent to common-lane collapse; (2) the one-pole identity `Delta_pole = nu4bar - 4 nu2bar^2 = -(3D2^2 + D0 D4)/D0^2`, equivalent to `D0 D4 + 3 D2^2 = 0`; (3) the forced one-parameter carrier `Y_Q^cons(omega) = 3/4 + 1/4 / (1 - omega^2/Omega_Q^2)` with `Omega_Q^2 = -D0/(4 D2)`, matching `1 + nu2bar omega^2 + nu4bar omega^4 + O(omega^6)`; (4) the Schur-complement firewall `D_{2,eff}(omega,chi) = D_2(omega) I_3 - chi^2 C(omega) D_0(omega)^{-1} C(omega)^T`, hence `partial_chi D_{2,eff}|_{chi=0} = 0` (no linear `l=0 <-> l=2` contamination of the grouped `l=2` carrier). The notes use the original-stage label "Stage 244" internally; the appendix maps that derivation to paper-stage 193, and the content is identical.

## What the script claims to verify

The SymPy script (banner "STAGE 176", ledger "Stage 193") exercises: in §I the grouped trace/anomaly forward+inverse map and vanishing of the conservative anomalies on the common-lane isotropic branch (claim 1); in §II the closed form `Delta_pole = -(3D2^2 + D0 D4)/D0^2` and the equivalence to `D0 D4 + 3 D2^2 = 0` via substitution `D4 = -3D2^2/D0` (claim 2); in §III the series of `3/4 + 1/4/(1 - omega^2/Omega_Q^2)` with `Omega_Q^2 = -D0/(4D2)` against `1 + nu2 omega^2 + nu4|_onepole omega^4` (claim 3); in §IV the firewall, building `Deff = D0 I_3 - chi^2 C^T C / K0` and asserting its chi^2 structure and that `d/dchi Deff|_{chi=0}=0` (claim 4).

## Paper <-> script cross-check

| Paper deliverable | Script check | Verdict |
|---|---|---|
| (1) `a2=b2=a4=b4=0` / common-lane collapse | §I lines 71–76 (`a2_iso`,`b2_iso`,`a4_iso`,`b4_iso`) + inverse-map round-trip 66–68 | match |
| (2) `Delta_pole=-(3D2^2+D0D4)/D0^2` and `D0D4+3D2^2=0` | §II lines 94–103 | match |
| (3) `Y=3/4+1/4/(1-omega^2/Omega_Q^2)`, `Omega_Q^2=-D0/(4D2)` | §III lines 110–120 | match |
| (4) Schur firewall `D_{2,eff}=...-chi^2 C D_0^{-1} C^T`, `partial_chi|_0=0` | §IV lines 133–145 | partial (form correct, but firewall is assumed not derived — see F1) |

`paper_alignment: aligned` — every paper deliverable maps to a script check; the only defect is the *strength* of the §IV checks (script-side), not a target/value mismatch with the paper. The math the paper states is correct.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 66–68 | `simplify(x_back - x) == 0` | claim 1 (exact inverse map) | yes |
| A2 | sympy | 73–76 | `a/b on (nu,nu,nu) == 0` | claim 1 (isotropic collapse) | yes |
| A3 | sympy | 94–97 | `Delta_pole + (3D2^2+D0D4)/D0^2 == 0` | claim 2 (closed form) | yes |
| A4 | sympy | 100–103 | `Delta_pole.subs(D4,-3D2^2/D0) == 0` | claim 2 (surface equivalence) | yes |
| A5 | sympy | 120 | `series(Y_pole) - Y_expected == 0` | claim 3 (one-parameter carrier) | yes |
| A6 | sympy | 140 | `Delta_geom/chi^2 + C^TC/K0 == 0` | claim 4 (firewall coeff) | no — tautological by construction |
| A7 | sympy | 141–144 | `d/dchi Deff|_{chi=0} == 0` | claim 4 (firewall, linear order vanishes) | no — tautological by construction |
| A8 | sympy | 145 | `Deff.subs(chi,0) - D0 I == 0` | claim 4 (chi=0 reduces to grouped block) | no — trivially true by construction |

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_sympy_audit.py:127-145`

**What's wrong:**
The paper's firewall claim (eqs. `app-part05-geometry-firewall-schur`/`-result`, notes §5.2) is a *Schur-complement* statement: the reduced block operator has a **linear** off-diagonal coupling,
```
D(omega,chi) = [[ D_0 , chi C^T ], [ chi C , D_2 I_3 ]],
```
and eliminating the scalar/geometry block by the exact Schur complement,
```
D_{2,eff} = D_2 I_3 - (chi C)(D_0)^{-1}(chi C)^T = D_2 I_3 - chi^2 C D_0^{-1} C^T,
```
**produces** the chi^2 factor. The non-trivial content is that a coupling that is *linear* in chi gets quadratured into chi^2 by the Schur complement, which is precisely why there is no O(chi) contamination.

The script never forms the block operator and never takes its Schur complement. It instead writes the *already-reduced, already-quadratic* answer by hand:
```python
Deff = sp.simplify(D0 * I3 - chi**2 * Cvec.T * Cvec / K0)   # line 133
Delta_geom = sp.simplify(Deff - D0 * I3)                     # line 134
```
and then "verifies" properties that are guaranteed by that very construction:
- line 140 `Delta_geom/chi**2 + (Cvec.T*Cvec)/K0`: since `Delta_geom = -chi^2 C^TC/K0` *by definition*, dividing by `chi^2` and adding `C^TC/K0` is identically `0` — it cannot fail for any choice of `C`, `K0`, or physics.
- lines 141–144 `Deff.diff(chi).subs(chi,0)`: because `Deff` was *written* with `chi**2`, its first chi-derivative is `-2 chi C^TC/K0`, which vanishes at `chi=0` for any chi^2 expression. This is true of any hand-written O(chi^2) form; it does not test that the Schur complement of a *linearly* coupled block is O(chi^2).
- line 145 `Deff.subs(chi,0) - D0 I3`: trivially `0` by construction.

So §IV asserts "a hand-written chi^2 expression is O(chi^2)" — true but vacuous. The load-bearing physics (linear coupling -> quadratic correction under Schur reduction) is assumed, not exercised. Secondarily, the leading block is written `D0 I3` (reusing the §II symbol `D0`, the l=2 conservative leading moment) where the paper's leading block is `D_2(omega) I_3` and the eliminated scalar block is `D_0(omega)` (here proxied by the scalar `K0`); the firewall conclusion is structurally independent of which scalar sits on the diagonal, so this is a naming/relabel issue, not a separate finding — but a faithful Schur derivation should make the block roles explicit.

**Why this matters:**
This is the only one of the four paper deliverables whose script "proof" cannot fail. A reader trusting the PASS transcript would believe the firewall theorem was verified, when in fact the script merely restated its conclusion. If the true Schur complement of `[[D_0, chi C^T],[chi C, D_2 I_3]]` did *not* reduce to a chi^2 correction (e.g., if the coupling were chi^0 or chi^1), this script would still PASS. The check provides zero adversarial protection for claim 4.

**Required change:**
Build the actual 4x4 (or 1+3 block) operator with *linear* chi coupling and compute the Schur complement, then assert the result equals the quadratic form. Concretely, in §IV:
1. Introduce the scalar block as a symbol `D0scalar` (= paper `D_0(omega)`, nonzero) and the l=2 leading block `D2blk` (= paper `D_2(omega)`); keep `Cvec` (1x3) and `chi`.
2. Form the block matrix
   ```python
   Dblock = sp.Matrix(sp.BlockMatrix([
       [sp.Matrix([[D0scalar]]), chi * Cvec],
       [chi * Cvec.T,            D2blk * I3]]))
   ```
   (so the off-diagonal is **linear** in chi, matching the paper's `chi C`).
3. Compute the Schur complement eliminating the scalar block:
   ```python
   Deff_schur = sp.simplify(D2blk * I3 - (chi * Cvec.T) * (sp.Matrix([[D0scalar]]).inv()) * (chi * Cvec))
   ```
4. Assert it equals the quadratic target and that its chi-linear part vanishes:
   ```python
   expect_zero("Schur complement matches D2 I - chi^2 C^T C / D0scalar",
               sp.Matrix(Deff_schur - (D2blk * I3 - chi**2 * Cvec.T * Cvec / D0scalar)))
   expect_zero("linear chi part of Schur complement vanishes",
               sp.Matrix(Deff_schur.diff(chi).subs(chi, 0)))
   ```
   The first assertion is now non-trivial: it can only pass because the off-diagonal was linear in chi and the Schur algebra genuinely produces chi^2. Keep lines 145's chi=0 reduction (it is fine as a sanity print). Use `D0scalar` (the eliminated scalar block) in the denominator, not the §II `D0`; do not reuse the conservative `D0` symbol for two distinct physical blocks.

**Verification:**
After Codex applies, the new §IV must contain (a) an explicit block matrix whose off-diagonal entries are `chi*c2j` (linear, not `chi**2`), and (b) an `expect_zero` that subtracts the *Schur-complement output* from the quadratic target. The refreshed output `.txt` must show both new `expect_zero` lines reporting `0`, and the script must still exit 0. The verifier confirms the off-diagonal block is linear in chi (grep for `chi * Cvec` / `chi*Cvec` without a `**2`) and that the Schur complement, not a hand-written quadratic, is what feeds the assertion.

## Independent-derivation check (Mathematica)

No `.wl` exists for this unit; `mathematica_transliteration` does not apply.

## Engine cross-check

Only one engine present; no cross-check. On the missing-Mathematica question: this stage is `is_status_only_candidate: False` and `is_checkpoint: False`, so it is non-status-only but SymPy-only. I do **not** flag `missing_mathematica`. §I, §II, §III contain substantive, non-tautological symbolic identities that SymPy settles completely and exactly (exact rational/series algebra over symbolic `D0,D2,D4,nu2,nu4,omega` — no numerics, no domain-fragile simplification), matching the established SymPy-only precedent (stages 121/122/123 and 56 pipeline-wide). A second engine would add nothing to the confidence in claims 1–3. The §IV weakness is a *content* defect (assumes the answer) that a Mathematica transliteration of the same hand-written `Deff` would reproduce verbatim, so a second engine is not the remedy; the remedy is the actual Schur derivation specified in F1. Hence no `missing_mathematica` finding.

## Verdict justification

`verdict: findings`. Attacks that failed (claims 1–3 hold up): I tried to break the §II surface-equivalence as a back-substitution tautology — it is not, because `D4_onepole = -3D2^2/D0` is the independent solution of `D0D4+3D2^2=0` for D4 and substituting it into the *independently defined* `Delta_pole = nu4 - 4 nu2^2` genuinely collapses to 0; that is a real implication, not `X-X`. I checked §III for a hidden hardcoded series — the `Y_expected` ω^4 coefficient uses `nu4_common.subs(D4, D4_onepole) = 4 D2^2/D0^2`, which I hand-verified matches the Taylor coefficient of `(D0+3D2 omega^2)/(D0+4D2 omega^2)`, so the series check is substantive and ties Omega_Q^2 to the carrier exactly. §I's anomaly-vanishing and inverse round-trip are genuine algebraic identities. The one defect is §IV: the firewall (claim 4) is asserted via a hand-written chi^2 form whose chi^2 structure and vanishing chi-linear part are guaranteed by construction (A6–A8 cannot fail), so the script never exercises the paper's actual Schur-complement-of-linear-coupling content. That is `insufficient_verification`, medium severity — the result is correct and fixable by deriving the Schur complement in-script, not UNFIXABLE, and not CRITICAL_DOWNSTREAM (the fix confirms the same result the paper already states; no downstream constant changes). I read the stage card, the notes, and the appendix rows; the script targets the correct paper claims (alignment is exact), and only the §IV verification strength is deficient. Output `.txt` (mtime 12:48:43) is newer than the script (11:56:50): no `stale_output`.

## Self-test notes

Trap 1 (variable independence in `Deff.diff(chi)`): chi IS wired into `Deff` (via `chi**2`), so the derivative is not identically-zero-from-an-unwired-symbol; it is zero only after `subs(chi,0)`. The defect is not a dead-symbol derivative but that the chi^2 form was hand-written rather than derived — captured as F1. Trap 3 (trivial-case pre-check): for the proposed F1 fix, with a concrete `D0scalar=1, c20=1,c21=c22=0`, the genuine Schur complement gives `D2blk I - chi^2 diag(1,0,0)`, whose chi-linear part is 0 and which equals the quadratic target — the new assertions are non-trivially satisfiable and would FAIL if the off-diagonal were written linear-in-chi-but-the-target-were-linear, confirming the fixed check has teeth. I also confirmed §II/§III are not back-substitution tautologies and that the appendix/notes claims (eqs. 1224–1260) match the script targets exactly.
