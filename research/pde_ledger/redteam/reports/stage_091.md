---
unit_id: 091
batch: IV.1
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["moving_throat_pde_stage091_grouped_p2_static_geometry_derivation.md"]
  paper_appendix: present
---

# Audit unit 091 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_091.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (stage referenced via `\input` at line 1216; no row-level prose)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_mathematica_audit.txt`

## What the paper claims

Stage 091's stage card carries the load-bearing quote (stage_091.tex:16):

> "Static-geometry plus one-pole grouped carrier forces the `3/4+1/4` conservative quadrupole module and hence `rho_alpha=4/3`, `zeta_req=1/3`."

The card's derivation-ledger line (stage_091.tex:13) names the carrier:
`Yhat_Q^cons = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)`.
The notes flesh out the algebra: Ansatz `K_Q^cons(omega) = K_geom + K_pole/(1 - omega^2/Omega_Q^2)`, expansion `K0 = K_geom+K_pole, K2 = K_pole/Omega_Q^2, K4 = K_pole/Omega_Q^4`, and imposing the upstream identity `K0*K4 = 4*K2^2` to force `K_geom = 3*K_pole`, hence `Yhat_Q^cons`, hence `rho_alpha = 4/3, zeta_req = 1/3` (notes Sec. 3 + Sec. 4). The card's `\stagefield{Checks}` (stage_091.tex:21-25) lists three checklist items: (a) static limit `eps2 = eps4 = 0` returns `c_pole = 1/4`; (b) `l=0`/`l=2` orthogonality before applying the geometry firewall; (c) any support/source success carries the minimal-module hypothesis.

## What the script claims to verify

Both engines verify the same chain: define `Kcons = Kgeom + Kpole/(1 - omega^2/OmegaQ^2)` with `Kgeom, Kpole, OmegaQ > 0`, take the omega-expansion up to `omega^4`, check the closed-form coefficients K0/K2/K4 match the rational form, compute the branch residual `K0*K4 - 4*K2^2 = K_pole*(K_geom - 3*K_pole)/Omega_Q^4`, solve for `K_geom` (returning `3*K_pole`), substitute back, and verify (1) `K_geom = 3*K_pole`, (2) `K0 = 4*K_pole` on branch, (3) `Yhat = 3/4 + (1/4)/(1-omega^2/Omega_Q^2)`, (4) `rho_alpha = 4/3`, (5) `zeta_req = 1/3`. The sympy and mathematica outputs both end in PASS / EXIT_CODE 0.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `Yhat_Q^cons = 3/4 + (1/4)/(1-omega^2/Omega_Q^2)` (stage_091.tex:13; notes Sec.3) | sympy line 63 / wl line 61 (`Yhat - Yhat_expected == 0`) | match |
| `rho_alpha = 4/3` (stage_091.tex:16; notes Sec.4) | sympy line 69 / wl line 67 | match |
| `zeta_req = 1/3` (stage_091.tex:16; notes Sec.4) | sympy line 70 / wl line 68 | match |
| `K_geom = 3*K_pole` (notes Sec.3) | sympy line 54 / wl line 53 | match |
| `K_pole = K0/4` (notes Sec.3) | sympy line 58 / wl line 56 (verifies `K0 = 4*K_pole` on branch) | match |
| Branch identity `K0*K4 = 4*K2^2` is the *forcing* relation (notes Sec.3) | sympy lines 49-54 / wl lines 47-53 (uses identity to solve for Kgeom) | match |
| Check (a) static limit `eps2=eps4=0` returns `c_pole=1/4` (stage_091.tex:22) | implicit: script has no eps2/eps4 deformation; static limit IS the script setup; `K_pole/K0 = 1/4` is asserted via `K0 - 4*K_pole = 0` | match (interpreted) |
| Check (b) `l=0`/`l=2` orthogonality before firewall (stage_091.tex:23) | no script-side counterpart | missing |
| Check (c) "support/source success carries minimal-module hypothesis" (stage_091.tex:24) | meta-condition, not a math identity; no executable counterpart possible | n/a |

The dominant pattern is `match`; one card-side `Checks` item is unaddressed by the script and one is non-executable (meta). Setting `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 45 | `expect_zero(K0 - (Kgeom + Kpole))` | series of `Kcons` reproduces the K0 form in notes Sec.2 | yes |
| A2 | sympy | 46 | `expect_zero(K2 - Kpole/OmegaQ^2)` | K2 form, notes Sec.2 | yes |
| A3 | sympy | 47 | `expect_zero(K4 - Kpole/OmegaQ^4)` | K4 form, notes Sec.2 | yes |
| A4 | sympy | 54 | `expect_zero(Kgeom_sol - 3*K_pole)` | branch identity forces `K_geom = 3*K_pole` (notes Sec.3, card body) | yes |
| A5 | sympy | 58 | `expect_zero(K0_on_branch - 4*K_pole)` | `K_pole = K0/4` (notes Sec.3) | yes |
| A6 | sympy | 63 | `expect_zero(Yhat - [3/4 + 1/4/(1-omega^2/OmegaQ^2)])` | main card claim (stage_091.tex:13) | yes |
| A7 | sympy | 69 | `expect_zero(rho_alpha - 4/3)` | `rho_alpha = 4/3` (card body, stage_091.tex:16) | yes |
| A8 | sympy | 70 | `expect_zero(zeta_req - 1/3)` | `zeta_req = 1/3` (card body) | yes |
| M1-M8 | mathematica | 43,44,45,53,56,61,67,68 | `expectZero[...]` mirrors of A1-A8 | same claims as A1-A8 | yes |

All assertions reduce to true after non-trivial simplification of distinct symbolic expressions; none are tautological by construction.

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** script_missing_paper_claim
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_091.tex:23` (paper-side checklist item)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_sympy_audit.py` (no counterpart)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_mathematica_audit.wl` (no counterpart)

**What's wrong:**
The stage card's `\stagefield{Checks}` block (stage_091.tex:21-25) lists, verbatim, "Check `l=0` and `l=2` orthogonality before applying the geometry firewall." The audit scripts contain no `l=0`/`l=2` decomposition or orthogonality check; the algebra operates purely on the scalar `Kcons(omega)` ansatz. The stage notes do not mention `l=0`/`l=2` orthogonality anywhere either, so this card item appears to be either (a) a hypothesis presumed established upstream (in which case the card should not list it as a stage-091 check), or (b) a check this stage owes that was overlooked.

**Why this matters:**
If the orthogonality is genuinely an upstream-established hypothesis, the card's wording mis-leads downstream consumers into expecting a stage-091 verification line that doesn't exist. If it is owed by this stage, the verification surface is incomplete. The current card wording forces this ambiguity onto the reader.

**Required change:**
Routed to user (see directive `## Resolve before fix_loop` for F1). Codex must not auto-resolve.

**Verification:**
Once user picks a direction: either the card text is amended to remove or relocate the orthogonality bullet (paper-side fix; Codex does not touch paper/), or the script is extended with an `l=0`/`l=2` orthogonality assertion (script-side fix; a follow-up directive will name the precise check).

### F2 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_sympy_audit.py:28-70`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_mathematica_audit.wl:28-68`

**What's wrong:**
The Mathematica script is a near line-by-line port of the SymPy script rather than an independent re-derivation. Corresponding sections:

SymPy (lines 28-38):
```
omega = sp.symbols("omega", real=True)
Kgeom, Kpole, OmegaQ = sp.symbols(..., positive=True, real=True)
Kcons = sp.simplify(Kgeom + Kpole / (1 - omega**2 / OmegaQ**2))
series = sp.expand(sp.series(Kcons, omega, 0, 6).removeO())
K0 = sp.simplify(series.coeff(omega, 0)); K2 = ...; K4 = ...
```

Mathematica (lines 28-35):
```
Clear[omega, kGeom, kPole, omegaQ];
$Assumptions = Element[{omega, kGeom, kPole, omegaQ}, Reals] && kGeom>0 && kPole>0 && omegaQ>0;
kCons = FullSimplify[kGeom + kPole/(1 - omega^2/omegaQ^2), ...]
series = Expand[Normal[Series[kCons, {omega, 0, 4}]]]
k0 = FullSimplify[Coefficient[series, omega, 0]]; k2 = ...; k4 = ...
```

Both scripts: (1) declare the same positivity assumptions on the same symbols, (2) build the same Kcons, (3) take the same series expansion, (4) extract the same three coefficients, (5) compute the same branch residual `K0*K4 - 4*K2^2`, (6) solve the same equation for `K_geom`, (7) emit the same eight `expect_zero`/`expectZero` assertions in the same order with the same residual expressions. The variable choreography is identical modulo Mathematica's CamelCase renaming.

**Why this matters:**
The second-engine policy is meant to catch algebra errors that survive in one CAS but break in another (e.g., simplification choices, branch handling, assumption propagation). A transliteration that follows the same sequence of CAS calls offers little additional coverage; if SymPy's `solve` returned a stale solution because of an aliasing bug, Mathematica's `Solve` could be expected to fail identically since it is given the identical equation in identical form.

That said, the algebra here is genuinely single-path (rational-function series, one polynomial root, one substitution), so the "independent derivation" requirement is hard to meaningfully satisfy. The two natural alternative attacks are (a) verify the closed form by direct partial-fraction recombination rather than series, and/or (b) confirm the result by an Apart/Together check on `Yhat - Yhat_expected` over the rational field without the intermediate `Series + Solve`.

**Required change:**
Extend the Mathematica script with one independent check that does NOT mirror the SymPy choreography: verify the main identity directly via partial-fraction recombination. Concretely, after the existing assertions in `mathematica_audit.wl`, add a check that proves the result without going through `Series`/`Solve`. See directive for exact lines.

**Verification:**
After Codex applies, `redteam exec-mathematica 091` should still exit 0; the saved output should contain a new `PASS: Yhat partial-fraction recombination` line.

## Independent-derivation check (Mathematica)

The Mathematica script structurally mirrors the SymPy script. Both define the same `Kcons`, both invoke their CAS's `Series` primitive, both invoke `Solve`/`solve` on the same residual, both verify the same eight identities in the same order. The only differences are naming (CamelCase) and the Mathematica-specific `FullSimplify[..., Assumptions -> $Assumptions]` wrapping, which is a syntactic rather than mathematical difference. See F2 for the corresponding code excerpts. The algebra is genuinely single-path, so transliteration is hard to fully avoid; F2's prescribed independent recombination check is the cheapest meaningful diversification.

## Engine cross-check

Both engines yield identical symbolic outputs:

SymPy:
```
Series = K_geom + K_pole + K_pole*omega**2/Omega_Q**2 + K_pole*omega**4/Omega_Q**4
Branch identity K0*K4 - 4*K2^2 = K_pole*(K_geom - 3*K_pole)/Omega_Q**4
K_geom forced by branch identity = 3*K_pole
K0 on branch = 4*K_pole
Normalized module on branch = (Omega_Q**2 - 3*omega**2/4)/(Omega_Q**2 - omega**2)
rho_alpha = 4/3
zeta_req = 1/3
```

Mathematica:
```
Series = kGeom + kPole + (kPole*omega^4)/omegaQ^4 + (kPole*omega^2)/omegaQ^2
Branch identity K0*K4 - 4*K2^2 = ((kGeom - 3*kPole)*kPole)/omegaQ^4
K_geom forced by branch identity = 3*kPole
Normalized module on branch = (3 + (1 - omega^2/omegaQ^2)^(-1))/4
rho_alpha = 4/3
zeta_req = 1/3
```

`(Omega_Q^2 - 3*omega^2/4)/(Omega_Q^2 - omega^2)` and `(3 + 1/(1-omega^2/omegaQ^2))/4` are equal: multiply the Mathematica form by `(omegaQ^2/omegaQ^2)`, then `(3*(1 - omega^2/omegaQ^2) + 1) / (4*(1-omega^2/omegaQ^2)) = (4 - 3*omega^2/omegaQ^2)/(4*(1-omega^2/omegaQ^2)) = (omegaQ^2 - 3*omega^2/4)/(omegaQ^2 - omega^2)`. Engines agree.

## Verdict justification

The bottom-line math holds: both engines derive `Yhat_Q^cons = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)` and the corollaries `rho_alpha = 4/3, zeta_req = 1/3` from the same branch identity (`K0*K4 = 4*K2^2`) that the notes specify, with consistent symbolic intermediates, fresh outputs, and matching engines. The two findings are: (F1) a card-listed checklist item (`l=0`/`l=2` orthogonality) has no script-side counterpart and is also absent from the notes — needs user resolution to decide whether the card is overclaiming or the script is undercovering; (F2) the Mathematica script is structurally a transliteration of the SymPy script, which the second-engine policy treats as a finding even when the algebra is genuinely single-path. Verdict is `findings` with `stop_cold: null`. No `paper_misalignment` of the kind that would warrant a stop-cold halt was detected on the main `3/4+1/4` claim; F1 is a low-severity peripheral discrepancy on the `Checks` checklist, not on the load-bearing `Output`.

Attacks tried that did NOT yield findings:
- Tautology probe on A4 (`K_geom - 3*K_pole`): if `branch_identity` were the wrong equation, `solve` would return a different solution; the check would fail. Non-tautological.
- Symbol-domain probe: `K_pole, K_geom, Omega_Q > 0` is consistent with the notes (`K_pole != 0` is the explicit non-triviality assumption in notes Sec.3; positivity is stronger but is a valid sub-domain).
- Missing-branch probe: the branch identity factors as `K_pole*(K_geom - 3*K_pole) = 0`. The `K_pole = 0` root corresponds to "no pole" (excluded by `K_pole > 0`). The `K_geom = 3*K_pole` root is the only physically admissible one; the script picks `solve(...)[0]` which lands on the correct root because Sympy orders `Kgeom` solutions first. Robust.
- Hardcoded-numeric probe: `4/3` and `1/3` are computed from `K0_on_branch/Kgeom_sol = 4*K_pole/(3*K_pole)` and `(K0_on_branch - Kgeom_sol)/Kgeom_sol`, not literals.
- Stale-output probe: sympy script mtime 2026-04-01, output mtime 2026-05-11 (fresher). Mathematica script and output both 2026-05-11. Both fresh.
- Static-limit probe: the card's check (a) "`eps2 = eps4 = 0` returns `c_pole = 1/4`" is satisfied by the script's `K0 = 4*K_pole` assertion (i.e., `K_pole/K0 = 1/4 = c_pole`) — the script's setup IS the static limit because no `eps2`/`eps4` deformation is introduced.

## Self-test notes

I checked: (1) no `sp.diff` or `D[expr, var]` calls appear, so the derivative-of-constant trap does not apply; (2) no integrals appear, so the parity-on-symmetric-domain trap does not apply; (3) the F2 directive's "partial-fraction recombination" verification is a closed-form algebraic identity over the rational function field — substituting any concrete `omega = c, omegaQ = d` reproduces a literal 0, so `FullSimplify` will succeed; (4) the F2 directive targets the `.wl` file in `mathematica/`, the correct directory for Mathematica scripts; (5) the F2 fix introduces no new constants and references the same `Yhat_expected` already present in the script, so no new `paper_misalignment` is created.
