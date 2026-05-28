---
unit_id: 127
batch: IV.4
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00-06:00
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
  notes_stage_files: [moving_throat_pde_stage127_penetration_families.md]
  paper_appendix: present
---

# Audit unit 127 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_127.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage127_penetration_families.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (single include line at L1288, no additional row text)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage127_penetration_families_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage127_penetration_families_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage127_penetration_families_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage127_penetration_families_mathematica_audit.txt`

## What the paper claims

The card's quoted Output is verbatim: "Slab and truncated-exponential families reach the lower branch at moderate penetration depths." Combined with the notes (which the card defers to), the deliverables are: (i) the uniform slab mouth source `sigma_x^slab(z) = 1/(xL)` on `[0, xL]` yields the closed-form mouth-bias factor `g_slab(x) = (2/(pi x)) sin(pi x / 2)`; (ii) the truncated-exponential source `sigma_x^exp(z) = e^{-z/(xL)} / (xL(1 - e^{-1/x}))` yields `g_exp(x) = 2(2 + pi x e^{-1/x}) / ((4 + pi^2 x^2)(1 - e^{-1/x}))`; (iii) solving `g_slab(x) = g_-^{F1}` gives a unique slab depth `x_*^slab ≈ 0.797839360904564`; (iv) solving `g_exp(x) = g_-^{F1}` gives a unique exponential depth `x_*^exp ≈ 0.662765402623160`; both penetration depths are "moderate" (0.66 and 0.80 of the throat span). The constant `g_-^{F1}` is the lower compensated branch fixed upstream (Stage 122 notes: `g_-^{F1} = (2 sqrt(4107 - 100 pi^2) - 37 sqrt(3)) / (20 pi) ≈ 0.758035078944663`).

## What the script claims to verify

The SymPy and Mathematica scripts construct the closed-form expressions `g_slab(x) = 2 sin(pi x / 2) / (pi x)` and `g_exp(x) = 2(2 + pi x e^{-1/x}) / ((4 + pi^2 x^2)(1 - e^{-1/x}))` as given, construct `g_-^{F1} = (2 sqrt(4107 - 100 pi^2) - 37 sqrt(3)) / (20 pi)`, then numerically root-find `g_slab(x*) = g_-^{F1}` near 0.8 and `g_exp(x*) = g_-^{F1}` near 0.66 at 80-digit precision. The final assertion checks the residuals `g_slab(x_*^slab) - g_-^{F1}` and `g_exp(x_*^exp) - g_-^{F1}` are below `1e-20` in magnitude. The scripts do NOT independently derive the closed forms from the underlying source integrals; they take the closed forms as given.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Slab closed form `g_slab(x) = (2/(pi x)) sin(pi x / 2)` | `g_slab = 2*sin(pi*x/2)/(pi*x)` in both scripts | match (closed form taken as given; integral not re-derived in-script — note below) |
| Exp closed form `g_exp(x) = 2(2 + pi x e^{-1/x}) / ((4 + pi^2 x^2)(1 - e^{-1/x}))` | `g_exp = 2*(2 + pi*x*exp(-1/x))/((4 + pi^2*x^2)*(1 - exp(-1/x)))` | match (closed form taken as given) |
| `x_*^slab ≈ 0.797839360904564` (slab Family-1 depth) | `nsolve(g_slab == gminus, 0.8)` → `0.79783936090456440508...`; SymPy residual `≈ -2.67e-83`; Mathematica residual `≈ 1.17e-58` | match |
| `x_*^exp ≈ 0.662765402623160` (exp Family-1 depth) | `nsolve(g_exp == gminus, 0.66)` → `0.66276540262316140251...`; SymPy residual `≈ -9.67e-82`; Mathematica residual `≈ -3.25e-58` | match |
| Constant `g_-^{F1} = (2 sqrt(4107 - 100 pi^2) - 37 sqrt(3)) / (20 pi)` (anchored Stage 122 notes) | `gminus = (2*R - 37*sqrt(3))/(20*pi)` with `R = sqrt(4107 - 100*pi^2)` | match (anchored to upstream) |
| Positivity / "moderate" framing | Implicitly satisfied by reported values (0.66 and 0.80, both in `(0, 1]`); no explicit positivity assertion | partial (no separate positivity assertion, but the reported roots are visibly in the "moderate" regime) |

`paper_alignment: aligned` — all four numerical deliverables and both closed forms match.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 38–39 | `if abs(res_slab) > 1e-20 or abs(res_exp) > 1e-20: raise AssertionError` | x_*^slab and x_*^exp are roots of `g_{slab,exp}(x) = g_-^{F1}` | yes (non-tautological: `gminus` is constructed independently from radicals; root is found by nsolve) |
| A2 | mathematica | 47 | `expectApprox["slab compensation root", gSlab /. x -> xStarSlab, gMinus, 1e-20]` | x_*^slab is a root of `g_slab(x) = g_-^{F1}` | yes |
| A3 | mathematica | 48 | `expectApprox["exponential compensation root", gExp /. x -> xStarExp, gMinus, 1e-20]` | x_*^exp is a root of `g_exp(x) = g_-^{F1}` | yes |

The assertions are non-tautological — they require the closed forms `g_slab`, `g_exp` and the algebraic constant `gminus` to be self-consistent under the numerical root finder. If, e.g., `gminus` were defined wrong (different branch), nsolve would either fail to converge or converge to a residual much larger than `1e-20`.

## Findings

### F1 — paper_misalignment

**Subtype:** notes_contradicts_script (cosmetic — banner label only)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage127_penetration_families_sympy_audit.py:12` quote: `banner("STAGE 110 — GEOMETRIC MOUTH-PENETRATION FAMILIES")`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage127_penetration_families_mathematica_audit.wl:26` quote: `banner["STAGE 110 — GEOMETRIC MOUTH-PENETRATION FAMILIES"];`
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_127.tex:1` quote: `\section[Stage 127]{Stage 127: Geometric Mouth-Penetration Families}`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage127_penetration_families.md:2` quote: `# Moving-Throat PDE — Stage 229: Geometric Mouth-Penetration Families`

**What's wrong:**
Three different stage numbers appear in this unit's documentation chain: the paper card and filename say "Stage 127", both script banners say "STAGE 110", and the source notes file title says "Stage 229". The transcript output also echoes "STAGE 110". The Mathematica end-of-run message `"Stage 127 Mathematica audit passed."` (L51) agrees with the paper card, and no other stage in the 12x range uses the "STAGE 110" banner — strongly suggesting "STAGE 110" is a copy-paste artifact, while "Stage 229" in the notes appears to be an old internal numbering pre-renumbering.

**Why this matters:**
Not load-bearing for the math, but the inconsistency confuses log-reading and any downstream tool that scrapes the banner line as a unit identifier. The discrepancy is also a tell that the script was duplicated from another file without full review. Per audit policy, paper-vs-script label disagreements are paper_misalignment and require user resolution before Codex acts.

**Required change:**
(For paper_misalignment, see the `## Resolve before fix_loop` block in the directive — Codex must not auto-resolve.)

**Verification:**
After user resolution and the chosen edits land, the next saved transcripts should print the canonical stage number consistently in the banner, the paper card title, the notes file title, and the Mathematica end-of-run message. No numerical change expected.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage127_penetration_families_mathematica_audit.wl:31–48`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage127_penetration_families_sympy_audit.py:14–39`

**What's wrong:**
The Mathematica script is a structural line-by-line port of the SymPy script — same variable choreography, same closed-form definitions (taken as given rather than re-derived), same numeric initial guesses, same final residual check. Side-by-side:

- SymPy L14: `R = sp.sqrt(4107 - 100*sp.pi**2)`; Mathematica L31: `rDisc = Sqrt[4107 - 100*Pi^2]`.
- SymPy L15: `gminus = sp.simplify((2*R - 37*sp.sqrt(3))/(20*sp.pi))`; Mathematica L32: `gMinus = N[(2*rDisc - 37*Sqrt[3])/(20*Pi), 80]`.
- SymPy L18: `g_slab = sp.simplify(2*sp.sin(sp.pi*x/2)/(sp.pi*x))`; Mathematica L34: `gSlab = FullSimplify[2*Sin[Pi*x/2]/(Pi*x), ...]`.
- SymPy L22: `g_exp = sp.simplify(2*(2 + sp.pi*x*sp.exp(-1/x))/((4 + sp.pi**2*x**2)*(1 - sp.exp(-1/x))))`; Mathematica L35: `gExp = FullSimplify[2*(2 + Pi*x*Exp[-1/x])/((4 + Pi^2*x^2)*(1 - Exp[-1/x])), ...]`.
- SymPy L26–27 use `nsolve(...)` with initial guesses 0.8 and 0.66; Mathematica L40–41 use `FindRoot[...]` with initial guesses 0.8 and 0.66.

Neither engine independently derives the closed forms `g_slab(x)` and `g_exp(x)` from the source integrals `int_0^{xL} (1/(xL)) cos(pi z/(2 L)) dz` and `int_0^L sigma_x^exp(z) cos(pi z/(2 L)) dz` defined in the notes. Both scripts copy the closed forms from the notes verbatim and then verify that those expressions, when evaluated at numerically-found roots, equal the upstream constant.

**Why this matters:**
The second-engine policy requires that both engines derive the load-bearing result independently from the physical premises. Here, neither script independently exercises the integral derivation — the only computation done in either engine is `FullSimplify` / `sp.simplify` on already-simplified closed forms, plus a numeric root find. A typo or sign error in the closed form `g_slab` or `g_exp` would propagate to both scripts identically and pass both audits.

**Required change:**
In the Mathematica script, add an independent symbolic integration step for both `g_slab` and `g_exp`, then compare against the closed forms used for root-finding. Specifically:

- After `Clear[x]` and before `gSlab = FullSimplify[...]`, introduce a fresh integration variable `z` and compute:
  ```
  Clear[x, z];
  $Assumptions = Element[x, Reals] && x > 0 && Element[z, Reals];
  gSlabFromIntegral = FullSimplify[Integrate[(1/x)*Cos[Pi*z/2], {z, 0, x}], Assumptions -> x > 0];
  gExpFromIntegral = FullSimplify[Integrate[(Exp[-z/x]/(x*(1 - Exp[-1/x])))*Cos[Pi*z/2], {z, 0, 1}], Assumptions -> x > 0];
  ```
  (Setting `L = 1` is acceptable since L drops out of the bias factor — all instances of `L` cancel between the source normalization and the integration variable change `u = z/L`.)
- Then assert `FullSimplify[gSlabFromIntegral - gSlab] === 0` and `FullSimplify[gExpFromIntegral - gExp] === 0` BEFORE the root-finding step, via the existing `pass`/`fail` helpers (e.g., `If[TrueQ[residual === 0], pass[...], fail[..., residual]]`).

**Verification:**
The new `.wl` script should contain explicit `Integrate[...]` calls for both source profiles, and the symbolic comparison `FullSimplify[gSlabFromIntegral - gSlab] == 0` and `FullSimplify[gExpFromIntegral - gExp] == 0` should appear as new PASS lines in the output. Mathematica's `Integrate` does evaluate the slab cosine integral and the exponential-cosine integral in closed form, so this is a genuine independent derivation rather than another transliteration.

## Independent-derivation check (Mathematica)

See F1. The Mathematica script is a structural line-by-line port of the SymPy script: same closed-form definitions, same numeric initial guesses (0.8 and 0.66), same final assertion shape. Neither engine performs the source-integral derivation. Quoted side-by-side excerpts are in F1.

## Engine cross-check

Both engines compute the same numerical roots to high precision:

| Quantity | SymPy | Mathematica |
|---|---|---|
| `g_-^{F1}` | `0.758035078944662826919680890414` | `0.7580350789446628269196808904141104577505...` |
| `x_*^slab` | `0.79783936090456440508490881279972941789601018958900271741142733512870884256665312` | `0.7978393609045644050849088127997294178960101895890027174113...` |
| `x_*^exp` | `0.66276540262316140251258875013529026617855113761165291470424652438625727951066625` | `0.6627654026231614025125887501352902661785511376116529147062...` |
| residual slab | `-2.67e-83` | `1.17e-58` |
| residual exp | `-9.67e-82` | `-3.25e-58` |

The numbers agree to all displayed digits. Both residuals are well below the asserted `1e-20` tolerance. Engines agree. (The Mathematica residuals are larger only because of `WorkingPrecision -> 80` vs SymPy's `prec=80` + `tol=1e-30` — not a disagreement.)

## Verdict justification

`verdict: findings`, `stop_cold: null`. The paper-side claim is fully exercised: both closed-form mouth-bias factors match the notes verbatim, the constant `g_-^{F1}` is correctly anchored to the Stage 122 algebraic value `(2 sqrt(4107 - 100 pi^2) - 37 sqrt(3)) / (20 pi)`, and the two compensation depths `x_*^slab ≈ 0.797839360904564` and `x_*^exp ≈ 0.662765402623160` reproduce the notes to all displayed digits. The assertions are non-tautological (a wrong branch for `gminus` would not converge), and both engines agree.

Two findings remain. F1 is a cosmetic `paper_misalignment` — three different stage numbers appear in this unit's documentation chain (paper card "Stage 127", script banners "STAGE 110", notes file title "Stage 229"); no numerical impact but user resolution required to pick the canonical number. F2 is a genuine `mathematica_transliteration`: neither engine derives the closed-form bias factors from the source integrals, so a copy-paste error in the closed form would slip through both audits.

Attacks tried that did NOT find a problem:
- Re-derived `g_slab(x) = (2/(pi x)) sin(pi x / 2)` from `int_0^{xL} (1/(xL)) cos(pi z / (2 L)) dz` by hand — matches.
- Re-derived `g_exp(x) = 2(2 + pi x e^{-1/x}) / ((4 + pi^2 x^2)(1 - e^{-1/x}))` from `int_0^L (e^{-z/(xL)} / (xL(1 - e^{-1/x}))) cos(pi z / (2 L)) dz` using the standard `int_0^1 e^{-au} cos(bu) du` evaluation — matches.
- Checked that `g_-^{F1} = (2 sqrt(4107 - 100 pi^2) - 37 sqrt(3)) / (20 pi) ≈ 0.7580351` is the same upstream constant established in Stage 122 notes — matches.
- Checked that the initial guesses for nsolve / FindRoot (0.8 and 0.66) are close to the true roots, so root-finding is well-posed; both engines converge to ~80-digit precision.
- Checked symbol assumptions: `x` declared `positive=True, real=True` in SymPy and `Element[x, Reals] && x > 0` in Mathematica — both match the notes' constraint `0 < x <= 1` for slab and `x > 0` for exp; no positivity is misused.
- Checked staleness: outputs (2026-05-11) are newer than scripts (sympy 2026-04-01; mathematica 2026-05-11), so `outputs_fresh: true`.

No `stop_cold` warranted. F1 is a real second-engine-policy violation but Codex can patch the Mathematica script to add the independent integration step without breaking anything. F2 is trivial. No downstream invalidation since this unit's numerical results are unchanged.

## Self-test notes

- For F1, mentally walked through the Mathematica `Integrate` of `Cos[Pi*z/2]` and `Exp[-z/x] Cos[Pi*z/2]` over their respective domains; both have closed forms that match the notes' boxed expressions, so the new assertions will pass on a correctly-coded patch.
- For F2, the change is a pure string edit; no derivative, integration, or simplification logic is touched, and no symbol redeclaration is needed.
- No `paper_misalignment` introduced by F1's proposed fix: the new symbolic integration verifies the SAME closed forms the paper notes give, just via an independent route.
