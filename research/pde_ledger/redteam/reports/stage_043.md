---
unit_id: 043
batch: III.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: false
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage043_support_direction_extraction.md"]
  paper_appendix: present
---

# Audit unit 043 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_043.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage043_support_direction_extraction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 64)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage043_support_direction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage043_support_direction_mathematica_audit.txt`

## What the paper claims

The paper card (`stage_043.tex`) is terse. Its `\stagefield{Output}` reads:
"A physical support-direction insertion rule \eqref{eq:app-stage043-replacements}", referring to the body equation that prescribes the continuum-selected replacements `q = t R_U, r = t R_phi, m = M_mix, n = M_supp` with `t = kappa_1/kappa_0`. The card asks the binary question `R_phi = R_U` vs `R_phi != R_U` (tracking vs. genuine rank-2 branch). The notes (Sections 1–7) elaborate seven deliverables: (i) rank-1 effective support vector `y = g_B v + g_U g_S D_U v`, (ii) `R_phi = [1 + sigma_0/(1+delta_U)]/(1 + sigma_0)`, (iii) source-tied iff `sigma_0 = 0` (or `delta_U = 0`) via `D_phi := kappa_0 y_1 - kappa_1 y_0`, (iv) `eps_phi^(split) = eps_phi [1 - (2/11) delta_U/(1+delta_U)]` with `kappa_1^2/sigma = 2/11`, (v) `M_supp = 8 Z_phi (1+sigma_0)^2 / [pi^2 (1-eps_eta)(1-eps_phi^(split))]`, (vi) tracking iff `g_B g_R = g_W g_S` i.e. `sigma_0 = rho_0` (via `D_(phi z) := y_0 z_1 - y_1 z_0`), (vii) mismatch `R_phi - R_U = delta_U (rho_0 - sigma_0)/[(1+delta_U)(1+rho_0)(1+sigma_0)]`. The appendix row marks the stage `Reduced/ExactClosure`.

## What the script claims to verify

Both scripts cover the same seven deliverables, in five sections matching the notes layout. The SymPy script's docstring enumerates exactly the same seven items, and the load-bearing `expect_zero(...)` / `expectZero[...]` assertions are: `R_phi - R_phi_expected`, `D_phi - D_phi_expected`, the overlap-contraction identity `v.D_U.v - (sigma/K_U)(1 - (2/11) delta_U/(1+delta_U))`, `A_phi^(eff)` and its minimal-limit reduction, the split-vs-minimal overlap ratio, `M_supp` independence from bare-mode masses and structural form (with free baseline `B` and at `B = 8/pi^2`), `D_(phi z) - expected`, the tracking specialization `g_B g_R = g_W g_S`, and the mismatch formula. The Mathematica script adds endpoint checks (`v.D_U.v` at `delta_U = 0` and `delta_U -> Infinity`) plus a leading-order Series check on the mismatch.

## Paper -> script cross-check

| paper-side deliverable | script-side check | status |
|---|---|---|
| Insertion rule `r = t R_phi` (R_phi formula) | `R_phi - R_phi_expected` (sympy L74, math L58) | match |
| `D_phi` and source-tied iff `sigma_0 = 0` or `delta_U = 0` | `D_phi - D_phi_expected` (sympy L81, math L59) | match (sign convention differs in Mathematica — see F2) |
| Overlap contraction and `eps_phi^(split)` | `v.D_U.v` and `A_phi^(eff)` checks (sympy L93,L101,L115; math L85–88,L98) | match |
| `M_supp` baseline `8 Z_phi (1+sigma_0)^2 / [pi^2(1-eps_eta)(1-eps_phi^split)]` | mu-independence + structural form + value at `B = 8/pi^2` (sympy L127–148, math L110–129) | match |
| Mixed-direction `D_(phi z)` and tracking iff `g_B g_R = g_W g_S` | `D_(phi z) - expected` and `R_phi = R_U` specialization (sympy L158,L162; math L144–145) | match |
| Mismatch formula `R_phi - R_U` | `mismatch - mismatch_expected` (sympy L170; math L160) | match |
| Insertion of `m = M_mix` (paper card eq. 2) | (none — `M_mix` is carried forward from Stage 22/041; paper card explicitly says so via `\stagefield{Inputs}`) | not-required-here |

`paper_alignment: aligned` — every paper-side deliverable that originates in this stage has a matching script-side check; the only paper symbol the script does not re-derive is `M_mix`, which the paper card itself lists as an Input rather than an Output.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 74 | `expect_zero("R_phi - expected", Rphi - Rphi_expected)` | R_phi formula (notes (ii)) | yes |
| A2 | sympy | 81 | `expect_zero("D_phi - expected", Dphi - Dphi_expected)` | direction-splitting invariant (notes (iii)) | yes |
| A3 | sympy | 93 | `expect_zero("support overlap contraction", ...)` | overlap contraction (notes (iv)) | yes |
| A4 | sympy | 101 | `expect_zero("A_phi^(eff) - expected", ...)` | effective pole (notes (iv)) | yes |
| A5 | sympy | 108 | `expect_zero("A_phi^(eff) at delta_U=0 (minimal)", ...)` | minimal-limit consistency | yes |
| A6 | sympy | 115 | `expect_zero("split-vs-minimal overlap ratio", ...)` | overlap-ratio identity | yes |
| A7 | sympy | 127–128 | `expect_zero("M_supp independent of mu_eta/mu_phi", ...)` | M_supp baseline (mass cancellation) | yes |
| A8 | sympy | 141 | `expect_zero("M_supp structural form (free baseline)", ...)` | M_supp structural form (notes (v)) | yes |
| A9 | sympy | 148 | `expect_zero("M_supp at baseline B = 8/pi^2", ...)` | M_supp at kappa_0^2 = 8/pi^2 | partial (literal 8/pi^2 imported from upstream sigma = 88/(9 pi^2), kappa_1^2/sigma = 2/11; consistent with notes Section 4) |
| A10 | sympy | 158 | `expect_zero("D_(phi z) - expected", ...)` | tracking invariant (notes (v),(vi)) | yes |
| A11 | sympy | 162 | `expect_zero("tracking condition via g_B g_R = g_W g_S", ...)` | tracking iff (notes (vi)) | yes |
| A12 | sympy | 170 | `expect_zero("mismatch formula", ...)` | mismatch formula (notes (vii)) | yes |
| B1 | mathematica | 58 | `expectZero["R_phi - expected", rPhi - rPhiExpected]` | R_phi formula | yes |
| B2 | mathematica | 59 | `expectZero["D_phi - expected", dPhi - dPhiExpected]` | direction-splitting invariant | yes (but `dPhi` uses opposite sign convention from notes — see F2) |
| B3 | mathematica | 85–86 | endpoint checks `delta_U = 0`, `delta_U -> Infinity` | overlap contraction (independent check) | yes |
| B4 | mathematica | 87 | `expectZero["support overlap contraction", ...]` | overlap contraction | yes |
| B5 | mathematica | 88 | `expectZero["A_phi^(eff) - expected", ...]` | effective pole | yes |
| B6 | mathematica | 93 | `expectZero["A_phi^(eff) at deltaU=0 (minimal)", ...]` | minimal-limit consistency | yes |
| B7 | mathematica | 98 | `expectZero["split-vs-minimal overlap ratio", ...]` | overlap-ratio identity | yes |
| B8 | mathematica | 110–111 | `expectZero["M_supp independent of muEta/muPhi", ...]` | M_supp mass independence | yes |
| B9 | mathematica | 123 | `expectZero["M_supp structural form (free baseline)", ...]` | M_supp structural form | yes |
| B10 | mathematica | 129 | `expectZero["M_supp at baseline B = 8/Pi^2", ...]` | M_supp value | partial (same 8/Pi^2 literal) |
| B11 | mathematica | 144 | `expectZero["D_(phi z) - expected", ...]` | tracking invariant | yes |
| B12 | mathematica | 145 | `expectZero["tracking condition via g_B g_R = g_W g_S", ...]` | tracking iff | yes |
| B13 | mathematica | 157 | `expectZero["mismatch leading-in-deltaU coefficient", ...]` | mismatch (leading-order check) | yes |
| B14 | mathematica | 160 | `expectZero["mismatch formula", ...]` | mismatch formula | yes |

All listed assertions are non-tautological in the strict sense (the LHS is computed by a different route than the RHS; equality is a real algebraic claim about the structure of `y`, `D_U`, and the substitutions). One assertion pair (A9 / B10) imports the value `kappa_0^2 = 8/pi^2` as a literal; this value is consistent with the notes' upstream identities `sigma = 88/(9 pi^2)`, `kappa_1^2/sigma = 2/11`, but is not re-derived in this stage. The script structures A8 / B9 around a free baseline `B`, then identifies `B = 8/pi^2` separately — the free-baseline check is independent of the literal, so the hardcoded number is acceptably isolated.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl:43-160`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py:59-170`

**What's wrong:**
The Mathematica script is structurally a near line-for-line port of the SymPy script: same variable choreography, same intermediate names (with case-only changes: `DU` -> `dU`, `Rphi` -> `rPhi`, etc.), same algebraic recipe in the same order. Both engines build `y = gB v + gU gS D_U v`, define `sigma_0 = gU gS/(KU gB)`, take the same ratio `(y1/k1)/(y0/k0)`, and then verify against the same closed-form expected. Both substitute `kappa_1^2 -> (2/11) sigma`, `kappa_0^2 -> (9/11) sigma` (sympy: L89; math: L65–67). Both isolate `M_supp` via the same free-baseline `B = kappa_0^2` substitution and the same `B = 8/pi^2` evaluation (sympy: L131–148; math: L116–129). Three corresponding excerpts:

Section 1 (effective support vector):

sympy:
```
DU = sp.diag(1 / KU, 1 / (KU * (1 + delta_U)))
v = sp.Matrix([kappa0, kappa1])
y = sp.simplify(gB * v + gU * gS * DU * v)
sigma0 = sp.simplify(gU * gS / (KU * gB))
Rphi = sp.simplify((y1 / y0) / (kappa1 / kappa0))
Rphi_expected = sp.simplify((1 + sigma0 / (1 + delta_U)) / (1 + sigma0))
```

mathematica:
```
dU = DiagonalMatrix[{1/kU, 1/(kU (1 + deltaU))}];
v = {kappa0, kappa1};
y = FullSimplify[gB v + gU gS (dU.v), Assumptions -> $Assumptions];
sigma0 = FullSimplify[gU gS/(kU gB), Assumptions -> $Assumptions];
rPhi = FullSimplify[(y[[2]]/kappa1) / (y[[1]]/kappa0), Assumptions -> $Assumptions];
rPhiExpected = FullSimplify[(1 + sigma0/(1 + deltaU))/(1 + sigma0), Assumptions -> $Assumptions];
```

Section 3 (M_supp structural form, free baseline):

sympy:
```
Msupp_cont_in_B = sp.simplify(Msupp_cont.subs(kappa0**2, B))
Msupp_struct_expected = sp.simplify(
    B
    * (cB**2 / (Keta_eff * Kphi_eff))
    * (1 + ceU * cUphi / (KU * cB))**2
    / ((1 - eps_eta) * (1 - eps_phi_split.subs(eps_phi, cUphi**2 * sigma / (KU * Kphi_eff))))
)
```

mathematica:
```
mSuppContInB = FullSimplify[mSuppCont /. kappa0^2 -> bBaseline, Assumptions -> $Assumptions];
mSuppStructExpected = FullSimplify[
  bBaseline (cB^2/(kEtaEff kPhiEff)) (1 + cEtaU cUphi/(kU cB))^2/
    ((1 - epsEta) (1 - epsPhiSplitPhys)),
  Assumptions -> $Assumptions
];
```

Section 5 (mismatch):

sympy:
```
mismatch = sp.simplify(Rphi_expected - RU)
mismatch_expected = sp.simplify(delta_U * (rho0 - sigma0) / ((1 + delta_U) * (1 + rho0) * (1 + sigma0)))
```

mathematica:
```
mismatch = FullSimplify[rPhiExpected - rU, Assumptions -> $Assumptions];
mismatchExpected = FullSimplify[
  deltaU (rho0 - sigma0)/((1 + deltaU) (1 + rho0) (1 + sigma0)),
  Assumptions -> $Assumptions
];
```

The Mathematica script does add two genuinely independent decorations (endpoint limits at `delta_U = 0` and `delta_U -> Infinity` on `v.D_U.v` at L79–86, and a `Series` leading-order coefficient check on the mismatch at L154–157). These are real second-engine value-adds. But the load-bearing assertions — R_phi, D_phi, A_phi^(eff), M_supp structural form, D_(phi z), mismatch formula — are computed by the same algebraic recipe in the same order using parallel variable names. The second engine is not independently rederiving the result; it is re-executing the SymPy script's algebra in Mathematica syntax.

**Why this matters:**
The second-engine policy exists to catch cases where the first engine has a hidden bug (mis-typed expected, wrong sign, sloppy `simplify` masking a residual). A transliterated second engine cannot catch any of those: it builds the same intermediate from the same recipe and finds the same answer. The endpoint and series checks are useful supplements but are not a substitute for an independent derivation.

**Required change:**
Augment the Mathematica script with at least one genuinely independent derivation route for the load-bearing claims. Concrete options Codex can pick from (any one of these would suffice; do not refactor the existing flow):
- For R_phi: solve a small linear system in the wall-basis components instead of pre-eliminating `U` symbolically. E.g. define `Eta`, `phi`, `U0`, `U1` as separate variables, write down the four reduced static equations from the notes (`L_(eta U) = -g_U Q.U`, `L_(eta phi) = -g_B (v.Q) Phi`, `L_(Uphi) = +g_S (v.U) Phi`), `Solve[]` the system for `(U0, U1)` in terms of `(Q, Phi)`, substitute back, and read off the residual `eta-phi` coupling. Compare the resulting `y` vector against the current `y` symbolically.
- For D_(phi z): instead of `Det[{{y0,y1},{z0,z1}}]`, write `z = gW v + gU gR D_U v` and derive the cross-product `y0 z1 - y1 z0` by expanding the bilinear in `gB, gW, gS, gR` and tracking the four resulting terms; verify the only surviving terms factor through `(gB gR - gW gS)`.
- For the mismatch: take `R_phi - R_U` via partial-fraction expansion in `(1+sigma_0)(1+rho_0)` rather than direct subtraction; check the numerator equals `delta_U (rho_0 - sigma_0)`.

Choose ONE of these; the goal is a second derivation route, not a refactor of all existing checks. The current parallel checks should remain (they confirm the simpler route still passes); the new route exercises the same identity by a structurally different path.

**Verification:**
After the patch, the Mathematica script must (a) still pass all current `expectZero[...]` checks, (b) contain a new `expectZero[...]` whose LHS is built by a derivation route absent from the SymPy script (e.g. a `Solve[]` call where SymPy uses direct substitution), and (c) the new check's variable names should NOT mirror the SymPy script's. The transcript should show the new check passing.

### F2 — paper_misalignment

**Subtype:** notes_contradicts_script

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage043_support_direction_extraction.md:119-121`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl:52-53`

**What's wrong:**
The notes explicitly define the support direction-splitting invariant as

> `D_phi := kappa_0 y_1 - kappa_1 y_0`
> `      = - kappa_0 kappa_1 g_B sigma_0 delta_U / (1 + delta_U).`

(notes lines 119–121). The SymPy script honors this convention: `Dphi = sp.factor(sp.expand(kappa0 * y1 - kappa1 * y0))` and `Dphi_expected = -kappa0 * kappa1 * gB * sigma0 * delta_U / (1 + delta_U)` (sympy L77–78). The printed SymPy output reads `D_phi = -delta_U*g_S*g_U*kappa0*kappa1 / (K_U*(delta_U+1))`, consistent with the notes' sign.

The Mathematica script reverses the sign:
```
dPhi = FullSimplify[Det[{{y0, y1}, {kappa0, kappa1}}], ...];
dPhiExpected = FullSimplify[kappa0 kappa1 gB sigma0 deltaU/(1 + deltaU), ...];
```
`Det[{{y0, y1}, {kappa0, kappa1}}] = y0*kappa1 - y1*kappa0 = -(kappa_0 y_1 - kappa_1 y_0)`, i.e. the negative of the notes' definition. The Mathematica `dPhiExpected` is also given without the leading minus sign, so the check passes internally — but the printed transcript reads `D_phi = (deltaU*gS*gU*kappa0*kappa1)/(kU + deltaU*kU)`, which has the OPPOSITE sign from the notes' definition and from the SymPy transcript.

The substantive claim `D_phi = 0 iff sigma_0 = 0 or delta_U = 0` is unaffected (the vanishing set is symmetric under sign), but the printed object labeled `D_phi` in the Mathematica transcript is `-D_phi` per the notes' convention, and the two engines' transcripts disagree on the sign of this named quantity.

**Why this matters:**
A reader cross-checking the Mathematica transcript against the notes will see a sign flip and may suspect a real error in either the notes or the script. The two engines do not agree on the printed sign of a load-bearing intermediate; this is exactly the kind of disagreement the second-engine policy is supposed to surface, not hide.

**Required change:**
This is a `paper_misalignment` (subtype `notes_contradicts_script`). Direction of resolution is the user's call — see `## Resolve before fix_loop` block in the directive. Codex must not silently flip either side.

**Verification:**
After user resolution, either (a) Mathematica L52 changes from `Det[{{y0, y1}, {kappa0, kappa1}}]` to `Det[{{kappa0, kappa1}, {y0, y1}}]` (equivalent to the notes' `kappa_0 y_1 - kappa_1 y_0`) and `dPhiExpected` gets a leading minus sign so the printed `D_phi` matches notes; or (b) the notes are updated to define `D_phi := y_0 kappa_1 - y_1 kappa_0` (sign-flipped). Either way, the two engines and the notes must agree on the sign of the printed `D_phi`.

## Independent-derivation check (Mathematica)

The Mathematica script is a near-transliteration of the SymPy script. The two endpoint checks (L79–86) and the Series check (L154–157) are independent, but the bulk of the load-bearing assertions follow the same recipe in the same order with parallel variable names. See F1 for the line-by-line excerpts.

## Engine cross-check

Both engines exit cleanly with `PASS:` lines on every assertion. The printed final residuals (after subtracting expected) are all `0`. However:
- `D_phi` printed values disagree in sign:
  - sympy: `D_phi = -delta_U*g_S*g_U*kappa0*kappa1 / (K_U*(delta_U+1))`
  - mathematica: `D_phi = (deltaU*gS*gU*kappa0*kappa1)/(kU + deltaU*kU)`

  This is F2 — sign-convention mismatch hidden by each engine's `expected` having the matching sign.
- `D_(phi z)`, `R_phi`, `A_phi^(eff)`, `M_supp` agree across engines.
- `mismatch (R_phi - R_U)` agrees (sympy: `K_U*delta_U*g_U*(g_B*g_R - g_S*g_W) / [(delta_U+1)*(K_U*g_B + g_S*g_U)*(K_U*g_W + g_R*g_U)]`; mathematica: `(deltaU*gU*(gB*gR - gS*gW)*kU)/((1+deltaU)*(gS*gU + gB*kU)*(gR*gU + gW*kU))` — algebraically identical).

`engines_agree: false` in the front-matter because of the sign convention on `D_phi`; the underlying math is the same, the printed/named quantity differs.

## Verdict justification

The stage substantively holds: every paper-side deliverable that originates in this stage has a corresponding non-tautological script check, and both engines reach `0` on every load-bearing identity. Attacks tried that failed: (a) tried to find a tautology — every `expected` is built from a different route than the LHS, none are algebraically identical by construction; (b) checked whether `kappa_0^2 -> B` substitution hides the structural claim — no, the script verifies the structural form with a free baseline `B` and only afterwards specializes `B = 8/pi^2`, so the structural identity is exercised symbolically; (c) checked whether `(2/11)`/`(9/11)` literal substitution is tautological — no, these come from upstream Stage 22 facts about `kappa_1^2/sigma`, not from this stage; (d) checked symbol assumptions — `kappa0, kappa1` are declared `real, nonzero`, `delta_U > 0`, all consistent with the physical setup; (e) checked sign conventions on `D_phi` and `D_(phi z)` — `D_phi` mismatches between engines (F2), `D_(phi z)` agrees.

The two findings are real but bounded: F1 (transliteration) is medium-severity because both engines reach the same answer but by the same recipe; F2 (sign of `D_phi`) is low-severity because the vanishing condition (the actual physical claim) is unaffected. The paper card and the script's load-bearing claims align.

Verdict: `findings`. Not `stop_cold` — neither finding propagates to invalidate downstream units. `paper_alignment: aligned` (the paper card's `\stagefield{Output}` insertion rule is exactly what the script verifies, with `M_mix` correctly excluded as an Input).

## Self-test notes

Traps checked: (1) Variable-independence trap — confirmed `M_supp` truly depends on `mu_eta, mu_phi` before cancellation (each appears in the numerator `1/(mu_eta mu_phi)` and in the denominators `Keta_eff/mu_eta` and `Kphi_eff/mu_phi`), so the `D[Msupp, mu_eta] = 0` check is non-trivial. (2) Sign / parity — confirmed `D_phi`'s sign convention differs between engines (-> F2); confirmed `D_(phi z)` signs agree. (3) Trivial-case pre-check — substituted `sigma_0 = 0` mentally into `R_phi = (1+sigma_0/(1+delta_U))/(1+sigma_0)` -> `1`, then into `D_phi` -> `0`, both consistent with the source-tied closure claim; substituted `sigma_0 = rho_0` into the mismatch formula -> `0`, consistent with tracking. (4) Path specifications — F1's directive change is in `mathematica/`, not `scripts/`. (5) Paper round-trip — F1's proposed change introduces a `Solve[]`-based route; this does not introduce new symbols beyond those already in the notes, and the paper card does not constrain the derivation method, only the final identities.
