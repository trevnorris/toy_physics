---
unit_id: 211
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction.md]
  paper_appendix: present
---

# Audit unit 211 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_211.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (read rows 49–65, 236, 909–916, 920–994, 1304–1345)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.txt`

## What the paper claims

The card's `\stagefield{Output}` is verbatim: *"Finite triple-interior candidate set, exact interior winner/non-improvement theorems, and special diagonal-isotropic and full-triple-symmetry reductions."* The stage solves the two-parameter interior optimizer for a primitive triple. The notes enumerate the deliverables in detail: (1) the exact two-component stationary numerator system `N_r,* = 2 M_r sqrt(Delta) + L_r`, `N_s,* = 2 M_s sqrt(Delta) + L_s` arising from `dPhi/dr`, `dPhi/ds`; (2) exact square-root elimination into one **quartic** cross-consistency polynomial `C = M_s L_r - M_r L_s` and two **sextic** square eliminants `S_r = L_r^2 - 4 M_r^2 Delta`, `S_s = L_s^2 - 4 M_s^2 Delta`; (3) the finite algebraic pre-candidate set with Bezout bound `4*6 = 24`; (4) two special reductions — diagonal-isotropic curvature (where `tau` depends only on `k_ijk(r,s)` and the optimizer is the gradient-optimal ray) and full triple-symmetry (where the equal-mix barycenter `(1,1,1)/sqrt(3)` is an exact stationary ray, i.e. `N_r(1,1)=N_s(1,1)=0`); (5) the interior winner / non-improvement theorems comparing against the imported Stage 209/210 pairwise brackets. The appendix row (line 53) summarizes: "Reduces the two-parameter triple interior optimizer to quartic/sextic algebraic candidate sets with finite interior brackets."

Crucially, the card's `\stagefield{Verification}` line states verbatim: *"Mathematica audit: none yet."*

## What the script claims to verify

The SymPy script verifies, by symbolic differentiation and algebra: (M1) that the closed-form `dPhi/dr`, `dPhi/ds` reduce to `N_r/[2(1+r^2+s^2)^{3/2} sqrt(Delta)]` etc. with `N_r = 2 M_r sqrt(Delta) + L_r`; (M2) the cross-elimination identity `M_s N_r - M_r N_s = C_cross` and that `C_cross` is total-degree 4; (M3) the square eliminant identities `N_r(N_r - 4 M_r sqrt(Delta)) = S_r` (and s) and that `S_r, S_s` are total-degree 6; (M4) the Bezout product `4*6 = 24`; (M5) the diagonal-isotropic reduction `Delta_iso = (k_i + r k_j + s k_k)^2 - 2 H0 u (1+r^2+s^2)` and `tau_iso = 2H0/(k_rs + sqrt(k_rs^2 - 2 H0 u))`; (M6) the equal-mix stationarity `N_r(1,1) = N_s(1,1) = 0` under the symmetric substitution. The Mathematica `.wl` (added in pass 1 as a retrofit) asserts the identical set M1–M6 with the same names. Neither script evaluates the winner / non-improvement theorems numerically (those are `\StatusNumerical` and require an actual PDE-derived packet, by the stage's own scope).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Stationary numerator system `N_r, N_s` from `dPhi/dr, dPhi/ds` | sympy M1 (L75–76), wl M1 (L96–99) | match |
| Quartic cross-consistency `C` deg 4 | sympy M2 (L94, L106), wl M2 (L106–108) | match |
| Sextic eliminants `S_r, S_s` deg 6 | sympy M3 (L95–96, L108–111), wl M3 (L117–122) | match |
| Bezout bound `4*6 = 24` | sympy V (L167–170), wl M4 (L126–128) | match |
| Diagonal-isotropic reduction → grad-optimal | sympy III (L131–135), wl M5 (L147–148) | match (reduction identity verified; "grad-optimal ray" conclusion is stated, not a separate assertion — same as paper's analytic argument) |
| Full-symmetry equal-mix stationarity `N_r(1,1)=N_s(1,1)=0` | sympy IV (L159–160), wl M6 (L166–167) | match |
| Interior winner / non-improvement theorems | (none — `\StatusNumerical`, requires real packet) | missing by design (in-scope deferral per notes §6) |
| Paper card: "Mathematica audit: none yet" | a `.wl` now EXISTS | mismatch (card stale post-retrofit) |

Both engines verify the same identity set faithfully against the analytic deliverables. The dominant pattern is `match`, but two defects pull alignment to `partial`: the stale `\stagefield{Verification}` line, and the `.wl` being a non-independent port (below).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 75 | `expect_zero(diff(Phi,r) - Nr/denom)` | stationary numerator law (M1) | yes |
| A2 | sympy | 76 | `expect_zero(diff(Phi,s) - Ns/denom)` | stationary numerator law | yes |
| A3 | sympy | 94 | `expect_zero(Ms*Nr - Mr*Ns - Ccross)` | quartic cross-elimination | yes |
| A4 | sympy | 95–96 | `expect_zero(Nr*(Nr-4 Mr sqrtDelta) - Sr)` (and s) | sextic eliminant | yes |
| A5 | sympy | 106–111 | degree asserts (4,6,6) | quartic/sextic claim | yes |
| A6 | sympy | 131–135 | `expect_zero(Delta_iso - ...)`, `tau_iso - ...` | diag-isotropic reduction | yes |
| A7 | sympy | 159–160 | `expect_zero(Nr_sym)`, `expect_zero(Ns_sym)` | equal-mix stationarity | yes |
| A8 | sympy | 167–170 | `bezout_bound == 24` | Bezout bound | yes |
| B1 | mathematica | 96–99 | `expectZero[D[Phi,r] numerator - Nr]`, `D[Phi,r]-Nr/denom` | stationary numerator law | yes (but port of A1/A2) |
| B2 | mathematica | 106–108 | `expectZero[Ms Nr - Mr Ns - Ccross]`, deg==4 | quartic | yes (port of A3/A5) |
| B3 | mathematica | 117–122 | `expectZero[Nr(Nr-4 Mr sqrtDelta) - Sr]`, deg==6 | sextic | yes (port of A4/A5) |
| B4 | mathematica | 126–128 | `bezoutBound == 24` | Bezout bound | yes (port of A8) |
| B5 | mathematica | 147–148 | `Delta_iso`, `tau_iso` reductions | diag-isotropic | yes (port of A6) |
| B6 | mathematica | 166–167 | `Nr(1,1)`, `Ns(1,1)` symmetric | equal-mix | yes (port of A7) |

All assertions are non-tautological and well-anchored to paper claims. The defect is not anchoring — it is that the B-rows replicate the A-rows operation-for-operation rather than re-deriving independently.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl:62-167`
- compared against `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py:42-170`

**What's wrong:**
The `.wl` is a line-for-line port of the `.py`, not an independent second engine. Both posit the identical objects from the same closed forms and then perform the identical extraction operation on the single load-bearing quantity (the stationary numerator), so the second engine adds no independent confirmation of the algebra.

Object definitions are 1:1 (only renames `C→Ccoef`, `D→Dcoef`, `E→Ecoef`, `F→Fcoef` to dodge Mathematica builtins):

- sympy L53–58:
  `Mr = (1+r**2+s**2)*kj - r*(ki+kj*r+kk*s)` … `Nr = 2*Mr*sqrtDelta + Lr`
- wl L87–92:
  `Mr = den kj - r linearK` … `Nr = 2 Mr sqrtDelta + Lr`

The load-bearing derivation — that `dPhi/dr` has numerator `N_r` — is extracted by the SAME operation (symbolic differentiation of the SAME `Phi`, compared to the SAME posited `Nr`):

- sympy L73–76:
  `expected_dPhi_dr = Nr / (2*(1+r**2+s**2)**Rational(3,2)*sqrtDelta)`
  `expect_zero("exact dPhi/dr compiler", sp.diff(Phi, r) - expected_dPhi_dr)`
- wl L82–99:
  `directNumerator[var_] := FullSimplify[Numerator[Together[D[Phi, var]]], ...]`
  `expectZero["M1 D[Phi,r] numerator minus paper N_r", directNumerator[r] - Nr]`
  `expectZero["M1 derivative law r", D[Phi, r] - Nr/stationaryDenominator]`

`sp.diff(Phi,r)` and `D[Phi,r]` are the same autodiff of the same expression. The elimination steps are likewise identical algebra on identical objects:

- sympy L83–96: `Ccross = Ms*Lr - Mr*Ls`; `Sr = Lr**2 - 4*Mr**2*Delta`; `expect_zero("cross-elimination identity", Ms*Nr - Mr*Ns - Ccross)`; `expect_zero("square elimination identity (r)", Nr*(Nr - 4*Mr*sqrtDelta) - Sr)`
- wl L103–118: `Ccross = Expand[Ms Lr - Mr Ls]`; `Sr = Expand[Lr^2 - 4 Mr^2 Delta]`; `expectZero["M2 square-root-free cross identity", Ms Nr - Mr Ns - Ccross]`; `expectZero["M3 square eliminant identity r", Nr (Nr - 4 Mr sqrtDelta) - Sr]`

And both apply byte-identical iso/symmetric substitutions (sympy L119–155 vs wl L132–162) and the same `4*6=24` product. There is no independent route in the `.wl` — Mathematica could have derived the eliminants via `Resultant`/`GroebnerBasis`, derived the candidate count from an actual `Solve`/`NSolve` on a concrete envelope, or derived `tau_iso` by an independent `Series`/`Reduce` rather than re-substituting the same posited form. None of that is done. The "each CAS runs its own FullSimplify" defense does not rescue independence: same premise + same operation = port. The discriminator (derive-vs-posit must differ between engines) fails on the one genuinely derived object (`N_r`), which both engines obtain by the same `diff` of the same `Phi`.

**Why this matters:**
The dual-engine policy exists so a transcription/algebra error in one engine is caught by an independent route in the other. A port catches only Mathematica-vs-SymPy implementation-of-`diff` discrepancies, not errors in the shared posited forms — both would carry the same mistake. The retrofit therefore does not deliver the independent confirmation the policy requires.

**Required change:**
Re-author the `.wl` so at least the load-bearing object is obtained by a route the `.py` does not use. Concretely, the `.wl` should NOT define `Nr`/`Ccross`/`Sr` by re-typing the paper formulas; instead it should (a) compute `D[Phi,r]`, take `Numerator[Together[...]]`, and then *independently rationalize the square root by Mathematica's own elimination* — e.g. obtain the quartic cross-consistency via `Resultant[numR, numS, Sqrt[Delta]]` or `GroebnerBasis`/`Eliminate` over the radical, and obtain the sextic via `Resultant[numR, RADICAL^2-Delta, RADICAL]` — then compare the *resulting* polynomial's total degree to 4 / 6 and its factor structure to the SymPy-posited `C`/`S_r`. The candidate-count `24` should fall out of the derived polynomials' degrees, not a hardcoded `4*6`. This makes the eliminants a Mathematica *derivation* checked against the SymPy *posit*, which is the independent cross-check. (Routes the new content through Codex per `feedback_claude_reviews_codex_codes`; this directive states requirement + acceptance, not a pre-designed script.)

**Verification:**
The re-authored `.wl` must (i) NOT contain the literal `Nr = 2 Mr sqrtDelta + Lr` / `Ccross = Ms Lr - Mr Ls` / `Sr = Lr^2 - 4 Mr^2 Delta` as the source of the eliminants, but derive them via Resultant/Eliminate; (ii) still exit 0 with all PASS lines; (iii) report degree 4 for the derived cross polynomial and 6 for the derived square eliminant; (iv) confirm the derived polynomials agree (up to nonzero factor) with the SymPy `C_cross`/`S_r`.

### F2 — paper_misalignment

**Subtype:** paper_missing_script_claim
**Severity:** low
**Files:**
- paper side: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_211.tex:11`
- script side: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl` (entire file)

**What's wrong:**
The stage card states verbatim (L11): *"SymPy audit: \StageFile{...sympy_audit.py}.  Mathematica audit: none yet."* But a Mathematica audit `.wl` now exists (retrofitted in pass 1, mtime 2026-06-02 10:46) and produces a full PASS transcript. The card is stale: it claims no Mathematica engine when one is present.

**Why this matters:**
The card understates the verification coverage and contradicts the actual repository state. A reader trusting the card would believe stage 211 is single-engine. This is a documentation-vs-script disagreement that the auditor must route to the user (Codex does not edit paper/).

**Required change:**
User decision (see directive `## Resolve before fix_loop`). Almost certainly: update L11 to "Mathematica audit: \StageFile{mathematica/...mathematica_audit.wl}." — but only after F1 is resolved (a port should arguably not be advertised as a second engine until re-authored). The direction and timing are the user's call.

**Verification:**
After user resolution, card L11 names the `.wl` path (and the `.wl` has been re-authored per F1, if the user chooses that order).

## Independent-derivation check (Mathematica)

PORT. See F1. Three corresponding sections, quoted above: (1) object definitions `Mr/Ms/Lr/Ls/Nr/Ns` are 1:1 identical formulas (sympy L53–58 ↔ wl L87–92); (2) the load-bearing stationary numerator is extracted by the same `diff(Phi,·)` autodiff compared to the same posited `Nr` (sympy L73–76 ↔ wl L82–99); (3) the eliminants `Ccross`, `Sr`, `Ss` are defined by identical algebra and checked by identical elimination identities (sympy L83–96 ↔ wl L103–118), with identical iso/symmetric substitutions and the same `4*6=24` product. The single discriminating operation: both posit `N_r = 2 M_r sqrt(Delta) + L_r` and obtain the eliminants by re-typing the paper's closed forms; neither engine *derives* the eliminants by an independent radical-elimination (Resultant/GroebnerBasis). The one derived object — `dPhi/dr`'s numerator — is obtained by the same symbolic-differentiation operation in both.

## Engine cross-check

Both engines agree: every assertion is `= 0` / `True` / PASS. SymPy output: "exact dPhi/dr compiler = 0", "cross-elimination identity = 0", "deg C_cross = 4", "deg S_r = 6", "Bezout bound ... = 24", "Delta_iso reduction = 0", "tau_iso ... = 0", "Nr(1,1) = 0", "Ns(1,1) = 0", "All Stage 211 identities verified." Mathematica output: matching PASS for M1–M6, "M4 Bezout product equals 24 = True", "M5 tau_iso reduction = 0", "M6 symmetric N_r(1,1) = 0", "All Stage 211 identities verified." Numerically/symbolically consistent — but the agreement is expected precisely because the `.wl` is a port (F1), so this consistency does not constitute independent corroboration.

## Verdict justification

The mathematics is correct and both transcripts pass: the stationary numerator law, the quartic/sextic elimination, the Bezout bound 24, and both special reductions all check out, and they faithfully exercise the paper's analytic deliverables (no tautologies — each `expect_zero` differentiates/eliminates real closed forms whose residual could be nonzero if a formula were wrong; the degree and Bezout literals match the documented 4/6/24). Attacks tried and failed: (a) is `dPhi/dr - Nr/denom` tautological? No — `Nr` is built from independent `Mr/Lr` definitions, not from `diff(Phi,r)`, so the subtraction genuinely tests the derivative formula. (b) Are the degree checks self-fulfilling? No — `total_degree()`/`totalDegreeRS` are computed from the expanded polynomial, not posited. (c) Does the symmetric substitution accidentally force `N_r(1,1)=0` for any inputs? No — it depends on the equal-`k` + permutation-symmetric envelope, which is the documented regime, and a wrong `M_r` would leave a residual. (d) Symbol domains: `Delta>0`, `k_*>0`, `H0>0`, `r,s>=0` are all justified by the physical setup. The verdict is `findings` (not clean) for two reasons: F1 — the `.wl` is a transliteration of the `.py` (medium; violates dual-engine independence on a retrofit batch), and F2 — the card's `\stagefield{Verification}` line still says "Mathematica audit: none yet" while a `.wl` exists (low; user-routed paper_misalignment). Neither finding breaks the math, so no stop-cold.

## Value Reconciliation (pass-2 augmentation)

Enumerating every RESULT/deliverable value the scripts emit (excluding pass/fail scaffolding, residual-zero checks, intermediate-only quantities):

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `N_r = 2 M_r sqrt(Delta) + L_r` (stationary numerator) | py L57, wl L91; sympy.txt L25–36 (M_r,L_r forms) | notes L196–198 (`\mathcal N_{r,\star}=2M_r\sqrt{\Delta}+L_r`); card "Output" implies | MATCH |
| `N_s = 2 M_s sqrt(Delta) + L_s` | py L58, wl L92 | notes L200–203 | MATCH |
| `M_r = k_j(1+s^2) - r(k_i + s k_k)` | py L53; sympy.txt L25–27 | notes L153–155 | MATCH |
| `M_s = k_k(1+r^2) - s(k_i + r k_j)` | py L54; sympy.txt L28–30 | notes L157–160 | MATCH |
| `C_cross = M_s L_r - M_r L_s` (quartic) | py L83, wl L103; sympy.txt L43–66 | notes L234–238; appendix L981 | MATCH |
| `deg C_cross = 4` | py L102/106, wl L107; sympy.txt L208; wl.txt L25 | notes L247 ("**quartic**"); appendix L981 | MATCH |
| `S_r = L_r^2 - 4 M_r^2 Delta` (sextic) | py L84, wl L112 | notes L259–262 | MATCH |
| `S_s = L_s^2 - 4 M_s^2 Delta` (sextic) | py L85, wl L113 | notes L266–269 | MATCH |
| `deg S_r = deg S_s = 6` | py L103–104/108–111, wl L119–122; sympy.txt L209–210 | notes L280 ("**sextic**"); appendix L981 | MATCH |
| Bezout bound `4*6 = 24` | py L167–170, wl L126–128; sympy.txt L227; wl.txt L46 | notes L303 (`4\cdot 6 = 24`); appendix L984 (`4\cdot 6=24`) | MATCH |
| `Delta_iso = (k_i + r k_j + s k_k)^2 - 2 H0 u (1+r^2+s^2)` | py L130, wl L141; sympy.txt L215 | notes L373–375 | MATCH |
| `tau_iso = 2H0/(k_rs + sqrt(k_rs^2 - 2 H0 u))` | py L134, wl L144; sympy.txt L216 | notes L379–383 | MATCH |
| equal-mix stationarity `N_r(1,1)=N_s(1,1)=0` | py L159–160, wl L166–167; sympy.txt L221–222 | notes L410–414 | MATCH |
| `tau = 2 H0 / Phi`, `Phi = (k_i+r k_j+s k_k+sqrt(Delta))/sqrt(1+r^2+s^2)` | py L50–51, wl L66–67; sympy.txt L9–24 | notes L113–122; appendix L974–979 | MATCH |

INTERNAL (scaffolding, no finding): `expected_dPhi_dr`/`expected_dPhi_ds` (driver expressions for the derivative-law assert), `den`/`linearK`/`stationaryDenominator` (intermediate), `cross-elimination identity`/`square elimination identity` residuals (assert residuals, all 0), `directNumerator[r/s]` printed forms (intermediate of the same assert), `k_rs` (intermediate).

reconciliation: complete; 14 deliverable values checked, 0 misaligned.

Note: every emitted *result value* reconciles with the card/notes/appendix. The two findings (F1 transliteration, F2 stale Verification line) are not value mismatches — F2 is a coverage/freshness misalignment in the card's Verification prose, not a numeric/symbolic value disagreement.

## Self-test notes

Variable-independence trap: checked the M1 derivative-law asserts — `Phi` genuinely depends on `r` and `s` (via `B r`, `D r^2`, `E r s`, `k_j r`, and `sqrt`), so `diff(Phi,r)` is not identically zero; the assert subtracts a separately-built `Nr`, so it is non-tautological and a wrong formula would leave a nonzero residual. Trivial-case check: under the symmetric substitution with `r=s=1`, `M_r(1,1)=k(1+1)-1(k+k)=0` analytically, so `N_r(1,1)=L_r(1,1)` and the symmetric `L_r` vanishes — the `=0` result is real, not forced. Degree/Bezout literals (4,6,24) are computed from expanded polynomials, not posited, so the `==4`/`==6`/`==24` asserts are substantive. Path check: F1 names the `.wl` under `mathematica/` and F2 names the card under `paper/stages/` — both correct directories; all paths are rooted at `/var/projects/toy_physics`.
