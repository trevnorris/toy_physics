---
unit_id: 052
batch: III.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage052_nontwin_asymmetry_threshold.md
  paper_appendix: present
---

# Audit unit 052 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_052.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage052_nontwin_asymmetry_threshold.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 82 inspected)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage052_nontwin_asymmetry_threshold_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage052_nontwin_asymmetry_threshold_mathematica_audit.txt`

## What the paper claims

Stage 052 asserts an exact three-regime classification of the lowest-lane support problem in terms of the branch product `Pi_tr` and the mixed product scale `C_mix`. The required support ratio is the closed form `zeta_req = (Pi_tr - C_mix)/[C_mix - eps(2 C_mix - Pi_tr)]` (eq:app-stage052-zeta-req), strictly monotone in `Pi_tr` with derivative `C_mix(1-eps)/[...]^2 > 0` (eq:app-stage052-zeta-derivative). The exact regime split (eq:app-stage052-regimes) is: `Pi_tr <= C_mix` → mixed-only enough; `C_mix < Pi_tr <= 2 C_mix` → symmetric lowest twin enough; `Pi_tr > 2 C_mix` → non-twin asymmetry required. The excess beyond the symmetric twin is `Delta_zeta = zeta_req - 1 = (1-eps)(Pi_tr - 2 C_mix)/[C_mix - eps(2 C_mix - Pi_tr)]` (eq:app-stage052-Dzeta), and the general lowest-lane resource factorization is `zeta_0^phys = (K_W^eff / K_phi,0^eff) Omega_0^2` (eq:app-stage052-zeta0-general). `\stagefield{Output}` enumerates: "Regime split (eq:app-stage052-regimes), excess formula (eq:app-stage052-Dzeta), and resource factorization (eq:app-stage052-zeta0-general)." Notes additionally enumerate diagnostics (`zeta_0^twin = 1`, pure-overlap rescue `Omega_0 >= sqrt(zeta_req)`, pure-softening rescue `K_phi,0 <= K_W/zeta_req`, softening fraction `(1-eps)(Pi_tr-2C_mix)/(Pi_tr-C_mix)`).

## What the script claims to verify

Both scripts assert exact symbolic identities for the non-twin asymmetry threshold: (a) the closed form `zeta_req` and its boundary values at `Pi_tr = C_mix` (zero) and `Pi_tr = 2 C_mix` (unity); (b) the closed form for `d zeta_req/d Pi_tr`, computed by SymPy/Mathematica differentiation and compared to the paper's expected closed form (Mathematica also adds an independent logarithmic-differentiation cross-check); (c) the excess `Delta_zeta = zeta_req - 1` against the paper's expected form (Mathematica uses `Together[zetaReq - 1]` as an independent derivation path); (d) the rescue-threshold inversions of `zeta_0^phys = zeta_req` for `Omega_0^2` and `K_phi,0` via `sp.solve`/`Solve`, compared against the closed-form `Omega_0_req^2 = zeta_req K_phi,0/K_W` and `K_phi,0_req = K_W Omega_0^2/zeta_req`; (e) the symmetric twin diagnostic `zeta_0^twin = 1`; and (f) the softening fraction `1 - K_phi,0_req/K_W` at equal overlap against the paper's expected `(1-eps)(Pi_tr - 2 C_mix)/(Pi_tr - C_mix)`, with Mathematica adding `Together[1 - 1/zetaReq]` as a second independent derivation path.

## Paper ↔ script cross-check

| Paper deliverable | Eq | Script-side check | Status |
|---|---|---|---|
| Closed form `zeta_req` | eq:app-stage052-zeta-req | boundary values at `Pi_tr = C_mix` (=0) and `Pi_tr = 2C_mix` (=1), construction from `Sreq` definition | match (boundary values + monotonicity uniquely identify the closed form on the support-needed branch) |
| Monotonicity `d zeta_req/d Pi_tr > 0` | eq:app-stage052-zeta-derivative | derivative computed and compared to the paper closed form (sympy line 61; mathematica line 62) | match |
| Three-regime split | eq:app-stage052-regimes | implicit via boundary values + strict monotonicity (no explicit per-regime sign assertion) | partial-match (implicit but mathematically sufficient) |
| Excess `Delta_zeta` | eq:app-stage052-Dzeta | `Delta_zeta - expected` zero check (sympy line 69; mathematica line 68) | match |
| Resource factorization `zeta_0^phys = (K_W/K_phi,0) Omega_0^2` | eq:app-stage052-zeta0-general | `zeta0_phys = KW * Omega0**2 / Kphi0` (sympy line 72; mathematica line 70) and `Solve` round-trip on `zeta_phys = zeta_req` | match (paper's `eff` superscripts treated as cosmetic labels by both engines, no physical operations distinguish them) |
| Symmetric twin diagnostic `zeta_0^twin = 1` (from notes §4) | notes:151-153 | `zeta_twin - 1 == 0` after substitution `Omega0=1, Kphi0=KW` (sympy line 106; mathematica line 94) | match |
| Softening-fraction diagnostic (notes §4) | notes:177-179 | `softFrac - softExpected` zero check (sympy line 122; mathematica line 98), Mathematica adds Together-based independent path | match |
| Pure-overlap and pure-softening thresholds (notes §3) | notes:124-133 | `Omega_req_equal_stiff = sqrt(zeta_req)`, `Kphi_req_equal_ov = KW/zeta_req` printed only, no `expect_zero` assertion | print-only (these are intermediates, not load-bearing claims; the load-bearing softening fraction is asserted) |

`paper_alignment: aligned`. The Output line of the paper card maps cleanly to script-side assertions; the regime split is mathematically implied by the verified boundary values plus the verified monotonicity, so an explicit per-regime sign check would be redundant.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 52 | `expect_zero("zeta_req at Pi_tr = C_mix", zeta_req.subs(Pi_tr, Cmix))` | zeta_req boundary at lower edge of support-needed branch | yes |
| A2 | sympy | 53 | `expect_zero("zeta_req at Pi_tr = 2 C_mix minus 1", zeta_req.subs(Pi_tr, 2*Cmix) - 1)` | zeta_req boundary at twin-failure edge | yes |
| A3 | sympy | 61 | `expect_zero("dzeta_req/dPi - expected", dz_dPi - dz_expected)` | monotonicity / derivative closed form | yes |
| A4 | sympy | 69 | `expect_zero("Delta_zeta - expected", Delta_zeta - Delta_expected)` | excess formula | yes |
| A5 | sympy | 90-93 | `expect_zero("solve(zeta_phys = zeta_req) for Omega0^2 - expected", sol[0] - Omega0_req_sq)` | rescue threshold via Omega_0^2 | yes (engine-derived solve vs closed form) |
| A6 | sympy | 99-102 | `expect_zero("solve(zeta_phys = zeta_req) for Kphi0 - expected", sol[0] - Kphi0_req)` | rescue threshold via K_phi,0 | yes (engine-derived solve vs closed form) |
| A7 | sympy | 106 | `expect_zero("zeta_0^(twin) - 1", zeta_twin - 1)` | twin diagnostic | yes |
| A8 | sympy | 122 | `expect_zero("softening fraction - expected", soft_frac - soft_expected)` | softening fraction closed form | yes |
| B1 | math  | 48 | `expectZero["zeta_req at Pi_tr = C_mix", zetaReq /. piTr -> cMix]` | zeta_req boundary at lower edge | yes |
| B2 | math  | 49 | `expectZero["zeta_req at Pi_tr = 2 C_mix minus 1", (zetaReq /. piTr -> 2 cMix) - 1]` | zeta_req boundary at twin-failure edge | yes |
| B3 | math  | 59 | `expectZero["dZdPi vs dZdPiAlt (independent path)", dZdPi - dZdPiAlt]` | derivative cross-engine independent path | yes |
| B4 | math  | 62 | `expectZero["dzeta_req/dPi - expected", dZdPi - dZExpected]` | derivative vs paper closed form | yes |
| B5 | math  | 68 | `expectZero["Delta_zeta - expected", deltaZetaDerived - deltaExpected]` | excess formula vs paper | yes |
| B6 | math  | 78-79 | `expectZero["solve(zeta_phys = zeta_req) for Omega0^2 - expected", omegaSqSol - omega0ReqSq]` | rescue threshold via Omega_0^2 | yes |
| B7 | math  | 82-83 | `expectZero["solve(zeta_phys = zeta_req) for Kphi0 - expected", kPhi0Sol - kPhi0Req]` | rescue threshold via K_phi,0 | yes |
| B8 | math  | 92 | `expectZero["softFrac vs Together[1 - 1/zetaReq] (independent path)", softFrac - softFracDerived]` | softening fraction independent path | yes |
| B9 | math  | 94 | `expectZero["zeta_0^(twin) - 1", zetaTwin - 1]` | twin diagnostic | yes |
| B10 | math | 98 | `expectZero["softening fraction - expected", softFrac - softExpected]` | softening fraction vs paper closed form | yes |

All rows trace to a paper-side deliverable or notes-side diagnostic; no orphan assertions.

## Findings

None. The script's claims align with the paper's `\stagefield{Output}` and the notes' diagnostic suite; every assertion is non-tautological under the v2 (post-v1-applied) state, the two engines verify the same identities by independent paths in three places (logarithmic differentiation of the derivative, Together-based delta, Together-based softening fraction), and both transcripts exit 0 with fresh mtimes.

## Independent-derivation check (Mathematica)

The `.wl` script (post-v1) introduces three engine-independent paths absent from the `.py`:

1. `zetaReq` is built via `Solve[(sReq - 1) - zetaSym (1 + eps (sReq - 2)) == 0, zetaSym]` (line 42) rather than direct substitution. The SymPy script types `(Sreq - 1)/(1 + eps(Sreq - 2))` literally (line 45).
2. `dZdPiAlt` is computed by logarithmic differentiation `zetaReq (D[numZ, piTr]/numZ - D[denZ, piTr]/denZ)` (lines 56-58) and cross-checked against `D[zetaReq, piTr]` at line 59. The SymPy script does no analogous independent-path check.
3. `softFracDerived` is computed as `Together[1 - 1/zetaReq]` (line 90) and cross-checked against the loop-driven `softFrac = 1 - Kphi_req_equal_ov/KW` at line 92.

The `dZExpected`, `deltaExpected`, and `softExpected` hand-typed closed forms remain in the `.wl`, but they are no longer the only Mathematica-side derivation: each is preceded by an engine-independent constructive path. The transliteration concern raised in v1 is resolved; the second engine now does carry information beyond echoing the first.

## Engine cross-check

Both engines produce identical canonical forms for `zeta_req` (`-(Cmix - Pi_tr)/(Cmix - 2 Cmix eps + eps Pi_tr)` / `-((cMix - piTr)/(cMix - 2*cMix*eps + eps*piTr))`), for `d zeta_req/d Pi_tr` (matching `Cmix(1-eps)/(Cmix - eps(2 Cmix - Pi_tr))^2` after sign-normalization), and for the softening fraction. Both transcripts print PASS on every `expectZero`/`expect_zero` and exit 0. No `engine_disagreement`.

`outputs_fresh: true`. sympy output (2026-05-22 17:35) is newer than the sympy script (2026-05-22 17:19); mathematica output (2026-05-22 17:24) is newer than the mathematica script (2026-05-22 17:22). No `stale_output`.

## Verdict justification

`clean`. The paper card's three-line `\stagefield{Output}` — regime split, excess `Delta_zeta`, and resource factorization — is exercised faithfully by the assertions in both engines: boundary values at `Pi_tr = C_mix` and `Pi_tr = 2 C_mix` plus strict monotonicity uniquely identify the closed-form `zeta_req` (and hence the three-regime split) on the support-needed branch; `Delta_zeta - expected = 0` is asserted directly; and the resource factorization is exercised via `Solve` round-trips on `zeta_0^phys = zeta_req` for both `Omega_0^2` and `K_phi,0`. The notes' additional diagnostics (twin lane `zeta_0^twin = 1`, softening fraction) are also asserted. Attacks I tried that failed: (i) probing for hidden sign errors at `Pi_tr < C_mix` where the formula technically yields `zeta_req < 0` — the script is explicit that the formula is only used on the support-needed branch and the boundary check at `Pi_tr = C_mix` pins the formula at the right value; (ii) attempting to construct a counterexample where `Omega0_req_sq` differs from the `sp.solve` output — solving the linear equation `K_W x/K_phi,0 = zeta_req` for `x` is one-step and the solver-vs-closed-form match holds identically (mild solver-sanity flavor, but the assertion is still non-tautological since `sp.solve` could in principle return a non-canonical form); (iii) checking for parity defects in the regime-classification chain — boundary values + monotonicity is mathematically airtight on the rational-function domain. The Mathematica transliteration concern from v1 is also resolved by the three independent-path checks at lines 59, 90-92, and 42-44. No `paper_misalignment`: every assertion traces to a paper or notes claim, and every paper Output item has a script-side counterpart.

## Self-test notes

- **Variable independence**: `zeta_req` is `(Pi_tr - C_mix)/(C_mix - eps(2 C_mix - Pi_tr))`, which has explicit `Pi_tr` dependence in both numerator and denominator, so `D[zeta_req, Pi_tr]` and the logarithmic-derivative alt path both produce nonzero residuals; `d zeta_req/d Pi_tr - expected` is a genuine zero-of-difference check, not a `0 == 0` trap.
- **Symmetry/parity**: No unbounded-domain integrals in this unit; trap not applicable.
- **Trivial-case pre-check**: With `Pi_tr=3, Cmix=1, eps=1/2`, `zeta_req = (3-1)/(1 - (1/2)(2-3)) = 2/(3/2) = 4/3`; `Delta_zeta = 1/3`; `softFrac = 1 - 1/(4/3) = 1/4`; `softExpected = (1/2)(3-2)/(3-1) = 1/4`. All match. At `Pi_tr = C_mix = 1`: `zeta_req = 0`. At `Pi_tr = 2 C_mix = 2`: `zeta_req = (2-1)/(1 - (1/2)(2-2)) = 1`. Boundary anchors confirmed by hand.
- **Path specifications**: No new scripts proposed; existing paths are correct (`.py` under `scripts/`, `.wl` under `mathematica/`).
- **Paper round-trip**: Re-checked the paper card `\stagefield{Output}` line and notes §5 ("What is exact now"). Every paper deliverable has a script-side check; no new misalignment introduced.
