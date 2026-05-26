---
unit_id: 022
batch: I.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md"]
  paper_appendix: present
---

# Audit unit 022 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_022.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 66: "022 & Grouped normalization bridge & ExactClosure / Open & Conversion from operator moments to normalized response moments and invariant outgoing-normalization product.")
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage022_grouped_p2_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.txt`

## What the paper claims

Per the `\stagefield{Output}` line of stage_022.tex:137: "Stage~022 outputs the normalized response formulas \eqref{eq:app-stage022-u2u4}, the grouped inverse map \eqref{eq:app-stage022-grouped-inverse}, the prefactor coefficients \eqref{eq:app-stage022-pref-coeffs}, the odd coefficient \eqref{eq:app-stage022-k-gamma}, and the invariant normalization test \eqref{eq:app-stage022-p0-target}." The five distinct deliverables are: (D1) u_2^{(A)} = -D_{A2}/D_{A0} and u_4^{(A)} = (D_{A2}^2 - D_{A0}D_{A4})/D_{A0}^2; (D2) grouped inverse map x_{20}=xbar+4a_x, x_{21}=xbar-a_x+b_x, x_{22}=xbar-a_x-b_x; (D3) P_{A0}, P_{A2}, P_{A4} coefficients of the internal prefactor D_{A0}N_A/(D_A^cons)^2; (D4) branch coefficients K_{A0}=P_{A0}, K_{A2}=P_{A2}+AP_{A0}, K_{A4}=P_{A4}+AP_{A2}+BP_{A0}, Gamma_5^{(A)}=G_5 P_{A0}; (D5) invariant normalization product mhat_0^2 P_0 = 54 G c_s^5/(5 a^5 c^5). The notes additionally include (D6, §6) the Stage-021 prototype N_0/N_2/N_4 in terms of (Delta0, S2, P0_proto, gW), and (D7, §8) constant-prefactor K_2/K_4 branch values. The compact-outgoing-fingerprint values A = a^2/(9c_s^2), B = 4a^4/(81c_s^4), G_5 = a^5/(27c_s^5) appear in the paper card eq:app-stage022-abg.

## What the script claims to verify

The SymPy script (Sections I-V) and Mathematica script (Sections I-VI) jointly assert: (a) the inverse series expansion of D0/Dcons gives u2, u4 with the boxed signs; (b) the prefactor series D0*N/Dcons^2 gives P0, P2, P4 with the boxed formulas; (c) multiplication by the outgoing fingerprint 1+A*omega^2+B*omega^4+iG5*omega^5 produces K0, K2, K4, Gamma5 with the claimed combinations; (d) the grouped inverse map (xbar, a_x, b_x) → (x20, x21, x22) recovers identity and collapses isotropically to a_2=b_2=0; (e) the Stage-021 one-port prototype's N0/N2/N4 reduce to the closed-form rational expressions in (Delta0, S2, P0_proto, gW); (f) the compact outgoing-fingerprint values A, B, G5 are reproduced from a spherical-Hankel-derived outgoing DtN expansion; (g) solving mhat^2*P0*Gamma5_port = 2G/(5c^5) for P0 and evaluating at mhat=1 yields P_0 = 54 G c_s^5 / (5 a^5 c^5), along with K2_target = 6 G c_s^3 / (5 a^3 c^5) and K4_target = 8 G c_s / (15 a c^5).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1 (u2, u4) | sympy:80-81; mathematica:47-48 | match |
| D2 (grouped inverse map) | sympy:181-183; mathematica:93-95 | match |
| D3 (P0, P2, P4) | sympy:118-123; mathematica:72-77 | match |
| D4 (K0, K2, K4, Gamma5) | sympy:133-136; mathematica:78-81 | match |
| D5 (mhat^2*P0 target) | sympy:298-300; mathematica:167 | match |
| D6 (notes §6: prototype N0/N2/N4) | sympy:217-222; mathematica:119-124 | match |
| D7 (notes §8: const-pref K2/K4) | sympy:314-320; mathematica:168-169 | match |
| paper eq:app-stage022-abg (A,B,G5 values) | sympy:275-277 (anchored via j2/y2 polynomial-rational form); mathematica:154-156 (anchored via SphericalHankelH1) | match |

Every paper-side deliverable has a corresponding script-side check, and every script-side assertion traces back to either an explicit paper equation or to the notes file. No orphan assertions.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 80 | `simplify(u2 + D2/D0) == 0` (u2 from series) | D1 | yes |
| A2 | sympy | 81 | `simplify(u4 - (D2^2-D0*D4)/D0^2) == 0` (u4 from series) | D1 | yes |
| A3 | sympy | 118 | `P0 - N0/D0 == 0` (P0 from series of D0*N/Dcons^2) | D3 | yes |
| A4 | sympy | 119 | `P2 - (D0*N2-2*D2*N0)/D0^2 == 0` | D3 | yes |
| A5 | sympy | 120-123 | `P4 - (D0^2*N4 - 2*D0*(D2*N2+D4*N0) + 3*D2^2*N0)/D0^3 == 0` | D3 | yes |
| A6 | sympy | 133 | `K0 - P0 == 0` | D4 | yes |
| A7 | sympy | 134 | `K2 - (P2+A*P0) == 0` | D4 | yes |
| A8 | sympy | 135 | `K4 - (P4+A*P2+B*P0) == 0` | D4 | yes |
| A9 | sympy | 136 | `Gamma5 - G5*P0 == 0` | D4 | yes |
| A10 | sympy | 181-183 | inverse map recovers x20, x21, x22 | D2 | yes (linear algebra check) |
| A11 | sympy | 193-195 | xbar_iso=xQ, ax_iso=0, bx_iso=0 | D2 | yes |
| A12 | sympy | 217-222 | prototype N0, N2, N4 match closed forms | D6 (notes §6) | yes |
| A13 | sympy | 275-277 | Stage-4 A, B, G5 from j2/y2 polynomial form match paper eq:app-stage022-abg | paper eq:app-stage022-abg | yes |
| A14 | sympy | 298-300 | `(NQ_target at mhat=1) - 54*G*c_s^5/(5*a^5*c^5) == 0` | D5 | yes |
| A15 | sympy | 313-320 | K2_target, K4_target at mhat=1 | D7 (notes §8) | yes |
| B1-B7 | mathematica | 47-48,72-81 | matches A1-A9 by inverse-relation Solve route (independent algorithm) | D1, D3, D4 | yes |
| B8-B13 | mathematica | 93-98 | inverse-map identity (parallel check; engines acknowledge no engine-independent route exists for 3x3 linear algebra) | D2 | yes |
| B14-B16 | mathematica | 119-124 | prototype N0, N2, N4 via inverse-relation Solve | D6 | yes |
| B17-B19 | mathematica | 154-156 | A, B, G5 from SphericalHankelH1[2,z] expansion (independent of SymPy's polynomial form) | paper eq:app-stage022-abg | yes |
| B20-B22 | mathematica | 167-169 | mhat=1 K0, K2, K4 targets | D5, D7 | yes |

All assertions are non-tautological in the sense that the LHS is computed independently of the RHS by series expansion or symbolic Solve, then the formula is checked against the paper-stated closed form. The one borderline case is the inverse-map identity (A10/B8-B10): the forward map (xbar, a_x, b_x) is a 3x3 linear transform of (x20, x21, x22), and any consistent inverse can be substituted back to recover identity. The Mathematica script explicitly flags this in a comment on line 90-92 ("Intentional parallel check: the 3x3 inverse-map identity admits no engine-independent route. Both engines verify the same algebra as a sanity cross-check."). I accept this disclosure; if the proposed inverse map formulas were wrong, the recovery would fail, so the check is substantive in the bounded sense.

## Findings

(No findings.)

## Independent-derivation check (Mathematica)

The Mathematica .wl is NOT a transliteration. Two specific structural differences confirm independent derivation:

1. **Section I/II inverse-relation route.** SymPy uses `sp.series(D0/Dcons, omega, 0, 6).removeO()` and `sp.series(D0*Nfac/Dcons**2, omega, 0, 6).removeO()` to extract u_n and P_n directly (forward series expansion). Mathematica instead posits a candidate `yRespCand = 1 + u2Sym*omega^2 + u4Sym*omega^4` (similarly `prefCand`), imposes `yRespCand*dCons - d0 = 0` coefficient-wise, and solves the resulting linear system via `Solve` (inverse-relation route). The two algorithms are mathematically equivalent but procedurally distinct.

   - SymPy `scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py:72,105`: `sp.series(D0 / Dcons, omega, 0, 6).removeO()` and `sp.series(D0 * Nfac / Dcons**2, omega, 0, 6).removeO()`.
   - Mathematica `mathematica/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.wl:36-40,53-57`: `prod = Expand[yRespCand*dCons - d0]; coeffEqs = Table[Coefficient[prod, omega, k] == 0, {k, 0, 4}]; sol = First[Solve[coeffEqs, {u2Sym, u4Sym}]]`.

2. **Section V outgoing fingerprint anchor.** SymPy constructs h_2^{(1)}(z) via the explicit polynomial-rational form `(3/z^3 - 1/z)*sin(z) - 3*cos(z)/z^2 + i*[-(3/z^3-1/z)*cos(z) - 3*sin(z)/z^2]`. Mathematica uses the built-in `SphericalHankelH1[2, z]`. The two representations are algebraically equivalent but pass through different simplification paths.

   - SymPy `scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py:262-264`: `j2 = (3/z**3 - 1/z) * sp.sin(z) - 3 * sp.cos(z) / z**2; y2 = -((3/z**3 - 1/z) * sp.cos(z) + 3 * sp.sin(z) / z**2); h2 = sp.simplify(j2 + I * y2)`.
   - Mathematica `mathematica/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.wl:143`: `h2 = SphericalHankelH1[2, z]`, with explicit comment "Use Mathematica's built-in spherical Hankel function instead of the explicit polynomial-rational form, so the derivation of A, B, G5 is independent of the SymPy script's choice of j2, y2 expressions."

The Section III inverse-map check is explicitly admitted as a parallel sanity check by both engines, with documented justification (linear algebra admits no engine-independent route). I accept this disclosure.

## Engine cross-check

Both engines produce identical final results, as verified in their saved outputs:

- u2 = -D2/D0, u4 = (D2^2 - D0*D4)/D0^2 (sympy:18-26; mathematica:9-10)
- P0 = N0/D0, P2 = (D0*N2 - 2*D2*N0)/D0^2, P4 = (D0^2*N4 - 2*D0*(D2*N2+D4*N0) + 3*D2^2*N0)/D0^3 (sympy:38-52; mathematica:19-24)
- Stage-4 prototype N0 = P0_proto^2/Delta0^2, N2 = 2*P0_proto*(P0_proto*S2 - Delta0*gW)/Delta0^3, etc. (sympy:114-130; mathematica:53-58)
- A = a^2/(9*c_s^2), B = 4*a^4/(81*c_s^4), G5 = a^5/(27*c_s^5) (sympy:155-172; mathematica:63-68)
- Required mhat^2 * P0 = 54*G*c_s^5 / (5*a^5*c^5*mhat^2) (sympy:188-193; mathematica:74)

All 27 SymPy `expect_zero` checks pass with residual 0; all 22 Mathematica `expectZero` checks emit PASS. No engine disagreement.

## Outputs freshness

Script mtimes: 2026-05-21 13:45 (both .py and .wl). Output mtimes: 2026-05-21 15:00 (sympy.txt) and 15:02 (mathematica.txt). Outputs are fresher than scripts. No `stale_output` concern.

## Verdict justification

Verdict: **clean**. I attacked the script along the following axes and failed to find a defect:

1. **Tautology probe.** Each `expect_zero` independently computes a series coefficient or solves an inverse relation, then compares to the paper-stated closed form. The closed forms are not built into the derivation; they could fail. The Section III inverse-map check borders on tautology (linear algebra round-trip), but both engines disclose this in code comments and the check still verifies that the proposed inverse formulas — not arbitrary ones — recover identity.
2. **Symbol-domain probe.** D0, D2, D4, Delta0 are nonzero (required to avoid division by zero in the series expansion of D0/Dcons and 1/Dcons^2). G, c, c_s, a, mhat, P0 are positive (required by the physical setup). No conflicting assumptions; no `simplify` masking branch issues.
3. **Hardcoded-constant probe.** The numbers 9, 27, 81, 54, 5, 2 that appear are not hardcoded — they emerge from the j2/y2 polynomial expansion (SymPy) and from the SphericalHankelH1 expansion (Mathematica). Both engines reproduce them independently. The 54 G c_s^5 / (5 a^5 c^5) target in eq:app-stage022-p0-target also derives from solving the invariant-product equation, not from a literal.
4. **Paper round-trip.** Every paper deliverable (D1-D5 from stage card Output line, plus D6-D7 from the notes) has a script-side check. The two "extra" script-side topics (prototype N0/N2/N4, spherical-Hankel anchor for A,B,G5) are both supported by the notes (§6) and by the paper card eq:app-stage022-abg respectively — they are not orphaned.
5. **Engine independence.** The two engines use procedurally distinct algorithms for both the series expansion (forward series vs inverse-relation Solve) and the spherical Hankel anchor (polynomial form vs SphericalHankelH1). This is a real second-engine check, not transliteration.

The only cosmetic anomaly I noted (not filed as a finding): the SymPy docstring (lines 3-23) and the Mathematica banner (line 26) carry stale "stage5"/"STAGE 005" labels from an earlier numbering. The actual assertions all target stage 022 content, the file names and paper card are correct, and no math is affected. This is a documentation-drift maintenance item, not an audit finding. It does NOT qualify as `paper_misalignment` because the script's bottom-line assertions (the highest-authority claim source) match the paper card exactly; only the lower-authority docstring/banner labels are stale.

## Self-test notes

I checked: (1) variable independence — every `sp.diff(h2, z)` operates on h2(z) which actually depends on z, so the derivative is nonzero (not a silent-pass trap); (2) series convergence — Lambda2 = omega * (dh2/dz)/h2 evaluated at z = omega*a/c_s has a finite static limit (h_2^{(1)} ~ -3i/z^3 has a pole, but the omega*dh/dz/h structure regularizes via cancellation, confirmed by output line 73 showing finite Gamma5_port = a^5/(27 c_s^5)); (3) inverse-map identity — substituting forward-map definitions into the proposed inverse formulas does recover x20, x21, x22 as required (independently verified algebraically: xbar+4a_x = (x20+2x21+2x22)/5 + 4*(2x20-x21-x22)/10 = x20 ✓); (4) K2/K4 target arithmetic — (54 G c_s^5/(5 a^5 c^5)) * a^2/(9 c_s^2) = 6 G c_s^3/(5 a^3 c^5) ✓, and (54 G c_s^5/(5 a^5 c^5)) * 4 a^4/(81 c_s^4) = 8 G c_s/(15 a c^5) ✓. No traps tripped.
