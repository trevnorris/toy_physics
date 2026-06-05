---
unit_id: 022
batch: I.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-04T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md]
  paper_appendix: present
---

# Audit unit 022 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_022.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage022_grouped_p2_normalization_bridge.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 66 + summary paragraph line 9)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage022_grouped_p2_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.txt`

## What the paper claims

Stage 022 is the grouped real \(P_2\) normalization bridge (checkpoint, Part I, anchors MTDC-T3/T4). It translates the Stages 004–021 projection-first Maxwell packet (one conservative grouped-lane operator \(D_A^{\rm cons}(\omega)\) plus the Stage-021 outgoing transfer factor \(N_A(\omega)\)) into the grouped real \(P_2\) response language. Per `\stagefield{Output}` (stage_022.tex:137), the stage outputs five deliverables: (1) the normalized-response coefficients \(u_2^{(A)}=-D_{A2}/D_{A0}\), \(u_4^{(A)}=(D_{A2}^2-D_{A0}D_{A4})/D_{A0}^2\) (eq. u2u4); (2) the grouped inverse map \(x_{20}=\bar x+4a_x\), \(x_{21}=\bar x-a_x+b_x\), \(x_{22}=\bar x-a_x-b_x\) (eq. grouped-inverse) with the isotropy gate \(a_2=b_2=a_4=b_4=0\); (3) the prefactor coefficients \(P_{A0}=N_{A0}/D_{A0}\), \(P_{A2}\), \(P_{A4}\) (eq. pref-coeffs); (4) the branch/odd coefficients \(K_{A0},K_{A2},K_{A4},\Gamma_5^{(A)}=G_5P_{A0}\) (eq. k-gamma), where only the static \(P_{A0}\) feeds the leading \(i\omega^5\) term; and (5) the invariant normalization test \(\widehat m_0^{\,2}P_0=54Gc_s^5/(5a^5c^5)\) (eq. p0-target), from \(\widehat m_0^{\,2}P_0\,a^5/(27c_s^5)=2G/(5c^5)\) against the GR target \(\gamma_{\rm GR}=2G/(5c^5)\). The notes additionally supply the Stage-021 one-lane prototype \(N_0/N_2/N_4\) closed forms (notes §6), the anisotropy norm \(A_2^2=4a_2^2+(4/5)b_2^2\) (notes:109), and the constant-prefactor even targets \(K_2=P_0a^2/(9c_s^2)\), \(K_4=4P_0a^4/(81c_s^4)\) (notes:287/289). The card flags the bridge formulas as exact-closure and the final universal equality as an open branch-realization test.

## What the script claims to verify

The SymPy script (5 sections) verifies, by independent series expansion: (I) \(u_2,u_4\) as coefficients of the series of \(D_0/D_{\rm cons}\), checked against the boxed forms; (II) \(P_0,P_2,P_4\) as the series of \(D_0 N/D_{\rm cons}^2\) and \(K_0,K_2,K_4,\Gamma_5\) from multiplying by the outgoing fingerprint \(1+A\omega^2+B\omega^4+iG_5\omega^5\); (III) the grouped inverse map round-trip (substitute \(\bar x,a_x,b_x\) definitions back, residual 0), the anisotropy norm, and the isotropic collapse; (IV) the Stage-021 prototype \(N_0/N_2/N_4\) from the rational \((P_0^{\rm proto}-g_W\omega^2)^2/(\Delta_0-S_2\omega^2+\omega^4)^2\), plus the dictionary back to \(\Omega_A,\Omega_W,R,g_A\); (V) the outgoing fingerprint coefficients \(A,B,G_5\) re-derived from the explicit spherical Bessel \(j_2,y_2\) Hankel ratio, then the invariant product \(\widehat m_0^2P_0=54Gc_s^5/(5a^5c^5)\) at \(\widehat m_0=1\) and the constant-prefactor \(K_2,K_4\) targets. The Mathematica script mirrors the deliverable list but derives each by an independent route (undetermined-coefficient `Solve` instead of `series`; built-in `SphericalHankelH1[2,z]` instead of explicit \(j_2/y_2\)). Every check is an `expect_zero`/`expectZero` of a residual.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) u2,u4 (eq. u2u4) | sympy I `expect_zero("u2 formula"…)`/`("u4 formula"…)` L80-81; math L47-48 | match |
| (2) grouped inverse map + isotropy gate (eq. grouped-inverse / isotropy-gate) | sympy III L181-195; math III L93-98 | match |
| (3) P0,P2,P4 (eq. pref-coeffs) | sympy II L118-123; math II L72-77 | match |
| (4) K0,K2,K4,Gamma5 (eq. k-gamma) | sympy II L133-136; math II L78-81 | match |
| (5) invariant normalization test (eq. p0-target / invariant-normalization) | sympy V L297-300; math VI L167 | match |
| notes §6 prototype N0/N2/N4 | sympy IV L217-222; math IV L119-124 | match |
| notes:109 anisotropy norm A_2^2=4a^2+(4/5)b^2 | sympy III L185 `A_sq` (printed, output L94) | match (printed, not asserted — see Assertion inventory note) |
| fingerprint A,B,G5 (eq. abg) | sympy V L275-277 (re-derived from j2/y2); math V L154-156 (re-derived from Hankel) | match |
| notes:287/289 constant-prefactor K2,K4 even targets | sympy V L313-320; math VI L168-169 | match |

Dominant pattern: every paper deliverable maps to a non-tautological script check. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 80 | `expect_zero(u2 - (-D2/D0))` from series of D0/Dcons | (1) u2 | yes |
| A2 | sympy | 81 | `expect_zero(u4 - (D2^2-D0 D4)/D0^2)` | (1) u4 | yes |
| A3 | sympy | 118-123 | `expect_zero(P0/P2/P4 - …)` from series of D0 N/Dcons^2 | (3) | yes |
| A4 | sympy | 133-136 | `expect_zero(K0/K2/K4/Gamma5 - …)` from Pref·fingerprint | (4) | yes |
| A5 | sympy | 181-183 | `expect_zero(x2j recovered)` inverse round-trip | (2) | yes |
| A6 | sympy | 193-195 | `expect_zero(xbar_iso-xQ, ax_iso, bx_iso)` | (2) gate | yes |
| A7 | sympy | 217-222 | `expect_zero(N0/N2/N4 - proto forms)` from series of Nproto | notes §6 | yes |
| A8 | sympy | 275-277 | `expect_zero(A/B/G5 - …)` re-derived from j2/y2 ratio | abg | yes |
| A9 | sympy | 297-300 | `expect_zero(NQ_target|mhat=1 - 54Gc_s^5/(5a^5c^5))` | (5) | yes |
| A10 | sympy | 313-320 | `expect_zero(K2/K4 target|mhat=1 - …)` | notes:287/289 | yes |
| M1 | math | 47-48 | `expectZero(u2/u4)` via undetermined-coeff Solve | (1) | yes |
| M2 | math | 72-81 | `expectZero(P0..Gamma5)` | (3),(4) | yes |
| M3 | math | 93-98 | `expectZero(inverse map + isotropy)` | (2) | yes |
| M4 | math | 119-124 | `expectZero(N0/N2/N4)` via Solve on nCand·dProto^2-pProto^2 | notes §6 | yes |
| M5 | math | 154-156 | `expectZero(A/B/G5)` re-derived from SphericalHankelH1 | abg | yes |
| M6 | math | 167-169 | `expectZero(mhat=1 K0/K2/K4 targets)` | (5),(287/289) | yes |

The anisotropy norm `A_sq` (sympy L185) is computed and printed (output L94) but not wrapped in an `expect_zero`; it is a display corollary of the inverse-map definitions, and its closed form matches notes:109 exactly (verified by hand below), so this is not a coverage gap — the load-bearing isotropy claim is exercised by A6/M3.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration; it derives each deliverable by a deliberately different route from the SymPy script. Three concrete examples:

1. **Series vs. undetermined coefficients.** SymPy gets \(u_2,u_4\) by `sp.series(D0/Dcons, omega, 0, 6).removeO()` then `.coeff` (py L72-74). Mathematica instead posits `yRespCand = 1 + u2Sym omega^2 + u4Sym omega^4`, forms `yRespCand*dCons - d0`, collects `Coefficient[...,omega,k]==0` for k=0..4, and `Solve`s (wl L36-42). Same for the prefactor (wl L53-57) and the prototype N's (wl L109-113). This is an algebraically independent inverse-relation route, not the SymPy series route.

2. **Fingerprint anchor — different special function.** SymPy builds the outgoing fingerprint from the explicit polynomial-rational spherical Bessels `j2 = (3/z^3 - 1/z) sin z - 3 cos z/z^2`, `y2 = -((3/z^3-1/z)cos z + 3 sin z/z^2)`, forms `h2=j2+I y2` (py L262-264). Mathematica uses the built-in `h2 = SphericalHankelH1[2, z]` (wl L143) and explicitly comments the intent (wl L140-142). Both then take the log-derivative \(\Lambda_2=\omega\,h_2'/h_2\), invert to \(Y_2\), normalize, and extract \(A,B,G_5\) — landing on the same \(a^2/9c_s^2\), \(4a^4/81c_s^4\), \(a^5/27c_s^5\). This is genuine cross-engine corroboration via independent function definitions.

3. **The one place that IS a deliberate parallel check is flagged as such.** The 3×3 grouped inverse-map identity (wl L90-92 comment) is acknowledged to "admit no engine-independent route"; both engines verify the same algebra as a sanity cross-check. This is honest and appropriate for a pure linear-algebra identity.

No `mathematica_transliteration` finding.

## Engine cross-check

Both engines assert the identical residual-zero set and both saved transcripts show all-zero / all-PASS:

| Quantity | sympy output | math output |
|---|---|---|
| u2,u4 | residual 0 (txt L15-16) | PASS (txt L12,14) |
| P0,P2,P4,K0,K2,K4,Gamma5 | residual 0 (txt L35-60) | PASS (txt L19-32) |
| inverse map + isotropy | residual 0 (txt L91-101) | PASS (txt L37-48) |
| N0,N2,N4 prototype | residual 0 (txt L111-113) | PASS (txt L53-58) |
| A,B,G5 fingerprint | residual 0 (txt L152-154) | PASS (txt L63-68) |
| mhat^2 P0 = 54Gc_s^5/(5a^5c^5) | residual 0 (txt L194) | PASS (txt L75-76) |
| K2,K4 const-prefactor targets | residual 0 (txt L217-218) | PASS (txt L77-80) |

Both transcripts print the same `Required (mhat^2) P0 = 54 G c_s^5/(5 a^5 c^5 mhat^2)` closed form (sympy txt L188-193; math txt L74). Engines agree.

## Verdict justification

Clean. I read the paper card, the full notes file, and the appendix row before opening the scripts, and built the five-deliverable model. Every deliverable is exercised by a non-tautological residual check in BOTH engines, and the two engines use deliberately independent derivation routes (series vs undetermined-coefficient Solve; explicit j2/y2 vs SphericalHankelH1). Attacks tried that failed: (a) tautology — the checks are `series → coeff → subtract claimed closed form`, where the closed form is independent of the series machinery, so they can genuinely fail; (b) hardcoded result — the load-bearing constants \(54/5\), \(2/5\), \(9/81/27\) are re-derived in §V from the GR target and the Hankel fingerprint, not pinned; (c) symbol-domain — `D0!=0`/`nonzero` is justified (it sits in denominators), positivity on \(G,c,c_s,a,\widehat m_0,P_0\) in §V is physically warranted and used only to let `solve` pick the single real branch; (d) value reconciliation — every emitted deliverable value matches the card/notes (see table below). The only blemish is a cosmetic stale-output mtime (the `.wl` was touched 2026-06-03 by a label-only banner fix `STAGE 005 → STAGE 022`, e2a4780, post-dating its 2026-05-21 transcript). The transcript still shows the old `STAGE 005` banner but is otherwise byte-faithful to the current math; a fresh run changes only the banner token. Per the freshness guard this does not warrant a finding (the math content does not disagree with what the current script produces). I tried to break the math and could not.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 14 deliverable values checked, 0 misaligned.

Every RESULT/deliverable value the scripts emit, and its location in the docs:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| u2 = -D2/D0 | py L80 / wl L47; out txt L9 | tex:31 (eq. u2u4); md:65 | MATCH |
| u4 = (D2^2 - D0 D4)/D0^2 | py L81 / wl L48; out txt L10 | tex:32-33; md:67 | MATCH |
| P0 = N0/D0 | py L118 / wl L72 | tex:80; md:133 | MATCH |
| P2 = (D0 N2 - 2 D2 N0)/D0^2 | py L119 / wl L73 | tex:81; md:135 | MATCH |
| P4 = (D0^2 N4 - 2 D0(D2 N2+D4 N0)+3 D2^2 N0)/D0^3 | py L120-123 / wl L74-77 | tex:82; md:137-138 | MATCH |
| K0=P0, K2=P2+A P0, K4=P4+A P2+B P0, Gamma5=G5 P0 | py L133-136 / wl L78-81 | tex:110-113 (eq. k-gamma); md:175-181 | MATCH |
| grouped inverse map x20=xbar+4ax, x21=xbar-ax+bx, x22=xbar-ax-bx | py L172-174 / wl L93-95 | tex:51-55 (eq. grouped-inverse); md:91-95 | MATCH |
| isotropy gate a2=b2=a4=b4=0 | py L193-195 / wl L96-98 | tex:60 (eq. isotropy-gate); md:115,245 | MATCH |
| anisotropy norm A^2 = 4 a^2 + (4/5) b^2 | py L185; out txt L94 | md:109 | MATCH |
| N0=P0proto^2/Delta0^2 | py L217 / wl L119 | md:209 | MATCH |
| N2=2 P0proto(P0proto S2 - Delta0 gW)/Delta0^3 | py L218 / wl L120 | md:211 | MATCH |
| N4=(Delta0^2 gW^2 - 2 Delta0 P0proto^2 - 4 Delta0 P0proto S2 gW + 3 P0proto^2 S2^2)/Delta0^4 | py L219-222 / wl L121-124 | md:213 | MATCH |
| fingerprint A=a^2/(9 c_s^2), B=4 a^4/(81 c_s^4), G5=a^5/(27 c_s^5) | py L275-277 / wl L154-156; out txt L155-172 | tex:97-100 (eq. abg); md:159-162 | MATCH |
| invariant product mhat^2 P0 = 54 G c_s^5/(5 a^5 c^5) | py L299 / wl L167; out txt L188-194 | tex:133 (eq. p0-target); md:270,281 | MATCH |

Notes on two emitted closed forms not given a stand-alone line in the docs (accounted, no finding):
- `K2_target|mhat=1 = 6 G c_s^3/(5 a^3 c^5)` (py L315) and `K4_target|mhat=1 = 8 G c_s/(15 a c^5)` (py L319) are NOT printed verbatim in tex/md as numbers, but they are the trivial products `NQ_target·A` and `NQ_target·B` of two values that ARE in the docs: the notes give the symbolic even targets `K_2=P_0 a^2/(9 c_s^2)` (md:287) and `K_4=4 P_0 a^4/(81 c_s^4)` (md:289) with `P_0=54 G c_s^5/(5 a^5 c^5)` (md:281). Substituting reproduces the script numbers exactly (hand-checked: 54/(5·9)=6/5, 54·4/(5·81)=8/15). These are corollary sanity substitutions, not stated deliverables → INTERNAL, no finding.
- `Required mhat^2 P0` with the explicit `mhat` denominator `54 G c_s^5/(5 a^5 c^5 mhat^2)` (out txt L188-193) is the general-mhat form of the boxed `mhat^2 P0 = 54 G c_s^5/(5 a^5 c^5)`; same content, MATCH above covers it.

INTERNAL-only items (scaffolding, no prose expected, no finding): `Y_resp` displayed series, `Pref`/`Ybranch` intermediate polynomials, the `Delta0/S2/P0_proto` back-dictionary to `Omega_A,Omega_W,R,g_A` (tex/md §6 carry these symbolically — md:202-204 — so even these reconcile), `Lambda2`/`Y2`/`Y2_hat` intermediate Hankel-ratio expressions, all `expect_zero`/`PASS` residual flags.

## Self-test notes

- **Variable independence:** the only derivative-bearing block is §V's `sp.diff(h2, z)` / `D[h2, z]` (py L265, wl L144); `h2(z)` genuinely depends on `z`, so the log-derivative is non-trivially nonzero — no identically-zero-derivative trap. No `sp.diff(EXPR, VAR)` where VAR is absent from EXPR anywhere.
- **Trivial-case / non-tautology:** each `expect_zero` subtracts an independently-stated closed form from a series/Solve-derived coefficient; substituting a concrete operator (e.g. D2=1,D0=1,D4=0) gives u2=-1 matching -D2/D0, and a wrong claimed sign would leave a nonzero residual — the checks can fail.
- **Constant round-trip:** re-derived 54/5, 9/81/27, 6/5, 8/15 by hand from γ_GR and the Hankel fingerprint; all consistent. No hardcoded-answer pin.
- **Freshness:** stale `.wl` output is a label-only banner mismatch (STAGE 005→022) with no math delta; intentionally not raised as a finding per the augmentation freshness guard.
