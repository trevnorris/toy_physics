---
unit_id: 095
batch: IV.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage095_second_order_geometry_contamination.md]
  paper_appendix: present
---

# Audit unit 095 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_095.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage095_second_order_geometry_contamination.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows: Prop. "Second-order onset of geometry contamination" L185-208; eps_2/eps_4 defs L130-152; isotropic decoupling L155-183; `\input{stages/stage_095}` L1224)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage095_second_order_geometry_contamination_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage095_second_order_geometry_contamination_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage095_second_order_geometry_contamination_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage095_second_order_geometry_contamination_mathematica_audit.txt`

## What the paper claims

Stage 095 is a geometry-lane firewall ledger step. The bottom-line claim, quoted from the card body block: "Weak \(l=0\leftrightarrow l=2\) mixing contaminates the grouped module only at \(O(\chi^2)\)." The notes derive this in a minimal two-mode reduced model: one grouped-`P2` quadrupole carrier `q` with conservative kernel `D_q(omega)=K_stat+K_pole/(1-omega^2/Omega_Q^2)`, one scalar/geometry mode `g` with `D_g(omega)=G_0+G_2 omega^2+G_4 omega^4+O(omega^6)`, coupled by a bilinear `chi M_0 q g`. Integrating out `g` yields the exact Schur complement `D_eff = D_q - chi^2 M_0^2/D_g`, whose low-frequency expansion gives the dynamic contamination coefficients `K_(g,2)^eff = chi^2 M_0^2 G_2/G_0^2` and `K_(g,4)^eff = chi^2 M_0^2 (G_0 G_4 - G_2^2)/G_0^3`. The obstruction numbers `eps_2 = Omega_Q^2 K_(g,2)^eff/K_pole`, `eps_4 = Omega_Q^4 K_(g,4)^eff/K_pole` are both `O(chi^2)`, so the pole fraction `c_pole=(1+eps_4)/[4(1+eps_2)^2]` deviates from `1/4` only at `O(chi^2)`; the static limit (`eps_2=eps_4=0`) returns `c_pole=1/4`. The deliverables are: (D1) the three closed-form `K_(g,k)^eff` coefficients; (D2) the eps_2, eps_4 definitions; (D3) the static-limit `c_pole=1/4`; (D4) the perturbative-stability statement that the first nonzero deviation is `O(chi^2)`.

## What the script claims to verify

The SymPy script symbolically series-expands `-chi^2 M0^2/Dg` to `O(w^4)`, extracts the `w^0/w^2/w^4` coefficients, and asserts each equals the notes' closed forms (`K0corr=-chi^2 M0^2/G0`, `K2corr=chi^2 M0^2 G2/G0^2`, `K4corr=chi^2 M0^2 (G0 G4-G2^2)/G0^3`). It then asserts each carries a `chi^2` factor, builds `eps2/eps4`, forms `c_pole=(1+eps4)/(4(1+eps2)^2)`, and verifies its chi-series: constant term `1/4` (`delta.subs(chi,0)==0`) and vanishing linear term (`diff(c_pole,chi)|0==0`), thus the leading deviation is `O(chi^2)`. The Mathematica script additionally and independently DERIVES the Schur complement: it writes the bilinear action `L=(1/2)q^2 dQ+(1/2)g^2 dG+chi m0 q g`, solves `D[L,g]==0` for `g*`, substitutes back, reads off `2*Coeff[LqEff,q,2]-dQ`, and asserts it equals `-chi^2 m0^2/dG` — then proceeds through the same coefficient and c_pole checks via `expectZero`.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1: K_(g,2)^eff, K_(g,4)^eff closed forms (+ static K_(g,0)) | sympy L39-41 asserts series coeffs == closed forms; wl L57-59 `expectZero` same; wl L43-44 independently derives the Schur form | match |
| D2: eps_2, eps_4 = Omega_Q^{2,4} K_eff/K_pole | sympy L51-52 builds & prints; wl L61-62 builds & prints | match |
| D3: static-limit c_pole = 1/4 | sympy L63 `delta.subs(chi,0)==0` (chi=0 ⇒ eps=0 ⇒ c_pole=1/4); wl L68/L72 | match |
| D4: first nonzero deviation is O(chi^2) | sympy L70 `chi1==0` + L57-59 c_pole series leading term chi^2; wl L72 `expectZero[d c_pole/dchi|0]` + L67 series | match |

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 39 | `simplify(K0corr - (-M0^2 chi^2/G0)) == 0` | D1 (static renorm) | yes |
| A2 | sympy | 40 | `simplify(K2corr - G2 M0^2 chi^2/G0^2) == 0` | D1 | yes |
| A3 | sympy | 41 | `simplify(K4corr - M0^2 chi^2 (G0 G4-G2^2)/G0^3) == 0` | D1 | yes |
| A4 | sympy | 46-48 | `factor(Kkcorr).has(chi**2)` | D4 (each coeff ∝ chi^2) | yes |
| A5 | sympy | 63 | `simplify(delta.subs(chi,0)) == 0` | D3 (static c_pole=1/4) | yes |
| A6 | sympy | 70 | `chi1 == 0` (linear chi term vanishes) | D4 | yes |
| A7 | mathematica | 43-44 | `expectZero[corrDerived - (-chi^2 m0^2/dGsym)]` | D1 (independent Schur derivation) | yes |
| A8 | mathematica | 57 | `expectZero[k0Corr/chi^2 + m0^2/g0]` | D1 | yes |
| A9 | mathematica | 58 | `expectZero[k2Corr/chi^2 - g2 m0^2/g0^2]` | D1 | yes |
| A10 | mathematica | 59 | `expectZero[k4Corr/chi^2 - m0^2(g0 g4-g2^2)/g0^3]` | D1 | yes |
| A11 | mathematica | 72 | `expectZero[D[cPole,chi]/.chi->0]` | D4 | yes |

## Findings

None. No script-side, paper-alignment, or reconciliation findings.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration of the `.py`. It is more independent than the SymPy script: the SymPy script takes the Schur form `-chi^2 M0^2/Dg` as a given premise from the notes (L25), whereas the Mathematica script (L37-44) re-derives that Schur complement from first principles. It constructs the bilinear Lagrangian `Lq = (1/2)*qSym^2*dQsym + (1/2)*gSym^2*dGsym + chi*m0*qSym*gSym`, eliminates `g` via `gStar = First[gSym /. Solve[D[Lq, gSym] == 0, gSym]]`, substitutes back, and identifies `dEffCoeff = 2*Coefficient[LqEff, qSym, 2]`, then asserts `corrDerived - (-chi^2*m0^2/dGsym) == 0`. The SymPy script has no analogue of this elimination step. The downstream coefficient checks are the same identities but reached via genuinely distinct machinery (`expectZero[FullSimplify[Together[Expand[...]]]]` vs `simplify(... ) == 0`), and the two engines' printed forms differ in term ordering (e.g. wl line 7 vs py line 1) while being algebraically equal — consistent with independent simplification, not echoing.

## Engine cross-check

Both engines agree on every emitted result:

| Quantity | SymPy output | Mathematica output | Agree |
|---|---|---|---|
| K0corr | `-M0**2*chi**2/G0` | `-((chi^2*m0^2)/g0)` | yes |
| K2corr | `G2*M0**2*chi**2/G0**2` | `(chi^2*g2*m0^2)/g0^2` | yes |
| K4corr | `M0**2*chi**2*(G0*G4-G2**2)/G0**3` | `(chi^2*(-g2^2+g0*g4)*m0^2)/g0^3` | yes |
| eps2 | `G2*M0**2*OQ**2*chi**2/(G0**2*Kpole)` | `(chi^2*g2*m0^2*oQ^2)/(g0^2*kPole)` | yes |
| eps4 | `M0**2*OQ**4*chi**2*(G0*G4-G2**2)/(G0**3*Kpole)` | `(chi^2*(-g2^2+g0*g4)*m0^2*oQ^4)/(g0^3*kPole)` | yes |
| c_pole leading dev | `1/4 - ... chi^2 ...` | `1/4 - ... chi^2 ...` | yes |
| d c_pole/dchi\|0 | `0` | `0` (PASS) | yes |
| d^2 c_pole/dchi^2\|0 | `M0^2*OQ^2*(-G0*G2+OQ^2*(G0*G4-G2^2)/2)/(G0^3*Kpole)` | `-1/2*(m0^2*(2*g0*g2*oQ^2+(g2^2-g0*g4)*oQ^4))/(g0^3*kPole)` | yes (sign-/factor-equal) |

The two `d^2 c_pole/dchi^2|0` forms are identical after expanding: SymPy's `(-G0*G2+OQ^2(G0*G4-G2^2)/2)/(G0^3 Kpole)` distributes to `-1/2(2 G0 G2 + OQ^2(G2^2-G0 G4))/(G0^3 Kpole) · OQ^2... ` matching Mathematica's grouped form. Outputs are fresh (sympy .txt 14:28 / .wl .txt 14:29, both newer than scripts at 11:14).

## Verdict justification

`clean`. The scripts faithfully verify every stated deliverable of stage 095. Attacks tried and failed: (1) **Tautology** — the K-coefficient asserts compare a series-expansion result (`series(-chi^2 M0^2/Dg)`) against independently-stated notes closed forms; they could disagree if the expansion were wrong, so they are non-tautological. (2) **Truncation error** — both engines expand to `w^4`/`w^5`; the K4 coefficient correctly captures the cross-term `G2^2/G0^3` from the geometric series, matching the notes; truncation order is sufficient since only `O(w^4)` is claimed. (3) **Trivial O(chi^2) check** — `chi1==0` alone would be weak, but it is backed by the printed c_pole series showing the leading correction is explicitly chi^2 and by each K_eff carrying a `chi^2` factor (A4), so D4 is genuinely exercised, not just asserted via a vanishing first derivative of an even function. (4) **Symbol domains** — `G0/Omega_Q/Kpole` declared `nonzero` (sympy L21-22; wl L31), which is exactly the physical setup (denominators must not vanish); no over-strong positivity assumption hides a branch. (5) **Engine independence** — confirmed the Mathematica script re-derives the Schur complement rather than transliterating. I read the paper card, the single stage-095 notes file, and the part04 appendix rows (Prop. + eps defs + isotropic decoupling), and the scripts' verified claims match all of them.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 6 deliverable values checked, 0 misaligned

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| K_(g,0)^eff = -chi^2 M0^2/G0 (static renorm) | py out L2; wl out L8 | notes L65-66 (`- chi^2 M_0^2/G_0`); appendix (absorbed into contact slot) | MATCH |
| K_(g,2)^eff = chi^2 M0^2 G2/G0^2 | py out L3; wl out L9 | notes L74; appendix L203 | MATCH |
| K_(g,4)^eff = chi^2 M0^2 (G0 G4-G2^2)/G0^3 | py out L4; wl out L10 | notes L76; appendix L205 | MATCH |
| eps_2 = OQ^2 K2corr/Kpole | py out L8; wl out L17 | notes L80; appendix L133 | MATCH |
| eps_4 = OQ^4 K4corr/Kpole | py out L9; wl out L18 | notes L82; appendix L135 | MATCH |
| c_pole = 1/4 + (deviation only at O(chi^2)); static limit = 1/4 | py out L10,L12; wl out L19,L21-22 | notes L98-104 (`c_pole=(1+eps_4)/[4(1+eps_2)^2]`, dev O(chi^2)); card Checks L22 (`c_pole=1/4`); appendix L144,L151 | MATCH |

INTERNAL (scaffolding, no finding expected in prose): `corr`/`corrDerived` (intermediate Schur expression), `delta c_pole`/`delta` (residual driving the static check), `d c_pole/dchi|0 = 0` and `d^2 c_pole/dchi^2|0` (derivative diagnostics used to assert the O(chi^2) onset), PASS flags, `gStar`/`dEffCoeff` (Mathematica elimination intermediates).

All six stage deliverables emitted by the scripts reconcile with the notes and/or appendix and the card; no MISMATCH or MISSING-DELIVERABLE.

## Self-test notes

Checked: (1) Variable independence — `diff(c_pole, chi)` is taken w.r.t. chi, on which c_pole genuinely depends (via eps2/eps4 ∝ chi^2), so the derivative is not identically zero by construction; its vanishing AT chi=0 is the real content. (2) Symmetry — c_pole is even in chi (eps ∝ chi^2), so the linear term legitimately vanishes; the script does not over-claim by stopping there, since the printed chi^2 leading term and the per-coefficient chi^2 factor checks establish the onset order. (3) Trivial-case — at chi=0, eps2=eps4=0 ⇒ c_pole=(1+0)/(4·1)=1/4, so `delta.subs(chi,0)==0` correctly reduces to 0; the K4 cross-term `G2^2/G0^3` is exactly what the geometric series produces, confirmed in both outputs. No directive written (zero findings).
