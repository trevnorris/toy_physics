---
title: Paper 4D 2Pn
type: source
status: current
sources:
- research/4d_2pn/paper/4d_2pn.tex
last_updated: '2026-08-31'
---

## Purpose and scope

### source-paper-4d-2pn--purpose-and-scope — Conservative two-body 2PN derivation

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper develops the complete near-zone conservative two-body ledger through order \(c^{-4}\) for a unified 4D toy model. Its target is a reduced two-body Lagrangian whose perturbative Legendre transform equals the standard generic-frame ADM Hamiltonian through 2PN, conditional on declared projection, worldtube, low-frequency response, and self/static closures.

Sources:

- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:intro}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:intro-claims}`

## Source-unit map

- Entry point: `research/4d_2pn/paper/4d_2pn.tex`
- Role: `primary`
- Read mode: `semantic`
- Unit shape: `file`

## Key statements

### source-paper-4d-2pn--one-body-denominator — Minimal one-body response closure

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

Within the paper’s minimal low-frequency DtN-invariant ansatz, preserving the frozen 1PN coefficient removes the term linear in \(G-1\), while matching the isotropic Schwarzschild test-mass target fixes
\[
D_{\mathrm{eff}}(u)=C_s(u)\left[1+\frac{32768}{3249}(G(u)-1)^2\right]
=1-4u+8u^2+O(u^3).
\]
The ansatz is a protocol closure, not a theorem of the fully dynamical moving-throat PDE.

Sources:

- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:one-body-ansatz}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:onebody-final-closure}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:onebody-final-target}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:status-protocol}`

### source-paper-4d-2pn--self-static-seed — Repaired self/static 2PN seed

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

Minimal Bernoulli exactification first produces raw one-body coefficients \(1/16\), \(7/8\), and \(6\) for \(v^6\), \(Uv^4\), and \(U^2v^2\). The DtN denominator repair removes the residual \(6-2=4\), and quadratic local mass scaling with the pure-static isotropic gate fixes \(\lambda_\rho=\tfrac12\). Within this combined closure hierarchy, the assembly seed is
\[
\begin{aligned}
L_{2,\mathrm{self}+\mathrm{static}}
={}&\frac{m_Av_A^6+m_Bv_B^6}{16}
+\frac{7Gm_Am_B}{8r}(v_A^4+v_B^4)\\
&+\frac{2G^2m_Am_B}{r^2}(m_Bv_A^2+m_Av_B^2)
+\frac{G^3m_Am_B(m_A^2+m_B^2)}{4r^3}.
\end{aligned}
\]

Sources:

- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:self-static-self-expansion}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:self-static-u2v2-residual}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:self-static-lambda-rho}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:self-static-final-two-body}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{clm:self-static-seed}`

### source-paper-4d-2pn--comparable-mass-solve — Unique comparable-mass residual solve

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

After freezing the lower-order ledger and repaired seed, the perturbative Legendre identity makes the remaining 2PN match linear in the new cross-sector coefficients. In the paper’s \(Q_i\), \(T_i\), and \(S_1\) invariant basis, the unique solution is
\[
(q_1,\ldots,q_7)=\left(-\frac74,-\frac14,\frac{11}{8},\frac14,-\frac58,\frac32,\frac38\right),
\]
\[
(t_1,\ldots,t_6)=\left(0,\frac{11}{8},-\frac{15}{4},0,0,\frac{15}{8}\right),
\qquad s_1=\frac54.
\]

Sources:

- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:adm-H1H2-formula}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:adm-new-block-linear}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:adm-Q-basis}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:adm-T-basis}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:adm-S-basis}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{clm:adm-solve}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:adm-q-solution}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:adm-t-solution}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:adm-s-solution}`

### source-paper-4d-2pn--charge-notation-firewall — Residual coefficients are not electric charge

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The residual-basis \(q_i\) are ADM-lift coefficients, not electric charges. The historical gravity-side scalar coefficient formerly written \(q\) is renamed \(\kappa_\rho\) and fixed to \(1\); electric charge instead follows the distinct ontology \(\eta_Q=\pm1\), \(q_*=\eta_Qe_*\), and \(q_{\mathrm{eff}}=q_*/\sqrt{Z_{\mathrm{int}}}\). Quantized circulation belongs to the magnetic/vortical sector.

Sources:

- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:intro-motivation}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:carry-forward}`

### source-paper-4d-2pn--full-adm-match — Closure-conditioned generic-frame ADM match

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

After imposing the paper’s one-body DtN denominator, Bernoulli, and self/static closure choices, the complete coefficient of \(c^{-4}\), with \(v_{AB}=\mathbf v_A\!\cdot\!\mathbf v_B\), \(v_{An}=\mathbf v_A\!\cdot\!\mathbf n\), and \(v_{Bn}=\mathbf v_B\!\cdot\!\mathbf n\), is
\[
\begin{aligned}
L_{2,\mathrm{full}}={}&\frac{m_Av_A^6+m_Bv_B^6}{16}
+\frac{Gm_Am_B}{r}\Big[\frac78(v_A^4+v_B^4)-\frac74v_{AB}(v_A^2+v_B^2)\\
&-\frac14v_{An}v_{Bn}(v_A^2+v_B^2)+\frac{11}{8}v_A^2v_B^2+\frac14v_{AB}^2
-\frac58(v_A^2v_{Bn}^2+v_B^2v_{An}^2)\\
&+\frac32v_{AB}v_{An}v_{Bn}+\frac38v_{An}^2v_{Bn}^2\Big]\\
&+\frac{G^2m_Am_B}{r^2}\Big[\left(2m_B+\frac{11}{8}m_A\right)v_A^2
+\left(2m_A+\frac{11}{8}m_B\right)v_B^2\\
&-\frac{15}{4}(m_A+m_B)v_{AB}
+\frac{15}{8}(m_Av_{An}^2+m_Bv_{Bn}^2)\Big]\\
&+\frac{G^3m_Am_B}{4r^3}(m_A^2+5m_Am_B+m_B^2).
\end{aligned}
\]
Conditional on those inherited controlled and protocol inputs, the paper reports that the subsequent perturbative Legendre-transform algebra gives \(H_2=H^{\mathrm{ADM}}_{2\mathrm{PN}}\) and leaves no free dimensionless coefficient in this conservative two-body 2PN ledger. The equality is exact at the algebraic matching stage, but it does not independently derive the physical closure inputs.

Sources:

- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:one-body-ansatz}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:self-static}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:full-assembly-L2full}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:full-assembly-adm}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{clm:full-2pn}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:full-assembly-exact-match}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:full-assembly-proof-chain}`

## Computed evidence represented by the source

### source-paper-4d-2pn--computed-evidence-availability — Prepared computed evidence

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

not supplied in the prepared unit

Sources:

- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:verification}`

## Assumptions, exclusions, and open questions

### source-paper-4d-2pn--scope-boundaries — Controlled reduction and finite-size regime

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The result is a conservative near-zone small-body/worldtube reduction, not an assumption-free particle-limit theorem from a solved moving-throat bulk background. It excludes radiation reaction, dissipative leakage, spin, non-adiabatic response, arbitrary finite-size structure, and strong-field completion. The source estimates
\[
\frac{F_{\mathrm{finite\text{-}size}}^{(2)}}{F_{\mathrm{universal}}^{(2)}}
\sim\left(\frac{ac^2}{GM}\right)^2,
\]
so its universal monopole 2PN interpretation requires \(a\ll GM/c^2\), not merely \(a\ll r\).

Sources:

- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:intro-nonclaims}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:self-static-finite-size-gate}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{app:static-finite-size}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:conclusion-open}`

### source-paper-4d-2pn--fixed-throat-response — Appendix-side zero-frequency packaging

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The appendix-side constructive packaging assigns the zero-frequency data
\[
R_{1\perp}=\frac72,\quad R_{10}=4,\quad
R_0=R_{20}=R_{21}=R_{22}=1,
\]
\[
J=\left(\frac4{\sqrt5},\frac54,0,0,0,0\right),\qquad
\Delta_{\mathrm{geom}}=\frac{281}{80}.
\]
Here \(J\) is a throat-response source-weight vector, not a Maxwell current. The paper explicitly treats this packaging as an interpretation and continuation of the conservative result, not as load-bearing evidence for the main 2PN theorem.

Sources:

- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:status-throat-fixed}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:status-fixed-odd-residues}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:status-fixed-even-residues}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:status-fixed-J}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:status-fixed-deltageom}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:intro-nonclaims}`

### source-paper-4d-2pn--dynamic-throat-observables — Dynamic DtN completion remains open

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

The dynamic pole scales
\[
\Omega_{1\perp},\Omega_{10},\Omega_0,\Omega_{20},\Omega_{21},\Omega_{22},\Omega_g
\]
remain unresolved. Absolute medium normalization, localization profiles, finite-size moments, damping, phase lag, multimode response, and a fully dynamical moving-throat PDE completion also remain outside the completed conservative coefficient ledger.

Sources:

- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:status-open}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:status-open-poles}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:status-symbolic}`

### source-paper-4d-2pn--exploratory-review-status — Author-declared review caveat

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The author describes the work as exploratory, states that the accompanying symbolic checks establish internal consistency only, and says that subject-matter-expert review is still needed.

Sources:

- `research/4d_2pn/paper/4d_2pn.tex` — heading `\section*{Acknowledgments}`

## Revision and supersession relationships

### source-paper-4d-2pn--lineage-position — Carry-forward and downstream scope

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper carries forward the earlier 4D projection/worldtube framework and the frozen Newtonian and 1PN ledger
\[
L_0=\frac12m_Av_A^2+\frac12m_Bv_B^2+\frac{Gm_Am_B}{r},
\]
together with the source’s stated full 1PN block, and extends that lineage to conservative 2PN order. It also adopts the corrected Papers 7/8 matter–gauge ontology and does not declare whole-document supersession of those earlier works.

Sources:

- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:intro-motivation}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:carry-forward-import}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:carry-forward-1pn-ledger}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:carry-forward-L0}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{eq:carry-forward-L1}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:intro}`

## Related topics and scripts

No related memory pages or script domains were supplied by the task.