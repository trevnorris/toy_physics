# Puncture-deflection electric-sign result

## Concise result (body; ≤2 pages)

**Scope.** This is the target-blind Tier-A far-field calculation for \(R\gg r_e\). G0 freezes the sleeve \(\Sigma\), U2 leaves its contact operator \(\mathcal B\) unresolved, and S_hold constrains only \(r_B-\tfrac12\). The committed bare model therefore does not select a mouth boundary class.

**Q-FIELD / Q-AMEND.** The geometric field is the normal embedding displacement

\[
\xi_w(x)=\ell h(x),\qquad [\xi_w]=L,
\]

with \(H=f_0h+H_\perp\), \(f_0=[\ell\cosh^2(w/\ell)]^{-1}\), \(N_0=8/(3\ell)\), and \(h=P_0H\). The preregistered held-datum candidate is \(h_A=\xi_w|_A/\ell=P_0H|_A\), dimensionless. The length-valued \(u_L\) remains distinct; \(C_{hu}\) mixes it with \(h\) but does not identify them. Since \(f_0(0)=1/\ell\), the live parent coupling \(-J_mQ_\chi H\) reduces to \(-g_{\chi h}Q_\chi h\).

REPLACE re-types the existing \(h\)-source BC as the candidate holding BC and adds zero rows. ADD retains that source and adds exactly one \(R_w\)-odd core-\(h\) holding row. All 13 other G0 [POSTULATE] zero-ledger rows remain zero; S_hold remains scoped to \(r_B-\tfrac12\). Fact: **internal_inconsistency = none**.

**Q-BC.** The exterior stationary solution is \(h(r)=h_Aa/r\), with positive holding curvature \(4\pi(D/B_{\rm eff})a\). A nonzero signed \(h_A\) is not stationary without a core term. The \(\pm w\) label is invariant, but G0 supplies neither the nonlinear core functional nor the barrier/map \(s\mapsto h_A\). One common classifier returns:

**UNDETERMINED_ANALYTICALLY(missing parent-throat/boundary closure)**.

Its controls return FIXED_SOURCE for a free relaxing mouth, DIRICHLET_VALUE for an imposed value, FIXED_MONOPOLE for imposed conormal flux, and MIXED for an imposed value/conormal relation. All four classes remain admissible.

**Q-FORCE / Q-COMBINE.** Put \(b=B_{\rm eff}\), \(k=K_h\), \(c=C_{hu}\), \(D=bk-c^2>0\). With the exact Robin escape factor \(z_g>0\),

\[
m_{gg}=\frac{b z_g^2}{D},\qquad \det m=\frac{z_g^2}{D},\qquad D_*=\frac74.
\]

For self response \(S_{gg}>0\), held value \(\phi\), core monopole/source magnitudes \(q,j>0\), committed source \(g=g_{\chi h}>0\), and \(0\le\lambda\le1\),

\[
U_X=\frac{s_1s_2A_X}{4\pi R},\qquad
F_{X,\mathrm{out}}=\frac{s_1s_2A_X}{4\pi R^2},
\]

\[
A_V=\frac{m_{gg}\phi^2}{S_{gg}^2},\quad
A_M=m_{gg}(q^2-g^2),\quad
A_J=-m_{gg}(j+g)^2,
\]

\[
A_{\rm MIXED}(\lambda)
=m_{gg}[(1-2\lambda)q^2-2\lambda qg-g^2].
\]

MIXED has endpoints \(A_M\) and \(-m_{gg}(q+g)^2\), and algebraic zero \(\lambda=(q-g)/(2q)\) when admissible. V is neutral-positive, M indefinite, J neutral-negative, and MIXED spans negative/null/positive over its full parameter-and-magnitude range. V and J flip under their bare-member mutations; M changes when the held-\(h\) work subtraction is mishandled.

REPLACE totals are \(\{m_{gg}\phi^2/S_{gg}^2,m_{gg}q^2,-m_{gg}j^2,m_{gg}(1-2\lambda)q^2\}\). ADD totals are the four \(A_X\) above. The realization fact is **variant_unresolved**. Across all admissible classes the full fact is **outcome_not_invariant**.

**Q-MAG / hooks.** Write \(a=c_ar_e\), \(r_e=k_ee^2/(m_ec^2)\), and \(\xi_A=c_\xi a\). Neither \(c_a\) nor \(c_\xi\) is fixed by the reduced action; core uncertainty is unbounded at this tier, while the far-field truncation is \(O(a/R)\). Facts: **magnitude_free_factor**, **tier_A_conditional**. Hooks: density **NO** (no_local_prediction; \(B_{\rm eff}=\rho_{B0}^2/\chi_c\) is background-only and local \(K_h,C_{hu},z_g\) laws are absent); radial monopole **UNDETERMINED**; universal quantization **NO** (continuous core modulus); close range **UNDETERMINED / out of scope**.

**§4 landing:** **R1_REQUIRED(bc_selection)**. The admissible classes do not share one outcome value, so the class gap fires before the later variant, tier, range, or magnitude-normalization gaps. This is neither an inconsistency nor a completed target mismatch.

**DECIDED:** field identity and live \(Q_\chi H\) reduction; both scoped amendments; positive coupled kernel and \(D_*=7/4\); every V/M/J/MIXED coefficient and \(1/R^2\) force; conditional REPLACE/ADD totals; free magnitude modulus; hooks; no internal inconsistency. **R1-deferred:** the nonlinear parent-throat boundary functional/barrier selecting the class, MIXED parameters if selected, REPLACE versus ADD, dimensionless throat/deflection normalization, and whether the signed core has a nonzero radial monopole. R1 may test branches; physical selection requires the specific throat solve/selection criterion.

---

## Appendix (evidence; exempt from the body cap)

### A. Field, action, kernel, and ledger

The term-for-term reduced action is

\[
S_{Lh}=\int dt\,d^3x[
\tfrac12A_{\rm eff}\dot u_L^2+\tfrac12M_h\dot h^2
-\tfrac12B_{\rm eff}|\nabla u_L|^2-\tfrac12K_h|\nabla h|^2
-C_{hu}\nabla u_L\cdot\nabla h].
\]

The parent mouth functional is

\[
\Omega_{\rm mouth}=\sum_i\int\eta_i[
\tfrac12K_mH(x,0)^2-J_mQ_\chi[r_\Sigma,s_i]H(x,0)]\,d^3x,
\]

reducing to \(\eta_i(k_mh-g_{\chi h}s_i)\). S_hold contains none of \(H,h,u_L\). REPLACE changes only that existing source-BC row. ADD activates only core_embedding_h_holding_row, with coefficient dimension \(E\), datum \(h\) dimension 1, and odd \(R_w\) parity.

Completing the square gives \(\kappa=D/b\) and

\[
m=Z^T\begin{pmatrix}b^{-1}&0\\0&\kappa^{-1}\end{pmatrix}Z,\qquad
Z=\begin{pmatrix}1&0\\-(c/b)z_b&z_g\end{pmatrix},
\]

\[
m_{uu}=\frac{D+c^2z_b^2}{bD},\quad
m_{ug}=-\frac{cz_bz_g}{D},\quad
m_{gg}=\frac{bz_g^2}{D},\quad
\det m=\frac{z_g^2}{D}>0.
\]

Here \(z_g=1-k_m\langle\eta,L_h^{-1}\eta\rangle\) and \(z_b=1-k_m\langle\eta,L_h^{-1}f_b\rangle\). At \(z_g=z_b=1\), the G0 witness is \(m=\frac17\left(\begin{smallmatrix}4&-2\\-2&8\end{smallmatrix}\right)\).

### B. Full ensemble definitions

Let \(y_i\) be the total \(h\)-channel source/conormal and \(E_0[y]\) its self-plus-pair stored energy. The one-body response is \(h_A=\partial E_0/\partial y\), with self coefficient \(S_{gg}\). Every row subtracts its self term.

| class | held data / total source | conjugate functional | pair coefficient | wrong-functional control |
|---|---|---|---|---|
| V / DIRICHLET_VALUE | \(h_{A,i}=s_i\phi\); \(y_i=s_i(g+Q_i)\), \(Q_i\) reacts | \(\Omega_V=E_0-gh-Qh=E_0-yh\) | \(+m_{gg}\phi^2/S_{gg}^2\) | bare \(E_0\) gives the negative |
| M / FIXED_MONOPOLE | core conormal \(s_iq\), committed \(s_ig\); \(y_i=s_i(q+g)\) | \(\Omega_M=E_0-gh\) | \(m_{gg}(q^2-g^2)\) | omitting held-\(h\) work gives \(m_{gg}(q+g)^2\) |
| J / FIXED_SOURCE | core source \(s_ij\), committed \(s_ig\); \(y_i=s_i(j+g)\) | \(\Omega_J=E_0-jh-gh\) | \(-m_{gg}(j+g)^2\) | bare \(E_0\) gives the positive |
| MIXED | core \(s_iq\), reservoir fraction \(0\le\lambda\le1\), committed \(s_ig\) | \(\Omega_\lambda=E_0-gh-\lambda qh\) | \(m_{gg}[(1-2\lambda)q^2-2\lambda qg-g^2]\) | omitting \(-\lambda qh\) collapses to M |

V is derived by exactly inverting the two-mouth response block
\(\left(\begin{smallmatrix}S_{gg}&\epsilon m_{gg}\\\epsilon m_{gg}&S_{gg}\end{smallmatrix}\right)\)
and extracting \(O(\epsilon s_1s_2)\). M/J/MIXED follow by subtracting the displayed reservoir work from the stored pair term \(m_{gg}y_1y_2\). This re-derives the ensembles for \(h=\xi_w/\ell\); it does not relabel the prior \(u_L\) datum.

| class | REPLACE total / fact | ADD total / fact |
|---|---|---|
| DIRICHLET_VALUE | \(m_{gg}\phi^2/S_{gg}^2\); POSITIVE_R2 | same; POSITIVE_R2 |
| FIXED_MONOPOLE | \(m_{gg}q^2\); POSITIVE_R2 | \(m_{gg}(q^2-g^2)\); outcome_not_invariant |
| FIXED_SOURCE | \(-m_{gg}j^2\); NEGATIVE_R2 | \(-m_{gg}(j+g)^2\); NEGATIVE_R2 |
| MIXED | \(m_{gg}(1-2\lambda)q^2\); outcome_not_invariant | \(A_{\rm MIXED}\); outcome_not_invariant |

### C. Q-BC evidence

The exterior functional gives

\[
E_{\rm ext}=2\pi\kappa a h_A^2,\quad
\partial_{h_A}E_{\rm ext}=4\pi\kappa a h_A,\quad
\partial_{h_A}^2E_{\rm ext}=4\pi\kappa a>0.
\]

| candidate | admissible variation | signed stationary | positive curvature | barrier computed | classifier output |
|---|---:|---:|---:|---:|---|
| actual core candidate | YES | NO | YES | NO | UNDETERMINED_ANALYTICALLY(missing parent-throat/boundary closure) |
| free-mouth control | YES | NO | YES | NO | FIXED_SOURCE |
| imposed-value control | NO | YES | YES | YES | DIRICHLET_VALUE |
| imposed-conormal control | YES | not required | YES | YES | FIXED_MONOPOLE |
| imposed mixed control | YES | relation-derived | YES | YES | MIXED |

Mutating the actual missing-closure premise returns the relaxing class. Mutating each control's decisive value, conormal, or mixed premise destroys only its own typed result.

### D. Exhaustive sealed §4 landing table

The full Cartesian domain is

\[
5\ {\rm QBC}\times6\ A_{\rm replace}\times6\ A_{\rm add}\times4\ {\rm variants}
\times2\ {\rm magnitude}\times2\ {\rm tier}\times2\ {\rm inconsistency}
\times2\ {\rm all\mbox{-}class\ invariance}\times2\ {\rm MIXED\ invariance}
=23040.
\]

QBC values are DIRICHLET_VALUE, FIXED_MONOPOLE, FIXED_SOURCE, MIXED, and UNDETERMINED_ANALYTICALLY. Each variant outcome is one of POSITIVE_R2, NEGATIVE_R2, NULL, POSITIVE_WRONG_RANGE, NEGATIVE_WRONG_RANGE, or outcome_not_invariant. Variant realization is replace, add, both, or variant_unresolved. The two invariance booleans are computed §4 summaries, not physical inputs.

This grouped table is lossless and exhaustive: predicates apply in order, are mutually exclusive after earlier matches, and counts sum to 23,040.

| first matching predicate | landing | cells |
|---|---|---:|
| internal inconsistency | NO_GO(sector) | 11,520 |
| unconditional; POSITIVE_R2; no free factor | SIGN_EARNED | 252 |
| unconditional; POSITIVE_R2; free factor | R1_REQUIRED(magnitude) | 252 |
| unconditional; NEGATIVE_R2 | MECHANISM_FALSIFIED(wrong_sign) | 504 |
| unconditional; either wrong-range outcome | MECHANISM_FALSIFIED(wrong_range) | 1,008 |
| unconditional; NULL | R1_REQUIRED(subleading) | 504 |
| QBC unresolved and outcomes disagree across classes/models | R1_REQUIRED(bc_selection) | 1,152 |
| MIXED parameters unfixed and range outcome not invariant | R1_REQUIRED(mixed_bc_parameters) | 1,152 |
| realization non-single and REPLACE/ADD differ | R1_REQUIRED(variant_selection) | 3,840 |
| Tier-A conditional or remaining outcome_not_invariant | R1_REQUIRED(sign_and_magnitude) | 2,856 |
| defensive fall-through | R1_REQUIRED(unclassified) | 0 |

“Unconditional” requires class resolution by a unique Tier-B-closed selection or full class/MIXED-range outcome invariance; variant resolution by a single realization or equal outcomes; and a constant sign/range or invariant null over the entire magnitude range.

Both engines emit exhaustive-row digest:

**7627417ace0f99a17187a90efa2523ca9b68df8b64f7960d38be2dc6f553ac49**.

The production tuple has unresolved QBC, unequal class outcomes, variant_unresolved, magnitude_free_factor, tier_A_conditional, and no inconsistency. Its first match is exactly **R1_REQUIRED(bc_selection)**.

### E. Production checks

Both scripts exit 0, agree term-by-term on every symbolic payload entry and the truth-table digest, and emit all 36 atomic teeth. The production landing is recomputed from the live typed Q-BC, Q-COMBINE, variant, Q-MAG, tier, and Q-AMEND facts; `LANDING_OWNERSHIP` rejects an upstream landing injection and requires the emitted landing to equal that recomputation. Units restored: \([A]=EL\), \([U]=E\), \([F]=EL^{-1}\). Separate range controls cover a sign flip, a zero touch without sign flip, and a derived subdominant constant outcome.
