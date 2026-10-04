I'll start by reading the method packet and evidence guide, then follow those files into the actual sources and saved operands.The packet guide points to source, contracts, and saved operands. I’ll inventory those files and start reading the kernel and original-method specs.I’ll now pull the kernel source, source-contracts, and the key saved physical/template operands the method depends on.The kernel identities look consistent so far. Next I’ll check complete templates, whole J/D factors, physical bindings, and the address census.Templates confirm opposite normal signs. I’ll now check envelopes, address census, complete templates, and the original (K,T) windows.CLEAR FOR THIS BOUNDED FINITE-WINDOW CONTRACTION METHOD

The proposed finite-domain reorderings of ordinary J and the three added direct densities are algebraically consistent with `kernel_components` and with the saved complete J/D densities, on the original continuous window \(k,l\in[-K,K]\), \(t\in[-T,T]\). The two substitutions, clipped windows, orientation, remaining resolvents, and the wrong-root mutant are the right objects for a bounded algebra/domain/readiness certificate. That certificate still has to be derived, residual-checked, and collision-transported. This is not clearance of an author, a numerical action, or an evaluator.

## Algebra of the four primitives

Source `kernel_components` (`original-source/S11c_d_defect_packet_inner_lib.py` lines 53–60; contracts fragment `kernel_components`) is

\[
\begin{aligned}
J&=a\mu^2\frac{WL}{4}A(t)A(l-k-t)\frac{(k+t)(2k+t)q(k)}{q(k+t)\,E(l)\,E(k+t)\,E(k)\,(q(k+t)+q(k))},\\
D_{\mathrm{ref}}&=\frac{WL}{4i}A(t)A(l-k-t)\frac{(-i\mu)}{E(k)E(l)}\frac{k(2l-t)}{q(l-t)+q(l)},\\
D_{\mathrm{ht}}&=\frac{WL}{4i}A(t)A(l-k-t)\frac{(-i\mu)}{E(k)E(l)}\frac{k(t+2k)q(k)}{q(k+t)(q(k+t)+q(k))},\\
D_{\mathrm{qd}}&=\frac{WL}{4i}A(t)A(l-k-t)\frac{(-i\mu)}{E(k)E(l)}\frac{q(k)^2}{q(k+t)},
\end{aligned}
\]

with \(E=q+\beta\), \(q_i=q(k)\), \(q_o=q(l)\), \(q_h=q(k+t)\), \(q_s=q(l-t)\), and the three direct braces added. The same expected densities appear in `original-source/S11c_d_defect_packet_preflight.py` lines 275–277, already written as residual targets against saved Jwhole/Dwhole densities.

Physical joins at the saved binding match the method’s \(a,\mu,\beta,W,L\):

- \(\omega=3\), \(\rho_m=1/10\) \(\Rightarrow\) \(\mu=3/10\)
- \(a=\Lambda_{A0}/(\rho_m^2(1-i\omega\tau_A))=1/(1-3i/10)\)
- \(\beta=a\mu=(30+9i)/109\), with strictly positive real and imaginary parts
- \(W=1\), \(L=10\), \(A(z)=Lz/[4\sinh(\pi L z/2)]\), agreeing with `profile` and `saved/pressure/profile-envelope.json`
- \(\kappa=\sqrt{595}/10\), \(c_s=\sqrt{6}/2\) from `saved/ends/left-match.json` and `saved/preflight/physical-plan.json`
- \(q(p)=\sqrt{\kappa^2-p^2}\) on the closed first quadrant (`q` in inner_lib; piecewise `common_outgoing_q` in `actual-whole-density-argument-maps.json`)

\(C_d=(-i\mu)WL/(4i)\) is the product of the two source phase factors \(WL/(4i)\) and \((-i\mu)\). After that join it equals \(-\mu WL/4\). The future certificate still has to persist both factors and the residual against the saved density.

With \(x(k)=X(k)/E(k)\) and \(y(l)=N_f(l)Y(l)/E(l)\), the four contracted formulas reproduce the source polynomials:

| Primitive | Map | Remaining factors | Polynomial |
|---|---|---|---|
| BJ | \(m=k+t\), \(dt=dm\) | \(m/[q(m)E(m)]\), \(k\) clipped to \(I_T(m)\), \(l\) full | \((k+t)(2k+t)=m(k+m)\) \(\to\) \(C_{1T}+m C_{0T}\) |
| BD_height | same | \(1/q(m)\), no extra \(E(m)\) | \(k(t+2k)=k(k+m)\) \(\to\) \(C_{2T}+m C_{1T}\) |
| BD_quadratic | same | \(1/q(m)\), no \((q(m)+q(k))\) | \(q(k)^2\) \(\to\) \(X_{02T}\) |
| BD_reflected | \(m=l-t\) | no \(1/q(m)\) in the outer factor | \(k(2l-t)=k(l+m)\) \(\to\) \(X_1\cdot(Y_{1CT}+m Y_{0CT})\) |

\(E(k)\) and \(E(l)\) are absorbed into \(x\) and \(y\). \(E(k+t)=E(m)\) remains only in BJ, which is the only density that carries three resolvents. Direct densities never pick up an extra \(E(m)\). The four depths stay distinct: after \(m=k+t\) one has \(q(m),q(k),q(l)\); after \(m=l-t\) one has \(q(m),q(l),q(k)\). Mapping both roots to \(q(m)\) would be a different (wrong) identity.

These are author substitutions. Exact old-versus-new residuals against restored kernel ASTs and all 20 complete templates remain a **build obligation**.

## Finite domain

\(M=K+T\). For \((K,T)=(27,122)\) this is \(M=149\), with declared enlargement \((29,124)\), matching `pressure-readiness-method.md` §3 and `geometry/journal-result.json` (`K: [27,29]`, `T: [122,124]`, `U: 75`).

- \(m=k+t\) restricts the **input** \(k\) to \(I_T(m)=[\max(-K,m-T),\min(K,m+T)]\), empty when the lower end exceeds the upper end. Output \(l\) stays on \([-K,K]\). That is why BJ/height/quadratic use clipped \(C_{rT},X_{02T}\) and unclipped \(Y_0\).
- \(m=l-t\) has Jacobian \(-1\). \(t:-T\to T\) sends \(m:l+T\to l-T\); reversing the limits restores a positive measure. This map restricts the **output** \(l\) to the same endpoint formula \(I_T(m)\), while \(k\) stays on \([-K,K]\). That is why reflected uses unclipped \(X_1\) and clipped \(Y_{rCT}\). Equal endpoint formulas do not identify the two windows; the integrated variable is part of the identity.
- For \(T>K\), \(I_T(m)=[-K,K]\) on the interior \(|m|\le T-K\). The clipped wings are the proper subsets for \(T-K<|m|\le T+K\). Replacing \(I_T\) by \([-K,K]\) on the whole of \([-M,M]\), or dropping those wings, would change the original finite-window set. \(m\) is an internal transfer coordinate on \([-K-T,K+T]\), not a new external momentum limited to \([-K,K]\).
- Linear maps have absolute Jacobian one and are bijections of the finite box up to boundaries. The continuous integral over the original \((k,l,t)\) set is the object being rewritten. Equality of two different finite quadrature sums is not claimed.

The separate height-PV domain \(k\in[-K,K]\), \(Q\in[0,75]\) is untouched. Finite-\(T\) H, contact, paired PV, slope, and flat remain original numerical obligations. Analytic-H comparison in `InnerEvaluator.H` stays a check; it is not a substitute for finite-\(T\) H.

## Integrability and Fubini

On the compact box the packet transforms are smooth and bounded (polynomial-in-tanh times Schwartz Gaussians; `fields.json` `productTransformedAsWhole`, Fourier product construction). \(A\) is continuous on the reals. \(\beta\) has positive real and imaginary parts, so \(|q+\beta|\ge|\beta|\) and the method’s weaker bound \(|q+\beta|\ge|\beta|/\sqrt{2}\) holds. For \(q_1,q_2\) in the closed first quadrant,

\[
|q_1+q_2|\ge\frac{|q_1|+|q_2|}{\sqrt{2}},\qquad \Bigl|\frac{q_1}{q_1+q_2}\Bigr|\le\sqrt{2}.
\]

The only non-integrable-looking factors on the finite box are simple inverse square roots at \(q(m)=0\) or \(q(k)=0\) or \(q(l)=0\). Those are locally integrable. Simultaneous vanishing of two depths is a measure-zero set and is excluded from pointwise identities. The original 3D ordinary J/D densities are therefore absolutely integrable on the finite window, so Fubini and the Jacobian-1 substitutions are available.

This argument is a compact, cutoff-dependent majorant. It is not an unsubtracted-PV interchange, not a uniform-in-\((K,T)\) error bound, and not a restored machine theorem (`saved/pressure/whole-envelopes.json`: `machineMeasureTheoryProof: false`). Printed first-quadrant inequalities are not an executed certificate. New compact bounds still need source/physical joins in the next instrument.

## Provenance, templates, census

Complete pressure templates already contain the J/D wholes and the mixed H resolvents. They do not multiply a closed Dwhole by further resolvents (`typed-direct.json`: `multiplyWholeTagByResolvents: false`; `whole-definitions.json`: `wholeDirectOnce`, `nativeIterationOnce`). Mixed actual factor (`accepted-units/complete-template-plus-pressure-NATIVE_MIXED_ITERATION.json`) is the H term plus one `Jwhole_plus(...)`. Direct actual factor is one `Dwhole_*` (pressure) or \(\pm i q(l)\) times that whole (normal). `x` and \(y\) only factor resolvents already inside those densities.

Normal multipliers are the saved slot assignments, not face-label inferences:

- pressure: \(N_f=1\)
- plus-normal: \(+i q(l)\) (`complete-template-plus-normal-NATIVE_MIXED_ITERATION.json`)
- minus-normal: \(-i q(l)\) (`complete-template-minus-normal-INHERITED_DIRECT_WHOLE_OFF_DIAGONAL.json`)

Normal depth stays at \(l\). Inner_lib `run` uses the same \(\pm i q_o\) on the already-complete pressure templates. Epsilon is one homogeneous factor (`epsilonCount: 1` on selected addresses; `Fourier-order.json` `nativeEpsilonOnce: true`). Jet is \((-3i)^{n_t}(i/5)^{n_2}(i/10)^{n_3}\partial_x^{n_1}\) with derivative-before-coefficient (`fourier_lib.product`; `Fourier-order.json`). \(X=\widehat{b D_j u}\), \(Y=2\pi\widehat{c v}(-l)\), pairing \(\int Y F X\,dl\,dk\).

Address census from saved metadata, not a numerical nonzero proof: 544 selected of 2652 THETA rows (`saved/selected/pressure-addresses.json`; geometry journal `addresses: 544`, `formalAddresses: 102`, `explicitZeroAddresses: 442`). Zeros remain in scope. 20 face/slot/component adapters are present in `numeric-factor-adapters.json`. Plus/minus pressure templates agree; Jwhole_plus/Jwhole_minus and Dwhole_plus/Dwhole_minus are the same expected kernels after the preflight `Jvalue`/`Dvalue` replacement. 16 local cells and 34 coefficient fields are recorded in geometry scope. Units of the original summands are inherited (`accepted-units/journal-result.json`). Transport of those units through \(x=X/E(k)\) and \(y=N_f Y/E(l)\) is new contraction algebra and remains a **runtime join**.

## Wrong-root mutant and controls

The original inner mutant replaces \(q_s\) by \(q_h\) in the reflected term only (`middle(..., mutate=True)`). That integrand still contains \(q(k+t)\) and \(2l-t\), so its contraction uses \(m=k+t\), not the correct reflected map \(m=l-t\). Expanding \(2l-t=2l-m+k\) on \(k\in I_T(m)\) and \(l\in[-K,K]\) gives exactly

\[
BD_{\mathrm{reflected,wrong}}=C_d\int_{[-M,M]}\bigl\{2 X_{1T}(m)Y_{1C}(m)+[X_{2T}(m)-m X_{1T}(m)]Y_{0C}(m)\bigr\}\,dm.
\]

\(X_{1T},X_{2T}\) are the clipped \(X_{rT}\); \(Y_{rC}\) are the unclipped \(1/(q(m)+q(l))\) contractions. This formula is a proposed original-mutant join, not a computed packet value.

Algebra/domain controls named in the method (reflected window applied to the wrong variable, omitted wing, wrong reflected numerator sign, omitted term of \(m(k+m)\)) hit the actual maps. Full numerical H-contact, Leibniz, and normal \(q(k)\) versus \(q(l)\) controls remain the original addressed obligations, with the original ten-times-summed-empirical-envelope rule, on an applicable off-diagonal normal address. Silence on flat support does not establish that control. Domain refusal is distinct from numerical responsiveness.

## Next instrument

The limited next worker — restore kernel ASTs, complete adapters, source/test/unit records, bounds, profile and finite-window operands; derive substitution, factorization, integrability, domain, and control statements; persist old/new expressions and residuals; transport collision labels per primitive; evaluate no packet integral — is a sound certificate step. Requirements before that build:

1. Exact residual of each of BJ, BD_height, BD_quadratic, BD_reflected against restored `kernel_components` plus all 20 complete templates, with \(a,\mu,\beta,W,L\) joined and \(C_d\)’s two phase factors persisted.
2. Forward and inverse maps, orientation of \(m=l-t\), \(I_T\) membership, empty-window predicate, and both clipped wings, separately for each primitive.
3. Per-primitive fate of the original collision set \(\{k=\pm\kappa,\,l=\pm\kappa,\,k+l=0,\,k+l=\pm 2\kappa,\,l=k\}\) and profile/packet cuts: retained as a singular boundary, or absent because that primitive’s actual depths do not depend on that root. Old eight-plan \(k,l\) arrangements are not the new \(m\)-geometry.
4. Compact Fubini/majorant certificate with actual envelopes and the saved binding; excluded pointwise zero sets; no \(0/0\) assignment.
5. Unit transport through \(x\) and \(y\).
6. Chosen symbolic/domain witnesses for the algebra/domain mutants, with full residual or domain-witness returns.
7. Own build assessment of this certificate worker. Guard, supervisor, no-deadline pool, and completion-hook pinning are execution policy for that later launch, not part of this method verdict.

A full numerical evaluator still needs a reviewed concrete discretization of the new nesting, finite error allocation that propagates inner empirical error through actual outer weights, two independent implementations, storage plan, and the original H/contact/height-PV/slope/flat routes. Those remain unapproved. No runtime, memory, or algorithm forecast is made here. The geometry occurrence figure 8,304,551,424 is inherited unshared A24+A48 product-occurrence metadata across eight plans; it is not a unique request count and was not re-summed in this review.

The two original real-\(\omega=3\) Gaussian weak pairings are unchanged. There is no leakage, current, inverse, or sweep claim.

## Classification

- **Fatal method error:** none found in the proposed identities, windows, mutant, or template/resolvent accounting.
- **Build obligations:** items 1–7 above; they change what the certificate may claim and must be discharged by the next instrument.
- **Optional suggestion:** when persisting operands, keep the names \(X_1\) (unclipped) and \(X_{1T}\) (clipped) in the same table as the mutant, so the two reflected windows cannot be swapped by filename.
- **Limits of coverage:** no CAS residual was computed; Fubini is a hand majorant; 544/102/442, 20 templates, 16 cells, 34 fields, and geometry inspection status are inherited metadata; H analytic product, inner 38-point bank, and old Fourier bank were not treated as new proofs; eight-plan occurrence totals were not re-added; field polynomials and factor residuals were sampled, not exhaustively re-read.

## Files and sections read

`input/method.md`; `input/evidence-guide.md`; `input/review-prompt.md`; `input/packet-index.json` (records through factors/preflight); `input/packet-action-method.md` §§1–5; `input/pressure-readiness-method.md`; `input/source-contracts.json`; `input/accepted-unit-source-contracts.json`.

`input/original-source/S11c_d_defect_packet_inner_lib.py` (full, especially `kernel_components`, `q`, `profile`, `middle`, `H`, `run`); `input/original-source/S11c_d_defect_packet_preflight.py` lines 220–277; Fourier product/request fragments in source-contracts.

`input/saved/physical-input.json`; `saved/preflight/physical-plan.json`; `saved/preflight/numeric-factor-adapters.json` (all 20 definitions); `saved/native/binding.json`; `saved/ends/left-match.json`; `saved/pressure/{whole-definitions,typed-direct,whole-tags,fields,global-parameter-domain,profile-envelope,whole-envelopes,normal-growth,J-numerator-envelope,D-height-envelope,D-reflected-envelope,shift-root-bound,Fourier-order}.json`; `saved/whole-origin/{Jwhole,Dwhole,H}.json`; `saved/inventory/{actual-whole-density-argument-maps,fourier-and-unit-provenance}.json`; `saved/selected/pressure-addresses.json` (header and first entry); `saved/inner/physical-H-and-PV-lemma.json` (native binding header); `saved/local/assembly.json`.

`input/accepted-units/journal-result.json`; `accepted-units/complete-summand-dimensions.json` (start); `accepted-units/complete-template-{plus-pressure-NATIVE_MIXED_ITERATION,plus-pressure-INHERITED_DIRECT_WHOLE_OFF_DIAGONAL,plus-normal-NATIVE_MIXED_ITERATION,minus-normal-INHERITED_DIRECT_WHOLE_OFF_DIAGONAL}.json`.

`input/geometry/journal-result.json` (plan header, closing totals/scope); `input/geometry/dependency-families.json` (start).

Not exhaustively read: all 544 address rows, all 34 field reconstruction triplets, all 17 factor operand sets, all eight full arrangements, tail JSON, and numeric-factor argument triplets as separate files (the 20 adapter definitions and four complete templates were used instead). Those omissions do not change the algebraic identities above; they remain joins the certificate worker must restore in full.