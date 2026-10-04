I will read the method proposal and evidence guide first, then the source and saved operands they name, and assess only this bounded numerical J/direct packet method.The proposal names a constant-field Fourier exception and a new panel plan. I will open the source index and packet index next and follow the operands they cite.CLEAR FOR THIS BOUNDED NUMERICAL J/D PACKET METHOD

The constant-field exception, the four contracted primitives, the six-slab window plan, the inherited tail expressions, and the finite-window controls are mathematically consistent with the saved operands and the original source. No fatal method error turned up. This reading did not execute a worker, certify a number, or clear a future build.

What was read: `input/method.md`, `input/evidence-guide.md`, the full applicable block and all 20 eligible records in `input/constant-subset-source-index.json`, `input/source-contracts.json`, `original-source/S11c_d_defect_packet_fourier_lib.py`, `original-source/S11c_d_defect_packet_inner_lib.py`, `original-source/S11c_d_defect_packet_preflight.py`, `original-source/S11c_d_defect_packet_evidence_store.py`, `runtime-input/saved/selected/pressure-addresses.json` (header), `runtime-input/saved/preflight/physical-plan.json`, `runtime-input/preflight/tail-plan.json` (header), `runtime-input/preflight/tail-bound-derivation.json` (through the exponential moments), `runtime-input/saved/inventory/fourier-and-unit-provenance.json` (opening contract), the `e_W_t` polynomial and its reconstruction return, `rules/extraction.json`, and `storage/preparation.json`. The 4 MB address file and the 18 MB inventory were not walked object by object.

## 1. Census and constant-field transforms

The selection header records `originalCount` 2652, `selectedCount` 544, and `formalCount` 102. The preflight checker requires the status split 102 nonzero, 106 source-jet zeros, and 336 consumer zeros. The index’s applicable block is exactly 64 addresses: 8 jets × {J, direct} × {plus, minus} × {pressure, normal}.

The 20 eligible entries are the pressure, grade 00/00/11, nonzero rows. Ten are `NATIVE_MIXED_ITERATION` and ten are `INHERITED_DIRECT_WHOLE_OFF_DIAGONAL`, five jets on each face: `e_W`, `e_W_d1d1`, `e_W_d2d2`, `e_W_d3d3`, `e_W_t`. The other 44 are inherited zeros: 12 pressure odd jets `EXACT_ZERO_SOURCE_JET`, and all 32 normal J/direct rows `EXACT_ZERO_CONSUMER`. A runtime census other than this is a refusal, as the proposal says.

Every eligible source and consumer record read has degree 0, denominator a positive integer, and reconstruction return `cancelled = 0`. The `e_W_t` source is numerator 1 over denominator 2, so the coefficient is the saved value `1/2`. Odd-jet and normal zeros stay proved zeros and do not become transform requests.

The displayed transforms match the original convention in `constant_reference` and `product`:

- \(G_u\) peaks at \(k=p_0\) with phase \(-i(k-p_0)x_u\), \(x_u=-5/2\).
- \(G_v\) is the dual at argument \(-l\), carrier \(-p_0\), center \(x_v=5/2\). Its \(\nu=p_0-l\), and \(\exp(-i\nu x_v)=\exp(+i(l-p_0)x_v)\), which is the displayed \(G_v\).
- \(X\) carries one factor \(1/(2\pi)\). \(Y\) carries none. That is the reduced one-dimensional factor after the two edge deltas in `fourier-and-unit-provenance.json`, together with the dual test having no extra power.
- Derivatives hit the packet before multiplication by the constant \(b\). The Route B recurrence \(Q\mapsto Q'+(ip_0-w/s^2)Q\) on \(w=x-x_u\), and the centered moments \(M_0\), \(M_1=-is^2\nu M_0\), \(M_{r+1}=-is^2\nu M_r+rs^2 M_{r-1}\), are the Gaussian moment identities for that derivative. For \(Y\), derivative order is 0 and the center is \(x_v\).
- The two routes are two writings of one analytic identity. The \(10^{-24}\) full-value gate and the Gaussian-stripped complex gate are roundoff checks. The stripped gate is on the complex difference, with phases left inside \(V\), and \(M_r/M_0\) is formed from the recurrence so the exponential is not divided out numerically. Actual discrepancies, not the thresholds, enter the error path.

Required runtime join, once: the preflight identity is

`waveMultiplier = (-3i)^{nt} (i p)^{n1} (i/5)^{n2} (i/10)^{n3}`.

Delta support `k=p` replaces that factor by \(P_j(ik)^{n_1}\). It is the same factor as the Fourier derivative, not a second copy. Checked examples: `e_W` gives `1`; `e_W_d1d1` gives `-composition_p**2`; `e_W_d2d2` gives `-1/25`; `e_W_d3d3` gives `-1/100`; `e_W_t` gives `-3*I`. The three second-derivative jets share coefficient id `69167940…` and field value `-1/109-3I/1090`. The jet factor is only in `waveMultiplier`. A field id, a constant flag, or the value at \(x=0\) does not separate them. `argumentDerivative` stays 0 because \(k^r\) and \(l^r\) sit outside \(X\) and \(Y\). `normalMultiplier` is `1` on these 20 pressure rows. Epsilon is already extracted once: `consumerOriginal` \(-I\epsilon(3-10I)/109\) equals the saved consumer `-10/109-3I/109`.

## 2. Four primitives

`kernel_components` and the contraction driver in `source-contracts.json` support the displayed formulas. \(E=q+\beta\), \(q\) is the positive-real / positive-imaginary square root, \(a=100/109+30i/109\), \(\mu=3/10\), \(\beta=a\mu=30/109+9i/109\), \(C_J=a\mu^2 WL/4\), and \(C_D=(-i\mu)WL/(4i)\). Both equal the prefactor \(WL(-i\mu)/(4i)\). \(A(z)=Lz/(4\sinh(\pi L z/2))\) with \(L=10\) matches `profile`, and \(A(0)=1/(2\pi)\).

\(J\), height, and quadratic clip input \(k\) on \(I(m)\). Reflected \(D_r\) clips output \(l\). \(X_1\) and \(Y_0\) stay on \([-K,K]\). The reflected root remains \(q(l-t)\). Substituting \(t=l-m\) produces \(k(l+m)\) and the displayed \(D_r\). The wrong-root substitution \(q(l-t)\to q(k+t)\), \(t=m-k\), clipped input and unclipped output, produces \(k(k-m+2l)\) and the displayed \(D_{r,\mathrm{wrong}}\). Wings run to \(M=K+T\), so \(M=149\) and \(M=153\), with kinks at \(|m|=95\).

\(J\) is the `Jwhole` addend only. Each mixed template still contains the separate `Hwhole` term. A normalized family \(x_n=k^n G_u/(2\pi E)\), \(y_0=G_v/E\), \(\alpha=bc P_j i^n\), with \(n\in\{0,2\}\), is exact linearity after the joins above. Face sharing needs the whole template and units. Restoring a small \(\alpha\) must not hide another address: every original address is checked after its own \(\alpha\) is put back, and unused primitive slots do not donate budget.

## 3. Panel geometry

For both windows, \(T-K=95>\kappa=\sqrt{595}/10\). The sorted set \(\{-M,-(T-K),-\kappa,0,\kappa,T-K,M\}\) cuts \((-M,M)\) into six affinity slabs before any branch selection. On those slabs the true window is \([-K,m+T]\), \([-K,K]\), and \([m-T,K]\) respectively, and adjacent formulas agree at \(\pm(T-K)\). The endpoints \(\pm M\) are single points and stay unsampled. \(d_g=\min(|m-\kappa|,|m+\kappa|)\) is \(-m-\kappa\), \(-m-\kappa\), \(m+\kappa\), \(\kappa-m\), \(m-\kappa\), \(m-\kappa\) on the six slabs, so both signs of \(\pm\kappa\pm d_g\) are affine there. One of those cuts is \(z=m\), already present as the \(d=0\) profile cut.

Fixed cuts, carrier offsets about \(p_0\), profile offsets about \(m\), grazing cuts, and clipped endpoints are intersected inside each slab, then clipped to the true interval. Coalescing is exact equality only. \(q(m)+q(z)=0\) on the common outgoing sheet only when both roots vanish, so those points stay endpoints and no \(0/0\) value is filled in. The opposite root is absent from each separated primitive.

The independent check against `max(-K,m-T)` and `min(K,m+T)`, with an interior witness, slope order, disjoint interiors, and separate \(k\) versus \(l\) coverage, is the right validator. Forcing \([-K,K]\) onto a wing must fail that same validator. Route A’s open squared map matches `gauss`: Jacobian \(2\cdot\mathrm{half}\cdot z\) and the extra weight \(1/2\), roots never nodes, no later panel deletion. Route B uses the same boundaries in the physical variable with G7/K15, global largest-error splits, and no A nodes. That pair is a suitable test of this bounded integral. A failing fixed A budget stops; it does not add orders.

## 4. Error transport and tails

The product rule \(|a|e_b+|b|e_a+e_a e_b\), absolute linear sums, and application before \(C_J\), \(C_D\), resolvents, powers, \(\alpha\), weights, and Jacobians are the right empirical propagation. They are not a forward enclosure.

The internal budget \(10^{-11}/80\) per eligible address and primitive, per carrier and window, matches a 20×4 split with no donation. The comparison tolerance \(10^{-9}+10^{-7}|A_{48}|\) stays the outer acceptance test. For B, \(w(m)=1+1/|q(m)|\) and \(W_M=3M+4\) are consistent: \(\int_{-M}^{M} dm/|q|=\pi+2\mathrm{acosh}(M/\kappa)\), and \(\kappa>2\), \(M>\kappa\) give \(\int w<W_M\). The \(1/|q|\) and the \(1\) in \(w\) are values in the original momentum unit and have to be nondimensionalized in the build. The pointwise target \(\epsilon w/(16 W_M)\) integrates to \(\epsilon/16\) under the lemma; the separate discrete check is \(\epsilon/4\). The lemma does not certify that discrete sum. Same-\(m\) inner GL24/GL48 differences stay on each route’s own nodes.

`tail-plan.json` is \(K=27\), \(U=75\), \(T=122\). `tail-bound-derivation.json` stores \(b_*=3000/11101\), and \(F_3=472/625\) matches `weighted_tail(3,0)`. The displayed formulas match `contributions` in the preflight source:

- outer \(=2 C_{\mathrm{ordinary}} C_X C_Y E_{30} F_3 \mathrm{WT}_3(K,5)\)
- middle J \(=4 C_X C_Y E_{30} F_3^2 (4\cdot 121/(5 b_*^3)) \mathrm{WT}_2(T,1)\)
- middle D \(=4 C_X C_Y E_{30} F_3^2 (36\cdot 121/b_*^2) \mathrm{WT}_1(T,1)\)

\(b_*\) is that positive beta lower bound, not the address coefficient \(b\). The extra mixed middle term \((4/b_*) C_X C_Y E_{30} F_2^2 (55/3) 2^{-T}\) is the H piece. Keeping the full mixed record is an overcount, not an H value. K29/T124 is the same expressions at the new integers, with the K27/T122 returns restored and the radius loop not replayed. There is no \(x\)-truncation share in the analytic Gaussian. The compact Fubini majorants stay off the all-real tail. No cancellation pays a tail or an error budget.

Before a direct primitive uses the single middle-D expression as its own upper bound, the build has to join the inherited absolute/triangle majorant of that primitive. A bound on a signed sum is not automatically a bound on each piece. Duplicating the shared D bound is overcount bookkeeping, not a budget transfer.

## 5. Controls

The wrong-root control belongs on the smallest eligible direct id, `8347`, separately for each carrier, on the distinct clipped-input / unclipped-output domain, on both windows. The derivative mutant belongs on the smallest \(n_1=2\) id, `8350`, which is a J address: \((ik)^2\) is replaced by \((ip_0)^2\) inside that J assembly. For \(p_0=0\) the mutant factor is zero. Movement must exceed ten times the summed finite-window empirical envelopes or the control is silent. Baseline all-real tails do not attach to either mutated kernel. H, normal-slot, and nonconstant Leibniz controls stay deferred with the pressure terms this subset does not evaluate. No responsiveness result is claimed.

## 6. Requests, storage, and the build gate

Reuse needs the full operand join and the original return receipt. A looser saved target does not satisfy a tighter one. A24, A48, and B50 do not share evaluated transforms, contractions, adaptive choices, or actions. Immutable rule payloads may be shared. `rules/extraction.json` only points at the saved GL24, GL48, and G7/K15 receipts.

`EvidenceStore` is an append-only writer with FULL transactions, route namespaces, and byte reconstruction. `storage/preparation.json` records 37 synthetic tests, `scientificComputationsRun: 0`, and `independentBuildClearance: false`. It has no scientific request index. The build still has to add immutable full-operand lookup, collision refusal, purpose in the namespace (`baseline`, `formula-check`, or a named mutant), MP-tuple records, failure prefixes, the 20 GiB disk reserve, the one-panel-at-a-time cap, and the 64 MiB route-local LRU. Cost and runtime are unknown. Passing this method assessment does not launch a worker. The concrete build needs its own review. Inferred gamma units and the accepted unit algebra remain inherited dependencies.

Mandatory build joins are the census and one-time wave/delta/denominator/epsilon joins, equality of the moment form to the displayed transform before quadrature, the six-slab validator and wing refusal, the per-address scaled budget and discrete \(w\)-sum, the K29/T124 tail expressions with the direct absolute-majorant join, the two finite-window mutants, and the request/storage semantics above. Wording that would help the build, without changing the mathematics, is to write \(Y\)’s center factor \(\exp(-i(p_0-l)x_v)\) in the same sentence as the \(X\) phase, and to write the family test as \(|\alpha_a|\) times the normalized indicator against each address budget.