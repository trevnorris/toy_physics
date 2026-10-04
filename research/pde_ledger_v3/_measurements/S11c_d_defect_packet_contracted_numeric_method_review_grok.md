I'll read the method proposal, the evidence guide, and the source and saved operands they name, then assess the bounded numerical J/direct packet on those merits.The proposal and evidence map are in. Next I’ll open the source index and the files it names so the census, joins, and formulas can be checked against the saved operands.The Fourier library and the constant-subset index are open. I’ll check the contraction identities, the census, and the panel and tail operands against the saved records.The assessment is still open. I am re-reading the proposal’s panel claim and the window arithmetic it depends on before giving the verdict.The proposal’s four-sector affinity sentence is the claim that can change the integration domain. I am checking it against the clipping rule and the saved kernel, Fourier, and census operands.NEEDS REVISION

The clipped-window endpoints are not affine on the four sectors the panel plan uses. Both routes are required to share those boundaries, so agreement between them would not catch a shared wrong wing. That has to be corrected before the fixed A and adaptive B prescriptions are a valid test of this bounded ordinary-J and three-direct target.

## Clipped endpoints

`input/method.md` lines 178–184 say that on each of the four sectors divided by −κ, 0, and κ, the inner-boundary functions are affine in m, and they include the clipped window endpoints in that list. The contraction window is

I(m) = [max(−K, m−T), min(K, m+T)],

with (K, T) = (27, 122) and (29, 124), so T−K = 95 and M = K+T is 149 or 153 (`input/contraction-method.md` lines 102–139; the same windows are restated in `input/method.md` lines 141–144). κ² = 595/100 = 5.95, so κ < 3. Thus ±95 lie strictly inside the outer sectors [−M, −κ] and [κ, M].

The right endpoint min(K, m+T) changes slope at m = −95. The left endpoint max(−K, m−T) changes slope at m = +95. On [−κ, 0] and [0, κ] both endpoints are the constants ±K, so those two sectors are affine. On the two outer sectors a single affine expression cannot represent the kink. An intersection sweep that uses one active endpoint formula per outer sector drops or extends the clipped wing. J, Dh, and Dq clip the input; Dr clips the output. Those wings are different objects (`input/contraction-method.md` lines 132–139).

The m-plan does list ±(T−K) as sample endpoints (`input/method.md` line 178). That list does not replace the four-sector affinity claim. Route B is required to use the same mathematical boundaries as Route A (line 197).

The grazing scale d_g(m) = min(|m−κ|, |m+κ|) is affine on those four sectors: m−κ, κ−m, κ+m, and −m−κ respectively. Fixed cuts and the shifts m+d are affine on the whole line. At m = κ, d_g = 0 and the resolution cuts coalesce onto ±κ; the method correctly keeps both labels and assigns no pointwise 0/0 value.

Required change: split the affinity slabs at ±(T−K), or carry both candidate graphs z = ±K and z = m±T on every slab and activate them by max/min. Then check order, disjoint open intervals, and exact coverage against the true I(m), separately for input clipping and output clipping. After that split, open squared GL24/GL48 with a fail-closed stop, and independent physical G7/K15 refinement, are a suitable pair for this empirical target. A coarse A panel may stop; that stop is part of the proposal.

## Constant-field transform

For exactly joined degree-zero fields, the displayed formulas match the saved convention. With hat f(s) = (1/(2π)) ∫ exp(−i s x) f(x) dx, X = hat[b D_j u] and Y = 2π hat[c v](−l) (`input/packet-action-method.md` lines 78–82; `input/original-source/S11c_d_defect_packet_fourier_lib.py` `product` and `constant_reference`).

G_u(k) = s √(2π) exp(−s²(k−p0)²/2 − i(k−p0)x_u) and X(k) = b P_j (i k)^n G_u(k)/(2π). G_v(l) = s √(2π) exp(−s²(l−p0)²/2 + i(l−p0)x_v) and Y(l) = c G_v(l). Derivatives act on u before multiplication by b. The carrier derivative is i k. Y uses Fourier argument −l, test carrier −p0, center x_v, and no extra 1/(2π). There is no extra conjugation. P_j = (−3i)^nt (i/5)^n2 (i/10)^n3, and the held frequency in `input/runtime-input/saved/preflight/physical-plan.json` is 3. The raw override recorded in `input/original-source/S11c_d_defect_packet_preflight.py` line 158 is old frequency 1, actual 3.

The moment formulas M0, M1 = −i s² ν M0, and M(r+1) = −i s² ν M_r + r s² M(r−1) agree with that transform for the centered Gaussian. Route A’s (i k)^n form and Route B’s moment form are one analytic identity. The absolute 1e−12 comparison is a roundoff check of that identity. The proposal says so (`input/method.md` lines 102–109). Old finite Fourier banks are not a uniform accuracy proof.

Route B writes Q_next = Q′ + (i p0 − w/s²) Q without defining w. The packet-action recurrence is in the physical coordinate, Q_{n+1}(z) = Q′ + (i p0 − (z−x_u)/s²) Q (`input/packet-action-method.md` lines 66–69). The Fourier library builds the same polynomial in w = z − center (`fourier_lib.py` lines 171–173). The build has to set w = x − x_u for X and w = x − x_v with carrier −p0 for Y, then require exact equality with the displayed (i k)^n form before quadrature. Route B must not also multiply by (i k)^n. A wrong absolute-x reading fails that equality check. That is a mandatory join of an already stated gate.

Nonconstant fields stay off the exception. `f588c943…-polynomial.json` is degree 1 in tanh(composition_x/10), denominator 6976. The zero field `9937a9eb…` in `fields.json` stays a proved zero. Constant flags are not enough: `08a5aec4…-polynomial.json` stores numerator 1 and denominator 2 for the value 1/2. Address 8346’s own source and consumer polynomials, as copied into `constant-subset-source-index.json`, are degree 0 with denominator 1, transfers k−p and r−l, and pressure factor 1. Its response is the H term plus Jwhole_plus. Ordinary J is one addend of native mixed.

## Census and the four contractions

`constant-subset-source-index.json` records 544, 64, and 20. The applicable-metadata array is 64 records of the same 29-line shape, from line 8 through line 1863. I read the plus-pressure block 8346–8361, the plus-normal end 9023–9024, the minus-pressure block 9672–9687, and the minus-normal end 10349–10350. Formal statuses on both pressure faces are the even and time jets e_W, e_W_d1d1, e_W_d2d2, e_W_d3d3, and e_W_t, each as J and direct: 20 addresses, grades 00 and 11. Odd spatial jets are EXACT_ZERO_SOURCE_JET. The normal-slot ends I read are EXACT_ZERO_CONSUMER. I did not read every normal-slot line or all 544 inventory rows.

`preflight.py` `validate_selection` requires 2652 inventory rows, 544 e_W rows, and the status counts 102, 106, and 336. The 102 formal addresses are the full selected formal set. The 20 are the formal rows inside the 64-entry J/direct index. Those are source and metadata requirements. I did not execute them, and they are not a new runtime census.

The four contracted primitives in `method.md` lines 126–139 match `contraction-method.md` lines 105–139. Cj = a μ² W L/4 and Cd = (−i μ) W L/(4 i) match `kernel_components` in `inner_lib.py` lines 53–60, with a = 100/109 + 30i/109, μ = 3/10, and β = 30/109 + 9i/109. q(p) = sqrt(κ² − p²) on the positive-real / positive-imaginary sheet matches `inner_lib.py` lines 101–103 and `contraction/original-outgoing-branch.json`. J, Dh, and Dq use m = k+t and clip input k. Dr uses m = l−t, absolute Jacobian 1, and clips output l. q(l−t) stays q(l−t). The four primitives stay separate. E(k) and E(l) are already inside x and y.

Normalized-family sharing is allowed only after exact template and unit identity, inside one route, precision, carrier, window, and request. α = b c P_j i^n stays in the product. A24, A48, and B50 do not share evaluated values. Face labels do not authorize sharing Jwhole_plus with a minus template.

## Error accounting

The product rule |a|e_b + |b|e_a + e_a e_b, applied before cancellation and including Cj, Cd, resolvents, powers of m, α, weights, and Jacobians, is the right empirical propagation. ε = 1e−11/80 per eligible address and primitive, with no donation from zeros and with the tightest budget of every address in a shared family, matches the stated bookkeeping. The comparison tolerance τ = 1e−9 + 1e−7 |A48| is a separate gate.

The lemma is correct in the fixed momentum unit. κ > 2 and M > κ give ∫_{−M}^{M} dm/|q(m)| = π + 2 acosh(M/κ) < M + 4, so ∫ w < W_M = 3M + 4 with w = 1 + 1/|q|. The proposal already says this continuous integral does not certify a discrete sum of w, and that weighted inner totals must be ≤ ε/4. The build must nondimensionalize w in that momentum unit. These indicators are empirical estimates. They are not a rigorous total-error enclosure, a flux uncertainty, or a loss.

The base positive allocations are in `input/runtime-input/preflight/tail-plan.json`: K = 27, U = 75, T = 122, budget 1/10^11, `noCancellationUsed` true, `quadratureReady` false. Address 8346 has a positive rational outer term and a positive rational middle term, with heightQ equal to 0. Address 8347 begins with the same positive outer term. `Fourier-envelope-constants.json` and the compact Fubini majorants are not substitutes. I found no `tail-K29.json`. The (29, 124) window still needs its own positive no-cancellation allocation, saved or derived by the same rule, before the enlargement comparison. Reordering the identical (27, 122) window adds no new tail. The analytic constant-field identity removes the x-radius truncation of that transform; the outer k/l and middle t tails remain.

## Controls, storage, and what this review did not establish

The wrong-root mutant uses m = k+t, clipped input, unclipped output, and Dr_wrong = Cd ∫ {2 X1T Y1C + (X2T − m X1T) Y0C} dm (`method.md` lines 284–288; `contraction-method.md` lines 210–216). `contraction/new-wrong-root-mutant-factorization-return.json` contains only the inherited residual 0. The smallest formal direct id in the index I read is 8347. The derivative mutant replaces (i k)² by (i p0)² on the smallest eligible n1 = 2 address before the same J assembly; that index id is 8350, the plus J jet e_W_d1d1. For p0 = 0 that mutant factor is zero while the baseline need not be. The proposal requires measured movement above ten times the summed empirical envelopes and states no numerical outcome. Silence is unestablished coverage. H, normal-slot q(k) versus q(l), and the nonconstant Leibniz control stay deferred. The preflight selectors for those last two controls require nonconstant fields, which this constant subset does not contain.

`input/storage/preparation.json` is a synthetic prototype: 37 tests, exit code 0, scientific payload unrestored, runtime not integrated. The proposed journal — separate route namespaces, exact MP tuples, commit input before evaluation, full failure prefix, one inner panel of at most 48 nodes, optional 64 MiB route-local cache, 20 GiB disk reserve, and the existing 4 GiB worker inside the 16 GiB pool — matches that tooling boundary and the standing guard. Total adaptive work and storage remain unknown. A future concrete build still needs its own review. The SQLite tests do not clear this method.

## Mandatory joins once the sectors are corrected

- Centered w for every actual n, then exact equality with the displayed X and Y before quadrature.
- Full numerator and denominator, reconstruction, jet, template, and unit for every eligible address. Degree ≥ 1 and exact zeros stay outside the exception.
- Square-map Jacobian from `inner_lib.py` `gauss`: n ∈ (−1, 1), z = (n+1)/2, stored Jacobian 2·length·z, sum divided by 2. Emit the Jacobian before any guard. Restore the saved MP tuples. In `rules/A-GL24.json` the positive and negative nodes are not bitwise mirrors. `rules/extraction.json` names A-GL48 and B-G7-K15; I did not open those two bodies.
- Original sinhc expansion in `inner_lib.py` `profile`: |5π t| ≤ 1e−6, series through z^10/39916800, stated bound ≤ 2|z|^12/13!.
- Per-address tails from `tail-plan.json`, plus the enlarged-window allocation. No cancellation against those budgets.
- Family sharing only after complete J or D template equality, including face-specific whole tags. I read the plus mixed template on address 8346. I did not open the minus template body.
- Control addresses by the stated minimum-id rule after the eligibility join. The index values are 8347 and 8350. One direct address does not move a different minus-face template.
- Request identity by full operands. A completed looser error target is not a tighter one. Unmatched old bank points are not misses.

## Coverage limits

I read `method.md`, `evidence-guide.md`, `contraction-method.md`, `packet-action-method.md` through its section 5, the Fourier library through `request`, the inner library through the start of `H`, all of `preflight.py`, the index samples above, `physical-plan.json`, `fields.json` through the first nonconstant entry, the two polynomials named above, `tail-plan.json` at its header, at addresses 8346 and 8347, and at its closing budget, `A-GL24.json`, `extraction.json`, `original-outgoing-branch.json`, the wrong-root residual, and `storage/preparation.json`. I did not execute a script or a CAS, and I did not read a peer report.

The 1e−12 check is absolute, so it is weak in the far Gaussian tail, and a bug shared by both formula routes would pass it. The two mutants test derivative routing and the reflected-root domain. They do not prove the Gaussian theorem. The closed-quadrant bound |q1/(q1+q2)| ≤ √2 does not assign a value at simultaneous zeros and does not give pointwise smoothness. Inferred gamma units and the accepted algebra remain inherited dependencies; `preflight.py` records that this packet adds no new post-binding dimension proof. This assessment is not a numerical value, a runtime forecast, an author clearance, or acceptance of the full packet.