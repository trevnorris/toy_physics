I'm drafting the verdict now. I've read the method, the evidence guide, the background weak-method file, the left-match record, the typed-direct and height-PV certificates, and parts of the packet index and local cells.

# Verdict: `NEEDS REVISION`

The method is mostly sound. I found one substantive blocker, which is a method-text amendment and not a tooling fix. I also list some non-blocking items.

## Files inspected
- **Read in full:** `input/method.md`, `input/evidence-guide.md`, `background/pressure-weak-method.md`, `ends/left-match.json`, `pressure/plus-height-PV.json`, `pressure/whole-definitions.json`, `pressure/typed-direct.json`.
- **Read in part:** `packet-index.json` (first 80 lines) and `selected/local-cells.json` (first 80 lines).
- **Grepped only:** `physical-input.json` and `local/assembly.json`.
- **Not inspected:**
  - `selected/pressure-addresses.json` (so the 544, 102, 106 and 336 counts are checked only for arithmetic: 102+106+336=544).
  - `pressure/fields.json` and `coefficient-certificates.json`.
  - The other `pressure/*` envelope files and `local/context.json`.
  - `full-weak-method.md` and the `source/*.py` files.
  - `background/numerical-route-inventory.json` and `review-prompt.md`.
- **Limit on the grade weights:** I confirmed eta=1/100 in `physical-input.json`. I did not find sigma=1/1000 there, so the 1/1000 and 1/100000 display weights are unverified.

## Checks that passed
- **Matching speed:** cs=√6/2 gives 9/cs²=6, so κ²=6−1/20=595/100. This matches `left-match.json` (sqrt(595)/10, sqrt(6)/2).
- **Height rule:** the positive half-line form `(W/4)f(0) + W/(2i)∫₀^∞ A[f(Q)−f(−Q)]/Q` follows from the inherited form (7). The χ terms cancel because A is even. The required Q-panel splits at |±κ−k| are correct.
- **Profiles:** h(x)=W(1+tanh(x/L))/4 and j(x)=sech²(x/L)/4 reproduce h(s)=(W/2)[δ/2+A/(is)] and j=LA/2 under the 1/(2π) convention. Both the contact W/4 and the PV prefactor match. The product-to-convolution factor is 1 in this convention. Method §3 correctly demands this as a runtime join.
- **Strip:** |tanh(x/10)|≤1 holds for |Im x|≤5, since 0.5<π/4. The nearest poles are at 5πi. The shifted Fourier factor is exp(25/(2s²)−5|ν|), and the carrier cancels in the offset ν.
- **J/GD endpoints:** the splits t=±κ−k and t=l±κ, with removable points t=0 and t=l−k, match `typed-direct.json` (qh=q(k+t), qs=q(l−t), sinh(5πt) and sinh(5π(l−k−t))).
- **Typed distinctness and scope:** `typed-direct.json` states that Rprod=qi·qo·E, that the closed density is not Rprod, and that the whole tag is not multiplied by resolvents. This matches the method's "once" rules. The scope is stated correctly as weak actions that are not power, loss or flux.

## Substantive blocker

**B1. The outer 2D (k,l) quadrature omits the diagonal non-smooth lines created by internal endpoint collisions.**

- **Evidence:** §3 splits the outer integrations only at k,l=±κ and names "k=l and k=−l" as collisions to record. The J/GD middle endpoints are t₁=s₁κ−k (from qm=q(k+t)) and t₂=l+s₂κ (from qs=q(l−t)). They coincide when k+l=(s₁−s₂)κ, which is k+l ∈ {0, +2κ, −2κ}.
  - At each of these lines two inverse-square-root endpoints pinch. For J this involves 1/qm and the qm+qi ratio. For GD it involves 1/qh and 1/(qs+qo). The result is a log-type, or at least non-Hölder-smooth, dependence of the inner whole integral on (k,l).
  - Lines of this kind are not axis-parallel to the ±κ splits.
  - The line k=l makes the two removable profile points t=0 and t=l−k coincide. That is a removable-point coincidence and not a branch singularity. The plan to run Gauss-Legendre 24/48 tensor panels with cuts only at ±κ, carrier ±1/s offsets and so on does not resolve the k+l=0, ±2κ lines. A tensor Gauss rule converges only algebraically across them. The base-versus-refined agreement, the 1e-9+1e-7 tolerance and the enlarged-box test could then fail for reasons unrelated to the physics, or could pass only by luck.
  - The "independent outer adaptive integration" is not described well enough to be an independent reference. A 2D adaptive scheme can also stall on diagonal log lines unless it is told where they are.
- **Correction:**
  1. Change the outer variables to (k, σ=k+l) or (k, Q=l−k), or add explicit diagonal cut lines.
  2. List the cut set for the outer integrals as k,l=±κ, k+l ∈ {0, ±2κ}, and l=k. Add the profile-induced lines too: A(l−k−t) singular points are at t=l−k, which is already covered by the l=±κ intersections.
  3. State the integrable-singularity treatment at those lines (for example, a square or logarithmic substitution on each side, with the inner integral evaluated at the line from both sides).
  4. Make the route-B outer rule use the same enumerated singular set but its own subdivision.
  5. Add a pre-declared outer test at nodes on and near k=−l, k+l=±2κ and k=l. The existing inner point tests list only "k=±l".

This is a modest method-text amendment and not a redesign. It does need to be settled before a worker is built, because it determines the outer panel structure and the empirical-envelope meaning.

## Non-blocking improvements
1. **Route A resolution:** Gauss-Hermite at orders 96/160 with s=8 has node spacing of about 2 in x near the centre. It cannot resolve e^{−iνx} at the large ν needed for K≈15–20, and the tanh factor makes it non-exact. The shift only lowers the integrand scale to e^{−5|ν|}. The quadrature error may therefore dominate the true (Gaussian-suppressed) transform at large |ν|. A trapezoid rule on the shifted line is exponentially accurate within the strip, and I recommend it as an alternative. At minimum the build should state the GH-versus-route-B transform tolerance as a stopping criterion at large |ν|.
2. **Grade weights:** the sigma weight 1/1000 should be joined to the saved input, since I did not find it.
3. **Middle truncation:** the claimed T≥K+4 with the loose envelope exp(−|t|) from (3) gives e^{−K−4} only after polynomial factors. This will probably force T in the 30s. The method already says to report infeasible limits as a method limitation, which is acceptable.
4. **Controls:** the controls are applied as "applicable-address-by-metadata", which is good. The coverage report should say explicitly that the carrier p0=0 setting makes some Leibniz controls silent when b is slowly varying, and that silence must stay "unestablished".

## Items that are correct and need no change
- The no-conjugation bilinear pairing, the epsilon-once rule, the independent grade triples and the derivative-before-coefficient order.
- The exclusion of a flux error bound and of a finite-delta extrapolation.
- The empirical-envelope language, and the separate treatment of the analytic tail bound and the numerical resolution error.

## Does any missing item need a larger change?
Only B1 is a method change, and it is a bounded one. The Gauss-Hermite resolution concern, the exact adapter joins (h/j transform, product factor, tail constants C_(field,jet)) and the numerical controls are runtime obligations of the concrete worker. They are not method defects.

This review does not clear a future build or the old finite benchmark.