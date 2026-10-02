**Verdict: CLEAR FOR THIS TRANSLATED WEAK-END METHOD**

I found no blocker. This clears the method only for translated weak end forms and their formal diagonal symbols. It does not clear a worker, a numerical result, or any of the out-of-scope items listed at the end.

**Files examined**
- Packet guides and method: `evidence-guide.md`, `method.md`, `reference/method.md`.
- Pressure certificates: `pressure/analytic-conclusion.json`, `global-parameter-domain.json`, `plus-global-height-PV-certificate.json`, `global-whole-kernel-envelopes.json`.
- Full-operator records: `full/complete-weak-assembly.json` (first 150 lines).
- Inventory records: `inventory/fourier-and-unit-provenance.json` (first 150 lines).
- Pressure records: `pressure/new-weak-duality-and-order.json` (first 100 lines).
- Reference records: `reference/plus-native-height-constant.json`, `minus-native-height-constant.json`, `minus-new-native-trace.json`.
- Targeted greps: `profile-jet-certificates.json` (names only), `physical-input.json`, `source/reference-worker.py` (lines around 164–227).

**What I checked and found sound**
- **Phase.** The saved profile phase is exp(-i(ko-ki)y). Translating the packets gives hat v_a(-l)·hat u_a(k) = e^{i(l-k)a}·hat v(-l)·hat u(k). That agrees with the method's exp(i(l-k)a) and with the saved Y(l)=2π·hat[cv](-l), X(k)=hat[bDu](k) pairing.
- **Height endpoints and Dirichlet limit.**
  - The Fourier convention requires h = (W/2)[δ/2 + PV A/(is)], and the saved profile has A(0)=1/2π.
  - Dirichlet gives PV∫e^{iQa}A·χ/Q → (i/2)sgn(a), so the PV term contributes ±W/4.
  - The contact contributes +W/4, so the total is 0 at a→−∞ and W/2 at a→+∞.
  - This matches the saved native height eta·w1/2, with w=(1+tanh)/2 giving endpoints 0 and 1/2. Removing the contact would give ±W/4, so the method's sensitivity claim is correct.
  - The reference worker maps w1 to 2h/W, which is consistent with H_+=W/2 at W=1. The instrument must still join this actual map, as the method says.
- **Height decomposition (1).** It correctly adds the oscillatory diagonal term that the unshifted subtraction formula lacks. The middle term is L1 by the saved Hölder-1/2 bound with the Pk^{3/2} weight. The Hölder constants also cover the normal coefficient, which is bounded directly without assuming qφ is C¹.
- **Both faces.**
  - The lower face has height −eta·w1/2 and normal jet −iq, so the trace product is +i·eta·q·h on both faces.
  - Equation (2) is consistent with F0 plus the reference height coefficient Fh = −iμq/(q+β)·h.
- **Uniform-in-a bounds.** The method avoids differentiating e^{iQa} and uses the sine-integral bound plus an L1 remainder. That is correct. X_a and Y_a converge in S, so their seminorms are uniformly bounded.
- **Ordinary mixed kernels.**
  - The kernels have the saved degree ≤3 polynomial envelope. A Schwartz-weighted pairing is therefore L1, and Riemann–Lebesgue applies to the linear phase.
  - Swapping X_a,Y_a for the limits costs only seminorm differences times a fixed L1 bound.
- **Uniformity over cs.**
  - The saved envelope certificate states that continuity comes from moving-endpoint uniform absolute continuity, compact tails and a.e. convergence. It explicitly does not assume a single pointwise dominant.
  - That is the right ingredient for Vitali-type L1 continuity, hence compactness of the kernel family and a uniform Riemann–Lebesgue limit. The method's fallback (restrict to fixed cs if this is not supplied) is explicit and correct.
- **Constants.** The β separation and the kappa, mu and q bounds in `global-parameter-domain.json` match the method's cited values. δ=0 is included in the domain.

**Nonblocking limitations**
1. Uniformity over cs rests on an inherited analytic certificate that its own text labels "not machine measure theory". If the later worker finds the a.e./uniform-integrability hypothesis unjoined, the result should be stated at fixed cs. It should not be silently widened.
2. The statement "Zh=0 on the diagonal" is pointwise only on the nongrazing domain. The method correctly avoids raw 0/0 at q=0 by using the closed reference form. The instrument should derive the end multiplier from the reference coefficient (Fh), not from the physical Zh.
3. The δ→0⁺ limit and the a→±∞ limit are taken on the already-defined limiting distribution. The method does not assert uniformity in δ jointly with a, and does not need it.
4. W appears in the native source as 2h/W and also as W_0 in the d_w coefficients. At W=1 these coincide. The instrument must join the actual map and not rely on W=1.
5. The lower-face endpoint, normal-sign and contact controls are formal sensitivity tests only. They cannot confirm the analytic argument, and the method says so.
6. Coverage of the 13,260 addresses, the 16 grade triples and the full source/consumer endpoint joins is deferred to the instrument.

**Material coverage gaps in my review**
- I did not read the 400-cell local polynomials, the address representatives, `all-local-cells.json`, the source workers beyond the greps above, `c1.py`, `c2.py`, or the retained-response census.
- I did not independently verify that the saved local endpoint values, address routing, epsilon normalization, Fourier 2π factor or source/consumer cross-grade endpoint values are correct. I only confirmed that the method requires those joins and that the instrument is positioned to test them.
- I read the complete-assembly, provenance and duality files only up to 150, 150 and 100 lines respectively.

**Scope boundary**
The clearance covers translated weak end forms and formal diagonal symbols E_±,g. It does not cover:
- old end-mode acceptance, roots or a mode census;
- plane-wave scattering, decay rate, a finite boundary or inverse, or loss;
- calibration or a draining model.

The method states these exclusions itself, and I found no inference that crosses them.