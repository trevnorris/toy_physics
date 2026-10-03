**Verdict: NEEDS REVISION**

I found one substantive blocker, in the Route A Fourier setting. Everything else I checked is sound or an optional improvement. The fix is a change to the declared settings in the method text. It does not need a larger redesign.

## Substantive blocker

**B1. Route A Fourier orders 96/160 cannot meet the method's own 1e-12 transform criterion, and a miss stops the action with no retry.** The relevant text is method.md §4 (Fourier routes, order 96/160, "A failed transform criterion stops the action; no automatic extra orders").

- **Resolution limit.**
  - With weight exp(-(x-x0)²/(2s²)) and t=(x-x0)/(s√2), the kernel e^{-iνx} becomes e^{-iωt} with ω=√2·s·ν≈11.3ν for s=8.
  - An n-point Gauss–Hermite rule resolves only ω≲√(2n). That is about 13.9 for n=96 (|ν|≲1.2) and about 17.9 for n=160 (|ν|≲1.6).
- **Error size.**
  - Beyond that limit the rule error is of the order of the shifted-integrand scale, e^{-5|ν|}·e^{25/128}·(moment constants). The analytic shift to |Im x|=5 only supplies that factor.
  - e^{-5|ν|} drops below 1e-12 only near |ν|≈5.5. So for roughly 1.2<|ν|<5.5 the A-96, A-160 and Route B values will disagree by far more than 1e-12. At |ν|=2 the shift factor alone is about 5e-5.
  - The true transform there is tiny (about e^{-5π|ν|}), so the error is pure quadrature error.
  - The tail rule (K chosen so e^{-5K} times the polynomial constants is below 1e-11) forces requests at |ν| of about 6 or more. The failing band is therefore inside the requested range.
- **Why this is a method problem.** The text says contour damping "does not prove" resolution, but it still fixes orders 96/160 and forbids more. The criterion would then fail predictably and halt the packet.
- **Correction.**
  - Replace the fixed 96/160 with a declared resolution rule. For example, require n ≥ ((1.3·√2·s·ν_*)²)/2. Here ν_* is the offset where the saved envelope C·e^{-c|ν|} falls below 1e-13, taken at the final c and the largest requested |ν|. Declare base n and refined 1.5n.
  - Alternatively, use a composite Gauss–Legendre rule on the shifted line, with panel width fixed by points per wavelength.
  - With c=5 the rule needs n of order 10³. A larger shift up to about 10–12 stays inside the tanh poles at 5π and reduces ν_*. Check |tanh(x/10)| on the larger strip before using it.
  - Write the rule into the method before build review.

## Checked and sound (no blocker)

- **Scope, speed and counts.**
  - cs=√6/2 and κ=√595/10 match `ends/left-match.json` (`weak_end_cs`, `weak_end_p`).
  - Ω²/cs²=6 and κ²=6−1/20=5.95, so q(p)=√(6−1/20−p²) is consistent.
  - The 544 addresses, with 102 formal and 442 zero, and the 2,652 total (13,260/5) are consistent with the address header and `local/assembly.json`.
  - The e_W jets in `pressure/source-jets.json` include time, d2 and d3 orders, which the D_j factor covers.
- **Grade weights.** `ends/left-source-binding.json`, key `origin`, gives eta_bg=1/100 and sigma_W=1/1000. The display B00+B10/100+B01/1000+B11/1e5 follows from that origin and is not a new calibration.
- **PV and contact.** The symmetric half-line form is equivalent to the saved χ subtraction. A is even, so the f(0)χ term is odd and integrates to zero. The Q cuts at |±κ−k| are the complete union of the f(Q) and f(−Q) roots.
- **H product reference.**
  - I derived the physical profiles from the saved h(s) and j(s=LA/2).
  - With A=5s/(2sinh 5πs) (`profile-envelope.json`) and the standard sine and cosine integrals, h(x)=W(1+tanh(x/L))/4 and j(x)=sech²(x/L)/4.
  - Under the 1/(2π) convention, hat(hj)=hat h * hat j with no extra factor. The proposed candidates and the join are correct.
- **Collision geometry.**
  - The middle roots are t=±κ−k, from qm or qh. The reflected route qs gives t=l∓κ. The profile points are t=0 and t=l−k.
  - All pairwise coincidences give exactly k=±κ, l=±κ, k+l=0, k+l=±2κ and l=k. The method's line set is complete.
  - T≥K+4 exceeds every |±κ−k| and |l±κ|, which are at most K+κ.
  - The four test pairs lie on k+l=0, k+l=2κ, k+l=−2κ and l=k, and avoid the external grazing lines.
- **Actual summed density.** `pressure/typed-direct.json` `closedDensity.Bc` is k(2l−t)/(qo+qs) + qi²/qh + qi·k(t+2k)/(qh(qh+qi)). It is a sum with no 1/(qh·qs) product.
  - Same-quadrant roots (real or positive imaginary) make qs+qo and qh+qi vanish only when both vanish.
  - |qi/(qh+qi)|≤1 keeps each term dominated by 1/|qh| or 1/|qs|, which is integrable.
  - So no logarithmic pinch is needed. The method correctly declines to claim differentiability.
- **Typed objects.** `pressure/typed-direct.json` (`barePlaceholderIsWholeTag: false`, `multiplyWholeTagByResolvents: false`) and `pressure/whole-definitions.json` (`wholeDirectOnce`, `nativeIterationOnce`) support the distinct-object and once-only statements.
- **Fourier strip.** The shift gives exp(25/(2s²))·exp(−5|ν|), and 5/10<π/4 gives |tanh|≤1. The carrier phase and the moment translation are as stated.
- **Honesty rules.** No-flux-error wording, empirical envelopes, no retry, and the e_W→THETA scope limit are properly stated.

## Optional improvements (not blockers)

1. **Route B singular treatment.** "Explicit integrable branch/PV subtraction" is not concrete. Plain adaptive Gauss–Kronrod converges poorly at inverse-square-root endpoints. Specify the transformation or subtraction before build review.
2. **Middle-integral tails.**
   - Cuts stop at ±8/L=0.8 from the profile points, so one panel runs out to T, which may be about 40 or more with the loose e^{-|t|} envelope.
   - On such a panel 24/48 Gauss nodes under-resolve the true e^{-5π|t|} tail. The tail is only about 1e-6 of the peak, so I did not rate it as a blocker, but a miss stops the run.
   - Add doubling cuts (16/L, 32/L, …) out to T.
3. **Absolute 1e-12 transform criterion.** It passes trivially for tiny large-|ν| values, and it gives no propagated bound on the action, since kernel growth P³ and the area (2K)² amplify errors. The method concedes this, and the complete-action comparison covers it. Record the large-|ν| requests, as already planned.
4. **Inner-comparison coverage.** Add the actual carrier points. For p0=κ the packet peak (κ,κ) is a four-line intersection. For p0=0 the peak (0,0) lies on l=k and k+l=0. Include them in the pre-selected inner comparisons.
5. **Missing control.** Add a control that removes the native height contact. Only the H contact removal is listed.
6. **Wording.** "Polynomial-Gaussian" Fourier products should say polynomial-in-tanh times polynomial times Gaussian. The carrier derivative factors (i·p0)ⁿ should be explicit.

## Does a missing item need a larger method change?

No. B1 is a revision of the declared numerical settings in the method document. The exact adapter joins, tail constants and numerical controls remain runtime obligations for the concrete worker, and I see no gap that needs more than tooling. Nothing here clears a future build or the old finite benchmark.

## Files inspected

- **Read fully:**
  - `input/method.md`
  - `input/evidence-guide.md`
  - `input/packet-index.json`
  - `input/pressure/typed-direct.json`
  - `input/pressure/whole-definitions.json`
  - `input/pressure/Fourier-order.json`
  - `input/pressure/source-jets.json`
  - `input/pressure/profile-envelope.json`
  - `input/local/assembly.json`
  - `input/ends/left-match.json`
  - `input/background/pressure-weak-method.md`
  - `input/background/numerical-route-inventory.json`
- **Partial:**
  - `input/ends/left-source-binding.json`: only the `origin` block, via grep. The full-file read was truncated and I did not use it.
  - `input/selected/pressure-addresses.json`: only the header counts, via grep. I did not check the individual addresses or zero statuses.
- **Not opened:**
  - `input/review-prompt.md`
  - `input/selected/local-cells.json`
  - `input/pressure/fields.json`
  - `input/pressure/whole-tags.json`
  - `input/pressure/coefficient-certificates.json`
  - `input/pressure/global-parameter-domain.json`
  - `input/pressure/whole-envelopes.json`
  - `input/pressure/normal-growth.json`
  - `input/pressure/H-bound.json`
  - `input/pressure/J-numerator-envelope.json`
  - `input/pressure/D-height-envelope.json`
  - `input/pressure/D-reflected-envelope.json`
  - `input/pressure/shift-root-bound.json`
  - `input/pressure/plus-height-PV.json`
  - `input/pressure/minus-height-PV.json`
  - `input/local/context.json`
  - `input/physical-input.json`
  - `input/ends/correspondence-summary.json`
  - `input/background/full-weak-method.md`
  - all of `input/source/*.py`

I did not verify the tanh-polynomial claims, the 34-field certificates, the bound constants (KJ, KD, Cq) or the 16 local cells. I relied on the background descriptions for those.