# Assessment of the closed direct mixed-kernel grazing method (S11c)

**(A) Source/evidence applicability: SUPPORTED WITH STATED LIMITS.**
**(B) Proposed method: CLEAR FOR BOUNDED CLOSED-GRAZING IMPLEMENTATION.** I found no blocker in the mathematics. The clearance holds only if the four acceptance amendments in item 6 go into the instrument specification. It clears the method, not any code, since none exists.

I checked the algebra by hand. I did not inspect the full THETA/E_W expanded rows, the three native source scripts, or the full 46 KB closed-increment expressions beyond their top-level keys and a few leading terms. I ran nothing, and everything below is source and algebra only.

## 1. Density, factors and rearrangement (1)

- **Factorization.** `raw/physical-factorization-return.json` gives C = tB. I checked it using only qⱼ² = κ² − pⱼ². The identities are t(t+2k) = qᵢ² − qₕ² and t(2l−t) = qₛ² − qₒ². Multiplying out B reproduces the numerator in `physical-factorization-input.json` exactly, with H = t. The identity is purely algebraic, so it also holds at complex Ω.
- **Rearrangement (1).** Multiplying B termwise by qᵢqₒ gives method.md (1) exactly. The three terms become k(2l−t)/(qₛ+qₒ), k(t+2k)qᵢ/[qₕ(qₕ+qᵢ)] and qᵢ²/qₕ.
- **External factor and iteration.** `raw/plus-closure-operands.json` has diagonal Z = 3/(10q), so β = aρω. The (0,2) closed entry splits as R·qᵢqₒ/((qᵢ+β)(qₒ+β)) plus a z01·z12 term. The z01·z12 term is the first-shape iteration, and `plus-iteration-unchanged-return.json` certifies the cancellation of the R part. The iteration is left once.
- **Limit of this evidence.** These certificates cover only the direct-R part. They do not check the iteration term's own closure or its grazing behavior.
- **Density and normalization.** The prefactor (WL/4i)·A·A matches `raw-ordered-before-cancel.json`. Its height numerator is 5t/(2 sinh 5πt) = A(t) and its jet is 5·A(u). No stray 2π appears.

Two conventions need to be pinned in any join:
- **Face sign.** G in method.md §2 carries no face sign. The face sign lives in the trace and jet entries: `plus-trace-domain.json` has +ηw/2 and +iqₒ, and `minus-trace-domain.json` has −ηw/2 and −iqₒ. I checked the minus-face jet against `minus-closed-raw-increment.json` and it equals −iqₒ × physical.
- **Mirror evidence.** `lower-upper-outward-mirror-return.json` is a residual between identical text in outward variables. It is not an independent derivation of the lower face. The real evidence is `lower-boundary-operands.json` and `native-lower-geometry-operands.json` (face = −1, outward slope σw′/2).

## 2. Continuation Ω = 3 + iδ

- **Beta formulas.** With x = 1+τδ and y = 3τ, I get Re β = 0.1·(3x − δy)/D = 0.3/D and Im β = 0.1(δx + 3y)/D. These match method.md. The minimum of Re β over δ ≤ 1/10 is 0.3/1.1101 = 3000/11101, as claimed.
- **Quadrants.** q² has imaginary part 6δ/cs² > 0, so the principal root lies in the first quadrant. At δ = 0 it reduces to the saved sheet (positive root, or i·√ for negative radicand). For u, v in the first quadrant, |u+v| ≥ max(|u|,|v|).
- **Prescription.** The prescription is the native outgoing one: e^{−iωt} with a = Λ/(ρ²(1−iωτ)) is analytic in the upper half-plane. No conjugation appears.
- **What the continuation covers.** It is applied to G only. The grade-zero rows are q-free constants frozen at ω = 3. `plus-source-zero-grade-operands.json` contains no q-symbols, and `effective-speed-reuse-domain.json` says `sourceSymbolsCs: false`. The result is therefore an L1 limit of G with the rows held fixed, not an analytic continuation of the whole term.
- **Frequency override.** `physical-input.json` says `omega: "1"`. The frequency 3 comes only from `plus-trace-domain.json` (`savedFrequency`/`newFrequency`). The instrument should cite this binding.

## 3. Domain, envelope (2), and the collision cases

- **Constants.** The constants are correct but loose. |Ω| ≤ 3.002, not 4. |q|² ≤ 18.06, not 25. The (18+3|t|) and (43+3|t|) coefficients follow from |k| ≤ 3, |2l−t| ≤ 6+|t| and |qᵢ| ≤ 5.
- **Inverse-root bound.** |q|² ≥ |Re q²| = |κδ² − p²| ≥ κ_min·dist(p, ±κδ) is valid. κ_min² = 879/400 is correct.
- **Local integrability.** Since ∫_E |t−a|^{−1/2} ≤ 2√(2m), the integrals over moving endpoints are uniformly integrable.
- **No product singularities.** (1) contains only separate 1/|qₛ| and 1/|qₕ|, never 1/(qₕqₛ). Collisions therefore add; they do not multiply. This holds for qₕ = qₛ = 0 at t = 0 when k = l = κ.
- **Case l = −k.** Here qₛ ≡ qₕ as functions, because q is even in p. The method's remark that "qₛ is not qₕ" fails in this case. The bound still holds, but the collision list must include it.
- **κ = 0 exclusion.** κ² = 5.95 (LEFT) and 6.01 (RIGHT); I verified both. κ = 0 needs cs² = 180, far outside [1,2]. The exclusion is sufficient and irrelevant to the actual matches. The endpoints −k±κ stay 2κ apart, so the simple-root bound applies.
- **Tails.** Tails are fine: A·A decays like e^{−5π(|t|+|Q−t|)}, against polynomial growth of at most (43+3|t|)/|qₕ|.
- **Limit existence.** Pointwise a.e. convergence plus uniform integrability plus tightness gives the L1 limit (Vitali). Equations (3) and (4) follow by direct substitution.

## 4. L1 limit, contact, and what is not claimed

- **Approaches.** The L1 limit covers radiating and evanescent approaches, δ → 0, and the joint δ/qᵢ/qₒ → 0 limits. Away from the finite endpoint set, convergence is by continuity of q into the closed first quadrant.
- **Contact.** For δ > 0 all q are nonzero, so B is finite at t = 0 and tδ(t)B = 0. Every G_δ is then a pure L1 function, and an L1 limit cannot gain a delta. A simpler argument also works at δ = 0. Near t = 0 (qᵢ = 0) we have |tBc| ≲ |t|^{1/2} → 0, so δ(t)·tBc = 0 directly.
- **Gap.** The saved contact-zero evidence (`upper-complete-contact-*`) is for real ω and uses symbols declared `omega, real=True`. The complex-Ω contact zero is new and needs its own certificate.
- **Not established:**
  - pointwise regularity;
  - differentiability in cs (the parameter-dependent endpoint singularities are only |t−a|^{−1/2});
  - continuity in (k,l) beyond L1 continuity at fixed test function 1;
  - any rate of convergence;
  - convergence of the iteration term.

## 5. Reference/jet and row/source reuse

- **Reference and jet.** The reference trace is 1 and the jet is ±iqₒ, consistent with the saved records. The jet is bounded and continuous in the domain.
- **Saved evidence limits.** `finish/applicability.json` states `exactMatch: UNRESOLVED`, `rawNongrazingOnly: true`. `minus-closed-raw-increment.json` states `uniformGrazingLimitEstablished: false`. The saved evidence is nongrazing only, and this is exactly what the proposed method must supply.
- **Execution guards.**
  - Join the full normal jet, not just iqₒG. `nativeNormal` includes slope components (−σw₁_dⱼ/2), and `reuse-native-slope-omission-response.json` shows the omission matters. The (1,1) direct term's jet is iqₒG, but the slope-times-reference pieces belong to the iteration.
  - Carry the face sign (item 1).
  - The rows are valid as ω = 3 constants for any cs and are bounded in k, l. Their finiteness at the exact match points should be an explicit certificate. Do not assume a finite-inverse theorem; `source/finite.py` and the guide warn against this.
- **Excluded.** Do not infer an incoming transverse excitation or finite-deficit protection.

## 6. Acceptance evidence and controls

**Adequate as designed:** the per-face join requirement, the β sign and bound certificate, root-quadrant certificates, the envelope with an explicit tail, the qᵢ-only/qₒ-only/simultaneous limits, the contact argument, and the ban on sampled grazing.

**Required amendments (execution guards, not changes to the proof):**

1. **Complex-Ω identity.** Re-derive C = tB and C(0) = 0 with unrestricted ω. Do not reuse the `real=True` saved residual.
2. **Collision list.** Include l = −k (qₛ ≡ qₕ), l = k, and qₕ = qₛ = 0 at t = 0.
3. **Face sign.** Join the ± height and ± jet sign to G for each face, from `lower-boundary-operands.json` and `native-lower-geometry-operands.json`. The outward-mirror equality alone is not enough.
4. **Predeclared control responses.**
   - **Missing external factor:** replacing Bc by B must show an infinite or nonintegrable density at qᵢ = 0.
   - **Wrong root sheet:** this will not necessarily break finiteness, since |q+β| ≥ Re β only needs Re q ≥ 0. The control must instead fail the quadrant certificate or change the value of G. Define its expected response as such.
   - **Wrong lower jet sign:** this must change the lower-face jet join residual to nonzero.

   Existing `finish/reuse-*` responses are sensitivity of the source-row reuse, not controls for the new density.

**Optional:** tighten |Ω| and |q| constants; replace the δ-device with the direct t·Bc → 0 argument.

## Strongest supportable claim

At real ω = 3 with the native a(3) and β = (30+9i)/109, in the declared compact domain (cs ∈ [1,2], |k|,|l| ≤ 3, κ > 0), including both selected matches and all combinations of qᵢ and qₒ grazing (qᵢ = 0, qₒ = 0, or both, i.e. l = ±k), the closed direct reference density G, joined face by face, is a well-defined a.e. function in L1(ℝ). It is the L1 limit along any parameter path of the nongrazing and δ-regularized densities. It has no concentrated contact term, and its t-integral is finite and continuous in those parameters. This claim holds only after the exact per-face joins.

**Excluded:**
- the first-shape iteration term and the full retained operator;
- the untruncated finite inverse, current balance, drain and loss;
- pure η² and σ_W² blocks;
- pointwise regularity, differentiability in cs, and any loss smoothness through speed matching;
- primitive speed calibration and κ = 0;
- an impermeable β = 0 limit;
- any incoming-excitation or finite-deficit statement;
- the defect sweep.
