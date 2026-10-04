CLEAR FOR THIS BOUNDED NUMERICAL J/D PACKET METHOD

I found no fatal mathematical error. The items below are mandatory obligations for the concrete build, not method defects. I did no execution, so the Gaussian identity and the lemma rest on hand derivation and source reading.

## What I read
- **Method and guide:** `input/method.md` and `input/evidence-guide.md`, in full.
- **Fourier and contraction sources:** `original-source/S11c_d_defect_packet_fourier_lib.py`, `_contraction.py`, `_contraction_lib.py` and `_evidence_store.py`, in full.
- **Kernel and tail sources:** the `kernel_components`, `q`, `profile` and `middle` parts of `_inner_lib.py`, and `exponential_moment`, `weighted_tail`, `contributions` and the tail-derivation block of `_preflight.py` (lines 120–150 and 290–358).
- **Saved operands:** `tooling-tests.log`, the first ~120 lines of `constant-subset-source-index.json`, targeted greps of it, and the fragment of `contraction/address-9990-input.json` around lines 640–700.
- **Not read:** the other 543 addresses, the 34 field files, the tail JSONs, the rule files, `AGENTS.md`, the guard and supervisor.
- **Census:** I did not independently recount 544/64/20. The index shows 40 "formal nonzero" status lines, which fits 20 entries listed twice, and the visible entries have grades 00/00 and target 11. The 10 J, 10 direct and 5-jets-per-face split is consistent with that but not rechecked entry by entry.

## 1. Constant-field Fourier exception
- **Gaussian transforms:** I re-derived both formulas.
  - X: `∫u e^{-ikx}dx = e^{-i(k-p0)x_u}·s√(2π)e^{-s²ν²/2}`. The Q_n recurrence `(ip0 − w/s²)` matches the lib's `q` loop at `fourier_lib.py:173`.
  - Y: ν = p0−l and phase `e^{+i(l-p0)x_v}`.
  - Moments: M1 = −is²νM0, and M_{r+1} = s²(rM_{r−1} − iνM_r) by integration by parts.
  - Routes A and B are algebraically consistent. Route B's moment sum already holds the derivative, as the method says.
- **Original analytic reference:** the old `constant_reference` at `fourier_lib.py:241-253` uses the same structure. So the exception specializes an existing reference rather than adding a new convention.
- **Wave multiplier:** the index shows `waveMultiplier` = 1 for n1=0 entries and `-composition_p**2` for n1=2. That matches (ik)^n at k=p, with n=2 only for `e_W_d1d1`. `preflight.py:207` has the same `(-3i)^nt (ip)^n1 (i/5)^n2 (i/10)^n3` structure.
- **Normal factor:** N=1 for pressure slots, per `contraction.py:229`.
- **Checks:** the two checks genuinely differ.
  - The full-value 1e-24 check is vacuous in the Gaussian wings, because the value is about e^{-7e5} at |ν|≈150.
  - The Gaussian-stripped scaled check carries the real gating. At 30 digits, phase argument error is about 4e-28, well inside 1e-24.
  - Both checks share one analytic identity, as stated. So the only independent evidence for the Fourier transform itself would be old bank matches, which are unlikely to match new nodes. That is a real limit of coverage.
- **Build obligations:**
  - The index lacks the `argumentDerivative` selector and any Y-interface carrier/argument for the 20 entries. In `address-9990` the `matchingYFamilyInterfaces` list is empty and the X interface has no carrier. The worker must construct and save the X and Y specs, with Y carrier −p0 and argument −l.
  - The worker must prove the native factors (`P_j`, n2/n3) appear exactly once relative to the template and factor proofs.
  - The worker must verify |X_analytic| ≤ CX e^{-5|ν|}, so the inherited tail envelope applies to the analytic X.

## 2. Contractions
- I matched `contraction.py:187-197` to the method's formulas for J, Dh, Dq and Dr.
  - J: `(k+t)(2k+t) = m(k+m)`, which gives `C1T + m·C0T`.
  - Dr: `k(2l−t) = k(l+m)`.
- The `Cj` and `Cd` constants match `kernel_components`.
- The input clipping for J/Dh/Dq and the output clipping for Dr are correct.
- The profile `5t/(2 sinh 5πt)` equals `A` with L=10. The sinhc cutoff bound is as stated.
- Family sharing through `alpha = b·c·P_j·i^n` is exact linearity. The method's all-address budget rule is the right safeguard.
- Because only n=0 or 2 and the primitives are shared, the numerical work is small relative to 20 addresses. This is an observation, not a claim about cost.

## 3. Panel geometry
- **Slabs:** the six slabs from {−M, −95, −κ, 0, κ, 95, M} are correct, and I verified the stated d_g per slab.
- **Clipping:** the I(m) endpoints on the left wing, centre and right wing are correct. The coincidence `z=−κ−d_g = m` on the left slabs is an exact coalescence, which the method allows.
- **Grazing scale:** the ±κ±d_g scale is the right one near a simple root. Near m→κ, |q(z)|≈|q(m)| gives |κ∓z|≈d_g.
- **Clustering:** the z² clustering turns the sqrt endpoint behaviour into smooth behaviour. Jacobians and weights are correct.
- **Missing detail:** the outer-m cuts say "carrier and profile offsets around p0" without an explicit list. The build must enumerate exactly which d values apply to m, per carrier, before any nodes exist. This is mandatory.

## 4. Error propagation
- **Lemma:** ∫_{−M}^{M} 1/|q| = π + 2·acosh(M/κ) is correct and below 4+M (about 13 versus M+4 ≥ 153). Then ∫w < 3M+4 = W_M.
- **Dimensions:** the "1" in w needs the stated nondimensionalization.
- **Budget:** ε = 1e-11/80 is conservative, since 10 J + 30 D = 40 primitives.
- **A-route indicators:** A24 and A48 inner differences are conservative (A48's indicator really measures the inner24 error). Empirical indicators are not rigorous bounds, and the method says so.
- **Tail expressions:** these match `preflight.py:331-336` exactly.
  - outer = 2·Cordinary·CX·CY·E30·F3·WT_3(K,5)
  - middle_J = 4·CX·CY·E30·F3²·(4·121/(5b³))·WT_2(T,1)
  - middle_D = 4·CX·CY·E30·F3²·(36·121/b²)·WT_1(T,1)
- `weighted_tail` needs an integer start, and K=29, T=124 are integers. T−K=95 ≥ 4. The mixed-middle H term is present in the original and correctly flagged as overcounting.
- **Build obligations:** join the ASTs, keep b_* = 3000/11101 distinct from source b, and require positive per-address, per-grade totals at both window pairs.

## 5. Controls
- **Wrong-root mutant:** I confirmed `Dr_wrong = Cd∫{2·X1T·Y1C + (X2T − m·X1T)·Y0C}` against `kernel_components` with `used_qs = q(k+t)`. The input is clipped and the output is unclipped, as stated.
- **Silent controls:** both routes are required and silence counts as unestablished coverage. That is the right posture. For p0=0, a symmetric cancellation could make a control silent. Nothing is asserted in advance.
- **Gap, Route B:** the derivative mutant "(ik)² → (ip0)²" is defined only for Route A's explicit `(ik)^n` factor. Route B embeds the derivative in Q_n and the moment sum. The build must define the Route B mutant (for example Q_n replaced by `(ip0)^n·M0` terms) or restrict the control to Route A and say so.
- **Gap, address choice:** "smallest eligible n1=2 address" should be the smallest NATIVE_MIXED (J) n1=2 entry, since it runs "the same complete J assembly".
- **Finite-window only:** baseline tails do not certify mutated kernels. H, normal and Leibniz controls stay deferred.

## 6. Identity, reuse and persistence
- **Store tests:** the 37 stdlib tests pass per `tooling-tests.log`. They cover immutability, chain tamper detection, route separation and rollback prefix.
- **What the store lacks:** the store at `evidence_store.py` has no full-operand lookup/index, no exact-return reader keyed by request identity, and no disk-reserve, record-size or LRU code. All of that is correctly flagged as a build obligation.
- **Namespace:** it requires `route ∈ {A24, A48, B50, mathematical-inputs}` plus a settings dict. Purpose and check precision therefore have to go in `settings`.
- **Record size:** the `json.iterencode` comment in `put` says a large string can be materialized. A ≤48-node record needs measured live-object RSS inside the 4 GiB cap.
- **Cost:** outer B leaves × joint inner adaptive solves × per-node durable records is unbounded and unknown. The method says so, and I do not forecast it.

## Limits
- Inferred gamma units and the old algebra remain inherited dependencies.
- This is not a numerical result or a runtime proof of any inherited zero count.
- No pairing, leakage or full action claim follows.
- A concrete build still needs its own review.