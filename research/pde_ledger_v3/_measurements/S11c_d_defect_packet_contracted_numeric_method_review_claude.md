**CLEAR FOR THIS BOUNDED NUMERICAL J/D PACKET METHOD**

I found no fatal error in the method as proposed. The Gaussian identities, the contraction identities and the allocation lemma all check out against the supplied source. The clearance covers the method only. The items below are obligations for the future concrete build, which still needs its own review. Nothing here is a numerical result or a runtime forecast.

## What I read
- **Method and navigation:** `input/method.md`, `input/evidence-guide.md` and `input/contraction-method.md`, all in full.
- **Source:** `original-source/S11c_d_defect_packet_fourier_lib.py` and `..._contraction_lib.py`, both in full.
- **Storage:** `..._evidence_store.py` in full, plus `storage/tooling-tests.log` (37 tests OK). I did not read `preparation.json` or the test source.
- **Index:** `constant-subset-source-index.json`, partly. I read the head of the 64-entry metadata and entries 8346 and 8347 in full. For the other 18 entries I read only grep extracts of the degree, wave multiplier and field text.
- **Not read:**
  - `packet-action-method.md`, `pressure-readiness-method.md`, the inner and geometry libs, `S11c_d_defect_packet_contraction.py`, and the runtime source (guard, supervisor).
  - Everything under `runtime-input/`, including the field, factor, accepted-unit, rule, tail and geometry files.
  - The `contraction/` return files, which were not in the listing I saw (it was cut off at 100 of 3165 files), and the rules directory.
- **Consequence:** I did not independently verify the 17 factor proofs, the 20 complete templates, the accepted units, the rule receipts or the tail allocations. Those remain inherited dependencies, as the method says.

## 1. Constant-field Fourier exception
- **Census:** the index agrees with the stated counts. It has 544 addresses, 64 applicable and 20 eligible. There are 40 `FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED` statuses, 20 in the metadata and 20 in the entries. Each of the 20 entries has two degree-0 polynomials, and the ones I opened have normal multiplier 1. This is a JSON count only. The index says `constantTransformIdentityProvedHere: false`, and I did not recompute the census from `pressure-addresses.json`.
- **Signs, centers and normalization:** I re-derived G_u and G_v from the original Fourier convention, and both match the method. For u, the phase is exp(−iν(x_u+w)) with ν=k−p0 and the 1/(2π) normalization. For v, the test argument is −l with carrier −p0, which gives exp(+i(l−p0)x_v) and no 1/(2π).
  - **Moments:** M1 and the recurrence M(r+1)=−i s²ν Mr + r s² M(r−1) are correct. I checked them by integration by parts, and Q_next=Q'+(ip0−w/s²)Q is the derivative of u.
  - **Derivative order:** the carrier derivative is ik, not ip0. It acts on u before b, which is irrelevant for constant b but harmless.
  - **Old code:** `fourier_lib.constant_reference` (lines 241–253) is the same analytic formula, so the method does not contradict the earlier code.
  - **Independence:** the old comparison of that formula against quadrature lives in banks that were not exported, so I have no independent corroboration of the identity from the packet. The method rightly says Routes A and B share it.
- **Wave-multiplier join, which the method omits.**
  - The index shows `waveMultiplier` values of 1, −p², −1/25, −1/100 and −3i for e_W, d1d1, d2d2, d3d3 and t. These equal P_j·(ik)^n exactly under the delta support k=p.
  - The build must join `waveMultiplier` as the same factor as P_j·(ik)^n. It must not apply both, and it must record that p→k comes from `deltaSupport`.
  - The source coefficient also differs per jet. For example, d1d1 has b = −1/109 − 3i/1090, which is not the e_W value. The build must take b per entry, not per field name.
- **Argument derivative:** the old `request` supports `argumentDerivative` 0 or 1. The J/D formulas use none, so the build must assert it is 0 and refuse otherwise.

## 2. Contractions
- I re-derived all four primitives from the saved J and D kernels.
  - **J and height:** m=k+t gives (k+m) for J and k(m+k) for height, hence C1T+m·C0T and C2T+m·C1T.
  - **Quadratic:** the quadratic term gives X02T.
  - **Reflected:** m=l−t gives k(l+m)/(q(m)+q(l)), hence X1·(Y1CT+m·Y0CT).
- The clipped variable is k for J, Dh and Dq, and l for Dr. The mutant formula Dr_wrong=Cd∫{2·X1T·Y1C+(X2T−m·X1T)·Y0C} also follows from k(2l−m+k).
- Normalized-family sharing is exact linearity, with n∈{0,2} and alpha=b·c·P_j·i^n. The n=2 family needs higher k powers inside the same CrT definitions, which is fine because x_n carries k^n. The build must record the full-expression equality, the per-address budget and alpha for every member.

## 3. Panel geometry
- **d_g(m) is affine on each sector.** On the four sectors cut at −κ, 0 and κ, d_g(m)=min(|m−κ|,|m+κ|) is affine (κ−m on (0,κ), m−κ above κ). The all-intersections plan is therefore well-posed. At κ=√595/10≈2.44, d_g is the right scale for |q(z)|≈|q(m)|.
- **Not a method error, but a risk to the test.** My own analysis, not something I read in a source, says two features of the fixed A rule may fail.
  - Near m=±κ the inner C_rT contain terms like δ·log δ, with δ=q(m). The squared map does not remove that, so A24 against A48 may converge slowly there.
  - A(z) decays like exp(−πL|z|/2) beyond only the 8/L cuts, and it is integrated over long intervals with one GL24 rule.
- The method already says a failure stops the run without added orders, so it is a legitimate test. An optional improvement is to add geometric decay cuts. Another is to state which endpoint each half's square map clusters toward.
- **Plan certificate.** The build must check order constancy, coverage and disjoint intervals per slab, as the method says.

## 4. Error accounting
- **Propagation and lemma.** The product and sum propagation formulas are acceptable as empirical bookkeeping. The w/W_M lemma holds, since ∫|1/q|=π+2acosh(M/κ)<4+M gives ∫w<3M+4. For M=149 and 153 it holds with large margin.
- **Check threshold is too loose.** The 1e-12 absolute check threshold is far larger than the per-address budget ε=1.25e-13 and the inner target ε·w/(16·W_M)≈1.7e-17. The actual discrepancy must be carried into the propagation, as the method already says. An absolute 1e-12 gate is also vacuous in the Gaussian tails, where |ν|≳0.9 gives |G|≲1e-12. I recommend a tighter threshold, a relative one, or both. The exact symbolic equality check stays mandatory.
- **A48 indicator is pessimistic.** The A48 inner indicator is the 24/48 difference, which measures the error of the coarser rule. It is conservative and could cause a stop. This is optional to address.
- **Tail budgets.** The old K27/U75/T122 tail budgets are inherited. Reordering the same windows adds no truncation, but the build must restate that no Fourier x-truncation share is claimed for the constant subset.

## 5. Controls
- **Wrong-root control.** The mutant uses an unclipped output and a clipped input, so its domain is distinct. Responsiveness is not asserted.
- **Adapter mutant.** For constant b it tests only k-power routing, not ordering or the center and sign conventions. That coverage limit should be stated. An optional addition is a center or sign mutant.
- **Coverage deferred, as the method says.** H, normal and Leibniz controls are deferred. Carrier p0=0 may give a zero mutant, which the method already states.

## 6. Request identity and storage
- **Missing request-identity index.** `EvidenceStore` (the write-once SQLite journal in `original-source/S11c_d_defect_packet_evidence_store.py`) has no lookup by request identity, only (namespace, name). The build must add an immutable index that stores the full operands, not just hashes, and retrieves the value without recomputation.
- **Per-record safeguards.** The `iterencode` comment admits a large record may materialize a large string, so the ≤48-node record limit must be enforced in the build. The store has no disk-space check, so the method's 20 GiB reserve and the LRU cap must be implemented and tested in the build.
- **Namespaces.** The store accepts only the namespaces A24, A48, B50 and mathematical-inputs. Controls, the check context and tails need declared homes, or an explicit extension reviewed with the build.
- **Cost.** Cost remains unknown. Nested adaptive B at 50 digits near the root singularities may be expensive. The method forecasts nothing, which is correct.

## Build obligations (mandatory)
1. Join the full constant-subset record per entry, as the method lists, plus the wave-multiplier equality and the p→k substitution above.
2. Assert argumentDerivative=0 and N=1 for all 20, and refuse any other census.
3. Persist both formula operands and the actual differences; use the tighter or relative check; propagate the differences.
4. Certify the exact plan with all intersections before nodes, and state the square-map orientation.
5. Implement the request-identity index and the safeguards above.
6. Keep inherited proofs labeled as inherited, and keep the empirical envelope unlabeled as a rigorous bound.