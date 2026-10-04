**CLEAR FOR THIS FIRST-ORDER TRANSVERSE SOURCE BUILD**

I found no substantive blockers. I read only `./input` and executed nothing, so every identity below was checked by hand against the supplied operands, not by SymPy.

## Checked against the concrete operands

**Incident lift, chart and wave sign**
- The saved lift columns are `i k×e` for k=(1/5, 1/10, p), so both columns satisfy k·U=0.
- Applying the weak-from-uniform chart gives U₀=(i/5, −ip, 0) and U₁=(0, i/10, −i/5). The chart sends uniform axes (3,1,2) to weak (x,e1,e2).
- `U[0,1]=0` is therefore true, and the worker's `geometric-polarization` zero is consistent with that.
- The wave is `exp(+ik·x − iωt)` for the sign=+1 LEFT leg. The worker uses `(−iω)^t ∏(ik)^n`. The saved contract, k²=6, ω²=cs²K²=9 and p²=119/20 all agree.

**S00 and S10 zeros**
- The saved S00 contains u only through 6i·div u, and S10 only through −6i·w₁·div u. Their other terms carry θ, e_W or `e_W_t`.
- With θ=e_W=0 and k·U=0, both vanish identically in x, not just on shell, for both columns and both faces.

**S01 hand identity**
- I expanded the complete saved S01 by hand using k·U=0 and the transverse profile jets being zero.
- The m₁ coefficient comes out as 400+600K², m₂ as 130p, w₁ as 1900+2100K², and w₂ as 430p, all times c·U_x. That is exactly the plan's bracket, with K²=6.
- Column 1 gives S01=0.

**Native profile L maps**
- Certificate `w1_profile_d1d1` has unscaled value (T³−T)/100 and applied value T³−T, matching the `L^n` rule in the saved contract. Likewise `w1_profile_d1` goes from (1−T²)/20 to (1−T²)/2.
- In the saved sources, the m₂ and w₂ coefficients already carry the extra 1/L. For example the chemical expression has 21/1010 for the m₂ term against 21/101 for the m₁ term.
- Jet fields are applied before profile coefficients.

**Chemical phase and face normalization**
- The chemical σ-coefficient reduces to −i·U_x·bracket/10100, so the worker's candidate identity holds.
- The saved face `chemicalCoefficient` at η=0 is (10+3i)/109, and c·10100·i equals the same number.
- So S01 = Cchem × chemical01/σ, which is the check at worker lines 189–193.
- Chemical at grades 00 and 10 is zero by div u and the θ/e_W factors.
- The velocity is `W₀·e_W,t·ε/2`, so it is zero for e_W=0, and no convection term is invented.

**Flat response and normal signs**
- β=3/(10−3i)=(30+9i)/109, and the saved flat factor 3/(10(q+30/109+9i/109)) equals R₀₀.
- The four route addresses carry normal multipliers 1, +iq, 1, −iq. The inherited cancel-proofs are literally identical operand pairs, so they prove nothing structural. The worker avoids relying on them by building `expected = R·(±iq)` independently and zero-checking it.
- E_W and THETA are the only rows with a nonzero C₀₀ (ε/2 on both faces for E_W; −i(3−10i)ε/109 on the minus face for THETA). Every consumer in the retained C₁₀ position is a profile function multiplying an identically zero S₀₀, so dropping it is legitimate.

**Local cells and rows**
- The baseline row U0 at grade 00 evaluates to 0 for both columns. For column 0 the three nonzero terms are −357i/200, +2023i/3000 and +833i/750, which sum to zero.
- Cell summands keep their own epsilon-power labels, and the worker emits them.
- The grade-10 U1 input shows the original unsummed child-form `left` next to the summed `right`. The worker checks the cell total against both, which is a real child-vs-cell join.

**Controls**
- Omitting the u₁ term moves S01 at T=0 by c·i·1900·w₁(0)·U_x ≠ 0.
- Reversing all spatial derivative signs moves it by −2cU_x·130p·m₂(0), with m₂(0)=−2/3, also ≠ 0.
- Both act through the same `action()` contraction as the real source.

**Gate, argv, helper, containment and hook**
- Gate and argv admission is exact, and the `main` order is gate, then output directory, then containment, then SymPy.
- The 30 s wait in the launcher is a startup handshake only, not a compute deadline.
- Resources match the manifest, and the unchanged guard asserts `RuntimeMaxUSec=infinity` and `Restart=no`.
- The launcher arms the hook (`state=waiting`) before the coordinator is released. The supervisor keeps stdout and stderr.
- On failure the worker writes `failure.json` with the active operation, the posthashes and `checks.json`, with no retry.
- The stdlib tests are admission and selector checks only, and the worker does not treat them as numerical evidence.

## Non-blocking notes

- **Unsupported branch label:** `receiving-sheet` hardcodes `'q>=0 real inside; q=+i sqrt(-q^2) outside'`. Nothing in the packet defines that branch, and no calculation uses the string. I would say "saved `common_outgoing_q` branch, not evaluated" instead.
- **Endpoint join:** stored local endpoints and `force.subs(T,±1)` are emitted separately, with no zero-check between them. A mismatch would be visible only on reading the artifacts.
- **Epsilon bookkeeping:** the cells mix `epsilonPower` 0 and 1 inside the same row-grade, and the worker sums them at ε=1 while the pressure side strips ε. This is consistent, and the powers are emitted per cell. The acceptance reader should keep that bookkeeping in view.
- **Minus-face proof identity:** the minus-face routes reuse proofs 0 and 12, and the worker compares the response-map content but not the minus `id` or `sha256`. The content is identical, so nothing is lost.

## Coverage I did not reach

- No hash was recomputed, so the worker and launcher pins are unverified here (the stated `651cb…` and `1b576…` are unchecked).
- I read about 3 of the 40 profile certificates, and the 400 cells only for U0 grade 00 and a sample of grade 10.
- I did not check the 17.7 MB original-address array, pointer equality, or runtime memory or time under the 4 GiB limit.
- The four inspected flat addresses all have an exactly zero consumer, since they come from U0. The nonzero-C₀₀ E_W and THETA pressure addresses are taken from the split files and not visited as flat addresses.
- `sp.cancel` on √595·I expressions is unexecuted. A false negative would fail closed and preserve evidence, but it would use the single authorized run.
- The inherited unit and algebra proofs were not re-derived.

## Still open, not claimed

Complete first-order face and work maps, matched-end forcing, grazing regularity, the receiving inverse, physical current and power, and leakage all remain open. The pressure forcing reduces to C₀₀R₀₀S₀₁ only for this finite record, and `k_out·f_U` is a coordinate diagnostic, not a loss projection.