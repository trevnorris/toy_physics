CLEAR FOR THIS NATIVE TRANSVERSE END-FLUX BUILD

I found no blocking error in the physics or source contract. This is a source and operand review only. I did not run anything and did not decode any pickle, so the RIGHT restoration is not certified. That remains a runtime obligation of the worker.

**Why the matched amplitude change is a linear survival deficit.** The transmitted flux is a†(I+λT1)†(G0+λG1)(I+λT1)a, which equals a†G0a + λ a†K1 a + O(λ²). Reflected amplitude is O(λ), so reflected flux enters at λ² as survival. Incident and transmitted flux are both weighted by the same native J, so the deficit −λ a†K1 a / a†G0a is the right quantity. The ε² factor and the native harmonic factor ¼ sit inside J at both ends and cancel in the ratio.

**Native current and conventions.**
- **Matrix orientation:** `frozen-engine.py:1833` defines the current matrix as `diff(polarized, Minus_i, Plus_j)`. Rows are bra and columns are ket, so G = U†JU is the correct orientation, with no J-versus-Jᵀ ambiguity.
- **Momentum legs:** The bra phase is exp(−i·km·z), the conjugate of a wave of momentum km. The worker's evaluation at (+p,+p) and (−p,−p) is therefore right, and not negating the bra momentum twice is correct.
- **G1 construction:** It has the five pieces the plan requires (explicit end dependence, both momentum legs, both basis factors). They are checked against a direct derivative of B(pR)†J_R(pR,pR)B(pR).
- **Truncation:** The native current is truncated to {1, η, σ, ησ}. That is exact for a first derivative along η=λ, σ=λ/10.
- **Chart:** J goes to the weak chart as S·J·Sᵀ, and U is already weak. The columns of the saved `incident-columns.json` match the weak rows of the uniform-chart object. B(l) is transverse (k·B=0) for every l, so evaluating the constraint and face checks at ±kp symbolically is legitimate.
- **Mass-rate correction:** The worker checks that slab equals conservative plus mass-rate correction at matrix level. It also checks that the correction vanishes on B(pR) for all λ. This is a real check and not an assertion.
- **Face maps:** The saved LEFT face outward velocity and amplitude depend on the θ and e_W amplitudes. Their vanishing on transverse waves is therefore a genuine check.
- **Sign detection:** G0 and Gref must be exactly Hermitian, and their principal minors must be positive. Baseline RIGHT(λ=0) must equal LEFT as full symbolic matrices. Unknown or nonpositive signs stop the run.
- **Gauge and origin:** The arbitrary complex outgoing change X gives G1→G1+X†G0+G0X and T1→T1−X, so K1 is invariant. Profile translation gives T1→T1−iδ_p·b·I, so K1 is again unchanged. Both are algebraic consistency checks and cannot detect an error in T1 itself.
- **Quadratic-route stop:** A nonzero K1 is persisted with an explicit stop of the assumed quadratic-deficit route. A zero K1 is recorded as removing one obstruction only.

**Missing physical arguments.** None prevents this finite prerequisite.

**Limitations the prerequisite retains (not blockers).**
- **Counter-propagating cross term:** The worker does not test the left-side incident×reflected cross term J(−p,+p). It must vanish for a conserved lossless transverse flux, and the faces are inert on transverse waves. That is physical support, not a check.
- **Reflection control:** The reflected-orientation control is vacuous (it only negates a positive number). The real sign protection is the positivity test on Gref.
- **Other channels:** The moment's θ and e_W rows are nonzero, so other-channel radiation is not covered. The same goes for coupled regularity, held-profile work and cross-flux with other channels. K1=0 does not accept any of these.
- **LEFT join:** It relies on the pinned file's saved `passed: true`. The file contains `actual` and `expected`, and the worker does not re-run `exact_structure` on them.
- **RIGHT opaque objects:** Their keys (`native`, `pairing.result.CLOSED_PENCIL_LEGS`, `profileBindings`, `slab`, `conservative`, `knownDimensions`, `acoustic`) are unverified from here. A mismatch fails closed at runtime.
- **Brittle zero tests:** `J.nonzero` and the `is_positive` tests may abort spuriously on radical expressions. K1 zero-detection is structural after `cancel`. These fail in the safe direction, but the authority covers only one execution.