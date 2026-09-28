# Literal Grok two-asymptote implementation review

Source: `/var/projects/toy_physics/_scratch/s11c/s11c-d-two-asymptote-20260927/build-review/grok.json`, field `text`; SHA-256 `77b38215449e781a3de4d5a60c0b163a57134ccc9e55f498425d8c7c38def439`.

I'll start by reading the review instructions and packet index, then stay strictly inside the indexed packet files.The packet is source-only, so I’ll read the review prompt and the two lead artifacts next.Next I’ll read the contract, response plan, and the lift implementation itself.The lift file continues after the native-grade join; I’ll read the rest of the construction.I’ll now read the manifests, source definitions, and evidence files that pin the saved schemas.I’ll inspect the remaining indexed schema, physics, and instrument files next, focusing on pencil conventions and native actions.I’ll read the schema constructors for pencils, grades, and native actions, plus the physical-spec excerpts that define those saved objects.**CLEAR FOR THIS BOUNDED INGREDIENT**

This covers the supplied worker’s construction of first-grade LEFT/RIGHT end-field lifts and native whole-line forcing/commutator expressions at unchanged `LAB_HELD/RHO4_CONSTANT`, inherited physical branch, real poles 16/17, full redundant rank-two residue frames, and independent grades `(1,0)` and `(0,1)`, together with the honest stop `END_LIFT_SAVED_FORCING_DOMAIN_UNRESOLVED`. It does not clear a forcing-domain pairing, an evaluated outgoing action, a solved response, Green/FORM, A11/A12, coincident-point extension, complex-frequency retarded equivalence, or a radiating witness.

This report is source-only. No CAS, pickle restore, or producer run was performed.

## 1. Local quotient, principal part, and end equation

The degree-three finite local quotient is the right object for a simple pole of `R0` with denominator order at most two. With Taylor coefficients of the saved rational chart,

`B = (n[d] - d[d+1] A) / d[d]`

is the next coefficient after the leading principal part (`S11c_d_two_asymptote_lift.py` 454–480; implementation.md 21–32). Forbidden higher pole coefficients and the residue join `n[d-1] - d[d] A` are required to vanish before `B` is used (494–495). The inherited branch is the saved positive radical: `r0 = q0/I` with `r0 > 0` and `r0^2 = minusSquare(k0)` (458–467), using the same quadratic `branchEvidence['minusSquare']` as the chart radicand. Both blocks have nullity two in the validated candidate, with a simple inverse pole (`saved-evidence.json` 227–289; outgoing construct excerpt in `source-definitions.md` 913–915).

The matrix identities `P A = A P = 0` and `P B + P' A = B P + A P' = I` are full 5×5 exact checks (501–504). Those are the correct simple-pole relations, so `B` is the holomorphic inverse coefficient on the full redundant frame.

On that frame,

`C2 = -A P_{e,h} A`  
`C1 = -(B P_{e,h} A + A P_{e,h} B + A P'_{e,h} A)`  
`E = exp(i k0 z) (C1 + i z C2)`

are the principal parts of the ordinary meromorphic product `-R0 P_{e,h} R0` (510–525; implementation.md 40–46). All five residue columns are retained; typed products use residue/kernel units, not flux-normalized channels (192–321).

The independent polynomial action `P f - i P' f' + P_h A`, with constant and linear coefficients required to vanish (518–523), is a real exact cancellation test of `P0 E_h + P_{e,h} E_{00} = 0`. The sign matches the plan’s `exp(+i k z)` convention (`S11c_d_two_asymptote_response_plan.md` 92–101). Pencil `j=0,1` tables already include `1/j!`; the `j=1` join is `mapped.diff(kn)` (428–439; `RectangularModeJets.__init__` in `source-definitions.md` 387–398). LEFT/RIGHT pencils are joined to the independent uniform end symbols after mapping the PLUS/right coordinates onto `kn` and `I sqrt(minusSquare)` (408–425).

That end-equation test does not uniquely pin kernel-valued pieces of `C1` (including `A P'_h A`). Selection of those pieces is the meromorphic product, which is what the code uses.

## 2. Native first-grade action

Grade selection is the intact original integral times the already extracted cell coefficient, after requiring the integral to be grade-free in `(ε, η, σ_W)` (346–358). That is the complete first-grade nonlocal action in the original assembly carrier. The factorization `termJoins` convolution is counted, not used as a whole-line integrand (`S11c_d_continuum_grades.py` 144–169; contract.md 40–45). Bounded cutoffs are not restored. No integral is evaluated.

Bound integral and `Subs` variables are alpha-renamed before insertion (532–555). Field substitution hits the single probe slot, then derivatives/`Subs` (557–574). Local jets are `c ∂_z^n` on the actual field, matching `ReducedActionAssembly.construct` (`source-definitions.md` 129–151). Definite three-part limits are required (350–351). The Abel regulator remains a live coordinate (`alpha` in `coordinate_symbols`, 203–204; `EdgeReduction.prescribe` in `source-definitions.md` 247–271). Half-line pieces that have already become local multipliers stay in `LOCAL`. Units are checked per local/nonlocal address (336–354). Auxiliary `chi` is independent of `η,σ` and of the physical profile binding (588–592).

Differentiating the accepted harmonic source-to-symbol identity is sufficient for **local** finite jets by Leibniz. For **nonlocal** regulated integrals it needs one extra premise: `k`-differentiation of the accepted harmonic identity passes through the existing Abel/half-line representation at the saved regulator, without evaluating the integrals. The worker states that as a method dependency and keeps both raw native plane-jet actions and symbol jets (649–657). It does not treat that identification as new physics evidence.

## 3. Forcing/commutator identity and the stop

Direct forcing is native `F = -L_h E_{00} - L0 T` (623–633). The decomposed form uses native `L0` on the cutoff fields `chi_e E_e` and the constant-coefficient end/reference symbol on the pure end fields. That is the identity that actually reconstructs native `F` from the end equation: `E_e` is an end kernel, not a whole-line homogeneous solution of the variable-profile operator. Reconstruction is linearity plus identical ordered limits, then an exact zero check after leftover `Integral`s are forbidden (642–648). That is a faithful formal regulated identity of the expressions they write, not an unbounded convergence proof.

The nonlocal mutation is saved with `integratedResponseNonzeroEstablished: False` (667–678). The exact end-lift mutation is labelled as a separate control and is not a nonlocal substitute (679–689). Local weighted tails are diagnostics only (690–693). Full-defect end limits stay unevaluated (721). The return status is `END_LIFT_SAVED_FORCING_DOMAIN_UNRESOLVED` (735–746).

That limited stop is useful and justified. No in-scope repair is required for this claim. The available source-specific tests the code already performs, and should keep performing, before stopping are:

- exact vanishing of both polynomial powers in `P f - i P' f' + P_h A`;
- formal reconstruction residual zero after linearity, with no leftover native `Integral`;
- persistence of every nonzero nonlocal operand, ordered limits, and Abel flag as `UNRESOLVED_UNBOUNDED_NONLOCAL_PAIRING`.

A syntactically nonzero `Integral` is not an established integrated response. Local tail limits do not settle the full defect. Those limitations are already recorded.

## 4. Instrument

Containment is checked before SymPy is imported (`main` 755–769; `containment` 54–73): 2 GiB, zero swap, 32 tasks, one CPU, nice ≥ 15, thread env `=1`, 840 s native alarm, `RLIMIT_AS` 2 GiB. The shared guard matches those bounds (`scripts/s11c_guarded_run.py` 30–47, 179–186). Inputs are rehashed into `prehashes.json` and required equal to the pinned manifest (115–117). Each journal op writes the full input pickle and started receipt, then the return pickle and completed receipt (96–108). Failure keeps traceback, `incompleteOperation`, partial index, and `failure-posthashes.json` (776–788). `open('x')` refuses overwrite. There is no retry/resume path. Gate identity covers worker, manifest, guard, 900/840 s, and `automaticRetry is False` (756–763). No producer is imported.

SavedCodec allows only sympy/numpy/builtins/collections (76–81). That is fail-closed if a consumed pickle still carries a producer class; the pinned outgoing/native packets are the same codec family already used for chart/residue restores.

## Substantive findings

None that change what this bounded ingredient computes or may claim.

## Optional observations

These do not require another review cycle for the limited claim.

- `A P'_h A` is kernel-valued under `P A = 0`, so the end-equation check cannot see a missing `P'_h` term. The meromorphic product is the selector (`lift.py` 510–523).
- `reference_plane` is written as `lifts/exp(...)` rather than the already stored polynomial (629–631). Equivalent if SymPy cancels; otherwise the reconstruction guard fails closed.
- Nonlocal-mutation pickle is taken from the first block/grade with a syntactic candidate (658–678). Later grades rely on that one control.
- `require(responsive)` (734) would fail the job if every first-grade `P_h A` were structurally the zero matrix. The contract’s structurally-absent-grade path is then a recorded zero, not a manufactured nonzero control (contract.md 70–71, 106–111). For this variable-profile first grade that is an unlikely structural case.
- Local native plane action versus local Leibniz expansion is not separately zero-checked. Reconstruction plus the end equation already constrain the used sign convention.

## Outside this verdict

Unaccepted remain: full Green/FORM, forcing-domain extension, two-asymptote response, A11/A12, coincident-point distributional extension, complex-frequency retarded equivalence, radiating witness, and any launch/method-acceptance reading of this review. Historical outgoing verdicts (Claude NEEDS REVISION; Grok CLEAR FOR THIS BOUNDED STAGE) apply only to that older artifact. The saved prescription still has `1/k` tails and domain `z ≠ z'`; the fixed slice is evanescent.
