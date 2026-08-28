# S11c-b Wolfram engine — build review, two legs (serialized), consolidated finding record

**Artifact:** `mathematica/S11c_b_brane_operator_mathematica_audit.wl` (Codex-built blind engine,
deliverable-verified: 73KB/1692 lines, 29 tags, exit 0, blindness + flush + no-stray-write confirmed).
**Legs (Codex-written deliverable → fresh Claude agent + Grok), serialized (Mathematica 2-seat), prompt:**
`_legs/S11c_b_wl_build_review.md`. **Raw:** Grok `~/.s11_build/S11c_b_wl_review_grok.txt`; Claude-agent
ablations `…/scratchpad/s11cb_leg/`; Grok ablations `/tmp/s11cb_review_grok/`. Each wrote an independent
derivation + FORM ablations under `timeout 600`, one kernel at a time (peak ~0.3GB; no orphans).

## Consolidated verdicts (both WL legs convergent)

**W1 — Admissibility operand identically zero (BLOCKING; CONVERGENT with SymPy B1 → 4-LEG shared-blind-spot,
SPEC GAP).** `backgroundBalanceFromModel` (L1121–1143) scales the perturbation (ε²) energy fields by
`backgroundOrder`, takes `D[·,backgroundOrder]/.backgroundOrder→0`, then `/.zeroWaveRules` — the ε→0 limit of
the §3b bilinear operator, ≡ 0. Grok ablation-3 and Claude-agent probe3 both show
`BODY_FORCE=<|U→{0,0,0},THETA→0,E_W→0|>` across ALL FOUR branch/density cases (`operand_reduced_zeroQ=True`);
face traction `Σ coeff·Inactive[Div][{0,0,0}]=0`. Claude-agent FORM ablation: injecting a genuine background
term `κ_W·∂W_bg²` leaves the operand UNCHANGED (`operand_UNCHANGED_by_background_injection=True`) — zero data
dependence on the background energy — and computed the correct object `EL wrt W_bg = −2κ_W·∂²W_bg`. Both WL
legs + the SymPy Claude leg independently computed the SAME correct object (profile's own κ_W bending force).
⇒ §3d is UNDER-SPECIFIED (both blind engines read it and produced the identical vacuous operand). **This is a
SPEC fix** (feedback: a one-engine — here two-engine — fix is a spec question; a closed spec is not a correct
spec). **Repair:** clarify §3d — the operand is the first variation of the FULL-field background energy at 𝔅⁰
with the profile's gradients retained, NOT `D[perturbation-energy(scaled fields),backgroundOrder]`.

**W2 — Off-diagonal kernel does not reduce to the decoupled zero in the uniform limit (both WL legs;
WL-specific bug; common root with SymPy B2).** `extractCoupling` (L897–925) applies the sector split
(`transverseInputRules`, `localSectorRules`) only to UNDIFFERENTIATED field occurrences, but the operator is
built from field GRADIENTS (`Derivative[…][uOne]`) — so the projection is INERT on the operator's actual
content. Grok/Claude-agent: the uniform-limit `TRANSVERSE_TO_THETA_EW_UL["THICKNESS_ROW"]` is a large nonzero
diagonal-thickness expression, IDENTICAL after substituting a solenoidal `u_T=curl(A)` (does not depend on the
transverse field). §3c requires the block's uniform limit to be S11b's decoupled zero; it is contaminated by
diagonal dynamics. ⚠ Claude-agent: the §5c smoke test won't catch it (both sides extracted with the same inert
projection cancel). **Repair (engine + spec §3c clarity):** the sector projection must act on the field
GRADIENTS (project `∇×u→u_T`, `∇·u→u_L`) before block extraction, so the block vanishes when the sectors
decouple. (SymPy's B2 is a different bug — parallel non-operator route + representational total-divergence
survivor — but the common root is that neither cleanly isolates the transverse↔thickness block; §3c should
name the extraction more concretely.)

**W3 — Adjointness residual tautological (both WL legs + both SymPy legs = 4-LEG shared; SPEC-CLARITY).**
`ADJOINTNESS_RESIDUAL = D[D[E,transverseProbe],thicknessProbe] − D[D[E,thicknessProbe],transverseProbe]`
(L945–957) ≡ 0 by Clairaut for ANY energy — Claude-agent verified a deliberately non-adjoint energy still
gives 0. §3c asks for the adjointness of the two off-diagonal OPERATOR blocks under the variational pairing;
both engines certified it tautologically via the scalar-energy Hessian symmetry. **Repair (spec §3c clarity):**
the adjointness residual is between the two off-diagonal OPERATOR blocks under the pairing
(`⟨u_T-block, thickness-test⟩ − ⟨thickness-block, u_T-test⟩`), NOT `∂ᵢ∂ⱼE − ∂ⱼ∂ᵢE`.

## Ablated-and-clean (both WL legs)
- **Energy CONSTRUCTED, not substituted (N15):** ablation-2 drops a `W_BG` spurion → WJET labels removed,
  count 26→18, `gammaWidth` kernel term →0; 16 genuine spurion invariants (e.g. `WJET_U_THETA ∝ σ_W·θ·(w1Jet·u)`).
  The undifferentiated-`u` couplings exist (a substitution would omit them). The D2 fix held in WL.
- **Kernel curl/div blocks ARE extracted from the operator:** ablation-4 zeroing the operator collapses the
  blocks to 0 (data-dependent). (Subject to W2's inert-projection contamination.)
- No VERDICT/PASS/FAIL; no assert-before-emit; relationals use `Inactive[Equal]` retaining operands.

## Cross-engine repair scope (all 4 legs consolidated)
**SPEC fixes (round-2 fold; rebuild both engines):**
- §3d admissibility — name the full-field background functional (W1, 4-leg). CONCRETE object known.
- §3c adjointness — operator-block relation under the pairing, not the scalar-Hessian Clairaut identity (W3,
  4-leg).
- §3c block extraction — the sector projection acts on the GRADIENT content; the block is extracted from the
  operator and is the genuine transverse↔thickness coupling (W2 + SymPy B2, both engines wrong differently).
- §5a — the one-sided corruption must PROPAGATE to the kernel and admissibility routes (SymPy B3 emitted A−A).
**ENGINE fixes (per-engine repair directives, against the corrected spec):**
- SymPy: extract kernel from operator (B2), non-tautological independence control (B3), fallback dimension (B4).
- WL: sector projection on gradients (W2). Both: admissibility per corrected §3d; adjointness per corrected §3c.
