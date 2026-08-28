# S11c-b SymPy engine — build review, two legs, consolidated finding record

**Artifact:** `scripts/S11c_b_brane_operator_sympy_audit.py` (Codex-built against `S11c_b_sympy_build_directive.md`,
deliverable-verified: 105KB script, 26MB exports, exit 0, two `s11c_b_` F9c writes).
**Legs (Codex-written deliverable → fresh Claude agent + Grok), prompt:** `_legs/S11c_b_sympy_build_review.md`.
**Raw:** Grok `~/.s11_build/S11c_b_sympy_review_grok.txt`; Claude-agent ablations under
`…/scratchpad/s11cb_review/`; Grok ablations `/tmp/s11cb_review_grok/`. Each leg wrote an independent
derivation + FORM ablations with literal diffs (rule 14). Verdicts below are the orchestrator's after
adjudication (rule 13). NOT yet committed (rule 9: engine commits after repair + re-legs).

## Consolidated verdicts

**B1 — Admissibility operand identically zero (CONVERGENT, BLOCKING).** `admissibility_operator` (L1944–1987)
builds the body force from `zero_wave(euler_derivative(wave_density, …))` — the ε→0 limit of the bilinear §3b
wave operator, which is identically zero. Both legs verified independently: Grok `BODY_IS_ZERO=True`,
`EL_ALL_ZERO=True`; Claude agent the structural identity `admissibility BODY_FORCE == zero_wave(§3b euler)/ε`,
`BODY_FORCE[THETA]=[E_W]=0`, `[U]=(0,0,0)`, `PER_FACE_TRACTION = 0`. The whole operand ≡ 0 ⇒ the can-fail N12
test is vacuous (`RESIDUAL = 0 − support = −support`, cannot detect a non-stationary profile). This is exactly
the ε→0-of-the-wave-operator construction spec §2b/§3d forbid — the D6 spec-review defect RECURRING in the
implementation. **Repair:** construct the background-order (ε⁰) functional with wave-INDEPENDENT profile
content (the profile's own full-field elastic/bending energy — the κ_W term at `W=W_bg`, `~ −κ W_bg''`) and
take its first variation w.r.t. the brane configuration at 𝔅⁰, in the `𝒮_hold⁰` pairing.
⚠ **SPEC-vs-ENGINE (resolve after WL legs):** §3d says "first variation of the total background functional"
but does not pin that the functional is the full-field energy at the background (profile gradients retained).
Under-specified ⇒ the blind WL engine may make the same error. If the WL legs show WL also vacuous ⇒ clarify
§3d in the spec (name the background functional) + rebuild both; if WL is correct ⇒ SymPy-only engine fix
(the spec is adequate and Codex violated the explicit `ε→0` prohibition).

**B2 — Coupling kernel not extracted from the operator (§3c), ships a representational total-divergence
survivor (legs DISAGREED; adjudicated).** `build_kernel` (L1834) uses `paired_kernel_from_density` (direct
mixed second variation of the energy density); `bulk_kernel_from_density` (L1725, the EL/curl-div route) and
`build_operator` are never used for the bulk (only μ_θ face binding). Ablation (Grok A4, Claude agent):
deleting the operator leaves the kernel byte-identical (`KERNEL_SURVIVED_OPERATOR_DELETION=True`). At jet→0
the direct route leaves a `μ_S` `u_T↔∇·u` (`DIV_U`) term (Grok F1: "forbidden uniform survivor"); the Claude
agent COMPUTED (`check_muS_divergence.py`) that this cross-energy is a **total divergence** — representational,
consistent with S11b bulk decoupling, and §1d/§5c route such representative differences to post-run
adjudication. **Adjudication:** the physical coupling is correct (total-divergence survivor is representational,
Claude leg right); but the engine violates spec §3c ("extract the block of the §3b operator") — the EL/operator
route (`bulk_kernel_from_density`) returns 0 there (Grok A1b), so extracting from the operator both satisfies
§3c and removes the representational pollution (Grok's repair right). **Repair (engine):** set the bulk kernel
blocks from the built operator's curl↔div EL structure, not the parallel `mixed_variation` route.

**B3 — Tautological controls (Grok F4 + Claude agent Secondary A).** (a) `task_independence` (L2140–2151):
for `COUPLING_KERNEL` and `ADMISSIBILITY_OPERATOR`, baseline and "corrupted" are IDENTICAL calls (no
corruption applied — only `SLAB_OPERATOR` gets `corrupt_material_constraint=True`), so the §5a/N6 residual is
`A−A ≡ 0`, a control that cannot fail. (b) `VARIATIONAL_ADJOINTNESS PAIRING_RESIDUAL` (L1762–1799): forward −
reciprocal ≡ 0 for ANY density (mixed partials of a scalar commute), verified on an arbitrary hand density —
non-failable ("operand theatre", rule-2 corollary 3). **Repair:** wire the named one-sided source corruption
to the kernel and admissibility routes (largely propagates once B2 is fixed and the kernel is extracted from
the corrupted operator); make adjointness a genuine second-route check or emit the objects honestly with "no
independent second route" rather than a structural zero.

**B4 — Dormant dimension crash (Claude agent Secondary C).** `uniform_coefficient` fallback (L1150–1159) mints
`a_s11cb_uniform_XX` without registering a dimension; if the quotient ever selects an uncatalogued uniform
invariant, `dimension_of` (L816) raises `KeyError` and kills the primary emit (observed under the μ_S FORM
removal). Dormant in the baseline (fallback never fires) but a latent crash on any future basis selection.
**Repair (engine):** register a dimension (or graceful report) for the fallback uniform invariants.

## Ablated-and-clean (both legs; the physics held)
- **Coupling is gradient-driven:** `transverse→THETA`/`→E_W` collapse to exactly 0 when the background first
  jet is removed (both legs, ablation 1). Genuinely built from `∂W_bg`/`∂μ_R,bg`.
- **Energy CONSTRUCTED, not substituted (N15 works — the D2 fix held):** dropping an undifferentiated-`u`
  spurion channel changes the new-invariant count (Grok 16→14; Claude agent 16→10) and removes the
  undifferentiated-`u`×θ kernel coupling — a `W₀→W_bg` substitution would omit these (both legs, ablation 2).
- **Kernel is a genuine mixed variation (not hand-typed):** it moves under every energy FORM change (ablation
  2 + 4). (B2 is about WHICH route, not hand-typing.)
- No assert-before-emit (every `raise` is an operational/type guard, physics residuals emitted unguarded);
  no VERDICT/PASS/FAIL; no native-boolean residual.
- ⚠ Triage tool `reduction/derived_or_declared.py` no longer exists; both legs traced §4 objects to their
  computing lines by hand instead (all COMPUTED except the admissibility operand, which computes to 0 = B1).

## Repair scope (to author after the WL legs resolve the §3d spec-vs-engine question)
B1 (admissibility — likely spec §3d clarification + engine), B2 (extract kernel from operator), B3
(non-tautological controls), B4 (fallback dimension). B1 is blocking; B2 is spec-compliance + removes
representational pollution; B3/B4 are correctness of the control/robustness.
