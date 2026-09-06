# S11c-c2 self-energy fold — independent physics review (fresh Claude leg)

Leg author: fresh general-purpose Claude agent (Codex/astra wrote the script → legs are fresh-Claude + Grok).
Identical prompt: `directives/_legs/S11c_c2_physics_review_prompt.md`. Sandbox: `/tmp/c2_ablate_agent`
(working tree never modified). Evidence: `/tmp/c2_ablate_agent/EVIDENCE_LOG.txt` + helper scripts/outputs there.

**Method.** Derived the physics from `S11c_c2_SHARED_PHYSICS.md` + the real imports, then tested the script by
(1) byte-offset extraction of the emitted `.out` objects, (2) loading the real folded model to inspect
bindings, (3) two mandatory FORM ablations on the TRIAGE core, (4) one-sided corruption checks.

## Findings by section

**A · Fold map — SOUND.** `build_face`@477 returns only `{delta_p_±, d_w_delta_p_±}` slot substitutions — no
`J_s` slot (no double-count). Jets computed by differentiation (`build_face`@506-508). Response elimination is
a genuine 2×2 matrix inverse (`kernel_bridge`@396-398), not scalar division. Load-bearing bulk object is the
two-momentum kernel `z_matrix=[[z0out,z1],[0,z0in]]`. Symbol map verified by model load: `face_velocity[LAB,+]
= W_0·e_W_t·ε/2` dim `[1,-1,0]` = `s11cc1_V` dim `[1,-1,0]` — `V_s` is the interfacial normal velocity;
`μ_θ`→`mu_theta_operator`.

**B · ε-normalization — SOUND.** Increment MULTIGRADE = `{(1,0,0),(1,0,1),(1,1,0),(1,1,1)}`. Per-term scan (LAB
54 / MATERIAL 318 terms): 0 terms without `epsilon_shape` (no O(1)), 0 terms with ε^≥2 (no spurious O(ε²)); the
O(εη) coupling present at (1,1,0). The stray `(0,0,0)` is the resolvent-denominator grade-support
over-approximation the code self-disclaims (`grades`@712-713). Strips exactly one ε.

**C · Close-then-extract — SOUND.** `ORDERING_COMMUTATOR` nonzero (0.40–1.79 MB, all 4). FORM-ablation form1
(freeze input leg `z0in=z0out`, the rule-17 one-leg freeze): increment sha `f4fba6b1→252d2552`. FORM-ablation
form2 (`z1=0`, kill off-diagonal): increment 2,575,195→1,899,256 bytes. Closure threading load-bearing.

**D1 · Density field-vs-field — SOUND.** `density[RHO4_CONSTANT] = rho_br·(1+eta_bg·w1_profile)` (live),
`density[RHOBR_CONSTANT] = rho_br` (frozen); `build_face`@483 rebinds the bare c1-frozen slot to the live
field. `DENSITY_LIVE_MINUS_FROZEN` carries x-dependent `w1_profile` + `eta_bg` — live-field difference, no ρ
jet, ρ_br-based (not ρ_m).

**D2 · Traction covector — SOUND.** `TRACTION_MECHANICAL_CONTRIB` = `Matrix(4,1)`, all 4 components nonzero;
`traction.dot(virtual_velocity)`@792. Native covector, not collapsed.

**D3 · DtN whole-form vs kernel — SOUND (non-vacuous).** `DTN_WHOLEFORM_DEPENDENCE` empty (all 4). Unfiltered
increment∩dtn_operator = `{W_0,c_s0,eta_bg,omega,rho_m,sigma_W}` (background constants only); zero
noncommutative whole-form symbols leak.

**D4 · Traction–slab power pairing — SOUND.** Emits SLAB_POWER / KINETIC_STORED_POWER / TRACTION_POWER /
FACE_GENERALIZED_POWER separately; residual = face_power − traction_power (two independent routes). One-sided
`t_s`-flip `TRACTION_SIGN_RESIDUAL` = 2·traction_power, nonzero (4.5–7.9 MB).

**D5 · Flat-resolvent leg-labeling — SOUND.** `FLAT_SYMBOL_USAGE`: REGRESSION == KERNEL_DIAGONAL = ω·ρ_m/q_out,
RESIDUAL=0. `dtn_flat_symbol` read only @917 as the uniform-limit diagonal regression operand; no MATERIAL
off-diagonal term consumes it.

**D6 · μ_R,bg FORM ablation — SOUND.** `modulus_form`@962 maps `m1_profile → m1_profile²` (plus every jet) — a
genuine FORM change; `MU_R_FORM_RESIDUAL` = 5.3 MB nonzero.

**E · N6 representation-invariance — NOT SETTLED (not a confirmed defect).** `REP_INVARIANCE_RESIDUAL` does not
syntactically vanish (279 genuine-coupling terms). BUT its MULTIGRADE = `{(1,0,1),(1,1,1)}` only — the leading
O(εη) coupling grade (1,1,0) cancels exactly, so leading-order rep-invariance holds; remnant confined to the
σ_W drain sector. Could not prove the σ_W remnant vanishes on this box (needs momentum-integral evaluation —
matches c1's rep-invariance family DEFERRED to ≥64 GB). One-sided `FLIP_FACE_SLOPE` moves each route.

**F · Uniform limit (astra's flagged concern) — SETTLED; not a defect.** `UNIFORM_LIMIT_OPERAND` not a
visibly-zero payload (~16.6 KB), but the genuine transverse-trial-dependent coupling term has Integral integrand
≡ 0 exactly in every block, all 4 cases + both densities. The only nonzero residue is the §3c bare open-slot
bookkeeping `coeff·(δp_minus+δp_plus)·Test`, `coeff = iΛ_A·ε/(ωρ_mτ_A+iρ_m) ≠ 0` — the `−extract(SLAB)` open
piece. So the genuine closure-induced coupling decouples (O(εη)→0). **Caveat:** the emitted uniform object is not
literally zero — a step record must say "the genuine coupling decouples," not "the increment vanishes."

**G · Tautology / adjointness — SOUND.** No `SELF_ENERGY_ADJOINTNESS_RESIDUAL` emitted; `CLOSED_COUPLING_KERNEL`
carries both off-diagonal blocks built by the same weak restriction — the omission is honest (§3b), not a
suppressed check. No emitted `*_RESIDUAL` has the `A−(A/B)·B` structural-zero form.

## Overall verdict
**The fold's physics is SOUND: zero confirmed physics defects.** Two items to carry, neither an error: (E) the
N6 σ_W-sector remnant could not be settled to zero on this box — leading-order invariance holds, the remnant is
consistent with the ≥64 GB deferral but a genuine σ_W pullback incompleteness cannot be ruled out here; (F) the
uniform-limit object is non-vanishing purely from the §3c bare-open-slot representation, so interpretation must
say "the genuine coupling decouples," not "the increment vanishes."
