# S11c-b spec — decision-list (spec) review, two legs, consolidated finding record

**Artifact reviewed:** `directives/S11c_b_SHARED_PHYSICS.md` (orchestrator-authored, 2026-08-27).
**Legs (orchestrator-written → Codex + Grok), rendered prompt:** `directives/_legs/S11c_b_spec_review.md`.
**Raw leg transcripts (outside repo, tree hygiene):** `~/.s11_build/S11c_b_spec_review_grok.txt`,
`~/.s11_build/S11c_b_spec_review_codex.txt` (Codex full transcript persisted at the session tool-results dir).
**Rule 7 TRIGGER:** these two legs are the spec's ONE pass before any engine runs; folded once below, not
iterated to green. Each finding was verified against source by the orchestrator (rule 13) before folding.

## Consolidated defects (8), verdicts, and folds

**D1 — face-DOF contradiction (Codex 1 = Grok secondary; CONVERGENT).** §1a wrote `ζ_± = ±δW/2` while also
defining `ζ_c` ⇒ forces `ζ_c=0`, contradicting S11c-a §3a `ζ_s=ζ_c+sδW/2` and its "no computation may set
`ζ_c=0`". Licenses both engines to drop centre-shift face terms the T-a..T-i substrate carries. VERIFIED real.
FOLD: §1a uses `ζ_s=ζ_c+sδW/2`; distinguishes internal slab fields `{u,δW,θ}` from the two independent face
variables `{ζ_c,δW}`; reserves `ζ_c=0` only for a named centre-fixed uniform regression.

**D2 — N15 channels suppressed by the uniform "u only through gradients" gloss (Grok 1; SERIOUS).** §1a's
"u_T enters the stored energy only through `∇×u`" and §1c's unqualified translation-invariance import state a
UNIFORM fact as if it governs the variable-coefficient energy. Under the broken translation invariance the
N15 spurion invariants (e.g. `(u_T·∇W_bg)θ`, symmetry-allowed: parity-even, w-even, first background jet)
have `u` UNDIFFERENTIATED — exactly the coupling channels the step exists to emit. Both engines would omit
them and agree the coupling is absent (a shared blind spot the residual can't catch). VERIFIED real (the
invariant is symmetry-admissible). FOLD: confine the "enters only through `∇×u`/`∇·u`" statement to the
UNIFORM quotient; state that the non-uniform background admits undifferentiated-`u` couplings contracted with
a background gradient, and §3a promotes `{∂W_bg,∂μ_R,bg}` to spurions so `u` may enter undifferentiated —
which channels appear is the computed §3a result.

**D3 — constitutive coefficient profiles not declared (Codex 2; N12).** The spec asks which uniform
invariants "acquire a variable coefficient" but never NAMES which quantities vary. The ansatz supplies only
`W_bg`, `μ_R,bg`, and the density factorizations; the moduli `B_ρ, C, k_W, κ_W, μ_W` are unaddressed.
Symmetry fixes allowed FORMS, not whether a modulus is constant/function ⇒ engines diverge or both adopt the
forbidden bare `W₀→W_bg` substitution. N12 requires naming which vary. VERIFIED real. FOLD: §2a declares the
varying background fields (`W_bg`, `μ_R,bg`, density reps) and declares `B_ρ, C, k_W, κ_W, μ_W, ρ_m, c_s0`
constant; the geometric `W₀→W_bg(y)` (physical thickness) is legitimate but is NOT the whole energy — §3a
builds the full basis (inherited terms with variable coefficients PLUS new gradient invariants).

**D4 — thickness sector over-counted (Codex 3).** §3a enumerated both `{e_W,∇e_W}` and `{δW,∇δW}`; but
`e_W≡δW/W₀` is one DOF and S11c-a supplies `e_W,bg=(W_0/W_bg)e_W`. Without imposing the map before the rank
test, engines mint duplicate invariants + spurious coefficients; "higher supplied jets" is also ungraded.
VERIFIED real. FOLD: enumerate one thickness coordinate `{e_W,∇e_W}` (as S11b), impose the `e_W,bg` map
before independence testing; drop "higher jets" (the divergence-form operator needs only first coefficient
jets, which the ansatz supplies).

**D5 — the total-divergence quotient does not lift to variable coefficients; comparator pre-registration
could MASK physics (Codex 4; subsumes Grok 3 citation-leak; DEEP).** S11b's uniform 10-dim quotient used
"equivalence modulo total in-plane divergences". For variable `c(x)`: `c∇·F ≡ −(∇c)·F` mod boundary, so the
formerly-equivalent representatives now differ by a first-background-jet invariant whose coefficient is
induced by `∇c` — that is PHYSICS in the coupling kernel, not a representational fold. My §1d cited the
step-record lines that literally print `μ_⊥=μ_R+μ_S/2` (Grok's leak-by-citation) and my §7 pre-registered
the representative identity as an inherited comparator fold — which could hide genuine first-jet coupling.
VERIFIED real. FOLD: §1d drops the leaky citation; states the uniform transverse-stiffness representative
difference is an S11b UNIFORM-LIMIT property, reconciled ONLY in the §5c uniform-limit regression, and warns
the uniform quotient does not lift trivially (IBP of a variable coefficient generates first-jet terms that
are part of the operator). §7 removes the pre-registered representative fold; points at the inherited T7
contract (N8) and states the comparator computes/decides nothing; any S11c-b representational identity is a
RESULT adjudicated after the run (as in S11c-a), never pre-registered.

**D6 — admissibility mis-defined as `ε→0` of the wave operator (Grok 2 = Codex 6; CONVERGENT, SERIOUS).**
§2b defined the operator operand as the §3 perturbation operator at `ε→0`; but §3 is bilinear in the
perturbations so its first-order operator is homogeneous and its `ε→0` limit is identically 0 ⇒ the can-fail
test collapses to a support-bookkeeping identity that never fails while looking like AGREE, and it is close
to the S11c-a-forbidden "perturbation equations as the admissibility test". Codex adds: `{f_hold⁰,t_hold,s⁰}`
was never mapped into the same ordered generalized-force vector as the operator, so the subtraction lacks
signs/measures/row-dimensions. VERIFIED real (aligns with S11c-a §2d's reserved "stationary operator on
𝔅⁰"). FOLD: define `S11CB_ADMISSIBILITY_OPERATOR_OPERAND` as the BACKGROUND-order (ε⁰, zero wave amplitude)
balance / first-variation of the variable-coefficient energy-and-geometry on 𝔅⁰ — the generalized body force
+ face traction the profile sources — expressed in the SAME ordered pairing as `𝒮_hold⁰`; residual is
dimension-matched. §3 gets an explicit background-order balance obligation distinct from the §3b O(ε) wave
operator. Keep the S11c-a §2d prohibition (no `W_bg−W_0` insertion into uniform perturbation eqns).

**D7 — Helmholtz split & adjoint under-specified on a variable background (Codex 5; N5).** `u=u_L+u_T` with
`∇×u_L=0,∇·u_T=0` plus a requested "adjointness relation" is not uniquely defined: a global Helmholtz
projection needs domain, boundary/decay conditions, inverse-Laplacian convention, and zero/harmonic-mode
handling, and N5 removed the plane-wave setting that could give algebraic k-projectors; adjointness needs a
field-space pairing. VERIFIED real. FOLD: define the sectors by the LOCAL differential structure (`u_T` = the
curl-acted part, `u_L` = the div-acted part) as an operator label, NOT a global spectral projection; define
the off-diagonal block operationally as the operator coupling the curl structure to the `{θ,e_W,∇·u}`
structure; supply the variational (stored-energy/kinetic) pairing for any adjointness statement and add a
domain/decay note; do not require a global projector or plane-wave projectors (N5).

**D8 — form ablation couples two independent profiles (Codex 7).** §5b set `∂ᵢW_bg` "and the induced
`∂ᵢμ_R,bg`" to zero TOGETHER; but `w₁` and `m₁` are independent with separate derivative maps — no induced
relation — so the control deletes two physical forms at once and can't separate thickness-slope from
modulus-gradient coupling (the very channels N6 wants distinguished). VERIFIED real. FOLD: §5b runs SEPARATE
one-source ablations of `∂ᵢW_bg` and `∂ᵢμ_R,bg`, holding the other fixed; only the density gradient
definitionally induced by `W_bg` (via §2b) co-varies with the `W_bg` ablation.

## Leg disagreement adjudicated (rule 13)
Codex reported "no transverse-stiffness-identity leak"; Grok reported a leak-by-citation (§1d → the
step-record lines that print `μ_⊥=μ_R+μ_S/2`). Not a contradiction: Codex checked for an INLINE formula (none)
while Grok checked the CITED line's content. Both resolve to the same action under D5 (drop the pointer + the
comparator pre-registration). Adopted.

## Checked-clean by both legs (ruled-out shared-blind-spot classes)
Inherited `U`/`T` fidelity; the uniform symmetry-group list; the background ansatz; the `(ε,η,σ_W)`
truncation-without-stated-grade (correct vs N12/N5); fresh-name/F9 reservations & `v_bulk_normal_0`≠`v_0`;
scope deferrals to S11c-c/d/e and the N11a standing limit; demotion of the uniform limit to a smoke test; the
by-reference inheritance leaves the blind Wolfram engine a complete substrate (aside from D6). No leaked
non-uniform coefficient/sign/grade, basis-count, or admissibility value.
