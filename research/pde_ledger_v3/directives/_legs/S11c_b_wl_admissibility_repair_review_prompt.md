# Independent decision-list review — S11c-b STEP 3 WL-engine repair directive

You are one of two independent review legs on an **orchestrator-written** directive that proposes to repair
two properties of a blind Mathematica engine. **No builder runs until this review clears.** Your review has
two jobs, and the FIRST is the more important.

Work in the repo `/var/projects/toy_physics`. All paths below are relative to it.

## Artifacts you are given
- The directive under review: `research/pde_ledger_v3/directives/S11c_b_wl_admissibility_repair_directive.md`
- Its measurements file: `research/pde_ledger_v3/directives/_measurements/S11c_b_wl_admissibility_repair_directive.md`
- The shared physics (the source of truth): `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md`
  — read §1c (L93–149), §2a/§2b (L172–235), §3a (L242–270), §3b (L272–288), §3d (L325–356), §1d (L150–171).
- The WL engine the directive would modify (you MAY read it — it is the target):
  `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`

## What you are NOT given, and must not substitute for a derivation
You are deliberately not handed the sibling (SymPy) engine's emitted values for the two contested objects,
nor the residual between engines. **Do not** open `scripts/S11c_b_brane_operator_sympy_audit.py` or any
committed run output to lift an answer. The point of this leg is an *independent* derivation from the shared
physics; a verdict that restates a value instead of deriving it is discarded. Your derivation **script and its
literal stdout are the evidence** — without them your conclusion carries nothing (a prose re-derivation is
worth nothing here, by the same standard this project applies to the engines).

---

## JOB 1 (primary) — settle the spec-question each repair item rests on

Each item repairs a WL output on the claim that WL *under-retains content the shared physics mandates*. That
claim is a spec-question, and it has three possible answers you must distinguish by derivation:
- **GENUINE UNDER-RETENTION** — the shared physics mandates content WL currently omits ⇒ the repair is
  justified in principle.
- **SPEC AMBIGUITY** — the shared physics does not determine the contested content either way ⇒ the spec must
  be fixed first; no engine repair yet.
- **OVER-PRODUCTION** — the shared physics mandates the content be ABSENT, so an engine that carries it is
  wrong ⇒ repairing WL toward it would corrupt the engine. Flag hard.

Derive, in SymPy/Python (save script + literal stdout to named absolute paths under `/tmp`; report the
paths). Prefer Python to avoid the two-seat Mathematica licence; if you must use Mathematica, wrap every
kernel run in `timeout 600` and run one kernel at a time.

### Item A — is a background-order body force conjugate to θ (the density DOF) mandated?
Per §3d, `S11CB_ADMISSIBILITY_OPERATOR_OPERAND` is the background-order (ε⁰) first variation, **with respect
to the full brane configuration**, of the §3a variable-coefficient **energy-and-geometry** functional written
in **full fields**, evaluated at 𝔅⁰. §3d states the full brane configuration contains `W=W_bg+δW`,
`ρ_4D=ρ_4D,bg⁰(1+θ)`, and the brane displacement, and that the gradient content must be that of the full
fields wherever a coefficient varies. §2a (`N12`) lists `ρ_4D,bg⁰`/`ρ_br,bg⁰` among the varying background
fields.

Derive: does taking the background-order first variation of the full-field energy-and-geometry functional
with respect to the density/θ sector yield a **nonzero** generalized body force at 𝔅⁰ (θ⁰=0, all
perturbations zero, profile and its jets retained)? In particular, where the full-field density
`ρ_4D,bg⁰(1+θ)` weights a background-thickness-dependent energy or the areal measure `Σ_E≡ρ_4D W`, does the
θ-linear piece survive with a pure-background (profile-gradient) coefficient at first background-jet order?
Compare with the analogous thickness sector (the E_W body force from `W=W_bg+δW`), which §3d treats
symmetrically. State whether the θ body force is mandated nonzero, mandated zero, or undetermined by the
spec, and show the CAS that decides it.

### Item B — is the thickness inertia coefficient built from the local W_bg(y) or the constant W₀?
Per §1c the supplied kinetic energy is `T = ½ ρ_br⁰ |∂_t u|² + ½ μ_W (∂_t δW)²`, with `δW` the **physical**
thickness perturbation. The normalized thickness strain relates to the physical thickness by the S11c-a §2a
exact map `e_W,bg=(W₀/W_bg)e_W`. §2a promotes every explicit `W₀` that is the physical background thickness
to `W_bg(y)`; §3b forbids freezing a coefficient at its constant binding before differentiation.

Derive: writing `δW` in terms of the normalized strain `e_W` on the non-uniform background, what is the
coefficient of `∂_t² e_W` in the thickness inertia — a function of the **local** `W_bg(y)`, or the
**constant** `W₀`? Show the CAS. State whether the local-thickness inertia is mandated, the constant is
mandated, or the spec is undetermined.

---

## JOB 2 (secondary) — review the directive itself
Report any of:
- **Over-specification (rule 3).** The directive must name the OBJECT to compute, not prescribe the recipe.
  Flag any instruction that dictates *how* to construct a term rather than *what* object must result.
- **Leaked value (rule 5).** The directive must not state, imply, or bound what either repaired object equals
  — no coefficient, sign, order, or residual. Flag any place a target value could be read off (including from
  quoted spec text). Confirm or refute the measurements file's leak-gate claim (§4) yourself.
- **Blindness.** The repair must keep the WL engine importing nothing and re-deriving independently. Flag any
  instruction that would make it copy or reference the sibling engine.
- **Buildability & object-naming.** Could a Mathematica engineer execute each item from the directive alone?
  Are the named spec objects (§3d admissibility operand; §1c/§2a/§3b kinetic term) the correct ones, and are
  the cited WL locations (L1334–1340, L1340; L565/L562/L560; L838/L923) accurate?
- **Scope.** The directive requires everything outside each item stay byte-identical. Is that boundary right,
  or would a faithful repair legitimately have to touch more (a regression risk to name now)?

## Physics filter
Report a finding only if it catches a way the *physics* or the *method* could go wrong — a mandated object
misidentified, a value leaked, blindness broken, a spec-question mis-answered. Do not report "the directive
would be wrong for a different spec."

## Output
A numbered list. Lead with the two JOB-1 verdicts, each as:
`ITEM A/B: <GENUINE UNDER-RETENTION | SPEC AMBIGUITY | OVER-PRODUCTION | CANNOT DETERMINE> — <one line> —
derivation: <abs script path> / <abs stdout path>`.
Then the JOB-2 directive findings (file:line — problem — suggested correction). If the directive is clean,
say so explicitly. Do not edit any file.
