# S11c-d SymPy engine — STAGED build plan (user decision 2026-09-10: "stage into sub-builds")

The full S11c-d distorted-wave construction is too large for one astra pass (2 monolith passes at `high`, ~2.4 h /
654K tokens, reached only the import + reduction + metadata + pencil-action + off-diagonal-extraction residual; the
scattering core + export were untouched, astra stopping at a runnable checkpoint each pass). Per the user's call
(rule 11 / [[project_v3_requirements_first]]), stage the construction into focused sub-builds on the **one** engine
script `scripts/S11c_d_mixing_scattering_sympy_audit.py` (preserving the single-engine pipeline shape → 2 build legs →
blind WL → T7 → reconcile → step record). Authority throughout = the cleared directive
`directives/S11c_d_sympy_build_directive.md` + spec `399a8516`.

**Heavy-CAS rule for every stage (governing, [[feedback_carrier_first_numeric_pit_for_heavy_cas]]):** emit heavy
objects (the full symbolic S-matrix over 4 cases, the pole solve) as a **carrier-first numeric-PIT fingerprint + SHA
digest**, ⛔ not a full-symbolic dump — so each `.out` stays bounded and reviewable and the comparator can join
fingerprints (the c2-N6 pattern). EXPORT stores the actual **transparent-factored** objects S11c-e binds (⛔ not mere
fingerprints); EMIT (`.out`) carries the object-or-fingerprint for review + the comparator.

## Stages (dependency order; each = a focused astra `high` pass extending the script; stop at a runnable checkpoint,
## ⛔ never fake/empty-substitute)
- **A — reduction substrate + metadata. ✅ DONE** (import wiring; `EdgeReduction`/`hat()`/action-integral reduction;
  ONM at 3 sites; 5-slot census; dimensions/grades; distributional handling; pencil-action + off-diagonal-extraction
  residual). Script 52 KB; `.out` ~97 MB (to be re-bounded by fingerprinting heavy objects in later stages).
- **B — §2 spectrum substrate.** Reduced full pencil `𝓛` from the reduced rows; computed `K₀,K₋,K₊` + reference/end
  baselines; the two full asymptotic block pencils `𝓛₋^full`,`𝓛₊^full`; left/right modes + spectral-projector
  classifiers; both-row reconstruction round-trips.
- **C — §3a scattering.** S11b modal energy current on the reduced operator; nonlinear-pencil normalization; the
  complete two-ended channel S-matrix (both incident ends, every open channel; fingerprint if heavy); conversion
  amplitude; continuum `T→H` flux functional.
- **D — §3b poles/survival + §3c/§3d bookkeeping.** Pole set + Riesz residues/projectors + sheet/normalizability/
  width/closure tests + spectral overlap (fingerprint the solve if heavy); transverse survival functional; amplitude
  components + multigrade; physical flux baseline/interference/quadratic slots; conversion fraction; induced-field
  quadratic form; N12 operands; weak Taylor coefficients + named strong-edge handoff.
- **E — §5 controls.** 5a coordinate-covariance regression + shape-sensitivity mutations (RHO4 absence = computed
  structural absence, ⛔ not `A−A`); 5b three uniform regressions; 5c profile-FORM ablation + edge-vs-bump + modulus
  discriminants (two separate operands); 5d computed flux-normalized dimensionless conversion FORM. Each = object +
  literal residual (both operands).
- **F — export.** `scripts/S11c_d_exports.py` own-rows delta (directive §5 membership: S-matrix / conversion amplitude
  / continuum flux + survival / bound pole set + Riesz + spectral overlap / §3d weak coefficients + recursive
  closure), `assert_delta_is_minimal` bind-closure guard, casewise compact-vs-expanded semantic check. Everything else
  EMIT-only.

## Review
Build all stages into the one engine, then **verify the complete deliverable** (script + bounded `.out` + export +
report) and run the **two build legs** (fresh Claude agent + Grok; Codex-written → that pairing) on the whole engine
(they ablate per-section; the engine is modular and the `.out` bounded). If a stage looks especially load-bearing or
a later stage exposes a defect in an earlier one, checkpoint-review that stage. ⛔ No commit before both legs report;
preserve the reviewed baseline before any repair overwrites it.
