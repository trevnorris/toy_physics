# S11c-d SymPy build — STAGE B (§2 spectrum substrate). Focused resume; extend the existing script.

⭐ **This is one STAGE of a staged build (plan: `_measurements/S11c_d_sympy_build_staging_plan.md`).** Do **only
Stage B** this pass. Your authority is unchanged: the cleared directive `directives/S11c_d_sympy_build_directive.md`
and, through it, the cleared spec `directives/S11c_d_SHARED_PHYSICS.md` (v10, `399a8516`). Read both. Model
`gpt-6-astra`. All directive obligations bind (three script clauses; input-driven — construction operands are the
§1c-reduced rows, ⛔ never a hand-typed payload; PRINT never assert; interpretation is the step record; builder lane =
build → run → report → **stop**, ⛔ no legs/comparator/WL/downstream/commits).

## Read and EXTEND (⛔ do NOT rewrite from scratch, ⛔ do NOT discard the reduction)
`scripts/S11c_d_mixing_scattering_sympy_audit.py` (Stage A: import wiring; `EdgeReduction`/`hat()`/action-integral
reduction; ONM at three sites; five-slot census; dimensions/grades; distributional handling; the pencil action +
off-diagonal-extraction residual). Preserve its correct discipline: `hat()` computes the engine's own reduction and
⛔ does NOT read `FOURIER_PROFILE_BINDINGS` as the answer; that slot is the **(i) 3-D operand** only.

## Stage B — implement ONLY the §2 spectrum substrate (spec §2, directive §1 HELD-PHYSICS)
On the engine's **§1c-reduced** closed operator + reduced coupling kernel (⛔ never unreduced 3-D content — the
reduced-representation rule governs), per `(α,ρ)` case:
1. The reduced full pencil `𝓛(y_n;ω,k_∥)` (the 2×2 sector block operator `[[L_TT,K_TH],[K_HT,L_HH]]`), built from
   the reduced rows; the reduced closed coupling kernel used as the canonical off-diagonal extraction. Emit the
   reduced-operator-off-diagonal-block-vs-reduced-kernel residual as a **surfaced finding** (⛔ do not silently
   override the canonical designation — Stage A's `REDUCED_OFF_DIAGONAL_EXTRACTION_RESIDUAL` may already do this;
   extend/confirm it).
2. `K₀` (the complete (η⁰,σ_W⁰) off-diagonal block at the uniform reference), `K₋ = lim_{ξ→−∞}K`, `K₊ = lim_{ξ→+∞}K`
   in both directions — **computed and emitted**; ⛔ never typed `=0` and ⛔ never inherit a decoupling value from the
   withdrawn F (compute it). Plus the reference/end operator baselines.
3. The two **full** asymptotic block pencils `𝓛₋^full = lim_{ξ→−∞}𝓛`, `𝓛₊^full = lim_{ξ→+∞}𝓛` (end limit =
   simultaneous translation of both normal arguments to that end, then the translation-invariant asymptotic symbol at
   fixed `(ω,k_∥,k_n)`).
4. Left/right modes (incoming, outgoing, evanescent, threshold) of each full pencil, and the transverse-like /
   thickness-like **classifiers** (computed spectral projectors or continuous continuation from the reference sector
   basis; emit the classification map; a degeneracy that makes it ill-defined is reported as a domain limitation, ⛔
   not silently replaced by bare-sector labels).
5. Both-row 3-D↔1-D **reconstruction round-trips** (if not already complete from Stage A's operand pairs).

Every emitted object carries its computed `(ε,η,σ_W)`/`λ` order **and** restored `[L,T,M]` dimension (⛔ neither
omitted). ⛔ Do NOT start §3a/§3b/§3c/§3d/§5/export this pass (later stages).

## Heavy-CAS discipline
Measure the process that runs. If a §2 object is impractical full-symbolic (a large ω-dependent nonlocal pencil /
mode over 4 cases), emit a **carrier-first numeric-PIT fingerprint + SHA digest** of it (⛔ not a full-symbolic dump),
so the `.out` stays reviewable and the comparator can join fingerprints. ⛔ Never run two memory-heavy CAS jobs
concurrently.

## Run + report + stop
Re-run the extended script (emit-before-guard). Produce the updated `.out`, and UPDATE
`_measurements/S11c_d_sympy_builder_report.md` (mark Stage B implemented; per-object which-line-computed-it; the §2
residuals literally; keep the outstanding-constructions list current). ⚠ If you cannot finish Stage B in this pass,
STOP at a runnable checkpoint and report exactly what of §2 remains — ⛔ do NOT fake, empty-substitute, or type any
unimplemented object. Then STOP (⛔ no legs, downstream, or commits).
