# Fresh Claude agent (Opus) review — S11c-d SHARED PHYSICS spec v4 (round 4, Codex-authored, DOCUMENT branch)

VERDICT: **SOUND** (1 non-material nit). Faithful transcription of the leg's returned review; the full JSONL agent
transcript is the out-of-repo task output file. This is the fresh-Claude-agent leg (Codex-authored → fresh Claude
agent + Grok); model = Opus.

---

## Method
Derived an own view of S11c-d from `S11c_decisions.md` (N5/N6/N7/N10–N15 + the `:52` row), the c2 step record, the c2
disposition, and the c2 house-format precedent BEFORE opening the spec. Ran two foreground symbolic checks (order
bookkeeping; strong-edge counterexample), verified the Fourier/δ²(Q_∥) convention by hand, and a mechanical
key-existence lookup against the real `scripts/S11c_c2_exports.py`.

## Verified physics checks (leg's own derivations)
- **Order bookkeeping:** converted amplitude O(εη), absolute converted flux O(ε²η²), incident O(ε²), fractional
  C=O(η²) (ε² cancels — a linear theory). Emitting BOTH the O(ε²η²) and O(η²) labels is correct (different objects).
- **Strong-edge counterexample:** F'(0)=1 nonzero, C=sin²(ηG)=0 at ηG=π → a nonzero Born coefficient gives no lower
  bound on strong-edge conversion. The spec's counterexample is exactly right.
- **Fourier/δ²(Q_∥) by hand:** with ξ=n̂·y/L_W, (2π)⁻³∫d³y e^{−iQ·y}f(ξ)=[L_W/(2π)]δ²(Q_∥)f̂_red(s), and the jet
  identity with n̂ᵢ(f′)̂_red(s). Both identities in §1c are correct; (f′)̂_red(0)=∫f′=Δf is the correct edge-vs-bump
  discriminant. [NOTE — this leg verified the integral identities but did NOT check whether the NAMED c2 carrier
  matches the normalized convention; the Grok round-4 leg found it does not (F1) — the c2 carrier is the unnormalised
  forward transform.]
- **Mechanical key-lookup:** all keys the spec names exist; s11cc2FieldTheta (capital T) = 0 hits — lowercase
  s11cc2Fieldtheta is the real key; the casing hazard the spec warns about is real.

## Eight scrutiny items — all SOUND
1. Profile class/regime (§1c–§1d): localized thickness interface (Δw₁≠0, W₋≠W₊), three genuinely independent grades
   (L_W free), forbiddances physically justified, w₁/m₁ independence is the M3-correct non-freezing choice, m₁ a real
   c2 export not invented.
2. Order bookkeeping (§3c): C not asserted O(η²) unconditionally — N12 labels attached to the physical T→H observable
   only conditionally on the computed baseline/interference disposition; asserting a₀=0 would leak the withdrawn F;
   the λ² term verified term-by-term; nonlinear-program guard present.
3. Strong-edge (§3d): honest — only weak coefficient; lab bounds C_strong(1); order-unity edge OUT of scope as a NEW
   construction whose weak limit must reproduce S11c-d; Born-misses flagged; conditional-constraint fallback; no forced
   factorization.
4. c2 import honesty (§1b): faithful — values per-engine SOUND only; cross-engine = N6 covariance thread only; (0)−(0)
   not operand agreement; DEBT 40/76/18 UNADJUDICATED + material; F/G WITHDRAWN; both R_N6=18/288 (raw, schema-
   unmatched) and R_cov no-nonzero preserved; B-trap guarded; no quiet upgrade.
5. N6 (§5a): honest — imports only the closed/extracted kernel, cannot re-derive shape content two ways; keeps
   kernel-level N6/N3/N4 as carried open debt; own control is a scattering-coordinate covariance regression (not a
   two-route independence residual); ∇W_bg→0/η→0 rejected; RHO4 structural absence no A−A; anchoring corruption
   excluded; tilt/advection probes match c2 §5c.
6. Two channels (§3b): confinement = transverse survival; continuum vs bound distinct; 1D-well rejection is a
   correctness strength (unequal asymptotes, multi-component/nonlocal/possibly non-Hermitian); existence
   profile-functional + computed conditional (pole set + Riesz residues, possibly empty) with truncated-operator/
   promotion caveat; not a Bloch band; confinement computed not asserted.
7. Answer/recipe (M2/M3): objects named; no value/sign/parity/grade beyond the (ε,η,σ_W) contract; O(1) reductio +
   numeric bound withheld; no varying quantity frozen; N14 reservations enforced.
8. Completeness/consistency: chain (positional load_model over base + c1 + c2), T7 comparator, blind-WL, §8
   supplied-vs-computed all present; N11a/N11b carried (density maps verified chain-rule-correct); comparator makes
   the honesty point that a projected d-agreement does not close any upstream operand family.

## Nit (does not change what is computed or may be claimed)
- §1a says S11c-d "consumes ONE already-built model, S11c-c2," while §7 correctly needs base + c1 delta + c2 delta
  (c2 is a bind-closure delta). §7 is authoritative and correct; no engine misled. Optional: add "(reconstructed from
  base + c1 + c2 deltas, §7)" to §1a for symmetry.

## Verdict
**SOUND** — nothing outstanding changes what an engine computes or what the spec may claim. [Adjudication note: the
Grok round-4 leg found one must-fix (F1, the c2 Fourier-carrier normalization) that this leg's by-hand Fourier check
did not test at the imported-carrier level; the orchestrator G4-verified F1 against the real c2 engine source and it
stands — so the round-4 gate verdict is NOT-SOUND on F1, folded in v5.]
