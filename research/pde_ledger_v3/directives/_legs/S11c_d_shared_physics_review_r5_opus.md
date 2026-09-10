# Fresh Claude agent (Opus) review — S11c-d SHARED PHYSICS spec v5 (round 5, DOCUMENT branch)

VERDICT: **SOUND** (1 provenance nit: stale v4 version label). Faithful transcription of the leg's returned review;
the full JSONL agent transcript is the out-of-repo task output file. Fresh-Claude-agent leg (model = Opus).

⚠ ADJUDICATION NOTE (orchestrator): this leg verified the MATH of the §1c Fourier identity and that the spec does not
literally print the `(2π)⁻³` map, but it did NOT test whether the prescribed COMPUTATION is EXECUTABLE against the
consumed objects. The Grok round-5 leg found (and the orchestrator G4-verified) that it is not — `DiracDelta` is absent
from `S11c_c2_exports.py` (c2 strips it; it lives on c1's `dtn_kernel`), the supplied unnormalised-forward sentence +
`f̂_red` + geometry still uniquely fixes the factor (leak), and §1c points WL at the SymPy audit script (blindness
violation). ⇒ the round-5 GATE verdict is NOT-SOUND on F1; refixed in v6 (Codex, per Grok's minimal fix).

---

## Method
Read the sources of truth first (S11c_decisions N5/N6/N7/N10–N15 + the row; c2 step record; c2 SHARED_PHYSICS; N6
disposition; the real c2 engine audit) and formed an own view BEFORE the spec. Verified order bookkeeping, the
strong-edge counterexample, and the ∫f'=Δf Fourier theorem with foreground SymPy; confirmed consumed write-keys exist
in the real 21 MB `S11c_c2_exports.py`.

## Symbolic checks (literal output, leg's own)
- Strong-edge: A_H=−iε sin(ηG), C=sin²(ηG); Born leading η² coeff = G²; C=0 at ηG=π,2π → nonzero Born coeff yet C
  returns to 0 at finite ηG ⇒ no lower bound. Confirmed.
- Flux expansion J_H^(0,1,2) coefficients match term-by-term; with a0=0, C starts at O(λ²)=O(η²), ε² cancels.
  Confirmed the spec's CONDITIONAL O(η²).
- ∫f'(ξ)dξ for f=tanh = 2 = Δf ⇒ (f')_red(0)=Δf. Confirmed.

## Eight scrutiny items — all SOUND (per this leg)
1. Profile class/regime (§1c–§1d): localized thickness interface Δw₁≠0; three grades handled as independent formal
   bookkeepers; numeric tie σ_W=ηW̄₀/L_W only along the homotopy λ≡η; do-not-re-expand framing correct; no freeze.
2. Order bookkeeping (§3c): C not asserted unconditionally — leading order conditioned on the computed
   baseline/interference; nonlinear misread blocked; truncation caveat present. Verified.
3. Strong-edge (§3d): counterexample verified; order-unity edge OUT of scope as a NEW construction; no forced
   factorization. Honest.
4. c2 import (§1b): values per-engine SOUND; cross-engine = N6 covariance thread only; (0)−(0) not operand agreement;
   DEBT 40/76/18 material; F↔S11b distinction handled (F is the c2 increment interpretation; S11b decoupling survives
   as M3 oracle; K_0 computed not typed); both R_N6 (census) + R_cov preserved; B-trap guarded.
5. N6 (§5a): kernel-level N6/N3/N4 kept as c2's unclosed debt; own control = scattering-coordinate covariance + honest
   shape-sensitivity probes; ∇W_bg→0/η→0 rejected; RHO4 structural absence no A−A; anchoring corruption excluded.
6. Two channels (§3b): confinement = transverse survival P_T,surv from reflected+transmitted; continuum vs bound
   distinct; 1D-well theorem correctly declined; profile-functional computed conditional; truncated-model promotion
   guard ‖𝓛_ret⁻¹R_{≥2}‖<1; not a Bloch band.
7. M2/M3 (§answer discipline): objects named; no value/sign/grade beyond the (ε,η,σ_W) contract; O(1) reductio +
   numeric bound withheld; representative "instance never the class"; no freeze.
8. §1c Fourier fold + completeness: [this leg judged it CORRECT AND COMPLETE — the orchestrator/Grok found otherwise,
   see the adjudication note above]. Chain/T7/blind-WL/§8 present; §7 blocks the over-claim that a d-projection closes
   an upstream operand family; all consumed keys exist (s11cc2FieldTheta capital-T does NOT exist — casing hazard real;
   no s11cc2SelfEnergy*/term-origin keys).

## Nit
- Header/lines 16/21 still read "v4"/"Spec v4" though the artifact is the v5 state (Codex v4 + the orchestrator §1c
  Fourier fold + claim-hygiene edits). Provenance/label only; nothing an engine reads changes. Bump to v5 + note the
  fold.

## Verdict
**SOUND** (this leg). [Gate verdict NOT-SOUND on Grok's F1 per the adjudication note.]
