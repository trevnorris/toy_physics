# Question-vet ONLY (F and G) — do NOT build, do NOT run, do NOT load any `.out`

You are Codex, consulted to **vet two questions** before diagnostics are built to answer them. ⛔ Do NOT write or
run a diagnostic; ⛔ do NOT load the ~499 MB `.out`. Read source/spec only and reason. Your job: for **each** of the
two questions below, judge whether it is **the right question at the retained order**, or a **convenient proxy**
(a proxy question already cost us the E over-clear), and propose the sharper question/diagnostic if mine is wrong.
⛔ I am giving you no answer.

## Context
STEP A's physics adjudication resolved F and G using **orchestrator-authored** CAS scripts (`verify_F.py`,
`verify_EG.py`) — now retired as biased. I am re-grounding F and G the corrected way (Codex writes the diagnostic,
G1 legs review, I adjudicate). First I need the QUESTIONS vetted. Both concern the emitted **self-energy
increment** `S11CC2_SELF_ENERGY_INCREMENT = extract(close(SLAB)) − extract(SLAB)` (per anchoring α, density ρ),
built at commit `8f3a017f`; its VALUES are unaffected by the recent §5c N6-control correction.

Read for grounding: c2 §3c (increment def, ≈203–222), §3b (adjointness: "emit both blocks; do not dress a
structural zero as a check"), §1d (routing), §5e (uniform limit is a secondary smoke test); the script
`scripts/S11c_c2_selfenergy_fold_sympy_audit.py` (`build_case`, `extract`, the uniform-limit control @≈1081,
`retained_shape`); and the prior adjudication `_measurements/S11c_c2_physics_review_adjudication.md` (F and G
sections) so you can see the two prior positions.

## QUESTION F (uniform-limit decoupling)
> In the uniform limit (`W_bg→W̄₀`, `η→0`, `σ_W→0`), does the **genuine closure-induced coupling** part of the
> increment decouple to zero — i.e. is the surviving `Integral(...)` integrand identically 0 after `.doit()` —
> leaving only the §3c `−extract(open)` **bare-slot bookkeeping** (free `δp_s` symbols)? (Prior split: one leg
> said the coupling decouples; Grok saw a surviving `Z₀·μ_θ Integral` and did not evaluate its integrand.)

Vet it: Is "integrand ≡ 0 after doit, leaving only the −extract(open) bare slots" the right decomposition, or a
proxy (e.g. a special-trial-field zero, an over-evaluated Fourier integral, a canonicalizer that hides a nonzero)?
Is the uniform limit even the right control to make a *decoupling* claim, or only a smoke test (§5e / decisions
N6: it "cannot validate coefficient/sign/parity")? What is the decisive diagnostic + the proxy traps to avoid?

## QUESTION G (directionality / adjointness)
> Is the increment **directional** — the thickness→transverse reverse blocks identically zero, only
> transverse→{θ,e_W} nonzero — and is emitting **both** off-diagonal blocks with **no** independent adjointness
> residual correct per §3b? (⛔ No dissipativity/passivity claim is to be made — directionality only. Prior split:
> one leg said directional/honest-omission; Grok said "blocks not adjoint, suppressed check hides asymmetry.")

Vet it: Is "reverse block ≡ 0 ⇒ directional" the right object, or a proxy (e.g. a basis/representation in which one
block is accidentally zero; a zero from a frozen field; comparing the wrong block pair)? Does §3b actually sanction
emitting both blocks without an adjointness residual, or does it require an independent second route where one
exists? What is the decisive diagnostic, and ⛔ what would make a "directional" conclusion an over-claim (e.g. a
representation-dependent zero, or a dissipativity reading the data does not support)?

## Output
For F and for G separately: is my question the right one at the retained order (yes / no + the corrected
question), the decisive diagnostic a build should implement, and the specific proxy-success traps to forbid.
⛔ Reasoning + source citations only — no script, no execution.
