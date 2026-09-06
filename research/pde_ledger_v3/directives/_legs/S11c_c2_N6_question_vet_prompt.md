# Question-vet ONLY — do NOT build, do NOT run, do NOT load any `.out`

You are Codex, consulted to **vet a question** before a diagnostic script is built to answer it. ⛔ Do NOT write
or run a diagnostic. ⛔ Do NOT import `S11c_c2_stdout_loader` or load the ~499 MB `.out`. Read source/spec only
and reason. Your entire job: judge whether the orchestrator's framed question is **the right question at the
retained order**, or whether it is itself a **convenient proxy** (a proxy question is exactly what caused the
defect below). Propose the sharper question/diagnostic if mine is wrong.

## The open finding
`S11CC2_REP_INVARIANCE_RESIDUAL` (N6 rep-invariance) is OPEN. Background:
`research/pde_ledger_v3/_measurements/S11c_c2_physics_review_adjudication.md` (the E section + VERDICT — read it).
A prior orchestrator-authored check (now retired as biased) computed only `residual.subs(σ_W→0) == 0` and read
that as "leading-order rep-invariance holds" — the **over-clear**. `subs(σ_W→0)==0` only proves the residual has
no η-only / constant part (it is ∝ σ_W); it does NOT establish N6 at the retained order.

## The spec claim and the truncation (read these)
- `directives/S11c_c2_SHARED_PHYSICS.md` §5c (≈ line 303): the two anchorings α∈{L,M} are the
  representation-invariance pair; map Eulerian↔material by the field redefinition `Δρ = δρ_E + u·∇ρ⁰`; the residual
  `S11CC2_REP_INVARIANCE_RESIDUAL` **must vanish** ("the same operator in two representations"). No scope
  qualification is stated there.
- `directives/S11c_c1_SHARED_PHYSICS.md` §2c (≈ line 223): `σ_W ≡ η·W̄₀/L_W`, varied independently by η and L_W;
  ⛔ no engine may replace σ_W by η or give them a common order; truncation = **first order in ε and first shape
  order in each background bookkeeper η and σ_W** (N5/N12).

## The code (read these; static)
- `scripts/S11c_c2_selfenergy_fold_sympy_audit.py`: `retained_shape`@≈571 keeps grades with η-power ≤ 1 **and**
  σ_W-power ≤ 1 **independently** (`shape_coefficients`@≈577, the Mul rule `a+c<=1 and b+d<=1`). So the code's
  own truncation RETAINS (η,σ_W) ∈ {(0,0),(1,0),(0,1),(1,1)}.
- `representation_pullback`@≈1136 builds `theta_shift = Σ_i wave_jet(u_i)·∂_i(ρ4)/ρ4` and
  `e_shift = Σ_i wave_jet(u_i)·W_bg_d{i}/W_0`, and shifts ONLY the fields `s11cc2Fieldtheta`, `s11cc2FieldeW`
  (θ, e_W) by those advective terms (composing through derivatives/integrals). The residual is emitted @≈1097.

## The orchestrator's framed question (VET THIS)
> Under the build's own `retained_shape` truncation (each of η, σ_W ≤ 1 independently), does the residual vanish
> at **every retained grade**? Decompose it into its `(η,σ_W)` grades; report which are nonzero — specifically
> whether the genuine first-shape-order grade `(0,1)=O(εσ_W)` (one background gradient) is nonzero, or only
> `(1,1)=O(εησ_W)`. Then: is the implemented `representation_pullback` the **complete** field redefinition
> `Δρ = δρ_E + u·∇ρ⁰`, or undercomplete in a way that leaves exactly that residual — i.e. a
> `representation_pullback` **build defect** vs a genuine order-of-retention scope decision?

## What I need from you (reason from the spec + code; quote lines)
1. **Is this the right question at the retained order, or a proxy?** Does the `(0,1)` vs `(1,1)` grade split
   actually discriminate "build defect" from "in-scope scope decision"? If not, say why and give the discriminator
   that does.
2. **The load-bearing scope subtlety:** since `σ_W ≡ η·W̄₀/L_W` (σ_W is *itself* O(η)), is `(1,1)=O(εησ_W)` a
   **retained** grade (as the code's each-axis-≤1 truncation treats it) or a **second physical shape order** term
   (η·σ_W ~ η², beyond "first shape order")? Does N6 "must vanish" bind at `(1,1)`, or only at `(1,0)`/`(0,1)`?
   This determines whether a nonzero `(1,1)`-only residual is a defect or in-scope. Ground the answer in §2c/§5c.
3. **Is a grade decomposition sufficient, or does the diagnostic also need a direct pullback-completeness test?**
   e.g. should it check whether *completing* the redefinition (shifting the other fields the increment depends on,
   or fully differentiating ρ4 = ρ_br/W_bg through both ρ_br and W_bg) drives the residual to zero at the retained
   grades — which would settle "build defect" affirmatively?
4. **Anything that would make the eventual built diagnostic itself a proxy** (a zero that proves nothing, a
   canonicalizer that hides a nonzero, an integral-dummy artifact, a freeze of a varying field).

## Output
A short verdict: is my question the right one at the retained order (yes / no + the corrected question), the
answer to the `(1,1)` scope question with spec citations, and the decisive diagnostic the build should implement.
⛔ Reasoning + source citations only — no script, no execution.
