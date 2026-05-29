---
unit_id: 122
batch: IV.3
created_at: 2026-05-29
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-29T12:55:03-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 122

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected script (`python3 <path>`) and iterate until it exits 0 with all in-file checks passing. Getting the script to run cleanly is your job; the orchestrator independently re-runs afterward. (Stage 122 is SymPy-only; there is no Mathematica script.)

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py:43-48` (definition of `T_ratio_minus`/`T_ratio_plus` and the four print lines that follow)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py:64-65` (the two `expect_zero("traction ratio ...")` assertions)

**Issue:**
The traction ratios are defined as bare reciprocals of the compensation roots and then "verified" against those same reciprocals.

Current lines 43-44:
```python
T_ratio_minus = sp.simplify(1/gminus)
T_ratio_plus  = sp.simplify(1/gplus)
```

Current lines 64-65:
```python
expect_zero("traction ratio (-) identity", gminus * T_ratio_minus - 1)
expect_zero("traction ratio (+) identity", gplus  * T_ratio_plus  - 1)
```

Because `T_ratio_minus = 1/gminus` and `T_ratio_plus = 1/gplus` are assigned immediately before the assertions, the checked expressions reduce to `gminus*(1/gminus) - 1` and `gplus*(1/gplus) - 1` by construction. These are pure reciprocal algebra and can never fail. A misquoted natural-branch normalization, a wrong proportionality constant, or a wrong branch normalization in the traction law would still pass. The assertion exercises no physical traction claim.

**Required change (SymPy):**

De-tautologize by deriving the traction ratio from the *independent* proportionality law established in the parent compensation family (stage 119 notes, the boxed normalized-mouth-ratio formula; this is the pre-renumber "Stage 221" the Codex review cites), and from the equal-normalized natural-branch ansatz `g_nat = 1` (stage 122 notes §1) — and only THEN comparing the result to the notes §5 closed form `1/g±`.

The physics: stage 119 §4 gives the boxed law
`g = sqrt(2 Z_q K_s) / (T_m * J_s * c_s * sqrt(mu_0 L_W))`,
i.e. `g = C / T_m` with `C = sqrt(2 Z_q K_s)/(J_s c_s sqrt(mu_0 L_W))` a *background-fixed, branch-independent* constant (it depends only on `Z_q, K_s, J_s, c_s, L_W`, none of which change between the natural and compensated branches). Hence the traction ratio is
`T_m^(±)/T_m^nat = (C/g_±)/(C/g_nat) = g_nat / g_±`,
and stage 122 §5 then quotes the closed form `T_m^(±)/T_m^nat = 1/g_±` *because* `g_nat = 1`. Encoding the ratio as `g_nat/g_±` (with `g_nat` taken from the line-34 ansatz) and asserting it equals `1/g_±` therefore exercises the `g_nat = 1` ansatz and the branch-independence of `C`, instead of inverting `g_±` against itself. This route does NOT use `1/g_±` as its own derivation primitive (no X−X tautology): the independent expression is built from `g_nat` and the cancellation of the symbolic constant `C`, and `1/g_±` appears only as the comparison target.

Replace lines 43-48.

Before (lines 43-48):
```python
T_ratio_minus = sp.simplify(1/gminus)
T_ratio_plus  = sp.simplify(1/gplus)
print("T_m(-)/T_m(nat) =", T_ratio_minus)
print("T_m(+)/T_m(nat) =", T_ratio_plus)
print("numeric T ratio (-) =", sp.N(T_ratio_minus, 20))
print("numeric T ratio (+) =", sp.N(T_ratio_plus, 20))
```

After (lines 43-48):
```python
# Independent traction ratio from stage 119 boxed law  g = C / T_m,
# with C = sqrt(2 Z_q K_s)/(J_s c_s sqrt(mu_0 L_W)) the SAME background-fixed
# constant on every branch. Then T_m^(±)/T_m^nat = (C/g_±)/(C/g_nat) = g_nat/g_±.
# C is carried symbolically so its cancellation is verified, not assumed; g_nat
# is the equal-normalized natural-branch ansatz value (line 34), not 1/g_±.
C = sp.symbols("C", positive=True)
Tm_nat   = C / g_nat
Tm_minus = C / gminus
Tm_plus  = C / gplus
T_ratio_minus = sp.simplify(Tm_minus / Tm_nat)
T_ratio_plus  = sp.simplify(Tm_plus  / Tm_nat)
print("T_m(-)/T_m(nat) =", T_ratio_minus)
print("T_m(+)/T_m(nat) =", T_ratio_plus)
print("numeric T ratio (-) =", sp.N(T_ratio_minus, 20))
print("numeric T ratio (+) =", sp.N(T_ratio_plus, 20))
```

Replace lines 64-65 (keep the assertion LABELS unchanged so output strings match), comparing the independently-derived ratio to the stage 122 §5 closed form `1/g±`.

Before (lines 64-65):
```python
expect_zero("traction ratio (-) identity", gminus * T_ratio_minus - 1)
expect_zero("traction ratio (+) identity", gplus  * T_ratio_plus  - 1)
```

After (lines 64-65):
```python
expect_zero("traction ratio (-) identity", T_ratio_minus - 1/gminus)
expect_zero("traction ratio (+) identity", T_ratio_plus  - 1/gplus)
```

Why this can now fail (non-tautology check): `T_ratio_minus` simplifies to `g_nat/gminus`; the assertion residual is `g_nat/gminus - 1/gminus = (g_nat - 1)/gminus`. This is `0` only because the line-34 ansatz sets `g_nat = 1`. If `g_nat` were changed to any other normalization (a mis-stated natural branch), the residual would be nonzero and the check would FAIL. The symbolic constant `C` must cancel for the ratio to be `g_nat/gminus` at all, so a branch-dependent traction normalization would also break it.

**Claim manifest:**
- **M1** — `expect_zero("traction ratio (-) identity", ...)` exercises the lower-branch traction-renormalization identity `T_m^(-)/T_m^nat = 1/g_-^{F1}`, derived from the proportionality `g ∝ 1/T_m` and the natural ansatz `g_nat = 1`. Notes anchors: `notes/stages/moving_throat_pde_stage122_mouth_source_compensation_test.md:118-129` (§5 `g ∝ 1/T_m`, boxed `T_m^(±)/T_m^nat = 1/g_±`) and `:24-28` (§1 boxed `g_nat = 1`); upstream proportionality `g = sqrt(2 Z_q K_s)/(T_m J_s c_s sqrt(mu_0 L_W))` at `notes/stages/moving_throat_pde_stage119_parent_balance_family.md:116-122`.
- **M2** — `expect_zero("traction ratio (+) identity", ...)` exercises the upper-branch identity `T_m^(+)/T_m^nat = 1/g_+^{F1}` from the same law and ansatz. Notes anchors: same as M1 (`stage122_...:118-129`, `:24-28`; `stage119_...:116-122`).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 122`. Expected output (labels unchanged):
- `T_m(-)/T_m(nat) = 1/g_minus_value` and `T_m(+)/T_m(nat) = 1/g_plus_value` print lines unchanged in numeric value from before (the simplified `g_nat/g±` equals the prior `1/g±` since `g_nat = 1`), so `numeric T ratio (-) = 1.319200163...` and `numeric T ratio (+) = 0.357404273...` remain as in the saved transcript.
- `traction ratio (-) identity = 0`
- `traction ratio (+) identity = 0`
- All other check lines unchanged; script exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py`
- summary: Replaced reciprocal-defined traction ratios with ratios derived from a branch-independent symbolic traction constant and compared them to the stage 122 closed forms.
- deviation: none
