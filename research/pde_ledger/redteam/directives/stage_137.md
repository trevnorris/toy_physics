---
unit_id: 137
batch: IV.4
created_at: 2026-05-29T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-29T22:42:50Z
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 137 (red-team v2, reconcile + fix)

This is a REWRITTEN directive. The ORIGINAL 2026-05-27 directive (F1–F4) was already
applied in a prior (tainted) pass, so the scripts ALREADY contain anchored assertions
(`Ms_paper`/`Mq_paper`, a `static_limit` check, an `S_q=0` outlet check, the banner is
ALREADY relabeled to `STAGE 137`). A read-only Codex REVIEW
(`redteam/codex_reviews/stage_137.md`, verdict FINDINGS, 3 findings) then flagged that the
applied edits are NOT adversarial:

- **R1** (tautological_check) — the static-limit check `static_limit - (rho_c - sigma_c)`
  is `X − X`: both sides are built from the SAME hand-assigned `rho_c, sigma_c`.
- **R2** (tautological_check) — the outlet check substitutes `S_q = 0`, which DELETES the
  entire `M_q` term; a flipped/missing/zeroed `M_q` would still pass.
- **R3** (insufficient_verification) — `rho_c, sigma_c` remain DIRECTLY ASSIGNED; there is
  no `M_core` / `sp.Matrix` / `Inverse[...]` / `rho_c_schur` / `sigma_c_schur` block. The
  prior directive's F3 matrix-Schur reconstruction was NEVER actually applied.

**R1 and R3 are coupled.** R3 (build `M_core`, invert it, DERIVE `rho_c_schur, sigma_c_schur`
from the physical core stiffness matrix) is what makes R1's static-limit comparison test real
content. Once `rho_c, sigma_c` are DERIVED from `M_core` rather than asserted-against-copies,
the static-limit / full-susceptibility comparison stops being `X − X`. Apply R3 first; R1's
new check consumes the matrix-derived objects.

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under
that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append
`## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line
ranges named. Do NOT introduce any NEW numeric literals — every symbol used
(`K_s, K_q, lam, g_s, g_q, kappa0, gamma0, L, Theta`) is already present in these scripts and
in the Stage 114 owner; no numbers are typed.

After editing, RUN the affected scripts (`timeout 600 python3 <path>` for SymPy,
`timeout 600 math -script <path>` for Mathematica) and iterate until they exit 0 with all
in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator
independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

---

## Reconcile note — PASSED anchors / banner (NO CHANGE)

The Codex review's verdict table marks the following as **PASS** (verdict-table rows 1, 2, 3).
They are present and correct in the CURRENT scripts; do NOT alter them:

- **`M_s` / `M_q` direct closed-form anchors.** SymPy lines 20–29 (`Ms_paper`, `Mq_paper`,
  the two `assert sp.simplify(...) == 0`, the two prints). Mathematica lines 45–49
  (`mSPaper`, `mQPaper`, the two `expectZero`). These compare the build-route `L*rho_c/Theta`
  / `-L*sigma_c/Theta` against the directly-typed paper-card forms — a genuine factor/sign
  catch. KEEP.
- **`sigma_c` r_c-form equivalence.** SymPy lines 48–51 (`expr_rc`, `assert sp.simplify(
  sigma_c - expr_rc) == 0`); Mathematica lines 61–63 (`sigmaCRc`, `expectZero[...]`). This is
  a genuinely independent denominator factorization (`K_s^2*K_q*(1+lam^2/(K_s*K_q))` vs
  `K_s*(K_s*K_q+lam^2)`), NOT X−X. KEEP.
- **Banner already STAGE 137.** Mathematica line 26 already reads
  `banner["STAGE 137 — EXPLICIT CORE-TO-MOUTH GAIN MAP"];` (the prior F4 hold was resolved —
  see `## RESOLVED` at the end). KEEP as-is; do NOT re-edit. There are NO residual `STAGE 120`
  / `Stage 120` strings in either script.

Do NOT re-edit any of the above. They are the reconciled-PASSED surface.

---

## F1 — insufficient_verification (R3): matrix-Schur reconstruction of `rho_c, sigma_c`

(Apply this finding FIRST — F2's de-tautologized static-limit check reuses its objects.)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py:9-10` (the directly-assigned `rho_c`, `sigma_c`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl:34-35` (the directly-assigned `rhoC`, `sigmaC`)

**Issue:**
`rho_c = g_s^2/K_s` and `sigma_c = (K_s g_q - lam g_s)^2 / [K_s (K_s K_q + lam^2)]` are hand-typed.
The notes (`notes/stages/moving_throat_pde_stage097_single_normalization_defect.md:8-89` and
`notes/stages/moving_throat_pde_stage114_concrete_core_schur.md`) DERIVE these as the
Schur complement of the concrete two-channel core stiffness matrix. The script never builds the
matrix, so the residue values behind the entire gain map are unverified — a wrong factor in
`sigma_c` is invisible to every downstream check (this is exactly why R1's static-limit check is
tautological: it compares `sigma_c` against itself).

**The independent source — the PHYSICAL core block (this is the crux, do not back-build it):**
Stage 97 (`...stage097...:16-33`) and the Stage 114 OWNER script
(`scripts/moving_throat_pde_stage114_concrete_core_schur_sympy_audit.py:25-49`) define the
concrete linear core system, NOT in terms of `rho_c, sigma_c`, but in terms of the physical
stiffnesses and couplings:
```
[[K_s,  lam     ]] (s)         (g_s)
[[lam, -K_q*D   ]] (q)  =  u * (g_q),     D = D_W_bare(z) = 1 - kappa0*z^2 - I*gamma0*z^5
```
with mouth feedback `delta_Lambda_core(z) * u = g_s*s + g_q*q`, i.e.
`delta_Lambda_core = v^T M_core^{-1} v` with `v = (g_s, g_q)`. Eliminating `(s,q)` (Schur
complement) gives the boxed Stage 97 result (`...stage097...:44-52`):
`delta_Lambda_core = g_s^2/K_s - (K_s g_q - lam g_s)^2 / [K_s (K_s K_q D + lam^2)]`.
The static residue (`D -> 1`, since `D_W_bare(0) = 1`) is `rho_c - sigma_c`; the `D -> oo`
limit (mixed channel frozen out) is `rho_c` alone. The matrix entries `K_s, K_q, lam` and the
coupling vector `(g_s, g_q)` are the PRIMITIVE physical core data — they are NOT functions of
`rho_c, sigma_c`. Inverting `M_core` is therefore a genuinely independent derivation primitive.

**Required change (SymPy):**

Insert, AFTER line 10 (the `sigma_c = ...` assignment) and BEFORE line 12 (`Ms = ...`), the
following matrix-Schur reconstruction block. Keep lines 9–10 as the hand-assigned reference
values; the new block DERIVES the same quantities from `M_core` and asserts they agree:

```python
# --- F1 (R3): independent matrix-Schur reconstruction of rho_c, sigma_c. ---
# The PHYSICAL core block (notes stage097 :16-33; owner script stage114 :25-49) is the
# two-channel stiffness matrix acting on the internal core coordinates (s, q); rho_c and
# sigma_c are NOT inputs to it. Inverting it is an independent derivation primitive.
kappa0, gamma0 = sp.symbols('kappa0 gamma0', positive=True, real=True)  # bare mixed pair, reused in F2
D_sch = sp.symbols('D_sch', positive=True)
M_core = sp.Matrix([[Ks, lam], [lam, -Kq * D_sch]])   # physical core stiffness matrix
v_coup = sp.Matrix([gs, gq])                          # mouth coupling vector (g_s, g_q)
# Mouth feedback delta_Lambda_core = v^T M_core^{-1} v  (Schur elimination of (s, q)).
delta_Lambda_schur = sp.apart((v_coup.T * M_core.inv() * v_coup)[0], D_sch)
# rho_c is the D -> oo limit (mixed side-channel frozen out); sigma_c is the residual that
# the finite-D term removes at the static point D = 1 (D_W_bare(0) = 1).
rho_c_schur = sp.simplify(sp.limit(delta_Lambda_schur, D_sch, sp.oo))
sigma_c_schur = sp.simplify(rho_c_schur - delta_Lambda_schur.subs(D_sch, 1))
assert sp.simplify(rho_c - rho_c_schur) == 0, "rho_c does not match the M_core Schur residue"
assert sp.simplify(sigma_c - sigma_c_schur) == 0, "sigma_c does not match the M_core Schur residue"
print('rho_c, sigma_c reproduced from explicit two-channel core Schur complement (M_core).')
```
(`(v_coup.T * M_core.inv() * v_coup)[0]` is the scalar quadratic form; `[0]` extracts the sole
entry of the resulting 1x1 Matrix. `sp.apart(..., D_sch)` puts it in the `rho_c - .../(D + r_c)`
partial-fraction form so the `D -> oo` and `D = 1` limits read off cleanly.)

**Anti-tautology guard (SymPy):** `rho_c_schur`, `sigma_c_schur` are computed by INVERTING the
physical stiffness matrix `M_core = [[K_s,lam],[lam,-K_q*D]]` whose entries are `K_s, K_q, lam`
and whose source is `(g_s, g_q)` — none of which are functions of the hand-assigned `rho_c,
sigma_c`. The asserted equalities therefore test that the hand-typed `g_s^2/K_s` and
`(K_s g_q-lam g_s)^2/[K_s(K_s K_q+lam^2)]` actually EQUAL the matrix-inversion result. A
wrong factor or sign in the hand-typed `sigma_c` makes `sp.simplify(sigma_c - sigma_c_schur)
!= 0` and FAILS. This is the same construction the Stage 114 OWNER uses
(`...stage114...:27-42`), confirming the primitive is the project-canonical independent route.

**Required change (Mathematica):**

Insert, AFTER line 35 (the `sigmaC = ...` assignment) and BEFORE line 37 (`mS = ...`), the
mirror block. Mathematica idiom for the same matrix inversion:

```mathematica
(* --- F1 (R3): independent matrix-Schur reconstruction of rhoC, sigmaC via Inverse. --- *)
(* Physical core stiffness matrix (notes stage097, owner script stage114); rhoC, sigmaC *)
(* are NOT inputs to it. Inverting it is the independent derivation primitive. *)
mCore = {{kS, lam}, {lam, -kQ*dSch}};
vCoup = {gS, gQ};
deltaLambdaSchur = Apart[(vCoup . Inverse[mCore] . vCoup), dSch];
rhoCSchur = FullSimplify[Limit[deltaLambdaSchur, dSch -> Infinity], Assumptions -> $Assumptions];
sigmaCSchur = FullSimplify[rhoCSchur - (deltaLambdaSchur /. dSch -> 1), Assumptions -> $Assumptions];
expectZero["rho_c equals M_core Schur residue (D -> Infinity)", rhoC - rhoCSchur];
expectZero["sigma_c equals M_core Schur residue (static D = 1)", sigmaC - sigmaCSchur];
```

Add `dSch, kappa0, gamma0, mCore, vCoup, deltaLambdaSchur, rhoCSchur, sigmaCSchur` to the
`Clear[...]` at line 28, and add `dSch, kappa0, gamma0` to the `Element[{...}, Reals]` /
positivity list in `$Assumptions` (lines 29–32) with `dSch > 0, kappa0 > 0, gamma0 > 0`.
(`zVar`, `piVar`, `sqVar` are ALREADY declared in the current line 30 `$Assumptions` — do not
re-add them.) The mirror uses `Inverse[mCore]` (a structurally distinct primitive from any
closed-form rewrite), satisfying the second-engine independence requirement.

**Anti-tautology guard (Mathematica):** same as SymPy — `rhoCSchur`/`sigmaCSchur` come from
`Inverse[mCore]`; `mCore`'s entries are the physical `kS, kQ, lam` and source `gS, gQ`, never
`rhoC`/`sigmaC`. Comment bodies above contain no `*)` substrings.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 137` / `redteam exec-mathematica 137`
and confirms: (a) both exit 0; (b) SymPy transcript gains the line "rho_c, sigma_c reproduced
from explicit two-channel core Schur complement (M_core)."; (c) Mathematica transcript gains
two new PASS lines: "PASS: rho_c equals M_core Schur residue (D -> Infinity)" and "PASS:
sigma_c equals M_core Schur residue (static D = 1)"; (d) grep confirms `sp.Matrix([[Ks, lam],
[lam, -Kq` (SymPy) and `Inverse[mCore]` with `mCore = {{kS, lam}, {lam, -kQ` (Mathematica) are
present; (e) the SymPy block uses `sp.limit(...)` for `rho_c_schur` while the Mathematica block
uses `Limit[...]` plus `Inverse` — two engines, structurally distinct primitives.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl`
- summary: Added independent matrix-Schur reconstruction of `rho_c`/`sigma_c` from the physical core matrix in both audit scripts.
- deviation: none

---

## F2 — tautological_check (R1): de-tautologize the static-limit / susceptibility check

(Apply AFTER F1 — this consumes `delta_Lambda_schur`, `rho_c_schur`, `sigma_c_schur`.)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py:31-38` (the current Schur static-limit block: `kappa_c, gamma_c, z_var` symbols, `delta_Lambda_core`, `static_limit = sp.limit(...)`, `assert sp.simplify(static_limit - (rho_c - sigma_c)) == 0`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl:51-55` (the current Schur static-limit block: `deltaLambdaCore`, `staticLimit = Normal[Series[...]]`, `expectZero["Schur static limit equals rho_c - sigma_c", staticLimit - (rhoC - sigmaC)]`)

**Issue:**
`static_limit` is the `z -> 0` limit of `delta_Lambda_core = rho_c - sigma_c/(1 - kappa_c z^2 -
I gamma_c z^5)`, which is `rho_c - sigma_c` by construction. The assert then compares it
against `(rho_c - sigma_c)` — the SAME two hand-assigned symbols. So
`static_limit - (rho_c - sigma_c) ≡ 0` regardless of whether `sigma_c` carries the right
factor; if `sigma_c` were misquoted, BOTH sides move together and the check still passes.
Codex R1 confirmed this is `X − X`. The Mathematica mirror (`staticLimit - (rhoC - sigmaC)`)
repeats it.

**Root cause:** the reduced envelope `rho_c - sigma_c/(1 - kappa_c z^2 - ...)` is a packaging of
`rho_c, sigma_c` themselves, so any check that extracts `rho_c, sigma_c` back out of it is
circular. The honest content is that the REDUCED envelope must equal the matrix-Schur source
`v^T M_core^{-1} v` evaluated on the bare denominator `D = D_W_bare(z)`, with the
Stage-97/114 coefficient relations `kappa_c = kappa0/(1+r_c)`, `gamma_c = gamma0/(1+r_c)`,
`sigma_c = sigma_tilde/(1+r_c)`, `r_c = lam^2/(K_s K_q)` (`...stage097...:55-89`;
`...stage114...:34-49`). That comparison is falsifiable: a wrong factor in the hand-typed
`sigma_c` (or a misquoted reduced envelope) leaves a NONZERO residual against the matrix route.

**Required change (SymPy):**

Replace SymPy lines 31–38 (the entire current Schur static-limit block) with the
matrix-anchored full-susceptibility identity:

Before (lines 31–38):
```python
# Schur-complement static-limit anchor (notes Sec. on Stage 97 form).
# delta_Lambda_core(z) = rho_c - sigma_c / (1 - kappa_c z^2 - I gamma_c z^5)
# at z -> 0 must reduce to rho_c - sigma_c.
kappa_c, gamma_c, z_var = sp.symbols('kappa_c gamma_c z_var', positive=True, real=True)
delta_Lambda_core = rho_c - sigma_c / (1 - kappa_c*z_var**2 - sp.I*gamma_c*z_var**5)
static_limit = sp.limit(delta_Lambda_core, z_var, 0)
assert sp.simplify(static_limit - (rho_c - sigma_c)) == 0
print('Schur-complement static limit matches rho_c - sigma_c.')
```

After (lines 31– …):
```python
# F2 (R1): full-susceptibility anchor against the matrix-Schur source (NOT X - X).
# The reduced envelope rho_c - sigma_c/(1 - kappa_c z^2 - I gamma_c z^5) must equal the
# matrix route v^T M_core^{-1} v evaluated on the bare denominator D_W_bare(z), using the
# Stage 97/114 coefficient relations (notes stage097 :55-89; stage114 :34-49). A wrong
# factor in the hand-typed sigma_c leaves a NONZERO residual against the matrix route.
z_var = sp.symbols('z_var', positive=True, real=True)
r_c = lam**2 / (Ks * Kq)
kappa_c = kappa0 / (1 + r_c)
gamma_c = gamma0 / (1 + r_c)
D_W_bare = 1 - kappa0*z_var**2 - sp.I*gamma0*z_var**5
# Independent source: the inverted physical matrix on the bare denominator.
delta_Lambda_matrix = sp.simplify(delta_Lambda_schur.subs(D_sch, D_W_bare))
# Reduced envelope built from the hand-assigned rho_c, sigma_c and Schur coefficient maps.
delta_Lambda_reduced = rho_c - sigma_c / (1 - kappa_c*z_var**2 - sp.I*gamma_c*z_var**5)
assert sp.simplify(delta_Lambda_matrix - delta_Lambda_reduced) == 0, (
    "reduced core susceptibility does not match the M_core Schur source on D_W_bare(z)"
)
print('Reduced core susceptibility matches the matrix-Schur source (full z dependence).')
# Static specialization (z -> 0) now tests DERIVED content, since both sides trace to M_core.
static_limit = sp.limit(delta_Lambda_matrix, z_var, 0)
assert sp.simplify(static_limit - (rho_c_schur - sigma_c_schur)) == 0, (
    "static core residue does not match rho_c_schur - sigma_c_schur from M_core"
)
print('Static core residue matches rho_c_schur - sigma_c_schur from M_core.')
```

**Anti-tautology guard (SymPy):** the LHS `delta_Lambda_matrix` is the inverted physical matrix
`v^T M_core^{-1} v` on the bare denominator; the RHS `delta_Lambda_reduced` is the reduced
envelope built from the hand-assigned `rho_c, sigma_c` and the Stage-97/114 coefficient maps.
The two sides reach the same expression through DIFFERENT routes (matrix inversion vs. reduced
rational form), so the residual is identically zero ONLY when the hand-typed `rho_c, sigma_c`
AND the coefficient maps `kappa_c, gamma_c, r_c` are all correct. The static assert compares
the matrix limit against `rho_c_schur - sigma_c_schur` (also matrix-derived in F1), not against
the hand-assigned pair — so a misquoted hand-typed `sigma_c` is caught by F1's assert and the
matrix side here is unaffected. This is NOT `X − X`.

**Required change (Mathematica):**

Replace Mathematica lines 51–55 (the entire current Schur static-limit block):

Before (lines 51–55):
```mathematica
(* Schur-complement static-limit anchor via Series (independent route from
   SymPy's sp.limit). *)
deltaLambdaCore = rhoC - sigmaC / (1 - kappaC*zVar^2 - I*gammaC*zVar^5);
staticLimit = Normal[Series[deltaLambdaCore, {zVar, 0, 0}]];
expectZero["Schur static limit equals rho_c - sigma_c", staticLimit - (rhoC - sigmaC)];
```

After (lines 51– …):
```mathematica
(* F2 (R1): full-susceptibility anchor against the matrix-Schur source (NOT X - X). *)
(* Reduced envelope must equal Inverse[mCore] source on the bare denominator, using *)
(* the Stage 97/114 coefficient relations (notes stage097, stage114). *)
rC = lam^2/(kS*kQ);
kappaC = kappa0/(1 + rC);
gammaC = gamma0/(1 + rC);
dWbare = 1 - kappa0*zVar^2 - I*gamma0*zVar^5;
deltaLambdaMatrix = FullSimplify[deltaLambdaSchur /. dSch -> dWbare, Assumptions -> $Assumptions];
deltaLambdaReduced = rhoC - sigmaC/(1 - kappaC*zVar^2 - I*gammaC*zVar^5);
expectZero["reduced core susceptibility equals matrix-Schur source", deltaLambdaMatrix - deltaLambdaReduced];
(* Static specialization via Series (independent route from SymPy's Limit). *)
staticLimit = Normal[Series[deltaLambdaMatrix, {zVar, 0, 0}]];
expectZero["static core residue equals rho_c_schur - sigma_c_schur", staticLimit - (rhoCSchur - sigmaCSchur)];
```

Add `rC, kappaC, gammaC` to the `Clear[...]` at line 28 (they are now assigned, not symbolic).
The Mathematica static specialization uses `Normal[Series[...]]` (Taylor extraction) while
SymPy uses `sp.limit(...)` — the prescribed independent-route divergence; do NOT switch them.

**Anti-tautology guard (Mathematica):** identical to SymPy — `deltaLambdaMatrix` comes from
`Inverse[mCore]`, `deltaLambdaReduced` from the reduced envelope; the static assert ties to
`rhoCSchur - sigmaCSchur` (matrix-derived in F1), not to the hand-assigned `rhoC - sigmaC`.
Comment bodies contain no `*)` substrings.

**Claim manifest:**
- M-R1 — the reduced core susceptibility `delta_Lambda_core(z) = rho_c - sigma_c/(1 - kappa_c
  z^2 - i gamma_c z^5)` is the exact Schur complement of the physical core matrix evaluated on
  `D_W_bare(z)`, with `kappa_c = kappa0/(1+r_c)`, `gamma_c = gamma0/(1+r_c)`,
  `r_c = lam^2/(K_s K_q)` (`notes/stages/moving_throat_pde_stage097_single_normalization_defect.md:44-89`;
  cross-checked by owner script `scripts/moving_throat_pde_stage114_concrete_core_schur_sympy_audit.py:30-49`).

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 137` / `redteam exec-mathematica 137`
and confirms: (a) both exit 0; (b) the OLD "Schur-complement static limit matches rho_c -
sigma_c." / "Schur static limit equals rho_c - sigma_c" PASS line is GONE from BOTH transcripts;
(c) SymPy transcript gains "Reduced core susceptibility matches the matrix-Schur source (full z
dependence)." and "Static core residue matches rho_c_schur - sigma_c_schur from M_core."; (d)
Mathematica transcript gains "PASS: reduced core susceptibility equals matrix-Schur source" and
"PASS: static core residue equals rho_c_schur - sigma_c_schur"; (e) grep confirms the SymPy block
uses `sp.limit(` and the Mathematica block uses `Series[` for the static specialization.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl`
- summary: Replaced the circular static-limit check with a full susceptibility comparison against the matrix-Schur source.
- deviation: none

---

## F3 — tautological_check (R2): outlet check with a NONZERO mixed channel

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py:40-46` (the current outlet block: `Pi_var, Sq_var` symbols, `outlet_residual`, `assert sp.simplify(outlet_residual.subs(Sq_var, 0).subs(Pi_var, Ms)) == 0`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl:57-59` (the current outlet block: `outletResidual = (piVar - (mS + mQ*sqVar)) /. sqVar -> 0 /. piVar -> mS`, `expectZero[...]`)

**Issue:**
Substituting `S_q = 0` into the fixed-point law `Pi = M_s + M_q * S_q(Pi)` DELETES the entire
`M_q` term, leaving `Pi - M_s` evaluated at `Pi = M_s`, i.e. `0`. A flipped sign of `M_q`, a
missing denominator in `M_q`, or even `M_q = 0` would still pass — the check never touches the
mixed channel. Codex R2 confirmed this does NOT test gain-pair outlet consistency.

**The independent reference for the mixed channel.** The fixed-point law (notes
`...stage137...:94-101`) is `Pi = M_s + M_q * S_q(Pi)`. To exercise `M_q`, evaluate the map at
a NONZERO susceptibility value `S_q` (keep it a free symbol — no numeric literal) and confirm
the mixed contribution `M_q * S_q` equals the INDEPENDENT matrix-Schur reconstruction
`-L * sigma_c_schur * S_q / Theta` (with `sigma_c_schur` derived from `M_core` in F1). This is
the M_q term reconstructed from the inverted physical matrix, so it pins both the NEGATIVE sign
and the Schur-derived `sigma_c` factor of `M_q`. A flipped/dropped sign or wrong factor in `M_q`
makes the residual NONZERO and FAILS even though `S_q != 0`.

**Required change (SymPy):**

Replace SymPy lines 40–46 (the entire current outlet block):

Before (lines 40–46):
```python
# Card Check 1 (downgraded to carry-forward at stage 134, but a self-contained
# Sq=0 limit of the gain pair fixed-point law is still well-defined here):
# Pi = M_s + M_q * S_q(Pi). At Sq = 0, recover Pi = M_s.
Pi_var, Sq_var = sp.symbols('Pi_var S_q_var', real=True)
outlet_residual = Pi_var - (Ms + Mq * Sq_var)
assert sp.simplify(outlet_residual.subs(Sq_var, 0).subs(Pi_var, Ms)) == 0
print('Outlet consistency reduces to Pi = M_s at S_q = 0.')
```

After (lines 40– …):
```python
# F3 (R2): outlet consistency with a NONZERO mixed channel (paper Checks item 1).
# Family-1 fixed-point law (notes stage137 :94-101): Pi = M_s + M_q * S_q(Pi). At a
# NONZERO susceptibility S_q the map value is Pi = M_s + M_q * S_q; the mixed contribution
# M_q * S_q must equal the matrix-Schur reconstruction -L*sigma_c_schur*S_q/Theta (M_q
# rebuilt from the inverted physical matrix of F1), pinning its sign AND its Schur factor.
Pi_var, Sq_var = sp.symbols('Pi_var S_q_var', real=True)
Pi_map = Ms + Mq * Sq_var                       # fixed-point map at a generic, NONZERO S_q
mixed_contribution = sp.simplify(Pi_map - Ms)   # isolates the M_q * S_q term (not deleted)
Mq_from_schur = -L * sigma_c_schur / Theta      # M_q rebuilt from the matrix-Schur sigma_c
assert sp.simplify(mixed_contribution - Mq_from_schur * Sq_var) == 0, (
    "M_q * S_q outlet term does not match -L*sigma_c_schur*S_q/Theta (sign/factor of M_q)"
)
print('Outlet mixed channel M_q*S_q matches the matrix-Schur reconstruction (S_q != 0).')
# Sanity: a flipped-sign M_q would NOT satisfy the above (guard is non-vacuous).
assert sp.simplify(mixed_contribution - (-Mq_from_schur) * Sq_var) != 0, (
    "outlet check is vacuous: +M_q and -M_q both pass"
)
print('Outlet consistency (paper Checks item 1) verified with nonzero S_q.')
```

**Anti-tautology guard (SymPy):** `mixed_contribution` isolates the `M_q * S_q` term with
`S_q` kept NONZERO (symbolic), so the term is NOT deleted. It is compared against
`Mq_from_schur * S_q`, where `Mq_from_schur = -L*sigma_c_schur/Theta` is rebuilt from the F1
matrix-Schur `sigma_c_schur` (which traces to `Inverse(M_core)`), NOT from `Mq` itself. A
flipped sign or wrong factor in `Mq` makes the residual NONZERO. The second assert
(`mixed_contribution - (-Mq_from_schur)*S_q != 0`) explicitly proves the check distinguishes
`+M_q` from `-M_q`, so it cannot be vacuously satisfied.

**Required change (Mathematica):**

Replace Mathematica lines 57–59 (the entire current outlet block):

Before (lines 57–59):
```mathematica
(* Outlet consistency at S_q = 0 (carry-forward downgrade from card Check 1). *)
outletResidual = (piVar - (mS + mQ*sqVar)) /. sqVar -> 0 /. piVar -> mS;
expectZero["Outlet consistency Pi = M_s at S_q = 0", outletResidual];
```

After (lines 57– …):
```mathematica
(* F3 (R2): outlet consistency with a NONZERO mixed channel (paper Checks item 1). *)
(* Family-1 fixed-point law (notes stage137): Pi = mS + mQ*sqVar; the mixed term *)
(* mQ*sqVar must equal the matrix-Schur reconstruction -lM*sigmaCSchur*sqVar/thetaSigma. *)
piMap = mS + mQ*sqVar;
mixedContribution = FullSimplify[piMap - mS, Assumptions -> $Assumptions];
mQFromSchur = -lM*sigmaCSchur/thetaSigma;
expectZero["outlet mixed channel equals matrix-Schur M_q (sq nonzero)", mixedContribution - mQFromSchur*sqVar];
(* Non-vacuity: +mQ and -mQ must NOT both pass; the flipped-sign residual is nonzero. *)
vacuityResidual = FullSimplify[mixedContribution - (-mQFromSchur)*sqVar, Assumptions -> $Assumptions];
If[TrueQ[vacuityResidual === 0],
  fail["outlet check vacuous: +mQ and -mQ both pass", vacuityResidual],
  pass["outlet consistency non-vacuous (sign of M_q is exercised)"]
];
```

Add `piMap, mixedContribution, mQFromSchur, vacuityResidual` to the `Clear[...]` at line 28.
(`piVar, sqVar` remain declared in `$Assumptions`; `piVar` is now unused but leave it — do not
chase unused-variable cleanup.)

**Anti-tautology guard (Mathematica):** identical to SymPy — `mixedContribution` keeps `sqVar`
symbolic/nonzero so the `mQ` term survives; it is compared against `mQFromSchur*sqVar` built
from the F1 matrix-Schur `sigmaCSchur`. The explicit non-vacuity `If[...]` proves `+mQ` and
`-mQ` give different residuals. Comment bodies contain no `*)` substrings.

**Claim manifest:**
- M-R2 — the Family-1 outlet fixed-point law `Pi = M_s + M_q * S_q(Pi)` with
  `M_q = -L*sigma_c/Theta` (`notes/stages/moving_throat_pde_stage137_core_to_mouth_gain_map.md:73-101`),
  where `sigma_c` is the matrix-Schur residue of F1.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 137` / `redteam exec-mathematica 137`
and confirms: (a) both exit 0; (b) the OLD "Outlet consistency reduces to Pi = M_s at S_q = 0."
/ "Outlet consistency Pi = M_s at S_q = 0" PASS line is GONE from BOTH transcripts; (c) SymPy
transcript gains "Outlet mixed channel M_q*S_q matches the matrix-Schur reconstruction (S_q !=
0)." and "Outlet consistency (paper Checks item 1) verified with nonzero S_q."; (d) Mathematica
transcript gains "PASS: outlet mixed channel equals matrix-Schur M_q (sq nonzero)" and "PASS:
outlet consistency non-vacuous (sign of M_q is exercised)"; (e) grep confirms NO `subs(Sq_var,
0)` (SymPy) and NO `sqVar -> 0` (Mathematica) remain in the outlet block.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl`
- summary: Replaced the `S_q = 0` outlet check with a nonzero mixed-channel check against the matrix-Schur reconstruction and a sign non-vacuity guard.
- deviation: none

---

## RESOLVED — F4 banner relabel (direction (a)), Claude+Codex math consult

The original directive's F4 was a `paper_misalignment` hold: the Mathematica banner at
`mathematica/moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl:26` read
`banner["STAGE 120 — EXPLICIT CORE-TO-MOUTH GAIN MAP"];` while the paper card
(`paper/stages/stage_137.tex:1-2`), the filename, and the notes all assign this unit to
**stage 137**. The `137 - 17 = 120` gap is the known stale-banner offset class from the linear
renumbering (`notes/LINEAR_STAGE_RENUMBERING_PLAN.md`).

**RESOLUTION: direction (a)** — the canonical internal stage number is **137**; the banner is
stale. The fix is to relabel the banner string to the canonical 137. This is a
script-side-only relabel of a display string (no algebra, no paper edit), settled
Claude+Codex per `feedback_claude_codex_resolve_math` (NOT escalated to the user; not a
conceptual change). It is already APPLIED in the current `.wl` — line 26 reads
`banner["STAGE 137 — EXPLICIT CORE-TO-MOUTH GAIN MAP"];`. Codex must:
- CONFIRM line 26 reads `STAGE 137 — EXPLICIT CORE-TO-MOUTH GAIN MAP` (ASCII-safe; the phrase
  contains no `*)` or `_*)` substring, so no comment-terminator hazard).
- `grep -n "120"` both scripts and confirm NO residual `STAGE 120` / `Stage 120` display string
  remains anywhere (banners, comments, or print strings). If any is found, relabel it to 137
  using the same ASCII-safe phrasing.
- Do NOT introduce a paper-side or notes-side change of any kind.
