---
unit_id: 242
batch: VII.2
created_at: 2026-06-02T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-06-03T08:51:18-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 242

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

---

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.wl` (whole file, keep the `_mathematica_audit.wl` filename suffix)

**Issue:** The current `.wl` is a line-by-line port of the SymPy script: identical section structure (I-VII), identical variable names, identical intermediate quantities (`eps`, `Cmix`, `PiTr=(4/3)Cmix`, `varrhoPhys`, `sigmaPhys`), the same artificial "above/below bound" term-split (`sigmaSel - (1/varrho-2) - 1/(3 varrho)`, wl:113 ≡ py:87), the same abstract-zeta-function support-blind substitution device (`RtrSbFn[zeta]` ≡ `sp.Function("Rtr_sb")(zeta)`), the same `Exp[t d]` first-order parametrization, the same orbit-compiler matrix and round-trip test vector. A transliterated second engine cannot independently catch an error in the first; the dual-engine guarantee on this checkpoint is void.

**Required change:**
Re-author the `.wl` so it verifies the SAME claims by GENUINELY INDEPENDENT Mathematica-native routes — a different decomposition and different primitives, not a syntactic rewrite of the Python choreography. Concretely, do NOT mirror the SymPy idioms; instead use routes such as:

- **Placement coordinate (M1, M2):** start from the defining chain `eps = epsW(1 - (2/11) deltaU/(1+deltaU))`, `Cmix = 8 Lambda (1-eps)/Pi^2`, `PiTr = (4/3) Cmix`, `varrhoPhys = Pi^2 PiTr/(16 Lambda)` and confirm `varrhoPhys == (2/3)(1-eps)` and `sigmaPhys == 2 eps/(1-eps)` via `Simplify`/`Together` (acceptable; these are direct identities), but verify `sigma_phys == 4/(3 varrhoPhys) - 2` as an independent consistency leg rather than re-deriving the same sequence.
- **Window inclusion (M3 — see also F2):** use `Reduce` or `Resolve` to certify `1 < PiTr/Cmix < 2` as a genuine inequality over `Lambda>0, 0<eps<1`, NOT the term-split equality. This is the independent route AND fixes F2 simultaneously.
- **Support-blindness (M5):** verify `D[eps, zeta] == 0`, `D[Rtr, zeta] == 0`, `D[Rtarget, zeta] == 0` directly on the closed forms (these are already independent and fine to keep), but REPLACE the abstract-function `RtrSbFn[zeta]` device for the q-packet with a direct argument: substitute the closed forms `Rtr`, `Rtarget`, `epsEta` (none of which contain `zeta`) into `qtr`, `qnt`, `qeta` and confirm `D[qtr, zeta] == D[qnt, zeta] == D[qeta, zeta] == 0` on the actual closed-form packet — no synthetic zeta-functions.
- **Infinitesimal compilers (M6):** instead of the `Exp[t d]` trick, compute the total logarithmic differential directly, e.g. for a quantity `X(chi0, deltaU, ...)` form `dlnX = Sum_v (v/X) D[X, v] dv` with `dv` the per-variable log-drift symbol (i.e. `d ln X = Sum_v D[Log[X], Log[v]] dlnv`, implemented as `v D[Log[X], v]` weighted by `dlnv`). Confirm the resulting `dlnEps`, `dlnRtr`, `dlnRtarget` match the boxed formulas. This is an independent route to the same first-order objects.
- **Orbit packet (M6):** keep the matrix identity but derive `Xi1`, `R1` from their definitions and confirm the compiler matrix maps the direct packet to the orbit packet; the matrix algebra in Mathematica (`LinearSolve`/`Inverse`) is acceptably native.
- **Two-packet split (M7):** confirm `Rtarget Mmix - Cmix == 0` and `D[Mtr, zeta] - Mmix (1-eps)/(1-zeta eps)^2 == 0` by `Simplify` (direct, fine).

Also fix the stale banner per F4 (`STAGE 225` → `STAGE 242`, wl:59) and keep the in-file footer `All Stage 242 ... passed.`

**Claim manifest** (the new `.wl` must independently verify each):
- M1: `varrho_phys = (2/3)(1 - eps)` with `eps = epsW(1 - (2/11) deltaU/(1+deltaU))`, `varrho_phys = pi^2 (4/3) Cmix/(16 Lambda)`, `Cmix = 8 Lambda (1-eps)/pi^2`.
- M2: `sigma_phys = 2 eps/(1-eps) = 4/(3 varrho_phys) - 2`.
- M3: threshold rewrites `eps_WLambda = 1/(2+beta^2)`, `eps_ULambda = beta/(1+beta+beta^2)` from `eps_* = 1 - (3/2) varrho_*`; AND the strict window inclusion `1 < PiTr/Cmix < 2` (i.e. `Cmix < PiTr < 2 Cmix`) certified as an inequality.
- M4: reduced-state bridge `R_target` invariant under `Zhat_W = Z_W Lambda_0/Lambda`.
- M5: `partial_zeta {eps, R_tr, R_target} = 0`, and `partial_zeta {q_tr, q_nt, q_eta} = 0` on the closed-form packet.
- M6: `dln eps = depsW - 2 deltaU/((1+deltaU)(11+9 deltaU)) ddeltaU`; `dln R_tr` per notes §5.2; `dln R_target` per notes §5.2; `Theta_1 = dln R_tr`, `Xi_1 = -dln R_target - (epsEta/(1-epsEta)) depsEta`, `R_1 = dln R_target`; orbit-compiler matrix `{{1,0,0},{0,-1,-c_eta},{0,1,0}}` with `c_eta = epsEta/(1-epsEta)` maps `(dln R_tr, dln R_target, depsEta)` to `(Theta_1, Xi_1, R_1)` and is invertible with `det = c_eta`.
- M7: `R_target M_mix = 8 Lambda (1-eps)/pi^2 = C_mix`; `partial_zeta M_tr = M_mix (1-eps)/(1-zeta eps)^2`.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 242` and confirm: (a) the script exits 0; (b) the banner reads `STAGE 242`; (c) the file is visibly NOT a port of the `.py` (no `RtrSbFn`-style abstract-zeta-function device; window inclusion via `Reduce`/`Resolve`; infinitesimal objects via direct total-log-differential rather than `Exp[t d]` where practical).

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.wl`
  - `mathematica/output/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.txt`
- summary: Re-authored the Mathematica audit around closed-form identities, `Resolve` window inclusion, closed-form q-packet zeta checks, and direct total-log-differential compilers.
- deviation: none

---

## F2 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.py:84-92` (and the mirrored `.wl:110-118` is handled by F1's re-author)

**Issue:** The paper deliverable is the STRICT inclusion `C_mix < Pi_tr < 2 C_mix` (appendix theorem `app-part07-selected-twin-placement`, notes §3). The two checks named "lies above mixed-only bound" / "lies below non-twin bound" are sign-free equalities — `sigma_sel - (1/varrho - 2) - 1/(3 varrho) == 0` (py:85-88) and `(2/varrho - 2) - sigma_sel - 2/(3 varrho) == 0` (py:89-92) — which only assert the gaps equal `1/(3 varrho)` and `2/(3 varrho)`, never that they are positive. The branch could be outside the window and these would still pass.

**Required change:**
Keep the two existing equality checks (they document the gap form) but ADD an explicit strict-window-membership check that exercises the inclusion. The cleanest, non-tautological route: the demand ratio is `Pi_tr/C_mix = 4/3` (a pure number — `C_mix` cancels). Assert membership strictly:

```python
# Strict lowest-symmetric-twin window inclusion: C_mix < Pi_tr < 2 C_mix.
ratio = sp.nsimplify(Pi_tr / C_mix)            # equals 4/3, C_mix cancels
assert ratio == sp.Rational(4, 3), f"demand ratio not 4/3: {ratio}"
assert ratio > 1 and ratio < 2, f"selected branch outside twin window: {ratio}"
print("[ok] selected branch strictly inside lowest symmetric twin window")
```

Place this immediately after the two existing bound checks (after py:92). Do NOT remove the existing two assertions.

(Note: `Pi_tr` is defined at py:60 as `sp.Rational(4,3)*C_mix`, so `Pi_tr/C_mix` simplifies to the literal `4/3`; the `> 1 and < 2` comparison is then a concrete Python bool, not a symbolic-domain question.)

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 242` and confirms a new `[ok] selected branch strictly inside lowest symmetric twin window` line appears and the script exits 0. The window membership must be exercised by a strict inequality that would fail if the ratio were 1 or 2.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.py`
  - `scripts/output/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.txt`
- summary: Added the explicit `Pi_tr/C_mix = 4/3` strict-window assertion while keeping the existing gap-form checks.
- deviation: none

---

## F3 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.py:208,219` (the mirrored `.wl:235,242` is handled by F1's re-author)

**Issue:** `Theta_1 = dln_Rtr` (py:208) followed by `assert_zero("Theta_1 direct-observable identity", Theta_1 - dln_Rtr)` (py:219) is `x - x == 0`, which cannot fail. The paper's substantive claim is `Theta_1 = d ln R_tr`, where `Theta_1 = -C_tr,* Sigma_tr` with `Sigma_tr = (1+deltaU) dln chi0 + (1+chi0) dln deltaU` and `C_tr,* = chi0 deltaU / ((1+chi0)(1+deltaU)(1+chi0+deltaU))` (notes §5.1). The check should confirm the packet-form equals the direct-observable drift.

**Required change:**
Build `Theta_1` from its packet definition and assert it equals `dln_Rtr` (which is itself derived independently at py:184). Replace the trivial check. Suggested edit near py:208-219:

```python
# Theta_1 from the packet form (notes §5.1): Theta_1 = -C_tr,* Sigma_tr.
Sigma_tr = (1 + deltaU) * dchi0 + (1 + chi0) * ddeltaU
C_tr_star = chi0 * deltaU / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU))
Theta_1 = sp.simplify(-C_tr_star * Sigma_tr)
...
assert_zero("Theta_1 packet form matches dln R_tr", Theta_1 - dln_Rtr)
```

This makes the LHS independent of `dln_Rtr` (it is built from `Sigma_tr` and `C_tr,*`), so the residual against `dln_Rtr` is a genuine (and from the algebra, vanishing) test. Keep the existing `Xi_1` and `R_1` checks; they already involve real cancellations. Verify `Theta_1` is still used consistently in the orbit-packet matrix block (py:224, `orbit_packet = sp.Matrix([Theta_1, Xi_1, R_1])`) — with the new definition `Theta_1` equals `dln_Rtr` so the matrix block remains valid.

**Self-test (do this before running):** confirm `-C_tr,* Sigma_tr` expands to `dln_Rtr` (py:184-192 already gives `dln_Rtr = -(chi0 deltaU/((1+chi0)(1+deltaU)(1+chi0+deltaU))) ((1+deltaU) dchi0 + (1+chi0) ddeltaU)`, which is exactly `-C_tr,* Sigma_tr`). So the new assertion is non-vacuous AND must pass.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 242` and confirms a check whose LHS is built from `C_tr_star`/`Sigma_tr` (not assigned `= dln_Rtr` directly) now appears and exits 0. The check must fail if `C_tr_star` is mis-stated (e.g., drop the `(1+chi0+deltaU)` factor).

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.py`
  - `scripts/output/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.txt`
- summary: Replaced the tautological `Theta_1 = dln_Rtr` assignment with the packet-form `-C_tr_star * Sigma_tr` check.
- deviation: none

---

## F4 — stale_output

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.wl:59`

**Issue:** Banner reads `banner["STAGE 225 — ..."]`; canonical numbering (filename, paper card, MANIFEST, in-file footer) is 242. Known stale-label artifact from the EM-extension renumber. Flag the single label only — no batch renumber.

**Required change:**
Change `wl:59` from `banner["STAGE 225 — ACTUAL TWIN-SUPPORT PLACEMENT AND COHERENT ORBIT-LOCK COMPILER"];` to `banner["STAGE 242 — ACTUAL TWIN-SUPPORT PLACEMENT AND COHERENT ORBIT-LOCK COMPILER"];`. (If F1's re-author rewrites the file, ensure the new banner reads `STAGE 242`.)

**Verification command:**
After Codex applies, the verifier confirms `wl:59` reads `STAGE 242` and the regenerated `.txt` header/banner shows `STAGE 242`.

## Applied: F4

- files_changed:
  - `mathematica/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.wl`
  - `mathematica/output/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.txt`
- summary: Corrected the Mathematica banner to `STAGE 242` and refreshed the saved output banner.
- deviation: none
