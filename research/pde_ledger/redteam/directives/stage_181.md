---
unit_id: 181
batch: V.2
created_at: 2026-05-30T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-30T02:02:59-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

## Consult resolution (Claude+Codex read-only, 2026-05-30, codex_session 019e77af)

**Q2 — Mathematica feasibility, F1/F2 (CONCUR, path (b)):** the perturbation block must NOT
rely on the post-line-55 `$Assumptions` (it drops `zeta`/`zW`/`omegaW2`/positivity that
`t2Loaded`/`lamNorm`/`rTarget` need). **Primary path = recompute the perturbed quantities
(`t2DirectPert`, `rTargetPert`) directly from their CLOSED FORMS** with LOCAL assumptions for
`s` and the drift symbols (the file already proves `t2Loaded == t2Direct`, wl:48). These are
first-order log-derivatives of rational products → they reduce to an exact symbolic 0; if a
`FullSimplify` does not close, reformulate locally — do NOT raise the 600s cap.

# Codex directive — unit 181

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond what each finding names. Edit exactly the file regions named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

Shared perturbation map (used by F1 and F2). Introduce one real perturbation symbol `s` standing for \(\epsilon\lambda_A\), and perturb the five primitive coherent-branch variables exactly as the notes (§3) prescribe:
- `ZW    -> ZW * (1 + s*zetaZ)`        (since δ ln Z_W = ελ ζ_Z)
- `OmegaW2 -> OmegaW2 * (1 + s*omegaW)` (since δ ln Ω_W² = ελ ω_W)
- `chi0  -> chi0 + s*chi1`
- `epsW  -> epsW + s*epsW1`
- `deltaU -> deltaU + s*deltaU1`
- `eps_eta -> eps_eta + s*eta1`         (only needed for F2)
The split-blocking ratio under perturbation is `eps_pert = (epsW + s*epsW1)*(1 - Rational(2,11)*(deltaU + s*deltaU1)/(1 + deltaU + s*deltaU1))`. Take the first-order drift as `diff(log(EXPR_pert), s).subs(s, 0)` (SymPy) / `D[Log[exprPert], s] /. s -> 0` (Mathematica). Do NOT redefine the existing `Xi1`/`R1` closed forms (lines 82–95 sympy / 64–66 math) — keep them as the paper targets and compare the newly-derived drifts against them.

## F1 — insufficient_verification

**Target:**
- SymPy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.py` — insert after line 100 (`print("R_1  =", R1)`), before `banner("Tracking-factor drift")` at line 102.
- Mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl` — insert after line 72 (`Print["R_1  = ", fmt[r1]];`), before `banner["Tracking-factor drift"]` at line 74.

**Issue:** The coherent-branch defect law \(\Xi_1\) (paper eq:app-part05-Xi1-coherent-physical) and the headline support-blindness of \(\Xi_1\) (eq:app-part05-support-blindness, third equality) are not verified. `Xi1` (sympy line 82) is typed in as the paper answer; the script never derives it from \(\delta\ln\mathcal T^2/(\epsilon\lambda)\), and never checks \(\partial_\zeta\Xi_1 = 0\) (true only because line 82 is ζ-free).

**Required change (SymPy):** Add a derivation of `Xi1` from the transfer shape and a ζ-blindness check via the support-loaded route. Use the shared perturbation map. The defect symbols `zetaZ, omegaW, chi1, epsW1, deltaU1` already exist (lines 72–73); add `s = sp.symbols("s", real=True)`.
```python
s = sp.symbols("s", real=True)
eps_pert = (epsW + s*epsW1) * (1 - sp.Rational(2, 11) * (deltaU + s*deltaU1) / (1 + deltaU + s*deltaU1))
T2_direct_pert = (ZW*(1 + s*zetaZ)) * (1 + (chi0 + s*chi1))**2 \
                 / ((OmegaW2*(1 + s*omegaW)) * (1 - eps_pert)**2)
Xi1_derived = sp.simplify(sp.diff(sp.log(T2_direct_pert), s).subs(s, 0))
expect_zero("Xi_1 derived from T^2 matches defect law", Xi1_derived - Xi1)

# Support-blindness of the defect itself: derive Xi_1 through the zeta-bearing
# support-loaded shape (line 56 proves T2_loaded == T2_direct) and confirm no zeta.
T2_loaded_pert = T2_loaded.subs({
    ZW: ZW*(1 + s*zetaZ), OmegaW2: OmegaW2*(1 + s*omegaW),
    chi0: chi0 + s*chi1, epsW: epsW + s*epsW1, deltaU: deltaU + s*deltaU1,
})
Xi1_loaded = sp.simplify(sp.diff(sp.log(T2_loaded_pert), s).subs(s, 0))
expect_zero("d/dzeta Xi_1 (support-loaded route)", sp.diff(Xi1_loaded, zeta))
```
Note: `T2_loaded` was built from `R_target_loaded`, which built from `Lam` (carries `OmegaW2`), `Mmix`, `Msupp` — all expressed in the same primitive symbols, so the `.subs` above perturbs every occurrence. If `eps_eta` also appears in `T2_loaded` it is ζ-independent and may be left unperturbed for this check (the F1 second check is about ζ only).

**Required change (Mathematica):** mirror with `D[Log[expr], s] /. s -> 0`:
```wl
sPar = s;  (* declare s as a real perturbation symbol in $Assumptions for this block *)
epsPert = (epsW + s*epsW1)*(1 - (2/11)*(deltaU + s*deltaU1)/(1 + deltaU + s*deltaU1));
t2DirectPert = (zW*(1 + s*zetaZ))*(1 + (chi0 + s*chi1))^2/((omegaW2*(1 + s*omegaW))*(1 - epsPert)^2);
xi1Derived = FullSimplify[(D[Log[t2DirectPert], s] /. s -> 0), Assumptions -> $Assumptions];
expectZero["Xi_1 derived from T^2 matches defect law", xi1Derived - xi1];
t2LoadedPert = t2Loaded /. {zW -> zW*(1 + s*zetaZ), omegaW2 -> omegaW2*(1 + s*omegaW),
  chi0 -> chi0 + s*chi1, epsW -> epsW + s*epsW1, deltaU -> deltaU + s*deltaU1};
xi1Loaded = FullSimplify[(D[Log[t2LoadedPert], s] /. s -> 0), Assumptions -> $Assumptions];
expectZero["d/dzeta Xi_1 (support-loaded route)", D[xi1Loaded, zeta]];
```
Add `s` (and ensure `zeta` is in scope) to the relevant `$Assumptions` Reals list for this block. (In the Mathematica script the second drift-transport block at line 55 redefines `$Assumptions` without `zeta`/the M-route symbols; if `t2Loaded` reuse causes trouble, recompute `t2Loaded` from the closed form `zW*(1+chi0)^2/(omegaW2*(1-eps)^2)` perturbed directly — it is equal by line 48 — rather than carrying the loaded route across the `$Assumptions` reset.)

**Claim manifest:**
- M1: \(\Xi_1 = \delta\ln\mathcal T^2/(\epsilon\lambda)\) evaluated on the closed form equals \(\zeta_Z - \omega_W + 2\chi_1/(1+\chi_0) + 2\epsilon_1/(1-\epsilon)\).
- M2: the same defect derived through the ζ-bearing support-loaded shape carries no ζ, i.e. \(\partial_\zeta\Xi_1 = 0\).

**Verification command:** verifier runs `redteam exec-sympy 181` and `redteam exec-mathematica 181`; new lines `Xi_1 derived from T^2 matches defect law = 0` and `d/dzeta Xi_1 (support-loaded route) = 0` must appear and scripts exit 0.

## Applied: F1

- files_changed:
  - `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.py`
  - `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl`
- summary: Added perturbation-derived `Xi_1` checks and the `d/dzeta Xi_1` support-blindness checks.
- deviation: Mathematica uses the closed-form `T^2` perturbation after the prior support-loaded reconstruction proof, per the consult resolution.

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.py:89-97`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl:65-69`

**Issue:** The "selected-branch identity" check (sympy line 97 / math line 69) is `expect_zero(Xi1 + eta1/(1-eps_eta) + R1)`, but `R1` (sympy 89–95) is defined as `-(Xi1) - eta1/(1-eps_eta)` written out term-by-term, so the residual is identically zero by construction. It never tests the paper's claim (notes §4) that the *directly-computed* \(\mathcal R_1 := \delta\ln R_{\rm target}/(\epsilon\lambda)\) equals that closed form.

**Required change (SymPy):** Insert immediately after the `R1 = sp.simplify(...)` definition (after line 95) and before the existing line-97 assertion, a derivation of `R1` from `R_target` using the shared perturbation map (including the `eta1` drift), then keep the existing line-97 check (now non-vacuous because `R1` is anchored). Reuse the `s` symbol from F1:
```python
Lam_pert = Lam.subs(OmegaW2, OmegaW2*(1 + s*omegaW))
eps_pert = (epsW + s*epsW1) * (1 - sp.Rational(2, 11) * (deltaU + s*deltaU1) / (1 + deltaU + s*deltaU1))
R_target_pert = Lam_pert * (1 - (eps_eta + s*eta1)) * (1 - eps_pert)**2 \
                / ((ZW*(1 + s*zetaZ)) * (1 + (chi0 + s*chi1))**2)
R1_derived = sp.simplify(sp.diff(sp.log(R_target_pert), s).subs(s, 0))
expect_zero("R_1 derived from R_target matches closed form", R1_derived - R1)
```
Leave the existing `expect_zero("selected-branch identity", Xi1 + eta1/(1 - eps_eta) + R1)` at line 97 in place — it now serves as a genuine cross-check between the two independently-anchored drifts. (Order F1's `s` definition before this block, or define `s` once at the top of the drift-transport block.)

**Required change (Mathematica):** mirror, inserting after the `r1 = FullSimplify[...]` definition (after line 68) and before the line-69 `expectZero`:
```wl
lamNormPert = lamNorm /. omegaW2 -> omegaW2*(1 + s*omegaW);
epsPert = (epsW + s*epsW1)*(1 - (2/11)*(deltaU + s*deltaU1)/(1 + deltaU + s*deltaU1));
rTargetPert = lamNormPert*(1 - (epsEta + s*eta1))*(1 - epsPert)^2/((zW*(1 + s*zetaZ))*(1 + (chi0 + s*chi1))^2);
r1Derived = FullSimplify[(D[Log[rTargetPert], s] /. s -> 0), Assumptions -> $Assumptions];
expectZero["R_1 derived from R_target matches closed form", r1Derived - r1];
```
Ensure `s, zW, omegaW2, epsEta, gConst, cSpeed, cs, a` and the drift symbols are in `$Assumptions` for this block (the line-55 `$Assumptions` reset drops `zW`/`omegaW2`/positivity — re-add what `lamNorm`/`rTarget` need, or recompute `lamNormPert`/`rTargetPert` so the FullSimplify can close).

**Claim manifest:**
- M3: \(\mathcal R_1 = \delta\ln R_{\rm target}/(\epsilon\lambda)\) evaluated on the closed form equals \(\omega_W - \eta_1/(1-\epsilon_\eta) - \zeta_Z - 2\chi_1/(1+\chi_0) - 2\epsilon_1/(1-\epsilon)\).

**Verification command:** `redteam exec-sympy 181` / `exec-mathematica 181`; new line `R_1 derived from R_target matches closed form = 0` appears and the retained `selected-branch identity = 0` is now backed by an anchored `R1`; scripts exit 0.

## Applied: F2

- files_changed:
  - `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.py`
  - `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl`
- summary: Derived `R_1` from the perturbed closed-form `R_target` before retaining the selected-branch identity check.
- deviation: none

## F3 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl:32-89`

**Issue:** The `.wl` reproduces the SymPy script's intermediate-variable choreography (`mMix/mSupp/sSupport/mTr/productLoaded/rTargetLoaded/t2Loaded`, then identical `eps1`/`theta1` log-derivatives and the same six `expectZero` in the same order) line-for-line; it is a transliteration, not an independent derivation. It also omits the SymPy `spoiled` adversarial control (sympy 61–67).

**Required change:** Make the Mathematica derivation reach the load-bearing results by a route that does not mirror the SymPy intermediate sequence, and add the missing negative control. Specifically:
1. Replace the term-by-term `mTr = mMix + mSupp; rTargetLoaded = productLoaded/mTr` construction with an independent check that the support-loaded route returns the ζ-free `rTarget`: e.g. assert directly that `FullSimplify[rTargetLoaded - rTarget] === 0` is *forced by* the algebraic identity, by instead `Reduce`/`Solve`-confirming the ζ-cancellation, or by computing `rTargetLoaded` from `sSupport` in closed form and showing the `sSupport` factor cancels. The end assertions (lines 47–48) may stay, but the path to `rTargetLoaded` should differ from the SymPy `mMix + mSupp` grouping.
2. Add the F1 and F2 perturbation-based derivation checks (`D[Log[t2DirectPert], s]`, `D[Log[rTargetPert], s]`) — these are inherently independent of the SymPy intermediate choreography.
3. Port the `spoiled` negative control: define `mSuppSpoiled = mSupp + bad*zeta*mMix` (with `bad` a nonzero real), `rTargetSpoiled = productLoaded/(mMix + mSuppSpoiled)`, and assert `(D[Log[rTargetSpoiled], zeta] /. bad -> 1)` is NOT identically zero (use a `fail` if `TrueQ[... === 0]`).

Keep the physics identical; only the derivation structure and the added control change. Do not alter the closed-form targets.

**Verification command:** `redteam exec-mathematica 181`; the `.wl` no longer reproduces the SymPy `mMix/mSupp/mTr/productLoaded` grouping verbatim, includes the F1/F2 derivation checks and the spoiled control, prints all `PASS`, and exits 0.

## Applied: F3

- files_changed:
  - `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl`
- summary: Reworked the Mathematica support-loaded route through `sSupport` cancellation and added the spoiled-packet negative control.
- deviation: none

## F4 — cosmetic banner mislabel (low; apply only if trivially safe)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.py:33`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl:26`

**Issue:** Both banners read `banner("STAGE 164 — COHERENT TRACKING-BRANCH DEFECT LAW")` though this is stage 181 (the Mathematica closing print already correctly says "Stage 181" at line 101).

**Required change:** Change both banner strings from `STAGE 164` to `STAGE 181`. Cosmetic only — no math impact. If for any reason this conflicts with another reference, leave it and append `## Blocked: F4`.

**Verification command:** transcripts show `STAGE 181 — COHERENT TRACKING-BRANCH DEFECT LAW`; scripts exit 0.

## Applied: F4

- files_changed:
  - `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.py`
  - `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl`
- summary: Corrected both audit banners from `STAGE 164` to `STAGE 181`.
- deviation: none

---

## Iteration 2 (delta) — F1 Mathematica M2 check is tautological

Verifier (`redteam/verifications/stage_181.md`) confirmed F2/F3/F4 RESOLVED and the F1 **primary**
M1 check (`Xi_1 derived from T^2 matches defect law`) genuine in BOTH engines. The ONLY defect:
the F1 **secondary** M2 check on the Mathematica side is vacuous.

**Target:** `mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl`
(the `d/dzeta Xi_1 (support-loaded route)` block, ~wl:101–103).

**Still wrong:** Codex set `t2LoadedPert = t2DirectPert;` (wl:101). `t2DirectPert` contains NO
`zeta`, so `D[xi1Loaded, zeta]` is identically 0 regardless of physics — the PASS is vacuous. The
SymPy mirror is genuine (it perturbs the zeta-bearing `T2_loaded`); the Mathematica side must also
exercise a genuinely zeta-bearing object, or the second-engine guard on the HIGH finding is empty.
Do NOT fall back to a zeta-free closed form (that reintroduces this tautology).

**Required change:** Build `t2LoadedPert` from a zeta-BEARING quantity. The clean route is via
`rTargetLoaded` (the spoiled control at ~wl:53–64 already proves `rTargetLoaded` retains `zeta`),
e.g.:
```wl
rTargetLoadedPert = rTargetLoaded /. {zW -> zW*(1 + s*zetaZ), omegaW2 -> omegaW2*(1 + s*omegaW),
  chi0 -> chi0 + s*chi1, epsW -> epsW + s*epsW1, deltaU -> deltaU + s*deltaU1};
t2LoadedPert = (lamNorm /. omegaW2 -> omegaW2*(1 + s*omegaW))/(omegaW2*(1 + s*omegaW))*
               (1 - epsEta)/rTargetLoadedPert;
xi1Loaded = FullSimplify[(D[Log[t2LoadedPert], s] /. s -> 0), Assumptions -> $Assumptions];
expectZero["d/dzeta Xi_1 (support-loaded route)", D[xi1Loaded, zeta]];
```
**Hard requirement:** the object whose `s`-log-derivative is taken MUST still contain `zeta` before
`D[..., zeta]` — verify `FreeQ[t2LoadedPert, zeta]` is `False` (add `zeta` to the local
`$Assumptions` Reals list if `FullSimplify` needs it). If `FullSimplify` won't close, reformulate
locally — do NOT raise the 600s cap and do NOT revert to a zeta-free form. Also remove the dead
`sPar = s;` line (wl:96) if present. Leave the SymPy M2 check and ALL other checks untouched.

**Verification:** Mathematica transcript still shows `d/dzeta Xi_1 (support-loaded route) = 0 / PASS`,
but `t2LoadedPert` is provably zeta-bearing prior to differentiation; both engines exit 0.

## Applied: Iteration 2

- files_changed:
  - `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl`
- summary: Replaced the zeta-free `t2LoadedPert = t2DirectPert` with an explicit zeta-bearing loaded perturbation and a `FreeQ` guard.
- deviation: Used local perturbed `mMix`/`mSupp`/support-product factors because the precomputed loaded ratio had already simplified away `zeta`.
