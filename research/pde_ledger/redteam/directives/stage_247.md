---
unit_id: 247
batch: VIII.1
created_at: 2026-06-03T00:00:00-06:00
findings_count: 6
stop_cold: null
applied: true
applied_at: 2026-06-03T10:16:18-06:00
findings_applied: 6
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 247

F1 (`paper_misalignment`) has been RESOLVED by the user (2026-06-03): correct the single NOTES-file typo to the script value — see the `## RESOLVED — F1` block below. Apply F1's authorized notes edit AND F2-F6 (all script-side) in order, then append an `## Applied: F<n>` block under each with `files_changed`, `summary` (one sentence), and `deviation` (or "none"). Do NOT edit the published paper card `paper/stages/stage_247.tex` (it is clean — it does NOT contain the disputed value) or any other prose document; the ONLY authorized prose edit is the single notes line named in `## RESOLVED — F1`. If a required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question and continue.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

Do NOT introduce new features, refactors, or stylistic changes beyond what each finding names.

## RESOLVED — F1 (paper_misalignment, USER-APPROVED 2026-06-03)

**Subtype:** value_mismatch — **NOTES-ONLY; the published paper card is UNAFFECTED.**

**Correction (apply this one line verbatim):**
- `notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md:406` — change `\Delta=210.17750000,` → `\Delta=142.17750000,`

**Why (verified by orchestrator):** The published card `paper/stages/stage_247.tex` is only 93 lines and contains NO occurrence of `210.1775` anywhere — the audit agent's `stage_247.tex:407` citation was a misattribution. The disputed value lives ONLY in the notes file (line 406). The notes' own formula `\Delta=\Omega_U^2\Omega_W^2-R_{\rm mix}^2 = 9·16 − 1.35² = 142.1775` and the adjacent notes value `D_0=3.76481862` (consistent only with Δ=142.1775) both confirm `210.17750000` is a notes transcription typo. The SymPy script asserts `142.1775` (py:235); the new `.wl` (F6/M6) independently computes `Δ=142.1775` from the same formula → cross-engine corroboration. User direction = (a): correct the notes to the script value.

**Scope:** Edit ONLY notes:406. Do NOT touch `paper/stages/stage_247.tex`, the script's `142.1775` literal, or any other line. This is a Codex-applied notes edit; Claude reviews.

## Applied: F1

- files_changed:
  - `notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md`
- summary: Corrected the notes-only benchmark typo from `\Delta=210.17750000` to `\Delta=142.17750000`.
- deviation: none

---

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py:122` (and add new numeric checks in Section 6, around py:228-238)

**Issue:** `Lvar_from_W` (py:113) is the exact algebraic inverse of the line defining `W_sess` (py:110), so `assert sp.simplify(Lvar_from_W - Lvar) == 0` (py:122) reduces to `sqrt(Lvar**2) == Lvar` for positive `Lvar` and can never fail. It does not verify the genuine Stage-244 inversion that recovers the *numeric* `Lvar` from the recorded `W_sess`.

**Required change:**
1. Keep the symbolic round-trip line py:122 only if you also harden it numerically; otherwise it is acceptable to leave py:122 as documentation. The real fix is to ADD, in Section 6 after `S_soft` is computed (after py:217), substantive numeric asserts:
   ```python
   # F2: the inverted Lvar reproduces the recorded benchmark work scalar, and
   #     matches the paper-stated Lvar(r_soft) = 20.01677473.
   Wsess_from_Lvar = sp.N(W_sess.subs({**subs_soft, Lvar: Lvar_soft}), 16)
   assert abs(float(Wsess_from_Lvar) - float(Wsess_obs)) < 1e-7
   assert abs(float(Lvar_soft) - 20.01677473) < 1e-6
   ```
   (`Wsess_obs` is already defined at py:205; `Lvar_soft` at py:213-216.)

**Self-test:** No `sp.diff` involved. `W_sess` depends on `Lvar`, `eta_leak`, `mu_w`, `q`, `rho0`, `lam`; substituting all of them yields the literal `1.51632107` (by construction of `Lvar_soft`, which was solved from `Wsess_obs`) — so the `Wsess_from_Lvar` check confirms the inversion is internally consistent, while the `Lvar_soft ~ 20.01677473` check pins it to the paper figure and IS falsifiable (it would fail if the inversion constants were wrong). Both satisfiable: printed `Lvar(session) = 20.01677472685125`.

**Verification command:** `redteam exec-sympy 247` — confirm the two new asserts appear in Section 6 and the script exits 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py`
- summary: Added numeric assertions that the inverted `Lvar_soft` reproduces `W_sess` and matches the paper benchmark value.
- deviation: none

## F3 — script_missing_paper_claim

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py` — Section 3 (after py:131) and Section 4 (after py:147), plus numeric checks in Section 6.

**Issue:** The paper makes `D_UV >= 0` (appendix `eq:app-part08-main-uv-drain`), `g(r) >= 2/pi` and `M_sigma >= 0` on the Session-I branch (notes 2.3) load-bearing for the lowering inequality, but the script only prints these and never asserts them.

**Required change:**
1. In Section 3, after `E_UV` is defined (after py:131), add a symbolic nonnegativity assert for the drain. `D_UV = a_V*chi_UV**2*f_U**2/Delta_UV**2`; with `a_V, chi_UV, f_U` positive and `Delta_UV` real-nonzero, this is manifestly `>= 0`. Use `sp.ask`/`is_nonnegative` style or an explicit structural check:
   ```python
   # F3: Stage-245 drain is nonnegative (square / square, a_V > 0).
   assert sp.ask(sp.Q.nonnegative(D_UV)) in (True, None)  # structural: numerator and denom are squares times a_V>0
   # Robust fallback (preferred): assert the simplified sign on a numeric probe.
   D_UV_probe = D_UV.subs({a_U: sp.Float("2.0"), a_V: sp.Float("1.5"),
                           chi_UV: sp.Float("0.7"), f_U: sp.Float("0.4")})
   assert float(D_UV_probe) >= 0
   ```
   If `sp.ask` returns `None` in your SymPy build, rely on the numeric probe assert (which is non-vacuous: it gives a positive literal).
2. In Section 6, after `M_sigma_num` is computed (py:211), add branch-positivity checks at `r_soft`:
   ```python
   # F3: source-response and g-bound positivity on the Session-I branch.
   g_soft = sp.N(g_r.subs(subs_soft), 16)
   assert float(g_soft) >= float(2 / sp.pi)          # g(r_soft) >= 2/pi
   assert float(g_soft) < float(subs_soft[rF1])      # g(r_soft) < rF1
   assert float(M_sigma_num) >= 0                     # M_sigma(r_soft) >= 0
   ```
   (`g_r` is defined at py:143; `rF1`, `r_soft` are in `subs_soft`.)

**Self-test:** No `sp.diff`. `g_r` depends on `r, r_sigma, a0, b0` (all in `subs_soft`), so `g_soft` is a concrete number — printed `g(r)` evaluates with `a0=2.2, b0=-0.6, r_sigma=0.8, r=0.18` to a value `>= 2/pi` (a0>0 raises g above 2/pi) and `< rF1=1.778` per the notes' stated branch fact; both asserts satisfiable. `M_sigma_num` printed = `0.1838612 >= 0`. `D_UV_probe` with the probe values = `1.5*0.49*0.16/(2*1.5-0.49)^2 = 0.1176/6.3001 = 0.01867 >= 0`, non-vacuous positive literal.

**Verification command:** `redteam exec-sympy 247` — confirm the drain-nonneg, g-bound, and M_sigma-nonneg asserts appear and the script exits 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py`
- summary: Added drain nonnegativity checks and Session-I `g(r)`/`M_sigma` positivity assertions.
- deviation: none

## F4 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py` — Section 5 (keep py:172) and add a numeric inequality check in Section 6.

**Issue:** The identity assert at py:172 only re-expands the definition of `V_eff` (py:166); it cannot fail. The paper's actual deliverable for claim (5) is the inequality `V_eff <= V_short` (appendix item 6), which is never checked.

**Required change:**
1. Leave py:172 as-is (valid structural check).
2. In Section 6, after the benchmark numbers are built (after the F3 block), add a numeric lowering-inequality assert on the Session-I slice using the paper's literal `lambda_L = 0.26971918`:
   ```python
   # F4: lowered potential is below the baseline on the benchmark slice (lowering theorem).
   lambda_L_paper = sp.Float("0.26971918")
   Veff_session = sp.N(Vshort_num - lambda_L_paper * S_soft - Wsess_obs - UVdrop_obs - M_sigma_num, 16)
   assert float(Vshort_num - Veff_session) >= 0     # V_eff <= V_short (gap >= 0)
   ```
   This depends on F5's `lambda_L_paper` literal; define it once and reuse.

**Self-test:** No `sp.diff`. `Vshort_num - Veff_session = lambda_L_paper*S_soft + Wsess_obs + UVdrop_obs + M_sigma_num`, a sum of positive numbers (`0.0838 + 1.5163 + 0.2106 + 0.1839 = 1.9946 >= 0`) — satisfiable and non-vacuous (would fail if any packet term were wrongly signed).

**Verification command:** `redteam exec-sympy 247` — confirm the `V_eff <= V_short` numeric assert appears and the script exits 0.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py`
- summary: Added the Session-I numeric lowering inequality using the paper literal `lambda_L = 0.26971918`.
- deviation: none

## F5 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py:219-220,237` plus new pinning asserts in Section 6.

**Issue:** `lambda_L_soft` is solved (py:219) precisely so that `Vrebuild_soft` (py:220) equals `Veff_obs`; the assert at py:237 (`abs(Vrebuild_soft - Veff_obs) < 1e-10`) is therefore true by construction for any inputs and cannot detect an error. The genuine benchmark claim is that the independently derived quantities plus the paper's `lambda_L = 0.26971918` reproduce `1.74701126`.

**Required change:**
1. Add asserts pinning each independently derived benchmark quantity to its paper-stated value (place after `S_soft` / `lambda_L_soft` are computed, py:217-219):
   ```python
   # F5: pin independently derived benchmark quantities to the paper figures.
   assert abs(float(Vshort_num) - 3.74163698) < 1e-6
   assert abs(float(M_sigma_num) - 0.18386120) < 1e-6
   assert abs(float(S_soft) - 0.31069599) < 1e-6
   assert abs(float(lambda_L_soft) - 0.26971918) < 1e-6   # solved lambda_L matches paper value
   ```
2. Replace the tautological closure check at py:237 with a FORWARD decomposition using the paper's literal `lambda_L` (define `lambda_L_paper = sp.Float("0.26971918")` once; reuse in F4):
   ```python
   # F5: forward benchmark decomposition with the paper's lambda_L (falsifiable closure).
   Veff_forward = sp.N(Vshort_num - lambda_L_paper * S_soft - Wsess_obs - UVdrop_obs - M_sigma_num, 16)
   assert abs(float(Veff_forward) - float(Veff_obs)) < 1e-6
   ```
   Keep `lambda_L_soft > 0` (py:238).

**Self-test:** No `sp.diff`. Hand-check the forward closure: `0.26971918*0.31069599 = 0.0838007`; `3.74163698 - 0.0838007 - 1.51632107 - 0.21064278 - 0.18386120 = 1.7470113`, within `1e-6` of `Veff_obs = 1.74701126`. The `lambda_L_soft ~ 0.26971918` pin is consistent because `lambda_L_soft` solves exactly to `(3.74163698 - 1.51632107 - 0.21064278 - 0.18386120 - 1.74701126)/0.31069599 = 0.0838007/0.31069599 = 0.2697192`. All asserts satisfiable and falsifiable (a wrong S_soft, M_sigma, or V_short breaks the pin or the forward closure).

**Verification command:** `redteam exec-sympy 247` — confirm the four pin asserts and the forward closure assert appear, py:237 no longer back-substitutes the solved lambda_L, and the script exits 0.

## Applied: F5

- files_changed:
  - `scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py`
- summary: Pinned independent benchmark quantities and replaced the solved-lambda closure assert with the forward paper-lambda decomposition.
- deviation: none

## F6 — missing_verification_script

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_mathematica_audit.wl` (NEW file; `.wl` lives in `mathematica/`)

**Issue:** Stage 247 is checkpoint:False, is_status_only_candidate:False → both engines required. No `.wl` exists. All claims are independently verifiable in native Mathematica, so a second-engine script is required. The `.wl` MUST be an independent re-derivation using native Mathematica primitives, NOT a line-by-line port of the `.py`.

**Required change:**
Create the `.wl` with `expectZero` / `expectTrue` / `expectApprox` helpers (each `Exit[1]` on failure) and an independent decomposition:
- Build the reduced stiffness matrix `Kred = {{Kstar,-GU,-GW},{-GU,OmU2,-Rmix},{-GW,-Rmix,OmW2}}` and take `Inv = Inverse[Kred]`. Extract the susceptibility entries from `Inv` directly and confirm against the paper's closed forms with `Simplify[Inv[[i,j]] - chiClosed]` (this is a DIFFERENT route than the SymPy `K_red.inv()` element-by-element compare only if you derive `chiClosed` from the paper formulas independently; do so).
- Use `Limit[g[r], r -> Infinity]` for `g_inf` (native limit; the SymPy script uses `sp.limit` — fine, but structure the `g[r]` definition from the moment formula `(2/Pi)(1 + a0 s/3 - b0 s/15)` rather than copying intermediate variable names).
- Verify `M_sigma`, the lowering identity `Simplify[(Vshort - Veff) - packetSum] == 0`, the drain nonnegativity, the `g`/`M_sigma` branch positivity, and the Session-I benchmark numerically.

**Claim manifest** (the new `.wl` must independently verify each):
- **M1** — Susceptibilities are the entries of `Inverse[Kred]`: `chi_qq = 1/D0`, `chi_qU = PU/(Delta D0)`, `chi_qW = P/(Delta D0)`, `chi_UU = (Kstar OmW2 - GW^2)/(Delta D0)`, `chi_UW = (Kstar Rmix + GU GW)/(Delta D0)`, `chi_WW = (Kstar OmU2 - GU^2)/(Delta D0)`, with `Delta = OmU2 OmW2 - Rmix^2`, `D0 = Kstar - Q/Delta`, `Q = GU^2 OmW2 + 2 GU GW Rmix + GW^2 OmU2`, `P = OmU2 GW + Rmix GU`, `PU = GU OmW2 + Rmix GW`. Verify via `Simplify[Inverse[Kred][[i,j]] - chiClosed] == 0`.
- **M2** — Leakage/work factoring: with `E0 = 16 eta_leak Lvar/Pi^2`, `S_leak = 8 Sqrt[2] eta_leak mu_w rho0 Lvar/(Pi^(5/2) lam^3)` equals `Sqrt[2] E0 mu_w rho0/(2 Sqrt[Pi] lam^3)`, and `W_sess = 512 eta_leak^2 mu_w q rho0 Lvar^2/(Pi^4 lam^2)` equals `2 E0^2 mu_w q rho0/lam^2`. Verify via `Simplify[... ] == 0`.
- **M3** — `Limit[g[r], r->Infinity] == 2/Pi` where `g[r] = (2/Pi)(1 + a0 s[r]/3 - b0 s[r]/15)`, `s[r] = r_sigma^2/(r^2 + r_sigma^2)`.
- **M4** — Drain nonnegativity: `D_UV = a_V chi_UV^2 f_U^2/(a_U a_V - chi_UV^2)^2 >= 0` for `a_V > 0` (square over square). Verify with `Reduce`/`Simplify` under positivity assumptions or a numeric probe giving a positive literal.
- **M5** — Lowering identity: with `Veff = Vshort - lambda_L S_leak - lambda_W W_sess - E_UV - M_sigma`, `Simplify[(Vshort - Veff) - (lambda_L S_leak + lambda_W W_sess + E_UV + M_sigma)] == 0`.
- **M6** — Session-I baseline numerics: `Delta = 142.1775` (NOTE: equals OmU2 OmW2 - Rmix^2 = 9*16 - 1.35^2; the paper card's 210.1775 is under user review per F1 — use the formula-correct value, do NOT hardcode 210.1775), `D0 ~ 3.76481862`, `Vshort(r_soft) ~ 3.74163698` at the recorded params.
- **M7** — Imported benchmark values: `M_sigma(r_soft) ~ 0.18386120` (and `>= 0`), `Lvar(r_soft) ~ 20.01677473`, `S_leak(r_soft) ~ 0.31069599`, `g(r_soft) >= 2/Pi` and `g(r_soft) < rF1`.
- **M8** — Forward benchmark decomposition with the paper's literal `lambda_L = 0.26971918` (forward, NOT solved): `Vshort - 0.26971918 S_leak - W_sess - DeltaE_UV - M_sigma ~ 1.74701126`, and `lambda_L = 0.26971918 > 0`.

**Self-test (new .wl):** M3 `Limit` is over `r`, on which `g[r]` genuinely depends (via `s[r]`) → not a vacuous limit. M4/M7 positivity checks are over variables the expressions contain. M8 forward arithmetic hand-checked above (1.7470113). No `D[expr, var]` with `var` absent from `expr` is prescribed. M6 deliberately uses the formula-correct `Delta=142.1775` so the `.wl` does not bake in the disputed paper value.

**Verification command:** After Codex applies, the verifier runs `redteam exec-mathematica 247` and confirms the new `.wl` exists at the `mathematica/` path, all `expect*` checks pass, and the script exits 0. The route must use `Inverse`/`Limit`/`Simplify` natively (independent decomposition), not echo the SymPy variable choreography.

## Applied: F6

- files_changed:
  - `mathematica/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_mathematica_audit.wl`
- summary: Created the independent Mathematica audit covering M1-M8 with native inverse, limit, positivity, identity, and benchmark checks.
- deviation: none
