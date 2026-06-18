# Directive pathA_06 — GATE-A freeze stamp (register frozen family + stamp candidate_freeze_hash)

**Owner:** Codex (codes). Claude reviews afterward (transliteration-fidelity + no-residual scan) and records the
freeze in the decision/pre-reg scaffolding + commits.
**Status gate:** USER authorized `frozen: YES` on `decisions/11_pathA_gate_a_freeze_sheet.md` (2026-06-18). This is
pre-reg §M Stage-2. **This directive is the freeze stamp ONLY — it must NOT run any calibration, closed solve, or
extraction, and must NOT compute `R_norm`, `D0`, `P0`, or any residual.** Freeze-before-solve discipline
(decision-07 / `m1c_physical_run.py` precedent): the frozen packet is byte-reproducible and contains NO target
residual or pass/fail field.

## Objective
Two deliverables, both target-blind:

### 1. Register the frozen `homogeneous_isotropic_hooke_v1` S_Σ provider
In `src/stage1_solver/patha_static_balance.py`, add a new registered family implementing decision-11 §1 EXACTLY
(strict harmonic wall, 1 material scale `τ`, `g=0`, preferred radius `R_*=a`). Follow the existing
`SSigmaProvider` interface (the 8 methods `mu, T_w, T_w_R, T_w_RR, T_Omega, U, U_R, U_RR`) and registry pattern
(`_S_SIGMA_REGISTRY`, an `SSigmaSpec` factory analogous to `smooth_positive_placeholder` / `patha_static_mms`).

Frozen forms (with `G=c=c_s=1`; parameters: `tau`, `a`, `w_min`, `w_max`; `R_star=a`):

| method | expression |
|---|---|
| `mu(R,w)`      | `τ`                  (constant) |
| `T_w(R,w)`     | `τ`                  (constant) |
| `T_w_R(R,w)`   | `0` |
| `T_w_RR(R,w)`  | `0` |
| `T_Omega(R,w)` | `τ/a²`               (constant; isotropic tie) |
| `U(R,w)`       | `½ (τ/a²) (R − a)²` |
| `U_R(R,w)`     | `(τ/a²) (R − a)` |
| `U_RR(R,w)`    | `τ/a²`               (constant) |

- `g=0`: NO quartic / anharmonic term anywhere. No `w`-dependence in any modulus (homogeneous). No
  sinusoidal-in-`w` placeholder modulations.
- Admissibility: require `τ>0` (and `a>0`). `U` is bounded below with a single well at `R=a`.
- Return tensors broadcasting over `R`,`w` exactly like the existing providers (e.g. constants via
  `torch.full_like`/`zeros_like` so shapes/devices/dtype match).
- Add a focused unit test (next to the existing patha_static_balance tests) that (a) checks the analytic
  derivatives `U_R, U_RR, T_w_R, T_w_RR` against finite differences of `U`/`T_w`, and (b) pins the harmonic values
  at a couple of `(R,w,τ,a)` points. Keep it within the repo's existing test layout.

### 2. Build the GATE-A freeze module + stamp the hash
Add a new module (e.g. `src/stage1_solver/patha_gate_a_freeze.py`) modeled on `m1c_physical_run.py`'s
`build_freeze_sheet()` / `write_freeze()` / `load_freeze()` (REUSE the same canonical-serialization + hash
convention verbatim: `json.dumps(..., sort_keys=True, separators=(",", ":"))` → `hashlib.sha256(...).hexdigest()`).
It assembles the canonical Path-A GATE-A freeze sheet covering **exactly** the `decisions/11` §8 hashed basis, then
stamps `candidate_freeze_hash` and writes the frozen packet.

Canonical freeze-sheet content (the §8 basis — encode each field; values as data/strings, not computed numerics):
- `schema`: `"pathA_gate_a_decision_11_freeze/v1"`
- `parent_action_status`: `"path_A_promoted_S_Sigma"`
- `branch_identity`: decision record path = `software/stage1_solver/decisions/11_pathA_gate_a_freeze_sheet.md`;
  `gate="GATE-A"`; `target_blind=true`; `no_post_residual_refit=true`;
  `firewalled_simulation_tree="not_touched"`; `physical_export_permitted_guard="not_touched"`.
- `family`: `"homogeneous_isotropic_hooke_v1"` + the 4 forms (as the expression strings above) + the constitutive
  derivatives + the ties (`mu_Σ=τ, T_w=τ, T_Ω=τ/a²`) + `g=0` + domains (`w∈[0,L]`, `R>0`) + units
  (`[μ]=[T_w]=1`, `[T_Ω]=L⁻²`, `[U]=1`) + `R_star="a"`.
- `geometry`: `a=1`, `L=1.85` (`L_exact="37/20"`), `boundary_class` (mouth Dirichlet `R0(0)=a`; natural
  zero-traction exit `T_w R0'(L)=0`, `Y_L=0`; no hard cap; `R0(w)>0`); modal BCs `η(0)=0`, `T_w η'(L)=0`.
- `source_port`: `mhat0=1`, `S_port=1`; `gauge` convention; constants `{G=1,c=1,c_s=1}`.
- `calibration_objective`: anchor `R_norm=0` ⇔ `P0 = N0/D0 = 54/5` under `mhat0²·S_port=1`; the **anchor constant
  string** `"54*G*c_s^5/(5*a^5*c^5)"` (this IS hashed — it is what we calibrate TO, known a priori, per §8
  "Hashed"); root-finder = deterministic 1-DOF root-find on `τ`; the **stable-side branch selection `D0>0`**;
  tolerances + mesh/grid-convergence ladder (carry the decision-11 §3/§7 + decision-10 settings); channel/mode
  selection (`l=2`, lowest-positive `L₂` eigenmode, `∫χ²=1`, `∫χ>0`).
- `extraction_formulas`: the §5b formula **strings** EXACTLY (`D0=K-B0-Z0`, `D2=-(M+B2+Z2)`, `D4=-(B4+Z4)`,
  `P0=N0/D0`, `P2=(D0*N2-2*D2*N0)/D0^2`, `P4=(D0^2*N4-2*D0*(D2*N2+D4*N0)+3*D2^2*N0)/D0^3`,
  `R_norm=mhat0^2*S_port*N0/D0-54*G*c_s^5/(5*a^5*c^5)`, `R_pole=D0*(B4+Z4)-3*(M+B2+Z2)^2`,
  `K=tau*kappahat`). These are SPEC strings; do NOT evaluate them.
- `source_revision`: `git_head`, `library_versions`, and `source_file_sha256` over the load-bearing files at
  freeze (at minimum: `patha_static_balance.py`, the new freeze module, the v2 fixtures
  `research/pde_audit/scripts/stage_v2_21_branch_extraction_fixture.py` and
  `…stage_v2_22a_profile_to_coefficient_adapter.py`; mirror the file-set m1c used where applicable).
- **EXCLUDED (must NOT appear):** any computed `R_norm`/`D0`/`P0`/`N0` numeric value, any target residual, any
  pass/fail field, any held-out (`R_pole`/`P2`/`P4`) target value.

Output: write `software/stage1_solver/frozen/pathA_gate_a/<freeze_hash>/freeze_sheet.json` and
`freeze_hash.txt` (the bare hash + newline). Provide a `load_freeze(hash)` that recomputes the hash from the
written sheet and asserts it matches (the m1c integrity check).

## Acceptance criteria (Codex iterates until ALL pass, exit 0)
1. The new family resolves via `resolve_s_sigma({"family":"homogeneous_isotropic_hooke_v1", "parameters":{...}})`
   and appears in `registered_s_sigma_families()`; its unit test passes.
2. The freeze module runs and prints the `candidate_freeze_hash` + the frozen packet path.
3. **Byte-reproducibility:** running the stamp twice yields the IDENTICAL hash and identical
   `freeze_sheet.json` bytes (HEAD unchanged between runs).
4. **No-residual scan:** a programmatic scan of the written `freeze_sheet.json` confirms it contains none of a
   computed `R_norm`/`D0`/`P0`/`N0` numeric value, nor any `pass`/`fail`/`residual`/target-value field. (The
   formula STRINGS and the anchor constant string are expected and allowed.)
5. **Firewall intact:** no import from or write under `research/pde_audit/simulation/`; the
   `physical_export_permitted` guard is not referenced/altered; all writes stay under `software/stage1_solver/`.
6. The existing patha_static_balance / chunk-1a–1c tests still pass (no regression to the placeholder families).

## Constraints
- Any script you run is wrapped `timeout 600`; a timeout (exit 124) is a failure → reformulate, never raise the
  cap. (This stamp is sub-second; the cap is just policy.)
- Do NOT solve, calibrate, or compute physical observables here. Freeze stamp + provider registration ONLY.
- Sandbox = workspace-write; git is the safety net. Do not `git add`/commit — the orchestrator commits the freeze.

## Report back
The `candidate_freeze_hash`, the frozen packet path, the byte-repro confirmation, the no-residual-scan result, the
list of files created/modified, and the unit-test + regression-test results.
