---
unit_id: 097
batch: IV.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage097_single_normalization_defect.md]
  paper_appendix: present
---

# Audit unit 097 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_097.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage097_single_normalization_defect.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows around §"Retarded 2.5PN collapse", eqs. 220–249, 264–289, and result anchor MTDC-T8.2 at line 1176)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage097_single_normalization_defect_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage097_single_normalization_defect_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage097_single_normalization_defect_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage097_single_normalization_defect_mathematica_audit.txt`

## What the paper claims

Stage 097 is a geometry-lane firewall ledger step that shows the actual isotropic passive/outgoing grouped-`P2` one-pole quadrupole branch collapses the entire reduced 2.5PN normalization problem to a single scalar defect `N_Q`. The card body fixes the conservative module `Yhat_Q^cons = 3/4 + (1/4)/(1-omega^2/Omega_Q^2)` and quotes the bottom line: "The isotropic passive/outgoing one-pole branch reduces to `N_Q = Kbar_0/Kbar_0^target`." The notes are the authoritative source and enumerate the concrete deliverables: even moments `Kbar_2 = Kbar_0/(4 Omega_Q^2)`, `Kbar_4 = Kbar_0/(4 Omega_Q^4)`; odd coefficient `Gammabar_5 = 9 Kbar_0/(32 Omega_Q^5)`; GR target `Kbar_0^target = 64 G Omega_Q^5/(45 c^5) = 54 G c_s^5/(5 a^5 c^5)` under `Omega_Q = 3 c_s/(2a)`; `Gammabar_5^target = 2G/(5c^5)`; and the single-defect conclusion `R_0 = R_2 = R_4 = R_5 = N_Q - 1`, i.e. the branch satisfies the full reduced 2.5PN theorem iff `N_Q = 1`. The appendix (eqs. \eqref{eq:app-part04-kbar-even-moments}, \eqref{eq:app-part04-kbar0-target}, \eqref{eq:app-part04-NQ-def}) reproduces these verbatim.

## What the script claims to verify

Both scripts verify the full deliverable chain symbolically over the positive-real domain `(G, c, c_s, a, K0, Omega, N_Q)`: (1) the odd coefficient `Gamma5 = 9 K2^(5/2)/K0^(3/2)` reduces to `9 K0/(32 Omega^5)`; (2) the GR target `K0_target = 64 G Omega^5/(45 c^5)` reduces under `Omega = 3 c_s/(2a)` to `54 G c_s^5/(5 a^5 c^5)`; (3) `Gamma5_target` reduces to `2G/(5c^5)`; (4) substituting `K0 := N_Q K0_target` makes each of the four defects `R0, R2, R4, R5` equal `N_Q - 1`. The SymPy script derives the moments by direct definition; the Mathematica script derives them independently by extracting `SeriesCoefficient` of the conservative module `k0*(3/4 + (1/4)/(1-w^2/omegaQ^2))`.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `Kbar_2 = Kbar_0/(4 Omega_Q^2)` (notes §1, app. eq. 231) | py L27 / wl L41 (series-derived, asserted `==0`) | match |
| `Kbar_4 = Kbar_0/(4 Omega_Q^4)` (notes §1, app. eq. 233) | py L28 / wl L42 (series-derived, asserted `==0`) | match |
| `Gammabar_5 = 9 Kbar_0/(32 Omega_Q^5)` (notes §2) | py L29–31 / wl L45–51 | match |
| `Kbar_0^target = 64 G Omega_Q^5/(45 c^5)` (notes §3, app. 240) | py L34 / wl L53 | match |
| geom reduction `= 54 G c_s^5/(5 a^5 c^5)`, `Omega_Q=3 c_s/(2a)` (notes §3, app. 242–244) | py L35–37 / wl L54–59 (asserted `==0`) | match |
| `Gammabar_5^target = 2G/(5c^5)` (notes §3) | py L41–42 / wl L65–66 (asserted `==0`) | match |
| `R_0=R_2=R_4=R_5=N_Q-1` (notes §4, app. 249–251) | py L46–54 / wl L69–87 (four asserted `==0`) | match |

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 31 | `simplify(Gamma5 - 9 K0/(32 Omega^5)) == 0` | odd coeff (notes §2) | yes |
| A2 | sympy | 37 | `simplify(K0_target_geom - 54 G c_s^5/(5 a^5 c^5)) == 0` | geom target (notes §3) | yes |
| A3 | sympy | 42 | `simplify(Gamma5_target - 2 G/(5 c^5)) == 0` | target odd coeff (notes §3) | yes |
| A4–A7 | sympy | 51–54 | `simplify(R{0,2,4,5} - (N_Q-1)) == 0` | single defect (notes §4) | yes |
| B1 | math | 41 | `expectZero[k2 series - k0/(4 omegaQ^2)]` | even moment K2 | yes |
| B2 | math | 42 | `expectZero[k4 series - k0/(4 omegaQ^4)]` | even moment K4 | yes |
| B3 | math | 51 | `expectZero[gamma5 - gamma5Expected]` | odd coeff | yes |
| B4 | math | 59 | `expectZero[k0TargetGeom - 54 ...]` | geom target | yes |
| B5 | math | 66 | `expectZero[gamma5Target - 2G/(5c^5)]` | target odd coeff | yes |
| B6–B9 | math | 84–87 | `expectZero[r{0,2,4,5} - (nQ-1)]` | single defect | yes |

Every script-side assertion traces to a specific notes/appendix deliverable. No orphaned checks; no checks unmentioned by the paper.

## Findings

None.

### Adversarial attacks attempted (all failed to break the scripts)

- **Tautology probe on R_i (A4–A7 / B6–B9):** The substitution `K0 := N_Q K0_target` is applied to the *actual-branch* coefficients only; `K2_target/K4_target/Gamma5_target` are built from the unsubstituted `K0_target`. So `R2 = K2(N_Q K0_target)/K2_target - 1 = N_Q - 1` is a genuine ratio that would expose a wrong proportionality (e.g. if `K2` had a different `Omega` power than `K2_target`, the `Omega` factors would not cancel and the residual would carry `Omega`). The checks can fail. Not tautological.
- **Independence of the two engines (transliteration probe):** SymPy defines `K2 = K0/(4 Omega^2)` directly; Mathematica derives `k2` by `SeriesCoefficient[k0*(3/4 + (1/4)/(1-w^2/omegaQ^2)), {w,0,2}]`. The Mathematica route reconstructs the moments from the pole module rather than re-stating the closed form, then asserts equality to `k0/(4 omegaQ^2)` (wl L41). That is an independent re-derivation, not a port. Similarly `Gamma5_target` uses the `9 k0/(32 omegaQ^5)` form in `.wl` (L65) but the `9 K2^(5/2)/K0^(3/2)` form in `.py` (L41) — different algebra reaching the same number.
- **Symbol-domain probe:** All symbols are `positive=True, real=True` (py L23–24, L45) / `> 0` Reals (wl L29–31). The half-integer powers `K2^(5/2)` and `K0^(3/2)` (py L29, L41) require a positive base to simplify to a single real branch; positivity is physically justified (`K0`, `Omega` are positive normalization/frequency scales) and matches the paper's setup. No domain error; positivity is needed and warranted, not a hidden over-assumption.
- **Series-coefficient parity (wl B1/B2):** `Yhat_Q^cons` is even in `w` (only `w^2` appears), so the order-2 and order-4 series coefficients exist and are nonzero — the `expectZero` differences against `k0/(4 omegaQ^2)` and `k0/(4 omegaQ^4)` are substantive, not trivially `0=0`.
- **Constant probe:** Every literal (`3/4`, `1/4`, `9`, `32`, `64`, `45`, `54`, `2/5`, `3/2`, `5/2`) matches the notes and appendix exactly. The `Gamma5 = 9 K0/(32 Omega^5)` "expected" form (py L30) is itself re-derived from `9 K2^(5/2)/K0^(3/2)` at A1, so it is not a hardcoded answer checked against itself.

## Independent-derivation check (Mathematica)

The `.wl` is an independent re-derivation, not a transliteration. The SymPy script writes the moments as bare closed forms (`K2 = K0/(4 Omega**2)`), whereas the Mathematica script builds them as Taylor coefficients of the physical conservative module (`SeriesCoefficient[kbarCons[w], {w,0,2}]`, wl L39–40) and only then asserts equality to the closed form (wl L41–42). The actual-branch R_i residuals are also built via the series route in `.wl` (k2Actual/k4Actual at wl L70–71) versus direct `.subs` in `.py` (L48–49). The two engines use different choreography (series extraction vs. direct substitution) and a different route for `Gamma5_target` (`9 k0/(32 omegaQ^5)` vs `9 K2^(5/2)/K0^(3/2)`), satisfying the second-engine independence policy.

## Engine cross-check

Both transcripts agree on every reported value:
- `K2 = K0/(4 Omega^2)` (sympy) ↔ `k0/(4*omegaQ^2)` (math)
- `K4 = K0/(4 Omega^4)` ↔ `k0/(4*omegaQ^4)`
- `Gamma5 = 9*K0/(32*Omega^5)` ↔ `(9*k0)/(32*omegaQ^5)`
- `K0_target = 64*G*Omega^5/(45*c^5)` ↔ `(64*gConst*omegaQ^5)/(45*cLight^5)`
- geom reduction `54*G*c_s^5/(5*a^5*c^5)` ↔ `(54*cSound^5*gConst)/(5*aRad^5*cLight^5)`
- `Gamma5_target = 2*G/(5*c^5)` ↔ confirmed `=0` residual
- `R0=R2=R4=R5 = N_Q - 1` ↔ `-1 + nQ` (all four)

All `expectZero`/`assert` residuals are exactly `0`. `engines_agree: true`.

## Verdict justification

Verdict is **clean**. I read the card, notes, and the part-IV appendix rows before opening the scripts; the script's verified claim (the K2/K4/Gamma5 reductions, the GR target, and the four-way R_i = N_Q - 1 collapse) matches the paper's stated deliverables exactly, with every constant agreeing. The R_i checks are non-tautological (target coefficients are built from the unsubstituted target K0, so a wrong power of Omega would survive), the positivity assumptions are physically warranted and required for the half-integer powers, the Mathematica engine independently re-derives the moments via series extraction, and both engines agree on every value with zero residuals. Output `.txt` files are newer than their scripts (see Value Reconciliation). Attacks on tautology, transliteration, symbol domain, parity, and constants all failed to find a defect.

## Value Reconciliation (pass-2 augmentation)

Output freshness: both `.txt` outputs are NEWER than their scripts — sympy `.py` mtime 2026-05-27 11:15:23 vs `.txt` 2026-05-27 14:28:41; math `.wl` mtime 2026-05-27 11:15:47 vs `.txt` 2026-05-27 14:30:31. `outputs_fresh: true`; no `stale_output` finding.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `K2 = K0/(4 Omega^2)` | py L27/L56, out L1; wl L41, out L9 | notes L34 `Kbar_2 = Kbar_0/(4 Omega_Q^2)`; app. eq. 231 | MATCH |
| `K4 = K0/(4 Omega^4)` | py L28/L57, out L2; wl L42, out L10 | notes L36 `Kbar_4 = Kbar_0/(4 Omega_Q^4)`; app. eq. 233 | MATCH |
| `Gamma5 = 9 K0/(32 Omega^5)` | py L30/L58, out L3; wl L46, out L11 | notes L50 `Gammabar_5 = 9 Kbar_0/(32 Omega_Q^5)`; app. eq. 281 | MATCH |
| `K0_target = 64 G Omega^5/(45 c^5)` | py L34/L59, out L4; wl L53, out L14 | notes L62 `Kbar_0^target = 64 G Omega_Q^5/(45 c^5)`; app. eq. 240 | MATCH |
| `K0_target geom = 54 G c_s^5/(5 a^5 c^5)` | py L36/L60, out L5; wl L55, out L15 | notes L68 `54 G c_s^5/(5 a^5 c^5)`; app. eq. 242 | MATCH |
| `Gamma5_target = 2 G/(5 c^5)` | py L41/L61, out L6; wl L65, out L18 | notes L77 `Gammabar_5^target = 2 G/(5 c^5)` | MATCH |
| `R0 = N_Q - 1` | py L47/L62, out L7; wl L74, out L20 | notes L104+L114 `R_0 = N_Q - 1` | MATCH |
| `R2 = N_Q - 1` | py L48/L63, out L8; wl L75, out L21 | notes L106+L114 | MATCH |
| `R4 = N_Q - 1` | py L49/L64, out L9; wl L76, out L22 | notes L108+L114 | MATCH |
| `R5 = N_Q - 1` | py L50/L65, out L10; wl L77, out L23 | notes L110+L114 | MATCH |
| `N_Q := Kbar_0/Kbar_0^target` (def) | implicit, py L46; wl L69 | notes L88; app. eq. 249; card body | MATCH |

INTERNAL (scaffolding, no finding expected): `Omega_geom = 3 c_s/(2 a)` (a carried upstream geometric pole — also stated in notes L66 and app. eq. 244, so it is in fact a documented MATCH, not orphaned); intermediate `K2_target/K4_target/K4Actual/k2Actual` etc. (built only to drive the R_i ratios); pass/fail flags and zero residuals.

reconciliation: complete; 11 deliverable values checked, 0 misaligned.

## Self-test notes

Checked: (1) variable-independence — no `diff`/`D` derivatives in either script, so the zero-derivative trap does not apply. (2) Series parity — `Yhat_Q^cons` is even in `w`, so the order-2/order-4 `SeriesCoefficient`s in the `.wl` are genuinely nonzero, making B1/B2 substantive. (3) Trivial-case — the R_i checks build target coefficients from unsubstituted `K0_target`, so substituting `N_Q=1` collapses each R_i to `0` while a wrong `Omega` power would leave a surviving factor; the checks can fail. (4) No missing-script findings, so no path-spec concern. (5) No fix prescribed (clean), so no paper round-trip risk. Conclusion: all traps clear; verdict clean.
