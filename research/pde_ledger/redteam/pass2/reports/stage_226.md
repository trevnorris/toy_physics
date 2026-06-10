---
unit_id: 226
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 226 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_226.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows/sections referencing 226: line 64, lines 681–745)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_mathematica_audit.txt`

## What the paper claims

The card (`\stagefield{Output}`) states: "Strict first-order packet: $\Xi_1$ is the prefactor slope, while $K_1$ and $H_{\rm even}$ are separate conservative gates; the survivor is a pure-transfer corridor rather than arbitrary grouped anisotropy." The derivation ledger states it "proves $\Xi_{\rm load}=\Xi_1=P_1/P_0$, shows that monomial/orbit rigidity kills the load defect but not the even gates, and narrows the viable same-charge corridor to a pure outgoing-transfer deformation with the conservative front end frozen." The notes (Section "Purpose", items 1–5) enumerate five deliverables: (1) the exact bridge $\Xi_{\rm load}=\Xi_1=P_1/P_0$; (2) the exact comparison of the Stage-225 compensation surface against the strict 5PN even gates (closed forms for $K_1$, $H_{\rm even}$ on that surface and on the one-pole branch); (3) the explicit mixed-sector-only strict-even-gate solve on the Stage-223 compatibility point (2×5 matrix, rank 2, nullity 3, null basis, $\Xi_1$ values, $\sigma_{\rm even}\approx2.67386816837173$); (4) the pure-transfer 3×5 intersection (rank 3, nullity 2, basis, $\Xi_1$/$N_{01}$ values, $\sigma_{\rm transfer}\approx2.31561904386057$); (5) the transported same-charge ceiling budgets on both corridors. The appendix (Theorem, line 740–743) and row 64 echo the same. The card is Mixed: ExactClosure + Reduced; the actual nonlinear PDE realization is explicitly Open.

## What the script claims to verify

The SymPy script (banner line 26) audits "strict 5PN even-gate package and pure-transfer subcorridor." Section 1 (lines 31–67) symbolically proves the bridge `Xi1 == Xi_load` and the three even-gate closed forms (compensation surface and one-pole branch). Section 2 (lines 69–188) builds the Stage-223 compatibility branch from physical primitives and pins D0/D2/D4/u2/u4/P0 plus the one-pole identity `u4 == 4 u2^2`. Section 3 (204–212) pins the D01-coefficients of $K_1$ and $H_{\rm even}$. Section 4 (217–265) builds the strict even-gate 2×5 matrix, checks rank 2 / nullity 3, the raw null basis, the $\Xi_1$ values, and `sigma_even`. Section 5 (270–318) does the pure-transfer 3×5 system, rank 3 / nullity 2, basis, $\Xi_1$/$N_{01}$, the explicit `D01=D21=D41=0` residual check on the basis, and `sigma_transfer`. Section 6 (323–360) divides the Stage-224 base budgets by each sigma. The Mathematica `.wl` mirrors the same seven blocks (M1–M7) but recomputes every quantity in its own kernel.

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| (1) bridge $\Xi_{\rm load}=\Xi_1=P_1/P_0$ | py L42–46 (`assert Xi1-Xi_load==0`); wl L83–87 (M1) | match |
| (2) even-gate closed forms ($K_1$,$H_{\rm even}$,one-pole) | py L48–61; wl L91–115 (M2) | match |
| compat branch D0/D2/D4/u2/u4/P0 + one-pole identity | py L160–179; wl L153–181 (M3) | match |
| D01-coefficients of $K_1$,$H_{\rm even}$ | py L204–208; wl L185–190 (M4) | match |
| (3) strict even-gate corridor (rank/nullity/basis/$\Xi_1$/$\sigma_{\rm even}$) | py L217–255; wl L223–238 (M5) | match (wl pins fewer literals — see Independence) |
| (4) pure-transfer subcorridor (rank/nullity/basis/$\Xi_1$/$N_{01}$/$\sigma_{\rm transfer}$/D=0) | py L270–310; wl L242–264 (M6) | match |
| (5) transported budgets (both corridors) | py L323–356; wl L266–310 (M7) | match |
| card Verification line "Mathematica audit: none yet" | `.wl` now exists | mismatch (stale doc pointer — F1) |

`paper_alignment: aligned` — every load-bearing physics deliverable is faithfully exercised by both engines; the single mismatch is a stale card-side verification pointer, not a math/value disagreement.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 46 | `simplify(Xi1 - Xi_load) == 0` | deliverable 1 | yes |
| A2 | sympy | 59–61 | 3× `simplify(... ) == 0` (gate closed forms) | deliverable 2 | yes |
| A3 | sympy | 171–178 | 8× `assert_close` branch data | compat branch | yes |
| A4 | sympy | 179 | `simplify(u4_s - 4 u2_s^2) == 0` | one-pole identity | yes |
| A5 | sympy | 207–208 | 2× `assert_close` D01-coeffs | gate coeffs | yes |
| A6 | sympy | 226–255 | matrix entries, rank==2, nullity==3, basis, $\Xi_1$, sigma_even | deliverable 3 | yes |
| A7 | sympy | 276–310 | rank==3, nullity==2, basis, $\Xi_1$/$N_{01}$, D=0 residuals, sigma_transfer | deliverable 4 | yes |
| A8 | sympy | 353–356 | 8× `assert_close` budgets | deliverable 5 | yes |
| M1 | math | 87 | `expectZero[xi1Series - xiLoad]` (Series route) | deliverable 1 | yes |
| M2 | math | 107–115 | 3× `expectZero` gate forms | deliverable 2 | yes |
| M3 | math | 174–181 | 7× `expectClose` + `expectZero[u4-4u2^2]` | compat branch | yes |
| M4 | math | 189–190 | 2× `expectClose` D01-coeffs | gate coeffs | yes |
| M5 | math | 232–238 | rank, nullity, sigma_even | deliverable 3 | yes |
| M6 | math | 243–264 | rank, nullity, residuals==0, sigma_transfer | deliverable 4 | yes |
| M7 | math | 292–310 | 8× `expectClose` budgets | deliverable 5 | yes |

All rows non-tautological and trace to a specific paper deliverable.

## Findings

### F1 — paper_misalignment

**Subtype:** paper_missing_script_claim (stale verification pointer)
**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_226.tex:11`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_mathematica_audit.wl` (whole file)

**What's wrong:**
The card states `\stagefield{Verification}{SymPy audit: ... Mathematica audit: none yet.}` (line 11), but a complete second-engine Mathematica audit now exists and passes all checks (`mathematica/output/...stage226...mathematica_audit.txt` ends "All Stage 226 Mathematica checks passed."). The notes line 11 (`\stagefield{Verification}`) is stale relative to the now-present `.wl`. This is a doc pointer the card carries about what is verified; the script side genuinely verifies more than the card admits.

**Why this matters:**
A reader trusting the card would believe stage 226 has no second-engine corroboration when in fact it has full SymPy+Mathematica agreement to ~20 digits. Pure doc drift; no math is wrong.

**Required change:**
Paper-side edit only — Codex must NOT auto-apply (paper card). Routed to user via `## Resolve before fix_loop`. Suggested fix: update line 11 to point at `mathematica/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_mathematica_audit.wl`.

**Verification:**
Card line 11 names the `.wl` path; no script change.

## Independent-derivation check (Mathematica)

INDEPENDENCE CALL: **independent**, with two noted cosmetic-method points.

The `.wl` recomputes every numeric quantity in its own kernel rather than echoing the `.py`'s numbers: M3 derives D0=24.2373099886222505…, M5/M6 build `evenMatrixExact`/`transferMatrixExact` from its own `coeffRow[…]` of `FullSimplify`-derived `d01Mixed`…`hEvenMixed`, and the nullspaces come from Mathematica's own `NullSpace`. Notably the `.wl` pins **fewer** hardcoded literals than the `.py` — it does NOT hardcode the 2×5/3×5 matrix entries or the null-basis vectors; it only checks the downstream invariants (rank, nullity, residuals==0, the sigmas, the budgets). That makes it strictly less of an echo than the `.py`.

The two deliberately-distinct operations confirmed:
- M1 uses `Series[p0A,{eps,0,1}]` + `Coefficient[…,eps]` (L83–84) vs the `.py`'s `sp.diff(P0A,eps).subs(eps,0)` (L42) — genuinely different expansion routes.
- The projector uses `Orthogonalize[N[basis,80]]` then `Transpose[qRows].qRows` (L48–52) vs the `.py`'s `QRdecomposition()` then `Q*Q.T` (L251–252, 306–307).

Honest caveat (not a finding): the projector swap is mathematically cosmetic — `Q Qᵀ` is the *unique* orthogonal projector onto the column space of the null basis regardless of whether the orthonormal basis came from Householder QR or Gram-Schmidt, so the two routes are guaranteed to produce the identical projector (and they do: `sigma_even` agrees to 2.67386816837172774819…, `sigma_transfer` to 2.31561904386055441584…). And the M3–M7 derivation *skeleton* (physical primitives Δ/S2/H/Q/P, B/Z bundles, dress-by-`Exp[eps·x]`, `D[…,eps]/.eps→0`, mixed-sector restriction, coefficient matrices) is the same algebra in both engines — but that is the one correct physics derivation, and each engine performs the computation independently in its own CAS, so this is corroboration, not policy-violating transliteration. One minor weakness: M2's `hOnePole` (L99–102) starts from the form `(u2var^2 - u4var)` (baking in `D4/D0 = u2²−u4`) and only verifies the `u4→4u2²` substitution, rather than deriving the one-pole form from `hEvenComp` the way the `.py` does (L55–57). The relation it bakes in is, however, independently confirmed numerically in M3 (`u4 == 4 u2²` and the actual `D4/D0` value), so the closed form is still genuinely corroborated. Verdict: independent; no `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree to full displayed precision:

| quantity | sympy out | mathematica out |
|---|---|---|
| D0 | 24.23730998862225 | 24.2373099886222505… |
| sigma_even | 2.67386816837172775 | 2.67386816837172774819… |
| sigma_transfer | 2.31561904386055442 | 2.31561904386055441584… |
| K1 coeff | 0.06219394707193093 | 0.06219394707193092742… |
| H_even coeff | −0.01160426115715844 | −0.01160426115715843674… |
| even budgets | [0.137602…,0.275862…,1.102857…,1.733464…] | identical to 18 digits |
| transfer budgets | [0.158890…,0.318540…,1.273480…,2.001648…] | identical to 18 digits |

No `engine_disagreement`.

## Verdict justification

`verdict: findings` solely because of the low-severity stale card-side verification pointer (F1, paper_misalignment routed to user). The math is sound: I attacked the bridge identity (M1/A1 — genuine symbolic, can fail), the even-gate closed forms (non-tautological), the one-pole identity `u4=4u2²` (a real residual==0, not assumed), the rank/nullity solves (computed by each kernel), the explicit `D01=D21=D41=0` residual check on the pure-transfer basis (M6 prints exact `{0,0,0}`), the projector-norm route (and confirmed the Orthogonalize-vs-QR swap yields the same unique projector, which is honest about the limited independence it buys), and the budget arithmetic. None broke. Both engines independently recompute and agree to ~20 digits. I read the card, notes, and appendix; the script's verified claims match the paper's deliverables. No `stale_output` (both `.txt` newer than both scripts). No stage-number-label drift in either script (both say 226 throughout — the deferred-class label drift is absent here). The `is_status_only_candidate: False` checkpoint bar is met: both engines present and substantive.

## Self-test notes

Variable-independence: every `D[…, eps]/.eps→0` / `sp.diff(…,eps)` acts on `…d/Exp[eps·x]`-dressed expressions that genuinely depend on `eps`, so no derivative is identically zero (the matrices are nonzero, as their printed entries show). Symmetry/parity: n/a — no unbounded-domain integrals; all checks are finite linear-algebra/series. Trivial-case: the pure-transfer residual check (M6, py L300–303) substitutes the actual null-basis vectors and the printed residuals are exactly `{0,0,0}`, confirming `assert_zero` is not vacuous; sigma checks return concrete nonzero literals. Paper round-trip: F1 is paper-side and routes to user; no script edit prescribed, so no new paper_misalignment introduced.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 23 deliverable values checked, 0 misaligned.

All emitted deliverable values reconcile to the notes `.md` (the natural carrier; the `.tex` card is intentionally terse and the appendix carries the symbolic forms but not the branch numerics — both legitimate per the augmentation guards). Outputs are fresh, so reconciliation is anchored in both source and committed `.txt`.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| $\Xi_{\rm load}=\Xi_1=N_{01}/N_0-D_{01}/D_0$ (bridge form) | py L44/out L4; wl L85/out L9 | notes L150–155 (boxed) | MATCH |
| $K_1$(comp) $=(1/9-u_2)D_{01}$ | py L59/out L5; wl L107/out L15 | notes L179 (boxed) | MATCH |
| $H_{\rm even}$(comp) $=(D_4/D_0+2u_2/3-1/27)D_{01}$ | py L60/out L6; wl L109/out L16 | notes L182 (boxed) | MATCH |
| $H_{\rm even}$(one-pole) $=(-3u_2^2+2u_2/3-1/27)D_{01}$ | py L61/out L7; wl L113/out L17 | notes L194 (boxed) | MATCH |
| D0 = 24.2373099886223 | py L172/out L11; wl L174/out L25 | notes L82 | MATCH |
| D2 = −1.18562046858190 | py L173/out L12; wl L175/out L26 | notes L83 | MATCH |
| D4 = −0.173991572849491 | py L174/out L13; wl L176/out L27 | notes L84 | MATCH |
| u2 = 0.0489171640391802 | py L175/out L14; wl L177/out L28 | notes L86 | MATCH |
| u4 = 0.00957155575054425 | py L176/out L15; wl L178/out L29 | notes L87 | MATCH |
| P0 = 0.00206979231806289 | py L178/out L16; wl L180/out L30 | notes L98 (0.002069792318062885) | MATCH |
| K1 coeff = 0.0621939470719309 | py L207/out L19; wl L189/out L43 | notes L202 | MATCH |
| H_even coeff = −0.0116042611571584 | py L208/out L20; wl L190/out L44 | notes L203 | MATCH |
| $\sigma_{\rm even}$ = 2.67386816837173 | py L255/out L36; wl L238/out L54 | notes L276 (boxed) | MATCH |
| $\Xi_1(w_1,w_2,w_3)$ = 1.33691841376792 / −13.9944400566810 / −5.02163500066813 | py L246/out L35; (wl computes, no literal) | notes L260–266 | MATCH |
| $\sigma_{\rm transfer}$ = 2.31561904386057 | py L310/out L44; wl L264/out L71 | notes L357 (boxed) | MATCH |
| $\Xi_1(t_1,t_2)$ = 11.0106276743889 / −5.66658382170817 | py L292/out L42; (wl computes) | notes L331–332 | MATCH |
| $N_{01}(t_1,t_2)$ = 0.552361328292489 / −0.284270966124842 | py L293/out L43; (wl computes) | notes L337–338 | MATCH |
| even budgets = 0.137602269567650 / 0.275862165676603 / 1.10285760977778 / 1.73346419189450 | py L341–346/out L47; wl L274–279/out L77 | notes L383,388,394,395 | MATCH |
| transfer budgets = 0.158890698998242 / 0.318540765855427 / 1.27348056877049 / 2.00164821411704 | py L347–352/out L48; wl L280–285/out L78 | notes L408,413,419,420 | MATCH |
| 2×5 even-gate matrix entries | py L222–225/out L24–30; (wl derives, no literal pin) | notes L228–231 | MATCH |
| null basis w1/w2/w3 | py L234–237/out L32–34 | notes L238–256 | MATCH |
| pure-transfer basis t1/t2 | py L280–283/out L40–41 | notes L311–322 | MATCH |
| base (Stage-224 import) budgets 0.367930328492646 / 0.737619063660757 / 2.94889585703134 / 4.63505472371892 | py L323–326; wl L268–273 | notes L126,130,135,137 | MATCH (carried-forward import) |

INTERNAL items (scaffolding, no prose expected, no finding): `K_compat = 24.4737548792910` (the stiffness K that *pins* the compatibility branch — a setup intermediate, asserted py L171 / wl L174 but not a stated deliverable; not carried in notes/.tex, correctly INTERNAL); pass/fail flags ("OK …", "All checks passed"); residual `{0,0,0}` near-zero verification values; tolerances; the `Series`/derivative intermediate `p0Series`; `Xi_coeff`/`N01_coeff` coefficient vectors.
