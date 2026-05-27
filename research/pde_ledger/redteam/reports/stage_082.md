---
unit_id: 082
batch: III.4
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage082_master_quadrupole_residual.md
  paper_appendix: present
---

# Audit unit 082 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_082.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage082_master_quadrupole_residual.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (read only the rows referencing stage 082; row 142 is the one-line summary, line 282 includes the stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.txt`

## What the paper claims

Paper card stage_082 boxes three formulas: the demand `zeta_req(Pi_tr, C_mix, eps_blk) = (Pi_tr - C_mix) / [ C_mix - eps_blk(2*C_mix - Pi_tr) ]` (eq. 082-zeta-req), the physical supply `zeta_phys(Pe, eta; kappa) = Omega_Pe^2 * (kappa + pi^2/4) / (kappa + y(eta)^2)` with `y tan y = eta` (eq. 082-zeta-phys), and the master residual `R_quad = zeta_req(Pi_tr,C_mix,eps_blk) - zeta_phys(Pe_*(Xi,eta,kappa), eta; kappa)` (eq. 082-Rquad, the `\stagefield{Output}`). The card also fixes the sign convention (`R_quad < 0` = excess support, `= 0` = saturation, `> 0` = failure) and specializes to Family-1 via `(eta, kappa) = (37, 12321/5)` and `Xi_F1 = W_wall = 1369*Upsilon_w = 136900*Theta_w` (eq. 082-XiF1). The notes section 4 also documents the inverse map `Pi_tr = C_mix * Q(zeta; eps_blk)` with `Q = (1 + (1-2*eps_blk)*zeta)/(1 - eps_blk*zeta)` and the product thresholds `Pi_suff = C_mix Q(zeta_-)`, `Pi_fail = C_mix Q(zeta_+)`, derived as consequences of `zeta_req`'s monotonicity in `Pi_tr`.

## What the script claims to verify

The sympy script docstring lists four claims: (1) exact inverse map between `zeta_req` and `Pi_tr = C_mix Q(zeta)`; (2) product thresholds `Pi_suff = C_mix Q(zeta_-)`, `Pi_fail = C_mix Q(zeta_+)`; (3) Family-1 strength identity `Xi_F1 = 1369 Upsilon_w = 136900 Theta_w`; (4) master residual definition. The assertions encode `zeta_req` in functional form, derive `Q` (independently in the .wl via `Solve`), then check `zeta_req(C_mix Q(z)) - z == 0`, the two threshold inversions, and `R_quad(Pi_suff, zeta_-) == 0`, `R_quad(Pi_fail, zeta_+) == 0`. Two derivative sanity checks `dR/dzeta_phys + 1 == 0` and `dR/dPi_tr - dzeta_req/dPi_tr == 0` are added. Family-1 constants are now display-only (v1 already converted the arithmetic `expect_zero` lines to `print` per directive F1). Critically, `zeta_phys` is treated as an opaque symbol throughout — the script never encodes the paper-side formula `Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2)`, and never instantiates the Family-1 numerical pair `(eta, kappa) = (37, 12321/5)`.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| eq. 082-zeta-req: `zeta_req = (Pi_tr - C_mix)/[C_mix - eps_blk(2 C_mix - Pi_tr)]` | sympy line 38; mathematica line 38 (re-derived via `Solve`) | match |
| eq. 082-zeta-phys: `zeta_phys = Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2)` | none — `zeta_phys` is left as a bare `sp.Symbol`/`Reals` element (sympy line 63; mathematica line 31) | missing |
| eq. 082-Rquad: `R_quad := zeta_req - zeta_phys(Pe_*(Xi,eta,kappa), eta; kappa)` | sympy line 64; mathematica line 59 — but only the algebraic skeleton `zeta_req - zeta_phys` is encoded; the `Pe_*` operator-selected bias dependence on `(Xi, eta, kappa)` is never exhibited | partial |
| eq. 082-sign convention (R_quad <0/=0/>0) | none — sign of R_quad is never checked at any concrete branch point | missing |
| eq. 082-XiF1: `Xi_F1 = W_wall = 1369*Upsilon_w = 136900*Theta_w` and `(eta,kappa) = (37, 12321/5)` | sympy lines 96-97; mathematica lines 88-90: arithmetic on hardcoded `Lambda_ell = 37` and `Upsilon_w = 100*Theta_w`. Now display-only after v1. `(eta, kappa) = (37, 12321/5)` numerical pair is not encoded anywhere. | partial |
| notes section 4: inverse map `Pi_tr = C_mix Q(zeta_req)` and threshold theorems | sympy lines 46-49, 51-52, 59-60; mathematica lines 45-57 | match |
| notes section 6: `R_quad` is the load-bearing residual everything else feeds | partial — R_quad is defined, but with `zeta_phys` opaque | partial |

Dominant pattern: the script verifies the inverse-map / threshold mechanics from notes section 4 faithfully, but two of the paper card's three boxed equations (eq. 082-zeta-phys and the Family-1 numerical specialization) are not exercised at all, and the sign-convention deliverable is absent. Set `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 46-49 | `zeta_req(C_mix*Q(zeta)) - zeta == 0` | notes section 4 inverse map | yes (non-tautological — exercises rational-function algebra) |
| A2 | sympy | 59 | `zeta_req(Pi_suff) - zeta_- == 0` | notes section 4 threshold | partial (A1 specialized at zeta -> zeta_-) |
| A3 | sympy | 60 | `zeta_req(Pi_fail) - zeta_+ == 0` | notes section 4 threshold | partial (A1 specialized at zeta -> zeta_+) |
| A4 | sympy | 68-71 | `R_quad(Pi_suff, zeta_-) == 0` | eq. 082-Rquad (definition only) | no (definition `R_quad = zeta_req - zeta_phys` + A2) |
| A5 | sympy | 72-75 | `R_quad(Pi_fail, zeta_+) == 0` | eq. 082-Rquad (definition only) | no (definition + A3) |
| A6 | sympy | 80 | `dR_quad/dzeta_phys + 1 == 0` | eq. 082-Rquad (definition only) | no (`d(a-b)/db = -1` by construction) |
| A7 | sympy | 83-87 | `dR_quad/dPi_tr - dzeta_req/dPi_tr == 0` | eq. 082-Rquad (definition only) | no (`d(zeta_req - 0)/dPi = dzeta_req/dPi` by construction since R_quad - zeta_req = -zeta_phys and zeta_phys is Pi_tr-independent) |
| A8 | sympy | 108-109 | display-only prints for Xi_F1 arithmetic | eq. 082-XiF1 | n/a (display-only after v1 fix) |
| A9 | mathematica | 38 | `Solve[Pi_tr = cMix*qMap, zetaSym]` recovers zeta_req | notes section 4 inverse map | yes (independent re-derivation via Solve) |
| A10 | mathematica | 45-48 | `zetaReq(cMix*qMap) - zeta == 0` | notes section 4 inverse map | partial (engine independence: SymPy uses pre-baked closed form, Mathematica solves) |
| A11 | mathematica | 56-57 | `zetaReq(piSuff/piFail) - zetaMinus/zetaPlus == 0` | notes section 4 threshold | partial |
| A12 | mathematica | 62-69 | `rQuad` evaluations at branch points | eq. 082-Rquad | no (mirrors A4/A5) |
| A13 | mathematica | 73-81 | derivative checks | eq. 082-Rquad | no (mirrors A6/A7) |
| A14 | mathematica | 99-102 | display-only prints for Xi_F1 arithmetic | eq. 082-XiF1 | n/a |

## Findings

### F1 — paper_misalignment

**Subtype:** script_missing_paper_claim

**Severity:** medium

**Files:**
- paper: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_082.tex:21-25` (eq. 082-zeta-phys)
- script: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py:63`
- script: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl:31`

**What's wrong:**

The paper card explicitly writes `zeta_phys` as a specific closed-form expression:

> `zeta_phys(Pe, eta; kappa) = Omega_Pe^2 * (kappa + pi^2/4) / (kappa + y(eta)^2)`  with  `y tan y = eta`.   (paper stage_082.tex, eq. 082-zeta-phys, lines 21-25)

The sympy script declares only

> `zeta_phys = sp.symbols("zeta_phys", real=True)`   (sympy line 63)

and the mathematica script declares only the type assumption

> `Element[{PiTr, epsBlk, zeta, zetaMinus, zetaPlus, zetaPhys}, Reals]`   (mathematica line 31).

`zeta_phys` is then carried through `R_quad = zeta_req - zeta_phys` as an opaque scalar. The closed-form expression from eq. 082-zeta-phys, the Robin-root condition `y tan y = eta`, and the transport-overlap factor `Omega_Pe` (notes section 1.4) are never instantiated. As a result, every assertion involving `zeta_phys` (A4-A7, A12-A13) reduces to a sanity check on linear algebra over a free symbol; the paper's specific functional form for the physical support ratio is not exercised by either engine. This is the paper card's second of three boxed equations, and the script does not touch it.

**Why this matters:**

The whole point of the master quadrupole residual is that `zeta_phys` carries a non-trivial functional dependence on `(Pe, eta, kappa)` through the Robin eigenvalue `y(eta)`. Verifying `R_quad := zeta_req - zeta_phys` at the symbol level confirms only the subtraction structure; it cannot catch sign errors, branch-cut errors, or normalization slips in the `Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2)` formula. The paper's Output equation (eq. 082-Rquad) cites both halves of the residual, but the script verifies only one half (zeta_req).

**Required change:**

`## Resolve before fix_loop` — see directive.

**Verification:**

Either the script grows checks that instantiate `zeta_phys = Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2)` and confirms at least one nontrivial property (e.g., that R_quad evaluates to a specific symbolic form at a representative Robin root), OR the paper card is updated to clarify that the script only verifies the algebraic skeleton of the residual and defers the zeta_phys functional verification to its upstream stage (stages 39-40 per notes section 1.4). The orchestrator routes this to the user.

### F2 — paper_misalignment

**Subtype:** script_missing_paper_claim

**Severity:** low

**Files:**
- paper: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_082.tex:43-47` (eq. 082-XiF1 and the `(eta,kappa) = (37, 12321/5)` line)
- script: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py:95-97`
- script: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl:88-90`

**What's wrong:**

The paper's Family-1 specialization in stage 082 has two parts:

> "The Family--1 specialization is obtained by setting `(eta, kappa) = (37, 12321/5)` and `Xi_F1 = W_wall = 1369*Upsilon_w = 136900*Theta_w`."   (paper stage_082.tex, lines 43-47)

The script only exercises the *Xi_F1* line of the specialization (and only as display-only arithmetic on hardcoded `Lambda_ell = 37` and `Upsilon_w = 100*Theta_w`). The numerical pair `(eta, kappa) = (37, 12321/5)` is never set as a variable or substituted into any expression in either script. The notes section 5 makes the same numerical specialization, but no script line references `12321/5` or pairs `eta = 37` against `kappa = 12321/5` in any check.

Additionally, the script's choice `Upsilon_w = 100*Theta_w` (sympy line 97, mathematica line 90) implicitly relies on the relation `136900/1369 = 100` from the paper card. But paper stage 075 line 7 states `Upsilon_w = 117*Theta_w`, not `100*Theta_w`. This is a paper-internal disagreement between stage 082 and stage 075 about the Upsilon_w/Theta_w relation — the script side is consistent with stage 082's own line but inconsistent with stage 075. Flag it for user attention; the script cannot resolve a paper-internal conflict.

**Why this matters:**

If the paper's intent in stage 082 is that the Family-1 specialization carries the full `(eta=37, kappa=12321/5, Xi_F1 = 136900*Theta_w)` triple, then the script should at least display `R_quad` with `(eta, kappa)` instantiated. The current state hides whether the numerical pair is even consistent with upstream stages 073-074 (which derive kappa = 9/5 * Lambda_ell^2 = 12321/5 from Lambda_ell = 37). The Upsilon_w 100-vs-117 conflict between stage 082 and stage 075 is a real paper-side inconsistency that may need the user to decide which one is canonical.

**Required change:**

`## Resolve before fix_loop` — see directive.

**Verification:**

User resolves whether (a) the script should add a numerical block instantiating `(eta, kappa) = (37, 12321/5)`, and (b) which paper stage's Upsilon_w/Theta_w relation is canonical (stage 082's implicit 100, or stage 075's explicit 117).

### F3 — insufficient_verification

**Severity:** low

**Files:**
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py:79-87`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl:73-81`

**What's wrong:**

The two derivative checks added by the v1 directive are mathematically trivial:

- `dR_quad/dzeta_phys + 1 == 0` is identically true by `R_quad = zeta_req - zeta_phys`: differentiating a linear difference w.r.t. its subtracted symbol always yields `-1`. The engine cannot fail this check for any well-formed `zeta_req` not containing `zeta_phys`.
- `dR_quad/dPi_tr - dzeta_req/dPi_tr == 0` is identically true because `zeta_phys` was declared as a free symbol with no `Pi_tr` dependence; so `d(zeta_req - zeta_phys)/dPi_tr = dzeta_req/dPi_tr` by construction.

These were added by v1 as "directional content of R_quad" but they don't exercise the directional content of R_quad against zeta_phys's actual functional form — they only confirm the engine's symbolic differentiation rules. Severity is low because they're cheap printouts and don't cost much, but they don't move the verification needle.

**Why this matters:**

The intent (per the inline comment "verify the partial derivatives that underwrite the 'guaranteed success / guaranteed failure' theorems") is to demonstrate that R_quad is strictly monotonic in zeta_phys and in Pi_tr. The actual checks verify only that the symbolic difference operator works as expected, not that the physical residual has the claimed monotonicity at any concrete point. The sign content of R_quad (paper eq. 082-sign) is not exercised anywhere.

**Required change:**

Either delete the two tautological derivative checks, or replace them with assertions that exercise the *sign* of `dzeta_req/dPi_tr` (the real physical content — that zeta_req is strictly increasing in Pi_tr on the allowed branch). For example, assert `sign(d zeta_req/d Pi_tr) > 0` under the branch positivity conditions, which is the actual content notes section 4 relies on ("Because zeta_req is strictly increasing in Pi_tr, the bounded residual can be translated back into exact product thresholds").

**Verification:**

After fix, the script should contain an assertion that fails if `dzeta_req/dPi_tr` were not strictly positive on the branch domain (e.g., evaluating it at a positive-branch sample point and checking the result simplifies to a positive ratio of squares). The current trivial-by-construction checks should either be removed or have a comment marking them as engine-correctness sanity rather than physics verification.

## Independent-derivation check (Mathematica)

The Mathematica script does one substantive thing differently from the SymPy script: it derives `zeta_req` by calling `Solve[PiTr == cMix * qMap(zetaSym), zetaSym]` (line 38), rather than starting from the closed-form `zeta_req = (Pi_tr - C_mix)/[C_mix - eps_blk(2 C_mix - Pi_tr)]` directly. This is a genuine independent re-derivation: SymPy writes the closed form and verifies the inverse; Mathematica writes the inverse map `Q` and solves for `zeta_req`, then verifies they agree. Quoting the two key sections:

- SymPy line 38: `zeta_req = sp.simplify((Pi_tr - Cmix) / (Cmix - eps_blk * (2 * Cmix - Pi_tr)))` — a direct write-down.
- Mathematica lines 38-40: `Solve[PiTr == cMix*((1 + (1 - 2*epsBlk)*zetaSym)/(1 - epsBlk*zetaSym)), zetaSym]; zetaReq = FullSimplify[(zetaSym /. First[zetaReqSolList]) /. ConditionalExpression[x_, _] :> x, ...]` — solves the inverse map and strips the conditional.

That's a non-trivial differentiation between the two engines (one writes the form, the other derives it from the inverse). No `mathematica_transliteration` finding. The rest of the script is parallel structure, but parallel structure of an algebraic identity check is hard to avoid; the load-bearing difference is at line 38.

## Engine cross-check

Both engines reach identical final forms:

- `zeta_req` form: identical up to sign of the denominator, `(-cMix + PiTr)/(cMix - 2*cMix*epsBlk + epsBlk*PiTr)` in both transcripts.
- `R_quad = (-cMix + PiTr)/(cMix - 2*cMix*epsBlk + epsBlk*PiTr) - zetaPhys` in both transcripts.
- All `expectZero` / `expect_zero` checks return `0` in both transcripts.
- Family-1 prints: both show `Xi_F1 from Upsilon_w = 1369*Upsilon_w`, `Xi_F1 from Theta_w = 136900*Theta_w`.

Engines agree. No `engine_disagreement` finding.

## Verdict justification

The script's notes-section-4 deliverables (inverse map and threshold theorems) hold up under attack — A1 is a genuine non-trivial rational-function identity, A2-A3 specialize it, and Mathematica's `Solve`-based re-derivation strengthens the independence claim. However, two of the paper card's three boxed equations (eq. 082-zeta-phys and the `(eta, kappa)` numerical specialization in eq. 082-XiF1) are not exercised by either engine, and the sign-convention deliverable is absent. The v1 derivative checks added at lines 79-87 / 73-81 are tautological by construction and don't compensate. Two `paper_misalignment` findings require user resolution before any further script edits; the third `insufficient_verification` finding can proceed independently. No `stop_cold` flag: F1 and F2 are paper-side-vs-script-side scope decisions, not mathematical inconsistencies, and the script's core mechanical content (inverse map, thresholds, residual definition) holds up.

## Self-test notes

Traps checked: (1) verified `dR/dzeta_phys = -1` is structurally identical to `d(a-b)/db = -1` since `zeta_req` contains no `zeta_phys` and vice versa — confirms F3's "trivial by construction" call; (2) checked `dR/dPi_tr - dzeta_req/dPi_tr` simplifies because `zeta_phys` is a free symbol with no `Pi_tr` dependence, again confirming F3; (3) verified that `136900/1369 = 100` matches stage 082 but `Upsilon_w = 117*Theta_w` in stage 075 line 7 is an independent paper-side fact, supporting F2's flag for user resolution; (4) confirmed F1 and F2 are paper_misalignment (script does not exercise paper claim) rather than insufficient_verification (which would require the script claim itself to fall short of its own docstring).
