---
unit_id: 089
batch: III.5
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage089_family1_minimal_isotropic_verdict.md]
  paper_appendix: present
---

# Audit unit 089 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_089.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage089_family1_minimal_isotropic_verdict.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 156)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_mathematica_audit.txt`

## What the paper claims

Stage 089 is a checkpoint stage whose `\stagefield{Output}` is the zero-bias Family-1 success theorem, with the boxed equations `\zeta_{\rm req}^{\rm min} = 1/3 < A_{\rm F1}` (eq. app-stage089-success) and `\mathrm{Pe}_{\rm req} = 0` (eq. app-stage089-Pe-zero). The card's `\stagefield{Inputs}` are `\zeta_{\rm req} = 1/3` and `A_{\rm F1} \simeq 1.00005192880220`. The derivation argues that at zero bias `\Omega_{\rm Pe=0} = 1`, so `\zeta_{\rm F1}(0) = A_{\rm F1}`, and the boxed inequality then gives `Pe_req = 0`. The notes add three further deliverables: (i) the loading-ratio margins `\Delta_{\rm suff}, \Delta_{\rm fail}, \Delta_{\rm max}` against the Stage-69/63 Family-1 window; (ii) the symmetric-lowest-twin regime check `0 < \zeta_{\rm req}^{\rm min} = 1/3 < 1`; (iii) the transport-map identity `\zeta_{\rm F1}(\mathrm{Pe}) = A_{\rm F1}\,\Omega_{\rm Pe}^{2}` from Stage 62. The appendix row (line 156) summarizes: "Minimal isotropic branch succeeds at zero transport bias."

## What the script claims to verify

The SymPy script defines `rho_min = 4/3`, `zeta_min = 1/3`, computes `A_F1` from `y tan(y) = 37`, `kappa_F1 = 12321/5`, builds `Omega(Pe)` and `zeta_F1(Pe) = A_F1 Omega^2`, evaluates `zeta_max = lim_{Pe->oo} zeta_F1`, then substitutes hardcoded `Pe_suff_chi = 96.5285...` and `Pe_fail_chi = 11220.544...` into `zeta_F1`, applies `Q(zeta; 0) = (1+(1-0)zeta)/(1-0·zeta) = 1+zeta` to get `rho_suff, rho_fail, rho_max`. It asserts: the algebraic anchor `zeta_max = A_F1 pi^2/4`, the Q-at-eps=0 identity for each rho, the ordering `rho_min < rho_suff < rho_fail < rho_max`, `zeta_min < 1`, `zeta_min < A_F1`, and `zeta_min < zeta_max`. The Mathematica script mirrors this structure verbatim. Neither script symbolically constructs `Pe_req = 0` nor verifies `Omega(Pe=0) = 1` (the 0/0 limit that connects `zeta_min < A_F1` to the boxed `Pe_req = 0`).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `A_F1 ≈ 1.00005192880220` (Input) | derived in script from `kappa_F1=12321/5`, `y tan(y)=37`; numeric match | match |
| `zeta_req^min = 1/3 < A_F1` (boxed eq.) | `if not (zeta_min < A_F1): raise` (sympy:91); `expectTrue[..., zetaMin < aF1]` (wl:83) | match |
| `Pe_req = 0` (boxed, Output) | not symbolically constructed; relies on unverified `Omega(Pe=0) = 1` link | partial |
| Notes: `Delta_suff, Delta_fail, Delta_max` against Stage-69 window | printed numerically (sympy:81-83; wl:75-77); ordering asserted | match |
| Notes: symmetric-lowest-twin regime `0 < 1/3 < 1` | asserted `zeta_min < 1` (sympy:89; wl:82) | match |
| Notes/Stage-62: `zeta_F1(Pe) = A_F1 Omega^2` form | built in code (sympy:45; wl:43) but not asserted | partial |

paper_alignment: `partial` — the headline boxed equation `Pe_req = 0` is not directly produced; its supporting link `Omega(Pe=0) = 1` is never checked.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 60 | `expect_zero(zeta_max - A_F1 pi^2/4)` | Stage-62 carry-forward identity (notes §3) | partial — algebraic anchor only |
| A2 | sympy | 61 | `expect_zero(Q - (1+zeta))` at eps=0 | none (algebraic substitution of eps=0) | no — tautological |
| A3 | sympy | 62-64 | `expect_zero(rho_X - (1+zeta_X))` | none | no — tautological by construction (Q.subs(eps=0)=1+zeta then subs zeta=zeta_X gives 1+zeta_X) |
| A4 | sympy | 87 | `rho_min < rho_suff < rho_fail < rho_max` | notes §1 loading-ratio margin (Delta_suff/fail/max ordering) | yes |
| A5 | sympy | 89 | `zeta_min < 1` | notes §2 symmetric-lowest-twin regime | yes |
| A6 | sympy | 91 | `zeta_min < A_F1` | boxed eq. app-stage089-success | yes |
| A7 | sympy | 93 | `zeta_min < zeta_max` | notes §1 support ceiling margin | yes |
| A8 | mathematica | 56 | `expectApprox[zetaMax, aF1 Pi^2/4, 1e-30]` | mirror of A1 | partial |
| A9 | mathematica | 57-59 | `expectApprox[rhoX, 1+zetaX, 1e-30]` | mirror of A2-A3 | no — tautological |
| A10 | mathematica | 81 | `rhoMin < rhoSuff < rhoFail < rhoMax` | mirror of A4 | yes |
| A11 | mathematica | 82 | `zetaMin < 1` | mirror of A5 | yes |
| A12 | mathematica | 83 | `zetaMin < aF1` | mirror of A6 | yes |
| A13 | mathematica | 84 | `zetaMin < zetaMax` | mirror of A7 | yes |

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Subtype:** script_missing_paper_claim
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_089.tex:27-29` quote: "\boxed{\mathrm{Pe}_{\rm req}=0}"
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py:91-92` quote: `if not (zeta_min < A_F1): raise AssertionError("Minimal isotropic branch no longer succeeds at zero transport bias.")`

**What's wrong:**
The paper's `\stagefield{Output}` is the zero-bias Family-1 success theorem `Pe_req = 0` (eq. app-stage089-Pe-zero), and the derivation hinges on `\Omega_{\rm Pe=0} = 1` so that `\zeta_{\rm F1}(0) = A_{\rm F1}`. The scripts never construct `Pe_req` as an object, never compute (or take the symbolic limit of) `Omega(Pe -> 0)`, and never assert `zeta_F1(0) = A_F1`. The Omega expression on `sympy:41-44` and `wl:42` is `pi Pe (2 Pe e^Pe + pi) / ((4 Pe^2 + pi^2)(e^Pe - 1))`, which is `0/0` at Pe = 0 — its limit is what closes the inference `zeta_min < A_F1 ⇒ Pe_req = 0`. The assertion `zeta_min < A_F1` is a necessary precondition for the boxed claim but does not itself establish `Pe_req = 0`.

**Why this matters:**
This is a checkpoint stage with `is_status_only_candidate: False`. The paper's bottom-line boxed equation is the unverified link. A reader could replace `Omega(Pe=0) = 1` with `Omega(Pe=0) = 2` and the script's PASS would not change.

**Required change:**
Add a one-line symbolic check `Limit[Omega(Pe), Pe -> 0] == 1` (Mathematica) and the sympy equivalent `sp.limit(Omega, Pe, 0) == 1`. Then assert `sp.simplify(zeta_F1.subs(Pe, 0_limit) - A_F1) == 0` (or take the limit explicitly), closing the chain to `Pe_req = 0`. This is the missing link the paper card relies on.

**Verification:**
After fix, sympy output should print a new line "Stage-62 Omega(Pe->0) = 1" and "zeta_F1(Pe->0) = A_F1" with residuals 0; Mathematica output should print the same checks; both should PASS.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_mathematica_audit.wl:34-54`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py:31-58`

**What's wrong:**
The `.wl` script is a line-by-line port of the `.py` script with identical variable choreography, identical constants, and identical intermediate names:

- `kappa_F1 = sp.Rational(12321, 5)` (sympy:37) ↔ `kappaF1 = 12321/5` (wl:38)
- `eta_F1 = sp.Integer(37)` (sympy:38) ↔ `etaF1 = 37` (wl:39)
- `Pe_suff_chi = sp.Float("96.5285247264386")` (sympy:49) ↔ `peSuffChi = SetPrecision[96.5285247264386, 40]` (wl:47)
- `Pe_fail_chi = sp.Float("11220.5441626259")` (sympy:50) ↔ `peFailChi = SetPrecision[11220.5441626259, 40]` (wl:48)
- `Q = (1 + (1 - 2 eps_blk) zeta)/(1 - eps_blk zeta)` (sympy:53) ↔ `q[zeta_, eps_] := (1 + (1 - 2 eps) zeta)/(1 - eps zeta)` (wl:49)
- `rho_suff = Q.subs(zeta, zeta_suff)` (sympy:56) ↔ `rhoSuff = N[q[zetaSuff, 0], 50]` (wl:52)

Even the saved `Pe_suff_chi`/`Pe_fail_chi` literals carry over identically, including the `0.0000004` digits — these are external Stage-63/69 numerics, not independently derived. Both engines also share the identical "anchor" sanity checks (A1↔A8, A2↔A9, A3↔A9). The second-engine policy requires independent derivation from physical premises; this is parallel transcription.

**Why this matters:**
Checkpoint-level audits require two independent engines. As written, an algebra error in `Omega(Pe)` or in the substitution of `Pe_suff_chi` would propagate identically to both engines, so engine agreement says nothing about correctness.

**Required change:**
Have the Mathematica script construct A_F1 and the Stage-69 thresholds from a different route — e.g., compute `Pe_suff_chi` and `Pe_fail_chi` by solving `zetaF1[Pe] == zeta_suff_target` directly with `Solve`/`FindRoot` rather than carrying the numeric Pe values forward as literals, OR introduce an alternative form of the Q identity (e.g., verify by expanding `Q` symbolically and matching coefficients), so the second engine does not echo the same numeric carry-forwards.

**Verification:**
The .wl script no longer contains the literals `96.5285247264386` and `11220.5441626259`; instead it derives them from Stage-69 inputs. The output residuals should still match the SymPy values within 1e-25.

### F3 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py:53,56-58,61-64`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_mathematica_audit.wl:49,52-54,57-59`

**What's wrong:**
At `eps_blk = 0` the symbolic identity `Q(zeta; 0) = (1 + (1 - 0)·zeta)/(1 - 0·zeta) = 1 + zeta` holds algebraically by construction. The script then defines `rho_suff = Q.subs(zeta, zeta_suff)` (sympy:56) and asserts `expect_zero("rho_suff anchor", rho_suff - (1 + zeta_suff))` (sympy:62). Because `Q.subs(zeta, zeta_suff)` evaluates to `1 + zeta_suff` literally — that is the same symbolic substitution sympy performs — the residual is zero by construction. The same pattern repeats for rho_fail and rho_max, and is mirrored in wl:57-59. The sympy output transcript confirms `Stage-69 Q(zeta;0) = 1 + zeta = 0`, `rho_suff anchor = 0`, `rho_fail anchor = 0`, `rho_max anchor = 0` (all four trivially zero).

**Why this matters:**
These three anchors look like Stage-69 carry-forward verification but verify nothing — they test the substitution `Q.subs(eps_blk=0) = 1 + zeta`, which is algebraic at the level of the symbolic definition. They give a false sense of "consistency with Stage 69."

**Required change:**
Replace the rho_X anchors with a non-trivial Stage-69 cross-check, e.g., assert that the numeric `rho_suff` matches the notes' literal `≈ 3.46622291347846` to a fixed tolerance — that does encode an external Stage-69 datum. (The Q-at-eps=0 identity check at sympy:61 / wl:57 may be retained as documentation but is also tautological.)

**Verification:**
The tautological anchors are replaced by an explicit numeric cross-check against the notes-quoted values, e.g., `expect_close(rho_suff, sp.Float("3.46622291347846"), tol=1e-12)`.

### F4 — hardcoded_result

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py:49-50`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_mathematica_audit.wl:47-48`

**What's wrong:**
`Pe_suff_chi = sp.Float("96.5285247264386")` and `Pe_fail_chi = sp.Float("11220.5441626259")` are hardcoded numeric literals. The script's comment says "Stage-63/69 thresholds evaluated at lambda_mu = 1" but there is no in-script derivation and no upstream-script citation pointing to the file that produced these values. The whole "loading-ratio ordering" assertion (`rho_min < rho_suff < rho_fail < rho_max`) depends on these literals — change `Pe_suff_chi` to `Pe_fail_chi` and back and the ordering still passes; change it to `0.0001` and the ordering fails. The constants are load-bearing but unprovenance-anchored.

**Why this matters:**
Checkpoint stages are required to be substantive. The numeric literals are inputs the user must trust as having been produced upstream — but the stage card's `\stagefield{Inputs}` lists only `\zeta_{\rm req} = 1/3` and `A_{\rm F1}`, not these Pe values. The notes §1 quote `rho_suff^(chi) ≈ 3.46622291347846` (a derived consequence) but do not quote `Pe_suff_chi`.

**Required change:**
Either (a) re-derive `Pe_suff_chi` and `Pe_fail_chi` in-script by solving the Stage-69 condition (e.g., `FindRoot[zetaF1[Pe] == zeta_suff_target_from_stage_63, ...]`); or (b) add a comment pointing to the exact upstream script that produces these numerics, and include a `# CARRY-FORWARD FROM stage_069_*.py output` annotation so the auditor can chase provenance.

**Verification:**
After fix, either the literals are absent and a derivation produces them, or each literal carries an in-line comment naming the producing script and the line of its output transcript.

## Independent-derivation check (Mathematica)

The `.wl` is **not** an independent re-derivation. Inspect the parallel:

- SymPy `kappa_F1 = sp.Rational(12321, 5)` (line 37) ↔ Mathematica `kappaF1 = 12321/5` (line 38)
- SymPy `y_F1 = sp.nsolve(y * sp.tan(y) - eta_F1, sp.Float("1.53", 80), ...)` (line 39) ↔ Mathematica `yF1 = y /. FindRoot[y Tan[y] == etaF1, {y, 1.53}, ...]` (line 40)
- SymPy `Omega = pi*Pe*(2*Pe*exp(Pe) + pi)/((4*Pe**2 + pi**2)*(exp(Pe)-1))` (line 41-44) ↔ Mathematica `omegaPe[pe_] := Pi pe (2 pe Exp[pe] + Pi)/((4 pe^2 + Pi^2) (Exp[pe] - 1))` (line 42)
- SymPy `Pe_suff_chi = 96.5285247264386` / `Pe_fail_chi = 11220.5441626259` (lines 49-50) ↔ Mathematica `peSuffChi = SetPrecision[96.5285247264386, 40]` / `peFailChi = SetPrecision[11220.5441626259, 40]` (lines 47-48)
- SymPy `Q = (1 + (1 - 2 eps_blk) zeta) / (1 - eps_blk zeta)` (line 53) ↔ Mathematica `q[zeta_, eps_] := (1 + (1 - 2 eps) zeta)/(1 - eps zeta)` (line 49)

This is a syntactic port. See F2.

## Engine cross-check

Both engines produce numerically matching outputs:

| Quantity | SymPy | Mathematica |
|---|---|---|
| rho_suff | 3.466222913478464001010330 | 3.466222913478464577913281378... |
| rho_fail | 3.467529132738703486267196 | 3.4675291327387033401575181822... |
| rho_max | 3.467529229456012233329585 | 3.46752922945601223332958450... |
| A_F1 | 1.000051928802195328659334 | 1.00005192880219532865933408... |
| zeta_max | 2.467529229456012233329585 | 2.46752922945601223332958450... |

Agreement to ~10^-15 (the limit set by Pe_suff_chi/Pe_fail_chi having only 12 significant figures). No engine_disagreement finding. But agreement is uninformative under F2 — both engines compute from identical literals.

## Verdict justification

The script holds up at the level of "minimal isotropic branch satisfies `zeta_min < A_F1` numerically with margin." It does **not** hold up at the level of the paper's boxed `Pe_req = 0` deliverable, because the link `Omega(Pe=0) = 1` is never verified (F1). The Mathematica script is not an independent derivation (F2). Three rho-anchor checks are tautological (F3), and two load-bearing Pe literals lack provenance (F4). The numeric ordering, the regime check (`1/3 < 1`), and the headline inequality (`1/3 < A_F1`) are substantive and pass. Verdict: `findings`, no stop_cold. Paper alignment is `partial` because four of six paper-side deliverables match and the headline boxed equation `Pe_req = 0` is only verified through a precondition.

## Self-test notes

- Variable independence: F1's proposed fix asks for `sp.limit(Omega, Pe, 0) == 1`. Omega does depend on Pe (the 0/0 limit at Pe=0 evaluates to 1 by L'Hopital); the limit is non-trivial, not identically zero. Check passes mentally.
- Symmetry/parity: not applicable — no integrals.
- Trivial-case pre-check: F3's proposed replacement uses notes-quoted numeric `3.46622291347846`; the existing output already prints `3.466222913478464...`, so the new check would pass.
- Paths: directive cites `scripts/...py` and `mathematica/...wl` correctly.
- Paper round-trip: F1's fix introduces a check (`Omega(0) = 1`) that exactly matches the paper card's prose ("\Omega_{\rm Pe=0}=1, so Family-1 supplies \zeta_{\rm F1}(0)=A_{\rm F1}"). No new misalignment introduced.
