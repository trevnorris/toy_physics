---
unit_id: 084
batch: III.4
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: missing
  mathematica: present
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [notes/stages/moving_throat_pde_stage084_full_reduced_pde_writeup.md]
  paper_appendix: present
---

# Audit unit 084 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_084.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage084_full_reduced_pde_writeup.md` (present)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 084 at line 146)
- sympy: (missing) — genuinely absent; see single-engine reasoning below
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.wl`
- sympy output: (missing)
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.txt`

## What the paper claims

Stage 084 is explicitly a **write-up skeleton / consolidation stage**, not a fresh derivation. The card's `\stagefield{Purpose}` reads: "Stage~084 organizes the completed reduced hierarchy into a write-up-ready theorem chain and isolates the final outgoing quadrupole product as the remaining gap." Its `\stagefield{Output}` reads verbatim: "A structured reduced-PDE theorem chain and handoff to the outgoing-branch loading-ratio problem." `\stagefield{Inputs}` is "Stages~037--083," and the card's `\claimstatus` is `\StatusExactClosure` for the reduced hierarchy organization but `\StatusOpen` for the final actual-branch realization. The notes elaborate the deliverables it consolidates (none new to this stage): (i) the explicit support family `zeta_phys(Pe,eta;kappa)=Omega_Pe^2(kappa+pi^2/4)/(kappa+y^2)` with `Omega_Pe = pi Pe(2 Pe e^Pe+pi)/[(4Pe^2+pi^2)(e^Pe-1)]` and `y tan y = eta` (§1.3); (ii) the selected-branch demand map `zeta_req=(Pi_tr-C_mix)/[C_mix-eps_blk(2C_mix-Pi_tr)]` with exact inverse `Pi_tr=C_mix Q(zeta;eps_blk)`, `Q=[1+(1-2eps_blk)zeta]/[1-eps_blk zeta]` (§1.4); (iii) the master residual `R_quad = zeta_req - zeta_phys` (§2); (iv) the Family-1 specialization `kappa_F1=12321/5`, `eta_F1=37`, `Xi_F1=1369 Upsilon_w = 136900 Theta_w` (§3); and (v) the numerical saturation picture — five carried-forward window/ceiling values (`zeta_-^chi≈2.46622291347846`, `zeta_+^chi≈2.46752913273870`, `zeta_-^J≈2.44257571477179`, `zeta_+^J≈2.46752736855058`, `zeta_max^F1≈2.46752922945601`) showing the explicit Family-1 branch is "effectively saturated." The remaining open datum (the outgoing quadrupole-normalization branch that fixes `Pi_tr`) is deliberately NOT proved here.

## What the script claims to verify

The Mathematica script defines the four symbolic objects (`zetaPhys`, `zetaReq`, `qMap`, `rQuad`) directly from the notes' closed forms, prints them, and runs four substantive checks: (A1) the demand map's inverse is genuine — `zeta_req(Pi_tr = C_mix·Q(zeta)) = zeta` identically; (A2) the Family-1 normalization double-equality is consistent — `1369·Upsilon_w` with `Upsilon_w→100 Theta_w` equals `136900·Theta_w`; (A3) the load-bearing one — independently recomputing `lim_{Pe→∞} zeta_phys(Pe, kappa_F1, eta_F1, y_F1)` (with `y_F1` solved from `y tan y = 37`) and matching it to the pinned ceiling `zeta_max^F1` to 1e-10; and (A4–A7) four ordering relations on the five pinned window/ceiling constants that encode the §3 saturation claim (`zeta_+^chi > zeta_-^chi`, `zeta_max^F1 > zeta_+^chi`, `zeta_+^J > zeta_-^J`, `zeta_+^J ≤ zeta_+^chi`). The script asserts no claim about the still-open outgoing-normalization branch, matching the card's `\StatusOpen` on that piece.

## Paper ↔ script cross-check

| paper-side deliverable | script-side check | status |
|---|---|---|
| `zeta_phys` closed form (§1.3) | defined lines 39–43, printed line 48; exercised by A3 limit | match |
| `zeta_req` demand map (§1.4) | defined line 44; exercised by A1 inverse check | match |
| `Q` inverse map (§1.4) | defined line 45; A1 verifies it is the true inverse of `zeta_req` | match |
| `R_quad = zeta_req - zeta_phys` (§2) | defined line 46, printed line 51 (symbolic only — no root claim, matching `\StatusOpen`) | match |
| Family-1 constants `kappa_F1, eta_F1, Xi_F1` double-eq (§3) | A2 (line 65) checks the `1369 Upsilon_w = 136900 Theta_w` consistency | match |
| `zeta_max^F1` ceiling reproduced from closed form (§3) | A3 (line 76) independent `Pe→∞` re-derivation matches pin to ~2e-15 | match |
| §3 saturation ordering of the 5 windows | A4–A7 (lines 84–87) | match |
| outgoing quadrupole-normalization branch (still open) | not tested (correct — card says `\StatusOpen`) | match (intentionally untested) |

`paper_alignment: aligned`. Every closed form and constant the script emits is present verbatim in the notes (and the card/appendix carry the same status), and the one item the paper marks open is correctly left unverified.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | mathematica | 53 | `expectZero[(zetaReq /. piTr->cMix*qMap) - zeta]` | §1.4 inverse demand map | yes |
| A2 | mathematica | 65 | `expectZero[(1369 upsilonW /. upsilonW->100 thetaW) - 136900 thetaW]` | §3 `Xi_F1` double-equality | yes (weak; arithmetic-level) |
| A3 | mathematica | 76 | `expectApprox[lim_{pe→∞} zetaPhys(F1), zetaMaxF1, 1e-10]` | §3 `zeta_max^F1` ceiling | yes (load-bearing) |
| A4 | mathematica | 84 | `expectZero[If[zetaPlusChi > zetaMinusChi, 0, 1]]` | §3 chi-window ordering | yes |
| A5 | mathematica | 85 | `expectZero[If[zetaMaxF1 > zetaPlusChi, 0, 1]]` | §3 ceiling above window | yes |
| A6 | mathematica | 86 | `expectZero[If[zetaPlusJ > zetaMinusJ, 0, 1]]` | §3 J-window ordering | yes |
| A7 | mathematica | 87 | `expectZero[If[zetaPlusJ <= zetaPlusChi, 0, 1]]` | §3 J below chi | yes |

## Findings

None. (clean)

## Independent-derivation check (Mathematica)

There is no SymPy script, so there is no transliteration to detect — this is a single-engine stage by design (`mathematica_transliteration` is N/A). The Mathematica script is not a port of any sibling script; its load-bearing check (A3) is a genuine independent re-derivation: it does not trust the pinned ceiling `zeta_max^F1 = 2.46752922945601`, it recomputes the `Pe→∞` limit of the §1.3 closed form at the Family-1 parameters — including independently root-solving `y tan y = 37` via `FindRoot` to get `y_F1 ≈ 1.5294824837…` — and only then compares. My own read-and-reason limit analysis (`Omega_Pe → pi/2` as `pe→∞`, so `zeta_phys → (pi^2/4)(kappa_F1+pi^2/4)/(kappa_F1+y_F1^2)`) reproduces `2.467529229456011`, matching the pinned ceiling to ~1.3e-15 and the script's reported diff of ~2.2e-15. This is the right way to anchor a carried-forward number.

## Engine cross-check

N/A — single engine. The stage card itself states "SymPy audit: none yet," the manifest records `sympy.exists: false` and `is_status_only_candidate: true`, and no `scripts/*stage084*` file exists. For a write-up/consolidation stage whose substantive content is checked through closed-form symbolic identities and numeric limit/ordering comparisons (rather than an independent algebraic re-derivation that two CAS engines could cross-validate differently), the single Mathematica engine is adequate; the only check that two engines would meaningfully duplicate is A3, and Mathematica already performs it from first principles. No `missing_sympy` finding is warranted: under the status-only rule, that finding requires the carry-forward to reference a result no upstream script verifies, which is not the case — `zeta_phys`/`zeta_req`/`Q`/`R_quad` and the Family-1 windows all originate in scripted upstream stages (037–083), and A3 re-anchors the one numeric ceiling here.

## Verdict justification

`clean`. I read the card, the notes, and the appendix row first, built the model that this is a consolidation/write-up skeleton (ExactClosure on organization, Open on the final outgoing branch), then attacked the script. Attacks tried and failed: (1) tautology hunt — A1 is a real inverse-map check (two independently written rationals shown to compose to the identity, would catch a sign/coefficient typo in either `zeta_req` or `Q`); A3 is a real re-derivation, not a number-against-itself; A4–A7 compare five distinct pinned constants and would fail on a mis-ordered carry-forward. (2) Hardcoded-value hunt — the five pinned windows ARE literals, but they are explicitly carried-forward upstream results, all five match the notes to every printed digit, the ceiling among them is independently re-derived by A3, and the rest are sanity-ordered; this is the correct treatment of carried constants, not a `hardcoded_result` defect. (3) Symbol-domain hunt — leaving `epsBlk` and `piTr` unrestricted-Real is correct (free demand parameters); the A1 simplification holds generically and the output confirms residual 0; `FindRoot` start `153/100` selects the correct first positive branch of `y tan y = 37`. (4) Paper-misalignment hunt — every symbolic form and constant the script emits is present verbatim in the notes; the one open item is correctly left untested. (5) Freshness — output mtime (02:19:50) is newer than script mtime (02:05:22), output content matches the current script line-by-line; not stale. The weakest check is A2 (arithmetic-level: it really just confirms `1369·100 = 136900` under the `Upsilon_w = 100 Theta_w` relation), but it is genuinely anchored to the notes' §3 double-equality and would catch a typo there, so it is a weak-but-valid check, not a finding.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every RESULT value the script emits (source: `.wl` source + committed `mathematica/output/...txt`; nothing executed):

| value | source (wl line / output line) | .tex/.md location | status |
|---|---|---|---|
| `zeta_phys = (pe^2 Pi^2 (2 e^pe pe+Pi)^2 (kappa+Pi^2/4))/((-1+e^pe)^2 (4pe^2+Pi^2)^2 (kappa+y^2))` | wl 43, out 5 | notes §1.3 lines 56–64 (same closed form, `Omega_Pe` expanded) | MATCH |
| `zeta_req = -(cMix-piTr)/(cMix-2cMix epsBlk+epsBlk piTr)` | wl 44, out 6 | notes §1.4 lines 73–74 (`(Pi_tr-C_mix)/[C_mix-eps_blk(2C_mix-Pi_tr)]`, sign-equivalent) | MATCH |
| `Q = (1+zeta-2 epsBlk zeta)/(1-epsBlk zeta)` | wl 45, out 7 | notes §1.4 line 80 | MATCH |
| `R_quad = zeta_req - zeta_phys` | wl 46, out 8 | notes §2 lines 90–92; card eq label `eq:app-stage084-chain` / appendix line 142 (`R_quad=zeta_req-zeta_phys`) | MATCH |
| `kappa_F1 = 12321/5` | wl 56, out 11 | notes §3 line 116 | MATCH |
| `eta_F1 = 37` | wl 56, out 12 | notes §3 line 117 | MATCH |
| `Xi_F1 = 1369 Upsilon_w` | wl 58, out 13 | notes §3 line 119 (`Xi_F1 = 1369 Upsilon_w`) | MATCH |
| `Xi_F1 = 136900 Theta_w` | wl 59, out 14 | notes §3 line 119 (`= 136900 Theta_w`) | MATCH |
| `zeta_phys(Pe→∞,F1) = 2.46752922945601…` | wl 75, out 19 | notes §3 line 144 (`zeta_max^F1 ≈ 2.46752922945601`) — same number, re-derived | MATCH |
| `zeta_-^chi = 2.46622291347846` | wl 67, out 22 | notes §3 line 132 | MATCH |
| `zeta_+^chi = 2.46752913273870` | wl 68, out 23 | notes §3 line 134 | MATCH |
| `zeta_-^J = 2.44257571477179` | wl 69, out 24 | notes §3 line 138 | MATCH |
| `zeta_+^J = 2.46752736855058` | wl 70, out 25 | notes §3 line 140 | MATCH |
| `zeta_max^F1 = 2.46752922945601` | wl 71, out 26 | notes §3 line 144 | MATCH |

INTERNAL scaffolding (accounted for, no finding): `lambdaEll = 37` (intermediate; `lambdaEll^2 = 1369`), `omegaPe` (intermediate building block of `zeta_phys`), `y_F1` (FindRoot intermediate ≈1.5294…), the per-check residual/diff printouts, and all PASS/FAIL flags and `Exit` codes.

`reconciliation: complete; 14 deliverable values checked, 0 misaligned`.

## Self-test notes

Checked: (1) variable independence — no spurious `D[...]`; the only derivative-like operation is the `Pe→∞` `Limit`, whose integrand/expression genuinely depends on `pe`, so the limit is non-trivial (confirmed it tends to `(pi^2/4)(kappa+pi^2/4)/(kappa+y^2)`, not 0). (2) Trivial-case / value pre-check — independently reproduced A3's `2.467529229456011` from the closed form and the `y tan y = 37` first root `y≈1.5295`, matching the pin to ~1e-15; confirmed A4–A7 orderings hold on the five pinned literals. (3) Branch choice — `FindRoot` start `153/100` lands on the correct first positive branch of `y tan y = 37`; no second-branch ambiguity. No directive written: zero findings, status-only single-engine stage by design, all 14 deliverable values reconcile.
