---
unit_id: 191
batch: V.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage191_minimal_pde_data_packet.md]
  paper_appendix: present
---

# Audit unit 191 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_191.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage191_minimal_pde_data_packet.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 113, 1139, 1216–1229 reference this stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage191_minimal_pde_data_packet_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage191_minimal_pde_data_packet_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

Stage 191 (`\StatusExactClosure`) is a finite-packet compiler card. Its `\stagefield{Output}` reads verbatim: "Defines Packet A, Packet B, \(\Delta_{\rm branch}\), \(\Delta_{\rm orbit}\), and the exact two-packet home-stretch theorem." The notes enumerate the deliverables: (1) the exact single-lane compiler from grouped-lane low-frequency data `(D0,D2,D4,N0,N2,N4)` to the response moments `nu_2 = -D_2/D_0`, `nu_4 = (D_2^2 - D_0 D_4)/D_0^2` and the prefactor moments `P_0=N_0/D_0`, `P_2=(D_0 N_2 - 2 D_2 N_0)/D_0^2`, `P_4=(D_0^2 N_4 - 2 D_0(D_2 N_2 + D_4 N_0)+3 D_2^2 N_0)/D_0^3`; (2) the grouped weighted trace/anomaly map `(xbar, a_x, b_x)` with its exact inverse `x20=xbar+4a_x`, `x21=xbar-a_x+b_x`, `x22=xbar-a_x-b_x`; (3) the one-pole defect `Delta_pole := bar_nu_4 - 4 bar_nu_2^2` and normalization defect `Delta_norm := mhat0^2 bar_P0 - 54 G c_s^5/(5 a^5 c^5)`; (4) the branch residual packet `Delta_branch = (a2,b2,a4,b4,aP0,bP0,Delta_pole,Delta_norm)`, which must vanish on the isotropic one-pole normalized branch (appendix eq. line 1233 gives the equivalent one-pole surface `D_0 D_4 + 3 D_2^2 = 0`); (5) Packet B orbit interconversion `R_tr=m_T^{1+chi}`, `R_nt=m_mu/(m_K m_T^F)`, `R_eta=1/m_K`, `q=log R`, with inverse `m_T=e^{q_tr/(1+chi)}`, `m_K=e^{-q_eta}`, `m_mu=e^{q_nt-q_eta+F q_tr/(1+chi)}`; and (6) the zero-set equivalence `Delta_orbit=0 ⟺ m=1 ⟺ R=1`. The stage is not a checkpoint and is not status-only; it imports forwards from Stages 239–241 but the identities it states are self-contained compiler algebra it verifies in-script.

## What the script claims to verify

The SymPy script verifies, symbolically: the single-lane response/prefactor compilers (Taylor expansion of `D0/D(ω)` and `D0 N(ω)/D(ω)^2` against the claimed closed-form moment coefficients); the one-lane `Delta_pole = -(3 D2^2 + D0 D4)/D0^2` identity; the grouped trace/anomaly inverse round-trip; the G-orthogonality of the three basis directions under metric `diag(1,2,2)` and the completeness/idempotency/orthogonality of the three projectors `Pbar,Pa,Pb`; isotropic collapse of all six anisotropy coordinates `(a2,b2,a4,b4,aP0,bP0)` and reduction of the grouped means to the single-lane moments; vanishing of the full `Delta_branch` on the isotropic one-pole (`D4c=-3 D2c^2/D0c`) normalized (`N0c=D0c·54Gc_s^5/(5a^5c^5)/mhat0^2`) branch; and the orbit-packet interconversion round-trips plus the orbit-lock zero-set. Every assertion routes through `expect_zero`, which simplifies and raises on nonzero.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| nu_2, nu_4 closed forms | A1 response series compiler (L87) | match |
| P_0,P_2,P_4 closed forms | A2 prefactor series compiler (L100) | match |
| one-pole defect form | A3 Delta_pole one-lane (L105–108) | match |
| trace/anomaly inverse map | A4–A6 x20/x21/x22 inverse (L125–127) | match |
| weighted-projector identities | A7–A16 G-orthogonality, projector algebra (L137–146) | match |
| isotropy gate (anisotropy → 0) | A17–A22 (L190–195) + A23–A25 mean→single-lane (L199–201) | match |
| Delta_branch=0 on iso one-pole normalized branch | A26 (L210) | match |
| orbit interconversion R(m), q(R), m(q) | A27–A29 round-trips (L257–259) | match |
| orbit-lock zero-set m=R=1, q=0 | A30–A32 (L260–262) | match |
| Packet A/B enumeration, home-stretch theorem | Section V print ledger (L271–292) | match (prose summary; substantive checks above) |

`paper_alignment: aligned`. Every paper deliverable has a faithful, non-tautological script-side check using the exact constants and forms the paper states; no orphaned assertions.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 87 | `Y_series - Y_expected == 0` | nu_2,nu_4 compiler | yes |
| A2 | sympy | 100 | `Pref_series - Pref_expected == 0` | P_0,P_2,P_4 compiler | yes |
| A3 | sympy | 105–108 | `Delta_pole_one + (3D2^2+D0D4)/D0^2 == 0` | one-pole defect form | yes |
| A4–A6 | sympy | 125–127 | `x_rec - x == 0` | trace/anomaly inverse | yes |
| A7–A9 | sympy | 137–139 | `e_i^T G e_j == 0` | weighted-projector orthogonality | yes |
| A10 | sympy | 140 | `Pbar+Pa+Pb-I == 0` | projector completeness | yes |
| A11–A13 | sympy | 141–143 | `P_i^2 - P_i == 0` | projector idempotency | yes |
| A14–A16 | sympy | 144–146 | `P_i P_j == 0` | projector orthogonality | yes |
| A17–A22 | sympy | 190–195 | anisotropy coords on iso branch `== 0` | isotropy gate | yes |
| A23–A25 | sympy | 199–201 | grouped mean − single-lane moment `== 0` | mean reduces to moment | yes |
| A26 | sympy | 210 | `Delta_branch` on iso one-pole normalized branch `== 0` | Delta_branch zero-set | yes |
| A27 | sympy | 257 | `R_from_m∘m_from_q − R_from_q == 0` | orbit interconversion (composition) | yes |
| A28 | sympy | 258 | `q_from_m∘m_from_q − q == 0` | inverse map consistency | yes |
| A29 | sympy | 259 | `q_from_R∘R_from_q − q == 0` | q=log R consistency | partial (log/exp inverse only) |
| A30–A32 | sympy | 260–262 | orbit-lock substitutions `== 0` | zero-set equivalence | yes |

No row is tautological in the harmful sense: each "claimed" closed form is constructed and then checked against an *independent* computation (Taylor series, log/exp inversion, projector algebra, or a guided substitution onto a defining surface). A29 is the weakest (just `log(exp(q))=q`) but it is harmless and corroborated by A27/A28.

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** target_mismatch (cosmetic stage-label mismatch; no verified identity differs)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage191_minimal_pde_data_packet_sympy_audit.py:65` quote: `banner("STAGE 174 — MINIMAL PDE DATA PACKET AND THE EXACT HOME-STRETCH THEOREM")`
- same file `:284` quote: `banner("STAGE 174 LEDGER")`
- paper side `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_191.tex:1` quote: `\section[Stage 191]{Stage 191: Minimal PDE data packet}`
- paper side `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex:113` quote: `191 & Minimal PDE data packet & ...`

**What's wrong:**
The script's two `banner(...)` print strings label the stage "STAGE 174", while the paper card and the part-05 appendix consistently number this unit Stage 191. The source notes additionally call the theorem the "Stage 242 home-stretch theorem" (notes line 311) and reference `..._stage242_...` supporting files (notes lines 421–422). These are legacy renumbering artifacts: the *mathematical content* of the script matches the Stage 191 deliverables exactly (every identity, constant `54 G c_s^5/(5 a^5 c^5)`, and form lines up). Only the human-readable banner labels carry a stale stage number. No verified assertion, constant, or deliverable differs between paper and script — the mismatch is confined to two print labels in the script's transcript header.

**Why this matters:**
The committed output (`...sympy_audit.txt`) is headed "STAGE 174" while the file is `stage191`. A future auditor cross-referencing the transcript against the paper card by stage number will hit a false signal that the transcript belongs to a different stage, undermining traceability of the audit trail. This is purely cosmetic — it changes no math — but the script-side label is the only place the audit trail can be made self-consistent, since Codex must not edit the notes prose.

**Required change:**
In the script, change the banner label string at line 65 from `STAGE 174 — MINIMAL PDE DATA PACKET AND THE EXACT HOME-STRETCH THEOREM` to `STAGE 191 — MINIMAL PDE DATA PACKET AND THE EXACT HOME-STRETCH THEOREM`, and at line 284 from `STAGE 174 LEDGER` to `STAGE 191 LEDGER`. Touch only these two string literals; do not alter any symbol, assertion, constant, or computation. Do NOT edit the notes file (the `242`/`174` legacy references there are prose and out of red-team scope).

**Verification:**
After the edit and a fresh run, the regenerated `scripts/output/moving_throat_pde_stage191_minimal_pde_data_packet_sympy_audit.txt` header banners read "STAGE 191 ..." and "STAGE 191 LEDGER", and the script still exits 0 with all `expect_zero` checks passing (the change is string-only and cannot affect any residual).

## Independent-derivation check (Mathematica)

No `.wl` exists for this unit, so transliteration is not in play.

## Engine cross-check

Only one engine (SymPy) is present; no cross-check applies. On the `missing_mathematica` question, I apply the prompt's line-114 judgment and do NOT flag it. Every claim this stage makes is a pure exact symbolic identity: a Taylor-series coefficient match, a closed-form log/exp inverse round-trip, projector algebra over a fixed `diag(1,2,2)` metric, and a guided substitution onto a defining surface. There is no numerical evaluation, no branch ambiguity (the symbol assumptions — `nonzero,real` on the D/N coefficients, `positive` on `chi0_star,F_star,m_*,R_*` — pin reality and positivity so `log(exp(real))=real` and `1+chi≠0` hold), and no transcendental special-function corner where an independent CAS would plausibly catch a SymPy defect. SymPy fully and genuinely settles each identity; the saved output shows every residual is exactly 0. This matches the established SymPy-only non-status-only precedent (stages 121/122/123). I cannot point to any claimed result that SymPy fails to genuinely verify, so a `missing_mathematica` finding would be unfounded.

## Verdict justification

Verdict is `findings` for one low-severity cosmetic `paper_misalignment` (stale "STAGE 174" banner labels in the script vs. paper-card Stage 191). The mathematics holds up under attack: I verified the prefactor `P_2` coefficient by hand against `d/d(ω^2)` of `D_0 N/D^2` at 0 (matches); confirmed `nu_4` as the ω^4 Taylor coefficient of `1/(1+u)` (matches); checked the `Delta_pole` reduction `nu_4 - 4 nu_2^2 = -(3 D_2^2 + D_0 D_4)/D_0^2` (matches and reduces to 0 exactly under the one-pole condition `D_4=-3 D_2^2/D_0`); confirmed `Delta_norm` cancels under the normalization substitution; confirmed the trace/anomaly inverse round-trip coefficients (the `4 a_x` factor is the genuine inverse, not a tautology); and walked the orbit interconversion composition `R_from_m∘m_from_q = R_from_q` by hand for the nontrivial `R_nt` component (exponents cancel to `e^{q_nt}`). None of the checks pass for the wrong reason: each "claimed" form is built once and tested against an independent computation, and the one guided-substitution check (`Delta_branch=0`) is the legitimate verification of the paper's stated zero-set surface, not a circular re-substitution of a value into its own defining equation. The anisotropy/orbit-lock substitution checks (A17–A22, A30–A32) are individually somewhat trivial but collectively exercise the stated isotropy and orbit-lock gates and would catch structural errors. I read the paper card, the notes, and the appendix rows, and the script's verified claims match the paper's deliverables exactly. Output is fresh (output 2026-05-11 newer than script 2026-04-11), so no `stale_output`. Not stop-cold: the lone finding is a string label with zero mathematical propagation, so no downstream result is affected.

## Self-test notes

Variable-independence trap: there are no `sp.diff(EXPR, VAR)` calls in this script (it uses `sp.series` and algebraic substitution), so the zero-derivative false-pass mode does not apply. Symmetry/parity trap: no unbounded integrals here. Trivial-case pre-check: I mentally substituted the iso one-pole normalized profile into `Delta_branch` and confirmed each of the eight components reduces to exactly 0 only because the defining conditions force it (not vacuously), and confirmed the orbit-lock substitutions give the literal identity `(1,1,1)`/`0`. F1's required change is a pure string edit and cannot introduce a new `paper_misalignment` or alter any residual; paper round-trip confirms the constant `54 G c_s^5/(5 a^5 c^5)` in the script already matches the notes verbatim, so no constant is touched.
