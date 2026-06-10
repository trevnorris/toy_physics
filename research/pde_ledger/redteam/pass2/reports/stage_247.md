---
unit_id: 247
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-10T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 247 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_247.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex` (row 92 + theorem block lines 117/138/242/249)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_mathematica_audit.txt`

## What the paper claims

Stage 247 assembles the first stationary lowered-barrier law from a carried one-port short-range baseline plus three imported "relaxed" packets (Stage-244 leakage/work, Stage-245 U/V drain, Stage-246 compensated source). The bottom-line deliverable (notes §6, appendix theorem item 6) is the compiler `V_eff^(247)(r) = V_short^(1p)(r) - lambda_L S_leak - lambda_W W_w^sess - DeltaE_UV - M_sigma` plus the exact lowering identity `V_short^(1p) - V_eff^(247) = lambda_L S_leak + lambda_W W_w^sess + DeltaE_UV + M_sigma` and the inequality `V_eff^(247) <= V_short^(1p)` when packet scalars/weights are nonnegative. Distinct deliverables: (1) one-port baseline `V_short^(1p)` with `A_6=3 alpha_6 + C_6/2`, `A_4=C_4`, `A_2=alpha_2 + C_2/2` and susceptibilities = entries of the inverse reduced stiffness matrix; (2) Stage-244 leakage/work scalars; (3) Stage-245 weighted drain `DeltaE_UV = eta_UV D_UV`, `D_UV>=0`; (4) Stage-246 source-response `M_sigma = xi_R[R_inf - R(r)]` with `g_inf=2/pi`, `g(r)>=2/pi`, `M_sigma>=0` on the Session-I branch; (5) lowering identity + inequality; (6) Session-I benchmark at `r_soft=0.18`: `3.74163698 - 0.26971918*0.31069599 - 1.51632107 - 0.21064278 - 0.18386120 = 1.74701126` (card eq:benchmark-decomposition lines 79-92). The **published card body carries only the abstract benchmark decomposition** (no `Delta` numeric); the benchmark intermediate `Delta=142.17750000` and `D0=3.76481862` live in the **notes** (md:406-409).

## What the script claims to verify

The SymPy docstring (py:5-14) enumerates six items mirroring the paper. The assertions test: (Sec 1) six susceptibilities equal the matching entries of `K_red.inv()` (py:60-65); (Sec 2) `S_leak`/`W_sess` match a factored form through `E0` (py:120-121) and `Lvar` is recoverable from `W_sess` (py:122); (Sec 3/4) `D_UV>=0` symbolically + numeric probe (py:133,140), `g_inf==2/pi` (py:167); (Sec 5) the lowering gap equals the symbolic packet sum (py:181); (Sec 6) numeric benchmark: `Delta=142.1775` (py:266), `D0~3.76481862` (py:267), the **forward** decomposition with the paper literal `lambda_L=0.26971918` reproduces `1.74701126` (py:248-250), plus pinned figures for `V_short`, `M_sigma`, `S_soft`, `Lvar_soft`, `lambda_L_soft` against paper values (py:235-243), branch positivity `g_soft`, `M_sigma>=0` (py:223-225), and the lowering inequality `V_short - V_eff_session >= 0` (py:247). The `.wl` mirrors all eight checkpoints M1-M8 with independent techniques.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) susceptibilities = inverse-matrix entries | py:60-65 / wl:128 (`Inverse[Kred]` vs closed forms) | match |
| (1) `A_6,A_4,A_2`, `V_short` | py:71-80 / wl:136-144 | match |
| (2) `S_leak`,`W_sess` factoring | py:120-121 / wl:152-159 | match |
| (3) `DeltaE_UV=eta_UV D_UV`, `D_UV>=0` | py:133,140 / wl:177-188 (Reduce) | match |
| (4) `g_inf=2/pi`, `g>=2/pi`, `M_sigma>=0` | py:167,223-225 / wl:170,261,264 | match |
| (5) lowering identity | py:181 / wl:198 | match |
| (5) inequality `V_eff<=V_short` | py:247 (benchmark slice) / wl (M8 numeric forward) | match (numeric on branch) |
| (6) benchmark `=1.74701126`, `Delta=142.1775` | py:240-250,266 / wl:237,276 | match |

`paper_alignment: aligned`. Card body is abstract (no `Delta` numeric); the `Delta=142.17750000` deliverable lives correctly in the notes and is reproduced by both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1-6 | sympy | 60-65 | `simplify(K_red_inv[i,j]-chi)==0` | claim 1 | yes |
| A7-8 | sympy | 120-121 | `simplify(S_leak-S_exp)==0`, `W_sess-W_exp` | claim 2 | yes |
| A9 | sympy | 122 | `simplify(Lvar_from_W - Lvar)==0` | claim 2 | no — tautological round-trip (non-load-bearing, see note) |
| A10-11 | sympy | 133,140 | `ask(nonnegative(D_UV))`, probe>=0 | claim 3 | yes |
| A12 | sympy | 167 | `simplify(g_inf-2/pi)==0` | claim 4 | yes |
| A13 | sympy | 181 | `expand(lowering_gap-lowering_expected)==0` | claim 5 | yes (definitional but correct) |
| A14-16 | sympy | 223-225 | `g_soft>=2/pi`, `g_soft<rF1`, `M_sigma>=0` | claim 4 | yes |
| A17-18 | sympy | 235-236 | `Wsess_from_Lvar≈Wsess_obs`, `Lvar_soft≈20.01677473` | claim 2/6 | yes (substantive, added pass-1) |
| A19-22 | sympy | 240-243 | pins `V_short,M_sigma,S_soft,lambda_L_soft` to paper | claim 6 | yes |
| A23 | sympy | 247 | `V_short - V_eff_session >= 0` | claim 5 (inequality) | yes |
| A24 | sympy | 250 | `Veff_forward(lambda_L_paper)≈Veff_obs` | claim 6 (forward closure) | yes — falsifiable |
| A25-26 | sympy | 266-267 | `Delta_num≈142.1775`, `D0≈3.76481862` | claim 6 | yes |
| M1-M8 | math | 128,152-198,237-276 | parallel `expectZero`/`expectApprox`/`expectTrue` | claims 1-6 | yes |

## Findings

None. The pass-1 findings are all resolved:
- **F1 (paper Delta=210.1775 typo)** — RESOLVED. `210.17750000` no longer appears anywhere (card/notes/scripts grep empty). Notes md:406 now reads `\Delta=142.17750000`, matching both engines. The published **card body is abstract** (lines 79-92, no `Delta` numeric) so it was never the carrier — the pass-1 misattribution to `stage_247.tex:407` was correctly overruled by the orchestrator; the card is clean and is NOT flagged.
- **F2 (Lvar inverse round-trip)** — DOWNGRADED to benign. The symbolic round-trip at py:122 still exists but is no longer load-bearing: the substantive numeric checks A17/A18 (py:235-236) were added and genuinely test that `Lvar_soft` reproduces `Wsess_obs` and matches `20.01677473`. Not a finding (a kept-but-harmless symbolic identity alongside real checks).
- **F5 (lambda_L self-solve resubstitution)** — RESOLVED. `Vrebuild_soft` (py:251) is computed but NEVER asserted; the load-bearing benchmark closure is A24 (py:248-250) using the **paper literal** `lambda_L_paper=0.26971918`, which is genuinely falsifiable.
- **F3/F4 (missing positivity + inequality)** — RESOLVED via A14-16, A23.
- **F6 (missing mathematica)** — RESOLVED; `.wl` present with independent M1-M8.

## Independent-derivation check (Mathematica)

The `.wl` is a genuine parallel re-derivation, not a transliteration. Evidence: (1) susceptibilities are checked by `Inverse[Kred]` minus closed forms (wl:104-128) — same physical premise, but the drain nonnegativity uses `Reduce[... Duv < 0, {aU,aV,chiUV,fU}, Reals]` returning `False` (wl:177-187), a logically distinct technique from SymPy's `sp.ask(Q.nonnegative)` + float probe (py:133-140). (2) The benchmark substitution uses **exact rationals** (`GW->5/4`, `Rmix->27/20`, `rF1->88899676773749/50000000000000`, wl:208-219) and `N[...,30]`, an independent numeric path vs SymPy `sp.Float`. (3) The `.wl` benchmark only performs the **forward** decomposition M8 with `lambdaLPaper=13485959/50000000` (=0.26971918, wl:269-276); it does NOT replicate the SymPy self-solve `lambda_L_soft`, so it cannot inherit that round-trip. Variable choreography and section ordering parallel the physics, not the SymPy line sequence.

## Engine cross-check

Both engines agree to high precision. SymPy out: `Delta=142.1775000000000`, `D0=3.764818624606566`, `V_short=3.741636979285524`, `M_sigma=0.1838612007190837`, `Lvar=20.01677472685125`, `S_leak=0.3106959887767410`, `lambda_L=0.2697191840048420`, `V_eff rebuilt=1.747011260000000`. Mathematica out: `Delta=142.1775`, `D0=3.76481862460656…`, `V_short=3.74163697928552…`, `M_sigma=0.18386120071908…`, `Lvar=20.0167747268512…`, `S_leak=0.31069598877674…`, `V_eff forward=1.74701126124428…` (diff 1.2e-9 vs target). All PASS in both transcripts; no engine_disagreement.

## Verdict justification

Clean. The pass-1 round-trip findings are genuinely cured: the load-bearing benchmark closure (py:248-250 / wl:269-276) now uses the paper-stated `lambda_L=0.26971918` rather than a self-solved root, so it is falsifiable; the surviving symbolic `Lvar_from_W` round-trip (py:122) is non-load-bearing scaffolding sitting beside real numeric checks (py:235-236). The paper Delta typo is fully purged (210.1775 absent everywhere); 142.17750000 holds in the notes and both engines, and the published card body is correctly abstract (not flagged). Both engines present and independent (Reduce vs ask; exact rationals vs floats). Attacks tried that failed: (a) searched for a surviving load-bearing tautology — the only tautology (py:122) is shadowed by substantive checks; (b) checked the forward `Veff` assert for a hidden self-solve — it uses the paper literal, confirmed by hand: `3.74163698 - 0.26971918*0.31069599 - 1.51632107 - 0.21064278 - 0.18386120 ≈ 1.7470113`; (c) probed `.wl` for line-by-line porting — the drain branch and numeric path diverge in method. I confirm I read the card, notes §1-6, and the part-08 appendix rows before opening the scripts.

## Value Reconciliation (pass-2 augmentation)

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Delta=142.17750000` | py:266 / wl:237 / out py L57, wl L46 | notes md:406 | MATCH |
| `D0=3.76481862` | py:267 / wl:238 / out py L58 | notes md:408; card eq has it only abstractly | MATCH (notes) |
| `V_short=3.74163698` | py:240 / wl:239 / out py L59 | notes md:412; card line 81 (`3.74163698`) | MATCH |
| `M_sigma=0.18386120` | py:241 / wl:260 / out py L60 | notes md:483; card line 89 | MATCH |
| `Lvar=20.01677473` | py:236 / wl:262 / out py L61 | notes md:457 | MATCH |
| `S_leak=0.31069599` | py:242 / wl:263 / out py L62 | notes md:463; card line 82 | MATCH |
| `lambda_L=0.26971918` | py:244 / wl:269 / out py L65 | notes md:518; card line 82 | MATCH |
| `V_eff=1.74701126` | py:250 / wl:276 / out py L66 | notes md:372,534; card line 91 | MATCH |
| `V_Coul=5.55555556` (context) | not asserted in script | notes md:374 | INTERNAL (notes-only context) |
| residual `2.01467313` | out py L63 | notes md:554 | MATCH |
| residual `1.83081193` | out py L64 | notes md:563 | MATCH |

INTERNAL scaffolding (no finding): pass/fail flags, `E0` intermediate, `g_soft=1.10521636`, `sigma_min`, `D_UV` probe `0.01866637`, ratio `0.31446203`/`0.67349466` (notes-only diagnostics), residual `0.08380067` (notes md:568).

reconciliation: complete; 11 deliverable values checked, 0 misaligned.

## Self-test notes

Variable-independence trap: no `sp.diff`/`D` against an absent variable anywhere; positivity/limit checks (`g_r`, `M_sigma`, `D_UV`, `g_inf`) operate on variables the expressions genuinely contain. Symmetry/parity: no unbounded-domain integrals in this stage. Trivial/forward pre-check: hand-evaluated the forward closure with the paper literal `lambda_L=0.26971918` → `1.7470113`, matching `Veff_obs` and confirming A24 is satisfiable and non-vacuous; verified `Delta=9*16-1.35^2=142.1775` and that `D0=4-33.4375/142.1775=3.76483` (matches notes 3.76481862), excluding any 210.1775 residue.
