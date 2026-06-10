---
unit_id: 238
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 238 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_238.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows at lines 88, 1067-1099)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}`: "Physical transfer-shape compiler: the finite rigid-mouth packet is completely described by tracking, $\mathcal T^2$, and $\epsilon_\eta$." The derivation ledger states the stage proves the coherent-branch identity `R_target T^2 = Lambda0(1-eps_eta)`, factorizes the finite packet, and splits the branch gates into tracking, transfer-shape, and dressing conditions. The notes enumerate the deliverables in detail: (1) the boxed coherent identity; (2) finite packet factorization `q_nt + (B*/C*) q_tr = ln(T^2/T_ref^2)`; (3) the rigid-mouth reduction `q_nt = ln(T^2/T_ref^2)`, `q_eta = ln(eps_eta/eps_eta_ref)`; (4) the exact first-order physical drift compiler for `delta ln R_tr`, `delta ln T^2`, `delta ln R_target`; (5) the first-order packet relation and its rigid-mouth slice; (6) full support-blindness of the direct packet w.r.t. `zeta` and `M_mix` (the support channel enters only through `M_tr = M_mix[1+zeta(1-eps)/(1-zeta eps)]`); (7) the three-gate post-static dressing-invariance theorem (tracking gate, transfer-shape gate, dressing gate). The appendix (part07, lines 1067-1099) restates the same identity, factorization, rigid reduction, and the three gate equivalences. Card says "Mathematica audit: none yet" — but a `.wl` exists (card text-lag; not a math defect).

## What the script claims to verify

The SymPy docstring and assertions verify exactly the seven deliverables above, symbolically and exactly (no floats). The WL script mirrors the same seven (labeled M1-M7) but derives the first-order compiler from a distinct log-form path. Crucially, both engines implement the orchestrator's pass-1 reframing of support-blindness: rather than relying solely on vacuous `diff(observable, zeta)==0` passes, they add (a) a NEGATIVE CONTROL asserting `∂_zeta M_tr ≠ 0` and `∂_{M_mix} M_tr ≠ 0`; (b) a LEAK DETECTOR asserting `∂_zeta(R_tr·M_tr/M_mix) ≠ 0`; (c) a STRUCTURAL EXCLUSION asserting `{zeta, M_mix}` are absent from the reduced observables and packet coordinates. The vacuous `≡0` derivatives are retained but are now load-bearing only in conjunction with these three guards.

## Paper ↔ script cross-check

| Paper deliverable | Script check (py / wl) | Status |
|---|---|---|
| Coherent identity `R_target T^2 = Lambda0(1-eps_eta)` | py:62 / wl:67 | match |
| Finite packet factorization `q_nt+(B/C)q_tr=ln(T^2/T_ref^2)` | py:78-82 / wl:113-117 | match |
| Rigid-mouth `q_nt=ln(T^2/T_ref^2)`, `q_eta=ln(eps_eta/eps_eta_ref)` | py:85-90 / wl:119-130 | match |
| First-order drift compiler (R_tr, T^2, R_target) | py:96-137 / wl:139-168 | match |
| First-order packet relation + rigid slice | py:140-150 / wl:170-183 | match |
| Support-blindness w.r.t. zeta, M_mix (+ neg control/leak/exclusion) | py:155-202 / wl:80-111 | match |
| Three-gate theorem (tracking/transfer-shape/dressing) | py:207-223 / wl:185-201 | match |

`paper_alignment: aligned`. Naming differences only: notes `B_*,C_*` → py `B,C` / wl `Bstar,Cstar` (cosmetic, same role).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 62 | `simplify(Rtarget*T2 - Lambda0(1-eps_eta))==0` | coherent identity | yes |
| A2 | sympy | 82 | exp-ratio factorization ==0 | packet factorization | yes |
| A3 | sympy | 89-90 | rigid-mouth q_nt, q_eta ==0 | rigid reduction | yes |
| A4 | sympy | 106,118,137 | three drift compilers ==0 | first-order compiler | yes |
| A5 | sympy | 144,148,150 | first-order packet + rigid ==0 | first-order packet | yes |
| A6 | sympy | 158-167 | neg control: dMtr/dzeta, dMtr/dMmix nonzero (raise if 0) | support neg control | yes |
| A7 | sympy | 171 | leak detector: diff(Rtr·Mtr/Mmix,zeta)!=0 (raise if 0) | support leak | yes |
| A8 | sympy | 184-187 | structural exclusion zeta,M_mix ∉ free_symbols | support exclusion | yes |
| A9 | sympy | 189-202 | support-blind diffs ==0 (vacuous, backed by A6-A8) | support-blindness | yes* |
| A10 | sympy | 209-219 | tracking gate factorization + locus ==0 | tracking gate | yes |
| M1-M7 | math | 67-201 | mirror of A1-A10 incl. expectNonZero neg control/leak | all claims | yes |

\* A9 is individually vacuous (Rtr/T2/eps_eta carry no zeta/M_mix) but is load-bearing as a set with A6 (control), A7 (leak), A8 (exclusion).

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is an independent re-derivation, not a transliteration. Three concrete divergences from the `.py`:
1. **First-order compiler path differs.** SymPy perturbs the raw rational expressions and differentiates `log(simplify(Rtr_pert))` (py:96-102, 108-116, 120-129). Mathematica instead pre-expands each observable into an explicit sum of `Log[...]` (wl:132-137 `logRtr/logT2/logRtarget`), applies a single combined `perturbCore` rule (wl:139-146), and differentiates the log-sum (wl:148-150). Different algebraic route to the same drift.
2. **Rigid q_nt built independently.** SymPy obtains `q_nt_rigid` by `.subs({Rtr:Rtr_ref, ...})` into the finite `q_nt` (py:85-87). Mathematica reconstructs `qNtRigid` directly as `Log[(1-epsEta)/(1-epsEtaRef)] - Log[Rtarget/(...)]` under `PowerExpand` (wl:119-122), not by substitution into `qNt`.
3. **Exclusion via `FreeQ` vs `free_symbols`**, and Mathematica uses `PowerExpand`/`Together`+`FullSimplify` machinery (wl:6-10,114,119,195,199) absent in the SymPy path.
Not a `mathematica_transliteration`.

## Engine cross-check

Both engines pass every check. SymPy output: 32 `[ok]` lines, "All Stage 238 symbolic checks passed." Mathematica output: all residuals `= 0`, all `expectNonZero` witnesses genuinely nonzero (e.g. `M2 negative control zeta live witness = -(((-1 + eps)*Mmix)/(-1 + eps*zeta)^2)`; `M2 leak detector ... witness = -(((1 + chi0 + deltaU)*(-1 + eps))/((1 + chi0)*(1 + deltaU)*(-1 + eps*zeta)^2))`), "All Stage 238 Mathematica checks passed." Same identities, same forms, agree.

## Verdict justification

Clean. I attacked the batch hotspot (variable-independence self-test trap) specifically. The bare `∂_zeta(observable)==0` / `∂_{M_mix}(observable)==0` checks are individually vacuous — Rtr, T2, eps_eta, and the q-coordinates genuinely carry no `zeta`/`M_mix` — but they are NOT the only support-blindness evidence. All three reframed guards are present and load-bearing in BOTH engines: (1) the negative control `∂_zeta M_tr` / `∂_{M_mix} M_tr` differentiates `M_tr = M_mix[1+zeta(1-eps)/(1-zeta eps)]`, which genuinely depends on both variables; the asserted witnesses are real nonzero rationals, and the `if ==0: raise` (py:166) / `expectNonZero` (wl:85-86) would FIRE if `M_tr` were faked constant; (2) the leak detector `∂_zeta(R_tr·M_tr/M_mix)` is genuinely nonzero (the `M_tr/M_mix` factor carries zeta) and would fire if support leaked into R_tr or if M_tr were constant; (3) the structural exclusion confirms `zeta`/`M_mix` are absent from the reduced set. Every `sp.diff(EXPR, h)` / `D[expr, h]` in the first-order block hits `h`, which is present (perturbations inject `Exp[h·dln...]`); every `diff(M_tr, zeta/M_mix)` hits a present variable. No tautology survives without an exercised counterpart, no hardcoded numerics (fully symbolic), assumptions are sound (SymPy `positive` is weaker/safer than WL `0<eps<1`; no false-pass risk), and paper↔script alignment is exact. Outputs are fresh (both `.txt` mtime 08:31 > script mtimes 08:22/08:24). The only non-math note is the card's "Mathematica audit: none yet" text-lag, which is documentation, not a script defect.

## Self-test notes

Checked: (1) Variable independence — every `diff(...,h)` hits `h` (present via Exp[h·dln]); `diff(M_tr,zeta)` and `diff(M_tr,M_mix)` hit present variables, yielding real nonzero witnesses; the bare-observable `diff==0` checks are vacuous-by-design but are exactly the cases the negative control + leak detector + structural exclusion are there to rescue, and all three are present and load-bearing in both engines. (2) Trivial-case — the negative control and leak detector would genuinely FAIL if support-blindness were faked (constant M_tr ⇒ derivative 0 ⇒ raise/Exit[1]); confirmed. (3) First-order rigid reduction `qNtFirstRigid - dlnT2` algebraically reduces to 0 by hand. No symmetry/parity integrals in this stage. Concluded: zero findings.

## Value Reconciliation (pass-2 augmentation)

This stage is fully symbolic; it emits no named numeric constants (no `Pi_star`, `gamma_0`, benchmark numbers). All deliverables are closed-form symbolic identities, each of which is asserted `==0` (or `!=0` for controls) and carried verbatim in the notes/appendix.

| value (symbolic deliverable) | source (py / wl + output) | .tex/.md location | status |
|---|---|---|---|
| `R_target T^2 = Lambda0(1-eps_eta)` | py:62 / wl:67-70 (out M1 res=0) | notes:94-96 (boxed); appendix:1068-1071; card:13 | MATCH |
| `q_nt+(B*/C*)q_tr = ln(T^2/T_ref^2)` | py:78-82 / wl:113-117 (out M3) | notes:136-141; appendix:1074-1078 | MATCH |
| `q_nt = ln(T^2/T_ref^2)` (rigid) | py:85-89 / wl:119-126 (out M3) | notes:152-155; appendix:1080-1085 | MATCH |
| `q_eta = ln(eps_eta/eps_eta_ref)` | py:90 / wl:127-130 (out M3) | notes:155; appendix:1084 | MATCH |
| `delta ln R_tr` compiler form | py:103-106 / wl:154,163 (out M4) | notes:176-178 | MATCH |
| `delta ln T^2` compiler form | py:117 / wl:155-157,164 (out M4) | notes:181-186 | MATCH |
| `delta ln R_target` compiler form | py:130-137 / wl:158-161,165 (out M4) | notes:188-194 | MATCH |
| `q_nt+(B*/C*)q_tr = delta ln T^2` | py:144 / wl:175-178 (out M5) | notes:200-206 | MATCH |
| rigid first-order `q_nt=delta ln T^2`, `q_eta=d ln eps_eta` | py:148,150 / wl:179-183 (out M5) | notes:214-218 | MATCH |
| support-blindness (zeta, M_mix) | py:155-202 / wl:80-111 (out M2 block) | notes:240-266; appendix:1061-1064 | MATCH |
| three-gate theorem (tracking/transfer/dressing) | py:207-223 / wl:185-201 (out M6,M7) | notes:276-313; appendix:1086-1099 | MATCH |

INTERNAL scaffolding (no finding): `M_tr` baseline factor, the negative-control witnesses, the leak-detector witness, `h`-perturbation symbols (`dlnchi`, `dlndelta`, `dlnZW`, `dlnOm2`, `dlneps`, `dlnepseta`), `prefactor`/`trackingCondition` intermediates, all residual-near-zero/PASS scaffolding.

reconciliation: complete; 11 deliverable values checked, 0 misaligned.
