---
unit_id: 127
batch: IV.4
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage127_penetration_families.md]
  paper_appendix: present
---

# Audit unit 127 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_127.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage127_penetration_families.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only an `\input{stages/stage_127}` at line 1288; no separate summary row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage127_penetration_families_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage127_penetration_families_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage127_penetration_families_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage127_penetration_families_mathematica_audit.txt`

## What the paper claims

Stage 127 translates the lower compensated mouth-source bias `g_-^F1` (imported from Stages 125–126) into concrete geometric penetration scales for two positive localized source families. The card's bottom-line `Output` (quote block) is: "*Slab and truncated-exponential families reach the lower branch at moderate penetration depths.*" The notes give the load-bearing detail and five distinct deliverables: (1) the slab bias factor `g_slab(x) = (2/(πx)) sin(πx/2)` derived from the uniform slab source `σ_slab(z)=1/(xL)` on `[0,xL]`; (2) the slab depth `x*_slab ≈ 0.797839360904564` solving `g_slab(x)=g_-^F1`; (3) the truncated-exponential bias factor `g_exp(x) = 2(2+πx e^{-1/x})/((4+π²x²)(1-e^{-1/x}))` derived from `σ_exp(z)=e^{-z/(xL)}/(xL(1-e^{-1/x}))`; (4) the exp depth `x*_exp ≈ 0.662765402623160` solving `g_exp(x)=g_-^F1`; (5) the qualitative conclusion that both families reach the same lower branch at *moderate* (non-sign-changing, non-oscillatory) depths. The card carries an explicit non-theorem status tag (`StatusExactClosure`) and is a derivation-ledger entry, not an unconditional actual-branch theorem.

## What the script claims to verify

Both scripts (a) pin the imported lower-branch target `g_-^F1 = (2√(4107-100π²) - 37√3)/(20π) ≈ 0.758035…`; (b) state the slab and exp closed-form bias factors; (c) numerically solve `g(x)=g_-^F1` for each family near the seeds 0.8 / 0.66; and (d) assert the residual `g(x*) − g_-^F1` is below 1e-20. The Mathematica script does materially more than the SymPy: it *independently derives* both closed forms by symbolically integrating the source profiles against `cos(πz/2)` (`Integrate[...]`) and `FullSimplify`-comparing the result to the boxed closed form (a real, non-tautological check). The SymPy script states the closed forms directly without re-deriving them. The audit verdict applies to: existence of a compensation root for each family near the reported depth, correctness of the two closed forms, and engine agreement on all three numeric values.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| (1) slab bias `g_slab=(2/πx)sin(πx/2)` | `.wl` integrates `∫_0^x (1/x)cos(πz/2)dz` and `FullSimplify`-matches it to the boxed form (`.wl:36–41`); `.py` states it (`.py:18`) | match |
| (2) slab depth `x*_slab≈0.797839360904564` | root solved + residual<1e-20 (`.py:26,33,38`; `.wl:54,61`); printed `0.7978393609045644…` | match (15 dp identical) |
| (3) exp bias `g_exp=2(2+πx e^{-1/x})/((4+π²x²)(1-e^{-1/x}))` | `.wl` integrates `∫_0^1 σ_exp cos(πz/2)dz` and `FullSimplify`-matches boxed form (`.wl:44–49`); `.py` states it (`.py:22`) | match (hand-verified, see below) |
| (4) exp depth `x*_exp≈0.662765402623160` | root solved + residual<1e-20 (`.py:27,34,38`; `.wl:55,62`); printed `0.66276540262316140251…` | **mismatch** — notes 15th decimal `0` vs computed `1` (see F1) |
| (5) both reach SAME lower branch at moderate depth | both roots solved against the SAME `g_-^F1` literal; depths 0.80/0.66 ∈ (0,1] | match (existence + magnitude); "moderate"/"unique-depth" are qualitative, not asserted |

`paper_alignment: partial` — four of five deliverables match exactly; deliverable (4)'s notes value disagrees with the engines' computed value in the 15th decimal place.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | mathematica | 39–41 | `FullSimplify[∫_0^x (1/x)cos(πz/2)dz − 2sin(πx/2)/(πx)] === 0` | claim 1 (slab closed form) | yes (independent integral) |
| A2 | mathematica | 47–49 | `FullSimplify[∫_0^1 σ_exp·cos(πz/2)dz − g_exp] === 0` | claim 3 (exp closed form) | yes (independent integral) |
| A3 | mathematica | 61 | `Abs[g_slab(x*_slab) − g_-^F1] < 1e-20` | claim 2 (slab root exists) | partial (root re-substitution; convergence check) |
| A4 | mathematica | 62 | `Abs[g_exp(x*_exp) − g_-^F1] < 1e-20` | claim 4 (exp root exists) | partial (root re-substitution) |
| A5 | sympy | 38–39 | `raise unless |g_slab(x*)−g_-^F1|≤1e-20 and |g_exp(x*)−g_-^F1|≤1e-20` | claims 2 & 4 (roots exist) | partial (root re-substitution) |

A3/A4/A5 are "plug the root back in" checks: a converged root-finder makes the residual small by construction, so these are convergence-sanity checks rather than independent verifications of the depth values. They are not *purely* tautological because the closed forms they evaluate are independently established by A1/A2 (Mathematica side), so a wrong closed form would surface. Within the card's bottom-line claim (existence of a compensation root at moderate depth) they are adequate; they do not assert uniqueness of the root, but uniqueness of the *depth* is not the card's `Output` claim (the uniqueness statement in the notes/global conclusion is between-branch uniqueness, established upstream by 125–126).

## Findings

### F1 — paper_misalignment (value_mismatch)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage127_penetration_families.md:78` (boxed `x_*^{\exp}\approx 0.662765402623160`)
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage127_penetration_families_sympy_audit.txt:10` (`x_*^exp = 0.66276540262316140251…`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage127_penetration_families_mathematica_audit.txt:11` (`x_*^exp = 0.6627654026231614025…`)

**What's wrong:**
The notes box the exponential compensation depth as `x_*^{\exp}\approx 0.662765402623160` (final shown digit `0`). Both engines independently compute `x*_exp = 0.6627654026231614025…`. Aligning digit-by-digit after the decimal point, positions 1–14 agree (`66276540262316`), but the 15th decimal differs: the notes show `0` where the true value is `1` (followed by `40251…`). Truncating or rounding the engine value to 15 decimals gives `…623161`, not `…623160`. The slab counterpart `x_*^{\rm slab}\approx 0.797839360904564` agrees with the engines to 15 dp (engine's 16th digit `4`, rounds/truncates correctly), confirming the author intended 15-dp precision — so the exp value's trailing digit is an off-by-one in the last shown decimal, not a deliberate 14-dp truncation. The notes use `≈`, and the value is correct to 14 dp, so this is low severity, but it is a stated boxed deliverable whose printed value disagrees with the verified computation.

**Why this matters:**
A downstream reader who copies the boxed 15-digit value from the notes would carry a wrong final digit. The exhaustive pass-2 reconciliation is meant to surface exactly this kind of digit-level drift between the verified script output and the prose carrier.

**Required change:**
None to apply mechanically — this is a paper/notes ↔ script disagreement. Routes to the user via the directive's `## Resolve before fix_loop` block. Almost certainly the notes digit is the typo (engines agree to full precision and the script is the authority on the computed root), but the direction is the user's call. Codex must not silently edit either side.

**Verification:**
After user resolution: if the notes are corrected, line 78 should read `x_*^{\exp}\approx 0.662765402623161` (or a clean 14-dp value `\approx 0.66276540262316` that does not assert a wrong 15th digit). No script change is expected (the script value is correct).

## Independent-derivation check (Mathematica)

The `.wl` is **not** a transliteration of the `.py`. It performs strictly more work: it derives both bias factors from the source profiles via symbolic integration and `FullSimplify`-compares them to the boxed closed forms (`.wl:36–49`), whereas the `.py` only writes the closed forms as given (`.py:18,22`). Corresponding sections:
- `.wl:36` `gSlabFromIntegral = FullSimplify[Integrate[(1/x)*Cos[Pi*z/2], {z, 0, x}], …]` — has no `.py` analogue (`.py` has no `integrate` call).
- `.wl:44` `gExpFromIntegral = FullSimplify[Integrate[(Exp[-z/x]/(x*(1 - Exp[-1/x])))*Cos[Pi*z/2], {z, 0, 1}], …]` — again no `.py` analogue.
- The root-finding stanza (`.wl:54–55` `FindRoot` vs `.py:26–27` `nsolve`) is parallel but uses each engine's native solver with independent seeds, which is the expected two-engine pattern, not echoed algebra.
The Mathematica side is the stronger, genuinely independent engine here.

## Engine cross-check

Both engines agree to full shown precision:

| value | SymPy | Mathematica |
|---|---|---|
| `g_-^F1` | `0.758035078944662826919680890414` | `0.7580350789446628269196808904141104577505296234775693214353` |
| `x*_slab` | `0.79783936090456440508490881279972941789601018958900271741142733512870884256665312` | `0.7978393609045644050849088127997294178960101895890027174113` |
| `x*_exp` | `0.66276540262316140251258875013529026617855113761165291470424652438625727951066625` | `0.6627654026231614025125887501352902661785511376116529147062` |

Each SymPy value contains the Mathematica value as a leading substring (engines agree to the full Mathematica-reported width). Closed forms: SymPy `g_exp = (2*pi*x + 4*exp(1/x))/((pi**2*x**2 + 4)*(exp(1/x) - 1))` equals Mathematica `(2*(2*E^(1/x) + Pi*x))/((-1 + E^(1/x))*(4 + Pi^2*x^2))` (SymPy multiplied numerator and denominator by `e^{1/x}`); `g_slab` forms are identical up to syntax. No `engine_disagreement`.

## Verdict justification

Attacks tried that failed: (a) hand-rederived both source integrals — the slab and exp closed forms are exactly correct (the exp closed form was verified term-by-term against the standard `∫e^{az}cos(bz)dz` formula and matches the boxed form). (b) Checked whether the residual assertions are tautological — they are convergence-sanity checks, but the closed forms they evaluate are independently established by the Mathematica integral derivations, so the chain is sound for the card's existence-claim. (c) Checked the `g_slab` small-`x` limit (→1, physically sensible). (d) Confirmed both engines agree to full precision and the `.wl` is an independent (stronger) derivation, not a transliteration. (e) Confirmed outputs are fresh (both `.txt` newer than their scripts). The single finding is a low-severity 15th-decimal value drift between the engines' computed `x*_exp` and the boxed value in the notes — a `paper_misalignment` (value_mismatch) that routes to the user, so the verdict is `findings` with `needs_user_resolution: true`. Everything else holds up against the paper. I read the paper card, notes, and appendix, and the scripts' verified claim matches the card's `Output` for four of five deliverables exactly.

## Self-test notes

Checked: (1) variable independence — no `diff`/`D`; the `.wl` integrals are over `z` of genuinely `z`-dependent integrands, results depend on `x`, not vacuous. (2) Parity — finite domains `[0,x]`,`[0,1]`, no symmetric-domain trap. (3) Trivial-case — `g_slab→1` as `x→0` (sensible), residual tolerances (1e-20) clear the reported residuals (~1e-58 / ~1e-83) by huge margins, and the closed forms are independently verified so a vacuous pass on a wrong value is excluded. (4) No missing-script findings (both engines present), so no path-spec trap. (5) Paper round-trip — the only finding prescribes NO Codex script edit (user-gated), so it cannot introduce a new misalignment; the script value is correct and engines agree.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every RESULT value the scripts emit, located in the `.tex` card and `.md` notes:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `g_slab(x) = 2 sin(πx/2)/(πx)` | py:18 / wl:37,51; sympy-out:5, math-out:7 | notes:29 (boxed `\frac{2}{\pi x}\sin(\pi x/2)`) | MATCH |
| `g_exp(x) = 2(2+πx e^{-1/x})/((4+π²x²)(1-e^{-1/x}))` | py:22 / wl:45,52; sympy-out:6, math-out:8 | notes:66–67 (boxed) | MATCH (hand-verified + `.wl` integral check) |
| `x*_slab ≈ 0.797839360904564` | py:26,30 / wl:54,58; sympy-out:9, math-out:10 | notes:41 (boxed `0.797839360904564`) | MATCH (15 dp identical) |
| `x*_exp ≈ 0.662765402623161…` | py:27,31 / wl:55,59; sympy-out:10 (`…623161 40251…`), math-out:11 | notes:78 (boxed `0.662765402623160`) | **MISMATCH** → F1 (15th decimal `0` vs `1`) |
| `g_-^F1 = (2√(4107-100π²)-37√3)/(20π) ≈ 0.758035078944662826919680890414` | py:14–15,29 / wl:31–32,57; sympy-out:8, math-out:9 | not stated in stage-127 `.tex` or `.md` (imported input from Stages 125–126; referenced only symbolically as `\mathfrak g_-^{F1}`) | IMPORTED-INPUT (not a stage-127 deliverable; per augmentation guard, not MISSING) |

INTERNAL scaffolding (accounted for, no finding): `R = √(4107-100π²)` / `rDisc` (intermediate of `g_-^F1`); root-finder seeds 0.8 / 0.66; residuals `res_slab ≈ -2.7e-83`, `res_exp ≈ -9.7e-82` (sympy), `1.17e-58`, `-3.25e-58` (mathematica); tolerance `1e-20`; `slabClosedFormResidual`/`expClosedFormResidual` (=0 flags); PASS/FAIL banner strings.

`reconciliation: complete; 5 deliverable values checked (4 MATCH + 1 IMPORTED-INPUT), 1 MISMATCH (x*_exp 15th decimal → F1)`.
