---
unit_id: 058
batch: III.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage058_coupled_support_source_operator.md]
  paper_appendix: present
---

# Audit unit 058 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_058.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage058_coupled_support_source_operator.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 94)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage058_coupled_support_source_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.txt`

## What the paper claims

Stage 058 builds the first explicitly **coupled** axial support/source operator so the transport bias `Pe` is no longer a free knob but the root of a fixed-point equation generated from the same support lane that carries `eta` and `kappa`. The card's `\stagefield{Output}` is verbatim: *"Fixed-point equation \eqref{eq:app-stage058-Pe-fixed} and bracket \eqref{eq:app-stage058-bracket}."* The boxed claims are the fixed point `Pe_* = Xi·Delta(Pe_*; kappa, eta)` (eq. app-stage058-Pe-fixed) and the operator-selected bracket `Xi·Delta_0 <= Pe_* <= Xi·Delta_inf` (eq. app-stage058-bracket), with the two boxed endpoint functions (eqs. app-stage058-Delta0 and Delta_inf, writing `alpha = sqrt(kappa)`). The notes add the full chain that supports these: (1) the minimal coupled operator and its nondimensionalization to `kappa, eta, Pe, Xi`; (2) the zero-flux source family `Sigma_Pe(x) = Pe·exp(Pe·x)/(exp(Pe)-1)`; (3) the support-drop kernel `K_(kappa,eta)(x)` and its closed-form `Delta` via the auxiliary integrals `Ic, Is`; (4) strict monotonicity of the kernel and of `Delta` in `Pe` (covariance identity); (5) endpoint values `Delta_0 = Delta(0)` and `Delta_inf = K(1) = lim_{Pe->oo} Delta`; (6) the IVT bracket-existence argument via `F(Pe) := Pe - Xi·Delta`; and (7) the weak-coupling branch law `Pe_*(Xi) = Xi·Delta_0 + O(Xi^2)`. The appendix row 94 summarizes it as the self-consistent equation plus bracket. There are no numeric deliverable constants — the entire stage is symbolic in `alpha, eta, Pe, Xi`.

## What the script claims to verify

Both scripts verify the full notes chain symbolically. They confirm: the kernel `K` and its derivative identity; the positivity of `dK/dx` (3×3×5 numeric sweep); normalization of `Sigma_Pe` to 1; the `Ic/Is` antiderivative regressions; the closed-form `Delta`; the endpoint formulas `Delta_0` (via `Pe->0` limit and via `integrate(K)`) and `Delta_inf` (via `K(1)` and via `Pe->oo` limit); that the closed `Delta` equals the kernel integral (SymPy: numeric sweep; Mathematica: a fully independent symbolic `Integrate[kernel*sigmaPe, ...]` compared to the `Ic/Is` combination form); the bracket-gap closed form and its positivity; the `Delta_0 <= Delta(Pe) <= Delta_inf` monotonicity sweep; the `F`-sign IVT bracket existence sweep; the weak-coupling constant and nonvanishing first-order series coefficient; and the IFT branch slope `dPe_*/dXi|_0 = Delta_0`. The bottom line tested is exactly the paper's fixed point + bracket + endpoint forms.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Endpoint `Delta_0 = eta(cosh a−1)/(a²·W)` (eq. Delta0) | `Delta0_expected` (py:84) vs `limit` and vs `integrate(K)` (py:86–87); Mma:91–97 | match |
| Endpoint `Delta_inf = (cosh a + (eta/a)sinh a − 1)/W` (eq. Deltainfty) | `Delta_inf_expected` (py:90) vs `K(1)` (py:92) and vs `Pe->oo` limit (py:201); Mma:100–105, 169–174 | match |
| Fixed point `Pe_* = Xi·Delta` (eq. Pe-fixed) | `Pe_lo/Pe_hi` (py:121–124); `F`-sign IVT existence sweep (py:169–196); Mma:114–117, 149–167 | match |
| Bracket `Xi·Delta_0 <= Pe_* <= Xi·Delta_inf` (eq. bracket) | bracket-gap closed form + positivity (py:127–140); monotonicity sweep (py:142–167); Mma:119–147 | match |
| Notes: kernel `K` + `dK/dx > 0` | Kprime identity + positivity sweep (py:42–61); Mma:38–60 | match |
| Notes: `Sigma_Pe` normalization | py:65; Mma:64 | match |
| Notes: closed `Delta` = `integral(K·Sigma_Pe)` | py:103–118 (numeric); Mma:76–84 (independent symbolic integral) | match |
| Notes: weak-coupling slope `= Delta_0` | py:219–227 (IFT); Mma:187–203 | match |

`paper_alignment: aligned` — every boxed paper equation and every numbered notes result has a corresponding non-tautological script check, with the engines computing the load-bearing `Delta`-equals-integral step by genuinely different routes.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 45–49 | `expect_zero(Kprime − closed form)` | notes kernel deriv | yes |
| A2 | sympy | 52–61 | numeric positivity sweep `> 0` | `dK/dx > 0` | yes |
| A3 | sympy | 65 | `integrate(Sigma)−1 == 0` | `Sigma_Pe` normalization | yes |
| A4 | sympy | 70–71 | `diff(Fc/Fs) − integrand == 0` | `Ic/Is` antiderivatives | yes |
| A5 | sympy | 86 | `Delta0 − Delta0_expected == 0` | endpoint `Delta_0` | yes |
| A6 | sympy | 87 | `Delta0 − integrate(K) == 0` | `Delta_0 = ∫K` | yes |
| A7 | sympy | 92 | `Delta_inf − expected == 0` | endpoint `Delta_inf` | yes |
| A8 | sympy | 109–118 | numeric `∫(K·Sigma) == Delta` sweep | closed `Delta` vs kernel integral | yes |
| A9 | sympy | 127–140 | bracket-gap form + positivity | bracket non-emptiness | yes |
| A10 | sympy | 142–167 | monotonicity sweep `Delta_0 ≤ Delta ≤ Delta_inf` | bracket | yes |
| A11 | sympy | 169–196 | `F`-sign IVT existence sweep | fixed-point root in bracket | yes |
| A12 | sympy | 199–201 | `Delta(Pe->oo) − Delta_inf == 0` | `Delta_inf` as limit | yes |
| A13 | sympy | 204–214 | weak-coupling const + nonvanishing 1st coeff | branch law | yes |
| A14 | sympy | 219–227 | `dPe_*/dXi|_0 − Delta_0 == 0` | weak-coupling slope | yes |
| B1 | mathematica | 46–49 | `expectZero(Kprime identity)` | kernel deriv | yes |
| B2 | mathematica | 64 | `expectZero(Sigma norm)` | normalization | yes |
| B3 | mathematica | 68–69 | `expectZero(Ic/Is antideriv)` | `Ic/Is` | yes |
| B4 | mathematica | 76–84 | independent `Integrate[K·Sigma]` − combination `== 0` | closed `Delta` vs kernel integral | yes |
| B5 | mathematica | 93–97 | `Delta0` formula + integral identity | endpoint `Delta_0` | yes |
| B6 | mathematica | 105 | `Delta_inf` substitution | endpoint `Delta_inf` | yes |
| B7 | mathematica | 124–147 | bracket gap + monotonicity sweeps | bracket | yes |
| B8 | mathematica | 150–167 | `F`-sign IVT sweep | fixed-point bracket | yes |
| B9 | mathematica | 174 | `Delta(Pe->oo) − Delta_inf == 0` | `Delta_inf` limit | yes |
| B10 | mathematica | 178–185, 203 | weak-coupling const/coeff + branch slope | branch law | yes |

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage058_coupled_support_source_sympy_audit.txt:3` and `:36`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.txt:3`
- (related stale self-labels in source) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py:2` and `:229`

**What's wrong:**
Both committed transcripts are older than their scripts (sympy `.txt` mtime 2026-05-26 11:00 vs `.py` 2026-06-03 15:59; mathematica `.txt` 2026-05-26 11:01 vs `.wl` 2026-06-03 15:59), and their banners predate the in-script banner fix. The sympy transcript line 3 reads `STAGE 41 — COUPLED SUPPORT/SOURCE OPERATOR` and line 36 reads `Stage 41 audit passed.`, while the current script banner (py:33) emits `STAGE 58` and the final print (py:229) emits `Stage 41 audit passed.` — note the script's own docstring (py:2, "Stage 41 SymPy audit") and closing line (py:229, "Stage 41 audit passed.") are themselves stale `Stage 41` self-labels. The mathematica transcript line 3 reads `STAGE 041` while the current `.wl` banner (wl:32) emits `STAGE 058`. All numeric/symbolic results in both transcripts otherwise show `= 0` / `PASS` and the closed forms agree with what the current scripts would produce; only the stage-number labels and the (already-correct in source) banner strings differ.

**Why this matters:**
The committed transcript no longer faithfully reflects the current script text, so a reader trusting the saved output sees a "Stage 41/041" provenance label on a Stage 058 artifact. Purely cosmetic for the math, but it is exactly the script/output-band label drift the project tracks.

**Required change:**
Re-run both scripts to refresh the transcripts so the banners read `STAGE 58`/`STAGE 058`. Separately, fix the stale self-labels in the sympy source: `scripts/...stage058...py:2` "Stage 41 SymPy audit" → "Stage 058 ..." and `:229` `print("\nStage 41 audit passed.")` → "Stage 058 audit passed." (The sympy banner at py:33 says "STAGE 58" — preferably "STAGE 058" for consistency with the Mma banner, but that is cosmetic.) Per project policy these self-label fixes are low-severity `stale_output`; the orchestrator resolves label scope.

**Verification:**
After a fresh run, sympy `.txt:3` reads `STAGE 058 ...` (or `STAGE 58`) and the closing line reads `Stage 058 audit passed.`; mathematica `.txt:3` reads `STAGE 058 ...`. Both transcripts' mtimes newer than the scripts.

## Independent-derivation check (Mathematica)

The `.wl` is **not** a transliteration. The decisive `Delta`-equals-kernel-integral step is computed by genuinely different routes: SymPy declares the closed `Delta` from the `Ic/Is` combination and validates it against the kernel integral *numerically* at four concrete `(alpha, eta, Pe)` points (py:103–118), explicitly noting the symbolic integral "does not terminate in sympy" (py:97–98). Mathematica instead performs the *symbolic* integral `delta = Integrate[kernel*sigmaPe, {x,0,1}, ...]` (wl:76–79) and asserts it equals the independently written combination form `deltaCombination` (wl:80–84) — a check SymPy cannot do at all. The kernel itself is re-derived (`FullSimplify` of the Green-difference, wl:38–41) rather than copied; the surface forms even differ between transcripts (sympy `K` line 5 vs Mma `K` line 5 are equivalent but not identical groupings). The two engines therefore exercise the same physics through different machinery.

## Engine cross-check

Final outputs agree. Endpoint forms are identical up to surface algebra: `Delta_0 = eta(cosh α−1)/(α²(α sinh α + η cosh α))` in both (sympy `.txt:16`, mma `.txt:22`); `Delta_inf` agrees (sympy `.txt:19` `(α(cosh α−1)+η sinh α)/(α·W)` ≡ mma `.txt:27` `(−1+cosh α+(η sinh α)/α)/W`). The weak-coupling first-order coefficients match term-for-term after clearing the factor-of-2 denominator convention: sympy `.txt:32` numerator ×2 = `2α²sinh + αη cosh + αη − 4α cosh + 4α − 2η sinh`, identical to mma `.txt:65` numerator `4α + αη − 4α cosh + αη cosh + 2α²sinh − 2η sinh`. The `Power::infy`/`N::meprec` warnings in the mma transcript (lines 37–58) arise from evaluating `delta` numerically near the removable `Pe = alpha` singularity during the IVT/monotonicity sweeps; they are benign — every gated check still reaches `PASS`. No `engine_disagreement`.

## Verdict justification

`verdict: findings` with a single low-severity `stale_output` item. I attacked every assertion: the endpoint checks compare a `limit`/substitution against an *independently written* closed form (not a definition echoed back), the `Delta_0 = integrate(K)` and the Mma symbolic `Delta = Integrate[K·Sigma]` checks are strong non-tautological cross-derivations, the positivity/monotonicity/IVT sweeps cover a real 3×3×(4–5) grid spanning small/moderate/large `alpha, eta, Pe`, and the removable `Pe = alpha` singularity is correctly skipped in the SymPy sweeps. The symbol assumptions (`alpha, eta, Pe, Xi > 0`, `0 ≤ x ≤ 1`) match both the notes' physical setup (`alpha>0, eta>0`, constructive branch `Pe≥0`) and are not strong enough to make a wrong candidate pass. No hardcoded numeric result exists (the stage is purely symbolic), no tautology, no symbol-domain error, no missing branch, no paper misalignment — the boxed fixed point, bracket, and both endpoint functions are each exercised by name. The only defect is that the committed transcripts predate the current scripts (and the sympy source still carries "Stage 41" self-labels), which is the tracked script/output-band label drift, not a math problem.

## Value Reconciliation (pass-2 augmentation)

This stage emits **no numeric deliverable constants** — every result is symbolic in `alpha, eta, Pe, Xi`. The deliverable-level reconciliation is therefore over the closed-form symbolic results, each of which must appear in the `.tex` card and/or `.md` notes.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `K_(α,η)(x) = (cosh αx + (η/α)sinh αx − cosh α(1−x))/W` | py:39–41, wl:38–41; sympy.txt:5, mma.txt:5 | notes:171–174 | MATCH |
| `dK/dx = (α sinh αx + η cosh αx + α sinh α(1−x))/W` | py:48, wl:48; sympy.txt:6 | notes:209–211 | MATCH |
| `Sigma_Pe(x) = Pe·exp(Pe x)/(exp(Pe)−1)` | py:63, wl:62; sympy.txt:9 | notes:167; tex (n/a) | MATCH (notes) |
| `Ic(Pe,α) = [exp(Pe)(Pe cosh α − α sinh α) − Pe]/(Pe²−α²)` | py:73; sympy.txt:13 | notes:186–188 | MATCH |
| `Is(Pe,α) = [exp(Pe)(Pe sinh α − α cosh α) + α]/(Pe²−α²)` | py:74; sympy.txt:14 | notes:190–192 | MATCH |
| closed `Delta(Pe;κ,η)` (Pe/(e^Pe−1)·[(1−cosh α)Ic+(η/α+sinh α)Is]/W) | py:78–79, wl:80–84; sympy.txt:15 | notes:196–199 | MATCH |
| `Delta_0 = η(cosh α−1)/(α²·W)` | py:84; sympy.txt:16 | tex:27–30 (eq Delta0); notes:231–235 | MATCH |
| `Delta_inf = (cosh α + (η/α)sinh α − 1)/W` | py:90; sympy.txt:19 | tex:31–36 (eq Deltainfty); notes:243–247 | MATCH |
| fixed point `Pe_* = Xi·Delta(Pe_*;κ,η)` | py:219, wl:188; (structural) | tex:14–18 (eq Pe-fixed); notes:261 | MATCH |
| bracket `Xi·Delta_0 ≤ Pe_* ≤ Xi·Delta_inf` | py:121–124; sympy.txt:22–23 | tex:20–24 (eq bracket); notes:269 | MATCH |
| bracket gap `((α²−η)(cosh α−1)+αη sinh α)/(α²W)` | py:128–131; sympy.txt:24 | (internal lemma; not a stated deliverable) | INTERNAL |
| weak-coupling slope `dPe_*/dXi|_0 = Delta_0` | py:226; sympy.txt:34 | notes:297 (`Pe_* = Xi·Delta_0 + O(Xi²)`) | MATCH |

INTERNAL (scaffolding, no prose expected): `W = α sinh α + η cosh α` (auxiliary), `Kprime` numerator positivity-sweep values, `Fc/Fs` antiderivatives, small-`Pe` series expansion, first-order series coefficient, all PASS flags / sweep residuals, the `bracket gap` closed form (an internal lemma proving `Delta_inf ≥ Delta_0`, not a card/notes deliverable).

reconciliation: complete; 11 values checked, 0 misaligned

## Self-test notes

I checked: (1) variable-independence — the only differentiations are `diff(K, x)`, `diff(Fc/Fs, x)`, and `diff(F, Pe/Xi)`, and in every case the differentiated expression genuinely depends on the variable (`K`, `Fc`, `Fs` all contain `x`; `F = Pe − Xi·Delta` depends on both `Pe` and `Xi`), so no derivative is identically zero. (2) Removable singularity — `Delta` has a 0/0 at `Pe = alpha`; both engines either skip that point in numeric sweeps (SymPy py:154, 182) or take symbolic limits, so no spurious infinities drive a false PASS (the Mma `Power::infy` warnings are benign and the gated checks still pass). (3) Trivial-case — `Delta_0` (Pe→0 limit) and `Delta_inf` (Pe→∞ limit and `K(1)`) each reduce to the independently-written closed forms, confirming the endpoint anchors are real cross-checks rather than echoes. No directive with a script-side math fix is warranted; the sole finding is the tracked low-severity stale-output/label-drift, which I report but do not block on.
</content>
</invoke>
