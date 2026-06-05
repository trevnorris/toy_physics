---
unit_id: 009
batch: I.1
auditor_model: claude-opus-4-8
audit_date: 2026-06-04T00:00:00Z
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
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 009 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_009.tex`
- notes: `(none)` — no file matching `notes/stages/moving_throat_pde_stage009_*.md` exists (verified by directory listing and recursive find)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 40, the stage-009 summary line)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage009_projected_maxwell_near_throat_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage009_projected_maxwell_near_throat_mathematica_audit.txt`

## What the paper claims

Stage 009 ("Near-throat projection scale") distinguishes two ways of locally observing the projected-Maxwell data near a throat. For a **symmetric interior kernel** the first finite-width correction to a smooth profile is even and begins at `O(σ²)`. For a **one-sided mouth kernel** anchored at the throat opening, `W_ℓ(s) = (1/ℓ) w(s/ℓ)` with `∫₀^∞ w(u) du = 1` (eq:stage009-mouth-kernel), the average expands as
`⟨X⟩_ℓ = X(0) + ℓ μ₁ X'(0) + O(ℓ²)`, with `μ₁ = ∫₀^∞ u w(u) du` (eq:stage009-mouth-expansion). The `\stagefield{Output}` reads verbatim: *"Stage 009 explains why mouth-local projected EM data enter at first order in the mouth scale, while symmetric interior thickness effects enter later."* Deliverables: (D1) symmetric-kernel correction is even/`O(σ²)`; (D2) mouth-kernel correction is first-order `O(ℓ)`, captured by the general first-moment formula `⟨X⟩_ℓ = X(0)+ℓμ₁X'(0)+O(ℓ²)` for a generic normalized one-sided kernel; (D3, from appendix row) the "near-throat projected scale factors and the leakage/source normalization ledger" (the `μ_eff`, `ξ_eff` effective parameters and the surviving mixed-sector derivative term). The card body equations carry the qualitative O(ℓ)-vs-O(σ²) message as the bottom line; the appendix tags the stage `\StatusExactClosure{}`.

## What the script claims to verify

The SymPy script extends the projection-first Maxwell law onto a finite/half-line throat domain and verifies, for a **specific exponential mouth kernel** `W_ℓ(w)=e^{-w/ℓ}/ℓ` and a **specific even Gaussian interior kernel**: (i) the exact integration-by-parts recombination on the half line so the apparent `1/ℓ` boundary singularity cancels (the `-q0/ℓ` boundary piece cancels the `+q0/ℓ` from `-⟨W'Q⟩`); (ii) the mouth average `⟨Q⟩_ℓ = q0+ℓq1+ℓ²q2+ℓ³q3+ℓ⁴q4` and derivative `⟨∂_wQ⟩_ℓ = q1+ℓq2+ℓ²q3+ℓ³q4`, giving a leading O(ℓ) correction; (iii) the symmetric even-kernel moment expansion `⟨Q⟩_σ = q0+m2σ²q2/2+m4σ⁴q4/24` (odd moments vanish), giving leading O(σ²); (iv) the first-order zero-mode effective-parameter series for `μ_eff=μ0⟨S⟩/⟨Z⟩` and `ξ_eff=ξ⟨Z⟩/⟨H⟩` in both regimes, plus perturbative profile-alignment corollaries (H=Z ⟹ ξ_eff=ξ; S=CZ ⟹ μ_eff correction vanishes); (v) Gaussian-localizer fingerprints `⟨Z⟩_sym = 1 − σ²/λ² + 3σ⁴/(2λ⁴)` and `⟨Z⟩_mouth = 1 − 2ℓ²/λ² + 12ℓ⁴/λ⁴ − 120ℓ⁶/λ⁶`. The Mathematica `.wl` independently re-verifies (i)–(v)'s core identities (M1–M5). The script's docstring frames the whole thing as making the near-throat structure auditable, especially the surviving mixed flux derivative `∂_w(Z F^{wν})`.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| D1: symmetric interior correction even, O(σ²) | sympy 77–87 / `.wl` M2,M4; Gaussian sym 229–232 / M5a | match |
| D2 (qualitative): mouth correction is first-order O(ℓ) | sympy 110–117 (`⟨Q⟩_ℓ=q0+ℓq1+…`) / `.wl` M1; mouth Gaussian 236–252 / M5b | match |
| D2 (general kernel): `⟨X⟩_ℓ=X(0)+ℓμ₁X'(0)+O(ℓ²)`, `μ₁=∫₀^∞ u w(u)du` for a **generic** one-sided kernel | script verifies ONLY the exponential kernel (`μ₁=∫₀^∞ u e^{-u}du = 1`); no free first-moment symbol `μ₁` is ever introduced or carried | partial |
| D3: near-throat scale factors / leakage / source-norm ledger | `μ_eff=μ0⟨S⟩/⟨Z⟩`, `ξ_eff=ξ⟨Z⟩/⟨H⟩` (sympy 148–214 / `.wl` M3,M4); mixed-sector term `⟨∂_w(ZF^{wν})⟩` (prose §1,§6) | match |
| (script extra) perturbative H=Z / S=CZ corollaries (sympy 190–214) | not mentioned in card/appendix; script-side elaboration of D3 | extra (benign — consistent with D3) |

Dominant pattern: the qualitative bottom-line `\stagefield{Output}` (O(ℓ) mouth vs O(σ²) interior) and the scale-factor ledger are faithfully exercised; the only gap is that the paper card displays a **generic-kernel** first-moment formula (with free `μ₁`) that the script only instantiates for the exponential representative. Hence `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 86 | `assert_zero(Gaussian even-kernel Q moments)` — integral vs literal moment formula | D1 | yes |
| A2 | sympy | 87 | `assert_zero(Gaussian even-kernel derivative moments)` | D1 | yes |
| A3 | sympy | 111 | `assert_zero(half-line Q expansion)` vs `q0+ℓq1+ℓ²q2+ℓ³q3+ℓ⁴q4` | D2 (exp kernel) | yes |
| A4 | sympy | 116 | `assert_zero(half-line boundary recombination)` (`1/ℓ` cancellation) | D2 / D3 | yes |
| A5 | sympy | 117 | `assert_zero(half-line derivative expansion)` vs `q1+ℓq2+ℓ²q3+ℓ³q4` | D2 (exp kernel) | yes |
| A6 | sympy | 151,154 | `assert_zero(symmetric μ_eff / ξ_eff series)` to O(σ²) | D3 | yes |
| A7 | sympy | 161,164 | `assert_zero(mouth μ_eff / ξ_eff series)` to O(ℓ) | D3 | yes |
| A8 | sympy | 203,204 | `assert_zero(H=Z / S=CZ eps=0 cancellation)` | D3 (corollary) | yes |
| A9 | sympy | 211,213 | `assert_zero(ξ_eff / μ_eff first-order in eps)` | D3 (corollary) | yes |
| A10 | sympy | 231 | `assert_zero(symmetric Gaussian asymptotic literal)` vs `1−σ²/λ²+3σ⁴/(2λ⁴)` | D1 | yes |
| A11 | sympy | 245 | `assert_zero(mouth Gaussian integral equals erfc closed form)` (guard) | D2 (exp kernel) | yes |
| A12 | sympy | 251 | `assert_zero(mouth Gaussian asymptotic from erfc)` vs `1−2ℓ²/λ²+12ℓ⁴/λ⁴−120ℓ⁶/λ⁶` | D2 (exp kernel) | yes |
| A13 | sympy | 252 | `assert_zero(mouth Gaussian asymptotic from Taylor integration)` — independent cross-check | D2 (exp kernel) | yes |
| A14 | math | 26 | `assertZero(M1a half-line IBP recombination)` | D2 / D3 | yes |
| A15 | math | 27,29 | `assertZero(M1b/M1c half-line dQ/Q expansion)` | D2 (exp kernel) | yes |
| A16 | math | 43,45 | `assertZero(M2a/M2b Gaussian even-kernel)` | D1 | yes |
| A17 | math | 60,62 | `assertZero(M3a/M3b mouth μ_eff/ξ_eff to O(ℓ))` | D3 | yes |
| A18 | math | 77,79 | `assertZero(M4a/M4b symmetric μ_eff/ξ_eff to O(σ²))` | D3 | yes |
| A19 | math | 91,98 | `assertZero(M5a/M5b Gaussian asymptotics)` | D1 / D2 | yes |

Every assertion traces to a paper deliverable and is non-tautological: each `assert_zero`/`assertZero` independently computes one side (symbolic integration, `series`, or limit) and compares against a hand-written literal closed form or against a second independent route. No assertion compares an expression to itself by construction.

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** script_missing_paper_claim
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_009.tex:23-29` (eq:stage009-mouth-expansion)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.py:107-117`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage009_projected_maxwell_near_throat_mathematica_audit.wl:20-30`

**What's wrong:**
The paper card states the mouth expansion for a **generic** normalized one-sided kernel and explicitly displays the first-moment functional:

> `\langle X\rangle_\ell = X(0)+\ell\mu_1X'(0)+O(\ell^2), \qquad \mu_1=\int_0^\infty u\,w(u)\,du`  (stage_009.tex:25-28)

with the kernel left general: `W_\ell(s)=\frac1\ell w(s/\ell), \int_0^\infty w(u)\,du=1` (stage_009.tex:18-21).

Both scripts only ever instantiate the **specific exponential kernel** `W_ℓ(w)=e^{-w/ℓ}/ℓ`:

> sympy:108 `Wexp = sp.exp(-w/ell)/ell` → asserts `⟨Q⟩_ℓ = q0+ℓq1+ℓ²q2+ℓ³q3+ℓ⁴q4` (line 111)
> `.wl`:21 `Wel = Exp[-w/ell]/ell` → M1c `avgQ - (q0 + ell q1 + …)` (line 30)

For the exponential kernel `w(u)=e^{-u}`, so `μ₁ = ∫₀^∞ u e^{-u} du = 1`, which is why the script's first-order coefficient is exactly `q1 = X'(0)` with an implicit factor `1`. The script's specific result is therefore *consistent with* the paper's general formula — but it never introduces a free `μ₁` symbol, and it never demonstrates the `μ₁` first-moment structure for any kernel other than the exponential. The general claim that the leading mouth coefficient is `μ₁ X'(0)` for an *arbitrary* normalized one-sided kernel is a stated paper deliverable (the displayed body equation) that no script-side check exercises. The qualitative bottom line (`\stagefield{Output}`: mouth data enter at first order in ℓ) IS verified; only the kernel-generality of the displayed `μ₁` formula is not.

**Why this matters:**
Low impact: the dominant `\stagefield{Output}` claim and the appendix-row ledger are fully verified, and the exponential-kernel result is a valid instance of the general formula. But strictly, the card displays a general-kernel identity (with the named functional `μ₁`) that the verification specializes to one kernel where `μ₁=1` collapses out of view. A reader checking "is eq:stage009-mouth-expansion verified?" will find only the special case. This is a scope/coverage mismatch between the displayed paper equation and the script, routed to the user because the direction of resolution is a paper-vs-script choice.

**Required change:**
Do NOT auto-edit. This routes to the user via the directive's `## Resolve before fix_loop` block. Two clean directions: (a) treat the displayed general `μ₁` formula as the load-bearing claim and add a script check with a generic kernel (e.g. a symbolic normalized kernel, or a second representative kernel such as a one-sided Gaussian/box, demonstrating the leading coefficient is `μ₁ X'(0)` with `μ₁` the explicit first moment); or (b) treat the card's `μ₁` line as an illustrative general statement whose verified instance is the exponential kernel (`μ₁=1`), and add one card-side sentence noting the exponential representative is what the audit script exercises. Either way the qualitative O(ℓ)-vs-O(σ²) claim is already verified, so this is non-blocking.

**Verification:**
If direction (a): a new generic-kernel assertion appears (sympy + `.wl`) producing leading coefficient `μ₁ X'(0)` with `μ₁ = ∫₀^∞ u w(u) du` as a free symbol, exiting 0. If direction (b): card line near stage_009.tex:23-29 gains an explicit note that the exponential kernel `e^{-w/ℓ}/ℓ` (`μ₁=1`) is the audited representative; no script change.

## Independent-derivation check (Mathematica)

The `.wl` is an **independent** re-derivation, not a transliteration. Evidence:
- The mouth asymptotic is obtained natively: `.wl`:95-97 uses `IMou /. ell -> 1/r` then `Series[IMouR, {r, Infinity, 7}]` and FullSimplify, whereas the `.py` first rewrites the integral via `erfc`, substitutes `ell→1/r`, takes `series(...,r,oo,8)`, AND adds two extra guards (erfc closed-form equality at sympy:245 and a Taylor-integration cross-check at sympy:252) that the `.wl` does not replicate.
- The `.wl` omits the `.py`'s prose-only sections (§1, §6) and the entire perturbative profile-alignment block (sympy:190-214); it does not echo the `.py`'s `assert_zero`/`sp.factor(sp.together(...))` choreography but uses `FullSimplify[..., Assumptions->$Assumptions]` with a fresh `$Assumptions` reset before each M-block.
- The shared literal RHS targets (`1 - 2ℓ²/λ² + 12ℓ⁴/λ⁴ - 120ℓ⁶/λ⁶`, `q0+ℓq1+ℓ²q2+...`) are the *paper's claimed answers* that both engines must legitimately reproduce — convergent verification, not algebra-echoing. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree. The SymPy output (`...sympy_audit.txt`) prints every section's closed form and ends `STATUS: PASS` with no `AssertionError`. The Mathematica output (`...mathematica_audit.txt`) reports `OK ... residual = 0` for all 11 M-checks (M1a–M5b) and `STATUS: PASS`. The overlapping deliverables match exactly: mouth `⟨Q⟩_ℓ` and `⟨∂_wQ⟩_ℓ` expansions (py §3 ↔ M1), even-kernel moments `q0+q2/2+q4/8` (py §2 ↔ M2), `μ_eff`/`ξ_eff` first-order series (py §4 ↔ M3,M4), and Gaussian fingerprints `1−σ²/λ²+3σ⁴/(2λ⁴)` and `1−2ℓ²/λ²+12ℓ⁴/λ⁴−120ℓ⁶/λ⁶` (py §5 ↔ M5). No sign/factor disagreement.

## Verdict justification

The scripts hold up under attack. I re-derived by hand: the exponential-kernel moment integrals `∫₀^∞ (e^{-w/ℓ}/ℓ) w^n dw = ℓⁿ n!` (giving the asserted `q0+ℓq1+ℓ²q2+...` and the `1/ℓ` boundary cancellation), the Gaussian even-moments (1,0,1,0,3,0 → `q0+q2/2+q4/8`), the symmetric Gaussian convolution `λ/√(λ²+2σ²)` with series `1−σ²/λ²+3σ⁴/(2λ⁴)`, and the mouth Gaussian asymptotic `1−2ℓ²/λ²+12ℓ⁴/λ⁴−120ℓ⁶/λ⁶` via the Taylor-integration route (`w^{2k}` term → `ℓ^{2k}(2k)!`). All match the script literals. Assertions are non-tautological (each side computed by an independent route) and the mouth Gaussian even has a triple cross-check (erfc closed form, asymptotic-of-erfc, and Taylor-integration). The single finding is a low-severity `paper_misalignment`: the card displays a *generic-kernel* first-moment formula (`⟨X⟩_ℓ=X(0)+ℓμ₁X'(0)`, `μ₁=∫₀^∞ u w(u)du`) that the verification only instantiates for the exponential kernel (`μ₁=1`). The qualitative `\stagefield{Output}` and the scale-factor/leakage ledger are fully verified, so the verdict is `findings` (not stop-cold), `paper_alignment: partial`, routed to the user for direction.

## Value Reconciliation (pass-2 augmentation)

Per `RECONCILIATION_AUGMENTATION.md`, every result/deliverable value the scripts emit is accounted for below. Outputs are FRESH (both scripts mtime 2026-05-21 11:11; sympy output 11:26, mathematica output 11:51 — both newer than their scripts), so the saved `.txt` files are authoritative for what the scripts produce.

Note on doc carriers: there are **no** `notes/stages/` files for stage 009 (the natural `.md` carrier is absent), so the only prose carrier is the terse `.tex` card. Per the augmentation guards, the card is allowed to omit intermediate quantities; only **stated deliverables absent from BOTH .tex and .md** count as MISSING. Most of the script's emitted closed forms are intermediate working steps that a terse card legitimately omits and are therefore INTERNAL.

| value | source (py / wl + output line) | .tex/.md location | status |
|---|---|---|---|
| Symmetric interior correction even, leading `O(σ²)`: `⟨Q⟩_σ = q0 + m2σ²q2/2 + m4σ⁴q4/24` | py:77 / out §2 L28; wl M2 / out L4 | `.tex:13-15` ("first finite-width correction … is even and begins at \(O(σ²)\)") | MATCH |
| Mouth average `⟨Q⟩_ℓ = q0 + ℓq1 + ℓ²q2 + ℓ³q3 + ℓ⁴q4` → leading `O(ℓ)` | py:111 / out §3 L47; wl M1c / out L3 | `.tex:23-29` (`⟨X⟩_ℓ = X(0)+ℓμ₁X'(0)+O(ℓ²)`) | MATCH (qualitative O(ℓ); general `μ₁` form is F1 partial) |
| General mouth first-moment functional `μ₁ = ∫₀^∞ u w(u) du` (generic kernel) | (not emitted as a free symbol; specialized to exp-kernel `μ₁=1`) | `.tex:28` (`\mu_1=\int_0^\infty u w(u) du`) | MISSING-DELIVERABLE → folded into F1 (paper_misalignment, low) |
| Effective coupling `μ_eff = μ0 ⟨S⟩/⟨Z⟩` and gauge `ξ_eff = ξ ⟨Z⟩/⟨H⟩` (the "scale factors / normalization ledger") | py:172 / out §4 L76; wl M3/M4 | appendix `part01.tex:40` ("near-throat projected scale factors and the leakage/source normalization ledger") | MATCH |
| Surviving mixed-sector ("leakage") term `⟨∂_w(Z F^{wν})⟩` | py prose §1/§6 / out L18,L127; wl (implicit via M1 IBP) | appendix `part01.tex:40` ("leakage … ledger"); `.tex` §"Interior versus mouth" theme | MATCH |
| Symmetric Gaussian fingerprint `⟨Z⟩ = 1 − σ²/λ² + 3σ⁴/(2λ⁴)` | py:232 / out §5 L108; wl M5a / out L10 | (illustrative example; not a card deliverable) | INTERNAL (example) |
| Mouth Gaussian fingerprint `⟨Z⟩ = 1 − 2ℓ²/λ² + 12ℓ⁴/λ⁴ − 120ℓ⁶/λ⁶` | py:251 / out §5 L116; wl M5b / out L11 | (illustrative example; not a card deliverable) | INTERNAL (example) |
| `μ_eff`/`ξ_eff` first-order series (sym O(σ²) & mouth O(ℓ), e.g. `μ0 s0/z0 + ℓ(μ0 s1/z0 − μ0 s0 z1/z0²)`) | py:152,162 / out §4 L79-84; wl M3a-M4b | (intermediate expansion; card states only the qualitative order) | INTERNAL (terse-card omission allowed) |
| Profile-alignment corollaries: `ξ_eff\|_{H=Z}=ξ`, `μ_eff\|_{S=CZ}=C μ0` | py:203,204 / out §4 L95-96 | (script-side corollary; not in card/appendix) | INTERNAL (script extra, benign) |
| Even-kernel Gaussian moment literal `q0 + q2/2 + q4/8` | py:86 / out (asserted); wl M2a | (intermediate; supports D1) | INTERNAL |

INTERNAL scaffolding (no finding): pass/fail flags `STATUS: PASS`; `assert_zero`/`assertZero` residual-zero outputs; the boundary-split intermediates `[WQ]_0^∞ = −q0/ℓ` and `−⟨W'Q⟩ = …+q0/ℓ` (these exist only to drive the `1/ℓ`-cancellation assertion); the erfc closed form `√π λ erfc(λ/2ℓ) e^{λ²/4ℓ²}/(2ℓ)` (intermediate guard, py:244); the `r=1/ℓ` asymptotic device; the symmetric convolution `λ/√(λ²+2σ²)` (intermediate, py:229).

reconciliation: 1 deliverable misaligned. 10 deliverable-level values checked; 9 reconcile (MATCH/illustrative), 1 MISSING-DELIVERABLE (the generic-kernel `μ₁` first-moment formula) — folded into finding F1 (low-severity paper_misalignment). The MISSING-DELIVERABLE does not introduce a *second* finding beyond F1; it is the same paper↔script scope gap already raised, so `findings_count` stays 1.

## Self-test notes

I checked: (1) Variable independence — every `sp.diff`/`D[...]` operates on `Q`/`Z` that genuinely depend on the differentiation variable `w` or `u`; no identically-zero-derivative trap. (2) Parity/symmetry — the symmetric kernel uses an even Gaussian over `(−∞,∞)` so odd moments correctly vanish (verified `⟨u⟩=⟨u³⟩=⟨u⁵⟩=0`), and the mouth integrals are over `[0,∞)` with no parity cancellation assumed. (3) Trivial-case substitution — set `q1=Q'(0)`: the mouth average's first-order coefficient is exactly `q1` (=`μ₁·q1` with `μ₁=1` for the exp kernel), confirming O(ℓ); the symmetric average's first nonzero correction is the `q2` term at O(σ²), confirming the even/later-order claim. (4) Mouth Gaussian coefficients re-derived two ways (Taylor-integration `w^{2k}→ℓ^{2k}(2k)!` and the script's erfc-asymptotic) both give `1−2ℓ²/λ²+12ℓ⁴/λ⁴−120ℓ⁶/λ⁶`. The only issue surfaced is the generic-`μ₁` coverage gap (F1), which is a paper↔script scope question, not a math error.
