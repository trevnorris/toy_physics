---
unit_id: 089
batch: III.5
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [notes/stages/moving_throat_pde_stage089_family1_minimal_isotropic_verdict.md]
  paper_appendix: present
---

# Audit unit 089 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_089.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage089_family1_minimal_isotropic_verdict.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 156)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_mathematica_audit.txt`

## What the paper claims

Stage 089 is a CHECKPOINT verdict: it evaluates the minimal isotropic passive/outgoing
quadrupole demand (`zeta_req = 1/3`, carried in from Stage 088 as `rho_alpha = 4/3`)
against the explicit Family-1 support/source branch. The card states `\stagefield{Inputs}`
as `zeta_req = 1/3` and the zero-bias Family-1 value `A_F1 ≈ 1.00005192880220`. The
derivation hinges on the 0/0 transport-map limit `Omega(Pe=0)=1`, which gives
`zeta_F1(0)=A_F1 ≈ 1.00005192880220` (eq `app-stage089-zero-bias-supply`). The boxed
success theorem is `zeta_req^min = 1/3 < A_F1` (eq `app-stage089-success`), and the boxed
`\stagefield{Output}` is the zero-bias success theorem `Pe_req = 0` (eq
`app-stage089-Pe-zero`). The notes add the supporting deliverables: the loading-ratio
window `rho_suff ≈ 3.46622291347846`, `rho_fail ≈ 3.46752913273870`,
`rho_max ≈ 3.46752922945601`; the support ceiling `zeta_max^(F1) ≈ 2.46752922945601`;
the margins `Delta_suff ≈ 2.13288958014513`, `Delta_fail ≈ 2.13419579940537`,
`Delta_max ≈ 2.13419589612268`; and the regime check `0 < zeta_req = 1/3 < 1`
(symmetric-lowest-twin window). Appendix row 156: "Minimal isotropic branch succeeds at
zero transport bias."

## What the script claims to verify

Both engines reconstruct the Family-1 transport map `zeta_F1(Pe) = A_F1 · Omega(Pe)^2`
from first ingredients: `A_F1 = (kappa_F1 + pi^2/4)/(kappa_F1 + y_F1^2)` with `kappa_F1 =
12321/5`, `y_F1` the root of `y tan(y) = 37` near `1.53` (SymPy `nsolve` seed 1.53;
Mathematica `FindRoot` seed 1.53), and the closed-form `Omega`. They assert the chain
closure `Omega(Pe→0)=1` and `zeta_F1(Pe→0)=A_F1` (l'Hopital / series of the genuine 0/0
form), the ceiling identity `zeta_max = A_F1·pi^2/4`, cross-check `rho_suff/rho_fail/rho_max`
against the Stage-082 literal quotes, assert the ordering `rho_min < rho_suff < rho_fail <
rho_max`, the regime `zeta_min < 1`, the zero-bias success condition `zeta_min < A_F1`, the
ceiling `zeta_min < zeta_max`, and finally `Pe_req = 0`. The Mathematica side independently
re-derives the upstream `Pe_suff/Pe_fail` via `FindRoot` (instead of SymPy's literal
carry-forwards), giving a genuinely independent second path.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `Omega(Pe=0)=1` (eq zero-bias-supply) | sympy L55 / wl L53 `expect[Omega(0)-1]` | match |
| `zeta_F1(0)=A_F1` (eq zero-bias-supply) | sympy L56 / wl L54 `expect[zeta_F1(0)-A_F1]` | match |
| `zeta_req^min = 1/3 < A_F1` (eq success, boxed) | sympy L119 `if not (zeta_min < A_F1)` / wl L105 `expectTrue` | match |
| `Pe_req = 0` (eq Pe-zero, boxed Output) | established by L55/56/119; restated by sympy L128-129 / wl L112-113 (tautological echo — see F2) | match (deliverable covered; echo cosmetic) |
| `A_F1 ≈ 1.00005192880220` (Inputs) | derived sympy L40 / wl L41, output 1.000051928802195 | match |
| `rho_alpha/rho_min = 4/3` (notes/088) | sympy L31 / wl L34 | match |
| `zeta_max^(F1) ≈ 2.46752922945601` (notes) | sympy L46/76 / wl L45/74 | match |
| `rho_suff/rho_fail/rho_max` (notes §1) | sympy L90-92 / wl L79-81 `expect_close` | match |
| `Delta_*` margins (notes §1) | sympy L94-97 / wl L83-87, printed | match |
| `0 < zeta_req < 1` regime (notes §2) | sympy L117 / wl L104 | match |

`paper_alignment: aligned` — every paper-side deliverable has a faithful script-side check.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 55 | `expect_zero[Omega(0)-1]` | Omega(Pe=0)=1 | yes (genuine 0/0 l'Hopital) |
| A2 | sympy | 56 | `expect_zero[zeta_F1(0)-A_F1]` | zeta_F1(0)=A_F1 | yes |
| A3 | sympy | 76 | `expect_zero[zeta_max - A_F1·pi^2/4]` | zeta_max ceiling | yes (limit derived independently L46) |
| A4 | sympy | 77 | `expect_zero[Q - (1+zeta)]` | (none — algebraic identity of Q at eps=0) | NO — tautological (F1) |
| A5 | sympy | 90-92 | `expect_close[rho_* vs literal]` | rho window | yes (A_F1/Omega/Q derived; vs upstream literal) |
| A6 | sympy | 115 | `rho_min<rho_suff<rho_fail<rho_max` | ordering | yes |
| A7 | sympy | 117 | `zeta_min < 1` | symmetric-twin regime | yes |
| A8 | sympy | 119 | `zeta_min < A_F1` | success theorem | yes (core can-fail check) |
| A9 | sympy | 121 | `zeta_min < zeta_max` | below ceiling | yes |
| A10 | sympy | 129 | `expect_zero[Pe_req]` with `Pe_req=0` | Pe_req=0 Output | NO — tautological echo (F2) |
| B1 | mathematica | 53 | `expectApprox[Omega(0),1]` | Omega(Pe=0)=1 | yes |
| B2 | mathematica | 54 | `expectApprox[zeta_F1(0),aF1]` | zeta_F1(0)=A_F1 | yes |
| B3 | mathematica | 74 | `expectApprox[zetaMax, aF1·pi^2/4]` | zeta_max ceiling | yes |
| B4 | mathematica | 79-81 | `expectApprox[rho_* vs literal]` | rho window | yes (Pe re-derived via FindRoot L65-66) |
| B5 | mathematica | 103 | `expectTrue[ordering]` | ordering | yes |
| B6 | mathematica | 104 | `expectTrue[zetaMin<1]` | regime | yes |
| B7 | mathematica | 105 | `expectTrue[zetaMin<aF1]` | success theorem | yes |
| B8 | mathematica | 106 | `expectTrue[zetaMin<zetaMax]` | below ceiling | yes |
| B9 | mathematica | 113 | `expectApprox[peReq,0]` with `peReq=0` | Pe_req=0 Output | NO — tautological echo (F2) |

## Findings

(No `paper_misalignment` findings. Two low-severity `tautological_check` findings.)

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py:77`

**What's wrong:**
Line 69 defines `Q = sp.simplify((1 + (1 - 2 * eps_blk) * zeta) / (1 - eps_blk * zeta))`
with `eps_blk = sp.Integer(0)` (line 67), so `Q` simplifies to exactly `1 + zeta`. Line 77
then asserts:
```
expect_zero("Stage-082 Q(zeta;0) = 1 + zeta", Q - (1 + zeta), tol=1e-30)
```
`Q - (1 + zeta)` is identically `0` by construction — the assertion cannot fail for any
physics. The script's own comment at lines 87-89 explicitly acknowledges this exact form is
tautological ("The previous `rho_X - (1 + zeta_X) == 0` form was tautological because
Q(zeta; eps=0) = 1 + zeta is the algebraic structure of Q…") and that is why the rho cross-
checks were upgraded to `expect_close` against literals (lines 90-92). The parallel
`Q - (1+zeta)` check at line 77 was left behind. The Mathematica side does not carry an
analogous check (it goes straight to the `expectApprox` rho checks at lines 79-81),
confirming this is a SymPy-only residual.

**Why this matters:**
On a checkpoint stage the bar is "assertions must be substantive." A check that is true by
construction adds no verification value and contradicts the script's own stated remediation
rationale two comments later. It does not affect any deliverable (the rho window is genuinely
checked by A5/B4), but it is dead scaffolding that should be removed for the checkpoint to be
clean.

**Required change:**
Delete line 77 entirely (`expect_zero("Stage-082 Q(zeta;0) = 1 + zeta", Q - (1 + zeta),
tol=1e-30)`). The non-tautological `expect_close` rho cross-checks at lines 90-92 already
exercise the `Q(zeta;0) = 1 + zeta` structure against external literals; no replacement
assertion is needed.

**Verification:**
After the fix, the sympy output no longer contains the line
`Stage-082 Q(zeta;0) = 1 + zeta = 0`; the script still exits 0 and the rho cross-checks
(lines now ~88-90) remain.

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py:128-129`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_mathematica_audit.wl:112-113`

**What's wrong:**
SymPy lines 128-129:
```
Pe_req = sp.Integer(0)
expect_zero("Pe_req (zero-bias bound from chain closure)", Pe_req)
```
Mathematica lines 112-113:
```
peReq = 0;
expectApprox["Pe_req (zero-bias bound from chain closure)", peReq, 0, 10^-30];
```
In both engines `Pe_req` is hardcoded to `0` and then checked to be `0`. The assertion is
`0 == 0` by construction and cannot fail regardless of the physics — it is a cosmetic echo
of the paper's boxed Output, not a verification of it.

Crucially, the deliverable `Pe_req = 0` IS genuinely established elsewhere by can-fail
checks: `Omega(Pe→0)=1` (A1/B1), `zeta_F1(0)=A_F1` (A2/B2), and `zeta_min < A_F1` (A8/B7).
Together these show the zero-bias supply `zeta_F1(0)=A_F1` already exceeds the demand
`zeta_req = 1/3`, so by the transport-map rule (notes line 98: "if zeta_req < A_F1, the
demand is already met at zero transport bias") `Pe_req = 0`. So the deliverable is covered;
the `Pe_req = sp.Integer(0)` line is a redundant tautological restatement.

**Why this matters:**
On a checkpoint the load-bearing Output (`Pe_req = 0`) must be asserted substantively. Here
it is — but the line that NAMES the Output is the one tautological line, while the genuine
support is spread across A1/A2/A8. A reader scanning for "where is Pe_req=0 verified?" lands
on a `0==0` check. The honest fix is to make the named assertion non-tautological by
deriving `Pe_req` from the success condition rather than hardcoding `0`.

**Required change:**
Replace the hardcoded `Pe_req = sp.Integer(0)` with a value derived from the established
zero-bias success condition so the named assertion can fail if the precondition breaks. For
example (SymPy):
```
# Pe_req = 0 iff the zero-bias supply already meets the demand; encode that
# implication so the named check fails if zeta_min >= zeta_F1(0).
Pe_req = sp.Piecewise((sp.Integer(0), zeta_min < zeta_F1_at_zero), (sp.nan, True))
expect_zero("Pe_req (zero-bias success => 0)", Pe_req)
```
and the Mathematica analogue:
```
peReq = If[zetaMin < zetaF1AtZero, 0, Indeterminate];
expectApprox["Pe_req (zero-bias success => 0)", peReq, 0, 10^-30];
```
This makes the named Output assertion contingent on the actual zero-bias success condition
(`zeta_min < zeta_F1(0)`) rather than a literal `0`. If Codex judges the Piecewise form
unsafe to apply mechanically, an acceptable minimal alternative is to drop the standalone
`Pe_req` echo and rely on A1/A2/A8 (which already establish the Output), leaving an inline
comment that points to those three lines as the substantive proof of `Pe_req = 0`.

**Verification:**
After the fix, the named `Pe_req` assertion is no longer `expect_zero(sp.Integer(0))` /
`expectApprox[0,0]` but is gated on `zeta_min < zeta_F1(0)`; both scripts still exit 0 and
the printed `Pe_req = 0` / Output line remains.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration. Key independence: the Mathematica script re-derives the
upstream `Pe_suff/Pe_fail` itself via `FindRoot[zetaF1[pe] == zetaTarget, {pe, 100/10000},
WorkingPrecision -> 40]` (lines 65-66), where SymPy keeps the literal carry-forwards
`Pe_suff_chi = 96.5285…`, `Pe_fail_chi = 11220.544…` (sympy lines 65-66) with an explicit
provenance comment. This is the robust root-finding route the checkpoint guidance calls for,
and it is a genuinely different computational path from SymPy's literals. `y_F1` is found by
`FindRoot` (wl L40) vs SymPy `nsolve` (py L39); both seed 1.53 on the principal branch of
`y tan(y) = 37` below the `pi/2` pole. The shared `omegaPe / zetaF1 / q` function definitions
mirror SymPy's `Omega / zeta_F1 / Q` — that is expected, since both must implement the same
physical map — but the Pe-value derivation differs, so this clears the second-engine bar.

## Engine cross-check

The two outputs agree to displayed precision:
- `A_F1`: sympy `1.000051928802195328659334`, wl `1.00005192880219532865933408…` — agree.
- `zeta_max`: both `2.46752922945601233329585…` — agree.
- `rho_suff/rho_fail/rho_max`: both `3.466222913478…`, `3.467529132738…`,
  `3.467529229456…` (all within `1e-12` of the Stage-082 quotes) — agree.
- `Delta_suff/fail/max/zeta/AF1`: both match (`2.132889580145…`, `2.134195799405…`,
  `2.134195896122…`, `2.134195896122…`, `0.666718595468…`). Note the internal consistency
  `Delta_max = Delta_zeta` (both `2.13419589612268`), since `rho = 1 + zeta` at `eps=0`.
- `Pe_req`: both `0`.
No `engine_disagreement`. Both transcripts end with all checks PASS / no `AssertionError`.

## Verdict justification

`findings`. The checkpoint's core deliverables — `Omega(Pe=0)=1`, `zeta_F1(0)=A_F1`, the
success theorem `zeta_req = 1/3 < A_F1`, the loading-ratio window, the regime check, and the
boxed Output `Pe_req = 0` — are all faithfully and (for the can-fail checks) substantively
verified, in two genuinely independent engines, and every emitted value reconciles exactly
with the card and notes. I attacked the `nsolve`-near-`tan`-pole pitfall (#10): the `y_F1`
seed 1.53 sits on the correct principal branch below `pi/2≈1.5708` and converges to the root
≈1.528 (confirmed by `A_F1` matching the paper to 14 digits); the fragile `Pe` re-derivation
is deliberately a literal carry-forward in SymPy (provenance-commented) and a robust
`FindRoot` in Mathematica — no fragile pole-adjacent `nsolve` is load-bearing. I confirmed
the chain-closure limits the first pass flagged are present AND asserted (not print-only) in
both engines. The only defects are two `tautological_check` echoes: line 77's
`Q - (1+zeta)` residual (true by construction; the script's own comment admits the form is
tautological) and the hardcoded `Pe_req = 0` checked against `0` in both engines. Both are
low-severity (the underlying deliverables are covered by genuine checks), but on a checkpoint
the named assertions should be substantive, so they are findings. No `paper_misalignment`,
no symbol-assumption error, no engine disagreement, no stale output, no stale self-label
(all "Stage 089" self-labels are correct; "Stage-62/082/075/074" are cross-references,
deferred per numbering policy).

## Self-test notes

I checked: (1) variable independence — no `sp.diff`/`D[]` in this stage, n/a. (2)
symmetry/parity — no unbounded integrals, n/a. (3) trivial-case pre-check — for F2's proposed
`Pe_req` derivation, substituting `zeta_min = 1/3 < zeta_F1(0) = A_F1 ≈ 1.00005` gives the
Piecewise true-branch `0`, so `expect_zero` passes (correct); if `zeta_min` ever exceeded
`zeta_F1(0)` the check would yield `nan` and fail (correct gating). (4) path specs — no
missing-script findings. (5) paper round-trip — F1 only deletes a tautology (introduces no
new constant); F2 keeps `Pe_req = 0` and only makes it contingent on the already-paper-stated
condition `zeta_req < A_F1`, so no new `paper_misalignment` is introduced. The `nsolve`/pole
and l'Hopital 0/0 concerns were checked and clear.

## Value Reconciliation (pass-2 augmentation)

Reconciliation based on script source + committed saved outputs (both fresh: sympy .txt
mtime 2026-05-27 10:34:33 > .py 10:34:25; mathematica .txt 10:26:31 > .wl 10:22:00). No
scripts were run.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `A_F1 = 1.000051928802195…` | py L40 / wl L41; sympy out L18, wl out L23 | tex:7 & tex:16 (`1.00005192880220`); notes:94, notes:103 | MATCH |
| `rho_min = 4/3 = 1.33333…` | py L31 / wl L34; sympy out L12 | notes:36/140 (`rho_alpha = 4/3`); tex via Stage 088 input | MATCH |
| `zeta_min = 1/3 = 0.33333…` | py L32 / wl L35; sympy out L16 | tex:7 & tex:22 (`zeta_req=1/3`); notes:63 | MATCH |
| `rho_suff = 3.46622291347846` | py L72 / wl L70; sympy out L13 | notes:32 | MATCH |
| `rho_fail = 3.46752913273870` | py L73 / wl L71; sympy out L14 | notes:33 | MATCH |
| `rho_max = 3.46752922945601` | py L74 / wl L72; sympy out L15 | notes:34 | MATCH |
| `zeta_max = 2.46752922945601` | py L46 / wl L45; sympy out L17 | notes:67 | MATCH |
| `Delta_suff = 2.13288958014513` | py L94; sympy out L21 | notes:43-44 | MATCH |
| `Delta_fail = 2.13419579940537` | py L95; sympy out L22 | notes:46 | MATCH |
| `Delta_max = 2.13419589612268` | py L96; sympy out L23 | notes:48-49 | MATCH |
| `Delta_zeta = 2.13419589612268` | py L97; sympy out L24 | notes:71-72 (`zeta_max - zeta_req ≈ 2.13419589612268`) | MATCH |
| `Omega(Pe→0) = 1` | py L53 / wl L51; sympy out L5 | tex:13 (`Omega_{Pe=0}=1`); notes:20 | MATCH |
| `zeta_F1(0) = A_F1` | py L54 / wl L52; sympy out L6 | tex:16; notes:90+94 | MATCH |
| `Pe_req = 0` (boxed Output) | py L128 / wl L112; sympy out L26 | tex:28 (boxed); notes:107, notes:141 | MATCH |

Internal scaffolding (accounted for, no finding expected in prose):
`Delta_AF1 = 0.666718595468…` (auxiliary margin `A_F1 - zeta_min`, not a stated deliverable);
`Pe_suff_chi = 96.5285247264386`, `Pe_fail_chi = 11220.5441626259` (upstream Stage-082
carry-forward Pe values, provenance-commented, not Stage-089 deliverables); `kappa_F1 =
12321/5`, `eta_F1 = 37`, `y_F1` (Stage-079 A_F1 ingredients); `eps_blk = 0` (blocking
parameter); residual/diff/PASS scaffolding.

**reconciliation: complete; 14 deliverable values checked, 0 misaligned.**

The two findings (F1, F2) are `tautological_check`, NOT `paper_misalignment` — no value
mismatch and no missing deliverable surfaced in the reconciliation. `needs_user_resolution:
false`.
