---
unit_id: 014
batch: I.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-04T00:00:00Z
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
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 014 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_014.tex`
- notes: `(none)` (no file matching `notes/stages/moving_throat_pde_stage014_*.md`)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 50)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_mathematica_audit.txt`

## What the paper claims

Stage 014 ("Mouth-Taylor gate bridge", Part I, anchor MTDC-T4) turns the
mouth-Taylor primitive map (from stage 013) into a gate-level compensation test.
The card states the conservative mouth shifts enter the 5PN even gates as
`K_1 = -z_2 - z_0/9` (eq:stage014-k1) and `H_even = -z_4 + (2/3) z_2 - z_0/27`
(eq:stage014-heven). Solving the even-gate pair against the primitive derivative
slots yields a compensation denominator "proportional to" `Q S_2 - Delta H_port`
(eq label `eq:stage007-compensation-denominator`, line 27). The `\stagefield{Output}`
reads verbatim: "Stage~014 exports the mechanism sieve: projected EM mouth data
can tune the even gates only away from the degeneracy" `Q S_2 - Delta H_port`.
The appendix row (line 50) summarizes: "Gate conditions for carrying mouth-local
projected data into the grouped bundle." Deliverables: (D1) the two even-gate
definitions in z-slot form; (D2) the compensation denominator `Q S_2 - Delta H_port`
and the fact that the compensation becomes singular on `Q S_2 = Delta H_port`
(equivalently the Z2 slot vanishes); (D3) the mechanism-sieve statement that no
pure single-sector correction closes the gates.

## What the script claims to verify

The SymPy script builds the primitive one-port slot expressions `z0, z2, z4`
(lines 45-47) in terms of primitive slips `q1, s1, h1, d1`, applies a one-sided
mouth-Taylor lift `slip = mu1 * X'` (line 50), divides out `mu1` to get the
mouth-derivative map `z0d, z2d, z4d`, and forms the even gates `K1 = -(z2 + z0/9)`
and `He = -z4 + (2/3) z2 - z0/27` (lines 57-58) — matching the card. It then
asserts: z2d/z4d are linear in each primitive slip (lines 80-82); K1/He are
independent of the spurious P', G_W' slips (lines 89-91); a "bundle pull-back"
rebuild of K1/He from z0d/z2d/z4d equals the directly-substituted forms (lines
94-97); the two sector Jacobian determinants are nonzero (lines 98-99); the
compensation surface denominators equal `9 Delta^2 (Delta H_port - Q S2)` and
`9 Delta (Delta H_port - Q S2)` (lines 100-101), restated as `-9 Delta^4 Z2`
and `-9 Delta^3 Z2` (lines 109-110); a sign-flipped denominator mutation is
nonzero (line 102, adversarial); and the pure single-sector solves yield only
the trivial solution (lines 104-107). The Mathematica script re-derives `z0d,
z2d, z4d` from a different starting point — compact closed-form Z-slots `Z2 =
(Q S2 - Hport Delta)/Delta^2`, `Z4 = (Q(S2^2 - Delta) - S2 Hport Delta)/Delta^3`
(lines 53-54) — via chain-rule `ell`-derivative, then mirrors the independence,
sieve-solve, Jacobian-determinant, and compensation-denominator checks (M1, M5-M10).

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| D1: `K_1 = -z_2 - z_0/9`, `H_even = -z_4 + (2/3) z_2 - z_0/27` | sympy 57-58 define exactly these; bundle pull-back consistency 94-97; Mathematica K1/He 82-83 | match |
| D2: compensation denom `Q S_2 - Delta H_port`, singular at `Z2 = 0` | sympy 100-101 + 109-110 assert denom = `9 Delta^{2,3}(Delta H_port - Q S2)` = `-9 Delta^{4,3} Z2`; Mathematica M9 119-123; output text 57-62 | match (sign: card "proportional to `Q S2 - Delta Hport`", script "proportional to `Delta Hport - Q S2`" — differ only by overall sign, identical degeneracy locus) |
| D3: no pure single-sector correction closes the gates | sympy qd_only/sh_only solves 61-62 asserted `[{Qx:0,Dx:0}]`/`[{Sx:0,Hx:0}]` at 104-107; Mathematica M5/M6 96/103; Jacobian dets M7/M8 nonzero | match |
| (extra) P'/G_W' independence of K1/He | sympy 89-91; Mathematica M1 86-88 | extra (supports the "only (Q',Delta',S2',H_port')" readout, line 145-147; consistent with card, not a misalignment) |

Dominant pattern: aligned.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 80-82 | `assert_zero(diff(z2d/z4d, slip, 2))` | D1 (linearity of lift) | yes |
| A2 | sympy | 89-91 | `assert_zero(diff(K1/He, Px/Gx))` | D3 readout (slip-set) | yes |
| A3 | sympy | 95 | `assert_zero(K1 - K1_bundle)` | D1 (pull-back consistency) | yes |
| A4 | sympy | 97 | `assert_zero(He - He_bundle)` | D1 | yes |
| A5 | sympy | 98-99 | `assert_nonzero(qd/sh det)` | D3 (sector solvability) | yes |
| A6 | sympy | 100-101 | `assert_zero(Hx_den/Sx_den - 9 Delta^k (Delta Hport - Q S2))` | D2 | yes |
| A7 | sympy | 102 | `assert_nonzero(Hx_den - 9 Delta^2 (Delta Hport + Q S2))` | D2 (sign mutation) | yes |
| A8 | sympy | 104-107 | `if qd_only/sh_only != trivial: raise` | D3 | yes |
| A9 | sympy | 109-110 | `assert_zero(Hx_den/Sx_den + 9 Delta^k Z2)` | D2 (Z2 restatement) | yes |
| M1 | mathematica | 86-87 | `expectZero(D[K1/He, Px/Gx])` | D3 readout | yes |
| M5 | mathematica | 96 | `expectEqual(qdSolve, {{Qx->0,Deltax->0}})` | D3 | yes |
| M6 | mathematica | 103 | `expectEqual(shSolve, {{Sx->0,Hx->0}})` | D3 | yes |
| M7/M8 | mathematica | 110-111 | `expectNonzero(Det jacobianQD/SH)` | D3 | yes |
| M9 | mathematica | 117-124 | `expectZero(hxDen/sxDen - 9 Delta^k (Delta Hport - Q S2))` | D2 | yes |
| M10 | mathematica | 125-128 | `expectNonzero(den - 9 Delta^2 (Delta Hport + Q S2))` | D2 (sign mutation) | yes |

No tautological rows. Every assertion traces to a paper deliverable.

## Findings

None. The script holds up against every attack tried (see Verdict justification).

## Independent-derivation check (Mathematica)

The Mathematica script is an **independent re-derivation**, not a transliteration.
The decisive evidence is the `z4d` cross-route check. SymPy obtains `z4d` by
substituting Taylor slips into the *expanded* primitive-slip numerator
(`scripts/...:47`: `-Delta^2 Hport s1 - Delta^2 S2 h1 - Delta^2 q1 + 2 Delta Hport S2 d1
+ 2 Delta Q S2 s1 + 2 Delta Q d1 + Delta S2^2 q1 - 3 Q S2^2 d1` over `Delta^4`).
Mathematica instead defines a *compact closed form* `Z4 = (Q(S2^2 - Delta) -
S2 Hport Delta)/Delta^3` (`mathematica/...:54`) and takes the chain-rule
`ell`-derivative of `Z4(Q + mu1 Qx ell, ...)`. I expanded the Mathematica chain
rule by hand: `∂Z4/∂Q = (S2^2 - Delta)/Delta^3`, `∂Z4/∂S2 = (2 Q S2 - Hport Delta)/Delta^3`,
`∂Z4/∂Hport = -S2/Delta^2`, `∂Z4/∂Delta = (2 Q Delta + 2 S2 Hport Delta - 3 Q S2^2)/Delta^4`;
combining gives exactly SymPy's printed `z4d` (output line 20) term-for-term.
Likewise `Z2` and the K1/He combinations reproduce SymPy's `z2d`/K1/He. The two
engines reach the same map from structurally different starting expressions (an
expanded slip numerator vs. a factored Z-slot), which is the hallmark of a true
second derivation. `Z0/Z2/Z4` use different variable choreography and a `D[...,ell]`
chain rule rather than SymPy's `.subs(...)/mu1` substitution — not a line-by-line port.

## Engine cross-check

Both engines agree at every comparable point:
- `z2d`, `z4d`, `K1`, `He` agree (hand-verified above and against output lines 17-28).
- Compensation denominators: SymPy A6/A9 and Mathematica M9 both reduce to 0 (sympy
  output is silent-on-pass; Mathematica output lines 18-21 print `residual = 0`, PASS).
- Sign-flip mutation: SymPy A7 residual is `-18 Delta^2 Q S2` (nonzero); Mathematica
  M10 output line 22 prints `residual = -18*Delta^2*Q*S2`. **Identical residual.**
- Sector solves: SymPy `[{Deltax:0,Qx:0}]`/`[{Hx:0,S2x:0}]` (output lines 40/43)
  ↔ Mathematica M5/M6 `True` (lines 10-13). Jacobian determinants: Mathematica M8
  prints `(-(Delta Hport) + Q S2)/Delta^4` (output line 16), the Z2-slot numerator
  up to the `Delta^4` weight — consistent with the SymPy `assert_nonzero` and with D2.

## Verdict justification

`clean`. I read the card, the (absent) notes, and the appendix row first, built
the model (even-gate definitions, compensation denominator `Q S2 - Delta Hport`,
mechanism-sieve), then attacked the scripts. Attacks tried and failed: (1) sign/factor
audit of `z2`,`z4` and the K1/He combinations against the card eqs — all match,
including the `2/3` and `1/9`,`1/27` weights; (2) tautology hunt — the bundle
pull-back checks (A3/A4) compare a `.subs`-substituted form against an independently
rebuilt-from-`z*d` form, so they can fail; the denominator checks (A6/A9) compare
a solved-surface denominator against a hand-written target, also non-trivial; (3)
the docstring's `D`-notation gates (`K1 = D21 + D01/9`) initially looked sign-flipped
vs. the code (`-(z2 + z0/9)`) and the card, but the code and card agree exactly and
the docstring `D`-symbols are a separate (negated) gate-slot naming convention, not
a paper claim — no misalignment; (4) symbol-domain check — SymPy `Q,S2,Hport,Delta`
`nonzero`, slips generic; Mathematica reals with `Delta,D0,mu1 != 0` — both consistent
with the physical mouth-Taylor setup and with the implicit `Delta Hport - Q S2 != 0`
the solves require; (5) engine-disagreement probe — the sign-flip mutation residual
`-18 Delta^2 Q S2` is byte-identical across engines, and the independently-routed
`z4d` matches. Adversarial mutation tests (A7/M10) are present and pass, which raises
confidence the denominator checks are real. The only cosmetic oddity (not a finding):
the equation label `eq:stage007-compensation-denominator` is defined and referenced
*entirely within* stage_014.tex (it resolves to local `B.93` in the `.aux`); the
`stage007` token in the label name is a legacy-naming artifact with no broken
reference and no content/value impact.

## Value Reconciliation (pass-2 augmentation)

All emitted deliverable values reconcile. The two `.tex` body equations match the
script's gate definitions exactly; the compensation denominator matches the card
up to the explicitly-stated "proportional to" overall sign (identical degeneracy
locus). There is no notes file; the terse card legitimately omits the intermediate
closed forms (`z2d`,`z4d`,`comp_surface`), which are scaffolding for the boxed
denominator/gate results and not stated deliverables.

`reconciliation: complete; 5 values checked, 0 misaligned`

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `K_1 = -z_2 - z_0/9` (even gate def) | py:57; wl:82; out (py) line 25-26 | `stage_014.tex:16-17` (eq:stage014-k1) | MATCH |
| `H_even = -z_4 + (2/3) z_2 - z_0/27` (even gate def) | py:58; wl:83; out (py) line 27-28 | `stage_014.tex:19-20` (eq:stage014-heven) | MATCH |
| compensation denominator `Delta H_port - Q S_2` (script sign) / `Q S_2 - Delta H_port` (card sign) | py:100-101,109-110; wl:119-123; out (py) line 57-58; out (wl) line 18-21 | `stage_014.tex:27` (eq label compensation-denominator) | MATCH (proportionality; card and script differ only by overall sign, same degeneracy `Q S2 = Delta H_port`) |
| degeneracy / `Z2 = (Q S2 - H_port Delta)/Delta^2` singular surface | py:108-110; wl:53; out (py) line 59-62 | `stage_014.tex:27,30-33` (Output: tune "only away from the degeneracy") | MATCH |
| mechanism sieve: no pure single-sector solution (`qd_only=[{Qx:0,Dx:0}]`, `sh_only=[{Sx:0,Hx:0}]`) | py:104-107; wl:96,103; out (py) lines 40,43; out (wl) lines 10-13 | `stage_014.tex:30-31` (Output: "exports the mechanism sieve") + appendix row 50 | MATCH |

INTERNAL (scaffolding, no prose expected): `z0d`, `z2d`, `z4d` (intermediate
mouth-derivative map); `K1_bundle`/`He_bundle` (pull-back rebuild); `qd_matrix`/
`sh_matrix` Jacobian determinants and M7/M8 residuals; `comp_surface[Hx]`/`comp_surface[Sx]`
full closed forms; sign-flip mutation residual `-18 Delta^2 Q S2`; `D0`/`mu1`
auxiliary symbols; PASS/FAIL flags and `residual = 0` lines.

## Self-test notes

Variable-independence: the `diff(K1/He, Px/Gx)` checks (A2/M1) target slips that
genuinely do not appear in K1/He, so the assert_zero is a real (and correct)
independence result, not a vacuous-derivative trap; the bundle pull-back compares
two independently constructed expressions (no `diff` involved). Symmetry/parity:
n/a — no unbounded-domain integrals; all checks are rational-function identities.
Trivial-case: the sign-flip mutation (A7/M10) correctly yields nonzero `-18 Delta^2 Q S2`,
and the denominator targets reduce to 0 against the solved surface; substituting a
generic nonzero `Delta` keeps every denominator finite. No fix prescribed (clean), so
no paper round-trip needed.
