---
unit_id: 114
batch: IV.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage114_concrete_core_schur.md]
  paper_appendix: present
---

# Audit unit 114 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_114.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage114_concrete_core_schur.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (line 28 audit-path summary; line 1262 input)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage114_concrete_core_schur_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage114_concrete_core_schur_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage114_concrete_core_schur_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage114_concrete_core_schur_mathematica_audit.txt`

## What the paper claims

The stage card (line 16) states: "Concrete shell/mixed core Schur complement reproduces the reduced Robin--mixed outlet." The notes (the authoritative source, since the .tex is terse) make this concrete with two boxed identities:

1. The exact Schur complement of the 2x2 core system `[[K_s, λ], [λ, −K_q D_W^bare(z)]] (s, q)^T = u(g_s, g_q)^T` yields
   `δΛ_core(z) = g_s^2/K_s − (K_s g_q − λ g_s)^2 / [K_s (K_s K_q D_W^bare(z) + λ^2)]`.
2. With `r_c := λ^2/(K_s K_q)` and `D_W^bare(z) = 1 − κ_0 z^2 − iγ_0 z^5 + O(z^6)`, this rearranges to the reduced Stage-95 form
   `δΛ_core(z) = ρ_c − σ_c / (1 − κ_c z^2 − iγ_c z^5) + O(z^6)`
   with `ρ_c = g_s^2/K_s`, `σ_c = (K_s g_q − λ g_s)^2 / (K_s^2 K_q (1 + r_c))`, `κ_c = κ_0/(1 + r_c)`, `γ_c = γ_0/(1 + r_c)`.

Part-IV audit-path map (appendix line 28) lists this as "two-channel core realization."

## What the script claims to verify

Both engines build the same 2x2 matrix `M = [[K_s, λ], [λ, −K_q D]]`, form `c^T M^{-1} c` with `c = (g_s, g_q)`, and reduce via `Apart` in D. They then assert two identities:

- **A1 ("Schur form identity"):** `δΛ(D) − [ρ_c − σ_tilde/(D + r_c)] == 0`, where `σ_tilde = (K_s g_q − λ g_s)^2 / (K_s^2 K_q)`. This verifies that the abstract Schur complement matches the boxed unreduced form (after pulling K_s out of the denominator).
- **A2 ("low-frequency normalized outlet identity"):** `δΛ(D)|_{D = D_W^bare} − [ρ_c − σ_c/(1 − κ_c z^2 − iγ_c z^5)] == 0`. This verifies the prefactor absorption: `(1 + r_c)` rolls into the `(κ_c, γ_c, σ_c)` definitions, so the substituted form matches the reduced Robin–mixed shape.

The script then prints the symbolic identifications of `ρ_c, σ_c, κ_c, γ_c`.

## Paper ↔ script cross-check

| Paper-side deliverable | Script check | Status |
|---|---|---|
| Exact Schur complement boxed form (notes eq. before `r_c`) | A1 | match |
| Reduced Robin–mixed shape with absorbed prefactors (notes second box) | A2 | match |
| Identifications `ρ_c, σ_c, κ_c, γ_c` (notes boxes) | Print statements + the σ_c, κ_c, γ_c definitions used in A2 | match (the definitions feed A2, so an off-value would break A2) |
| Stage card heuristic "checks" (Schur signs, L/a normalization, parent overlap) | A1 implicitly exercises signs (the −K_q D entry); L/a and parent overlap are not in scope of this stage's formal Output | n/a (heuristic guidance, not Output) |

The paper's `Output` line (the boxed claim in the stage card) is "Concrete shell/mixed core Schur complement reproduces the reduced Robin–mixed outlet." A1+A2 jointly cover this. Paper alignment is `aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 42 | `expect_zero("Schur form identity", delta_D - target_D)` | Exact Schur complement boxed form | yes |
| A2 | sympy | 49 | `expect_zero("low-frequency normalized outlet identity", delta_z - target_z)` | Prefactor absorption to reduced Robin–mixed form | yes |
| A1' | mathematica | 47 | `expectZero["Schur form identity", deltaD - targetD]` | Same as A1 | yes |
| A2' | mathematica | 52 | `expectZero["low-frequency normalized outlet identity", deltaZ - targetZ]` | Same as A2 | yes |

Each assertion compares an independently-computed quantity (matrix-inverse Schur complement) against an independently-constructed target from the notes' boxed equations. Each could fail under a sign flip, factor change, or off-by-one in any of `K_s, K_q, λ, g_s, g_q, κ_0, γ_0`. Non-tautological.

## Findings

(none)

## Independent-derivation check (Mathematica)

The Mathematica script follows the same algebraic choreography as the SymPy script: same matrix, same `Inverse[m]`, same `Apart` in `dSym`, same `targetD`, same substitution `dSym -> dBare`, same `targetZ`. Variable names differ trivially (`K_s` ↔ `kS`, `lam` ↔ `lam`, `g_s` ↔ `gS`, etc.). Helper functions `banner`/`expectZero` are pure mirrors of the SymPy `banner`/`expect_zero`.

I considered filing this as `mathematica_transliteration`. However, the algebraic content of this stage is essentially a single linear-algebra fact — the Schur complement of a 2x2 matrix — followed by a substitution and prefactor absorption. The set of "natural derivation paths" is therefore narrow: any engine starting from the 2x2 matrix must build the same inverse and obtain the same rational function in D. The Mathematica engine does *not* re-derive via, e.g., a continued-fraction expansion or a different elimination order, but for a 2x2 matrix this is not a meaningfully separable choice. The check that ultimately matters — does both `Apart[c.Inverse[m].c, dSym]` (Mathematica) and `apart((c.T M.inv() c)[0], D)` (SymPy) agree with the boxed `ρ_c − σ_tilde/(D + r_c)` form? — is exercised independently by each engine's symbolic kernel. I am not filing this as a finding; flagging here for the verifier's awareness.

## Engine cross-check

Both engines reach the same final identifications:

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| `δΛ(D)` | `g_s^2/K_s − (K_s g_q − g_s lam)^2 / [K_s (D K_q K_s + lam^2)]` | `(dSym gS^2 kQ − gQ^2 kS + 2 gQ gS lam)/(dSym kQ kS + lam^2)` |
| `ρ_c` | `g_s^2/K_s` | `gS^2/kS` |
| `σ_c` | `(K_s g_q − g_s lam)^2 / (K_s (K_q K_s + lam^2))` | `(gQ kS − gS lam)^2 / (kS (kQ kS + lam^2))` |
| `κ_c` | `K_q K_s κ_0 / (K_q K_s + lam^2)` | `(κ_0 kQ kS) / (kQ kS + lam^2)` |
| `γ_c` | `K_q K_s γ_0 / (K_q K_s + lam^2)` | `(γ_0 kQ kS) / (kQ kS + lam^2)` |

The `δΛ(D)` forms appear different but are algebraically identical:
`g_s^2/K_s − (K_s g_q − g_s lam)^2 / [K_s (D K_q K_s + lam^2)]`
`= [g_s^2 (D K_q K_s + lam^2) − (K_s g_q − g_s lam)^2] / [K_s (D K_q K_s + lam^2)]`
`= [D g_s^2 K_q K_s − K_s^2 g_q^2 + 2 K_s g_q g_s lam] / [K_s (D K_q K_s + lam^2)]`
`= [D g_s^2 K_q − K_s g_q^2 + 2 g_q g_s lam] / (D K_q K_s + lam^2)`

which matches the Mathematica combined-fraction form. Both A1 and A2 simplify to literal `0` in both engines. `engines_agree: true`.

Output mtimes: sympy script Apr 1 → sympy output May 11 (fresh); mathematica script May 11 11:56 → output May 11 13:08 (fresh). `outputs_fresh: true`.

## Verdict justification

The script holds up. Both engines independently form the Schur complement of the 2x2 core matrix, reduce it via `apart`/`Apart` in D, and assert it equals the boxed unreduced form from the notes (A1). Both then substitute `D_W^bare(z) = 1 − κ_0 z^2 − iγ_0 z^5` and assert that the result equals the reduced Robin–mixed form `ρ_c − σ_c/(1 − κ_c z^2 − iγ_c z^5)` once the prefactor `(1 + r_c)` is absorbed into the constants (A2). Both assertions are non-tautological — a sign flip on the `−K_q D` entry, a factor flip on any of `ρ_c, σ_c, κ_c, γ_c`, or an off-by-one in the σ_tilde numerator would cause failure. Both engines agree on every printed identification; the apparently different forms of `δΛ(D)` are algebraically identical (verified by hand-combining the SymPy form). Outputs are fresh. Paper card claim, notes boxed equations, and script assertions all line up.

Attacks tried that failed: (a) checked whether A2 is a tautology after A1 — it isn't, because A2 exercises the `(1+r_c)` absorption that A1 leaves implicit (A1 uses `σ_tilde = σ_c (1+r_c)` and divisor `D + r_c`, while A2 uses `σ_c` and divisor `1 − κ_c z^2 − iγ_c z^5`); (b) checked symbol assumptions — `z, kappa0, gamma0` are declared `positive=True, real=True`, which is more restrictive than needed but does not invalidate the algebraic identity (the identity is rational and holds for any z where the denominator is nonzero); (c) checked the stage-numbering inconsistency (banner says "STAGE 97/097", filename and paper card say "Stage 114/131") — this is a relabeling artifact, not a math finding; (d) checked whether the paper's `+ O(z^6)` truncation is honored — the script defines `D_W^bare = 1 − κ_0 z^2 − iγ_0 z^5` exactly (no extra terms), so the identity becomes exact (not just asymptotic), which is consistent with the notes' boxed form.

## Self-test notes

I checked: (1) variable independence is not relevant here (no derivatives are taken). (2) Symmetry/parity is not relevant (no unbounded integrals). (3) Trivial-case pre-check: at `λ = 0`, `D = 1`, the matrix decouples; the Schur complement gives `g_s^2/K_s − g_q^2/K_q` for the unreduced form and equivalent values for the reduced form — both checks reduce to consistent algebraic identities. (4) Path specs — no missing-script findings. (5) Paper round-trip — no prescribed fix; nothing to round-trip.
