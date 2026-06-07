---
unit_id: 138
batch: IV.4
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00-06:00
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage138_normalized_mouth_gain_family.md]
  paper_appendix: present
---

# Audit unit 138 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_138.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage138_normalized_mouth_gain_family.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (read anchor block at lines 19, 1175-1179, and the `\input{stages/stage_138}` row at line 1310; no per-stage status row beyond the card inclusion)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage138_normalized_mouth_gain_family_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage138_normalized_mouth_gain_family_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage138_normalized_mouth_gain_family_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage138_normalized_mouth_gain_family_mathematica_audit.txt`

## What the paper claims

The stage rewrites the explicit mouth-gain map of Stage 137 in the core-normalized parent variables (`\mathfrak r := λ/√(K_sK_q)`, `\mathfrak g_c := g_q√(K_s)/(g_s√(K_q))`, `Σ₀ := L g_s²/(K_s Θ_σ)`). The card's load-bearing line (line 16) states verbatim: "Gain ratio \(R_q=(\mathfrak g_c-\mathfrak r)^2/(1+\mathfrak r^2)\); exact compensation gives \(R_q=1/4\)." The notes expand this into the full deliverable set: `M_s = Σ₀` (notes:25); `M_q = -Σ₀ (g_c-r)²/(1+r²)` (notes:35-40); the exact mixed-to-shell ratio `R_q := -M_q/M_s = (g_c-r)²/(1+r²)` (notes:45/99); insertion of the Stage-115 core-balance family `g_c = r ± ½√(1+r²)` collapsing the ratio to `R_q = 1/4` (notes:69/104); and the corollaries on the compensated branch `M_s = Σ₀, M_q = -Σ₀/4` (notes:75-77) equal to the Stage-135 closure `M_s = 4Σ_m, M_q = -Σ_m` with `Σ_m = Σ₀/4` (notes:82-86). The distinct script-checkable deliverables are: (D1) the normalized `M_s = Σ₀`; (D2) the normalized `M_q = -Σ₀(g_c-r)²/(1+r²)`; (D3) the exact ratio formula `R_q = (g_c-r)²/(1+r²)`; (D4) `R_q = 1/4` on the `+` compensation branch; (D5) `R_q = 1/4` on the `-` compensation branch. The Stage-135 closure / `Σ_m` corollaries are downstream re-labelings that follow trivially from `R_q=1/4` and are not independent script targets.

## What the script claims to verify

Both scripts start from the dimensional gain laws `M_s = L g_s²/(K_s Θ)` and `M_q = -L(K_s g_q - λ g_s)²/(K_s(K_s K_q + λ²)Θ)`, apply the normalization substitution `{λ → r√(K_sK_q), g_q → g_c g_s √(K_q/K_s)}`, simplify, and form `R_q = -M_q/M_s`. They PRINT `M_s normalized = Σ₀`, `M_q normalized = -Σ₀(g_c-r)²/(1+r²)`, and `R_q = (g_c-r)²/(1+r²)`. The SymPy script then asserts (load-bearing) that substituting both branches of the compensation family `g_c = r ± ½√(1+r²)` reduces `R_q` to `1/4` (py:31-32). The Mathematica script asserts the same two branch identities (wl:51-52) and ADDITIONALLY asserts that the derived `rQ` equals the closed form `(g_c-r)²/(1+r²)` (wl:45, `expectZero["R_q exact formula", rQ - (gCore-rHat)^2/(1+rHat^2)]`) — a check the SymPy version only prints rather than asserts.

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| D1 `M_s = Σ₀` | py:22 / wl:42 print; derived as `Ms_norm = simplify(Ms.subs(subs))` then displayed as `Σ₀` | match (printed; the `Σ₀` label is the definition of `L g_s²/(K_s Θ)`) |
| D2 `M_q = -Σ₀(g_c-r)²/(1+r²)` | py:23 / wl:43 print of `Mq_norm` with `Ms_norm → Σ₀`; outputs show `-Σ₀(g_c-r)²/(r²+1)` | match |
| D3 `R_q = (g_c-r)²/(1+r²)` | py:24 print; wl:44 print + wl:45 `expectZero` assertion | match (mathematica asserts; sympy prints) |
| D4 `R_q = 1/4` on `+` branch | py:31 `assert`, wl:51 `expectZero` | match |
| D5 `R_q = 1/4` on `-` branch | py:32 `assert`, wl:52 `expectZero` | match |

`paper_alignment: aligned`. Every card/notes deliverable maps to a script-side print or assertion; the load-bearing compensation result (`R_q=1/4`) is asserted by both engines, on both branches.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 31 | `assert simplify(sol_plus - 1/4) == 0` | D4 (`R_q=1/4`, `+` branch) | yes |
| A2 | sympy | 32 | `assert simplify(sol_minus - 1/4) == 0` | D5 (`R_q=1/4`, `-` branch) | yes |
| A3 | mathematica | 45 | `expectZero["R_q exact formula", rQ - (gCore-rHat)^2/(1+rHat^2)]` | D3 (exact ratio formula) | yes |
| A4 | mathematica | 51 | `expectZero["R_q on + branch - 1/4", solPlus - 1/4]` | D4 | yes |
| A5 | mathematica | 52 | `expectZero["R_q on - branch - 1/4", solMinus - 1/4]` | D5 | yes |

Note: D1/D2 (the normalized `M_s`, `M_q` forms) are PRINTED but never `assert`-guarded in either engine; they are exercised only by the derivation chain that feeds the asserted `R_q`. This is acceptable here because `R_q = -M_q/M_s` and the asserted `R_q` identities (A3-A5) cannot hold unless the upstream `M_q`/`M_s` normalization is correct — i.e. the printed forms are load-bearing inputs to the asserted checks, not orphans. I considered raising `insufficient_verification` for the un-asserted `M_s`/`M_q` prints but rejected it: the closed-form ratio `R_q=(g_c-r)²/(1+r²)` is the card's single stated bottom line and IS asserted (A3 in Mathematica), and the branch-collapse to `1/4` is the load-bearing physics and IS asserted in both engines.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` shares the algebraic choreography of the `.py`: identical starting expressions (`mS`, `mQ`), identical substitution dictionary (`{lam -> rHat*Sqrt[kS*kQ], gQ -> gCore*gS*Sqrt[kQ/kS]}`), and the identical no-op tail substitution `/. mSNorm -> sigma0` (which, after `rQRaw` is already fully simplified to `(gCore-rHat)²/(1+rHat²)`, does nothing — `mSNorm` no longer appears). That copied no-op is a mild transliteration tell. However, the derivation here is a single substitute-and-simplify whose structure is essentially forced by the physics, the `.wl` performs its own independent `FullSimplify[Together[Expand[...]]]` reductions rather than echoing a SymPy-side simplified literal, and the `.wl` adds a genuinely independent assertion (`R_q exact formula`, wl:45) that the `.py` does not contain. Side-by-side example: sympy `Rq_raw = sp.simplify(-Mq_norm / Ms_norm)` (py:19) vs. mathematica `rQRaw = FullSimplify[-mQNorm/mSNorm, ...]` (wl:39) — same target, distinct CAS engine. I judged this short of a `mathematica_transliteration` finding: both engines independently reduce the same physical premises and the second engine carries an extra check; the structural overlap is the minimum the math permits, not evidence of a non-independent port. (Noted, not filed.)

## Engine cross-check

Both engines agree exactly:

| quantity | SymPy output | Mathematica output |
|---|---|---|
| `M_s normalized` | `Sigma_0` (line 9) | `sigma0` (line 5) |
| `M_q normalized` | `-Sigma_0*(g_c - r)**2/(r**2 + 1)` (line 10) | `-(((gCore - rHat)^2*sigma0)/(1 + rHat^2))` (line 6) |
| `R_q` | `(g_c - r)**2/(r**2 + 1)` (line 11) | `(gCore - rHat)^2/(1 + rHat^2)` (line 7) |
| `R_q` exact-formula residual | (printed only) | `0` → PASS (line 8-9) |
| `R_q` + branch | `1/4` (line 12) | `1/4` → PASS (lines 11-13) |
| `R_q` − branch | `1/4` (line 13) | `1/4` → PASS (lines 14-15) |

`engines_agree: true`.

## Verdict justification

Clean. I attacked the substitution algebra by hand: `K_s g_q - λ g_s` under the subs becomes `g_s√(K_sK_q)(g_c-r)`, squared gives `g_s²K_sK_q(g_c-r)²`, and the denominator `K_s(K_sK_q+λ²)Θ = K_s²K_q(1+r²)Θ`, so `M_q/M_s = -(g_c-r)²/(1+r²)` — matching both engines and the notes exactly. The branch assertions are non-tautological: `R_q` is derived symbolically (not pinned to a literal) and only collapses to `1/4` because `(±½√(1+r²))²/(1+r²) = ¼`; substituting any other `g_c` would not give `1/4`, so A1/A2/A4/A5 can genuinely fail. I checked the no-op `subs(Ms_norm, Sigma0)` (py:20 / wl:40) — harmless, it changes nothing in an already-simplified `R_q`, and the corresponding `M_q normalized` print correctly relabels the surviving `L g_s²/(K_s Θ)` factor as `Σ₀`. Symbol domains (`r, g_c` real; `1+r²` never zero) introduce no branch hazard. Outputs are fresh (both `.txt` mtimes post-date their scripts). I read the card, the notes, and the appendix anchor block; the script's verified claim (`R_q = (g_c-r)²/(1+r²)`, collapsing to `1/4` on the compensation family) matches the paper's stated `\StatusExactClosure` claim exactly. The only un-asserted deliverables (`M_s`, `M_q` prints) are load-bearing inputs to the asserted `R_q` checks, so the coverage gap is cosmetic, not substantive.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 5 deliverable values checked, 0 misaligned

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `M_s = Σ₀` | py:22, wl:42; sympy out:9 (`M_s normalized = Sigma_0`), math out:5 | notes:25 `\boxed{M_s=\Sigma_0.}`; also notes:75 | MATCH |
| `M_q = -Σ₀(g_c-r)²/(1+r²)` | py:23, wl:43; sympy out:10, math out:6 | notes:35-40 `M_q = -\Sigma_0 (\mathfrak g_c-\mathfrak r)^2/(1+\mathfrak r^2)` | MATCH |
| `R_q = (g_c-r)²/(1+r²)` | py:24, wl:44-45; sympy out:11, math out:7 | card:16 `R_q=(\mathfrak g_c-\mathfrak r)^2/(1+\mathfrak r^2)`; notes:45/99 | MATCH |
| `R_q = 1/4` (`+` branch) | py:27/29/31, wl:47/49/51; sympy out:12, math out:11 | card:16 `exact compensation gives R_q=1/4`; notes:69/104 | MATCH |
| `R_q = 1/4` (`-` branch) | py:28/30/32, wl:48/50/52; sympy out:13, math out:11 | card:16; notes:69 | MATCH |

Notes-only corollaries NOT separately emitted by the scripts (and therefore not reconciliation targets, but cross-checked as consistent): `M_q = -Σ₀/4` on the compensated branch (notes:77 — follows from `R_q=1/4`); Stage-135 closure `M_s = 4Σ_m, M_q = -Σ_m` with `Σ_m = Σ₀/4` (notes:82-86). The scripts emit `R_q=1/4`, from which these follow by definition; no script-emitted value is missing from the docs.

INTERNAL scaffolding (accounted for, no finding): pass/fail flags (`PASS:`/`FAIL:`), `Normalized mouth-gain family verified.` banner, `expectZero` residual `= 0` lines, intermediate symbols `Ms_norm`/`mSNorm`/`Mq_norm`/`mQNorm`/`Rq_raw`/`rQRaw`/`subs` (exist only to drive the printed/asserted results), banner string, exit codes.

Every emitted deliverable value reconciles with the card and/or notes. No `value_mismatch` and no `script_missing_paper_claim` raised.

## Self-test notes

(1) Variable independence: no `sp.diff`/`D[]` in either script — the only operations are substitution and simplify, so the zero-derivative trap does not apply. (2) Symmetry/parity: no integrals over unbounded domains — N/A. (3) Trivial-case pre-check: I substituted `g_c = r + ½√(1+r²)` into `(g_c-r)²/(1+r²)` by hand and confirmed it reduces to `¼(1+r²)/(1+r²) = 1/4`; the `-` branch gives the same since the term is squared; both `assert`/`expectZero` are genuinely satisfiable-or-falsifiable, not trivially true. (4) Path specs: N/A (no missing-script directive). (5) Paper round-trip: no fix prescribed, so no new misalignment introduced. Conclusion: no findings; both engines agree; reconciliation complete.
