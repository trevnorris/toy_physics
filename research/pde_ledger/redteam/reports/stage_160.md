---
unit_id: 160
batch: IV.6
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage160_bare_mixed_port_slippage.md]
  paper_appendix: present
---

# Audit unit 160 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_160.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage160_bare_mixed_port_slippage.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only `\input{stages/stage_160}` row at line 1354 — no per-stage prose summary)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage160_bare_mixed_port_slippage_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage160_bare_mixed_port_slippage_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage160_bare_mixed_port_slippage_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage160_bare_mixed_port_slippage_mathematica_audit.txt`

## What the paper claims

The stage card quotes a single boxed result:
> `\delta\gamma_W=(\delta\gamma_0-\delta\kappa_0/3)/(1+\mathfrak r_{F1}^2)`.

The supporting notes derive this via three substantive steps that together constitute the stage's claim: (1) the exact compensated-branch identity `\delta\gamma_W - (1/3)\delta\kappa_W = (\delta\gamma_0 - (1/3)\delta\kappa_0)/(1+r_{c,*})` from linearizing `\kappa_W=\kappa_0/(1+r_c)` and `\gamma_W=\gamma_0/(1+r_c)` around the canonical branch `\kappa_0_*=(1+r_{c,*})/3, \gamma_0_*=(1+r_{c,*})/9`; (2) collapse to the boxed card identity under the Stage 244 canonical-even gate `\delta\kappa_W=0`; (3) pure-scale harmlessness, `\delta\gamma_0=\delta\kappa_0/3 \Rightarrow \delta\gamma_W=0`. Sections 4-5 of the notes chain this into a final defect law (`\Delta_Q`, `N_Q-1`) using the Stage 159 transport `\delta\Pi_{tan} = 0.832409\,\delta\Sigma_0 - 1.16276\,\delta\mathcal S` and the susceptibility `\Upsilon_\Pi`; the card itself does not box these chained forms, so they are downstream narrative rather than core stage claims.

## What the script claims to verify

Both scripts assert exactly two zero-residuals: (A1) the exact compensated-branch identity `dγ_W - (1/3)dκ_W - (dγ_0 - (1/3)dκ_0)/(1+r_c) == 0`, derived in-script by series expansion of `(k0_star + ε dk0)/(1+rc+ε drc)` and `(g0_star + ε dg0)/(1+rc+ε drc)` to first order in ε; (A2) pure-scale harmlessness, `[(dγ_0 - (1/3)dκ_0)/(1+r_c)] |_{dγ_0=dκ_0/3} == 0`. Section 2 (canonical-even gate) is rendered as a `print` of the form `dγ_W` collapses to, not an assertion (this is appropriate, since the gate is upstream Stage 159's result, not Stage 160's). The "tangential DtN susceptibility and final defect law" block is print-only — it composes symbolic expressions and displays them without asserting any equality. The carry-forward formulas block is pure print narrative.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Exact compensated-branch identity `dγ_W - (1/3)dκ_W = (dγ_0 - (1/3)dκ_0)/(1+r_{c,*})` (Section 1, notes) | A1 (sympy:50-51; math:41-42) | match |
| Boxed card identity `dγ_W = (dγ_0 - dκ_0/3)/(1+r_F1²)` under `dκ_W=0` gate (card + Section 2) | sympy:54-55, math:44-45 (print-only) | partial (derived by substitution from A1; not separately asserted but algebraically forced by A1 + the externally-supplied gate dκ_W=0) |
| Pure-scale harmlessness (Section 3) | A2 (sympy:58-59; math:46) | match |
| Final defect law (Sections 4-5, notes, not boxed in card) | sympy:63-74, math:50-61 (print-only) | extra-scaffolding (paper card does not require it; downstream-narrative only) |

Dominant pattern: aligned. The card's single boxed deliverable is structurally forced by A1 plus the upstream `dκ_W=0` gate; A1 is asserted directly, A2 covers the pure-scale corollary.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 50-51 | `expect_zero(... dgW - (1/3)dkW - (dg0 - (1/3)dk0)/(1+rc))` | Section-1 exact identity (and by gate substitution, card identity) | yes |
| A2 | sympy | 58-59 | `expect_zero("pure-scale harmlessness", dgW_purescale)` | Section-3 pure-scale | yes |
| B1 | mathematica | 41-42 | `expectZero["exact compensated-branch slippage identity", identity]` | Section-1 exact identity | yes |
| B2 | mathematica | 46 | `expectZero["pure-scale harmlessness", dgWGate /. dg0 -> dk0/3]` | Section-3 pure-scale | yes |

Both A1/B1 are non-tautological: the residual is a non-trivial combination of four independent symbols (rc, drc, dk0, dg0) and the linearization step has real algebraic content (the `dr_c` dependence must cancel between `dκ_W` and `(1/3)dγ_W` — verified by the script outputs `dκ_W = (-dr_c/3 + dκ0)/(rc+1)` and `dγ_W = (-dr_c/9 + dγ0)/(rc+1)`, and the cross-cancellation is enforced precisely because `κ_0_*/3 = γ_0_*/1` here gives `(1+rc)/9` for both `γ_0_*` and `κ_0_*/3 = (1+rc)/9` — i.e., A1 is the algebraic content of choosing those specific star values). A2 is the corollary substitution, asserted; passes only because A1 + substitution work.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage160_bare_mixed_port_slippage_sympy_audit.py:34-59`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage160_bare_mixed_port_slippage_mathematica_audit.wl:28-46`

**What's wrong:**
The Mathematica script is a structural line-by-line port of the SymPy script. Corresponding excerpts:

SymPy (lines 34-46):
```
eps = sp.symbols("eps", real=True)
rc, drc = sp.symbols("r_c dr_c", real=True)
dk0, dg0 = sp.symbols("dκ0 dγ0", real=True)
k0_star = (1 + rc) / 3
g0_star = (1 + rc) / 9
kW = ((k0_star + eps * dk0) / (1 + rc + eps * drc)).series(eps, 0, 2).removeO()
gW = ((g0_star + eps * dg0) / (1 + rc + eps * drc)).series(eps, 0, 2).removeO()
dkW = sp.expand(kW.coeff(eps, 1))
dgW = sp.expand(gW.coeff(eps, 1))
```

Mathematica (lines 28-36):
```
Clear[eps, rc, drc, dk0, dg0];
$Assumptions = Element[{eps, rc, drc, dk0, dg0}, Reals];
k0Star = (1 + rc)/3;
g0Star = (1 + rc)/9;
kW = Normal[Series[(k0Star + eps*dk0)/(1 + rc + eps*drc), {eps, 0, 1}]];
gW = Normal[Series[(g0Star + eps*dg0)/(1 + rc + eps*drc), {eps, 0, 1}]];
dkW = Expand[Coefficient[kW, eps, 1]];
dgW = Expand[Coefficient[gW, eps, 1]];
```

Same symbol names (`eps, rc, drc, dk0, dg0`), same intermediate names (`k0Star/k0_star`, `g0Star/g0_star`, `kW`, `gW`, `dkW`, `dgW`), same algorithmic shape (`Series` to first order then coefficient extraction), same identity assembly (`dgW - 1/3*dkW - (dg0 - 1/3*dk0)/(1 + rc)`), same gate substitution (`dgWGate /. dg0 -> dk0/3` vs `dgW_gate.subs(dg0, dk0 / 3)`), and the tangential-block hardcoded numerics are byte-identical (`0.832409471081635` and `1.16275838754222`). This violates the second-engine policy that both engines must derive the claim independently from the physical premises.

**Why this matters:**
The cross-engine sanity check exists to catch a systematic error that one engine's simplification machinery might silently introduce (e.g., branch-cut convention, series convention). A line-by-line port re-executes the same algebra in a different toolchain, which only confirms that both engines can faithfully simulate each other on the same algorithm — not that two independent derivations agree on the physics.

**Required change:**
Reshape the Mathematica derivation so it follows a structurally different path to the same identity. A natural independent route: instead of `Series` + `Coefficient`, write the linearized perturbations symbolically as

```
kWFull = (k0Star + dk0)/(1 + rc + drc);
gWFull = (g0Star + dg0)/(1 + rc + drc);
dkWind = Normal[Series[kWFull, {dk0, 0, 1}, {drc, 0, 1}]] - k0Star/(1+rc);
```

then form the residual `dgWind - (1/3) dkWind - (dg0 - (1/3) dk0)/(1+rc)` and `FullSimplify` it to zero. Or, equivalently, derive `dgW` and `dkW` analytically from `D[...]/.{dk0->...}` (treating the variation as a differential), keeping the algebra distinct from SymPy's `series(...).coeff(eps, 1)` recipe. Symbol names should also drift away from the SymPy script's exact spellings (e.g., `dKappa0`, `dGamma0`, `deltaR`, `kappaWperp` etc.) so a future audit can tell at a glance that the two scripts were authored independently. Codex: do not touch the SymPy script; only modify the `.wl` file.

**Verification:**
The new Mathematica script should (a) still pass `expectZero` for the compensated-branch identity and pure-scale harmlessness; (b) not contain the literal substrings `k0Star`, `g0Star`, `Series[..., {eps, 0, 1}]`, `Coefficient[..., eps, 1]`; (c) reach the same residual `0` via a `D[...]`-based or direct-symbolic-expansion route. Verifier confirms by diffing the structural shape of the `.wl` derivation against the `.py` and by re-running `redteam exec-mathematica 160` to confirm exit 0.

## Independent-derivation check (Mathematica)

See F1. The Mathematica script is a transliteration: identical variable choreography, identical `Series` + `Coefficient` recipe, identical substitution chain, identical hardcoded numerics in the print-only tangential block.

## Engine cross-check

Both engines produce the same first-order pieces:
- SymPy: `dκ_W = (-dr_c/3 + dκ0)/(r_c + 1)`, `dγ_W = (-dr_c/9 + dγ0)/(r_c + 1)`
- Mathematica: `dκ_W = (3*dk0 - drc)/(3 + 3*rc)`, `dγ_W = (9*dg0 - drc)/(9 + 9*rc)`

These are algebraically identical (factor 1/3 vs. 3/(3+3rc) is just rewritten). Both `expectZero` pass on the identity and pure-scale residuals. The hardcoded transport coefficients in the print-only block are byte-identical between scripts. Engines agree — but per F1, agreement on a transliterated algorithm carries less weight than agreement reached by two independent derivations.

## Verdict justification

The math holds: the script's two assertions (exact compensated-branch identity and pure-scale harmlessness) are non-tautological, traceable to the paper card's boxed identity (modulo the upstream `dκ_W=0` gate), and produce the correct residual `0`. The card's boxed deliverable is exercised. The single finding is structural: the Mathematica script does not provide an independent second engine — it re-runs SymPy's algorithm through Mathematica's syntax. No `paper_misalignment`, no `tautological_check`, no `symbol_assumption_error`, no `missing_branch`, no `engine_disagreement` (both produce the same identity), no `hardcoded_result` in the assertions (the print-only numerics in the tangential block are not load-bearing for any assertion and match the notes' Stage-244-derived values). Outputs are fresh (sympy mtime 12:47, script mtime 11:58; mathematica output 13:17, script 11:56). Banner labels in both scripts say "STAGE 143" rather than "STAGE 160" — a stale-stage-number cosmetic artifact in a `Print` line, not load-bearing for any check, and the notes header similarly says "Stage 245"; this is a renumbering inheritance across the codebase and not flagged as a script defect under the audit categories. Verdict: `findings` (one medium-severity `mathematica_transliteration`).

## Self-test notes

(1) Variable independence: no `sp.diff`/`D[]` in the script; the linearization is via `Series` + `Coefficient`, which by construction extracts the coefficient of `eps` from a known expansion — not subject to the silent-zero-derivative trap. (2) Symmetry/parity: no integrals on unbounded domains in this stage. (3) Trivial-case pre-check: substituting `dk0=0, dg0=0, drc=0` makes both `dkW` and `dgW` vanish, and the identity reduces to `0 - 0 - 0 = 0` — fine; with `dk0=0, dg0=1, drc=0` the identity becomes `1/(1+rc) - (1/3)·0 - (1 - 0)/(1+rc) = 0` — non-trivially correct, confirms the assertion can distinguish wrong star-values. (4) Path specification: F1 modifies the existing `.wl` file at the given absolute path; no missing-script subtype. (5) Paper round-trip: F1's required change keeps the same paper-side identity; only the derivation mechanism inside Mathematica changes; no new `paper_misalignment` introduced.
