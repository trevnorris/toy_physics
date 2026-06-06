---
unit_id: 106
batch: IV.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage106_canonical_outgoing_reduced_closure.md]
  paper_appendix: present
---

# Audit unit 106 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_106.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage106_canonical_outgoing_reduced_closure.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows: compression chain l.44–62; theorem l.65–94; one-scalar-defect l.215–258; retarded factorization l.260–295; canonical DtN branch l.297–343; invariant-coeff box l.328–338; result-anchor row l.1176)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.txt`

## What the paper claims

Stage 106 is the reduced 2.5PN closure on the canonical outgoing DtN branch. The card's Derivation ledger states the load-bearing claim: "The computation isolates the reduced product \(\widehat m_0^{\,2}\chi_Q N_Q=1\) and the canonical condition \(\chi_Q=1\)." `Inputs` declares \(\chi_Q\) (the outgoing \(l=2\) DtN coefficient), \(\widehat m_0\), and \(N_Q\) as **imports** — i.e. \(\chi_Q=1\) is a carried-in anchored value, not a quantity 106 re-derives. The notes make the deliverables explicit: imposing the carried \(\chi_Q=1\) into the frozen factorization \(\widehat m_0^{\,2}\chi_Q N_Q=1\), in the strict point-particle limit \(\widehat m_0=1\), yields the boxed \(\boxed{N_Q=1}\); the canonical invariant coefficients then fix to \(\overline K_0=54Gc_s^5/(5a^5c^5)\), \(\overline K_2=6Gc_s^3/(5a^3c^5)\), \(\overline K_4=8Gc_s/(15ac^5)\), \(\overline\Gamma_5=2G/(5c^5)\), and the effective odd coefficient \(\gamma_{\rm quad}^{\rm eff}=\widehat m_0^{\,2}\Gamma_5=2G/(5c^5)\). The appendix corroborates every literal (eq. `app-part04-canonical-invariant-coeffs` l.328–338, eq. `app-part04-factorized-defect-again` l.286–288, eq. `app-part04-sigmaQcan` l.275). The card's `Checks` further name (ii) higher odd terms beginning beyond the 2.5PN \(\omega^5\) coefficient and (iii) the outgoing \(l=2\) DtN fingerprint vs. the normalized \(z=\omega a/c_s\) expansion; the SymPy docstring (l.4–13) explicitly delegates (ii)→Stage 102, (iii)→Stage 104, and the \(\chi_Q=1\) fix→Stage 105, with 106 owning item (i) — factorization separability and the closure to \(N_Q=1\).

## What the script claims to verify

Both engines verify: (a) solving \(\widehat m_0^{\,2}\chi_Q N_Q=1\) gives \(N_Q=1/(\chi_Q\widehat m_0^{\,2})\), collapsing to \(N_Q=1\) at \(\chi_Q=1,\widehat m_0=1\); (b) the four carried target literals \(\overline K_0,\overline K_2,\overline K_4,\overline\Gamma_5\) satisfy two mutual algebraic identities, \(K_0 K_4-4K_2^2=0\) and \(\Gamma_5-9\sqrt{K_2^5/K_0^3}=0\) (these cross-relate the independent literals, so a typo in any one would trip them); (c) the canonical \(\gamma_{\rm eff}=\widehat m_0^{\,2}N_Q\Gamma_5\) equals the target \(2G/(5c^5)\); (d) a first-order \(\Delta_Q\) sensitivity slope of \(-2G/(5c^5)\) around \(\chi_Q=1+\Delta_Q\). The Mathematica engine ADDS an independent route covering Check (iii)-adjacent structure: it builds the retarded one-pole transfer function \(\widehat Y_Q^{\rm ret}(\omega)\) (appendix eq. `retarded-one-pole`), series-expands to \(\omega^7\), checks \(\sigma_Q^{\rm can}=4a^5/(27c_s^5)\), confirms the next odd term begins at \(\omega^7\) (Check (ii)), and re-fixes \(\chi_Q=1\) by matching the \(\omega^5\) amplitude to \(\overline\Gamma_5^{\rm target}\) (appendix eq. `chiQ-equals-one` logic).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| \(\widehat m_0^{\,2}\chi_Q N_Q=1\) separable factorization (card ledger) | py l.40–44 constraint + solve; wl l.98 `sourceNormalizer` | match |
| \(N_Q=1\) in strict point-particle limit (notes box; card) | py l.46–48 `expect_zero`; wl l.99–102 | match |
| \(\chi_Q=1\) carried-in / canonical condition (card Inputs+ledger) | py: consumed as `.subs(chi_Q,1)` carry-in (l.46,74,86); wl l.80–81 re-derives via \(\omega^5\) match (appendix-faithful) | match |
| \(\overline K_0,\overline K_2,\overline K_4,\overline\Gamma_5\) literals (notes box; appendix l.328–338) | py l.50–53 + identity asserts l.65–72; wl l.70–73 + l.86–93 | match |
| \(\gamma_{\rm quad}^{\rm eff}=2G/(5c^5)\) (notes box) | py l.74–78; wl l.104–107 | match |
| Check (ii) higher odd terms begin beyond \(\omega^5\) | delegated to Stage 102 (py docstring); wl l.55–67 also confirms next odd at \(\omega^7\) | match (wl exceeds) |
| Check (iii) outgoing \(l=2\) DtN fingerprint vs \(z\)-expansion | delegated to Stage 104/105 (py docstring); wl l.45–67 reconstructs the one-pole \(\omega^5\) amplitude route | match (delegated; wl partial-extra) |
| \(\sigma_Q^{\rm can}=4a^5/(27c_s^5)\) (appendix eq. sigmaQcan) | wl l.39–41 only | match (single-engine, by design) |

`paper_alignment: aligned`. Every boxed/deliverable value the scripts emit is present and agreeing in the notes and/or appendix. The delegation of Checks (ii)/(iii) to upstream stages is stated openly in both scripts' headers and matches the card's `Inputs` import structure; the Mathematica engine voluntarily re-exercises both, which strengthens rather than weakens coverage.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 48 | `expect_zero(NQ_canonical.subs(m0hat,1)-1)` | N_Q=1 closure | yes |
| A2 | sympy | 65–68 | `K0_t*K4_t - 4*K2_t^2 == 0` | even-coeff literals mutual consistency | yes (cross-relates 3 literals) |
| A3 | sympy | 69–72 | `Gamma5_t - 9*sqrt(K2_t^5/K0_t^3) == 0` | odd↔even literal consistency | yes |
| A4 | sympy | 74–78 | `gamma_eff_canonical - Gamma5_target == 0` | gamma_eff=2G/(5c^5) | yes |
| A5 | sympy | 89–92 | `linear_coeff - (-2G/(5c^5)) == 0` | Delta_Q slope (sensitivity) | yes (non-trivial series coeff) |
| A6 | sympy | 93–97 | `zeroth_coeff - Gamma5_target == 0` | Delta_Q zeroth = target | yes |
| B1 | math | 41 | `sigmaQcan - 4a^5/(27 cs^5) == 0` | sigma_Q^can literal (appendix) | yes |
| B2 | math | 63 | `omega5Coeff - I*chiQ*sigmaQcan/4 == 0` | one-pole omega^5 form | yes (series-derived) |
| B3 | math | 64–67 | `omega7Coeff - I*sigmaQcan/(2 OmegaQ^2) == 0` | next odd at omega^7 (Check ii) | yes |
| B4 | math | 76–79 | `oddScaleFromSeries - chiQ*gamma5Target == 0` | omega^5 amplitude = chi_Q*Gamma5_t | yes |
| B5 | math | 80–81 | `chiFromOmega5Match - 1 == 0` | chi_Q=1 fixed by canonical match | yes (see note in §Independent-derivation) |
| B6 | math | 86–93 | two target-literal identities | even/odd literal consistency | yes |
| B7 | math | 99–102 | `(sourceNormalizer at m0=1,chi=1)-1 == 0` | N_Q=1 closure | yes |
| B8 | math | 107 | `gammaEffCanonical - gamma5Target == 0` | gamma_eff=2G/(5c^5) | yes |
| B9 | math | 116–124 | Delta_Q slope and zeroth | Delta_Q sensitivity | yes |

No row is tautological-by-construction: A2/A3/B6 cross-relate three independent hardcoded literals (a single mistyped literal trips them); A1/A4/B7/B8 substitute into a solved/constructed expression and would change under a sign or factor error in the factorization; A5/B9 extract a genuine series coefficient (a sign flip in \(N_Q=1/\chi_Q\) flips the slope). B1–B5 derive from the transfer-function series, a route absent from SymPy.

## Findings

None. (Zero script-side findings; zero `paper_misalignment`.)

## Independent-derivation check (Mathematica)

The `.wl` is a **genuinely independent** route, not a transliteration. Decisive evidence:

1. **Different core mechanism.** SymPy never constructs a transfer function or performs any series expansion; its closure is a one-line algebraic `solve(m0hat^2*chi_Q*N_Q=1, N_Q)` (py l.40–43). Mathematica instead *builds* the retarded one-pole module `Yret[om,chi] = 3/4 + (1/4)/(1 - om^2/OmegaQ^2 - I*chi*sigmaQcan*om^5)` (wl l.45) and runs `Normal[Series[Yret,{omegaSym,0,7}]]` (wl l.49), extracting \(\omega^5/\omega^6/\omega^7\) coefficients (wl l.55–60). There is no corresponding `series()`/`Series[]` call, no `Yret`, and no `sigmaQcan` anywhere in the SymPy file.

2. **Mathematica re-derives chi_Q from first structure; SymPy carries it in.** SymPy treats \(\chi_Q=1\) purely as a substitution (`.subs(chi_Q,1)`, py l.46). Mathematica derives it: `oddScaleFromSeries = k0Target*omega5Coeff/I`, then `Solve[oddScaleFromSeries == gamma5Target, chiQ]` → `chiQ -> 1` (wl l.75–81). Opposite handling of the same quantity is the antithesis of a line-by-line port.

3. **Shared portions are only the trivial substitution closure and the two target-literal identities** (py l.65–72 ↔ wl l.86–93; py l.48 ↔ wl l.99–102; py l.85–97 ↔ wl l.109–124). These are short, are literal-consistency / point-substitution checks dictated by the shared appendix box, and constitute a minority of the Mathematica script; they do not make the engine an echo of the SymPy algebra.

One observation (not a finding): the Mathematica header comment (wl l.30–32) says "This engine uses chi_Q = 1 as carry-in," yet the body at l.80–81 actively re-derives \(\chi_Q=1\) from the \(\omega^5\) match. The re-derivation is sound and faithful to appendix eq. `app-part04-chiQ-equals-one` ("Comparison with the retarded one-pole fixes \(\chi_Q=1\)"), it is NOT a tautological re-pin of a free parameter (it cross-constrains the three independent literals `k0Target`, `sigmaQcan`, `gamma5Target` — a wrong literal would solve to \(\chi_Q\ne1\)), and it does not contradict the card, which lists \(\chi_Q\) as an import while also naming "the canonical condition \(\chi_Q=1\)." The comment is merely slightly looser than the code. This is a cosmetic prose nit inside a script comment, below the threshold for a finding; I leave it rather than manufacture a finding.

## Engine cross-check

Both saved outputs are present and fresh (outputs newer than scripts; see Self-test notes). Every shared deliverable agrees:

- N_Q closure: sympy `point-particle canonical branch gives N_Q = 1 = 0` (txt l.8) ↔ math `N_Q on natural branch at m0hat=1, chi_Q=1 = 0 / PASS` (txt l.24–25).
- Target-literal identities: sympy txt l.13–14 ↔ math txt l.20–23, both `= 0`.
- gamma_eff: sympy `canonical gamma_eff - target = 0` (txt l.15) ↔ math (txt l.26–27).
- Delta_Q slope and zeroth: sympy txt l.16–17 ↔ math txt l.28–31, both `= 0`.
- Mathematica-only series data is internally consistent: `omega^5 coefficient (general) = I/27 a^5 chiQ/cS^5`, and `chi_Q fixed by canonical omega^5 match = 0 / PASS` (txt l.9, l.18–19), agreeing with the SymPy carry-in \(\chi_Q=1\).

No residual sign, factor-of-two, or magnitude disagreement. Engines agree.

## Verdict justification

`clean`. I attacked the four standard vectors and each held: (1) **Tautology** — the two target-literal identities cross-relate three independently typed constants (a single literal typo trips them), and the \(N_Q\)/gamma_eff/Delta_Q checks substitute into solved/constructed expressions rather than into a value defined to equal itself; none is guaranteed by construction. (2) **Hardcoded-not-derived** — the four \(\overline K\)/\(\Gamma_5\) literals and `gamma5Target` are carried-in appendix-box values (legitimate carry-ins for a ledger/closure stage), and the appendix derives them (eq. l.328–338); they are anchored on the paper side and mutually constrained on the script side. (3) **Single-point / insufficient** — the \(\Delta_Q\) series check and the Mathematica \(\omega^7\)-order expansion exercise the algebra beyond the static point, and the \(\chi_Q=1\) re-derivation is amplitude-matched, not asserted. (4) **Transliteration** — the Mathematica engine uses a structurally distinct transfer-function/series route absent from SymPy and re-derives \(\chi_Q\) where SymPy carries it in; not a port. Paper ↔ script alignment is exact: every boxed deliverable maps to a faithful check, the delegation of Checks (ii)/(iii) to upstream stages is openly declared and matches the card's `Inputs` import structure, and 106 correctly **consumes** \(\chi_Q=1\) as an upstream anchored input rather than silently re-pinning a free parameter. I confirm I read the card, the notes, and the Part IV appendix rows before judging the scripts.

## Value Reconciliation (pass-2 augmentation)

Deliverable values emitted by the scripts/outputs, located in the card (`.tex`), notes (`.md`), and Part IV appendix:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| \(N_Q=1\) | py l.48 / txt l.8; wl l.99–102 / txt l.24 | notes box l.35 `\boxed{N_Q=1}`; card l.13 (canonical condition); appendix l.249,288 | MATCH |
| \(\widehat m_0^{\,2}\chi_Q N_Q=1\) | py l.40 / txt l.5; wl l.98 | card l.13; notes l.16; appendix eq. `factorized-defect-again` l.288 | MATCH |
| \(\chi_Q=1\) (consumed/re-confirmed) | py `.subs(chi_Q,1)` l.46; wl l.80–81 / txt l.18 | card l.9,13; notes l.8,31,73; appendix eq. `chiQ-equals-one` l.326 | MATCH |
| \(\overline K_0=54Gc_s^5/(5a^5c^5)\) | py l.50 / txt l.9; wl l.70 | notes box l.43; appendix l.331 | MATCH |
| \(\overline K_2=6Gc_s^3/(5a^3c^5)\) | py l.51 / txt l.10; wl l.71 | notes box l.44; appendix l.333 | MATCH |
| \(\overline K_4=8Gc_s/(15ac^5)\) | py l.52 / txt l.11; wl l.72 | notes box l.45; appendix l.335 | MATCH |
| \(\overline\Gamma_5=2G/(5c^5)\) | py l.53 / txt l.12; wl l.73 | notes box l.49; appendix l.337 | MATCH |
| \(\gamma_{\rm quad}^{\rm eff}=2G/(5c^5)\) | py l.74–78 / txt l.15; wl l.104–107 / txt l.26 | notes box l.55–61; appendix l.286,337 | MATCH |
| \(\sigma_Q^{\rm can}=4a^5/(27c_s^5)\) | wl l.39–41 / txt l.5 (Mathematica only) | appendix eq. `sigmaQcan` l.275 | MATCH |
| \(\omega^5\) coeff \(=\tfrac{i}{27}a^5\chi_Q/c_s^5\), next odd at \(\omega^7\) | wl l.55–67 / txt l.7–11 | appendix eq. `retarded-one-pole` l.264–270, `Yout-dtn` l.321 (\(z^5/27\)) | MATCH |

INTERNAL scaffolding (no finding): `NQ_general = 1/(chi_Q*m0hat^2)` (general pre-substitution form), printed `K0/K2/K4/Gamma5` general expressions carrying the `1/(chi_Q*m0hat^2)` factor (txt l.9–12; intermediate, the canonical/box values are the deliverables), `\Delta_Q` first-order slope \(-2G/(5c^5)\) and zeroth coeff (sensitivity probe, not a stage deliverable — \(\Delta_Q:=\chi_Q-1\) itself is the appendix object, its slope is internal), `oddScaleFromSeries`, `chiFromOmega5Match`, pass/fail flags.

`reconciliation: complete; 11 deliverable values checked, 0 misaligned`

## Self-test notes

Variable-independence: the only `series()`/`Series[]` are over \(\Delta_Q\) (py l.87, on `gamma_eff_off` which genuinely depends on \(\Delta_Q\)) and over `omegaSym` (wl l.49, on `Yret` which genuinely depends on `omegaSym`); neither is a degenerate zero-derivative. Symmetry/parity: no unbounded-domain integrals in either script. Trivial-case pre-check: substituting \(\chi_Q=1,\widehat m_0=1\) into `N_Q=1/(chi_Q m0hat^2)` gives 1 (assert_zero residual 0 ✓), and the \(\Delta_Q\) linear coefficient of \(1/(1+\Delta_Q)\) is \(-1\Rightarrow -\)Gamma5_target (assert_zero residual 0 ✓, nonzero literal as required). Paper round-trip: the four \(\overline K\)/\(\Gamma_5\) literals and `gamma5Target` used in-script byte-match the appendix box (l.328–338) and notes box (l.42–61); no fix prescribed, so no risk of introducing a new misalignment. Output freshness confirmed by mtime (both `.txt` newer than their scripts).
</content>
</invoke>
