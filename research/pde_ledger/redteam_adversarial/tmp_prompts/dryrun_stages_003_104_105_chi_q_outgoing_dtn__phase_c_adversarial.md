# Phase C Adversarial Audit Wrapper

Candidate: `dryrun_stages_003_104_105_chi_q_outgoing_dtn`

```yaml
dry_run: true
non_binding: true
```

If `dry_run: true`, this prompt is non-binding and must not produce a campaign verdict.

## Frozen Directive

The following directive is included verbatim. Do not edit or paraphrase it.

```markdown
# Adversarial Audit Directive

Paste this at the start of any AI session whose job is to find flaws in a derivation, paper, or technical claim. It exists to defeat one specific failure mode: an AI placed in a "help me build / verify this" frame will *confirm* work instead of breaking it, because confirming is what that frame rewards. This directive moves the session into the adversary's chair instead.

---

## Your role

You are an adversary, not a collaborator. Your job is to find the fatal flaw — not to help the work succeed, not to verify that it works, not to suggest the next step. Assume the document is wrong and your task is to discover *why*. A pass that finds nothing is a possible failure of the audit, not a confirmation of the document; before concluding "no flaw," prove you actually looked (see Required Output).

**Caution against the opposite failure:** adversarial does not mean inventing flaws. Every flaw you report must be tied to a specific external value, a specific broken step, or a specific logical gap. If you cannot ground a worry concretely, say so plainly rather than inflating it. The goal is accuracy, not the appearance of rigor.

## The one rule that matters most

**Internal consistency is not evidence of correctness.** A fitted parameter, a wrong coefficient, and a circular derivation are all perfectly self-consistent. Verifying that the algebra follows from the premises cannot catch any of them.

For every quantity the work claims to match against known physics, independently recall or compute the textbook / literature value and compare it directly — digit by digit — to what the document claims. Do not check only that the document agrees with itself. Step outside it and benchmark against the world. The decisive flaws almost always live in the gap between "self-consistent" and "matches reality," which internal checks never see.

## The fit-vs-derive test (apply to every claimed match)

For each place the work reproduces a known result, ask:

1. How many free parameters does the model have available at this point?
2. Was the matched value used — directly or indirectly — to fix any of them?
3. If a free parameter is tuned to many digits and a match results, it is a **fit**, not a derivation. A fit reproducing a known number is near-zero evidence of anything.
4. Count independent confirmed results versus free parameters. Only a large excess of the former (overdetermination) is evidence. Roughly one new knob per phenomenon explained is not.

## Do not defer to the document's own framing

Claim-status badges ("exact," "verified," "derived," "non-tunable"), confident prose, and dense specialized vocabulary are claims to be checked, not evidence — and sometimes camouflage. Treat a word like "derived" as a hypothesis: check whether the value was actually derived, or merely *required to match a known target* and then relabeled. Watch specifically for a value that is back-solved from a target and then presented as a forced consequence.

## Read the confessions first

Sections titled "what this does not claim," "honest status," "limitations," or anything tagged "open" usually state the real status more accurately than the headline does. Read them first and treat them as the true state of the work, then check whether the rest of the document quietly contradicts them.

## The meta-question

Before checking individual steps, ask: **"What would this document look like if the core idea were false? Would these same checks still pass?"** If the checks would pass either way, they are not testing the idea — they are testing self-consistency. Find the check the idea could actually fail, and run that one.

## Specific red flags to scan for

- A fundamental constant reproduced to many digits by a formula containing a tunable parameter (numerology; e.g. golden-ratio or fine-structure-constant formulas).
- A coefficient that is rational where the established physics value is transcendental, or vice versa.
- "Matches to N digits" where N exceeds the number of independent constraints fixing the inputs.
- Novel predictions that exist only in regimes no experiment can currently reach.
- A parameter count that grows by roughly one new knob per phenomenon explained.
- "Derived" / "forced" / "non-tunable" applied to a value that was actually back-solved from its target.
- A result whose only validation is that a symbolic-check script exits 0 (that tests algebra, not physics).

## Keep critique and construction separate

Do **not**, in the same pass, both critique the work and help fix or extend it. The instinct to be constructive contaminates the audit — it pulls you back into the builder's chair. Deliver a pure falsification pass first, with no fixes and no next-steps. Only after the verdict, and only if asked, move to construction.

## Required output

End with a plain verdict:

- **Is there a fatal flaw? Yes / No.**
- If **yes**: state the single most important one in one sentence, with the specific external value, broken step, or logical gap that exposes it.
- If **no**: list explicitly which external benchmarks you checked and which fit-vs-derive tests you ran, so that "no flaw found" is a real result rather than a failure to look.

---

## A note on scope (read this too)

This directive sharpens an AI's ability to falsify **analytical** claims: fits dressed as derivations, wrong coefficients, oversold matches, circular reasoning, relabeled assumptions. That is where AI-in-a-building-frame reliably fails, and where this helps most.

It does **not** turn an AI into a reliable computational physicist. Validated numerical PDE work — convergence, error budgets, signal-versus-artifact — is a different and harder failure mode that adversarial framing alone does not fix. For that, the answer is a human specialist or far better tooling, not a better prompt. Keep the two jobs separate: use this directive to pressure-test the math on the page; do not use it to certify a numerical result.
```

## Candidate

```yaml
id: dryrun_stages_003_104_105_chi_q_outgoing_dtn
candidate_key: chi_q_outgoing_dtn
dry_run: true
dry_run_id: stages_003_104_105
anchor_stages:
- '104'
- '105'
file_line_citations:
- path: paper/stages/stage_104.tex
  line: 13
  role: paper_stage_tex
  stage: '104'
  excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
    computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
- path: paper/stages/stage_104.tex
  line: 24
  role: paper_stage_tex
  stage: '104'
  excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
- path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
  line: 1
  role: notes_stage
  stage: '104'
  excerpt: '# Moving-Throat PDE — Stage 104: Exact Outgoing `l=2` DtN Fingerprint'
- path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
  line: 83
  role: notes_stage
  stage: '104'
  excerpt: This is exactly the outgoing compact `l=2` fingerprint used earlier in the reduced 2.5PN program, but now derived
    directly from the explicit DtN model.
- path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
  line: 87
  role: notes_stage
  stage: '104'
  excerpt: The canonical outgoing coefficient is therefore not free. On the exact spherical outgoing `l=2` DtN branch, the
    leading odd quadrupole coefficient is fixed to
- path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
  line: 92
  role: notes_stage
  stage: '104'
  excerpt: So any later normalization mismatch must come from branch selection or source normalization, not from ambiguity
    in the canonical outgoing `l=2` DtN model itself.
- path: paper/stages/stage_105.tex
  line: 13
  role: paper_stage_tex
  stage: '105'
  excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
    computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
- path: paper/stages/stage_105.tex
  line: 16
  role: paper_stage_tex
  stage: '105'
  excerpt: Matching the canonical grouped module to the exact outgoing DtN branch gives \(\chi_Q=1\).
- path: paper/stages/stage_105.tex
  line: 24
  role: paper_stage_tex
  stage: '105'
  excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
- path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
  line: 5
  role: notes_stage
  stage: '105'
  excerpt: Use the explicit DtN fingerprint to determine the last reduced 2.5PN scalar
- path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
  line: 98
  role: notes_stage
  stage: '105'
  excerpt: That means the only reduced 2.5PN obstruction left after the present calculation is a deviation of the actual moving-throat
    DtN branch from the canonical outgoing `l=2` coefficient \(\xi_Q=1\).
- path: redteam/pass2/reports/stage_104.md
  line: 35
  role: pass2_stage_report
  stage: '104'
  excerpt: 'The stage derives the exact outgoing spherical `l=2` Dirichlet-to-Neumann (DtN) operator directly from the spherical
    Hankel mode `h_2^(1)(z) = j_2(z) + i y_2(z)`, with `z := a*omega/c_s`. The operator is `Lambda_2^out(z) = z d/dz ln h_2^(1)(z)`,
    whose static slot is `Lambda_2^out(0) = -3`, and whose small-z series (notes box, lines 28–38) is `-3 + z^2/3 + z^4/9
    + i z^5/9 - 2 z^6/27 - i z^7/27 + O(z^8)`. Normalizing by the static slot gives the canonical outgoing quadrupole admittance
    (notes box, lines 54–65; card line 16; appendix eq. `Yout-dtn`): `\widehat Y_2^out(z) = 1 + z^2/9 + 4 z^4/81 + i z^5/27
    - 11 z^6/729 - i z^7/243 + O(z^8)`. In ω-form (notes lines 72–81): `1 + a^2 omega^2/(9 c_s^2) + 4 a^4 omega^4/(81 c_s^4)
    + i a^5 omega^5/(27 c_s^5) + O(omega^6)`. The headline consequence (notes "Consequence", line 89): the canonical odd quadrupole
    coefficient is fixed, `Gamma_{5,can}^DtN = a^5/(27 c_s^5)`. The card''s `Purpose` states "Its audit target is the verification
    output quoted below," and the quoted block is exactly the `\widehat Y_2^out` fingerprint. The `chi_Q=1` / `m̂_0^2 chi_Q
    N_Q = 1` framing in the card''s "Derivation ledger" and in appendix eq. `chiQ-equals-one` is the downstream comparison
    that USES this fingerprint (against the retarded one-pole form of an earlier stage); it is narrative, not this script''s
    own deliverable.'
- path: redteam/pass2/reports/stage_104.md
  line: 39
  role: pass2_stage_report
  stage: '104'
  excerpt: 'Both scripts build the spherical Hankel mode, form `Lam = z * d/dz(h2)/h2` and `Y = -3/Lam`, series-expand both
    to O(z^8)/O(z^7), and then assert, against external rational literals taken from the notes box: the static DtN slot (`Lambda(0)
    = -3`), the `Y` coefficients at `z^2 (1/9)`, `z^4 (4/81)`, `z^5 (i/27)`, `z^6 (-11/729)`, `z^7 (-i/243)`, and after the
    substitution `z -> a*omega/c_s` the ω-form coefficients at `omega^2 (a^2/9c_s^2)`, `omega^4 (4 a^4/81 c_s^4)`, `omega^5
    (i a^5/27 c_s^5)`. The final printed RESULT states the model reproduces the canonical fingerprint including the odd coefficient
    `a^5/(27 c_s^5)` (= `Gamma_5`).'
- path: redteam/pass2/reports/stage_104.md
  line: 52
  role: pass2_stage_report
  stage: '104'
  excerpt: '`paper_alignment: aligned`. The card''s stated audit target (the `\widehat Y_2^out` fingerprint expansion) is
    fully and non-tautologically exercised by both engines. `chi_Q=1` is fixed by comparing this fingerprint to an earlier
    stage''s retarded one-pole form in the appendix narrative; it is not an in-script deliverable of stage 104, so its absence
    from the script is not `script_missing_paper_claim`.'
- path: redteam/pass2/reports/stage_104.md
  line: 59
  role: pass2_stage_report
  stage: '104'
  excerpt: '| A2 | sympy | 38 | `Y[z^2] − 1/9 == 0` | Y fingerprint | yes |'
- path: redteam/pass2/reports/stage_104.md
  line: 60
  role: pass2_stage_report
  stage: '104'
  excerpt: '| A3 | sympy | 39 | `Y[z^4] − 4/81 == 0` | Y fingerprint | yes |'
- path: redteam/pass2/reports/stage_104.md
  line: 61
  role: pass2_stage_report
  stage: '104'
  excerpt: '| A4 | sympy | 40 | `Y[z^5]/I − 1/27 == 0` | Y fingerprint (odd) | yes |'
- path: redteam/pass2/reports/stage_104.md
  line: 62
  role: pass2_stage_report
  stage: '104'
  excerpt: '| A5 | sympy | 41 | `Y[z^6] + 11/729 == 0` | Y fingerprint | yes |'
- path: redteam/pass2/reports/stage_104.md
  line: 63
  role: pass2_stage_report
  stage: '104'
  excerpt: '| A6 | sympy | 42 | `Y[z^7]/I + 1/243 == 0` | Y fingerprint (odd) | yes |'
- path: redteam/pass2/reports/stage_104.md
  line: 68
  role: pass2_stage_report
  stage: '104'
  excerpt: '| B2 | mathematica | 43 | `Y[z^2] − 1/9 == 0` | Y fingerprint | yes |'
- path: redteam/pass2/reports/stage_104.md
  line: 69
  role: pass2_stage_report
  stage: '104'
  excerpt: '| B3 | mathematica | 44 | `Y[z^4] − 4/81 == 0` | Y fingerprint | yes |'
- path: redteam/pass2/reports/stage_104.md
  line: 70
  role: pass2_stage_report
  stage: '104'
  excerpt: '| B4 | mathematica | 45 | `Y[z^5]/I − 1/27 == 0` | Y fingerprint (odd) | yes |'
- path: redteam/pass2/reports/stage_104.md
  line: 71
  role: pass2_stage_report
  stage: '104'
  excerpt: '| B5 | mathematica | 46 | `Y[z^6] + 11/729 == 0` | Y fingerprint | yes |'
- path: redteam/pass2/reports/stage_104.md
  line: 72
  role: pass2_stage_report
  stage: '104'
  excerpt: '| B6 | mathematica | 47 | `Y[z^7]/I + 1/243 == 0` | Y fingerprint (odd) | yes |'
- path: redteam/pass2/reports/stage_104.md
  line: 90
  role: pass2_stage_report
  stage: '104'
  excerpt: 'The static-slot check is also done by **different mechanisms**: SymPy uses `sp.limit(Lam, z, 0) + 3` (L37) whereas
    Mathematica reads the `z^0` series coefficient `Coefficient[lamSeries, z, 0] + 3` (L42). The downstream choreography (`Lam
    = z*D[h2,z]/h2`, `Y = -3/Lam`, series, coefficient extraction, compare-to-external-literal) is necessarily parallel because
    that IS the single mathematical definition of the DtN operator stated in the notes — there is no second "route" to a DtN
    operator, and the parallelism here reflects the math, not a transliteration. Crucially, the two series engines (`sp.series`
    operating on an explicit trig closed form vs. `Series[]` operating on a built-in special function) and the differing static-slot
    mechanisms mean the `.wl` is NOT a line-by-line re-typing of the `.py`''s numeric algebra. Verdict: **independent** (the
    primitive and the static-slot mechanism differ; anchors are external literals, not cross-engine echoes). I considered
    `partial` because the mid-pipeline steps are textually similar, but the divergence at the two ends (primitive construction
    + slot mechanism) and the external anchoring place it on the independent side of the line.'
- path: redteam/pass2/reports/stage_104.md
  line: 98
  role: pass2_stage_report
  stage: '104'
  excerpt: '`clean`. I read the card, the notes box, and the appendix "Canonical outgoing DtN branch" subsection before opening
    the scripts, and the script''s verified identity (the exact `\widehat Y_2^out` fingerprint and its ω-form, including `Gamma_5
    = a^5/(27 c_s^5)`) is precisely the card''s stated audit target. Attacks tried and failed: (1) tautology — assertions
    anchor to external rational literals from the notes box, not to in-script re-substitutions, so each can fail; (2) imag-part
    extraction — dividing by `I` is correct and would expose any spurious real part; (3) symbol traps — `a,omega,c_s` declared
    `positive,real` and `z` real are physically justified (throat radius, sound speed, normalized frequency) and do not over-constrain
    the series; (4) transliteration — the Hankel primitive and the static-slot mechanism are constructed by genuinely different
    routes across engines; (5) staleness — both `.txt` outputs are newer than their scripts (sympy out 15:18 vs script 15:08;
    mma out 15:24 vs script 15:08); (6) paper misalignment — the `chi_Q=1` narrative is a cross-stage comparison, not this
    script''s deliverable, so its absence is correct.'
- path: redteam/pass2/reports/stage_105.md
  line: 35
  role: pass2_stage_report
  stage: '105'
  excerpt: 'Stage 105 (checkpoint, anchor MTDC-T8) fixes the last reduced 2.5PN scalar `chi_Q` of the canonical retarded grouped-`P2`
    one-pole module. The card states the result as `\widehat m_0^2 chi_Q N_Q = 1` together with the canonical condition `\chi_Q=1`,
    and quotes: *"Matching the canonical grouped module to the exact outgoing DtN branch gives \(\chi_Q=1\)."* The notes/appendix
    give the mechanism: the retarded module `\widehat Y_Q^{ret}(\omega)=3/4+(1/4)/(1-\omega^2/\Omega_Q^2-i\chi_Q\sigma_Q^{can}\omega^5)`
    expands to `1+a^2\omega^2/(9c_s^2)+4a^4\omega^4/(81c_s^4)+i\chi_Q a^5\omega^5/(27c_s^5)`, with `\Omega_Q=3c_s/(2a)`, `\sigma_Q^{can}=9/(8\Omega_Q^5)=4a^5/(27c_s^5)`.
    Crucially, the **outgoing target** coefficient `i z^5/27` (i.e. `i a^5\omega^5/(27c_s^5)`) is asserted by the paper to
    come from the exact spherical Hankel DtN operator `\Lambda_2^{out}(z)=z\,d/dz\ln h_2^{(1)}(z)=-3+z^2/3+z^4/9+i z^5/9`,
    normalized as `\widehat Y_2^{out}=-3/\Lambda_2^{out}=1+z^2/9+4z^4/81+i z^5/27`. Matching the `O(\omega^5)` coefficients
    fixes `\chi_Q=1`. The card''s third explicit Check requires: *"Check the outgoing \(l=2\) DtN fingerprint against the
    normalized \(z=\omega a/c_s\) expansion."* A second deliverable is the general deformed result `\chi_Q=\xi_Q` (the only
    remaining 2.5PN obstruction is the deviation of the actual branch from the canonical `\xi_Q=1`). Deliverables: (D1) `\sigma_Q^{can}=4a^5/(27c_s^5)`;
    (D2) the retarded series coefficients (z²,z⁴,z⁵); (D3) `\chi_Q=1` from the outgoing DtN fingerprint match; (D4) `\chi_Q=\xi_Q`
    general.'
- path: redteam/pass2/reports/stage_105.md
  line: 47
  role: pass2_stage_report
  stage: '105'
  excerpt: '| D3 `\chi_Q=1` **from outgoing DtN fingerprint** | solve/Reduce retarded ω⁵ coeff `= a^5/(27c_s^5)` → `\chi_Q-1=0`
    (py:49–55; wl:52–59) | **partial** — pins `\chi_Q=1` but matches against a HARDCODED literal `a^5/(27c_s^5)`, NOT the
    Hankel DtN value; the DtN fingerprint identity (`z d/dz ln h_2^{(1)}`) is never exercised (card Check #3 not performed)
    |'
- path: redteam/pass2/reports/stage_105.md
  line: 69
  role: pass2_stage_report
  stage: '105'
  excerpt: A5/B3 are the load-bearing `\chi_Q=1` pin and are the rows that fail the "Anchored?" test.
- path: redteam/pass2/reports/stage_105.md
  line: 83
  role: pass2_stage_report
  stage: '105'
  excerpt: '105 owns the canonical pin `\chi_Q=1`. The paper grounds this pin in a genuine identity: the *outgoing target*
    coefficient `i a^5\omega^5/(27c_s^5)` is the normalized exact-Hankel DtN fingerprint `\widehat Y_2^{out}=-3/\Lambda_2^{out}`,
    `\Lambda_2^{out}(z)=z\,d/dz\ln h_2^{(1)}(z)` (appendix `eq:app-part04-Lambda-out-dtn`/`eq:app-part04-Yout-dtn`), and the
    card''s Check #3 explicitly is "Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion."'
- path: redteam/pass2/reports/stage_105.md
  line: 96
  role: pass2_stage_report
  stage: '105'
  excerpt: 'The left side, `Yret_series.coeff(omega,5)/i` (= `imagPart5`), is **by construction** exactly `chi_Q * a^5/(27c_s^5)`:
    it was generated from `sigma_can`, which line 33/36 *defined* as `(9/8)/Omega_Q^5 = 4a^5/(27c_s^5)`. The right side is
    the same literal `a^5/(27c_s^5)`, typed in directly. So the solved equation is identically `chi_Q · K = K` with `K = a^5/(27c_s^5)`,
    whose root is `chi_Q = 1` for **any** nonzero `K` — the equation cannot detect whether `K` is the true outgoing DtN coefficient.
    The DtN operator `z d/dz ln h_2^{(1)}(z)` is never evaluated, never series-expanded, never normalized; the "fingerprint
    match" the card promises is replaced by self-matching a hardcoded copy of the answer. Re-derived independently: yes, `chi_Q=1`
    is the *correct* value (the retarded ω⁵ coefficient is `i·chi_Q·a^5/(27c_s^5)` and the genuine normalized Hankel value
    is `i·a^5/(27c_s^5)`), but the script''s assertion does not *test* that the two agree — it asserts agreement by typing
    the same number on both sides.'
- path: redteam/pass2/reports/stage_105.md
  line: 98
  role: pass2_stage_report
  stage: '105'
  excerpt: 'Note the asymmetry: the deformed half (A6–A8 / B4) *does* derive its target coefficients from the operator (`-3/\Lambda`
    in SymPy; polynomial inversion of `\Lambda\cdot Y=-3` in Mathematica), so `\chi_Q=\xi_Q` is genuinely exercised. Only
    the canonical `\chi_Q=1` pin — the checkpoint''s headline result — is left self-referential.'
- path: redteam/pass2/reports/stage_105.md
  line: 102
  role: pass2_stage_report
  stage: '105'
  excerpt: 'This is a checkpoint, and 105 is the unit that establishes `\chi_Q=1` as the DtN l=2 fingerprint pin carried forward
    by stages 097/100/106 and the 2.5PN/4PN bridge (107–113). A reader auditing the chain would expect the pin to be forced
    by the Hankel DtN identity, per Check #3. As written, a transcription error or sign error in the *true* DtN coefficient
    could not be caught here: the script would still "pass" because both sides carry the same literal. The pin is correct,
    but it is asserted, not verified.'
- path: redteam/pass2/reports/stage_105.md
  line: 106
  role: pass2_stage_report
  stage: '105'
  excerpt: 'Make the canonical match exercise the actual outgoing DtN fingerprint rather than a hardcoded literal, in BOTH
    engines. Construct the exact `l=2` outgoing DtN operator from the spherical Hankel function and take its series, then
    normalize and read off the imaginary z⁵ coefficient as the *target*, instead of typing `a^5/(27c_s^5)`:'
- path: redteam/pass2/reports/stage_105.md
  line: 110
  role: pass2_stage_report
  stage: '105'
  excerpt: '- assert `Lam_out_series == -3 + z^2/3 + z^4/9 + i z^5/9` (this is the genuine fingerprint exercise — Check #3);'
- path: redteam/pass2/reports/stage_105.md
  line: 114
  role: pass2_stage_report
  stage: '105'
  excerpt: Then the `chi_Q=1` pin is forced by the genuine Hankel DtN coefficient. (If 104 `outgoing_dtn_fingerprint` already
    verifies the Hankel `→ -3+z^2/3+z^4/9+i z^5/9` expansion and 105 is intended to consume that as a carry-forward, the alternative
    acceptable fix is to make the match symbolic against a single named DtN-coefficient symbol `T_dtn`, derive `chi_Q = T_dtn
    / (a^5/(27c_s^5)) · (a^5/(27c_s^5))`-style so the literal appears once, and add an in-script `expect_zero` tying `T_dtn`
    to the 104-verified value — but the self-test prefers deriving the Hankel coefficient locally so the pin is not self-referential.)
- path: redteam/pass2/reports/stage_105.md
  line: 145
  role: pass2_stage_report
  stage: '105'
  excerpt: 'For a checkpoint owning the canonical pin, the retarded-half port means both engines share the same self-referential
    weakness in F1 with no cross-check value — Mathematica cannot catch a transcription error in the SymPy retarded route
    because it re-types the same literal and the same module. Fixing F1 (deriving the DtN coefficient from `h_2^{(1)}` in
    each engine) is the natural opportunity to also de-transliterate: have one engine reach the outgoing coefficient by `Series[z
    D[Log[SphericalHankelH1[2,z]],z]]` and the other by the recurrence/explicit-Hankel-rational form, so the retarded-half
    match no longer mirrors.'
- path: redteam/pass2/reports/stage_105.md
  line: 167
  role: pass2_stage_report
  stage: '105'
  excerpt: '**Checkpoint higher-bar verdict: not-cleared.** chi_Q=1 derivation result: the value is *correct* (independently
    re-derived), but the in-script assertion does not *exercise* the DtN-fingerprint identity — it self-matches a hardcoded
    copy of the target coefficient, so the canonical pin is asserted rather than forced. The higher bar (substantive, non-tautological,
    fingerprint actually exercised) is not met until F1 is fixed.'
- path: redteam/pass2/reports/stage_105.md
  line: 184
  role: pass2_stage_report
  stage: '105'
  excerpt: '| `\chi_Q = 1` (canonical pin) | py:53–55 / wl:58–59; sympy txt:19, math txt:17–18 | notes:69 (boxed); card line
    16, appendix:326 `eq:app-part04-chiQ-equals-one` | MATCH (value); see F1 re anchoring |'
- path: notes/CHECKPOINT_CONSTANT_PROVENANCE.md
  line: 86
  role: checkpoint_constant_provenance
  stage: '104'
  excerpt: fingerprint coefficients from stage 104).
parameter_names:
- Gamma_5
- T_dtn
- chi_Q
- matched_fingerprint_value
- outgoing_dtn_fingerprint
- sigma_Q
- xi_Q
batch_id: null
status: provenance_built
status_timestamps:
  pending: '2026-06-11T07:20:03Z'
  scanned: '2026-06-11T07:22:51Z'
  provenance_built: '2026-06-11T07:23:21Z'
  audit_pending: '2026-06-11T07:21:01Z'
codex_session:
  by_parameter:
    Gamma_5:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    T_dtn:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    chi_Q:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    matched_fingerprint_value:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    outgoing_dtn_fingerprint:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    sigma_Q:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
    xi_Q:
      session_id: null
      log_paths: []
      last_iter: null
      last_exit: null
      last_run: null
paths:
  report: redteam_adversarial/reports/dryrun_stages_003_104_105_chi_q_outgoing_dtn__dry_run_prompt_render.md
  defense: null
  verdict: null
  provenance:
  - redteam_adversarial/provenance/dryrun_stages_003_104_105_chi_q_outgoing_dtn__gamma_5.yaml
  - redteam_adversarial/provenance/dryrun_stages_003_104_105_chi_q_outgoing_dtn__t_dtn.yaml
  - redteam_adversarial/provenance/dryrun_stages_003_104_105_chi_q_outgoing_dtn__chi_q.yaml
  - redteam_adversarial/provenance/dryrun_stages_003_104_105_chi_q_outgoing_dtn__matched_fingerprint_value.yaml
  - redteam_adversarial/provenance/dryrun_stages_003_104_105_chi_q_outgoing_dtn__outgoing_dtn_fingerprint.yaml
  - redteam_adversarial/provenance/dryrun_stages_003_104_105_chi_q_outgoing_dtn__sigma_q.yaml
  - redteam_adversarial/provenance/dryrun_stages_003_104_105_chi_q_outgoing_dtn__xi_q.yaml
  phase_c_prompt: redteam_adversarial/tmp_prompts/dryrun_stages_003_104_105_chi_q_outgoing_dtn__phase_c_adversarial.md
verdict:
  adversarial: null
  adjudication: null
modality_attribution:
- claim_label
- existing_provenance
- numeric_literal
modality_fragments:
- candidate_key: stage_104_chi_q
  modality: numeric_literal
  anchor_stage: '104'
  parameter_names:
  - chi_Q
  reason: target-related numerical literal or closed-form coefficient
  citation:
    path: paper/stages/stage_104.tex
    line: 13
    role: paper_stage_tex
    stage: '104'
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
- candidate_key: stage_104_matched_fingerprint_value
  modality: numeric_literal
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: target-related numerical literal or closed-form coefficient
  citation:
    path: paper/stages/stage_104.tex
    line: 24
    role: paper_stage_tex
    stage: '104'
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
- candidate_key: stage_104_matched_fingerprint_value
  modality: numeric_literal
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: target-related numerical literal or closed-form coefficient
  citation:
    path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 1
    role: notes_stage
    stage: '104'
    excerpt: '# Moving-Throat PDE — Stage 104: Exact Outgoing `l=2` DtN Fingerprint'
- candidate_key: stage_104_matched_fingerprint_value
  modality: numeric_literal
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: target-related numerical literal or closed-form coefficient
  citation:
    path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 83
    role: notes_stage
    stage: '104'
    excerpt: This is exactly the outgoing compact `l=2` fingerprint used earlier in the reduced 2.5PN program, but now derived
      directly from the explicit DtN model.
- candidate_key: stage_104_matched_fingerprint_value
  modality: numeric_literal
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: target-related numerical literal or closed-form coefficient
  citation:
    path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 87
    role: notes_stage
    stage: '104'
    excerpt: The canonical outgoing coefficient is therefore not free. On the exact spherical outgoing `l=2` DtN branch, the
      leading odd quadrupole coefficient is fixed to
- candidate_key: stage_104_matched_fingerprint_value
  modality: numeric_literal
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: target-related numerical literal or closed-form coefficient
  citation:
    path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 92
    role: notes_stage
    stage: '104'
    excerpt: So any later normalization mismatch must come from branch selection or source normalization, not from ambiguity
      in the canonical outgoing `l=2` DtN model itself.
- candidate_key: stage_105_chi_q
  modality: numeric_literal
  anchor_stage: '105'
  parameter_names:
  - chi_Q
  reason: target-related numerical literal or closed-form coefficient
  citation:
    path: paper/stages/stage_105.tex
    line: 13
    role: paper_stage_tex
    stage: '105'
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
- candidate_key: stage_105_chi_q
  modality: numeric_literal
  anchor_stage: '105'
  parameter_names:
  - chi_Q
  - matched_fingerprint_value
  reason: target-related numerical literal or closed-form coefficient
  citation:
    path: paper/stages/stage_105.tex
    line: 16
    role: paper_stage_tex
    stage: '105'
    excerpt: Matching the canonical grouped module to the exact outgoing DtN branch gives \(\chi_Q=1\).
- candidate_key: stage_105_matched_fingerprint_value
  modality: numeric_literal
  anchor_stage: '105'
  parameter_names:
  - matched_fingerprint_value
  reason: target-related numerical literal or closed-form coefficient
  citation:
    path: paper/stages/stage_105.tex
    line: 24
    role: paper_stage_tex
    stage: '105'
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
- candidate_key: stage_105_matched_fingerprint_value
  modality: numeric_literal
  anchor_stage: '105'
  parameter_names:
  - matched_fingerprint_value
  reason: target-related numerical literal or closed-form coefficient
  citation:
    path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 5
    role: notes_stage
    stage: '105'
    excerpt: Use the explicit DtN fingerprint to determine the last reduced 2.5PN scalar
- candidate_key: stage_105_matched_fingerprint_value
  modality: numeric_literal
  anchor_stage: '105'
  parameter_names:
  - matched_fingerprint_value
  - xi_Q
  reason: target-related numerical literal or closed-form coefficient
  citation:
    path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 98
    role: notes_stage
    stage: '105'
    excerpt: That means the only reduced 2.5PN obstruction left after the present calculation is a deviation of the actual
      moving-throat DtN branch from the canonical outgoing `l=2` coefficient \(\xi_Q=1\).
- candidate_key: stage_104_chi_q
  modality: claim_label
  anchor_stage: '104'
  parameter_names:
  - chi_Q
  reason: claim label or status wording near target-related parameter
  citation:
    path: paper/stages/stage_104.tex
    line: 13
    role: paper_stage_tex
    stage: '104'
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
- candidate_key: stage_104_matched_fingerprint_value
  modality: claim_label
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: claim label or status wording near target-related parameter
  citation:
    path: paper/stages/stage_104.tex
    line: 24
    role: paper_stage_tex
    stage: '104'
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
- candidate_key: stage_104_matched_fingerprint_value
  modality: claim_label
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: claim label or status wording near target-related parameter
  citation:
    path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 1
    role: notes_stage
    stage: '104'
    excerpt: '# Moving-Throat PDE — Stage 104: Exact Outgoing `l=2` DtN Fingerprint'
- candidate_key: stage_104_matched_fingerprint_value
  modality: claim_label
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: claim label or status wording near target-related parameter
  citation:
    path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 83
    role: notes_stage
    stage: '104'
    excerpt: This is exactly the outgoing compact `l=2` fingerprint used earlier in the reduced 2.5PN program, but now derived
      directly from the explicit DtN model.
- candidate_key: stage_104_matched_fingerprint_value
  modality: claim_label
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: claim label or status wording near target-related parameter
  citation:
    path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 87
    role: notes_stage
    stage: '104'
    excerpt: The canonical outgoing coefficient is therefore not free. On the exact spherical outgoing `l=2` DtN branch, the
      leading odd quadrupole coefficient is fixed to
- candidate_key: stage_104_matched_fingerprint_value
  modality: claim_label
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: claim label or status wording near target-related parameter
  citation:
    path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 92
    role: notes_stage
    stage: '104'
    excerpt: So any later normalization mismatch must come from branch selection or source normalization, not from ambiguity
      in the canonical outgoing `l=2` DtN model itself.
- candidate_key: stage_105_chi_q
  modality: claim_label
  anchor_stage: '105'
  parameter_names:
  - chi_Q
  reason: claim label or status wording near target-related parameter
  citation:
    path: paper/stages/stage_105.tex
    line: 13
    role: paper_stage_tex
    stage: '105'
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
- candidate_key: stage_105_chi_q
  modality: claim_label
  anchor_stage: '105'
  parameter_names:
  - chi_Q
  - matched_fingerprint_value
  reason: claim label or status wording near target-related parameter
  citation:
    path: paper/stages/stage_105.tex
    line: 16
    role: paper_stage_tex
    stage: '105'
    excerpt: Matching the canonical grouped module to the exact outgoing DtN branch gives \(\chi_Q=1\).
- candidate_key: stage_105_matched_fingerprint_value
  modality: claim_label
  anchor_stage: '105'
  parameter_names:
  - matched_fingerprint_value
  reason: claim label or status wording near target-related parameter
  citation:
    path: paper/stages/stage_105.tex
    line: 24
    role: paper_stage_tex
    stage: '105'
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
- candidate_key: stage_105_matched_fingerprint_value
  modality: claim_label
  anchor_stage: '105'
  parameter_names:
  - matched_fingerprint_value
  reason: claim label or status wording near target-related parameter
  citation:
    path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 5
    role: notes_stage
    stage: '105'
    excerpt: Use the explicit DtN fingerprint to determine the last reduced 2.5PN scalar
- candidate_key: stage_105_matched_fingerprint_value
  modality: claim_label
  anchor_stage: '105'
  parameter_names:
  - matched_fingerprint_value
  - xi_Q
  reason: claim label or status wording near target-related parameter
  citation:
    path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 98
    role: notes_stage
    stage: '105'
    excerpt: That means the only reduced 2.5PN obstruction left after the present calculation is a deviation of the actual
      moving-throat DtN branch from the canonical outgoing `l=2` coefficient \(\xi_Q=1\).
- candidate_key: stage_104_chi_q
  modality: existing_provenance
  anchor_stage: '104'
  parameter_names:
  - chi_Q
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_104.md
    line: 35
    role: pass2_stage_report
    stage: '104'
    excerpt: 'The stage derives the exact outgoing spherical `l=2` Dirichlet-to-Neumann (DtN) operator directly from the spherical
      Hankel mode `h_2^(1)(z) = j_2(z) + i y_2(z)`, with `z := a*omega/c_s`. The operator is `Lambda_2^out(z) = z d/dz ln
      h_2^(1)(z)`, whose static slot is `Lambda_2^out(0) = -3`, and whose small-z series (notes box, lines 28–38) is `-3 +
      z^2/3 + z^4/9 + i z^5/9 - 2 z^6/27 - i z^7/27 + O(z^8)`. Normalizing by the static slot gives the canonical outgoing
      quadrupole admittance (notes box, lines 54–65; card line 16; appendix eq. `Yout-dtn`): `\widehat Y_2^out(z) = 1 + z^2/9
      + 4 z^4/81 + i z^5/27 - 11 z^6/729 - i z^7/243 + O(z^8)`. In ω-form (notes lines 72–81): `1 + a^2 omega^2/(9 c_s^2)
      + 4 a^4 omega^4/(81 c_s^4) + i a^5 omega^5/(27 c_s^5) + O(omega^6)`. The headline consequence (notes "Consequence",
      line 89): the canonical odd quadrupole coefficient is fixed, `Gamma_{5,can}^DtN = a^5/(27 c_s^5)`. The card''s `Purpose`
      states "Its audit target is the verification output quoted below," and the quoted block is exactly the `\widehat Y_2^out`
      fingerprint. The `chi_Q=1` / `m̂_0^2 chi_Q N_Q = 1` framing in the card''s "Derivation ledger" and in appendix eq. `chiQ-equals-one`
      is the downstream comparison that USES this fingerprint (against the retarded one-pole form of an earlier stage); it
      is narrative, not this script''s own deliverable.'
- candidate_key: stage_104_gamma_5
  modality: existing_provenance
  anchor_stage: '104'
  parameter_names:
  - Gamma_5
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_104.md
    line: 39
    role: pass2_stage_report
    stage: '104'
    excerpt: 'Both scripts build the spherical Hankel mode, form `Lam = z * d/dz(h2)/h2` and `Y = -3/Lam`, series-expand both
      to O(z^8)/O(z^7), and then assert, against external rational literals taken from the notes box: the static DtN slot
      (`Lambda(0) = -3`), the `Y` coefficients at `z^2 (1/9)`, `z^4 (4/81)`, `z^5 (i/27)`, `z^6 (-11/729)`, `z^7 (-i/243)`,
      and after the substitution `z -> a*omega/c_s` the ω-form coefficients at `omega^2 (a^2/9c_s^2)`, `omega^4 (4 a^4/81
      c_s^4)`, `omega^5 (i a^5/27 c_s^5)`. The final printed RESULT states the model reproduces the canonical fingerprint
      including the odd coefficient `a^5/(27 c_s^5)` (= `Gamma_5`).'
- candidate_key: stage_104_chi_q
  modality: existing_provenance
  anchor_stage: '104'
  parameter_names:
  - chi_Q
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_104.md
    line: 52
    role: pass2_stage_report
    stage: '104'
    excerpt: '`paper_alignment: aligned`. The card''s stated audit target (the `\widehat Y_2^out` fingerprint expansion) is
      fully and non-tautologically exercised by both engines. `chi_Q=1` is fixed by comparing this fingerprint to an earlier
      stage''s retarded one-pole form in the appendix narrative; it is not an in-script deliverable of stage 104, so its absence
      from the script is not `script_missing_paper_claim`.'
- candidate_key: stage_104_matched_fingerprint_value
  modality: existing_provenance
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_104.md
    line: 59
    role: pass2_stage_report
    stage: '104'
    excerpt: '| A2 | sympy | 38 | `Y[z^2] − 1/9 == 0` | Y fingerprint | yes |'
- candidate_key: stage_104_matched_fingerprint_value
  modality: existing_provenance
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_104.md
    line: 60
    role: pass2_stage_report
    stage: '104'
    excerpt: '| A3 | sympy | 39 | `Y[z^4] − 4/81 == 0` | Y fingerprint | yes |'
- candidate_key: stage_104_matched_fingerprint_value
  modality: existing_provenance
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_104.md
    line: 61
    role: pass2_stage_report
    stage: '104'
    excerpt: '| A4 | sympy | 40 | `Y[z^5]/I − 1/27 == 0` | Y fingerprint (odd) | yes |'
- candidate_key: stage_104_matched_fingerprint_value
  modality: existing_provenance
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_104.md
    line: 62
    role: pass2_stage_report
    stage: '104'
    excerpt: '| A5 | sympy | 41 | `Y[z^6] + 11/729 == 0` | Y fingerprint | yes |'
- candidate_key: stage_104_matched_fingerprint_value
  modality: existing_provenance
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_104.md
    line: 63
    role: pass2_stage_report
    stage: '104'
    excerpt: '| A6 | sympy | 42 | `Y[z^7]/I + 1/243 == 0` | Y fingerprint (odd) | yes |'
- candidate_key: stage_104_matched_fingerprint_value
  modality: existing_provenance
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_104.md
    line: 68
    role: pass2_stage_report
    stage: '104'
    excerpt: '| B2 | mathematica | 43 | `Y[z^2] − 1/9 == 0` | Y fingerprint | yes |'
- candidate_key: stage_104_matched_fingerprint_value
  modality: existing_provenance
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_104.md
    line: 69
    role: pass2_stage_report
    stage: '104'
    excerpt: '| B3 | mathematica | 44 | `Y[z^4] − 4/81 == 0` | Y fingerprint | yes |'
- candidate_key: stage_104_matched_fingerprint_value
  modality: existing_provenance
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_104.md
    line: 70
    role: pass2_stage_report
    stage: '104'
    excerpt: '| B4 | mathematica | 45 | `Y[z^5]/I − 1/27 == 0` | Y fingerprint (odd) | yes |'
- candidate_key: stage_104_matched_fingerprint_value
  modality: existing_provenance
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_104.md
    line: 71
    role: pass2_stage_report
    stage: '104'
    excerpt: '| B5 | mathematica | 46 | `Y[z^6] + 11/729 == 0` | Y fingerprint | yes |'
- candidate_key: stage_104_matched_fingerprint_value
  modality: existing_provenance
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_104.md
    line: 72
    role: pass2_stage_report
    stage: '104'
    excerpt: '| B6 | mathematica | 47 | `Y[z^7]/I + 1/243 == 0` | Y fingerprint (odd) | yes |'
- candidate_key: stage_104_matched_fingerprint_value
  modality: existing_provenance
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_104.md
    line: 90
    role: pass2_stage_report
    stage: '104'
    excerpt: 'The static-slot check is also done by **different mechanisms**: SymPy uses `sp.limit(Lam, z, 0) + 3` (L37) whereas
      Mathematica reads the `z^0` series coefficient `Coefficient[lamSeries, z, 0] + 3` (L42). The downstream choreography
      (`Lam = z*D[h2,z]/h2`, `Y = -3/Lam`, series, coefficient extraction, compare-to-external-literal) is necessarily parallel
      because that IS the single mathematical definition of the DtN operator stated in the notes — there is no second "route"
      to a DtN operator, and the parallelism here reflects the math, not a transliteration. Crucially, the two series engines
      (`sp.series` operating on an explicit trig closed form vs. `Series[]` operating on a built-in special function) and
      the differing static-slot mechanisms mean the `.wl` is NOT a line-by-line re-typing of the `.py`''s numeric algebra.
      Verdict: **independent** (the primitive and the static-slot mechanism differ; anchors are external literals, not cross-engine
      echoes). I considered `partial` because the mid-pipeline steps are textually similar, but the divergence at the two
      ends (primitive construction + slot mechanism) and the external anchoring place it on the independent side of the line.'
- candidate_key: stage_104_gamma_5
  modality: existing_provenance
  anchor_stage: '104'
  parameter_names:
  - Gamma_5
  - chi_Q
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_104.md
    line: 98
    role: pass2_stage_report
    stage: '104'
    excerpt: '`clean`. I read the card, the notes box, and the appendix "Canonical outgoing DtN branch" subsection before
      opening the scripts, and the script''s verified identity (the exact `\widehat Y_2^out` fingerprint and its ω-form, including
      `Gamma_5 = a^5/(27 c_s^5)`) is precisely the card''s stated audit target. Attacks tried and failed: (1) tautology —
      assertions anchor to external rational literals from the notes box, not to in-script re-substitutions, so each can fail;
      (2) imag-part extraction — dividing by `I` is correct and would expose any spurious real part; (3) symbol traps — `a,omega,c_s`
      declared `positive,real` and `z` real are physically justified (throat radius, sound speed, normalized frequency) and
      do not over-constrain the series; (4) transliteration — the Hankel primitive and the static-slot mechanism are constructed
      by genuinely different routes across engines; (5) staleness — both `.txt` outputs are newer than their scripts (sympy
      out 15:18 vs script 15:08; mma out 15:24 vs script 15:08); (6) paper misalignment — the `chi_Q=1` narrative is a cross-stage
      comparison, not this script''s deliverable, so its absence is correct.'
- candidate_key: stage_105_chi_q
  modality: existing_provenance
  anchor_stage: '105'
  parameter_names:
  - chi_Q
  - matched_fingerprint_value
  - sigma_Q
  - xi_Q
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_105.md
    line: 35
    role: pass2_stage_report
    stage: '105'
    excerpt: 'Stage 105 (checkpoint, anchor MTDC-T8) fixes the last reduced 2.5PN scalar `chi_Q` of the canonical retarded
      grouped-`P2` one-pole module. The card states the result as `\widehat m_0^2 chi_Q N_Q = 1` together with the canonical
      condition `\chi_Q=1`, and quotes: *"Matching the canonical grouped module to the exact outgoing DtN branch gives \(\chi_Q=1\)."*
      The notes/appendix give the mechanism: the retarded module `\widehat Y_Q^{ret}(\omega)=3/4+(1/4)/(1-\omega^2/\Omega_Q^2-i\chi_Q\sigma_Q^{can}\omega^5)`
      expands to `1+a^2\omega^2/(9c_s^2)+4a^4\omega^4/(81c_s^4)+i\chi_Q a^5\omega^5/(27c_s^5)`, with `\Omega_Q=3c_s/(2a)`,
      `\sigma_Q^{can}=9/(8\Omega_Q^5)=4a^5/(27c_s^5)`. Crucially, the **outgoing target** coefficient `i z^5/27` (i.e. `i
      a^5\omega^5/(27c_s^5)`) is asserted by the paper to come from the exact spherical Hankel DtN operator `\Lambda_2^{out}(z)=z\,d/dz\ln
      h_2^{(1)}(z)=-3+z^2/3+z^4/9+i z^5/9`, normalized as `\widehat Y_2^{out}=-3/\Lambda_2^{out}=1+z^2/9+4z^4/81+i z^5/27`.
      Matching the `O(\omega^5)` coefficients fixes `\chi_Q=1`. The card''s third explicit Check requires: *"Check the outgoing
      \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion."* A second deliverable is the general deformed
      result `\chi_Q=\xi_Q` (the only remaining 2.5PN obstruction is the deviation of the actual branch from the canonical
      `\xi_Q=1`). Deliverables: (D1) `\sigma_Q^{can}=4a^5/(27c_s^5)`; (D2) the retarded series coefficients (z²,z⁴,z⁵); (D3)
      `\chi_Q=1` from the outgoing DtN fingerprint match; (D4) `\chi_Q=\xi_Q` general.'
- candidate_key: stage_105_chi_q
  modality: existing_provenance
  anchor_stage: '105'
  parameter_names:
  - chi_Q
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_105.md
    line: 47
    role: pass2_stage_report
    stage: '105'
    excerpt: '| D3 `\chi_Q=1` **from outgoing DtN fingerprint** | solve/Reduce retarded ω⁵ coeff `= a^5/(27c_s^5)` → `\chi_Q-1=0`
      (py:49–55; wl:52–59) | **partial** — pins `\chi_Q=1` but matches against a HARDCODED literal `a^5/(27c_s^5)`, NOT the
      Hankel DtN value; the DtN fingerprint identity (`z d/dz ln h_2^{(1)}`) is never exercised (card Check #3 not performed)
      |'
- candidate_key: stage_105_chi_q
  modality: existing_provenance
  anchor_stage: '105'
  parameter_names:
  - chi_Q
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_105.md
    line: 69
    role: pass2_stage_report
    stage: '105'
    excerpt: A5/B3 are the load-bearing `\chi_Q=1` pin and are the rows that fail the "Anchored?" test.
- candidate_key: stage_105_chi_q
  modality: existing_provenance
  anchor_stage: '105'
  parameter_names:
  - chi_Q
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_105.md
    line: 83
    role: pass2_stage_report
    stage: '105'
    excerpt: '105 owns the canonical pin `\chi_Q=1`. The paper grounds this pin in a genuine identity: the *outgoing target*
      coefficient `i a^5\omega^5/(27c_s^5)` is the normalized exact-Hankel DtN fingerprint `\widehat Y_2^{out}=-3/\Lambda_2^{out}`,
      `\Lambda_2^{out}(z)=z\,d/dz\ln h_2^{(1)}(z)` (appendix `eq:app-part04-Lambda-out-dtn`/`eq:app-part04-Yout-dtn`), and
      the card''s Check #3 explicitly is "Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\)
      expansion."'
- candidate_key: stage_105_chi_q
  modality: existing_provenance
  anchor_stage: '105'
  parameter_names:
  - chi_Q
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_105.md
    line: 96
    role: pass2_stage_report
    stage: '105'
    excerpt: 'The left side, `Yret_series.coeff(omega,5)/i` (= `imagPart5`), is **by construction** exactly `chi_Q * a^5/(27c_s^5)`:
      it was generated from `sigma_can`, which line 33/36 *defined* as `(9/8)/Omega_Q^5 = 4a^5/(27c_s^5)`. The right side
      is the same literal `a^5/(27c_s^5)`, typed in directly. So the solved equation is identically `chi_Q · K = K` with `K
      = a^5/(27c_s^5)`, whose root is `chi_Q = 1` for **any** nonzero `K` — the equation cannot detect whether `K` is the
      true outgoing DtN coefficient. The DtN operator `z d/dz ln h_2^{(1)}(z)` is never evaluated, never series-expanded,
      never normalized; the "fingerprint match" the card promises is replaced by self-matching a hardcoded copy of the answer.
      Re-derived independently: yes, `chi_Q=1` is the *correct* value (the retarded ω⁵ coefficient is `i·chi_Q·a^5/(27c_s^5)`
      and the genuine normalized Hankel value is `i·a^5/(27c_s^5)`), but the script''s assertion does not *test* that the
      two agree — it asserts agreement by typing the same number on both sides.'
- candidate_key: stage_105_chi_q
  modality: existing_provenance
  anchor_stage: '105'
  parameter_names:
  - chi_Q
  - xi_Q
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_105.md
    line: 98
    role: pass2_stage_report
    stage: '105'
    excerpt: 'Note the asymmetry: the deformed half (A6–A8 / B4) *does* derive its target coefficients from the operator (`-3/\Lambda`
      in SymPy; polynomial inversion of `\Lambda\cdot Y=-3` in Mathematica), so `\chi_Q=\xi_Q` is genuinely exercised. Only
      the canonical `\chi_Q=1` pin — the checkpoint''s headline result — is left self-referential.'
- candidate_key: stage_105_chi_q
  modality: existing_provenance
  anchor_stage: '105'
  parameter_names:
  - chi_Q
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_105.md
    line: 102
    role: pass2_stage_report
    stage: '105'
    excerpt: 'This is a checkpoint, and 105 is the unit that establishes `\chi_Q=1` as the DtN l=2 fingerprint pin carried
      forward by stages 097/100/106 and the 2.5PN/4PN bridge (107–113). A reader auditing the chain would expect the pin to
      be forced by the Hankel DtN identity, per Check #3. As written, a transcription error or sign error in the *true* DtN
      coefficient could not be caught here: the script would still "pass" because both sides carry the same literal. The pin
      is correct, but it is asserted, not verified.'
- candidate_key: stage_105_matched_fingerprint_value
  modality: existing_provenance
  anchor_stage: '105'
  parameter_names:
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_105.md
    line: 106
    role: pass2_stage_report
    stage: '105'
    excerpt: 'Make the canonical match exercise the actual outgoing DtN fingerprint rather than a hardcoded literal, in BOTH
      engines. Construct the exact `l=2` outgoing DtN operator from the spherical Hankel function and take its series, then
      normalize and read off the imaginary z⁵ coefficient as the *target*, instead of typing `a^5/(27c_s^5)`:'
- candidate_key: stage_105_matched_fingerprint_value
  modality: existing_provenance
  anchor_stage: '105'
  parameter_names:
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_105.md
    line: 110
    role: pass2_stage_report
    stage: '105'
    excerpt: '- assert `Lam_out_series == -3 + z^2/3 + z^4/9 + i z^5/9` (this is the genuine fingerprint exercise — Check
      #3);'
- candidate_key: stage_105_t_dtn
  modality: existing_provenance
  anchor_stage: '105'
  parameter_names:
  - T_dtn
  - chi_Q
  - matched_fingerprint_value
  - outgoing_dtn_fingerprint
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_105.md
    line: 114
    role: pass2_stage_report
    stage: '105'
    excerpt: Then the `chi_Q=1` pin is forced by the genuine Hankel DtN coefficient. (If 104 `outgoing_dtn_fingerprint` already
      verifies the Hankel `→ -3+z^2/3+z^4/9+i z^5/9` expansion and 105 is intended to consume that as a carry-forward, the
      alternative acceptable fix is to make the match symbolic against a single named DtN-coefficient symbol `T_dtn`, derive
      `chi_Q = T_dtn / (a^5/(27c_s^5)) · (a^5/(27c_s^5))`-style so the literal appears once, and add an in-script `expect_zero`
      tying `T_dtn` to the 104-verified value — but the self-test prefers deriving the Hankel coefficient locally so the pin
      is not self-referential.)
- candidate_key: stage_105_matched_fingerprint_value
  modality: existing_provenance
  anchor_stage: '105'
  parameter_names:
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_105.md
    line: 145
    role: pass2_stage_report
    stage: '105'
    excerpt: 'For a checkpoint owning the canonical pin, the retarded-half port means both engines share the same self-referential
      weakness in F1 with no cross-check value — Mathematica cannot catch a transcription error in the SymPy retarded route
      because it re-types the same literal and the same module. Fixing F1 (deriving the DtN coefficient from `h_2^{(1)}` in
      each engine) is the natural opportunity to also de-transliterate: have one engine reach the outgoing coefficient by
      `Series[z D[Log[SphericalHankelH1[2,z]],z]]` and the other by the recurrence/explicit-Hankel-rational form, so the retarded-half
      match no longer mirrors.'
- candidate_key: stage_105_chi_q
  modality: existing_provenance
  anchor_stage: '105'
  parameter_names:
  - chi_Q
  - matched_fingerprint_value
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_105.md
    line: 167
    role: pass2_stage_report
    stage: '105'
    excerpt: '**Checkpoint higher-bar verdict: not-cleared.** chi_Q=1 derivation result: the value is *correct* (independently
      re-derived), but the in-script assertion does not *exercise* the DtN-fingerprint identity — it self-matches a hardcoded
      copy of the target coefficient, so the canonical pin is asserted rather than forced. The higher bar (substantive, non-tautological,
      fingerprint actually exercised) is not met until F1 is fixed.'
- candidate_key: stage_105_chi_q
  modality: existing_provenance
  anchor_stage: '105'
  parameter_names:
  - chi_Q
  reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
  citation:
    path: redteam/pass2/reports/stage_105.md
    line: 184
    role: pass2_stage_report
    stage: '105'
    excerpt: '| `\chi_Q = 1` (canonical pin) | py:53–55 / wl:58–59; sympy txt:19, math txt:17–18 | notes:69 (boxed); card
      line 16, appendix:326 `eq:app-part04-chiQ-equals-one` | MATCH (value); see F1 re anchoring |'
- candidate_key: stage_104_matched_fingerprint_value
  modality: existing_provenance
  anchor_stage: '104'
  parameter_names:
  - matched_fingerprint_value
  reason: checkpoint provenance seed mentions a candidate parameter
  citation:
    path: notes/CHECKPOINT_CONSTANT_PROVENANCE.md
    line: 86
    role: checkpoint_constant_provenance
    stage: '104'
    excerpt: fingerprint coefficients from stage 104).
```

## Primary Sources

- notes/CHECKPOINT_CONSTANT_PROVENANCE.md
- notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
- notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
- paper/stages/stage_104.tex
- paper/stages/stage_105.tex
- redteam/pass2/reports/stage_104.md
- redteam/pass2/reports/stage_105.md

## Provenance Slice

```yaml
- schema_version: 1
  candidate_id: dryrun_stages_003_104_105_chi_q_outgoing_dtn
  parameter_name: Gamma_5
  dry_run: true
  dry_run_id: stages_003_104_105
  non_binding: true
  generated_content_kind: mechanical_evidence_bundle
  agent_prompt_path: redteam_adversarial/tmp_prompts/dryrun_stages_003_104_105_chi_q_outgoing_dtn__phase_b__gamma_5.md
  agent_synthesis_required: true
  dry_run_fixture: true
  fixture_note: Dry-run provenance stores only opened source evidence and graph context. Origin claims, constraints, and contradiction
    findings must be produced by the Phase B agent path.
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '104'
  - '105'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_104.tex
    line: 1
    role: paper_stage_tex
    excerpt: '\section[Stage~104]{Stage~104: Exact Outgoing l=2 DtN Fingerprint}'
  - path: paper/stages/stage_104.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{This stage imports the cleaned conservative isotropic branch, the source-map factor \(\widehat
      m_0\), the conservative normalization \(N_Q\), and the outgoing \(l=2\) DtN coefficient \(\chi_Q\).  It does not change
      the parent action or the grouped-lane ontology; it works inside the declared reduced branch and tests the next algebraic
      or realization condition.}
  - path: paper/stages/stage_104.tex
    line: 11
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py};
      transcript: \StageFile{scripts/output/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.txt}.  Mathematica
      audit: \StageFile{mathematica/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.wl}; transcript:
      \StageFile{mathematica/output/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.txt}.}'
  - path: paper/stages/stage_104.tex
    line: 13
    role: paper_stage_tex
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
  - path: paper/stages/stage_104.tex
    line: 16
    role: paper_stage_tex
    excerpt: Spherical Hankel DtN expansion gives \(\widehat Y_2^{\rm out}=1+z^2/9+4z^4/81+iz^5/27+\cdots\).
  - path: paper/stages/stage_104.tex
    line: 22
    role: paper_stage_tex
    excerpt: \item Check the product \(\widehat m_0^{\,2}\chi_QN_Q\) keeps source, conservative, and outgoing factors separate.
  - path: paper/stages/stage_104.tex
    line: 24
    role: paper_stage_tex
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
  - path: paper/stages/stage_104.tex
    line: 27
    role: paper_stage_tex
    excerpt: \stagefield{Downstream use}{This stage feeds Stages~107--113, the \(2.5\)PN/4PN outgoing bridge, and later same-charge
      prefactor audits.  When cited downstream, the status tag above must be carried with the result; the card is a derivation
      ledger entry, not an unconditional actual-branch theorem.}
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 104: Exact Outgoing `l=2` DtN Fingerprint'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 5
    role: notes_stage
    excerpt: Compute the canonical compact passive/outgoing quadrupole normalization directly from an explicit outgoing `l=2`
      Dirichlet-to-Neumann model instead of carrying the coefficient `chi_Q` symbolically.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 7
    role: notes_stage
    excerpt: '## Exact outgoing DtN model'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 13
    role: notes_stage
    excerpt: and let the outgoing partial wave be the spherical Hankel mode
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 18
    role: notes_stage
    excerpt: The exact `l=2` outgoing DtN operator is
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 42
    role: notes_stage
    excerpt: '## Normalized outgoing quadrupole admittance'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 44
    role: notes_stage
    excerpt: Define the normalized outgoing branch by
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 83
    role: notes_stage
    excerpt: This is exactly the outgoing compact `l=2` fingerprint used earlier in the reduced 2.5PN program, but now derived
      directly from the explicit DtN model.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 87
    role: notes_stage
    excerpt: The canonical outgoing coefficient is therefore not free. On the exact spherical outgoing `l=2` DtN branch, the
      leading odd quadrupole coefficient is fixed to
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 89
    role: notes_stage
    excerpt: \Gamma_{5,\rm can}^{\rm DtN}=\frac{a^5}{27c_s^5}.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 92
    role: notes_stage
    excerpt: So any later normalization mismatch must come from branch selection or source normalization, not from ambiguity
      in the canonical outgoing `l=2` DtN model itself.
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 18
    role: sympy_script
    excerpt: banner("STAGE 104 — EXACT OUTGOING l=2 DtN FINGERPRINT")
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 37
    role: sympy_script
    excerpt: expect_zero("static DtN slot", sp.limit(Lam, z, 0) + 3)
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 55
    role: sympy_script
    excerpt: print("  The explicit outgoing spherical l=2 DtN model reproduces the canonical compact")
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 56
    role: sympy_script
    excerpt: print("  quadrupole fingerprint exactly, including the odd coefficient a^5/(27 c_s^5).")
  - path: paper/stages/stage_105.tex
    line: 1
    role: paper_stage_tex
    excerpt: '\section[Stage~105]{Stage~105: Exact Fixing of chi-Q}'
  - path: paper/stages/stage_105.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{This stage imports the cleaned conservative isotropic branch, the source-map factor \(\widehat
      m_0\), the conservative normalization \(N_Q\), and the outgoing \(l=2\) DtN coefficient \(\chi_Q\).  It does not change
      the parent action or the grouped-lane ontology; it works inside the declared reduced branch and tests the next algebraic
      or realization condition.}
  - path: paper/stages/stage_105.tex
    line: 11
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py};
      transcript: \StageFile{scripts/output/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.txt}.  Mathematica
      audit: \StageFile{mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl}; transcript:
      \StageFile{mathematica/output/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.txt}.}'
  - path: paper/stages/stage_105.tex
    line: 13
    role: paper_stage_tex
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
  - path: paper/stages/stage_105.tex
    line: 16
    role: paper_stage_tex
    excerpt: Matching the canonical grouped module to the exact outgoing DtN branch gives \(\chi_Q=1\).
  - path: paper/stages/stage_105.tex
    line: 22
    role: paper_stage_tex
    excerpt: \item Check the product \(\widehat m_0^{\,2}\chi_QN_Q\) keeps source, conservative, and outgoing factors separate.
  - path: paper/stages/stage_105.tex
    line: 24
    role: paper_stage_tex
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
  - path: paper/stages/stage_105.tex
    line: 27
    role: paper_stage_tex
    excerpt: \stagefield{Downstream use}{This stage feeds Stages~107--113, the \(2.5\)PN/4PN outgoing bridge, and later same-charge
      prefactor audits.  When cited downstream, the status tag above must be carried with the result; the card is a derivation
      ledger entry, not an unconditional actual-branch theorem.}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 105: Exact Fixing of `chi_Q`'
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 5
    role: notes_stage
    excerpt: Use the explicit DtN fingerprint to determine the last reduced 2.5PN scalar
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 7
    role: notes_stage
    excerpt: \chi_Q
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 22
    role: notes_stage
    excerpt: -i\,\chi_Q\,\sigma_Q^{\rm can}\,\omega^5
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 34
    role: notes_stage
    excerpt: \sigma_Q^{\rm can}:=\frac{9}{8\Omega_Q^5}.
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 39
    role: notes_stage
    excerpt: \sigma_Q^{\rm can}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 53
    role: notes_stage
    excerpt: +i\,\chi_Q\,\frac{a^5\omega^5}{27c_s^5}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 57
    role: notes_stage
    excerpt: From Stage 87, the explicit outgoing DtN branch gives
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 67
    role: notes_stage
    excerpt: Matching the \(O(\omega^5)\) coefficient yields
  origin_claims: []
  constraints: []
  graph_context:
    contexts: []
    graph_gaps:
    - source: paper/stages/stage_104.tex
      attempted_sources:
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: paper/stages/stage_105.tex
      attempted_sources:
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
  provenance_findings:
  - type: graph_gap
    severity: dry_run_non_binding
    summary: The atlas graph wrapper returned no source nodes for one or more primary sources.
    gaps:
    - source: paper/stages/stage_104.tex
      attempted_sources:
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: paper/stages/stage_105.tex
      attempted_sources:
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
- schema_version: 1
  candidate_id: dryrun_stages_003_104_105_chi_q_outgoing_dtn
  parameter_name: T_dtn
  dry_run: true
  dry_run_id: stages_003_104_105
  non_binding: true
  generated_content_kind: mechanical_evidence_bundle
  agent_prompt_path: redteam_adversarial/tmp_prompts/dryrun_stages_003_104_105_chi_q_outgoing_dtn__phase_b__t_dtn.md
  agent_synthesis_required: true
  dry_run_fixture: true
  fixture_note: Dry-run provenance stores only opened source evidence and graph context. Origin claims, constraints, and contradiction
    findings must be produced by the Phase B agent path.
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '104'
  - '105'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_104.tex
    line: 1
    role: paper_stage_tex
    excerpt: '\section[Stage~104]{Stage~104: Exact Outgoing l=2 DtN Fingerprint}'
  - path: paper/stages/stage_104.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{This stage imports the cleaned conservative isotropic branch, the source-map factor \(\widehat
      m_0\), the conservative normalization \(N_Q\), and the outgoing \(l=2\) DtN coefficient \(\chi_Q\).  It does not change
      the parent action or the grouped-lane ontology; it works inside the declared reduced branch and tests the next algebraic
      or realization condition.}
  - path: paper/stages/stage_104.tex
    line: 11
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py};
      transcript: \StageFile{scripts/output/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.txt}.  Mathematica
      audit: \StageFile{mathematica/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.wl}; transcript:
      \StageFile{mathematica/output/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.txt}.}'
  - path: paper/stages/stage_104.tex
    line: 13
    role: paper_stage_tex
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
  - path: paper/stages/stage_104.tex
    line: 16
    role: paper_stage_tex
    excerpt: Spherical Hankel DtN expansion gives \(\widehat Y_2^{\rm out}=1+z^2/9+4z^4/81+iz^5/27+\cdots\).
  - path: paper/stages/stage_104.tex
    line: 22
    role: paper_stage_tex
    excerpt: \item Check the product \(\widehat m_0^{\,2}\chi_QN_Q\) keeps source, conservative, and outgoing factors separate.
  - path: paper/stages/stage_104.tex
    line: 24
    role: paper_stage_tex
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
  - path: paper/stages/stage_104.tex
    line: 27
    role: paper_stage_tex
    excerpt: \stagefield{Downstream use}{This stage feeds Stages~107--113, the \(2.5\)PN/4PN outgoing bridge, and later same-charge
      prefactor audits.  When cited downstream, the status tag above must be carried with the result; the card is a derivation
      ledger entry, not an unconditional actual-branch theorem.}
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 104: Exact Outgoing `l=2` DtN Fingerprint'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 5
    role: notes_stage
    excerpt: Compute the canonical compact passive/outgoing quadrupole normalization directly from an explicit outgoing `l=2`
      Dirichlet-to-Neumann model instead of carrying the coefficient `chi_Q` symbolically.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 7
    role: notes_stage
    excerpt: '## Exact outgoing DtN model'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 13
    role: notes_stage
    excerpt: and let the outgoing partial wave be the spherical Hankel mode
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 18
    role: notes_stage
    excerpt: The exact `l=2` outgoing DtN operator is
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 42
    role: notes_stage
    excerpt: '## Normalized outgoing quadrupole admittance'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 44
    role: notes_stage
    excerpt: Define the normalized outgoing branch by
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 83
    role: notes_stage
    excerpt: This is exactly the outgoing compact `l=2` fingerprint used earlier in the reduced 2.5PN program, but now derived
      directly from the explicit DtN model.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 87
    role: notes_stage
    excerpt: The canonical outgoing coefficient is therefore not free. On the exact spherical outgoing `l=2` DtN branch, the
      leading odd quadrupole coefficient is fixed to
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 89
    role: notes_stage
    excerpt: \Gamma_{5,\rm can}^{\rm DtN}=\frac{a^5}{27c_s^5}.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 92
    role: notes_stage
    excerpt: So any later normalization mismatch must come from branch selection or source normalization, not from ambiguity
      in the canonical outgoing `l=2` DtN model itself.
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 18
    role: sympy_script
    excerpt: banner("STAGE 104 — EXACT OUTGOING l=2 DtN FINGERPRINT")
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 37
    role: sympy_script
    excerpt: expect_zero("static DtN slot", sp.limit(Lam, z, 0) + 3)
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 55
    role: sympy_script
    excerpt: print("  The explicit outgoing spherical l=2 DtN model reproduces the canonical compact")
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 56
    role: sympy_script
    excerpt: print("  quadrupole fingerprint exactly, including the odd coefficient a^5/(27 c_s^5).")
  - path: paper/stages/stage_105.tex
    line: 1
    role: paper_stage_tex
    excerpt: '\section[Stage~105]{Stage~105: Exact Fixing of chi-Q}'
  - path: paper/stages/stage_105.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{This stage imports the cleaned conservative isotropic branch, the source-map factor \(\widehat
      m_0\), the conservative normalization \(N_Q\), and the outgoing \(l=2\) DtN coefficient \(\chi_Q\).  It does not change
      the parent action or the grouped-lane ontology; it works inside the declared reduced branch and tests the next algebraic
      or realization condition.}
  - path: paper/stages/stage_105.tex
    line: 11
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py};
      transcript: \StageFile{scripts/output/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.txt}.  Mathematica
      audit: \StageFile{mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl}; transcript:
      \StageFile{mathematica/output/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.txt}.}'
  - path: paper/stages/stage_105.tex
    line: 13
    role: paper_stage_tex
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
  - path: paper/stages/stage_105.tex
    line: 16
    role: paper_stage_tex
    excerpt: Matching the canonical grouped module to the exact outgoing DtN branch gives \(\chi_Q=1\).
  - path: paper/stages/stage_105.tex
    line: 22
    role: paper_stage_tex
    excerpt: \item Check the product \(\widehat m_0^{\,2}\chi_QN_Q\) keeps source, conservative, and outgoing factors separate.
  - path: paper/stages/stage_105.tex
    line: 24
    role: paper_stage_tex
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
  - path: paper/stages/stage_105.tex
    line: 27
    role: paper_stage_tex
    excerpt: \stagefield{Downstream use}{This stage feeds Stages~107--113, the \(2.5\)PN/4PN outgoing bridge, and later same-charge
      prefactor audits.  When cited downstream, the status tag above must be carried with the result; the card is a derivation
      ledger entry, not an unconditional actual-branch theorem.}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 105: Exact Fixing of `chi_Q`'
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 5
    role: notes_stage
    excerpt: Use the explicit DtN fingerprint to determine the last reduced 2.5PN scalar
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 7
    role: notes_stage
    excerpt: \chi_Q
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 22
    role: notes_stage
    excerpt: -i\,\chi_Q\,\sigma_Q^{\rm can}\,\omega^5
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 34
    role: notes_stage
    excerpt: \sigma_Q^{\rm can}:=\frac{9}{8\Omega_Q^5}.
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 39
    role: notes_stage
    excerpt: \sigma_Q^{\rm can}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 53
    role: notes_stage
    excerpt: +i\,\chi_Q\,\frac{a^5\omega^5}{27c_s^5}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 57
    role: notes_stage
    excerpt: From Stage 87, the explicit outgoing DtN branch gives
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 67
    role: notes_stage
    excerpt: Matching the \(O(\omega^5)\) coefficient yields
  origin_claims: []
  constraints: []
  graph_context:
    contexts: []
    graph_gaps:
    - source: paper/stages/stage_104.tex
      attempted_sources:
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: paper/stages/stage_105.tex
      attempted_sources:
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
  provenance_findings:
  - type: graph_gap
    severity: dry_run_non_binding
    summary: The atlas graph wrapper returned no source nodes for one or more primary sources.
    gaps:
    - source: paper/stages/stage_104.tex
      attempted_sources:
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: paper/stages/stage_105.tex
      attempted_sources:
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
- schema_version: 1
  candidate_id: dryrun_stages_003_104_105_chi_q_outgoing_dtn
  parameter_name: chi_Q
  dry_run: true
  dry_run_id: stages_003_104_105
  non_binding: true
  generated_content_kind: mechanical_evidence_bundle
  agent_prompt_path: redteam_adversarial/tmp_prompts/dryrun_stages_003_104_105_chi_q_outgoing_dtn__phase_b__chi_q.md
  agent_synthesis_required: true
  dry_run_fixture: true
  fixture_note: Dry-run provenance stores only opened source evidence and graph context. Origin claims, constraints, and contradiction
    findings must be produced by the Phase B agent path.
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '104'
  - '105'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_104.tex
    line: 1
    role: paper_stage_tex
    excerpt: '\section[Stage~104]{Stage~104: Exact Outgoing l=2 DtN Fingerprint}'
  - path: paper/stages/stage_104.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{This stage imports the cleaned conservative isotropic branch, the source-map factor \(\widehat
      m_0\), the conservative normalization \(N_Q\), and the outgoing \(l=2\) DtN coefficient \(\chi_Q\).  It does not change
      the parent action or the grouped-lane ontology; it works inside the declared reduced branch and tests the next algebraic
      or realization condition.}
  - path: paper/stages/stage_104.tex
    line: 11
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py};
      transcript: \StageFile{scripts/output/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.txt}.  Mathematica
      audit: \StageFile{mathematica/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.wl}; transcript:
      \StageFile{mathematica/output/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.txt}.}'
  - path: paper/stages/stage_104.tex
    line: 13
    role: paper_stage_tex
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
  - path: paper/stages/stage_104.tex
    line: 16
    role: paper_stage_tex
    excerpt: Spherical Hankel DtN expansion gives \(\widehat Y_2^{\rm out}=1+z^2/9+4z^4/81+iz^5/27+\cdots\).
  - path: paper/stages/stage_104.tex
    line: 22
    role: paper_stage_tex
    excerpt: \item Check the product \(\widehat m_0^{\,2}\chi_QN_Q\) keeps source, conservative, and outgoing factors separate.
  - path: paper/stages/stage_104.tex
    line: 24
    role: paper_stage_tex
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
  - path: paper/stages/stage_104.tex
    line: 27
    role: paper_stage_tex
    excerpt: \stagefield{Downstream use}{This stage feeds Stages~107--113, the \(2.5\)PN/4PN outgoing bridge, and later same-charge
      prefactor audits.  When cited downstream, the status tag above must be carried with the result; the card is a derivation
      ledger entry, not an unconditional actual-branch theorem.}
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 104: Exact Outgoing `l=2` DtN Fingerprint'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 5
    role: notes_stage
    excerpt: Compute the canonical compact passive/outgoing quadrupole normalization directly from an explicit outgoing `l=2`
      Dirichlet-to-Neumann model instead of carrying the coefficient `chi_Q` symbolically.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 7
    role: notes_stage
    excerpt: '## Exact outgoing DtN model'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 13
    role: notes_stage
    excerpt: and let the outgoing partial wave be the spherical Hankel mode
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 18
    role: notes_stage
    excerpt: The exact `l=2` outgoing DtN operator is
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 42
    role: notes_stage
    excerpt: '## Normalized outgoing quadrupole admittance'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 44
    role: notes_stage
    excerpt: Define the normalized outgoing branch by
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 83
    role: notes_stage
    excerpt: This is exactly the outgoing compact `l=2` fingerprint used earlier in the reduced 2.5PN program, but now derived
      directly from the explicit DtN model.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 87
    role: notes_stage
    excerpt: The canonical outgoing coefficient is therefore not free. On the exact spherical outgoing `l=2` DtN branch, the
      leading odd quadrupole coefficient is fixed to
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 89
    role: notes_stage
    excerpt: \Gamma_{5,\rm can}^{\rm DtN}=\frac{a^5}{27c_s^5}.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 92
    role: notes_stage
    excerpt: So any later normalization mismatch must come from branch selection or source normalization, not from ambiguity
      in the canonical outgoing `l=2` DtN model itself.
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 18
    role: sympy_script
    excerpt: banner("STAGE 104 — EXACT OUTGOING l=2 DtN FINGERPRINT")
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 37
    role: sympy_script
    excerpt: expect_zero("static DtN slot", sp.limit(Lam, z, 0) + 3)
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 55
    role: sympy_script
    excerpt: print("  The explicit outgoing spherical l=2 DtN model reproduces the canonical compact")
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 56
    role: sympy_script
    excerpt: print("  quadrupole fingerprint exactly, including the odd coefficient a^5/(27 c_s^5).")
  - path: paper/stages/stage_105.tex
    line: 1
    role: paper_stage_tex
    excerpt: '\section[Stage~105]{Stage~105: Exact Fixing of chi-Q}'
  - path: paper/stages/stage_105.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{This stage imports the cleaned conservative isotropic branch, the source-map factor \(\widehat
      m_0\), the conservative normalization \(N_Q\), and the outgoing \(l=2\) DtN coefficient \(\chi_Q\).  It does not change
      the parent action or the grouped-lane ontology; it works inside the declared reduced branch and tests the next algebraic
      or realization condition.}
  - path: paper/stages/stage_105.tex
    line: 11
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py};
      transcript: \StageFile{scripts/output/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.txt}.  Mathematica
      audit: \StageFile{mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl}; transcript:
      \StageFile{mathematica/output/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.txt}.}'
  - path: paper/stages/stage_105.tex
    line: 13
    role: paper_stage_tex
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
  - path: paper/stages/stage_105.tex
    line: 16
    role: paper_stage_tex
    excerpt: Matching the canonical grouped module to the exact outgoing DtN branch gives \(\chi_Q=1\).
  - path: paper/stages/stage_105.tex
    line: 22
    role: paper_stage_tex
    excerpt: \item Check the product \(\widehat m_0^{\,2}\chi_QN_Q\) keeps source, conservative, and outgoing factors separate.
  - path: paper/stages/stage_105.tex
    line: 24
    role: paper_stage_tex
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
  - path: paper/stages/stage_105.tex
    line: 27
    role: paper_stage_tex
    excerpt: \stagefield{Downstream use}{This stage feeds Stages~107--113, the \(2.5\)PN/4PN outgoing bridge, and later same-charge
      prefactor audits.  When cited downstream, the status tag above must be carried with the result; the card is a derivation
      ledger entry, not an unconditional actual-branch theorem.}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 105: Exact Fixing of `chi_Q`'
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 5
    role: notes_stage
    excerpt: Use the explicit DtN fingerprint to determine the last reduced 2.5PN scalar
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 7
    role: notes_stage
    excerpt: \chi_Q
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 22
    role: notes_stage
    excerpt: -i\,\chi_Q\,\sigma_Q^{\rm can}\,\omega^5
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 34
    role: notes_stage
    excerpt: \sigma_Q^{\rm can}:=\frac{9}{8\Omega_Q^5}.
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 39
    role: notes_stage
    excerpt: \sigma_Q^{\rm can}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 53
    role: notes_stage
    excerpt: +i\,\chi_Q\,\frac{a^5\omega^5}{27c_s^5}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 57
    role: notes_stage
    excerpt: From Stage 87, the explicit outgoing DtN branch gives
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 67
    role: notes_stage
    excerpt: Matching the \(O(\omega^5)\) coefficient yields
  origin_claims: []
  constraints: []
  graph_context:
    contexts: []
    graph_gaps:
    - source: paper/stages/stage_104.tex
      attempted_sources:
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: paper/stages/stage_105.tex
      attempted_sources:
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
  provenance_findings:
  - type: graph_gap
    severity: dry_run_non_binding
    summary: The atlas graph wrapper returned no source nodes for one or more primary sources.
    gaps:
    - source: paper/stages/stage_104.tex
      attempted_sources:
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: paper/stages/stage_105.tex
      attempted_sources:
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
- schema_version: 1
  candidate_id: dryrun_stages_003_104_105_chi_q_outgoing_dtn
  parameter_name: matched_fingerprint_value
  dry_run: true
  dry_run_id: stages_003_104_105
  non_binding: true
  generated_content_kind: mechanical_evidence_bundle
  agent_prompt_path: redteam_adversarial/tmp_prompts/dryrun_stages_003_104_105_chi_q_outgoing_dtn__phase_b__matched_fingerprint_value.md
  agent_synthesis_required: true
  dry_run_fixture: true
  fixture_note: Dry-run provenance stores only opened source evidence and graph context. Origin claims, constraints, and contradiction
    findings must be produced by the Phase B agent path.
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '104'
  - '105'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_104.tex
    line: 1
    role: paper_stage_tex
    excerpt: '\section[Stage~104]{Stage~104: Exact Outgoing l=2 DtN Fingerprint}'
  - path: paper/stages/stage_104.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{This stage imports the cleaned conservative isotropic branch, the source-map factor \(\widehat
      m_0\), the conservative normalization \(N_Q\), and the outgoing \(l=2\) DtN coefficient \(\chi_Q\).  It does not change
      the parent action or the grouped-lane ontology; it works inside the declared reduced branch and tests the next algebraic
      or realization condition.}
  - path: paper/stages/stage_104.tex
    line: 11
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py};
      transcript: \StageFile{scripts/output/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.txt}.  Mathematica
      audit: \StageFile{mathematica/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.wl}; transcript:
      \StageFile{mathematica/output/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.txt}.}'
  - path: paper/stages/stage_104.tex
    line: 13
    role: paper_stage_tex
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
  - path: paper/stages/stage_104.tex
    line: 16
    role: paper_stage_tex
    excerpt: Spherical Hankel DtN expansion gives \(\widehat Y_2^{\rm out}=1+z^2/9+4z^4/81+iz^5/27+\cdots\).
  - path: paper/stages/stage_104.tex
    line: 22
    role: paper_stage_tex
    excerpt: \item Check the product \(\widehat m_0^{\,2}\chi_QN_Q\) keeps source, conservative, and outgoing factors separate.
  - path: paper/stages/stage_104.tex
    line: 24
    role: paper_stage_tex
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
  - path: paper/stages/stage_104.tex
    line: 27
    role: paper_stage_tex
    excerpt: \stagefield{Downstream use}{This stage feeds Stages~107--113, the \(2.5\)PN/4PN outgoing bridge, and later same-charge
      prefactor audits.  When cited downstream, the status tag above must be carried with the result; the card is a derivation
      ledger entry, not an unconditional actual-branch theorem.}
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 104: Exact Outgoing `l=2` DtN Fingerprint'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 5
    role: notes_stage
    excerpt: Compute the canonical compact passive/outgoing quadrupole normalization directly from an explicit outgoing `l=2`
      Dirichlet-to-Neumann model instead of carrying the coefficient `chi_Q` symbolically.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 7
    role: notes_stage
    excerpt: '## Exact outgoing DtN model'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 13
    role: notes_stage
    excerpt: and let the outgoing partial wave be the spherical Hankel mode
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 18
    role: notes_stage
    excerpt: The exact `l=2` outgoing DtN operator is
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 42
    role: notes_stage
    excerpt: '## Normalized outgoing quadrupole admittance'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 44
    role: notes_stage
    excerpt: Define the normalized outgoing branch by
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 83
    role: notes_stage
    excerpt: This is exactly the outgoing compact `l=2` fingerprint used earlier in the reduced 2.5PN program, but now derived
      directly from the explicit DtN model.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 87
    role: notes_stage
    excerpt: The canonical outgoing coefficient is therefore not free. On the exact spherical outgoing `l=2` DtN branch, the
      leading odd quadrupole coefficient is fixed to
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 89
    role: notes_stage
    excerpt: \Gamma_{5,\rm can}^{\rm DtN}=\frac{a^5}{27c_s^5}.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 92
    role: notes_stage
    excerpt: So any later normalization mismatch must come from branch selection or source normalization, not from ambiguity
      in the canonical outgoing `l=2` DtN model itself.
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 18
    role: sympy_script
    excerpt: banner("STAGE 104 — EXACT OUTGOING l=2 DtN FINGERPRINT")
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 37
    role: sympy_script
    excerpt: expect_zero("static DtN slot", sp.limit(Lam, z, 0) + 3)
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 55
    role: sympy_script
    excerpt: print("  The explicit outgoing spherical l=2 DtN model reproduces the canonical compact")
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 56
    role: sympy_script
    excerpt: print("  quadrupole fingerprint exactly, including the odd coefficient a^5/(27 c_s^5).")
  - path: paper/stages/stage_105.tex
    line: 1
    role: paper_stage_tex
    excerpt: '\section[Stage~105]{Stage~105: Exact Fixing of chi-Q}'
  - path: paper/stages/stage_105.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{This stage imports the cleaned conservative isotropic branch, the source-map factor \(\widehat
      m_0\), the conservative normalization \(N_Q\), and the outgoing \(l=2\) DtN coefficient \(\chi_Q\).  It does not change
      the parent action or the grouped-lane ontology; it works inside the declared reduced branch and tests the next algebraic
      or realization condition.}
  - path: paper/stages/stage_105.tex
    line: 11
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py};
      transcript: \StageFile{scripts/output/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.txt}.  Mathematica
      audit: \StageFile{mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl}; transcript:
      \StageFile{mathematica/output/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.txt}.}'
  - path: paper/stages/stage_105.tex
    line: 13
    role: paper_stage_tex
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
  - path: paper/stages/stage_105.tex
    line: 16
    role: paper_stage_tex
    excerpt: Matching the canonical grouped module to the exact outgoing DtN branch gives \(\chi_Q=1\).
  - path: paper/stages/stage_105.tex
    line: 22
    role: paper_stage_tex
    excerpt: \item Check the product \(\widehat m_0^{\,2}\chi_QN_Q\) keeps source, conservative, and outgoing factors separate.
  - path: paper/stages/stage_105.tex
    line: 24
    role: paper_stage_tex
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
  - path: paper/stages/stage_105.tex
    line: 27
    role: paper_stage_tex
    excerpt: \stagefield{Downstream use}{This stage feeds Stages~107--113, the \(2.5\)PN/4PN outgoing bridge, and later same-charge
      prefactor audits.  When cited downstream, the status tag above must be carried with the result; the card is a derivation
      ledger entry, not an unconditional actual-branch theorem.}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 105: Exact Fixing of `chi_Q`'
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 5
    role: notes_stage
    excerpt: Use the explicit DtN fingerprint to determine the last reduced 2.5PN scalar
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 7
    role: notes_stage
    excerpt: \chi_Q
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 22
    role: notes_stage
    excerpt: -i\,\chi_Q\,\sigma_Q^{\rm can}\,\omega^5
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 34
    role: notes_stage
    excerpt: \sigma_Q^{\rm can}:=\frac{9}{8\Omega_Q^5}.
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 39
    role: notes_stage
    excerpt: \sigma_Q^{\rm can}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 53
    role: notes_stage
    excerpt: +i\,\chi_Q\,\frac{a^5\omega^5}{27c_s^5}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 57
    role: notes_stage
    excerpt: From Stage 87, the explicit outgoing DtN branch gives
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 67
    role: notes_stage
    excerpt: Matching the \(O(\omega^5)\) coefficient yields
  origin_claims: []
  constraints: []
  graph_context:
    contexts: []
    graph_gaps:
    - source: paper/stages/stage_104.tex
      attempted_sources:
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: paper/stages/stage_105.tex
      attempted_sources:
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
  provenance_findings:
  - type: graph_gap
    severity: dry_run_non_binding
    summary: The atlas graph wrapper returned no source nodes for one or more primary sources.
    gaps:
    - source: paper/stages/stage_104.tex
      attempted_sources:
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: paper/stages/stage_105.tex
      attempted_sources:
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
- schema_version: 1
  candidate_id: dryrun_stages_003_104_105_chi_q_outgoing_dtn
  parameter_name: outgoing_dtn_fingerprint
  dry_run: true
  dry_run_id: stages_003_104_105
  non_binding: true
  generated_content_kind: mechanical_evidence_bundle
  agent_prompt_path: redteam_adversarial/tmp_prompts/dryrun_stages_003_104_105_chi_q_outgoing_dtn__phase_b__outgoing_dtn_fingerprint.md
  agent_synthesis_required: true
  dry_run_fixture: true
  fixture_note: Dry-run provenance stores only opened source evidence and graph context. Origin claims, constraints, and contradiction
    findings must be produced by the Phase B agent path.
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '104'
  - '105'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_104.tex
    line: 1
    role: paper_stage_tex
    excerpt: '\section[Stage~104]{Stage~104: Exact Outgoing l=2 DtN Fingerprint}'
  - path: paper/stages/stage_104.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{This stage imports the cleaned conservative isotropic branch, the source-map factor \(\widehat
      m_0\), the conservative normalization \(N_Q\), and the outgoing \(l=2\) DtN coefficient \(\chi_Q\).  It does not change
      the parent action or the grouped-lane ontology; it works inside the declared reduced branch and tests the next algebraic
      or realization condition.}
  - path: paper/stages/stage_104.tex
    line: 11
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py};
      transcript: \StageFile{scripts/output/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.txt}.  Mathematica
      audit: \StageFile{mathematica/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.wl}; transcript:
      \StageFile{mathematica/output/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.txt}.}'
  - path: paper/stages/stage_104.tex
    line: 13
    role: paper_stage_tex
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
  - path: paper/stages/stage_104.tex
    line: 16
    role: paper_stage_tex
    excerpt: Spherical Hankel DtN expansion gives \(\widehat Y_2^{\rm out}=1+z^2/9+4z^4/81+iz^5/27+\cdots\).
  - path: paper/stages/stage_104.tex
    line: 22
    role: paper_stage_tex
    excerpt: \item Check the product \(\widehat m_0^{\,2}\chi_QN_Q\) keeps source, conservative, and outgoing factors separate.
  - path: paper/stages/stage_104.tex
    line: 24
    role: paper_stage_tex
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
  - path: paper/stages/stage_104.tex
    line: 27
    role: paper_stage_tex
    excerpt: \stagefield{Downstream use}{This stage feeds Stages~107--113, the \(2.5\)PN/4PN outgoing bridge, and later same-charge
      prefactor audits.  When cited downstream, the status tag above must be carried with the result; the card is a derivation
      ledger entry, not an unconditional actual-branch theorem.}
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 104: Exact Outgoing `l=2` DtN Fingerprint'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 5
    role: notes_stage
    excerpt: Compute the canonical compact passive/outgoing quadrupole normalization directly from an explicit outgoing `l=2`
      Dirichlet-to-Neumann model instead of carrying the coefficient `chi_Q` symbolically.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 7
    role: notes_stage
    excerpt: '## Exact outgoing DtN model'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 13
    role: notes_stage
    excerpt: and let the outgoing partial wave be the spherical Hankel mode
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 18
    role: notes_stage
    excerpt: The exact `l=2` outgoing DtN operator is
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 42
    role: notes_stage
    excerpt: '## Normalized outgoing quadrupole admittance'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 44
    role: notes_stage
    excerpt: Define the normalized outgoing branch by
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 83
    role: notes_stage
    excerpt: This is exactly the outgoing compact `l=2` fingerprint used earlier in the reduced 2.5PN program, but now derived
      directly from the explicit DtN model.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 87
    role: notes_stage
    excerpt: The canonical outgoing coefficient is therefore not free. On the exact spherical outgoing `l=2` DtN branch, the
      leading odd quadrupole coefficient is fixed to
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 89
    role: notes_stage
    excerpt: \Gamma_{5,\rm can}^{\rm DtN}=\frac{a^5}{27c_s^5}.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 92
    role: notes_stage
    excerpt: So any later normalization mismatch must come from branch selection or source normalization, not from ambiguity
      in the canonical outgoing `l=2` DtN model itself.
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 18
    role: sympy_script
    excerpt: banner("STAGE 104 — EXACT OUTGOING l=2 DtN FINGERPRINT")
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 37
    role: sympy_script
    excerpt: expect_zero("static DtN slot", sp.limit(Lam, z, 0) + 3)
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 55
    role: sympy_script
    excerpt: print("  The explicit outgoing spherical l=2 DtN model reproduces the canonical compact")
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 56
    role: sympy_script
    excerpt: print("  quadrupole fingerprint exactly, including the odd coefficient a^5/(27 c_s^5).")
  - path: paper/stages/stage_105.tex
    line: 1
    role: paper_stage_tex
    excerpt: '\section[Stage~105]{Stage~105: Exact Fixing of chi-Q}'
  - path: paper/stages/stage_105.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{This stage imports the cleaned conservative isotropic branch, the source-map factor \(\widehat
      m_0\), the conservative normalization \(N_Q\), and the outgoing \(l=2\) DtN coefficient \(\chi_Q\).  It does not change
      the parent action or the grouped-lane ontology; it works inside the declared reduced branch and tests the next algebraic
      or realization condition.}
  - path: paper/stages/stage_105.tex
    line: 11
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py};
      transcript: \StageFile{scripts/output/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.txt}.  Mathematica
      audit: \StageFile{mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl}; transcript:
      \StageFile{mathematica/output/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.txt}.}'
  - path: paper/stages/stage_105.tex
    line: 13
    role: paper_stage_tex
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
  - path: paper/stages/stage_105.tex
    line: 16
    role: paper_stage_tex
    excerpt: Matching the canonical grouped module to the exact outgoing DtN branch gives \(\chi_Q=1\).
  - path: paper/stages/stage_105.tex
    line: 22
    role: paper_stage_tex
    excerpt: \item Check the product \(\widehat m_0^{\,2}\chi_QN_Q\) keeps source, conservative, and outgoing factors separate.
  - path: paper/stages/stage_105.tex
    line: 24
    role: paper_stage_tex
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
  - path: paper/stages/stage_105.tex
    line: 27
    role: paper_stage_tex
    excerpt: \stagefield{Downstream use}{This stage feeds Stages~107--113, the \(2.5\)PN/4PN outgoing bridge, and later same-charge
      prefactor audits.  When cited downstream, the status tag above must be carried with the result; the card is a derivation
      ledger entry, not an unconditional actual-branch theorem.}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 105: Exact Fixing of `chi_Q`'
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 5
    role: notes_stage
    excerpt: Use the explicit DtN fingerprint to determine the last reduced 2.5PN scalar
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 7
    role: notes_stage
    excerpt: \chi_Q
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 22
    role: notes_stage
    excerpt: -i\,\chi_Q\,\sigma_Q^{\rm can}\,\omega^5
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 34
    role: notes_stage
    excerpt: \sigma_Q^{\rm can}:=\frac{9}{8\Omega_Q^5}.
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 39
    role: notes_stage
    excerpt: \sigma_Q^{\rm can}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 53
    role: notes_stage
    excerpt: +i\,\chi_Q\,\frac{a^5\omega^5}{27c_s^5}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 57
    role: notes_stage
    excerpt: From Stage 87, the explicit outgoing DtN branch gives
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 67
    role: notes_stage
    excerpt: Matching the \(O(\omega^5)\) coefficient yields
  origin_claims: []
  constraints: []
  graph_context:
    contexts: []
    graph_gaps:
    - source: paper/stages/stage_104.tex
      attempted_sources:
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: paper/stages/stage_105.tex
      attempted_sources:
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
  provenance_findings:
  - type: graph_gap
    severity: dry_run_non_binding
    summary: The atlas graph wrapper returned no source nodes for one or more primary sources.
    gaps:
    - source: paper/stages/stage_104.tex
      attempted_sources:
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: paper/stages/stage_105.tex
      attempted_sources:
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
- schema_version: 1
  candidate_id: dryrun_stages_003_104_105_chi_q_outgoing_dtn
  parameter_name: sigma_Q
  dry_run: true
  dry_run_id: stages_003_104_105
  non_binding: true
  generated_content_kind: mechanical_evidence_bundle
  agent_prompt_path: redteam_adversarial/tmp_prompts/dryrun_stages_003_104_105_chi_q_outgoing_dtn__phase_b__sigma_q.md
  agent_synthesis_required: true
  dry_run_fixture: true
  fixture_note: Dry-run provenance stores only opened source evidence and graph context. Origin claims, constraints, and contradiction
    findings must be produced by the Phase B agent path.
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '104'
  - '105'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_104.tex
    line: 1
    role: paper_stage_tex
    excerpt: '\section[Stage~104]{Stage~104: Exact Outgoing l=2 DtN Fingerprint}'
  - path: paper/stages/stage_104.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{This stage imports the cleaned conservative isotropic branch, the source-map factor \(\widehat
      m_0\), the conservative normalization \(N_Q\), and the outgoing \(l=2\) DtN coefficient \(\chi_Q\).  It does not change
      the parent action or the grouped-lane ontology; it works inside the declared reduced branch and tests the next algebraic
      or realization condition.}
  - path: paper/stages/stage_104.tex
    line: 11
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py};
      transcript: \StageFile{scripts/output/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.txt}.  Mathematica
      audit: \StageFile{mathematica/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.wl}; transcript:
      \StageFile{mathematica/output/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.txt}.}'
  - path: paper/stages/stage_104.tex
    line: 13
    role: paper_stage_tex
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
  - path: paper/stages/stage_104.tex
    line: 16
    role: paper_stage_tex
    excerpt: Spherical Hankel DtN expansion gives \(\widehat Y_2^{\rm out}=1+z^2/9+4z^4/81+iz^5/27+\cdots\).
  - path: paper/stages/stage_104.tex
    line: 22
    role: paper_stage_tex
    excerpt: \item Check the product \(\widehat m_0^{\,2}\chi_QN_Q\) keeps source, conservative, and outgoing factors separate.
  - path: paper/stages/stage_104.tex
    line: 24
    role: paper_stage_tex
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
  - path: paper/stages/stage_104.tex
    line: 27
    role: paper_stage_tex
    excerpt: \stagefield{Downstream use}{This stage feeds Stages~107--113, the \(2.5\)PN/4PN outgoing bridge, and later same-charge
      prefactor audits.  When cited downstream, the status tag above must be carried with the result; the card is a derivation
      ledger entry, not an unconditional actual-branch theorem.}
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 104: Exact Outgoing `l=2` DtN Fingerprint'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 5
    role: notes_stage
    excerpt: Compute the canonical compact passive/outgoing quadrupole normalization directly from an explicit outgoing `l=2`
      Dirichlet-to-Neumann model instead of carrying the coefficient `chi_Q` symbolically.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 7
    role: notes_stage
    excerpt: '## Exact outgoing DtN model'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 13
    role: notes_stage
    excerpt: and let the outgoing partial wave be the spherical Hankel mode
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 18
    role: notes_stage
    excerpt: The exact `l=2` outgoing DtN operator is
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 42
    role: notes_stage
    excerpt: '## Normalized outgoing quadrupole admittance'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 44
    role: notes_stage
    excerpt: Define the normalized outgoing branch by
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 83
    role: notes_stage
    excerpt: This is exactly the outgoing compact `l=2` fingerprint used earlier in the reduced 2.5PN program, but now derived
      directly from the explicit DtN model.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 87
    role: notes_stage
    excerpt: The canonical outgoing coefficient is therefore not free. On the exact spherical outgoing `l=2` DtN branch, the
      leading odd quadrupole coefficient is fixed to
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 89
    role: notes_stage
    excerpt: \Gamma_{5,\rm can}^{\rm DtN}=\frac{a^5}{27c_s^5}.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 92
    role: notes_stage
    excerpt: So any later normalization mismatch must come from branch selection or source normalization, not from ambiguity
      in the canonical outgoing `l=2` DtN model itself.
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 18
    role: sympy_script
    excerpt: banner("STAGE 104 — EXACT OUTGOING l=2 DtN FINGERPRINT")
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 37
    role: sympy_script
    excerpt: expect_zero("static DtN slot", sp.limit(Lam, z, 0) + 3)
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 55
    role: sympy_script
    excerpt: print("  The explicit outgoing spherical l=2 DtN model reproduces the canonical compact")
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 56
    role: sympy_script
    excerpt: print("  quadrupole fingerprint exactly, including the odd coefficient a^5/(27 c_s^5).")
  - path: paper/stages/stage_105.tex
    line: 1
    role: paper_stage_tex
    excerpt: '\section[Stage~105]{Stage~105: Exact Fixing of chi-Q}'
  - path: paper/stages/stage_105.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{This stage imports the cleaned conservative isotropic branch, the source-map factor \(\widehat
      m_0\), the conservative normalization \(N_Q\), and the outgoing \(l=2\) DtN coefficient \(\chi_Q\).  It does not change
      the parent action or the grouped-lane ontology; it works inside the declared reduced branch and tests the next algebraic
      or realization condition.}
  - path: paper/stages/stage_105.tex
    line: 11
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py};
      transcript: \StageFile{scripts/output/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.txt}.  Mathematica
      audit: \StageFile{mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl}; transcript:
      \StageFile{mathematica/output/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.txt}.}'
  - path: paper/stages/stage_105.tex
    line: 13
    role: paper_stage_tex
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
  - path: paper/stages/stage_105.tex
    line: 16
    role: paper_stage_tex
    excerpt: Matching the canonical grouped module to the exact outgoing DtN branch gives \(\chi_Q=1\).
  - path: paper/stages/stage_105.tex
    line: 22
    role: paper_stage_tex
    excerpt: \item Check the product \(\widehat m_0^{\,2}\chi_QN_Q\) keeps source, conservative, and outgoing factors separate.
  - path: paper/stages/stage_105.tex
    line: 24
    role: paper_stage_tex
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
  - path: paper/stages/stage_105.tex
    line: 27
    role: paper_stage_tex
    excerpt: \stagefield{Downstream use}{This stage feeds Stages~107--113, the \(2.5\)PN/4PN outgoing bridge, and later same-charge
      prefactor audits.  When cited downstream, the status tag above must be carried with the result; the card is a derivation
      ledger entry, not an unconditional actual-branch theorem.}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 105: Exact Fixing of `chi_Q`'
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 5
    role: notes_stage
    excerpt: Use the explicit DtN fingerprint to determine the last reduced 2.5PN scalar
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 7
    role: notes_stage
    excerpt: \chi_Q
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 22
    role: notes_stage
    excerpt: -i\,\chi_Q\,\sigma_Q^{\rm can}\,\omega^5
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 34
    role: notes_stage
    excerpt: \sigma_Q^{\rm can}:=\frac{9}{8\Omega_Q^5}.
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 39
    role: notes_stage
    excerpt: \sigma_Q^{\rm can}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 53
    role: notes_stage
    excerpt: +i\,\chi_Q\,\frac{a^5\omega^5}{27c_s^5}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 57
    role: notes_stage
    excerpt: From Stage 87, the explicit outgoing DtN branch gives
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 67
    role: notes_stage
    excerpt: Matching the \(O(\omega^5)\) coefficient yields
  origin_claims: []
  constraints: []
  graph_context:
    contexts: []
    graph_gaps:
    - source: paper/stages/stage_104.tex
      attempted_sources:
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: paper/stages/stage_105.tex
      attempted_sources:
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
  provenance_findings:
  - type: graph_gap
    severity: dry_run_non_binding
    summary: The atlas graph wrapper returned no source nodes for one or more primary sources.
    gaps:
    - source: paper/stages/stage_104.tex
      attempted_sources:
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: paper/stages/stage_105.tex
      attempted_sources:
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
- schema_version: 1
  candidate_id: dryrun_stages_003_104_105_chi_q_outgoing_dtn
  parameter_name: xi_Q
  dry_run: true
  dry_run_id: stages_003_104_105
  non_binding: true
  generated_content_kind: mechanical_evidence_bundle
  agent_prompt_path: redteam_adversarial/tmp_prompts/dryrun_stages_003_104_105_chi_q_outgoing_dtn__phase_b__xi_q.md
  agent_synthesis_required: true
  dry_run_fixture: true
  fixture_note: Dry-run provenance stores only opened source evidence and graph context. Origin claims, constraints, and contradiction
    findings must be produced by the Phase B agent path.
  taxonomy:
    ledger_scope: parameter-value provenance only
    excluded_layers:
    - file/source provenance
    - result/stage provenance without value genealogy
  anchor_stages:
  - '104'
  - '105'
  notes_source_opened: true
  source_evidence:
  - path: paper/stages/stage_104.tex
    line: 1
    role: paper_stage_tex
    excerpt: '\section[Stage~104]{Stage~104: Exact Outgoing l=2 DtN Fingerprint}'
  - path: paper/stages/stage_104.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{This stage imports the cleaned conservative isotropic branch, the source-map factor \(\widehat
      m_0\), the conservative normalization \(N_Q\), and the outgoing \(l=2\) DtN coefficient \(\chi_Q\).  It does not change
      the parent action or the grouped-lane ontology; it works inside the declared reduced branch and tests the next algebraic
      or realization condition.}
  - path: paper/stages/stage_104.tex
    line: 11
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py};
      transcript: \StageFile{scripts/output/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.txt}.  Mathematica
      audit: \StageFile{mathematica/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.wl}; transcript:
      \StageFile{mathematica/output/moving_throat_pde_stage104_outgoing_dtn_fingerprint_mathematica_audit.txt}.}'
  - path: paper/stages/stage_104.tex
    line: 13
    role: paper_stage_tex
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
  - path: paper/stages/stage_104.tex
    line: 16
    role: paper_stage_tex
    excerpt: Spherical Hankel DtN expansion gives \(\widehat Y_2^{\rm out}=1+z^2/9+4z^4/81+iz^5/27+\cdots\).
  - path: paper/stages/stage_104.tex
    line: 22
    role: paper_stage_tex
    excerpt: \item Check the product \(\widehat m_0^{\,2}\chi_QN_Q\) keeps source, conservative, and outgoing factors separate.
  - path: paper/stages/stage_104.tex
    line: 24
    role: paper_stage_tex
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
  - path: paper/stages/stage_104.tex
    line: 27
    role: paper_stage_tex
    excerpt: \stagefield{Downstream use}{This stage feeds Stages~107--113, the \(2.5\)PN/4PN outgoing bridge, and later same-charge
      prefactor audits.  When cited downstream, the status tag above must be carried with the result; the card is a derivation
      ledger entry, not an unconditional actual-branch theorem.}
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 104: Exact Outgoing `l=2` DtN Fingerprint'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 5
    role: notes_stage
    excerpt: Compute the canonical compact passive/outgoing quadrupole normalization directly from an explicit outgoing `l=2`
      Dirichlet-to-Neumann model instead of carrying the coefficient `chi_Q` symbolically.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 7
    role: notes_stage
    excerpt: '## Exact outgoing DtN model'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 13
    role: notes_stage
    excerpt: and let the outgoing partial wave be the spherical Hankel mode
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 18
    role: notes_stage
    excerpt: The exact `l=2` outgoing DtN operator is
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 42
    role: notes_stage
    excerpt: '## Normalized outgoing quadrupole admittance'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 44
    role: notes_stage
    excerpt: Define the normalized outgoing branch by
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 83
    role: notes_stage
    excerpt: This is exactly the outgoing compact `l=2` fingerprint used earlier in the reduced 2.5PN program, but now derived
      directly from the explicit DtN model.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 87
    role: notes_stage
    excerpt: The canonical outgoing coefficient is therefore not free. On the exact spherical outgoing `l=2` DtN branch, the
      leading odd quadrupole coefficient is fixed to
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 89
    role: notes_stage
    excerpt: \Gamma_{5,\rm can}^{\rm DtN}=\frac{a^5}{27c_s^5}.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 92
    role: notes_stage
    excerpt: So any later normalization mismatch must come from branch selection or source normalization, not from ambiguity
      in the canonical outgoing `l=2` DtN model itself.
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 18
    role: sympy_script
    excerpt: banner("STAGE 104 — EXACT OUTGOING l=2 DtN FINGERPRINT")
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 37
    role: sympy_script
    excerpt: expect_zero("static DtN slot", sp.limit(Lam, z, 0) + 3)
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 55
    role: sympy_script
    excerpt: print("  The explicit outgoing spherical l=2 DtN model reproduces the canonical compact")
  - path: scripts/moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py
    line: 56
    role: sympy_script
    excerpt: print("  quadrupole fingerprint exactly, including the odd coefficient a^5/(27 c_s^5).")
  - path: paper/stages/stage_105.tex
    line: 1
    role: paper_stage_tex
    excerpt: '\section[Stage~105]{Stage~105: Exact Fixing of chi-Q}'
  - path: paper/stages/stage_105.tex
    line: 9
    role: paper_stage_tex
    excerpt: \stagefield{Inputs}{This stage imports the cleaned conservative isotropic branch, the source-map factor \(\widehat
      m_0\), the conservative normalization \(N_Q\), and the outgoing \(l=2\) DtN coefficient \(\chi_Q\).  It does not change
      the parent action or the grouped-lane ontology; it works inside the declared reduced branch and tests the next algebraic
      or realization condition.}
  - path: paper/stages/stage_105.tex
    line: 11
    role: paper_stage_tex
    excerpt: '\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py};
      transcript: \StageFile{scripts/output/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.txt}.  Mathematica
      audit: \StageFile{mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl}; transcript:
      \StageFile{mathematica/output/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.txt}.}'
  - path: paper/stages/stage_105.tex
    line: 13
    role: paper_stage_tex
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
  - path: paper/stages/stage_105.tex
    line: 16
    role: paper_stage_tex
    excerpt: Matching the canonical grouped module to the exact outgoing DtN branch gives \(\chi_Q=1\).
  - path: paper/stages/stage_105.tex
    line: 22
    role: paper_stage_tex
    excerpt: \item Check the product \(\widehat m_0^{\,2}\chi_QN_Q\) keeps source, conservative, and outgoing factors separate.
  - path: paper/stages/stage_105.tex
    line: 24
    role: paper_stage_tex
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
  - path: paper/stages/stage_105.tex
    line: 27
    role: paper_stage_tex
    excerpt: \stagefield{Downstream use}{This stage feeds Stages~107--113, the \(2.5\)PN/4PN outgoing bridge, and later same-charge
      prefactor audits.  When cited downstream, the status tag above must be carried with the result; the card is a derivation
      ledger entry, not an unconditional actual-branch theorem.}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 1
    role: notes_stage
    excerpt: '# Moving-Throat PDE — Stage 105: Exact Fixing of `chi_Q`'
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 5
    role: notes_stage
    excerpt: Use the explicit DtN fingerprint to determine the last reduced 2.5PN scalar
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 7
    role: notes_stage
    excerpt: \chi_Q
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 22
    role: notes_stage
    excerpt: -i\,\chi_Q\,\sigma_Q^{\rm can}\,\omega^5
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 34
    role: notes_stage
    excerpt: \sigma_Q^{\rm can}:=\frac{9}{8\Omega_Q^5}.
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 39
    role: notes_stage
    excerpt: \sigma_Q^{\rm can}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 53
    role: notes_stage
    excerpt: +i\,\chi_Q\,\frac{a^5\omega^5}{27c_s^5}
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 57
    role: notes_stage
    excerpt: From Stage 87, the explicit outgoing DtN branch gives
  - path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 67
    role: notes_stage
    excerpt: Matching the \(O(\omega^5)\) coefficient yields
  origin_claims: []
  constraints: []
  graph_context:
    contexts: []
    graph_gaps:
    - source: paper/stages/stage_104.tex
      attempted_sources:
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: paper/stages/stage_105.tex
      attempted_sources:
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
  provenance_findings:
  - type: graph_gap
    severity: dry_run_non_binding
    summary: The atlas graph wrapper returned no source nodes for one or more primary sources.
    gaps:
    - source: paper/stages/stage_104.tex
      attempted_sources:
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
      - paper/stages/stage_104.tex
      - research/pde_ledger/paper/stages/stage_104.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: paper/stages/stage_105.tex
      attempted_sources:
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
    - source: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      attempted_sources:
      - notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
      - paper/stages/stage_105.tex
      - research/pde_ledger/paper/stages/stage_105.tex
      graph_gap: true
      reason: no atlas node tied to this source path
```

## Benchmarks

External-match checks must use these sourced benchmark entries. Do not adjudicate a match from model memory.

```yaml
- id: dryrun_stages_003_104_105_chi_q_outgoing_dtn__dry_run_benchmark_placeholder
  candidate_id: dryrun_stages_003_104_105_chi_q_outgoing_dtn
  dry_run: true
  dry_run_id: stages_003_104_105
  non_binding: true
  fixture_type: dry-run source-anchor placeholder
  claim: 'DRY-RUN PLACEHOLDER: Phase C must use sourced benchmark entries, not model memory.'
  value: null
  source_type: dry-run project-source anchor list
  source_citations:
  - path: paper/stages/stage_104.tex
    line: 13
    role: paper_stage_tex
    stage: '104'
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
  - path: paper/stages/stage_104.tex
    line: 24
    role: paper_stage_tex
    stage: '104'
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 1
    role: notes_stage
    stage: '104'
    excerpt: '# Moving-Throat PDE — Stage 104: Exact Outgoing `l=2` DtN Fingerprint'
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 83
    role: notes_stage
    stage: '104'
    excerpt: This is exactly the outgoing compact `l=2` fingerprint used earlier in the reduced 2.5PN program, but now derived
      directly from the explicit DtN model.
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 87
    role: notes_stage
    stage: '104'
    excerpt: The canonical outgoing coefficient is therefore not free. On the exact spherical outgoing `l=2` DtN branch, the
      leading odd quadrupole coefficient is fixed to
  - path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 92
    role: notes_stage
    stage: '104'
    excerpt: So any later normalization mismatch must come from branch selection or source normalization, not from ambiguity
      in the canonical outgoing `l=2` DtN model itself.
  - path: paper/stages/stage_105.tex
    line: 13
    role: paper_stage_tex
    stage: '105'
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
  - path: paper/stages/stage_105.tex
    line: 16
    role: paper_stage_tex
    stage: '105'
    excerpt: Matching the canonical grouped module to the exact outgoing DtN branch gives \(\chi_Q=1\).
  requires_agent_sourcing: true
  external_match_policy: For real campaign use, replace this placeholder with an agent-built sourced benchmark entry.
```

## Graph Context

```yaml
contexts: []
graph_gaps:
- source: paper/stages/stage_104.tex
  attempted_sources:
  - paper/stages/stage_104.tex
  - research/pde_ledger/paper/stages/stage_104.tex
  graph_gap: true
  reason: no atlas node tied to this source path
- source: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
  attempted_sources:
  - notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
  - research/pde_ledger/notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
  - paper/stages/stage_104.tex
  - research/pde_ledger/paper/stages/stage_104.tex
  graph_gap: true
  reason: no atlas node tied to this source path
- source: paper/stages/stage_105.tex
  attempted_sources:
  - paper/stages/stage_105.tex
  - research/pde_ledger/paper/stages/stage_105.tex
  graph_gap: true
  reason: no atlas node tied to this source path
- source: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
  attempted_sources:
  - notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
  - research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
  - paper/stages/stage_105.tex
  - research/pde_ledger/paper/stages/stage_105.tex
  graph_gap: true
  reason: no atlas node tied to this source path
```

Output a markdown adversarial report only when this is not a dry run. For dry runs, stop after confirming the prompt has enough context and record no verdict.
