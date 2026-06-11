# Phase A Completeness Critic

You are the completeness critic. This pass runs only after the four blind modalities have emitted independent YAML fragments.

Stages: `003, 104, 105`

Modality fragment paths:
- redteam_adversarial/phase_a_fragments/stages_003_104_105_numeric_literal.yaml
- redteam_adversarial/phase_a_fragments/stages_003_104_105_claim_label.yaml
- redteam_adversarial/phase_a_fragments/stages_003_104_105_graph.yaml
- redteam_adversarial/phase_a_fragments/stages_003_104_105_existing_provenance.yaml

Union output:

```yaml
candidates:
- id: dryrun_stages_003_104_105_chi_q_outgoing_dtn
  candidate_key: chi_q_outgoing_dtn
  dry_run: true
  dry_run_id: stages_003_104_105
  anchor_stages:
  - '104'
  - '105'
  parameter_names:
  - Gamma_5
  - T_dtn
  - chi_Q
  - matched_fingerprint_value
  - outgoing_dtn_fingerprint
  - sigma_Q
  - xi_Q
  citations:
  - &id001
    path: paper/stages/stage_104.tex
    line: 13
    role: paper_stage_tex
    stage: '104'
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
  - &id002
    path: paper/stages/stage_104.tex
    line: 24
    role: paper_stage_tex
    stage: '104'
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
  - &id003
    path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 1
    role: notes_stage
    stage: '104'
    excerpt: '# Moving-Throat PDE — Stage 104: Exact Outgoing `l=2` DtN Fingerprint'
  - &id004
    path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 83
    role: notes_stage
    stage: '104'
    excerpt: This is exactly the outgoing compact `l=2` fingerprint used earlier in the reduced 2.5PN program, but now derived
      directly from the explicit DtN model.
  - &id005
    path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 87
    role: notes_stage
    stage: '104'
    excerpt: The canonical outgoing coefficient is therefore not free. On the exact spherical outgoing `l=2` DtN branch, the
      leading odd quadrupole coefficient is fixed to
  - &id006
    path: notes/stages/moving_throat_pde_stage104_outgoing_dtn_fingerprint.md
    line: 92
    role: notes_stage
    stage: '104'
    excerpt: So any later normalization mismatch must come from branch selection or source normalization, not from ambiguity
      in the canonical outgoing `l=2` DtN model itself.
  - &id007
    path: paper/stages/stage_105.tex
    line: 13
    role: paper_stage_tex
    stage: '105'
    excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization} block.  The
      computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\).}
  - &id008
    path: paper/stages/stage_105.tex
    line: 16
    role: paper_stage_tex
    stage: '105'
    excerpt: Matching the canonical grouped module to the exact outgoing DtN branch gives \(\chi_Q=1\).
  - &id009
    path: paper/stages/stage_105.tex
    line: 24
    role: paper_stage_tex
    stage: '105'
    excerpt: \item Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\) expansion.
  - &id010
    path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 5
    role: notes_stage
    stage: '105'
    excerpt: Use the explicit DtN fingerprint to determine the last reduced 2.5PN scalar
  - &id011
    path: notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
    line: 98
    role: notes_stage
    stage: '105'
    excerpt: That means the only reduced 2.5PN obstruction left after the present calculation is a deviation of the actual
      moving-throat DtN branch from the canonical outgoing `l=2` coefficient \(\xi_Q=1\).
  - &id012
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
  - &id013
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
  - &id014
    path: redteam/pass2/reports/stage_104.md
    line: 52
    role: pass2_stage_report
    stage: '104'
    excerpt: '`paper_alignment: aligned`. The card''s stated audit target (the `\widehat Y_2^out` fingerprint expansion) is
      fully and non-tautologically exercised by both engines. `chi_Q=1` is fixed by comparing this fingerprint to an earlier
      stage''s retarded one-pole form in the appendix narrative; it is not an in-script deliverable of stage 104, so its absence
      from the script is not `script_missing_paper_claim`.'
  - &id015
    path: redteam/pass2/reports/stage_104.md
    line: 59
    role: pass2_stage_report
    stage: '104'
    excerpt: '| A2 | sympy | 38 | `Y[z^2] − 1/9 == 0` | Y fingerprint | yes |'
  - &id016
    path: redteam/pass2/reports/stage_104.md
    line: 60
    role: pass2_stage_report
    stage: '104'
    excerpt: '| A3 | sympy | 39 | `Y[z^4] − 4/81 == 0` | Y fingerprint | yes |'
  - &id017
    path: redteam/pass2/reports/stage_104.md
    line: 61
    role: pass2_stage_report
    stage: '104'
    excerpt: '| A4 | sympy | 40 | `Y[z^5]/I − 1/27 == 0` | Y fingerprint (odd) | yes |'
  - &id018
    path: redteam/pass2/reports/stage_104.md
    line: 62
    role: pass2_stage_report
    stage: '104'
    excerpt: '| A5 | sympy | 41 | `Y[z^6] + 11/729 == 0` | Y fingerprint | yes |'
  - &id019
    path: redteam/pass2/reports/stage_104.md
    line: 63
    role: pass2_stage_report
    stage: '104'
    excerpt: '| A6 | sympy | 42 | `Y[z^7]/I + 1/243 == 0` | Y fingerprint (odd) | yes |'
  - &id020
    path: redteam/pass2/reports/stage_104.md
    line: 68
    role: pass2_stage_report
    stage: '104'
    excerpt: '| B2 | mathematica | 43 | `Y[z^2] − 1/9 == 0` | Y fingerprint | yes |'
  - &id021
    path: redteam/pass2/reports/stage_104.md
    line: 69
    role: pass2_stage_report
    stage: '104'
    excerpt: '| B3 | mathematica | 44 | `Y[z^4] − 4/81 == 0` | Y fingerprint | yes |'
  - &id022
    path: redteam/pass2/reports/stage_104.md
    line: 70
    role: pass2_stage_report
    stage: '104'
    excerpt: '| B4 | mathematica | 45 | `Y[z^5]/I − 1/27 == 0` | Y fingerprint (odd) | yes |'
  - &id023
    path: redteam/pass2/reports/stage_104.md
    line: 71
    role: pass2_stage_report
    stage: '104'
    excerpt: '| B5 | mathematica | 46 | `Y[z^6] + 11/729 == 0` | Y fingerprint | yes |'
  - &id024
    path: redteam/pass2/reports/stage_104.md
    line: 72
    role: pass2_stage_report
    stage: '104'
    excerpt: '| B6 | mathematica | 47 | `Y[z^7]/I + 1/243 == 0` | Y fingerprint (odd) | yes |'
  - &id025
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
  - &id026
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
  - &id027
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
  - &id028
    path: redteam/pass2/reports/stage_105.md
    line: 47
    role: pass2_stage_report
    stage: '105'
    excerpt: '| D3 `\chi_Q=1` **from outgoing DtN fingerprint** | solve/Reduce retarded ω⁵ coeff `= a^5/(27c_s^5)` → `\chi_Q-1=0`
      (py:49–55; wl:52–59) | **partial** — pins `\chi_Q=1` but matches against a HARDCODED literal `a^5/(27c_s^5)`, NOT the
      Hankel DtN value; the DtN fingerprint identity (`z d/dz ln h_2^{(1)}`) is never exercised (card Check #3 not performed)
      |'
  - &id029
    path: redteam/pass2/reports/stage_105.md
    line: 69
    role: pass2_stage_report
    stage: '105'
    excerpt: A5/B3 are the load-bearing `\chi_Q=1` pin and are the rows that fail the "Anchored?" test.
  - &id030
    path: redteam/pass2/reports/stage_105.md
    line: 83
    role: pass2_stage_report
    stage: '105'
    excerpt: '105 owns the canonical pin `\chi_Q=1`. The paper grounds this pin in a genuine identity: the *outgoing target*
      coefficient `i a^5\omega^5/(27c_s^5)` is the normalized exact-Hankel DtN fingerprint `\widehat Y_2^{out}=-3/\Lambda_2^{out}`,
      `\Lambda_2^{out}(z)=z\,d/dz\ln h_2^{(1)}(z)` (appendix `eq:app-part04-Lambda-out-dtn`/`eq:app-part04-Yout-dtn`), and
      the card''s Check #3 explicitly is "Check the outgoing \(l=2\) DtN fingerprint against the normalized \(z=\omega a/c_s\)
      expansion."'
  - &id031
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
  - &id032
    path: redteam/pass2/reports/stage_105.md
    line: 98
    role: pass2_stage_report
    stage: '105'
    excerpt: 'Note the asymmetry: the deformed half (A6–A8 / B4) *does* derive its target coefficients from the operator (`-3/\Lambda`
      in SymPy; polynomial inversion of `\Lambda\cdot Y=-3` in Mathematica), so `\chi_Q=\xi_Q` is genuinely exercised. Only
      the canonical `\chi_Q=1` pin — the checkpoint''s headline result — is left self-referential.'
  - &id033
    path: redteam/pass2/reports/stage_105.md
    line: 102
    role: pass2_stage_report
    stage: '105'
    excerpt: 'This is a checkpoint, and 105 is the unit that establishes `\chi_Q=1` as the DtN l=2 fingerprint pin carried
      forward by stages 097/100/106 and the 2.5PN/4PN bridge (107–113). A reader auditing the chain would expect the pin to
      be forced by the Hankel DtN identity, per Check #3. As written, a transcription error or sign error in the *true* DtN
      coefficient could not be caught here: the script would still "pass" because both sides carry the same literal. The pin
      is correct, but it is asserted, not verified.'
  - &id034
    path: redteam/pass2/reports/stage_105.md
    line: 106
    role: pass2_stage_report
    stage: '105'
    excerpt: 'Make the canonical match exercise the actual outgoing DtN fingerprint rather than a hardcoded literal, in BOTH
      engines. Construct the exact `l=2` outgoing DtN operator from the spherical Hankel function and take its series, then
      normalize and read off the imaginary z⁵ coefficient as the *target*, instead of typing `a^5/(27c_s^5)`:'
  - &id035
    path: redteam/pass2/reports/stage_105.md
    line: 110
    role: pass2_stage_report
    stage: '105'
    excerpt: '- assert `Lam_out_series == -3 + z^2/3 + z^4/9 + i z^5/9` (this is the genuine fingerprint exercise — Check
      #3);'
  - &id036
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
  - &id037
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
  - &id038
    path: redteam/pass2/reports/stage_105.md
    line: 167
    role: pass2_stage_report
    stage: '105'
    excerpt: '**Checkpoint higher-bar verdict: not-cleared.** chi_Q=1 derivation result: the value is *correct* (independently
      re-derived), but the in-script assertion does not *exercise* the DtN-fingerprint identity — it self-matches a hardcoded
      copy of the target coefficient, so the canonical pin is asserted rather than forced. The higher bar (substantive, non-tautological,
      fingerprint actually exercised) is not met until F1 is fixed.'
  - &id039
    path: redteam/pass2/reports/stage_105.md
    line: 184
    role: pass2_stage_report
    stage: '105'
    excerpt: '| `\chi_Q = 1` (canonical pin) | py:53–55 / wl:58–59; sympy txt:19, math txt:17–18 | notes:69 (boxed); card
      line 16, appendix:326 `eq:app-part04-chiQ-equals-one` | MATCH (value); see F1 re anchoring |'
  - &id040
    path: notes/CHECKPOINT_CONSTANT_PROVENANCE.md
    line: 86
    role: checkpoint_constant_provenance
    stage: '104'
    excerpt: fingerprint coefficients from stage 104).
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
    citation: *id001
  - candidate_key: stage_104_matched_fingerprint_value
    modality: numeric_literal
    anchor_stage: '104'
    parameter_names:
    - matched_fingerprint_value
    reason: target-related numerical literal or closed-form coefficient
    citation: *id002
  - candidate_key: stage_104_matched_fingerprint_value
    modality: numeric_literal
    anchor_stage: '104'
    parameter_names:
    - matched_fingerprint_value
    reason: target-related numerical literal or closed-form coefficient
    citation: *id003
  - candidate_key: stage_104_matched_fingerprint_value
    modality: numeric_literal
    anchor_stage: '104'
    parameter_names:
    - matched_fingerprint_value
    reason: target-related numerical literal or closed-form coefficient
    citation: *id004
  - candidate_key: stage_104_matched_fingerprint_value
    modality: numeric_literal
    anchor_stage: '104'
    parameter_names:
    - matched_fingerprint_value
    reason: target-related numerical literal or closed-form coefficient
    citation: *id005
  - candidate_key: stage_104_matched_fingerprint_value
    modality: numeric_literal
    anchor_stage: '104'
    parameter_names:
    - matched_fingerprint_value
    reason: target-related numerical literal or closed-form coefficient
    citation: *id006
  - candidate_key: stage_105_chi_q
    modality: numeric_literal
    anchor_stage: '105'
    parameter_names:
    - chi_Q
    reason: target-related numerical literal or closed-form coefficient
    citation: *id007
  - candidate_key: stage_105_chi_q
    modality: numeric_literal
    anchor_stage: '105'
    parameter_names:
    - chi_Q
    - matched_fingerprint_value
    reason: target-related numerical literal or closed-form coefficient
    citation: *id008
  - candidate_key: stage_105_matched_fingerprint_value
    modality: numeric_literal
    anchor_stage: '105'
    parameter_names:
    - matched_fingerprint_value
    reason: target-related numerical literal or closed-form coefficient
    citation: *id009
  - candidate_key: stage_105_matched_fingerprint_value
    modality: numeric_literal
    anchor_stage: '105'
    parameter_names:
    - matched_fingerprint_value
    reason: target-related numerical literal or closed-form coefficient
    citation: *id010
  - candidate_key: stage_105_matched_fingerprint_value
    modality: numeric_literal
    anchor_stage: '105'
    parameter_names:
    - matched_fingerprint_value
    - xi_Q
    reason: target-related numerical literal or closed-form coefficient
    citation: *id011
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
      excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization}
        block.  The computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition
        \(\chi_Q=1\).}
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
      excerpt: The canonical outgoing coefficient is therefore not free. On the exact spherical outgoing `l=2` DtN branch,
        the leading odd quadrupole coefficient is fixed to
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
      excerpt: \stagefield{Derivation ledger}{The verification step belongs to the \emph{Retarded \(2.5\)PN factorization}
        block.  The computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition
        \(\chi_Q=1\).}
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
    citation: *id012
  - candidate_key: stage_104_gamma_5
    modality: existing_provenance
    anchor_stage: '104'
    parameter_names:
    - Gamma_5
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id013
  - candidate_key: stage_104_chi_q
    modality: existing_provenance
    anchor_stage: '104'
    parameter_names:
    - chi_Q
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id014
  - candidate_key: stage_104_matched_fingerprint_value
    modality: existing_provenance
    anchor_stage: '104'
    parameter_names:
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id015
  - candidate_key: stage_104_matched_fingerprint_value
    modality: existing_provenance
    anchor_stage: '104'
    parameter_names:
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id016
  - candidate_key: stage_104_matched_fingerprint_value
    modality: existing_provenance
    anchor_stage: '104'
    parameter_names:
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id017
  - candidate_key: stage_104_matched_fingerprint_value
    modality: existing_provenance
    anchor_stage: '104'
    parameter_names:
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id018
  - candidate_key: stage_104_matched_fingerprint_value
    modality: existing_provenance
    anchor_stage: '104'
    parameter_names:
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id019
  - candidate_key: stage_104_matched_fingerprint_value
    modality: existing_provenance
    anchor_stage: '104'
    parameter_names:
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id020
  - candidate_key: stage_104_matched_fingerprint_value
    modality: existing_provenance
    anchor_stage: '104'
    parameter_names:
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id021
  - candidate_key: stage_104_matched_fingerprint_value
    modality: existing_provenance
    anchor_stage: '104'
    parameter_names:
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id022
  - candidate_key: stage_104_matched_fingerprint_value
    modality: existing_provenance
    anchor_stage: '104'
    parameter_names:
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id023
  - candidate_key: stage_104_matched_fingerprint_value
    modality: existing_provenance
    anchor_stage: '104'
    parameter_names:
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id024
  - candidate_key: stage_104_matched_fingerprint_value
    modality: existing_provenance
    anchor_stage: '104'
    parameter_names:
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id025
  - candidate_key: stage_104_gamma_5
    modality: existing_provenance
    anchor_stage: '104'
    parameter_names:
    - Gamma_5
    - chi_Q
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id026
  - candidate_key: stage_105_chi_q
    modality: existing_provenance
    anchor_stage: '105'
    parameter_names:
    - chi_Q
    - matched_fingerprint_value
    - sigma_Q
    - xi_Q
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id027
  - candidate_key: stage_105_chi_q
    modality: existing_provenance
    anchor_stage: '105'
    parameter_names:
    - chi_Q
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id028
  - candidate_key: stage_105_chi_q
    modality: existing_provenance
    anchor_stage: '105'
    parameter_names:
    - chi_Q
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id029
  - candidate_key: stage_105_chi_q
    modality: existing_provenance
    anchor_stage: '105'
    parameter_names:
    - chi_Q
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id030
  - candidate_key: stage_105_chi_q
    modality: existing_provenance
    anchor_stage: '105'
    parameter_names:
    - chi_Q
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id031
  - candidate_key: stage_105_chi_q
    modality: existing_provenance
    anchor_stage: '105'
    parameter_names:
    - chi_Q
    - xi_Q
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id032
  - candidate_key: stage_105_chi_q
    modality: existing_provenance
    anchor_stage: '105'
    parameter_names:
    - chi_Q
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id033
  - candidate_key: stage_105_matched_fingerprint_value
    modality: existing_provenance
    anchor_stage: '105'
    parameter_names:
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id034
  - candidate_key: stage_105_matched_fingerprint_value
    modality: existing_provenance
    anchor_stage: '105'
    parameter_names:
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id035
  - candidate_key: stage_105_t_dtn
    modality: existing_provenance
    anchor_stage: '105'
    parameter_names:
    - T_dtn
    - chi_Q
    - matched_fingerprint_value
    - outgoing_dtn_fingerprint
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id036
  - candidate_key: stage_105_matched_fingerprint_value
    modality: existing_provenance
    anchor_stage: '105'
    parameter_names:
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id037
  - candidate_key: stage_105_chi_q
    modality: existing_provenance
    anchor_stage: '105'
    parameter_names:
    - chi_Q
    - matched_fingerprint_value
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id038
  - candidate_key: stage_105_chi_q
    modality: existing_provenance
    anchor_stage: '105'
    parameter_names:
    - chi_Q
    reason: pass-2 reconciliation or red-team provenance seed mentions a candidate value
    citation: *id039
  - candidate_key: stage_104_matched_fingerprint_value
    modality: existing_provenance
    anchor_stage: '104'
    parameter_names:
    - matched_fingerprint_value
    reason: checkpoint provenance seed mentions a candidate parameter
    citation: *id040
stage_results:
  '003':
    candidate_count: 0
    structurally_vacuous: true
    dry_run: true
  '104':
    candidate_count: 1
    structurally_vacuous: false
    dry_run: true
  '105':
    candidate_count: 1
    structurally_vacuous: false
    dry_run: true
```

Task:
Ask one question: which stage, value, parameter, or claim class did no modality cover? If you identify a missing class, emit a fifth-scan request. Do not re-score candidates already covered.

Emit only YAML:

```yaml
critic:
  uncovered_items: []
  fifth_scan_requests: []
  notes:
```
