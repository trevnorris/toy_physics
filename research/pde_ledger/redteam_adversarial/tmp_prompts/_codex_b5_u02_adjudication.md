# Codex adjudication request — batch-5 u02 (γ₀ overclaim?)

You are the INDEPENDENT co-adjudicator in a layer-2 fit-vs-derive audit of a toy-physics derivation
ledger. Read-only. Do NOT edit any file. Do NOT run scripts. Read and reason, then return a verdict.

## The question (narrow, layer-2)
The parameter γ₀ = (1+𝔯²)/9 is classified `free_choice` in its provenance bundle. An adversarial auditor
returned a YES verdict — a "posit-dressed-as-derived" FATAL overclaim — at **stage 162**, on the grounds
that the notes call it an "exact parent-family formula" and say "Stages 115–116 **fix** the bare odd
normalization by γ₀=(1+𝔯²)/9", under a `\StatusExactClosure{}` card badge with no "ansatz" tag.

A defense auditor REFUTED the YES as fatal, arguing it is at most a non-fatal disclosure/scoping nit.

**Your job: independently decide FATAL vs NON-FATAL.** "Fatal" = the program ASSERTS γ₀ is a derived/forced
consequence in a load-bearing way that would be false. "Non-fatal" = a presentation/scoping imperfection
(e.g. an "exact" badge that should be scoped, or a missing ansatz hedge) that does not make the program
claim something false. We are NOT asking you to decide whether γ₀ is "really" derived or posited — that is a
separate, already-user-owned conceptual call (see ANSATZ_LEDGER §5 item 1).

## Files to read (all under /var/projects/toy_physics/research/pde_ledger/)
- paper/stages/stage_162.tex  (esp. :5 status badge, :16 quote, :27 caveat)
- notes/stages/moving_throat_pde_stage162_parent_compensation_rigidity.md  (§1 lines ~38-64, §2 lines ~68-104)
- paper/stages/stage_115.tex , paper/stages/stage_116.tex
- notes/stages/moving_throat_pde_stage116_dn_mixed_tube_realization.md  (the "A simple concrete realization is to take…" line ~62)
- notes/CHECKPOINT_CONSTANT_PROVENANCE.md  (around line 1210, the γ₀ "POSTULATED pure-scale ANSATZ, NOT derived" entry)
- redteam_adversarial/ANSATZ_LEDGER.md  (§5 item 1 — the γ₀ borderline entry)

## Decisive sub-questions
1. Does the word "exact" in stage 162 attach to the algebraic FAMILY LAW / relation (γ₀ written as a
   function of the single parent variable 𝔯), or does it assert the γ₀ VALUE is derived? Compare to how the
   genuinely-derived L_W/a relation is presented in the same section.
2. Does stage 162's actual RESULT (the rigidity identity / δ⊥=0 on the family) only USE γ₀'s functional
   form, so the conclusion is an exact consequence OF the posited family — not a claim that γ₀ is forced?
3. Is "Stages 115–116 fix …" a back-reference to where γ₀ was set upstream, and do 115/116 POSIT γ₀ (origin
   116:62 "A simple concrete realization is to take…") rather than derive it from an external target?
4. Is γ₀=(1+r_c)/9 already logged in ANSATZ_LEDGER §5 as a borderline (derived-vs-postulated) item awaiting
   a human gate — i.e. is the derived-vs-posited status an already-known open call, not a newly-discovered
   contradiction?

## REQUIRED OUTPUT (concise)
1. VERDICT: FATAL or NON-FATAL (and if NON-FATAL, the minimal fix class).
2. One line per sub-question (1-4): your finding + the file:line that grounds it.
3. Whether you AGREE or DISAGREE with the defense's "non-fatal, contingent on the user's existing
   ANSATZ_LEDGER §5 borderline call" disposition, and why.
