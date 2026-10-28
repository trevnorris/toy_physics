# S15 reading decision list, review pass 2: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S15_reading_decision_list_review_p2_lookups.md`, beside this
file (generator `_scratch/s15/gen/s15_reading_p2_lookups.sh`). "M-n" cites a section of that file. Every computed
result quoted here is a leg's saved stdout; the legs' scripts were not re-run.

## The pass

- **Artifact.** `directives/S15_reading_decision_list.md`, revision 2, committed at `48d06300`. It folds pass 1
  (`S15_reading_decision_list_review_p1_disposition.md`, F1–F4) into revision 1, which is preserved at `1c6695df`.
  The list, its setup lookups and the pass-1 disposition and lookups matched the baseline after the pass (M0).
- **Gate.** The list's spec content is reviewed until clear. Revision 2 changed it, so this is the second Codex + Grok
  pass.
- **Prompt and legs.** `_scratch/s15/reading_dl/reading_dl_review_prompt_p2.md`, sha256 `db29c2f7…` (M0). Both legs
  received it identically, in parallel:
  - Codex (`gpt-6.1-sol`, effort xhigh): `/var/projects/_review_scratch/s15-readdl-p2-codex/REVIEW.md`;
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s15-readdl-p2-grok/REVIEW.md`.
- **Verdicts.** Codex FINDINGS, one finding, MEDIUM. Grok CLEAR (M1).
- **The tree.** The status matched the pre-pass snapshot, and no file in the repository root is newer than the
  prompt (M0).

## Findings

**N1 (Codex 1): the fold's named pressure conditions are stated as conditions without which the derivative does
not exist. They are not that, and they are not a complete set of conditions under which it does. Confirmed.**
- *The list* (`:50–53`, M2): "For instance, the separation derivative of the whole pressure, and of the whole minus
  the expansion, needs `p_br` differentiable at the total density. The response part of the expansion needs `p_br′`
  differentiable at `ρ_2 + ϱ`." Item 2 asks for the parts "exact wherever they exist on the domain". Line 84 reads
  "Every object is exact on S15's domain wherever it exists", and the intended claim reads "wherever it exists".
- *Not necessary (Codex, M2):*
  - **The expansion.** Take `p(ρ) = ρ + (ρ − 1)|ρ − 1|/2` and the base `ρ_2 + ϱ = 1`. Then `p′(1) = 1` from both
    sides, and the one-sided derivatives of `p′` at the base are `(−1, 1)`. The base's change is zero, and drain 1's
    amplitude change `B′` is live. Codex's stdout prints `expansion_Eresp: 8*pi*Bprime*R**3*d/3`, which exists and
    is nonzero.
  - **The whole pressure.** Take `p(ρ) = |ρ − ρ_0|` with the total density `ρ_0 + (z − z_0 − εd)²`. That density
    reaches the kink on the surface. The composite pressure is `(d*eps - z + z0)**2`, and the separation derivative
    is `8*pi*R**3*d/3`.
- *Not sufficient (Grok's stdout, M3; Grok did not file it).*
  - Take `p = (ρ − 1)|ρ − 1|/2`, base `1` and amplitude `1`. The total density is `2`, where `p` is differentiable.
  - With the base moving, the whole minus the expansion has one-sided derivatives `2` and `0`.
  - Grok's report: "A whole-minus response whose base moves has to carry differentiability of `p′` at `ρ_2 + ϱ` as
    well as differentiability of `p` at the total density". Grok judged that the list's general sentence covers this,
    since "The two pressure sentences are marked "For instance"" (M3).
- *Effect.*
  - A build following the named sentences can print a condition under the label of one the object needs. The record
    would then read that the object does not exist where it does: for instance, for a pressure law with a kink.
  - Read as the conditions under which the objects exist, the named sentences omit one: the whole minus the
    expansion with its base moving.
  - "Wherever it exists" promises the objects on a larger set than any printed condition establishes.
- *What must be true:*
  - nothing in the list presents a condition as one without which a separation derivative does not exist, unless it
    is shown to be one;
  - nothing presents a set of conditions as complete for an object unless it is;
  - each departure object carries a condition under which its separation derivative exists, labelled where it is
    printed as a condition under which it exists;
  - the claim and the record claim each departure object where its printed condition holds, and nothing about it
    where that condition does not hold.

## Examined and not filed

- **The fold.** Grok checked each "What must be true" of F1–F4 against revision 2, with the fold diff saved (M4).
  Codex: "F1/F2 are folded; F3/F4 now have a universal existence clause".
- **Codex's other probes** (`physics_checks_corrected/stdout`, hash in M5).
  - The operator: zero departure for a scalar inverse square and a radial vector; nonzero for a turning vector and a
    one-point zero; `2*Q + Qprime` for an object without `X_1` whose change is live.
  - The split: `0` for the same named fields inside or outside the geometric part.
  - The two F1 cases again: `(0, 2, 2)` and `(3, −3, 0)`.
  - The instant jet; a moving shell's transport residuals `0`; held radii unchanged along the family; both velocity
    parts; held motion still in the rate; dimension residuals `(0, 0, 0)`.
  - Changing the labelled split moves `−A*d` between the parts with the total unchanged. That is pass 1's examined
    item, and revision 2 still conditions the claim on the labelled choices.
- **Grok's other probes** (`01`–`08`, hashes in M5, empty stderr).
  - The operator: homogeneity of degree −2 on every ray, per Cartesian component.
  - The parts as objects: one pair for the unsplit integral and the product.
  - The held quantities kept in the objects. Holding `R` keeps a finite-size term that scaling it with `d` would
    remove.
  - The instant jet differs from a time path, and the list chooses the instant.
  - The profile family's closed forms, with special exponents `α = −2`, `−3`, `−4`, and two writings that agree on
    the domain. A truncation's departure is zero where the exact one is not.
  - The velocity parts linear in each velocity; the domain with room for a two-sided step; S15's tags found in its
    transcript, with the three named sums built from them.
  - A carried profile not differentiable in its radial argument. Grok reads it as covered by the general sentence, as
    is the case in M3. N1's "What must be true" keeps that general obligation and removes the named instances' claim
    to be more than examples.

## The repair

N1 changes what may be claimed: a printed label can read as non-existence where the object exists, and the claim's
coverage exceeds what the printed conditions establish. Revision 2 is preserved as committed at `48d06300`, beside this
disposition. Revision 3 makes N1 true. It is the third pass on this spec, so how to run it goes to the user first.
