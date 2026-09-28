# Two-asymptote implementation review disposition

Both approved independent reports finished before inspection, adjudication or
edits. Their literal verdicts are **Claude: NEEDS REVISION** and **Grok: CLEAR
FOR THIS BOUNDED INGREDIENT**. The local corrections below do not constitute
a fresh independent CLEAR. No reviewer rerun, peer sharing or science ran.

The exact 28-file / 611,908-byte packet remains in
`_scratch/s11c/s11c-d-two-asymptote-20260927/build-review/packet`, with SHA-256
`0bfabfac59801cd96b780b1d878095dbdcbad9fdd15237803426e21631d44e9c`.
The archive, prompt, approval and all source hashes were verified before edits;
the original worker/manifest/implementation description were also copied to
the sibling `reviewed-source` directory. Complete literal reports are retained
in the original Claude JSON `result` and Grok JSON `text` fields and in the
[Claude](S11c_d_two_asymptote_claude_review.md) and
[Grok](S11c_d_two_asymptote_grok_review.md) text copies.

Claude completed in 612.752 seconds with empty stderr, a completed response
and no permission denials. Grok completed in 790.522 seconds, after Claude,
with `end_turn` and 1,431 bytes of startup configuration warnings concerning
plugin precedence and unsupported hook settings. These warnings are retained;
they are not a scientific stderr stream. The coordinator and watcher stderr
were empty. The hook finished for session
`01a0e01b-ef84-7192-817f-584cda5d339b`. Literal content, rather than exit codes,
supplies the verdicts. The historical state field saying submission approval
was pending is superseded by `explicit-packet-approval.json`, which records
the user's “You may submit” and exact packet/archive identities.

## Substantive findings and corrections

**S1 — formal integral coefficient collection: accepted and corrected.**
The reviewed reducer split Add terms inside an Integral but did not pull out
coefficients independent of its bound variables. Equivalent combined and
separate end terms could therefore leave different Integral carriers. The
corrected `linear_form` uses `as_independent` on all bound variables, preserves
the original ordered limits, and performs no integral or limit evaluation.
`formal_coefficients` encodes intact maximal Integral carriers with temporary
symbols and records every row/column, carrier power and scalar coefficient.
The baseline guard reduces those actual scalar coefficients to exact zero.
No evaluated-integral identity, integration interchange or convergence claim
is inferred from this formal polynomial algebra.

**S2 — mutation of a real decomposition term: accepted and corrected.**
The reviewed `side - chi*plain` operand was not the actual addend used in the
decomposition. The corrected `decompose` builds each side from its local part
and full addressed native term list. The mutation omits exactly one addressed
nonlocal side term and reassembles the decomposition, holding the independently
computed direct forcing fixed. The same carrier reducer processes the mutated
residual. A surviving carrier must have an exactly nonzero reduced coefficient
before the formal control is marked responsive. Complete inputs, omitted term,
changed decomposition, raw residual and coefficient reductions are saved.
The job stops if no such certificate is established. This closes a **formal
instrument control** only: `integratedResponseNonzeroEstablished` remains false.
The algebraic end-lift mutation remains a separate control.

**S3 — persistence before reconstruction guard: accepted and corrected.**
Both raw native plane actions and the symbol-plane comparison operands are
now saved before constructing or guarding the formal reconstruction residual.
A later failed carrier check cannot erase this method evidence.

**S4 — nonlocal plane-jet identification: accepted and corrected.**
`UNVERIFIED_NONLOCAL_PLANE_JET_EXTENSION` replaces the established-sounding
flag. The worker records whether the actual zero-grade strong action contains
the source Abel regulator or profile functions. These absence checks are
necessary diagnostics; they do not prove translation invariance or that
momentum differentiation passes through the ordered integrals and regulator
limit. Both obligations remain explicit. The saved decomposition is labelled
conditional on that identification; its algebraic reconstruction can be
verified without claiming the identification is an evaluated operator result.

**S5 — grade-free integrals: confirmed by existing inventory, guarded join added.**
The source-generated `record-inventory.json` has 375 records: 100 local, 160
cell, 80 factor and 35 source. Every factor/source record lists only grade
`(0,0,0)`. The 115 individual saved record bytes match the inventory hashes,
and the combined packet matches the accepted checkpoint hash. See the
[metadata evidence](S11c_d_two_asymptote_grade_inventory_evidence.json).
No scientific object was restored to obtain this evidence; the inventory is
not represented as a newly validated mathematical result. The corrected
worker must join those exact addresses to the restored combined packet,
require its actual `COEFFICIENTS` keys and original factors/sources to be
grade-free, and still require the intact original nonlocal integrals to be
grade-free. Thus no missing internal grade convolution is silently omitted.
An unexpected mismatch stops with actual operands; no derivative fallback
or producer replay is added.

Both reviewers support the degree-three local quotient, first-grade principal
parts, physical branch, factorial convention and exact end-symbol equations.
Grok did not identify S1–S4 as blockers, but its agreement does not override
Claude's concrete instrumentation findings. The local fixes directly implement
the requested bounded repairs. Both reviewers accept a useful ingredient stop
that leaves the forcing domain unresolved; Claude conditions that acceptance
on the substantive repairs above.

## Other bounded dispositions

The strict unit checker now includes the source's exact dimensionless
`[wm]1_profile` symbol rule. Origin selection uses `Subs(...).doit(deep=False)`
on integral/limit-free operands to align the native substitution convention;
unsupported nonlocal cases stop instead of evaluating an integral. The
accepted checkpoint chain, operation inventory and domain hash are joined
explicitly. Success/failure inventories now hash every saved artifact,
including non-journal bundles.

The 200 optional local summand limit evaluations are removed from this first
bounded run, as Claude suggested. Every local weighted limit operand and
full-defect limit operand is still saved and marked **unevaluated**. A local
summand limit would not establish the full nonlocal forcing domain; running
all of them before later core records risked losing the useful deliverable to
the native deadline. No decay result or weighted bound is now claimed for
them. This supersedes the local-limit-evaluation paragraph in the preserved
reviewed implementation description.

Optional right-sided principal-part checks, algebraic optimization, codec
tightening and further gauge tests are deferred. The exact local construction
already checks both inverse constant identities and both residue identities.
The end equation does not uniquely test homogeneous kernel-valued pieces of
the field correction; their selection remains the saved meromorphic product
formula. No new comparison/optimization campaign is introduced to address
optional observations. The existing restrictive hash-pinned saved codec is
unchanged. No optional wording review cycle is requested.

## Actual readiness and authorization

The corrected worker is checked by source/AST, byte comparisons, manifest
joins and standard-library persistence checks only. Its mathematical body and
new carrier reducer have **not** been executed. Runtime checks must establish
their actual source-bound outputs under containment. The local disposition
closes the implementation findings; it does not accept any scientific result.

The proposed single stage still stops at
`END_LIFT_SAVED_FORCING_DOMAIN_UNRESOLVED`. Remaining dependencies include the
nonlocal plane-jet extension, combined weighted tails and prescription pairing,
and any integrated nonlocal responsiveness claim. No full Green, response,
FORM, A11/A12, coincident-point extension, retarded equivalence or radiation
witness is accepted. Physical inputs and the evanescent slice are unchanged.

**Science approval is pending.** The explicit authorization received was to
submit the fixed review packet. The completion event expressly says: “No
science launch is authorized by this review-completion event.” A pending
gate and exact launcher can be prepared, but only a subsequent approval can
activate the one 900s outer/840s native stage, under the shared guard and
normalization supervisor with 2GiB/zero swap/one CPU/nice15/32 tasks/one native
thread. Preserve global locking, all prior results/failures and native limits.
No unlimited-duration inheritance, overlap, retry, replay or new review is
authorized. The local hook must be armed first for the existing session.
