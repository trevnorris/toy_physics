# Required instructions for ledger Lean work

Read [FORMALIZATION_POLICY.md](FORMALIZATION_POLICY.md) before planning or
changing a formalization in this directory. Read it again when resuming after
a handoff or context reset, and before proposing any expansion of proof scope.
It also governs work on the associated Lean generators, checks and reports
outside this directory. Do not rely on an older checkpoint's next-step list.

- In the first substantive progress update for Lean work, identify the policy
  and the specific theorem, fidelity or coverage-contract obligation being
  addressed. A document-only policy task needs document checks, not a proof build.
- The deliverable is the proof, a compact action/operator identification, an
  exhaustive coverage contract, and meaningful mutation controls, with statement
  fidelity reviewed. Match the effort to the claim and reuse existing proofs.
- Do not extend a systematic per-output CAS bridge. Repeated fields, solver
  metadata, every matrix entry/minor, and every engine/dimension transcript do
  not each become Lean obligations merely because they exist.
- Keep the committed S10 bridge as evidence. Its former expansion plan is
  superseded. The current S10 task is to consolidate and review the compact
  contract described in the policy and [s10/COVERAGE.md](s10/COVERAGE.md).
- Stop when the declared Lean completion criteria are met. CAS production and
  comparator/export work retain their own requirements and status; do not
  silently turn them into more Lean work or claim they are complete.

Follow the user's current instructions. A request to continue advances the
agreed scope; it does not by itself enlarge that scope. Apply this policy without
repeatedly asking for approval for work already authorized within the contract.
