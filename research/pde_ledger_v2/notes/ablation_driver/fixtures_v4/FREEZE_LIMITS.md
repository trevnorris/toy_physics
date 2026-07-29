# Freeze scope and trust boundary

This freeze governs bytes, not the execution environment. It binds the committed
fixture suite, its normative contract and requirements, the legacy comparison
inputs, and the fixed stage023 programs and data used by the retrofit. It requires
the authority's exact committed `HEAD` bytes, rejects suite-local Python bytecode,
and rejects worktree drift in every governed path.

Import-time substitution supplied from outside the suite can defeat the verifier
or conformance process while leaving every frozen byte intact. Examples include
startup hooks supplied through `PYTHONPATH`, user-site configuration, `.pth`
processing, or substituted executables found through `PATH`. A passing in-process
verification is therefore not evidence that the intended interpreter, standard
library behavior, or executable resolution actually ran.

`verify_freeze.git()` shells out to the name `git` with the verifier's inherited
environment. Executable resolution and Git behavior can consequently be
substituted without changing frozen bytes. Git also reads repository-local
`.git/config`; that file is not in the freeze inventory. The conformance
caller-tree guard includes `.git` when detecting changes made during an item, but
that later snapshot does not authenticate the configuration used by the freeze
verifier itself.

Closing this execution-environment boundary requires an out-of-process control
that is governed independently of this repository and of the process it checks.
That control must launch or attest the interpreter and executable resolution from
a controlled environment, suppress or independently govern startup/import hooks
and Git configuration, and compare the governed bytes against an independently
owned authority. Repository-local in-process code cannot supply that assurance
about itself.

Separately, the worktree verifier is its own decision function. An isolated edit
to that verifier can approve itself, and an edit to the conformance bootstrap can
skip it. A coordinated change can weaken the suite or verifier, regenerate the
authority, and commit both. The authority is a current-`HEAD` manifest, not an
independently controlled historical or signed trust anchor, and it does not
govern adjacent non-normative notes or the driver under test.

Closing that byte-authority gap requires an authority outside the fixture
author's and driver builder's write control, such as a protected external
commit/digest or a signature checked by independently owned policy, plus review
enforcement that those parties cannot update.
