# Preflight build evidence map

Start with `build.md`, `worker.py` and `manifest.json`. Every runtime JSON file is
supplied in full under `saved/` at its manifest alias. The private packet is a
fixed source snapshot. Host paths in the manifest are provenance, not permission
to open the live repository. `packet-index.json` maps every private file to the
actual source path/hash and distinguishes reviewed projections from originals.
No old reviewer reports or live credentials are included.

## Execution paths

- `worker.py`: the new adapter/tail preflight. `run_science` contains all new
  mathematics; it is called only after gate and actual containment checks.
- `launcher.py`: the future hook-first guarded launch. No READY gate exists yet.
- `execution-authority.json`: standing authority for one bounded preflight.
- `runtime-source/raw-helper.py`: only the named inert helper definitions are
  extracted, including containment, JSON decode and durable Journal methods.
  Its old scientific body is NOT executed.
- `runtime-source/shared-guard.py`, `supervisor.py`, `completion-hook.py`: unchanged
  actual execution route. `resources` in the manifest joins 4GiB native/cgroup,
  16GiB pool, zero swap, one CPU/thread,32 tasks and no deadline.
- `tooling-tests.py`: stdlib metadata/source/synthetic tests only. No science is
  restored or computed by these tests. `tooling-test-record.json` is a receipt,
  not a source-review verdict or runtime scientific result.

## Restore versus new adapter obligations

- `saved/inventory/THETA_BALANCE-ordered-addresses.json`: original 2,652-entry row.
  `saved/selected/pressure-addresses.json` is the COMPLETE 544-entry projection.
  `saved/local/all-local-cells.json` has the original400; selected/local-cells
  has all16 selected cells, their summands and complete inherited proof records.
- `saved/pressure/fields.json` and `coefficient-certificates.json`:34 field IDs,
  original expressions and saved quotient polynomial certificates.
  `saved/field/<ID>-polynomial.json` stores numerator polynomial/coefficients and
  a constant denominator; the reconstruction input and literal zero return are
  alongside it. Local cell `polynomial.coefficients` already stores quotients.
  A useful nonconstant field is
  f588c94387ca5afe8f5285403a543afc513deab9465dd1108aebab14058e7126.
- `saved/factors/address-full-factor-<N>-operands.json`, full-mapped-residual input
  and zero return, normal-source-join input and zero return: all17 used proof
  sets. Proof reuse joins actual source/normal/response maps and complete factors.
- `saved/weak/weak-address-coverage.json`: full13,260-entry completed route record;
  the worker joins selected addresses by actual ID/face/slot/grade/jet/field.
  It does not regenerate the original inventory.
- `saved/local/context.json`, `physical-input.json`, inventory/native-profile-scale
  and prebinding-native-speed-inventory: actual held omega3 and original physical
  source, native L10 and source-speed independence. `ends/left-match.json` holds
  the saved cs/kappa point; left-source-binding holds optional finite origin.
- `saved/inventory/fourier-and-unit-provenance.json`: native source normalization,
  edge reduction and units. These are inherited contracts; the preflight does
  not claim newly reconstructed dimensions for every numerical action summand.

## Whole kernels and new bounds

- `saved/pressure/whole-tags.json` includes actual full definitions and byte
  receipts. Their originals are `saved/saved/reference/left-height-subtracted-PV.json`,
  `right-height-PV-operands.json` and `saved/saved/direct/closed-density.json`.
  The double `saved/saved` is just the manifest alias, not a second integration.
- `saved/pressure/typed-direct.json`, `whole-definitions.json` distinguish bare
  Rprod, factored external denominator and complete closed D density. The worker
  uses the density, not an isolated raw factor. No whole value is evaluated.
- `saved/pressure/global-parameter-domain.json`, `profile-envelope.json`,
  `shift-root-bound.json`, `whole-envelopes.json`, `H-bound.json`, the two
  height-PV files, normal-growth and J/D numerator files: full inherited global
  certificate operands. The new rational tail bounds in `build.md` are derived
  from these constants and the exact numerical adapters; they are NOT old
  existence constants relabelled as truncation errors.
- `method.md`: the unchanged wavelength-resolved packet-action method.
  `background/pressure-weak-method.md` details the source-bound global envelopes.
  `background/full-weak-method.md` gives the completed local complement.

The new worker stops after exact adapters, conditional analytic pressure-tail
constants/K/U/T and metadata eligibility of the four controls. It does not carry
out the method's independent Fourier routes, collision quadrature, local action
or numerical controls. Their absence is explicit, and neither source clearance
nor preflight completion is called numerical-action acceptance.
