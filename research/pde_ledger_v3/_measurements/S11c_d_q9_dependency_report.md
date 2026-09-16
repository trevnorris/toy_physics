# Q9 orientation repair: S11c-d dependency disposition

The current S11c-d momentum-domain job can continue. The reported Q9-only
correction has no identified value dependency into its closed slab/kernel,
energy-current input or bound numerical action. No S11c-d numerical rerun is
indicated by this trace. This does not validate the repair, Q9 in D3–D5, or
other unrelated properties of the downstream calculations.

## Evidence

- The working producer differs from HEAD only in `compute_q9`; its patch is
  owned by the Lean session and remains untouched. The producer and export
  baseline hashes match the reported D2 counterexample. Wolfram's source
  already uses `Transpose[actionRows]` before its connected nullspace.
- `emit_q9` and `PD_TERM` produce a 72-row family across D2–D5, including
  unchanged ordering/reflection bookkeeping. All 72 values are copied exactly
  from S11 through S11b, S11c-a and S11c-b accumulated exports. This identifies
  affected containers, not 72 independently disproved values.
- `package_build` uses the computed `PD_DENSITY_PLACEHOLDER` in the
  `XFORM_EXTRA` action only. MAIN uses the separately constructed curl/divergence
  densities. `q7_objects` accepts a Q9 argument but never reads it. The Q9
  bases, their reflected/Euler–Lagrange descendants and parity-odd construction
  require revalidation; whether `PD_TERM` and extra-action outputs change must
  be measured, not presumed from the basis defect.
- S11b reads S11's `c_s0`, `mu_R`, `rho_br` and their dimension rows for its
  physical construction. Its energy basis is independently constructed from
  gradient contractions (`construct_energy_basis`, line 601). Whole-ledger
  iteration is export carry/roundtrip work. S11c-a binds named S11b identities;
  S11c-b constructs its uniform and gradient-energy basis from contractions and
  Euler signatures (`enumerate_uniform_candidates`, line 1477), with named
  shape/balance inputs. Their inherited-access inventory is saved in the trace.
- S11c-c1/c2/d manifests have 44/35/73 keys and no Q9-family keys. The c1/c2
  delta exports contain none of the 72 Q9 rows. Executing only the actual pure
  S11c-d `bind` function on opaque serialized rows gives exactly the same 73
  inputs after deleting all 72 Q9 rows. Mutating an actual consumed closed-row
  value is detected. This is a selector dependency control, not a new physics
  equality check.
- The accepted reduction transcript's decoded import lookup list is exactly
  the current 73-key manifest; its recursive closure has no Q9-family key.
  The reduced-action, action-probe and bound-momentum packet hashes match their
  accepted checkpoints. A nonexecuting pickle-string census finds no Q9 row
  names in those packets, consistent with the source-level data path.
- All 78 current/frozen source hashes of the active momentum-domain run match.
  It does not pin or execute `S11_stray_longitudinal_sympy_audit.py`. Its
  ancestor export inputs are the frozen S11c-b/c1/c2 modules. The host process
  census showed this supervisor/coordinator, its watcher and four wider-box
  numerical workers; no concurrent S11 Q9 or Wolfram CAS job was found.

The audit used ASTs, hashes, the existing lossless text decoder and pickle
opcodes. It imported no CAS engine and repeated no numerical integration.
The initial instrument compared the encoded lookup wrapper against the manifest;
that guard correctly rejected the extra codec marker. Decoding with the unchanged
codec resolves the comparison. Both audit versions and the instrument note are
preserved in `_scratch/s11c/s11c-q9-impact-20260916/`.

## What to rebuild or refresh after the upstream repair is accepted

| Object | Required disposition |
| --- | --- |
| S11 Q9 results, D2–D5 | Recompute and validate the invariant spans and descendants, including exceptional dimensions; counts alone are insufficient. Independent fidelity reviews remain with the repair session. |
| `PD_TERM` / `XFORM_EXTRA`, D2–D5 | Compare corrected parity-odd densities first. Recompute dependent extra-action results where the density changes; retain explicit equality evidence where it does not. |
| S11 MAIN non-Q9 dynamics | Source trace shows independence from the Q9 value. Preserve them with a before/after dependency/value check rather than infer a need for a full dynamics rebuild. |
| S11/S11b/S11c-a/S11c-b accumulated exports | Refresh corrected carried Q9 rows and provenance through a controlled, reproducible process. Preserve original files and compare all other row payloads. Full original run transcripts remain historical; do not silently relabel them as repaired. |
| S11c-c1/c2 delta exports and source manifests | Compare the actually consumed roots and recursive symbol/dimension closure. Their numerical payloads have no identified Q9 dependence, but ancestor whole-file hash changes need explicit provenance joins. |
| S11c-d normalization, channels, reduction/factorization and numerical caches | Retain the accepted operands. Before rebasing to repaired exports, prove consumed-row/closure identity and record an explicit old/new source join. No expensive rerun is indicated for a Q9-only carried-row update. |
| Running momentum-domain production | Continue against its immutable current sources; perform normal completion validation and preserve all partials/results. |

Do not overwrite the pinned S11c-b/c1/c2 exports while production runs. Changing
an unused carried row still changes a whole-file SHA and will make the strict
source guard reject publication. After completion, use explicit source-rebase
records; never patch hashes alone or rewrite old evidence in place. If the
upstream repair changes any consumed physical root, stop and reassess the
minimal dependent recomputation before continuing.

This trace does not authorize adopting the pending upstream patch. No source,
export, transcript, cache or Lean file was changed by the dependency audit.
Only this report, its instrument/checkpoint and the S11c-d tracking records are
new work.
