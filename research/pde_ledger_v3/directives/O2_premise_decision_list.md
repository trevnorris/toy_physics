# O2: premise decision list (orchestrator)

**Author:** Claude (orchestrator), 2026-10-06. **Status:** folded once after one Codex + Grok pass (`CLAUDE.md` G2;
both NEEDS REVISION).
It implements sub-step 1 of the cleared inventory `directives/O2_steady_brane_balance_scoping.md` (`59855a38`),
and points at that inventory for every operand named here.

## Filing
O2 is a bounded sub-step of its own. It declares interfaces to S12, S14a/S14, S16 and S21 (inventory §4) and owns
none of them. This filing was the orchestrator's recommendation, and the user did not object (2026-10-06).

## User-selected premises (2026-10-06)
Each is an **adopted substrate input to a conditional model**, for S21's sort (inventory §6). None is a derived
result.

1. **Material reference: viscoelastic (the "glacier" choice).** The shear carrier responds elastically in the
   optical shear regime and relaxes under steady load. The user chose this on 2026-10-06 over a separately offered
   option: first resuming the paused retained-reference/relaxed-reference comparison (inventory §2). This
   supersedes the order of that earlier steer.
   - The choice revisits S9's no-dissipation and frequency-independent-moduli limits
     (`steps/S9_light_requires_shear.md:349–350`).
   - The reference and strain evolution law is the relaxation response. It is a general unknown, not an
     engine-chosen family. The input contract (sub-step 2) supplies its form, or keeps it as a named operand.
   - The input contract carries the steady state's energy balance explicitly, including the relaxation
     response's power. Its sign is left to the computation, and any net power names its supplier.
   - The relaxation response's consequences for light in the optical regime are a later light-compatibility
     question. They are not an O2 deliverable.
2. **Drive: the drain itself.** No separate external body force. The drive is the committed order-conversion
   drain with its source and boundary/return data, acting through O2's stress, traction and O3 operands. Either
   `F_drive` is identified with those terms or it is absent. The far-field link to `GM` stays with the S14a
   bridge and S16, and is open in O2.
3. **Exchange momentum: the brane's local velocity.** Converted material carries the brane material's local
   velocity `V`. This closes the transported momentum only.
   - Any additional non-variational momentum partners, and the system that carries their reaction, stay OPEN
     with S12.
   - The drive, the face/support tractions `T_hold,s` and the O3 transport are declared as separate terms.
4. **Bulk: normal loading only.** Keep the postulated shear-free scalar bulk (`steps/S9_light_requires_shear.md:185`;
   `directives/S11b_SHARED_PHYSICS.md:164`). Bulk loading on the brane is face-normal, and momentum also moves
   through O3. No tangential bulk stress.

## Retained OPEN operands
Each stays named, as a general unknown:
- the material branch (A13);
- the live brane stress, inertia and normal material response. Premise 1 fixes its character, not its form.
- the stiffness and density responses `ℳ_⊥` and `ℛ_br`;
- the live embedding relation `ℰ_h^live` (O4). L3's field identity `ξ_w = ℓh` is retained, and O4's identity with
  O2's normal part stays unsettled.
- `T_hold,s`. Its bulk part is face-normal under premise 4. A declared support would be a held input.
- the core holder and mouth data (O5, owned by Q2/S22);
- the sheet/slab map `𝒥_map`;
- every missing grade. v9's counting is kept, and the deliverable is the untruncated named balance.

## Intended claim
A conditional O2 balance law: the steady in-plane and normal momentum and support relation with `V` and `j_n`
live, conditional on premises 1–4, with every OPEN operand named. It is not a derived brane law. Its record
carries register entries for S21.

## Next
Sub-step 2 of the inventory: an input contract, authored by Codex and reviewed until clear by a fresh Claude
agent and Grok.
