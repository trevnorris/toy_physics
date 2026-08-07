# `_asbuilt/` — frozen predecessors that the registry still cites

⛔ **Nothing in this directory is live code.** ⛔ Do not import it, do not run it, do not adapt it, and
⛔ do not read it when building a replacement engine.

## Why it exists

`reduction/quantities.yaml` and `reduction/relations.yaml` carry `source_locus` / `execution_locus` entries
that are **file paths with line ranges**. The registry reader validates them on every load
(`registry_read.py::_validate_loci`): the file must exist and the range must lie inside it.
`source_locus` is **schema-required** — it cannot be blanked.

⇒ When a step's engine is rebuilt in place, a locus keyed into it **can break the registry** (if the new
file is shorter than the cited range) or **can come to point at unrelated code** (if it is long enough to
validate but the cited lines now hold something else). ⚠ The second is the dangerous one: the registry keeps
loading and the provenance is quietly false.

⛔ **Neither outcome is inevitable** — a replacement that happened to keep the cited content at the same
line numbers would stay correct. ⚠ But **nothing checks that**, and `_validate_loci` cannot: it verifies a
range *exists*, ⛔ never that the lines still say what the row claims.

⭐ So the artifact that actually produced a declared dimension is frozen here, and the locus points at the
frozen copy. ⚠ The declared dimension's provenance then stays **true** across the rebuild.

## What this is NOT

⛔ This is **not** a claim that the frozen artifact is correct, and ⛔ **not** a second opinion. It is the
predecessor, kept addressable. ⭐ A step record comparing a rebuilt engine against a registry row whose
locus points in here must say plainly that the declared dimension **came from the artifact being
replaced**, so a zero residual means *the new engine reproduces its predecessor* — ⛔ not that two
independent routes agreed. Shared-physics §Q6r exists to make exactly that visible, which is why it emits
the locus.

## Contents

| file | frozen from | cited by |
|---|---|---|
| `S11_stray_longitudinal_sympy_audit.py` | `scripts/S11_stray_longitudinal_sympy_audit.py`, byte-identical, at tag `s11-as-built` | `Q.brane.B_comp`, `Q.brane.c_L` (dimension provenance and `source_loci`), relation `R5` |

⚠ **A correction was made while repointing, and it is a real defect that predates the rebuild.** The
dimension provenance for `B_comp` and `c_L` pointed at lines `632-664` — which spans the tail of the `A2`
cross-derivative check and the start of the `A3` degeneracy block (`A3` begins at `638`), and ends
mid-condition at `664`. ⛔ **No part of it derives a dimension.** The dimensions are derived in the **`A4`**
block, lines `676-715` (`modulus_dimension`, `phase_speed_dimension`). ⇒ the dimension loci now point at
`A4`. The other loci keep their original ranges and changed **path only**.

⚠⚠ **And the two rows are not alike, which the repoint made visible.** At general `D` the `A4` block gives
`modulus_dimension = (2−D, −2, 1)` and `phase_speed_dimension = (1, −1, 0)`. ⇒ `Q.brane.c_L`'s declared
`(1, −1, 0)` matches at **every** `D`, while `Q.brane.B_comp`'s declared `(−1, −2, 1)` matches **only at
`D = 3`** — the general-`D` difference is `(3−D, 0, 0)`. ⭐ `B_comp`'s declared row is a `D = 3`
specialisation, ⛔ not a closed-form-in-`D` declaration, and a `Q6r` residual against it is a `D = 3`
statement.
