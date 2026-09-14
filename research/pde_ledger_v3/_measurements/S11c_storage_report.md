# Durable S11c working data

The user required important calculation files to live in the repository, not
in `/tmp`. The active endpoint controller and its RIGHT pairing child were
suspended while their data and retained dependencies were moved on the same
filesystem into `/var/projects/toy_physics/_scratch/s11c/`.

The migration covers **35 run directories and 80 adjacent files**, containing
**10,718 regular files and 13,891,191,404 bytes**. Every regular file's SHA256,
size and inode was compared before and after the move; all matched. There
were no symlinks inside these runs. Dependency discovery found no missing
referenced S11c run directories. The tracked
[migration checkpoint](S11c_storage_migration_checkpoint.json) records each
old/durable path and tree digest. Its complete per-file inventory is durably
stored with the migration records, with SHA256
`68665cc3e34dc431ddfe04d43292467b057dda943c831d99cdd4f58353f094de`.

The active run is now physically located at:

`/var/projects/toy_physics/_scratch/s11c/s11c-thickness-coordinate-20260914`

Old `/tmp` names are compatibility symlinks only. Keeping those names preserves
absolute paths embedded in frozen manifests, source signatures and scalar
cache provenance without rewriting historical evidence. New runs must use
repository paths directly. If cleanup removes the compatibility links, run
from the repository root:

```bash
python research/pde_ledger_v3/_measurements/S11c_storage_paths.py --restore
```

The helper preflights the whole mapping, refuses conflicts, and creates only
missing links. Recovery fixtures verified link recreation after removal,
preservation of the durable payload, and rejection of a conflicting file.
All 115 actual compatibility paths and the active controller's 32 source/input
pins passed checks after relocation. The approved recurring follow-up now
discovers the controller through the durable repository location.

Working data remains local and Git-ignored under the existing scratch policy;
validated published `.out` files remain versioned through DataLad/git-annex.
Code, reports, inventories and recovery tooling are in Git. The migration does
not change any physical construction, published result or source pin. The
storage pause is recorded separately because the running calculation's elapsed
wall time includes it.
