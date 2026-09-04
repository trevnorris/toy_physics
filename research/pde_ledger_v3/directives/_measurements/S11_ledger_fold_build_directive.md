# _measurements — S11_ledger_fold_build_directive.md

The fold+guard module directive is physics-free infrastructure implementing the committed design's §D2/§D3;
each provenance claim carries its command (rule 2).

## authority — the design is committed c04e071f
```
$ git -C /var/projects/toy_physics log --oneline -1 -- research/pde_ledger_v3/directives/export_ledger_bind_closure_design.md
c04e071f Export architecture redesign: the bind-closure LEDGER (design, 2-leg-gated)
```

## §D2 (fold: write-time F9, last-wins, no re-apply) and §D3 (guard) are the authority
```
$ grep -nE '^## D2 ·|^## D3 ·' research/pde_ledger_v3/directives/export_ledger_bind_closure_design.md
64:## D2 · TOPOLOGY — generate over a frozen base; F9 at write time; fold = last-wins
95:## D3 · UNDER-EXPORT GUARD — a manifest + edge-wise route/closure check (NOT D3-round-trip)
```

## the predecessor trap the guard closes (base has S11b flat face_response; c1 curved ⇒ F9c write-key)
```
$ grep -nE "^    'face_response':     \{" research/pde_ledger_v3/scripts/S11c_b_exports.py
3626:    'face_response':     {
```

## the unresolved routed-key contract the (c) honesty-check raises on
```
$ sed -n '205,206p' research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md
- ⛔⛔ **WHAT "BARE" MEANS TO A LATER READER.** ⚠ Measured: ⛔ no presently instantiated broken consumer, ⭐ and
  a latent one — after `F9c` the bare key still answers with the **predecessor**. ⇒ ⭐ **a routed-key
```

## the module writes only its two products; imports no engine export (physics-free)
```
$ grep -nE 'ledger_fold.py.*in full|only writes|imports no export|physics-free|infrastructure, not physics' research/pde_ledger_v3/directives/S11_ledger_fold_build_directive.md | head
5:Write `research/pde_ledger_v3/scripts/ledger_fold.py` in full, plus its self-test
6:`research/pde_ledger_v3/scripts/test_ledger_fold.py`. Those two files are the only writes.
8:`CLAUDE.md` binds. The physics-free **authority** is the committed design
11:write-time key routing; **F3** — a row carries its own evidence). This is **infrastructure, not physics**: it
104:itself imports no export and pins nothing.
```
