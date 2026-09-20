# Preserve repeated complete current checkpoints

The original current production stopped after 24.77 seconds. Its coordinator
saved the accepted reference acoustic packet, then attempted to save it again
at the common case-completion step. The strict create-once writer rejected that
second write. No new acoustic or pairing construction had started. Preserve the
original 93 copied inputs, that complete acoustic packet and all original logs.

Keep both existing current helpers and all native constructors unchanged. The
recovery uses the identical coordinator AST/bytecode and replaces only its
packet-writer binding and saved-input loader. All other file/validation functions
retain their original identities. For `acoustic.pickle` and `pairing.pickle`
only, an existing packet must equal the entire requested payload with the same
literal type/content comparator. A mismatch rejects without touching the file.
All other packets keep the original strict create-once behavior. Save the actual
requested/existing pair and both hashes for every repeated current write.

Copy the completed original acoustic packet byte-for-byte after its full join
to the original accepted result. Reproduce every original input copy and exact
source pin. No baseline current, source, energy, factor, quadrature or mode is
reconstructed. A bounded check runs the complete original LAB_HELD/RHO4 packet
validation at its three backgrounds, including full current/pairing dictionaries
and unit joins. It computes no new closure or pairing. Actual changed acoustic
surface density and pairing matrix controls must reject, leaving saved bytes
unchanged. Retain this limited baseline scope explicitly.

After that check passes, run all twelve addresses with the accepted source
routing and compute only the one genuinely new RHOBR-right acoustic/pairing
family. Retain full rational operands, denominators, cross-product certificates,
face/radical/retained-balance checks and 5-by-5 independent-frequency maps.
Require clean supervisor exit, empty stderr, checks/stdout identity and all
original/copied/source/pre-post hashes. The final one-new-family guard must pass
after checks are emitted. Preserve any later incomplete packet and finish only
remaining work. Use one native thread, 2 GiB, 900 seconds, one supervisor and a
silent local completion/error hook. Keep approved inputs, unchanged engine,
pinned exports, retained builder suffix and practical scope restrictions.
