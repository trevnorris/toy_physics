# Reduction measurements

- `q7_index_map_discrimination.py` compares the declared and transposed gradient maps over every emitted D=3 §Q7 pair.
- `q7_payload_invariance_group.py` exhaustively measures the gradient-symbol permutation group admitted by every emitted D=3 §Q7 density.
- `q7_map_diagonal_gap.py` measures whether cycling only the diagonal gradient declarations survives the existing gates.
- `q7_map_declaration_ablation.py` reports §Q7 verdicts after transposing, deranging, dropping, or removing the gradient map.
- `declaration_load_ablation.py` measures the cross-engine rows and counters carried by each naming or identity declaration.
- `stiffness_coefficient_recovery.py` recovers each package's stiffness coefficient from its emitted Lagrangian and §Q7 density.
- `q8_stratum_manual_comparison.py` pairs and prints the engines' otherwise-undeclared Q8 stratum count payloads.
- `bare_word_payload_scan.py` inventories emitted payloads that are bare words rather than CAS objects.
- `q3_sign_adjudication.py` asks SymPy for each disputed Q3 root sign under the declared premises.
- `route_row_information_content.py` measures the downstream-route row's verdict under representative replacement payloads.
