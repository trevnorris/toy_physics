"""Demonstrate the three baseline holes the fix round-1 directive targets."""
import importlib.util, sys
from pathlib import Path
SRC = Path("research/pde_ledger_v3/scripts/S11b_cross_engine_comparator.py")
spec = importlib.util.spec_from_file_location("cmp", SRC); m = importlib.util.module_from_spec(spec)
sys.modules["cmp"] = m; spec.loader.exec_module(m)
def chk(label, py, wl):
    r = m.compare_records("S11B_PROBE", m.Record("PY","S11B_PROBE","PY_S11B_PROBE",py,Path("/tmp/x"),1),
                          m.Record("WL","S11B_PROBE","WL_S11B_PROBE",wl,Path("/tmp/x"),1))
    print(f"{label}: {r.classification} (reason={r.reason})")
chk("D1 function-head collision", "(f_a(x), fA(x))", "{fA[x], fA[x]}")
chk("D2 status buries sibling", '((Str(\"STATUS_TOKEN\"), UNDECIDED), (Str(\"RESULT\"), a))', '<|\"STATUS_TOKEN\" -> UNDECIDED, \"RESULT\" -> b|>')
chk("D3 symbol-pair tuple promoted", "((p, 1), (q, 2))", '<|\"p\" -> 1, \"q\" -> 2|>')
