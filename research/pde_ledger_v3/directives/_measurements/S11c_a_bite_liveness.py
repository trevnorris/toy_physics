# Correct bite-liveness: a genuine bite = operand_A is a real nonzero sympy expression carrying a Symbol
# (the ablation/corruption moved PY's operand). Sentinels (TextAtom/Mismatch/UNDEFINED/<MISSING>) and
# literal zeros do NOT count. Run against the FINAL committed-state comparator run.
import re, collections
RUN = "/home/trevnorris/.s11_build/comparator_final_cccb4f9e.out"
def classify(fam):
    per = collections.defaultdict(lambda: {"BITE":0,"NULL":0,"SENTINEL":0})
    obj=None; want=False
    zero_re = re.compile(r'^\(?Integer\(0\)\)?(, Integer\(0\))*\)?$|^Integer\(0\)$')
    with open(RUN) as f:
        for line in f:
            if line.startswith(f"CASE family={fam} "):
                m=re.search(r'OBJECT=([A-Z_0-9]+)', line)
                obj = m.group(1) if m else "(noobj)"; want=True; continue
            if want and line.startswith("operand_A = "):
                v=line[len("operand_A = "):].strip()
                if v.startswith(("TextAtom","Mismatch")) or v=="<MISSING>" or "UNDEFINED" in v:
                    per[obj]["SENTINEL"]+=1
                elif zero_re.match(v):
                    per[obj]["NULL"]+=1
                elif "Symbol(" in v:
                    per[obj]["BITE"]+=1
                else:
                    per[obj]["NULL"]+=1
                want=False
    return per
for fam in ("CONTROL_FORM_RESIDUAL","CONTROL_INDEPENDENCE_RESIDUAL"):
    per=classify(fam)
    print(f"\n=== {fam} (genuine BITE = nonzero sympy w/ Symbol; SENTINEL excluded) ===")
    for obj in sorted(per):
        d=per[obj]
        verdict = "BITES" if d["BITE"]>0 else ("DEAD(only sentinels/null)" )
        print(f"  {obj:28s} BITE={d['BITE']:3d} NULL={d['NULL']:3d} SENTINEL={d['SENTINEL']:3d}  -> {verdict}")
