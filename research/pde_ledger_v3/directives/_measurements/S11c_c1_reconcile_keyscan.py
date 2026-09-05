import re, sys, collections
path=sys.argv[1]; want=set(sys.argv[2:])
# Walk CASE blocks: family, key, which side present, case_note
fam=key=None; opA=opB=None; note=None
rows=[]
def flush():
    if fam and (not want or fam in want):
        side = "JOIN" if (opA!="<MISSING>" and opB!="<MISSING>") else ("PY" if opB=="<MISSING>" else "WL")
        rows.append((fam,key,side,note))
with open(path) as f:
    for line in f:
        if line.startswith("CASE family="):
            flush()
            m=re.match(r"CASE family=(\S+) key=(.*)",line.rstrip("\n"))
            fam=m.group(1); key=m.group(2); opA=opB=note=None
        elif line.startswith("operand_A = "):
            opA="<MISSING>" if line.rstrip("\n")=="operand_A = <MISSING>" else "P"
        elif line.startswith("operand_B = "):
            opB="<MISSING>" if line.rstrip("\n")=="operand_B = <MISSING>" else "W"
        elif line.startswith("case_note = "):
            note=line.rstrip("\n")[12:]
    flush()
# Summarize per family: distinct key axis-signatures per side
for target in sorted(want) if want else sorted({r[0] for r in rows}):
    sub=[r for r in rows if r[0]==target]
    print(f"\n######## {target}  (cases={len(sub)}) ########")
    byside=collections.defaultdict(list)
    for fam,key,side,note in sub:
        byside[side].append((key,note))
    for side in ("PY","WL","JOIN"):
        items=byside.get(side,[])
        if not items: continue
        print(f"  --- {side} ({len(items)}) ---")
        seen=set()
        for key,note in items:
            # collapse the key to its axis names + short values
            sig=re.sub(r"\s+","",key)
            if sig in seen: continue
            seen.add(sig)
            if len(seen)<=6:
                print(f"    key={key}")
                if note: print(f"       note={note[:200]}")
