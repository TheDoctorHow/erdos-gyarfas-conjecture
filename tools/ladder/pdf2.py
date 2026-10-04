import re, zlib, sys
path = sys.argv[1]
data = open(path,"rb").read()
out=[]
for m in re.finditer(rb"stream\r?\n", data):
    s=m.end(); e=data.find(b"endstream", s)
    if e<0: continue
    try: dec=zlib.decompress(data[s:e])
    except Exception: continue
    if b"TJ" not in dec and b"Tj" not in dec: continue
    txt=[]
    for tm in re.finditer(rb"\[(.*?)\]\s*TJ|\((.*?)\)\s*Tj", dec, re.S):
        if tm.group(1) is not None:
            txt.append(b"".join(re.findall(rb"\((.*?)(?<!\\)\)", tm.group(1), re.S)))
        else: txt.append(tm.group(2))
    s2=b" ".join(txt).replace(b"\\(",b"(").replace(b"\\)",b")").replace(b"\\\\",b"\\")
    try: out.append(s2.decode("latin-1"))
    except Exception: pass
full=re.sub(r"[ \t]+"," ","\n".join(out))
print("chars:",len(full))
kw=[k.strip() for k in sys.argv[2].split(",")] if len(sys.argv)>2 else []
if kw:
    for i,line in enumerate(full.split("\n")):
        if any(k.lower() in line.lower() for k in kw):
            print("---",i); print(line[:2500])
else:
    print(full[:4000])
