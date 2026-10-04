import re,zlib,sys
p=sys.argv[1]
raw=open(p,'rb').read()
out=[]
for m in re.finditer(rb'stream\r?\n(.*?)endstream',raw,re.S):
    s=m.group(1)
    try: s=zlib.decompress(s)
    except Exception: continue
    txt=[]
    for t in re.finditer(rb'\[(.*?)\]\s*TJ|\((.*?)\)\s*Tj',s,re.S):
        seg=t.group(1) or t.group(2) or b''
        for q in re.finditer(rb'\((.*?)(?<!\\)\)',seg) if t.group(1) else [t]:
            try: txt.append((q.group(1) if t.group(1) else seg).decode('latin-1'))
            except Exception: pass
    if txt: out.append(''.join(txt))
doc='\n'.join(out)
doc=doc.replace('\\(','(').replace('\\)',')')
kw=[k.lower() for k in sys.argv[2:]] or ['arkin','modular','path']
lines=doc.split('\n')
for i,l in enumerate(lines):
    if any(k in l.lower() for k in kw):
        print('---',i)
        print('\n'.join(lines[max(0,i-2):i+3])[:1800])
