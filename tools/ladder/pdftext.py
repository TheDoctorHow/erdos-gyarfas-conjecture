import re, zlib, io, sys

path = sys.argv[1]
data = io.open(path, 'rb').read()
texts = []
for m in re.finditer(b'stream\r?\n(.*?)endstream', data, re.S):
    try:
        texts.append(zlib.decompress(m.group(1)))
    except Exception:
        pass

PAREN = re.compile(r'\((?:\\.|[^()\\])*\)', re.S)

def extract(b):
    s = b.decode('latin-1')
    return ' '.join(x[1:-1] for x in PAREN.findall(s))

full = ' '.join(extract(t) for t in texts)
full = re.sub(r'\\([()])', r'\1', full)
full = re.sub(r'\s+', ' ', full)
print('streams', len(texts), 'chars', len(full))
for kw in sys.argv[2:]:
    idxs = [m.start() for m in re.finditer(re.escape(kw), full, re.I)][:3]
    print('=== %s : %d hits' % (kw, len(idxs)))
    for i in idxs:
        print('   ...' + full[max(0, i-450):i+450] + '...')
        print()
