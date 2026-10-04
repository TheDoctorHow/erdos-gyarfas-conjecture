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
full = re.sub(r'\\[0-9]{3}', '~', full)
full = re.sub(r'\s+', ' ', full)

# space-insensitive index
compact = re.sub(r'[^A-Za-z0-9]', '', full)
mapping = [i for i, c in enumerate(full) if re.match(r'[A-Za-z0-9]', c)]
print('chars', len(full), 'compact', len(compact))

for kw in sys.argv[2:]:
    k = re.sub(r'[^A-Za-z0-9]', '', kw).lower()
    hits = [m.start() for m in re.finditer(re.escape(k), compact.lower())][:4]
    print('=== %s : %d hits' % (kw, len(hits)))
    for h in hits:
        a = mapping[max(0, h-500)]
        b = mapping[min(len(mapping)-1, h+len(k)+500)]
        print('   ...' + full[a:b] + '...')
        print()
