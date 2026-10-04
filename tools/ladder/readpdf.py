"""Print a window of extracted PDF text around a keyword."""
import re, zlib, io, sys

path, kw = sys.argv[1], sys.argv[2]
count = int(sys.argv[3]) if len(sys.argv) > 3 else 3
width = int(sys.argv[4]) if len(sys.argv) > 4 else 1200
data = io.open(path, 'rb').read()
texts = []
for m in re.finditer(b'stream\r?\n(.*?)endstream', data, re.S):
    try:
        texts.append(zlib.decompress(m.group(1)))
    except Exception:
        pass
PAREN = re.compile(r'\((?:\\.|[^()\\])*\)', re.S)
full = ' '.join(' '.join(x[1:-1] for x in PAREN.findall(t.decode('latin-1')))
                for t in texts)
full = re.sub(r'\\([()])', r'\1', full)
full = re.sub(r'\\[0-9]{3}', '~', full)
full = re.sub(r'\s+', ' ', full)
compact = re.sub(r'[^A-Za-z0-9]', '', full).lower()
mapping = [i for i, c in enumerate(full) if c.isalnum()]
k = re.sub(r'[^A-Za-z0-9]', '', kw).lower()
hits = [m.start() for m in re.finditer(re.escape(k), compact)][:count]
print('hits for %r: %d shown' % (kw, len(hits)))
for h in hits:
    a = mapping[max(0, h - width // 2)]
    b = mapping[min(len(mapping) - 1, h + len(k) + width)]
    print('...' + full[a:b].encode('ascii', 'replace').decode() + '...')
    print('-----')
