"""Fetch only CI's ten pinned models and verify the committed corpus manifest."""
import hashlib
import urllib.request
from pathlib import Path
pin = '56257eea85b433ce6aa67d26156b36385318fd6f'
names = 'afiro sc50a sc50b sc105 adlittle blend kb2 share2b stocfor1 scagr7'.split()
manifest = dict(reversed(line.split()) for line in Path('reproduction/MANIFEST_corpus_sha256.txt').read_text().splitlines() if line.strip())
root = Path('out/ci/corpus'); root.mkdir(parents=True, exist_ok=True)
for name in names:
    filename = name + '.mps'
    path = root / filename
    if not path.exists():
        path.write_bytes(urllib.request.urlopen('https://raw.githubusercontent.com/ozy4dm/lp-data-netlib/' + pin + '/mps_files/' + filename, timeout=30).read())
    assert hashlib.sha256(path.read_bytes()).hexdigest() == manifest[filename], filename + ': corpus hash mismatch'
print('Ten pinned corpus hashes verified')
