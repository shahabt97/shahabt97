"""Summary statistics of the shipped certificates (sizes of pieces, magnitude of entries) and fingerprints of the data."""
import collections
import gzip
import hashlib
import json
import os

DATA = os.path.join(os.path.dirname(__file__), '..', 'data')
for form in ('L', 'U'):
    path = os.path.join(DATA, f'cert_{form}.json.gz')
    raw = json.load(gzip.open(path, 'rt'))
    types = raw['types']
    sizes = collections.Counter(len(t['N']) for t in types)
    mx = max(abs(v) for t in types for row in t['N'] for v in row)
    print(f"certificate {form}: type radius Rt = {raw['Rt']}, entries are integers over 2^{raw['bits']}, {len(types)} pieces")
    print(f"  piece sizes (size: number of pieces): {dict(sorted(sizes.items()))}")
    if form == 'L':
        print(f"  pieces by child index j of the centre: {dict(sorted(collections.Counter(t['j'] for t in types).items()))}")
    print(f"  largest |entry| = {mx / 2 ** raw['bits']:.6g}")
    print(f"  sha256 {hashlib.sha256(open(path, 'rb').read()).hexdigest()}  (gzip file, {os.path.getsize(path)} bytes)")
    small = min(types, key=lambda t: (len(t['N']), t['j']))
    big = max(types, key=lambda t: len(t['N']))
    for name, t in (('smallest', small), ('largest', big)):
        S = 2 ** raw['bits']
        print(f"  example {name} piece: j = {t['j']}, {len(t['N'])}x{len(t['N'])}, offsets {t['offs'][:4]}{'...' if len(t['offs']) > 4 else ''}")
        for row in t['N'][:3]:
            print("     ", [round(v / S, 5) for v in row[:6]], '...' if len(row) > 6 else '')
