"""Census of the finite 'interior' clusters of the level-K patch: components of the same-species nearest-neighbour graph
all of whose sites have all six lattice neighbours present in the patch.  Such a component is a genuine cluster of the
infinite tiling (no later site can occupy a position that is already occupied).  The sizes found are exactly those of the
two chains.   usage: python -m tests.census [K]"""
import collections
import sys

from ww.core import D, build_patch

K = int(sys.argv[1]) if len(sys.argv) > 1 else 9
patch = build_patch(K)
seen = set()
sizes = collections.Counter()
for p in patch:
    if p in seen:
        continue
    t = patch[p][1]
    comp, stack = [p], [p]
    seen.add(p)
    while stack:
        x = stack.pop()
        for d in D:
            y = (x[0] + d[0], x[1] + d[1])
            if y in patch and y not in seen and patch[y][1] == t:
                seen.add(y)
                stack.append(y)
                comp.append(y)
    if all((s[0] + d[0], s[1] + d[1]) in patch for s in comp for d in D):
        sizes[len(comp)] += 1
A = [1, 4, 16, 63, 245, 941, 3583, 13574, 51293, 193605]
B = [2, 10, 42, 165, 632, 2396, 9049, 34135]
print(f"patch level {K}: {len(patch)} sites; interior finite clusters by size (size: number):")
print(dict(sorted(sizes.items())))
extra = set(sizes) - set(A) - set(B)
print("sizes outside the two chains A, B:", sorted(extra))
assert not extra
