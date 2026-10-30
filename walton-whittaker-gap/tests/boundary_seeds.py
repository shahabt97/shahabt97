"""Chains that start next to a hole of the tiling (singleton clusters with a lattice neighbour that is not a site): sizes
|C_k| from the substitution of ww.counts, using the same cells as ww.verify --all-seeds.  These chains are "on the boundary"
and are not part of the statement; they are listed only to show how the interior chains differ from them.
usage: python -m tests.boundary_seeds"""
import collections

from ww.closure import cells as build_cells
from ww.counts import children_types, sizes

F, plain, seeds, wins = build_cells(7, interior_only=False)
_, _, seeds_in, _ = build_cells(7, interior_only=True)
interior = {v.tobytes() for _, v in seeds_in}
type_index = {v.tobytes(): i for i, v in enumerate(wins)}
trans = [children_types(F, F.to_dict(v), type_index) for v in wins]
seq = collections.Counter()
for c, v in seeds:
    if c == 'A' and v.tobytes() not in interior:
        seq[tuple(sizes(trans, {type_index[v.tobytes()]: 1}, 7))] += 1
print("singleton seed windows adjacent to a hole:", sum(seq.values()))
for s, n in sorted(seq.items(), reverse=True):
    print(f"  {n:3d} window classes: |C_0..C_7| = {list(s)}")
