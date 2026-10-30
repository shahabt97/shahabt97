"""How the extras (flipped fourth children) are attached to the rest of the actual clusters.  The statement of the problem
says "attached as a matching"; the proof of FULL_PROOF.md does not use that phrase (only H0-H2).  For every extra of C_k the
number of lattice neighbours of the same species that are triangle vertices, and the number that are extras.
usage: python -m tests.attachments"""
import collections

from tests.cluster_numpy import chain_levels
from tests.crosscheck import cluster_data

for chain, L in (('A', 4), ('B', 3)):
    for k, prev, cur in chain_levels(chain, L):
        cs, adj, tri, parent, child, plab, ppos, n_old = cluster_data(prev, cur)
        stat = collections.Counter()
        for i, p in enumerate(cs):
            if child[p] == 3:
                nt = sum(1 for j in adj[i] if child[cs[j]] != 3)
                ne = sum(1 for j in adj[i] if child[cs[j]] == 3)
                stat[(nt, ne)] += 1
        print(f"chain {chain} k={k}: |C_k| = {len(cs)}; extras by (#triangle-vertex neighbours, #extra neighbours): "
              f"{dict(sorted(stat.items()))}")
