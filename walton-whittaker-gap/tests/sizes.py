"""Brute-force sizes |C_k| of both chains (numpy inflation of the actual clusters), for comparison with the exact
substitution of ww.counts.   usage: python -m tests.sizes [LA LB]"""
import sys
from tests.cluster_numpy import chain_levels

LA = int(sys.argv[1]) if len(sys.argv) > 1 else 9
LB = int(sys.argv[2]) if len(sys.argv) > 2 else 8
for chain, L in (('A', LA), ('B', LB)):
    sizes = []
    for k, prev, cur in chain_levels(chain, L):
        if k == 1:
            sizes.append(int(prev[4].sum()))
        sizes.append(int(cur[4].sum()))
    print(f"chain {chain}: |C_0..C_{L}| =", sizes, flush=True)
