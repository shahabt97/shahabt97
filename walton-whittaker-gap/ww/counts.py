"""Cluster sizes |C_k| from the substitution on inflation cells, their linear recurrence, and the gap labels.

Every site of C_{k-1} contributes the 3 vertices of its triangle; every extra (flipped fourth child of a neighbouring site
w of the other species) is attributed to the neighbour z* of w, inside the same cluster, with the smallest
(dir(w -> z') - r_w) mod 6.  The cell of z determines the cells of all attributed new members, so the sizes obey a
linear recurrence; the Perron root gives lim |C_{k-1}|/|C_k|, i.e. the integrated density of states above the gap.
"""
import collections
import sys
from fractions import Fraction as Fr

import numpy as np

from .core import D, Frame, child_pos, inflate, enc
from .closure import cells as build_cells, same_species_nbrs, component_from

REC = [8, -22, 25, -7, -5, 2]       # s_{k+6} = 8 s_{k+5} - 22 s_{k+4} + 25 s_{k+3} - 7 s_{k+2} - 5 s_{k+1} + 2 s_k
CHARPOLY = "x^6 - 8x^5 + 22x^4 - 25x^3 + 7x^2 + 5x - 2 = (x-1)(x^5 - 7x^4 + 15x^3 - 10x^2 - 3x + 2)"


def children_types(F, d, type_index):
    """types (indices into the cell list) of the new members attributed to the centre of the cell d"""
    fine = inflate(d)
    z = (0, 0)
    r0, t0, _ = d[z]
    tri = [child_pos(z, r0, j) for j in range(3)]
    comp = component_from(fine, tri)
    cands = list(tri)
    for q in D:
        if q in d and d[q][1] != t0:
            f = child_pos(q, d[q][0], 3)
            if f not in comp:
                continue
            rw = d[q][0]
            members = []
            for dl in D:
                zz = (q[0] + dl[0], q[1] + dl[1])
                if zz in d and d[zz][1] == t0 and all(child_pos(zz, d[zz][0], j) in comp for j in range(3)):
                    members.append(zz)
            assert z in members
            zstar = min(members, key=lambda zz: (D.index((zz[0] - q[0], zz[1] - q[1])) - rw) % 6)
            if zstar == z:
                cands.append(f)
    out = []
    for f in cands:
        w = {}
        for o in F.offs:
            g = (f[0] + o[0], f[1] + o[1])
            if g in fine:
                w[o] = (fine[g][0], fine[g][1], 0)
        k, _, _ = F.canon(F.from_dict(w))
        out.append(type_index[k])
    return out


def sizes(trans, start, K):
    vec = dict(start)
    s = [sum(vec.values())]
    for _ in range(K):
        nv = collections.Counter()
        for t, c in vec.items():
            for ch in trans[t]:
                nv[ch] += c
        vec = nv
        s.append(sum(vec.values()))
    return s


def krylov_dimension(trans, start, maxdim=60):
    """exact dimension of span{N_0, m N_0, m^2 N_0, ...}: the vectors N_k of cell counts live in this space, so every
    linear functional of N_k (e.g. |C_k|) obeys a linear recurrence of at most this order"""
    vec = {t: Fr(c) for t, c in start.items()}
    basis = []
    for _ in range(maxdim):
        v = dict(vec)
        for piv, b in basis:
            c = v.get(piv, 0)
            if c:
                for k, x in b.items():
                    v[k] = v.get(k, 0) - c * x
                    if v[k] == 0:
                        del v[k]
        if not v:
            return len(basis)
        piv = min(v)
        c = v[piv]
        basis.append((piv, {k: x / c for k, x in v.items()}))
        nv = {}
        for t, c in vec.items():
            for ch in trans[t]:
                nv[ch] = nv.get(ch, 0) + c
        vec = nv
    return None


def main():
    F, plain, seeds, wins = build_cells(7, interior_only=True)
    type_index = {v.tobytes(): i for i, v in enumerate(wins)}
    trans = [children_types(F, F.to_dict(v), type_index) for v in wins]
    A_seeds = [v for c, v in seeds if c == 'A']
    seqA = {tuple(sizes(trans, {type_index[v.tobytes()]: 1}, 60)) for v in A_seeds}
    assert len(seqA) == 1, "all singleton seeds give the same sizes"
    sA = list(seqA.pop())
    kdA = {krylov_dimension(trans, {type_index[v.tobytes()]: 1}, maxdim=200) for v in A_seeds}
    print(f"exact Krylov dimension of the substitution started from a chain-A seed (all {len(A_seeds)} seeds): {sorted(kdA)}  "
          f"(=> |C_k| obeys a linear recurrence of order <= {max(kdA)} for ALL k; the minimal one has order 6)")
    # chain B: a pair seed = two member windows that see each other (matched through their radius-6 windows)
    F6 = Frame(6)

    def sub6(d, o):
        v = np.zeros(len(F6.offs), dtype=np.int8)
        for i, q in enumerate(F6.offs):
            l = d.get((o[0] + q[0], o[1] + q[1]))
            if l is not None:
                v[i] = enc(l[0], l[1], 0)
        return F6.canon(v)[0]

    B = [v for c, v in seeds if c == 'B']
    info = []
    for v in B:
        d = F.to_dict(v)
        info.append((v, sub6(d, (0, 0)), sub6(d, same_species_nbrs(d, (0, 0))[0])))
    seqB, kdB = set(), set()
    for i, (v, c6, p6) in enumerate(info):
        for j, (w, c6j, p6j) in enumerate(info):
            if c6j == p6 and p6j == c6:
                start = collections.Counter([type_index[v.tobytes()], type_index[w.tobytes()]])
                seqB.add(tuple(sizes(trans, start, 60)))
                kdB.add(krylov_dimension(trans, start, maxdim=200))
    assert len(seqB) == 1
    sB = list(seqB.pop())
    print(f"exact Krylov dimension of the substitution started from a chain-B seed (all matched pairs): {sorted(kdB)}")
    assert max(kdA | kdB) + 6 <= len(sA), "too few terms to conclude that the recurrence holds for all k"
    for name, s in (('A', sA), ('B', sB)):
        assert all(s[k + 6] == sum(c * s[k + 5 - i] for i, c in enumerate(REC)) for k in range(len(s) - 6))
        print(f"chain {name}: |C_k| = {s[:10]} ...   (satisfies the recurrence with characteristic polynomial {CHARPOLY})")
    quintic = [1, -7, 15, -10, -3, 2]
    theta = max(r.real for r in np.roots(quintic) if abs(r.imag) < 1e-12)
    print(f"Perron root theta = {theta:.12f},  lim |C_(k-1)|/|C_k| = 1/theta = {1 / theta:.12f},  1 - 1/theta = {1 - 1 / theta:.12f}")
    print("finite-k labels 1 - |C_(k-1)|/|C_k|  (A):", [str(1 - Fr(sA[k - 1], sA[k])) for k in range(1, 8)])
    print("finite-k labels 1 - |C_(k-1)|/|C_k|  (B):", [str(1 - Fr(sB[k - 1], sB[k])) for k in range(1, 8)])
    try:
        import sympy as sp
        x = sp.symbols('x')
        print("quintic irreducible over Q:", sp.Poly(sum(c * x ** (5 - i) for i, c in enumerate(quintic)), x).is_irreducible)
    except ImportError:
        pass


if __name__ == '__main__':
    main()
