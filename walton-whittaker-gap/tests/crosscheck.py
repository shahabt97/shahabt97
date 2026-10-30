"""Independent numerical / exact cross-checks on the actual clusters C_k (brute-force numpy construction).

For chain A (levels 1..LA) and chain B (levels 1..LB) it
  * builds the cluster graph by brute-force inflation (tests/cluster_numpy.py, no cell machinery),
  * counts the eigenvalues of A(C_k) above 9/8 and in [9/8, 23/20] (dense float eigensolver; a float is NOT the proof),
  * assembles the forms Q1 and 20*M globally and checks their smallest eigenvalues against the certified margins,
  * assembles the certificate pieces globally (exact Fractions) and recomputes the Gershgorin margin of the residual
    R = F - sum_c iota_c^T G_c iota_c on the whole cluster: it must reproduce (be >=) the margins certified on the cells.
  * checks the structural hypotheses H1, H2 of FULL_PROOF.md on the actual clusters: the triangle vertices of C_k are exactly
    the triangles of the sites of C_{k-1}, and every other vertex is the flipped fourth child of a lattice neighbour of a site of C_{k-1}.
Run:  python -m tests.crosscheck [LA LB]              (from the package root)
      python -m tests.crosscheck LA LB structure      (only the structural checks, much faster, deeper levels)
"""
import collections
import os
import sys
import time
from fractions import Fraction as Fr

import numpy as np

from ww.core import D, child_pos, enc, rotk, Frame
from ww.cells import C20
from ww.verify import Certificate, DATA
from tests.cluster_numpy import chain_levels

RHO_L = Fr(85844911, 1073741824)        # certified by ww.verify for the shipped certificates
RHO_U = Fr(696063391, 8388608)


def cluster_data(prev, cur):
    pa, pb, pr, pt, pmask = prev
    a, b, r, t, mask, par, ch = cur
    cs = [(int(x), int(y)) for x, y, mk in zip(a, b, mask) if mk]
    idx = {p: i for i, p in enumerate(cs)}
    parent = {(int(x), int(y)): int(p) for x, y, mk, p in zip(a, b, mask, par) if mk}
    child = {(int(x), int(y)): int(c) for x, y, mk, c in zip(a, b, mask, ch) if mk}
    plab = {(int(x), int(y)): (int(rr), int(tt)) for x, y, rr, tt in zip(pa, pb, pr, pt)}
    ppos = {i: (int(x), int(y)) for i, (x, y) in enumerate(zip(pa, pb))}
    n_old = int(pmask.sum())
    # structural hypotheses H1, H2 (FULL_PROOF.md, section 1.3)
    old = {(int(x), int(y)) for x, y, mk in zip(pa, pb, pmask) if mk}
    for p in cs:
        w = ppos[parent[p]]
        if child[p] != 3:
            assert w in old, "H1 violated: triangle vertex whose parent is not in C_(k-1)"
        else:
            assert any((w[0] + d[0], w[1] + d[1]) in old for d in D), "H2 violated: extra whose parent is not adjacent to C_(k-1)"
    for v in old:
        assert all(child_pos(v, plab[v][0], j) in idx for j in range(3)), "H1 violated: triangle of a site of C_(k-1) not in C_k"
    adj = [[idx[(p[0] + dl[0], p[1] + dl[1])] for dl in D if (p[0] + dl[0], p[1] + dl[1]) in idx] for p in cs]
    tri = collections.defaultdict(list)
    for p in cs:
        if child[p] != 3:
            tri[parent[p]].append(idx[p])
    return cs, adj, tri, parent, child, plab, ppos, n_old


def dense_forms(adj, tri):
    n = len(adj)
    A = np.zeros((n, n))
    for i, nb in enumerate(adj):
        A[i, nb] = 1
    tids = sorted(tri)
    U = np.zeros((n, len(tids)))
    for k, t in enumerate(tids):
        U[tri[t], k] = 1
    AU = A @ U
    Q1 = 9 / 8 * np.eye(n) - A + 2 * AU @ AU.T
    P = sum(c * np.linalg.matrix_power(A, i) for i, c in enumerate(C20))
    M20 = U.T @ P @ U
    return A, Q1, np.rint(M20), tids


def canon_plain(FR, plab, p):
    v = np.zeros(len(FR.offs), dtype=np.int8)
    for i, o in enumerate(FR.offs):
        l = plab.get((p[0] + o[0], p[1] + o[1]))
        if l is not None:
            v[i] = enc(l[0], l[1], 0)
    return FR.canon(v)


def ball_dist(nbrs, x, rad):
    dist = {x: 0}
    fr = [x]
    for s in range(rad):
        nf = []
        for y in fr:
            for g in nbrs[y]:
                if g not in dist:
                    dist[g] = s + 1
                    nf.append(g)
        fr = nf
    return dist


def margin_L(cert, cs, adj, tri, parent, child, plab, ppos):
    n = len(cs)
    Q = [collections.defaultdict(Fr) for _ in range(n)]
    for i in range(n):
        Q[i][i] += Fr(9, 8)
        for j in adj[i]:
            Q[i][j] -= 1
    for v, members in tri.items():
        Au = collections.Counter()
        for t in members:
            for j in adj[t]:
                Au[j] += 1
        for i, ci in Au.items():
            for j, cj in Au.items():
                Q[i][j] += 2 * ci * cj
    S = Fr(1 << cert.bits)
    Fh = [collections.defaultdict(Fr) for _ in range(n)]
    for c, p in enumerate(cs):
        key, k, s = canon_plain(cert.FR, plab, ppos[parent[p]])
        piece = cert.pieces[(child[p], key)]
        dist = ball_dist(adj, c, 2)
        coffs = {y: rotk((cs[y][0] - p[0], cs[y][1] - p[1]), k) for y in dist}
        order = sorted(dist, key=lambda y: coffs[y])
        assert tuple(coffs[y] for y in order) == piece['offs']
        for a_, y in enumerate(order):
            for b_, w in enumerate(order):
                Fh[y][w] += Fr(piece['N'][a_][b_]) / S
    worst = None
    for i in range(n):
        R = {w: Q[i].get(w, Fr(0)) - Fh[i].get(w, Fr(0)) for w in set(Q[i]) | set(Fh[i])}
        m = R[i] - sum(abs(v) for w, v in R.items() if w != i)
        worst = m if worst is None or m < worst else worst
    return worst


def margin_U(cert, cs, adj, tri, plab, ppos, M20, tids):
    m = len(tids)
    site2tri = {s: k for k, t in enumerate(tids) for s in tri[t]}
    H = [set() for _ in range(m)]
    for k, t in enumerate(tids):
        for s in tri[t]:
            for g in adj[s]:
                if g in site2tri and site2tri[g] != k:
                    H[k].add(site2tri[g])
    S = Fr(1 << cert.bits)
    Mh = [collections.defaultdict(Fr) for _ in range(m)]
    for k, t in enumerate(tids):
        p = ppos[t]
        key, rk, s = canon_plain(cert.FR, plab, p)
        piece = cert.pieces[(-1, key)]
        dist = ball_dist(H, k, 2)
        coffs = {y: rotk((ppos[tids[y]][0] - p[0], ppos[tids[y]][1] - p[1]), rk) for y in dist}
        order = sorted(dist, key=lambda y: coffs[y])
        assert tuple(coffs[y] for y in order) == piece['offs']
        for a_, y in enumerate(order):
            for b_, w in enumerate(order):
                Mh[y][w] += Fr(piece['N'][a_][b_]) / S
    worst = None
    for i in range(m):
        nz = np.nonzero(M20[i])[0]
        R = {w: Fr(int(M20[i, w])) - Mh[i].get(w, Fr(0)) for w in set(nz.tolist()) | set(Mh[i])}
        mm = R[i] - sum(abs(v) for w, v in R.items() if w != i)
        worst = mm if worst is None or mm < worst else worst
    return worst


def main():
    LA = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    LB = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    if len(sys.argv) > 3 and sys.argv[3] == 'structure':
        for chain, Lmax in (('A', LA), ('B', LB)):
            for k, prev, cur in chain_levels(chain, Lmax):
                cs, adj, tri, parent, child, plab, ppos, n_old = cluster_data(prev, cur)
                assert len(tri) == n_old
                print(f"chain {chain} k={k}: |C_k|={len(cs)}, |C_(k-1)|={n_old}: H1, H2 hold", flush=True)
        return
    L = Certificate(os.path.join(DATA, 'cert_L.json.gz'))
    U = Certificate(os.path.join(DATA, 'cert_U.json.gz'))
    lam1, lam2 = 9 / 8, 23 / 20
    for chain, Lmax in (('A', LA), ('B', LB)):
        sizes = []
        for k, prev, cur in chain_levels(chain, Lmax):
            t0 = time.time()
            cs, adj, tri, parent, child, plab, ppos, n_old = cluster_data(prev, cur)
            n, m = len(cs), len(tri)
            assert m == n_old, "triangles of C_k <-> sites of C_{k-1}"
            A, Q1, M20, tids = dense_forms(adj, tri)
            w = np.linalg.eigvalsh(A)
            above = int((w >= lam1).sum())
            inwin = int(((w >= lam1) & (w <= lam2)).sum())
            top = w[w > lam2].min()
            bulk = w[w < lam1].max()
            eq1 = np.linalg.eigvalsh(Q1).min()
            em = np.linalg.eigvalsh((M20 + M20.T) / 2).min()
            ml = margin_L(L, cs, adj, tri, parent, child, plab, ppos)
            mu = margin_U(U, cs, adj, tri, plab, ppos, M20, tids)
            ok = (above == m and inwin == 0 and ml >= RHO_L and mu >= RHO_U and eq1 >= float(RHO_L) and em >= float(RHO_U))
            print(f"chain {chain} k={k}: |C_k|={n}, |C_(k-1)|={m}; #eig>=9/8: {above}; #eig in [9/8,23/20]: {inwin}; "
                  f"max bulk {bulk:.4f}, min band {top:.4f}; lmin(Q1)={eq1:.4f} (>= {float(RHO_L):.4f}); "
                  f"lmin(20M)={em:.2f} (>= {float(RHO_U):.2f}); global residual margins {float(ml):.4f}, {float(mu):.3f}  "
                  f"{'OK' if ok else 'FAIL'}  [{time.time() - t0:.0f}s]", flush=True)
            assert ok
    print("cross-checks passed")


if __name__ == '__main__':
    main()
