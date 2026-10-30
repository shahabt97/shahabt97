"""Exact-arithmetic verification of the two local certificates over all inflation cells.

For a quadratic form F on the vertices of a cluster (F = Q1 on fine vertices, F = 20*M on coarse sites) the certificate is
a family of positive definite matrices G_tau (tau = local type of a centre c) such that

        F = sum_c  iota_c^T G_{tau(c)} iota_c  +  R,        R_xx >= sum_{y != x} |R_xy| + rho   for every row x,

where iota_c embeds the ball W_c around the centre c.  Since every term of the sum is positive semidefinite and R is
diagonally dominant (Gershgorin), F >= rho * I.  The row inequality only involves the cell of x, so it is checked
for the finitely many cells.  Everything below is computed with python integers / Fractions.
"""
import argparse
import collections
import gzip
import json
import os
import time
from fractions import Fraction as Fr

import numpy as np

from .core import Frame, enc, euclid, hexd, rotk
from .cells import Cell, q1_entry, coarse_data, hball, m_row
from .closure import cells as build_cells, ZONE

DATA = os.path.join(os.path.dirname(__file__), '..', 'data')


def is_pd(M):
    """exact Gaussian elimination without pivoting: (True, smallest pivot) iff the symmetric matrix M is positive definite"""
    n = len(M)
    A = [row[:] for row in M]
    minp = None
    for k in range(n):
        p = A[k][k]
        if p <= 0:
            return False, p
        minp = p if minp is None or p < minp else minp
        for i in range(k + 1, n):
            if A[i][k] != 0:
                f = A[i][k] / p
                for j in range(k, n):
                    A[i][j] -= f * A[k][j]
    return True, minp


class Certificate:
    def __init__(self, path):
        with gzip.open(path, 'rt') as fh:
            raw = json.load(fh)
        self.form, self.Rt, self.bits = raw['form'], raw['Rt'], raw['bits']
        self.FR = Frame(self.Rt)
        self.pieces = {}
        for t in raw['types']:
            key = (t['j'], bytes.fromhex(t['key']))
            self.pieces[key] = dict(offs=tuple(tuple(o) for o in t['offs']), N=t['N'])

    def check_pieces(self):
        """every G_tau symmetric and positive definite (exact); returns the smallest pivot met"""
        S = Fr(1 << self.bits)
        minpiv = None
        for key, p in self.pieces.items():
            n = len(p['N'])
            assert len(p['offs']) == n and all(len(row) == n for row in p['N']), "piece shape"
            assert all(p['N'][i][j] == p['N'][j][i] for i in range(n) for j in range(i)), "piece not symmetric"
            ok, mp = is_pd([[Fr(v) / S for v in row] for row in p['N']])
            assert ok, "piece not positive definite"
            minpiv = mp if minpiv is None or mp < minpiv else minpiv
        return minpiv

    def type_of(self, d, o, j):
        """canonical type of the centre with child index j whose parent sits at offset o of the window d"""
        FR, Rt = self.FR, self.Rt
        assert hexd(*o) + Rt <= 7
        v = np.zeros(len(FR.offs), dtype=np.int8)
        for i, q in enumerate(FR.offs):
            l = d.get((o[0] + q[0], o[1] + q[1]))
            if l is not None:
                v[i] = enc(l[0], l[1], 0)
        key, k, s = FR.canon(v)
        return (j, key), k

    def piece(self, d, o, j):
        """the piece matrix (integers over 2^bits, canonical frame) of the centre, and the rotation k to that frame"""
        tkey, k = self.type_of(d, o, j)
        return self.pieces[tkey], k


def row_margin_L(cert, cell, x, pw=2):
    S = Fr(1 << cert.bits)
    cache = {}

    def centre(c):
        if c not in cache:
            v = cell.fine[c]
            p, k = cert.piece(cell.d, v[3], v[2])
            ball_ = cell.dist_ball(c, 2)
            coffs = {y: rotk((y[0] - c[0], y[1] - c[1]), k) for y in ball_}
            order = sorted(ball_, key=lambda y: coffs[y])
            assert tuple(coffs[y] for y in order) == p['offs'], "piece window mismatch"
            cache[c] = (p, {y: i for i, y in enumerate(order)})
        return cache[c]

    acc = collections.defaultdict(Fr)
    for c in cell.dist_ball(x, 2):
        p, idx = centre(c)
        ix = idx[x]
        for y, iy in idx.items():
            acc[y] += Fr(p['N'][ix][iy]) / S
    cols = set(acc) | set(cell.dist_ball(x, 3))
    R = {y: q1_entry(cell, x, y, pw) - acc.get(y, Fr(0)) for y in cols}
    return R[x] - sum(abs(v) for y, v in R.items() if y != x)


def row_margin_U(cert, cell, Wrad=2):
    S = Fr(1 << cert.bits)
    tri, site2tri, H = coarse_data(cell)
    z = (0, 0)
    cache = {}
    inside = lambda po: all(euclid(f) <= ZONE - 1 for f in tri[po])

    def centre(c):
        if c not in cache:
            p, k = cert.piece(cell.d, c, -1)
            ball_ = hball(H, c, Wrad, inside)
            coffs = {w: rotk((w[0] - c[0], w[1] - c[1]), k) for w in ball_}
            order = sorted(ball_, key=lambda w: coffs[w])
            assert tuple(coffs[w] for w in order) == p['offs'], "piece window mismatch"
            cache[c] = (p, {w: i for i, w in enumerate(order)})
        return cache[c]

    row = m_row(cell, tri, site2tri, z)
    acc = collections.defaultdict(Fr)
    for c in hball(H, z, Wrad, inside):
        p, idx = centre(c)
        ix = idx[z]
        for w, iw in idx.items():
            acc[w] += Fr(p['N'][ix][iw]) / S
    cols = set(acc) | set(row)
    R = {w: Fr(row.get(w, 0)) - acc.get(w, Fr(0)) for w in cols}
    return R[z] - sum(abs(v) for w, v in R.items() if w != z)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--all-seeds', action='store_true', help='also include the chains that start next to holes of the tiling')
    ap.add_argument('--data', default=DATA)
    args = ap.parse_args(argv)
    t0 = time.time()
    F, plain, seeds, wins = build_cells(7, interior_only=not args.all_seeds)
    nA = sum(1 for c, _ in seeds if c == 'A')
    print(f"all windows of the tiling (radius 7, mod G): {len(plain)}")
    print(f"seeds: {nA} singleton windows (chain A), {len(seeds) - nA} pair windows (chain B)")
    print(f"inflation cells (closure under inflation, mod G): {len(wins)}   [{time.time() - t0:.0f}s]")
    L = Certificate(os.path.join(args.data, 'cert_L.json.gz'))
    U = Certificate(os.path.join(args.data, 'cert_U.json.gz'))
    pL, pU = L.check_pieces(), U.check_pieces()
    print(f"lower edge: {len(L.pieces)} pieces, all symmetric and positive definite (exact), smallest pivot {float(pL):.3g}")
    print(f"upper edge: {len(U.pieces)} pieces, all symmetric and positive definite (exact), smallest pivot {float(pU):.3g}")
    rhoL = rhoU = None
    nrows = 0
    for n, v in enumerate(wins):
        cell = Cell(F.to_dict(v))
        for x in cell.owned:
            r = row_margin_L(L, cell, x)
            rhoL = r if rhoL is None or r < rhoL else rhoL
            nrows += 1
        r = row_margin_U(U, cell)
        rhoU = r if rhoU is None or r < rhoU else rhoU
        if n % 500 == 0:
            print(f"  cell {n}/{len(wins)}  [{time.time() - t0:.0f}s]", flush=True)
    print(f"lower edge: {nrows} rows checked; min row margin rho_L = {rhoL} = {float(rhoL):.6f} > 0")
    print(f"upper edge: {len(wins)} rows checked; min row margin rho_U = {rhoU} = {float(rhoU):.6f} > 0 (units of 20*M)")
    assert rhoL > 0 and rhoU > 0
    print(f"=> Q1 = 9/8 - A + 2 sum_v (A u_v)(A u_v)^T  >=  {float(rhoL):.4f} * I")
    print(f"=> M = U^T p5(A) U  >=  {float(rhoU) / 20:.4f} * I     (u_v = indicator of the triangle of v, not normalised)")
    print(f"VERIFIED in {time.time() - t0:.0f}s")


if __name__ == '__main__':
    main()
