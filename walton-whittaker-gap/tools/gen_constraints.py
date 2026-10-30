"""Generate the linear constraints of the piece SDP over all inflation cells.

Lower edge (L):  pieces G_tau sit on the ball of graph radius 2 around a fine vertex c; type tau = (child index j of c,
                 canonical radius-Rt window of the parent of c).
Upper edge (U):  pieces sit on the ball of radius 2 (triangle adjacency graph H) around a coarse site c; type tau =
                 canonical radius-Rt window of c.
For every row x of every cell the equation  sum_c G_{tau(c)}[x, y] = F_xy  (minus a slack eps on the diagonal) is one
constraint; constraints with the same left-hand side and right-hand side are merged.

usage: python -m tools.gen_constraints L 2 out.pkl [--all-seeds]       (Rt = 2 for L, 4 for U)
"""
import argparse
import collections
import pickle
import sys
from fractions import Fraction as Fr

import numpy as np

from ww.core import Frame, enc, hexd, rotk
from ww.cells import Cell, q1_entry, coarse_data, hball, m_row
from ww.closure import cells as build_cells


class TypeBook:
    def __init__(self, Rt):
        self.Rt, self.FR = Rt, Frame(Rt)
        self.ids, self.keys, self.windows = {}, [], []

    def type_of(self, d, o, j):
        FR = self.FR
        assert hexd(*o) + self.Rt <= 7
        v = np.zeros(len(FR.offs), dtype=np.int8)
        for i, q in enumerate(FR.offs):
            l = d.get((o[0] + q[0], o[1] + q[1]))
            if l is not None:
                v[i] = enc(l[0], l[1], 0)
        key, k, s = FR.canon(v)
        full = (j, key)
        if full not in self.ids:
            self.ids[full] = len(self.ids)
            self.keys.append(full)
            self.windows.append(None)
        return self.ids[full], k

    def set_window(self, tid, win):
        if self.windows[tid] is None:
            self.windows[tid] = win
        assert self.windows[tid] == win, "the type does not determine the piece window: increase Rt"


def rows_L(tb, cell):
    cinfo = {}

    def centre(c):
        if c not in cinfo:
            v = cell.fine[c]
            tid, k = tb.type_of(cell.d, v[3], v[2])
            ball_ = cell.dist_ball(c, 2)
            coffs = {y: rotk((y[0] - c[0], y[1] - c[1]), k) for y in ball_}
            order = sorted(ball_, key=lambda y: coffs[y])
            tb.set_window(tid, tuple(coffs[y] for y in order))
            cinfo[c] = (tid, {y: i for i, y in enumerate(order)})
        return cinfo[c]

    for x in cell.owned:
        contrib = collections.defaultdict(list)
        for c in cell.dist_ball(x, 2):
            tid, idx = centre(c)
            for y, iy in idx.items():
                contrib[y].append((tid, idx[x], iy))
        for y in set(contrib) | set(cell.dist_ball(x, 3)):
            yield q1_entry(cell, x, y), tuple(sorted(contrib.get(y, []))), y == x


def rows_U(tb, cell, Wrad=2):
    tri, site2tri, H = coarse_data(cell)
    z = (0, 0)
    cinfo = {}

    def centre(c):
        if c not in cinfo:
            tid, k = tb.type_of(cell.d, c, -1)
            ball_ = hball(H, c, Wrad)
            coffs = {w: rotk((w[0] - c[0], w[1] - c[1]), k) for w in ball_}
            order = sorted(ball_, key=lambda w: coffs[w])
            tb.set_window(tid, tuple(coffs[w] for w in order))
            cinfo[c] = (tid, {w: i for i, w in enumerate(order)})
        return cinfo[c]

    row = m_row(cell, tri, site2tri, z)
    contrib = collections.defaultdict(list)
    for c in hball(H, z, Wrad):
        tid, idx = centre(c)
        for w, iw in idx.items():
            contrib[w].append((tid, idx[z], iw))
    for w in set(contrib) | set(row):
        assert not (row.get(w, 0) != 0 and w not in contrib), "uncovered nonzero entry"
        yield Fr(row.get(w, 0)), tuple(sorted(contrib.get(w, []))), w == z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('form', choices=['L', 'U'])
    ap.add_argument('Rt', type=int)
    ap.add_argument('out')
    ap.add_argument('--all-seeds', action='store_true')
    a = ap.parse_args()
    F, plain, seeds, wins = build_cells(7, interior_only=not a.all_seeds)
    tb = TypeBook(a.Rt)
    sigs = {}
    for n, v in enumerate(wins):
        cell = Cell(F.to_dict(v))
        for rhs, lst, isdiag in (rows_L(tb, cell) if a.form == 'L' else rows_U(tb, cell)):
            sigs[(lst, rhs, isdiag)] = 1
    print(f"{a.form}: {len(wins)} cells, {len(tb.ids)} types, {len(sigs)} constraints")
    pickle.dump(dict(Rt=a.Rt, keys=tb.keys, windows=tb.windows, sigs=list(sigs)), open(a.out, 'wb'))


if __name__ == '__main__':
    main()
